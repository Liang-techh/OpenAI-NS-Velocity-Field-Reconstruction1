"""Controlled radial swirl loops: shear gain versus full viscous residual cost."""
import json
import numpy as np
from numpy.polynomial.legendre import leggauss,legfit,legint,legval,legvander
from wide_modes import CachedWidth,WideJointModes
from joined_field import ROOT,coordinates
from affine_momentum import jets,momentum

class SwirlLoop:
 def __init__(self,base,n,amplitude):
  self.base=base;self.inner=base.inner;self.nu=base.nu;self.ratio=base.ratio;self.n=n;self.amplitude=amplitude
 def fields(self,points,tau):
  pts=np.asarray(points,float);u,p=self.base.fields(pts,tau);r=np.hypot(pts[:,0],pts[:,1]);co=coordinates(r/np.sqrt(self.nu),pts[:,2]/np.sqrt(self.nu),np.broadcast_to(tau,r.shape),self.inner.h)
  ri=np.sqrt(2*self.nu*co['q']*3/64);raw=(r/ri-1)/5;y=np.clip(raw,0,1);active=(raw>0)&(raw<1)
  bubble=256*y**4*(1-y)**4;delta=np.where(active,np.sqrt(self.nu)*co['q']**(-self.inner.A)*self.amplitude*bubble*np.sin(2*np.pi*self.n*y),0)
  safe=np.where(r>0,r,1);u[:,0]-=delta*pts[:,1]/safe;u[:,1]+=delta*pts[:,0]/safe
  return u,p

def evaluate(f,k,eta,g,Q):
 tau=.5*2**(-k);ip=f.inner.from_similarity([3/64],[eta],tau)[0];ri=ip[0];width=5*ri;r=ri+width*(g+1)/2;pts=np.column_stack((r,np.zeros(len(r)),np.full(len(r),ip[2])))
 u,J,L=jets(f,pts,tau,.0005*np.sqrt(f.nu*tau),.00025*tau);R=momentum((u,J,L))
 T=np.column_stack((-width/2*(Q@(r*r*R[:,1]))/(r*r),-width/2*(Q@(r*R[:,2]))/r))
 F=u[:,1]/r;shear=np.column_stack((J[:,1,0]-F,J[:,2,0]));sn=np.linalg.norm(shear,axis=1);N=shear/np.maximum(sn[:,None],1e-15);lam2=-2*F*N[:,0]*(2*F*N[:,0]+sn)
 lam=np.sqrt(np.maximum(lam2,0));tn=np.sum(T*N,axis=1);tk=T[:,1]*N[:,0]-T[:,0]*N[:,1];q=2*F*N[:,0]
 direction=(lam2>0)&(tn<0)&(lam*np.abs(tk)<np.abs(q)*(-tn));growth=lam>f.nu/r**2;middle=((g+1)/2>=.12)&((g+1)/2<=.88)
 return dict(max_momentum=float(np.max(np.linalg.norm(R,axis=1))),rms_momentum=float(np.sqrt(np.mean(np.sum(R*R,axis=1)))),middle_direction=int(np.count_nonzero(direction&middle)),middle_joint=int(np.count_nonzero(direction&growth&middle)),all_joint=int(np.count_nonzero(direction&growth)),middle_joint_y=((g+1)/2)[direction&growth&middle].tolist(),tangential_moment=float(T[-1,0]),axial_moment=float(T[-1,1]))

def run():
 a=json.loads((ROOT/'wide_collocation/report.json').read_text())['amplitudes'];base=WideJointModes(CachedWidth(6),a);g,w=leggauss(16);v=legvander(g,15);A=np.column_stack([legval(g,legint([0]*j+[1]))-legval(-1,legint([0]*j+[1])) for j in range(16)]);Q=A@np.linalg.inv(v)
 configs=[(0,0)]+[(n,amp) for n in (2,4) for amp in (-.3,-.1,.1,.3)];rows=[]
 for n,amp in configs:
  f=base if n==0 else SwirlLoop(base,n,amp)
  for eta in (-.3,0,.3):
   row=dict(radial_cycles=n,amplitude=amp,eta=eta,k=5.5,**evaluate(f,5.5,eta,g,Q));rows.append(row)
  print('config',n,amp,'joint passes',[row['middle_joint'] for row in rows[-3:]],'max',[row['max_momentum'] for row in rows[-3:]],flush=True)
 report=dict(rows=rows,construction='Axisymmetric pure swirl sqrt(nu) q^-A amplitude 256 y^4(1-y)^4 sin(2pi*n*y), y=(r-ri)/(5ri), n=2 or4. Exact divergence-free, preserves value and spatial jets through order2 at both radial interfaces.',
 scope='Finite late-time radial/axial slice screen. Actual full Cartesian FD momentum includes viscous cost. Approximate physical frozen direction/growth conditions, not the paper normalized theorem. No physical-volume L2, axial closure or globally supported wave.',pde_validated=False,global_field_ready=False)
 out=ROOT/'radial_shear_loop';out.mkdir(exist_ok=True);(out/'report.json').write_bytes((json.dumps(report,indent=2)+'\n').encode())
if __name__=='__main__':run()
