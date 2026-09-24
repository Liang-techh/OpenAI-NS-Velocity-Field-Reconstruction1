"""Physical-volume audit of the two-cycle radial shear modulation."""
import json
import numpy as np
from numpy.polynomial.legendre import leggauss,legvander,legint,legval
from radial_shear_loop import SwirlLoop,evaluate
from wide_modes import CachedWidth,WideJointModes
from joined_field import ROOT,independent_fd

def run():
 a=json.loads((ROOT/'wide_collocation/report.json').read_text())['amplitudes'];base=WideJointModes(CachedWidth(6),a);mod=SwirlLoop(base,2,-.3);g,w=leggauss(16);v=legvander(g,15);A=np.column_stack([legval(g,legint([0]*j+[1]))-legval(-1,legint([0]*j+[1])) for j in range(16)]);Q=A@np.linalg.inv(v)
 slices=[]
 for eta in (-.3,0,.3):
  for label,f in [('prior',base),('modulated',mod)]:
   item=dict(eta=eta,k=5.5,field=label,**evaluate(f,5.5,eta,g,Q));slices.append(item);print(item,flush=True)
 volume=[]
 for k in (2.5,5.5):
  tau=.5*2**(-k);gy,wy=leggauss(12);ge,we=leggauss(8);eta=np.repeat(.4*ge,12);y=np.tile((gy+1)/2,8);q=tau/(1-eta**2)
  ri=np.sqrt(2*base.nu*q*3/64);r=ri*(1+5*y);z=np.sqrt(base.nu)*q**base.inner.D*eta;ze=np.sqrt(base.nu)*q**base.inner.D*(1+2*base.inner.D*eta**2/(1-eta**2));wv=np.repeat(.4*we,12)*np.tile(wy/2,8)*2*np.pi*r*5*ri*ze
  pts=np.column_stack((r,np.zeros_like(r),z))
  for label,f in [('prior',base),('modulated',mod)]:
   R,div=independent_fd(f,pts,tau,.0005*np.sqrt(base.nu*tau),.00025*tau)
   item=dict(k=k,field=label,volume=float(wv.sum()),momentum_max=float(np.max(np.linalg.norm(R,axis=1))),physical_volume_L2=float(np.sqrt(np.sum(wv*np.sum(R*R,axis=1)))),divergence_max=float(np.max(np.abs(div))))
   volume.append(item);print(item,flush=True)
 report=dict(slices=slices,volume=volume,scope='Unfitted two-cycle amplitude -.3 mean swirl loop on corrected width6 field, exact axisymmetric solenoidality, full unforced Cartesian FD residual. Approximate directional screen at finite nodes; no PDE acceptance or globally supported waves.',pde_validated=False,global_field_ready=False)
 out=ROOT/'radial_shear_loop';(out/'candidate_audit.json').write_bytes((json.dumps(report,indent=2)+'\n').encode())
if __name__=='__main__':run()
