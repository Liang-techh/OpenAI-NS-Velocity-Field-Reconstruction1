"""Coupled solenoidal streamfunction / pressure annular collocation."""
import json
import numpy as np
from numpy.polynomial import Polynomial as Poly
from numpy.polynomial.legendre import Legendre,leggauss
from scipy.optimize import least_squares
from angular_collocation import CachedBase,sample,metrics
from joined_field import ROOT,coordinates
from affine_momentum import jets,momentum,combine

class PoloidalModes:
 def __init__(self,base,amplitudes):
  self.base=base;self.inner=base.inner;self.nu=base.nu;self.a=np.asarray(amplitudes).reshape(2,3,2)
 def fields(self,points,tau):
  pts=np.asarray(points,float);u,p=self.base.fields(pts,tau)
  r=np.hypot(pts[:,0],pts[:,1]);sn=np.sqrt(self.nu);safe=np.where(r>0,r,1)
  co=coordinates(r/sn,pts[:,2]/sn,np.broadcast_to(tau,r.shape),self.inner.h)
  q=np.asarray(co['q']);eta=np.asarray(co['eta']);qz=np.asarray(co['q_z'])/sn;ez=np.asarray(co['eta_z'])/sn
  ri=np.sqrt(2*self.nu*q*3/64);raw_y=r/ri-1;active=(raw_y>0)&(raw_y<1);y=np.clip(raw_y,0,1);yz=-(1+y)*qz/(2*q)
  psi_scale=self.nu**1.5*(3/32)*q**(1-self.inner.A)
  psi_r=np.zeros_like(r);psi_z=np.zeros_like(r);dp=np.zeros_like(r)
  for j in range(3):
   leg=Legendre.basis(j).convert(kind=Poly)(Poly([-1,2]))
   bubble=256*Poly([0,1])**4*Poly([1,-1])**4*leg
   pressure=64*Poly([0,1])**3*Poly([1,-1])**3*leg
   B=bubble(y);By=bubble.deriv()(y)
   for m in range(2):
    axial=(eta/.3)**m;axial_z=np.zeros_like(eta) if m==0 else ez/.3
    psi_r+=self.a[0,j,m]*psi_scale*By*axial/ri
    psi_z+=self.a[0,j,m]*psi_scale*((1-self.inner.A)*qz/q*B*axial+By*yz*axial+B*axial_z)
    dp+=self.a[1,j,m]*self.nu*q**(-2*self.inner.A)*pressure(y)*axial
  ur=np.where(active,-psi_z/safe,0);uz=np.where(active,psi_r/safe,0);dp=np.where(active,dp,0)
  u[:,0]+=ur*pts[:,0]/safe;u[:,1]+=ur*pts[:,1]/safe;u[:,2]+=uz
  return u,p+dp

def run():
 base=CachedBase();zero=PoloidalModes(base,np.zeros(12));train=[(k,e) for k in (1.,4.) for e in (-.3,0,.3)]
 data=[]
 for k,e in train:
  tau=.5*2**(-k);g,w=leggauss(12);ip=base.inner.from_similarity([3/64],[e],tau)[0]
  pts=np.column_stack((ip[0]*(1+(g+1)/2),np.zeros(12),np.full(12,ip[2])))
  args=(pts,tau,.001*np.sqrt(base.nu*tau),.00025*tau)
  jb=jets(zero,*args);mode=[]
  for j in range(12):
   a=np.zeros(12);a[j]=1;ju=jets(PoloidalModes(base,a),*args)
   mode.append(tuple(x-y for x,y in zip(ju,jb)))
  modes=tuple(np.stack([v[i] for v in mode]) for i in range(3));scale=np.max(np.linalg.norm(momentum(jb),axis=1))
  data.append((jb,modes,scale));print('assembled',k,e,flush=True)
 def objective(a):
  return np.r_[np.concatenate([(momentum(combine(b,m,a))/s).ravel() for b,m,s in data]),.001*a]
 fit=least_squares(objective,np.zeros(12),bounds=(-4,4),max_nfev=120,ftol=1e-10,xtol=1e-10,gtol=1e-10)
 corrected=PoloidalModes(base,fit.x);rows=[]
 for split,points,order in [('training',train,12),('holdout',[(k,e) for k in (2.5,5.5) for e in (-.2,.2)],18)]:
  for k,e in points:
   rb,db,mw=sample(zero,k,e,order);ra,da,_=sample(corrected,k,e,order)
   row=dict(split=split,k=k,eta=e,order=order,before=metrics(rb,db,mw),after=metrics(ra,da,mw));rows.append(row);print(json.dumps(row),flush=True)
 rr,dd,mw=sample(corrected,5.5,.2,18,.0005)
 # Check affine-jet prediction against a fresh complete-field FD evaluation.
 rb,db,_=sample(corrected,4.,.3,12)
 predicted=momentum(combine(data[-1][0],data[-1][1],fit.x))
 report=dict(amplitudes=fit.x.reshape(2,3,2).tolist(),bounds=[-4,4],fit_nfev=fit.nfev,fit_success=bool(fit.success),fit_cost=float(fit.cost),rows=rows,refined_holdout=metrics(rr,dd,mw),surrogate_full_fd_max_difference=float(np.max(np.abs(rb-predicted))),
 streamfunction='nu^1.5 (3/32) q^(1-A) 256 y^4(1-y)^4 P_j(2y-1)(eta/.3)^m, j=0..2,m=0..1; ur=-psi_z/r, uz=psi_r/r, physical z derivatives include moving ri',
 pressure='nu q^(-2A) 64 y^3(1-y)^3 P_j(2y-1)(eta/.3)^m',
 scope='Autonomous compact bubble basis. Full nonlinear unforced momentum objective, each training slice normalized by its baseline full maximum. No global axial closure, finite-energy proof, physical-volume L2 or supremum certificate. Frozen core and original swirl retained.',pde_validated=False,global_field_ready=False)
 out=ROOT/'poloidal_collocation';out.mkdir(exist_ok=True);(out/'report.json').write_bytes((json.dumps(report,indent=2)+'\n').encode())
if __name__=='__main__':run()
