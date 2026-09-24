"""Multi-mode angular PDE collocation; no full momentum admission implied."""
import json
import numpy as np
from numpy.polynomial.legendre import leggauss, legval
from scipy.optimize import lsq_linear
from joined_field import ROOT, JoinedField, coordinates, independent_fd

class CachedBase(JoinedField):
 def __init__(self):super().__init__();self.cache={}
 def fields(self,points,tau):
  pts=np.asarray(points,float);ts=np.broadcast_to(tau,(len(pts),))
  key=(pts.shape,pts.tobytes(),ts.tobytes())
  if key not in self.cache:self.cache[key]=super().fields(pts,ts)
  u,p=self.cache[key];return u.copy(),p.copy()

class AngularModes:
 def __init__(self,base,amplitudes):
  self.base=base;self.inner=base.inner;self.nu=base.nu;self.amplitudes=np.asarray(amplitudes).reshape(4,3)
 def fields(self,points,tau):
  pts=np.asarray(points,float);u,p=self.base.fields(pts,tau)
  r=np.hypot(pts[:,0],pts[:,1]);sn=np.sqrt(self.nu)
  co=coordinates(r/sn,pts[:,2]/sn,np.broadcast_to(tau,r.shape),self.inner.h)
  q=np.asarray(co['q']);eta=np.asarray(co['eta']);ri=np.sqrt(2*self.nu*q*3/64)
  y=np.clip(r/ri-1,0,1);amp=np.zeros_like(r)
  for j in range(4):
   radial=legval(2*y-1,[0]*j+[1])
   for m in range(3):amp+=self.amplitudes[j,m]*radial*(eta/.3)**m
  delta=sn*q**(-self.inner.A)*64*y**3*(1-y)**3*amp
  safe=np.where(r>0,r,1);u[:,0]-=delta*pts[:,1]/safe;u[:,1]+=delta*pts[:,0]/safe
  return u,p

def sample(f,k,eta,order,step=.001):
 tau=.5*2**(-k);g,w=leggauss(order);ip=f.inner.from_similarity([3/64],[eta],tau)[0]
 ri=ip[0];r=ri*(1+(g+1)/2);pts=np.column_stack((r,np.zeros(order),np.full(order,ip[2])))
 res,div=independent_fd(f,pts,tau,step*np.sqrt(f.nu*tau),.00025*tau)
 moment_weights=w*ri/2*r*r/(2*ri)**2
 return res,div,moment_weights

def metrics(res,div,mw):
 return dict(momentum_max=float(np.max(np.linalg.norm(res,axis=1))),theta_max=float(np.max(np.abs(res[:,1]))),terminal_rtheta=float(-mw@res[:,1]),divergence_max=float(np.max(np.abs(div))))

def run():
 base=CachedBase();zero=AngularModes(base,np.zeros(12));train=[(k,e) for k in (1.,4.) for e in (-.3,0,.3)]
 matrices=[];targets=[]
 for k,e in train:
  rb,db,mw=sample(zero,k,e,12);scale=np.max(np.abs(rb[:,1]));cols=[]
  for j in range(12):
   unit=np.zeros(12);unit[j]=1
   r,_,_=sample(AngularModes(base,unit),k,e,12);cols.append(r[:,1]-rb[:,1])
  matrix=np.column_stack(cols);average=mw/np.sum(mw)
  matrices.append(np.vstack((matrix/scale,3*(average@matrix)[None,:]/scale)))
  targets.append(np.r_[-rb[:,1]/scale,-3*average@rb[:,1]/scale])
  print('assembled',k,e,flush=True)
 M=np.vstack(matrices);b=np.concatenate(targets)
 fit=lsq_linear(np.vstack((M,.001*np.eye(12))),np.r_[b,np.zeros(12)],bounds=(-4,4),tol=1e-10)
 corrected=AngularModes(base,fit.x);rows=[]
 for split,points,order in [('training',train,12),('holdout',[(k,e) for k in (2.5,5.5) for e in (-.2,.2)],18)]:
  for k,e in points:
   rb,db,mw=sample(zero,k,e,order);ra,da,_=sample(corrected,k,e,order)
   row=dict(split=split,k=k,eta=e,order=order,before=metrics(rb,db,mw),after=metrics(ra,da,mw));rows.append(row);print(json.dumps(row),flush=True)
 rr,dd,mw=sample(corrected,5.5,.2,18,.0005)
 report=dict(amplitudes=fit.x.reshape(4,3).tolist(),bounds=[-4,4],rows=rows,refined_holdout=metrics(rr,dd,mw),fit_cost=float(fit.cost),
 basis='sqrt(nu) q^-A 64 y^3(1-y)^3 P_j(2y-1) (eta/.3)^m; j=0..3,m=0..2',
 objective='Per slice normalize theta pointwise residual by baseline theta max; append 3 times r^2-weighted average theta residual; coefficient regularization .001. Soft moment penalty, not exact constraint.',
 scope='Same fixed core, poloidal field, pressure and heat exterior. Pure axisymmetric swirl preserves divergence and interface jets. No force. Full momentum evaluated on disjoint scale/axial holdouts; samples do not certify supremum or physical-volume L2.',pde_validated=False,global_field_ready=False)
 out=ROOT/'angular_collocation';out.mkdir(exist_ok=True);(out/'report.json').write_bytes((json.dumps(report,indent=2)+'\n').encode())
if __name__=='__main__':run()
