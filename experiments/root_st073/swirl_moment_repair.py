"""Interface-preserving swirl moment fit, with disjoint time/axial holdouts."""
import json
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import lsq_linear
from joined_field import ROOT, JoinedField, independent_fd, coordinates

class SwirlRepair(JoinedField):
 def __init__(self, amplitudes=(0.,0.)):
  super().__init__(); self.amplitudes=np.asarray(amplitudes,float)
 def fields(self, points, tau):
  pts=np.asarray(points,float); u,p=super().fields(pts,tau)
  r=np.hypot(pts[:,0],pts[:,1]); sn=np.sqrt(self.nu)
  co=coordinates(r/sn,pts[:,2]/sn,np.broadcast_to(tau,r.shape),self.inner.h)
  q=np.asarray(co['q']); eta=np.asarray(co['eta']); ri=np.sqrt(2*self.nu*q*3/64)
  y=np.clip((r-ri)/ri,0,1)
  delta=sn*q**(-self.inner.A)*64*y**3*(1-y)**3*(self.amplitudes[0]+self.amplitudes[1]*eta**2)
  safe=np.where(r>0,r,1)
  u[:,0]-=delta*pts[:,1]/safe; u[:,1]+=delta*pts[:,0]/safe
  return u,p

def measure(f,k,eta,order=8,step=.001):
 tau=.5*2**(-k); g,w=leggauss(order)
 inner=f.inner.from_similarity([3/64],[eta],tau)[0]; ri=inner[0]; r=ri*(1+(g+1)/2)
 pts=np.column_stack((r,np.zeros(order),np.full(order,inner[2])))
 residual,div=independent_fd(f,pts,tau,step*np.sqrt(f.nu*tau),.00025*tau)
 return dict(k=k,eta=eta,order=order,step=step,
  terminal_rtheta=float(-np.sum(w*ri/2*r**2*residual[:,1])/(2*ri)**2),
  terminal_rz=float(-np.sum(w*ri/2*r*residual[:,2])/(2*ri)),
  momentum_max=float(np.max(np.linalg.norm(residual,axis=1))),
  component_max=np.max(np.abs(residual),axis=0).tolist(),
  divergence_max=float(np.max(np.abs(div))))

def run():
 training=[(3.,e) for e in (-.3,0,.3)]
 base=SwirlRepair(); units=[SwirlRepair((1,0)),SwirlRepair((0,1))]
 b=np.array([measure(base,*pt)['terminal_rtheta'] for pt in training])
 matrix=np.column_stack([np.array([measure(f,*pt)['terminal_rtheta'] for pt in training])-b for f in units])
 fit=lsq_linear(matrix,-b,bounds=(-4,4),tol=1e-12)
 corrected=SwirlRepair(fit.x)
 rows=[]
 for split,points in [('training',training),('holdout',[(.7,-.2),(.7,.2),(5.5,-.2),(5.5,.2)])]:
  for pt in points:
   before=measure(base,*pt,order=12); after=measure(corrected,*pt,order=12)
   rows.append(dict(split=split,before=before,after=after))
   print(json.dumps(rows[-1]),flush=True)
 refined=measure(corrected,5.5,.2,order=12,step=.0005)
 report=dict(amplitudes=fit.x.tolist(),bounds=[-4,4],training_matrix=matrix.tolist(),training_baseline=b.tolist(),linear_fit_defect=(matrix@fit.x+b).tolist(),rows=rows,refined_holdout=refined,
  correction='delta u_theta=sqrt(nu) q^(-A) 64 y^3(1-y)^3 (a0+a1 eta^2), y=(r-ri)/ri, supported on 0<y<1',
  invariants='Axisymmetric pure swirl is solenoidal. Cubic endpoint zeros preserve velocity and first/second spatial jets at both moving interfaces. Core and pressure unchanged.',
  limitations='Unforced full momentum finite differences. Fit only angular radial moment, not full pointwise PDE. Axial slab only, no global finite energy or continuum certificate.',pde_validated=False,global_field_ready=False)
 out=ROOT/'swirl_moment_repair';out.mkdir(exist_ok=True)
 (out/'report.json').write_bytes((json.dumps(report,indent=2)+'\n').encode())
 print('AMPLITUDES',fit.x,flush=True)
if __name__=='__main__':run()
