"""Radial inverse of actual full tangential residual; not full tensor realization."""
import json,numpy as np
from numpy.polynomial.legendre import leggauss
from joined_field import ROOT,JoinedField,independent_fd

def run():
 f=JoinedField();rows=[]
 for order in (8,12):
  g,w=leggauss(order)
  for k in (.4,5.5):
   tau=.5*2**(-k)
   for eta in (-.3,0,.3):
    inner=f.inner.from_similarity([3/64],[eta],tau)[0];ri=inner[0];ro=2*ri;r=ri+(g+1)*ri/2
    points=np.column_stack((r,np.zeros(order),np.full(order,inner[2])))
    residual,div=independent_fd(f,points,tau,.001*np.sqrt(f.nu*tau),.00025*tau)
    theta_moment=float(np.sum(w*ri/2*r*r*residual[:,1]));axial_moment=float(np.sum(w*ri/2*r*residual[:,2]))
    terminal=np.array([-theta_moment/ro**2,-axial_moment/ro])
    rows.append(dict(order=order,k=k,eta=eta,inner_radius=ri,outer_radius=ro,theta_weighted_moment=theta_moment,axial_weighted_moment=axial_moment,terminal_rtheta_stress=float(terminal[0]),terminal_rz_stress=float(terminal[1]),psd_trace_lower_bound=float(2*np.linalg.norm(terminal)),sampled_momentum_max=float(np.max(np.linalg.norm(residual,axis=1))),sampled_divergence_max=float(np.max(np.abs(div)))))
 report=dict(rows=rows,definition='T_rtheta(r)=-r^-2 integral_ri^r s² Rtheta(s) ds; T_rz(r)=-r^-1 integral_ri^r s Rz(s) ds at fixed physical z,t, with zero inner stress.',limitations='Tangential radial inverse only, assuming no theta-z or z-z stress contribution. A symmetric r-z stress also affects radial momentum through z derivative. PSD completion adds diagonal stresses and is NOT a dynamically realizable wave. Terminal zero requires these two moments vanish for this restricted inverse; not an impossibility theorem for all stress tensors.',forcing='No external force fitted or added',global_field_ready=False,pde_validated=False)
 out=ROOT/'transition_stress';out.mkdir(exist_ok=True);(out/'moments.json').write_bytes((json.dumps(report,indent=2)+'\n').encode());print(json.dumps(rows[-3:],indent=2))
if __name__=='__main__':run()
