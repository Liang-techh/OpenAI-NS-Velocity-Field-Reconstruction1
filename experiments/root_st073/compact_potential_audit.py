"""Check plateau matching, collar momentum, and finite-volume kinetic energy."""
import json
import numpy as np
from numpy.polynomial.legendre import leggauss
from compact_potential import CompactPotentialField
from joined_field import ROOT,independent_fd

def run():
 f=CompactPotentialField();rows=[]
 for k in (.4,5.5):
  tau=.5*2**(-k);q=tau;sn=np.sqrt(f.nu);ri=sn*np.sqrt(2*q*3/64);rflat,rsupp,zflat,zsupp=f.support(tau)
  plateau=np.array([[ri*.5,0,0],[ri*3,0,0],[ri*5,0,0],[ri*3,0,zflat*.5]])
  u,p=f.fields(plateau,tau);ub,pb=f.field.fields(plateau,tau)
  # Physical cylindrical Gauss quadrature; finite only, no convergence claim.
  g,w=leggauss(6);rr=(g+1)*rsupp/2;zz=g*zsupp;pts=np.array([[rad,0,z] for z in zz for rad in rr]);v,_=f.fields(pts,tau)
  weight=np.repeat(w*zsupp,6)*np.tile(w*rsupp/2,6)*2*np.pi*np.tile(rr,6)
  energy=float(.5*np.sum(weight*np.sum(v*v,axis=1)))
  zcollar=(zflat+zsupp)/2;rcollar=(rflat+rsupp)/2
  sample=np.array([[ri*3,0,zcollar],[ri*.5,0,zcollar],[rcollar,0,0],[rcollar,0,zcollar]])
  residual,div=independent_fd(f,sample,tau,.0005*np.sqrt(f.nu*tau),.0001*tau)
  row=dict(k=k,tau=tau,support=dict(rflat=float(rflat),rsupp=float(rsupp),zflat=float(zflat),zsupp=float(zsupp)),plateau_velocity_max_difference=float(np.max(np.abs(u-ub))),plateau_pressure_max_difference=float(np.max(np.abs(p-pb))),sampled_kinetic_energy=energy,collar_points=sample.tolist(),collar_momentum_norms=np.linalg.norm(residual,axis=1).tolist(),collar_divergence=div.tolist())
  rows.append(row);print(json.dumps(row),flush=True)
 report=dict(rows=rows,construction='V=curl(chi*A), A_theta=psi/r, A_z=-integral_0^r u_theta ds; chi compactly supported in registered eta slab. Therefore div V=0 identically and V has finite spatial energy for each registered time if smooth.',
 limitations='Only tau in[.5/64,.5]; no continuation to critical tau=0. Kinetic energy uses6x6 axisymmetric Gauss quadrature, not a uniform bound or convergence certificate. Collar residuals are four samples, not global max or L2. C2 radial joins; cutoff-induced momentum defect expected. The pressure is chi times background pressure.',pde_validated=False,global_field_ready=False)
 out=ROOT/'compact_potential';out.mkdir(exist_ok=True);(out/'report.json').write_bytes((json.dumps(report,indent=2)+'\n').encode())
if __name__=='__main__':run()
