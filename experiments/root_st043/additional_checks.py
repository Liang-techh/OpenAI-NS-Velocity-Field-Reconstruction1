"""Additional fixed-field precision and a pressure-direction necessary condition."""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import bootstrap
from spacetime import Family,force,TIMES
from validate import cartesian_residual,norms
from structure_audit import jets

def run(path,out,seed=9174301):
 f,raw=Family.load(path);field=lambda x,t:f.fields(raw,x,t);forcing=lambda x,t:force(x,t,*raw[-2:]);rng=np.random.default_rng(seed);x=rng.uniform(-2,2,(4096,3));rows=[]
 for h in [.0025,.00125]:
  for t in TIMES:
   R,div=cartesian_residual(field,forcing,x,float(t),h,.0025);idx=np.argmax(np.linalg.norm(R,axis=1));rows.append(dict(time=float(t),space_step=h,time_step=.0025,**norms(R,div),maximum_at=x[idx].tolist()))
  print('completed additional fixed-sample h',h,flush=True)
 # If sign(z)*p_z >=0, then |R_z| >= max(sign(z)*M_z,0),
 # M_z=u_t,z+(u.grad)u_z-nu*lap(u_z)-f_z; this bound fixes u and f.
 R,Z=np.meshgrid(np.linspace(.04,.2,11),np.r_[-np.linspace(.04,.2,11)[::-1],np.linspace(.04,.2,11)],indexing='ij');bounds=[]
 for t in TIMES:
  tau=1-t;r=R.ravel()*np.sqrt(tau);z=Z.ravel()*tau**.495;p=np.c_[r,np.zeros(len(r)),z]
  j=jets(f,raw,r,z,float(t));res=f.analytic_residual(raw,p,float(t));M=res[:,2]-j[:,7];lower=np.maximum(np.sign(z)*M,0);idx=int(np.argmax(lower))
  # Independent FD crosscheck at this same maximizing point.
  fd,_=cartesian_residual(field,forcing,p[idx:idx+1],float(t),.00125,.000625)
  e=np.array([0,0,1e-5]);pg=(field(p[idx:idx+1]+e,float(t))[1]-field(p[idx:idx+1]-e,float(t))[1])/(2e-5)
  fdM=float(fd[0,2]-pg[0]);bounds.append(dict(time=float(t),sampled_pressure_direction_lower_bound=float(lower[idx]),point=p[idx].tolist(),frozen_axial_pressure_gradient=float(j[idx,7]),frozen_axial_residual=float(res[idx,2]),analytic_material_defect=float(M[idx]),independent_FD_material_defect=fdM,derivative_discrepancy=abs(fdM-float(M[idx]))))
 d=dict(candidate_sha256=hashlib.sha256(Path(path).read_bytes()).hexdigest(),validation_seed=seed,uniform_points=4096,spatial_refinement=rows,pressure_direction_necessary_condition=bounds,scope='Additional derivative resolution does not overwrite original gate results. Inequality is pointwise for fixed velocity/force; numerical lower-bound values are floating estimates at stated points, not interval certificates or a family-wide no-go.',pde_validated=False)
 Path(out).write_text(json.dumps(d,indent=2)+'\n')

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('candidate');p.add_argument('--out',required=True);p.add_argument('--seed',type=int,default=9174301);a=p.parse_args();run(a.candidate,a.out,a.seed)
