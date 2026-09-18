"""Supplemental fresh random space-time sample and independent FD at its peak."""
from pathlib import Path
import argparse,hashlib,json
import numpy as np
import aligned_continuation
from spacetime import Family,force
from validate import cartesian_residual
from pressure_morph import atomic_json

def compute(candidate,out,seed=9175194):
    f,raw=Family.load(candidate);rng=np.random.default_rng(seed)
    x=rng.uniform(-2,2,(4096,3));t=rng.uniform(.25,.75,4096)
    res=np.concatenate([f.analytic_residual(raw,x[i:i+256],t[i:i+256]) for i in range(0,len(x),256)])
    norms=np.linalg.norm(res,axis=1);k=int(norms.argmax())
    field=lambda x,t:f.fields(raw,x,t);forcing=lambda x,t:force(x,t,*raw[-2:])
    fd,div=cartesian_residual(field,forcing,x[k:k+1],float(t[k]),.00125,.000625)
    report=dict(candidate_sha256=hashlib.sha256(Path(candidate).read_bytes()).hexdigest(),seed=seed,points=len(x),sampled_spacetime_max=float(norms[k]),spacetime_RMS=float(np.sqrt(np.mean(norms**2))),point=x[k].tolist(),time=float(t[k]),residual_vector=res[k].tolist(),independent_fd_vector=fd[0].tolist(),fd_difference=float(np.linalg.norm(fd[0]-res[k])),scope='Fresh random space AND time, analytic residual cross-checked with independent Cartesian FD at sampled peak. Not a continuum maximum, not fixed-time volume L2, not a replacement acceptance gate.',pde_validated=False)
    atomic_json(out,report);print(json.dumps(report),flush=True);return report
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('candidate');p.add_argument('--out',required=True);a=p.parse_args();compute(a.candidate,a.out)
