"""Supplementary frozen random-sample peak refinement; never changes coefficients."""
from __future__ import annotations
import argparse,hashlib
from pathlib import Path
import numpy as np
import aligned_continuation
from spacetime import Family,force,TIMES
from validate import cartesian_residual
from pressure_morph import atomic_json

def compute(candidate,out,seed):
    f,raw=Family.load(candidate);x=np.random.default_rng(seed).uniform(-2,2,(4096,3))
    field=lambda x,t:f.fields(raw,x,t);forcing=lambda x,t:force(x,t,*raw[-2:])
    rows=[]
    for t in TIMES:
        r=np.concatenate([f.analytic_residual(raw,x[i:i+256],t) for i in range(0,len(x),256)])
        norms=np.linalg.norm(r,axis=1);k=int(norms.argmax());pt=x[k:k+1];levels=[]
        for h in (.005,.0025,.00125):
            fd,div=cartesian_residual(field,forcing,pt,float(t),h,.0025)
            levels.append(dict(space_step=h,time_step=.0025,full_norm=float(np.linalg.norm(fd[0])),divergence=float(div[0]),analytic_fd_vector_difference=float(np.linalg.norm(r[k]-fd[0]))))
        rows.append(dict(time=float(t),point=x[k].tolist(),sampled_analytic_max=float(norms[k]),levels=levels))
    result=dict(candidate_sha256=hashlib.sha256(Path(candidate).read_bytes()).hexdigest(),seed=seed,rows=rows,scope='Same frozen 4096-point validation locations. Find analytic peak at each fixed time, then vary only spatial FD step there. Not a new domain-wide maximization or continuum bound; original reports remain unchanged.',pde_validated=False)
    atomic_json(out,result);print(out,max(row['levels'][-1]['full_norm'] for row in rows),flush=True);return result
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('candidate');p.add_argument('--out',required=True);p.add_argument('--seed',type=int,required=True);a=p.parse_args();compute(a.candidate,a.out,a.seed)
