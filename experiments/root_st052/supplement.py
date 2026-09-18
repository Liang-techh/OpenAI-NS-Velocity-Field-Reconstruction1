"""Post-freeze random space-time and fixed-peak FD diagnostics, no fitting."""
from pathlib import Path
import argparse,json,hashlib
import numpy as np
import minimax_exchange
from spacetime import Family,force,TIMES
from validate import cartesian_residual
from minimax_exchange import atomic_json

def compute(candidate,out):
    f,raw=Family.load(candidate);rng=np.random.default_rng(9175294)
    x=rng.uniform(-2,2,(4096,3));t=rng.uniform(.25,.75,4096)
    R=np.concatenate([f.analytic_residual(raw,x[k:k+256],t[k:k+256]) for k in range(0,len(x),256)])
    sizes=np.linalg.norm(R,axis=1);k=int(sizes.argmax());rows=[]
    field=lambda p,t:f.fields(raw,p,t);forcing=lambda p,t:force(p,t,*raw[-2:])
    fd,_=cartesian_residual(field,forcing,x[k:k+1],float(t[k]),.00125,.000625)
    result=dict(candidate_sha256=hashlib.sha256(Path(candidate).read_bytes()).hexdigest(),random_spacetime=dict(seed=9175294,points=len(x),maximum=float(sizes[k]),RMS=float(np.sqrt(np.mean(sizes**2))),point=x[k].tolist(),time=float(t[k]),independent_FD_difference=float(np.linalg.norm(fd[0]-R[k]))),fixed_peak_refinements=[])
    for seed in (9175291,9175292):
        x=np.random.default_rng(seed).uniform(-2,2,(4096,3))
        for t in TIMES:
            r=np.concatenate([f.analytic_residual(raw,x[k:k+256],t) for k in range(0,len(x),256)])
            k=int(np.linalg.norm(r,axis=1).argmax());pt=x[k:k+1];levels=[]
            for h in (.005,.0025,.00125):
                v,div=cartesian_residual(field,forcing,pt,float(t),h,.0025)
                levels.append(dict(space_step=h,time_step=.0025,norm=float(np.linalg.norm(v[0])),divergence=float(div[0])))
            result['fixed_peak_refinements'].append(dict(seed=seed,time=float(t),point=pt[0].tolist(),analytic_norm=float(np.linalg.norm(r[k])),levels=levels))
    result['scope']='Random-space-time RMS is NOT fixed-time volume L2; fixed-peak refinement is NOT a continuous maximum search. Frozen fields unchanged.'
    atomic_json(out,result);return result
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('candidate');p.add_argument('--out',required=True);a=p.parse_args();compute(a.candidate,a.out)
