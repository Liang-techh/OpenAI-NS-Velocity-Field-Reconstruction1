"""Deterministic diagnostic max search, fresh times and separate FD refinement."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import acceleration_fit
from spacetime import Family,TIMES,force
from validate import cartesian_residual

def check(path):
    f,raw=Family.load(path);times=TIMES;rows=[]
    # Fixed before any resulting evaluations, not used in optimizing coefficients.
    r=np.linspace(0,1.995,81);z=np.linspace(-1.995,1.995,121);R,Z=np.meshgrid(r,z,indexing='ij');x=np.c_[R.ravel(),np.zeros(R.size),Z.ravel()]
    for t in times:
        vals=np.concatenate([f.analytic_residual(raw,x[i:i+256],t) for i in range(0,len(x),256)]);n=np.linalg.norm(vals,axis=1);i=int(np.argmax(n));p=x[i:i+1]
        fd,div=cartesian_residual(lambda x,t:f.fields(raw,x,t),lambda x,t:force(x,t,*raw[-2:]),p,t,.00125,.00125)
        rows.append(dict(time=float(t),analytic_grid_max=float(n[i]),point=p[0].tolist(),FD_at_grid_peak=float(np.linalg.norm(fd[0])),analytic_FD_difference=float(np.linalg.norm(fd[0]-vals[i]))))
    rng=np.random.default_rng(9174692);x=rng.uniform(-2,2,(1024,3));newtimes=np.sort(rng.uniform(.25,.75,11));temporal=[]
    for t in newtimes:
        vals=np.concatenate([f.analytic_residual(raw,x[i:i+256],t) for i in range(0,len(x),256)]);norm=np.linalg.norm(vals,axis=1)
        temporal.append(dict(time=float(t),max=float(norm.max()),volume_L2=float(np.sqrt(64*np.mean(norm*norm)))))
    return dict(candidate_sha256=hashlib.sha256(Path(path).read_bytes()).hexdigest(),deterministic_grid_shape=[81,121],dense_grid=rows,additional_times=temporal,scope='Supplementary diagnostic samples, not rigorous maximum/interval bounds. Analytic operator with independent Cartesian FD at each grid peak; no fitting on these samples.',pde_validated=False)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('candidate');p.add_argument('--out',required=True);a=p.parse_args();d=check(a.candidate);Path(a.out).write_text(json.dumps(d,indent=2)+'\n');print(json.dumps({'dense_max':max(x['analytic_grid_max'] for x in d['dense_grid']),'fresh_time_max':max(x['max'] for x in d['additional_times'])},indent=2))
