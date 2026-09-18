"""Supplementary fixed-field cylinder-grid maxima with independent FD checks."""
import argparse,hashlib,json,sys
from pathlib import Path
import numpy as np
import continuation
from spacetime import Family,force,TIMES
from validate import cartesian_residual

def run(candidate,out):
 f,raw=Family.load(candidate);r,z=np.meshgrid(np.linspace(0,1.995,81),np.linspace(-1.995,1.995,121),indexing='ij')
 x=np.c_[r.ravel(),np.zeros(r.size),z.ravel()];rows=[]
 field=lambda p,t:f.fields(raw,p,t);forcing=lambda p,t:force(p,t,*raw[-2:])
 for t in TIMES:
  vals=np.concatenate([f.analytic_residual(raw,x[k:k+256],t) for k in range(0,len(x),256)])
  n=np.linalg.norm(vals,axis=1);i=int(np.argmax(n));fd,div=cartesian_residual(field,forcing,x[i:i+1],float(t),.00125,.000625)
  rows.append(dict(time=float(t),sampled_cylinder_grid_max=float(n[i]),point=x[i].tolist(),component_vector=vals[i].tolist(),independent_FD_vector=fd[0].tolist(),FD_vector_difference=float(np.linalg.norm(fd[0]-vals[i]))))
  print(json.dumps(rows[-1]),flush=True)
 result=dict(candidate_sha256=hashlib.sha256(Path(candidate).read_bytes()).hexdigest(),grid=[81,121],time_slices=len(TIMES),grid_includes_axis=True,rows=rows,scope='Supplementary sampled maxima, NOT a continuous supremum or an L2 quadrature. No fitting after observing this grid.',pde_validated=False)
 Path(out).write_text(json.dumps(result,indent=2));return result
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('candidate');p.add_argument('--out',required=True);a=p.parse_args();run(a.candidate,a.out)
