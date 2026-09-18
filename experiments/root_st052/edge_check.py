"""Independent finer axial/radial edge grid after freezing; no optimizer calls."""
from pathlib import Path
import argparse,hashlib
import numpy as np
import minimax_exchange
from spacetime import Family,force
from validate import cartesian_residual
from minimax_exchange import atomic_json

def compute(candidate,out):
    f,raw=Family.load(candidate);times=np.linspace(.25,.75,13);rows=[]
    grids=[('axial_edges',np.linspace(.012,1.4,41),np.r_[-np.linspace(1.7,1.997,51),np.linspace(1.7,1.997,51)]),('radial_edge',np.linspace(1.7,1.997,31),np.linspace(-1.997,1.997,31))]
    field=lambda x,t:f.fields(raw,x,t);forcing=lambda x,t:force(x,t,*raw[-2:])
    for name,r,z in grids:
        R,Z=np.meshgrid(r,z,indexing='ij');x=np.c_[R.ravel(),np.zeros(R.size),Z.ravel()]
        for t in times:
            rr=np.concatenate([f.analytic_residual(raw,x[i:i+384],t) for i in range(0,len(x),384)])
            norm=np.linalg.norm(rr,axis=1);k=int(norm.argmax());fd,_=cartesian_residual(field,forcing,x[k:k+1],float(t),.00125,.000625)
            rows.append(dict(grid=name,time=float(t),points=len(x),maximum=float(norm[k]),point=x[k].tolist(),FD_difference=float(np.linalg.norm(fd[0]-rr[k]))))
    result=dict(candidate_sha256=hashlib.sha256(Path(candidate).read_bytes()).hexdigest(),rows=rows,scope='Independent axial/radial edge grids at13times after freezing. Not domain-complete, not a continuum supremum or a replacement for original max/L2 gates.')
    atomic_json(out,result);print('EDGE',out,max(r['maximum'] for r in rows),flush=True);return result
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('candidate');p.add_argument('--out',required=True);a=p.parse_args();compute(a.candidate,a.out)
