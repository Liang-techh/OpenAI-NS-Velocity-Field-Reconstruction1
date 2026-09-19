"""Fresh deterministic spatial grid at intermediate times, plus pressure-only limits."""
from __future__ import annotations
import argparse,json,hashlib
from pathlib import Path
import numpy as np
from pressure_completion import Family,residual,save
from spacetime import force
from validate import cartesian_residual

def check(parent,child,out):
    f,b=Family.load(parent);fc,c=Family.load(child)
    rr=np.unique(np.r_[np.linspace(0,1.995,41),np.linspace(.025,.65,13)])
    zz=np.unique(np.r_[np.linspace(-1.995,1.995,81),np.linspace(-1.98,-1.7,17),np.linspace(1.7,1.98,17)])
    R,Z=np.meshgrid(rr,zz,indexing='ij');p=np.c_[R.ravel(),np.zeros(R.size),Z.ravel()]
    ts=np.linspace(.25,.75,13);rows=[]
    for t in ts:
        r0=residual(f,b,p,t);r1=residual(f,c,p,t)
        n0=np.linalg.norm(r0,axis=1);n1=np.linalg.norm(r1,axis=1);i0=int(n0.argmax());i1=int(n1.argmax());it=int(abs(r0[:,1]).argmax())
        fd,_=cartesian_residual(lambda x,t:f.fields(c,x,t),lambda x,t:force(x,t,*c[-2:]),p[i1:i1+1],float(t),.00125,.000625)
        rows.append(dict(time=float(t),parent_max=float(n0[i0]),child_max=float(n1[i1]),parent_peak=p[i0].tolist(),child_peak=p[i1].tolist(),child_peak_fd=float(np.linalg.norm(fd[0])),fd_vector_difference=float(np.linalg.norm(fd[0]-r1[i1])),azimuthal_sampled_floor=float(abs(r0[it,1])),azimuthal_peak=p[it].tolist(),azimuthal_change=float(np.max(abs(r0[:,1]-r1[:,1])))))
        print('GRID',t,n0[i0],n1[i1],flush=True)
    result=dict(parent_sha256=hashlib.sha256(Path(parent).read_bytes()).hexdigest(),child_sha256=hashlib.sha256(Path(child).read_bytes()).hexdigest(),shape=[len(rr),len(zz)],times=len(ts),rows=rows,parent_max=max(r['parent_max'] for r in rows),child_max=max(r['child_max'] for r in rows),pressure_only_azimuthal_floor=max(r['azimuthal_sampled_floor'] for r in rows),scope='Frozen fields, independent grid and extra times. Finite samples not continuum upper bounds. At fixed u/f every axisymmetric pressure has exactly the same azimuthal residual; numerical floor is not an interval certificate. Independent Cartesian FD at every child grid peak.',pde_validated=False)
    save(out,result);return result
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--parent',required=True);p.add_argument('--child',required=True);p.add_argument('--out',required=True);a=p.parse_args();check(a.parent,a.child,a.out)
