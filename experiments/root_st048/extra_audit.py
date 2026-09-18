"""New off-grid shear and random space-time FD checks, never used for fitting."""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import boundary_shear
from spacetime import Family,force,NU
from structure_audit import jets

def velocity_pressure(f,raw,points,times):
    times=np.broadcast_to(np.asarray(times),(len(points),));u=np.empty_like(points);p=np.empty(len(points))
    for k in range(0,len(points),512):
        sl=slice(k,k+512);u[sl],p[sl]=f.fields(raw,points[sl],times[sl])
    return u,p

def spacetime_fd(f,raw,x,t,h=.0025,ht=.000625):
    """Independent Cartesian FD4, with central time stencil strictly inside window."""
    x=np.asarray(x,float);t=np.asarray(t,float)
    if x.ndim!=2 or x.shape[1]!=3 or t.shape!=(len(x),) or not np.isfinite(x).all() or not np.isfinite(t).all():raise ValueError('Invalid finite points/times')
    if h<=0 or ht<=0 or np.any((t-2*ht<.25)|(t+2*ht>.75)):raise ValueError('Invalid central time window')
    field=lambda X,T:velocity_pressure(f,raw,X,T)
    u,p=field(x,t);J=np.empty((len(x),3,3));lap=np.zeros_like(u);gp=np.zeros_like(u)
    for j in range(3):
        e=np.eye(3)[j]*h
        m2,qm2=field(x-2*e,t);m1,qm1=field(x-e,t);p1,qp1=field(x+e,t);p2,qp2=field(x+2*e,t)
        J[:,:,j]=(m2-8*m1+8*p1-p2)/(12*h)
        lap+=(-m2+16*m1-30*u+16*p1-p2)/(12*h*h)
        gp[:,j]=(qm2-8*qm1+8*qp1-qp2)/(12*h)
    ut=(field(x,t-2*ht)[0]-8*field(x,t-ht)[0]+8*field(x,t+ht)[0]-field(x,t+2*ht)[0])/(12*ht)
    R=ut+np.einsum('nij,nj->ni',J,u)+gp-NU*lap-force(x,t,*raw[-2:])
    return R,np.trace(J,axis1=1,axis2=2)

def audit(parent,candidate,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    f0,x0=Family.load(parent);f,x=Family.load(candidate)
    rng=np.random.default_rng(9174881);r,t=np.meshgrid(rng.uniform(.04,.6,31),np.r_[.25,np.sort(rng.uniform(.25,.75,19)),.75],indexing='ij')
    pjets=jets(f0,x0,r,0,t);cjets=jets(f,x,r,0,t)
    a,b=pjets[...,8],cjets[...,8];ratio=b/a;err=(b-a)
    shear=dict(seed=9174881,points=r.size,minimum_signed_ratio=float(ratio.min()),maximum_signed_ratio=float(ratio.max()),median_signed_ratio=float(np.median(ratio)),ratio_below_training_target_fraction=float(np.mean(ratio<.995)),same_shear_direction_fraction=float(np.mean(a*b>0)),parent_max_abs=float(abs(a).max()),child_max_abs=float(abs(b).max()),shear_difference_L2_sample=float(np.sqrt(np.mean(err**2))),note='New finite points; training target does not imply a continuum signed shear bound')
    rng=np.random.default_rng(9174882);points=rng.uniform(-2,2,(4096,3));times=rng.uniform(.255,.745,4096)
    rows=[]
    for ident,F,raw in [('parent',f0,x0),('child',f,x)]:
        R,d=spacetime_fd(F,raw,points,times);norm=np.linalg.norm(R,axis=1);i=int(np.argmax(norm))
        exact=np.concatenate([F.analytic_residual(raw,points[k:k+256],times[k:k+256]) for k in range(0,len(points),256)])
        rows.append(dict(id=ident,seed=9174882,random_interior_spacetime_points=4096,sampled_spacetime_max=float(norm.max()),volume_scaled_spacetime_RMS=float(np.sqrt(64*np.mean(norm**2))),divergence_max=float(abs(d).max()),FD_analytic_max_difference=float(np.linalg.norm(R-exact,axis=1).max()),max_point=points[i].tolist(),max_time=float(times[i]),norm_scope='Interior space-time average, NOT fixed-time spatial L2 and not a continuous supremum'))
    result=dict(parent_sha256=hashlib.sha256(Path(parent).read_bytes()).hexdigest(),candidate_sha256=hashlib.sha256(Path(candidate).read_bytes()).hexdigest(),offgrid_shear=shear,random_spacetime=rows,step=.0025,time_step=.000625,pde_validated=False)
    (out/'extra_audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)
    return result
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--parent',required=True);p.add_argument('--candidate',required=True);p.add_argument('--out',required=True);a=p.parse_args();audit(a.parent,a.candidate,a.out)
