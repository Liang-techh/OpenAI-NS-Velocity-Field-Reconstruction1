"""Pre-holdout bounded minimum-norm feasibility restoration and training-pool exchange."""
from __future__ import annotations
import argparse,time,json,hashlib
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from moment_step import MomentObjective,pool_residual
from aligned_continuation import EdgeModel
from pressure_morph import atomic_json

def run(parent,reference,source,out,max_steps=8,seconds=180):
    source,out=Path(source),Path(out)
    if out.exists():raise ValueError('Output exists')
    out.mkdir(parents=True)
    atomic_json(out/'registration.json',dict(utc_before_restore=datetime.now(timezone.utc).isoformat(),parent_sha256=hashlib.sha256(Path(parent).read_bytes()).hexdigest(),modifier_sha256=hashlib.sha256((source/'modifiers.npy').read_bytes()).hexdigest(),max_steps=max_steps,seconds=seconds,feasibility_tolerance=1e-7,selection='restore all original run-Q constraints; optionally add up to256 violated TRAININGpool points; no holdout used',original_physical_gates='UNCHANGED'))
    start=time.monotonic();m=EdgeModel(parent,edge_modes=True)
    o=MomentObjective(m,shear_reference=reference,robust_grid=True,poisson_weight=.003,morph_ratio=.9999,axis_ratio=1.75,shear_ratio=.999,edge_weight=.08,pressure_target=0.,profile_ratio=1.,pressure_ratio=1.,acceleration_ratio=1.,axis_weight=.02,space_order=(32,48),anchor_tolerance=.08,time_cap=1.,peak_weight=0.,softmax_weight=.03)
    c=np.load(source/'modifiers.npy');T=np.load(source/'solver_transform.npy');lo,hi=np.array(m.bounds).T;bs=1/np.maximum(hi-lo,1e-4)
    res=pool_residual(m,c,o.pool);extra=np.flatnonzero(res>o.cap*(1+1e-8));extra=extra[np.argsort(res[extra])[-256:]]
    active=np.union1d(o.active,extra);o.peakD=m.cache(*o.pool[active].T);o.last_r=o.last_c=None
    def con(c):
        v,J=o.constraints(c)
        return np.r_[v,(c-lo)*bs,(hi-c)*bs],np.vstack((J,np.diag(bs),-np.diag(bs)))
    history=[]
    for k in range(max_steps):
        v,J=con(c);old=float(v.min())
        if old>=-1e-7:break
        A=J@T;scale=np.clip(1/np.maximum(np.linalg.norm(A,axis=1),1e-4),.01,1e4)
        b=v*scale;A=A*scale[:,None]
        opt=minimize(lambda y:(.5*y@y,y),np.zeros(m.dim),jac=True,method='SLSQP',constraints=[dict(type='ineq',fun=lambda y:b+A@y,jac=lambda y:A)],options=dict(maxiter=80,ftol=1e-13))
        step=T@opt.x;alpha=1.;new=old;chosen=None
        for _ in range(10):
            trial=c+alpha*step;new=float(con(trial)[0].min())
            if new>old+1e-10 or new>=-1e-7:chosen=trial;break
            alpha*=.5
        row=dict(step=k+1,old_min=old,new_min=new,step_fraction=alpha,QP_success=bool(opt.success),QP_iterations=int(opt.nit),QP_message=str(opt.message),elapsed=time.monotonic()-start);history.append(row)
        if chosen is None:break
        c=chosen;np.save(out/'checkpoint.npy',c);m.f.save(m.candidate(c),out/'checkpoint.json',row);atomic_json(out/'history.json',history);print('REPAIR',row,flush=True)
        if time.monotonic()-start>seconds:break
    v,_=con(c);res=pool_residual(m,c,o.pool)
    valid=bool(v.min()>=-1e-7 and np.all((c>=lo-1e-10)&(c<=hi+1e-10)))
    summary=dict(active_constraint_feasible=valid,min_constraint=float(v.min()),steps=len(history),active_count=len(active),added_points=len(extra),pool_max=float(res.max()),parent_pool_max=o.cap,whole_pool_peak_pass=bool(res.max()<=o.cap*(1+1e-7)),loss=o.fun(c)[0],elapsed=time.monotonic()-start,pde_validated=False,source_correspondence_verified=False)
    np.save(out/'modifiers.npy',c);np.save(out/'pool_residual.npy',res);m.f.save(m.candidate(c),out/'candidate.json',summary)
    atomic_json(out/'summary.json',summary);atomic_json(out/'history.json',history);print('SUMMARY',json.dumps(summary),flush=True)
    return valid
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--parent',required=True);p.add_argument('--reference',required=True);p.add_argument('--source',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();run(a.parent,a.reference,a.source,a.out)
