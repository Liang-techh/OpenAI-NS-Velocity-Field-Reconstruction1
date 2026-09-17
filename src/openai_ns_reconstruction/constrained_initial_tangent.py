"""Equation-driven initial tangent fit, not a full-time solution."""
from dataclasses import replace
import json
from pathlib import Path
import numpy as np
from scipy.optimize import lsq_linear
from .constrained_outer_momentum import AngularMomentumCandidate


def initial_residual(c,x,step=0.0001):
    t=.25;u=c.velocity(x,t)
    ut=(-3*u+4*c.velocity(x,t+step)-c.velocity(x,t+2*step))/(2*step)
    conv=np.zeros_like(u);lap=np.zeros_like(u);gp=np.zeros_like(u)
    for j in range(3):
        d=np.eye(3)[j]*step;up=c.velocity(x+d,t);um=c.velocity(x-d,t)
        conv+=u[:,j,None]*(up-um)/(2*step);lap+=(up-2*u+um)/step**2
        gp[:,j]=(c.pressure(x+d,t)-c.pressure(x-d,t))/(2*step)
    return ut+conv+gp-.01*lap-c.force(x,t)


def run():
    old=AngularMomentumCandidate.load('artifacts/constrained/outer_shape/candidate.json')
    indices=[n for n,(_,_,k) in enumerate(old.base.BASIS) if k==1]
    def candidate(v):
        pol=list(old.base.poloidal_coefficients);sw=list(old.base.swirl_coefficients)
        for n,i in enumerate(indices):pol[i]=float(v[n]);sw[i]=float(v[n+len(indices)])
        base=replace(old.base,poloidal_coefficients=tuple(pol),swirl_coefficients=tuple(sw))
        return AngularMomentumCandidate(base,old.force,96,old.outer_shape,old.pressure_coefficients)
    rng=np.random.default_rng(20260916);x=rng.uniform(-2,2,(2048,3))
    c0=candidate(np.zeros(18));r0=initial_residual(c0,x)
    columns=[]
    for n in range(18):
        e=np.zeros(18);e[n]=1
        columns.append((initial_residual(candidate(e),x)-r0).ravel())
    matrix=np.column_stack(columns)
    fit=lsq_linear(matrix,-r0.ravel(),bounds=(-1,1),tol=1e-10,max_iter=200)
    if not fit.success:raise RuntimeError(fit.message)
    c=candidate(fit.x);r=initial_residual(c,x)
    # Verify actual residual against the affine model used by the linear fit.
    mismatch=float(np.max(np.abs(r.ravel()-(r0.ravel()+matrix@fit.x))))
    if mismatch>1e-5:raise RuntimeError('initial residual is not sufficiently affine')
    out=Path('artifacts/constrained/initial_tangent');out.mkdir(exist_ok=True);c.save(out/'candidate.json')
    report={'status':'initial_tangent_only_not_validated','seed':20260916,'time':.25,
        'force':{'a':c.force.a,'c':c.force.c},'active_indices':indices,'bounds':[-1,1],
        'old_initial_sample_max':float(np.max(np.linalg.norm(initial_residual(old,x),axis=1))),
        'new_initial_sample_max':float(np.max(np.linalg.norm(r,axis=1))),
        'affine_model_error':mismatch,'initial_velocity_max_change':float(np.max(np.abs(c.velocity(x,.25)-old.velocity(x,.25))))}
    (out/'training.json').write_text(json.dumps(report,indent=2)+'\n');print(report)

if __name__=='__main__':run()
