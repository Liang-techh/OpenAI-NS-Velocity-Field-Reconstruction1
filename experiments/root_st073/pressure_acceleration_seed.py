"""Peak-constrained pressure-only endpoint seed; velocity/shape unchanged."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):
    os.environ[key]='1'
import json
import hashlib
from pathlib import Path
import numpy as np
from scipy.optimize import minimize

ROOT=Path(__file__).resolve().parent


def run():
    path=ROOT/'acceleration_momentum_cache.npz'
    with np.load(path,allow_pickle=False) as cache:
        R=cache['base_residual'];w=cache['weights'];D=cache['design']
        convention=str(cache['design_convention'])
        source_projection=str(cache['projection_sha256'])
    if convention != 'real-minus-imag':
        raise ValueError('Regenerate acceleration cache with corrected complex packing')
    if D.shape != (R.size, 324) or w.shape != (len(R),):
        raise ValueError('Unexpected acceleration cache dimensions')
    indices=np.r_[np.arange(27,36),np.arange(90,108),np.arange(162,180),np.arange(234,252),np.arange(306,324)]
    D=D[:,indices]
    sw=np.repeat(np.sqrt(w),3);L2=float(np.linalg.norm(R.ravel()*sw))
    scales=np.linalg.norm(D*sw[:,None],axis=0)
    active=scales>1e-30;indices=indices[active];D=D[:,active];scales=scales[active]
    B=D/scales[None,:];A=B*sw[:,None]
    b=R.ravel()*sw/L2
    H=A.T@A;g=A.T@b;ridge=1e-8
    def objective(y):return .5*float(y@H@y)+float(g@y)+.5*ridge*float(y@y)
    def gradient(y):return H@y+g+ridge*y
    norm0=np.linalg.norm(R,axis=1);cap=float(norm0.max())*(1+1e-6)
    pool=set(np.argsort(norm0)[-20:].tolist());best=np.zeros(len(indices));best_value=objective(best)
    rows=[]
    for round_index in range(5):
        ids=np.array(sorted(pool));P=B.reshape(-1,3,len(indices))[ids]*L2/cap;offset=R[ids]/cap
        def constraint(y):
            r=offset+np.einsum('nij,j->ni',P,y)
            return 1-np.sum(r*r,axis=1)
        def jac(y):
            r=offset+np.einsum('nij,j->ni',P,y)
            return -2*np.einsum('ni,nij->nj',r,P)
        fit=minimize(objective,best,jac=gradient,method='SLSQP',
            constraints=[dict(type='ineq',fun=constraint,jac=jac)],
            options=dict(maxiter=150,ftol=1e-12))
        residual=R+(B@fit.x*L2).reshape(-1,3);norms=np.linalg.norm(residual,axis=1)
        violations=np.flatnonzero(norms>cap*(1+1e-10))
        feasible=len(violations)==0
        if feasible and objective(fit.x)<best_value:
            best=fit.x.copy();best_value=objective(best)
        rows.append(dict(round=round_index+1,optimizer_success=bool(fit.success),message=str(fit.message),
            peak=float(norms.max()),L2=float(np.sqrt(w@np.sum(residual**2,axis=1))),violations=len(violations)))
        if feasible:break
        pool.update(np.argsort(norms)[-20:].tolist())
    control=np.zeros(324);control[indices]=best*L2/scales
    residual=R+(B@best*L2).reshape(-1,3)
    report=dict(status='completed',accepted=False,pde_validated=False,scale_recursion_established=False,
        scope='Pressure time slope only; endpoint velocity and all velocity-based shapes unchanged. Sampled linear pressure response, no independent momentum acceptance.',
        cache_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),design_convention=convention,
        projection_sha256=source_projection,coefficients=control.tolist(),pressure_columns=indices.tolist(),
        base_L2=L2,base_max=float(norm0.max()),selected_L2=float(np.sqrt(w@np.sum(residual**2,axis=1))),
        selected_max=float(np.max(np.linalg.norm(residual,axis=1))),peak_cap=cap,rounds=rows,
        endpoint_velocity_unchanged=True)
    (ROOT/'pressure_acceleration_seed.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('coefficients','pressure_columns')}),flush=True)


if __name__=='__main__':run()
