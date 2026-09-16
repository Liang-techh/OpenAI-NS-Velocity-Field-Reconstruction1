"""Fit compact pressure corrections at fixed velocity and fixed force."""
import json
from pathlib import Path
import numpy as np
from scipy.optimize import lsq_linear, least_squares
from .constrained_outer_momentum import AngularMomentumCandidate
from .constrained_optimize import training_residual


def run(*, initial="artifacts/constrained/outer_shape/candidate.json", output="artifacts/constrained/outer_pressure", dense_cylindrical=False, fourth_power=False):
    family=json.loads(Path(initial).read_text())["family"]
    if family=="inner_swirl_v1":
        from .constrained_inner_swirl import InnerSwirlCandidate
        c=InnerSwirlCandidate.load(initial);pressure_owner=c.parent
    else:
        c=AngularMomentumCandidate.load(initial);pressure_owner=c
    rng=np.random.default_rng(20260916);x=rng.uniform(-2,2,(2048,3));t=rng.uniform(.252,.748,len(x))
    if dense_cylindrical:
        rr,zz,tt=np.meshgrid(np.linspace(.03,1.97,48),np.linspace(0,1.97,48),np.linspace(.252,.748,9),indexing="ij")
        x=np.vstack((x,np.column_stack((rr.ravel(),np.zeros(rr.size),zz.ravel()))));t=np.r_[t,tt.ravel()]
    h=.001;r=training_residual(c,c.force,x,t,.01,h)
    matrix=np.empty((len(x),3,18))
    for j in range(3):
        d=np.eye(3)[j]*h
        matrix[:,j,:]=(pressure_owner.pressure_basis(x+d,t)-pressure_owner.pressure_basis(x-d,t))/(2*h)
    A=matrix.reshape(-1,18)
    previous=np.asarray(pressure_owner.pressure_coefficients)
    fit=lsq_linear(A,-r.ravel(),bounds=(-100-previous,100-previous),tol=1e-10,max_iter=200)
    if fourth_power:
        def fun(a):
            value=r+matrix@a
            return (value*np.linalg.norm(value,axis=1)[:,None]).ravel()/np.sqrt(len(x))
        def jac(a):
            value=r+matrix@a;norms=np.linalg.norm(value,axis=1)
            inner=np.einsum('nc,nck->nk',value,matrix)
            inverse=np.divide(1.,norms,out=np.zeros_like(norms),where=norms>0)
            J=norms[:,None,None]*matrix+value[:,:,None]*inner[:,None,:]*inverse[:,None,None]
            return J.reshape(-1,18)/np.sqrt(len(x))
        fit=least_squares(fun,fit.x,jac=jac,bounds=(-100-previous,100-previous),max_nfev=150)
    if not fit.success:raise RuntimeError(fit.message)
    pressure_owner.pressure_coefficients=tuple(previous+fit.x)
    out=Path(output);out.mkdir(exist_ok=True);c.save(out/'candidate.json')
    after=training_residual(c,c.force,x,t,.01,h)
    report={'status':'pressure_fit_only_not_validated','force':{'a':c.force.a,'c':c.force.c},
        'dense_cylindrical':dense_cylindrical,'fourth_power':fourth_power,'training_points':len(x),'initial_artifact':str(initial),'seed':20260916,'basis':'B(r²/4)B(z²/4)(r²/4)^i(z²/4)^j[2(t-.25)]^k; i,j=0..2,k=0..1',
        'coefficient_bounds':[-100,100],'initial_L2_sample_mean':float(np.mean(np.sum(r*r,axis=1))),
        'final_L2_sample_mean':float(np.mean(np.sum(after*after,axis=1))),'velocity_changed':False}
    (out/'training.json').write_text(json.dumps(report,indent=2)+'\n');print(report)

if __name__=='__main__':run()
