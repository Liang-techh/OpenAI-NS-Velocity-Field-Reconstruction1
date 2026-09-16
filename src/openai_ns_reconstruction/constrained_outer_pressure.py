"""Fit compact pressure corrections at fixed velocity and fixed force."""
import json
from pathlib import Path
import numpy as np
from scipy.optimize import lsq_linear
from .constrained_outer_momentum import AngularMomentumCandidate
from .constrained_optimize import training_residual


def run(*, initial="artifacts/constrained/outer_shape/candidate.json", output="artifacts/constrained/outer_pressure"):
    c=AngularMomentumCandidate.load(initial)
    rng=np.random.default_rng(20260916);x=rng.uniform(-2,2,(2048,3));t=rng.uniform(.252,.748,len(x))
    h=.001;r=training_residual(c,c.force,x,t,.01,h)
    matrix=np.empty((len(x),3,18))
    for j in range(3):
        d=np.eye(3)[j]*h
        matrix[:,j,:]=(c.pressure_basis(x+d,t)-c.pressure_basis(x-d,t))/(2*h)
    A=matrix.reshape(-1,18)
    previous=np.asarray(c.pressure_coefficients)
    fit=lsq_linear(A,-r.ravel(),bounds=(-100-previous,100-previous),tol=1e-10,max_iter=200)
    if not fit.success:raise RuntimeError(fit.message)
    c.pressure_coefficients=tuple(previous+fit.x)
    out=Path(output);out.mkdir(exist_ok=True);c.save(out/'candidate.json')
    after=training_residual(c,c.force,x,t,.01,h)
    report={'status':'pressure_fit_only_not_validated','force':{'a':c.force.a,'c':c.force.c},
        'initial_artifact':str(initial),'seed':20260916,'basis':'B(r²/4)B(z²/4)(r²/4)^i(z²/4)^j[2(t-.25)]^k; i,j=0..2,k=0..1',
        'coefficient_bounds':[-100,100],'initial_L2_sample_mean':float(np.mean(np.sum(r*r,axis=1))),
        'final_L2_sample_mean':float(np.mean(np.sum(after*after,axis=1))),'velocity_changed':False}
    (out/'training.json').write_text(json.dumps(report,indent=2)+'\n');print(report)

if __name__=='__main__':run()
