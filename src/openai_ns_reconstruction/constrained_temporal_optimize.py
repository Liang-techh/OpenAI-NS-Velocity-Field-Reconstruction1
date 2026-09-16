"""Bounded zero-moment temporal swirl fit; original force and core preserved."""
import json
from dataclasses import asdict
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from .constrained_outer_momentum import AngularMomentumCandidate
from .constrained_temporal_swirl import TemporalSwirlCandidate
from .constrained_optimize import training_residual
from .constrained_temporal_cache import cached_residual


def run(initial='artifacts/constrained/continued_pressure/candidate.json', output='artifacts/constrained/temporal_swirl', call_budget=500, tail_adaptive=False):
    parent=AngularMomentumCandidate.load(initial)
    rng=np.random.default_rng(20260916)
    x=rng.uniform(-2,2,(2048,3));t=rng.uniform(.252,.748,len(x))
    if tail_adaptive:
        pool_rng=np.random.default_rng(20261017)
        pool=pool_rng.uniform(-2,2,(8192,3));pool_t=np.full(len(pool),.748)
        magnitudes=np.linalg.norm(training_residual(parent,parent.force,pool,pool_t,.01),axis=1)
        chosen=np.argsort(magnitudes)[-256:]
        x=np.vstack((x,pool[chosen]));t=np.r_[t,pool_t[chosen]]
    candidate=TemporalSwirlCandidate(parent)
    evaluate=cached_residual(candidate,x,t)
    probe=np.linspace(-.2,.2,15)
    direct=training_residual(TemporalSwirlCandidate(parent,tuple(probe)),parent.force,x,t,.01)
    cache_error=float(np.max(np.abs(evaluate(probe)-direct)))
    if cache_error>1e-7:raise RuntimeError(f"residual cache discrepancy {cache_error}")
    calls=0;best={'loss':float('inf')};history=[]
    class Budget(Exception):pass
    def fun(v):
        nonlocal calls
        if calls>=call_budget:raise Budget()
        calls+=1
        c=TemporalSwirlCandidate(parent,tuple(v-2))
        r=evaluate(v-2).ravel()/np.sqrt(len(x))
        energies=np.array([c.energy(float(tt),24) for tt in np.linspace(.25,.75,9)])
        violation=np.r_[np.maximum(.1005-energies,0),np.maximum(energies-9.999,0)]
        loss=float(r@r)
        if not np.any(violation) and loss<best['loss']:
            best.update(loss=loss,coefficients=tuple(v-2));history.append({'call':calls,'loss':loss})
        return np.r_[r,100*violation]
    try:
        fit=least_squares(fun,np.full(15,2.),bounds=(np.ones(15),np.full(15,3.)),diff_step=1e-4,max_nfev=50,ftol=1e-8)
        reason=fit.message
    except Budget:reason=f'actual {call_budget}-call budget reached'
    if 'coefficients' not in best:raise RuntimeError('no energy-feasible candidate')
    result=TemporalSwirlCandidate(parent,best['coefficients'])
    out=Path(output);out.mkdir(parents=True,exist_ok=True);result.save(out/'candidate.json')
    report={'status':'training_only_not_validated','initial_artifact':initial,'force':asdict(parent.force),
        'tail_adaptive':tail_adaptive,'adaptive_seed':20261017 if tail_adaptive else None,'calls':calls,'call_budget':call_budget,'termination':str(reason),'seed':20260916,
        'cache_direct_max_discrepancy':cache_error,'training_points':len(x),'training_loss':best['loss'],'history':history,
        'coefficient_bounds':[-1,1],'basis':'five zero-moment spatial modes times three cubic Bernstein temporal modes',
        'initial_velocity_max_change':float(np.max(np.abs(result.velocity(x,.25)-parent.velocity(x,.25))))}
    (out/'training.json').write_text(json.dumps(report,indent=2)+'\n');print(report)

if __name__=='__main__':run()
