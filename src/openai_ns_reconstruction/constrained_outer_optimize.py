"""Fit two outer shape coefficients; force and global torque constraint fixed."""
import json
from pathlib import Path
from scipy.optimize import least_squares
import numpy as np
from .constrained_outer_momentum import AngularMomentumCandidate
from .constrained_optimize import training_residual


def run():
    c=AngularMomentumCandidate.load('artifacts/constrained/outer_momentum/candidate.json')
    rng=np.random.default_rng(20260916)
    x=rng.uniform(-2,2,(2048,3));t=rng.uniform(.252,.748,len(x))
    calls=0;history=[]
    def fun(v):
        nonlocal calls
        calls+=1
        q=AngularMomentumCandidate(c.base,c.force,96,tuple(v))
        r=training_residual(q,c.force,x,t,.01).ravel()/np.sqrt(len(x))
        history.append({'call':calls,'loss':float(r@r),'shape':v.tolist()})
        return r
    result=least_squares(fun,[0.,0.],bounds=([-2,-2],[2,2]),max_nfev=40,diff_step=1e-4)
    q=AngularMomentumCandidate(c.base,c.force,96,tuple(result.x))
    out=Path('artifacts/constrained/outer_shape');out.mkdir(exist_ok=True);q.save(out/'candidate.json')
    (out/'training.json').write_text(json.dumps({'status':'training_only_not_validated',
        'force':{'a':c.force.a,'c':c.force.c},'calls':calls,'history':history,
        'termination':result.message,'shape_bounds':[-2,2],'seed':20260916},indent=2)+'\n')
    print(result.x,history[-1]['loss'],calls)

if __name__=='__main__':run()
