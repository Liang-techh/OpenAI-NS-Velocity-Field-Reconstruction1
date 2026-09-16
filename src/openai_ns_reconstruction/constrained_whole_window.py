"""Whole-window time-coefficient fit with initial field and prescribed torque fixed."""
from dataclasses import replace
import json
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares, minimize, LinearConstraint
from .constrained_outer_momentum import AngularMomentumCandidate
from .constrained_optimize import training_residual
from .constrained_validation import structure_metrics


def run(core_equalities=False, *, initial='artifacts/constrained/outer_shape/candidate.json', output=None, call_budget=600, max_iterations=20):
    old=AngularMomentumCandidate.load(initial)
    config=json.loads(Path('configs/constraints_tensor.json').read_text())
    indices=[n for n,(i,j,k) in enumerate(old.base.BASIS) if i+j<=2 and k>=1]
    size=len(indices)
    v0=np.r_[[old.base.poloidal_coefficients[i]+2 for i in indices],
              [old.base.swirl_coefficients[i]+2 for i in indices]]
    rng=np.random.default_rng(20260916);x=rng.uniform(-2,2,(2048,3));t=rng.uniform(.252,.748,len(x))
    calls=0;best={'loss':float('inf')};history=[]
    def candidate(v):
        pol=list(old.base.poloidal_coefficients);sw=list(old.base.swirl_coefficients)
        for n,i in enumerate(indices):pol[i]=float(v[n]-2);sw[i]=float(v[n+size]-2)
        base=replace(old.base,poloidal_coefficients=tuple(pol),swirl_coefficients=tuple(sw))
        return AngularMomentumCandidate(base,old.force,96,old.outer_shape,old.pressure_coefficients)
    class Budget(Exception):pass
    def fun(v):
        nonlocal calls
        if calls>=call_budget:raise Budget()
        calls+=1;c=candidate(v)
        r=training_residual(c,c.force,x,t,.01).ravel()/np.sqrt(len(x))
        ts=np.linspace(.25,.75,9)
        energy=np.array([c.energy(float(tt),24) for tt in ts])
        points=np.column_stack((.1*np.sqrt(1-ts),np.zeros(9),.1*(1-ts)**.495))
        u=c.velocity(points,ts);scaled=u*np.column_stack((np.sqrt(1-ts),(1-ts)**.505,(1-ts)**.505))
        drift=np.linalg.norm(scaled-scaled[0],axis=1)/np.linalg.norm(scaled[0])
        violation=np.r_[np.maximum(.1005-energy,0),np.maximum(energy-9.999,0),np.maximum(drift-.0495,0),np.maximum(u[:,0],0),np.maximum(-u[:,1:],0).ravel()]
        loss=float(r@r)
        feasible=np.max(violation)==0 and np.all(u[:,0]<0) and np.all(u[:,1:]>0)
        if feasible and loss<best['loss']:
            best.update(loss=loss,v=v.copy());history.append({'call':calls,'loss':loss})
        return np.r_[r,100*violation]
    try:
        if core_equalities:
            ts=np.array([.4,.65]);p=np.column_stack((.1*np.sqrt(1-ts),np.zeros(2),.1*(1-ts)**.495))
            initial_core=candidate(v0).velocity(p,ts).ravel();columns=[]
            for j in range(len(v0)):
                h=.001 if v0[j]<2.999 else -.001
                trial=v0.copy();trial[j]+=h
                columns.append((candidate(trial).velocity(p,ts).ravel()-initial_core)/h)
            C=np.column_stack(columns)
            C=C/np.linalg.norm(C,axis=1)[:,None]
            def objective(v):
                r=fun(v);return float(r@r)
            result=minimize(objective,v0,method='SLSQP',bounds=[(1.,3.)]*len(v0),
                constraints=[LinearConstraint(C,C@v0,C@v0)],
                options={'maxiter':max_iterations,'eps':1e-4,'ftol':1e-9})
        else:
            result=least_squares(fun,v0,bounds=(np.ones(2*size),np.full(2*size,3.)),max_nfev=80,diff_step=1e-4)
        reason=result.message
    except Budget:reason=f'actual {call_budget}-call budget reached'
    if 'v' not in best:raise RuntimeError('no feasible iterate')
    c=candidate(best['v']);out=Path(output or ('artifacts/constrained/whole_window_equalities' if core_equalities else 'artifacts/constrained/whole_window'));out.mkdir(exist_ok=True);c.save(out/'candidate.json')
    report={'status':'training_only_not_validated','force':{'a':c.force.a,'c':c.force.c},
        'initial_artifact':initial,'call_budget':call_budget,'max_iterations':max_iterations,'core_equalities':core_equalities,'calls':calls,'termination':reason,'training_loss':best['loss'],'history':history,
        'active_time_indices':indices,'seed':20260916,'initial_velocity_max_change':float(np.max(np.abs(c.velocity(x,.25)-old.velocity(x,.25))))}
    (out/'training.json').write_text(json.dumps(report,indent=2)+'\n');print(report['calls'],report['training_loss'])

if __name__=='__main__':run()
