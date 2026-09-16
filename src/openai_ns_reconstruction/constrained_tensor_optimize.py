"""Bounded tensor-correction fit; keeps base geometry and force family fixed."""
import json
from pathlib import Path
from dataclasses import asdict
import numpy as np
from scipy.optimize import least_squares
from .constrained_candidate import CompactCandidate
from .constrained_force import RestrictedForce
from .constrained_optimize import training_residual


def run(output='artifacts/constrained/tensor_stage1', *, feasible_selection=False):
    from .constrained_tensor_candidate import TensorCandidate
    cfg=json.loads(Path('configs/constraints_tensor.json').read_text())
    out=Path(output)
    base=CompactCandidate.load(cfg['tensor_experiment']['initial_candidate'])
    old=json.loads(Path(cfg['tensor_experiment']['initial_candidate']).with_name('training.json').read_text())
    rng=np.random.default_rng(cfg['optimization']['seed'])
    x=rng.uniform(-2,2,(2048,3));t=rng.uniform(.252,.748,len(x))
    extra=rng.uniform(-.3,.3,(256,3));extra[:64,:2]=0
    x=np.vstack((x,extra));t=np.r_[t,rng.uniform(.252,.748,len(extra))]
    indices=[n for n,(i,j,k) in enumerate((i,j,k) for i in range(3) for j in range(3) for k in range(3)) if i+j<=2 and k<=1]
    size=len(indices);pn=('pressure_constant','pressure_radial','pressure_axial')
    # Shift coefficient variables by two, so relative finite differences at
    # zero corrections still use a meaningful physical parameter step.
    v0=np.r_[np.full(2*size,2.),[getattr(base,n) for n in pn],old['force']['a'],old['force']['c']]
    lower=np.r_[np.ones(2*size),[-100]*3,[0,0]]
    upper=np.r_[np.full(2*size,3.),[100]*3,[10,10]]
    calls=0;best={'loss':float('inf')};feasible_best={'loss':float('inf')};history=[]
    class Budget(Exception):pass
    def decode(v):
        pol=np.zeros(27);sw=np.zeros(27)
        pol[indices]=v[:size]-2;sw[indices]=v[size:2*size]-2
        args=asdict(base);args.update(dict(zip(pn,v[2*size:2*size+3])))
        c=TensorCandidate(**args,poloidal_coefficients=tuple(pol),swirl_coefficients=tuple(sw)).normalized()
        return c,RestrictedForce(*v[-2:])
    def fun(v):
        nonlocal calls
        if calls>=cfg['optimization']['maximum_function_evaluations']:raise Budget()
        calls+=1;c,f=decode(v)
        rr=training_residual(c,f,x,t,cfg['nu']).ravel()/np.sqrt(len(x))
        ts=np.linspace(.25,.75,9)
        energy=np.array([c.energy(float(tt),24) for tt in ts])
        ep=np.r_[np.maximum(.1-energy,0),np.maximum(energy-10,0)]/3
        probes=np.column_stack((.1*np.sqrt(1-ts),np.zeros(9),.1*(1-ts)**.495))
        u=c.velocity(probes,ts)
        scaled=u*np.column_stack((np.sqrt(1-ts),(1-ts)**.505,(1-ts)**.505))
        drift=np.linalg.norm(scaled-scaled[0],axis=1)/np.linalg.norm(scaled[0])
        constraints=np.r_[ep,np.maximum(u[:,0],0),np.maximum(-u[:,1:],0).ravel(),np.maximum(drift-.05,0)]
        result=np.r_[rr,(100 if feasible_selection else 1)*constraints];loss=float(result@result)
        pde_loss=float(rr@rr)
        admissible=bool(np.min(energy)>=.1005 and np.max(energy)<=9.999
                        and np.max(drift)<=.0495 and np.all(u[:,0]<0) and np.all(u[:,1:]>0))
        if admissible and pde_loss<feasible_best['loss']:
            feasible_best.update(loss=pde_loss,parameters=v.copy(),call=calls)
        if loss<best['loss']:
            best.update(loss=loss,parameters=v.copy());history.append({'call':calls,'loss':loss})
        return result
    try:
        fit=least_squares(fun,v0,bounds=(lower,upper),diff_step=1e-4,max_nfev=100,ftol=1e-8,xtol=1e-8,gtol=1e-8)
        reason=fit.message
    except Budget:reason='actual call budget reached'
    selected=feasible_best if feasible_selection else best
    if 'parameters' not in selected:
        raise RuntimeError('No structurally feasible sampled candidate found; no artifact accepted')
    c,f=decode(selected['parameters']);out.mkdir(parents=True,exist_ok=True);c.save(out/'candidate.json')
    result={'status':'training_only_not_validated','calls':calls,'termination':reason,'loss':selected['loss'],'weighted_best_loss':best['loss'],
            'feasible_selection':feasible_selection,'structure_multiplier':100 if feasible_selection else 1,
            'selection_limits':{'energy_min':.1005,'energy_max':9.999,'drift_max':.0495},
            'active_tensor_indices':indices,'force':asdict(f),'history':history,'config':'configs/constraints_tensor.json'}
    (out/'training.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k] for k in ('calls','termination','loss')}))

if __name__=='__main__':run()
