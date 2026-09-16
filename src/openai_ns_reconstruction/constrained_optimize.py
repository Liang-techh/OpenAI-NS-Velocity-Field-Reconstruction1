"""Bounded training run using second-order differences, separate from validation."""
import json
from pathlib import Path
from dataclasses import asdict
import numpy as np
from scipy.optimize import least_squares
from .constrained_candidate import CompactCandidate
from .constrained_force import RestrictedForce


def training_residual(c, f, x, t, nu, h=0.001):
    u=c.velocity(x,t)
    ut=(c.velocity(x,t+h)-c.velocity(x,t-h))/(2*h)
    conv=np.zeros_like(u);lap=np.zeros_like(u);gp=np.zeros_like(u)
    for j in range(3):
        d=np.eye(3)[j]*h
        up,um=c.velocity(x+d,t),c.velocity(x-d,t)
        conv+=u[:,j,None]*(up-um)/(2*h)
        lap+=(up-2*u+um)/h**2
        gp[:,j]=(c.pressure(x+d,t)-c.pressure(x-d,t))/(2*h)
    return ut+conv+gp-nu*lap-f(x,t)


def run(config='configs/constraints.json', output='artifacts/constrained/optimized'):
    cfg=json.loads(Path(config).read_text()); opt=cfg['optimization']
    rng=np.random.default_rng(opt['seed']);n=opt['interior_points']
    box=np.asarray(cfg['domain']['evaluation_box']);lo,hi=cfg['domain']['time_interval']
    x=rng.uniform(box[:,0],box[:,1],(n,3));t=rng.uniform(lo+0.002,hi-0.002,n)
    # Independent core/axis strata supplement the uniform training samples.
    extra=rng.uniform(-0.3,0.3,(256,3));extra[:64,:2]=0
    x=np.vstack((x,extra));t=np.r_[t,rng.uniform(lo+0.002,hi-0.002,len(extra))]
    names=['swirl_ratio','radial_width','axial_width','radial_shape','axial_shape',
           'pressure_constant','pressure_radial','pressure_axial']
    bounds=opt['candidate_parameter_bounds']
    lows=[bounds[k][0] for k in names]+[0,0];highs=[bounds[k][1] for k in names]+[10,10]
    initial=CompactCandidate();v0=np.array([getattr(initial,k) for k in names]+[1,1])
    best={'loss':float('inf')};history=[];calls=0
    class BudgetReached(Exception): pass
    def decode(v):
        return CompactCandidate(**dict(zip(names,v[:8]))).normalized(),RestrictedForce(*v[8:])
    def fun(v):
        nonlocal calls
        if calls>=opt['maximum_function_evaluations']: raise BudgetReached()
        calls+=1;c,f=decode(v)
        pde=training_residual(c,f,x,t,cfg['nu']).ravel()/np.sqrt(len(x))
        times=np.linspace(lo,hi,9)
        energy=np.array([c.energy(float(tt),24) for tt in times])
        emin=cfg['nontriviality']['minimum_energy_each_validation_time']
        emax=cfg['nontriviality']['maximum_energy_each_validation_time']
        ep=np.r_[np.maximum(emin-energy,0),np.maximum(energy-emax,0)]/np.sqrt(len(times))
        probes=np.column_stack((0.1*np.sqrt(1-times),np.zeros(9),0.1*(1-times)**0.495))
        u=c.velocity(probes,times)
        signs=np.column_stack((np.maximum(u[:,0],0),np.maximum(-u[:,1],0),np.maximum(-u[:,2],0))).ravel()
        w=opt['loss_weights'];r=np.r_[np.sqrt(w['pde'])*pde,np.sqrt(w['energy'])*ep,np.sqrt(w['structure'])*signs]
        loss=float(r@r)
        if loss<best['loss']:
            best.update(loss=loss,parameters=v.copy())
            history.append({'call':calls,'loss':loss})
        return r
    try:
        result=least_squares(fun,v0,bounds=(lows,highs),max_nfev=150,ftol=1e-8,xtol=1e-8,gtol=1e-8,diff_step=1e-4)
        termination=result.message
    except BudgetReached:
        termination='actual function evaluation budget reached'
    c,f=decode(best['parameters']);out=Path(output);out.mkdir(parents=True,exist_ok=True)
    c.save(out/'candidate.json')
    report={'status':'training_only_not_validated','seed':opt['seed'],'uniform_samples':n,
        'additional_core_axis_samples':256,'calls':calls,'termination':termination,
        'best_training_loss':best['loss'],'force':asdict(f),'parameters':asdict(c),
        'history':history,'training_spatial_time_step':0.001,'normalization_order':96,
        'solver':'scipy.optimize.least_squares','max_nfev':150,'actual_call_budget':opt['maximum_function_evaluations']}
    (out/'training.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ('calls','termination','best_training_loss','force')}))

if __name__=='__main__':
    run()
