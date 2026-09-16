"""Bounded zero-moment temporal swirl fit; original force and core preserved."""
import json
from dataclasses import asdict
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from .constrained_outer_momentum import AngularMomentumCandidate
from .constrained_temporal_swirl import TemporalSwirlCandidate
from .constrained_optimize import training_residual
from .constrained_temporal_cache import cached_residual, cached_energy


def run(initial='artifacts/constrained/continued_pressure/candidate.json', output='artifacts/constrained/temporal_swirl', call_budget=500, tail_adaptive=False, localized=False, fourth_power=False, axial=False, dense_cylindrical=False, analytic_jacobian=False, coefficient_initial=None, grid_resolution=24, quintic=False, azimuthal_only=False):
    parent=AngularMomentumCandidate.load(initial)
    cls=TemporalSwirlCandidate
    if localized:
        from .constrained_localized_swirl import LocalizedSwirlCandidate
        cls=LocalizedSwirlCandidate
    if axial:
        from .constrained_axial_swirl import AxialSwirlCandidate
        from .constrained_localized_swirl import LocalizedSwirlCandidate
        cls=AxialSwirlCandidate
    if quintic:
        from .constrained_quintic_swirl import QuinticSwirlCandidate, elevate_cubic
        from .constrained_axial_swirl import AxialSwirlCandidate
        cls=QuinticSwirlCandidate
    count=cls.COEFFICIENT_COUNT
    rng=np.random.default_rng(20260916)
    x=rng.uniform(-2,2,(2048,3));t=rng.uniform(.252,.748,len(x))
    if dense_cylindrical:
        rr,zz,tt=np.meshgrid(np.linspace(.03,1.97,grid_resolution),np.linspace(0,1.97,grid_resolution),np.linspace(.252,.748,9),indexing="ij")
        extra=np.column_stack((rr.ravel(),np.zeros(rr.size),zz.ravel()))
        x=np.vstack((x,extra));t=np.r_[t,tt.ravel()]
    if tail_adaptive:
        pool_rng=np.random.default_rng(20261017)
        pool=pool_rng.uniform(-2,2,(8192,3));pool_t=np.full(len(pool),.748)
        magnitudes=np.linalg.norm(training_residual(parent,parent.force,pool,pool_t,.01),axis=1)
        chosen=np.argsort(magnitudes)[-256:]
        x=np.vstack((x,pool[chosen]));t=np.r_[t,pool_t[chosen]]
    radius=np.linalg.norm(x[:,:2],axis=1)
    tangents=np.column_stack((-x[:,1],x[:,0],np.zeros(len(x))))/np.maximum(radius[:,None],1e-30)
    candidate=cls(parent)
    evaluate=cached_residual(candidate,x,t)
    energy=cached_energy(candidate,np.linspace(.25,.75,9))
    probe=np.linspace(-.2,.2,count)
    direct=training_residual(cls(parent,tuple(probe)),parent.force,x,t,.01)
    cache_error=float(np.max(np.abs(evaluate(probe)-direct)))
    if cache_error>1e-7:raise RuntimeError(f"residual cache discrepancy {cache_error}")
    start=np.full(count,2.)
    if coefficient_initial:
        family=json.loads(Path(coefficient_initial).read_text())["family"]
        previous=AxialSwirlCandidate.load(coefficient_initial) if quintic and family=="axial_swirl_v1" else cls.load(coefficient_initial)
        from .constrained_temporal_swirl import _parent_payload
        if _parent_payload(previous.parent)!=_parent_payload(parent):raise ValueError("warm start parent mismatch")
        start+=np.asarray(elevate_cubic(previous.coefficients) if quintic and len(previous.coefficients)==36 else previous.coefficients)
    elif axial:
        previous=LocalizedSwirlCandidate.load('artifacts/constrained/localized_swirl/candidate.json')
        from .constrained_temporal_swirl import _parent_payload
        if _parent_payload(previous.parent)!=_parent_payload(parent):
            raise ValueError('axial warm start parent differs from requested parent')
        start[:18]+=np.asarray(previous.coefficients)
    calls=0;best={'loss':float('inf')};history=[]
    class Budget(Exception):pass
    def fun(v):
        nonlocal calls
        if calls>=call_budget:raise Budget()
        calls+=1
        raw=evaluate(v-2)
        if azimuthal_only:raw=np.sum(raw*tangents,axis=1)[:,None]
        if fourth_power:raw=raw*np.linalg.norm(raw,axis=1)[:,None]
        r=raw.ravel()/np.sqrt(len(x))
        energies=energy(v-2)
        violation=np.r_[np.maximum(.1005-energies,0),np.maximum(energies-9.999,0)]
        loss=float(r@r)
        if not np.any(violation) and loss<best['loss']:
            best.update(loss=loss,coefficients=tuple(v-2));history.append({'call':calls,'loss':loss})
        return np.r_[r,100*violation]
    def jac(v):
        a=v-2;raw=evaluate(a);J=evaluate.jacobian(a)
        if azimuthal_only:
            raw=np.sum(raw*tangents,axis=1)[:,None]
            J=np.einsum("nc,nck->nk",tangents,J)[:,None,:]
        if fourth_power:
            norms=np.linalg.norm(raw,axis=1)
            projected=np.einsum('nc,nck->nk',raw,J)
            inverse=np.divide(1.,norms,out=np.zeros_like(norms),where=norms>0)
            J=norms[:,None,None]*J+raw[:,:,None]*projected[:,None,:]*inverse[:,None,None]
        e=energy(a);de=energy.jacobian(a)
        penalty=np.vstack((-de*(e<.1005)[:,None],de*(e>9.999)[:,None]))
        return np.vstack((J.reshape(-1,count)/np.sqrt(len(x)),100*penalty))
    try:
        fit=least_squares(fun,start,bounds=(np.ones(count),np.full(count,3.)),jac=jac if analytic_jacobian else "2-point",diff_step=1e-4,max_nfev=150 if analytic_jacobian else 50,ftol=1e-8)
        reason=fit.message
    except Budget:reason=f'actual {call_budget}-call budget reached'
    if 'coefficients' not in best:raise RuntimeError('no energy-feasible candidate')
    result=cls(parent,best['coefficients'])
    out=Path(output);out.mkdir(parents=True,exist_ok=True);result.save(out/'candidate.json')
    report={'status':'training_only_not_validated','initial_artifact':initial,'force':asdict(parent.force),
        'azimuthal_only':azimuthal_only,'quintic':quintic,'grid_resolution':grid_resolution if dense_cylindrical else None,'analytic_parameter_jacobian':analytic_jacobian,'dense_cylindrical':dense_cylindrical,'coefficient_initial_artifact':coefficient_initial or ('artifacts/constrained/localized_swirl/candidate.json' if axial else None),'axial':axial,'localized':localized,'fourth_power':fourth_power,'tail_adaptive':tail_adaptive,'adaptive_seed':20261017 if tail_adaptive else None,'calls':calls,'call_budget':call_budget,'termination':str(reason),'seed':20260916,
        'cache_direct_max_discrepancy':cache_error,'training_points':len(x),'training_loss':best['loss'],'history':history,
        'coefficient_bounds':[-1,1],'basis':('twelve axial zero-moment rings' if (axial or quintic) else 'six localized zero-moment rings' if localized else 'five zero-moment polynomial modes')+(' times five quintic Bernstein modes' if quintic else ' times three cubic Bernstein temporal modes'),
        'initial_velocity_max_change':float(np.max(np.abs(result.velocity(x,.25)-parent.velocity(x,.25))))}
    (out/'training.json').write_text(json.dumps(report,indent=2)+'\n');print(report)

if __name__=='__main__':run()
