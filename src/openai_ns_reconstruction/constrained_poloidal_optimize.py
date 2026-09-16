"""Joint bounded poloidal/pressure fit with analytic parameter derivatives."""
from dataclasses import asdict,replace
import json
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from .constrained_local_pressure import LocalPressureCandidate
from .constrained_poloidal import PoloidalCandidate
from .constrained_temporal_cache import cached_residual,cached_energy
from .constrained_optimize import training_residual


def run(output='artifacts/constrained/poloidal_joint',anchor_core=False,coupled=False):
    initial='artifacts/constrained/poloidal_anchor/candidate.json' if coupled else 'artifacts/constrained/local_pressure/candidate.json'
    if coupled:
        from .constrained_coupled import CoupledCandidate
        old=PoloidalCandidate.load(initial);cls=CoupledCandidate;anchor_core=old.anchor_core
        base=replace(old.base,coefficients=(0.,)*90)
        start_velocity=np.r_[old.base.coefficients,old.coefficients]
    else:
        old=LocalPressureCandidate.load(initial);cls=PoloidalCandidate;base=old.base;start_velocity=np.zeros(27)
    count=cls.COEFFICIENT_COUNT;total=count+27
    zero=cls(base,old.pressure_coefficients,anchor_core=anchor_core)
    rng=np.random.default_rng(20260916);x=rng.uniform(-2,2,(2048,3));t=rng.uniform(.252,.748,len(x))
    rr,zz,tt=np.meshgrid(np.linspace(.03,1.97,48),np.linspace(0,1.97,48),np.linspace(.252,.748,9),indexing='ij')
    x=np.vstack((x,np.column_stack((rr.ravel(),np.zeros(rr.size),zz.ravel()))));t=np.r_[t,tt.ravel()]
    evaluate=cached_residual(zero,x,t);energy=cached_energy(zero,np.linspace(.25,.75,9))
    G=np.empty((len(x),3,27));h=.001
    for j in range(3):
        d=np.eye(3)[j]*h;G[:,j,:]=(zero.pressure_basis(x+d,t)-zero.pressure_basis(x-d,t))/(2*h)
    p0=np.asarray(old.pressure_coefficients)
    probe=np.linspace(-.03,.03,count)
    direct=training_residual(cls(base,old.pressure_coefficients,tuple(probe),anchor_core),old.force,x,t,.01)
    discrepancy=float(np.max(np.abs(evaluate(probe)-direct)))
    if discrepancy>1e-7:raise RuntimeError('poloidal residual cache mismatch')
    calls=0;best={'loss':float('inf')};history=[]
    def raw(v):return evaluate(v[:count])+G@(v[count:]-p0)
    def fun(v):
        nonlocal calls
        calls+=1;R=raw(v);norm=np.linalg.norm(R,axis=1)
        e=energy(v[:count]);violation=np.r_[np.maximum(.1005-e,0),np.maximum(e-9.999,0)]
        loss=float(np.mean(norm**4))
        if not np.any(violation) and loss<best['loss']:
            best.update(loss=loss,v=v.copy());history.append({'call':calls,'loss':loss})
        return np.r_[(R*norm[:,None]).ravel()/np.sqrt(len(x)),100*violation]
    def jac(v):
        R=raw(v);J=np.concatenate((evaluate.jacobian(v[:count]),G),axis=-1)
        norm=np.linalg.norm(R,axis=1);inverse=np.divide(1.,norm,out=np.zeros_like(norm),where=norm>0)
        projected=np.einsum('nc,nck->nk',R,J)
        J=norm[:,None,None]*J+R[:,:,None]*projected[:,None,:]*inverse[:,None,None]
        e=energy(v[:count]);de=np.column_stack((energy.jacobian(v[:count]),np.zeros((9,27))))
        penalty=np.vstack((-de*(e<.1005)[:,None],de*(e>9.999)[:,None]))
        return np.vstack((J.reshape(-1,total)/np.sqrt(len(x)),100*penalty))
    fit=least_squares(fun,np.r_[start_velocity,p0],jac=jac,bounds=(-np.ones(total),np.ones(total)),max_nfev=80,ftol=1e-8)
    if 'v' not in best:raise RuntimeError('no energy-feasible iterate')
    v=best['v'];candidate=cls(base,tuple(v[count:]),tuple(v[:count]),anchor_core)
    out=Path(output);out.mkdir(parents=True,exist_ok=True);candidate.save(out/'candidate.json')
    report={'coupled':coupled,'anchor_core':anchor_core,'status':'joint_fit_not_validated','force':asdict(old.force),'calls':calls,'max_nfev':80,'termination':fit.message,
            'initial_artifact':initial,'seed':20260916,'grid':[48,48,9],
            'training_points':len(x),'loss':'mean fourth power full residual','training_loss':best['loss'],'history':history,
            'coefficient_bounds':[-1,1],'velocity_coefficients':count,'pressure_coefficients':27,'cache_direct_discrepancy':discrepancy}
    (out/'training.json').write_text(json.dumps(report,indent=2)+'\n');print(report)

if __name__=='__main__':run()
