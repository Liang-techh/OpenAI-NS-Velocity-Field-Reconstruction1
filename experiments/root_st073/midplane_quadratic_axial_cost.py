"""Reduce full sampled momentum inside the enriched moment-feasible family."""
import json
import numpy as np
from scipy.optimize import minimize
from adaptive_bridge_recursive_defect import build_fields
from adaptive_bridge_moment_fit import moment_slices, outer_moments
from affine_momentum import combine, momentum, jets
from midplane_resolved_feasibility import ZeroBackground, RADIAL_BREAKS, _integrated_moment_coefficients, _evaluate, _jacobian, _cone_replay
from separated_moment_modes import SeparatedMomentModes, RADIAL_WINDOWS_THREE
from radial_continuation import ROOT


def run():
    inner,fields=build_fields();base=fields['two_sided_cone']
    old=np.array(json.loads((ROOT/'midplane_axial_cone_all_knots_repair.json').read_text())['amplitudes']).reshape(3,2,3,2)
    start=np.zeros((3,2,3,3));start[:,:,:,:2]=old;start=start.ravel()
    previous=json.loads((ROOT/'midplane_quadratic_axial_trial.json').read_text())
    def field(a,background=base):
        return SeparatedMomentModes(background,a,windows=RADIAL_WINDOWS_THREE,knots=(11.,15.,19.),axial_powers=(0,1,2))
    current=field(start)
    units=[field(np.eye(54)[i],ZeroBackground(base)) for i in range(18)]
    data=moment_slices(inner,base,current,orders=(11.,),n=48,unit_fields=units,
        unit_fields_are_deltas=True,radial_breaks=RADIAL_BREAKS)[0]
    c=_integrated_moment_coefficients(data)
    scale=np.maximum(abs(c[0]),1.)
    constraints=[dict(type='eq',fun=lambda x:_evaluate(c,x)/scale,
        jac=lambda x:_jacobian(c,x)/scale[:,None])]
    before=momentum(data['baseline']);rscale=float(max(np.linalg.norm(before,axis=1)))
    weight=np.concatenate([p['weights']*p['radii'] for p in data['panels']]);weight/=sum(weight)
    def objective(x):
        jet=combine(data['baseline'],data['modes'],x)
        r=momentum(jet)/rscale
        mu,mg,ml=data['modes']
        dr=(ml+np.einsum('pnij,nj->pni',mg,jet[0])+np.einsum('nij,pnj->pni',jet[1],mu))/rscale
        norm2=np.sum(r*r,axis=1)
        value=np.dot(weight,norm2+.1*norm2**2)
        grad=2*np.einsum('n,ni,pni->p',weight*(1+.2*norm2),r,dr)
        return float(value),grad
    x0=np.array(previous['correction'])
    direction=np.random.default_rng(73).normal(size=18);direction/=np.linalg.norm(direction)
    error=abs((objective(x0+1e-5*direction)[0]-objective(x0-1e-5*direction)[0])/2e-5-objective(x0)[1]@direction)
    assert error < 1e-5*max(1,abs(objective(x0)[1]@direction)),error
    normfit=minimize(lambda x:(float(x@x),2*x),x0,jac=True,method='SLSQP',bounds=[(-40.,40.)]*18,
        constraints=constraints,options=dict(maxiter=300,ftol=1e-11))
    fit=minimize(objective,normfit.x,jac=True,method='SLSQP',bounds=[(-40.,40.)]*18,
        constraints=constraints,options=dict(maxiter=400,ftol=1e-11))
    stages=[]
    for label,x,success in [('previous',x0,False),('minimum_norm',normfit.x,normfit.success),('momentum',fit.x,fit.success)]:
        rr=momentum(combine(data['baseline'],data['modes'],x))
        stages.append(dict(label=label,success=bool(success),moment_max=float(max(abs(_evaluate(c,x)))),
            correction_norm=float(np.linalg.norm(x)),momentum_peak=float(max(np.linalg.norm(rr,axis=1))),objective=objective(x)[0]))
    print(json.dumps(stages),flush=True)
    candidate=start.copy();candidate[:18]+=fit.x;f=field(candidate)
    replay=[]
    for n in (96,128):
        d=moment_slices(inner,base,f,orders=(11.,),n=n,unit_fields=[f],radial_breaks=RADIAL_BREAKS)[0]
        values=outer_moments(d,np.zeros(1));rr=momentum(d['baseline'])
        replay.append(dict(order=n,moments=values.tolist(),absolute_max=float(max(abs(values))),momentum_peak=float(max(np.linalg.norm(rr,axis=1)))))
        print(json.dumps(replay[-1]),flush=True)
    tau=.5*2**-11
    points=np.concatenate([inner.from_similarity(inner.p.X_max*(1+15*np.linspace(.01,.99,61))**2,np.full(61,eta),tau) for eta in (-.3,-.1,0.,.1,.3)])
    args=(points,tau,.0005*np.sqrt(inner.nu*tau),.0001*tau)
    heldout={}
    for label,fld in [('baseline',current),('candidate',f)]:
        rr=momentum(jets(fld,*args));heldout[label]=float(max(np.linalg.norm(rr,axis=1)))
    supports=json.loads((ROOT/'midplane_physical_covariance_pairs.json').read_text())
    cone=_cone_replay(inner,f,supports,11,order=96)
    report=dict(k=11,fit_order=48,stages=stages,optimizer_message=str(fit.message),correction=fit.x.tolist(),amplitudes=candidate.tolist(),
        gradient_directional_error=float(error),replay=replay,heldout_momentum_peaks=heldout,cone=cone,
        scope='Single-scale moment-constrained optimization of sampled full momentum; heldout radial/axial grid. No global maximum, spatial-volume L2, finite energy, or scale recursion established.',accepted=False,scale_recursion_established=False)
    (ROOT/'midplane_quadratic_axial_cost.json').write_bytes((json.dumps(report,indent=2)+'\n').encode())
    print(json.dumps(dict(heldout=heldout,cone_pass=cone['pass_count'])),flush=True)

if __name__=='__main__':
    run()
