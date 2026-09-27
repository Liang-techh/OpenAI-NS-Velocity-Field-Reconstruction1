"""Single-scale constructive trial with an additional quadratic axial factor."""
import json
import numpy as np
from scipy.optimize import least_squares
from adaptive_bridge_recursive_defect import build_fields
from adaptive_bridge_moment_fit import moment_slices, outer_moments
from affine_momentum import combine, momentum
from midplane_resolved_feasibility import ZeroBackground, RADIAL_BREAKS, _integrated_moment_coefficients, _evaluate, _jacobian, _cone_replay
from separated_moment_modes import SeparatedMomentModes, RADIAL_WINDOWS_THREE
from radial_continuation import ROOT


def run():
    inner,fields=build_fields();base=fields['two_sided_cone']
    old=np.array(json.loads((ROOT/'midplane_axial_cone_all_knots_repair.json').read_text())['amplitudes']).reshape(3,2,3,2)
    start=np.zeros((3,2,3,3));start[:,:,:,:2]=old;start=start.ravel()
    def field(a,background=base):
        return SeparatedMomentModes(background,a,windows=RADIAL_WINDOWS_THREE,knots=(11.,15.,19.),axial_powers=(0,1,2))
    current=field(start)
    units=[field(np.eye(54)[i],ZeroBackground(base)) for i in range(18)]
    data=moment_slices(inner,base,current,orders=(11.,),n=24,unit_fields=units,
        unit_fields_are_deltas=True,radial_breaks=RADIAL_BREAKS)[0]
    c=_integrated_moment_coefficients(data)
    probe=np.random.default_rng(73).normal(size=18)*.01
    np.testing.assert_allclose(_evaluate(c,probe),outer_moments(data,probe),rtol=1e-10,atol=1e-7)
    scale=np.maximum(abs(c[0]),1.)
    fit=least_squares(lambda x:_evaluate(c,x)/scale,np.zeros(18),
        jac=lambda x:_jacobian(c,x)/scale[:,None],bounds=(-40.,40.),max_nfev=500,
        ftol=1e-12,xtol=1e-12,gtol=1e-12,x_scale='jac')
    candidate=start.copy();candidate[:18]+=fit.x;f=field(candidate)
    print(json.dumps(dict(stage='fit',success=bool(fit.success),moments=_evaluate(c,fit.x).tolist(),max_delta=float(max(abs(fit.x))))),flush=True)
    replay=[]
    for n in (48,96):
        d=moment_slices(inner,base,f,orders=(11.,),n=n,unit_fields=[f],radial_breaks=RADIAL_BREAKS)[0]
        values=outer_moments(d,np.zeros(1))
        replay.append(dict(order=n,moments=values.tolist(),absolute_max=float(max(abs(values)))))
        print(json.dumps(replay[-1]),flush=True)
    before=momentum(data['baseline']);after=momentum(combine(data['baseline'],data['modes'],fit.x))
    supports=json.loads((ROOT/'midplane_physical_covariance_pairs.json').read_text())
    cone=_cone_replay(inner,f,supports,11,order=48)
    report=dict(k=11,axial_powers=[0,1,2],optimizer_success=bool(fit.success),message=fit.message,nfev=fit.nfev,
        correction=fit.x.tolist(),amplitudes=candidate.tolist(),fit_moments=_evaluate(c,fit.x).tolist(),replay=replay,
        sampled_momentum_max_before=float(max(np.linalg.norm(before,axis=1))),
        sampled_momentum_max_after=float(max(np.linalg.norm(after,axis=1))),cone=cone,
        scope='Single-scale sampled moment repair with fixed pressure. Momentum peaks are on radial quadrature nodes only, not a full-domain max or volume L2. Polynomial axial factors do not establish finite energy or temporal recursion.',
        accepted=False,scale_recursion_established=False)
    (ROOT/'midplane_quadratic_axial_trial.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(momentum_ratio=report['sampled_momentum_max_after']/report['sampled_momentum_max_before'],cone_pass=cone['pass_count'])),flush=True)

if __name__=='__main__':
    run()
