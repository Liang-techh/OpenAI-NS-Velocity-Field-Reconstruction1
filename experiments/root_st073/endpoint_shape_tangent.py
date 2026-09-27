"""Optimize full tangent momentum with direct finite-endpoint geometry bounds."""
import hashlib
import json
import os
from pathlib import Path
for name in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):
    os.environ[name] = '1'
import numpy as np
from scipy.optimize import LinearConstraint, minimize
from constrained_tangent_projection import ConstrainedTangent, quadratic_moment_target
from wave_momentum_projection import residual_and_jacobian, _metric

ROOT = Path(__file__).resolve().parent


def run():
    from wave_endpoint_shape_cache import EndpointShape
    oracle = EndpointShape()
    paths = dict(seed=ROOT/'shape_direction_margin_tangent.json',
        moments=ROOT/'wave_dynamics_mean_compatibility.json',
        cones=ROOT/'wave_mean_cone_projection.json',
        endpoint=ROOT/'wave_endpoint_shape_cache.json')
    raw = {key:path.read_bytes() for key,path in paths.items()}
    data = {key:json.loads(value) for key,value in raw.items()}
    if any(value.get('status') != 'completed' for value in data.values()):
        raise ValueError('Wait for completed frozen inputs')
    seed = data['seed']
    if not np.array_equal(np.asarray(oracle.coefficients_original),
                          np.asarray(seed['selected']['coefficients_original'])):
        raise ValueError('Endpoint oracle initial wave differs from candidate')
    moment = data['moments']['reusable_moment_linearization']
    cone = data['cones'].get('reusable_cone_linearization',data['cones'])
    with np.load(ROOT/'wave_momentum_projection.npz',allow_pickle=False) as loaded:
        cache = {key:loaded[key] for key in loaded.files}
    c = np.array(seed['selected']['coefficients_original'])
    x = np.r_[c[:,0],c[:,1]]
    residual,_ = residual_and_jacobian(cache,x)
    sqrtw = np.repeat(np.sqrt(cache['weights']),3)
    D = cache['tangent_design']*sqrtw[:,None]
    E = np.array(moment['moment_rows'])
    target,_ = quadratic_moment_target(moment['baseline_moments'],moment['wave_moment_forms'],x)
    projector = ConstrainedTangent(D,E)
    particular = (projector.particular@target)/projector.scales
    u,s,vh = np.linalg.svd(projector.design@projector.nullspace,full_matrices=False)
    keep = s>1e-10*s[0]
    mapping = ((projector.nullspace@vh[keep].T)/s[keep])/projector.scales[:,None]
    shifted = residual.reshape(-1)*sqrtw+D@particular
    scale = max(np.linalg.norm(shifted),1.)
    mapping *= scale
    z0 = -u[:,keep].T@shifted/scale
    initial = np.linalg.lstsq(mapping,np.array(seed['selected']['tangent_coefficients'])-particular,rcond=1e-12)[0]
    C = np.array(cone['cone_control_rows'])
    cone_base = np.array(cone['cone_baseline'])+np.einsum('i,kij,j->k',x,np.array(cone['cone_wave_forms']),x)
    cone_lower = np.array(cone['cone_lower'])
    cone_matrix = C@mapping
    norms = np.linalg.norm(cone_matrix,axis=1)
    responsive = norms>max(norms.max()*1e-12,1e-14)
    rhs = cone_lower+1e-4*responsive-cone_base-C@particular
    if np.any(rhs[~responsive]>1e-7):
        raise ValueError('A fixed cone row is infeasible')
    matrix_n = cone_matrix[responsive]/norms[responsive,None]
    rhs_n = rhs[responsive]/norms[responsive]
    reference = np.asarray(oracle.reference,float)
    signs = np.array([-1.,1.,np.sign(reference[2])])
    reference_scale = np.maximum(abs(reference),1e-30)
    # Direct fractional endpoint improvement, not instantaneous rate per k.
    requested_change = np.array([1e-6,1e-6,1e-3])
    def geometry(z):
        values,jac = oracle.evaluate(particular+mapping@z)
        change = signs*(values-reference)/reference_scale
        derivative = signs[:,None]*(jac@mapping)/reference_scale[:,None]
        return change-requested_change, derivative
    initial_geometry,_ = geometry(initial)
    if np.min(initial_geometry)<-1e-8 or np.min(matrix_n@initial-rhs_n)<-1e-8:
        raise ValueError('Expected the known direction-margin seed to be endpoint feasible')
    fit = minimize(lambda z:.5*np.sum((z-z0)**2),initial,
        jac=lambda z:z-z0,method='SLSQP',
        constraints=[LinearConstraint(matrix_n,rhs_n,np.inf),
            dict(type='ineq',fun=lambda z:geometry(z)[0],jac=lambda z:geometry(z)[1])],
        options=dict(maxiter=100,ftol=1e-12))
    selected = initial
    if np.min(geometry(fit.x)[0])>=-1e-10 and np.min(matrix_n@fit.x-rhs_n)>=-1e-9 and np.sum((fit.x-z0)**2)<np.sum((initial-z0)**2):
        selected = fit.x
    control = particular+mapping@selected
    corrected = residual+(cache['tangent_design']@control).reshape(-1,3)
    values,_ = oracle.evaluate(control)
    margins = cone_base+C@control-cone_lower
    report = dict(status='completed',accepted=False,pde_validated=False,
        scale_recursion_established=False,source=seed['source'],
        sources={key:dict(path=paths[key].name,sha256=hashlib.sha256(value).hexdigest()) for key,value in raw.items()},
        scope='Instantaneous momentum fitted with sampled finite-endpoint geometry constraints. Endpoint velocity is an affine tangent candidate; endpoint momentum, interval matching and NS evolution are not certified.',
        requested_signed_fractional_endpoint_change=requested_change.tolist(),
        initial_endpoint_constraint_margin=initial_geometry.tolist(),
        optimizer=dict(success=bool(fit.success),message=str(fit.message),iterations=int(fit.nit)),
        selected=dict(coefficients_original=seed['selected']['coefficients_original'],
            coefficients_whitened=seed['selected']['coefficients_whitened'],
            tangent_coefficients=control.tolist(),growth_lambda=seed['selected']['growth_lambda'],
            training_momentum=_metric(corrected,cache['weights']),
            assembled_moment_max_abs=float(np.max(abs(E@control-target))),
            assembled_cone_min_margin=float(np.min(margins)),
            assembled_cone_location_pass_count=int(np.sum(np.all(margins.reshape(-1,3)>=-1e-7,axis=1))),
            endpoint_observables=values.tolist(),
            signed_fractional_endpoint_change=(signs*(values-reference)/reference_scale).tolist()))
    (ROOT/'endpoint_shape_tangent.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({key:value for key,value in report['selected'].items() if 'coefficients' not in key}),flush=True)
    return report


if __name__=='__main__':
    run()
