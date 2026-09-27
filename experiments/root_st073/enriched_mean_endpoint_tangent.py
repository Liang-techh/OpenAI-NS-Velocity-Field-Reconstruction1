"""Degree-3 mean/mode-2 tangent with matching and direct endpoint geometry."""
import hashlib
import argparse
import json
import os
from pathlib import Path
for name in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):
    os.environ[name] = '1'
import numpy as np
from scipy.optimize import LinearConstraint, minimize
from constrained_tangent_projection import quadratic_moment_target
from enriched_shape_tangent import _weighted_constraint_setup
from wave_higher_harmonic_tangent import _mode_block_columns
from wave_momentum_projection import residual_and_jacobian, _metric

ROOT = Path(__file__).resolve().parent


def embed_seed(control):
    control = np.asarray(control,float)
    if len(control)!=236:
        raise ValueError('Expected the degree-2 mean / degree-3 mode-2 seed')
    new = np.zeros(264)
    for family in range(4):
        for a in range(3):
            for b in range(3):
                new[family*16+a*4+b] = control[family*9+a*3+b]
    new[64:] = control[36:]
    return new


def run(refined_cache_path=None,output_path=None):
    from enriched_endpoint_shape_cache import EnrichedEndpointShape
    oracle = EnrichedEndpointShape()
    paths = dict(seed=ROOT/'enriched_shape_candidate_snapshot.json',
        mean_rows=ROOT/'enriched_mean_constraint_rows.json',
        moment=ROOT/'wave_dynamics_mean_compatibility.json',
        cone=ROOT/'wave_mean_cone_projection.json',
        endpoint=ROOT/'enriched_endpoint_shape_cache.json')
    raw = {key:path.read_bytes() for key,path in paths.items()}
    data = {key:json.loads(value) for key,value in raw.items()}
    if any(value.get('status')!='completed' for value in data.values()):
        raise ValueError('Wait for frozen completed sources')
    seed = data['seed']
    if not np.array_equal(np.asarray(oracle.coefficients_original),
                          np.asarray(seed['selected']['coefficients_original'])):
        raise ValueError('Endpoint initial wave differs from seed')
    rows = data['mean_rows']
    moment = data['moment']['reusable_moment_linearization']
    cone = data['cone'].get('reusable_cone_linearization',data['cone'])
    snapshot = json.loads((ROOT/'full_wave_frozen_cache.json').read_text())
    geometry = snapshot['inputs']['wave']
    carrier = np.array(geometry['carrier'])
    cache_path=Path(refined_cache_path) if refined_cache_path else ROOT/'wave_momentum_projection.npz'
    with np.load(cache_path,allow_pickle=False) as loaded:
        cache = {key:loaded[key] for key in loaded.files}
    points,weights = cache['points'],cache['weights']
    if refined_cache_path:
        if not np.array_equal(np.asarray(cache['coefficients_original']),
                              np.asarray(seed['selected']['coefficients_original'])):
            raise ValueError('Refined cache initial wave differs from seed')
        D=np.asarray(cache['tangent_design'])
    else:
        mode0 = _mode_block_columns(points,geometry['center'],geometry['widths'],0,3,np.zeros(2))[:,::2]
        mode2 = _mode_block_columns(points,geometry['center'],geometry['widths'],2,3,2*carrier)
        D = np.column_stack((mode0,cache['tangent_design'][:,36:108],mode2))
    E = np.column_stack((np.asarray(rows['moment_rows']),np.zeros((4,200))))
    C = np.column_stack((np.asarray(rows['cone_control_rows']),np.zeros((81,200))))
    if D.shape[1]!=264 or E.shape!=(4,264) or C.shape!=(81,264):
        raise ValueError('Unexpected enriched mean layout')
    c = np.array(seed['selected']['coefficients_original'])
    x = np.r_[c[:,0],c[:,1]]
    if refined_cache_path:
        residual=np.asarray(cache['residual'])
    else:
        residual,_ = residual_and_jacobian(cache,x)
    target,_ = quadratic_moment_target(moment['baseline_moments'],moment['wave_moment_forms'],x)
    setup = _weighted_constraint_setup(D,E,target,residual,weights)
    particular = setup['particular']
    mapping = setup['mapping']*setup['scale']
    z0 = setup['q0']/setup['scale']
    old_control = embed_seed(seed['selected']['tangent_coefficients'])
    # Recover objective coordinates through the orthonormal weighted design;
    # least-squares inversion of physical control coefficients is ill-scaled.
    initial = setup['Q'].T@(setup['weighted_design']@(old_control-particular))/setup['scale']
    cone_base = np.array(cone['cone_baseline'])+np.einsum('i,kij,j->k',x,np.array(cone['cone_wave_forms']),x)
    lower = np.array(cone['cone_lower'])
    matrix = C@mapping
    norms = np.linalg.norm(matrix,axis=1)
    responsive = norms>max(norms.max()*1e-12,1e-14)
    rhs = lower+1e-4*responsive-cone_base-C@particular
    if np.any(rhs[~responsive]>1e-7):
        raise ValueError('Fixed cone row failure')
    matrix_n = matrix[responsive]/norms[responsive,None]
    rhs_n = rhs[responsive]/norms[responsive]
    reference = np.asarray(oracle.reference,float)
    signs = np.array([-1.,1.,np.sign(reference[2])])
    ref_scale = np.maximum(abs(reference),1e-30)
    requested = np.array([1e-6,1e-6,1e-3])
    def endpoint(z):
        values,jac = oracle.evaluate(particular+mapping@z)
        return signs*(values-reference)/ref_scale-requested, signs[:,None]*(jac@mapping)/ref_scale[:,None]
    report = dict(status='initialized',accepted=False,pde_validated=False,
        scale_recursion_established=False,source='wave_stress_growth_codesign.json',control_count=264,
        degrees=dict(mode0=3,mode1=2,mode2=3),
        sources={key:dict(path=paths[key].name,sha256=hashlib.sha256(value).hexdigest()) for key,value in raw.items()},
        scope='Reference-time full momentum with integral moments, sampled cones and nonlinear fixed-cylinder endpoint geometry. No endpoint momentum or continuous trajectory acceptance.',
        momentum_cache=dict(path=cache_path.name,sha256=hashlib.sha256(cache_path.read_bytes()).hexdigest()),
        objective_rank=setup['keep_rank'],requested_signed_fractional_endpoint_change=requested.tolist(),
        mapped_seed_control_max_error=float(np.max(abs(particular+mapping@initial-old_control))),
        initial_endpoint_constraint_margin=endpoint(initial)[0].tolist(),
        initial_cone_normalized_min_margin=float(np.min(matrix_n@initial-rhs_n)))
    output = Path(output_path) if output_path else ROOT/'enriched_mean_endpoint_tangent.json'
    def save():
        output.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    if np.min(endpoint(initial)[0])<-1e-8 or np.min(matrix_n@initial-rhs_n)<-1e-8:
        report['status']='mapped_seed_infeasible_or_unresolved'
        save()
        return report
    fit = minimize(lambda z:.5*np.sum((z-z0)**2),initial,jac=lambda z:z-z0,
        method='SLSQP',constraints=[LinearConstraint(matrix_n,rhs_n,np.inf),
            dict(type='ineq',fun=lambda z:endpoint(z)[0],jac=lambda z:endpoint(z)[1])],
        options=dict(maxiter=100,ftol=1e-12))
    selected = initial
    if np.min(endpoint(fit.x)[0])>=-1e-10 and np.min(matrix_n@fit.x-rhs_n)>=-1e-9 and np.sum((fit.x-z0)**2)<np.sum((initial-z0)**2):
        selected=fit.x
    control=particular+mapping@selected
    corrected=residual+(D@control).reshape(-1,3)
    values,_=oracle.evaluate(control)
    margins=cone_base+C@control-lower
    report['optimizer']=dict(success=bool(fit.success),message=str(fit.message),iterations=int(fit.nit))
    report['selected']=dict(coefficients_original=seed['selected']['coefficients_original'],
        coefficients_whitened=seed['selected']['coefficients_whitened'],
        tangent_coefficients=control.tolist(),growth_lambda=seed['selected']['growth_lambda'],
        training_momentum=_metric(corrected,weights),
        assembled_moment_max_abs=float(np.max(abs(E@control-target))),
        assembled_cone_min_margin=float(np.min(margins)),
        assembled_cone_location_pass_count=int(np.sum(np.all(margins.reshape(-1,3)>=-1e-7,axis=1))),
        endpoint_observables=values.tolist(),
        signed_fractional_endpoint_change=(signs*(values-reference)/ref_scale).tolist(),
        legacy_180_replay_compatible=False)
    report['assembled_feasible']=bool(
        report['selected']['assembled_moment_max_abs']<1e-5 and
        np.min(margins)>=-1e-7 and np.min(endpoint(selected)[0])>=-1e-10)
    report['status']='completed'
    save()
    print(json.dumps({key:value for key,value in report['selected'].items() if 'coefficients' not in key}),flush=True)
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--refined-cache',type=Path)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    run(args.refined_cache,args.output)
