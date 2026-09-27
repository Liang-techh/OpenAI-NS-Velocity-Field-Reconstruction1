"""Full tangent QP preserving moments, cones and local shape-rate directions.

Shape rows are sampled instantaneous derivatives, not a time trajectory.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path

for name in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):
    os.environ[name] = '1'
import numpy as np
from scipy.optimize import LinearConstraint, linprog, minimize
from constrained_tangent_projection import ConstrainedTangent, quadratic_moment_target
from wave_momentum_projection import residual_and_jacobian, _metric

ROOT = Path(__file__).resolve().parent


def run(relative_rate=0.01, output_path=None):
    paths = dict(seed=ROOT/'wave_moment_cone_tangent.json',
        shape=ROOT/'wave_shape_tangent_rows.json',
        moments=ROOT/'wave_dynamics_mean_compatibility.json',
        cones=ROOT/'wave_mean_cone_projection.json')
    raw = {key:path.read_bytes() for key,path in paths.items()}
    data = {key:json.loads(value) for key,value in raw.items()}
    if any(value.get('status') != 'completed' for value in data.values()):
        raise ValueError('Wait for frozen completed source reports')
    seed, shape = data['seed'], data['shape']
    moments = data['moments']['reusable_moment_linearization']
    cones = data['cones'].get('reusable_cone_linearization', data['cones'])
    if shape['source_sha256'] != hashlib.sha256(raw['seed']).hexdigest():
        raise ValueError('Shape rows must describe this exact frozen seed')
    with np.load(ROOT/'wave_momentum_projection.npz',allow_pickle=False) as loaded:
        cache = {key:loaded[key] for key in loaded.files}
    c = np.array(seed['selected']['coefficients_original'])
    x = np.r_[c[:,0],c[:,1]]
    residual, _ = residual_and_jacobian(cache,x)
    sqrtw = np.repeat(np.sqrt(cache['weights']),3)
    D = cache['tangent_design']*sqrtw[:,None]
    E = np.array(moments['moment_rows'])
    target, _ = quadratic_moment_target(moments['baseline_moments'],moments['wave_moment_forms'],x)
    projector = ConstrainedTangent(D,E)
    particular = (projector.particular@target)/projector.scales
    reduced = projector.design@projector.nullspace
    u,s,vh = np.linalg.svd(reduced,full_matrices=False)
    keep = s > 1e-10*s[0]
    Q = u[:,keep]
    mapping = ((projector.nullspace@vh[keep].T)/s[keep])/projector.scales[:,None]
    shifted = residual.reshape(-1)*sqrtw+D@particular
    q0 = -Q.T@shifted
    scale = max(float(np.linalg.norm(shifted)),1.)
    z0 = q0/scale
    C = np.array(cones['cone_control_rows'])
    cone_base = np.array(cones['cone_baseline'])+np.einsum(
        'i,kij,j->k',x,np.array(cones['cone_wave_forms']),x)
    cone_lower = np.array(cones['cone_lower'])
    S = np.array(shape['shape_control_rows'])
    shape_base = np.array(shape['shape_baseline'])
    signs = np.array(shape['desired_signs'])
    reference_data = shape['reference_observables']
    if isinstance(reference_data,dict):
        reference_data = [reference_data[key] for key in
            ('enstrophy_radial_rms','enstrophy_aspect_ratio','enstrophy_weighted_angular_speed')]
    reference = np.array(reference_data)
    shape_scale = np.maximum(abs(reference),1e-30)
    # Normalize geometric rates to fractional change per k. A small positive
    # margin separates a desired direction from numerical equality at zero.
    C_all = np.vstack((C,signs[:,None]*S/shape_scale[:,None]))
    base_all = np.r_[cone_base,signs*shape_base/shape_scale]
    lower_all = np.r_[cone_lower,np.full(3,relative_rate)]
    matrix = C_all@mapping*scale
    norms = np.linalg.norm(matrix,axis=1)
    responsive = norms > max(norms.max()*1e-12,1e-14)
    safety = np.r_[np.full(len(C),1e-4),np.zeros(3)]*responsive
    rhs = lower_all+safety-base_all-C_all@particular
    report = dict(status='assembling',accepted=False,pde_validated=False,
        scale_recursion_established=False,source=seed['source'],
        sources={key:dict(path=paths[key].name,sha256=hashlib.sha256(value).hexdigest()) for key,value in raw.items()},
        relative_shape_rate_floor=relative_rate, objective_rank=int(keep.sum()),
        scope='Fixed wave; all 180 time/pressure controls. Four integral moments, 81 sampled cones and three instantaneous sampled shape directions. No finite-time, full-domain or PDE acceptance.')
    output = Path(output_path) if output_path else ROOT/'shape_constrained_wave_tangent.json'
    def save():
        output.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    if np.any(rhs[~responsive] > 1e-7):
        report.update(status='fixed_row_failure',fixed_row_rhs=rhs[~responsive].tolist())
        save()
        return report
    matrix_n = matrix[responsive]/norms[responsive,None]
    rhs_n = rhs[responsive]/norms[responsive]
    phase = linprog(np.r_[np.zeros(len(z0)),1.],
        A_ub=np.column_stack((-matrix_n,-np.ones(len(rhs_n)))),
        b_ub=-rhs_n,bounds=[(None,None)]*len(z0)+[(0,None)],method='highs')
    report['phase_one'] = dict(success=bool(phase.success),message=phase.message,
        slack=float(phase.x[-1]) if phase.x is not None else None)
    if not phase.success or phase.x[-1] > 1e-8:
        report.update(status='restricted_basis_infeasible_or_unresolved')
        save()
        return report
    fit = minimize(lambda z:.5*np.sum((z-z0)**2),phase.x[:-1],
        jac=lambda z:z-z0,method='SLSQP',
        constraints=[LinearConstraint(matrix_n,rhs_n,np.inf)],
        options=dict(maxiter=500,ftol=1e-12))
    z = phase.x[:-1]
    if np.min(matrix_n@fit.x-rhs_n) >= -1e-9 and np.sum((fit.x-z0)**2) <= np.sum((z-z0)**2):
        z = fit.x
    control = particular+mapping@(scale*z)
    corrected = residual+(cache['tangent_design']@control).reshape(-1,3)
    cone_margins = cone_base+C@control-cone_lower
    shape_rates = shape_base+S@control
    report['optimizer'] = dict(success=bool(fit.success),message=str(fit.message))
    report['selected'] = dict(coefficients_original=seed['selected']['coefficients_original'],
        coefficients_whitened=seed['selected']['coefficients_whitened'],
        tangent_coefficients=control.tolist(),growth_lambda=seed['selected']['growth_lambda'],
        training_momentum=_metric(corrected,cache['weights']),
        assembled_moment_max_abs=float(np.max(abs(E@control-target))),
        assembled_cone_min_margin=float(np.min(cone_margins)),
        assembled_cone_location_pass_count=int(np.sum(np.all(cone_margins.reshape(-1,3)>=-1e-7,axis=1))),
        assembled_shape_rates=shape_rates.tolist(),
        signed_relative_shape_rates=(signs*shape_rates/shape_scale).tolist())
    report['assembled_feasible'] = bool(np.max(abs(E@control-target))<1e-5 and
        np.min(cone_margins)>=-1e-7 and np.min(signs*shape_rates/shape_scale-relative_rate)>=-1e-7)
    report['status'] = 'completed'
    save()
    print(json.dumps({k:v for k,v in report['selected'].items() if 'coefficients' not in k}),flush=True)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--relative-rate',type=float,default=.01)
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    run(args.relative_rate,args.output)
