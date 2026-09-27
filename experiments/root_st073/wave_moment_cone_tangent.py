"""Fit mode-0 time/pressure controls subject to moments and sampled cones.

Wave shape and nonzero-harmonic controls remain fixed at the saved moment
candidate. A phase-I LP tests this bounded basis before the convex QP.
No spatial/temporal acceptance follows from assembled constraints alone.
"""
import hashlib
import json
import os
from pathlib import Path

for _name in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):
    os.environ[_name] = '1'

import numpy as np
from scipy.optimize import LinearConstraint,linprog,minimize

from constrained_tangent_projection import ConstrainedTangent,quadratic_moment_target
from wave_momentum_projection import residual_and_jacobian,_metric

ROOT = Path(__file__).resolve().parent


def run():
    paths = {name:ROOT/filename for name,filename in dict(
        seed='wave_dynamics_moment_codesign.json',
        moments='wave_dynamics_mean_compatibility.json',
        cones='wave_mean_cone_projection.json').items()}
    raw = {name:path.read_bytes() for name,path in paths.items()}
    sources = {name:json.loads(value) for name,value in raw.items()}
    if any(report.get('status','completed') != 'completed' for report in sources.values()):
        raise ValueError('Wait for completed source reports before solving')
    seed = sources['seed']
    moment = sources['moments']['reusable_moment_linearization']
    cones = sources['cones']
    if 'cone_control_rows' not in cones:
        cones = cones['reusable_cone_linearization']
    with np.load(ROOT/'wave_momentum_projection.npz',allow_pickle=False) as data:
        cache = {key:data[key] for key in data.files}
    packed = np.array(seed['selected']['coefficients_original'])
    x = np.r_[packed[:,0],packed[:,1]]
    residual,_ = residual_and_jacobian(cache,x)
    weights = cache['weights']
    sqrtw = np.repeat(np.sqrt(weights),3)
    design = cache['tangent_design']
    weighted = design*sqrtw[:,None]
    control = np.array(seed['selected']['tangent_coefficients'])
    E = np.array(moment['moment_rows'])
    C = np.array(cones['cone_control_rows'])
    if max(np.max(abs(E[:,36:])),np.max(abs(C[:,36:]))) > 1e-12:
        raise ValueError('This restricted solve expects only mode0 constraint responses')
    target,_ = quadratic_moment_target(moment['baseline_moments'],moment['wave_moment_forms'],x)
    projector = ConstrainedTangent(weighted[:,:36],E[:,:36])
    particular = (projector.particular@target)/projector.scales
    reduced = projector.design@projector.nullspace
    u,s,vh = np.linalg.svd(reduced,full_matrices=False)
    keep = s > 1e-10*s[0]
    Q = u[:,keep]
    mapping = ((projector.nullspace@vh[keep].T)/s[keep])/projector.scales[:,None]
    shifted = residual.reshape(-1)*sqrtw+weighted[:,36:]@control[36:]+weighted[:,:36]@particular
    q0 = -Q.T@shifted
    scale = max(float(np.linalg.norm(q0)),float(np.linalg.norm(shifted)),1.)
    z0 = q0/scale
    forms = np.array(cones['cone_wave_forms'])
    wave_rows = np.einsum('i,kij,j->k',x,forms,x)
    baseline = np.array(cones['cone_baseline'])
    lower = np.array(cones['cone_lower'])
    safety_margin = 1e-4
    matrix = C[:,:36]@mapping*scale
    row_norms = np.linalg.norm(matrix,axis=1)
    responsive = row_norms > max(row_norms.max()*1e-10,1e-14)
    rhs = lower+safety_margin*responsive-baseline-wave_rows-C[:,:36]@particular
    if np.any(rhs[~responsive] > 1e-7):
        report = dict(status='restricted_basis_fixed_row_failure',accepted=False,
            pde_validated=False,scale_recursion_established=False,
            reason='Moment-preserving retained basis cannot change an already violated cone row',
            fixed_rows=np.flatnonzero(~responsive).tolist(),
            fixed_row_violations=rhs[~responsive].tolist(),
            responsiveness_cutoff=max(row_norms.max()*1e-10,1e-14),
            sources={name:dict(path=paths[name].name,sha256=hashlib.sha256(raw[name]).hexdigest()) for name in paths})
        (ROOT/'wave_moment_cone_tangent.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
        print(json.dumps(report),flush=True)
        return report
    normalized = matrix[responsive]/row_norms[responsive,None]
    normalized_rhs = rhs[responsive]/row_norms[responsive]
    phase = linprog(np.r_[np.zeros(len(z0)),1.],
        A_ub=np.column_stack((-normalized,-np.ones(len(normalized_rhs)))),
        b_ub=-normalized_rhs,bounds=[(None,None)]*len(z0)+[(0,None)],method='highs')
    report = dict(status='phase_one_complete',accepted=False,pde_validated=False,
        constraints_maintained=False,scale_recursion_established=False,
        source=seed['source'],sources={name:dict(path=paths[name].name,
          sha256=hashlib.sha256(raw[name]).hexdigest()) for name in paths},
        mode0_control_count=36,objective_rank=int(keep.sum()),
        moment_rank=projector.metadata['moment_rank'],cone_count=len(lower),
        safety_margin=safety_margin,
        responsive_cone_rows=np.flatnonzero(responsive).tolist(),
        fixed_cone_rows=np.flatnonzero(~responsive).tolist(),
        safety_margin_scope='Applied only to responsive rows; fixed satisfied rows retain the original target.',
        phase_one=dict(success=bool(phase.success),message=phase.message,
            slack=float(phase.x[-1]) if phase.x is not None else None),
        scope='Fixed wave and nonzero-harmonic controls; only mode0 tangent controls fitted. Assembled moments and cone inequalities; no actual-field replay or time/recursion acceptance.')
    output = ROOT/'wave_moment_cone_tangent.json'

    def save():
        output.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')

    save()
    if not phase.success or phase.x[-1] > 1e-8:
        report['status'] = 'restricted_basis_infeasible_or_unresolved'
        save()
        print(json.dumps(report['phase_one']),flush=True)
        return report
    feasible = phase.x[:-1]
    fit = minimize(lambda z:.5*np.sum((z-z0)**2),feasible,
        jac=lambda z:z-z0,method='SLSQP',
        constraints=[LinearConstraint(normalized,normalized_rhs,np.inf)],
        options=dict(maxiter=500,ftol=1e-12,disp=False))
    selected = feasible
    if min(matrix@fit.x-rhs) >= -1e-6 and np.sum((fit.x-z0)**2) <= np.sum((feasible-z0)**2):
        selected = fit.x
    control[:36] = particular+mapping@(scale*selected)
    corrected = residual+(design@control).reshape(-1,3)
    margins = baseline+wave_rows+C@control-lower
    report['optimizer'] = dict(success=bool(fit.success),status=int(fit.status),message=str(fit.message))
    report['selected'] = dict(coefficients_original=seed['selected']['coefficients_original'],
        coefficients_whitened=seed['selected']['coefficients_whitened'],
        tangent_coefficients=control.tolist(),growth_lambda=seed['selected']['growth_lambda'],
        training_momentum=_metric(corrected,weights),
        assembled_moment_max_abs=float(max(abs(E@control-target))),
        assembled_cone_min_margin=float(min(margins)),
        assembled_cone_location_pass_count=int(np.sum(np.all(margins.reshape(-1,3)>=-1e-7,axis=1))),
        assembled_cone_inequality_pass_count=int(np.sum(margins>=-1e-7)))
    report['status'] = 'completed'
    save()
    print(json.dumps({key:value for key,value in report['selected'].items() if 'coefficients' not in key}),flush=True)
    return report


if __name__ == '__main__':
    run()
