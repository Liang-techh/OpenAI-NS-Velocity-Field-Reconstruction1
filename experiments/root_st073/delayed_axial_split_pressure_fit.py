"""Joint 15-node cone and momentum screen after axial-split swirl repair.

The two split modes are supported only for X>1.02 and in disjoint eta
bands. Previously computed stress primitives remain valid at the other
thirteen nodes; the cached exact quadratic response updates the two new
outer nodes. Full velocity momentum is then recomputed at all centers
and a nearby-time grid before pressure fitting.
"""

import json

import numpy as np

from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from delayed_axial_split_swirl_response import (AXIAL_INTERVALS,
                                                  CACHE_NAME as RESPONSE_CACHE,
                                                  RADIAL, load_rows)
from delayed_midband_moment_patch import SimilaritySwirlMode, mode_profile
from delayed_momentum_tangent_screen import nodes, stats
from delayed_multimode_cone_fit import BASE_NAME, load_cache
from delayed_multimode_pressure_admission import geometry as response_geometry
from delayed_outer_radial_extension import (CACHE_NAME as OUTER_CACHE,
                                             make_extended_modes, signature)
from delayed_plateau_pressure_screen import (PlateauPressurePatch,
                                              minmax_pressure,
                                              plateau_stress_columns)
from delayed_pressure_cone_geometry import line_interval
from delayed_pressure_cone_window import pressure_stress_columns
from delayed_shifted_swirl_cone_fit import (AMPLITUDE,
                                             NEW_AXIAL, NEW_RADIAL,
                                             geometry as shear_geometry)
from delayed_similarity_curl_screen import CurlPatchedLift
from delayed_similarity_pressure_screen import SimilarityPressurePatch
from joint_collar_fit import kinematics
from moment_shear_slice_repair import moment_vector
from radial_continuation import ROOT


def run(tuned=False):
    mean_source = json.loads((ROOT/'delayed_outer_admission_search.json').read_text())
    shifted = json.loads((ROOT/'delayed_shifted_swirl_cone_fit.json').read_text())
    split = json.loads((ROOT/'delayed_axial_split_swirl_response.json').read_text())
    split_rows = load_rows()
    if tuned:
        tuning = json.loads((ROOT/'delayed_axial_split_amplitude_lp.json').read_text())
        amplitudes = tuning['best']['amplitudes']
    else:
        amplitudes = [row['selected']['amplitude'] for row in split['rows']]
    base = CoupledMomentPhysicalLift(slice_filename=BASE_NAME)
    broad = SimilaritySwirlMode(base, NEW_RADIAL,
                                eta_interval=NEW_AXIAL)
    split_modes = [SimilaritySwirlMode(base, RADIAL,
                                       eta_interval=axial)
                   for axial in AXIAL_INTERVALS]
    mean = CurlPatchedLift(
        base, make_extended_modes(base)+[broad]+split_modes,
        [*mean_source['coefficients'], AMPLITUDE, *amplitudes])
    old_pressure = SimilarityPressurePatch(base, np.zeros(6))
    plateau = PlateauPressurePatch(base)
    tau = shifted['tau']
    geometry_rows = shifted['geometry']
    X = np.array([row['X'] for row in geometry_rows])
    eta = np.array([row['eta'] for row in geometry_rows])
    points = base.compact.joined.inner.from_similarity(X, eta, tau)
    hs, ht = .0005*np.sqrt(base.nu*tau), .0001*tau
    u, grad, part = kinematics(mean, points, tau, hs, ht)
    train_residual = part+np.einsum('pab,pb->pa', grad, u)
    target_by_eta = {}
    for i, row in enumerate(split_rows):
        c = np.zeros(2)
        c[i] = amplitudes[i]
        target_by_eta[row['eta']] = response_geometry(row, c)[0]
    A, rhs, pointwise = [], [], []
    for i, row in enumerate(geometry_rows):
        target = (target_by_eta[row['eta']]
                  if row['X'] == 1.03 and row['eta'] in target_by_eta
                  else np.asarray(row['target']))
        N, K, lam2, multiplier = shear_geometry(
            u[i], grad[i], points[i, 0])
        if lam2 <= 0 or multiplier is None:
            raise ValueError(f'Invalid lambda at {row["X"]},{row["eta"]}')
        columns = np.column_stack((
            pressure_stress_columns(old_pressure, points[i], tau),
            plateau_stress_columns(plateau, points[i], tau)))
        directions = [.8*N+sign*multiplier*K
                      for sign in (-1., 1.)]
        for direction in directions:
            A.append(direction@columns)
            rhs.append(-.1-direction@target)
        pointwise.append(dict(X=row['X'], eta=row['eta'],
                              lambda_squared=lam2,
                              target=target.tolist(),
                              free_axial=bool(line_interval(
                                  directions, target, .1, 1) is not None),
                              residual_norm=float(np.linalg.norm(
                                  train_residual[i]))))
    A, rhs = np.asarray(A), np.asarray(rhs)
    nearby_tau = .5*2**(-5.4)
    nearby_points, _, _ = nodes(base, (1.01, 1.015, 1.025),
                                 (.22, .28, .32), nearby_tau)
    nearby_hs, nearby_ht = (.001*np.sqrt(base.nu*nearby_tau),
                            .00025*nearby_tau)
    nearby_u, nearby_grad, nearby_part = kinematics(
        mean, nearby_points, nearby_tau, nearby_hs, nearby_ht)
    nearby_residual = (nearby_part+
                       np.einsum('pab,pb->pa', nearby_grad, nearby_u))
    _, old_train = old_pressure.basis(points, tau)
    _, new_train = plateau.basis(points, tau)
    _, old_near = old_pressure.basis(nearby_points, nearby_tau)
    _, new_near = plateau.basis(nearby_points, nearby_tau)
    train_gradient = np.concatenate((old_train, new_train), axis=2)
    near_gradient = np.concatenate((old_near, new_near), axis=2)
    matrix = np.concatenate((train_gradient.reshape(-1, 18),
                             near_gradient.reshape(-1, 18)))
    residual = np.r_[train_residual.ravel(), nearby_residual.ravel()]
    pressure_fits = {}
    for name, active, limit in (
            ('old_six', list(range(6)), None),
            ('all_eighteen', list(range(18)), None),
            ('all_cap_1e3', list(range(18)), 1e3)):
        fit = minmax_pressure(A, rhs, matrix, residual,
                              active, limit)
        if fit['feasible']:
            pc = np.asarray(fit['coefficients'])
            fit.update(train_max=stats(train_residual+
                                       train_gradient@pc)['max'],
                       nearby_max=stats(nearby_residual+
                                        near_gradient@pc)['max'])
        pressure_fits[name] = fit
        print(json.dumps(dict(name=name, feasible=fit['feasible'],
                              component_max=fit.get('actual_component_max'),
                              train_max=fit.get('train_max'),
                              nearby_max=fit.get('nearby_max'))), flush=True)
    _, moment_rows = load_cache(ROOT/OUTER_CACHE, signature())
    original_c = np.asarray(mean_source['coefficients'])
    moments = []
    for row in moment_rows:
        U = row['U']+row['Bu']@original_c
        E = row['E']+row['Be']@original_c
        for mode, amplitude in zip([broad]+split_modes,
                                   [AMPLITUDE]+amplitudes):
            du, de = mode_profile(mode, base, row['X'],
                                  row['eta'], tau)
            U += amplitude*du
            E += amplitude*de
        defect = moment_vector(U, E, row['X'], row['weights'])-row['target']
        moments.append(dict(eta=row['eta'], defect=defect.tolist(),
                            min_relative_E=float(np.min(E/row['E0']))))
    report = dict(source=('delayed_axial_split_amplitude_lp.json' if tuned
                          else 'delayed_axial_split_swirl_response.json'),
                  response_cache=RESPONSE_CACHE,
                  split_amplitudes=amplitudes,
                  tau=tau, X=X.tolist(), eta=eta.tolist(),
                  pointwise=pointwise,
                  free_axial_count=sum(row['free_axial'] for row in pointwise),
                  train_before=stats(train_residual),
                  nearby_before=stats(nearby_residual),
                  moments=moments, pressure_fits=pressure_fits,
                  scope='Two axial-split solenoidal swirl corrections '
                        'with updated physical stress at formerly '
                        'blocked nodes; sampled pressure cone LP on '
                        '15 nodes, full center/nearby momentum and '
                        'three-slice moments. No continuous cone, '
                        'supported wave, volume L2, or PDE acceptance.',
                  accepted=False)
    output_name = ('delayed_axial_split_pressure_fit_tuned.json' if tuned
                   else 'delayed_axial_split_pressure_fit.json')
    (ROOT/output_name).write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps(dict(free_axial_count=report['free_axial_count'],
                          train_before=report['train_before'],
                          nearby_before=report['nearby_before'],
                          largest_moment_defect=max(
                              abs(value) for row in moments
                              for value in row['defect']))), flush=True)


if __name__ == '__main__':
    import sys
    run(tuned='--tuned' in sys.argv[1:])
