"""Jointly tune two split-swirl amplitudes for smaller cone pressure.

The exact quadratic stress response changes only the X=1.03 low/high
eta nodes on this 15-node set. Other stress targets and pressure columns
are frozen. Each amplitude pair solves a linear program minimizing the
largest absolute coefficient among eighteen compact pressure modes.
"""

import json

import numpy as np
from scipy.optimize import linprog

from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from delayed_axial_split_swirl_response import load_rows
from delayed_multimode_cone_fit import BASE_NAME
from delayed_multimode_pressure_admission import geometry as response_geometry
from delayed_outer_radial_extension import make_extended_modes
from delayed_plateau_pressure_screen import (PlateauPressurePatch,
                                              plateau_stress_columns)
from delayed_pressure_cone_window import pressure_stress_columns
from delayed_shifted_swirl_cone_fit import (AMPLITUDE, NEW_AXIAL,
                                             NEW_RADIAL,
                                             geometry as shear_geometry)
from delayed_midband_moment_patch import SimilaritySwirlMode
from delayed_similarity_curl_screen import CurlPatchedLift
from delayed_similarity_pressure_screen import SimilarityPressurePatch
from joint_collar_fit import kinematics
from radial_continuation import ROOT


def pressure_lp(A, rhs):
    n = A.shape[1]
    scale = 1000.
    fit = linprog(
        np.r_[np.zeros(n), 1.],
        A_ub=np.vstack((np.column_stack((scale*A,
                                         np.zeros(len(rhs)))),
                        np.column_stack((np.eye(n), -np.ones(n))),
                        np.column_stack((-np.eye(n), -np.ones(n))))),
        b_ub=np.r_[rhs, np.zeros(2*n)],
        bounds=[(None, None)]*n+[(0, None)], method='highs')
    if not fit.success:
        return dict(feasible=False, message=str(fit.message))
    coefficients = scale*fit.x[:n]
    return dict(feasible=True,
                max_abs_coefficient=float(np.max(np.abs(coefficients))),
                coefficients=coefficients.tolist(),
                min_cone_slack=float(np.min(rhs-A@coefficients)))


def run():
    shifted = json.loads((ROOT/'delayed_shifted_swirl_cone_fit.json').read_text())
    split = json.loads((ROOT/'delayed_axial_split_swirl_response.json').read_text())
    mean_source = json.loads((ROOT/'delayed_outer_admission_search.json').read_text())
    response_rows = load_rows()
    base = CoupledMomentPhysicalLift(slice_filename=BASE_NAME)
    broad = SimilaritySwirlMode(base, NEW_RADIAL,
                                eta_interval=NEW_AXIAL)
    shifted_mean = CurlPatchedLift(
        base, make_extended_modes(base)+[broad],
        [*mean_source['coefficients'], AMPLITUDE])
    old_pressure = SimilarityPressurePatch(base, np.zeros(6))
    plateau = PlateauPressurePatch(base)
    rows = shifted['geometry']
    tau = shifted['tau']
    points = base.compact.joined.inner.from_similarity(
        np.array([row['X'] for row in rows]),
        np.array([row['eta'] for row in rows]), tau)
    hs, ht = .0005*np.sqrt(base.nu*tau), .0001*tau
    u, grad, _ = kinematics(shifted_mean, points, tau, hs, ht)
    fixed = []
    for i, row in enumerate(rows):
        columns = np.column_stack((
            pressure_stress_columns(old_pressure, points[i], tau),
            plateau_stress_columns(plateau, points[i], tau)))
        if row['X'] == 1.03 and row['eta'] in (.2, .3):
            fixed.append(None)
            continue
        target = np.asarray(row['target'])
        N, K, lam2, multiplier = shear_geometry(
            u[i], grad[i], points[i, 0])
        if lam2 <= 0 or multiplier is None:
            raise ValueError('A fixed node has invalid shear')
        directions = [.8*N+sign*multiplier*K for sign in (-1., 1.)]
        fixed.append(dict(columns=columns, target=target,
                          directions=directions))
    intervals = [row['amplitude_range'] for row in split['rows']]
    trials = []
    for low in np.linspace(*intervals[0], 11):
        for high in np.linspace(*intervals[1], 11):
            amplitudes = (float(low), float(high))
            A, rhs = [], []
            for i, row in enumerate(rows):
                if fixed[i] is None:
                    k = 0 if row['eta'] == .2 else 1
                    c = np.zeros(2)
                    c[k] = amplitudes[k]
                    target, N, K, multiplier, lam2 = response_geometry(
                        response_rows[k], c)
                    if lam2 <= 0 or multiplier is None:
                        A = []
                        break
                    directions = [.8*N+sign*multiplier*K
                                  for sign in (-1., 1.)]
                    columns = np.column_stack((
                        pressure_stress_columns(
                            old_pressure, points[i], tau),
                        plateau_stress_columns(
                            plateau, points[i], tau)))
                else:
                    target = fixed[i]['target']
                    directions = fixed[i]['directions']
                    columns = fixed[i]['columns']
                for direction in directions:
                    A.append(direction@columns)
                    rhs.append(-.1-direction@target)
            fit = (pressure_lp(np.asarray(A), np.asarray(rhs))
                   if A else dict(feasible=False))
            trials.append(dict(amplitudes=amplitudes,
                               pressure=fit))
    feasible = [row for row in trials if row['pressure']['feasible']]
    best = min(feasible, key=lambda row:
               row['pressure']['max_abs_coefficient']) if feasible else None
    original_amplitudes = [row['selected']['amplitude']
                           for row in split['rows']]
    nearest = min(trials, key=lambda row:
                  sum((a-b)**2 for a, b in zip(row['amplitudes'],
                                                 original_amplitudes)))
    report = dict(source='delayed_axial_split_swirl_response.json',
                  X=1.03, etas=(.2, .3),
                  amplitude_intervals=intervals,
                  grid_shape=(11, 11),
                  original_amplitudes=original_amplitudes,
                  original_nearest_grid=nearest,
                  feasible_count=len(feasible),
                  best=best,
                  trials=[dict(amplitudes=row['amplitudes'],
                               feasible=row['pressure']['feasible'],
                               max_abs_pressure_coefficient=(
                                   row['pressure'].get(
                                       'max_abs_coefficient')))
                          for row in trials],
                  scope='Finite two-amplitude grid minimizing the '
                        'sampled pressure coefficient infinity norm '
                        'under 15-node cone inequalities. No momentum '
                        'objective, continuous cone, exact moments, '
                        'wave, or PDE acceptance.', accepted=False)
    (ROOT/'delayed_axial_split_amplitude_lp.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps(dict(original_nearest=dict(
                              amplitudes=nearest['amplitudes'],
                              max_abs_pressure_coefficient=(
                                  nearest['pressure'].get(
                                      'max_abs_coefficient'))),
                          best=dict(amplitudes=best['amplitudes'],
                                    max_abs_pressure_coefficient=(
                                        best['pressure'][
                                            'max_abs_coefficient']))
                          if best else None,
                          feasible_count=len(feasible))), flush=True)


if __name__ == '__main__':
    run()
