"""Joint pressure/cone optimization across scales of robust annular velocity."""
import json

import numpy as np

from annular_cone_scale_screen import CorrectionOnly
from annular_pressure_joint_screen import AnnularPressureCandidate
from annular_robust_fit import load_candidate as load_robust
from axial_pressure_cone_screen import AxialPressureMode, make_blocks
from joined_field import ROOT
from joint_collar_fit import kinematics, nodes
from radial_peak_cone import current_field
from transition_poloidal_pressure_screen import solve_pressure


def load_candidate():
    report = json.loads((ROOT/'compact_potential'/'annular_pressure_scale_screen.json').read_text())
    best = report['best']
    robust = load_robust()
    return AnnularPressureCandidate(robust.base, CorrectionOnly(robust),
                                    best['velocity_scale'], best['pressure_amplitudes'])


def kinematic_parts(base, correction, pressure_modes, points, tau, hs, ht):
    u0, J0, part0 = kinematics(base, points, tau, hs, ht)
    du, dJ, dpart = kinematics(correction, points, tau, hs, ht)
    pressure = np.stack([kinematics(mode, points, tau, hs, ht)[2]
                         for mode in pressure_modes], axis=-1)
    return u0, J0, part0, du, dJ, dpart, pressure


def recombine(data, scale):
    u0, J0, part0, du, dJ, dpart, pressure = data
    u = u0+scale*du
    J = J0+scale*dJ
    part = part0+scale*dpart
    residual = part+np.einsum('nij,nj->ni', J, u)
    return u, J, residual, pressure


def run():
    robust = load_robust()
    base = robust.base
    correction = CorrectionOnly(robust)
    pressure_modes = [AxialPressureMode(base.base, degree) for degree in range(4)]
    tau = .5/64
    radius = .0056890761915166545
    heights = tuple(np.linspace(.0025, .003432627453438968, 13))
    blocks = make_blocks(base, tau, radius, heights)
    line_points = np.vstack([block[3] for block in blocks])
    hs = .0005*np.sqrt(base.nu*tau)
    ht = .0001*tau
    line = kinematic_parts(base, correction, pressure_modes,
                           line_points, tau, hs, ht)
    volume = []
    for order in (6, 8, 10):
        points, weights = nodes(base, tau, order)
        data = kinematic_parts(base, correction, pressure_modes,
                               points, tau, hs, ht)
        _, _, baseline, _ = recombine(data, 0.)
        baseline_l2 = float(np.sqrt(weights @ np.linalg.norm(baseline, axis=1)**2))
        volume.append({'order': order, 'weights': weights, 'data': data,
                       'baseline_l2': baseline_l2})
        print(json.dumps({'loaded_order': order, 'baseline_l2': baseline_l2}),
              flush=True)
    scales = (0., .1, .2, .3, .35, .38, .4, .42, .44, .46, .48)
    candidates = []
    for scale in scales:
        lu, lJ, lR, lP = recombine(line, scale)
        states = [recombine(item['data'], scale) for item in volume]
        objective_residual = np.vstack([state[2] for state in states])
        objective_columns = np.concatenate([state[3] for state in states], axis=0)
        objective_weights = np.concatenate([
            item['weights']/item['baseline_l2']**2 for item in volume])
        fit = solve_pressure(lu, lJ, lR, lP, blocks, radius, ratio=.8,
                             objective_residual=objective_residual,
                             objective_columns=objective_columns,
                             objective_weights=objective_weights)
        if fit is None:
            print(json.dumps({'velocity_scale': scale, 'feasible': False}),
                  flush=True)
            continue
        volume_rows = []
        for item, state in zip(volume, states):
            residual = state[2]+np.einsum('nik,k->ni', state[3],
                                          fit['pressure_amplitudes'])
            norm = np.linalg.norm(residual, axis=1)
            l2 = float(np.sqrt(item['weights'] @ norm**2))
            volume_rows.append({'order': item['order'],
                                'baseline_l2': item['baseline_l2'],
                                'joint_l2': l2,
                                'ratio': l2/item['baseline_l2'],
                                'joint_max': float(np.max(norm))})
        candidate = {'velocity_scale': scale,
                     'pressure_amplitudes': fit['pressure_amplitudes'],
                     'pass_count': fit['pass_count'],
                     'max_cone_ratio': fit['max_cone_ratio'],
                     'volume_rows': volume_rows,
                     'worst_l2_ratio': max(row['ratio'] for row in volume_rows)}
        candidates.append(candidate)
        print(json.dumps({'velocity_scale': scale,
                          'worst_l2_ratio': candidate['worst_l2_ratio']}),
              flush=True)
    candidates.sort(key=lambda item: item['worst_l2_ratio'])
    report = {'tau': tau, 'radius': radius, 'heights': heights,
              'screened_velocity_scales': scales,
              'candidates': candidates,
              'best': candidates[0] if candidates else None,
              'scope': 'Exact nonlinear velocity recombination and linear pressure fit at one time; local cone on 13 axial points and Gauss6/8/10 physical-volume objective. No midpoint or time holdout for selected scale, radial patch, or continuum certificate.',
              'accepted': False}
    out = ROOT/'compact_potential'/'annular_pressure_scale_screen.json'
    out.write_bytes((json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps({'best': report['best']}), flush=True)


if __name__ == '__main__':
    run()
