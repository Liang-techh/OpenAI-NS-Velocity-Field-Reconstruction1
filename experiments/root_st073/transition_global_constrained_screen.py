"""Minimize full-support sampled residual while retaining the local cone."""
import json

import numpy as np

from axial_pressure_cone_screen import AxialPressureMode, make_blocks
from joined_field import ROOT
from joint_collar_fit import TransitionPoloidalMode, kinematics, nodes
from radial_peak_cone import current_field
from transition_poloidal_pressure_screen import solve_pressure


def kinematic_parts(field, mode, pressure_modes, points, tau, hs, ht):
    u0, J0, part0 = kinematics(field, points, tau, hs, ht)
    du, dJ, dpart = kinematics(mode, points, tau, hs, ht)
    pressure = np.stack([kinematics(p, points, tau, hs, ht)[2]
                         for p in pressure_modes], axis=-1)
    return u0, J0, part0, du, dJ, dpart, pressure


def recombine(data, amplitude):
    u0, J0, part0, du, dJ, dpart, pressure = data
    u = u0+amplitude*du
    J = J0+amplitude*dJ
    part = part0+amplitude*dpart
    residual = part+np.einsum('nij,nj->ni', J, u)
    return u, J, residual, pressure


def run():
    field = current_field()
    tau = .5/64
    radius = .0056890761915166545
    heights = tuple(np.linspace(.0025, .003432627453438968, 13))
    blocks = make_blocks(field, tau, radius, heights)
    line_points = np.vstack([block[3] for block in blocks])
    volume_points, volume_weights = nodes(field, tau, 6)
    hs = .0005*np.sqrt(field.nu*tau)
    ht = .0001*tau
    mode = TransitionPoloidalMode(field.base, 2., tau)
    pressure_modes = [AxialPressureMode(field.base, i) for i in range(4)]
    line = kinematic_parts(field, mode, pressure_modes,
                           line_points, tau, hs, ht)
    volume = kinematic_parts(field, mode, pressure_modes,
                             volume_points, tau, hs, ht)
    amplitudes = (-.010, -.008, -.006, -.004, -.002,
                  0., .002, .004, .006, .008, .010)
    rows = []
    for amplitude in amplitudes:
        lu, lJ, lR, lP = recombine(line, amplitude)
        _, _, vR, vP = recombine(volume, amplitude)
        result = solve_pressure(lu, lJ, lR, lP, blocks, radius,
                                objective_residual=vR,
                                objective_columns=vP,
                                objective_weights=volume_weights)
        if result is not None:
            result['poloidal_amplitude'] = amplitude
            result['global_gauss6_max'] = float(np.max(np.linalg.norm(
                vR+np.einsum('nik,k->ni', vP,
                              result['pressure_amplitudes']), axis=1)))
            rows.append(result)
    rows.sort(key=lambda row: row['objective_weighted_l2'])
    _, _, baseline, _ = recombine(volume, 0.)
    report = {'tau': tau, 'radius': radius, 'heights': heights,
              'screened_poloidal_amplitudes': amplitudes,
              'baseline_gauss6_l2': float(np.linalg.norm(
                  np.sqrt(volume_weights)[:, None]*baseline)),
              'results': rows, 'best': rows[0] if rows else None,
              'scope': 'Local cone constraints on 13 axial points, objective Gauss6 physical-volume residual at one time. Full nonlinear velocity recombination, linear pressure. No spacetime/global certificate.',
              'accepted': False}
    out = ROOT/'compact_potential'/'transition_global_constrained_screen.json'
    out.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'baseline_gauss6_l2': report['baseline_gauss6_l2'],
                      'feasible': len(rows),
                      'best': {k: rows[0][k] for k in
                               ('poloidal_amplitude', 'pressure_amplitudes',
                                'pass_count', 'max_cone_ratio',
                                'objective_weighted_l2', 'global_gauss6_max')}
                      if rows else None}), flush=True)


if __name__ == '__main__':
    run()
