"""Unfitted volume and time checks for the refitted meridional correction."""
import json

import numpy as np

from axial_swirl_multitime_screen import load_dense_candidate
from dynamic_poloidal_multitime_screen import full_residual, make_modes
from dynamic_poloidal_volume_refit import build_block, evaluate, metrics
from joined_field import ROOT
from joint_collar_fit import kinematics


def relative_collar_points(base, tau):
    _, radius, zflat, zsupport = base.support(tau)
    return np.array([[rf*radius, 0., sign*(zflat+s*(zsupport-zflat))]
                     for sign in (-1, 1) for rf in (.3, .5)
                     for s in (.275, .425, .575, .725)])


def point_metrics(residual):
    norms = np.linalg.norm(residual, axis=1)
    return {'max': float(np.max(norms)),
            'radial_max': float(np.max(np.abs(residual[:, 0]))),
            'angular_max': float(np.max(np.abs(residual[:, 1]))),
            'axial_max': float(np.max(np.abs(residual[:, 2])))}


def run():
    report = json.loads((ROOT/'compact_potential'/'dynamic_poloidal_volume_refit.json').read_text())
    amplitudes = np.array(report['amplitudes'])
    base = load_dense_candidate()
    modes = make_modes(base)
    tau = .0084
    volume = build_block(base, modes, tau, 10)
    volume_row = {'tau': tau, 'order': 10, 'point_count': len(volume['points']),
                  'baseline': metrics(volume, np.zeros(len(modes))),
                  'refit': metrics(volume, amplitudes)}
    print(json.dumps({'volume_order': 10, 'baseline_max': volume_row['baseline']['max'],
                      'refit_max': volume_row['refit']['max'],
                      'baseline_l2': volume_row['baseline']['physical_volume_l2'],
                      'refit_l2': volume_row['refit']['physical_volume_l2']}), flush=True)
    rows = []
    for tau in (.0084, .012, .024):
        points = relative_collar_points(base, tau)
        hs = .0005*np.sqrt(base.nu*tau)
        ht = .0001*tau
        u0, J0, part0 = kinematics(base, points, tau, hs, ht)
        baseline = full_residual(u0, J0, part0)
        for mode, amplitude in zip(modes, amplitudes):
            if amplitude:
                du, dJ, dpart = kinematics(mode, points, tau, hs, ht)
                u0 += amplitude*du
                J0 += amplitude*dJ
                part0 += amplitude*dpart
        corrected = full_residual(u0, J0, part0)
        row = {'tau': tau, 'point_count': len(points),
               'baseline': point_metrics(baseline),
               'refit': point_metrics(corrected)}
        rows.append(row)
        print(json.dumps(row), flush=True)
    result = {'volume_order10': volume_row, 'collar_holdout': rows,
              'scope': 'Gauss10 whole-support volume is disjoint from Gauss6/8 fitting; relative collar points are also disjoint. Finite sampled checks, not continuum max or converged L2 proof. Other times were not fit.',
              'accepted': False}
    (ROOT/'compact_potential'/'dynamic_poloidal_volume_refit_holdout.json').write_bytes(
        (json.dumps(result, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
