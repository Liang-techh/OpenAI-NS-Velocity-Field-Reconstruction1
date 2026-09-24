"""Find how much of robust annular correction preserves positive cone lambda²."""
import json

import numpy as np

from annular_robust_fit import load_candidate
from joined_field import ROOT
from joint_collar_fit import kinematics
from radial_peak_cone import current_field


class CorrectionOnly:
    def __init__(self, candidate):
        self.candidate = candidate
        self.nu = candidate.nu

    def fields(self, points, tau):
        velocity = np.zeros_like(np.asarray(points, float))
        for amplitude, mode in zip(self.candidate.amplitudes,
                                   self.candidate.modes):
            if amplitude:
                velocity += amplitude*mode.fields(points, tau)[0]
        return velocity, np.zeros(len(points))


def lambda_squared(u, J, radius):
    F = u[:, 1]/radius
    shear = np.column_stack((J[:, 1, 0]-F, J[:, 2, 0]))
    shear_norm = np.linalg.norm(shear, axis=1)
    N0 = shear[:, 0]/shear_norm
    return -2*F*N0*(2*F*N0+shear_norm)


def run():
    base = current_field()
    candidate = load_candidate()
    tau = .5/64
    radius = .0056890761915166545
    heights = np.linspace(.0025, .003432627453438968, 13)
    points = np.column_stack((np.full(len(heights), radius),
                              np.zeros(len(heights)), heights))
    hs = .0005*np.sqrt(base.nu*tau)
    ht = .0001*tau
    u0, J0, _ = kinematics(base, points, tau, hs, ht)
    du, dJ, _ = kinematics(CorrectionOnly(candidate), points, tau, hs, ht)
    scales = np.linspace(0., 1., 101)
    rows = []
    for scale in scales:
        lam2 = lambda_squared(u0+scale*du, J0+scale*dJ, radius)
        rows.append({'scale': float(scale),
                     'positive_count': int(np.sum(lam2 > 0)),
                     'minimum_lambda_squared': float(np.min(lam2))})
    allowed = [row['scale'] for row in rows if row['positive_count'] == len(heights)]
    report = {'tau': tau, 'radius': radius, 'heights': heights.tolist(),
              'largest_grid_scale_with_all_positive_lambda_squared': max(allowed) if allowed else None,
              'rows': rows,
              'scope': 'Positive lambda-squared only at thirteen centerline points, no stress target or cone-ratio condition. Not a wave certificate.',
              'accepted': False}
    out = ROOT/'compact_potential'/'annular_cone_scale_screen.json'
    out.write_bytes((json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps({'largest_grid_scale_with_all_positive_lambda_squared':
                      report['largest_grid_scale_with_all_positive_lambda_squared'],
                      'selected': [rows[i] for i in (0, 25, 50, 75, 100)]}),
          flush=True)


if __name__ == '__main__':
    run()
