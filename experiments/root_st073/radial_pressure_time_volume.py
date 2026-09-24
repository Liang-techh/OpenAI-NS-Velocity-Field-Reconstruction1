"""Physical-volume residual comparison for the fixed time-box pressure fit."""
import json

import numpy as np

from joined_field import ROOT
from joint_collar_fit import nodes
from radial_peak_cone import operator
from radial_pressure_patch_screen import load_dense_candidate
from radial_pressure_time_validate import load_candidate


def run():
    before, after = load_dense_candidate(), load_candidate()
    rows = []
    for tau in (.00825, .0084, .00855):
        for order in (6, 8):
            points, weights = nodes(before, tau, order)
            left = operator(before, points, tau)[2]
            right = operator(after, points, tau)[2]
            def measure(residual):
                norm = np.linalg.norm(residual, axis=1)
                return {'max': float(np.max(norm)),
                        'physical_volume_l2': float(np.sqrt(weights@norm**2))}
            row = {'tau': tau, 'order': order,
                   'before': measure(left), 'after': measure(right)}
            rows.append(row)
            print(json.dumps(row), flush=True)
    report = {'rows': rows,
              'scope': 'Three time slices, Gauss6/8 spatial quadrature. These orders disagree and do not establish convergence, a spacetime bound, or momentum accuracy.',
              'accepted': False}
    path = ROOT/'compact_potential'/'radial_pressure_time_volume.json'
    path.write_bytes((json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
