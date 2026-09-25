"""Independent full-momentum check for the time-scaled collar candidate."""
import json

import numpy as np

from axial_swirl_multitime_screen import load_dense_candidate, momentum
from dynamic_poloidal_multitime_screen import load_dynamic_candidate
from joined_field import ROOT


def run():
    base = load_dense_candidate()
    candidate = load_dynamic_candidate()
    rows = []
    for tau in (.0084, .012, .024):
        _, radius, zflat, zsupport = base.support(tau)
        points = np.array([[rf*radius, 0., sign*(zflat+s*(zsupport-zflat))]
                           for sign in (-1, 1) for rf in (.3, .5)
                           for s in (.275, .425, .575, .725)])
        u0, _, residual0 = momentum(base, points, tau)
        u1, _, residual1 = momentum(candidate, points, tau)
        row = {'tau': tau, 'point_count': len(points),
               'baseline_full_max': float(np.max(np.linalg.norm(residual0, axis=1))),
               'corrected_full_max': float(np.max(np.linalg.norm(residual1, axis=1))),
               'baseline_radial_max': float(np.max(np.abs(residual0[:, 0]))),
               'corrected_radial_max': float(np.max(np.abs(residual1[:, 0]))),
               'baseline_angular_max': float(np.max(np.abs(residual0[:, 1]))),
               'corrected_angular_max': float(np.max(np.abs(residual1[:, 1]))),
               'baseline_axial_max': float(np.max(np.abs(residual0[:, 2]))),
               'corrected_axial_max': float(np.max(np.abs(residual1[:, 2]))),
               'correction_speed_max': float(np.max(np.linalg.norm(u1-u0, axis=1)))}
        rows.append(row)
        print(json.dumps(row), flush=True)
    result = {'rows': rows,
              'scope': 'Sixteen relative axial-collar points per time, disjoint from poloidal training coordinates, at three registered times. Full nonlinear physical momentum. This is not a global max or volume L2 certificate.',
              'accepted': False}
    (ROOT/'compact_potential'/'dynamic_poloidal_holdout.json').write_bytes(
        (json.dumps(result, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
