"""Independent collar points for the multi-time swirl correction."""
import json
import sys

import numpy as np

from axial_swirl_multitime_screen import AxialSwirlMode, SwirlCandidate, momentum
from joined_field import ROOT
from local_poloidal_basis_screen import load_robust_candidate


def run():
    low_order = '--low-order' in sys.argv
    medium_order = '--medium-order' in sys.argv
    dense = '--dense' in sys.argv
    if sum((low_order, medium_order, dense)) > 1:
        raise ValueError('Choose one basis order')
    source_name = ('axial_swirl_low_order.json' if low_order else
                   'axial_swirl_medium_order.json' if medium_order else
                   'axial_swirl_dense.json' if dense else
                   'axial_swirl_multitime_screen.json')
    report = json.loads((ROOT/'compact_potential'/source_name).read_text())
    base = load_robust_candidate(.1)
    modes = [AxialSwirlMode(base, item['radial_degree'], item['axial_degree'],
                            report['reference_tau'], report['temporal_power'])
             for item in report['mode_ids']]
    candidate = SwirlCandidate(base, modes, report['amplitudes'])
    rows = []
    for tau in (.0084, .012, .024):
        _, radius, zflat, zsupport = base.support(tau)
        radial_fractions = (.3, .5) if dense else (.25, .55)
        axial_fractions = (.275, .425, .575, .725) if dense else (.2, .5, .8)
        points = np.array([[rf*radius, 0., sign*(zflat+s*(zsupport-zflat))]
                           for sign in (-1, 1) for rf in radial_fractions
                           for s in axial_fractions])
        u0, _, residual0 = momentum(base, points, tau)
        u1, _, residual1 = momentum(candidate, points, tau)
        row = {'tau': tau, 'point_count': len(points),
               'baseline_angular_max': float(np.max(np.abs(residual0[:, 1]))),
               'corrected_angular_max': float(np.max(np.abs(residual1[:, 1]))),
               'baseline_full_max': float(np.max(np.linalg.norm(residual0, axis=1))),
               'corrected_full_max': float(np.max(np.linalg.norm(residual1, axis=1))),
               'correction_speed_max': float(np.max(np.linalg.norm(u1-u0, axis=1)))}
        rows.append(row)
        print(json.dumps(row), flush=True)
    result = {'rows': rows, 'low_order': low_order,
              'medium_order': medium_order, 'dense': dense,
              'scope': 'Relative axial-collar points disjoint from the fitting coordinates at three registered times. This is a local holdout, not a full-volume error certificate.',
              'accepted': False}
    filename = ('axial_swirl_low_order_holdout.json' if low_order else
                'axial_swirl_medium_order_holdout.json' if medium_order else
                'axial_swirl_dense_holdout.json' if dense else
                'axial_swirl_holdout.json')
    (ROOT/'compact_potential'/filename).write_bytes(
        (json.dumps(result, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
