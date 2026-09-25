"""Finite-difference sanity check for the exact-solenoidal streamfunctions."""
import json

import numpy as np

from local_poloidal_basis_screen import LocalPoloidalMode, compact_base
from radial_pressure_volume_constrained import load_candidate


def run():
    base = compact_base(load_candidate())
    points = np.array([[.005, .004, .00345],
                       [.007, .003, .0037],
                       [.005, .004, -.00345]])
    tau = .0084
    rows = []
    for parity in ('even', 'odd'):
        for radial_degree in (0, 1):
            for axial_degree in (0, 1):
                mode = LocalPoloidalMode(base, radial_degree,
                                         axial_degree, parity)
                estimates = []
                for h in (2e-6, 1e-6):
                    divergence = np.zeros(len(points))
                    for axis in range(3):
                        step = h*np.eye(3)[axis]
                        up = mode.fields(points+step, tau)[0]
                        um = mode.fields(points-step, tau)[0]
                        divergence += (up[:, axis]-um[:, axis])/(2*h)
                    estimates.append(float(np.max(np.abs(divergence))))
                rows.append({'parity': parity,
                             'radial_degree': radial_degree,
                             'axial_degree': axial_degree,
                             'max_abs_divergence_fd': estimates})
    report = {'tau': tau, 'points': points.tolist(),
              'steps': [2e-6, 1e-6], 'rows': rows,
              'scope': 'Finite-difference sanity check at three interior points. Exact divergence-freeness follows analytically from the streamfunction representation, not from these samples.'}
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    run()
