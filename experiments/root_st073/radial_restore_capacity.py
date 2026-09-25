"""Screen radial step return location against delayed U moment capacity."""

import json

import numpy as np

from azimuthal_capacity_optimize import INTERVALS
from delayed_taper_capacity_screen import capacity, grid
from extended_relaxed_cone_screen import profile
from high_frequency_shear_screen import make_field
from moment_shear_slice_repair import moment_vector
from paper_moment_bridge import bump
from radial_continuation import ROOT
from radial_moment_step import RadialMomentStep


def run():
    tau = .5*2**(-5.5)
    X, weights = grid(order=48)
    base = make_field(16, 2.)
    source = json.loads((ROOT/'delayed_e_capacity_optimize.json').read_text())
    e_coeff = {row['eta']: np.asarray(row['coefficients'])
               for row in source['rows']}
    eb = np.array([bump(X, *interval)[0] for interval in INTERVALS])
    targets = {}
    for eta in (.2, .3):
        U0, E0 = profile(base, X, eta, tau)
        targets[eta] = dict(E=E0+e_coeff[eta]@eb,
                            target=moment_vector(U0, E0, X, weights))
    rows = []
    for restore_start in (1.55, 1.75, 2., 2.2, 2.4):
        changed = RadialMomentStep(base, -.5,
                                   restore_start=restore_start,
                                   restore_end=3.)
        for eta in (.2, .3):
            U, _ = profile(changed, X, eta, tau)
            for degree in (19, 23, 27, 31):
                result = capacity(X, weights, U, targets[eta]['E'],
                                  targets[eta]['target'],
                                  1.005, .4, degree)
                row = dict(restore_start=restore_start, eta=eta,
                           **result)
                rows.append(row)
            print(json.dumps(dict(restore_start=restore_start, eta=eta,
                                  slacks={row['degree']: row['S_slack']
                                          for row in rows[-4:]})),
                  flush=True)
    report = dict(tau=tau, step_amplitude=-.5,
                  restore_end=3., width=.4, start_X=1.005,
                  quadrature_order=48, E_source='delayed_e_capacity_optimize.json',
                  rows=rows,
                  scope='Finite-basis fixed-slice M/J/S capacity as the '
                        'outer radial streamfunction restoration starts '
                        'at different X. Entrance step and E profile '
                        'fixed. No physical complete-momentum or cone '
                        'acceptance.', accepted=False)
    (ROOT/'radial_restore_capacity.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
