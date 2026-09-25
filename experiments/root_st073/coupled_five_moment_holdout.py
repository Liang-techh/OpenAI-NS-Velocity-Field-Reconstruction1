"""Reintegrate fixed-slice five-moment solutions on independent nodes."""
import json

import numpy as np

from azimuthal_capacity_optimize import INTERVALS, grid
from coupled_five_moment_slice import correction_modes
from extended_relaxed_cone_screen import profile
from high_frequency_shear_screen import make_field
from moment_shear_slice_repair import moment_vector
from paper_moment_bridge import bump
from radial_continuation import ROOT
from radial_moment_step import RadialMomentStep


def run():
    tau = .5*2**(-5.5)
    base = make_field(16, 2.)
    changed = RadialMomentStep(base, -.5)
    source = json.loads((ROOT/'coupled_five_moment_slice.json').read_text())
    rows = []
    for order in (128, 160):
        X, weights = grid(order)
        eb = np.array([bump(X, *interval)[0] for interval in INTERVALS])
        ub = correction_modes(X)
        for item in source['rows']:
            if not item.get('five_moments_restored'):
                continue
            eta = item['eta']
            U0, E0 = profile(base, X, eta, tau)
            U, _ = profile(changed, X, eta, tau)
            E = E0+np.asarray(item['e_coefficients'])@eb
            U = U+np.asarray(item['u_coefficients'])@ub
            delta = (moment_vector(U, E, X, weights)
                     -moment_vector(U0, E0, X, weights))
            row = dict(order=order, eta=eta, variant=item['variant'],
                       delta=delta.tolist(),
                       max_abs_defect=float(np.max(np.abs(delta))),
                       min_relative_E=float(np.min(E/E0)))
            rows.append(row)
            print(json.dumps(row), flush=True)
    report = dict(tau=tau, rows=rows,
                  scope='Independent Gauss128/160 quadrature of fixed-slice coefficients fitted at Gauss96. This only tests radial integration, not eta/time smoothness, solenoidal lifting, cone width, or full NS residual.',
                  accepted=False)
    (ROOT/'coupled_five_moment_holdout.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
