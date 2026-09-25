"""See whether more smooth radial modes rescue a wider moment repair."""
import json

import numpy as np

from azimuthal_capacity_optimize import INTERVALS
from coupled_five_moment_slice import correction_modes
from extended_relaxed_cone_screen import profile
from high_frequency_shear_screen import make_field
from moment_shear_slice_repair import moment_vector
from paper_moment_bridge import bump
from radial_continuation import ROOT
from radial_moment_step import RadialMomentStep
from taper_width_capacity import grid


def minimum_slack(U0, E0, U, E, X, weights, width, degree):
    B = correction_modes(X, width, degree)
    orth, factor = np.linalg.qr((B*np.sqrt(weights)).T)
    Q = orth.T/np.sqrt(weights)
    target = moment_vector(U0, E0, X, weights)
    H = np.sqrt(2*X)*E
    zero = -(Q@(weights*U))
    constraints = np.array([Q@weights, Q@(weights*H)])
    needed = np.array([target[0]-weights@U,
                       target[2]-weights@(U*H)])
    coefficients = zero+constraints.T@np.linalg.solve(
        constraints@constraints.T, needed-constraints@zero)
    U_min = U+coefficients@Q
    achieved = moment_vector(U_min, E, X, weights)
    return dict(S_slack=float(target[3]-achieved[3]),
                max_abs_U=float(np.max(np.abs(U_min))),
                condition=float(np.linalg.cond(factor)))


def run():
    tau = .5*2**(-5.5)
    base = make_field(16, 2.)
    changed = RadialMomentStep(base, -.5)
    source = json.loads((ROOT/'taper_width_capacity.json').read_text())
    rows = []
    for prior in source['rows']:
        width, eta = prior['width'], prior['eta']
        if width == .02:
            continue
        X, weights = grid(width)
        U0, E0 = profile(base, X, eta, tau)
        U, _ = profile(changed, X, eta, tau)
        eb = np.array([bump(X, *interval)[0] for interval in INTERVALS])
        E = E0+np.asarray(prior['coefficients'])@eb
        for degree in (11, 19, 31):
            row = dict(width=width, eta=eta, degree=degree,
                       **minimum_slack(U0, E0, U, E, X, weights,
                                       width, degree))
            rows.append(row)
            print(json.dumps(row), flush=True)
    report = dict(tau=tau, rows=rows,
                  scope='Finite-basis lower S moment at fixed E coefficients from width optimization. Higher polynomial degree may reduce L2 cost but can worsen derivatives. No exact physical field, pressure, cone or momentum acceptance.',
                  accepted=False)
    (ROOT/'taper_degree_capacity.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
