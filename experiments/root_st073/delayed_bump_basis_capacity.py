"""Screen local C4 radial U bumps for exact five-moment capacity."""

import json

import numpy as np

from azimuthal_capacity_optimize import INTERVALS
from coupled_five_moment_slice import correction_modes
from delayed_taper_capacity_screen import grid
from extended_relaxed_cone_screen import profile
from high_frequency_shear_screen import make_field
from moment_shear_slice_repair import moment_vector
from paper_moment_bridge import bump
from radial_continuation import ROOT
from radial_moment_step import RadialMomentStep


def bump_modes(X, count, support_width, start_X=1.005, end_X=3.):
    starts = np.linspace(start_X, end_X-support_width, count)
    return np.array([bump(X, float(a), float(a+support_width))[0]
                     for a in starts])


def capacity(B, X, weights, U, E, target):
    orth, factor = np.linalg.qr((B*np.sqrt(weights)).T)
    Q = orth.T/np.sqrt(weights)
    H = np.sqrt(2*X)*E
    linear = np.array([Q@weights, Q@(weights*H)])
    required = np.array([target[0]-weights@U,
                         target[2]-weights@(U*H)])
    dual = linear@linear.T
    minimum = -(Q@(weights*U))
    minimum += linear.T@np.linalg.solve(
        dual, required-linear@minimum)
    U_min = U+minimum@Q
    S_min = float(weights@(U_min**2-E**2/2))
    return dict(S_slack=float(target[3]-S_min),
                basis_condition=float(np.linalg.cond(factor)),
                max_abs_minimum_U=float(np.max(np.abs(U_min))))


def run():
    tau = .5*2**(-5.5)
    X, weights = grid(order=48)
    base = make_field(16, 2.)
    changed = RadialMomentStep(base, -.5)
    e_data = json.loads((ROOT/'delayed_e_capacity_optimize.json').read_text())
    e_coeff = {row['eta']: np.asarray(row['coefficients'])
               for row in e_data['rows']}
    eb = np.array([bump(X, *interval)[0] for interval in INTERVALS])
    profiles = {}
    for eta in (.2, .3):
        U0, E0 = profile(base, X, eta, tau)
        U, _ = profile(changed, X, eta, tau)
        E = E0+e_coeff[eta]@eb
        profiles[eta] = (U, E, moment_vector(U0, E0, X, weights))
    rows = []
    global_modes = correction_modes(X, .4, 11, 1.005)
    for support_width in (.3, .5, .8):
        for count in (16, 24, 32):
            local = bump_modes(X, count, support_width)
            for kind, B in (('local', local),
                            ('hybrid_degree11', np.vstack((global_modes,
                                                           local)))):
                row = dict(kind=kind, support_width=support_width,
                           count=count)
                for eta in (.2, .3):
                    U, E, target = profiles[eta]
                    row[str(eta)] = capacity(B, X, weights, U, E, target)
                rows.append(row)
                print(json.dumps({key: row[key] for key in
                                  ('kind', 'support_width', 'count',
                                   '0.2', '0.3')}), flush=True)
    report = dict(tau=tau, start_X=1.005, end_X=3.,
                  quadrature_order=48, rows=rows,
                  scope='Finite-basis fixed-eta M/J/S capacity with E '
                        'fixed to the reoptimized five-bump profile. '
                        'Each U mode is compact C4 and endpoint-flat. '
                        'No physical lift, stress cone, or momentum '
                        'claim.', accepted=False)
    (ROOT/'delayed_bump_basis_capacity.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
