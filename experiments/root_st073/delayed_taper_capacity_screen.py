"""Check whether the five moments permit delaying the U repair from X=1.

The physical stress cone survives at X=1 but fails immediately outside it
for the current lift. This is a slice-level capacity test, not a physical
momentum or continuous-cone construction.
"""
import json

import numpy as np
from numpy.polynomial.legendre import leggauss

from azimuthal_capacity_optimize import INTERVALS
from coupled_five_moment_slice import correction_modes
from extended_relaxed_cone_screen import profile
from high_frequency_shear_screen import make_field
from moment_shear_slice_repair import moment_vector
from paper_moment_bridge import bump
from radial_continuation import ROOT
from radial_moment_step import RadialMomentStep


STARTS = (1., 1.005, 1.01, 1.015, 1.02, 1.05, 1.1)
WIDTHS = (.05, .1, .2, .3, .4, .6, .8)


def grid(order=32):
    edges = sorted(set((0., 3.5, 3/32, 1., 1.5, 1.75, 3.,
                        *STARTS, *(s+w for s in STARTS for w in WIDTHS),
                        *(3.-w for w in WIDTHS),
                        *(v for pair in INTERVALS for v in pair))))
    g, w = leggauss(order)
    X = np.concatenate([(a+b)/2+(b-a)*g/2
                        for a, b in zip(edges[:-1], edges[1:])])
    weights = np.concatenate([(b-a)*w/2
                              for a, b in zip(edges[:-1], edges[1:])])
    return X, weights


def capacity(X, weights, U, E, target, start, width, degree):
    B = correction_modes(X, width, degree, start)
    orth, factor = np.linalg.qr((B*np.sqrt(weights)).T)
    Q = orth.T/np.sqrt(weights)
    H = np.sqrt(2*X)*E
    gram = (Q*weights)@Q.T
    linear = np.array([Q@weights, Q@(weights*H)])
    required = np.array([target[0]-weights@U,
                         target[2]-weights@(U*H)])
    inverse = np.linalg.solve(gram, linear.T)
    dual = linear@inverse
    unconstrained = -np.linalg.solve(gram, Q@(weights*U))
    minimum = unconstrained+inverse@np.linalg.solve(
        dual, required-linear@unconstrained)
    u_min = U+minimum@Q
    S_min = float(weights@(u_min**2-E**2/2))
    return dict(start_X=start, width=width, degree=degree,
                S_slack=float(target[3]-S_min),
                basis_condition=float(np.linalg.cond(factor)),
                max_abs_minimum_U=float(np.max(np.abs(u_min))))


def run():
    tau = .5*2**(-5.5)
    X, weights = grid()
    base = make_field(16, 2.)
    changed = RadialMomentStep(base, -.5)
    prior = json.loads((ROOT/'wide_taper_five_moment_slice.json').read_text())
    rows = []
    for source in prior['rows']:
        eta = source['eta']
        U0, E0 = profile(base, X, eta, tau)
        U, _ = profile(changed, X, eta, tau)
        eb = np.array([bump(X, *interval)[0] for interval in INTERVALS])
        E = E0+np.asarray(source['e_coefficients'])@eb
        target = moment_vector(U0, E0, X, weights)
        for start in STARTS:
            for width in WIDTHS:
                for degree in (19, 31):
                    result = capacity(X, weights, U, E, target,
                                      start, width, degree)
                    result.update(eta=eta,
                                  E_I_defect=float(moment_vector(U, E, X, weights)[1]
                                                   -target[1]),
                                  E_Cp_defect=float(moment_vector(U, E, X, weights)[4]
                                                    -target[4]))
                    rows.append(result)
                    print(json.dumps(result), flush=True)
    report = dict(tau=tau, quadrature_per_piece=32, rows=rows,
                  scope='Two fixed-eta slice capacity tests using the '
                        'existing five E coefficients. Positive S slack '
                        'means an M/J/S U repair exists in this finite '
                        'basis; no selected physical field, continuous '
                        'cone, or PDE admission.', accepted=False)
    (ROOT/'delayed_taper_capacity_screen.json').write_text(
        json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    run()
