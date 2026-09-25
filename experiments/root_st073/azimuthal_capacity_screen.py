"""Test whether positive E redistribution can reopen U-only moment capacity.

For each outer E bump amplitude, two other disjoint bumps restore I and
Cp. An exact Hilbert projection then gives the minimum U kinetic moment
needed to restore M and J after X=1. Positive slack would make all five
moments feasible at a fixed eta slice, before smooth physical lifting.
"""
import json

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import least_squares

from extended_relaxed_cone_screen import profile
from high_frequency_shear_screen import make_field
from moment_shear_slice_repair import moment_vector
from paper_moment_bridge import bump
from radial_continuation import ROOT
from radial_moment_step import RadialMomentStep


INTERVALS = ((1.005, 1.40), (1.45, 2.15), (2.20, 2.98))
OUTER_AMPLITUDES = (-.14, -.10, -.05, 0., .05, .10, .20, .40, .80, 1.20)


def grid(order=96):
    edges = sorted(set((0., 3.5, 3/32, 1., 1.5, 1.75, 3.,
                        *(v for interval in INTERVALS for v in interval))))
    g, w = leggauss(order)
    X = np.concatenate([(a+b)/2+(b-a)*g/2
                        for a, b in zip(edges[:-1], edges[1:])])
    weights = np.concatenate([(b-a)*w/2
                              for a, b in zip(edges[:-1], edges[1:])])
    return X, weights


def solve_slice(base, changed, X, weights, eta, tau):
    U0, E0 = profile(base, X, eta, tau)
    U, _ = profile(changed, X, eta, tau)
    target = moment_vector(U0, E0, X, weights)
    eb = np.array([bump(X, *interval)[0] for interval in INTERVALS])
    lower = [float(np.max(-E0[b > 1e-10]/b[b > 1e-10])+1e-6)
             for b in eb]
    inside = (X > 1.) & (X < 3.)
    outside = ~inside
    wi = weights[inside]
    rows = []
    previous = np.maximum(np.zeros(2), np.array(lower[:2])+1e-5)
    for outer in OUTER_AMPLITUDES:
        if outer <= lower[2]:
            rows.append(dict(eta=eta, outer=outer,
                             failure='outer bump violates positive E bound'))
            continue

        def restore(two):
            E = E0+two@eb[:2]+outer*eb[2]
            data = moment_vector(U0, E, X, weights)
            return (data-target)[[1, 4]]

        fit = least_squares(restore, previous,
                            bounds=(lower[:2], [np.inf, np.inf]),
                            max_nfev=500, xtol=1e-12, ftol=1e-12,
                            gtol=1e-12)
        previous = fit.x
        E = E0+fit.x@eb[:2]+outer*eb[2]
        H = np.sqrt(2*X)*E
        Hi = H[inside]
        gram = np.array([[sum(wi), wi@Hi],
                         [wi@Hi, wi@(Hi**2)]])
        need = np.array([target[0]-weights[outside]@U[outside],
                         target[2]-weights[outside]@(
                             U[outside]*H[outside])])
        min_U2 = (float(need@np.linalg.solve(gram, need))
                  +float(weights[outside]@(U[outside]**2)))
        required_U2 = float(target[3]+.5*weights@(E**2))
        row = dict(eta=eta, outer=outer,
                   inner_coefficients=fit.x.tolist(),
                   I_Cp_defect=restore(fit.x).tolist(),
                   min_E=float(np.min(E)),
                   required_U2=required_U2,
                   minimum_U2=min_U2,
                   energy_slack=required_U2-min_U2,
                   candidate_feasible_by_projection=bool(
                       np.max(np.abs(restore(fit.x))) < 1e-8
                       and np.min(E) > 0 and required_U2 >= min_U2))
        rows.append(row)
    return rows


def run():
    tau = .5*2**(-5.5)
    base = make_field(16, 2.)
    changed = RadialMomentStep(base, -.5)
    X, weights = grid()
    rows = []
    for eta in (.2, .3):
        rows.extend(solve_slice(base, changed, X, weights, eta, tau))
    for row in rows:
        print(json.dumps(row), flush=True)
    report = dict(tau=tau, intervals=INTERVALS,
                  quadrature_per_piece=96, rows=rows,
                  scope='Fixed-eta E redistribution with I/Cp restoration and variational U feasibility; no smooth eta/time lift, solenoidal 3D field, continuous cone or full NS acceptance.',
                  accepted=False)
    (ROOT/'azimuthal_capacity_screen.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
