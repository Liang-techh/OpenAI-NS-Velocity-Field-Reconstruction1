"""Optimize positive E redistribution while restoring I and Cp exactly.

The objective is the variational slack for an arbitrary downstream U
correction matching M and J. This is a fixed-eta profile feasibility
search and does not itself construct the downstream U or a 3D field.
"""
import json

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import minimize

from extended_relaxed_cone_screen import profile
from high_frequency_shear_screen import make_field
from moment_shear_slice_repair import moment_vector
from paper_moment_bridge import bump
from radial_continuation import ROOT
from radial_moment_step import RadialMomentStep


INTERVALS = ((1.005, 1.18), (1.20, 1.40), (1.45, 1.75),
             (1.80, 2.30), (2.35, 2.98))
RELATIVE_SWIRL_FLOOR = .11
SLACK_RESERVE = .01


def grid(order=96):
    edges = sorted(set((0., 3.5, 3/32, 1., 1.02, 1.5, 1.75, 3.,
                        *(v for interval in INTERVALS for v in interval))))
    g, w = leggauss(order)
    X = np.concatenate([(a+b)/2+(b-a)*g/2
                        for a, b in zip(edges[:-1], edges[1:])])
    weights = np.concatenate([(b-a)*w/2
                              for a, b in zip(edges[:-1], edges[1:])])
    return X, weights


def optimize_slice(base, changed, X, weights, eta, tau):
    U0, E0 = profile(base, X, eta, tau)
    U, _ = profile(changed, X, eta, tau)
    target = moment_vector(U0, E0, X, weights)
    eb = np.array([bump(X, *interval)[0] for interval in INTERVALS])
    lower = np.array([np.max((RELATIVE_SWIRL_FLOOR-1)*E0[b > 1e-10]
                             /b[b > 1e-10])+1e-5
                      for b in eb])
    inside = (X > 1.) & (X < 3.)
    outside = ~inside
    wi = weights[inside]
    outer_U2 = float(weights[outside]@(U[outside]**2))
    outer_M = float(weights[outside]@U[outside])
    outer_J = float(weights[outside]@(
        U[outside]*np.sqrt(2*X[outside])*E0[outside]))

    def metrics(c):
        E = E0+c@eb
        H = np.sqrt(2*X)*E
        I_Cp = np.array([weights@H-target[1],
                         weights@(E**2/(2*X))-target[4]])
        Hi = H[inside]
        gram = np.array([[sum(wi), wi@Hi],
                         [wi@Hi, wi@(Hi**2)]])
        need = np.array([target[0]-outer_M,target[2]-outer_J])
        minimum_U2 = float(need@np.linalg.solve(gram, need)+outer_U2)
        required_U2 = float(target[3]+.5*weights@(E**2))
        return I_Cp, required_U2-minimum_U2, float(np.min(E))

    constraint = dict(type='eq', fun=lambda c: metrics(c)[0])
    rng = np.random.default_rng(73075)
    starts = [np.zeros(5)]
    starts += [np.maximum(rng.normal(size=5)*s, lower+1e-4)
               for s in (.05, .15, .30) for _ in range(3)]
    fits = [minimize(lambda c: -metrics(c)[1]+.001*(c@c),
                     start, method='SLSQP', bounds=list(zip(lower, [1.5]*5)),
                     constraints=[constraint],
                     options=dict(maxiter=500, ftol=1e-12))
            for start in starts]
    eligible = [f for f in fits if np.max(np.abs(metrics(f.x)[0])) < 1e-7]
    fit = max(eligible, key=lambda f: metrics(f.x)[1]) if eligible else min(
        fits, key=lambda f: np.linalg.norm(metrics(f.x)[0]))
    reserve = dict(type='ineq', fun=lambda c: metrics(c)[1]-SLACK_RESERVE)
    small = minimize(lambda c: float(c@c), fit.x, method='SLSQP',
                     bounds=list(zip(lower, [1.5]*5)),
                     constraints=[constraint, reserve],
                     options=dict(maxiter=500, ftol=1e-12))
    selected = (small if small.success and
                np.max(np.abs(metrics(small.x)[0])) < 1e-7 and
                metrics(small.x)[1] >= SLACK_RESERVE-1e-8 else fit)
    errors, slack, minimum_E = metrics(selected.x)
    E_new = E0+selected.x@eb
    return dict(eta=eta, starts_tried=len(starts),
                maximum_slack_found=float(metrics(fit.x)[1]),
                maximum_slack_coefficients=fit.x.tolist(),
                coefficients=selected.x.tolist(),
                I_Cp_defect=errors.tolist(), energy_slack=slack,
                min_E=minimum_E,
                min_relative_E=float(np.min(E_new/E0)),
                optimizer_success=bool(selected.success),
                minimum_norm_with_reserve=bool(selected is small),
                fixed_slice_feasible=bool(np.max(np.abs(errors)) < 1e-7
                                          and slack > 0 and
                                          np.min(E_new/E0) > RELATIVE_SWIRL_FLOOR))


def run():
    tau = .5*2**(-5.5)
    base = make_field(16, 2.)
    changed = RadialMomentStep(base, -.5)
    X, weights = grid()
    rows = [optimize_slice(base, changed, X, weights, eta, tau)
            for eta in (.2, .3)]
    for row in rows:
        print(json.dumps(row), flush=True)
    report = dict(tau=tau, intervals=INTERVALS,
                  relative_swirl_floor=RELATIVE_SWIRL_FLOOR,
                  slack_reserve=SLACK_RESERVE,
                  quadrature_per_piece=96, rows=rows,
                  scope='Bounded five-E-bump fixed-slice optimization with I/Cp equalities and variational M/J/S capacity. Positive slack is not an actual smooth U correction, eta/time lift, continuous cone, or full NS validation.',
                  accepted=False)
    (ROOT/'azimuthal_capacity_optimize.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
