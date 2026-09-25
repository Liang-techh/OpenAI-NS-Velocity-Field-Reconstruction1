"""Bounded fixed-slice five-moment repair using the active bridge as capacity.

Two derivative U bumps preserve total M, while three E bumps are constrained
to keep sampled swirl positive. This is a profile-level feasibility solve;
eta/time interpolation and a solenoidal physical lift are separate work.
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


U_INTERVALS = ((1.08, 1.45), (1.76, 2.98))
E_INTERVALS = ((1.08, 1.45), (1.76, 2.20), (2.25, 2.98))


def grid(order=64):
    edges = sorted(set((0., 3.5, 3/32, 1.5, 1.75, 3.,
                        *(v for interval in U_INTERVALS+E_INTERVALS
                          for v in interval))))
    g, w = leggauss(order)
    xs = np.concatenate([(a+b)/2+(b-a)*g/2
                         for a, b in zip(edges[:-1], edges[1:])])
    ws = np.concatenate([(b-a)*w/2
                         for a, b in zip(edges[:-1], edges[1:])])
    return xs, ws


def solve_slice(base, field, X, weights, eta, tau):
    original_U, original_E = profile(base, X, eta, tau)
    U, E = profile(field, X, eta, tau)
    target = moment_vector(original_U, original_E, X, weights)
    ub = np.array([bump(X, *interval)[1] for interval in U_INTERVALS])
    eb = np.array([bump(X, *interval)[0] for interval in E_INTERVALS])
    lower_e = [float(np.max(-E[b > 1e-10]/b[b > 1e-10])+1e-5)
               for b in eb]
    lower = np.array([-np.inf, -np.inf, *lower_e])
    upper = np.full(5, np.inf)

    def evaluate(coefficients):
        new_U = U+coefficients[:2]@ub
        new_E = E+coefficients[2:]@eb
        return new_U, new_E, moment_vector(new_U, new_E, X, weights)

    def objective(coefficients):
        return (evaluate(coefficients)[2]-target)[1:]

    rng = np.random.default_rng(73074)
    starts = [np.maximum(np.zeros(5), lower+1e-6)]
    starts += [np.maximum(rng.normal(size=5)*scale, lower+1e-6)
               for scale in (.05, .2, .5, 1.) for _ in range(5)]
    fits = [least_squares(objective, start, bounds=(lower, upper),
                          max_nfev=500, xtol=1e-12, ftol=1e-12,
                          gtol=1e-12) for start in starts]
    fit = min(fits, key=lambda result: np.linalg.norm(result.fun))
    new_U, new_E, achieved = evaluate(fit.x)
    delta = achieved-target
    return dict(eta=eta, starts_tried=len(starts),
                lower_e_bounds=lower_e,
                coefficients=fit.x.tolist(),
                initial_delta=objective(starts[0]).tolist(),
                final_delta=delta.tolist(),
                max_abs_four_moment_defect=float(np.max(np.abs(delta[1:]))),
                min_corrected_E=float(np.min(new_E)),
                max_abs_corrected_U=float(np.max(np.abs(new_U))),
                moment_repaired=bool(np.max(np.abs(delta[1:])) < 1e-6),
                positive_swirl=bool(np.min(new_E) > 0),
                optimizer_success=bool(fit.success))


def run():
    tau = .5*2**(-5.5)
    base = make_field(16, 2.)
    field = RadialMomentStep(base, -.5)
    X, weights = grid()
    rows = [solve_slice(base, field, X, weights, eta, tau)
            for eta in (.2, .3)]
    for row in rows:
        print(json.dumps(row), flush=True)
    report = dict(tau=tau, x_out=3.5, quadrature_per_piece=64,
                  u_intervals=U_INTERVALS, e_intervals=E_INTERVALS,
                  rows=rows,
                  scope='Fixed-slice bounded nonlinear profile solve. Positive swirl is sampled on Gauss nodes only. No eta/time coefficient field, solenoidal lift, continuous cone or PDE acceptance.',
                  accepted=False)
    (ROOT/'moment_shear_bridge_repair.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
