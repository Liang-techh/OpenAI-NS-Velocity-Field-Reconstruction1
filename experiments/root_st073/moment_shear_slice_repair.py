"""Feasibility screen for four outgoing moments after the M-step.

This operates on normalized fixed-eta profiles only. It must not be used
as a three-dimensional velocity field until the eta-dependent amplitudes
are lifted through a solenoidal streamfunction and retested.
"""
import json

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import least_squares

from extended_relaxed_cone_screen import profile
from high_frequency_shear_screen import make_field
from paper_moment_bridge import bump
from radial_continuation import ROOT
from radial_moment_step import RadialMomentStep


E_INTERVALS = ((1.76, 2.10), (2.12, 2.48), (2.50, 2.98))
U_INTERVAL = (1.76, 2.98)


def grid(order=64):
    edges = sorted(set((0., 3.5, 1.5, 1.75, 3.,
                        3/32, *U_INTERVAL,
                        *(edge for interval in E_INTERVALS
                          for edge in interval))))
    g, w = leggauss(order)
    xs = np.concatenate([(a+b)/2+(b-a)*g/2
                         for a, b in zip(edges[:-1], edges[1:])])
    ws = np.concatenate([(b-a)*w/2
                         for a, b in zip(edges[:-1], edges[1:])])
    return xs, ws


def moment_vector(U, E, X, weights):
    H = np.sqrt(2*X)*E
    return np.array([weights@U, weights@H, weights@(U*H),
                     weights@(U**2-E**2/2),
                     weights@(E**2/(2*X))])


def solve_slice(base, field, X, weights, eta, tau):
    original_U, original_E = profile(base, X, eta, tau)
    U, E = profile(field, X, eta, tau)
    target = moment_vector(original_U, original_E, X, weights)
    ub = bump(X, *U_INTERVAL)[1]
    eb = np.array([bump(X, *interval)[0] for interval in E_INTERVALS])

    def evaluate(coefficients):
        new_U = U+coefficients[0]*ub
        new_E = E+coefficients[1:]@eb
        return new_U, new_E, moment_vector(new_U, new_E, X, weights)

    def objective(coefficients):
        return (evaluate(coefficients)[2]-target)[1:]

    initial = np.zeros(4)
    initial_residual = objective(initial)
    jac = np.column_stack([(objective(initial+np.eye(4)[i]*1e-5)
                            -objective(initial-np.eye(4)[i]*1e-5))/2e-5
                           for i in range(4)])
    singular_values = np.linalg.svd(jac, compute_uv=False)
    linear = np.linalg.solve(jac, -initial_residual)
    rng = np.random.default_rng(73073)
    starts = [initial, linear, -linear]
    starts += [linear+rng.normal(size=4)*s
               for s in (.1, .3, 1.) for _ in range(4)]
    candidates = [least_squares(objective, start, max_nfev=500,
                                xtol=1e-12, ftol=1e-12, gtol=1e-12)
                  for start in starts]
    fit = min(candidates, key=lambda result: np.linalg.norm(result.fun))
    new_U, new_E, achieved = evaluate(fit.x)
    return dict(eta=eta,
                initial_delta=(moment_vector(U, E, X, weights)-target).tolist(),
                initial_jacobian_singular_values=singular_values.tolist(),
                linear_coefficients=linear.tolist(),
                starts_tried=len(starts),
                coefficients=fit.x.tolist(),
                final_delta=(achieved-target).tolist(),
                max_abs_four_moment_defect=float(np.max(np.abs(
                    (achieved-target)[1:]))),
                min_corrected_E=float(np.min(new_E)),
                max_abs_corrected_U=float(np.max(np.abs(new_U))),
                optimizer_success=bool(fit.success),
                moment_repaired=bool(np.max(np.abs(
                    (achieved-target)[1:])) < 1e-6),
                positive_swirl=bool(np.min(new_E) > 0))


def run():
    tau = .5*2**(-5.5)
    base = make_field(16, 2.)
    field = RadialMomentStep(base, -.5)
    X, weights = grid()
    rows = [solve_slice(base, field, X, weights, eta, tau)
            for eta in (.2, .3)]
    for row in rows:
        print(json.dumps(row), flush=True)
    report = dict(tau=tau, x_out=3.5, u_interval=U_INTERVAL,
                  e_intervals=E_INTERVALS, quadrature_per_piece=64,
                  rows=rows,
                  scope='Fixed-eta normalized nonlinear moment feasibility only. M is preserved by the derivative U bump. A slice solve does not define a smooth eta/time coefficient field, physical solenoidal correction, continuous stress cone, or PDE-valid solution.',
                  accepted=False)
    (ROOT/'moment_shear_slice_repair.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
