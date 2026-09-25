"""Solve J/S outgoing moments with streamfunction derivative modes only.

Each U mode is d/dX of a compact bump after the physical-cone point, so
its radial integral vanishes. E is unchanged. The J constraint is linear
and S is quadratic; a constrained quadratic minimum gives a direct
feasibility test for this fixed family before any eta/time-dependent lift.
"""
import json

import numpy as np
from numpy.polynomial.legendre import leggauss

from extended_relaxed_cone_screen import profile
from high_frequency_shear_screen import make_field
from moment_shear_slice_repair import moment_vector
from paper_moment_bridge import bump
from radial_continuation import ROOT
from radial_moment_step import RadialMomentStep


U_INTERVALS = ((1.025, 1.25), (1.26, 1.49),
               (1.60, 2.15), (2.20, 2.98))


def grid(order=96):
    edges = sorted(set((0., 3.5, 3/32, 1.5, 1.75, 3.,
                        *(v for interval in U_INTERVALS for v in interval))))
    g, w = leggauss(order)
    xs = np.concatenate([(a+b)/2+(b-a)*g/2
                         for a, b in zip(edges[:-1], edges[1:])])
    ws = np.concatenate([(b-a)*w/2
                         for a, b in zip(edges[:-1], edges[1:])])
    return xs, ws


def solve_slice(base, field, X, weights, eta, tau):
    U0, E = profile(base, X, eta, tau)
    U, _ = profile(field, X, eta, tau)
    target = moment_vector(U0, E, X, weights)
    current = moment_vector(U, E, X, weights)
    modes = np.array([bump(X, *interval)[1] for interval in U_INTERVALS])
    H = np.sqrt(2*X)*E
    j = modes@(weights*H)
    u = modes@(weights*U)
    gram = (modes*weights)@modes.T
    required_j = target[2]-current[2]
    unconstrained = -np.linalg.solve(gram, u)
    dual = np.linalg.solve(gram, j)
    dual_norm = float(j@dual)
    minimum = unconstrained + ((required_j-j@unconstrained)
                               /dual_norm)*dual
    U_min = U+minimum@modes
    achieved_min = moment_vector(U_min, E, X, weights)
    gap = target[3]-achieved_min[3]
    result = dict(eta=eta, initial_delta=(current-target).tolist(),
                  J_response=j.tolist(), constrained_min_coefficients=minimum.tolist(),
                  constrained_min_S=float(achieved_min[3]),
                  target_S=float(target[3]), S_above_minimum=float(gap),
                  fixed_basis_feasible=bool(gap >= -1e-8),
                  min_solution_delta=(achieved_min-target).tolist())
    if gap >= 0:
        # A null-J direction adds exactly the required quadratic S amount.
        trial = np.zeros(len(j))
        trial[np.argmin(np.diag(gram))] = 1.
        null = trial-(j@trial/dual_norm)*dual
        norm = float(null@gram@null)
        coeff = minimum+np.sqrt(gap/norm)*null
        achieved = moment_vector(U+coeff@modes, E, X, weights)
        result.update(coefficients=coeff.tolist(),
                      final_delta=(achieved-target).tolist(),
                      max_abs_five_moment_defect=float(np.max(np.abs(
                          achieved-target))),
                      max_abs_corrected_U=float(np.max(np.abs(U+coeff@modes))))
    return result


def run():
    tau = .5*2**(-5.5)
    base = make_field(16, 2.)
    field = RadialMomentStep(base, -.5)
    X, weights = grid()
    rows = [solve_slice(base, field, X, weights, eta, tau)
            for eta in (.2, .3)]
    for row in rows:
        print(json.dumps(row), flush=True)
    report = dict(tau=tau, x_out=3.5, quadrature_per_piece=96,
                  u_intervals=U_INTERVALS, rows=rows,
                  scope='Exact finite-dimensional J/S feasibility on two normalized eta slices. Every added U basis has zero M integral and E is unchanged. This is not a smooth eta/time coefficient field, a physical solenoidal lift, a continuous cone, or PDE acceptance.',
                  accepted=False)
    (ROOT/'meridional_moment_repair.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
