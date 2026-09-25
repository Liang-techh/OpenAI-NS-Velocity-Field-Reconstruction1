"""Hilbert-space lower bound for any U-only repair on a radial interval.

With E fixed, matching M and J imposes two linear constraints on U.
The least possible integral of U**2 follows from projection onto the
span of 1 and H=sqrt(2X)E. If this minimum exceeds the baseline U**2
budget, no meridional-only repair supported on that interval can also
match S. The bound is evaluated by finite quadrature, not interval proof.
"""
import json

import numpy as np
from numpy.polynomial.legendre import leggauss

from extended_relaxed_cone_screen import profile
from high_frequency_shear_screen import make_field
from radial_continuation import ROOT
from radial_moment_step import RadialMomentStep


STARTS = (.1, .3, .5, .75, .95, .975, .99, 1., 1.005,
          1.01, 1.025, 1.05)
END = 3.


def grid(order=96):
    edges = sorted(set((0., 3.5, 3/32, 1.5, 1.75, 3., *STARTS)))
    g, w = leggauss(order)
    xs = np.concatenate([(a+b)/2+(b-a)*g/2
                         for a, b in zip(edges[:-1], edges[1:])])
    ws = np.concatenate([(b-a)*w/2
                         for a, b in zip(edges[:-1], edges[1:])])
    return xs, ws


def run():
    tau = .5*2**(-5.5)
    base = make_field(16, 2.)
    changed = RadialMomentStep(base, -.5)
    X, weights = grid()
    rows = []
    for eta in (.2, .3):
        U0, E = profile(base, X, eta, tau)
        U, _ = profile(changed, X, eta, tau)
        H = np.sqrt(2*X)*E
        target_M = float(weights@U0)
        target_J = float(weights@(U0*H))
        target_U2 = float(weights@(U0**2))
        for start in STARTS:
            inside = (X > start) & (X < END)
            outside = ~inside
            wi = weights[inside]
            Hi = H[inside]
            gram = np.array([[sum(wi), wi@Hi],
                             [wi@Hi, wi@(Hi**2)]])
            need = np.array([target_M-weights[outside]@U[outside],
                             target_J-weights[outside]@(
                                 U[outside]*H[outside])])
            min_inside_U2 = float(need@np.linalg.solve(gram, need))
            min_total_U2 = (min_inside_U2
                            +float(weights[outside]@(U[outside]**2)))
            row = dict(eta=eta, support_start=start, support_end=END,
                       gram_condition=float(np.linalg.cond(gram)),
                       required_M_J=need.tolist(),
                       target_U2=target_U2,
                       min_total_U2=min_total_U2,
                       energy_slack=target_U2-min_total_U2,
                       U_only_possible_by_bound=bool(
                           target_U2-min_total_U2 >= -1e-8))
            rows.append(row)
            print(json.dumps(row), flush=True)
    report = dict(tau=tau, quadrature_per_piece=96, rows=rows,
                  scope='Finite-quadrature projection bound with E fixed, prescribed outgoing M/J/S and U change supported on indicated X interval. Negative slack rules out U-only repair within numerical accuracy, not E/profile redesign. Positive slack is necessary but not a smooth physical construction.',
                  accepted=False)
    (ROOT/'meridional_moment_lower_bound.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
