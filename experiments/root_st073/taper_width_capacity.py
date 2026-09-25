"""Optimize exact finite-basis moment capacity against U taper width."""
import json

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import minimize

from azimuthal_capacity_optimize import INTERVALS, RELATIVE_SWIRL_FLOOR
from coupled_five_moment_slice import correction_modes
from extended_relaxed_cone_screen import profile
from high_frequency_shear_screen import make_field
from moment_shear_slice_repair import moment_vector
from paper_moment_bridge import bump
from radial_continuation import ROOT
from radial_moment_step import RadialMomentStep


WIDTHS = (.02, .05, .10, .20)


def grid(width, order=96):
    edges = sorted(set((0., 3.5, 3/32, 1., 1.+width,
                        1.5, 1.75, 3.-width, 3.,
                        *(v for interval in INTERVALS for v in interval))))
    g, w = leggauss(order)
    X = np.concatenate([(a+b)/2+(b-a)*g/2
                        for a, b in zip(edges[:-1], edges[1:])])
    weights = np.concatenate([(b-a)*w/2
                              for a, b in zip(edges[:-1], edges[1:])])
    return X, weights


def optimize_slice(base, changed, X, weights, eta, tau, width, prior):
    U0, E0 = profile(base, X, eta, tau)
    U, _ = profile(changed, X, eta, tau)
    target = moment_vector(U0, E0, X, weights)
    eb = np.array([bump(X, *interval)[0] for interval in INTERVALS])
    lower = np.array([np.max((RELATIVE_SWIRL_FLOOR-1)*E0[b > 1e-10]
                             /b[b > 1e-10])+1e-5 for b in eb])
    B = correction_modes(X, width)
    orth, factor = np.linalg.qr((B*np.sqrt(weights)).T)
    Q = orth.T/np.sqrt(weights)
    linear_u = Q@(weights*U)
    zero = -linear_u
    mass_row = Q@weights
    M_current = float(weights@U)

    def metrics(c):
        E = E0+c@eb
        H = np.sqrt(2*X)*E
        I_Cp = np.array([weights@H-target[1],
                         weights@(E**2/(2*X))-target[4]])
        J_current = float(weights@(U*H))
        D = np.array([mass_row, Q@(weights*H)])
        dual = D@D.T
        need = np.array([target[0]-M_current, target[2]-J_current])
        coeff = zero+D.T@np.linalg.solve(dual, need-D@zero)
        U_min = U+coeff@Q
        S_min = float(weights@(U_min**2-E**2/2))
        return I_Cp, float(target[3]-S_min), float(np.min(E/E0))

    starts = [np.zeros(5), np.asarray(prior, float)]
    rng = np.random.default_rng(73076)
    starts += [np.maximum(np.asarray(prior)+rng.normal(size=5)*.15,
                          lower+1e-4) for _ in range(3)]
    bounds = list(zip(lower, [3.]*5))
    equality = dict(type='eq', fun=lambda c: metrics(c)[0])
    fits = [minimize(lambda c: -metrics(c)[1]+.0005*(c@c),
                     np.clip(start, lower+1e-6, 3.),
                     method='SLSQP', bounds=bounds,
                     constraints=[equality],
                     options=dict(maxiter=300, ftol=1e-11))
            for start in starts]
    eligible = [f for f in fits if np.max(np.abs(metrics(f.x)[0])) < 1e-7]
    fit = max(eligible, key=lambda f: metrics(f.x)[1]) if eligible else min(
        fits, key=lambda f: np.linalg.norm(metrics(f.x)[0]))
    defects, slack, min_ratio = metrics(fit.x)
    return dict(width=width, eta=eta, coefficients=fit.x.tolist(),
                I_Cp_defect=defects.tolist(), finite_basis_S_slack=slack,
                min_relative_E=min_ratio,
                original_basis_condition=float(np.linalg.cond(factor)),
                optimizer_success=bool(fit.success),
                feasible=bool(np.max(np.abs(defects)) < 1e-7
                              and slack > 0 and min_ratio > .10))


def run():
    tau = .5*2**(-5.5)
    base = make_field(16, 2.)
    changed = RadialMomentStep(base, -.5)
    prior = json.loads((ROOT/'azimuthal_capacity_optimize.json').read_text())
    prior_by_eta = {row['eta']: row['maximum_slack_coefficients']
                    for row in prior['rows']}
    rows = []
    for width in WIDTHS:
        X, weights = grid(width)
        for eta in (.2, .3):
            row = optimize_slice(base, changed, X, weights, eta, tau,
                                 width, prior_by_eta[eta])
            rows.append(row)
            print(json.dumps(row), flush=True)
    report = dict(tau=tau, widths=WIDTHS, quadrature_per_piece=96,
                  relative_swirl_floor_training=RELATIVE_SWIRL_FLOOR,
                  upper_e_coefficient=3., rows=rows,
                  scope='Fixed-slice finite-basis S capacity jointly optimized with positive E and exact I/Cp for each taper width. A positive slack permits an exact U moment correction in this basis, but says nothing about physical momentum, eta/time continuity, pressure or full-domain gates.',
                  accepted=False)
    (ROOT/'taper_width_capacity.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
