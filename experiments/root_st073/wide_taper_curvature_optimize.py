"""Minimize radial U curvature within the exact five-moment slice manifold."""
import json

import numpy as np
from scipy.optimize import minimize

from azimuthal_capacity_optimize import INTERVALS
from coupled_five_moment_slice import correction_modes
from extended_relaxed_cone_screen import profile
from high_frequency_shear_screen import make_field
from moment_shear_slice_repair import moment_vector
from paper_moment_bridge import bump
from radial_continuation import ROOT
from radial_moment_step import RadialMomentStep
from taper_width_capacity import grid


def optimize_slice(base, changed, X, weights, eta, tau, row):
    width, degree = row['taper_width'], row['u_degree']
    U0, E0 = profile(base, X, eta, tau)
    U, _ = profile(changed, X, eta, tau)
    eb = np.array([bump(X, *interval)[0] for interval in INTERVALS])
    E = E0+np.asarray(row['e_coefficients'])@eb
    target = moment_vector(U0, E0, X, weights)
    H = np.sqrt(2*X)*E
    B = correction_modes(X, width, degree)
    dense_X = np.linspace(1., 3., 2001)
    dx = dense_X[1]-dense_X[0]
    dense_B = correction_modes(dense_X, width, degree)
    Bxx = np.gradient(np.gradient(dense_B, dx, axis=1), dx, axis=1)
    stiffness = (Bxx*dx)@Bxx.T
    stiffness += np.eye(len(B))*1e-7
    scale = float(np.linalg.eigvalsh(stiffness)[-1])
    K = stiffness/scale
    start = np.asarray(row['u_coefficients'])
    target_linear = np.array([target[0]-weights@U,
                              target[2]-weights@(U*H)])
    linear = np.array([B@weights, B@(weights*H)])

    def constraints(c):
        u_new = U+c@B
        return np.r_[linear@c-target_linear,
                     weights@(u_new**2-E**2/2)-target[3]]

    def jac(c):
        u_new = U+c@B
        return np.vstack((linear, 2*B@(weights*u_new)))

    fit = minimize(lambda c: float(c@K@c), start,
                   jac=lambda c: 2*K@c, method='SLSQP',
                   constraints=[dict(type='eq', fun=constraints, jac=jac)],
                   options=dict(maxiter=1000, ftol=1e-12))
    candidates = [start, fit.x]
    feasible = [c for c in candidates if np.max(np.abs(constraints(c))) < 1e-7]
    selected = min(feasible, key=lambda c: c@K@c) if feasible else start
    final = moment_vector(U+selected@B, E, X, weights)-target
    return dict(variant='curvature_optimized', eta=eta,
                taper_width=width, u_degree=degree,
                e_coefficients=row['e_coefficients'],
                u_coefficients=selected.tolist(),
                original_curvature_proxy=float(start@stiffness@start),
                optimized_curvature_proxy=float(selected@stiffness@selected),
                curvature_ratio=float(selected@stiffness@selected
                                      /(start@stiffness@start)),
                final_delta=final.tolist(),
                max_abs_five_moment_defect=float(np.max(np.abs(final))),
                optimizer_success=bool(fit.success),
                five_moments_restored=bool(np.max(np.abs(final)) < 1e-6))


def run():
    tau = .5*2**(-5.5)
    base = make_field(16, 2.)
    changed = RadialMomentStep(base, -.5)
    X, weights = grid(.05)
    source = json.loads((ROOT/'wide_taper_five_moment_slice.json').read_text())
    rows = [optimize_slice(base, changed, X, weights, row['eta'], tau, row)
            for row in source['rows']]
    for row in rows:
        print(json.dumps(row), flush=True)
    report = dict(tau=tau, x_out=3.5, taper_width=.05, u_degree=19,
                  quadrature_per_piece=96, rows=rows,
                  scope='Fixed-slice radial second-derivative minimization subject to exact M/J/S, with E and I/Cp fixed. Proxy improvement does not imply physical momentum improvement; no eta/time or pressure repair.',
                  accepted=False)
    (ROOT/'wide_taper_curvature_optimize.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
