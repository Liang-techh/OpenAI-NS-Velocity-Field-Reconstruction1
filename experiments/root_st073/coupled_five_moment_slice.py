"""Construct fixed-eta E/U profiles that restore all five outgoing moments.

The E coefficients come from the positive-swirl I/Cp solve. A C3 radial
window times twelve Legendre modes supplies U corrections after X=1. The
constrained quadratic solve enforces M/J and then S exactly when its
finite-dimensional kinetic minimum lies below the target. No physical
eta/time-dependent solenoidal field is asserted by this slice result.
"""
import json

import numpy as np
from numpy.polynomial.legendre import legval

from azimuthal_capacity_optimize import INTERVALS, grid
from extended_relaxed_cone_screen import profile
from high_frequency_shear_screen import make_field
from moment_shear_slice_repair import moment_vector
from paper_moment_bridge import bump
from radial_continuation import ROOT
from radial_moment_step import RadialMomentStep, septic_step


def correction_modes(X):
    up = septic_step((X-1.)/.02)[0]
    down = septic_step((X-2.98)/.02)[0]
    window = up*(1-down)
    argument = np.clip(X-2., -1., 1.)
    return np.array([window*legval(argument, [0.]*degree+[1.])
                     for degree in range(12)])


def construct(base, changed, X, weights, eta, tau, e_coeff):
    U0, E0 = profile(base, X, eta, tau)
    U, _ = profile(changed, X, eta, tau)
    eb = np.array([bump(X, *interval)[0] for interval in INTERVALS])
    E = E0+np.asarray(e_coeff)@eb
    target = moment_vector(U0, E0, X, weights)
    current = moment_vector(U, E, X, weights)
    B = correction_modes(X)
    H = np.sqrt(2*X)*E
    # Orthonormalize in the quadrature L2 norm before the KKT solve.
    orth, factor = np.linalg.qr((B*np.sqrt(weights)).T)
    Q = (orth.T/np.sqrt(weights))
    gram = (Q*weights)@Q.T
    linear_u = Q@(weights*U)
    constraints = np.array([Q@weights, Q@(weights*H)])
    required = np.array([target[0]-current[0],
                         target[2]-weights@(U*H)])
    inverse_constraint = np.linalg.solve(gram, constraints.T)
    dual = constraints@inverse_constraint
    unconstrained = -np.linalg.solve(gram, linear_u)
    minimum = (unconstrained+inverse_constraint@np.linalg.solve(
        dual, required-constraints@unconstrained))
    min_data = moment_vector(U+minimum@Q, E, X, weights)
    gap = float(target[3]-min_data[3])
    result = dict(eta=eta, e_coefficients=list(e_coeff),
                  minimum_S=float(min_data[3]), target_S=float(target[3]),
                  S_slack=gap,
                  original_basis_condition=float(np.linalg.cond(factor)),
                  min_relative_E=float(np.min(E/E0)),
                  finite_basis_feasible=bool(gap >= -1e-8))
    if gap < 0:
        return result
    options = []
    for index in range(len(Q)):
        trial = np.eye(len(Q))[index]
        null = trial-inverse_constraint@np.linalg.solve(
            dual, constraints@trial)
        norm = float(null@gram@null)
        if norm < 1e-12:
            continue
        for sign in (-1., 1.):
            coeff = minimum+sign*np.sqrt(gap/norm)*null
            U_new = U+coeff@Q
            options.append((float(np.max(np.abs(U_new))), coeff, U_new))
    maximum, coefficients, U_new = min(options, key=lambda row: row[0])
    achieved = moment_vector(U_new, E, X, weights)
    result.update(u_coefficients=np.linalg.solve(factor,
                  coefficients).tolist(),
                  final_delta=(achieved-target).tolist(),
                  max_abs_five_moment_defect=float(np.max(np.abs(
                      achieved-target))),
                  max_abs_corrected_U=maximum,
                  five_moments_restored=bool(np.max(np.abs(
                      achieved-target)) < 1e-6))
    return result


def run():
    tau = .5*2**(-5.5)
    base = make_field(16, 2.)
    changed = RadialMomentStep(base, -.5)
    X, weights = grid()
    source = json.loads((ROOT/'azimuthal_capacity_optimize.json').read_text())
    rows = []
    for prior in source['rows']:
        eta = prior['eta']
        for variant, key in (('minimum_norm', 'coefficients'),
                             ('maximum_slack', 'maximum_slack_coefficients')):
            row = dict(variant=variant, **construct(
                base, changed, X, weights, eta, tau, prior[key]))
            rows.append(row)
            print(json.dumps(row), flush=True)
    report = dict(tau=tau, x_out=3.5, quadrature_per_piece=96,
                  e_intervals=INTERVALS, u_degree=11, rows=rows,
                  scope='Fixed-eta constructive five-moment solve after the X=1 cone point. Does not establish coefficient smoothness in eta/time, a physical divergence-free lift, a continuous cone, or full NS momentum acceptance.',
                  accepted=False)
    (ROOT/'coupled_five_moment_slice.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
