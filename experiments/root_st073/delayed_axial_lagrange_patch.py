"""Couple three exact-curl/toroidal slice repairs through an axial basis.

The Lagrange profiles vanish at eta=.2,.3 and isolate three interior
collocation slices. All poloidal modes are curls of physical streamfunctions.
The intervening eta and independent momentum samples decide usefulness.
"""

import json

import numpy as np
from scipy.optimize import minimize

from azimuthal_capacity_optimize import (
    INTERVALS as SWIRL_INTERVALS, grid, optimize_slice,
)
from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from delayed_midband_moment_patch import (
    BASE_NAME, ETA_INTERVAL,
    SimilaritySwirlMode, mode_profile,
)
from delayed_momentum_tangent_screen import nodes, residual, stats
from delayed_similarity_curl_screen import CurlPatchedLift, SimilarityCurlMode
from extended_relaxed_cone_screen import profile
from high_frequency_shear_screen import make_field
from moment_shear_slice_repair import moment_vector
from paper_moment_bridge import bump
from radial_continuation import ROOT


COLLOCATION_ETA = (.225, .25, .275)
LAGRANGE = ((0., -1., 2.), (1., 0., -4.), (0., 1., 2.))
POLOIDAL_INTERVALS = ((1.005, 1.75), (1.1, 2.5), (1.5, 3.),
                      (1.005, 1.35), (1.3, 2.), (1.8, 2.5),
                      (2.35, 3.), (1.005, 1.5), (1.2, 1.7),
                      (1.4, 1.9), (1.6, 2.1), (1.8, 2.3),
                      (2., 2.5), (2.2, 2.7), (2.4, 3.),
                      (1.005, 2.), (1.2, 2.2), (1.4, 2.4),
                      (1.6, 2.6), (1.8, 3.))


def axial_shape(index):
    polynomial = np.asarray(LAGRANGE[index])
    node = COLLOCATION_ETA[index]
    normalization = bump(np.array([node]), *ETA_INTERVAL)[0][0]

    def shape(eta):
        eta = np.asarray(eta, float)
        s = (eta-.25)/.05
        value = np.polynomial.polynomial.polyval(s, polynomial)
        derivative = np.polynomial.polynomial.polyval(
            s, np.polynomial.polynomial.polyder(polynomial))/.05
        envelope, envelope_d = bump(eta, *ETA_INTERVAL)
        return (envelope*value/normalization,
                (envelope_d*value+envelope*derivative)/normalization)

    return shape


def fit_slice(base, target, modes, X, weights, eta, tau):
    U, E = profile(base, X, eta, tau)
    U0, E0 = profile(target, X, eta, tau)
    desired = moment_vector(U0, E0, X, weights)
    bu, be = [np.column_stack(items) for items in zip(*[
        mode_profile(mode, base, X, eta, tau) for mode in modes])]
    eb = be[:, len(POLOIDAL_INTERVALS):]
    optimized = optimize_slice(target, base.base, X, weights, eta, tau)
    current_e = base.coefficients(eta)[0]
    e_options = [np.asarray(optimized[key])-current_e
                 for key in ('maximum_slack_coefficients', 'coefficients')]
    def moments(c):
        return moment_vector(U+bu@c, E+be@c, X, weights)
    candidates = []
    ub = bu[:, :len(POLOIDAL_INTERVALS)]
    gram_u = ub.T@(weights[:, None]*ub)
    inverse_u = np.linalg.pinv(gram_u, rcond=1e-11)

    def capacity(e_coeff):
        E_new = E+eb@e_coeff
        H_new = np.sqrt(2*X)*E_new
        a_u = ub.T@(weights*H_new)
        g_u = ub.T@(weights*U)
        need_J = desired[2]-weights@(U*H_new)
        multiplier_u = (-a_u@inverse_u@g_u-need_J)/(
            a_u@inverse_u@a_u)
        minimum_u = -inverse_u@(g_u+multiplier_u*a_u)
        required_u2 = desired[3]+weights@(E_new**2/2)
        achieved_min_u2 = weights@((U+ub@minimum_u)**2)
        return required_u2-achieved_min_u2, minimum_u, a_u, E_new

    def e_equalities(e_coeff):
        E_new = E+eb@e_coeff
        return np.array([weights@(np.sqrt(2*X)*E_new)-desired[1],
                         weights@(E_new**2/(2*X))-desired[4]])/.01

    fit_e = minimize(lambda c: -capacity(c)[0]+.001*(c@c),
                     e_options[0], method='SLSQP',
                     constraints=[dict(type='eq', fun=e_equalities),
                                  dict(type='ineq', fun=lambda c:
                                       np.min((E+eb@c)/E0)-.110001)],
                     options=dict(maxiter=500, ftol=1e-12))
    if (np.all(np.isfinite(fit_e.x)) and
            np.max(np.abs(e_equalities(fit_e.x))) < 1e-7 and
            np.min((E+eb@fit_e.x)/E0) >= .11):
        e_options.append(fit_e.x)
    print(json.dumps(dict(eta=eta,
                          optimized_e_success=bool(fit_e.success),
                          optimized_e_slack=float(capacity(fit_e.x)[0]),
                          optimized_e_constraint=e_equalities(fit_e.x).tolist())),
          flush=True)
    for e_coeff in e_options:
        u_slack, minimum_u, a_u, E_new = capacity(e_coeff)
        print(json.dumps(dict(eta=eta,
                              optimized_energy_slack=optimized['energy_slack'],
                              e_coefficients=e_coeff.tolist(),
                              u_slack=float(u_slack),
                              min_relative_E=float(np.min(E_new/E0)))),
              flush=True)
        if u_slack < -1e-8 or np.min(E_new/E0) < .11:
            continue
        _, _, vh = np.linalg.svd(a_u[None, :])
        for direction in vh[1:]:
            step_u = np.sqrt(max(0., u_slack)/(
                direction@gram_u@direction))
            for sign in (-1., 1.):
                u_coeff = minimum_u+sign*step_u*direction
                c = np.r_[u_coeff, e_coeff]
                defect = moments(c)-desired
                score = (np.max(np.abs(ub@u_coeff))
                         +np.max(np.abs(eb@e_coeff)))
                candidates.append((score, c, defect, u_slack,
                                   np.min(E_new/E0)))
    if not candidates:
        raise RuntimeError(f'eta={eta} no J/S solution with optimized E '
                           f'{optimized["energy_slack"]}')
    selected = min(candidates, key=lambda item: item[0])
    _, c, defect, u_slack, relative_e = selected
    success = bool(np.max(np.abs(defect[1:])) < 1e-7 and relative_e >= .11)
    return c, dict(eta=eta, success=success,
                   optimized_energy_slack=optimized['energy_slack'],
                   u_slack=float(u_slack),
                   min_relative_E=float(relative_e),
                   coefficient_norm=float(np.linalg.norm(c)),
                   before=(moments(np.zeros(len(modes)))-desired).tolist(),
                   after=defect.tolist())


def run():
    tau = .5*2**(-5.5)
    base = CoupledMomentPhysicalLift(slice_filename=BASE_NAME)
    target = make_field(16, 2.)
    X, weights = grid(order=48)
    mode_groups = []
    for index in range(len(COLLOCATION_ETA)):
        shape = axial_shape(index)
        mode_groups.append(
            [SimilarityCurlMode(base, radial, ETA_INTERVAL, shape)
             for radial in POLOIDAL_INTERVALS]
            +[SimilaritySwirlMode(base, radial, shape)
              for radial in SWIRL_INTERVALS])
    slices = [fit_slice(base, target, group, X, weights, eta, tau)
              for eta, group in zip(COLLOCATION_ETA, mode_groups)]
    if any(not item[1]['success'] or
           not np.all(np.isfinite(item[0])) for item in slices):
        raise RuntimeError('Axial slice solve failed: '+json.dumps(
            [item[1] for item in slices]))
    modes = [mode for group in mode_groups for mode in group]
    coefficients = np.concatenate([item[0] for item in slices])
    patched = CurlPatchedLift(base, modes, coefficients)

    eta_samples = (.2, .2125, .225, .2375, .25,
                   .2625, .275, .2875, .3)
    rows = []
    for eta in eta_samples:
        U0, E0 = profile(target, X, eta, tau)
        U_base, E_base = profile(base, X, eta, tau)
        U, E = profile(patched, X, eta, tau)
        desired = moment_vector(U0, E0, X, weights)
        rows.append(dict(eta=eta,
                         before=(moment_vector(U_base, E_base, X, weights)
                                 -desired).tolist(),
                         after=(moment_vector(U, E, X, weights)
                                -desired).tolist(),
                         min_relative_E=float(np.min(E/E0))))
    tau_holdout = .5*2**(-5.4)
    points, Xh, etah = nodes(base, (1.015, 1.08, 1.3),
                             (.225, .2375, .25, .2625, .275), tau_holdout)
    before, _ = residual(base, points, tau_holdout)
    after, divergence = residual(patched, points, tau_holdout)
    report = dict(base_slice=BASE_NAME, tau=tau,
                  eta_interval=ETA_INTERVAL,
                  collocation_eta=COLLOCATION_ETA,
                  poloidal_intervals=POLOIDAL_INTERVALS,
                  swirl_intervals=SWIRL_INTERVALS,
                  coefficients=coefficients.tolist(),
                  slice_fits=[item[1] for item in slices], rows=rows,
                  holdout=dict(tau=tau_holdout, X=Xh.tolist(),
                               eta=etah.tolist(), before=stats(before),
                               after=stats(after),
                               max_abs_divergence=float(np.max(np.abs(
                                   divergence)))),
                  scope='Three interior axial moment collocation slices '
                        'with independent intervening eta and full '
                        'momentum holdout. Analytically divergence-free; '
                        'no interval, cone, volume norm, or PDE acceptance.',
                  accepted=False)
    (ROOT/'delayed_axial_lagrange_patch.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps(dict(slice_fits=report['slice_fits'],
                          rows=rows, holdout=report['holdout'])), flush=True)


if __name__ == '__main__':
    run()
