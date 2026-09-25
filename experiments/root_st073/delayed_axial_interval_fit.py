"""Fit the physical axial patch against an interval rather than three slices.

Fixed-slice U/E moments have analytic derivatives in the compact mode
coefficients. The fit includes interleaved eta points and the relative
swirl floor; separate eta and time points remain holdouts. This is still
only a mean-flow screen, not a continuous cone or PDE certificate.
"""

import json

import numpy as np
from scipy.optimize import least_squares

from azimuthal_capacity_optimize import grid
from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from delayed_axial_lagrange_patch import (
    COLLOCATION_ETA, POLOIDAL_INTERVALS, SWIRL_INTERVALS, axial_shape,
)
from delayed_midband_moment_patch import (
    BASE_NAME, ETA_INTERVAL, SimilaritySwirlMode,
)
from delayed_momentum_tangent_screen import nodes, residual, stats
from delayed_similarity_curl_screen import CurlPatchedLift, SimilarityCurlMode
from extended_relaxed_cone_screen import profile
from high_frequency_shear_screen import make_field
from moment_shear_slice_repair import moment_vector
from paper_moment_bridge import bump
from radial_continuation import ROOT


TRAIN_ETA = (.2125, .225, .2375, .25, .2625, .275, .2875)
HOLDOUT_ETA = (.20625, .21875, .24375, .26875, .29375)
MOMENT_SCALES = np.array([1e-5, 5e-4, 5e-4, 1e-5])
SWIRL_FLOOR = .11005


def basis_matrices(X, eta):
    u_radial = np.column_stack([bump(X, *interval)[1]
                                for interval in POLOIDAL_INTERVALS])
    e_radial = np.column_stack([bump(X, *interval)[0]
                                for interval in SWIRL_INTERVALS])
    pieces_u, pieces_e = [], []
    for index in range(len(COLLOCATION_ETA)):
        amplitude = axial_shape(index)(np.array([eta]))[0][0]
        pieces_u.append(np.column_stack((amplitude*u_radial,
                                         np.zeros((len(X), e_radial.shape[1])))))
        pieces_e.append(np.column_stack((np.zeros((len(X), u_radial.shape[1])),
                                         amplitude*e_radial)))
    return np.column_stack(pieces_u), np.column_stack(pieces_e)


def physical_modes(base):
    modes = []
    for index in range(len(COLLOCATION_ETA)):
        shape = axial_shape(index)
        modes.extend(SimilarityCurlMode(base, radial, ETA_INTERVAL, shape)
                     for radial in POLOIDAL_INTERVALS)
        modes.extend(SimilaritySwirlMode(base, radial, shape)
                     for radial in SWIRL_INTERVALS)
    return modes


def run():
    tau = .5*2**(-5.5)
    base = CoupledMomentPhysicalLift(slice_filename=BASE_NAME)
    target = make_field(16, 2.)
    X, weights = grid(order=48)
    prior = json.loads((ROOT/'delayed_axial_lagrange_patch.json').read_text())
    initial = np.asarray(prior['coefficients'])
    rows = []
    for eta in TRAIN_ETA:
        U0, E0 = profile(target, X, eta, tau)
        U, E = profile(base, X, eta, tau)
        bu, be = basis_matrices(X, eta)
        rows.append(dict(eta=eta, U=U, E=E, U0=U0, E0=E0,
                         target=moment_vector(U0, E0, X, weights),
                         bu=bu, be=be))
    print('Precomputed interval profiles', flush=True)

    def evaluate(c, include_jac=False):
        values, derivatives = [], []
        for row in rows:
            U = row['U']+row['bu']@c
            E = row['E']+row['be']@c
            defect = (moment_vector(U, E, X, weights)
                      -row['target'])[1:]/MOMENT_SCALES
            values.extend(defect)
            if include_jac:
                H_weight = weights*np.sqrt(2*X)
                derivatives.extend((
                    H_weight@row['be']/MOMENT_SCALES[0],
                    (H_weight*E)@row['bu']/MOMENT_SCALES[1]
                    +(H_weight*U)@row['be']/MOMENT_SCALES[1],
                    (2*weights*U)@row['bu']/MOMENT_SCALES[2]
                    -(weights*E)@row['be']/MOMENT_SCALES[2],
                    (weights*E/X)@row['be']/MOMENT_SCALES[3]))
        for row in rows:
            E = row['E']+row['be']@c
            ratio = E/row['E0']
            index = int(np.argmin(ratio))
            gap = SWIRL_FLOOR-ratio[index]
            values.append(max(0., gap)/1e-4)
            if include_jac:
                derivatives.append((-row['be'][index]/row['E0'][index]/1e-4
                                    if gap > 0 else np.zeros(len(c))))
        values.extend(.001*c)
        if include_jac:
            derivatives.extend(.001*np.eye(len(c)))
            return np.asarray(derivatives)
        return np.asarray(values)

    fit = least_squares(evaluate, initial, jac=lambda c: evaluate(c, True),
                        method='trf', max_nfev=200, ftol=1e-11,
                        xtol=1e-11, gtol=1e-11)
    c = fit.x
    train = []
    for row in rows:
        U, E = row['U']+row['bu']@c, row['E']+row['be']@c
        train.append(dict(eta=row['eta'],
                          before=(moment_vector(row['U'], row['E'], X, weights)
                                  -row['target']).tolist(),
                          after=(moment_vector(U, E, X, weights)
                                 -row['target']).tolist(),
                          min_relative_E=float(np.min(E/row['E0']))))
    print(json.dumps(dict(fit_success=bool(fit.success),
                          fit_cost=float(fit.cost), train=train)), flush=True)

    modes = physical_modes(base)
    patched = CurlPatchedLift(base, modes, c)
    holdout = []
    for eta in HOLDOUT_ETA:
        U0, E0 = profile(target, X, eta, tau)
        Ubase, Ebase = profile(base, X, eta, tau)
        U, E = profile(patched, X, eta, tau)
        desired = moment_vector(U0, E0, X, weights)
        holdout.append(dict(eta=eta,
                            before=(moment_vector(Ubase, Ebase, X, weights)
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
                  source='delayed_axial_lagrange_patch.json',
                  train_eta=TRAIN_ETA, holdout_eta=HOLDOUT_ETA,
                  coefficients=c.tolist(), fit_success=bool(fit.success),
                  fit_message=fit.message, fit_cost=float(fit.cost),
                  train=train, holdout=holdout,
                  momentum_holdout=dict(tau=tau_holdout, X=Xh.tolist(),
                                        eta=etah.tolist(),
                                        before=stats(before), after=stats(after),
                                        max_abs_divergence=float(np.max(np.abs(
                                            divergence)))),
                  scope='Seven-eta joint moment/floor least squares with '
                        'five distinct eta holdouts and independent '
                        'full-momentum space/time holdout. No continuous '
                        'cone, volume norm, or PDE acceptance.',
                  accepted=False)
    (ROOT/'delayed_axial_interval_fit.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps(dict(holdout=holdout,
                          momentum_holdout=report['momentum_holdout'])),
          flush=True)


if __name__ == '__main__':
    run()
