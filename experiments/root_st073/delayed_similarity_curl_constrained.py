"""Fit exact-curl similarity modes with exact two-slice moment constraints.

This is a local momentum screen, not a relaxed-cone or PDE acceptance test.
"""

import json

import numpy as np
from scipy.optimize import minimize

from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from delayed_momentum_tangent_screen import nodes, residual, stats
from delayed_similarity_curl_screen import (
    CurlPatchedLift, INTERVALS, SimilarityCurlMode, nonlinear_residual,
    screen_data,
)
from delayed_taper_capacity_screen import grid
from extended_relaxed_cone_screen import profile
from moment_shear_slice_repair import moment_vector
from radial_continuation import ROOT


ETA_INTERVALS = ((.15, .275), (.225, .38))
BASE_NAME = 'delayed005_rise146_degree31.json'


def moment_model(base, modes, tau):
    X, weights = grid(order=48)
    rows = []
    for eta in (.2, .3):
        U, E = profile(base, X, eta, tau)
        points = base.compact.joined.inner.from_similarity(
            X, np.full(len(X), eta), tau)
        q = tau/(1-eta**2)
        factor = q**(.5+base.heat.h)/np.sqrt(base.nu)
        B = np.column_stack([factor*mode.fields(points, tau)[0][:, 2]
                             for mode in modes])
        H = np.sqrt(2*X)*E
        rows.append(dict(eta=eta, J=weights@(H[:, None]*B),
                         S=2*weights@(U[:, None]*B),
                         Q=B.T@(weights[:, None]*B),
                         base=moment_vector(U, E, X, weights)))
    return rows


def run():
    base = CoupledMomentPhysicalLift(slice_filename=BASE_NAME)
    modes = [SimilarityCurlMode(base, radial, axial)
             for axial in ETA_INTERVALS for radial in INTERVALS]
    taus = (.5*2**(-5.5), .5*2**(-5.25))
    train = [screen_data(base, modes, (1.01, 1.015, 1.02),
                         (.2, .3), tau) for tau in taus]
    R0 = np.concatenate([item['base_residual'] for item in train])
    C = np.concatenate([item['columns'] for item in train])
    R_scale = np.linalg.norm(R0)
    coef_scale = R_scale/np.linalg.norm(C.reshape(-1, len(modes)), axis=0)
    moments = moment_model(base, modes, taus[0])

    def coefficients(z):
        return coef_scale*z

    def objective(z):
        c = coefficients(z)
        R = np.concatenate([nonlinear_residual(item, c) for item in train])
        return np.sum(R**2)/R_scale**2 + .0005*np.sum(z**2)

    max_scale = stats(R0)['max']

    def peak_objective(z):
        c = coefficients(z)
        R = np.concatenate([nonlinear_residual(item, c) for item in train])
        norms = np.linalg.norm(R, axis=1)/max_scale
        return np.mean(norms**8)**.25 + .0005*np.sum(z**2)

    def constraints(z):
        c = coefficients(z)
        return np.array([value for row in moments for value in (
            row['J']@c, row['S']@c+c@row['Q']@c)])

    equality = [dict(type='eq', fun=constraints)]
    rms_result = minimize(objective, np.zeros(len(modes)), method='SLSQP',
                          constraints=equality,
                          options=dict(ftol=1e-12, maxiter=500))
    peak_result = minimize(peak_objective, rms_result.x, method='SLSQP',
                           constraints=equality,
                           options=dict(ftol=1e-12, maxiter=500))
    fits = []
    for name, fit in (('rms', rms_result), ('peak', peak_result)):
        fit_residual = np.concatenate([
            nonlinear_residual(item, coefficients(fit.x)) for item in train])
        fits.append(dict(name=name, result=fit, metrics=stats(fit_residual)))
    chosen = min(fits, key=lambda item: item['metrics']['max'])
    result = chosen['result']
    c = coefficients(result.x)
    patched = CurlPatchedLift(base, modes, c)
    tau_holdout = .5*2**(-5.4)
    points, Xh, etah = nodes(base, (1.0125, 1.0175, 1.025),
                             (.22, .28, .32), tau_holdout)
    before, _ = residual(base, points, tau_holdout)
    after, divergence = residual(patched, points, tau_holdout)
    at_train = np.concatenate([nonlinear_residual(item, c) for item in train])
    report = dict(base_slice=BASE_NAME,
                  radial_intervals=INTERVALS, eta_intervals=ETA_INTERVALS,
                  selected_objective=chosen['name'],
                  candidate_train_metrics=[dict(name=item['name'],
                                                metrics=item['metrics'],
                                                success=bool(item['result'].success))
                                           for item in fits],
                  coefficients=c.tolist(), optimizer_success=bool(result.success),
                  optimizer_message=result.message,
                  objective=float(result.fun),
                  moment_constraint_residual=constraints(result.x).tolist(),
                  max_abs_moment_constraint=float(np.max(np.abs(
                      constraints(result.x)))),
                  train_before=stats(R0), train_after=stats(at_train),
                  holdout=dict(tau=tau_holdout, X=Xh.tolist(),
                               eta=etah.tolist(), before=stats(before),
                               after=stats(after),
                               max_abs_divergence=float(np.max(np.abs(
                                   divergence)))),
                  scope='Six compact exact-curl similarity modes. Four '
                        'equalities enforce J and S at two reference '
                        'slices; M is exact by compact streamfunction, '
                        'and H/Cp are unchanged because E is unchanged. '
                        'Local full-momentum train and disjoint space/time '
                        'holdout only; no global cone or PDE acceptance.',
                  accepted=False)
    (ROOT/'delayed_similarity_curl_constrained.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    run()
