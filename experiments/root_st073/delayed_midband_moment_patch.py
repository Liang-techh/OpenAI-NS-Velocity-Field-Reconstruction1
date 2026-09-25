"""Screen a compact physical patch that closes the missing middle slice.

The poloidal modes are Stokes-streamfunction curls. Axisymmetric swirl modes
also have zero divergence. This is a three-slice experiment, not continuous
moment closure or a Navier--Stokes acceptance test.
"""

import json

import numpy as np
from scipy.optimize import minimize

from azimuthal_capacity_optimize import grid
from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from delayed_momentum_tangent_screen import nodes, residual, stats
from delayed_similarity_curl_screen import CurlPatchedLift, SimilarityCurlMode
from extended_relaxed_cone_screen import profile
from high_frequency_shear_screen import make_field
from joined_field import coordinates
from moment_shear_slice_repair import moment_vector
from paper_moment_bridge import bump
from radial_continuation import ROOT


BASE_NAME = 'delayed005_rise146_degree31.json'
ETA_INTERVAL = (.2, .3)
POLOIDAL_INTERVALS = ((1.005, 1.75), (1.1, 2.5), (1.5, 3.))
SWIRL_INTERVALS = ((1.005, 1.75), (1.8, 2.98))


class SimilaritySwirlMode:
    def __init__(self, base, radial_interval, eta_shape=None,
                 eta_interval=ETA_INTERVAL):
        self.base = base
        self.nu = base.nu
        self.radial_interval = radial_interval
        self.eta_shape = eta_shape
        self.eta_interval = eta_interval

    def fields(self, points, tau):
        pts = np.asarray(points, float)
        ts = np.broadcast_to(np.asarray(tau, float), (len(pts),))
        radius = np.hypot(pts[:, 0], pts[:, 1])
        sn = np.sqrt(self.nu)
        co = coordinates(radius/sn, pts[:, 2]/sn, ts, self.base.heat.h)
        q, X, eta = (np.asarray(co[key]) for key in ('q', 'X', 'eta'))
        bx, _ = bump(X, *self.radial_interval)
        be, _ = (self.eta_shape(eta) if self.eta_shape is not None
                 else bump(eta, *self.eta_interval))
        ut = sn*q**(-self.base.A)*bx*be
        ca = np.divide(pts[:, 0], radius, out=np.ones_like(radius),
                       where=radius > 0)
        sa = np.divide(pts[:, 1], radius, out=np.zeros_like(radius),
                       where=radius > 0)
        return np.column_stack((-ut*sa, ut*ca, np.zeros(len(pts)))), \
            np.zeros(len(pts))


def mode_profile(mode, base, X, eta, tau):
    points = base.compact.joined.inner.from_similarity(
        X, np.full(len(X), eta), tau)
    velocity, _ = mode.fields(points, tau)
    q = tau/(1-eta**2)
    scale = q**(.5+base.heat.h)/np.sqrt(base.nu)
    return scale*velocity[:, 2], scale*velocity[:, 1]


def run():
    tau = .5*2**(-5.5)
    base = CoupledMomentPhysicalLift(slice_filename=BASE_NAME)
    target = make_field(16, 2.)
    modes = ([SimilarityCurlMode(base, interval, ETA_INTERVAL)
              for interval in POLOIDAL_INTERVALS]
             +[SimilaritySwirlMode(base, interval)
               for interval in SWIRL_INTERVALS])
    X, weights = grid(order=48)
    eta = .25
    U, E = profile(base, X, eta, tau)
    U0, E0 = profile(target, X, eta, tau)
    target_moments = moment_vector(U0, E0, X, weights)
    bu, be = [np.column_stack(items) for items in zip(*[
        mode_profile(mode, base, X, eta, tau) for mode in modes])]

    def moments(c):
        return moment_vector(U+bu@c, E+be@c, X, weights)

    def constraints(z):
        return (moments(scales*z)-target_moments)[1:]/.01

    baseline = moments(np.zeros(len(modes)))
    scales = .05/np.maximum(np.max(np.abs(bu)+np.abs(be), axis=0), 1e-10)
    jac = np.column_stack([(constraints(1e-4*np.eye(len(modes))[j])
                            -constraints(-1e-4*np.eye(len(modes))[j]))/2e-4
                           for j in range(len(modes))])
    initial = np.linalg.lstsq(jac, -constraints(np.zeros(len(modes))),
                              rcond=None)[0]
    fit = minimize(lambda z: float(z@z), initial, method='SLSQP',
                   constraints=[dict(type='eq', fun=constraints)],
                   options=dict(maxiter=500, ftol=1e-12))
    coefficients = scales*fit.x
    patched = CurlPatchedLift(base, modes, coefficients)
    tau_holdout = .5*2**(-5.4)
    points, Xh, etah = nodes(base, (1.015, 1.08, 1.3),
                             (.225, .25, .275), tau_holdout)
    before, _ = residual(base, points, tau_holdout)
    after, divergence = residual(patched, points, tau_holdout)
    rows = []
    for slice_eta in (.2, .225, .25, .275, .3):
        ub, eb = profile(base, X, slice_eta, tau)
        ut, et = profile(target, X, slice_eta, tau)
        us = np.column_stack([mode_profile(mode, base, X, slice_eta,
                                           tau)[0] for mode in modes])
        es = np.column_stack([mode_profile(mode, base, X, slice_eta,
                                           tau)[1] for mode in modes])
        desired = moment_vector(ut, et, X, weights)
        rows.append(dict(eta=slice_eta,
                         before=(moment_vector(ub, eb, X, weights)
                                 -desired).tolist(),
                         after=(moment_vector(ub+us@coefficients,
                                              eb+es@coefficients, X, weights)
                                -desired).tolist(),
                         min_relative_E=float(np.min((eb+es@coefficients)/et))))
    report = dict(base_slice=BASE_NAME, tau=tau,
                  eta_interval=ETA_INTERVAL,
                  poloidal_intervals=POLOIDAL_INTERVALS,
                  swirl_intervals=SWIRL_INTERVALS,
                  coefficients=coefficients.tolist(),
                  optimizer_success=bool(fit.success),
                  optimizer_message=fit.message,
                  midpoint_before=(baseline-target_moments).tolist(),
                  midpoint_after=(moments(coefficients)-target_moments).tolist(),
                  rows=rows,
                  holdout=dict(tau=tau_holdout, X=Xh.tolist(),
                               eta=etah.tolist(), before=stats(before),
                               after=stats(after),
                               max_abs_divergence=float(np.max(np.abs(
                                   divergence)))),
                  scope='Three-slice compact physical moment screen only. '
                        'Analytic solenoidality, but no continuous moment '
                        'identity, cone, volume norm, or PDE acceptance.',
                  accepted=False)
    (ROOT/'delayed_midband_moment_patch.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    run()
