"""Fit a small full-momentum step inside the exact fixed-slice moment manifold.

Only the eta=.3 U row moves. M and J are linear null constraints; a scalar
quadratic repair enforces S exactly for each trial. E, I, and Cp stay fixed.
The physical objective uses all three Cartesian momentum components.
"""
import json

import numpy as np

from azimuthal_capacity_optimize import INTERVALS
from coupled_five_moment_slice import correction_modes
from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from delayed_taper_capacity_screen import grid
from extended_relaxed_cone_screen import profile
from high_frequency_shear_screen import make_field
from joined_field import independent_fd
from moment_shear_slice_repair import moment_vector
from paper_moment_bridge import bump
from radial_continuation import ROOT


SLICE = 'delayed005_wide04_reoptimized_E.json'
INTERVALS_TRIAL = ((1.005, 1.06), (1.005, 1.12), (1.005, 1.24))


def fit_geometry(field, tau):
    X, weights = grid(order=48)
    eta = .3
    U_changed, _ = profile(field.base, X, eta, tau)
    U_actual, E = profile(field, X, eta, tau)
    B = correction_modes(X, field.width, field.degree, field.start_X)
    orth, factor = np.linalg.qr((B*np.sqrt(weights)).T)
    Q = orth.T/np.sqrt(weights)
    c0 = field.u_rows[1].copy()
    U = U_changed+c0@B
    if np.max(np.abs(U-U_actual)) > 1e-7:
        raise ValueError('Physical lift and coefficient representation differ')
    H = np.sqrt(2*X)*E
    linear = np.array([Q@weights, Q@(weights*H)])
    energy = 2*Q@(weights*U)
    constraints = np.vstack((linear, energy))
    wq = energy-linear.T@np.linalg.solve(
        linear@linear.T, linear@energy)
    if np.linalg.norm(wq) < 1e-9:
        raise ValueError('No quadratic S repair direction')
    dense = np.linspace(1., 3., 2001)
    dx = dense[1]-dense[0]
    Bdense = correction_modes(dense, field.width,
                              field.degree, field.start_X)
    Bxx = np.gradient(np.gradient(Bdense, dx, axis=1), dx, axis=1)
    baseline_curvature = float(np.linalg.norm(c0@Bxx))
    directions = []
    for lo, hi in INTERVALS_TRIAL:
        shape = bump(X, lo, hi)[0]
        vq = Q@(weights*shape)
        vq -= constraints.T@np.linalg.solve(
            constraints@constraints.T, constraints@vq)
        v_raw = np.linalg.solve(factor, vq)
        v_profile = v_raw@B
        scale = float(np.max(np.abs(v_profile)))
        if scale < 1e-8:
            raise ValueError('Projected direction vanished')
        vq /= scale
        v_raw /= scale
        curvature = float(np.linalg.norm(v_raw@Bxx))
        bound = min(.75, .5*baseline_curvature/max(curvature, 1e-30))
        directions.append(dict(interval=(lo, hi), vq=vq,
                               curvature=curvature, bound=bound))
    return X, weights, Q, factor, U, E, c0, wq, directions


def repaired_coefficients(a, vq, wq, Q, factor, U, weights, c0):
    vf, wf = vq@Q, wq@Q
    quadratic = float(weights@(wf**2))
    linear = float(2*weights@(U*wf)+2*a*weights@(vf*wf))
    constant = float(2*a*weights@(U*vf)+a*a*weights@(vf**2))
    roots = np.roots([quadratic, linear, constant])
    real = [float(np.real(root)) for root in roots
            if abs(np.imag(root)) < 1e-8]
    if not real:
        return None
    b = min(real, key=abs)
    raw_delta = np.linalg.solve(factor, a*vq+b*wq)
    return c0+raw_delta, b


def nodes(field, xs, etas, tau):
    X, eta = np.meshgrid(xs, etas, indexing='ij')
    points = field.compact.joined.inner.from_similarity(
        X.ravel(), eta.ravel(), tau)
    return points, X.ravel(), eta.ravel()


def residual(field, points, tau):
    return independent_fd(field, points, tau,
                          .001*np.sqrt(field.nu*tau), .00025*tau)


def stats(R):
    return dict(max=float(np.max(np.linalg.norm(R, axis=1))),
                rms=float(np.sqrt(np.mean(np.sum(R**2, axis=1)))))


def run():
    tau = .5*2**(-5.5)
    field = CoupledMomentPhysicalLift(slice_filename=SLICE)
    Xq, weights, Q, factor, U, E, c0, wq, directions = fit_geometry(
        field, tau)
    train, train_x, train_eta = nodes(
        field, (1.01, 1.0125, 1.015, 1.0175, 1.02, 1.025), (.3,), tau)
    baseline, _ = residual(field, train, tau)
    candidates = []
    eps = .05
    for direction in directions:
        samples = []
        for sign in (-1., 1.):
            repaired = repaired_coefficients(
                sign*eps, direction['vq'], wq, Q, factor, U, weights, c0)
            if repaired is None:
                raise ValueError('Small tangent trial left S manifold')
            field.u_rows[1] = repaired[0]
            R, _ = residual(field, train, tau)
            samples.append(R)
        field.u_rows[1] = c0
        slope = (samples[1]-samples[0])/(2*eps)
        suggested = float(-np.sum(baseline*slope)
                          /max(np.sum(slope**2), 1e-30))
        a = float(np.clip(suggested, -direction['bound'],
                          direction['bound']))
        for strength in (1., .5):
            trial_a = a*strength
            repaired = repaired_coefficients(
                trial_a, direction['vq'], wq, Q, factor, U, weights, c0)
            if repaired is None:
                continue
            field.u_rows[1] = repaired[0]
            R, _ = residual(field, train, tau)
            candidates.append(dict(interval=direction['interval'],
                                   a=trial_a, b=repaired[1],
                                   coefficients=repaired[0].tolist(),
                                   train=stats(R),
                                   curvature_response=direction['curvature'],
                                   trust_bound=direction['bound']))
        field.u_rows[1] = c0
    selected = min(candidates, key=lambda item: item['train']['max'])
    field.u_rows[1] = np.asarray(selected['coefficients'])
    U_trial = U+(field.u_rows[1]-c0)@correction_modes(
        Xq, field.width, field.degree, field.start_X)
    target_moments = moment_vector(U, E, Xq, weights)
    defect = moment_vector(U_trial, E, Xq, weights)-target_moments
    if np.max(np.abs(defect)) > 1e-7:
        raise ValueError('Selected trial lost five-moment constraints')
    screens = []
    for name, xs, etas, t in (
        ('train', (1.01, 1.0125, 1.015, 1.0175, 1.02, 1.025),
         (.3,), tau),
        ('space_holdout', (1.01125, 1.01375, 1.01625,
                           1.01875, 1.0225), (.29, .31), tau),
        ('time_holdout', (1.0115, 1.0145, 1.0175,
                          1.0205, 1.0235), (.292, .308),
         .5*2**(-5.4))):
        points, xx, ee = nodes(field, xs, etas, t)
        after, div = residual(field, points, t)
        field.u_rows[1] = c0
        before, _ = residual(field, points, t)
        field.u_rows[1] = np.asarray(selected['coefficients'])
        screen = dict(name=name, node_count=len(points),
                      before=stats(before), after=stats(after),
                      max_abs_divergence=float(np.max(np.abs(div))),
                      nodes=[dict(X=float(x), eta=float(e),
                                  before=float(np.linalg.norm(rb)),
                                  after=float(np.linalg.norm(ra)))
                             for x, e, rb, ra in zip(xx, ee, before, after)])
        screens.append(screen)
        print(json.dumps({key: screen[key] for key in
                          ('name', 'node_count', 'before', 'after')}), flush=True)
    report = dict(slice_filename=SLICE, tau=tau,
                  directions=[dict(interval=d['interval'],
                                   curvature_response=d['curvature'],
                                   trust_bound=d['bound'])
                              for d in directions],
                  baseline_train=stats(baseline),
                  candidates=candidates, selected=selected,
                  selected_slice_moment_defect=defect.tolist(),
                  screens=screens,
                  scope='Three projected low-dimensional U directions on '
                        'the fixed-eta=.3 M/J/S manifold, screened with '
                        'full Cartesian residual on disjoint nodes and '
                        'nearby time. E fixed, eta=.2 row fixed; no '
                        'continuous cone or global admission.',
                  accepted=False)
    (ROOT/'delayed_momentum_tangent_screen.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
