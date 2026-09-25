"""Pressure-only diagnostic for the curvature-optimized physical lift.

The compact pressure modes preserve the velocity and its five slice moments.
They cannot alter the azimuthal component of an axisymmetric residual.
"""
import json

import numpy as np

from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from joined_field import coordinates, independent_fd
from paper_moment_bridge import bump
from radial_continuation import ROOT
from radial_moment_step import septic_step


RADIAL_INTERVALS = ((1., 1.08), (1., 1.2), (1., 1.5),
                    (1.03, 1.3), (1.1, 2.), (1.4, 2.6),
                    (1.8, 3.), (2.5, 3.), (2.9, 3.))


def pressure_basis(points, tau, field):
    """Return compact pressures and their Cartesian physical gradients."""
    points = np.asarray(points, float)
    sn = np.sqrt(field.nu)
    radius = np.hypot(points[:, 0], points[:, 1])
    co = coordinates(radius/sn, points[:, 2]/sn, tau, field.heat.h)
    X, eta, q = (np.asarray(co[key]) for key in ('X', 'eta', 'q'))
    rise, rise_d = septic_step((eta-.1)/.1)
    fall, fall_d = septic_step((eta-.3)/.1)
    cutoff = rise*(1-fall)
    cutoff_d = rise_d*(1-fall)/.1-rise*fall_d/.1
    axial = ((cutoff, cutoff_d),
             (cutoff*(eta-.25)/.05,
              cutoff_d*(eta-.25)/.05+cutoff/.05))
    scale = field.nu*q**(-2*field.A)
    xr = np.asarray(co['X_r'])/sn
    xz = np.asarray(co['X_z'])/sn
    ez = np.asarray(co['eta_z'])/sn
    qz = np.asarray(co['q_z'])/sn
    cos = points[:, 0]/np.maximum(radius, 1e-300)
    sin = points[:, 1]/np.maximum(radius, 1e-300)
    values, gradients = [], []
    for lo, hi in RADIAL_INTERVALS:
        radial, radial_d = bump(X, lo, hi)
        for axial_value, axial_d in axial:
            shape = radial*axial_value
            shape_x = radial_d*axial_value
            shape_eta = radial*axial_d
            gr = scale*shape_x*xr
            gz = scale*(shape_x*xz+shape_eta*ez
                        -2*field.A*qz/q*shape)
            values.append(scale*shape)
            gradients.append(np.column_stack((gr*cos, gr*sin, gz)))
    return np.stack(values, axis=1), np.stack(gradients, axis=-1)


class CurvaturePressureField:
    def __init__(self, base, amplitudes):
        self.base = base
        self.nu = base.nu
        self.amplitudes = np.asarray(amplitudes, float)

    def fields(self, points, tau):
        velocity, pressure = self.base.fields(points, tau)
        basis, _ = pressure_basis(points, tau, self.base)
        return velocity, pressure+basis@self.amplitudes


def sample(field, xs, etas, tau):
    X, eta = np.meshgrid(xs, etas, indexing='ij')
    points = field.compact.joined.inner.from_similarity(
        X.ravel(), eta.ravel(), tau)
    return points, X.ravel(), eta.ravel()


def components(residual, points):
    radius = np.hypot(points[:, 0], points[:, 1])
    radial = (residual[:, 0]*points[:, 0]
              +residual[:, 1]*points[:, 1])/radius
    angular = (-residual[:, 0]*points[:, 1]
               +residual[:, 1]*points[:, 0])/radius
    return radial, angular, residual[:, 2]


def run():
    base = CoupledMomentPhysicalLift(
        slice_filename='wide_taper_curvature_optimize.json')
    tau = .5*2**(-5.5)
    train_points, _, _ = sample(
        base, (1.01, 1.025, 1.05, 1.1, 1.25, 1.5,
               2., 2.5, 2.975, 2.99), (.2, .25, .3), tau)
    step = .001*np.sqrt(base.nu*tau)
    baseline, _ = independent_fd(base, train_points, tau,
                                 step, .00025*tau)
    _, grad = pressure_basis(train_points, tau, base)
    matrix = grad.reshape(-1, grad.shape[-1])
    norms = np.maximum(np.linalg.norm(matrix, axis=0), 1e-300)
    ridge = .01
    augmented = np.vstack((matrix/norms,
                           ridge*np.eye(matrix.shape[1])))
    rhs = np.r_[-baseline.ravel(), np.zeros(matrix.shape[1])]
    coefficients = np.linalg.lstsq(augmented, rhs, rcond=None)[0]/norms
    candidate = CurvaturePressureField(base, coefficients)
    rows = []
    for label, xs, etas, t in (
        ('train', (), (), tau),
        ('space_holdout', (1.015, 1.075, 1.4, 2.75, 2.985),
         (.225, .275), tau),
        ('time_holdout', (1.015, 1.075, 1.4, 2.75, 2.985),
         (.225, .275), .5*2**(-5.25))):
        if label == 'train':
            points, before = train_points, baseline
        else:
            points, _, _ = sample(base, xs, etas, t)
            before, _ = independent_fd(
                base, points, t, .001*np.sqrt(base.nu*t), .00025*t)
        _, gradients = pressure_basis(points, t, base)
        after = before+np.einsum('nik,k->ni', gradients, coefficients)
        angular = components(before, points)[1]
        row = dict(name=label, node_count=len(points),
                   before_max=float(np.max(np.linalg.norm(before, axis=1))),
                   after_max=float(np.max(np.linalg.norm(after, axis=1))),
                   angular_floor_max=float(np.max(np.abs(angular))))
        rows.append(row)
        print(json.dumps(row), flush=True)
    # A direct finite-difference check keeps the analytic pressure derivative honest.
    checked, _ = independent_fd(candidate, train_points[:1], tau,
                                step, .00025*tau)
    predicted = baseline[:1]+np.einsum(
        'nik,k->ni', grad[:1], coefficients)
    report = dict(slice_filename='wide_taper_curvature_optimize.json',
                  radial_intervals=RADIAL_INTERVALS,
                  axial_support=[.1, .4], ridge=ridge,
                  pressure_amplitudes=coefficients.tolist(), rows=rows,
                  gradient_fd_difference=float(np.linalg.norm(
                      checked-predicted)),
                  scope='Pressure-only compact basis; sampled nodes only. '
                        'No global residual or volume L2 admission.',
                  accepted=False)
    (ROOT/'curvature_pressure_screen.json').write_text(
        json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    run()
