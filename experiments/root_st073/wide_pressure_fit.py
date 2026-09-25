"""Fit a smooth pressure increment on the wider ST073-V transition.

This is a pressure-only first stage. It cannot remove pressure-independent
curl or angular momentum defects; those require a velocity correction.
"""
import json

import numpy as np

from joined_field import JoinedField, independent_fd, coordinates
from radial_continuation import ROOT, FullRadialField


def basis_and_gradient(points, tau, field):
    points = np.asarray(points, float)
    sn = np.sqrt(field.nu)
    r = np.hypot(points[:, 0], points[:, 1])
    co = coordinates(r/sn, points[:, 2]/sn, tau, field.inner.h)
    X, eta, q = (np.asarray(co[name]) for name in ('X', 'eta', 'q'))
    ratio = field.outer_ratio
    y = (np.sqrt(X/field.join_X)-1)/(ratio-1)
    active = (y > 0) & (y < 1)
    yy = np.clip(y, 0, 1)
    B = 256*yy**4*(1-yy)**4
    By = 1024*yy**3*(1-yy)**3*(1-2*yy)
    B *= active
    By *= active
    safe_X = np.maximum(X, 1e-300)
    safe_r = np.maximum(r, 1e-300)
    yr = np.asarray(co['X_r'])/(sn*2*(ratio-1)*np.sqrt(safe_X*field.join_X))
    yz = np.asarray(co['X_z'])/(sn*2*(ratio-1)*np.sqrt(safe_X*field.join_X))
    ez = np.asarray(co['eta_z'])/sn
    qz = np.asarray(co['q_z'])/sn
    scale = field.nu*q**(-1-2*field.inner.h)
    radial = [(np.ones_like(y), np.zeros_like(y)), (2*y-1, 2*np.ones_like(y))]
    axial = [(np.ones_like(eta), np.zeros_like(eta)),
             (eta/.4, np.ones_like(eta)/.4),
             ((eta/.4)**2, 2*eta/.4**2)]
    values, gradients = [], []
    for R, Ry in radial:
        for A, Ae in axial:
            shape = B*R*A
            sy = (By*R+B*Ry)*A
            se = B*R*Ae
            val = scale*shape
            gr = scale*sy*yr
            gz = scale*(sy*yz+se*ez-(1+2*field.inner.h)*qz/q*shape)
            grad = np.column_stack((gr*points[:, 0]/safe_r,
                                    gr*points[:, 1]/safe_r, gz))
            values.append(val)
            gradients.append(grad)
    return np.stack(values, axis=1), np.stack(gradients, axis=-1)


class PressureBubbleField:
    def __init__(self, base, amplitudes):
        self.base = base
        self.nu = base.nu
        self.amplitudes = np.asarray(amplitudes, float)

    def fields(self, points, tau):
        velocity, pressure = self.base.fields(points, tau)
        basis, _ = basis_and_gradient(points, tau, self.base)
        return velocity, pressure + basis@self.amplitudes


def load_pressure_candidate():
    report = json.loads((ROOT/'wide_pressure_fit.json').read_text())
    inner = FullRadialField.load(
        ROOT/'NS_ST073_Full_Local_Recurrence'/'data'/'ST073-V-wide14.json')
    base = JoinedField(inner=inner, join_X=report['join_X'],
                       outer_ratio=report['outer_ratio'])
    return PressureBubbleField(base, report['pressure_amplitudes'])


def points_at(field, k, ys, etas):
    tau = .5*2**(-k)
    y, eta = np.meshgrid(ys, etas, indexing='ij')
    X = field.join_X*(1+(field.outer_ratio-1)*y.ravel())**2
    return field.inner.from_similarity(X, eta.ravel(), tau), tau


def residual(field, points, tau):
    hs = .001*np.sqrt(field.base.nu*tau) if isinstance(field, PressureBubbleField) else .001*np.sqrt(field.nu*tau)
    return independent_fd(field, points, tau, hs, .00025*tau)


def run():
    inner = FullRadialField.load(
        ROOT/'NS_ST073_Full_Local_Recurrence'/'data'/'ST073-V-wide14.json')
    base = JoinedField(inner=inner, join_X=3/32, outer_ratio=4.)
    train_points, train_tau = points_at(
        base, 5.5, [.2, .5, .8], [-.3, -.15, 0, .15, .3])
    R0, div0 = residual(base, train_points, train_tau)
    _, G = basis_and_gradient(train_points, train_tau, base)
    matrix = G.reshape(-1, G.shape[-1])
    column_scale = np.maximum(np.linalg.norm(matrix, axis=0), 1e-30)
    ridge = .001
    system = np.vstack((matrix/column_scale, ridge*np.eye(matrix.shape[1])))
    rhs = np.r_[-R0.ravel(), np.zeros(matrix.shape[1])]
    scaled_coefficients = np.linalg.lstsq(system, rhs, rcond=None)[0]
    amplitudes = scaled_coefficients/column_scale
    candidate = PressureBubbleField(base, amplitudes)
    rows = []
    for name, k, ys, etas in (
            ('train', 5.5, [.2, .5, .8], [-.3, -.15, 0, .15, .3]),
            ('space_holdout', 5.5, [.3, .7], [-.25, .1, .25]),
            ('time_holdout', 5.25, [.3, .7], [-.25, .1, .25])):
        points, tau = points_at(base, k, ys, etas)
        before = R0 if name == 'train' else residual(base, points, tau)[0]
        after, divergence = residual(candidate, points, tau)
        predicted = before + np.einsum(
            'nik,k->ni', basis_and_gradient(points, tau, base)[1], amplitudes)
        row = dict(name=name, k=k, point_count=len(points),
                   baseline_max=float(np.max(np.linalg.norm(before, axis=1))),
                   corrected_max=float(np.max(np.linalg.norm(after, axis=1))),
                   baseline_rms=float(np.sqrt(np.mean(before**2))),
                   corrected_rms=float(np.sqrt(np.mean(after**2))),
                   corrected_angular_max=float(np.max(np.abs(after[:, 1]))),
                   divergence_max=float(np.max(np.abs(divergence))),
                   analytic_gradient_fd_discrepancy=float(
                       np.max(np.linalg.norm(after-predicted, axis=1))))
        rows.append(row)
        print(json.dumps(row), flush=True)
    report = dict(join_X=base.join_X, outer_ratio=base.outer_ratio,
                  pressure_mode_count=len(amplitudes),
                  pressure_amplitudes=amplitudes.tolist(),
                  ridge=ridge, rows=rows,
                  scope='Six smooth pressure modes supported in the radial '
                        'bridge, fitted at one late time and checked at disjoint '
                        'space/time points. No velocity correction, axial '
                        'closure, five-moment match or full-domain acceptance.',
                  accepted=False)
    (ROOT/'wide_pressure_fit.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
