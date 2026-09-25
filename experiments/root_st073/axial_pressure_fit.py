"""Compact axial-pressure correction for the finite-energy ST073 field."""
import json

import numpy as np

from axial_compact_join import load_compact_candidate
from joined_field import coordinates, independent_fd
from radial_continuation import ROOT


def pressure_basis(points, tau, field):
    pts = np.asarray(points, float)
    sn = np.sqrt(field.nu)
    r = np.hypot(pts[:, 0], pts[:, 1])
    co = coordinates(r/sn, pts[:, 2]/sn, tau, field.joined.inner.h)
    X, eta, q = (np.asarray(co[name]) for name in ('X', 'eta', 'q'))
    s = (np.abs(eta)-field.eta_flat)/(field.eta_outer-field.eta_flat)
    active = (s > 0) & (s < 1)
    ss = np.clip(s, 0, 1)
    B = 256*ss**4*(1-ss)**4 * active
    Bs = 1024*ss**3*(1-ss)**3*(1-2*ss) * active
    Xr = np.asarray(co['X_r'])/sn
    Xz = np.asarray(co['X_z'])/sn
    ez = np.asarray(co['eta_z'])/sn
    qz = np.asarray(co['q_z'])/sn
    sz = np.sign(eta)*ez/(field.eta_outer-field.eta_flat)
    A = .5+field.joined.inner.h
    scale = field.nu*q**(-2*A)
    safe_r = np.maximum(r, 1e-300)
    values, gradients = [], []
    for radial_scale in (.1, 1., 10.):
        denominator = 1+X/radial_scale
        R = 1/denominator
        Rx = -1/(radial_scale*denominator**2)
        for axial_degree in (0, 1):
            T = np.ones_like(s) if axial_degree == 0 else 2*s-1
            Ts = np.zeros_like(s) if axial_degree == 0 else 2*np.ones_like(s)
            shape = B*T
            shapes = Bs*T+B*Ts
            val = scale*R*shape
            gr = scale*Rx*Xr*shape
            gz = scale*((Rx*Xz-2*A*qz/q*R)*shape+R*shapes*sz)
            grad = np.column_stack((gr*pts[:, 0]/safe_r,
                                    gr*pts[:, 1]/safe_r, gz))
            values.append(val)
            gradients.append(grad)
    return np.stack(values, axis=1), np.stack(gradients, axis=-1)


class AxialPressureField:
    def __init__(self, base, amplitudes):
        self.base = base
        self.nu = base.nu
        self.amplitudes = np.asarray(amplitudes, float)

    def fields(self, points, tau):
        velocity, pressure = self.base.fields(points, tau)
        basis, _ = pressure_basis(points, tau, self.base)
        return velocity, pressure+basis@self.amplitudes


def load_axial_pressure_candidate():
    report = json.loads((ROOT/'axial_pressure_fit.json').read_text())
    return AxialPressureField(load_compact_candidate(),
                              report['amplitudes'])


def sample_points(field, k, X_values, fractions):
    tau = .5*2**(-k)
    X, s = np.meshgrid(X_values, fractions, indexing='ij')
    eta = field.eta_flat+(field.eta_outer-field.eta_flat)*s
    points = field.joined.inner.from_similarity(X.ravel(), eta.ravel(), tau)
    return points, tau


def full_residual(field, points, tau):
    return independent_fd(field, points, tau,
                          .001*np.sqrt(field.nu*tau), .00025*tau)


def run():
    base = load_compact_candidate()
    training, train_tau = sample_points(
        base, 5.5, [.03, .5859375, 2.], [.25, .5, .75])
    baseline, _ = full_residual(base, training, train_tau)
    _, G = pressure_basis(training, train_tau, base)
    matrix = G.reshape(-1, G.shape[-1])
    scales = np.maximum(np.linalg.norm(matrix, axis=0), 1e-30)
    ridge = .001
    system = np.vstack((matrix/scales, ridge*np.eye(matrix.shape[1])))
    rhs = np.r_[-baseline.ravel(), np.zeros(matrix.shape[1])]
    amplitudes = np.linalg.lstsq(system, rhs, rcond=None)[0]/scales
    field = AxialPressureField(base, amplitudes)
    rows = []
    for name, k, Xs, fractions in (
            ('train', 5.5, [.03, .5859375, 2.], [.25, .5, .75]),
            ('space_holdout', 5.5, [.06, .3, 1.], [.375, .625]),
            ('time_holdout', 5.25, [.06, .3, 1.], [.375, .625])):
        points, tau = sample_points(base, k, Xs, fractions)
        before = baseline if name == 'train' else full_residual(
            base, points, tau)[0]
        after, div = full_residual(field, points, tau)
        predicted = before+np.einsum(
            'nik,k->ni', pressure_basis(points, tau, base)[1], amplitudes)
        row = dict(name=name, k=k, point_count=len(points),
                   baseline_max=float(np.max(np.linalg.norm(before, axis=1))),
                   corrected_max=float(np.max(np.linalg.norm(after, axis=1))),
                   baseline_vector_rms=float(np.sqrt(np.mean(before**2))),
                   corrected_vector_rms=float(np.sqrt(np.mean(after**2))),
                   baseline_angular_max=float(np.max(np.abs(before[:, 1]))),
                   corrected_angular_max=float(np.max(np.abs(after[:, 1]))),
                   divergence_max=float(np.max(np.abs(div))),
                   pressure_gradient_fd_discrepancy=float(
                       np.max(np.linalg.norm(after-predicted, axis=1))))
        rows.append(row)
        print(json.dumps(row), flush=True)
    report = dict(mode_count=len(amplitudes),
                  amplitudes=amplitudes.tolist(), ridge=ridge,
                  eta_flat=base.eta_flat, eta_outer=base.eta_outer,
                  rows=rows, scope='Six compact axial-pressure modes on the '
                        'full-space finite-energy field, fitted at one late '
                        'time and checked at disjoint space/time points. '
                        'Velocity and spatial divergence unchanged. '
                        'Angular residual and dynamical return flow remain.',
                  accepted=False, pde_validated=False)
    (ROOT/'axial_pressure_fit.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
