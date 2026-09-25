"""Project the balanced bridge's complete residual onto six pressure modes.

This diagnoses the pressure-removable portion of the physical momentum
defect. The pressure increment is localized by the existing axial cutoff
and vanishes with four radial endpoint jets.
"""
import json

import numpy as np

from bridge_swirl_moment_balance import make_field
from joined_field import coordinates, independent_fd
from radial_continuation import ROOT
from wide_pressure_fit import basis_and_gradient


class AdditionalBridgePressure:
    def __init__(self, base, amplitudes):
        self.base = base
        self.compact = base.compact
        self.heat = base.heat
        self.nu = base.nu
        self.amplitudes = np.asarray(amplitudes, float)

    def attachment_radius(self, tau):
        return self.base.attachment_radius(tau)

    def fields(self, points, tau):
        pts = np.asarray(points, float)
        ts = np.broadcast_to(np.asarray(tau, float), (len(pts),))
        velocity, pressure = self.base.fields(pts, ts)
        co = coordinates(np.hypot(pts[:, 0], pts[:, 1])/np.sqrt(self.nu),
                         pts[:, 2]/np.sqrt(self.nu), ts,
                         self.compact.joined.inner.h)
        chi, _ = self.compact.cutoff(co['eta'])
        active = chi > 0
        if np.any(active):
            values, _ = basis_and_gradient(
                pts[active], ts[active], self.compact.joined)
            pressure[active] += chi[active]*(values@self.amplitudes)
        return velocity, pressure


def points_at(field, k, xs, etas):
    tau = .5*2**(-k)
    X, eta = np.meshgrid(xs, etas, indexing='ij')
    return field.compact.joined.inner.from_similarity(
        X.ravel(), eta.ravel(), tau), tau


def residual(field, points, tau):
    return independent_fd(field, points, tau,
                          .001*np.sqrt(field.nu*tau), .00025*tau)


def fit_pressure(base, points, tau, ridge=.003):
    before, _ = residual(base, points, tau)
    _, gradients = basis_and_gradient(points, tau, base.compact.joined)
    matrix = gradients.reshape(-1, gradients.shape[-1])
    scales = np.maximum(np.linalg.norm(matrix, axis=0), 1e-30)
    system = np.vstack((matrix/scales,
                        ridge*np.eye(matrix.shape[1])))
    rhs = np.r_[-before.ravel(), np.zeros(matrix.shape[1])]
    amplitudes = np.linalg.lstsq(system, rhs, rcond=None)[0]/scales
    return amplitudes, before, gradients


def run():
    base = make_field(-.27311627313724157)
    train, tau = points_at(base, 5.5, [.5859375, .75, 1., 1.25],
                           [.2, .3])
    amplitudes, train_before, gradients = fit_pressure(base, train, tau)
    corrected = AdditionalBridgePressure(base, amplitudes)
    rows = []
    for name, k, xs, etas in (
            ('train', 5.5, [.5859375, .75, 1., 1.25], [.2, .3]),
            ('space_holdout', 5.5, [.68, .9, 1.1, 1.34], [.25, .35]),
            ('time_holdout', 5.25, [.68, .9, 1.1, 1.34], [.25, .35])):
        points, time = points_at(base, k, xs, etas)
        before = train_before if name == 'train' else residual(base, points, time)[0]
        after, divergence = residual(corrected, points, time)
        _, G = basis_and_gradient(points, time, base.compact.joined)
        predicted = before+np.einsum('nik,k->ni', G, amplitudes)
        row = dict(name=name, k=k, X=xs, eta=etas,
                   before_max=float(np.max(np.linalg.norm(before, axis=1))),
                   after_max=float(np.max(np.linalg.norm(after, axis=1))),
                   after_angular_max=float(np.max(np.abs(after[:, 1]))),
                   after_divergence_max=float(np.max(np.abs(divergence))),
                   analytic_fd_gradient_discrepancy=float(
                       np.max(np.linalg.norm(after-predicted, axis=1))))
        rows.append(row)
        print(json.dumps(row), flush=True)
    report = dict(bridge_coefficients=dict(symmetric=-.27311627313724157,
                                           outer=1., minimum=-2.75,
                                           moment=4.5),
                  pressure_amplitudes=amplitudes.tolist(),
                  ridge=.003, rows=rows,
                  scope='Six pressure modes fitted to physical complete momentum on eight points, with disjoint space/time holdouts. This is not a five-moment match or global PDE acceptance.',
                  accepted=False)
    (ROOT/'balanced_bridge_pressure_screen.json').write_text(
        json.dumps(report, indent=2)+'\n')
    print(json.dumps(dict(pressure_amplitudes=report['pressure_amplitudes'],
                          rows=rows)))


if __name__ == '__main__':
    run()
