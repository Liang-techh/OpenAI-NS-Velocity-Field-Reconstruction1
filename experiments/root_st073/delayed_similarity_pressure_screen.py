"""Fit compact similarity-scaled pressure to the full ST073 momentum defect.

Pressure leaves velocity, divergence, and the fixed-slice moments unchanged.
The fit is a local diagnostic; it cannot remove solenoidal momentum error.
"""

import json

import numpy as np

from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from delayed_momentum_tangent_screen import nodes, residual
from joined_field import coordinates
from paper_moment_bridge import bump
from radial_continuation import ROOT


RADIAL_INTERVALS = ((1.005, 1.05), (1.005, 1.12), (1.02, 1.2))
AXIAL_INTERVALS = ((.14, .36), (.2, .4))


class SimilarityPressurePatch:
    def __init__(self, base, coefficients):
        self.base = base
        self.nu = base.nu
        self.compact = base.compact
        self.coefficients = np.asarray(coefficients, float)
        self.beta = 2*base.A

    def basis(self, points, tau):
        pts = np.asarray(points, float)
        ts = np.broadcast_to(np.asarray(tau, float), (len(pts),))
        radius = np.hypot(pts[:, 0], pts[:, 1])
        sn = np.sqrt(self.nu)
        co = coordinates(radius/sn, pts[:, 2]/sn, ts,
                         self.base.heat.h)
        q, X, eta = (np.asarray(co[key]) for key in ('q', 'X', 'eta'))
        factor = q**(-self.beta)
        Xr = np.asarray(co['X_r'])/sn
        Xz = np.asarray(co['X_z'])/sn
        etaz = np.asarray(co['eta_z'])/sn
        qz = np.asarray(co['q_z'])/sn
        ca = np.divide(pts[:, 0], radius, out=np.ones_like(radius),
                       where=radius > 0)
        sa = np.divide(pts[:, 1], radius, out=np.zeros_like(radius),
                       where=radius > 0)
        values = []
        gradients = []
        for ra, rb in RADIAL_INTERVALS:
            bx, bxd = bump(X, ra, rb)
            for ea, eb in AXIAL_INTERVALS:
                be, bed = bump(eta, ea, eb)
                value = factor*bx*be
                pr = factor*bxd*Xr*be
                pz = factor*(bxd*Xz*be+bx*bed*etaz
                             -self.beta*bx*be*qz/q)
                values.append(value)
                gradients.append(np.column_stack((pr*ca, pr*sa, pz)))
        return np.column_stack(values), np.stack(gradients, axis=2)

    def fields(self, points, tau):
        velocity, pressure = self.base.fields(points, tau)
        values, _ = self.basis(points, tau)
        return velocity, pressure+values@self.coefficients


def stats(R):
    norms = np.linalg.norm(R, axis=1)
    return dict(max=float(np.max(norms)),
                rms=float(np.sqrt(np.mean(norms**2))),
                component_max=np.max(np.abs(R), axis=0).tolist())


def run():
    field_name = 'delayed005_rise146_degree31.json'
    base = CoupledMomentPhysicalLift(slice_filename=field_name)
    patch = SimilarityPressurePatch(base, np.zeros(6))
    train = []
    for tau in (.5*2**(-5.5), .5*2**(-5.25)):
        points, X, eta = nodes(
            base, (1.01, 1.015, 1.02, 1.03, 1.05), (.2, .3), tau)
        R, _ = residual(base, points, tau)
        _, G = patch.basis(points, tau)
        train.append((R, G))
    R_train = np.concatenate([item[0] for item in train], axis=0)
    G_train = np.concatenate([item[1] for item in train], axis=0)
    matrix = G_train.reshape(-1, 6)
    columns = np.linalg.norm(matrix, axis=0)
    if np.min(columns) <= 0:
        raise ValueError('Pressure basis has a zero column')
    normalized = matrix/columns
    ridge = .03
    fitted = np.linalg.solve(
        normalized.T@normalized+ridge**2*np.eye(6),
        -normalized.T@R_train.reshape(-1))
    coefficients = fitted/columns
    patch.coefficients = coefficients
    after_train = R_train+G_train@coefficients
    holdout = []
    tau_holdout = .5*2**(-5.4)
    points, Xh, etah = nodes(
        base, (1.0125, 1.0175, 1.025, 1.045, 1.07),
        (.22, .28, .32), tau_holdout)
    R_holdout, _ = residual(base, points, tau_holdout)
    _, G_holdout = patch.basis(points, tau_holdout)
    after_holdout = R_holdout+G_holdout@coefficients
    R_fd, _ = residual(patch, points[:3], tau_holdout)
    gradient_check = float(np.max(np.abs(R_fd-after_holdout[:3])))
    report = dict(slice_filename=field_name, beta=patch.beta,
                  radial_intervals=RADIAL_INTERVALS,
                  axial_intervals=AXIAL_INTERVALS,
                  ridge=ridge, coefficients=coefficients.tolist(),
                  train=dict(before=stats(R_train),
                             after=stats(after_train)),
                  holdout=dict(tau=tau_holdout, X=Xh.tolist(),
                               eta=etah.tolist(),
                               before=stats(R_holdout),
                               after=stats(after_holdout),
                               node_before=np.linalg.norm(R_holdout,
                                                          axis=1).tolist(),
                               node_after=np.linalg.norm(after_holdout,
                                                         axis=1).tolist()),
                  analytic_vs_independent_fd_max=gradient_check,
                  scope='Six compact similarity-scaled pressure modes '
                        'fitted to full Cartesian momentum at two times; '
                        'disjoint space/time holdout. Velocity and five '
                        'moments unchanged. No global pressure matching, '
                        'spatial L2, or PDE acceptance.', accepted=False)
    (ROOT/'delayed_similarity_pressure_screen.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps({key: report[key] for key in
                      ('train', 'holdout',
                       'analytic_vs_independent_fd_max')}), flush=True)


if __name__ == '__main__':
    run()
