"""Paper-inspired radial heat swirl with a 1D axial heat envelope.

Exterior-only: the source radial profile is singular on the axis. This
module supplies an exact angular heat target for a later global join.
"""
import json

import numpy as np

from joined_field import independent_fd
from radial_continuation import ROOT
from heat_exterior import physical
from wide_pressure_fit import load_pressure_candidate


class AxiallyHeatedExterior:
    def __init__(self, nu=.01, h=.005, amplitude=None, heat_age=2048.):
        if heat_age <= 0:
            raise ValueError('Heat age must be positive')
        self.nu = float(nu)
        self.h = float(h)
        self.heat_age = float(heat_age)
        self.amplitude = (float(amplitude) if amplitude is not None else
                          load_pressure_candidate().base.c)

    def envelope(self, z, tau):
        age = self.heat_age+.5-np.asarray(tau, float)
        z = np.asarray(z, float)
        G = np.sqrt(self.heat_age/age)*np.exp(
            -z*z/(4*self.nu*age))
        Gz = -z/(2*self.nu*age)*G
        return G, Gz

    def fields(self, points, tau):
        pts = np.asarray(points, float)
        if pts.ndim != 2 or pts.shape[1] != 3:
            raise ValueError('Cartesian (n,3) points required')
        if np.any(np.hypot(pts[:, 0], pts[:, 1]) <= 0):
            raise ValueError('Exterior excludes the singular axis')
        radial = physical(pts, tau, c=self.amplitude,
                          h=self.h, nu=self.nu)
        G, _ = self.envelope(pts[:, 2], tau)
        return G[:, None]*radial['velocity'], G*G*radial['pressure']

    def analytic_residual(self, points, tau):
        pts = np.asarray(points, float)
        radial = physical(pts, tau, c=self.amplitude,
                          h=self.h, nu=self.nu)
        G, Gz = self.envelope(pts[:, 2], tau)
        residual = G[:, None]*radial['residual']
        residual[:, 2] = 2*G*Gz*radial['pressure']
        return residual

    def axial_energy_factor(self, tau):
        """Integral of G(z,tau)^2 dz over the whole axial line."""
        age = self.heat_age+.5-tau
        return self.heat_age*np.sqrt(2*np.pi*self.nu/age)


def run():
    field = AxiallyHeatedExterior()
    inner = load_pressure_candidate().base.inner
    rows = []
    for k in (3., 5.5):
        tau = .5*2**(-k)
        etas = np.array([.2725, .345, .4175])
        points = inner.from_similarity(np.full(3, 2.), etas, tau)
        direct, divergence = independent_fd(
            field, points, tau, .001*np.sqrt(field.nu*tau), .00025*tau)
        analytic = field.analytic_residual(points, tau)
        row = dict(k=k, tau=tau, X=2., eta=etas.tolist(),
                   direct_angular_residual=direct[:, 1].tolist(),
                   direct_axial_residual=direct[:, 2].tolist(),
                   analytic_axial_residual=analytic[:, 2].tolist(),
                   direct_max=float(np.max(np.linalg.norm(direct, axis=1))),
                   analytic_direct_vector_discrepancy=float(np.max(
                       np.linalg.norm(direct-analytic, axis=1))),
                   finite_difference_divergence_max=float(np.max(
                       np.abs(divergence))),
                   axial_energy_factor=field.axial_energy_factor(tau))
        rows.append(row)
        print(json.dumps(row), flush=True)
    report = dict(rows=rows, heat_age=field.heat_age,
                  radial_heat_amplitude=field.amplitude,
                  scope='Exterior r>0 target only. Radial OpenAI-inspired '
                        'heat swirl multiplied by an axial Gaussian solving '
                        'the one-dimensional heat equation, so angular '
                        'viscous/time residual cancels analytically. '
                        'Axial pressure gradient remains; singular axis '
                        'requires inner matching before a global field.',
                  accepted=False, pde_validated=False)
    (ROOT/'axially_heated_exterior.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
