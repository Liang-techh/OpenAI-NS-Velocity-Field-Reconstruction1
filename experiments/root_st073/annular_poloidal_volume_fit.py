"""Fit two solenoidal radial/axial collar bubbles to physical-volume momentum."""
import json

import numpy as np
from numpy.polynomial.legendre import Legendre
from scipy.optimize import minimize

from joined_field import ROOT
from joint_collar_fit import kinematics, nodes
from radial_peak_cone import current_field


class AnnularPoloidalMode:
    """Odd streamfunction in z gives even radial velocity at axial hotspots."""

    def __init__(self, base, lower, upper, reference_tau=.5/64,
                 temporal_power=2., axial_degree=0):
        self.base = base
        self.nu = base.nu
        self.lower = lower
        self.upper = upper
        self.reference_tau = reference_tau
        self.temporal_power = temporal_power
        self.axial_degree = axial_degree
        self.axial_poly = Legendre.basis(axial_degree)
        self.axial_poly_derivative = self.axial_poly.deriv()

    def fields(self, points, tau):
        pts = np.asarray(points, float)
        times = np.broadcast_to(tau, (len(pts),))
        velocity = np.zeros_like(pts)
        for i, (point, t) in enumerate(zip(pts, times)):
            x, ycart, z = point
            r = np.hypot(x, ycart)
            _, radius, zflat, zsupp = self.base.support(t)
            lo, hi = self.lower*radius, self.upper*radius
            if r <= lo or r >= hi or abs(z) <= zflat or abs(z) >= zsupp:
                continue
            rho = (r-lo)/(hi-lo)
            s = (abs(z)-zflat)/(zsupp-zflat)
            radial = 1024*rho**5*(1-rho)**5
            radial_r = (5120*rho**4*(1-rho)**4*(1-2*rho)
                        /(hi-lo))
            axial_bubble = 1024*s**5*(1-s)**5
            axial_bubble_s = 5120*s**4*(1-s)**4*(1-2*s)
            axial = axial_bubble*self.axial_poly(2*s-1)
            axial_z = ((axial_bubble_s*self.axial_poly(2*s-1)
                        +2*axial_bubble*self.axial_poly_derivative(2*s-1))
                       /(zsupp-zflat))
            ur = -radius**2*radial*axial_z/r
            uz = radius**2*radial_r*axial*np.sign(z)/r
            velocity[i] = [ur*x/r, ur*ycart/r, uz]
        velocity *= ((self.reference_tau/times)**self.temporal_power)[:, None]
        return velocity, np.zeros(len(pts))


class AnnularPoloidalCandidate:
    windows = ((.12, .36), (.28, .55))

    def __init__(self, base, amplitudes, reference_tau=.5/64):
        self.base = base
        self.nu = base.nu
        self.amplitudes = np.asarray(amplitudes)
        self.modes = [AnnularPoloidalMode(base.base, lo, hi, reference_tau)
                      for lo, hi in self.windows]

    def support(self, tau):
        return self.base.support(tau)

    def fields(self, points, tau):
        u, p = self.base.fields(points, tau)
        for amplitude, mode in zip(self.amplitudes, self.modes):
            if amplitude:
                u += amplitude*mode.fields(points, tau)[0]
        return u, p


def load_candidate():
    report = json.loads((ROOT/'compact_potential'/'annular_poloidal_volume_fit.json').read_text())
    return AnnularPoloidalCandidate(current_field(), report['amplitudes'])


def run():
    field = current_field()
    tau = .5/64
    points, weights = nodes(field, tau, 8)
    hs = .0005*np.sqrt(field.nu*tau)
    ht = .0001*tau
    u0, J0, part0 = kinematics(field, points, tau, hs, ht)
    modes = [AnnularPoloidalMode(field.base, lo, hi, tau)
             for lo, hi in AnnularPoloidalCandidate.windows]
    data = [kinematics(mode, points, tau, hs, ht) for mode in modes]
    du = np.stack([row[0] for row in data], axis=-1)
    dJ = np.stack([row[1] for row in data], axis=-1)
    dpart = np.stack([row[2] for row in data], axis=-1)
    rootw = np.sqrt(weights)[:, None]
    baseline = part0+np.einsum('nij,nj->ni', J0, u0)
    linear = (dpart+np.einsum('nijk,nj->nik', dJ, u0)
              +np.einsum('nij,njk->nik', J0, du))
    matrix = (rootw[:, :, None]*linear).reshape(-1, len(modes))
    target = -(rootw*baseline).reshape(-1)
    initial, *_ = np.linalg.lstsq(matrix, target, rcond=None)
    def recombine(a):
        u = u0+np.einsum('nik,k->ni', du, a)
        J = J0+np.einsum('nijk,k->nij', dJ, a)
        part = part0+np.einsum('nik,k->ni', dpart, a)
        return u, J, part+np.einsum('nij,nj->ni', J, u)
    scale = max(np.linalg.norm(rootw*baseline), 1.)
    def objective(a):
        u, J, R = recombine(a)
        tangent = (dpart+np.einsum('nijk,nj->nik', dJ, u)
                   +np.einsum('nij,njk->nik', J, du))
        value = .5*np.sum(weights[:, None]*R**2)/scale**2
        gradient = np.einsum('n,ni,nik->k', weights, R, tangent)/scale**2
        return value, gradient
    fit = minimize(objective, initial, jac=True, method='L-BFGS-B',
                   bounds=[(-1., 1.)]*len(modes),
                   options={'maxiter': 300, 'ftol': 1e-13})
    _, _, R = recombine(fit.x)
    report = {'tau': tau, 'order': 8,
              'radial_windows_fraction_of_support': AnnularPoloidalCandidate.windows,
              'temporal_power': 2., 'reference_tau': tau,
              'linear_amplitudes': initial.tolist(),
              'amplitudes': fit.x.tolist(),
              'optimizer_success': bool(fit.success),
              'optimizer_message': fit.message,
              'baseline': {'max': float(np.max(np.linalg.norm(baseline, axis=1))),
                           'physical_volume_l2': float(np.linalg.norm(rootw*baseline))},
              'candidate': {'max': float(np.max(np.linalg.norm(R, axis=1))),
                            'physical_volume_l2': float(np.linalg.norm(rootw*R))},
              'scope': 'Two compact solenoidal axial-collar modes fitted to full-support Gauss8 physical-volume residual at one time. No local cone, held-out time, or global supremum certificate.',
              'accepted': False}
    out = ROOT/'compact_potential'/'annular_poloidal_volume_fit.json'
    out.write_bytes((json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    run()
