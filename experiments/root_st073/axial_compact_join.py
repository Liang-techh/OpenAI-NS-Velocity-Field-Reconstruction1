"""Finite-energy axial localization of the pressure-corrected ST073 join.

The meridional velocity is localized through its streamfunction, not by
componentwise multiplication, so spatial incompressibility is retained.
The resulting full-space field is a PDE diagnostic candidate, not accepted.
"""
import json

import numpy as np
from scipy.special import expit

from joined_field import coordinates, independent_fd
from poloidal_bridge import coefficients
from radial_continuation import ROOT
from wide_pressure_fit import load_pressure_candidate


class AxiallyCompactField:
    def __init__(self, base, eta_flat=.2, eta_outer=.49):
        if not 0 < eta_flat < eta_outer < base.base.inner.p.eta_max:
            raise ValueError('Axial cutoff must lie inside the registered slab')
        self.base = base
        self.joined = base.base
        self.nu = self.joined.nu
        self.eta_flat = float(eta_flat)
        self.eta_outer = float(eta_outer)

    def cutoff(self, eta):
        eta = np.asarray(eta, float)
        abs_eta = np.abs(eta)
        chi = np.zeros_like(eta)
        chi_e = np.zeros_like(eta)
        chi[abs_eta <= self.eta_flat] = 1.
        middle = (abs_eta > self.eta_flat) & (abs_eta < self.eta_outer)
        s = (abs_eta[middle]-self.eta_flat)/(self.eta_outer-self.eta_flat)
        logistic = expit(1/s-1/(1-s))
        chi[middle] = logistic
        chi_e[middle] = (logistic*(1-logistic)
                         *(-1/s**2-1/(1-s)**2)
                         *np.sign(eta[middle])
                         /(self.eta_outer-self.eta_flat))
        return chi, chi_e

    def streamfunction(self, points, tau, co=None):
        pts = np.asarray(points, float)
        ts = np.broadcast_to(np.asarray(tau, float), (len(pts),))
        sn = np.sqrt(self.nu)
        r = np.hypot(pts[:, 0], pts[:, 1])
        if co is None:
            co = coordinates(r/sn, pts[:, 2]/sn, ts, self.joined.inner.h)
        psi = np.zeros(len(pts))
        for i, (radius, t, X, eta, q) in enumerate(
                zip(r, ts, co['X'], co['eta'], co['q'])):
            ri, ro, _, bridge = coefficients(
                self.joined.inner, float(eta), float(t),
                self.joined.join_X, self.joined.outer_ratio)
            if radius <= ri:
                c = self.joined.inner.coefficients(float(eta), float(q))[2, :, 0]
                primitive = np.r_[0, c/np.arange(1, len(c)+1)]
                psi[i] = self.nu**1.5*q*np.polynomial.polynomial.polyval(
                    X, primitive)
            elif radius < ro:
                y = (radius-ri)/(ro-ri)
                psi[i] = np.polynomial.polynomial.polyval(y, bridge)
        return psi

    def fields(self, points, tau):
        pts = np.asarray(points, float)
        if pts.ndim != 2 or pts.shape[1] != 3 or not np.isfinite(pts).all():
            raise ValueError('Finite (n,3) Cartesian points required')
        ts = np.broadcast_to(np.asarray(tau, float), (len(pts),))
        if np.any(ts < .5/64) or np.any(ts > .5):
            raise ValueError('Outside the registered time interval')
        velocity = np.zeros_like(pts)
        pressure = np.zeros(len(pts))
        sn = np.sqrt(self.nu)
        r = np.hypot(pts[:, 0], pts[:, 1])
        co = coordinates(r/sn, pts[:, 2]/sn, ts, self.joined.inner.h)
        chi, chi_e = self.cutoff(co['eta'])
        active = chi > 0
        if not np.any(active):
            return velocity, pressure
        chosen = pts[active]
        sub_co = {key: np.asarray(value)[active] for key, value in co.items()}
        base_u, base_p = self.base.fields(chosen, ts[active])
        psi = self.streamfunction(chosen, ts[active], sub_co)
        radius = r[active]
        chi_z = chi_e[active]*sub_co['eta_z']/sn
        correction = -chi_z*psi/np.maximum(radius, 1e-300)
        radial = (np.column_stack((chosen[:, 0], chosen[:, 1])) /
                  np.maximum(radius[:, None], 1e-300))
        velocity[active] = chi[active, None]*base_u
        velocity[active, :2] += correction[:, None]*radial
        pressure[active] = chi[active]*base_p
        return velocity, pressure

    def exterior_energy_tail_bound(self, tau, radius=1.):
        """Upper bound on 1/2 integral |u|^2 for r>=radius at fixed tau."""
        qmax = tau/(1-self.eta_outer**2)
        max_join_radius = self.joined.outer_ratio*np.sqrt(
            2*self.nu*qmax*self.joined.join_X)
        if radius < max_join_radius:
            raise ValueError('Tail radius must be beyond the radial join')
        h = self.joined.inner.h
        zcap = np.sqrt(self.nu)*qmax**(.5-h)*self.eta_outer
        return (np.pi*2*zcap*self.nu*self.joined.c**2
                *(2*self.nu)**(1+2*h)*radius**(-4*h)/(4*h))


def load_compact_candidate():
    return AxiallyCompactField(load_pressure_candidate())


def run():
    field = load_compact_candidate()
    rows = []
    for k in (.4, 5.5):
        tau = .5*2**(-k)
        X = np.array([.03, .5859375, 2.])
        eta = np.full(3, .375)
        pts = field.joined.inner.from_similarity(X, eta, tau)
        residual, div = independent_fd(
            field, pts, tau, .001*np.sqrt(field.nu*tau), .00025*tau)
        rows.append(dict(k=k, tau=tau, X=X.tolist(), eta=eta.tolist(),
                         complete_residual_norms=np.linalg.norm(
                             residual, axis=1).tolist(),
                         finite_difference_divergence=div.tolist(),
                         exterior_energy_tail_bound_r_ge_1=
                         field.exterior_energy_tail_bound(tau)))
        print(json.dumps(rows[-1]), flush=True)
    tau=.5*2**(-5.5)
    probe = field.joined.inner.from_similarity([.03], [.375], tau)
    divergence_refinement = []
    for factor in (.001, .0005, .00025):
        _, divergence = independent_fd(
            field, probe, tau, factor*np.sqrt(field.nu*tau), .00025*tau)
        divergence_refinement.append(dict(space_step_factor=factor,
                                          divergence=float(divergence[0])))
    tau=.5*2**(-5.5)
    axis = field.joined.inner.from_similarity([0.], [0.], tau)
    outside = field.joined.inner.from_similarity([.03], [.6], tau)
    axis_u, axis_p = field.fields(axis, tau)
    outside_u, outside_p = field.fields(outside, tau)
    report = dict(rows=rows, divergence_refinement=divergence_refinement,
                  tau_window=[.5/64, .5], viscosity=field.nu,
                  forcing='zero for the reported residual diagnostic',
                  global_velocity_callable=True,
                  finite_energy_each_registered_time=True,
                  solenoidal_by_construction=True,
                  eta_flat=field.eta_flat,
                  eta_outer=field.eta_outer,
                  axis_velocity=axis_u[0].tolist(),
                  axis_pressure=float(axis_p[0]),
                  outside_velocity=outside_u[0].tolist(),
                  outside_pressure=float(outside_p[0]),
                  scope='Globally callable, axisymmetric, solenoidal and '
                        'finite-energy at each registered positive time. '
                        'The heat radial tail decays slowly but is integrable '
                        'for h>0. Sampled unforced residual in the axial '
                        'collar is diagnostic only; no global PDE acceptance, '
                        'critical-time extension, five-moment restoration or '
                        'nonaxisymmetric stress correction.',
                  pde_validated=False, accepted=False)
    (ROOT/'axial_compact_join.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
