"""Solenoidal two-stage streamfunction step for the radial M moment.

The streamfunction rises inside the old bridge, remains constant on a
reserved annulus, and returns to zero farther out. In the flat axial
region its plateau has no velocity, while its value changes cumulative
axial flux M at the old bridge exit. The outer field is unchanged beyond
the restoration radius. This is a diagnostic mean mode, not a PDE patch.
"""
import numpy as np

from high_frequency_shear_screen import make_field
from joined_field import coordinates


def septic_step(x):
    y = np.clip(np.asarray(x, float), 0., 1.)
    value = y**4*(35-84*y+70*y**2-20*y**3)
    derivative = 140*y**3*(1-y)**3
    return value, derivative


class RadialMomentStep:
    def __init__(self, base=None, amplitude=0., restore_start=1.75,
                 restore_end=3.):
        self.base = base if base is not None else make_field(16, 0.)
        self.compact = self.base.compact
        self.heat = self.base.heat
        self.nu = self.base.nu
        self.amplitude = float(amplitude)
        radial = self.compact.joined
        self.rise_start = radial.join_X
        self.rise_end = radial.join_X*radial.outer_ratio**2
        self.restore_start = float(restore_start)
        self.restore_end = float(restore_end)
        if not self.rise_start < self.restore_start < self.restore_end:
            raise ValueError('Restoration must follow the rise start')

    def attachment_radius(self, tau):
        return self.base.attachment_radius(tau)

    def shape(self, X):
        X = np.asarray(X, float)
        up, up_prime = septic_step(
            (X-self.rise_start)/(self.rise_end-self.rise_start))
        down, down_prime = septic_step(
            (X-self.restore_start)/(self.restore_end-self.restore_start))
        S = up*(1-down)
        S_X = (up_prime*(1-down)/(self.rise_end-self.rise_start)
               -up*down_prime/(self.restore_end-self.restore_start))
        return S, S_X

    def fields(self, points, tau):
        pts = np.asarray(points, float)
        ts = np.broadcast_to(np.asarray(tau, float), (len(pts),))
        velocity, pressure = self.base.fields(pts, ts)
        if self.amplitude == 0:
            return velocity, pressure
        r = np.hypot(pts[:, 0], pts[:, 1])
        co = coordinates(r/np.sqrt(self.nu),
                         pts[:, 2]/np.sqrt(self.nu), ts, self.heat.h)
        X = np.asarray(co['X'])
        chi, chi_eta = self.compact.cutoff(co['eta'])
        active = ((X > self.rise_start) & (X < self.restore_end)
                  & (chi > 0) & (r > 0))
        if not np.any(active):
            return velocity, pressure
        xa = X[active]
        qa = np.asarray(co['q'])[active]
        ra = r[active]
        S, S_X = self.shape(xa)
        q_z = np.asarray(co['q_z'])[active]/np.sqrt(self.nu)
        eta_z = np.asarray(co['eta_z'])[active]/np.sqrt(self.nu)
        X_z = -xa*q_z/qa
        A = .5+self.heat.h
        scale = (self.amplitude*self.nu**1.5
                 *ts[active]**(1-A))
        psi_z = scale*(chi_eta[active]*eta_z*S
                       +chi[active]*S_X*X_z)
        ur = -psi_z/ra
        uz = scale*chi[active]*S_X/(self.nu*qa)
        velocity[active, 0] += ur*pts[active, 0]/ra
        velocity[active, 1] += ur*pts[active, 1]/ra
        velocity[active, 2] += uz
        return velocity, pressure
