"""Solenoidal poloidal bubble in the moving ST073 radial bridge.

This experimental streamfunction mode can change meridional shear without
changing the inner/outer C2 radial traces. It is not a solved momentum patch.
"""
import numpy as np

from extended_compact_join import load_extended_heated_candidate
from joined_field import coordinates


class BridgePoloidalMode:
    def __init__(self, base=None, amplitude=0., shape='standard'):
        self.base = base if base is not None else load_extended_heated_candidate()
        self.nu = self.base.nu
        self.compact = self.base.compact
        self.heat = self.base.heat
        self.amplitude = float(amplitude)
        if shape not in ('standard', 'curvature'):
            raise ValueError('Unknown poloidal bubble shape')
        self.shape = shape

    def attachment_radius(self, tau):
        return self.base.attachment_radius(tau)

    def fields(self, points, tau):
        pts = np.asarray(points, float)
        ts = np.broadcast_to(np.asarray(tau, float), (len(pts),))
        velocity, pressure = self.base.fields(pts, ts)
        if self.amplitude == 0:
            return velocity, pressure
        r = np.hypot(pts[:, 0], pts[:, 1])
        co = coordinates(r/np.sqrt(self.nu),
                         pts[:, 2]/np.sqrt(self.nu), ts, self.heat.h)
        q, eta = co['q'], co['eta']
        radial = self.compact.joined
        ri = np.sqrt(2*self.nu*q*radial.join_X)
        ratio = radial.outer_ratio
        y = (r/ri-1)/(ratio-1)
        chi, chi_eta = self.compact.cutoff(eta)
        active = (y > 0) & (y < 1) & (chi > 0)
        if not np.any(active):
            return velocity, pressure
        yy = y[active]
        bubble = 1024*yy**5*(1-yy)**5
        bubble_y = 5120*yy**4*(1-yy)**4*(1-2*yy)
        if self.shape == 'curvature':
            y0 = (np.sqrt(1/radial.join_X)-1)/(ratio-1)
            factor = 8*(yy-y0)**2
            bubble_y = bubble_y*factor+16*bubble*(yy-y0)
            bubble *= factor
        qa, ra = q[active], r[active]
        qz = co['q_z'][active]/np.sqrt(self.nu)
        etaz = co['eta_z'][active]/np.sqrt(self.nu)
        yz = -(1+(ratio-1)*yy)*qz/(2*qa*(ratio-1))
        A = .5+self.heat.h
        scale = self.nu**1.5*self.amplitude*qa**(1-A)
        psi_z = scale*(chi_eta[active]*etaz*bubble
                       +(1-A)*qz/qa*chi[active]*bubble
                       +chi[active]*bubble_y*yz)
        psi_r = scale*chi[active]*bubble_y/((ratio-1)*ri[active])
        ur, uz = -psi_z/ra, psi_r/ra
        velocity[active, 0] += ur*pts[active, 0]/ra
        velocity[active, 1] += ur*pts[active, 1]/ra
        velocity[active, 2] += uz
        return velocity, pressure
