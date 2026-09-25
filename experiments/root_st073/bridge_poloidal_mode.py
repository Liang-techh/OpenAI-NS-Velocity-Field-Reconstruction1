"""Solenoidal poloidal bubble in the moving ST073 radial bridge.

This experimental streamfunction mode can change meridional shear without
changing the inner/outer C2 radial traces. It is not a solved momentum patch.
"""
import numpy as np
from functools import lru_cache
from numpy.polynomial import Polynomial
from numpy.polynomial.legendre import leggauss

from extended_compact_join import load_extended_heated_candidate
from joined_field import coordinates


@lru_cache(maxsize=16)
def minimum_curvature_coefficients(join_X, ratio):
    """Minimize integral |phi'''|^2 with fixed phi'' at similarity X=1."""
    y0 = (np.sqrt(1/join_X)-1)/(ratio-1)
    if not 0 < y0 < 1:
        raise ValueError('Reference X=1 must lie inside the radial bridge')
    y_poly = Polynomial([0., 1.])
    common = y_poly**4*(1-y_poly)**4*(y_poly-y0)**2
    basis = [common*(2*y_poly-1)**k for k in range(5)]
    nodes, weights = leggauss(64)
    nodes, weights = (nodes+1)/2, weights/2
    third = np.column_stack([b.deriv(3)(nodes) for b in basis])
    gram = (third*weights[:, None]).T@third
    constraint = np.array([b.deriv(2)(y0) for b in basis])
    old_second = 16*1024*y0**5*(1-y0)**5
    inverse_constraint = np.linalg.solve(gram, constraint)
    coefficients = (inverse_constraint*old_second
                    /(constraint@inverse_constraint))
    return y0, coefficients


class BridgePoloidalMode:
    def __init__(self, base=None, amplitude=0., shape='standard',
                 curvature_weight=.5):
        self.base = base if base is not None else load_extended_heated_candidate()
        self.nu = self.base.nu
        self.compact = self.base.compact
        self.heat = self.base.heat
        self.amplitude = float(amplitude)
        if shape not in ('standard', 'curvature', 'minimum_curvature', 'mixed',
                         'moment'):
            raise ValueError('Unknown poloidal bubble shape')
        self.shape = shape
        if not 0 <= curvature_weight <= 1:
            raise ValueError('Curvature weight must lie in [0,1]')
        self.curvature_weight = float(curvature_weight)
        if shape in ('minimum_curvature', 'mixed'):
            radial = self.compact.joined
            self.y0, self.curvature_coefficients = minimum_curvature_coefficients(
                radial.join_X, radial.outer_ratio)
        elif shape == 'moment':
            radial = self.compact.joined
            self.y0 = (np.sqrt(1/radial.join_X)-1)/(radial.outer_ratio-1)

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
        if self.shape in ('curvature', 'mixed'):
            y0 = (np.sqrt(1/radial.join_X)-1)/(ratio-1)
            factor = 8*(yy-y0)**2
            curvature_bubble_y = bubble_y*factor+16*bubble*(yy-y0)
            curvature_bubble = bubble*factor
            if self.shape == 'curvature':
                bubble, bubble_y = curvature_bubble, curvature_bubble_y
        if self.shape in ('minimum_curvature', 'mixed'):
            offset = yy-self.y0
            common = yy**4*(1-yy)**4
            common_y = 4*yy**3*(1-yy)**3*(1-2*yy)
            polynomial = np.polynomial.polynomial.polyval(
                2*yy-1, self.curvature_coefficients)
            polynomial_y = 2*np.polynomial.polynomial.polyval(
                2*yy-1, np.polynomial.polynomial.polyder(
                    self.curvature_coefficients))
            minimum_bubble = common*offset**2*polynomial
            minimum_bubble_y = ((common_y*offset**2+2*common*offset)*polynomial
                                +common*offset**2*polynomial_y)
            if self.shape == 'minimum_curvature':
                bubble, bubble_y = minimum_bubble, minimum_bubble_y
            else:
                weight = self.curvature_weight
                bubble = weight*curvature_bubble+(1-weight)*minimum_bubble
                bubble_y = (weight*curvature_bubble_y
                            +(1-weight)*minimum_bubble_y)
        elif self.shape == 'moment':
            offset = yy-self.y0
            common = yy**4*(1-yy)**4
            common_y = 4*yy**3*(1-yy)**3*(1-2*yy)
            bubble = 1024*common*offset**3
            bubble_y = 1024*(common_y*offset**3+3*common*offset**2)
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
