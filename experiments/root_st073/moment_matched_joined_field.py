"""Diagnostic two-moment field; pressure and five-moment interface remain open."""
import json

import numpy as np
from numpy.polynomial.chebyshev import chebder, chebfit, chebval

from joined_field import JoinedField, ROOT, coordinates
from paper_moment_bridge import bump, null_bump, even_bump, slice_data


class MomentMatchedJoinedField:
    def __init__(self, base=None, null_amplitude=0.,
                 meridional_null_amplitude=0.,
                 meridional_even_amplitude=0., patch=None):
        self.base = base if base is not None else JoinedField()
        self.nu = self.base.nu
        self.inner = self.base.inner
        self.c = self.base.c
        self.null_amplitude = float(null_amplitude)
        self.meridional_null_amplitude = float(meridional_null_amplitude)
        self.meridional_even_amplitude = float(meridional_even_amplitude)
        report = json.loads((ROOT/'paper_moment_coefficients.json').read_text())
        self.patch_start, self.patch_end = (
            tuple(map(float, patch)) if patch is not None
            else (report['patch_start'], report['patch_end']))
        if self.patch_start < 4. or self.patch_end <= self.patch_start:
            raise ValueError('remote patch must lie in the heat exterior')
        self.eta_scale = report['eta_chebyshev_scale']
        self.swirl_coefficients = np.array(report['swirl_chebyshev_coefficients'])
        self.meridional_coefficients = np.array(report['meridional_chebyshev_coefficients'])
        if (self.null_amplitude or self.meridional_null_amplitude
                or self.meridional_even_amplitude or patch is not None):
            training = [slice_data(self.base, float(eta), tau, 32,
                                   self.patch_start, self.patch_end,
                                   self.null_amplitude,
                                   self.meridional_null_amplitude,
                                   self.meridional_even_amplitude)
                        for tau in (.5/64, .032, .128)
                        for eta in np.linspace(-.4, .4, 9)]
            if any(row['meridional_amplitude'] is None for row in training):
                raise ValueError('moment-null shape has no real meridional root')
            eta = np.array([row['eta'] for row in training])
            fit_degree = 6 if self.meridional_even_amplitude else 4
            self.swirl_coefficients = chebfit(
                eta/self.eta_scale,
                [row['swirl_amplitude'] for row in training], fit_degree)
            self.meridional_coefficients = chebfit(
                eta/self.eta_scale,
                [row['meridional_amplitude'] for row in training], fit_degree)
        self.meridional_derivative = chebder(self.meridional_coefficients)

    def fields(self, points, tau):
        pts = np.asarray(points, float)
        velocity, pressure = self.base.fields(pts, tau)
        source = pts/np.sqrt(self.nu)
        radius = np.hypot(pts[:, 0], pts[:, 1])
        source_radius = radius/np.sqrt(self.nu)
        co = coordinates(source_radius, source[:, 2], tau, self.inner.h)
        X = np.asarray(co['X'])
        active = (X > self.patch_start) & (X < self.patch_end)
        if not np.any(active):
            return velocity, pressure
        q = np.asarray(co['q'])[active]
        eta = np.asarray(co['eta'])[active]
        b, bx = bump(X[active], self.patch_start, self.patch_end)
        null, null_x = null_bump(X[active], self.patch_start, self.patch_end)
        even, even_x = even_bump(X[active], self.patch_start, self.patch_end)
        aE = chebval(eta/self.eta_scale, self.swirl_coefficients)
        aU = chebval(eta/self.eta_scale, self.meridional_coefficients)
        aU_eta = (chebval(eta/self.eta_scale,
                           self.meridional_derivative)/self.eta_scale)
        A = .5+self.inner.h
        delta_theta = np.sqrt(self.nu)*q**(-A)*(aE*b+self.null_amplitude*null)
        delta_axial = np.sqrt(self.nu)*q**(-A)*(
            aU*bx+self.meridional_null_amplitude*null_x
            +self.meridional_even_amplitude*even_x)
        streamfunction_shape = (aU*b+self.meridional_null_amplitude*null
                                +self.meridional_even_amplitude*even)
        axial_shape = (aU*bx+self.meridional_null_amplitude*null_x
                       +self.meridional_even_amplitude*even_x)
        bracket = ((1-A)*np.asarray(co['q_z'])[active]/q*streamfunction_shape
                   + aU_eta*np.asarray(co['eta_z'])[active]*b
                   + axial_shape*np.asarray(co['X_z'])[active])
        delta_radial = (-self.nu*q**(1-A)/radius[active])*bracket
        ca = pts[active, 0]/radius[active]
        sa = pts[active, 1]/radius[active]
        velocity[active, 0] += delta_radial*ca-delta_theta*sa
        velocity[active, 1] += delta_radial*sa+delta_theta*ca
        velocity[active, 2] += delta_axial
        return velocity, pressure
