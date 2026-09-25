"""Remote-moment patched joined field, before stress realization and localization."""
import json

import numpy as np
from numpy.polynomial.chebyshev import chebder, chebval

from joined_field import JoinedField, ROOT, coordinates
from paper_moment_bridge import bump


class MomentMatchedJoinedField:
    def __init__(self, base=None):
        self.base = base if base is not None else JoinedField()
        self.nu = self.base.nu
        self.inner = self.base.inner
        self.c = self.base.c
        report = json.loads((ROOT/'paper_moment_coefficients.json').read_text())
        self.patch_start = report['patch_start']
        self.patch_end = report['patch_end']
        self.eta_scale = report['eta_chebyshev_scale']
        self.swirl_coefficients = np.array(report['swirl_chebyshev_coefficients'])
        self.meridional_coefficients = np.array(report['meridional_chebyshev_coefficients'])
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
        aE = chebval(eta/self.eta_scale, self.swirl_coefficients)
        aU = chebval(eta/self.eta_scale, self.meridional_coefficients)
        aU_eta = (chebval(eta/self.eta_scale,
                           self.meridional_derivative)/self.eta_scale)
        A = .5+self.inner.h
        delta_theta = np.sqrt(self.nu)*q**(-A)*aE*b
        delta_axial = np.sqrt(self.nu)*q**(-A)*aU*bx
        bracket = ((1-A)*np.asarray(co['q_z'])[active]/q*aU*b
                   + aU_eta*np.asarray(co['eta_z'])[active]*b
                   + aU*bx*np.asarray(co['X_z'])[active])
        delta_radial = (-self.nu*q**(1-A)/radius[active])*bracket
        ca = pts[active, 0]/radius[active]
        sa = pts[active, 1]/radius[active]
        velocity[active, 0] += delta_radial*ca-delta_theta*sa
        velocity[active, 1] += delta_radial*sa+delta_theta*ca
        velocity[active, 2] += delta_axial
        return velocity, pressure

