"""Compact solenoidal physical lift of the two fixed-slice moment repairs.

The max-slack coefficients at eta=.2,.3 are smoothly interpolated on a
positive axial band and tapered to zero outside it. The azimuthal change
is divergence-free by axisymmetry; the meridional change is the curl of
a Stokes streamfunction. Exact five-moment matching is claimed only at
the two fitted slices and reference time, not throughout the band.
"""
import json

import numpy as np
from numpy.polynomial.legendre import leggauss

from azimuthal_capacity_optimize import INTERVALS
from coupled_five_moment_slice import correction_modes
from high_frequency_shear_screen import make_field
from joined_field import coordinates
from paper_moment_bridge import bump
from radial_continuation import ROOT
from radial_moment_step import RadialMomentStep, septic_step


class CoupledMomentPhysicalLift:
    def __init__(self, base=None, quadrature_order=24):
        self.base = base if base is not None else RadialMomentStep(
            make_field(16, 2.), -.5)
        self.compact = self.base.compact
        self.heat = self.base.heat
        self.nu = self.base.nu
        self.A = .5+self.heat.h
        data = json.loads((ROOT/'coupled_five_moment_slice.json').read_text())
        chosen = sorted((row for row in data['rows']
                         if row['variant'] == 'maximum_slack'
                         and row.get('five_moments_restored')),
                        key=lambda row: row['eta'])
        if len(chosen) != 2 or [row['eta'] for row in chosen] != [.2, .3]:
            raise ValueError('Two fitted maximum-slack slices required')
        self.e_rows = np.array([row['e_coefficients'] for row in chosen])
        self.u_rows = np.array([row['u_coefficients'] for row in chosen])
        self.nodes, self.weights = leggauss(quadrature_order)
        self.mass = self.primitive_basis(3.)

    def attachment_radius(self, tau):
        return self.base.attachment_radius(tau)

    def primitive_basis(self, X):
        X = float(np.clip(X, 1., 3.))
        total = np.zeros(12)
        for lo, hi in ((1., 1.02), (1.02, 2.98), (2.98, 3.)):
            end = min(X, hi)
            if end <= lo:
                continue
            nodes = (lo+end)/2+(end-lo)*self.nodes/2
            weights = (end-lo)*self.weights/2
            total += correction_modes(nodes)@weights
        return total

    def coefficients(self, eta):
        eta = float(eta)
        rise, rise_d = septic_step((eta-.1)/.1)
        fall, fall_d = septic_step((eta-.3)/.1)
        rho = float(rise*(1-fall))
        rho_d = float(rise_d*(1-fall)/.1-rise*fall_d/.1)
        blend, blend_d = septic_step((eta-.2)/.1)
        blend = float(blend)
        blend_d = float(blend_d)/.1
        E = (1-blend)*self.e_rows[0]+blend*self.e_rows[1]
        U = (1-blend)*self.u_rows[0]+blend*self.u_rows[1]
        Ed = blend_d*(self.e_rows[1]-self.e_rows[0])
        Ud = blend_d*(self.u_rows[1]-self.u_rows[0])
        E, Ed = rho*E, rho_d*E+rho*Ed
        U, Ud = rho*U, rho_d*U+rho*Ud
        # Enforce exact zero total U correction for every interpolated eta.
        mass_norm = float(self.mass@self.mass)
        projection = float(U@self.mass)/mass_norm
        projection_d = float(Ud@self.mass)/mass_norm
        U = U-projection*self.mass
        Ud = Ud-projection_d*self.mass
        return E, Ed, U, Ud

    def fields(self, points, tau):
        pts = np.asarray(points, float)
        ts = np.broadcast_to(np.asarray(tau, float), (len(pts),))
        velocity, pressure = self.base.fields(pts, ts)
        r = np.hypot(pts[:, 0], pts[:, 1])
        sn = np.sqrt(self.nu)
        co = coordinates(r/sn, pts[:, 2]/sn, ts, self.heat.h)
        for i, (radius, X, eta, q) in enumerate(zip(
                r, co['X'], co['eta'], co['q'])):
            if not (radius > 0 and 1. < X < 3. and .1 < eta < .4):
                continue
            eco, _, uco, uco_d = self.coefficients(eta)
            e_delta = sum(eco[j]*bump(np.array([X]), *interval)[0][0]
                          for j, interval in enumerate(INTERVALS))
            basis = correction_modes(np.array([X]))[:, 0]
            primitive = self.primitive_basis(X)
            G = float(uco@primitive)
            G_x = float(uco@basis)
            G_eta = float(uco_d@primitive)
            q_z = float(co['q_z'][i])/sn
            X_z = float(co['X_z'][i])/sn
            eta_z = float(co['eta_z'][i])/sn
            psi_z = (self.nu**1.5*((1-self.A)*q**(-self.A)*q_z*G
                     +q**(1-self.A)*(G_x*X_z+G_eta*eta_z)))
            ur = -psi_z/radius
            uz = sn*q**(-self.A)*G_x
            ut = sn*q**(-self.A)*e_delta
            ca, sa = pts[i, 0]/radius, pts[i, 1]/radius
            velocity[i] += np.array([ur*ca-ut*sa, ur*sa+ut*ca, uz])
        return velocity, pressure
