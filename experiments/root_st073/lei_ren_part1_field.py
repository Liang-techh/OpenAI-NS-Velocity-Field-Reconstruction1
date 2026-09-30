"""Physical 3D background seed with constrained axial/radial matching.

This is a computable geometry/background candidate, not an accepted NS
solution. The axial cutoff is applied to the meridional streamfunction;
swirl is localized by the same cutoff. The radial heat tail is retained.
"""
import math
import sys
from functools import lru_cache
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'src'))
from openai_ns_reconstruction.coordinates import similarity_coordinates_from_tau
from openai_ns_reconstruction.paper_compact_field import cutoff
from lei_ren_part1_pressure_core import load_core
from lei_ren_part1_axial_match import AxialMatchedProfile


class PartIBackgroundField:
    def __init__(self, nu=.01, z_inner=.25, z_outer=.5):
        if not all(math.isfinite(v) for v in (nu, z_inner, z_outer)) or not (nu > 0 and 0 < z_inner < z_outer):
            raise ValueError('Require nu>0 and 0<z_inner<z_outer')
        self.nu, self.z_inner, self.z_outer = nu, z_inner, z_outer
        self.core, self.joined = load_core()
        self.matched = AxialMatchedProfile(self.core, self.joined)
        self.profile = self.matched.profile
        self.h = self.core.reference.h

    def z_cutoff(self, z):
        value, derivative = cutoff(z*z, self.z_inner**2, self.z_outer**2)
        return value, 2*z*derivative

    def from_similarity(self, R, Z, tau):
        if R < 0 or abs(Z) >= 1 or tau <= 0:
            raise ValueError('Require R>=0, |Z|<1, tau>0')
        q = tau/(1-Z*Z)
        return math.sqrt(2*self.nu*q*R), math.sqrt(self.nu)*q**(.5-self.h)*Z

    def cylindrical_from_tau(self, r, z, tau):
        if not all(math.isfinite(v) for v in (r, z, tau)) or r < 0 or tau <= 0:
            raise ValueError('Require finite r>=0, z and tau>0')
        B, _ = self.z_cutoff(z)
        if B == 0:
            return np.zeros(3)
        root_nu = math.sqrt(self.nu)
        s = similarity_coordinates_from_tau(r/root_nu, z/root_nu, tau, self.h)
        return self.cylindrical_at_chart(s.X, s.eta, s.q, z, d=s.d)

    def cylindrical_at_chart(self, R, Z, q, physical_z, *, d=None):
        """Evaluate a known chart without an inverse solve (quadrature API)."""
        B, Bz = self.z_cutoff(physical_z)
        if B == 0:
            return np.zeros(3)
        d = 1-Z*Z if d is None else d
        average = self.profile.radial_average_U(R, Z)
        flux = self.profile.radial_flux_factor(R, Z, self.h, d=d, L=1-2*self.h*Z*Z)
        ur = math.sqrt(self.nu/q)*math.sqrt(R/2)*flux
        ut = math.sqrt(self.nu)*q**(-.5-self.h)*math.sqrt(2*R)*self.profile.F(R, Z)
        uz = math.sqrt(self.nu)*q**(-.5-self.h)*self.profile.U(R, Z)
        # psi/r = nu*q^(-h)*sqrt(R/2)*average(U), regular at R=0.
        psi_over_r = self.nu*q**(-self.h)*math.sqrt(R/2)*average
        return np.array([B*ur-Bz*psi_over_r, B*ut, B*uz])

    def velocity_from_tau(self, x, y, z, tau):
        if not all(math.isfinite(v) for v in (x, y)):
            raise ValueError('Coordinates must be finite')
        r = math.hypot(x, y)
        ur, ut, uz = self.cylindrical_from_tau(r, z, tau)
        if r == 0:
            return np.array([0., 0., uz])
        return np.array([(x*ur-y*ut)/r, (y*ur+x*ut)/r, uz])

    def velocity(self, x, y, z, t):
        if not math.isfinite(t) or not 0 <= t < 1:
            raise ValueError('Require 0<=t<1; use tau interface near t=1')
        return self.velocity_from_tau(x, y, z, 1-t)

    def pressure_from_tau(self, x, y, z, tau):
        if not all(math.isfinite(v) for v in (x, y, z, tau)) or tau <= 0:
            raise ValueError('Require finite coordinates and tau>0')
        B, _ = self.z_cutoff(z)
        if B == 0:
            return 0.
        s = similarity_coordinates_from_tau(math.hypot(x, y)/math.sqrt(self.nu),
                                            z/math.sqrt(self.nu), tau, self.h)
        return self.nu*B*B*s.q**(-1-2*self.h)*self.joined.pressure(s.X, s.eta)

    def metadata(self):
        return dict(candidate_id='lr1-axial-matched-heat-background-001', nu=self.nu,
                    h=self.h, delta=2*self.h, time='T=1, tau=T-t>0',
                    support=f'|z|<{self.z_outer}; algebraically decaying radial heat tail',
                    cutoff_plateau=f'|z|<={self.z_inner}',
                    dependent_radial_component='Recovered from axial primitive and localized streamfunction',
                    pressure='Same joined swirl/heat pressure; axial localization p=B(z)^2 p_profile',
                    diagnostic_forcing='zero for reporting background residual; final smooth forcing unresolved',
                    acceptance_stage='geometry/background candidate; stress and momentum acceptance open',
                    pde_validated=False, scale_recursion_established=False,
                    full_five_moment_matching=False)


@lru_cache(maxsize=1)
def default_field():
    return PartIBackgroundField()


def velocity(x, y, z, t):
    return tuple(default_field().velocity(x, y, z, t))
