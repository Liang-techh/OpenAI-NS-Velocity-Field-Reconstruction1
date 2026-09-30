"""Physical extended Part I candidate with shared pressure and axial repair.

Uses the same implicit physical chart and streamfunction localization as
the retained short-join seed. The whole stress cone, recursive corrections
and smooth final forcing remain open.
"""
import math
from functools import lru_cache
from pathlib import Path

from lei_ren_part1_field import PartIBackgroundField


class ExtendedPartIBackgroundField(PartIBackgroundField):
    def __init__(self, nu=.01, z_inner=.25, z_outer=.5, *, matched=None):
        if not all(math.isfinite(v) for v in (nu,z_inner,z_outer)) or not (nu>0 and 0<z_inner<z_outer):
            raise ValueError('Require nu>0 and 0<z_inner<z_outer')
        if matched is None:
            from lei_ren_part1_smooth_extended_axial import SmoothExtendedAxial
            seed=Path(__file__).with_name('lei_ren_part1_smooth_extended_axial_seed.json')
            matched=(SmoothExtendedAxial.from_seed(seed) if seed.is_file()
                     else SmoothExtendedAxial())
        self.nu,self.z_inner,self.z_outer=nu,z_inner,z_outer
        self.matched=matched
        self.joined=matched.swirl
        self.core=self.joined.core
        self.profile=matched.profile
        self.h=self.core.reference.h

    def metadata(self):
        return {
            'candidate_id':'lr1-extended-common-pressure-axial-background-001',
            'nu':self.nu,'h':self.h,'delta':2*self.h,'time':'T=1, tau=T-t>0',
            'support':f'|z|<{self.z_outer}; infinite radial heat tail',
            'core_radius':self.joined.R_core,
            'collar_inner_radius':self.joined.collar.collar_inner_radius,
            'heat_radius':self.joined.collar.R_b,
            'dependent_radial_component':'Recovered from same repaired axial primitive; localized streamfunction',
            'pressure':'Integral of same extended swirl, normalized by actual heat tail; physical p=B(z)^2 p_profile',
            'diagnostic_forcing':'zero; final smooth forcing unresolved',
            'acceptance_stage':'extended background candidate; whole cone and recursion remain open',
            'profile':self.matched.metadata(),
            'pde_validated':False,'scale_recursion_established':False,
            'whole_stress_cone_validated':False,
        }


@lru_cache(maxsize=1)
def default_field():
    """Reuse one assembled profile for repeated physical velocity calls."""
    return ExtendedPartIBackgroundField()


def velocity(x,y,z,t):
    return tuple(default_field().velocity(x,y,z,t))


def velocity_from_tau(x,y,z,tau):
    return tuple(default_field().velocity_from_tau(x,y,z,tau))
