"""Nested localized swirl basis with six additional axially shifted rings."""
from dataclasses import dataclass
from functools import lru_cache
from typing import ClassVar
import numpy as np
from .constrained_localized_swirl import LocalizedSwirlCandidate, SPATIAL_BASIS
from .constrained_outer_momentum import AngularMomentumCandidate
from .constrained_momentum_budget import angular_moment
from .constrained_temporal_swirl import _quadrature_order, _broadcast_physical_coordinates, _bernstein_nonzero

RINGS=tuple((r,w,0.,a) for r,w,a in SPATIAL_BASIS)+tuple(
    (r,w,.8,a) for r,w in ((.16,.12),(.35,.2),(.7,.35)) for a in (.25,.5))


def factors(points):
    x,y,z=np.moveaxis(np.asarray(points),-1,0)
    return np.stack([np.exp(-((x*x+y*y-r)/w)**2-((z*z-center)/a)**2)
                     for r,w,center,a in RINGS],axis=-1)


@lru_cache(maxsize=32)
def ratios(parent,order):
    denominator=angular_moment(parent.outer_basis,.5,order)
    return tuple(angular_moment(lambda p,t,i=i:parent.outer_basis(p,t)*factors(p)[...,i,None],.5,order)/denominator for i in range(len(RINGS)))


@dataclass(frozen=True,eq=False)
class AxialSwirlCandidate(LocalizedSwirlCandidate):
    coefficients: tuple[float,...]=(0.,)*36
    FAMILY_ID: ClassVar[str]='axial_swirl_v1'
    COEFFICIENT_COUNT: ClassVar[int]=36

    def __post_init__(self):
        if not isinstance(self.parent,AngularMomentumCandidate):raise TypeError('expected angular momentum parent')
        order=_quadrature_order(self.parent.order if self.order is None else self.order)
        a=np.asarray(self.coefficients,dtype=float)
        if a.shape!=(self.COEFFICIENT_COUNT,) or not np.all(np.isfinite(a)) or np.any(np.abs(a)>1):
            raise ValueError(f'{self.COEFFICIENT_COUNT} finite coefficients in [-1,1] required')
        object.__setattr__(self,'coefficients',tuple(float(v) for v in a))
        object.__setattr__(self,'order',order)
        object.__setattr__(self,'_moment_ratios',ratios(self.parent,order))

    def correction_basis(self,points,time):
        x,y,z,t=_broadcast_physical_coordinates(points,time)
        p=np.stack((x,y,z),axis=-1)
        spatial=self.parent.outer_basis(p,t)[..., :,None]*(factors(p)-self._moment_ratios)[...,None,:]
        temporal=_bernstein_nonzero(2*(t-.25))
        return (spatial[..., :, :,None]*temporal[...,None,None,:]).reshape(x.shape+(3,36))
