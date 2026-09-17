"""Additional zero-moment swirl modes inside the former frozen collar."""
from dataclasses import dataclass
from functools import lru_cache
from math import comb
from typing import ClassVar
import numpy as np
from .constrained_quintic_swirl import QuinticSwirlCandidate
from .constrained_temporal_swirl import _broadcast_physical_coordinates
from .constrained_force import compact_bump
from .constrained_momentum_budget import angular_moment
from .cutoffs import standard_cutoff

INNER_RINGS=tuple((r,w,z,.5) for r,w in ((.04,.035),(.1,.06),(.2,.1)) for z in (0.,.8))


def inner_raw(points):
    x,y,z=np.moveaxis(np.asarray(points),-1,0);r2=x*x+y*y
    guard=1-np.fromiter((standard_cutoff(float(v)) for v in (r2/.02).ravel()),float,count=r2.size).reshape(r2.shape)
    envelope=compact_bump(r2/4)[0]*compact_bump(z*z/4)[0]*guard
    swirl=np.stack((-y*envelope,x*envelope,np.zeros_like(x)),axis=-1)
    factors=np.stack([np.exp(-((r2-r)/w)**2-((z*z-center)/a)**2) for r,w,center,a in INNER_RINGS],axis=-1)
    return swirl[..., :,None]*factors[...,None,:]


@lru_cache(maxsize=32)
def inner_ratios(parent,order):
    denominator=angular_moment(parent.outer_basis,.5,order)
    return tuple(angular_moment(lambda p,t,i=i:inner_raw(p)[...,i],.5,order)/denominator for i in range(6))


@dataclass(frozen=True,eq=False)
class InnerSwirlCandidate(QuinticSwirlCandidate):
    coefficients: tuple[float,...]=(0.,)*90
    FAMILY_ID: ClassVar[str]='inner_swirl_v1'
    COEFFICIENT_COUNT: ClassVar[int]=90

    def __post_init__(self):
        super().__post_init__()
        object.__setattr__(self,'_inner_ratios',inner_ratios(self.parent,self.order))

    def correction_basis(self,points,time):
        original=super().correction_basis(points,time)
        x,y,z,t=_broadcast_physical_coordinates(points,time);p=np.stack((x,y,z),axis=-1)
        modes=inner_raw(p)-self.parent.outer_basis(p,t)[..., :,None]*np.asarray(self._inner_ratios)
        w=2*(t-.25);temporal=np.stack([comb(5,k)*w**k*(1-w)**(5-k) for k in range(1,6)],axis=-1)
        additional=(modes[..., :, :,None]*temporal[...,None,None,:]).reshape(x.shape+(3,30))
        return np.concatenate((original,additional),axis=-1)
