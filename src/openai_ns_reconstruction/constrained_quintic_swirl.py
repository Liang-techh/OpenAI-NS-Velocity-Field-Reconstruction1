"""Degree-five temporal extension containing the cubic swirl family exactly."""
from dataclasses import dataclass
from math import comb
from typing import ClassVar
import numpy as np
from .constrained_axial_swirl import AxialSwirlCandidate,factors
from .constrained_temporal_swirl import _broadcast_physical_coordinates


def elevate_cubic(coefficients):
    old=np.column_stack((np.zeros(12),np.asarray(coefficients).reshape(12,3)))
    new=np.zeros((12,5))
    for k in range(1,6):
        for i in range(max(0,k-2),min(3,k)+1):
            new[:,k-1]+=old[:,i]*comb(3,i)*comb(2,k-i)/comb(5,k)
    return tuple(np.clip(new.ravel(),-1,1))


@dataclass(frozen=True,eq=False)
class QuinticSwirlCandidate(AxialSwirlCandidate):
    coefficients: tuple[float,...]=(0.,)*60
    FAMILY_ID: ClassVar[str]='quintic_swirl_v1'
    COEFFICIENT_COUNT: ClassVar[int]=60

    def correction_basis(self,points,time):
        x,y,z,t=_broadcast_physical_coordinates(points,time)
        p=np.stack((x,y,z),axis=-1)
        spatial=self.parent.outer_basis(p,t)[..., :,None]*(factors(p)-self._moment_ratios)[...,None,:]
        w=2*(t-.25)
        temporal=np.stack([comb(5,k)*w**k*(1-w)**(5-k) for k in range(1,6)],axis=-1)
        return (spatial[..., :, :,None]*temporal[...,None,None,:]).reshape(x.shape+(3,60))
