"""A single linear velocity basis exposing swirl and poloidal coefficients."""
from dataclasses import dataclass
from typing import ClassVar
import numpy as np
from .constrained_poloidal import PoloidalCandidate


@dataclass(frozen=True)
class CoupledCandidate(PoloidalCandidate):
    coefficients: tuple[float,...]=(0.,)*117
    COEFFICIENT_COUNT: ClassVar[int]=117
    FAMILY_ID: ClassVar[str]='coupled_velocity_v1'

    def __post_init__(self):
        super().__post_init__()
        if np.any(self.base.coefficients):raise ValueError("coupled base swirl coefficients must be zero; use the leading90 coefficients")

    def correction_basis(self,points,time):
        return np.concatenate((self.base.correction_basis(points,time),super().correction_basis(points,time)),axis=-1)
