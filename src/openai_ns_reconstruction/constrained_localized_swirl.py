"""Bounded, zero-moment localized pure-swirl corrections.

``LocalizedSwirlCandidate`` freezes an ``AngularMomentumCandidate`` parent and
adds eighteen bounded coefficients.  The spatial modes are six fixed
Gaussian rings, each multiplied by the parent's compact outer swirl basis and
orthogonalized against that basis at the selected angular-moment quadrature.
The three nonzero cubic Bernstein factors preserve the initial field.

The moment cancellation is numerical evidence for the selected quadrature;
the Cartesian pure-swirl form is divergence free by construction.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
from functools import lru_cache
import json
from numbers import Integral
from pathlib import Path
from typing import ClassVar

import numpy as np

from .constrained_force import RestrictedForce
from .constrained_momentum_budget import angular_moment
from .constrained_outer_momentum import AngularMomentumCandidate
from .constrained_tensor_candidate import TensorCandidate
from .constrained_temporal_swirl import (
    TEMPORAL_BASIS,
    _bernstein_nonzero,
    _broadcast_physical_coordinates,
    _parent_from_payload,
    _parent_payload,
    _quadrature_order,
)
from .constrained_candidate import CompactCandidate
from .quadrature import unit_rule


SPATIAL_BASIS = (
    (0.22, 0.14, 0.25),
    (0.22, 0.14, 0.8),
    (0.5, 0.3, 0.25),
    (0.5, 0.3, 0.8),
    (1.0, 0.5, 0.25),
    (1.0, 0.5, 0.8),
)
"""Fixed ``(center, width, axial_width)`` parameters for the six rings."""

# Descriptive aliases make the fixed geometry easy to discover without
# changing the ordering used by the optimizer or serialized artifacts.
RING_PARAMETERS = SPATIAL_BASIS
LOCALIZED_BASIS = SPATIAL_BASIS

COEFFICIENT_COUNT = len(SPATIAL_BASIS) * len(TEMPORAL_BASIS)
ZERO_COEFFICIENTS = (0.0,) * COEFFICIENT_COUNT
FAMILY_ID = "localized_swirl_v1"


def _coefficient_tuple(value) -> tuple[float, ...]:
    """Validate and freeze the bounded length-eighteen vector."""
    try:
        array = np.asarray(value, dtype=float)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(
            f"coefficients must be a finite length-{COEFFICIENT_COUNT} sequence"
        ) from exc
    if array.ndim != 1 or array.shape[0] != COEFFICIENT_COUNT:
        raise ValueError(f"coefficients must have length {COEFFICIENT_COUNT}")
    if not np.all(np.isfinite(array)) or np.any((array < -1.0) | (array > 1.0)):
        raise ValueError("coefficients entries must be finite and in [-1, 1]")
    return tuple(float(entry) for entry in array)


@lru_cache(maxsize=32)
def _localized_moment_ratios(
    parent: AngularMomentumCandidate, order: int
) -> tuple[float, ...]:
    """Return raw-ring/parent angular-moment ratios.

    ``AngularMomentumCandidate`` uses identity hashing, so this cache is per
    parent instance and does not modify the parent or its derived moment data.
    """
    parent_moment = angular_moment(parent.outer_basis, 0.5, order)
    if not np.isfinite(parent_moment) or parent_moment <= 0.0:
        raise ValueError("parent outer basis has a nonpositive angular moment")

    ratios = []
    for center, width, axial_width in SPATIAL_BASIS:

        def raw_mode(points, time, c=center, w=width, aw=axial_width):
            points = np.asarray(points, dtype=float)
            x, y, z = np.moveaxis(points, -1, 0)
            radius_sq = x * x + y * y
            factor = np.exp(
                -((radius_sq - c) / w) ** 2 - ((z * z) / aw) ** 2
            )
            return parent.outer_basis(points, time) * factor[..., None]

        raw_moment = angular_moment(raw_mode, 0.5, order)
        if not np.isfinite(raw_moment):
            raise ValueError("nonfinite angular moment for localized swirl mode")
        ratios.append(float(raw_moment / parent_moment))
    return tuple(ratios)


def _base_family(base) -> str:
    """Return the supported serialized family for a parent base."""
    if isinstance(base, TensorCandidate):
        return TensorCandidate.FAMILY_ID
    if isinstance(base, CompactCandidate):
        return "compact_axisymmetric_window_v4"
    raise ValueError("localized swirl parent base must be a compact candidate")


@dataclass(frozen=True, eq=False)
class LocalizedSwirlCandidate:
    """Angular-momentum-preserving bounded localized swirl wrapper."""

    parent: AngularMomentumCandidate
    coefficients: tuple[float, ...] = ZERO_COEFFICIENTS
    order: int | None = None

    SPATIAL_BASIS: ClassVar[tuple[tuple[float, float, float], ...]] = SPATIAL_BASIS
    RING_PARAMETERS: ClassVar[tuple[tuple[float, float, float], ...]] = RING_PARAMETERS
    LOCALIZED_BASIS: ClassVar[tuple[tuple[float, float, float], ...]] = LOCALIZED_BASIS
    TEMPORAL_BASIS: ClassVar[tuple[int, ...]] = TEMPORAL_BASIS
    COEFFICIENT_COUNT: ClassVar[int] = COEFFICIENT_COUNT
    FAMILY_ID: ClassVar[str] = FAMILY_ID

    _moment_ratios: tuple[float, ...] = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        if not isinstance(self.parent, AngularMomentumCandidate):
            raise TypeError("parent must be an AngularMomentumCandidate")
        selected_order = (
            int(self.parent.order) if self.order is None else _quadrature_order(self.order)
        )
        selected_order = _quadrature_order(selected_order)
        object.__setattr__(self, "order", selected_order)
        object.__setattr__(self, "coefficients", _coefficient_tuple(self.coefficients))
        object.__setattr__(
            self,
            "_moment_ratios",
            _localized_moment_ratios(self.parent, selected_order),
        )

    def __eq__(self, other) -> bool:
        if not isinstance(other, LocalizedSwirlCandidate):
            return NotImplemented
        return (
            self.coefficients == other.coefficients
            and self.order == other.order
            and _parent_payload(self.parent) == _parent_payload(other.parent)
        )

    @property
    def force(self):
        """The prescribed force delegated from the parent candidate."""
        return self.parent.force

    def correction_basis(self, points, time) -> np.ndarray:
        """Return the physical correction basis with shape ``(...,3,18)``."""
        x, y, z, time = _broadcast_physical_coordinates(points, time)
        broadcast_points = np.stack((x, y, z), axis=-1)
        outer = self.parent.outer_basis(broadcast_points, time)
        radius_sq = x * x + y * y
        axial_sq = z * z
        raw_factors = np.stack(
            [
                np.exp(
                    -((radius_sq - center) / width) ** 2
                    - (axial_sq / axial_width) ** 2
                )
                for center, width, axial_width in self.SPATIAL_BASIS
            ],
            axis=-1,
        )
        modes = outer[..., :, None] * (
            raw_factors - np.asarray(self._moment_ratios, dtype=float)
        )[..., None, :]
        temporal = _bernstein_nonzero(2.0 * (time - 0.25))
        combined = modes[..., :, :, None] * temporal[..., None, None, :]
        return combined.reshape(x.shape + (3, COEFFICIENT_COUNT))

    def velocity(self, points, time):
        """Evaluate parent velocity plus the bounded localized correction."""
        parent_velocity = self.parent.velocity(points, time)
        basis = self.correction_basis(points, time)
        return parent_velocity + np.einsum(
            "...cm,m->...c", basis, np.asarray(self.coefficients, dtype=float)
        )

    def pressure(self, points, time):
        """Delegate pressure and its gauge to the parent candidate."""
        return self.parent.pressure(points, time)

    def energy(self, time=0.25, order=96):
        """Numerically integrate kinetic energy with the shared unit rule."""
        time_array = np.asarray(time, dtype=float)
        if time_array.ndim != 0 or not np.isfinite(time_array):
            raise ValueError("energy time must be a finite scalar")
        if not 0.25 <= float(time_array) <= 0.75:
            raise ValueError("candidate time domain is [0.25,0.75]")
        order = _quadrature_order(order)
        nodes, weights = unit_rule(order)
        radius = 2.0 * nodes
        axial = 4.0 * nodes - 2.0
        rr, zz = np.meshgrid(radius, axial, indexing="ij")
        points = np.stack((rr, np.zeros_like(rr), zz), axis=-1)
        speed2 = np.sum(self.velocity(points, float(time_array)) ** 2, axis=-1)
        return float(
            8.0
            * np.pi
            * np.sum(rr * speed2 * weights[:, None] * weights[None, :])
        )

    def save(self, path) -> None:
        """Save coefficients and a complete serialized parent candidate."""
        payload = {
            "schema_version": 1,
            "status": "candidate",
            "paper_exact": False,
            "family": self.FAMILY_ID,
            "parameters": {
                "coefficients": list(self.coefficients),
                "quadrature_order": self.order,
                "parent": _parent_payload(self.parent),
            },
        }
        Path(path).write_text(
            json.dumps(payload, indent=2) + "\n", encoding="utf-8"
        )

    @classmethod
    def load(cls, path):
        """Load a self-contained localized swirl artifact."""
        try:
            data = json.loads(Path(path).read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise ValueError("invalid localized swirl artifact") from exc
        if data.get("schema_version") != 1 or data.get("family") != cls.FAMILY_ID:
            raise ValueError("unsupported localized swirl candidate artifact")
        parameters = data.get("parameters")
        if not isinstance(parameters, dict):
            raise ValueError("localized swirl parameters must be an object")
        parent_data = parameters.get("parent", data.get("parent"))
        if parent_data is None:
            parent_data = parameters.get("parent_parameters")
        parent = _parent_from_payload(parent_data)
        coefficients = parameters.get("coefficients", data.get("coefficients"))
        if coefficients is None:
            raise ValueError("localized swirl coefficients are missing")
        order = parameters.get("quadrature_order", parent.order)
        return cls(parent, coefficients, order)


__all__ = [
    "COEFFICIENT_COUNT",
    "FAMILY_ID",
    "LOCALIZED_BASIS",
    "RING_PARAMETERS",
    "SPATIAL_BASIS",
    "TEMPORAL_BASIS",
    "ZERO_COEFFICIENTS",
    "LocalizedSwirlCandidate",
]
