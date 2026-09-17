"""Bounded, zero-moment temporal pure-swirl corrections.

``TemporalSwirlCandidate`` wraps :class:`AngularMomentumCandidate` and adds
fifteen bounded coefficients.  The five spatial factors are fixed in
physical coordinates and the three nonzero cubic Bernstein factors provide
the temporal degrees of freedom.  Each spatial factor is orthogonal to the
parent outer basis for the numerical angular-moment quadrature used by the
candidate, so the correction does not change the parent torque budget.

The quadrature orthogonalization is numerical evidence for the selected
quadrature order; it is not an interval-certified identity.  The spatial
factor itself is a smooth compactly supported Cartesian pure-swirl field,
which makes its divergence vanish identically away from the numerical
evaluation details.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
import json
from numbers import Integral
from pathlib import Path
from typing import ClassVar

import numpy as np

from .constrained_candidate import CompactCandidate
from .constrained_force import RestrictedForce
from .constrained_momentum_budget import angular_moment
from .constrained_outer_momentum import AngularMomentumCandidate
from .constrained_tensor_candidate import TensorCandidate
from .quadrature import unit_rule


SPATIAL_BASIS = ((1, 0), (0, 1), (2, 0), (1, 1), (0, 2))
"""The five ``(r^2/4)^i (z^2/4)^j`` spatial factors."""

TEMPORAL_BASIS = (1, 2, 3)
"""The nonzero cubic Bernstein indices; index zero preserves initial data."""

COEFFICIENT_COUNT = len(SPATIAL_BASIS) * len(TEMPORAL_BASIS)
ZERO_COEFFICIENTS = (0.0,) * COEFFICIENT_COUNT
FAMILY_ID = "temporal_swirl_v1"


def _coefficient_tuple(value) -> tuple[float, ...]:
    """Validate and freeze the bounded length-fifteen coefficient vector."""
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


def _quadrature_order(value: int) -> int:
    """Validate the order accepted by :func:`unit_rule`."""
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise ValueError("quadrature order must be an integer in [2, 2048]")
    value = int(value)
    if not 2 <= value <= 2048:
        raise ValueError("quadrature order must be an integer in [2, 2048]")
    # This also keeps the validation tied to the actual shared rule.
    unit_rule(value)
    return value


def _zero_moment_ratios(
    parent: AngularMomentumCandidate, order: int
) -> tuple[float, ...]:
    """Return raw-mode/parent-moment ratios at the selected quadrature order."""
    parent_moment = angular_moment(parent.outer_basis, 0.5, order)
    if not np.isfinite(parent_moment) or parent_moment <= 0.0:
        raise ValueError("parent outer basis has a nonpositive angular moment")

    ratios = []
    for radial_degree, axial_degree in SPATIAL_BASIS:
        def raw_mode(points, time, i=radial_degree, j=axial_degree):
            points = np.asarray(points, dtype=float)
            x, y, z = np.moveaxis(points, -1, 0)
            radial = (x * x + y * y) / 4.0
            axial = (z * z) / 4.0
            factor = radial**i * axial**j
            return parent.outer_basis(points, time) * factor[..., None]

        raw_moment = angular_moment(raw_mode, 0.5, order)
        if not np.isfinite(raw_moment):
            raise ValueError("nonfinite angular moment for temporal swirl mode")
        ratios.append(float(raw_moment / parent_moment))
    return tuple(ratios)


def _broadcast_physical_coordinates(points, time):
    """Validate inputs and broadcast physical coordinates with pointwise time."""
    points = np.asarray(points, dtype=float)
    time = np.asarray(time, dtype=float)
    if points.ndim < 1 or points.shape[-1] != 3:
        raise ValueError("points must have shape (...,3)")
    if not np.all(np.isfinite(points)) or not np.all(np.isfinite(time)):
        raise ValueError("coordinates and time must be finite")
    if np.any((time < 0.25) | (time > 0.75)):
        raise ValueError("candidate time domain is [0.25,0.75]")
    x, y, z, time = np.broadcast_arrays(*np.moveaxis(points, -1, 0), time)
    return x, y, z, time


def _bernstein_nonzero(w: np.ndarray) -> np.ndarray:
    """Evaluate cubic Bernstein factors ``k=1,2,3`` on ``w``."""
    one_minus_w = 1.0 - w
    return np.stack(
        (
            3.0 * w * one_minus_w**2,
            3.0 * w**2 * one_minus_w,
            w**3,
        ),
        axis=-1,
    )


def _base_family(base) -> str:
    if isinstance(base, TensorCandidate):
        return TensorCandidate.FAMILY_ID
    if isinstance(base, CompactCandidate):
        return "compact_axisymmetric_window_v4"
    raise ValueError("temporal swirl parent base must be a compact candidate")


def _parent_payload(parent: AngularMomentumCandidate) -> dict:
    """Build a JSON-safe, self-contained representation of the parent."""
    return {
        "schema_version": 1,
        "family": "angular_momentum_outer_v1",
        "base_family": _base_family(parent.base),
        "base_parameters": asdict(parent.base),
        "force": asdict(parent.force),
        "quadrature_order": int(parent.order),
        "outer_shape": list(parent.outer_shape),
        "pressure_coefficients": list(parent.pressure_coefficients),
    }


def _parent_from_payload(payload: dict) -> AngularMomentumCandidate:
    """Restore an AngularMomentumCandidate from the nested save payload."""
    if not isinstance(payload, dict):
        raise ValueError("temporal swirl parent must be an object")
    # Accept the parent artifact's own top-level shape when it is embedded by
    # a caller, while preferring the explicit nested representation we save.
    if isinstance(payload.get("parameters"), dict):
        payload = {**payload, **payload["parameters"]}
    base_parameters = payload.get("base_parameters")
    force_parameters = payload.get("force")
    if not isinstance(base_parameters, dict) or not isinstance(force_parameters, dict):
        raise ValueError("temporal swirl parent parameters are incomplete")

    base_family = payload.get("base_family")
    if base_family == TensorCandidate.FAMILY_ID or "poloidal_coefficients" in base_parameters:
        base = TensorCandidate(**base_parameters)
    elif base_family in (None, "compact_axisymmetric_window_v4"):
        base = CompactCandidate(**base_parameters)
    else:
        raise ValueError("unsupported temporal swirl parent base artifact")
    force = RestrictedForce(**force_parameters)
    try:
        order = payload.get("quadrature_order", 96)
        outer_shape = payload.get("outer_shape", (0.0, 0.0))
        pressure = payload.get("pressure_coefficients", (0.0,) * 18)
        return AngularMomentumCandidate(base, force, order, outer_shape, pressure)
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("invalid temporal swirl parent artifact") from exc


@dataclass(frozen=True, eq=False)
class TemporalSwirlCandidate:
    """Angular-momentum-preserving bounded temporal swirl wrapper."""

    parent: AngularMomentumCandidate
    coefficients: tuple[float, ...] = ZERO_COEFFICIENTS
    order: int | None = None

    SPATIAL_BASIS: ClassVar[tuple[tuple[int, int], ...]] = SPATIAL_BASIS
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
            self, "_moment_ratios", _zero_moment_ratios(self.parent, selected_order)
        )

    def __eq__(self, other) -> bool:
        if not isinstance(other, TemporalSwirlCandidate):
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
        """Return the physical correction basis with shape ``(...,3,15)``."""
        x, y, z, time = _broadcast_physical_coordinates(points, time)
        broadcast_points = np.stack((x, y, z), axis=-1)
        outer = self.parent.outer_basis(broadcast_points, time)
        radial = (x * x + y * y) / 4.0
        axial = (z * z) / 4.0
        raw_factors = np.stack(
            [radial**i * axial**j for i, j in self.SPATIAL_BASIS], axis=-1
        )
        modes = outer[..., :, None] * (
            raw_factors - np.asarray(self._moment_ratios, dtype=float)
        )[..., None, :]
        temporal = _bernstein_nonzero(2.0 * (time - 0.25))
        combined = modes[..., :, :, None] * temporal[..., None, None, :]
        return combined.reshape(x.shape + (3, self.COEFFICIENT_COUNT))

    def velocity(self, points, time):
        """Evaluate parent velocity plus the bounded temporal swirl correction."""
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
        parent = _parent_payload(self.parent)
        payload = {
            "schema_version": 1,
            "status": "candidate",
            "paper_exact": False,
            "family": self.FAMILY_ID,
            "parameters": {
                "coefficients": list(self.coefficients),
                "quadrature_order": self.order,
                "parent": parent,
            },
        }
        Path(path).write_text(
            json.dumps(payload, indent=2) + "\n", encoding="utf-8"
        )

    @classmethod
    def load(cls, path):
        """Load a self-contained temporal swirl artifact."""
        try:
            data = json.loads(Path(path).read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise ValueError("invalid temporal swirl artifact") from exc
        if data.get("schema_version") != 1 or data.get("family") != cls.FAMILY_ID:
            raise ValueError("unsupported temporal swirl candidate artifact")
        parameters = data.get("parameters")
        if not isinstance(parameters, dict):
            raise ValueError("temporal swirl parameters must be an object")
        parent_data = parameters.get("parent", data.get("parent"))
        if parent_data is None:
            parent_data = parameters.get("parent_parameters")
        parent = _parent_from_payload(parent_data)
        coefficients = parameters.get("coefficients", data.get("coefficients"))
        if coefficients is None:
            raise ValueError("temporal swirl coefficients are missing")
        order = parameters.get("quadrature_order", parent.order)
        return cls(parent, coefficients, order)


__all__ = [
    "COEFFICIENT_COUNT",
    "FAMILY_ID",
    "SPATIAL_BASIS",
    "TEMPORAL_BASIS",
    "ZERO_COEFFICIENTS",
    "TemporalSwirlCandidate",
]
