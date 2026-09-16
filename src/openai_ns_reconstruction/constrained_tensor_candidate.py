"""Bounded tensor-polynomial corrections to :class:`CompactCandidate`.

The compact candidate already provides the axisymmetric chart and its smooth
compact bump.  ``TensorCandidate`` keeps that geometry and adds two bounded
polynomials in the fixed coordinates

``s = R**2 / 4``, ``v = Z**2 / 4`` and ``w = 2 * (t - .25)``.

The poloidal polynomial is added to the streamfunction factor ``q``.  The
swirl polynomial is added to the existing swirl multiplier.  Both corrections
are evaluated only through the compact chart, so the bump still supplies the
zero extension.  Spatial derivatives of the poloidal polynomial are included
explicitly in the Cartesian curl formula; this keeps the removable axis
singularity out of the direct velocity evaluation.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import json
from pathlib import Path
from typing import ClassVar

import numpy as np

from .constrained_candidate import CompactCandidate


TENSOR_BASIS = tuple(
    (i, j, k)
    for i in range(3)
    for j in range(3)
    for k in range(3)
)
COEFFICIENT_COUNT = len(TENSOR_BASIS)
ZERO_COEFFICIENTS = (0.0,) * COEFFICIENT_COUNT
FAMILY_ID = "compact_axisymmetric_tensor_v1"


def _coefficient_tuple(value, name: str) -> tuple[float, ...]:
    """Validate and freeze one 27 coefficient vector."""
    try:
        array = np.asarray(value, dtype=float)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(
            f"{name} must be a finite length-{COEFFICIENT_COUNT} sequence"
        ) from exc
    if array.ndim != 1 or array.shape[0] != COEFFICIENT_COUNT:
        raise ValueError(
            f"{name} must have length {COEFFICIENT_COUNT}"
        )
    if not np.all(np.isfinite(array)) or np.any((array < -1.0) | (array > 1.0)):
        raise ValueError(f"{name} entries must be finite and in [-1, 1]")
    return tuple(float(entry) for entry in array)


def _polynomial_and_spatial_derivatives(
    s: np.ndarray,
    v: np.ndarray,
    w: np.ndarray,
    coefficients: tuple[float, ...],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return ``P``, ``partial_s P`` and ``partial_v P`` for the tensor basis."""
    # The chart has already broadcast these arrays.  Explicit powers avoid a
    # negative exponent at the axis when differentiating the i=1 or j=1 rows.
    s_powers = (np.ones_like(s, dtype=float), s, s * s)
    v_powers = (np.ones_like(v, dtype=float), v, v * v)
    w_powers = (np.ones_like(w, dtype=float), w, w * w)
    value = np.zeros_like(s, dtype=float)
    derivative_s = np.zeros_like(s, dtype=float)
    derivative_v = np.zeros_like(s, dtype=float)
    for coefficient, (i, j, k) in zip(coefficients, TENSOR_BASIS):
        spatial_term = v_powers[j] * w_powers[k]
        value += coefficient * s_powers[i] * spatial_term
        if i:
            derivative_s += coefficient * i * s_powers[i - 1] * spatial_term
        if j:
            derivative_v += coefficient * j * s_powers[i] * v_powers[j - 1] * w_powers[k]
    return value, derivative_s, derivative_v


@dataclass(frozen=True)
class TensorCandidate(CompactCandidate):
    """Compact candidate with two bounded 3 by 3 by 3 tensor corrections."""

    poloidal_coefficients: tuple[float, ...] = ZERO_COEFFICIENTS
    swirl_coefficients: tuple[float, ...] = ZERO_COEFFICIENTS

    BASIS: ClassVar[tuple[tuple[int, int, int], ...]] = TENSOR_BASIS
    COEFFICIENT_COUNT: ClassVar[int] = COEFFICIENT_COUNT
    FAMILY_ID: ClassVar[str] = FAMILY_ID

    def __post_init__(self) -> None:
        # CompactCandidate validates all inherited geometric, time and
        # pressure parameters before the tensor-specific checks below.
        super().__post_init__()
        object.__setattr__(
            self,
            "poloidal_coefficients",
            _coefficient_tuple(self.poloidal_coefficients, "poloidal_coefficients"),
        )
        object.__setattr__(
            self,
            "swirl_coefficients",
            _coefficient_tuple(self.swirl_coefficients, "swirl_coefficients"),
        )

    def _corrected_chart(self, points, time):
        """Return the inherited chart plus corrected factors and derivatives."""
        (
            x,
            y,
            tau,
            lr,
            lz,
            R2,
            Z,
            bump,
            bump_R_over_R,
            bump_Z,
            base_q,
        ) = super()._chart(points, time)

        s = R2 / 4.0
        v = Z * Z / 4.0
        # The inherited chart stores tau=1-t, so this is exactly
        # w=2*(t-.25), including for broadcasted point/time arrays.
        w = 2.0 * (0.75 - tau)
        poloidal, poloidal_s, poloidal_v = _polynomial_and_spatial_derivatives(
            s, v, w, self.poloidal_coefficients
        )
        swirl, _, _ = _polynomial_and_spatial_derivatives(
            s, v, w, self.swirl_coefficients
        )

        q = base_q + poloidal
        q_s = self.radial_shape + (0.75 - tau) * self.poloidal_radial_time
        q_s = q_s + poloidal_s
        q_v = self.axial_shape + (0.75 - tau) * self.poloidal_axial_time
        q_v = q_v + poloidal_v

        base_swirl = (
            1.0
            + self.swirl_radial_shape * s
            + self.swirl_axial_shape * v
            + (0.75 - tau)
            * (
                self.swirl_radial_time * (s - 0.0025)
                + self.swirl_axial_time * (v - 0.0025)
            )
        )
        swirl_multiplier = base_swirl + swirl
        return (
            x,
            y,
            tau,
            lr,
            lz,
            R2,
            Z,
            bump,
            bump_R_over_R,
            bump_Z,
            q,
            q_s,
            q_v,
            swirl_multiplier,
        )

    def velocity(self, points, time):
        if (
            self.poloidal_coefficients == ZERO_COEFFICIENTS
            and self.swirl_coefficients == ZERO_COEFFICIENTS
        ):
            return super().velocity(points, time)
        (
            x,
            y,
            tau,
            lr,
            lz,
            R2,
            Z,
            bump,
            bump_R_over_R,
            bump_Z,
            q,
            q_s,
            q_v,
            swirl_multiplier,
        ) = self._corrected_chart(points, time)
        scale = self.amplitude * tau ** (-0.505)

        # For H=scale*Z*bump*q and A=(-y*H,x*H,0),
        # curl(A)=(-x*H_Z,-y*H_Z,2*H+R*H_R).  Since
        # q_R/R=q_s/2 and q_Z=q_v*Z/2, these expressions remain Cartesian
        # and regular at R=0.
        radial = -scale / lz * (
            bump * q + Z * (bump_Z * q + bump * q_v * Z / 2.0)
        )
        swirl = scale * self.swirl_ratio * bump / lr * swirl_multiplier
        axial = scale * Z * (
            2.0 * bump * q
            + R2 * (bump_R_over_R * q + bump * q_s / 2.0)
        )
        return np.stack(
            (x * radial - y * swirl, y * radial + x * swirl, axial), axis=-1
        )

    def vector_potential(self, points, time):
        """Axis-regular vector potential for the inherited LocalField bridge."""
        if (
            self.poloidal_coefficients == ZERO_COEFFICIENTS
            and self.swirl_coefficients == ZERO_COEFFICIENTS
        ):
            return super().vector_potential(points, time)
        (
            x,
            y,
            tau,
            _,
            _,
            _,
            Z,
            bump,
            _,
            _,
            q,
            _,
            _,
            _,
        ) = self._corrected_chart(points, time)
        H = self.amplitude * tau ** (-0.505) * Z * bump * q
        return np.stack((-y * H, x * H, np.zeros_like(H)), axis=-1)

    def save(self, path) -> None:
        """Save a versioned tensor candidate artifact."""
        Path(path).write_text(
            json.dumps(
                {
                    "schema_version": 1,
                    "status": "candidate",
                    "paper_exact": False,
                    "family": self.FAMILY_ID,
                    "parameters": asdict(self),
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )

    @classmethod
    def load(cls, path):
        """Load a tensor candidate and restore immutable coefficient tuples."""
        data = json.loads(Path(path).read_text(encoding="utf-8"))
        if data.get("schema_version") != 1 or data.get("family") != cls.FAMILY_ID:
            raise ValueError("unsupported tensor candidate artifact")
        parameters = data.get("parameters")
        if not isinstance(parameters, dict):
            raise ValueError("tensor candidate parameters must be an object")
        return cls(**parameters)
