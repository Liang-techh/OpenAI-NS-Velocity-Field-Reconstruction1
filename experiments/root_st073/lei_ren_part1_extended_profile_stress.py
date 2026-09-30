"""Stress adapter for the actual ``SmoothExtendedAxial`` candidate.

This module binds the saved smooth axial repair to the profile-derived stress
formulas in :mod:`lei_ren_part1_profile_stress`.  Every moment is recomputed
from the same actual ``F`` and ``U`` supplied by the matched profile.  In
particular, no heat or collar target is substituted for the inner cumulative
moments, so a genuine closure defect remains visible in ``I`` and ``T``.

The axial cumulative moment ``M^z`` and its axial derivative use the matched
profile's ``R * average_U`` and ``R * average_dU_deta`` identities.  The other
moments use split Gauss--Legendre quadrature at the core, cross, bump, collar,
and geometric breakpoints.  ``F_eta``, ``P_Z``, and the noncore ``U_R`` use
declared finite differences or the analytic method available from the source
profile; their numerical errors are reported by the checks module.
"""

from __future__ import annotations

from functools import lru_cache
import math
from pathlib import Path
import sys
from typing import Any

import numpy as np
from numpy.polynomial.legendre import leggauss


ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
LOCAL = ROOT / "experiments" / "root_st073"
for path in (SRC, LOCAL):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from lei_ren_part1_profile_stress import ProfileStress  # noqa: E402
from lei_ren_part1_smooth_extended_axial import SmoothExtendedAxial  # noqa: E402


SOURCE = "https://arxiv.org/html/2609.35406v1"
SOURCE_VERSION = "2609.35406v1"
SOURCE_SECTIONS = "Sections 2.5 and 3.1--3.4, equations (2.21)--(2.22), (3.5), (3.9), (3.16)--(3.18)"
DEFAULT_RADIAL_QUADRATURE = 64
DEFAULT_Z_DIFFERENCE = 2.0e-4
DEFAULT_U_RADIAL_DIFFERENCE = 2.0e-5
R_GEOMETRIC_BREAK_COUNT = 20


def _finite(value: Any, name: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be finite") from exc
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def _validate_z(z: Any) -> float:
    value = _finite(z, "Z")
    if abs(value) > 1.0:
        raise ValueError("require |Z|<=1")
    return value


class ExtendedProfileStressAdapter:
    """Supply actual SmoothExtendedAxial data to ``ProfileStress``."""

    def __init__(
        self,
        matched: SmoothExtendedAxial | None = None,
        *,
        radial_quadrature_order: int = DEFAULT_RADIAL_QUADRATURE,
        z_difference: float = DEFAULT_Z_DIFFERENCE,
        u_radial_difference: float = DEFAULT_U_RADIAL_DIFFERENCE,
    ) -> None:
        if matched is None:
            matched = SmoothExtendedAxial()
        required = ("F", "U", "dU_deta", "average_U", "average_dU_deta", "swirl")
        if any(not callable(getattr(matched, name, None)) for name in required[:-1]):
            raise TypeError("matched must expose F, U, dU_deta, average_U and average_dU_deta")
        if not hasattr(matched, "swirl") or not hasattr(matched.swirl, "F_R"):
            raise TypeError("matched.swirl must expose F_R")
        radial_quadrature_order = int(radial_quadrature_order)
        if radial_quadrature_order not in (32, 64, 96, 128):
            raise ValueError("radial_quadrature_order must be 32, 64, 96, or 128")
        z_difference = _finite(z_difference, "z_difference")
        u_radial_difference = _finite(u_radial_difference, "u_radial_difference")
        if z_difference <= 0.0 or z_difference >= 0.01:
            raise ValueError("z_difference must lie in (0,.01)")
        if u_radial_difference <= 0.0 or u_radial_difference >= 0.01:
            raise ValueError("u_radial_difference must lie in (0,.01)")
        self.matched = matched
        self.swirl = matched.swirl
        self.core = matched.core
        self.collar = self.swirl.collar
        self.delta = float(getattr(matched, "delta", 2.0 * self.core.reference.h))
        self.radial_quadrature_order = radial_quadrature_order
        self.z_difference = z_difference
        self.u_radial_difference = u_radial_difference
        self.R_core = float(getattr(matched, "R_core", self.swirl.R_core))
        self.R_cutoff = float(getattr(matched, "R_cutoff", 0.02))
        self.R_anchor = float(self.swirl.R_anchor)
        self.R_a = float(self.collar.collar_inner_radius)
        self.R_b = float(self.collar.R_b)
        self._break_geometry_cache: dict[float, tuple[float, ...]] = {}

    @property
    def profile(self) -> ProfileStress:
        return ProfileStress(self, delta=self.delta)

    def _validate_R(self, R: Any) -> float:
        value = _finite(R, "R")
        if value < 0.0:
            raise ValueError("require R>=0")
        return value

    def F(self, R: Any, Z: Any) -> float:
        return float(self.swirl.F(self._validate_R(R), _validate_z(Z)))

    def U(self, R: Any, Z: Any) -> float:
        return float(self.matched.U(self._validate_R(R), _validate_z(Z)))

    def F_R(self, R: Any, Z: Any) -> float:
        return float(self.swirl.F_R(self._validate_R(R), _validate_z(Z)))

    def F_eta(self, R: Any, Z: Any) -> float:
        radius = self._validate_R(R)
        z = _validate_z(Z)
        if radius <= self.R_core:
            return float(self.core.F_eta(radius, z))
        if radius >= self.R_b:
            return float(self.collar.F_heat_Z(radius, z))
        if radius >= self.R_a:
            return float(self.collar.F_Z(radius, z))
        # The local anchor and smooth cross blend do not expose an analytic
        # Z derivative; use a fourth-order interior finite difference.
        step = self.z_difference
        if abs(z) < 2.0 * step:
            return float(
                (
                    self.F(radius, z + step)
                    - self.F(radius, z - step)
                )
                / (2.0 * step)
            )
        return float(
            (
                -self.F(radius, z + 2.0 * step)
                + 8.0 * self.F(radius, z + step)
                - 8.0 * self.F(radius, z - step)
                + self.F(radius, z - 2.0 * step)
            )
            / (12.0 * step)
        )

    def U_R(self, R: Any, Z: Any) -> float:
        radius = self._validate_R(R)
        z = _validate_z(Z)
        if radius <= self.R_core:
            return float(self.core.U_radial_derivative(radius, z))
        if radius >= self.R_a:
            return 0.0
        step = min(self.u_radial_difference, 0.2 * max(radius, self.R_core))
        left = max(0.0, radius - step)
        right = radius + step
        if left == 0.0:
            return float((self.U(right, z) - self.U(radius, z)) / (right - radius))
        return float((self.U(right, z) - self.U(left, z)) / (right - left))

    def U_eta(self, R: Any, Z: Any) -> float:
        return float(self.matched.dU_deta(self._validate_R(R), _validate_z(Z)))

    def P(self, R: Any, Z: Any) -> float:
        return float(self.swirl.pressure(self._validate_R(R), _validate_z(Z)))

    def P_Z(self, R: Any, Z: Any) -> float:
        radius = self._validate_R(R)
        z = _validate_z(Z)
        step = self.z_difference
        if abs(z) < 2.0 * step:
            return float(
                (self.P(radius, z + step) - self.P(radius, z - step))
                / (2.0 * step)
            )
        return float(
            (
                -self.P(radius, z + 2.0 * step)
                + 8.0 * self.P(radius, z + step)
                - 8.0 * self.P(radius, z - step)
                + self.P(radius, z - 2.0 * step)
            )
            / (12.0 * step)
        )

    def _breakpoints(self, Z: float, upper: float) -> tuple[float, ...]:
        z = _validate_z(Z)
        upper = self._validate_R(upper)
        cached = self._break_geometry_cache.get(z)
        if cached is None:
            cross = float(self.swirl._cross(z))
            values = [
                0.0,
                self.R_core,
                self.R_anchor,
                self.R_cutoff,
                cross,
                *[edge for pair in self.matched.bump_supports for edge in pair],
                self.R_a,
                self.R_b,
            ]
            if self.R_cutoff < self.R_a:
                values.extend(np.geomspace(self.R_cutoff, self.R_a, R_GEOMETRIC_BREAK_COUNT).tolist())
            cached = tuple(sorted({float(value) for value in values if value >= 0.0}))
            self._break_geometry_cache[z] = cached
        return tuple(value for value in cached if 0.0 <= value <= upper)

    def _radial_integrals(self, R: float, Z: float, order: int) -> dict[str, float]:
        radius = self._validate_R(R)
        z = _validate_z(Z)
        if radius == 0.0:
            return {"theta": 0.0, "theta_z": 0.0, "z_theta": 0.0, "p": 0.0}
        nodes, weights = leggauss(order)
        edges = self._breakpoints(z, radius)
        if not edges or edges[0] != 0.0:
            edges = (0.0, *edges)
        if edges[-1] < radius:
            edges = (*edges, radius)
        totals = {"theta": 0.0, "theta_z": 0.0, "z_theta": 0.0, "p": 0.0}
        for left, right in zip(edges[:-1], edges[1:]):
            if right <= left:
                continue
            points = left + (right - left) * (nodes + 1.0) / 2.0
            radial_weights = (right - left) * weights / 2.0
            f_values = np.asarray([self.F(float(point), z) for point in points], dtype=float)
            u_values = np.asarray([self.U(float(point), z) for point in points], dtype=float)
            totals["theta"] += float(radial_weights @ (2.0 * points * f_values))
            totals["theta_z"] += float(radial_weights @ (2.0 * points * f_values * u_values))
            totals["z_theta"] += float(radial_weights @ (u_values**2 - points * f_values**2))
            totals["p"] += float(radial_weights @ (f_values**2))
        if not all(math.isfinite(value) for value in totals.values()):
            raise ArithmeticError("actual profile moment quadrature returned non-finite value")
        return totals

    @lru_cache(maxsize=8192)
    def _moments_cached(self, R: float, Z: float, order: int) -> dict[str, float]:
        radius = self._validate_R(R)
        z = _validate_z(Z)
        totals = self._radial_integrals(radius, z, int(order))
        # Keep the matched primitive gauge exactly as used by the physical
        # field wrapper, rather than replacing it with an independently
        # quadratured U integral.
        totals["z"] = float(radius * self.matched.average_U(radius, z))
        return totals

    def moments(self, R: Any, Z: Any) -> dict[str, float]:
        return dict(self._moments_cached(self._validate_R(R), _validate_z(Z), self.radial_quadrature_order))

    def moments_Z(self, R: Any, Z: Any) -> dict[str, float]:
        radius = self._validate_R(R)
        z = _validate_z(Z)
        step = self.z_difference
        if abs(z) >= 2.0 * step:
            lower2 = self.moments(radius, z - 2.0 * step)
            lower = self.moments(radius, z - step)
            upper = self.moments(radius, z + step)
            upper2 = self.moments(radius, z + 2.0 * step)
            result = {
                key: float(
                    (
                        lower2[key]
                        - 8.0 * lower[key]
                        + 8.0 * upper[key]
                        - upper2[key]
                    )
                    / (12.0 * step)
                )
                for key in ("theta", "theta_z", "z_theta", "p")
            }
        else:
            lower = self.moments(radius, z - step)
            upper = self.moments(radius, z + step)
            result = {
                key: float((upper[key] - lower[key]) / (2.0 * step))
                for key in ("theta", "theta_z", "z_theta", "p")
            }
        # This is the prescribed matched identity for M^z_Z.
        result["z"] = float(radius * self.matched.average_dU_deta(radius, z))
        return result

    def stress(self, R: Any, Z: Any) -> tuple[float, float]:
        return self.profile.stress(self._validate_R(R), _validate_z(Z))

    def metadata(self) -> dict[str, Any]:
        matched_metadata_method = getattr(self.matched, "metadata", None)
        if callable(matched_metadata_method):
            matched_metadata = matched_metadata_method()
            if not isinstance(matched_metadata, dict):
                matched_metadata = {
                    "metadata_unavailable": True,
                    "metadata_type": type(matched_metadata).__name__,
                }
        else:
            matched_metadata = {"metadata_unavailable": True}
        return {
            "source": SOURCE,
            "source_version": SOURCE_VERSION,
            "source_sections": SOURCE_SECTIONS,
            "class": type(self).__name__,
            "matched_class": type(self.matched).__name__,
            "matched_metadata": matched_metadata,
            "matched_metadata_embedded": True,
            "pressure_source": matched_metadata.get(
                "profile_source", "unrecorded"
            ),
            "delta": self.delta,
            "R_core": self.R_core,
            "R_cutoff": self.R_cutoff,
            "R_anchor": self.R_anchor,
            "R_a": self.R_a,
            "R_b": self.R_b,
            "radial_quadrature_order": self.radial_quadrature_order,
            "z_difference": self.z_difference,
            "u_radial_difference": self.u_radial_difference,
            "moment_quadrature": (
                "actual F/U cumulative quadrature split at core, anchor, cutoff, "
                "cross, bump edges, geometric points, collar and heat edges"
            ),
            "Mz_gauge": "R * matched.average_U; Mz_Z=R * matched.average_dU_deta",
            "F_eta_gauge": "core.F_eta and collar F_Z where available; fourth-order Z FD in connection",
            "U_R_gauge": "core.U_radial_derivative in core; local centered FD elsewhere before R_a; exact zero for R>=R_a",
            "P_Z_gauge": "fourth-order centered Z finite difference",
            "heat_or_collar_targets_substituted": False,
            "whole_cone_validated": False,
            "full_five_moment_closure": False,
            "pde_validated": False,
        }


def load_extended_profile_stress(**kwargs: Any) -> ExtendedProfileStressAdapter:
    """Build the adapter on the saved SmoothExtendedAxial candidate."""

    return ExtendedProfileStressAdapter(**kwargs)


__all__ = [
    "DEFAULT_RADIAL_QUADRATURE",
    "ExtendedProfileStressAdapter",
    "SOURCE",
    "SOURCE_SECTIONS",
    "SOURCE_VERSION",
    "load_extended_profile_stress",
]
