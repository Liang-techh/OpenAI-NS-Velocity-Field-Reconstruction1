"""Finite Part I axial-moment repair on the saved pressure-bound core.

This module owns one deliberately narrow bridge.  It keeps the saved scalar
core swirl ``F`` and the pressure supplied by ``JoinedOuterPressure``.  The
core axial profile is multiplied by a flat radial cutoff and two normalized
positive C-infinity bumps are added so that the two axial moment identities

    integral U dR = 0,
    integral 2 R F U dR = 0

hold on the Chebyshev eta grid used by the saved degree-four core.  The first
coefficient is interpolated spectrally; the second is always defined by the
first axial identity, ``c2=-B-c1``.  This keeps the first identity exact in
the finite quadrature representation at off-grid eta values as well.

The object is a finite numerical adapter, not the Lei--Ren outer construction:
it does not implement the O.1--O.8 profile, all five moments, stress, cone,
energy, forcing, or recursive flatness.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property, lru_cache
from pathlib import Path
import math
import sys
from typing import Any, Callable

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from openai_ns_reconstruction.profiles import LeadingProfile  # noqa: E402
from openai_ns_reconstruction.quadrature import unit_rule  # noqa: E402

from lei_ren_part1_pressure_core import (  # noqa: E402
    JoinedOuterPressure,
    load_core,
)


SOURCE = "https://arxiv.org/html/2609.35406v1"
SOURCE_VERSION = "2609.35406v1"

R_CORE = 0.05
R_JOIN = 0.20
BUMP_SUPPORTS = ((0.055, 0.095), (0.115, 0.180))
DEFAULT_QUADRATURE_ORDER = 64


def _finite(value: Any, name: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be finite") from exc
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def _smooth_step(s: float) -> float:
    """A scalar C-infinity step, flat at both endpoints."""

    s = _finite(s, "s")
    if s <= 0.0:
        return 0.0
    if s >= 1.0:
        return 1.0
    left = math.exp(-1.0 / s)
    right = math.exp(-1.0 / (1.0 - s))
    return left / (left + right)


def _chi(radius: float) -> float:
    """One through ``R_CORE`` and flat to zero at ``R_JOIN``."""

    radius = _finite(radius, "R")
    if radius <= R_CORE:
        return 1.0
    if radius >= R_JOIN:
        return 0.0
    return 1.0 - _smooth_step((radius - R_CORE) / (R_JOIN - R_CORE))


def _raw_bump(radius: float, left: float, right: float) -> float:
    """Positive flat bump on the open interval ``(left,right)``."""

    radius = _finite(radius, "R")
    if radius <= left or radius >= right:
        return 0.0
    s = (radius - left) / (right - left)
    return math.exp(-1.0 / s - 1.0 / (1.0 - s))


def _split_edges(*extra: float) -> tuple[float, ...]:
    values = {0.0, R_CORE, R_JOIN}
    for left, right in BUMP_SUPPORTS:
        values.update((left, right))
    values.update(float(item) for item in extra)
    return tuple(sorted(values))


def _quadrature_segments(
    function: Callable[[float], float],
    left: float,
    right: float,
    *,
    order: int,
    edges: tuple[float, ...] = (),
) -> float:
    """Integrate a scalar function with Gauss rules split at flat joins."""

    left = _finite(left, "left")
    right = _finite(right, "right")
    if right < left:
        raise ValueError("right must be at least left")
    if right == left:
        return 0.0
    breaks = [left]
    breaks.extend(value for value in edges if left < value < right)
    breaks.append(right)
    nodes, weights = unit_rule(order)
    total = 0.0
    for a, b in zip(breaks[:-1], breaks[1:]):
        width = b - a
        points = a + width * nodes
        values = np.asarray([function(float(point)) for point in points], dtype=float)
        if values.shape != nodes.shape or not np.all(np.isfinite(values)):
            raise ValueError("quadrature integrand returned non-finite values")
        total += width * float(weights @ values)
    if not math.isfinite(total):
        raise ArithmeticError("quadrature returned a non-finite value")
    return float(total)


@dataclass(frozen=True)
class _RadialGrid:
    nodes: np.ndarray
    weights: np.ndarray
    chi: np.ndarray
    bump1: np.ndarray
    bump2: np.ndarray


class AxialMatchedProfile:
    """Saved core with a finite two-bump axial moment repair.

    Parameters are the tuple returned by
    :func:`lei_ren_part1_pressure_core.load_core`.  The public ``profile``
    property is a repository ``LeadingProfile`` with the original joined
    ``F`` and pressure and the repaired ``U``/radial averages.
    """

    def __init__(
        self,
        core: Any,
        joined: JoinedOuterPressure,
        *,
        quadrature_order: int = DEFAULT_QUADRATURE_ORDER,
    ) -> None:
        if not hasattr(core, "u_coefficients") or not hasattr(core, "grid"):
            raise TypeError("core must expose u_coefficients and grid")
        if not isinstance(joined, JoinedOuterPressure):
            raise TypeError("joined must be a JoinedOuterPressure")
        if isinstance(quadrature_order, bool) or not isinstance(quadrature_order, int):
            raise TypeError("quadrature_order must be an integer")
        if not 2 <= quadrature_order <= 2048:
            raise ValueError("quadrature_order must be in [2,2048]")
        if abs(float(joined.R_core) - R_CORE) > 1e-14 or abs(float(joined.R_join) - R_JOIN) > 1e-14:
            raise ValueError("axial repair requires the saved R_core=.05, R_join=.20 geometry")

        self.core = core
        self.joined = joined
        self.quadrature_order = int(quadrature_order)
        self.R_core = R_CORE
        self.R_join = R_JOIN
        self.bump_supports = BUMP_SUPPORTS
        coefficients = np.asarray(core.u_coefficients, dtype=float)
        if coefficients.ndim != 2 or coefficients.shape[1] != core.grid.nodes:
            raise ValueError("core.u_coefficients must be (degree+1, eta_nodes)")
        if not np.all(np.isfinite(coefficients)):
            raise ValueError("core.u_coefficients must be finite")
        self._u_coefficients = coefficients.copy()
        self._u_coefficients.setflags(write=False)
        derivative_coefficients = np.asarray(core.grid.differentiate(self._u_coefficients), dtype=float)
        derivative_coefficients.setflags(write=False)
        self._d_u_coefficients = derivative_coefficients
        self._degree = coefficients.shape[0] - 1

        self._radial = self._build_radial_grid()
        self._base_integral_coefficients = self._base_integrals_to_join()
        self._bump_normalizations = tuple(
            self._bump_normalization(index) for index in range(len(BUMP_SUPPORTS))
        )
        if not all(math.isfinite(value) and value > 0.0 for value in self._bump_normalizations):
            raise ArithmeticError("bump normalization is nonpositive or non-finite")
        self._bump1_grid = self._normalized_bump_values(0)
        self._bump2_grid = self._normalized_bump_values(1)

        # Precompute the joined F values only once.  HeatExterior's adaptive
        # quadrature is comparatively expensive and the same nodes are needed
        # for the full Chebyshev eta grid and both bump moments.
        self._F_grid = self._joined_F_grid()
        self._base_U_grid = self._base_U_values()
        self._B_grid = self._base_integral_coefficients @ self._u_coefficients
        self._B_derivative_grid = self.core.grid.differentiate(self._B_grid)
        self._C_grid = self._integral_over_grid(
            2.0 * self._radial.nodes[:, None] * self._F_grid * self._base_U_grid
        )
        self._w1_grid = self._integral_over_grid(
            2.0 * self._radial.nodes[:, None] * self._F_grid * self._bump1_grid[:, None]
        )
        self._w2_grid = self._integral_over_grid(
            2.0 * self._radial.nodes[:, None] * self._F_grid * self._bump2_grid[:, None]
        )
        denominator = self._w1_grid - self._w2_grid
        if not np.all(np.isfinite(denominator)) or np.any(denominator == 0.0):
            raise ArithmeticError("axial moment matrix is singular")
        self._denominator_grid = denominator
        self._c1_grid = (-self._C_grid + self._w2_grid * self._B_grid) / denominator
        self._c2_grid = -self._B_grid - self._c1_grid
        self._c1_derivative_grid = self.core.grid.differentiate(self._c1_grid)
        self._c2_derivative_grid = -self._B_derivative_grid - self._c1_derivative_grid

        for name, values in (
            ("B_grid", self._B_grid),
            ("C_grid", self._C_grid),
            ("w1_grid", self._w1_grid),
            ("w2_grid", self._w2_grid),
            ("c1_grid", self._c1_grid),
            ("c2_grid", self._c2_grid),
        ):
            values = np.asarray(values, dtype=float)
            if not np.all(np.isfinite(values)):
                raise ArithmeticError(f"{name} contains a non-finite value")
            values.setflags(write=False)
            setattr(self, "_" + name, values)
        self._profile_cache: LeadingProfile | None = None

    def _build_radial_grid(self) -> _RadialGrid:
        nodes, weights = unit_rule(self.quadrature_order)
        radial_nodes: list[float] = []
        radial_weights: list[float] = []
        edges = _split_edges()
        for left, right in zip(edges[:-1], edges[1:]):
            width = right - left
            radial_nodes.extend((left + width * nodes).tolist())
            radial_weights.extend((width * weights).tolist())
        radii = np.asarray(radial_nodes, dtype=float)
        weights_all = np.asarray(radial_weights, dtype=float)
        chi = np.asarray([_chi(float(radius)) for radius in radii], dtype=float)
        bump1 = np.asarray(
            [
                _raw_bump(float(radius), *BUMP_SUPPORTS[0])
                for radius in radii
            ],
            dtype=float,
        )
        bump2 = np.asarray(
            [
                _raw_bump(float(radius), *BUMP_SUPPORTS[1])
                for radius in radii
            ],
            dtype=float,
        )
        return _RadialGrid(radii, weights_all, chi, bump1, bump2)

    def _base_integrals_to_join(self) -> np.ndarray:
        degree = self._degree
        result = np.empty(degree + 1, dtype=float)
        for power in range(degree + 1):
            result[power] = float(
                self._radial.weights
                @ (self._radial.chi * self._radial.nodes**power)
            )
        return result

    def _bump_normalization(self, index: int) -> float:
        left, right = BUMP_SUPPORTS[index]
        return _quadrature_segments(
            lambda radius: _raw_bump(radius, left, right),
            left,
            right,
            order=self.quadrature_order,
            edges=(left, right),
        )

    def _normalized_bump_values(self, index: int) -> np.ndarray:
        left, right = BUMP_SUPPORTS[index]
        normalization = self._bump_normalizations[index]
        return np.asarray(
            [
                _raw_bump(float(radius), left, right) / normalization
                for radius in self._radial.nodes
            ],
            dtype=float,
        )

    def _joined_F_grid(self) -> np.ndarray:
        return np.asarray(
            [
                [float(self.joined.F(float(radius), float(eta))) for eta in self.core.grid.eta]
                for radius in self._radial.nodes
            ],
            dtype=float,
        )

    def _base_U_values(self) -> np.ndarray:
        values = np.asarray(
            self.core.U(self._radial.nodes[:, None], self.core.grid.eta[None, :]),
            dtype=float,
        )
        if values.shape != (self._radial.nodes.size, self.core.grid.nodes):
            raise ValueError("core.U returned an unexpected broadcast shape")
        return self._radial.chi[:, None] * values

    def _integral_over_grid(self, values: np.ndarray) -> np.ndarray:
        values = np.asarray(values, dtype=float)
        if values.shape != (self._radial.nodes.size, self.core.grid.nodes):
            raise ValueError("radial grid integrand has an unexpected shape")
        return self._radial.weights @ values

    def _validate_point(self, R: Any, Z: Any) -> tuple[float, float]:
        radius, eta = _finite(R, "R"), _finite(Z, "Z")
        if radius < 0.0 or abs(eta) > 1.0:
            raise ValueError("require R>=0 and |Z|<=1")
        return radius, eta

    def _interpolate(self, values: np.ndarray, eta: float) -> float:
        return float(self.core.grid.interpolate(values, eta))

    def _coefficients_at(self, eta: float) -> tuple[float, float, float, float, float, float]:
        B = self._interpolate(self._B_grid, eta)
        B_eta = self._interpolate(self._B_derivative_grid, eta)
        c1 = self._interpolate(self._c1_grid, eta)
        c1_eta = self._interpolate(self._c1_derivative_grid, eta)
        c2 = -B - c1
        c2_eta = -B_eta - c1_eta
        return B, B_eta, c1, c1_eta, c2, c2_eta

    def _normalized_bump(self, radius: float, index: int) -> float:
        left, right = BUMP_SUPPORTS[index]
        return _raw_bump(radius, left, right) / self._bump_normalizations[index]

    @lru_cache(maxsize=512)
    def _base_primitive_coefficients(self, radius: float) -> tuple[float, ...]:
        """Return ``integral_0^R chi(rho) rho**n d rho`` for every core row."""

        radius = _finite(radius, "R")
        if radius <= 0.0:
            return (0.0,) * (self._degree + 1)
        if radius <= R_CORE:
            return tuple(radius ** (power + 1) / (power + 1) for power in range(self._degree + 1))
        if radius >= R_JOIN:
            return tuple(float(value) for value in self._base_integral_coefficients)
        values = []
        for power in range(self._degree + 1):
            values.append(
                radius**0 * (
                    R_CORE ** (power + 1) / (power + 1)
                    + _quadrature_segments(
                        lambda point, p=power: _chi(point) * point**p,
                        R_CORE,
                        radius,
                        order=self.quadrature_order,
                        edges=_split_edges(radius),
                    )
                )
            )
        return tuple(float(value) for value in values)

    @lru_cache(maxsize=512)
    def _bump_primitive(self, radius: float, index: int) -> float:
        radius = _finite(radius, "R")
        left, right = BUMP_SUPPORTS[index]
        if radius <= left:
            return 0.0
        if radius >= right:
            return 1.0
        normalization = self._bump_normalizations[index]
        return _quadrature_segments(
            lambda point: _raw_bump(point, left, right) / normalization,
            left,
            radius,
            order=self.quadrature_order,
            edges=(left, right),
        )

    def F(self, R: Any, Z: Any) -> float:
        radius, eta = self._validate_point(R, Z)
        return float(self.joined.F(radius, eta))

    def E(self, R: Any, Z: Any) -> float:
        radius, eta = self._validate_point(R, Z)
        if radius == 0.0:
            return 0.0
        return math.sqrt(2.0 * radius) * self.F(radius, eta)

    def U(self, R: Any, Z: Any) -> float:
        radius, eta = self._validate_point(R, Z)
        if radius >= R_JOIN:
            return 0.0
        core_value = float(self.core.U(radius, eta))
        if radius <= R_CORE:
            return core_value
        B, _B_eta, c1, _c1_eta, c2, _c2_eta = self._coefficients_at(eta)
        del B
        return (
            _chi(radius) * core_value
            + c1 * self._normalized_bump(radius, 0)
            + c2 * self._normalized_bump(radius, 1)
        )

    def dU_deta(self, R: Any, Z: Any) -> float:
        radius, eta = self._validate_point(R, Z)
        if radius >= R_JOIN:
            return 0.0
        core_value = float(self.core.dU_deta(radius, eta))
        if radius <= R_CORE:
            return core_value
        _B, _B_eta, _c1, c1_eta, _c2, c2_eta = self._coefficients_at(eta)
        return (
            _chi(radius) * core_value
            + c1_eta * self._normalized_bump(radius, 0)
            + c2_eta * self._normalized_bump(radius, 1)
        )

    def average_U(self, R: Any, Z: Any) -> float:
        radius, eta = self._validate_point(R, Z)
        if radius >= R_JOIN:
            return 0.0
        if radius == 0.0:
            return float(self.core.U(0.0, eta))
        coefficients = np.asarray(
            self.core.grid.interpolate(self._u_coefficients, eta), dtype=float
        )
        primitive = float(np.dot(self._base_primitive_coefficients(radius), coefficients))
        _B, _B_eta, c1, _c1_eta, c2, _c2_eta = self._coefficients_at(eta)
        del _B
        primitive += c1 * self._bump_primitive(radius, 0)
        primitive += c2 * self._bump_primitive(radius, 1)
        return float(primitive / radius)

    def average_dU_deta(self, R: Any, Z: Any) -> float:
        radius, eta = self._validate_point(R, Z)
        if radius >= R_JOIN:
            return 0.0
        if radius == 0.0:
            return float(self.core.dU_deta(0.0, eta))
        coefficients = np.asarray(
            self.core.grid.interpolate(self._d_u_coefficients, eta), dtype=float
        )
        primitive = float(np.dot(self._base_primitive_coefficients(radius), coefficients))
        _B, _B_eta, _c1, c1_eta, _c2, c2_eta = self._coefficients_at(eta)
        primitive += c1_eta * self._bump_primitive(radius, 0)
        primitive += c2_eta * self._bump_primitive(radius, 1)
        return float(primitive / radius)

    def radial_flux_factor(self, R: Any, Z: Any, *, h: float | None = None) -> float:
        """Convenience wrapper for the regular incompressibility formula."""

        radius, eta = self._validate_point(R, Z)
        profile = self.profile
        exponent = self.core.reference.h if h is None else _finite(h, "h")
        return float(profile.radial_flux_factor(radius, eta, exponent, d=1.0 - eta**2,
                                               L=1.0 - 2.0 * exponent * eta**2))

    def V0(self, R: Any, Z: Any) -> float:
        radius, eta = self._validate_point(R, Z)
        return radius * self.radial_flux_factor(radius, eta)

    @property
    def profile(self) -> LeadingProfile:
        if self._profile_cache is None:
            self._profile_cache = LeadingProfile(
                E=self.E,
                U=self.U,
                dU_deta=self.dU_deta,
                Pi=self.joined.pressure,
                F=self.F,
                average_U=self.average_U,
                average_dU_deta=self.average_dU_deta,
                name="lei-ren-finite-axial-moment-match",
                paper_exact=False,
                provenance=(
                    "Finite degree-four saved core with flat radial cutoff and "
                    "two normalized positive axial bumps; source arXiv:2609.35406v1, "
                    "Section 2.5 moment identities."
                ),
            )
        return self._profile_cache

    def metadata(self) -> dict[str, Any]:
        denominator = np.asarray(self._denominator_grid, dtype=float)
        scale = max(float(np.max(np.abs(self._w1_grid))), float(np.max(np.abs(self._w2_grid))), np.finfo(float).tiny)
        nodal_mz = self._B_grid + self._c1_grid + self._c2_grid
        nodal_mixed = self._C_grid + self._c1_grid * self._w1_grid + self._c2_grid * self._w2_grid
        return {
            "source": SOURCE,
            "source_version": SOURCE_VERSION,
            "source_section": "Section 2.5, equations (2.21)-(2.22)",
            "status": "finite two-axial-moment repair",
            "core_degree": self._degree,
            "eta_nodes": int(self.core.grid.nodes),
            "quadrature_order": self.quadrature_order,
            "R_core": self.R_core,
            "R_join": self.R_join,
            "bump_supports": [list(item) for item in self.bump_supports],
            "bump_normalizations": list(self._bump_normalizations),
            "base_integral_coefficients": self._base_integral_coefficients.tolist(),
            "matrix_conditioning": {
                "w1_min": float(np.min(self._w1_grid)),
                "w1_max": float(np.max(self._w1_grid)),
                "w2_min": float(np.min(self._w2_grid)),
                "w2_max": float(np.max(self._w2_grid)),
                "denominator_min_abs": float(np.min(np.abs(denominator))),
                "denominator_max_abs": float(np.max(np.abs(denominator))),
                "relative_min_abs": float(np.min(np.abs(denominator)) / scale),
                "condition_proxy": float(np.max(np.abs(denominator)) / np.min(np.abs(denominator))),
            },
            "nodal_moment_residuals": {
                "max_abs_integral_U": float(np.max(np.abs(nodal_mz))),
                "max_abs_integral_2RFU": float(np.max(np.abs(nodal_mixed))),
            },
            "pressure_owner": "joined.pressure unchanged",
            "paper_exact": False,
            "all_five_moment_matching": False,
            "admissible_stress": False,
            "pde_validated": False,
            "scale_recursion_established": False,
            "scope": (
                "Finite radial axial-moment repair on the saved core; U is zero "
                "outside R_join and the original joined F/pressure are retained. "
                "No Lei-Ren O.1-O.8 outer construction, collar, cone, energy, "
                "forcing, all-five-moment closure, PDE, or recursion claim."
            ),
        }


__all__ = ["AxialMatchedProfile", "SOURCE", "SOURCE_VERSION", "load_core"]
