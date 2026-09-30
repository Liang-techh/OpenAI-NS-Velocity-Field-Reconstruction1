"""Smooth axial moment repair for the shared-pressure extended profile.

This module binds a finite, Z-dependent axial profile to the actual
``ExtendedSwirl`` returned by
``lei_ren_part1_extended_pressure_core.load_extended_profile()``.  The core
``U`` is retained through ``R=.005`` and multiplied by a flat cutoff which
vanishes at ``R=.02``.  Four fixed compact C-infinity bumps, supported inside
the collar radius ``R_a``, supply the axial repair.

At each Chebyshev-Lobatto Z node, the two linear moments are solved with a
deterministic orthogonal projector.  The projector acts on the fixed physical
reference coefficient vector ``(1,-1,1,-1)``; its positive reference inner
product fixes the branch orientation without choosing signed SVD vectors.
The nonzero quadratic target is then reached along that oriented direction.
Moment arrays are interpolated on the same Chebyshev grid, and the two linear
constraints plus oriented quadratic root are solved algebraically at every
evaluation.  The last bump coefficient is recovered from the mass identity.
This is a continuous finite repair gauge; it deliberately does not claim the
minimum-L2 choice over the whole nullspace.

The public methods are scalar-callable and accept broadcast NumPy inputs:
``F``, ``U``, ``dU_deta``, ``average_U`` and ``average_dU_deta``.  The axial
primitive is defined to be exactly zero for ``R >= R_a``.  The returned
``LeadingProfile`` is an adapter for the repository's local field helpers;
this module does not certify a PDE, stress cone, global energy, or complete
five-moment Part I construction.

The Z-dependent moment arrays can be saved as a replay seed.  A seed records
the shared-pressure receipt hash, core polynomial/primitive inputs, geometry,
and all precomputed radial moments.  Replay rejects a changed pressure
receipt or incompatible core/geometry before constructing the adapter.

Source pinned to Z. Lei and X. Ren, arXiv:2609.35406v1, Section 2.5.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import hashlib
import json
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

from openai_ns_reconstruction.paper_core_series import ChebyshevEtaGrid  # noqa: E402
from openai_ns_reconstruction.profiles import LeadingProfile  # noqa: E402
from lei_ren_part1_exterior_targets import exterior_targets  # noqa: E402
from lei_ren_part1_extended_pressure_core import load_extended_profile  # noqa: E402


SOURCE = "https://arxiv.org/html/2609.35406v1"
SOURCE_VERSION = "2609.35406v1"
SOURCE_SECTION = "Section 2.5, equations (2.21)--(2.22)"

R_CORE = 0.005
R_CUTOFF = 0.02
DEFAULT_Z_NODES = 257
DEFAULT_RADIAL_QUADRATURE = 32
DEFAULT_TARGET_QUADRATURE = 64
LOG_BREAK_COUNT = 24
REFERENCE_COEFFICIENTS = np.asarray((1.0, -1.0, 1.0, -1.0), dtype=float)
BUMP_FRACTIONS = ((0.08, 0.20), (0.30, 0.42), (0.54, 0.66), (0.78, 0.90))
SEED_SCHEMA_VERSION = 1
PRESSURE_RECEIPT_NAME = "lei_ren_part1_extended_pressure_core.json"


def _finite(value: Any, name: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be finite") from exc
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def _positive_integer(value: Any, name: str) -> int:
    if isinstance(value, (bool, np.bool_)):
        raise TypeError(f"{name} must be a positive integer")
    try:
        result = int(value)
    except (TypeError, ValueError) as exc:
        raise TypeError(f"{name} must be a positive integer") from exc
    if result <= 0 or result != value:
        raise ValueError(f"{name} must be a positive integer")
    return result


def _flat_step(value: float) -> float:
    """C-infinity step from zero to one on ``[0,1]``."""

    value = _finite(value, "step coordinate")
    if value <= 0.0:
        return 0.0
    if value >= 1.0:
        return 1.0
    left = math.exp(-1.0 / value)
    right = math.exp(-1.0 / (1.0 - value))
    return left / (left + right)


def _flat_bump(radius: float, left: float, right: float) -> float:
    """Unit-height compact C-infinity bump on ``(left,right)``."""

    radius = _finite(radius, "R")
    left = _finite(left, "left")
    right = _finite(right, "right")
    if not left < right:
        raise ValueError("require left < right")
    if radius <= left or radius >= right:
        return 0.0
    coordinate = (radius - left) / (right - left)
    return float(math.exp(4.0 - 1.0 / coordinate - 1.0 / (1.0 - coordinate)))


def _broadcast_apply(function: Any, R: Any, Z: Any) -> Any:
    """Apply a scalar profile method to scalar or broadcast array inputs."""

    radii, eta = np.broadcast_arrays(
        np.asarray(R, dtype=float), np.asarray(Z, dtype=float)
    )
    if not np.all(np.isfinite(radii)) or not np.all(np.isfinite(eta)):
        raise ValueError("R and Z must contain only finite values")
    if radii.ndim == 0:
        return float(function(float(radii), float(eta)))
    vectorized = np.vectorize(function, otypes=[float])
    return vectorized(radii, eta)


def _pressure_receipt_path() -> Path:
    return LOCAL / PRESSURE_RECEIPT_NAME


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _pressure_receipt_identity() -> dict[str, str]:
    path = _pressure_receipt_path()
    if not path.is_file():
        raise FileNotFoundError(f"shared-pressure receipt is missing: {path}")
    return {"name": path.name, "sha256": _file_sha256(path)}


def _read_seed(path: Any) -> dict[str, Any]:
    seed_path = Path(path)
    if not seed_path.is_file():
        raise FileNotFoundError(f"smooth axial seed is missing: {seed_path}")
    try:
        data = json.loads(seed_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"invalid smooth axial seed: {seed_path}") from exc
    if not isinstance(data, dict) or data.get("schema_version") != SEED_SCHEMA_VERSION:
        raise ValueError("unsupported smooth axial seed schema")
    return data


@dataclass(frozen=True)
class _SliceData:
    """Radial moments used to solve one Chebyshev slice."""

    base_mass: float
    base_mixed: float
    base_energy: float
    base_bump_cross: np.ndarray
    bump_mass: np.ndarray
    bump_mixed: np.ndarray
    bump_gram: np.ndarray
    swirl_energy: float
    quadratic_target: float
    cross_radius: float


class SmoothExtendedAxial:
    """Smooth finite axial repair on the actual shared-pressure swirl.

    The coefficient functions are Chebyshev interpolants in ``Z``.  Their
    branch is fixed by projection of ``REFERENCE_COEFFICIENTS`` into the
    deterministic nullspace of the two moment rows.  The class intentionally
    exposes the profile adapter but does not assert that the resulting field
    is a complete Part I solution.
    """

    def __init__(
        self,
        swirl: Any | None = None,
        *,
        z_nodes: int = DEFAULT_Z_NODES,
        radial_quadrature_order: int = DEFAULT_RADIAL_QUADRATURE,
        target_quadrature_order: int = DEFAULT_TARGET_QUADRATURE,
        _seed_data: dict[str, Any] | None = None,
    ) -> None:
        if swirl is None:
            swirl = load_extended_profile()
        for name in ("F", "core", "collar", "_cross"):
            if not hasattr(swirl, name):
                raise TypeError("swirl must be the ExtendedSwirl returned by load_extended_profile()")
        self.swirl = swirl
        self.core = swirl.core
        self.R_core = R_CORE
        self.R_anchor = float(swirl.R_anchor)
        self.R_cutoff = R_CUTOFF
        self.Rmax = float(swirl.collar.collar_inner_radius)
        self.bump_supports = tuple(
            (self.Rmax * left, self.Rmax * right)
            for left, right in BUMP_FRACTIONS
        )
        self.z_nodes = _positive_integer(z_nodes, "z_nodes")
        if self.z_nodes < 3:
            raise ValueError("z_nodes must be at least 3")
        self.radial_quadrature_order = _positive_integer(
            radial_quadrature_order, "radial_quadrature_order"
        )
        self.target_quadrature_order = _positive_integer(
            target_quadrature_order, "target_quadrature_order"
        )
        self.z_grid = ChebyshevEtaGrid(self.z_nodes)
        self.reference_coefficients = REFERENCE_COEFFICIENTS.copy()
        self.reference_coefficients.setflags(write=False)

        coefficients = np.asarray(self.core.u_coefficients, dtype=float)
        if coefficients.ndim != 2 or coefficients.shape[1] != self.core.grid.nodes:
            raise ValueError("core.u_coefficients has an unexpected shape")
        if not np.all(np.isfinite(coefficients)):
            raise ValueError("core.u_coefficients must be finite")
        self._core_u_coefficients = coefficients.copy()
        self._core_u_eta_coefficients = np.asarray(
            self.core.grid.differentiate(coefficients), dtype=float
        )
        self._degree = coefficients.shape[0] - 1
        self._base_primitive_full = self._base_primitive_coefficients(self.R_cutoff)
        self._base_primitive_full.setflags(write=False)
        self._base_primitive_eta_full = self._base_primitive_eta_coefficients(
            self.R_cutoff
        )
        self._base_primitive_eta_full.setflags(write=False)

        self._bump_mass = self._bump_mass_integrals()
        self._bump_mass.setflags(write=False)

        if _seed_data is None:
            rows = [self._compute_slice(float(z)) for z in self.z_grid.eta]
        else:
            self._validate_seed(_seed_data)
            rows = self._rows_from_seed(_seed_data)
        self._slice_data = rows
        self._base_mass_grid = np.asarray([row.base_mass for row in rows], dtype=float)
        self._base_mixed_grid = np.asarray([row.base_mixed for row in rows], dtype=float)
        self._base_energy_grid = np.asarray([row.base_energy for row in rows], dtype=float)
        self._base_bump_cross_grid = np.asarray(
            [row.base_bump_cross for row in rows], dtype=float
        ).T
        self._bump_mixed_grid = np.asarray(
            [row.bump_mixed for row in rows], dtype=float
        ).T
        self._bump_gram_grid = np.asarray(
            [row.bump_gram for row in rows], dtype=float
        )
        self._swirl_energy_grid = np.asarray(
            [row.swirl_energy for row in rows], dtype=float
        )
        self._quadratic_target_grid = np.asarray(
            [row.quadratic_target for row in rows], dtype=float
        )
        self._cross_radius_grid = np.asarray(
            [row.cross_radius for row in rows], dtype=float
        )
        # The bump Gram matrix is independent of Z.  Recompute it once on a
        # high-order fixed-support rule instead of interpolating tiny changes
        # caused only by the per-slice geometry splitting.
        self._bump_gram = self._bump_gram_integrals()
        self._bump_gram.setflags(write=False)
        self._base_mixed_eta_grid = self.z_grid.differentiate(self._base_mixed_grid)
        self._base_energy_eta_grid = self.z_grid.differentiate(self._base_energy_grid)
        self._bump_mixed_eta_grid = self.z_grid.differentiate(self._bump_mixed_grid)
        self._swirl_energy_eta_grid = self.z_grid.differentiate(self._swirl_energy_grid)
        self._quadratic_target_eta_grid = self.z_grid.differentiate(
            self._quadratic_target_grid
        )
        self._moment_row_scales = np.asarray(
            (
                np.max(np.abs(self._bump_mass)),
                np.max(np.abs(self._bump_mixed_grid)),
            ),
            dtype=float,
        )

        solved = [
            self._solve_interpolated(float(z), with_derivative=True)
            for z in self.z_grid.eta
        ]
        self._coefficient_grid = np.asarray(
            [item["coefficients"] for item in solved], dtype=float
        ).T
        self._coefficient_eta_grid = np.asarray(
            self.z_grid.differentiate(self._coefficient_grid), dtype=float
        )
        self._direction_grid = np.asarray(
            [item["direction"] for item in solved], dtype=float
        ).T
        self._orientation_dot_grid = np.asarray(
            [item["orientation_dot"] for item in solved], dtype=float
        )
        self._conditioning = {
            "scaled_linear_condition_min": float(
                min(item["linear_condition"] for item in solved)
            ),
            "scaled_linear_condition_max": float(
                max(item["linear_condition"] for item in solved)
            ),
            "direction_norm_min": float(
                min(item["direction_norm"] for item in solved)
            ),
            "direction_norm_max": float(
                max(item["direction_norm"] for item in solved)
            ),
            "quadratic_discriminant_min": float(
                min(item["discriminant"] for item in solved)
            ),
            "quadratic_parameter_min": float(min(item["parameter"] for item in solved)),
            "quadratic_parameter_max": float(max(item["parameter"] for item in solved)),
            "moment_row_scale_mass": float(self._moment_row_scales[0]),
            "moment_row_scale_mixed": float(self._moment_row_scales[1]),
        }
        for value in (
            self._base_mass_grid,
            self._base_mixed_grid,
            self._base_energy_grid,
            self._bump_mixed_grid,
            self._bump_gram_grid,
            self._swirl_energy_grid,
            self._quadratic_target_grid,
            self._cross_radius_grid,
            self._base_mixed_eta_grid,
            self._base_energy_eta_grid,
            self._bump_mixed_eta_grid,
            self._swirl_energy_eta_grid,
            self._quadratic_target_eta_grid,
            self._coefficient_grid,
            self._coefficient_eta_grid,
            self._direction_grid,
            self._orientation_dot_grid,
        ):
            if not np.all(np.isfinite(value)):
                raise ArithmeticError("smooth axial precomputation is non-finite")
            value.setflags(write=False)
        self._profile_cache: LeadingProfile | None = None

    @staticmethod
    def _seed_array(
        seed: dict[str, Any], key: str, shape: tuple[int, ...]
    ) -> np.ndarray:
        arrays = seed.get("moments")
        if not isinstance(arrays, dict) or key not in arrays:
            raise ValueError(f"smooth axial seed is missing moments.{key}")
        try:
            result = np.asarray(arrays[key], dtype=float)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"smooth axial seed moments.{key} is not numeric") from exc
        if result.shape != shape or not np.all(np.isfinite(result)):
            raise ValueError(f"smooth axial seed moments.{key} has the wrong shape or non-finite values")
        return result

    @staticmethod
    def _require_close(
        actual: np.ndarray, expected: Any, name: str, *, tolerance: float = 2.0e-13
    ) -> None:
        candidate = np.asarray(expected, dtype=float)
        if actual.shape != candidate.shape or not np.allclose(
            actual, candidate, rtol=0.0, atol=tolerance
        ):
            raise ValueError(f"smooth axial seed {name} does not match the active profile")

    def _validate_seed(self, seed: dict[str, Any]) -> None:
        """Reject a seed produced for another pressure/core/geometry receipt."""

        expected_receipt = _pressure_receipt_identity()
        if seed.get("pressure_receipt") != expected_receipt:
            raise ValueError(
                "smooth axial seed pressure receipt hash does not match "
                f"{expected_receipt['name']}"
            )
        if int(seed.get("z_nodes", -1)) != self.z_nodes:
            raise ValueError("smooth axial seed z_nodes does not match the requested grid")
        if int(seed.get("radial_quadrature_order", -1)) != self.radial_quadrature_order:
            raise ValueError("smooth axial seed radial quadrature order does not match")
        if int(seed.get("target_quadrature_order", -1)) != self.target_quadrature_order:
            raise ValueError("smooth axial seed target quadrature order does not match")
        geometry = seed.get("geometry")
        if not isinstance(geometry, dict):
            raise ValueError("smooth axial seed geometry metadata is missing")
        for name, actual in (
            ("R_core", self.R_core),
            ("R_anchor", self.R_anchor),
            ("R_cutoff", self.R_cutoff),
            ("Rmax_Ra", self.Rmax),
        ):
            if abs(float(geometry.get(name, float("nan"))) - actual) > 2.0e-13:
                raise ValueError(f"smooth axial seed geometry {name} does not match")
        self._require_close(
            np.asarray(self.bump_supports, dtype=float),
            geometry.get("bump_supports"),
            "bump_supports",
        )
        self._require_close(self.z_grid.eta, seed.get("z_grid"), "z_grid")
        core_inputs = seed.get("core_inputs")
        if not isinstance(core_inputs, dict):
            raise ValueError("smooth axial seed core inputs are missing")
        self._require_close(
            self._core_u_coefficients,
            core_inputs.get("u_coefficients"),
            "core u_coefficients",
        )
        self._require_close(
            self._core_u_eta_coefficients,
            core_inputs.get("u_eta_coefficients"),
            "core u_eta_coefficients",
        )
        self._require_close(
            self._base_primitive_full,
            core_inputs.get("base_primitive_full"),
            "base primitive",
        )
        self._require_close(
            self._base_primitive_eta_full,
            core_inputs.get("base_primitive_eta_full"),
            "base eta primitive",
        )
        n = self.z_nodes
        for key, shape in (
            ("base_mass", (n,)),
            ("base_mixed", (n,)),
            ("base_energy", (n,)),
            ("base_bump_cross", (4, n)),
            ("bump_mass", (4, n)),
            ("bump_mixed", (4, n)),
            ("bump_gram", (n, 4, 4)),
            ("swirl_energy", (n,)),
            ("quadratic_target", (n,)),
            ("cross_radius", (n,)),
        ):
            self._seed_array(seed, key, shape)

    def _rows_from_seed(self, seed: dict[str, Any]) -> list[_SliceData]:
        n = self.z_nodes
        arrays = {
            key: self._seed_array(seed, key, shape)
            for key, shape in (
                ("base_mass", (n,)),
                ("base_mixed", (n,)),
                ("base_energy", (n,)),
                ("base_bump_cross", (4, n)),
                ("bump_mass", (4, n)),
                ("bump_mixed", (4, n)),
                ("bump_gram", (n, 4, 4)),
                ("swirl_energy", (n,)),
                ("quadratic_target", (n,)),
                ("cross_radius", (n,)),
            )
        }
        return [
            _SliceData(
                base_mass=float(arrays["base_mass"][index]),
                base_mixed=float(arrays["base_mixed"][index]),
                base_energy=float(arrays["base_energy"][index]),
                base_bump_cross=np.asarray(
                    arrays["base_bump_cross"][:, index], dtype=float
                ),
                bump_mass=np.asarray(arrays["bump_mass"][:, index], dtype=float),
                bump_mixed=np.asarray(arrays["bump_mixed"][:, index], dtype=float),
                bump_gram=np.asarray(arrays["bump_gram"][index], dtype=float),
                swirl_energy=float(arrays["swirl_energy"][index]),
                quadratic_target=float(arrays["quadratic_target"][index]),
                cross_radius=float(arrays["cross_radius"][index]),
            )
            for index in range(n)
        ]

    def _base_chi(self, radius: float) -> float:
        if radius <= self.R_core:
            return 1.0
        if radius >= self.R_cutoff:
            return 0.0
        return 1.0 - _flat_step(
            (radius - self.R_core) / (self.R_cutoff - self.R_core)
        )

    def _base_U_scalar(self, radius: float, eta: float) -> float:
        if radius <= self.R_core:
            return float(self.core.U(radius, eta))
        if radius >= self.R_cutoff:
            return 0.0
        return float(self._base_chi(radius) * self.core.U(radius, eta))

    def _base_U_eta_scalar(self, radius: float, eta: float) -> float:
        if radius <= self.R_core:
            return float(self.core.dU_deta(radius, eta))
        if radius >= self.R_cutoff:
            return 0.0
        return float(self._base_chi(radius) * self.core.dU_deta(radius, eta))

    def _breakpoints(self, eta: float) -> tuple[float, ...]:
        cross = float(self.swirl._cross(float(eta)))
        geometric = np.geomspace(self.R_cutoff, self.Rmax, LOG_BREAK_COUNT)
        values = [self.R_core, self.R_cutoff, cross]
        values.extend(edge for pair in self.bump_supports for edge in pair)
        values.extend(float(value) for value in geometric[1:-1])
        return tuple(sorted({value for value in values if 0.0 < value < self.Rmax}))

    def _radial_rule(
        self, eta: float, order: int
    ) -> tuple[np.ndarray, np.ndarray]:
        nodes, weights = leggauss(order)
        edges = (0.0, *self._breakpoints(eta), self.Rmax)
        radii: list[float] = []
        radial_weights: list[float] = []
        for left, right in zip(edges[:-1], edges[1:]):
            width = right - left
            radii.extend((left + width * (nodes + 1.0) / 2.0).tolist())
            radial_weights.extend((width * weights / 2.0).tolist())
        return np.asarray(radii, dtype=float), np.asarray(radial_weights, dtype=float)

    def _bump_values(self, radii: np.ndarray) -> np.ndarray:
        return np.asarray(
            [
                [
                    _flat_bump(float(radius), left, right)
                    for left, right in self.bump_supports
                ]
                for radius in radii
            ],
            dtype=float,
        )

    def _compute_slice(self, eta: float) -> _SliceData:
        radii, weights = self._radial_rule(eta, self.radial_quadrature_order)
        f_values = np.asarray(
            [float(self.swirl.F(float(radius), eta)) for radius in radii],
            dtype=float,
        )
        base_values = np.asarray(
            [self._base_U_scalar(float(radius), eta) for radius in radii],
            dtype=float,
        )
        bumps = self._bump_values(radii)
        base_mass = float(weights @ base_values)
        base_mixed = float(weights @ (2.0 * radii * f_values * base_values))
        base_energy = float(weights @ (base_values * base_values))
        base_bump_cross = bumps.T @ (weights * base_values)
        bump_mass = weights @ bumps
        bump_mixed = (weights * 2.0 * radii * f_values) @ bumps
        bump_gram = bumps.T @ (weights[:, None] * bumps)
        swirl_energy = float(weights @ (radii * f_values * f_values))
        target = exterior_targets(
            self.swirl.collar.F,
            self.swirl.collar.heat,
            self.Rmax,
            self.swirl.collar.R_b,
            eta,
            n=self.target_quadrature_order,
        )["quadratic_target"]
        values = (
            base_mass,
            base_mixed,
            base_energy,
            base_bump_cross,
            bump_mass,
            bump_mixed,
            bump_gram,
            swirl_energy,
            float(target),
            float(self.swirl._cross(eta)),
        )
        if not all(np.all(np.isfinite(value)) for value in values):
            raise ArithmeticError("radial moment precomputation returned non-finite values")
        return _SliceData(*values)

    def _solve_slice(
        self, row: _SliceData, *, base_mass_override: float | None = None
    ) -> dict[str, Any]:
        """Solve the two linear moments and oriented nonzero quadratic root."""

        base_mass = row.base_mass if base_mass_override is None else float(base_mass_override)
        linear = np.vstack((self._bump_mass, row.bump_mixed))
        rhs = -np.asarray((base_mass, row.base_mixed), dtype=float)
        row_scale = np.maximum(np.max(np.abs(linear), axis=1), np.finfo(float).tiny)
        scaled = linear / row_scale[:, None]
        scaled_rhs = rhs / row_scale
        gram_linear = scaled @ scaled.T
        particular = scaled.T @ np.linalg.solve(gram_linear, scaled_rhs)
        projector = np.eye(4) - scaled.T @ np.linalg.solve(gram_linear, scaled)
        direction = projector @ self.reference_coefficients
        orientation_dot = float(self.reference_coefficients @ direction)
        if orientation_dot < 0.0:
            direction = -direction
            orientation_dot = -orientation_dot
        direction_norm = math.sqrt(max(0.0, float(direction @ row.bump_gram @ direction)))
        if not math.isfinite(direction_norm) or direction_norm <= 1.0e-15:
            raise ArithmeticError("fixed reference projection has a degenerate direction")

        # q(p+t*d) is evaluated from the same radial rule that supplied the
        # linear moments.  ``base_bump_cross`` is integral(base_U*bump_i).
        q0 = float(
            row.base_energy
            + 2.0 * (particular @ row.base_bump_cross)
            + particular @ row.bump_gram @ particular
            - row.swirl_energy
        )
        cross = float(
            direction
            @ (row.base_bump_cross + row.bump_gram @ particular)
        )
        level = float(row.quadratic_target - q0)
        discriminant = float(cross * cross + direction_norm**2 * level)
        tolerance = 1.0e-12 * max(1.0, abs(cross * cross), abs(direction_norm**2 * level))
        if discriminant < -tolerance:
            raise ArithmeticError(
                "quadratic target is unreachable on the oriented reference branch"
            )
        discriminant = max(0.0, discriminant)
        parameter = float(
            (-cross + math.sqrt(discriminant)) / (direction_norm**2)
        )
        coefficients = particular + parameter * direction

        # Recover the final coefficient from the mass identity.  This is
        # intentionally explicit: the same recovery is used by off-grid
        # evaluation, so the terminal primitive does not acquire a numerical
        # jump when it is set to zero at R=Rmax.
        mass_defect = float(base_mass + self._bump_mass @ coefficients)
        coefficients[-1] -= mass_defect / self._bump_mass[-1]
        if not np.all(np.isfinite(coefficients)):
            raise ArithmeticError("oriented quadratic coefficients are non-finite")
        singular_values = np.linalg.svd(scaled, compute_uv=False)
        linear_condition = float(singular_values[0] / singular_values[-1])
        return {
            "coefficients": coefficients,
            "particular": particular,
            "direction": direction,
            "orientation_dot": orientation_dot,
            "direction_norm": direction_norm,
            "linear_condition": linear_condition,
            "discriminant": discriminant,
            "parameter": parameter,
            "q0": q0,
            "mass_defect_before_recovery": mass_defect,
        }

    @lru_cache(maxsize=256)
    def _base_primitive_coefficients(self, radius: float) -> np.ndarray:
        """Integrate ``chi(R) R**n`` for the retained core polynomial."""

        radius = _finite(radius, "R")
        if radius <= 0.0:
            return np.zeros(self._degree + 1, dtype=float)
        upper = min(radius, self.R_cutoff)
        result = np.empty(self._degree + 1, dtype=float)
        if upper <= self.R_core:
            for power in range(self._degree + 1):
                result[power] = upper ** (power + 1) / float(power + 1)
            return result
        nodes, weights = leggauss(96)
        result[:] = self.R_core ** (np.arange(self._degree + 1) + 1) / (
            np.arange(self._degree + 1) + 1.0
        )
        left, right = self.R_core, upper
        points = left + (right - left) * (nodes + 1.0) / 2.0
        radial_weights = (right - left) * weights / 2.0
        for power in range(self._degree + 1):
            result[power] += float(
                radial_weights
                @ np.asarray(
                    [self._base_chi(float(point)) * point**power for point in points],
                    dtype=float,
                )
            )
        return result

    def _base_primitive_eta_coefficients(self, radius: float) -> np.ndarray:
        """Same primitive weights for the core's analytic eta derivative."""

        return self._base_primitive_coefficients(radius)

    def _bump_mass_integrals(self) -> np.ndarray:
        nodes, weights = leggauss(128)
        result = []
        for left, right in self.bump_supports:
            points = left + (right - left) * (nodes + 1.0) / 2.0
            radial_weights = (right - left) * weights / 2.0
            result.append(
                float(
                    radial_weights
                    @ np.asarray(
                        [_flat_bump(float(point), left, right) for point in points],
                        dtype=float,
                    )
                )
            )
        return np.asarray(result, dtype=float)

    def _bump_gram_integrals(self) -> np.ndarray:
        """Return the fixed-support Gram matrix ``integral bump_i bump_j``."""

        nodes, weights = leggauss(160)
        gram = np.zeros((len(self.bump_supports), len(self.bump_supports)), dtype=float)
        for left, right in self.bump_supports:
            points = left + (right - left) * (nodes + 1.0) / 2.0
            radial_weights = (right - left) * weights / 2.0
            values = np.asarray(
                [
                    [
                        _flat_bump(float(point), support_left, support_right)
                        for support_left, support_right in self.bump_supports
                    ]
                    for point in points
                ],
                dtype=float,
            )
            # The supports are disjoint, so only the active interval contributes
            # to this local rule.  Adding the block keeps the implementation
            # explicit and avoids relying on a dense global radial mesh.
            gram += values.T @ (radial_weights[:, None] * values)
        return gram

    @staticmethod
    def _validate_point(radius: Any, eta: Any) -> tuple[float, float]:
        radius = _finite(radius, "R")
        eta = _finite(eta, "Z")
        if radius < 0.0 or abs(eta) > 1.0:
            raise ValueError("require R>=0 and |Z|<=1")
        return radius, eta

    def _interpolate(self, values: np.ndarray, eta: float) -> float:
        return float(np.asarray(self.z_grid.interpolate(values, eta), dtype=float))

    def _interpolate_rows(self, values: np.ndarray, eta: float) -> np.ndarray:
        result = np.asarray(self.z_grid.interpolate(values, eta), dtype=float)
        return result.reshape(result.shape[:-1] if result.ndim > 1 else result.shape)

    @lru_cache(maxsize=8192)
    def _core_coefficients(self, eta: float, *, eta_derivative: bool = False) -> np.ndarray:
        source = (
            self._core_u_eta_coefficients
            if eta_derivative
            else self._core_u_coefficients
        )
        result = np.asarray(
            self.core.grid.interpolate(source, eta), dtype=float
        ).reshape(-1)
        result.setflags(write=False)
        return result

    def _base_mass_at(self, eta: float) -> float:
        return float(self._base_primitive_full @ self._core_coefficients(eta))

    def _base_mass_eta_at(self, eta: float) -> float:
        return float(
            self._base_primitive_eta_full
            @ self._core_coefficients(eta, eta_derivative=True)
        )

    def _interpolated_moments(
        self, eta: float, *, with_derivative: bool
    ) -> dict[str, Any]:
        """Evaluate the same Chebyshev moment interpolants used by the gauge."""

        base_mass = self._base_mass_at(eta)
        base_mass_eta = self._base_mass_eta_at(eta) if with_derivative else 0.0
        base_mixed = self._interpolate(self._base_mixed_grid, eta)
        bump_mixed = np.asarray(
            self.z_grid.interpolate(self._bump_mixed_grid, eta), dtype=float
        ).reshape(4)
        base_energy = self._interpolate(self._base_energy_grid, eta)
        swirl_energy = self._interpolate(self._swirl_energy_grid, eta)
        target = self._interpolate(self._quadratic_target_grid, eta)
        if not with_derivative:
            return {
                "base_mass": base_mass,
                "base_mixed": base_mixed,
                "base_energy": base_energy,
                "bump_mixed": bump_mixed,
                "swirl_energy": swirl_energy,
                "target": target,
            }
        return {
            "base_mass": base_mass,
            "base_mass_eta": base_mass_eta,
            "base_mixed": base_mixed,
            "base_mixed_eta": self._interpolate(self._base_mixed_eta_grid, eta),
            "base_energy": base_energy,
            "base_energy_eta": self._interpolate(self._base_energy_eta_grid, eta),
            "bump_mixed": bump_mixed,
            "bump_mixed_eta": np.asarray(
                self.z_grid.interpolate(self._bump_mixed_eta_grid, eta), dtype=float
            ).reshape(4),
            "swirl_energy": swirl_energy,
            "swirl_energy_eta": self._interpolate(self._swirl_energy_eta_grid, eta),
            "target": target,
            "target_eta": self._interpolate(self._quadratic_target_eta_grid, eta),
        }

    def _solve_interpolated(
        self, eta: float, *, with_derivative: bool
    ) -> dict[str, Any]:
        """Solve the algebraic moment system on interpolated Z-moment data.

        The row scaling is constant in Z, so the projector and its analytic
        derivative retain the exact nullspace of the unscaled constraints.
        """

        moments = self._interpolated_moments(eta, with_derivative=with_derivative)
        linear = np.vstack((self._bump_mass, moments["bump_mixed"]))
        rhs = -np.asarray((moments["base_mass"], moments["base_mixed"]), dtype=float)
        scales = self._moment_row_scales
        scaled = linear / scales[:, None]
        scaled_rhs = rhs / scales
        gram_linear = scaled @ scaled.T
        y = np.linalg.solve(gram_linear, scaled_rhs)
        particular = scaled.T @ y
        projector_factor = np.linalg.solve(gram_linear, scaled)
        projector = np.eye(4) - scaled.T @ projector_factor
        direction = projector @ self.reference_coefficients
        orientation_dot = float(self.reference_coefficients @ direction)
        if orientation_dot < 0.0:
            direction = -direction
            orientation_dot = -orientation_dot

        if with_derivative:
            linear_eta = np.vstack(
                (
                    np.zeros(4, dtype=float),
                    moments["bump_mixed_eta"],
                )
            )
            rhs_eta = -np.asarray(
                (moments["base_mass_eta"], moments["base_mixed_eta"]), dtype=float
            )
            scaled_eta = linear_eta / scales[:, None]
            scaled_rhs_eta = rhs_eta / scales
            gram_linear_eta = scaled_eta @ scaled.T + scaled @ scaled_eta.T
            y_eta = np.linalg.solve(
                gram_linear, scaled_rhs_eta - gram_linear_eta @ y
            )
            particular_eta = scaled_eta.T @ y + scaled.T @ y_eta
            projector_factor_eta = np.linalg.solve(
                gram_linear,
                scaled_eta - gram_linear_eta @ projector_factor,
            )
            projector_eta = -(
                scaled_eta.T @ projector_factor
                + scaled.T @ projector_factor_eta
            )
            direction_eta = projector_eta @ self.reference_coefficients
            if orientation_dot < 0.0:  # defensive; orientation is fixed above
                direction_eta = -direction_eta
        else:
            particular_eta = np.zeros(4, dtype=float)
            direction_eta = np.zeros(4, dtype=float)

        gram = self._bump_gram
        direction_norm_sq = float(direction @ gram @ direction)
        if not math.isfinite(direction_norm_sq) or direction_norm_sq <= 1.0e-15:
            raise ArithmeticError("fixed reference projection has a degenerate direction")
        direction_norm = math.sqrt(direction_norm_sq)
        p_gram_d = float(particular @ gram @ direction)
        q0 = float(
            moments["base_energy"]
            - moments["swirl_energy"]
            + particular @ gram @ particular
        )
        cross = p_gram_d
        level = float(moments["target"] - q0)
        discriminant = float(cross * cross + direction_norm_sq * level)
        tolerance = 1.0e-12 * max(
            1.0, abs(cross * cross), abs(direction_norm_sq * level)
        )
        if discriminant < -tolerance:
            raise ArithmeticError(
                "quadratic target is unreachable on the oriented interpolated branch"
            )
        discriminant = max(0.0, discriminant)
        root = math.sqrt(discriminant)
        parameter = float((-cross + root) / direction_norm_sq)
        coefficients = particular + parameter * direction

        if with_derivative:
            q0_eta = float(
                moments["base_energy_eta"]
                - moments["swirl_energy_eta"]
                + 2.0 * (particular_eta @ gram @ particular)
            )
            cross_eta = float(
                particular_eta @ gram @ direction
                + particular @ gram @ direction_eta
            )
            direction_norm_sq_eta = float(2.0 * (direction_eta @ gram @ direction))
            level_eta = float(moments["target_eta"] - q0_eta)
            denominator = 2.0 * direction_norm_sq * parameter + 2.0 * cross
            if abs(denominator) <= 1.0e-14:
                raise ArithmeticError("quadratic branch derivative has a pole")
            parameter_eta = float(
                -(
                    direction_norm_sq_eta * parameter * parameter
                    + 2.0 * cross_eta * parameter
                    - level_eta
                )
                / denominator
            )
            coefficients_eta = (
                particular_eta + parameter_eta * direction + parameter * direction_eta
            )
        else:
            parameter_eta = 0.0
            coefficients_eta = np.zeros(4, dtype=float)

        # Enforce the mass identity from the exact 257-node core primitive at
        # every Z query.  Applying the differentiated identity gives the same
        # analytic eta derivative for the recovered fourth coefficient.
        mass_defect = float(moments["base_mass"] + self._bump_mass @ coefficients)
        coefficients = coefficients.copy()
        coefficients[-1] -= mass_defect / self._bump_mass[-1]
        if with_derivative:
            mass_eta_defect = float(
                moments["base_mass_eta"] + self._bump_mass @ coefficients_eta
            )
            coefficients_eta = coefficients_eta.copy()
            coefficients_eta[-1] -= mass_eta_defect / self._bump_mass[-1]
        singular_values = np.linalg.svd(scaled, compute_uv=False)
        linear_condition = float(singular_values[0] / singular_values[-1])
        return {
            "coefficients": coefficients,
            "coefficients_eta": coefficients_eta,
            "particular": particular,
            "direction": direction,
            "direction_eta": direction_eta,
            "orientation_dot": orientation_dot,
            "direction_norm": direction_norm,
            "linear_condition": linear_condition,
            "discriminant": discriminant,
            "parameter": parameter,
            "parameter_eta": parameter_eta,
            "q0": q0,
            "moments": moments,
        }

    @lru_cache(maxsize=8192)
    def _coefficients_at(self, eta: float) -> tuple[np.ndarray, np.ndarray]:
        """Return coefficients and derivatives from the interpolated moment solve."""

        solved = self._solve_interpolated(eta, with_derivative=True)
        coefficients = np.asarray(solved["coefficients"], dtype=float)
        eta_derivative = np.asarray(solved["coefficients_eta"], dtype=float)
        coefficients.setflags(write=False)
        eta_derivative.setflags(write=False)
        return coefficients, eta_derivative

    @lru_cache(maxsize=4096)
    def _bump_primitive(self, radius: float, index: int) -> float:
        radius = _finite(radius, "R")
        if not 0 <= index < len(self.bump_supports):
            raise ValueError("invalid bump index")
        left, right = self.bump_supports[index]
        if radius <= left:
            return 0.0
        if radius >= right:
            return float(self._bump_mass[index])
        nodes, weights = leggauss(96)
        points = left + (radius - left) * (nodes + 1.0) / 2.0
        radial_weights = (radius - left) * weights / 2.0
        return float(
            radial_weights
            @ np.asarray(
                [_flat_bump(float(point), left, right) for point in points],
                dtype=float,
            )
        )

    def F(self, R: Any, Z: Any) -> Any:
        """Return the supplied actual extended swirl scalar ``F(R,Z)``."""

        return _broadcast_apply(self._F_scalar, R, Z)

    def _F_scalar(self, radius: float, eta: float) -> float:
        radius, eta = self._validate_point(radius, eta)
        return float(self.swirl.F(radius, eta))

    def U(self, R: Any, Z: Any) -> Any:
        """Return the finite repaired axial profile, zero for ``R>=R_a``."""

        return _broadcast_apply(self._U_scalar, R, Z)

    def _U_scalar(self, radius: float, eta: float) -> float:
        radius, eta = self._validate_point(radius, eta)
        if radius >= self.Rmax:
            return 0.0
        if radius <= self.R_core:
            return float(self.core.U(radius, eta))
        coefficients, _ = self._coefficients_at(eta)
        correction = sum(
            coefficient * _flat_bump(radius, *support)
            for coefficient, support in zip(coefficients, self.bump_supports)
        )
        return float(self._base_U_scalar(radius, eta) + correction)

    def dU_deta(self, R: Any, Z: Any) -> Any:
        """Return the analytic eta derivative of the same interpolated U."""

        return _broadcast_apply(self._dU_deta_scalar, R, Z)

    def _dU_deta_scalar(self, radius: float, eta: float) -> float:
        radius, eta = self._validate_point(radius, eta)
        if radius >= self.Rmax:
            return 0.0
        if radius <= self.R_core:
            return float(self.core.dU_deta(radius, eta))
        _, eta_derivative = self._coefficients_at(eta)
        correction = sum(
            coefficient * _flat_bump(radius, *support)
            for coefficient, support in zip(eta_derivative, self.bump_supports)
        )
        return float(self._base_U_eta_scalar(radius, eta) + correction)

    def average_U(self, R: Any, Z: Any) -> Any:
        """Return ``R**-1 * integral_0^R U`` with an exact zero tail."""

        return _broadcast_apply(self._average_U_scalar, R, Z)

    def _average_U_scalar(self, radius: float, eta: float) -> float:
        radius, eta = self._validate_point(radius, eta)
        if radius >= self.Rmax:
            return 0.0
        if radius == 0.0:
            return float(self.core.U(0.0, eta))
        if radius <= self.R_core:
            return float(self.core.radial_average_U(radius, eta))
        coefficients, _ = self._coefficients_at(eta)
        primitive = float(
            self._base_primitive_coefficients(radius)
            @ self._core_coefficients(eta)
        )
        primitive += sum(
            coefficient * self._bump_primitive(radius, index)
            for index, coefficient in enumerate(coefficients)
        )
        return float(primitive / radius)

    def average_dU_deta(self, R: Any, Z: Any) -> Any:
        """Return the eta derivative of ``average_U`` using the same gauge."""

        return _broadcast_apply(self._average_dU_deta_scalar, R, Z)

    def _average_dU_deta_scalar(self, radius: float, eta: float) -> float:
        radius, eta = self._validate_point(radius, eta)
        if radius >= self.Rmax:
            return 0.0
        if radius == 0.0:
            return float(self.core.dU_deta(0.0, eta))
        if radius <= self.R_core:
            return float(self.core.radial_average_dU_deta(radius, eta))
        _, eta_derivative = self._coefficients_at(eta)
        primitive = float(
            self._base_primitive_coefficients(radius)
            @ self._core_coefficients(eta, eta_derivative=True)
        )
        primitive += sum(
            coefficient * self._bump_primitive(radius, index)
            for index, coefficient in enumerate(eta_derivative)
        )
        return float(primitive / radius)

    def quadratic_target(self, Z: Any) -> Any:
        """Return the Chebyshev-interpolated inward quadratic target."""

        eta = np.asarray(Z, dtype=float)
        if not np.all(np.isfinite(eta)) or np.any(np.abs(eta) > 1.0):
            raise ValueError("Z must be finite and satisfy |Z|<=1")
        result = self.z_grid.interpolate(self._quadratic_target_grid, eta)
        return float(result) if eta.ndim == 0 else np.asarray(result, dtype=float)

    @property
    def profile(self) -> LeadingProfile:
        """Return the repository ``LeadingProfile`` adapter."""

        if self._profile_cache is None:
            self._profile_cache = LeadingProfile(
                E=lambda R, Z: 0.0
                if float(R) == 0.0
                else math.sqrt(2.0 * float(R)) * self._F_scalar(float(R), float(Z)),
                U=self.U,
                dU_deta=self.dU_deta,
                Pi=None,
                F=self.F,
                average_U=self.average_U,
                average_dU_deta=self.average_dU_deta,
                name="lei-ren-smooth-extended-axial-repair",
                paper_exact=False,
                provenance=(
                    "Finite Z-Chebyshev shared-pressure ExtendedSwirl with four "
                    "fixed radial C-infinity bumps and oriented quadratic repair; "
                    "source arXiv:2609.35406v1, Section 2.5."
                ),
            )
        return self._profile_cache

    def metadata(self) -> dict[str, Any]:
        """Return JSON-ready construction and branch diagnostics."""

        return {
            "source": SOURCE,
            "source_version": SOURCE_VERSION,
            "source_section": SOURCE_SECTION,
            "profile_source": "lei_ren_part1_extended_pressure_core.load_extended_profile()",
            "class": type(self).__name__,
            "R_core": self.R_core,
            "R_anchor": self.R_anchor,
            "R_cutoff": self.R_cutoff,
            "Rmax_Ra": self.Rmax,
            "bump_supports": [list(item) for item in self.bump_supports],
            "bump_fractions": [list(item) for item in BUMP_FRACTIONS],
            "z_nodes": self.z_nodes,
            "z_range": [-1.0, 1.0],
            "radial_quadrature_order": self.radial_quadrature_order,
            "target_quadrature_order": self.target_quadrature_order,
            "log_break_count": LOG_BREAK_COUNT,
            "reference_coefficients": self.reference_coefficients.tolist(),
            "branch_orientation": {
                "rule": "positive inner product with fixed reference projection",
                "orientation_dot_min": float(np.min(self._orientation_dot_grid)),
                "orientation_dot_max": float(np.max(self._orientation_dot_grid)),
                "all_positive_at_nodes": bool(np.all(self._orientation_dot_grid > 0.0)),
                "quadratic_root": "positive parameter along the oriented direction",
                "nonzero_quadratic_target": True,
            },
            "conditioning": self._conditioning,
            "mass_recovery": (
                "last bump coefficient is recovered from base mass and the first "
                "three coefficients using the same 257-node core primitive"
            ),
            "coefficient_gauge": (
                "Chebyshev interpolation of moment functions followed by an "
                "oriented fixed-reference algebraic solve; not the full-nullspace "
                "minimum-L2 gauge"
            ),
            "replay_seed": {
                "schema_version": SEED_SCHEMA_VERSION,
                "save_method": "SmoothExtendedAxial.save_seed(path)",
                "load_factory": "load_smooth_extended_axial_from_seed(path)",
                "pressure_receipt_hash_bound": True,
            },
            "moment_interpolation": {
                "base_mixed": True,
                "base_energy": True,
                "bump_mixed": True,
                "swirl_energy": True,
                "quadratic_target": True,
                "linear_and_quadratic_re_solve_at_query": True,
                "analytic_eta_derivative": True,
            },
            "primitive_tail": "U, dU_deta, average_U, average_dU_deta are exactly zero for R>=R_a",
            "five_moment_closure": False,
            "pde_validated": False,
            "whole_cone_validated": False,
            "scope": (
                "Finite smooth axial repair of the shared-pressure ExtendedSwirl "
                "candidate. No complete Part I outer construction, PDE, stress "
                "cone, global energy, or five-moment certificate is claimed."
            ),
        }

    def save_seed(self, path: Any) -> Path:
        """Save a replayable JSON seed for this exact pressure/core receipt."""

        destination = Path(path)
        moments = {
            "base_mass": np.asarray(self._base_mass_grid, dtype=float).tolist(),
            "base_mixed": np.asarray(self._base_mixed_grid, dtype=float).tolist(),
            "base_energy": np.asarray(self._base_energy_grid, dtype=float).tolist(),
            "base_bump_cross": np.asarray(
                self._base_bump_cross_grid, dtype=float
            ).tolist(),
            "bump_mass": np.asarray(
                [row.bump_mass for row in self._slice_data], dtype=float
            ).T.tolist(),
            "bump_mixed": np.asarray(self._bump_mixed_grid, dtype=float).tolist(),
            "bump_gram": np.asarray(self._bump_gram_grid, dtype=float).tolist(),
            "swirl_energy": np.asarray(self._swirl_energy_grid, dtype=float).tolist(),
            "quadratic_target": np.asarray(
                self._quadratic_target_grid, dtype=float
            ).tolist(),
            "cross_radius": np.asarray(self._cross_radius_grid, dtype=float).tolist(),
        }
        seed = {
            "schema_version": SEED_SCHEMA_VERSION,
            "source": SOURCE,
            "source_version": SOURCE_VERSION,
            "algorithm": "Chebyshev moment interpolation plus oriented algebraic solve",
            "pressure_receipt": _pressure_receipt_identity(),
            "profile_source": "lei_ren_part1_extended_pressure_core.load_extended_profile()",
            "z_nodes": self.z_nodes,
            "radial_quadrature_order": self.radial_quadrature_order,
            "target_quadrature_order": self.target_quadrature_order,
            "geometry": {
                "R_core": self.R_core,
                "R_anchor": self.R_anchor,
                "R_cutoff": self.R_cutoff,
                "Rmax_Ra": self.Rmax,
                "bump_fractions": [list(item) for item in BUMP_FRACTIONS],
                "bump_supports": [list(item) for item in self.bump_supports],
            },
            "reference_coefficients": self.reference_coefficients.tolist(),
            "z_grid": np.asarray(self.z_grid.eta, dtype=float).tolist(),
            "core_inputs": {
                "u_coefficients": np.asarray(
                    self._core_u_coefficients, dtype=float
                ).tolist(),
                "u_eta_coefficients": np.asarray(
                    self._core_u_eta_coefficients, dtype=float
                ).tolist(),
                "base_primitive_full": np.asarray(
                    self._base_primitive_full, dtype=float
                ).tolist(),
                "base_primitive_eta_full": np.asarray(
                    self._base_primitive_eta_full, dtype=float
                ).tolist(),
            },
            "moments": moments,
            "metadata": self.metadata(),
        }
        destination.write_text(json.dumps(seed, indent=2) + "\n", encoding="utf-8")
        return destination

    @classmethod
    def from_seed(cls, path: Any, *, swirl: Any | None = None) -> "SmoothExtendedAxial":
        """Construct from a validated replay seed without radial precomputation."""

        seed = _read_seed(path)
        return cls(
            swirl=swirl,
            z_nodes=int(seed["z_nodes"]),
            radial_quadrature_order=int(seed["radial_quadrature_order"]),
            target_quadrature_order=int(seed["target_quadrature_order"]),
            _seed_data=seed,
        )


def load_smooth_extended_axial_profile(**kwargs: Any) -> SmoothExtendedAxial:
    """Construct a smooth finite axial adapter on the saved extended profile."""

    return SmoothExtendedAxial(**kwargs)


def load_smooth_extended_axial_from_seed(
    path: Any, *, swirl: Any | None = None
) -> SmoothExtendedAxial:
    """Explicitly replay a saved seed after validating its pressure receipt."""

    return SmoothExtendedAxial.from_seed(path, swirl=swirl)


def build_smooth_extended_axial(**kwargs: Any) -> SmoothExtendedAxial:
    """Compatibility factory for the public ``SmoothExtendedAxial`` class."""

    return SmoothExtendedAxial(**kwargs)


__all__ = [
    "BUMP_FRACTIONS",
    "R_CORE",
    "R_CUTOFF",
    "SOURCE",
    "SOURCE_SECTION",
    "SOURCE_VERSION",
    "SmoothExtendedAxial",
    "build_smooth_extended_axial",
    "load_smooth_extended_axial_from_seed",
    "load_smooth_extended_axial_profile",
]
