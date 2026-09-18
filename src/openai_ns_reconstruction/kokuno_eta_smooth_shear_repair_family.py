"""Eta-smooth PA.17 repair family for the selected Kokuno-form modulation.

This module advances only the leading/base-profile lane.  It consumes the
source-form finite-frequency radial insertion exposed by
``KokunoRadialModulationDiscrepancy`` and the later PA.17 I1 inverse exposed by
``KokunoShearMomentRepair``.  The public corrected reconstruction requires the
five modulation-induced prefix-moment discrepancies to be restored by a smooth
coefficient family ``c(eta)`` on I1.

The corrected public source does not publish one hidden numerical admissible
loop.  Therefore the upstream modulation remains the explicitly autonomous,
source-form realization already recorded by that module.  This file does not
rename it recovered Kokuno data.  It only turns the already selected
provenance-labelled realization into a differentiable coefficient family that
can be used directly by the I1 velocity-correction evaluator.

For the selected realization the target depends on eta only through
``f(eta)=1/(1+eta^2)``.  We therefore interpolate the five solved coefficients
as global Chebyshev polynomials in ``s=eta^2`` on ``[0,1]``.  This makes the
family globally smooth and exactly even by construction; its analytic eta
Derivative is obtained by the chain rule.  Interpolation nodes and degree are
repository-autonomous numerical choices, not source constants.

The resulting object exposes

    coefficients(eta), coefficient_eta(eta),
    profile_correction_logX(log_X, eta),
    velocity_correction(x,y,z,t) -> [...,3].

It is a local I1 repair contribution, not yet a global matched leading field,
not independent PDE validation, and not paper-exact.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from numpy.polynomial import chebyshev as cheb

from .kokuno_radial_modulation_discrepancy import KokunoRadialModulationDiscrepancy
from .kokuno_shear_moment_repair import KokunoShearMomentRepair
from .kokuno_similarity_coordinates import KokunoNativeSimilarityCoordinates


SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-eta-smooth-shear-repair-family-v1"
MOMENT_ORDER = ("M", "J", "I", "S", "C_p")

_SOURCE_FORMULAS = {
    "radial_insertion": "E_N=E exp(A/N); U_N=U+B/N with phase N log X",
    "repair_location": "the five O_m(N^-1) prefix-moment discrepancies are restored on reserved I1",
    "PA17_background": "U=0; E=K(eta) X^(-1/2-lambda)",
    "PA17_rows": "row order (M,J,I,S,C_p)",
    "smooth_family_requirement": "repair coefficients are selected smoothly in eta",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "selected_autonomous_modulation_consumed": True,
    "eta_smooth_repair_coefficient_family_reconstructed": True,
    "I1_repair_applied_to_selected_autonomous_modulation": True,
    "native_I1_velocity_correction_executable": True,
    "autonomous_eta_interpolation": True,
    "source_hidden_loop_parameters_recovered": False,
    "source_admissible_loop_reconstructed": False,
    "actual_source_admissible_loop_discrepancy_supplied": False,
    "I1_cone_repair_applied_to_actual_source_modulation": False,
    "cone_modulation_completed": False,
    "global_pressure_matched": False,
    "global_leading_profile_reconstructed": False,
    "complete_kokuno_composite_velocity": False,
    "formal_full_domain_pde_gate_assessed": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _finite_array(value: Any, name: str) -> np.ndarray:
    array = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain only finite values")
    return array


def _eta_array(value: Any) -> np.ndarray:
    eta = _finite_array(value, "eta")
    if np.any(np.abs(eta) > 1.0):
        raise ValueError("eta must lie in [-1,1]")
    return eta


def _evaluate_chebyshev(coefficients: np.ndarray, x: np.ndarray) -> np.ndarray:
    values = cheb.chebval(np.asarray(x, dtype=float), coefficients)
    return np.moveaxis(np.asarray(values, dtype=float), 0, -1)


@dataclass(frozen=True)
class KokunoEtaSmoothShearRepairFamily:
    """Smooth PA.17 coefficient family for the selected autonomous insertion."""

    modulation: KokunoRadialModulationDiscrepancy = field(
        default_factory=KokunoRadialModulationDiscrepancy
    )
    interpolation_nodes: int = 11
    closure_tolerance: float = 5.0e-7

    _repair: KokunoShearMomentRepair = field(init=False, repr=False, compare=False)
    _coordinates: KokunoNativeSimilarityCoordinates = field(init=False, repr=False, compare=False)
    _s_nodes: np.ndarray = field(init=False, repr=False, compare=False)
    _eta_nodes: np.ndarray = field(init=False, repr=False, compare=False)
    _nodal_coefficients: np.ndarray = field(init=False, repr=False, compare=False)
    _poly_coefficients: np.ndarray = field(init=False, repr=False, compare=False)
    _poly_derivative: np.ndarray = field(init=False, repr=False, compare=False)
    _nodal_max_residual: float = field(init=False, repr=False, compare=False)
    _holdout_max_residual: float = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        if not isinstance(self.modulation, KokunoRadialModulationDiscrepancy):
            raise TypeError("modulation must be a KokunoRadialModulationDiscrepancy")
        if isinstance(self.interpolation_nodes, bool) or not isinstance(
            self.interpolation_nodes, (int, np.integer)
        ):
            raise TypeError("interpolation_nodes must be an integer")
        count = int(self.interpolation_nodes)
        if not 7 <= count <= 21:
            raise ValueError("interpolation_nodes must lie in [7,21]")
        tolerance = float(self.closure_tolerance)
        if not math.isfinite(tolerance) or not (0.0 < tolerance <= 1.0e-5):
            raise ValueError("closure_tolerance must lie in (0,1e-5]")

        repair = KokunoShearMomentRepair(outer_schedule=self.modulation.outer_schedule)
        coordinates = KokunoNativeSimilarityCoordinates(h=float(self.modulation.outer_schedule.h))

        # Chebyshev-Lobatto nodes in s=eta^2, ordered from 0 to 1.
        j = np.arange(count, dtype=float)
        s_nodes = 0.5 * (1.0 - np.cos(math.pi * j / float(count - 1)))
        eta_nodes = np.sqrt(s_nodes)
        solved: list[np.ndarray] = []
        residuals: list[float] = []
        for eta in eta_nodes:
            result = self.modulation.solve_repair_at_eta(float(eta))
            if not result.success:
                raise RuntimeError("pointwise PA.17 solve failed at an interpolation node")
            coefficients = np.asarray(result.coefficients, dtype=float)
            if np.any(~np.isfinite(coefficients)):
                raise RuntimeError("pointwise PA.17 solve returned non-finite coefficients")
            if np.any(np.abs(coefficients) > repair.coefficient_limit + 1.0e-14):
                raise RuntimeError("pointwise PA.17 solve exceeded the declared coefficient bound")
            solved.append(coefficients)
            residuals.append(float(result.max_abs_residual))
        nodal_coefficients = np.asarray(solved, dtype=float)

        # Interpolate all five columns at once in the Chebyshev basis.
        x_nodes = 2.0 * s_nodes - 1.0
        vandermonde = cheb.chebvander(x_nodes, count - 1)
        poly_coefficients = np.linalg.solve(vandermonde, nodal_coefficients)
        poly_derivative = cheb.chebder(poly_coefficients, axis=0)

        for array in (
            s_nodes,
            eta_nodes,
            nodal_coefficients,
            poly_coefficients,
            poly_derivative,
        ):
            array.setflags(write=False)

        object.__setattr__(self, "interpolation_nodes", count)
        object.__setattr__(self, "closure_tolerance", tolerance)
        object.__setattr__(self, "_repair", repair)
        object.__setattr__(self, "_coordinates", coordinates)
        object.__setattr__(self, "_s_nodes", s_nodes)
        object.__setattr__(self, "_eta_nodes", eta_nodes)
        object.__setattr__(self, "_nodal_coefficients", nodal_coefficients)
        object.__setattr__(self, "_poly_coefficients", poly_coefficients)
        object.__setattr__(self, "_poly_derivative", poly_derivative)
        object.__setattr__(self, "_nodal_max_residual", max(residuals, default=0.0))

        # Midpoints in s are not interpolation nodes and provide a deterministic
        # fail-closed guard against a polynomial that only replays the training nodes.
        s_mid = 0.5 * (s_nodes[:-1] + s_nodes[1:])
        eta_mid = np.sqrt(s_mid)
        holdout_residual = self._closure_residual(eta_mid)
        holdout_max = float(np.max(np.abs(holdout_residual))) if holdout_residual.size else 0.0
        object.__setattr__(self, "_holdout_max_residual", holdout_max)
        if holdout_max > tolerance:
            raise RuntimeError(
                "eta-smooth PA.17 family misses its preregistered held-out moment-closure tolerance"
            )

    @property
    def repair(self) -> KokunoShearMomentRepair:
        return self._repair

    @property
    def eta_nodes(self) -> np.ndarray:
        return self._eta_nodes.copy()

    @property
    def nodal_coefficients(self) -> np.ndarray:
        return self._nodal_coefficients.copy()

    def coefficients(self, eta: Any) -> np.ndarray:
        eta_array = _eta_array(eta)
        x = 2.0 * eta_array * eta_array - 1.0
        return _evaluate_chebyshev(self._poly_coefficients, x)

    def coefficient_eta(self, eta: Any) -> np.ndarray:
        eta_array = _eta_array(eta)
        x = 2.0 * eta_array * eta_array - 1.0
        dcdx = _evaluate_chebyshev(self._poly_derivative, x)
        return dcdx * (4.0 * eta_array)[..., None]

    def _closure_residual(self, eta: Any) -> np.ndarray:
        eta_array = _eta_array(eta)
        shape = eta_array.shape
        flat_eta = eta_array.reshape(-1)
        flat_coefficients = self.coefficients(flat_eta).reshape((-1, 5))
        target = np.asarray(self.modulation.repair_target_normalized(flat_eta), dtype=float).reshape(
            (-1, 5)
        )
        achieved = np.empty_like(target)
        f_values = np.asarray(self.modulation.outer_schedule.source_f(flat_eta), dtype=float).reshape(-1)
        for index, (coefficients, f_eta) in enumerate(zip(flat_coefficients, f_values)):
            achieved[index] = self._repair.normalized_increment(
                coefficients,
                f_eta=float(f_eta),
            )
        return (achieved - target).reshape(shape + (5,))

    def closure_report(self, eta: Any) -> dict[str, Any]:
        eta_array = _eta_array(eta)
        residual = self._closure_residual(eta_array)
        return {
            "eta": np.asarray(eta_array, dtype=float).tolist(),
            "moment_order": list(MOMENT_ORDER),
            "residual_normalized": np.asarray(residual, dtype=float).tolist(),
            "max_abs_residual": float(np.max(np.abs(residual))) if residual.size else 0.0,
            "closure_tolerance": self.closure_tolerance,
            "passed": bool(np.all(np.abs(residual) <= self.closure_tolerance)),
        }

    def family_report(self) -> dict[str, Any]:
        return {
            "interpolation_variable": "s=eta^2",
            "interpolation_basis": "global Chebyshev polynomial on [0,1]",
            "interpolation_nodes": self.interpolation_nodes,
            "polynomial_degree": self.interpolation_nodes - 1,
            "eta_even_by_construction": True,
            "analytic_eta_derivative": True,
            "nodal_max_abs_moment_residual": self._nodal_max_residual,
            "heldout_midpoint_max_abs_moment_residual": self._holdout_max_residual,
            "closure_tolerance": self.closure_tolerance,
            "max_abs_coefficient": float(np.max(np.abs(self._nodal_coefficients))),
            "coefficient_limit": self._repair.coefficient_limit,
            "source_admissible_loop_reconstructed": False,
            "selected_autonomous_modulation_consumed": True,
        }

    def profile_correction_logX(self, log_X: Any, eta: Any) -> dict[str, np.ndarray]:
        """Evaluate the I1 repair profile using the smooth c(eta) family."""

        log_array, eta_array = np.broadcast_arrays(
            _finite_array(log_X, "log_X"), _eta_array(eta)
        )
        shape = log_array.shape
        flat_log = log_array.reshape(-1)
        flat_eta = eta_array.reshape(-1)
        coefficients = self.coefficients(flat_eta).reshape((-1, 5))
        derivatives = self.coefficient_eta(flat_eta).reshape((-1, 5))
        keys = (
            "xi",
            "delta_E",
            "delta_E_X",
            "delta_E_eta",
            "delta_F",
            "delta_F_X",
            "delta_F_eta",
            "delta_U",
            "delta_U_X",
            "delta_U_eta",
            "delta_M",
            "delta_M_eta",
            "delta_v0",
        )
        collected = {key: [] for key in keys}
        for log_value, eta_value, coeff, coeff_eta in zip(
            flat_log, flat_eta, coefficients, derivatives
        ):
            values = self._repair.profile_correction_logX(
                float(log_value),
                float(eta_value),
                coefficients=coeff,
                coefficient_eta=coeff_eta,
            )
            for key in keys:
                collected[key].append(float(np.asarray(values[key]).item()))
        return {
            key: np.asarray(values, dtype=float).reshape(shape)
            for key, values in collected.items()
        }

    def velocity_correction(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        """Return the selected-modulation I1 repair contribution in Cartesian form."""

        coordinates = self._coordinates.evaluate(x, y, z, t)
        profiles = self.profile_correction_logX(
            np.log(np.asarray(coordinates["X"], dtype=float)),
            np.asarray(coordinates["eta"], dtype=float),
        )
        x_array, y_array, q = np.broadcast_arrays(
            _finite_array(x, "x"),
            _finite_array(y, "y"),
            np.asarray(coordinates["q"], dtype=float),
        )
        x_array = np.broadcast_to(x_array, np.shape(profiles["delta_F"]))
        y_array = np.broadcast_to(y_array, np.shape(profiles["delta_F"]))
        q = np.broadcast_to(q, np.shape(profiles["delta_F"]))
        h = float(self.modulation.outer_schedule.h)
        swirl = np.power(q, -1.0 - h) * profiles["delta_F"]
        radial = profiles["delta_v0"] / (2.0 * q)
        axial = np.power(q, -0.5 - h) * profiles["delta_U"]
        return np.stack(
            (
                radial * x_array - swirl * y_array,
                radial * y_array + swirl * x_array,
                axial,
            ),
            axis=-1,
        )

    def _unsigned_payload(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "corrected_release": CORRECTED_RELEASE,
                "corrected_release_date": CORRECTED_RELEASE_DATE,
            },
            "source_formulas": dict(_SOURCE_FORMULAS),
            "modulation": self.modulation.to_payload(),
            "interpolation_nodes": self.interpolation_nodes,
            "closure_tolerance": self.closure_tolerance,
            "s_nodes": [float(value) for value in self._s_nodes],
            "nodal_coefficients": self._nodal_coefficients.tolist(),
            "chebyshev_coefficients": self._poly_coefficients.tolist(),
            "family_report": self.family_report(),
            "truth_boundary": dict(_TRUTH_BOUNDARY),
        }

    @property
    def sha256(self) -> str:
        return hashlib.sha256(
            _canonical_json(self._unsigned_payload()).encode("utf-8")
        ).hexdigest()

    def to_payload(self) -> dict[str, Any]:
        payload = self._unsigned_payload()
        payload["sha256"] = self.sha256
        return payload

    @classmethod
    def from_payload(cls, payload: dict[str, Any]) -> "KokunoEtaSmoothShearRepairFamily":
        if not isinstance(payload, dict) or payload.get("schema") != SCHEMA:
            raise ValueError("unexpected Kokuno eta-smooth shear-repair schema")
        modulation_payload = payload.get("modulation")
        if not isinstance(modulation_payload, dict):
            raise ValueError("modulation payload is missing")
        obj = cls(
            modulation=KokunoRadialModulationDiscrepancy.from_payload(modulation_payload),
            interpolation_nodes=int(payload.get("interpolation_nodes")),
            closure_tolerance=float(payload.get("closure_tolerance")),
        )
        if obj.to_payload() != payload:
            raise ValueError("eta-smooth shear-repair payload hash or content mismatch")
        return obj

    def save_json(self, path: str | Path) -> Path:
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(
            json.dumps(self.to_payload(), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        return target

    @classmethod
    def load_json(cls, path: str | Path) -> "KokunoEtaSmoothShearRepairFamily":
        return cls.from_payload(json.loads(Path(path).read_text(encoding="utf-8")))
