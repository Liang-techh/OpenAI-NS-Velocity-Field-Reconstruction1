"""Executable five-moment repair primitive from the corrected Kokuno reader.

This module implements the *local repair map*, not the still-missing final
outer/heat assembly.  The pinned 2026-09-09 reconstruction keeps five prefix
moments and, on the inner-join scaled coordinate ``x``, perturbs the ideal pair
with two compact ``U`` bumps and three compact ``E`` bumps.  The exact nonlinear
increment integrands are the source's PA.14 and the physical/scaled moment map
is PA.15.

The source proves local solvability using ordered nonnegative bumps and the
five-by-five block structure PA.16.  It does not prescribe numerical bump
centres for this repository.  The concrete C-infinity bump locations and the
finite coefficient bound below are therefore explicit autonomous numerical
choices.  A solved repair is evidence that the local five-moment map is
executable; it is not a completed Kokuno outer profile or an NS validation.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import least_squares


SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-five-moment-repair-v1"

# Source inner-join interval in the scaled coordinate x=X/X_R.
X_SCALED_MIN = math.exp(-6.0)
X_SCALED_MAX = math.exp(-5.0)

# Autonomous deterministic C-infinity bump placement.  The source requires two
# ordered nonnegative U bumps and three ordered nonnegative E bumps but does not
# supply repository-ready numerical centres/widths.
_U_LOG_CENTERS = (-5.82, -5.18)
_E_LOG_CENTERS = (-5.86, -5.50, -5.14)
_LOG_HALF_WIDTH = 0.105

_SOURCE_FORMULAS = {
    "five_moments": (
        "M=int U dx; I=int H dx; J=int U H dx; "
        "S=int(U^2-E^2/2)dx; C_p=int E^2/(2x)dx"
    ),
    "PA14_exact_increment": (
        "(u, sqrt(2x)e, H u + U sqrt(2x)e + sqrt(2x)u e, "
        "2Uu+u^2-Ee-e^2/2, Ee/x+e^2/(2x))"
    ),
    "PA15_scaling": (
        "(M,I,J,S,C_p)=(X_R*Mhat,X_R^(3/2)*Ihat,"
        "X_R^(3/2)*Jhat,X_R*Shat,Cphat)"
    ),
    "inner_base_pair": "U=4 eta; E=P_* f x^(1/10) on exp(-6)<x<exp(-5)",
    "PA16_row_operations": "J->J-4 eta I; S->S-8 eta M",
    "PA16_linear_blocks": (
        "U: (1,sqrt(2) P_* f x^(3/5)); "
        "E: (sqrt(2)x^(1/2),-P_*f x^(1/10),P_*f x^(-9/10))"
    ),
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "source_exact_five_moment_integrands_implemented": True,
    "source_inner_join_scaled_coordinate_implemented": True,
    "autonomous_bump_shapes_and_locations": True,
    "autonomous_coefficient_bound": True,
    "source_outer_target_discrepancy_supplied": False,
    "eta_smooth_coefficient_family_reconstructed": False,
    "cone_modulation_completed": False,
    "heat_compensation_completed": False,
    "core_to_heat_matching_completed": False,
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


def _coefficient_array(value: Any, limit: float) -> np.ndarray:
    coefficients = _finite_array(value, "coefficients")
    if coefficients.shape != (5,):
        raise ValueError("coefficients must have shape (5,)")
    if np.any(np.abs(coefficients) > limit + 1.0e-15):
        raise ValueError("repair coefficient exceeds the declared bound")
    return coefficients


def _cinfty_log_bump(x: np.ndarray, center: float) -> np.ndarray:
    """Peak-one C-infinity bump compactly supported in log(x)."""

    values = np.asarray(x, dtype=float)
    out = np.zeros_like(values)
    positive = values > 0.0
    if not np.any(positive):
        return out
    local = values[positive]
    s = (np.log(local) - center) / _LOG_HALF_WIDTH
    active = np.abs(s) < 1.0
    tmp = np.zeros_like(local)
    if np.any(active):
        sa = s[active]
        tmp[active] = np.exp(1.0 - 1.0 / (1.0 - sa * sa))
    out[positive] = tmp
    return out


def _cinfty_log_bump_dx(x: np.ndarray, center: float) -> np.ndarray:
    values = np.asarray(x, dtype=float)
    out = np.zeros_like(values)
    positive = values > 0.0
    if not np.any(positive):
        return out
    local = values[positive]
    s = (np.log(local) - center) / _LOG_HALF_WIDTH
    active = np.abs(s) < 1.0
    tmp = np.zeros_like(local)
    if np.any(active):
        sa = s[active]
        bump = np.exp(1.0 - 1.0 / (1.0 - sa * sa))
        tmp[active] = (
            bump
            * (-2.0 * sa / (1.0 - sa * sa) ** 2)
            / (_LOG_HALF_WIDTH * local[active])
        )
    out[positive] = tmp
    return out


@dataclass(frozen=True)
class FiveMomentSolveResult:
    coefficients: tuple[float, float, float, float, float]
    target_scaled_moments: tuple[float, float, float, float, float]
    achieved_scaled_moments: tuple[float, float, float, float, float]
    residual_scaled_moments: tuple[float, float, float, float, float]
    max_abs_residual: float
    nfev: int
    success: bool


@dataclass(frozen=True)
class KokunoFiveMomentRepair:
    """Bounded numerical realization of the source five-moment local inverse.

    Coefficients are ordered ``(u1,u2,e1,e2,e3)``.  The first two multiply
    compact perturbations of ``U`` and the last three compact perturbations of
    ``E``.  ``P_star`` and ``f_eta`` belong to the source ideal pair and are
    supplied by the future outer-profile assembly rather than guessed here.
    """

    quadrature_points: int = 192
    coefficient_limit: float = 0.05

    def __post_init__(self) -> None:
        if isinstance(self.quadrature_points, bool) or not isinstance(
            self.quadrature_points, (int, np.integer)
        ):
            raise TypeError("quadrature_points must be an integer")
        order = int(self.quadrature_points)
        if not 64 <= order <= 256:
            raise ValueError("quadrature_points must lie in [64,256]")
        limit = float(self.coefficient_limit)
        if not math.isfinite(limit) or not (0.0 < limit <= 0.2):
            raise ValueError("coefficient_limit must lie in (0,0.2]")
        object.__setattr__(self, "quadrature_points", order)
        object.__setattr__(self, "coefficient_limit", limit)

    @staticmethod
    def _validate_source_inputs(eta: float, P_star: float, f_eta: float) -> tuple[float, float, float]:
        eta = float(eta)
        P_star = float(P_star)
        f_eta = float(f_eta)
        if not all(math.isfinite(v) for v in (eta, P_star, f_eta)):
            raise ValueError("eta, P_star and f_eta must be finite")
        if not -1.0 <= eta <= 1.0:
            raise ValueError("eta must lie in [-1,1]")
        if P_star <= 0.0 or f_eta <= 0.0:
            raise ValueError("P_star and f_eta must be positive")
        return eta, P_star, f_eta

    def _quadrature(self) -> tuple[np.ndarray, np.ndarray]:
        nodes, weights = leggauss(self.quadrature_points)
        half = 0.5 * (X_SCALED_MAX - X_SCALED_MIN)
        center = 0.5 * (X_SCALED_MAX + X_SCALED_MIN)
        x = center + half * nodes
        w = half * weights
        return np.asarray(x, dtype=float), np.asarray(w, dtype=float)

    @staticmethod
    def basis_values(x: Any) -> np.ndarray:
        values = _finite_array(x, "x")
        return np.stack(
            [
                *(_cinfty_log_bump(values, center) for center in _U_LOG_CENTERS),
                *(_cinfty_log_bump(values, center) for center in _E_LOG_CENTERS),
            ],
            axis=0,
        )

    @staticmethod
    def basis_derivatives(x: Any) -> np.ndarray:
        values = _finite_array(x, "x")
        return np.stack(
            [
                *(_cinfty_log_bump_dx(values, center) for center in _U_LOG_CENTERS),
                *(_cinfty_log_bump_dx(values, center) for center in _E_LOG_CENTERS),
            ],
            axis=0,
        )

    def perturbation(self, x: Any, coefficients: Any) -> tuple[np.ndarray, np.ndarray]:
        values = _finite_array(x, "x")
        coeff = _coefficient_array(coefficients, self.coefficient_limit)
        basis = self.basis_values(values)
        u = np.tensordot(coeff[:2], basis[:2], axes=(0, 0))
        e = np.tensordot(coeff[2:], basis[2:], axes=(0, 0))
        return np.asarray(u, dtype=float), np.asarray(e, dtype=float)

    def perturbation_dx(self, x: Any, coefficients: Any) -> tuple[np.ndarray, np.ndarray]:
        values = _finite_array(x, "x")
        coeff = _coefficient_array(coefficients, self.coefficient_limit)
        basis = self.basis_derivatives(values)
        u_x = np.tensordot(coeff[:2], basis[:2], axes=(0, 0))
        e_x = np.tensordot(coeff[2:], basis[2:], axes=(0, 0))
        return np.asarray(u_x, dtype=float), np.asarray(e_x, dtype=float)

    def corrected_scaled_profiles(
        self,
        x: Any,
        eta: float,
        P_star: float,
        f_eta: float,
        coefficients: Any,
    ) -> dict[str, np.ndarray]:
        eta, P_star, f_eta = self._validate_source_inputs(eta, P_star, f_eta)
        values = _finite_array(x, "x")
        if np.any(values <= 0.0):
            raise ValueError("scaled x must be positive")
        u, e = self.perturbation(values, coefficients)
        u_x, e_x = self.perturbation_dx(values, coefficients)
        U0 = np.full_like(values, 4.0 * eta, dtype=float)
        E0 = P_star * f_eta * values ** 0.1
        U = U0 + u
        E = E0 + e
        if np.any(E <= 0.0):
            raise ValueError("repair would violate positivity of the source E profile")
        return {
            "U": U,
            "E": E,
            "U_x": u_x,
            "E_x": 0.1 * E0 / values + e_x,
            "delta_U": u,
            "delta_E": e,
        }

    def physical_profile_correction(
        self, X: Any, X_R: float, coefficients: Any
    ) -> dict[str, np.ndarray]:
        X_values = _finite_array(X, "X")
        X_R = float(X_R)
        if not math.isfinite(X_R) or X_R <= 0.0:
            raise ValueError("X_R must be positive and finite")
        x = X_values / X_R
        u, e = self.perturbation(x, coefficients)
        u_x, e_x = self.perturbation_dx(x, coefficients)
        return {
            "delta_U": u,
            "delta_E": e,
            "delta_U_X": u_x / X_R,
            "delta_E_X": e_x / X_R,
        }

    def moment_increment(
        self,
        coefficients: Any,
        *,
        eta: float,
        P_star: float,
        f_eta: float,
    ) -> np.ndarray:
        """Return exact PA.14 increments in scaled-moment order M,I,J,S,C_p."""

        eta, P_star, f_eta = self._validate_source_inputs(eta, P_star, f_eta)
        coeff = _coefficient_array(coefficients, self.coefficient_limit)
        x, weights = self._quadrature()
        basis = self.basis_values(x)
        u = np.tensordot(coeff[:2], basis[:2], axes=(0, 0))
        e = np.tensordot(coeff[2:], basis[2:], axes=(0, 0))
        U0 = 4.0 * eta
        E0 = P_star * f_eta * x ** 0.1
        root = np.sqrt(2.0 * x)
        H0 = root * E0
        integrands = np.stack(
            [
                u,
                root * e,
                H0 * u + U0 * root * e + root * u * e,
                2.0 * U0 * u + u * u - E0 * e - 0.5 * e * e,
                E0 * e / x + 0.5 * e * e / x,
            ],
            axis=0,
        )
        return np.asarray(integrands @ weights, dtype=float)

    def moment_jacobian(
        self,
        coefficients: Any,
        *,
        eta: float,
        P_star: float,
        f_eta: float,
    ) -> np.ndarray:
        """Analytic coefficient Jacobian of the exact PA.14 moment map."""

        eta, P_star, f_eta = self._validate_source_inputs(eta, P_star, f_eta)
        coeff = _coefficient_array(coefficients, self.coefficient_limit)
        x, weights = self._quadrature()
        basis = self.basis_values(x)
        u = np.tensordot(coeff[:2], basis[:2], axes=(0, 0))
        e = np.tensordot(coeff[2:], basis[2:], axes=(0, 0))
        U0 = 4.0 * eta
        E0 = P_star * f_eta * x ** 0.1
        root = np.sqrt(2.0 * x)
        H0 = root * E0
        jac = np.zeros((5, 5), dtype=float)
        zeros = np.zeros_like(x)
        for column, bump in enumerate(basis[:2]):
            integrands = np.stack(
                [
                    bump,
                    zeros,
                    (H0 + root * e) * bump,
                    (2.0 * U0 + 2.0 * u) * bump,
                    zeros,
                ],
                axis=0,
            )
            jac[:, column] = integrands @ weights
        for local, bump in enumerate(basis[2:]):
            integrands = np.stack(
                [
                    zeros,
                    root * bump,
                    root * (U0 + u) * bump,
                    -(E0 + e) * bump,
                    (E0 + e) * bump / x,
                ],
                axis=0,
            )
            jac[:, 2 + local] = integrands @ weights
        return jac

    @staticmethod
    def source_row_transform(values: Any, eta: float) -> np.ndarray:
        """Apply the PA.16 row operations and order rows by its two blocks."""

        vector = _finite_array(values, "moments")
        if vector.shape != (5,):
            raise ValueError("moments must have shape (5,)")
        eta = float(eta)
        if not math.isfinite(eta) or not -1.0 <= eta <= 1.0:
            raise ValueError("eta must lie in [-1,1]")
        M, I, J, S, C_p = vector
        U0 = 4.0 * eta
        return np.asarray(
            [M, J - U0 * I, I, S - 2.0 * U0 * M, C_p], dtype=float
        )

    @staticmethod
    def _source_row_matrix(eta: float) -> np.ndarray:
        matrix = np.zeros((5, 5), dtype=float)
        # Output order: M, J-U0 I, I, S-2U0 M, C_p.
        matrix[0, 0] = 1.0
        matrix[1, 1] = -4.0 * eta
        matrix[1, 2] = 1.0
        matrix[2, 1] = 1.0
        matrix[3, 0] = -8.0 * eta
        matrix[3, 3] = 1.0
        matrix[4, 4] = 1.0
        return matrix

    def transformed_jacobian(
        self,
        coefficients: Any,
        *,
        eta: float,
        P_star: float,
        f_eta: float,
    ) -> np.ndarray:
        return self._source_row_matrix(float(eta)) @ self.moment_jacobian(
            coefficients, eta=eta, P_star=P_star, f_eta=f_eta
        )

    @staticmethod
    def physical_moment_increment(scaled_increment: Any, X_R: float) -> np.ndarray:
        values = _finite_array(scaled_increment, "scaled_increment")
        if values.shape != (5,):
            raise ValueError("scaled_increment must have shape (5,)")
        X_R = float(X_R)
        if not math.isfinite(X_R) or X_R <= 0.0:
            raise ValueError("X_R must be positive and finite")
        return values * np.asarray(
            [X_R, X_R ** 1.5, X_R ** 1.5, X_R, 1.0], dtype=float
        )

    def solve(
        self,
        target_scaled_moments: Any,
        *,
        eta: float,
        P_star: float,
        f_eta: float,
        initial_coefficients: Any | None = None,
        max_nfev: int = 200,
        absolute_tolerance: float = 2.0e-11,
    ) -> FiveMomentSolveResult:
        """Solve a supplied small five-moment discrepancy inside fixed bounds.

        The target is supplied by the caller because the source outer/cone/heat
        assembly that determines the *actual* discrepancy is not reconstructed
        yet.  No target is inferred from residuals or images here.
        """

        eta, P_star, f_eta = self._validate_source_inputs(eta, P_star, f_eta)
        target = _finite_array(target_scaled_moments, "target_scaled_moments")
        if target.shape != (5,):
            raise ValueError("target_scaled_moments must have shape (5,)")
        if isinstance(max_nfev, bool) or int(max_nfev) < 1 or int(max_nfev) > 2000:
            raise ValueError("max_nfev must be an integer in [1,2000]")
        tolerance = float(absolute_tolerance)
        if not math.isfinite(tolerance) or not (0.0 < tolerance <= 1.0e-6):
            raise ValueError("absolute_tolerance must lie in (0,1e-6]")

        zero = np.zeros(5, dtype=float)
        source_matrix = self._source_row_matrix(eta)
        target_t = source_matrix @ target
        jac0_t = self.transformed_jacobian(
            zero, eta=eta, P_star=P_star, f_eta=f_eta
        )
        if np.linalg.matrix_rank(jac0_t) != 5:
            raise RuntimeError("source five-moment linearized map lost rank")

        if initial_coefficients is None:
            initial = np.linalg.solve(jac0_t, target_t)
            initial = np.clip(
                initial,
                -0.9 * self.coefficient_limit,
                0.9 * self.coefficient_limit,
            )
        else:
            initial = _coefficient_array(
                initial_coefficients, self.coefficient_limit
            ).copy()

        row_scale = np.maximum(
            np.linalg.norm(jac0_t, axis=1) * self.coefficient_limit,
            np.maximum(np.abs(target_t), 1.0e-14),
        )

        def residual(coefficients: np.ndarray) -> np.ndarray:
            achieved = self.moment_increment(
                coefficients, eta=eta, P_star=P_star, f_eta=f_eta
            )
            return (source_matrix @ (achieved - target)) / row_scale

        def jacobian(coefficients: np.ndarray) -> np.ndarray:
            return self.transformed_jacobian(
                coefficients, eta=eta, P_star=P_star, f_eta=f_eta
            ) / row_scale[:, None]

        result = least_squares(
            residual,
            initial,
            jac=jacobian,
            bounds=(-self.coefficient_limit, self.coefficient_limit),
            method="trf",
            ftol=1.0e-13,
            xtol=1.0e-13,
            gtol=1.0e-13,
            max_nfev=int(max_nfev),
        )
        coefficients = np.asarray(result.x, dtype=float)
        achieved = self.moment_increment(
            coefficients, eta=eta, P_star=P_star, f_eta=f_eta
        )
        raw_residual = achieved - target
        max_abs = float(np.max(np.abs(raw_residual)))

        # The source construction also requires E>0.  Check the numerical
        # repaired pair on the integration nodes rather than assuming a small
        # coefficient is automatically harmless.
        x, _ = self._quadrature()
        repaired = self.corrected_scaled_profiles(
            x, eta, P_star, f_eta, coefficients
        )
        if np.any(repaired["E"] <= 0.0):
            raise RuntimeError("five-moment solve violated E positivity")
        success = bool(result.success and max_abs <= tolerance)
        if not success:
            raise RuntimeError(
                f"bounded five-moment solve did not close: max residual {max_abs:.3e}"
            )
        return FiveMomentSolveResult(
            coefficients=tuple(float(v) for v in coefficients),
            target_scaled_moments=tuple(float(v) for v in target),
            achieved_scaled_moments=tuple(float(v) for v in achieved),
            residual_scaled_moments=tuple(float(v) for v in raw_residual),
            max_abs_residual=max_abs,
            nfev=int(result.nfev),
            success=True,
        )

    def to_payload(self) -> dict[str, Any]:
        payload = {
            "schema": SCHEMA,
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "corrected_release": CORRECTED_RELEASE,
                "corrected_release_date": CORRECTED_RELEASE_DATE,
            },
            "source_formulas": dict(_SOURCE_FORMULAS),
            "autonomous_choices": {
                "U_log_centers": list(_U_LOG_CENTERS),
                "E_log_centers": list(_E_LOG_CENTERS),
                "log_half_width": _LOG_HALF_WIDTH,
                "bump": "exp(1-1/(1-s^2)) for |s|<1, else 0",
                "finite_coefficient_bound": self.coefficient_limit,
                "solver": "bounded scipy.optimize.least_squares with analytic Jacobian",
            },
            "parameters": {
                "quadrature_points": self.quadrature_points,
                "coefficient_limit": self.coefficient_limit,
                "scaled_interval": [X_SCALED_MIN, X_SCALED_MAX],
            },
            "truth_boundary": dict(_TRUTH_BOUNDARY),
        }
        payload["sha256"] = hashlib.sha256(
            _canonical_json(payload).encode("utf-8")
        ).hexdigest()
        return payload

    @property
    def sha256(self) -> str:
        return str(self.to_payload()["sha256"])

    @classmethod
    def from_payload(cls, payload: dict[str, Any]) -> "KokunoFiveMomentRepair":
        if not isinstance(payload, dict):
            raise ValueError("payload must be a mapping")
        claimed = payload.get("sha256")
        unsigned = {key: value for key, value in payload.items() if key != "sha256"}
        expected = hashlib.sha256(_canonical_json(unsigned).encode("utf-8")).hexdigest()
        if claimed != expected:
            raise ValueError("payload SHA mismatch")
        expected_source = {
            "repository": SOURCE_REPOSITORY,
            "commit": SOURCE_COMMIT,
            "path": SOURCE_PATH,
            "corrected_release": CORRECTED_RELEASE,
            "corrected_release_date": CORRECTED_RELEASE_DATE,
        }
        if payload.get("schema") != SCHEMA or payload.get("source") != expected_source:
            raise ValueError("source/schema metadata mismatch")
        if payload.get("source_formulas") != _SOURCE_FORMULAS:
            raise ValueError("source formula metadata mismatch")
        if payload.get("truth_boundary") != _TRUTH_BOUNDARY:
            raise ValueError("truth-boundary metadata mismatch")
        choices = payload.get("autonomous_choices", {})
        if (
            choices.get("U_log_centers") != list(_U_LOG_CENTERS)
            or choices.get("E_log_centers") != list(_E_LOG_CENTERS)
            or choices.get("log_half_width") != _LOG_HALF_WIDTH
        ):
            raise ValueError("autonomous bump metadata mismatch")
        parameters = payload.get("parameters", {})
        if parameters.get("scaled_interval") != [X_SCALED_MIN, X_SCALED_MAX]:
            raise ValueError("scaled interval metadata mismatch")
        candidate = cls(
            quadrature_points=parameters.get("quadrature_points"),
            coefficient_limit=parameters.get("coefficient_limit"),
        )
        if candidate.to_payload() != payload:
            raise ValueError("payload does not replay the exact repair primitive")
        return candidate

    def save_json(self, path: str | Path) -> Path:
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(
            json.dumps(self.to_payload(), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        return target

    @classmethod
    def load_json(cls, path: str | Path) -> "KokunoFiveMomentRepair":
        payload = json.loads(Path(path).read_text(encoding="utf-8"))
        return cls.from_payload(payload)
