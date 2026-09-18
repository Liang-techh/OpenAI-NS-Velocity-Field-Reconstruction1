"""Executable PA.17 five-moment shear repair on Kokuno reserved interval I1.

The corrected 2026-09-09 Kokuno reconstruction uses two different five-row
repairs. ``kokuno_five_moment_repair`` implements the earlier PA.16 inner join
on ``U=4 eta, E=P_* f x^(1/10)``. After the admissible shear loop is inserted
through the source high-frequency radial modulation, the resulting O(1/N)
five-moment discrepancy is repaired on a compact subinterval of the first
reserved power-law patch I1. There the source background is

    U = 0,  E = K(eta) X^(-1/2-lambda),

and PA.17 gives the exact linear blocks, in row order ``(M,J ; I,S,C_p)``,

    U columns: (1, sqrt(2) K X^(-lambda)),
    E columns: (sqrt(2X), -K X^(-1/2-lambda), K X^(-3/2-lambda)).

This module implements that later local inverse rather than reusing PA.16. For
numerical conditioning it chooses an autonomous interior scale X_1 in I1,
sets ``xi=X/X_1`` and ``e_1=c_patch X_1^(-1/2-lambda)``, and scales both
perturbations by e_1. With ``b=f xi^(-1/2-lambda)`` the exact normalized
increments are

    dM /(X_1 e_1)         = int u dxi,
    dJ /(X_1^(3/2)e_1^2)  = int sqrt(2xi) u (b+e) dxi,
    dI /(X_1^(3/2)e_1)    = int sqrt(2xi) e dxi,
    dS /(X_1 e_1^2)       = int (u^2-b e-e^2/2) dxi,
    dCp/e_1^2             = int (b e+e^2/2)/xi dxi.

Their zero-correction Jacobian is exactly PA.17 after removal of nonzero
physical scale factors. The source only requires ordered nonnegative compact
bumps; the concrete bump locations, I1 centre, quadrature and finite coefficient
bound below are repository-autonomous numerical choices.

The class also exposes the repair-only profile/velocity correction on I1. A
caller must supply both c(eta) and dc/deta so the incompressibility-derived
radial correction is not silently evaluated with a false eta derivative. The
actual modulation discrepancy and its smooth coefficient family remain future
inputs. This local source-map replay is not independent Navier--Stokes
validation and does not make the global Kokuno field paper-exact.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import least_squares

from .kokuno_outer_reserved_patch_schedule import KokunoOuterReservedPatchSchedule
from .kokuno_similarity_coordinates import KokunoNativeSimilarityCoordinates


SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-pa17-shear-moment-repair-v1"

_I1_LOG_OFFSET = -22.5
_XI_MIN = 0.80
_XI_MAX = 1.70
_U_CENTERS = (0.94, 1.43)
_E_CENTERS = (0.87, 1.18, 1.57)
_BUMP_HALF_WIDTH = 0.065

_SOURCE_FORMULAS = {
    "radial_modulation": "E_N=E exp(A(X,eta,N log X)/N); U_N=U+B(X,eta,N log X)/N",
    "discrepancy_size": "the five prefix-moment differences before I1 are O_m(N^-1)",
    "I1_background": "U=0; E=K(eta) X^(-1/2-lambda) on a compact subinterval of I1",
    "PA17_rows": "row order (M,J;I,S,C_p)",
    "PA17_U_block": "(1,sqrt(2) K X^(-lambda))",
    "PA17_E_block": "(sqrt(2X),-K X^(-1/2-lambda),K X^(-3/2-lambda))",
    "exact_nonlinearity": "M,I linear; J bilinear; S,C_p quadratic; retain u*e,u^2,e^2 terms",
    "endpoint_map": "exact five-moment restoration plus compact profile support preserves all later moments, pressure, velocity and stress",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "PA17_later_shear_repair_executable": True,
    "PA16_inner_join_reused_as_PA17": False,
    "source_exact_nonlinear_increment_rows_implemented": True,
    "source_I1_power_background_implemented": True,
    "native_velocity_correction_executable": True,
    "autonomous_I1_repair_center": True,
    "autonomous_bump_shapes_and_locations": True,
    "autonomous_quadrature": True,
    "autonomous_coefficient_bound": True,
    "actual_radial_modulation_discrepancy_supplied": False,
    "eta_smooth_repair_coefficient_family_reconstructed": False,
    "cone_modulation_completed": False,
    "I1_cone_repair_applied_to_actual_modulation": False,
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


def _coefficient_array(value: Any, limit: float, name: str = "coefficients") -> np.ndarray:
    coefficients = _finite_array(value, name)
    if coefficients.shape != (5,):
        raise ValueError(f"{name} must have shape (5,)")
    if name == "coefficients" and np.any(np.abs(coefficients) > limit + 1.0e-15):
        raise ValueError("repair coefficient exceeds the declared bound")
    return coefficients


def _cinfty_bump(x: np.ndarray, center: float) -> np.ndarray:
    values = np.asarray(x, dtype=float)
    s = (values - center) / _BUMP_HALF_WIDTH
    out = np.zeros_like(values)
    active = np.abs(s) < 1.0
    if np.any(active):
        sa = s[active]
        out[active] = np.exp(1.0 - 1.0 / (1.0 - sa * sa))
    return out


def _cinfty_bump_dx(x: np.ndarray, center: float) -> np.ndarray:
    values = np.asarray(x, dtype=float)
    s = (values - center) / _BUMP_HALF_WIDTH
    out = np.zeros_like(values)
    active = np.abs(s) < 1.0
    if np.any(active):
        sa = s[active]
        bump = np.exp(1.0 - 1.0 / (1.0 - sa * sa))
        out[active] = bump * (-2.0 * sa / (1.0 - sa * sa) ** 2) / _BUMP_HALF_WIDTH
    return out


def _safe_exp(log_value: float, name: str) -> float:
    maximum = math.log(np.finfo(float).max)
    minimum = math.log(np.nextafter(0.0, 1.0))
    if not minimum <= log_value <= maximum:
        raise OverflowError(f"{name} is outside float64 exponent range; use log_scale_report()")
    value = math.exp(log_value)
    if not math.isfinite(value) or value <= 0.0:
        raise OverflowError(f"{name} could not be materialized as a positive finite float")
    return value


@dataclass(frozen=True)
class ShearMomentSolveResult:
    coefficients: tuple[float, float, float, float, float]
    target_normalized: tuple[float, float, float, float, float]
    achieved_normalized: tuple[float, float, float, float, float]
    residual_normalized: tuple[float, float, float, float, float]
    max_abs_residual: float
    nfev: int
    success: bool


@dataclass(frozen=True)
class KokunoShearMomentRepair:
    """Numerical realization of the corrected reader's PA.17 local inverse."""

    outer_schedule: KokunoOuterReservedPatchSchedule = field(default_factory=KokunoOuterReservedPatchSchedule)
    quadrature_points: int = 192
    coefficient_limit: float = 0.05
    i1_log_offset: float = _I1_LOG_OFFSET

    _nodes: np.ndarray = field(init=False, repr=False, compare=False)
    _weights: np.ndarray = field(init=False, repr=False, compare=False)
    _coordinates: KokunoNativeSimilarityCoordinates = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        if not isinstance(self.outer_schedule, KokunoOuterReservedPatchSchedule):
            raise TypeError("outer_schedule must be a KokunoOuterReservedPatchSchedule")
        if isinstance(self.quadrature_points, bool) or not isinstance(self.quadrature_points, (int, np.integer)):
            raise TypeError("quadrature_points must be an integer")
        order = int(self.quadrature_points)
        if not 64 <= order <= 256:
            raise ValueError("quadrature_points must lie in [64,256]")
        limit = float(self.coefficient_limit)
        if not math.isfinite(limit) or not (0.0 < limit <= 0.2):
            raise ValueError("coefficient_limit must lie in (0,0.2]")
        offset = float(self.i1_log_offset)
        if not math.isfinite(offset) or not (-25.0 < offset < -20.0):
            raise ValueError("i1_log_offset must lie strictly inside source I1 offsets (-25,-20)")
        if offset + math.log(_XI_MIN) <= -25.0:
            raise ValueError("autonomous repair support crosses the lower I1 edge")
        if offset + math.log(_XI_MAX) >= -20.0:
            raise ValueError("autonomous repair support crosses the upper I1 edge")
        nodes, weights = leggauss(order)
        nodes = np.asarray(nodes, dtype=float)
        weights = np.asarray(weights, dtype=float)
        nodes.setflags(write=False)
        weights.setflags(write=False)
        object.__setattr__(self, "quadrature_points", order)
        object.__setattr__(self, "coefficient_limit", limit)
        object.__setattr__(self, "i1_log_offset", offset)
        object.__setattr__(self, "_nodes", nodes)
        object.__setattr__(self, "_weights", weights)
        object.__setattr__(self, "_coordinates", KokunoNativeSimilarityCoordinates(h=float(self.outer_schedule.h)))

    @property
    def lambda_outer(self) -> float:
        return float(self.outer_schedule.lambda_outer)

    @property
    def h(self) -> float:
        return float(self.outer_schedule.h)

    @property
    def log_X_1(self) -> float:
        return float(self.outer_schedule.log_X_w + self.outer_schedule.T_w + self.i1_log_offset)

    @property
    def log_e_1(self) -> float:
        return float(self.outer_schedule.log_c_patch - (0.5 + self.lambda_outer) * self.log_X_1)

    @property
    def X_1(self) -> float:
        return _safe_exp(self.log_X_1, "X_1")

    @property
    def e_1(self) -> float:
        return _safe_exp(self.log_e_1, "e_1")

    @staticmethod
    def source_f(eta: Any) -> np.ndarray:
        values = _finite_array(eta, "eta")
        if np.any(np.abs(values) > 1.0):
            raise ValueError("eta must lie in [-1,1]")
        return 1.0 / (1.0 + values * values)

    def _quadrature(self) -> tuple[np.ndarray, np.ndarray]:
        half = 0.5 * (_XI_MAX - _XI_MIN)
        centre = 0.5 * (_XI_MAX + _XI_MIN)
        xi = centre + half * self._nodes
        return np.asarray(xi), np.asarray(half * self._weights)

    @staticmethod
    def basis_values(xi: Any) -> np.ndarray:
        values = _finite_array(xi, "xi")
        return np.stack([*(_cinfty_bump(values, c) for c in _U_CENTERS), *(_cinfty_bump(values, c) for c in _E_CENTERS)], axis=0)

    @staticmethod
    def basis_derivatives(xi: Any) -> np.ndarray:
        values = _finite_array(xi, "xi")
        return np.stack([*(_cinfty_bump_dx(values, c) for c in _U_CENTERS), *(_cinfty_bump_dx(values, c) for c in _E_CENTERS)], axis=0)

    def _basis_primitives(self, xi: np.ndarray) -> np.ndarray:
        values = np.asarray(xi, dtype=float)
        flat = values.reshape(-1)
        out = np.zeros((5, flat.size), dtype=float)
        for j, upper in enumerate(flat):
            if upper <= _XI_MIN:
                continue
            end = min(float(upper), _XI_MAX)
            half = 0.5 * (end - _XI_MIN)
            centre = 0.5 * (end + _XI_MIN)
            points = centre + half * self._nodes
            out[:, j] = half * (self.basis_values(points) @ self._weights)
        return out.reshape((5,) + values.shape)

    def background_over_e1(self, xi: Any, f_eta: float) -> np.ndarray:
        values = _finite_array(xi, "xi")
        if np.any(values <= 0.0):
            raise ValueError("xi must be positive")
        f = float(f_eta)
        if not math.isfinite(f) or f <= 0.0:
            raise ValueError("f_eta must be positive and finite")
        return f * np.power(values, -0.5 - self.lambda_outer)

    def normalized_increment(self, coefficients: Any, *, f_eta: float) -> np.ndarray:
        coeff = _coefficient_array(coefficients, self.coefficient_limit)
        xi, weights = self._quadrature()
        basis = self.basis_values(xi)
        u = np.tensordot(coeff[:2], basis[:2], axes=(0, 0))
        e = np.tensordot(coeff[2:], basis[2:], axes=(0, 0))
        b = self.background_over_e1(xi, f_eta)
        root = np.sqrt(2.0 * xi)
        integrands = np.stack([u, root * u * (b + e), root * e, u * u - b * e - 0.5 * e * e, (b * e + 0.5 * e * e) / xi], axis=0)
        return np.asarray(integrands @ weights, dtype=float)

    def coefficient_jacobian(self, coefficients: Any, *, f_eta: float) -> np.ndarray:
        coeff = _coefficient_array(coefficients, self.coefficient_limit)
        xi, weights = self._quadrature()
        basis = self.basis_values(xi)
        u = np.tensordot(coeff[:2], basis[:2], axes=(0, 0))
        e = np.tensordot(coeff[2:], basis[2:], axes=(0, 0))
        b = self.background_over_e1(xi, f_eta)
        root = np.sqrt(2.0 * xi)
        jac = np.zeros((5, 5), dtype=float)
        for column, beta in enumerate(basis):
            if column < 2:
                rows = np.stack([beta, root * beta * (b + e), np.zeros_like(beta), 2.0 * u * beta, np.zeros_like(beta)], axis=0)
            else:
                rows = np.stack([np.zeros_like(beta), root * u * beta, root * beta, -(b + e) * beta, (b + e) * beta / xi], axis=0)
            jac[:, column] = rows @ weights
        return jac

    def linear_weights(self, xi: Any, *, f_eta: float) -> dict[str, np.ndarray]:
        values = _finite_array(xi, "xi")
        if np.any(values <= 0.0):
            raise ValueError("xi must be positive")
        f = float(f_eta)
        if not math.isfinite(f) or f <= 0.0:
            raise ValueError("f_eta must be positive and finite")
        return {"U": np.stack([np.ones_like(values), np.sqrt(2.0) * f * np.power(values, -self.lambda_outer)], axis=0), "E": np.stack([np.sqrt(2.0 * values), -f * np.power(values, -0.5 - self.lambda_outer), f * np.power(values, -1.5 - self.lambda_outer)], axis=0)}

    def solve(self, target_normalized: Any, *, f_eta: float) -> ShearMomentSolveResult:
        target = _finite_array(target_normalized, "target_normalized")
        if target.shape != (5,):
            raise ValueError("target_normalized must have shape (5,)")
        f = float(f_eta)
        if not math.isfinite(f) or f <= 0.0:
            raise ValueError("f_eta must be positive and finite")
        zero = np.zeros(5, dtype=float)
        jac0 = self.coefficient_jacobian(zero, f_eta=f)
        row_scale = np.linalg.norm(jac0, axis=1)
        singular = np.linalg.svd(jac0, compute_uv=False)
        if np.any(~np.isfinite(row_scale)) or np.any(row_scale <= 0.0) or singular[-1] <= 1.0e-12 * singular[0]:
            raise RuntimeError("PA.17 normalized Jacobian is numerically rank deficient")

        def residual(c: np.ndarray) -> np.ndarray:
            return (self.normalized_increment(c, f_eta=f) - target) / row_scale

        def jacobian(c: np.ndarray) -> np.ndarray:
            return self.coefficient_jacobian(c, f_eta=f) / row_scale[:, None]

        result = least_squares(residual, zero, jac=jacobian, bounds=(-self.coefficient_limit, self.coefficient_limit), method="trf", x_scale="jac", ftol=1.0e-13, xtol=1.0e-13, gtol=1.0e-13, max_nfev=200)
        coefficients = np.asarray(result.x, dtype=float)
        achieved = self.normalized_increment(coefficients, f_eta=f)
        raw_residual = achieved - target
        return ShearMomentSolveResult(tuple(float(v) for v in coefficients), tuple(float(v) for v in target), tuple(float(v) for v in achieved), tuple(float(v) for v in raw_residual), float(np.max(np.abs(raw_residual))), int(result.nfev), bool(result.success))

    def log_scale_report(self) -> dict[str, Any]:
        i1 = self.outer_schedule.reserved_log_intervals()["I1"]
        scale_logs = {"M": self.log_X_1 + self.log_e_1, "J": 1.5 * self.log_X_1 + 2.0 * self.log_e_1, "I": 1.5 * self.log_X_1 + self.log_e_1, "S": self.log_X_1 + 2.0 * self.log_e_1, "C_p": 2.0 * self.log_e_1}
        return {"I1_log_interval": list(i1), "i1_log_offset": self.i1_log_offset, "xi_support": [_XI_MIN, _XI_MAX], "log_X_1": self.log_X_1, "log_e_1": self.log_e_1, "X_1": self.X_1, "e_1": self.e_1, "physical_moment_scale_logs": scale_logs}

    def physical_increment(self, normalized: Any) -> np.ndarray:
        values = _finite_array(normalized, "normalized")
        if values.shape != (5,):
            raise ValueError("normalized must have shape (5,)")
        logs = self.log_scale_report()["physical_moment_scale_logs"]
        scales = np.asarray([_safe_exp(float(logs[key]), f"{key} moment scale") for key in ("M", "J", "I", "S", "C_p")], dtype=float)
        return values * scales

    def profile_correction_logX(self, log_X: Any, eta: Any, *, coefficients: Any, coefficient_eta: Any) -> dict[str, np.ndarray]:
        log_array, eta_array = np.broadcast_arrays(_finite_array(log_X, "log_X"), _finite_array(eta, "eta"))
        if np.any(np.abs(eta_array) > 1.0):
            raise ValueError("eta must lie in [-1,1]")
        low, high = self.outer_schedule.reserved_log_intervals()["I1"]
        if np.any(log_array <= low) or np.any(log_array >= high):
            raise ValueError("PA.17 repair-only profile evaluation is restricted to open I1")
        coeff = _coefficient_array(coefficients, self.coefficient_limit)
        coeff_eta = _coefficient_array(coefficient_eta, float("inf"), name="coefficient_eta")
        xi = np.exp(log_array - self.log_X_1)
        basis = self.basis_values(xi)
        basis_dx = self.basis_derivatives(xi)
        primitives = self._basis_primitives(xi)
        u = np.tensordot(coeff[:2], basis[:2], axes=(0, 0))
        e = np.tensordot(coeff[2:], basis[2:], axes=(0, 0))
        u_xi = np.tensordot(coeff[:2], basis_dx[:2], axes=(0, 0))
        e_xi = np.tensordot(coeff[2:], basis_dx[2:], axes=(0, 0))
        u_eta = np.tensordot(coeff_eta[:2], basis[:2], axes=(0, 0))
        e_eta = np.tensordot(coeff_eta[2:], basis[2:], axes=(0, 0))
        m = np.tensordot(coeff[:2], primitives[:2], axes=(0, 0))
        m_eta = np.tensordot(coeff_eta[:2], primitives[:2], axes=(0, 0))
        X = np.exp(log_array)
        e1 = self.e_1
        X1 = self.X_1
        delta_U = e1 * u
        delta_E = e1 * e
        delta_U_X = e1 * u_xi / X1
        delta_E_X = e1 * e_xi / X1
        delta_U_eta = e1 * u_eta
        delta_E_eta = e1 * e_eta
        delta_M = X1 * e1 * m
        delta_M_eta = X1 * e1 * m_eta
        d = 1.0 - eta_array * eta_array
        D = 0.5 - self.h
        L = 1.0 - 2.0 * self.h * eta_array * eta_array
        delta_V0 = (2.0 * eta_array * X * delta_U - 2.0 * D * eta_array * delta_M - d * delta_M_eta) / L
        delta_v0 = delta_V0 / X
        root = np.sqrt(2.0 * X)
        delta_F = delta_E / root
        # Algebraically identical to delta_E/(2*X*root), but avoids forming
        # X*sqrt(2X), which overflows on the source existence-scale I1.
        delta_F_X = delta_E_X / root - 0.5 * delta_F / X
        delta_F_eta = delta_E_eta / root
        return {"xi": xi, "delta_E": delta_E, "delta_E_X": delta_E_X, "delta_E_eta": delta_E_eta, "delta_F": delta_F, "delta_F_X": delta_F_X, "delta_F_eta": delta_F_eta, "delta_U": delta_U, "delta_U_X": delta_U_X, "delta_U_eta": delta_U_eta, "delta_M": delta_M, "delta_M_eta": delta_M_eta, "delta_v0": delta_v0}

    def velocity_correction(self, x: Any, y: Any, z: Any, t: Any, *, coefficients: Any, coefficient_eta: Any) -> np.ndarray:
        coordinates = self._coordinates.evaluate(x, y, z, t)
        profiles = self.profile_correction_logX(np.log(coordinates["X"]), coordinates["eta"], coefficients=coefficients, coefficient_eta=coefficient_eta)
        x_array, y_array = np.broadcast_arrays(_finite_array(x, "x"), _finite_array(y, "y"))
        x_array = np.broadcast_to(x_array, np.shape(profiles["delta_F"]))
        y_array = np.broadcast_to(y_array, np.shape(profiles["delta_F"]))
        q = np.asarray(coordinates["q"], dtype=float)
        swirl = np.power(q, -1.0 - self.h) * profiles["delta_F"]
        radial = profiles["delta_v0"] / (2.0 * q)
        axial = np.power(q, -0.5 - self.h) * profiles["delta_U"]
        return np.stack((radial * x_array - swirl * y_array, radial * y_array + swirl * x_array, axial), axis=-1)

    def _unsigned_payload(self) -> dict[str, Any]:
        return {"schema": SCHEMA, "source": {"repository": SOURCE_REPOSITORY, "commit": SOURCE_COMMIT, "path": SOURCE_PATH, "corrected_release": CORRECTED_RELEASE, "corrected_release_date": CORRECTED_RELEASE_DATE}, "source_formulas": dict(_SOURCE_FORMULAS), "outer_schedule": self.outer_schedule.to_payload(), "quadrature_points": self.quadrature_points, "coefficient_limit": self.coefficient_limit, "i1_log_offset": self.i1_log_offset, "autonomous_numerics": {"xi_support": [_XI_MIN, _XI_MAX], "u_centers": list(_U_CENTERS), "e_centers": list(_E_CENTERS), "bump_half_width": _BUMP_HALF_WIDTH, "quadrature": "fixed Gauss-Legendre on autonomous xi support", "both_delta_U_and_delta_E_scaled_by": "e_1"}, "truth_boundary": dict(_TRUTH_BOUNDARY)}

    @property
    def sha256(self) -> str:
        return hashlib.sha256(_canonical_json(self._unsigned_payload()).encode("utf-8")).hexdigest()

    def to_payload(self) -> dict[str, Any]:
        payload = self._unsigned_payload()
        payload["sha256"] = self.sha256
        return payload

    @classmethod
    def from_payload(cls, payload: dict[str, Any]) -> "KokunoShearMomentRepair":
        if not isinstance(payload, dict) or payload.get("schema") != SCHEMA:
            raise ValueError("unexpected Kokuno PA.17 shear-repair schema")
        schedule_payload = payload.get("outer_schedule")
        if not isinstance(schedule_payload, dict):
            raise ValueError("outer_schedule is missing")
        obj = cls(outer_schedule=KokunoOuterReservedPatchSchedule.from_payload(schedule_payload), quadrature_points=int(payload.get("quadrature_points")), coefficient_limit=float(payload.get("coefficient_limit")), i1_log_offset=float(payload.get("i1_log_offset")))
        if obj.to_payload() != payload:
            raise ValueError("PA.17 shear-repair payload hash or content mismatch")
        return obj

    def save_json(self, path: str | Path) -> Path:
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(self.to_payload(), indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return target

    @classmethod
    def load_json(cls, path: str | Path) -> "KokunoShearMomentRepair":
        return cls.from_payload(json.loads(Path(path).read_text(encoding="utf-8")))
