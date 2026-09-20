"""Executable five-moment repair mechanics for the Kokuno final-join annulus.

Pinned public provenance: KokunoYumeto/yang-mills-interacting-workbench at
143f6773feb424ad9ed3a8d116653200f20346b7, corrected 2026-09-09 reader.

The corrected reader states that on

    -6 < log x < -5,       x = X/X_R,

one chooses two fixed ordered nonnegative bumps for U and three for E.  For a
reference pair (U,E), with u=delta U and e=delta E, it gives the complete
five moment-increment integrands (PA.14), in order (M,I,J,S,C_p),

    u,
    sqrt(2X)e,
    H u + U sqrt(2X)e + sqrt(2X)u e,
    2Uu + u^2 - Ee - e^2/2,
    Ee/X + e^2/(2X).

After X=X_R x and the public moment normalization, the same exact increment
map is evaluated here in x coordinates.

The reader fixes the bump count/order/support class and proves local
invertibility, but the material audited by Agent 1 does not specify numerical
bump shapes.  The compact C-infinity bump intervals below are therefore a
repository-autonomous realization, chosen before any residual calculation.
They are not OpenAI/paper-exact coefficients.  This module also does not know
the real upstream five-moment discrepancy, so it exposes a local solver for a
*supplied small discrepancy* and a manufactured replay only.  It is not the
missing global inner-to-outer bridge and is not PDE validation.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass, field
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

from .kokuno_pa10_ideal_join_profile_contract import (
    KokunoPA10IdealJoinProfileContract,
)
from .kokuno_pa10_physical_center_profile_contract import (
    CORRECTED_RELEASE,
    CORRECTED_RELEASE_DATE,
    SOURCE_COMMIT,
    SOURCE_PATH,
    SOURCE_REPOSITORY,
)


SCHEMA = "kokuno-pa10-five-moment-repair-system-v1"
LOG_X_SUPPORT = (-6.0, -5.0)
QUADRATURE_ORDER = 256
NEWTON_MAX_ITERATIONS = 10
NEWTON_ABS_TOL = 2.0e-13
COEFFICIENT_TRUST_RADIUS = 5.0e-5

# Fixed before any PDE/residual calculation. Each tuple is an open interval
# in log x. The final support ends before -5, leaving an autonomous
# ideal-profile neighborhood at the public join x=e^-5.
AUTONOMOUS_U_BUMP_INTERVALS = (
    (-5.94, -5.68),
    (-5.48, -5.22),
)
AUTONOMOUS_E_BUMP_INTERVALS = (
    (-5.92, -5.74),
    (-5.61, -5.43),
    (-5.30, -5.12),
)

_SOURCE_FORMULAS = {
    "repair_annulus": "-6<log x<-5, x=X/X_R",
    "bump_count": "two fixed ordered nonnegative U bumps plus three E bumps",
    "delta_M": "integral u dX",
    "delta_I": "integral sqrt(2X)e dX",
    "delta_J": "integral [H u + U sqrt(2X)e + sqrt(2X)u e] dX",
    "delta_S": "integral [2Uu+u^2-Ee-e^2/2] dX",
    "delta_Cp": "integral [Ee/X+e^2/(2X)] dX",
    "normalized_scaling": (
        "(M,I,J,S,C_p)="
        "(X_R hatM,X_R^(3/2) hatI,X_R^(3/2) hatJ,X_R hatS,hatC_p)"
    ),
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "public_pa14_exact_moment_increment_map_executable": True,
    "public_two_u_three_e_support_class_respected": True,
    "repository_autonomous_bump_shapes_materialized": True,
    "local_small_discrepancy_solver_executable": True,
    "source_bump_shapes_numerically_identified": False,
    "source_repair_coefficients_identified": False,
    "actual_upstream_five_moment_discrepancy_materialized": False,
    "eta_dependent_repair_coefficient_functions_materialized": False,
    "actual_inner_to_outer_bridge_materialized": False,
    "source_join_neighborhood_width_known": False,
    "outer_global_leading_velocity_materialized": False,
    "source_center_is_final_corrected_fixed_point": False,
    "matched_global_pressure_materialized": False,
    "restricted_forcing_materialized": False,
    "heldout_ns_residual_assessed": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _canonical_json(payload: Mapping[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _finite(value: Any, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if np.any(~np.isfinite(out)):
        raise ValueError(f"{name} must contain only finite values")
    return out


def _check_eta_scalar(eta: Any) -> float:
    eta_arr = _finite(eta, "eta")
    if eta_arr.ndim != 0:
        raise ValueError("eta must be scalar for a five-coefficient local solve")
    value = float(eta_arr)
    if abs(value) > 1.0:
        raise ValueError("eta must lie in the public profile interval [-1,1]")
    return value


def _check_coefficients(coefficients: Sequence[float]) -> np.ndarray:
    coeff = _finite(coefficients, "coefficients")
    if coeff.shape != (5,):
        raise ValueError("coefficients must have shape (5,)")
    if float(np.max(np.abs(coeff))) > COEFFICIENT_TRUST_RADIUS:
        raise ValueError(
            "repair coefficients exceed the preregistered local coefficient trust radius "
            f"{COEFFICIENT_TRUST_RADIUS:g}"
        )
    return coeff


def _raw_cinf_bump(log_x: np.ndarray, interval: tuple[float, float]) -> np.ndarray:
    a, b = interval
    z = (2.0 * log_x - (a + b)) / (b - a)
    out = np.zeros_like(log_x, dtype=float)
    mask = np.abs(z) < 1.0
    z_in = z[mask]
    out[mask] = np.exp(-1.0 / (1.0 - z_in * z_in))
    return out


@dataclass(frozen=True)
class KokunoPA10FiveMomentRepairSystem:
    """Exact PA.14 moment mechanics with one autonomous compact bump basis."""

    ideal: KokunoPA10IdealJoinProfileContract = field(
        default_factory=KokunoPA10IdealJoinProfileContract,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        if not isinstance(self.ideal, KokunoPA10IdealJoinProfileContract):
            raise TypeError("ideal must be KokunoPA10IdealJoinProfileContract")

    @property
    def P_star(self) -> float:
        return float(self.ideal.P_star)

    @property
    def X_R(self) -> float:
        return float(self.ideal.X_R)

    @property
    def X_h(self) -> float:
        return float(self.ideal.X_h)

    @staticmethod
    def _quadrature(order: int = QUADRATURE_ORDER) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        if not isinstance(order, int) or order < 32:
            raise ValueError("quadrature order must be an integer >=32")
        nodes, weights = np.polynomial.legendre.leggauss(order)
        y0, y1 = LOG_X_SUPPORT
        log_x = 0.5 * (y1 - y0) * nodes + 0.5 * (y0 + y1)
        dy_weight = 0.5 * (y1 - y0) * weights
        x = np.exp(log_x)
        # dx = x d(log x)
        dx_weight = dy_weight * x
        return log_x, x, dx_weight

    @classmethod
    def _normalizers(
        cls, intervals: Sequence[tuple[float, float]], order: int = QUADRATURE_ORDER
    ) -> np.ndarray:
        log_x, _, dx_weight = cls._quadrature(order)
        values = []
        for interval in intervals:
            integral = float(np.sum(dx_weight * _raw_cinf_bump(log_x, interval)))
            if not math.isfinite(integral) or integral <= 0.0:
                raise RuntimeError("autonomous bump normalization failed")
            values.append(integral)
        return np.asarray(values, dtype=float)

    @classmethod
    def _basis_family(
        cls,
        x: Any,
        intervals: Sequence[tuple[float, float]],
        *,
        order: int = QUADRATURE_ORDER,
    ) -> np.ndarray:
        x_arr = _finite(x, "x")
        if np.any(x_arr <= 0.0):
            raise ValueError("normalized x must be strictly positive")
        log_x = np.log(x_arr)
        normalizers = cls._normalizers(intervals, order=order)
        pieces = [
            _raw_cinf_bump(log_x, interval) / normalizer
            for interval, normalizer in zip(intervals, normalizers, strict=True)
        ]
        return np.stack(pieces, axis=-1)

    def basis_values(self, x: Any) -> dict[str, np.ndarray]:
        return {
            "U_bumps": self._basis_family(x, AUTONOMOUS_U_BUMP_INTERVALS),
            "E_bumps": self._basis_family(x, AUTONOMOUS_E_BUMP_INTERVALS),
        }

    def quadrature_snapshot(self, *, order: int = QUADRATURE_ORDER) -> dict[str, np.ndarray]:
        log_x, x, dx_weight = self._quadrature(order)
        return {
            "log_x": log_x,
            "x": x,
            "dx_weight": dx_weight,
            "U_bumps": self._basis_family(x, AUTONOMOUS_U_BUMP_INTERVALS),
            "E_bumps": self._basis_family(x, AUTONOMOUS_E_BUMP_INTERVALS),
        }

    def _ideal_x_values(self, x: np.ndarray, eta: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        f = 1.0 / (1.0 + eta * eta)
        U0 = np.full_like(x, 4.0 * eta)
        E0 = self.P_star * f * np.power(x, 0.1)
        H0_normalized = np.sqrt(2.0 * x) * E0
        return U0, E0, H0_normalized

    def perturbations(
        self, x: Any, coefficients: Sequence[float]
    ) -> dict[str, np.ndarray]:
        coeff = _check_coefficients(coefficients)
        basis = self.basis_values(x)
        u = np.einsum("...j,j->...", basis["U_bumps"], coeff[:2])
        e = np.einsum("...j,j->...", basis["E_bumps"], coeff[2:])
        return {"delta_U": u, "delta_E": e}

    def repaired_values(
        self, x: Any, eta: Any, coefficients: Sequence[float]
    ) -> dict[str, np.ndarray]:
        coeff = _check_coefficients(coefficients)
        x_arr, eta_arr = np.broadcast_arrays(_finite(x, "x"), _finite(eta, "eta"))
        if np.any(x_arr <= 0.0):
            raise ValueError("normalized x must be strictly positive")
        if np.any(np.abs(eta_arr) > 1.0):
            raise ValueError("eta must lie in the public profile interval [-1,1]")
        log_x = np.log(x_arr)
        lo, hi = LOG_X_SUPPORT
        if np.any(log_x < lo) or np.any(log_x > hi):
            raise ValueError(
                "repaired_values is intentionally restricted to the public repair annulus "
                "-6<=log x<=-5; this module is not a global outer profile"
            )
        perturb = self.perturbations(x_arr, coeff)
        f = 1.0 / (1.0 + eta_arr * eta_arr)
        U0 = 4.0 * eta_arr
        E0 = self.P_star * f * np.power(x_arr, 0.1)
        return {
            "x": x_arr,
            "eta": eta_arr,
            "U_ideal": U0,
            "E_ideal": E0,
            "delta_U": perturb["delta_U"],
            "delta_E": perturb["delta_E"],
            "U_repaired": U0 + perturb["delta_U"],
            "E_repaired": E0 + perturb["delta_E"],
        }

    def moment_increment(
        self,
        coefficients: Sequence[float],
        eta: Any,
        *,
        order: int = QUADRATURE_ORDER,
    ) -> np.ndarray:
        coeff = _check_coefficients(coefficients)
        eta_value = _check_eta_scalar(eta)
        q = self.quadrature_snapshot(order=order)
        x = q["x"]
        w = q["dx_weight"]
        U_bumps = q["U_bumps"]
        E_bumps = q["E_bumps"]
        u = U_bumps @ coeff[:2]
        e = E_bumps @ coeff[2:]
        U0, E0, H0 = self._ideal_x_values(x, eta_value)
        root = np.sqrt(2.0 * x)
        return np.asarray(
            [
                np.sum(w * u),
                np.sum(w * root * e),
                np.sum(w * (H0 * u + U0 * root * e + root * u * e)),
                np.sum(w * (2.0 * U0 * u + u * u - E0 * e - 0.5 * e * e)),
                np.sum(w * (E0 * e / x + 0.5 * e * e / x)),
            ],
            dtype=float,
        )

    def coefficient_jacobian(
        self,
        coefficients: Sequence[float],
        eta: Any,
        *,
        order: int = QUADRATURE_ORDER,
    ) -> np.ndarray:
        coeff = _check_coefficients(coefficients)
        eta_value = _check_eta_scalar(eta)
        q = self.quadrature_snapshot(order=order)
        x = q["x"]
        w = q["dx_weight"]
        U_bumps = q["U_bumps"]
        E_bumps = q["E_bumps"]
        u = U_bumps @ coeff[:2]
        e = E_bumps @ coeff[2:]
        U0, E0, H0 = self._ideal_x_values(x, eta_value)
        root = np.sqrt(2.0 * x)
        jac = np.zeros((5, 5), dtype=float)

        for j in range(2):
            du = U_bumps[:, j]
            jac[:, j] = [
                np.sum(w * du),
                0.0,
                np.sum(w * (H0 * du + root * du * e)),
                np.sum(w * (2.0 * U0 + 2.0 * u) * du),
                0.0,
            ]

        for k in range(3):
            de = E_bumps[:, k]
            j = 2 + k
            jac[:, j] = [
                0.0,
                np.sum(w * root * de),
                np.sum(w * (U0 * root + root * u) * de),
                np.sum(w * (-E0 - e) * de),
                np.sum(w * (E0 / x + e / x) * de),
            ]
        return jac

    def linearized_system(self, eta: Any) -> dict[str, Any]:
        eta_value = _check_eta_scalar(eta)
        jac = self.coefficient_jacobian(np.zeros(5), eta_value)
        singular = np.linalg.svd(jac, compute_uv=False)
        rank = int(np.linalg.matrix_rank(jac))
        return {
            "eta": eta_value,
            "jacobian": jac,
            "rank": rank,
            "determinant": float(np.linalg.det(jac)),
            "condition_number_2": float(singular[0] / singular[-1]),
            "singular_values": singular,
        }

    def solve_moment_delta(self, target_delta: Sequence[float], eta: Any) -> dict[str, Any]:
        target = _finite(target_delta, "target_delta")
        if target.shape != (5,):
            raise ValueError("target_delta must have shape (5,)")
        eta_value = _check_eta_scalar(eta)

        jac0 = self.coefficient_jacobian(np.zeros(5), eta_value)
        if np.linalg.matrix_rank(jac0) != 5:
            raise RuntimeError("autonomous five-bump linearization is not full rank")
        coeff = np.linalg.solve(jac0, target)
        if float(np.max(np.abs(coeff))) > COEFFICIENT_TRUST_RADIUS:
            raise ValueError("linearized repair exits the preregistered local coefficient trust region")

        converged = False
        iterations = 0
        for iteration in range(NEWTON_MAX_ITERATIONS + 1):
            residual = self.moment_increment(coeff, eta_value) - target
            residual_inf = float(np.max(np.abs(residual)))
            if residual_inf <= NEWTON_ABS_TOL:
                converged = True
                iterations = iteration
                break
            if iteration == NEWTON_MAX_ITERATIONS:
                break
            jac = self.coefficient_jacobian(coeff, eta_value)
            step = np.linalg.solve(jac, residual)

            # Deterministic fail-closed backtracking within the local trust region.
            accepted = False
            base_norm = residual_inf
            damping = 1.0
            for _ in range(12):
                proposal = coeff - damping * step
                if float(np.max(np.abs(proposal))) <= COEFFICIENT_TRUST_RADIUS:
                    proposal_residual = self.moment_increment(proposal, eta_value) - target
                    if float(np.max(np.abs(proposal_residual))) < base_norm:
                        coeff = proposal
                        accepted = True
                        break
                damping *= 0.5
            if not accepted:
                raise RuntimeError("Newton repair could not decrease the exact PA.14 moment residual")

        final_residual = self.moment_increment(coeff, eta_value) - target
        final_inf = float(np.max(np.abs(final_residual)))
        if not converged:
            raise RuntimeError(
                "five-moment repair did not converge within the preregistered iteration budget; "
                f"residual_inf={final_inf:.3e}"
            )
        return {
            "coefficients": coeff,
            "eta": eta_value,
            "target_delta": target,
            "achieved_delta": self.moment_increment(coeff, eta_value),
            "residual": final_residual,
            "residual_inf": final_inf,
            "iterations": iterations,
            "converged": True,
        }

    def configuration(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "ideal_configuration": self.ideal.configuration(),
            "log_x_support": list(LOG_X_SUPPORT),
            "quadrature_order": QUADRATURE_ORDER,
            "newton_max_iterations": NEWTON_MAX_ITERATIONS,
            "newton_abs_tol": NEWTON_ABS_TOL,
            "coefficient_trust_radius": COEFFICIENT_TRUST_RADIUS,
            "autonomous_U_bump_intervals": [list(v) for v in AUTONOMOUS_U_BUMP_INTERVALS],
            "autonomous_E_bump_intervals": [list(v) for v in AUTONOMOUS_E_BUMP_INTERVALS],
            "bump_semantics": "repository-autonomous-cinf-unit-dx-integral-not-source-identified",
        }

    @classmethod
    def from_configuration(cls, payload: Mapping[str, Any]) -> "KokunoPA10FiveMomentRepairSystem":
        if not isinstance(payload, Mapping):
            raise TypeError("configuration must be a mapping")
        expected = cls().configuration()
        if payload.get("schema") != SCHEMA:
            raise ValueError("configuration schema mismatch")
        for key in (
            "log_x_support",
            "quadrature_order",
            "newton_max_iterations",
            "newton_abs_tol",
            "coefficient_trust_radius",
            "autonomous_U_bump_intervals",
            "autonomous_E_bump_intervals",
            "bump_semantics",
        ):
            if payload.get(key) != expected[key]:
                raise ValueError(f"serialized fixed repair setting mismatch: {key}")
        ideal = KokunoPA10IdealJoinProfileContract.from_configuration(payload["ideal_configuration"])
        return cls(ideal=ideal)

    def save_configuration(self, path: str | Path) -> dict[str, Any]:
        payload = self.configuration()
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        return payload

    @classmethod
    def load_configuration(cls, path: str | Path) -> "KokunoPA10FiveMomentRepairSystem":
        return cls.from_configuration(json.loads(Path(path).read_text()))

    @property
    def semantic_sha256(self) -> str:
        payload = {
            "schema": SCHEMA,
            "source_commit": SOURCE_COMMIT,
            "source_formulas": _SOURCE_FORMULAS,
            "ideal_semantic_sha256": self.ideal.semantic_sha256,
            "configuration": self.configuration(),
        }
        return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()

    @property
    def truth_boundary(self) -> dict[str, bool]:
        return copy.deepcopy(_TRUTH_BOUNDARY)

    def report(self) -> dict[str, Any]:
        linear = self.linearized_system(0.3)
        return {
            "schema": SCHEMA,
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "release": CORRECTED_RELEASE,
                "release_date": CORRECTED_RELEASE_DATE,
            },
            "source_formulas": copy.deepcopy(_SOURCE_FORMULAS),
            "configuration": self.configuration(),
            "semantic_sha256": self.semantic_sha256,
            "linearized_probe_eta_0p3": {
                "rank": linear["rank"],
                "determinant": linear["determinant"],
                "condition_number_2": linear["condition_number_2"],
                "singular_values": linear["singular_values"].tolist(),
            },
            "truth_boundary": self.truth_boundary,
        }
