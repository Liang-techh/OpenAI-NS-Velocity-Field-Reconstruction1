"""Executable later half of the Kokuno PA.16 inner-to-outer join.

Pinned provenance: KokunoYumeto/yang-mills-interacting-workbench at
143f6773feb424ad9ed3a8d116653200f20346b7, corrected 2026-09-09 reader.

The public reconstruction fixes X_i=110, X_R=110(C P_*)^10 and x=X/X_R.
Writing ell_i=log(C E(X_i,eta)), G_i=U(X_i,eta), f=(1+eta^2)^(-1), it
prescribes, for y_i=log(X/X_i),

  log E = -log C + y_i/10 + (1-sigma(y_i/T_sh)) ell_i
          + sigma(y_i/T_sh) log f,
  U = G_i.

After X_sep=110 exp(T_sh), E is exactly the ideal outer angular profile.
The source restores G_i to 4 eta on -8<log x<-7, then applies the PA.16
2-U/3-E nonlinear five-moment repair on -6<log x<-5.  At
X_h=X_R exp(-5), equality of the fields and the five prefix moments gives the
exact exterior preservation map (for the common axis pressure datum).

This file implements only that later join map.  T_sh, ell_i, G_i and the
incoming PA.15-scaled five-moment discrepancy are caller supplied because they
are outputs of earlier Appendix-B existence choices not yet executable in this
repository.  They are never inferred from plots or residuals.

A numerical subtlety matters at the repository's existence-style P_*~3e13:
the PA.16 Jacobian is mathematically full rank but its raw rows span many
orders of magnitude, so an unscaled float64 matrix_rank test falsely reports
rank loss.  The solve below uses only invertible row scaling before rank/linear
solve; the exact nonlinear PA.16 moment map and target are unchanged.
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

from .kokuno_five_moment_repair import KokunoFiveMomentRepair, FiveMomentSolveResult
from .kokuno_outer_reserved_patch_schedule import KokunoOuterReservedPatchSchedule
from .kokuno_reference_continuation import source_smooth_step

SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-inner-join-exit-v2"

X_I = 110.0
LOG_RESTORE_START = -8.0
LOG_RESTORE_END = -7.0
LOG_REPAIR_START = -6.0
LOG_JOIN_EXIT = -5.0

_SOURCE_FORMULAS = {
    "scales": "X_i=110; X_R=110*(C*P_*)^10; x=X/X_R; y_i=log(X/X_i)",
    "angular_transition": (
        "log E=-log C+y_i/10+(1-sigma(y_i/T_sh))*ell_i+"
        "sigma(y_i/T_sh)*log f; U=G_i"
    ),
    "ideal_pair": "U_0=4 eta; E_0=P_* f x^(1/10); f=(1+eta^2)^(-1)",
    "axial_restoration": "restore G_i to 4 eta with a fixed flat step on -8<log x<-7",
    "PA16_repair": "two U bumps and three E bumps on -6<log x<-5 restore M,I,J,S,C_p",
    "join_exit": "X_h=X_R*exp(-5); matched fields+five moments preserve the exterior",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "source_later_inner_join_formula_executable": True,
    "source_PA16_inverse_reused": True,
    "scale_conditioned_PA16_solve": True,
    "stable_log_x_discrepancy_propagation": True,
    "caller_supplied_T_sh": True,
    "source_T_sh_lower_bound_verified": False,
    "caller_supplied_upstream_boundary_data": True,
    "actual_upstream_appendix_B_boundary_data_bound": False,
    "actual_inner_join_target_from_full_appendix_B": False,
    "inner_to_outer_join_completed": False,
    "global_pressure_matched": False,
    "global_leading_profile_reconstructed": False,
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


def _step_derivative(s: Any) -> np.ndarray:
    values = _finite_array(s, "s")
    out = np.zeros_like(values)
    mask = (values > 0.0) & (values < 1.0)
    if np.any(mask):
        local = values[mask]
        sigma = source_smooth_step(local)
        logit_prime = 2.0 / local**3 + 2.0 / (1.0 - local) ** 3
        out[mask] = sigma * (1.0 - sigma) * logit_prime
    return out


@dataclass(frozen=True)
class InnerJoinSolveResult:
    eta: float
    incoming_scaled_discrepancy: tuple[float, float, float, float, float]
    pre_repair_scaled_discrepancy: tuple[float, float, float, float, float]
    repair: FiveMomentSolveResult
    exit_scaled_discrepancy: tuple[float, float, float, float, float]
    max_abs_exit_discrepancy: float
    scaled_jacobian_condition: float


@dataclass(frozen=True)
class KokunoInnerJoinExit:
    """Later PA.16 join map with explicit upstream inputs and fail-closed scope."""

    T_sh: float
    outer_schedule: KokunoOuterReservedPatchSchedule = field(
        default_factory=KokunoOuterReservedPatchSchedule
    )
    repair: KokunoFiveMomentRepair = field(default_factory=KokunoFiveMomentRepair)
    quadrature_points: int = 96

    _nodes: np.ndarray = field(init=False, repr=False, compare=False)
    _weights: np.ndarray = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        T_sh = float(self.T_sh)
        if not math.isfinite(T_sh) or T_sh <= 0.0:
            raise ValueError("T_sh must be positive and finite")
        if not isinstance(self.outer_schedule, KokunoOuterReservedPatchSchedule):
            raise TypeError("outer_schedule must be a KokunoOuterReservedPatchSchedule")
        if not isinstance(self.repair, KokunoFiveMomentRepair):
            raise TypeError("repair must be a KokunoFiveMomentRepair")
        if isinstance(self.quadrature_points, bool) or not isinstance(
            self.quadrature_points, (int, np.integer)
        ):
            raise TypeError("quadrature_points must be an integer")
        order = int(self.quadrature_points)
        if not 32 <= order <= 256:
            raise ValueError("quadrature_points must lie in [32,256]")
        object.__setattr__(self, "T_sh", T_sh)
        object.__setattr__(self, "quadrature_points", order)
        if not self.log_x_sep < LOG_RESTORE_START:
            raise ValueError(
                "source join requires x_sep<exp(-8); increase C/P_* or use a shorter caller-supplied T_sh"
            )
        nodes, weights = leggauss(order)
        nodes = np.asarray(nodes, dtype=float)
        weights = np.asarray(weights, dtype=float)
        nodes.setflags(write=False)
        weights.setflags(write=False)
        object.__setattr__(self, "_nodes", nodes)
        object.__setattr__(self, "_weights", weights)

    @property
    def log_x_i(self) -> float:
        return -10.0 * (
            math.log(float(self.outer_schedule.C)) + float(self.outer_schedule.log_P_star)
        )

    @property
    def log_x_sep(self) -> float:
        return self.log_x_i + self.T_sh

    @property
    def log_X_i(self) -> float:
        return math.log(X_I)

    @property
    def log_X_sep(self) -> float:
        return self.log_X_i + self.T_sh

    @property
    def log_X_h(self) -> float:
        return float(self.outer_schedule.log_X_R) + LOG_JOIN_EXIT

    @property
    def P_star(self) -> float:
        return math.exp(float(self.outer_schedule.log_P_star))

    @staticmethod
    def source_f(eta: float) -> float:
        eta = float(eta)
        if not math.isfinite(eta) or not -1.0 <= eta <= 1.0:
            raise ValueError("eta must lie in [-1,1]")
        return 1.0 / (1.0 + eta * eta)

    def geometry_report(self) -> dict[str, Any]:
        return {
            "X_i": X_I,
            "log_X_i": self.log_X_i,
            "T_sh": self.T_sh,
            "log_X_sep": self.log_X_sep,
            "log_x_i": self.log_x_i,
            "log_x_sep": self.log_x_sep,
            "axial_restore_log_x": [LOG_RESTORE_START, LOG_RESTORE_END],
            "PA16_repair_log_x": [LOG_REPAIR_START, LOG_JOIN_EXIT],
            "log_X_h": self.log_X_h,
            "source_T_sh_lower_bound_verified": False,
        }

    def _ideal(self, log_x: Any, eta: float) -> tuple[np.ndarray, np.ndarray]:
        log_values = np.asarray(log_x, dtype=float)
        f = self.source_f(eta)
        E = np.exp(float(self.outer_schedule.log_P_star) + math.log(f) + 0.1 * log_values)
        U = np.full_like(log_values, 4.0 * float(eta), dtype=float)
        return U, E

    def profile_values_scaled(
        self,
        x: Any,
        *,
        eta: float,
        ell_i: float,
        G_i: float,
        coefficients: Any | None = None,
    ) -> dict[str, np.ndarray]:
        raw = _finite_array(x, "x")
        shape = raw.shape
        values = np.atleast_1d(raw).astype(float, copy=False)
        if np.any(values <= 0.0):
            raise ValueError("x must be positive")
        eta = float(eta)
        f = self.source_f(eta)
        ell_i = float(ell_i)
        G_i = float(G_i)
        if not math.isfinite(ell_i) or not math.isfinite(G_i):
            raise ValueError("ell_i and G_i must be finite")
        log_x = np.log(values)
        if np.any(log_x < self.log_x_i - 2.0e-14) or np.any(
            log_x > LOG_JOIN_EXIT + 2.0e-14
        ):
            raise ValueError("scaled x is outside the executable inner-join interval")

        U0, E0 = self._ideal(log_x, eta)
        U = U0.copy()
        E = E0.copy()
        U_x = np.zeros_like(values)
        E_x = 0.1 * E0 / values

        mask = log_x <= self.log_x_sep
        if np.any(mask):
            y = log_x[mask] - self.log_x_i
            s = y / self.T_sh
            sigma = source_smooth_step(s)
            dsigma = _step_derivative(s)
            # Algebraically identical source formula, expressed relative to the
            # ideal pair.  This makes exact ideal boundary data replay exactly
            # instead of leaving cancellation-sized quadrature noise.
            mismatch = ell_i - math.log(f)
            log_local_E = np.log(E0[mask]) + (1.0 - sigma) * mismatch
            local_E = np.exp(log_local_E)
            D_log_E = 0.1 - dsigma * mismatch / self.T_sh
            U[mask] = G_i
            E[mask] = local_E
            U_x[mask] = 0.0
            E_x[mask] = local_E * D_log_E / values[mask]

        mask = (log_x > self.log_x_sep) & (log_x <= LOG_RESTORE_START)
        if np.any(mask):
            U[mask] = G_i

        mask = (log_x > LOG_RESTORE_START) & (log_x < LOG_RESTORE_END)
        if np.any(mask):
            s = log_x[mask] - LOG_RESTORE_START
            sigma = source_smooth_step(s)
            dsigma = _step_derivative(s)
            U[mask] = G_i + sigma * (4.0 * eta - G_i)
            U_x[mask] = dsigma * (4.0 * eta - G_i) / values[mask]

        mask = (log_x > LOG_REPAIR_START) & (log_x < LOG_JOIN_EXIT)
        if np.any(mask):
            if coefficients is None:
                raise ValueError("PA.16 coefficients are required inside the repair interval")
            repaired = self.repair.corrected_scaled_profiles(
                values[mask], eta, self.P_star, f, coefficients
            )
            U[mask] = repaired["U"]
            E[mask] = repaired["E"]
            U_x[mask] = repaired["U_x"]
            E_x[mask] = repaired["E_x"]

        result = {
            "U": U,
            "E": E,
            "H": np.sqrt(2.0 * values) * E,
            "F_scaled": E / np.sqrt(2.0 * values),
            "U_x": U_x,
            "E_x": E_x,
        }
        return {key: np.asarray(value).reshape(shape) for key, value in result.items()}

    @staticmethod
    def _moment_density_dx(
        x: np.ndarray,
        U: np.ndarray,
        E: np.ndarray,
        U0: np.ndarray,
        E0: np.ndarray,
    ) -> np.ndarray:
        root = np.sqrt(2.0 * x)
        H = root * E
        H0 = root * E0
        return np.stack(
            [
                U - U0,
                H - H0,
                U * H - U0 * H0,
                U * U - U0 * U0 - 0.5 * (E * E - E0 * E0),
                (E * E - E0 * E0) / (2.0 * x),
            ],
            axis=0,
        )

    def _integrate_log_interval(
        self,
        lo: float,
        hi: float,
        *,
        eta: float,
        ell_i: float,
        G_i: float,
    ) -> np.ndarray:
        if hi <= lo:
            return np.zeros(5, dtype=float)
        y = 0.5 * (hi - lo) * self._nodes + 0.5 * (hi + lo)
        x = np.exp(y)
        profile = self.profile_values_scaled(
            x, eta=eta, ell_i=ell_i, G_i=G_i, coefficients=None
        )
        U0, E0 = self._ideal(y, eta)
        density_dx = self._moment_density_dx(
            x, np.asarray(profile["U"]), np.asarray(profile["E"]), U0, E0
        )
        return 0.5 * (hi - lo) * ((density_dx * x[np.newaxis, :]) @ self._weights)

    def pre_repair_scaled_discrepancy(
        self,
        *,
        eta: float,
        ell_i: float,
        G_i: float,
        incoming_scaled_discrepancy: Any,
    ) -> np.ndarray:
        incoming = _finite_array(incoming_scaled_discrepancy, "incoming_scaled_discrepancy")
        if incoming.shape != (5,):
            raise ValueError("incoming_scaled_discrepancy must have shape (5,)")
        eta = float(eta)
        self.source_f(eta)
        ell_i = float(ell_i)
        G_i = float(G_i)
        if not math.isfinite(ell_i) or not math.isfinite(G_i):
            raise ValueError("ell_i and G_i must be finite")

        delta = self._integrate_log_interval(
            self.log_x_i,
            self.log_x_sep,
            eta=eta,
            ell_i=ell_i,
            G_i=G_i,
        )

        a = math.exp(self.log_x_sep)
        b = math.exp(LOG_RESTORE_START)
        dU = G_i - 4.0 * eta
        if b > a and dU != 0.0:
            delta[0] += dU * (b - a)
            h_prefactor = math.sqrt(2.0) * self.P_star * self.source_f(eta)
            delta[2] += dU * h_prefactor * (b**1.6 - a**1.6) / 1.6
            delta[3] += (G_i * G_i - (4.0 * eta) ** 2) * (b - a)

        delta += self._integrate_log_interval(
            LOG_RESTORE_START,
            LOG_RESTORE_END,
            eta=eta,
            ell_i=ell_i,
            G_i=G_i,
        )
        return incoming + delta

    @staticmethod
    def _row_matrix(eta: float) -> np.ndarray:
        matrix = np.zeros((5, 5), dtype=float)
        matrix[0, 0] = 1.0
        matrix[1, 1] = -4.0 * eta
        matrix[1, 2] = 1.0
        matrix[2, 1] = 1.0
        matrix[3, 0] = -8.0 * eta
        matrix[3, 3] = 1.0
        matrix[4, 4] = 1.0
        return matrix

    def _conditioned_pa16_solve(
        self,
        target: np.ndarray,
        *,
        eta: float,
        absolute_tolerance: float,
    ) -> tuple[FiveMomentSolveResult, float]:
        target = np.asarray(target, dtype=float)
        f = self.source_f(eta)
        zero = np.zeros(5, dtype=float)
        row_matrix = self._row_matrix(eta)
        target_t = row_matrix @ target
        jac0_t = row_matrix @ self.repair.moment_jacobian(
            zero, eta=eta, P_star=self.P_star, f_eta=f
        )
        rank_scale = np.linalg.norm(jac0_t, axis=1)
        if np.any(~np.isfinite(rank_scale)) or np.any(rank_scale <= 0.0):
            raise RuntimeError("PA.16 Jacobian has a zero/nonfinite row")
        jac_rank = jac0_t / rank_scale[:, None]
        target_rank = target_t / rank_scale
        singular = np.linalg.svd(jac_rank, compute_uv=False)
        if singular[-1] <= 1.0e-12 * singular[0]:
            raise RuntimeError("source five-moment linearized map lost rank after row scaling")
        scaled_condition = float(singular[0] / singular[-1])
        initial = np.linalg.solve(jac_rank, target_rank)
        initial = np.clip(
            initial,
            -0.9 * self.repair.coefficient_limit,
            0.9 * self.repair.coefficient_limit,
        )
        row_scale = np.maximum(
            rank_scale * self.repair.coefficient_limit,
            np.maximum(np.abs(target_t), 1.0e-14),
        )

        def residual(coeff: np.ndarray) -> np.ndarray:
            achieved = self.repair.moment_increment(
                coeff, eta=eta, P_star=self.P_star, f_eta=f
            )
            return (row_matrix @ (achieved - target)) / row_scale

        def jacobian(coeff: np.ndarray) -> np.ndarray:
            jac = self.repair.moment_jacobian(
                coeff, eta=eta, P_star=self.P_star, f_eta=f
            )
            return (row_matrix @ jac) / row_scale[:, None]

        result = least_squares(
            residual,
            initial,
            jac=jacobian,
            bounds=(-self.repair.coefficient_limit, self.repair.coefficient_limit),
            method="trf",
            ftol=1.0e-13,
            xtol=1.0e-13,
            gtol=1.0e-13,
            max_nfev=400,
        )
        coefficients = np.asarray(result.x, dtype=float)
        achieved = self.repair.moment_increment(
            coefficients, eta=eta, P_star=self.P_star, f_eta=f
        )
        raw = achieved - target
        max_abs = float(np.max(np.abs(raw)))
        sample_x = np.exp(np.linspace(LOG_REPAIR_START, LOG_JOIN_EXIT, 257))
        repaired = self.repair.corrected_scaled_profiles(
            sample_x, eta, self.P_star, f, coefficients
        )
        if np.any(np.asarray(repaired["E"]) <= 0.0):
            raise RuntimeError("PA.16 conditioned solve violated E positivity")
        if not bool(result.success) or max_abs > float(absolute_tolerance):
            raise RuntimeError(
                f"conditioned PA.16 solve did not close: max residual {max_abs:.3e}"
            )
        receipt = FiveMomentSolveResult(
            coefficients=tuple(float(v) for v in coefficients),
            target_scaled_moments=tuple(float(v) for v in target),
            achieved_scaled_moments=tuple(float(v) for v in achieved),
            residual_scaled_moments=tuple(float(v) for v in raw),
            max_abs_residual=max_abs,
            nfev=int(result.nfev),
            success=True,
        )
        return receipt, scaled_condition

    def solve_at_eta(
        self,
        *,
        eta: float,
        ell_i: float,
        G_i: float,
        incoming_scaled_discrepancy: Any,
        absolute_tolerance: float = 2.0e-11,
    ) -> InnerJoinSolveResult:
        tolerance = float(absolute_tolerance)
        if not math.isfinite(tolerance) or not 0.0 < tolerance <= 1.0e-6:
            raise ValueError("absolute_tolerance must lie in (0,1e-6]")
        incoming = _finite_array(incoming_scaled_discrepancy, "incoming_scaled_discrepancy")
        if incoming.shape != (5,):
            raise ValueError("incoming_scaled_discrepancy must have shape (5,)")
        pre = self.pre_repair_scaled_discrepancy(
            eta=eta,
            ell_i=ell_i,
            G_i=G_i,
            incoming_scaled_discrepancy=incoming,
        )
        repair_result, condition = self._conditioned_pa16_solve(
            -pre, eta=float(eta), absolute_tolerance=tolerance
        )
        exit_error = pre + np.asarray(repair_result.achieved_scaled_moments, dtype=float)
        return InnerJoinSolveResult(
            eta=float(eta),
            incoming_scaled_discrepancy=tuple(float(v) for v in incoming),
            pre_repair_scaled_discrepancy=tuple(float(v) for v in pre),
            repair=repair_result,
            exit_scaled_discrepancy=tuple(float(v) for v in exit_error),
            max_abs_exit_discrepancy=float(np.max(np.abs(exit_error))),
            scaled_jacobian_condition=condition,
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
            "parameters": {
                "T_sh": self.T_sh,
                "quadrature_points": self.quadrature_points,
                "X_i": X_I,
                "outer_schedule": self.outer_schedule.to_payload(),
                "repair": self.repair.to_payload(),
            },
            "geometry": self.geometry_report(),
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
    def from_payload(cls, payload: dict[str, Any]) -> "KokunoInnerJoinExit":
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
        parameters = payload.get("parameters", {})
        if parameters.get("X_i") != X_I:
            raise ValueError("X_i metadata mismatch")
        candidate = cls(
            T_sh=parameters.get("T_sh"),
            quadrature_points=parameters.get("quadrature_points"),
            outer_schedule=KokunoOuterReservedPatchSchedule.from_payload(
                parameters.get("outer_schedule")
            ),
            repair=KokunoFiveMomentRepair.from_payload(parameters.get("repair")),
        )
        if candidate.to_payload() != payload:
            raise ValueError("payload does not replay the exact inner-join exit")
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
    def load_json(cls, path: str | Path) -> "KokunoInnerJoinExit":
        return cls.from_payload(json.loads(Path(path).read_text(encoding="utf-8")))
