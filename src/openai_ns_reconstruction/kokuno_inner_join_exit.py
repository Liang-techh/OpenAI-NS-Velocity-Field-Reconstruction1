"""Executable later half of the corrected-reader PA.16 inner-to-outer join.

The corrected Kokuno 2026-09-09 reconstruction fixes the *form* of the
continuation from the Appendix-B regular field into the temporary RF40 outer
pair, but its left-boundary data are outputs of earlier existence choices.  In
particular, at the fixed radius ``X_i=110`` it writes

    ell_i(eta) = log(C E(X_i,eta)),     G_i(eta) = U(X_i,eta),

then, with ``x=X/X_R`` and ``X_R=110(C P_*)^10``, uses a logarithmic angular
transition of length ``T_sh``, becomes exactly

    U_0 = 4 eta,     E_0 = P_* f(eta) x^(1/10),
    f(eta)=(1+eta^2)^(-1),

after a fixed axial restoration on ``-8 < log x < -7``, and finally restores
all five prefix moments with the PA.16 two-U/three-E repair on
``-6 < log x < -5``.  At ``X_h=X_R exp(-5)`` the source matching theorem says
that equality of the two fields and all five moments propagates to the entire
outer field (for the same axis pressure datum).

This module makes that *later join map* executable.  It deliberately accepts
``ell_i``, ``G_i`` and the incoming five-moment discrepancy as caller inputs;
it does not substitute the currently available temporary frozen reference
continuation for the still-missing earlier Appendix-B continuation.  Likewise
``T_sh`` is caller supplied.  This preserves the source dependency graph
instead of inventing hidden numerical data.

The discrepancy is propagated in the exact PA.15 scaled moment order
``(M,I,J,S,C_p)``.  The very long scale separation before ``log x=-8`` is
handled in logarithmic coordinates (plus an analytic constant-U segment), so
no enormous absolute ``X`` or tiny ``x`` moment intermediate is required.
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

from .kokuno_five_moment_repair import KokunoFiveMomentRepair, FiveMomentSolveResult
from .kokuno_outer_reserved_patch_schedule import KokunoOuterReservedPatchSchedule
from .kokuno_reference_continuation import source_smooth_step


SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-inner-join-exit-v1"

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
    "axial_restoration": (
        "after X_sep=110*exp(T_sh), E is ideal; restore G_i to 4 eta with a "
        "fixed flat step on -8<log x<-7"
    ),
    "PA16_repair": "two U bumps and three E bumps on -6<log x<-5 restore M,I,J,S,C_p",
    "join_exit": (
        "X_h=X_R*exp(-5); equality of fields plus all five moments at X_h "
        "propagates to the exterior field for common Pi_ax"
    ),
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "source_later_inner_join_formula_executable": True,
    "source_PA16_inverse_reused": True,
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


@dataclass(frozen=True)
class KokunoInnerJoinExit:
    """Source-form transition/restoration/PA.16 map from ``X_i`` to ``X_h``.

    ``T_sh`` is intentionally required: the public reconstruction chooses it
    from bounds supplied by earlier continuation data, so this component may
    not silently invent that hidden dependency.  The constructor checks only
    the later geometric requirement ``x_sep < exp(-8)``.  Satisfying the
    reader's earlier lower-bound proof for ``T_sh`` remains an upstream duty.
    """

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
        if not self.log_x_sep < LOG_RESTORE_START:
            raise ValueError(
                "source join requires x_sep<exp(-8); increase C/P_* or use a shorter caller-supplied T_sh"
            )
        nodes, weights = leggauss(order)
        nodes = np.asarray(nodes, dtype=float)
        weights = np.asarray(weights, dtype=float)
        nodes.setflags(write=False)
        weights.setflags(write=False)
        object.__setattr__(self, "T_sh", T_sh)
        object.__setattr__(self, "quadrature_points", order)
        object.__setattr__(self, "_nodes", nodes)
        object.__setattr__(self, "_weights", weights)

    @property
    def log_x_i(self) -> float:
        # log(110/X_R) = -10(log C + log P_*) exactly under the source X_R law.
        return -10.0 * (math.log(float(self.outer_schedule.C)) + float(self.outer_schedule.log_P_star))

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

    def _ideal(self, log_x: np.ndarray, eta: float) -> tuple[np.ndarray, np.ndarray]:
        f = self.source_f(eta)
        E = np.exp(float(self.outer_schedule.log_P_star) + math.log(f) + 0.1 * log_x)
        U = np.full_like(log_x, 4.0 * eta, dtype=float)
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
        """Evaluate the source join profiles and radial derivatives in scaled ``x``.

        The API covers ``x_i <= x <= exp(-5)``.  PA.16 coefficients are only
        needed inside the repair interval; before it, the transition is fully
        determined by ``ell_i,G_i``.
        """

        values = _finite_array(x, "x")
        if np.any(values <= 0.0):
            raise ValueError("x must be positive")
        eta = float(eta)
        f = self.source_f(eta)
        ell_i = float(ell_i)
        G_i = float(G_i)
        if not math.isfinite(ell_i) or not math.isfinite(G_i):
            raise ValueError("ell_i and G_i must be finite")
        log_x = np.log(values)
        if np.any(log_x < self.log_x_i - 2.0e-14) or np.any(log_x > LOG_JOIN_EXIT + 2.0e-14):
            raise ValueError("scaled x is outside the executable inner-join interval")

        U0, E0 = self._ideal(log_x, eta)
        U = U0.copy()
        E = E0.copy()
        U_x = np.zeros_like(values)
        E_x = 0.1 * E0 / values

        # X_i -> X_sep: exact source logarithmic angular transition, U=G_i.
        mask = log_x <= self.log_x_sep
        if np.any(mask):
            y = log_x[mask] - self.log_x_i
            s = y / self.T_sh
            sigma = source_smooth_step(s)
            dsigma = _step_derivative(s)
            log_E = (
                -math.log(float(self.outer_schedule.C))
                + 0.1 * y
                + (1.0 - sigma) * ell_i
                + sigma * math.log(f)
            )
            local_E = np.exp(log_E)
            D_log_E = 0.1 + dsigma * (math.log(f) - ell_i) / self.T_sh
            U[mask] = G_i
            E[mask] = local_E
            U_x[mask] = 0.0
            E_x[mask] = local_E * D_log_E / values[mask]

        # X_sep -> exp(-8): angular field ideal, U still G_i.
        mask = (log_x > self.log_x_sep) & (log_x <= LOG_RESTORE_START)
        if np.any(mask):
            U[mask] = G_i

        # -8 < log x < -7: fixed flat axial restoration at ideal E.
        mask = (log_x > LOG_RESTORE_START) & (log_x < LOG_RESTORE_END)
        if np.any(mask):
            s = log_x[mask] - LOG_RESTORE_START
            sigma = source_smooth_step(s)
            dsigma = _step_derivative(s)
            U[mask] = G_i + sigma * (4.0 * eta - G_i)
            U_x[mask] = dsigma * (4.0 * eta - G_i) / values[mask]

        # -6 < log x < -5: exact PA.16 local repair on the ideal pair.
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

        H = np.sqrt(2.0 * values) * E
        F_scaled = E / np.sqrt(2.0 * values)
        return {
            "U": U,
            "E": E,
            "H": H,
            "F_scaled": F_scaled,
            "U_x": U_x,
            "E_x": E_x,
        }

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
        density_dx = self._moment_density_dx(x, profile["U"], profile["E"], U0, E0)
        # dx = x d(log x).  This avoids materializing huge physical X scales.
        density_dlogx = density_dx * x[np.newaxis, :]
        return 0.5 * (hi - lo) * (density_dlogx @ self._weights)

    def pre_repair_scaled_discrepancy(
        self,
        *,
        eta: float,
        ell_i: float,
        G_i: float,
        incoming_scaled_discrepancy: Any,
    ) -> np.ndarray:
        """Propagate the exact five-moment discrepancy to ``log x=-6``.

        ``incoming_scaled_discrepancy`` is the difference (actual minus ideal)
        of the PA.15 normalized moments at ``X_i``.  Keeping it explicit is the
        fail-closed boundary to the not-yet-executable earlier Appendix-B
        continuation.
        """

        incoming = _finite_array(incoming_scaled_discrepancy, "incoming_scaled_discrepancy")
        if incoming.shape != (5,):
            raise ValueError("incoming_scaled_discrepancy must have shape (5,)")
        eta = float(eta)
        self.source_f(eta)
        ell_i = float(ell_i)
        G_i = float(G_i)
        if not math.isfinite(ell_i) or not math.isfinite(G_i):
            raise ValueError("ell_i and G_i must be finite")

        # Source angular transition: log x_i -> log x_sep.
        delta = self._integrate_log_interval(
            self.log_x_i,
            self.log_x_sep,
            eta=eta,
            ell_i=ell_i,
            G_i=G_i,
        )

        # Long ideal-E / constant-U segment: x_sep -> exp(-8).  Integrate the
        # exact powers analytically; numerical quadrature over hundreds of
        # logarithmic units would be unnecessarily ill-conditioned.
        a = math.exp(self.log_x_sep)
        b = math.exp(LOG_RESTORE_START)
        dU = G_i - 4.0 * eta
        f = self.source_f(eta)
        if b > a and dU != 0.0:
            delta[0] += dU * (b - a)
            h_prefactor = math.sqrt(2.0) * self.P_star * f
            delta[2] += dU * h_prefactor * (b**1.6 - a**1.6) / 1.6
            delta[3] += (G_i * G_i - (4.0 * eta) ** 2) * (b - a)

        # Fixed flat U restoration on -8 < log x < -7.
        delta += self._integrate_log_interval(
            LOG_RESTORE_START,
            LOG_RESTORE_END,
            eta=eta,
            ell_i=ell_i,
            G_i=G_i,
        )
        # From -7 to -6 the pair is exactly ideal, so no further discrepancy.
        return incoming + delta

    def solve_at_eta(
        self,
        *,
        eta: float,
        ell_i: float,
        G_i: float,
        incoming_scaled_discrepancy: Any,
        absolute_tolerance: float = 2.0e-11,
    ) -> InnerJoinSolveResult:
        incoming = _finite_array(incoming_scaled_discrepancy, "incoming_scaled_discrepancy")
        if incoming.shape != (5,):
            raise ValueError("incoming_scaled_discrepancy must have shape (5,)")
        pre = self.pre_repair_scaled_discrepancy(
            eta=eta,
            ell_i=ell_i,
            G_i=G_i,
            incoming_scaled_discrepancy=incoming,
        )
        repair_result = self.repair.solve(
            -pre,
            eta=float(eta),
            P_star=self.P_star,
            f_eta=self.source_f(float(eta)),
            absolute_tolerance=absolute_tolerance,
        )
        exit_error = pre + np.asarray(repair_result.achieved_scaled_moments, dtype=float)
        return InnerJoinSolveResult(
            eta=float(eta),
            incoming_scaled_discrepancy=tuple(float(v) for v in incoming),
            pre_repair_scaled_discrepancy=tuple(float(v) for v in pre),
            repair=repair_result,
            exit_scaled_discrepancy=tuple(float(v) for v in exit_error),
            max_abs_exit_discrepancy=float(np.max(np.abs(exit_error))),
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
        payload["sha256"] = hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()
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
        if payload.get("schema") != SCHEMA:
            raise ValueError("schema mismatch")
        expected_source = {
            "repository": SOURCE_REPOSITORY,
            "commit": SOURCE_COMMIT,
            "path": SOURCE_PATH,
            "corrected_release": CORRECTED_RELEASE,
            "corrected_release_date": CORRECTED_RELEASE_DATE,
        }
        if payload.get("source") != expected_source:
            raise ValueError("source metadata mismatch")
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
        target.write_text(json.dumps(self.to_payload(), indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return target

    @classmethod
    def load_json(cls, path: str | Path) -> "KokunoInnerJoinExit":
        payload = json.loads(Path(path).read_text(encoding="utf-8"))
        return cls.from_payload(payload)
