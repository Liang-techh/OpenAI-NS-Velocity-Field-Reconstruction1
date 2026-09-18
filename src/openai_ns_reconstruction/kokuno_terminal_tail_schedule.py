"""Executable terminal release/tail scale binding for the Kokuno leading profile.

This module implements the *public* post-pulse scale chain in the corrected
2026-09-09 reconstruction.  It is intentionally stacked on
``KokunoOuterReservedPatchSchedule``: the earlier ``X_R/X_w`` and reserved
patch scales remain the same object, while this module continues the schedule

    X_p = X_w exp(T_w),
    pulse length = 13/lambda,
    interpolation length = T_f,
    uniform hold length = 30 log(1/lambda),
    X_rel = endpoint of that hold,

then executes the source exterior release

    l: -lambda -> -1  (one unit),
    l = -1            (4 log(1/h) units),
    l: -1 -> -h       (one unit),
    l = -h            (until Q_s reaches Q_p),

before the terminal interval ``0 <= y=log(X/X_tail) <= 3``.  The source gives

    Q_s(X_rel) = (lambda-h)/(1-lambda),
    f_o(y)=1-rho_o psi_o(y), rho_o=c_o h,
    l=-h+f_o'/f_o,
    Q_p=int_0^3 exp(int_0^v(1+l)) f_o'/f_o dv,

and the hold length is ``log(Q_in/Q_p)/(1-h)``.  Continuity of the terminal
power law ``E=c_inf X^(-A) f_o`` then determines ``c_inf``.

The reconstruction does not publish one hidden numerical ``T_f`` or ``c_o``.
They remain explicit autonomous choices here and are guarded by the source
inequalities.  In particular this object does *not* claim paper-exact hidden
parameters.  All large radial scales are retained in log space.  The default
diagnostic schedule has ``X_tail`` beyond IEEE-754 float range; materialization
therefore fails closed instead of silently overflowing.

This binds the terminal-tail geometry and amplitude needed by the corrected
heat-replacement discrepancy.  It does not apply the three-bump compensation,
does not assemble the global leading velocity, and is not a PDE validation.
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
from scipy.integrate import solve_ivp

from .kokuno_outer_reserved_patch_schedule import KokunoOuterReservedPatchSchedule


SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-terminal-tail-schedule-v1"

_SOURCE_FORMULAS = {
    "post_reserved_chain": (
        "X_p=X_w*exp(T_w); pulse length=13/lambda; after a source-free fixed "
        "T_f interpolation retain the eta-independent exponential for "
        "30*log(1/lambda); its endpoint is X_rel"
    ),
    "release_q_initial": "Q_s(X_rel)=(lambda-h)/(1-lambda)",
    "release_slopes": (
        "l:-lambda->-1 on one unit; l=-1 for 4*log(1/h); "
        "l:-1->-h on one unit; then l=-h until Q_s=Q_p"
    ),
    "terminal_multiplier": (
        "f_o(y)=1-rho_o*psi_o(y), psi_o=1-sigma((y-1)/2), "
        "rho_o=c_o*h, 0<=f_o'/f_o<h/4 on 0<=y<=3"
    ),
    "terminal_q_target": (
        "Q_p=int_0^3 exp(int_0^v(1+l(s))ds)*(f_o'(v)/f_o(v))dv, "
        "l=-h+f_o'/f_o"
    ),
    "hold_length": "L_hold=log(Q_in/Q_p)/(1-h)",
    "terminal_scale": "X_tail=X_hold_start*exp(L_hold)",
    "terminal_power": "E=c_inf*X^(-A)*f_o(y), A=1/2+h",
    "terminal_continuity": (
        "c_inf=E_tail_start*X_tail^A/f_o(0); this is a continuity identity, "
        "not a recovered hidden parameter"
    ),
    "heat_geometry": "X_K=exp(.2)*X_tail; X_b=exp(3)*X_tail",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "source_post_pulse_scale_chain_executable": True,
    "source_release_Q_equation_executable": True,
    "source_X_tail_relation_executable": True,
    "source_c_inf_continuity_relation_executable": True,
    "log_space_large_scale_tracking": True,
    "autonomous_T_f": True,
    "autonomous_c_o": True,
    "source_hidden_numeric_choices_recovered": False,
    "actual_heat_discrepancy_applied_to_patch": False,
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


def _exp_if_float(log_value: float) -> float | None:
    value = float(log_value)
    if not math.isfinite(value):
        return None
    max_log = math.log(np.finfo(float).max)
    min_log = math.log(np.nextafter(0.0, 1.0))
    if value > max_log or value < min_log:
        return None
    result = math.exp(value)
    return result if math.isfinite(result) and result > 0.0 else None


def _flat_step(value: Any) -> np.ndarray:
    s = _finite_array(value, "step_argument")
    out = np.zeros_like(s)
    out[s >= 1.0] = 1.0
    active = (s > 0.0) & (s < 1.0)
    if np.any(active):
        a = s[active]
        left = np.exp(-1.0 / (a * a))
        right = np.exp(-1.0 / ((1.0 - a) ** 2))
        out[active] = left / (left + right)
    return out


def _flat_step_derivative(value: Any) -> np.ndarray:
    s = _finite_array(value, "step_argument")
    out = np.zeros_like(s)
    active = (s > 0.0) & (s < 1.0)
    if np.any(active):
        a = s[active]
        left = np.exp(-1.0 / (a * a))
        right = np.exp(-1.0 / ((1.0 - a) ** 2))
        sigma = left / (left + right)
        out[active] = sigma * (1.0 - sigma) * (
            2.0 / (a**3) + 2.0 / ((1.0 - a) ** 3)
        )
    return out


@dataclass(frozen=True)
class KokunoTerminalTailSchedule:
    """Continue a source-aligned outer schedule through the terminal power tail.

    ``T_f`` and ``c_o`` are explicit autonomous finite choices because the
    public existence proof only requires a sufficiently large fixed ``T_f``
    and a fixed positive ``c_o`` satisfying ``f_o'/f_o<h/4``.
    """

    outer_schedule: KokunoOuterReservedPatchSchedule = field(
        default_factory=KokunoOuterReservedPatchSchedule
    )
    T_f: float = 60.0
    c_o: float = 0.05
    terminal_quadrature_order: int = 192
    transition_rtol: float = 2.0e-12
    transition_atol: float = 2.0e-14

    def __post_init__(self) -> None:
        if not isinstance(self.outer_schedule, KokunoOuterReservedPatchSchedule):
            raise TypeError("outer_schedule must be a KokunoOuterReservedPatchSchedule")
        T_f = float(self.T_f)
        c_o = float(self.c_o)
        rtol = float(self.transition_rtol)
        atol = float(self.transition_atol)
        if not math.isfinite(T_f) or not (1.0 <= T_f <= 1.0e4):
            raise ValueError("T_f must be finite and lie in [1,1e4]")
        if not math.isfinite(c_o) or not (0.0 < c_o < 0.2):
            raise ValueError("c_o must be finite and lie in (0,0.2)")
        if isinstance(self.terminal_quadrature_order, bool) or not isinstance(
            self.terminal_quadrature_order, (int, np.integer)
        ):
            raise TypeError("terminal_quadrature_order must be an integer")
        order = int(self.terminal_quadrature_order)
        if not 64 <= order <= 512:
            raise ValueError("terminal_quadrature_order must lie in [64,512]")
        if not math.isfinite(rtol) or not (0.0 < rtol <= 1.0e-8):
            raise ValueError("transition_rtol must lie in (0,1e-8]")
        if not math.isfinite(atol) or not (0.0 < atol <= 1.0e-10):
            raise ValueError("transition_atol must lie in (0,1e-10]")
        if not self.h < self.lambda_outer:
            raise ValueError("terminal schedule requires h<lambda_outer")

        object.__setattr__(self, "T_f", T_f)
        object.__setattr__(self, "c_o", c_o)
        object.__setattr__(self, "terminal_quadrature_order", order)
        object.__setattr__(self, "transition_rtol", rtol)
        object.__setattr__(self, "transition_atol", atol)

        lo, hi = self.interpolation_l_bounds()
        if lo < -self.lambda_outer - 0.1 - 2.0e-8 or hi > -self.lambda_outer + 2.0e-8:
            raise ValueError("T_f does not satisfy the source interpolation slope band")
        if self.max_terminal_log_slope() >= 0.25 * self.h:
            raise ValueError("c_o violates the source terminal bound f_o'/f_o<h/4")
        if not self.Q_in > self.Q_p > 0.0:
            raise ValueError("release schedule requires Q_in>Q_p>0")

    @property
    def h(self) -> float:
        return self.outer_schedule.h

    @property
    def lambda_outer(self) -> float:
        return self.outer_schedule.lambda_outer

    @property
    def A(self) -> float:
        return 0.5 + self.h

    @property
    def rho_o(self) -> float:
        return self.c_o * self.h

    @property
    def pulse_length(self) -> float:
        return 13.0 / self.lambda_outer

    @property
    def uniform_hold_length(self) -> float:
        return 30.0 * math.log(1.0 / self.lambda_outer)

    @property
    def steep_hold_length(self) -> float:
        return 4.0 * math.log(1.0 / self.h)

    @property
    def Q_release(self) -> float:
        return (self.lambda_outer - self.h) / (1.0 - self.lambda_outer)

    @staticmethod
    def sigma(value: Any) -> np.ndarray:
        return _flat_step(value)

    @staticmethod
    def sigma_derivative(value: Any) -> np.ndarray:
        return _flat_step_derivative(value)

    def interpolation_l_bounds(self) -> tuple[float, float]:
        """Screen the source interpolation band on its compact parameter box."""

        s = np.linspace(0.0, 1.0, 8193)
        theta_prime = -self.sigma_derivative(s) / self.T_f
        l_min = self.lambda_outer * -1.0 + float(np.min(theta_prime)) * math.log(2.0)
        l_max = -self.lambda_outer
        return l_min, l_max

    def terminal_multiplier(self, y: Any) -> np.ndarray:
        values = _finite_array(y, "terminal_y")
        if np.any(values < 0.0) or np.any(values > 3.0):
            raise ValueError("terminal_y must lie in [0,3]")
        psi = 1.0 - self.sigma((values - 1.0) / 2.0)
        return 1.0 - self.rho_o * psi

    def terminal_multiplier_derivative(self, y: Any) -> np.ndarray:
        values = _finite_array(y, "terminal_y")
        if np.any(values < 0.0) or np.any(values > 3.0):
            raise ValueError("terminal_y must lie in [0,3]")
        return 0.5 * self.rho_o * self.sigma_derivative((values - 1.0) / 2.0)

    def terminal_log_slope(self, y: Any) -> np.ndarray:
        f = self.terminal_multiplier(y)
        return self.terminal_multiplier_derivative(y) / f

    def max_terminal_log_slope(self) -> float:
        y = np.linspace(0.0, 3.0, 16385)
        return float(np.max(self.terminal_log_slope(y)))

    def release_l_first(self, y: float) -> float:
        value = float(y)
        if not 0.0 <= value <= 1.0:
            raise ValueError("release transition coordinate must lie in [0,1]")
        return -self.lambda_outer - (1.0 - self.lambda_outer) * float(
            self.sigma(np.asarray(value))
        )

    def release_l_second(self, y: float) -> float:
        value = float(y)
        if not 0.0 <= value <= 1.0:
            raise ValueError("release transition coordinate must lie in [0,1]")
        return -1.0 + (1.0 - self.h) * float(self.sigma(np.asarray(value)))

    def _evolve_Q(self, initial: float, l_law: Any) -> float:
        def rhs(y: float, q: np.ndarray) -> list[float]:
            l_value = float(l_law(y))
            return [(-l_value - self.h) - (1.0 + l_value) * float(q[0])]

        result = solve_ivp(
            rhs,
            (0.0, 1.0),
            [float(initial)],
            method="DOP853",
            rtol=self.transition_rtol,
            atol=self.transition_atol,
        )
        if not result.success:
            raise RuntimeError("release transition ODE solve failed")
        value = float(result.y[0, -1])
        if not math.isfinite(value) or value <= 0.0:
            raise RuntimeError("release transition produced a nonpositive Q_s")
        return value

    @property
    def Q_after_first_transition(self) -> float:
        return self._evolve_Q(self.Q_release, self.release_l_first)

    @property
    def Q_after_steep_hold(self) -> float:
        return self.Q_after_first_transition + (1.0 - self.h) * self.steep_hold_length

    @property
    def Q_in(self) -> float:
        return self._evolve_Q(self.Q_after_steep_hold, self.release_l_second)

    @property
    def Q_p(self) -> float:
        nodes, weights = leggauss(self.terminal_quadrature_order)
        y = 1.5 * (nodes + 1.0)
        weights = 1.5 * weights
        f0 = float(self.terminal_multiplier(np.asarray(0.0)))
        fp = self.terminal_multiplier_derivative(y)
        value = np.sum(weights * np.exp((1.0 - self.h) * y) * fp / f0)
        return float(value)

    @property
    def release_hold_length(self) -> float:
        value = math.log(self.Q_in / self.Q_p) / (1.0 - self.h)
        if not math.isfinite(value) or value <= 0.0:
            raise RuntimeError("source l=-h release hold must have positive finite length")
        return value

    def terminal_Q(self, y: Any) -> np.ndarray:
        """Evaluate the exact terminal Q_s representation on 0<=y<=3."""

        values = _finite_array(y, "terminal_y")
        if np.any(values < 0.0) or np.any(values > 3.0):
            raise ValueError("terminal_y must lie in [0,3]")
        flat = values.reshape(-1)
        nodes, weights = leggauss(self.terminal_quadrature_order)
        out = np.empty_like(flat)
        for index, y0 in enumerate(flat):
            if y0 == 3.0:
                out[index] = 0.0
                continue
            half = 0.5 * (3.0 - float(y0))
            mid = 0.5 * (3.0 + float(y0))
            v = mid + half * nodes
            f0 = float(self.terminal_multiplier(np.asarray(float(y0))))
            fp = self.terminal_multiplier_derivative(v)
            out[index] = half * np.sum(
                weights * np.exp((1.0 - self.h) * (v - float(y0))) * fp / f0
            )
        return out.reshape(values.shape)

    @property
    def log_X_p(self) -> float:
        return self.outer_schedule.log_X_w + self.outer_schedule.T_w

    @property
    def log_e_b(self) -> float:
        return self.outer_schedule.log_e_w - (0.5 + self.lambda_outer) * self.outer_schedule.T_w

    @property
    def log_e_end(self) -> float:
        return self.log_e_b - (0.5 + self.lambda_outer) * self.pulse_length

    @property
    def log_X_rel(self) -> float:
        return self.log_X_p + self.pulse_length + self.T_f + self.uniform_hold_length

    @property
    def log_e_rel(self) -> float:
        return (
            self.log_e_end
            - (0.5 + self.lambda_outer) * self.T_f
            - math.log(2.0)
            - (0.5 + self.lambda_outer) * self.uniform_hold_length
        )

    @property
    def log_X_hold_start(self) -> float:
        return self.log_X_rel + 2.0 + self.steep_hold_length

    @property
    def log_e_hold_start(self) -> float:
        return (
            self.log_e_rel
            - 1.0
            - 0.5 * self.lambda_outer
            - 1.5 * self.steep_hold_length
            - 1.0
            - 0.5 * self.h
        )

    @property
    def log_X_tail(self) -> float:
        return self.log_X_hold_start + self.release_hold_length

    @property
    def log_e_tail_start(self) -> float:
        return self.log_e_hold_start - self.A * self.release_hold_length

    @property
    def log_c_inf(self) -> float:
        f0 = float(self.terminal_multiplier(np.asarray(0.0)))
        return self.log_e_tail_start + self.A * self.log_X_tail - math.log(f0)

    @property
    def log_X_K(self) -> float:
        return self.log_X_tail + 0.2

    @property
    def log_X_b(self) -> float:
        return self.log_X_tail + 3.0

    @property
    def log_e_K(self) -> float:
        f_k = float(self.terminal_multiplier(np.asarray(0.2)))
        return self.log_c_inf - self.A * self.log_X_K + math.log(f_k)

    def log_ratio_report(self) -> dict[str, float | None]:
        log_xstar_over_xk = self.outer_schedule.log_X_star - self.log_X_K
        log_estar_over_ek = self.outer_schedule.log_e_star - self.log_e_K
        return {
            "log_X_star_over_X_K": log_xstar_over_xk,
            "log_e_star_over_e_K": log_estar_over_ek,
            "X_star_over_X_K": _exp_if_float(log_xstar_over_xk),
            "e_star_over_e_K": _exp_if_float(log_estar_over_ek),
            "log_nominal_X_K_inverse": -self.log_X_K,
            "nominal_X_K_inverse": _exp_if_float(-self.log_X_K),
        }

    def materialize_heat_inputs(self) -> dict[str, float]:
        """Return V2 heat-discrepancy inputs when all absolute scales fit floats."""

        X_tail = _exp_if_float(self.log_X_tail)
        c_inf = _exp_if_float(self.log_c_inf)
        if X_tail is None or c_inf is None:
            raise OverflowError(
                "terminal tail scales exceed finite float range; use log_scale_report()"
            )
        return {"h": self.h, "X_tail": X_tail, "c_inf": c_inf, "rho_o": self.rho_o}

    def log_scale_report(self) -> dict[str, Any]:
        l_lo, l_hi = self.interpolation_l_bounds()
        return {
            "log_X_p": self.log_X_p,
            "log_e_b": self.log_e_b,
            "pulse_length": self.pulse_length,
            "T_f": self.T_f,
            "uniform_hold_length": self.uniform_hold_length,
            "log_X_rel": self.log_X_rel,
            "log_e_rel": self.log_e_rel,
            "Q_release": self.Q_release,
            "Q_after_first_transition": self.Q_after_first_transition,
            "steep_hold_length": self.steep_hold_length,
            "Q_after_steep_hold": self.Q_after_steep_hold,
            "Q_in": self.Q_in,
            "Q_p": self.Q_p,
            "release_hold_length": self.release_hold_length,
            "rho_o": self.rho_o,
            "terminal_log_slope_max": self.max_terminal_log_slope(),
            "interpolation_l_bounds": [l_lo, l_hi],
            "log_X_tail": self.log_X_tail,
            "log_e_tail_start": self.log_e_tail_start,
            "log_c_inf": self.log_c_inf,
            "log_X_K": self.log_X_K,
            "log_e_K": self.log_e_K,
            "log_X_b": self.log_X_b,
            "X_tail": _exp_if_float(self.log_X_tail),
            "c_inf": _exp_if_float(self.log_c_inf),
            "X_K": _exp_if_float(self.log_X_K),
            "e_K": _exp_if_float(self.log_e_K),
            "X_b": _exp_if_float(self.log_X_b),
            "repair_ratios": self.log_ratio_report(),
        }

    def _outer_parameters(self) -> dict[str, float]:
        outer = self.outer_schedule
        return {
            "M_d": outer.M_d,
            "log_p_star_margin": outer.log_p_star_margin,
            "C": outer.C,
            "lambda_outer": outer.lambda_outer,
            "h": outer.h,
            "repair_log_offset": outer.repair_log_offset,
        }

    def _unsigned_payload(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "provenance": {
                "source_repository": SOURCE_REPOSITORY,
                "source_commit": SOURCE_COMMIT,
                "source_path": SOURCE_PATH,
                "corrected_release": CORRECTED_RELEASE,
                "corrected_release_date": CORRECTED_RELEASE_DATE,
            },
            "parameters": {
                "outer_schedule": self._outer_parameters(),
                "T_f": self.T_f,
                "c_o": self.c_o,
                "terminal_quadrature_order": self.terminal_quadrature_order,
                "transition_rtol": self.transition_rtol,
                "transition_atol": self.transition_atol,
            },
            "source_formulas": dict(_SOURCE_FORMULAS),
            "truth_boundary": dict(_TRUTH_BOUNDARY),
        }

    @property
    def sha256(self) -> str:
        return hashlib.sha256(_canonical_json(self._unsigned_payload()).encode("utf-8")).hexdigest()

    def to_payload(self) -> dict[str, Any]:
        payload = self._unsigned_payload()
        payload["sha256"] = self.sha256
        return payload

    @classmethod
    def from_payload(cls, payload: dict[str, Any]) -> "KokunoTerminalTailSchedule":
        if not isinstance(payload, dict) or payload.get("schema") != SCHEMA:
            raise ValueError("unexpected Kokuno terminal-tail schedule schema")
        params = payload.get("parameters")
        if not isinstance(params, dict):
            raise ValueError("parameters are missing")
        outer_params = params.get("outer_schedule")
        if not isinstance(outer_params, dict):
            raise ValueError("outer_schedule parameters are missing")
        obj = cls(
            outer_schedule=KokunoOuterReservedPatchSchedule(**outer_params),
            T_f=params.get("T_f"),
            c_o=params.get("c_o"),
            terminal_quadrature_order=params.get("terminal_quadrature_order"),
            transition_rtol=params.get("transition_rtol"),
            transition_atol=params.get("transition_atol"),
        )
        if obj.to_payload() != payload:
            raise ValueError("terminal-tail schedule payload hash or content mismatch")
        return obj

    def save_json(self, path: str | Path) -> None:
        Path(path).write_text(
            json.dumps(self.to_payload(), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )

    @classmethod
    def load_json(cls, path: str | Path) -> "KokunoTerminalTailSchedule":
        payload = json.loads(Path(path).read_text(encoding="utf-8"))
        return cls.from_payload(payload)
