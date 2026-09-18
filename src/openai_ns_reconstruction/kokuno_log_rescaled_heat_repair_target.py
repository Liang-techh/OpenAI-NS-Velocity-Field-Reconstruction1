"""Underflow-safe normalized heat-repair target for the Kokuno outer schedule.

This module composes the public corrected heat-replacement formulas with the
Agent-1 outer/terminal scale chain without materializing enormous absolute
``X_tail`` values.  The three target channels are returned as signed logarithms,

    -(Delta C_p/e_*^2,
      Delta S/(X_* e_*^2),
      Delta I_sub/(X_*^(3/2) e_*)),

so the corrected V2 target remains executable even when IEEE-754 would round a
nonzero channel to zero.

The transition and far-tail integrals are algebraically factored into scale
logs and O(1) dimensionless quadratures.  This is a numerical representation
change only: it uses the same public V2 heat replacement, the same autonomous
chi_K transition, and the same source X_*/e_* normalization already present in
the repository.

It deliberately does not claim that the existing float64 three-bump inverse can
resolve the resulting multi-hundred-decade channel hierarchy.  No correction is
applied to the velocity here, and no PDE validation or paper-exact claim is
made.
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
from scipy.special import roots_laguerre

from .kokuno_heat_exterior_profile import KokunoHeatExteriorProfile
from .kokuno_terminal_tail_schedule import KokunoTerminalTailSchedule


SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-log-rescaled-heat-repair-target-v1"

_SOURCE_FORMULAS = {
    "repair_target": (
        "target=-(Delta C_p/e_*^2, Delta S/(X_*e_*^2), "
        "Delta I_sub/(X_*^(3/2)e_*))"
    ),
    "heat_replacement": (
        "E_new=E_cl*(1+chi_K(y)*(H(2*(1-eta^2)/X)-1)); "
        "y=log(X/X_tail), chi_K=0 for y<=.2 and 1 for y>=.5"
    ),
    "v2_cp_tail": (
        "the y>=3 first-order Delta C_p tail carries the corrected factor 2"
    ),
    "log_factoring": (
        "X=X_tail*exp(y); factor all powers of X_tail,c_inf,X_*,e_* "
        "before quadrature and retain channel sign/log-absolute-value"
    ),
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "corrected_v2_heat_target_formula_executable": True,
    "source_outer_and_terminal_scales_bound": True,
    "signed_log_underflow_safe_target": True,
    "autonomous_chi_transition_shape": True,
    "float64_three_bump_inverse_certified": False,
    "actual_heat_discrepancy_applied_to_patch": False,
    "eta_smooth_repair_coefficients_reconstructed": False,
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


def _safe_exp(log_value: np.ndarray) -> np.ndarray:
    values = np.asarray(log_value, dtype=float)
    out = np.zeros_like(values)
    lower = math.log(np.nextafter(0.0, 1.0))
    upper = math.log(np.finfo(float).max)
    active = (values >= lower) & (values <= upper)
    out[active] = np.exp(values[active])
    out[values > upper] = np.inf
    return out


def _signed_logsum(terms: list[tuple[float, float]]) -> tuple[int, float]:
    active = [
        (float(sign), float(log_abs))
        for sign, log_abs in terms
        if sign != 0.0 and math.isfinite(log_abs)
    ]
    if not active:
        return 0, -math.inf
    reference = max(log_abs for _, log_abs in active)
    total = sum(sign * math.exp(log_abs - reference) for sign, log_abs in active)
    if total == 0.0:
        return 0, -math.inf
    return (1 if total > 0.0 else -1), reference + math.log(abs(total))


@dataclass(frozen=True)
class KokunoLogRescaledHeatRepairTarget:
    """Evaluate the corrected source heat-repair target in signed-log form."""

    terminal_schedule: KokunoTerminalTailSchedule = field(
        default_factory=KokunoTerminalTailSchedule
    )
    heat_quadrature_order: int = 96
    transition_quadrature_order: int = 192
    tail_quadrature_order: int = 96

    _heat: KokunoHeatExteriorProfile = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        if not isinstance(self.terminal_schedule, KokunoTerminalTailSchedule):
            raise TypeError("terminal_schedule must be a KokunoTerminalTailSchedule")
        for name, value, lo, hi in (
            ("heat_quadrature_order", self.heat_quadrature_order, 32, 256),
            ("transition_quadrature_order", self.transition_quadrature_order, 64, 512),
            ("tail_quadrature_order", self.tail_quadrature_order, 32, 256),
        ):
            if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
                raise TypeError(f"{name} must be an integer")
            if not lo <= int(value) <= hi:
                raise ValueError(f"{name} must lie in [{lo},{hi}]")
            object.__setattr__(self, name, int(value))
        object.__setattr__(
            self,
            "_heat",
            KokunoHeatExteriorProfile(
                h=self.h,
                c_inf=1.0,
                quadrature_order=self.heat_quadrature_order,
            ),
        )

    @property
    def h(self) -> float:
        return self.terminal_schedule.h

    @property
    def A(self) -> float:
        return 0.5 + self.h

    @property
    def outer_schedule(self):
        return self.terminal_schedule.outer_schedule

    def _heat_minus_one_over_z(self, z: Any) -> np.ndarray:
        values = _finite_array(z, "z")
        if np.any(values < 0.0):
            raise ValueError("heat argument must be nonnegative")
        out = np.empty_like(values)
        small = values <= 1.0e-4
        if np.any(~small):
            z_big = values[~small]
            out[~small] = (self._heat.heat_factor(z_big) - 1.0) / z_big
        if np.any(small):
            z_small = values[small]
            derivatives = [
                float(self._heat.heat_factor(np.asarray(0.0), derivative=k))
                for k in range(1, 5)
            ]
            out[small] = (
                derivatives[0]
                + 0.5 * derivatives[1] * z_small
                + derivatives[2] * z_small * z_small / 6.0
                + derivatives[3] * z_small * z_small * z_small / 24.0
            )
        return out

    def _z_from_log(self, log_z: Any) -> np.ndarray:
        values = _finite_array(log_z, "log_z")
        z = _safe_exp(values)
        if np.any(~np.isfinite(z)):
            raise OverflowError("heat argument overflowed unexpectedly")
        return z

    def _scalar_signed_log_target(self, eta: float) -> tuple[np.ndarray, np.ndarray]:
        value = float(eta)
        if not math.isfinite(value) or abs(value) > 1.0:
            raise ValueError("eta must lie in [-1,1]")
        d = 1.0 - value * value
        if d == 0.0:
            return np.zeros(3, dtype=int), np.full(3, -math.inf, dtype=float)

        schedule = self.terminal_schedule
        outer = self.outer_schedule
        Lx = schedule.log_X_tail
        Lc = schedule.log_c_inf
        Lxs = outer.log_X_star
        Les = outer.log_e_star

        nodes, weights = leggauss(self.transition_quadrature_order)
        y = 1.5 * (nodes + 1.0)
        weights = 1.5 * weights
        chi = _flat_step((y - 0.2) / 0.3)
        f = schedule.terminal_multiplier(y)
        log_two_d = math.log(2.0 * d)
        z = self._z_from_log(log_two_d - Lx - y)
        R = self._heat_minus_one_over_z(z)

        lag_nodes, lag_weights = roots_laguerre(self.tail_quadrature_order)
        log_z0 = log_two_d - Lx - 3.0

        def tail_R(rate: float) -> np.ndarray:
            z_tail = self._z_from_log(log_z0 - lag_nodes / rate)
            return self._heat_minus_one_over_z(z_tail)

        k_cp_1 = 2.0 * self.A + 1.0
        k_cp_2 = 2.0 * self.A + 2.0
        b_cp_1 = 2.0 * d * (
            np.sum(weights * f * f * np.exp(-k_cp_1 * y) * chi * R)
            + math.exp(-3.0 * k_cp_1)
            / k_cp_1
            * np.sum(lag_weights * tail_R(k_cp_1))
        )
        b_cp_2 = 2.0 * d * d * (
            np.sum(weights * f * f * np.exp(-k_cp_2 * y) * chi * chi * R * R)
            + math.exp(-3.0 * k_cp_2)
            / k_cp_2
            * np.sum(lag_weights * tail_R(k_cp_2) ** 2)
        )
        cp_sign, cp_log = _signed_logsum(
            [
                (
                    math.copysign(1.0, b_cp_1),
                    2.0 * Lc - k_cp_1 * Lx - 2.0 * Les + math.log(abs(b_cp_1)),
                ),
                (
                    math.copysign(1.0, b_cp_2),
                    2.0 * Lc - k_cp_2 * Lx - 2.0 * Les + math.log(abs(b_cp_2)),
                ),
            ]
        )

        k_s_1 = 1.0 + 2.0 * self.h
        k_s_2 = 2.0 + 2.0 * self.h
        b_s_1 = -2.0 * d * (
            np.sum(weights * f * f * np.exp(-k_s_1 * y) * chi * R)
            + math.exp(-3.0 * k_s_1)
            / k_s_1
            * np.sum(lag_weights * tail_R(k_s_1))
        )
        b_s_2 = -2.0 * d * d * (
            np.sum(weights * f * f * np.exp(-k_s_2 * y) * chi * chi * R * R)
            + math.exp(-3.0 * k_s_2)
            / k_s_2
            * np.sum(lag_weights * tail_R(k_s_2) ** 2)
        )
        s_sign, s_log = _signed_logsum(
            [
                (
                    math.copysign(1.0, b_s_1),
                    2.0 * Lc
                    - k_s_1 * Lx
                    - Lxs
                    - 2.0 * Les
                    + math.log(abs(b_s_1)),
                ),
                (
                    math.copysign(1.0, b_s_2),
                    2.0 * Lc
                    - k_s_2 * Lx
                    - Lxs
                    - 2.0 * Les
                    + math.log(abs(b_s_2)),
                ),
            ]
        )

        b_i = 2.0 * math.sqrt(2.0) * d * (
            np.sum(weights * f * np.exp(-self.h * y) * chi * R)
            + math.exp(-3.0 * self.h)
            / self.h
            * np.sum(lag_weights * tail_R(self.h))
        )
        i_sign = 1 if b_i > 0.0 else -1
        i_log = Lc - self.h * Lx - 1.5 * Lxs - Les + math.log(abs(b_i))

        signs = -np.asarray([cp_sign, s_sign, i_sign], dtype=int)
        logs = np.asarray([cp_log, s_log, i_log], dtype=float)
        return signs, logs

    def signed_log_target(self, eta: Any) -> dict[str, np.ndarray]:
        """Return target channel signs and natural-log absolute values."""

        values = _finite_array(eta, "eta")
        if np.any(np.abs(values) > 1.0):
            raise ValueError("eta must lie in [-1,1]")
        flat = values.reshape(-1)
        signs = np.empty((flat.size, 3), dtype=int)
        logs = np.empty((flat.size, 3), dtype=float)
        for index, entry in enumerate(flat):
            signs[index], logs[index] = self._scalar_signed_log_target(float(entry))
        return {
            "sign": signs.reshape(values.shape + (3,)),
            "log_abs": logs.reshape(values.shape + (3,)),
        }

    def materialize_target(self, eta: Any) -> np.ndarray:
        """Materialize float64 target only when every nonzero channel is representable."""

        encoded = self.signed_log_target(eta)
        signs = encoded["sign"]
        logs = encoded["log_abs"]
        lower = math.log(np.nextafter(0.0, 1.0))
        upper = math.log(np.finfo(float).max)
        nonzero = signs != 0
        if np.any(nonzero & (logs < lower)):
            raise OverflowError(
                "nonzero normalized repair target underflows float64; use signed_log_target()"
            )
        if np.any(nonzero & (logs > upper)):
            raise OverflowError(
                "normalized repair target overflows float64; use signed_log_target()"
            )
        result = np.zeros_like(logs)
        result[nonzero] = signs[nonzero] * np.exp(logs[nonzero])
        return result

    def rescaled_target(self, eta: Any) -> dict[str, np.ndarray]:
        """Return an O(1)-referenced representation without discarding channel logs."""

        encoded = self.signed_log_target(eta)
        signs = encoded["sign"]
        logs = encoded["log_abs"]
        masked = np.where(signs != 0, logs, -np.inf)
        reference = np.max(masked, axis=-1)
        finite_reference = np.isfinite(reference)
        reference_expanded = np.expand_dims(reference, axis=-1)
        active = (signs != 0) & np.expand_dims(finite_reference, axis=-1)
        shifted = np.full_like(logs, -np.inf)
        np.subtract(logs, reference_expanded, out=shifted, where=active)
        mantissa = np.zeros_like(logs)
        mantissa[active] = signs[active] * _safe_exp(shifted[active])
        return {
            "mantissa": mantissa,
            "log_reference": reference,
            "sign": signs,
            "log_abs": logs,
        }

    def precision_report(self, eta: Any) -> dict[str, Any]:
        """Report the channel hierarchy relevant to the existing float64 inverse."""

        encoded = self.signed_log_target(eta)
        signs = np.asarray(encoded["sign"])
        logs = np.asarray(encoded["log_abs"])
        if signs.shape != (3,):
            raise ValueError("precision_report requires a scalar eta")
        active = logs[signs != 0]
        if active.size == 0:
            span_decades = 0.0
        else:
            span_decades = float((np.max(active) - np.min(active)) / math.log(10.0))
        required_digits = int(math.ceil(span_decades)) + 4
        float64_digits = int(math.floor(-math.log10(np.finfo(float).eps)))
        lower = math.log(np.nextafter(0.0, 1.0))
        underflow_channels = [
            name
            for name, sign, log_abs in zip(
                ("C_p", "S", "I_sub"), signs.tolist(), logs.tolist()
            )
            if sign != 0 and log_abs < lower
        ]
        return {
            "eta": float(np.asarray(eta)),
            "sign": signs.tolist(),
            "log_abs": logs.tolist(),
            "log10_abs": [
                value / math.log(10.0) if math.isfinite(value) else -math.inf
                for value in logs.tolist()
            ],
            "span_decades": span_decades,
            "estimated_required_decimal_digits": required_digits,
            "float64_decimal_digits": float64_digits,
            "underflow_channels": underflow_channels,
            "float64_direct_repair_safe": bool(
                not underflow_channels and required_digits <= float64_digits - 2
            ),
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
            "terminal_schedule": self.terminal_schedule.to_payload(),
            "parameters": {
                "heat_quadrature_order": self.heat_quadrature_order,
                "transition_quadrature_order": self.transition_quadrature_order,
                "tail_quadrature_order": self.tail_quadrature_order,
            },
            "source_formulas": dict(_SOURCE_FORMULAS),
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
    def from_payload(cls, payload: dict[str, Any]) -> "KokunoLogRescaledHeatRepairTarget":
        if not isinstance(payload, dict) or payload.get("schema") != SCHEMA:
            raise ValueError("unexpected log-rescaled heat-repair target schema")
        if payload.get("source_formulas") != _SOURCE_FORMULAS:
            raise ValueError("source formula contract mismatch")
        if payload.get("truth_boundary") != _TRUTH_BOUNDARY:
            raise ValueError("truth boundary mismatch")
        schedule_payload = payload.get("terminal_schedule")
        params = payload.get("parameters")
        if not isinstance(schedule_payload, dict) or not isinstance(params, dict):
            raise ValueError("serialized schedule/parameters are missing")
        obj = cls(
            terminal_schedule=KokunoTerminalTailSchedule.from_payload(schedule_payload),
            **params,
        )
        if obj.to_payload() != payload:
            raise ValueError("log-rescaled heat-repair target payload hash or content mismatch")
        return obj

    def save_json(self, path: str | Path) -> None:
        Path(path).write_text(
            json.dumps(self.to_payload(), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )

    @classmethod
    def load_json(cls, path: str | Path) -> "KokunoLogRescaledHeatRepairTarget":
        return cls.from_payload(json.loads(Path(path).read_text(encoding="utf-8")))
