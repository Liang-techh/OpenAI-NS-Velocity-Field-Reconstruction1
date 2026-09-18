"""High-precision continuous-moment backend for the Kokuno I2 heat repair.

Agent-4's independent audit of the first hierarchical implementation showed that
Decimal Newton arithmetic was not enough: the three-bump moment tensor itself
had first been assembled with float64 Gauss--Legendre quadrature, and its
``I_sub`` row missed an independently refined continuous integral by about
5.5e-5 relative.  For the default existence-style schedule the signed-log
repair target spans more than 400 decades, so converting those float64 tensor
entries to Decimal cannot preserve the source moment constraints.

This module keeps the public Kokuno three-row formulas and the repository's
explicit autonomous bump geometry, but rebuilds the five required moments on
each disjoint bump support directly with arbitrary-precision tanh--sinh
quadrature.  The resulting Decimal tensor is then consumed by the existing
hierarchical solve and local I2 velocity evaluator.

The tanh--sinh rule and concrete bump locations are autonomous numerical
choices, not parameters recovered from OpenAI.  This is still a local I2
construction.  Independent high-precision alternate-quadrature certification,
global core-to-heat assembly and full-domain Navier--Stokes validation remain
separate tasks.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, getcontext
from functools import lru_cache
import json
import math
from pathlib import Path
from typing import Any

import mpmath as mp

from .kokuno_hierarchical_heat_repair import KokunoHierarchicalHeatRepair
from .kokuno_log_rescaled_heat_repair_target import KokunoLogRescaledHeatRepairTarget


SCHEMA = "kokuno-high-precision-hierarchical-heat-repair-v1"
_BUMP_CENTERS = ("0.96", "1.24", "1.53")
_BUMP_HALF_WIDTH = "0.105"
_TANH_SINH_SUBDIVISIONS = 2048
_EXTRA_WORK_DIGITS = 48


@lru_cache(maxsize=16)
def _base_moment_strings(
    lambda_text: str, precision_digits: int
) -> tuple[tuple[str, str, str, str, str], ...]:
    """Return unscaled bump moments as high-precision decimal strings.

    For each autonomous bump the five entries are

    ``(L_Cp/f, L_S/f, L_I, Q_Cp, Q_S)``

    with the source normalized weights.  Cross-bump quadratic terms are exactly
    zero because the declared compact supports are disjoint.

    A fixed tanh--sinh mesh is used in the transformed support coordinate
    ``s in [-1,1]``.  The endpoint tail is stopped only once the C-infinity
    bump itself is smaller than the requested working precision.
    """

    digits = int(precision_digits)
    if digits < 80 or digits > 900:
        raise ValueError("precision_digits must lie in [80,900]")
    lam_float = float(lambda_text)
    if not math.isfinite(lam_float) or not (0.0 < lam_float <= 0.5):
        raise ValueError("lambda_outer must lie in (0,0.5]")

    work_digits = digits + _EXTRA_WORK_DIGITS
    with mp.workdps(work_digits):
        lam = mp.mpf(lambda_text)
        half_width = mp.mpf(_BUMP_HALF_WIDTH)
        h = mp.mpf(1) / _TANH_SINH_SUBDIVISIONS
        log10 = mp.log(10)
        cutoff = -(digits + 32) * log10
        pi_over_2 = mp.pi / 2
        rows: list[tuple[str, str, str, str, str]] = []

        for center_text in _BUMP_CENTERS:
            center = mp.mpf(center_text)
            sums = [mp.mpf(0) for _ in range(5)]
            k = 0
            while True:
                t = k * h
                sinh_t = mp.sinh(t)
                cosh_t = mp.cosh(t)
                a = pi_over_2 * sinh_t
                s = mp.tanh(a)
                one_minus_s2 = 1 - s * s
                if one_minus_s2 <= 0:
                    break
                bump_exponent = 1 - 1 / one_minus_s2
                if k > 0 and bump_exponent < cutoff:
                    break

                ds_dt = pi_over_2 * cosh_t / (mp.cosh(a) ** 2)
                bump = mp.exp(bump_exponent)

                def accumulate(s_value: mp.mpf) -> None:
                    x = center + half_width * s_value
                    beta = bump
                    dx_weight = half_width * ds_dt
                    values = (
                        mp.power(x, -mp.mpf("1.5") - lam) * beta,
                        -mp.power(x, -mp.mpf("0.5") - lam) * beta,
                        mp.sqrt(2 * x) * beta,
                        mp.mpf("0.5") * beta * beta / x,
                        -mp.mpf("0.5") * beta * beta,
                    )
                    for index, value in enumerate(values):
                        sums[index] += value * dx_weight

                accumulate(s)
                if k:
                    accumulate(-s)
                k += 1
                if k > 8192:
                    raise RuntimeError("tanh-sinh support integration did not terminate")

            moments = [value * h for value in sums]
            rows.append(
                tuple(
                    mp.nstr(value, n=digits + 8, strip_zeros=False)
                    for value in moments
                )
            )

    return tuple(rows)


def _decimal(value: float | str | Decimal) -> Decimal:
    if isinstance(value, Decimal):
        return value
    if isinstance(value, str):
        return Decimal(value)
    number = float(value)
    if not math.isfinite(number):
        raise ValueError("non-finite value cannot be converted to Decimal")
    return Decimal(repr(number))


@dataclass(frozen=True)
class KokunoHighPrecisionHierarchicalHeatRepair(KokunoHierarchicalHeatRepair):
    """Hierarchical I2 repair using a high-precision continuous moment tensor."""

    def _discrete_tensors(
        self, f_eta: float
    ) -> tuple[list[list[Decimal]], list[list[list[Decimal]]]]:
        # ``solve`` enters a Decimal localcontext at the target-dependent
        # precision before calling this method, so the current context is the
        # precision that must survive the >400-decade hierarchy.
        precision_digits = int(getcontext().prec)
        lambda_text = repr(float(self.outer_schedule.lambda_outer))
        base = _base_moment_strings(lambda_text, precision_digits)
        f_value = _decimal(f_eta)

        linear = [[Decimal(0) for _ in range(3)] for _ in range(3)]
        quadratic = [
            [[Decimal(0) for _ in range(3)] for _ in range(3)]
            for _ in range(3)
        ]
        for column, moments in enumerate(base):
            cp_linear, s_linear, i_linear, cp_quad, s_quad = map(
                Decimal, moments
            )
            linear[0][column] = f_value * cp_linear
            linear[1][column] = f_value * s_linear
            linear[2][column] = i_linear
            quadratic[0][column][column] = cp_quad
            quadratic[1][column][column] = s_quad
        return linear, quadratic

    def precision_report(self, eta: float) -> dict[str, Any]:
        report = dict(super().precision_report(float(eta)))
        report.update(
            {
                "moment_tensor_backend": "arbitrary-precision tanh-sinh",
                "moment_tensor_subdivisions_per_t_unit": _TANH_SINH_SUBDIVISIONS,
                "moment_tensor_extra_work_digits": _EXTRA_WORK_DIGITS,
                "float64_moment_tensor_used": False,
                "independent_high_precision_moment_audit_completed": False,
                "continuous_moment_compensation_certified": False,
            }
        )
        return report

    def _unsigned_payload(self) -> dict[str, Any]:
        payload = super()._unsigned_payload()
        payload["schema"] = SCHEMA
        payload["autonomous_numerics"] = {
            "high_precision_engine": "mpmath arbitrary precision",
            "moment_tensor": (
                "per-support tanh-sinh quadrature in the normalized bump "
                "coordinate; no float64 moment coefficient is promoted to Decimal"
            ),
            "tanh_sinh_subdivisions_per_t_unit": _TANH_SINH_SUBDIVISIONS,
            "extra_work_digits": _EXTRA_WORK_DIGITS,
            "bump_centers": list(_BUMP_CENTERS),
            "bump_half_width": _BUMP_HALF_WIDTH,
            "hierarchy": (
                "existing per-target-channel linear pieces followed by Decimal "
                "Newton counterterms against the high-precision tensor"
            ),
        }
        truth = dict(payload["truth_boundary"])
        truth.update(
            {
                "repository_float64_moment_tensor_used": False,
                "high_precision_continuous_moment_tensor_used": True,
                "independent_high_precision_moment_audit_completed": False,
                "continuous_source_moment_compensation_certified": False,
                "heat_compensation_completed": False,
                "global_leading_profile_reconstructed": False,
                "formal_full_domain_pde_gate_assessed": False,
                "pde_validated": False,
                "paper_exact": False,
            }
        )
        payload["truth_boundary"] = truth
        return payload

    @classmethod
    def from_payload(
        cls, payload: dict[str, Any]
    ) -> "KokunoHighPrecisionHierarchicalHeatRepair":
        if not isinstance(payload, dict) or payload.get("schema") != SCHEMA:
            raise ValueError("unexpected high-precision heat-repair schema")
        target_payload = payload.get("target")
        params = payload.get("parameters")
        if not isinstance(target_payload, dict) or not isinstance(params, dict):
            raise ValueError("serialized target/parameters are missing")
        obj = cls(
            target=KokunoLogRescaledHeatRepairTarget.from_payload(target_payload),
            **params,
        )
        if obj.to_payload() != payload:
            raise ValueError("high-precision heat-repair payload hash or content mismatch")
        return obj

    def save_json(self, path: str | Path) -> None:
        Path(path).write_text(
            json.dumps(self.to_payload(), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )

    @classmethod
    def load_json(
        cls, path: str | Path
    ) -> "KokunoHighPrecisionHierarchicalHeatRepair":
        return cls.from_payload(json.loads(Path(path).read_text(encoding="utf-8")))
