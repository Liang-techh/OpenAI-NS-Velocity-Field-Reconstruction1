"""Independent audit of the Kokuno signed-log heat-repair target.

This is Kokuno Agent 4 validation code for the Agent-1
``KokunoLogRescaledHeatRepairTarget`` handoff.  The production implementation
uses Gauss--Legendre transition quadrature, Gauss--Laguerre tail quadrature and
the repository heat-profile evaluator.  This module deliberately follows a
different numerical route:

* it consumes only the public terminal schedule scale logs;
* it uses the exact small-Z heat limit ``(H(Z)-1)/Z -> -h(1+h)``;
* it independently checks that every tested source schedule is so deep in the
  far tail that the first omitted heat-series correction is negligible;
* it integrates the compact transition with adaptive Gauss--Kronrod ``quad``;
* it evaluates the infinite power-law tail analytically in the heat-limit
  approximation; and
* it combines signed logarithms with ``math.fsum`` rather than the production
  helper.

The public scale chain used as input was independently preflighted by Kokuno
Agent 4 in PR #314.  This audit is intentionally scoped to the new signed-log
heat-target algebra in Agent-1 PR #320.  It does not apply a velocity
correction, does not fit pressure/forcing, and does not evaluate the formal
full-domain Navier--Stokes gate.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from scipy.integrate import quad

from .kokuno_log_rescaled_heat_repair_target import (
    KokunoLogRescaledHeatRepairTarget,
)
from .kokuno_outer_reserved_patch_schedule import (
    KokunoOuterReservedPatchSchedule,
)
from .kokuno_terminal_tail_schedule import KokunoTerminalTailSchedule


SCHEMA = "kokuno-agent4-independent-signed-log-heat-target-audit-v1"
FORMAL_MOMENTUM_GATE = 1.0e-3
FORMAL_DIVERGENCE_GATE = 1.0e-5
ST006_MOMENTUM_MAX = 0.1082289305112118
ST006_VOLUME_L2 = 0.10758432876230622

# Frozen before exact-head CI.  These are local implementation-consistency
# guards and do not replace the preregistered full-domain PDE thresholds.
LOCAL_GUARDS = {
    "max_public_log_abs_error": 2.0e-7,
    "max_medium_to_fine_log_change": 5.0e-9,
    "max_materialized_relative_error": 2.0e-7,
    "max_even_symmetry_log_error": 2.0e-12,
    "maximum_tested_log_z": -80.0,
    "minimum_old_tail_bug_log_shift": 2.0e-3,
    "minimum_default_span_decades": 400.0,
}

# Three independent adaptive-quadrature accuracy levels.
QUAD_TOLERANCES = (2.0e-7, 2.0e-9, 2.0e-11)
ETA_PROBES = (-0.8, -0.2, 0.0, 0.2, 0.8)
CHANNELS = ("C_p", "S", "I_sub")


def _flat_step_scalar(value: float) -> float:
    s = float(value)
    if not math.isfinite(s):
        raise ValueError("flat-step input must be finite")
    if s <= 0.0:
        return 0.0
    if s >= 1.0:
        return 1.0
    left = math.exp(-1.0 / (s * s))
    right = math.exp(-1.0 / ((1.0 - s) ** 2))
    return left / (left + right)


def _chi(y: float) -> float:
    return _flat_step_scalar((float(y) - 0.2) / 0.3)


def _terminal_multiplier(y: float, *, h: float, c_o: float) -> float:
    value = float(y)
    if value < 0.0:
        raise ValueError("terminal coordinate must be nonnegative")
    if value >= 3.0:
        return 1.0
    rho_o = c_o * h
    return 1.0 - rho_o * (1.0 - _flat_step_scalar((value - 1.0) / 2.0))


def _transition_quad(
    integrand: Any,
    *,
    tolerance: float,
) -> float:
    """Integrate the compact transition with an independent adaptive rule."""

    total = 0.0
    # Split at every formula-transition point so the adaptive integrator does
    # not have to discover the flat-step seams.
    for left, right in ((0.2, 0.5), (0.5, 1.0), (1.0, 3.0)):
        value, _ = quad(
            integrand,
            left,
            right,
            epsabs=max(2.0e-14, tolerance * 1.0e-4),
            epsrel=tolerance,
            limit=200,
        )
        total += float(value)
    return total


def _signed_log_sum(terms: list[tuple[int, float]]) -> tuple[int, float]:
    active = [
        (int(sign), float(log_abs))
        for sign, log_abs in terms
        if int(sign) != 0 and math.isfinite(float(log_abs))
    ]
    if not active:
        return 0, -math.inf
    reference = max(log_abs for _, log_abs in active)
    scaled = math.fsum(
        float(sign) * math.exp(log_abs - reference)
        for sign, log_abs in active
    )
    if scaled == 0.0:
        return 0, -math.inf
    return (1 if scaled > 0.0 else -1), reference + math.log(abs(scaled))


def _case_matrix() -> list[tuple[str, KokunoTerminalTailSchedule]]:
    return [
        ("default_extreme", KokunoTerminalTailSchedule()),
        (
            "moderate_materializable",
            KokunoTerminalTailSchedule(
                outer_schedule=KokunoOuterReservedPatchSchedule(
                    M_d=0.1,
                    lambda_outer=0.65,
                    h=0.005,
                )
            ),
        ),
        (
            "perturbed_outer",
            KokunoTerminalTailSchedule(
                outer_schedule=KokunoOuterReservedPatchSchedule(
                    M_d=2.0,
                    log_p_star_margin=1.5,
                    C=2.5,
                    lambda_outer=0.08,
                    h=0.0035,
                    repair_log_offset=-17.0,
                ),
                T_f=72.0,
                c_o=0.04,
            ),
        ),
    ]


def _independent_scalar_target(
    schedule: KokunoTerminalTailSchedule,
    eta: float,
    *,
    tolerance: float,
    old_tail_bug: bool = False,
) -> dict[str, Any]:
    """Return an independent signed-log target for one eta.

    The source schedules tested here have ``Z=2(1-eta^2)/X`` astronomically
    close to zero throughout the terminal/tail region.  The heat ratio is
    therefore evaluated from its exact Z=0 limit rather than through the
    production heat quadrature.  ``maximum_log_z`` and the first omitted series
    term scale are returned explicitly so that this approximation cannot be
    silently used outside its regime.
    """

    value = float(eta)
    if not math.isfinite(value) or abs(value) > 1.0:
        raise ValueError("eta must lie in [-1,1]")
    d = 1.0 - value * value
    if d == 0.0:
        return {
            "sign": [0, 0, 0],
            "log_abs": [-math.inf, -math.inf, -math.inf],
            "maximum_log_z": -math.inf,
            "first_omitted_relative_log_scale": -math.inf,
        }

    h = float(schedule.h)
    c_o = float(schedule.c_o)
    A = 0.5 + h
    outer = schedule.outer_schedule
    Lx = float(schedule.log_X_tail)
    Lc = float(schedule.log_c_inf)
    Lxs = float(outer.log_X_star)
    Les = float(outer.log_e_star)

    maximum_log_z = math.log(2.0 * d) - Lx
    # H(Z)=1-h(h+1)Z + 1/2*h(h+1)^2(h+2)Z^2+...
    # so the first relative correction to R=(H-1)/Z is
    # 0.5*(h+1)*(h+2)*Z.
    first_omitted_relative_log_scale = (
        math.log(0.5 * (h + 1.0) * (h + 2.0)) + maximum_log_z
    )
    if maximum_log_z > LOCAL_GUARDS["maximum_tested_log_z"]:
        raise RuntimeError(
            "heat-limit independent oracle left its preregistered small-Z regime"
        )

    R0 = -h * (1.0 + h)

    def f(y: float) -> float:
        return _terminal_multiplier(y, h=h, c_o=c_o)

    def transition(
        *,
        rate: float,
        chi_power: int,
        f_power: int,
        r_power: int,
    ) -> float:
        return _transition_quad(
            lambda y: (
                (f(y) ** f_power)
                * math.exp(-rate * y)
                * (_chi(y) ** chi_power)
                * (R0 ** r_power)
            ),
            tolerance=tolerance,
        )

    def tail(*, rate: float, r_power: int) -> float:
        return (R0 ** r_power) * math.exp(-3.0 * rate) / rate

    k_cp_1 = 2.0 * A + 1.0
    k_cp_2 = 2.0 * A + 2.0
    cp_tail_factor = 0.5 if old_tail_bug else 1.0
    b_cp_1 = 2.0 * d * (
        transition(
            rate=k_cp_1,
            chi_power=1,
            f_power=2,
            r_power=1,
        )
        + cp_tail_factor * tail(rate=k_cp_1, r_power=1)
    )
    b_cp_2 = 2.0 * d * d * (
        transition(
            rate=k_cp_2,
            chi_power=2,
            f_power=2,
            r_power=2,
        )
        + tail(rate=k_cp_2, r_power=2)
    )
    cp_sign, cp_log = _signed_log_sum(
        [
            (
                1 if b_cp_1 > 0.0 else -1,
                2.0 * Lc
                - k_cp_1 * Lx
                - 2.0 * Les
                + math.log(abs(b_cp_1)),
            ),
            (
                1 if b_cp_2 > 0.0 else -1,
                2.0 * Lc
                - k_cp_2 * Lx
                - 2.0 * Les
                + math.log(abs(b_cp_2)),
            ),
        ]
    )

    k_s_1 = 1.0 + 2.0 * h
    k_s_2 = 2.0 + 2.0 * h
    b_s_1 = -2.0 * d * (
        transition(
            rate=k_s_1,
            chi_power=1,
            f_power=2,
            r_power=1,
        )
        + tail(rate=k_s_1, r_power=1)
    )
    b_s_2 = -2.0 * d * d * (
        transition(
            rate=k_s_2,
            chi_power=2,
            f_power=2,
            r_power=2,
        )
        + tail(rate=k_s_2, r_power=2)
    )
    s_sign, s_log = _signed_log_sum(
        [
            (
                1 if b_s_1 > 0.0 else -1,
                2.0 * Lc
                - k_s_1 * Lx
                - Lxs
                - 2.0 * Les
                + math.log(abs(b_s_1)),
            ),
            (
                1 if b_s_2 > 0.0 else -1,
                2.0 * Lc
                - k_s_2 * Lx
                - Lxs
                - 2.0 * Les
                + math.log(abs(b_s_2)),
            ),
        ]
    )

    b_i = 2.0 * math.sqrt(2.0) * d * (
        transition(
            rate=h,
            chi_power=1,
            f_power=1,
            r_power=1,
        )
        + tail(rate=h, r_power=1)
    )
    i_sign = 1 if b_i > 0.0 else -1
    i_log = (
        Lc
        - h * Lx
        - 1.5 * Lxs
        - Les
        + math.log(abs(b_i))
    )

    return {
        "sign": [-cp_sign, -s_sign, -i_sign],
        "log_abs": [cp_log, s_log, i_log],
        "maximum_log_z": maximum_log_z,
        "first_omitted_relative_log_scale": first_omitted_relative_log_scale,
    }


def _underflow_channels(sign: np.ndarray, logs: np.ndarray) -> list[str]:
    lower = math.log(np.nextafter(0.0, 1.0))
    return [
        CHANNELS[index]
        for index in range(3)
        if int(sign[index]) != 0 and float(logs[index]) < lower
    ]


def _materialized_from_signed_log(sign: np.ndarray, logs: np.ndarray) -> np.ndarray:
    lower = math.log(np.nextafter(0.0, 1.0))
    upper = math.log(np.finfo(float).max)
    if np.any((sign != 0) & ((logs < lower) | (logs > upper))):
        raise OverflowError("independent target is not fully float64 materializable")
    values = np.zeros(3, dtype=float)
    active = sign != 0
    values[active] = sign[active] * np.exp(logs[active])
    return values


def _relative_error(actual: np.ndarray, expected: np.ndarray) -> float:
    scale = np.maximum(np.abs(expected), np.finfo(float).tiny)
    return float(np.max(np.abs(actual - expected) / scale))


def run_audit() -> dict[str, Any]:
    cases: list[dict[str, Any]] = []
    max_public_log_abs_error = 0.0
    max_medium_to_fine_log_change = 0.0
    max_materialized_relative_error = 0.0
    max_even_symmetry_log_error = 0.0
    maximum_tested_log_z = -math.inf
    minimum_old_tail_bug_log_shift = math.inf
    sign_match_all = True
    underflow_classification_matches = True
    endpoint_zero_match = True
    default_span_decades = 0.0

    for name, schedule in _case_matrix():
        public = KokunoLogRescaledHeatRepairTarget(terminal_schedule=schedule)
        level_records: list[dict[str, Any]] = []

        for tolerance in QUAD_TOLERANCES:
            independent_rows: list[dict[str, Any]] = []
            for eta in ETA_PROBES:
                reference = _independent_scalar_target(
                    schedule,
                    eta,
                    tolerance=tolerance,
                )
                encoded = public.signed_log_target(eta)
                public_sign = np.asarray(encoded["sign"], dtype=int)
                public_log = np.asarray(encoded["log_abs"], dtype=float)
                reference_sign = np.asarray(reference["sign"], dtype=int)
                reference_log = np.asarray(reference["log_abs"], dtype=float)

                sign_match = bool(np.array_equal(public_sign, reference_sign))
                sign_match_all = sign_match_all and sign_match
                active = reference_sign != 0
                log_error = (
                    float(np.max(np.abs(public_log[active] - reference_log[active])))
                    if np.any(active)
                    else 0.0
                )
                max_public_log_abs_error = max(
                    max_public_log_abs_error,
                    log_error,
                )
                maximum_tested_log_z = max(
                    maximum_tested_log_z,
                    float(reference["maximum_log_z"]),
                )

                independent_rows.append(
                    {
                        "eta": eta,
                        "sign": reference_sign.tolist(),
                        "log_abs": reference_log.tolist(),
                        "public_log_abs_error": log_error,
                        "maximum_log_z": float(reference["maximum_log_z"]),
                        "first_omitted_relative_log_scale": float(
                            reference["first_omitted_relative_log_scale"]
                        ),
                    }
                )

            level_records.append(
                {
                    "tolerance": tolerance,
                    "rows": independent_rows,
                }
            )

        medium = level_records[1]["rows"]
        fine = level_records[2]["rows"]
        for medium_row, fine_row in zip(medium, fine, strict=True):
            medium_log = np.asarray(medium_row["log_abs"], dtype=float)
            fine_log = np.asarray(fine_row["log_abs"], dtype=float)
            active = np.isfinite(fine_log)
            if np.any(active):
                max_medium_to_fine_log_change = max(
                    max_medium_to_fine_log_change,
                    float(np.max(np.abs(medium_log[active] - fine_log[active]))),
                )

        # Even-in-eta structure is checked on the independent oracle, not inferred
        # from the production object.
        finest_by_eta = {
            float(row["eta"]): row
            for row in level_records[-1]["rows"]
        }
        for eta in (0.2, 0.8):
            left = np.asarray(finest_by_eta[-eta]["log_abs"], dtype=float)
            right = np.asarray(finest_by_eta[eta]["log_abs"], dtype=float)
            max_even_symmetry_log_error = max(
                max_even_symmetry_log_error,
                float(np.max(np.abs(left - right))),
            )

        finest_eta02 = finest_by_eta[0.2]
        fine_sign = np.asarray(finest_eta02["sign"], dtype=int)
        fine_log = np.asarray(finest_eta02["log_abs"], dtype=float)
        public_eta02 = public.signed_log_target(0.2)
        public_sign = np.asarray(public_eta02["sign"], dtype=int)
        public_log = np.asarray(public_eta02["log_abs"], dtype=float)
        independent_underflow = _underflow_channels(fine_sign, fine_log)
        public_underflow = _underflow_channels(public_sign, public_log)
        report_underflow = list(public.precision_report(0.2)["underflow_channels"])
        underflow_match = (
            independent_underflow == public_underflow == report_underflow
        )
        underflow_classification_matches = (
            underflow_classification_matches and underflow_match
        )

        span_decades = float((np.max(fine_log) - np.min(fine_log)) / math.log(10.0))
        if name == "default_extreme":
            default_span_decades = span_decades

        materialized_relative_error = None
        if not independent_underflow:
            independent_values = _materialized_from_signed_log(fine_sign, fine_log)
            public_values = np.asarray(public.materialize_target(0.2), dtype=float)
            materialized_relative_error = _relative_error(
                public_values,
                independent_values,
            )
            max_materialized_relative_error = max(
                max_materialized_relative_error,
                materialized_relative_error,
            )
        else:
            try:
                public.materialize_target(0.2)
            except OverflowError:
                pass
            else:
                underflow_classification_matches = False

        # Recreate the old factor-of-two terminal C_p bug.  A useful validator
        # must react measurably rather than agreeing with both formulas.
        mutated = _independent_scalar_target(
            schedule,
            0.2,
            tolerance=QUAD_TOLERANCES[-1],
            old_tail_bug=True,
        )
        mutation_shift = abs(float(mutated["log_abs"][0]) - float(fine_log[0]))
        minimum_old_tail_bug_log_shift = min(
            minimum_old_tail_bug_log_shift,
            mutation_shift,
        )

        # eta=+-1 are exact zero channels in both paths.
        for endpoint in (-1.0, 1.0):
            reference = _independent_scalar_target(
                schedule,
                endpoint,
                tolerance=QUAD_TOLERANCES[-1],
            )
            encoded = public.signed_log_target(endpoint)
            endpoint_ok = (
                reference["sign"] == [0, 0, 0]
                and all(math.isinf(v) and v < 0.0 for v in reference["log_abs"])
                and np.all(np.asarray(encoded["sign"], dtype=int) == 0)
                and np.all(np.isneginf(np.asarray(encoded["log_abs"], dtype=float)))
            )
            endpoint_zero_match = endpoint_zero_match and bool(endpoint_ok)

        cases.append(
            {
                "name": name,
                "schedule": {
                    "h": float(schedule.h),
                    "c_o": float(schedule.c_o),
                    "log_X_tail": float(schedule.log_X_tail),
                    "log_c_inf": float(schedule.log_c_inf),
                    "log_X_star": float(schedule.outer_schedule.log_X_star),
                    "log_e_star": float(schedule.outer_schedule.log_e_star),
                },
                "finest_eta_0p2": {
                    "sign": fine_sign.tolist(),
                    "log_abs": fine_log.tolist(),
                    "log10_abs": (fine_log / math.log(10.0)).tolist(),
                    "span_decades": span_decades,
                    "underflow_channels": independent_underflow,
                },
                "materialized_relative_error": materialized_relative_error,
                "old_tail_bug_cp_log_shift": mutation_shift,
                "levels": level_records,
            }
        )

    guards = {
        "public_log_accuracy": (
            max_public_log_abs_error
            <= LOCAL_GUARDS["max_public_log_abs_error"]
        ),
        "independent_refinement": (
            max_medium_to_fine_log_change
            <= LOCAL_GUARDS["max_medium_to_fine_log_change"]
        ),
        "materialized_accuracy": (
            max_materialized_relative_error
            <= LOCAL_GUARDS["max_materialized_relative_error"]
        ),
        "even_symmetry": (
            max_even_symmetry_log_error
            <= LOCAL_GUARDS["max_even_symmetry_log_error"]
        ),
        "small_z_regime": (
            maximum_tested_log_z
            <= LOCAL_GUARDS["maximum_tested_log_z"]
        ),
        "old_tail_bug_detected": (
            minimum_old_tail_bug_log_shift
            >= LOCAL_GUARDS["minimum_old_tail_bug_log_shift"]
        ),
        "default_span_preserved": (
            default_span_decades
            >= LOCAL_GUARDS["minimum_default_span_decades"]
        ),
        "signs_match": sign_match_all,
        "underflow_classification_matches": underflow_classification_matches,
        "eta_endpoints_match": endpoint_zero_match,
    }
    structural_preflight_passed = all(guards.values())

    return {
        "schema": SCHEMA,
        "task": "KOKUNO-A4-SIGNED-LOG-HEAT-TARGET-INDEPENDENT-013",
        "source_agent1_pr": 320,
        "source_agent1_head": "2a7d51a05bcf46a2071ed7e767e8f9f467a79f86",
        "previous_scale_chain_audit_pr": 314,
        "quadrature_tolerances": list(QUAD_TOLERANCES),
        "eta_probes": list(ETA_PROBES),
        "local_guards": dict(LOCAL_GUARDS),
        "metrics": {
            "max_public_log_abs_error": max_public_log_abs_error,
            "max_medium_to_fine_log_change": max_medium_to_fine_log_change,
            "max_materialized_relative_error": max_materialized_relative_error,
            "max_even_symmetry_log_error": max_even_symmetry_log_error,
            "maximum_tested_log_z": maximum_tested_log_z,
            "minimum_old_tail_bug_log_shift": minimum_old_tail_bug_log_shift,
            "default_span_decades": default_span_decades,
        },
        "guards": guards,
        "cases": cases,
        "structural_preflight_passed": structural_preflight_passed,
        "formal_full_domain_pde_gate_assessed": False,
        "formal_momentum_gate": FORMAL_MOMENTUM_GATE,
        "formal_divergence_gate": FORMAL_DIVERGENCE_GATE,
        "normalized_ns_residual_le_1e-3_claimed": False,
        "pde_validated": False,
        "st006_context": {
            "momentum_max": ST006_MOMENTUM_MAX,
            "volume_l2": ST006_VOLUME_L2,
            "directly_comparable_to_this_local_audit": False,
        },
        "limitations": [
            "This audit validates the signed-log heat-repair target algebra, not a global velocity field.",
            "No three-bump repair is applied to I2 in Agent-1 PR #320.",
            "No leading+oscillatory+correction composite velocity/pressure contract exists yet.",
            "No pressure or forcing is fitted here, and no free f=R residual cancellation is used.",
            "The exact heat factor is replaced only after checking the tested schedules lie deep in the small-Z regime.",
        ],
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Independently audit the Kokuno signed-log heat-repair target"
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "artifacts/kokuno_agent4/independent_signed_log_heat_target_audit.json"
        ),
    )
    args = parser.parse_args(argv)

    report = run_audit()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(report["metrics"], indent=2, sort_keys=True))
    print(
        "structural_preflight_passed="
        f"{report['structural_preflight_passed']}"
    )
    return 0 if report["structural_preflight_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
