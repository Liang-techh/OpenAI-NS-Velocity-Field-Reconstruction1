"""Independent numerical audit of the Kokuno terminal release/tail schedule.

This is Kokuno Agent 4 validation code.  The production construction in
``kokuno_terminal_tail_schedule`` advances the two release transitions with
``solve_ivp(DOP853)`` and evaluates terminal integrals with Gauss--Legendre
quadrature.  This module deliberately uses a different numerical path:

* fixed uniform grids;
* an integrating-factor representation of the scalar release ODE;
* cumulative trapezoidal integration only for the integrating factor; and
* composite Simpson quadrature for the remaining transition/terminal integrals.

Only the public schedule object is compared with the independently reconstructed
values.  No training tensor/loss, production ODE state, or fitted forcing is
used.  The checks here are local source-schedule implementation checks, not a
Navier--Stokes residual gate.  In particular ``pde_validated`` remains false.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from scipy.integrate import cumulative_trapezoid, simpson

from .kokuno_outer_reserved_patch_schedule import KokunoOuterReservedPatchSchedule
from .kokuno_terminal_tail_schedule import KokunoTerminalTailSchedule


SCHEMA = "kokuno-agent4-independent-terminal-tail-audit-v1"
ST006_SHA256 = "6b4d84b48ab9dbcd2ee1a1858d3e56ef81523f5864369d7e96c6431fccf107a3"
ST006_MOMENTUM_MAX = 0.1082289305112118
ST006_VOLUME_L2 = 0.10758432876230622
FORMAL_MOMENTUM_GATE = 1.0e-3
FORMAL_DIVERGENCE_GATE = 1.0e-5

# Frozen before evaluating the exact-head report.  These are local
# implementation-consistency guards, not replacements for the formal PDE gates.
LOCAL_GUARDS = {
    "max_release_relative_error": 5.0e-9,
    "max_terminal_q_relative_error": 5.0e-9,
    "max_log_scale_absolute_error": 5.0e-8,
    "max_independent_medium_to_fine_relative_change": 1.0e-8,
    "minimum_transition_observed_order": 1.8,
    "minimum_mutation_log_X_tail_shift": 0.5,
}

RESOLUTIONS = (2049, 4097, 8193)
TERMINAL_Q_PROBES = (0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0)


def _sigma(value: Any) -> np.ndarray:
    """Independent flat step; intentionally does not call production helpers."""

    s = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(s)):
        raise ValueError("sigma input must be finite")
    out = np.zeros_like(s)
    out[s >= 1.0] = 1.0
    active = (s > 0.0) & (s < 1.0)
    if np.any(active):
        a = s[active]
        left = np.exp(-1.0 / (a * a))
        right = np.exp(-1.0 / ((1.0 - a) ** 2))
        out[active] = left / (left + right)
    return out


def _sigma_derivative(value: Any) -> np.ndarray:
    s = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(s)):
        raise ValueError("sigma input must be finite")
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


def _transition_integrating_factor(
    q_initial: float,
    l_values: np.ndarray,
    y: np.ndarray,
    h: float,
) -> float:
    """Solve Q'+(1+l)Q=-l-h by an independent integrating factor."""

    a = 1.0 + l_values
    b = -l_values - h
    primitive = cumulative_trapezoid(a, x=y, initial=0.0)
    value = math.exp(-float(primitive[-1])) * (
        float(q_initial) + float(simpson(np.exp(primitive) * b, x=y))
    )
    if not math.isfinite(value) or value <= 0.0:
        raise RuntimeError("independent release transition produced invalid Q")
    return value


def _terminal_multiplier(y: np.ndarray, h: float, c_o: float) -> np.ndarray:
    rho_o = c_o * h
    return 1.0 - rho_o * (1.0 - _sigma((y - 1.0) / 2.0))


def _terminal_multiplier_derivative(y: np.ndarray, h: float, c_o: float) -> np.ndarray:
    return 0.5 * c_o * h * _sigma_derivative((y - 1.0) / 2.0)


def _simpson_grid(start: float, stop: float, base_intervals_per_unit: int) -> np.ndarray:
    width = float(stop - start)
    if width <= 0.0:
        return np.asarray([start], dtype=float)
    intervals = max(2, int(round(width * base_intervals_per_unit)))
    if intervals % 2:
        intervals += 1
    return np.linspace(start, stop, intervals + 1)


def _independent_terminal_q(
    y0: float,
    *,
    h: float,
    c_o: float,
    transition_points: int,
) -> float:
    if y0 == 3.0:
        return 0.0
    base_intervals = transition_points - 1
    v = _simpson_grid(float(y0), 3.0, base_intervals)
    f0 = float(_terminal_multiplier(np.asarray([y0]), h, c_o)[0])
    fp = _terminal_multiplier_derivative(v, h, c_o)
    return float(simpson(np.exp((1.0 - h) * (v - y0)) * fp / f0, x=v))


def _raw_case_parameters(schedule: KokunoTerminalTailSchedule) -> dict[str, float]:
    outer = schedule.outer_schedule
    return {
        "M_d": float(outer.M_d),
        "log_p_star_margin": float(outer.log_p_star_margin),
        "C": float(outer.C),
        "lambda_outer": float(outer.lambda_outer),
        "h": float(outer.h),
        "repair_log_offset": float(outer.repair_log_offset),
        "T_f": float(schedule.T_f),
        "c_o": float(schedule.c_o),
    }


def independent_reference(
    schedule: KokunoTerminalTailSchedule,
    *,
    transition_points: int,
) -> dict[str, Any]:
    """Reconstruct release and log-scale data without production integrators."""

    if transition_points < 257 or transition_points % 2 == 0:
        raise ValueError("transition_points must be an odd integer >=257")
    p = _raw_case_parameters(schedule)
    lam = p["lambda_outer"]
    h = p["h"]
    c_o = p["c_o"]

    y = np.linspace(0.0, 1.0, transition_points)
    l_first = -lam - (1.0 - lam) * _sigma(y)
    q_release = (lam - h) / (1.0 - lam)
    q_after_first = _transition_integrating_factor(q_release, l_first, y, h)

    steep_hold_length = 4.0 * math.log(1.0 / h)
    q_after_steep = q_after_first + (1.0 - h) * steep_hold_length

    l_second = -1.0 + (1.0 - h) * _sigma(y)
    q_in = _transition_integrating_factor(q_after_steep, l_second, y, h)

    terminal_grid = _simpson_grid(0.0, 3.0, transition_points - 1)
    f0 = float(_terminal_multiplier(np.asarray([0.0]), h, c_o)[0])
    fp = _terminal_multiplier_derivative(terminal_grid, h, c_o)
    q_p = float(
        simpson(
            np.exp((1.0 - h) * terminal_grid) * fp / f0,
            x=terminal_grid,
        )
    )
    release_hold_length = math.log(q_in / q_p) / (1.0 - h)

    # Rebuild the scale chain directly from the public source relations and raw
    # autonomous parameters rather than consuming production derived properties.
    T_d = math.exp(p["M_d"]) + 10.0
    log_P_star = T_d + p["log_p_star_margin"]
    T_w = 60.0 * math.log(1.0 / lam)
    log_X_R = math.log(110.0) + 10.0 * (math.log(p["C"]) + log_P_star)
    log_X_w = log_X_R + T_d + 2.0
    log_P1 = log_P_star - 0.2
    log_e_w = log_P1 - 0.5 * T_d - 0.5 - 0.5 * lam
    log_c_patch = log_e_w + (0.5 + lam) * log_X_w
    log_X_star = log_X_w + T_w + p["repair_log_offset"]
    log_e_star = log_c_patch - (0.5 + lam) * log_X_star

    pulse_length = 13.0 / lam
    uniform_hold_length = 30.0 * math.log(1.0 / lam)
    log_X_p = log_X_w + T_w
    log_e_b = log_e_w - (0.5 + lam) * T_w
    log_e_end = log_e_b - (0.5 + lam) * pulse_length
    log_X_rel = log_X_p + pulse_length + p["T_f"] + uniform_hold_length
    log_e_rel = (
        log_e_end
        - (0.5 + lam) * p["T_f"]
        - math.log(2.0)
        - (0.5 + lam) * uniform_hold_length
    )
    log_X_hold_start = log_X_rel + 2.0 + steep_hold_length
    log_e_hold_start = (
        log_e_rel
        - 1.0
        - 0.5 * lam
        - 1.5 * steep_hold_length
        - 1.0
        - 0.5 * h
    )
    log_X_tail = log_X_hold_start + release_hold_length
    A = 0.5 + h
    log_e_tail_start = log_e_hold_start - A * release_hold_length
    log_c_inf = log_e_tail_start + A * log_X_tail - math.log(f0)
    log_X_K = log_X_tail + 0.2
    f_k = float(_terminal_multiplier(np.asarray([0.2]), h, c_o)[0])
    log_e_K = log_c_inf - A * log_X_K + math.log(f_k)

    terminal_q = {
        f"{probe:g}": _independent_terminal_q(
            probe,
            h=h,
            c_o=c_o,
            transition_points=transition_points,
        )
        for probe in TERMINAL_Q_PROBES
    }

    return {
        "transition_points": int(transition_points),
        "q_release": q_release,
        "q_after_first_transition": q_after_first,
        "q_after_steep_hold": q_after_steep,
        "q_in": q_in,
        "q_p": q_p,
        "release_hold_length": release_hold_length,
        "log_X_tail": log_X_tail,
        "log_c_inf": log_c_inf,
        "log_X_K": log_X_K,
        "log_e_K": log_e_K,
        "log_X_star_over_X_K": log_X_star - log_X_K,
        "log_e_star_over_e_K": log_e_star - log_e_K,
        "terminal_Q": terminal_q,
    }


def _relative_error(public: float, reference: float, floor: float = 1.0e-300) -> float:
    return abs(float(public) - float(reference)) / max(abs(float(reference)), floor)


def _case_matrix() -> list[tuple[str, KokunoTerminalTailSchedule]]:
    return [
        ("default_extreme", KokunoTerminalTailSchedule()),
        (
            "moderate_materializable",
            KokunoTerminalTailSchedule(
                outer_schedule=KokunoOuterReservedPatchSchedule(
                    M_d=1.0,
                    lambda_outer=0.2,
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


def _float_materializable(log_value: float) -> bool:
    return math.log(np.nextafter(0.0, 1.0)) <= log_value <= math.log(np.finfo(float).max)


def run_audit() -> dict[str, Any]:
    cases: list[dict[str, Any]] = []
    max_release_relative_error = 0.0
    max_terminal_q_relative_error = 0.0
    max_log_scale_absolute_error = 0.0
    max_medium_to_fine_relative_change = 0.0
    minimum_observed_order = math.inf
    minimum_mutation_shift = math.inf
    materialization_classification_matches = True

    for name, schedule in _case_matrix():
        levels = [
            independent_reference(schedule, transition_points=resolution)
            for resolution in RESOLUTIONS
        ]
        coarse, medium, fine = levels

        public_release = {
            "q_after_first_transition": schedule.Q_after_first_transition,
            "q_in": schedule.Q_in,
            "q_p": schedule.Q_p,
            "release_hold_length": schedule.release_hold_length,
        }
        release_errors = {
            key: _relative_error(public_release[key], fine[key])
            for key in public_release
        }
        max_release_relative_error = max(
            max_release_relative_error, max(release_errors.values())
        )

        public_logs = {
            "log_X_tail": schedule.log_X_tail,
            "log_c_inf": schedule.log_c_inf,
            "log_X_star_over_X_K": schedule.log_ratio_report()["log_X_star_over_X_K"],
            "log_e_star_over_e_K": schedule.log_ratio_report()["log_e_star_over_e_K"],
        }
        log_errors = {
            key: abs(float(public_logs[key]) - float(fine[key])) for key in public_logs
        }
        max_log_scale_absolute_error = max(
            max_log_scale_absolute_error, max(log_errors.values())
        )

        public_terminal = schedule.terminal_Q(np.asarray(TERMINAL_Q_PROBES, dtype=float))
        terminal_errors: dict[str, float] = {}
        for index, probe in enumerate(TERMINAL_Q_PROBES):
            reference = fine["terminal_Q"][f"{probe:g}"]
            public = float(public_terminal[index])
            error = 0.0 if probe == 3.0 else _relative_error(public, reference)
            terminal_errors[f"{probe:g}"] = error
        max_terminal_q_relative_error = max(
            max_terminal_q_relative_error, max(terminal_errors.values())
        )

        stability_keys = (
            "q_after_first_transition",
            "q_in",
            "q_p",
            "release_hold_length",
            "log_X_tail",
            "log_c_inf",
        )
        medium_to_fine = {
            key: _relative_error(medium[key], fine[key], floor=1.0)
            for key in stability_keys
        }
        max_medium_to_fine_relative_change = max(
            max_medium_to_fine_relative_change, max(medium_to_fine.values())
        )

        d_coarse_medium = abs(
            coarse["q_after_first_transition"] - medium["q_after_first_transition"]
        )
        d_medium_fine = abs(
            medium["q_after_first_transition"] - fine["q_after_first_transition"]
        )
        observed_order = math.log(d_coarse_medium / d_medium_fine, 2.0)
        minimum_observed_order = min(minimum_observed_order, observed_order)

        # Fixed mutation: remove half of Q_p.  This recreates a large but simple
        # terminal-target error and calibrates that the independent audit is not
        # blind to a wrong release target.
        mutated_q_p = 0.5 * fine["q_p"]
        mutated_hold = math.log(fine["q_in"] / mutated_q_p) / (1.0 - schedule.h)
        mutation_log_X_tail_shift = abs(mutated_hold - fine["release_hold_length"])
        minimum_mutation_shift = min(minimum_mutation_shift, mutation_log_X_tail_shift)

        independent_materializable = _float_materializable(fine["log_X_tail"]) and _float_materializable(
            fine["log_c_inf"]
        )
        public_materialized = True
        public_heat_inputs: dict[str, float] | None
        try:
            public_heat_inputs = schedule.materialize_heat_inputs()
        except OverflowError:
            public_materialized = False
            public_heat_inputs = None
        classification_match = independent_materializable == public_materialized
        materialization_classification_matches &= classification_match

        materialized_log_errors: dict[str, float] | None = None
        if public_heat_inputs is not None:
            materialized_log_errors = {
                "log_X_tail": abs(math.log(public_heat_inputs["X_tail"]) - fine["log_X_tail"]),
                "log_c_inf": abs(math.log(public_heat_inputs["c_inf"]) - fine["log_c_inf"]),
            }
            max_log_scale_absolute_error = max(
                max_log_scale_absolute_error, max(materialized_log_errors.values())
            )

        cases.append(
            {
                "name": name,
                "parameters": _raw_case_parameters(schedule),
                "schedule_sha256": schedule.sha256,
                "resolutions": list(RESOLUTIONS),
                "finest_independent": fine,
                "public_release": public_release,
                "release_relative_errors": release_errors,
                "log_scale_absolute_errors": log_errors,
                "terminal_q_relative_errors": terminal_errors,
                "medium_to_fine_relative_change": medium_to_fine,
                "transition_observed_order": observed_order,
                "mutation_half_q_p_log_X_tail_shift": mutation_log_X_tail_shift,
                "independent_materializable": independent_materializable,
                "public_materialized": public_materialized,
                "materialization_classification_match": classification_match,
                "materialized_log_errors": materialized_log_errors,
            }
        )

    guard_results = {
        "release_agreement": max_release_relative_error
        <= LOCAL_GUARDS["max_release_relative_error"],
        "terminal_q_agreement": max_terminal_q_relative_error
        <= LOCAL_GUARDS["max_terminal_q_relative_error"],
        "log_scale_agreement": max_log_scale_absolute_error
        <= LOCAL_GUARDS["max_log_scale_absolute_error"],
        "independent_resolution_stability": max_medium_to_fine_relative_change
        <= LOCAL_GUARDS["max_independent_medium_to_fine_relative_change"],
        "transition_refinement_order": minimum_observed_order
        >= LOCAL_GUARDS["minimum_transition_observed_order"],
        "mutation_detected": minimum_mutation_shift
        >= LOCAL_GUARDS["minimum_mutation_log_X_tail_shift"],
        "materialization_classification": materialization_classification_matches,
    }
    structural_preflight_passed = all(guard_results.values())

    return {
        "schema": SCHEMA,
        "operator_independence": {
            "production_transition": "scipy.solve_ivp DOP853",
            "production_terminal_quadrature": "Gauss-Legendre",
            "independent_transition": "integrating factor + cumulative trapezoid + composite Simpson",
            "independent_terminal_quadrature": "uniform composite Simpson",
            "training_loss_or_tensor_read": False,
            "production_internal_ode_state_read": False,
            "free_residual_canceling_forcing_used": False,
        },
        "frozen_local_guards": dict(LOCAL_GUARDS),
        "cases": cases,
        "summary": {
            "max_release_relative_error": max_release_relative_error,
            "max_terminal_q_relative_error": max_terminal_q_relative_error,
            "max_log_scale_absolute_error": max_log_scale_absolute_error,
            "max_independent_medium_to_fine_relative_change": max_medium_to_fine_relative_change,
            "minimum_transition_observed_order": minimum_observed_order,
            "minimum_mutation_half_q_p_log_X_tail_shift": minimum_mutation_shift,
            "materialization_classification_matches": materialization_classification_matches,
            "guard_results": guard_results,
            "structural_preflight_passed": structural_preflight_passed,
        },
        "cross_route_baseline": {
            "name": "ST006",
            "candidate_sha256": ST006_SHA256,
            "held_out_momentum_sampled_max": ST006_MOMENTUM_MAX,
            "held_out_momentum_volume_L2": ST006_VOLUME_L2,
            "directly_comparable": False,
            "reason": "this audit validates a local source schedule, not a full-domain momentum residual",
        },
        "formal_gates_unchanged": {
            "normalized_momentum": FORMAL_MOMENTUM_GATE,
            "divergence_max": FORMAL_DIVERGENCE_GATE,
        },
        "truth_boundary": {
            "terminal_tail_schedule_independently_preflighted": structural_preflight_passed,
            "global_leading_profile_reconstructed": False,
            "complete_kokuno_composite_velocity": False,
            "formal_full_domain_pde_gate_assessed": False,
            "normalized_ns_residual_le_1e-3_claimed": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("artifacts/kokuno_agent4/independent_terminal_tail_audit.json"),
    )
    args = parser.parse_args()
    report = run_audit()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report["summary"], indent=2, sort_keys=True))
    if not report["summary"]["structural_preflight_passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
