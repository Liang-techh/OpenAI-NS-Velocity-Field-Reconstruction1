"""One bounded Agent-3 mean-correction cycle on the actual Kokuno leading core.

This increment moves the Agent-3 correction machinery off the temporary capped-
bipolar bridge.  It consumes Agent-1's public source-native leading-core velocity
and Agent-2's existing complete-curl oscillatory correction, measures the real
phase-mean nonlinear defect with :class:`PhaseMeanDefectContract`, reconstructs
the compact theta/e=2 radial stress, and screens one bounded mean-velocity
correction on held-in data before freezing it for disjoint held-out evaluation.

KokunoYumeto's corrected 2026-09-09 reader supplies the weighted moment-
complement radial inverse used here.  The current mapping from reconstructed
radial stress to a mean swirl velocity is still the explicit autonomous
engineering adapter introduced by Agent 3 in PR #232; the reader does not state
``delta u_mean = stress``.  A previous bridge screen reduced the aggregate mean
defect while worsening its theta component, so this core screen adds a
pre-registered theta non-worsening guard.  No pressure or force is fitted and no
project acceptance threshold is changed.
"""
from __future__ import annotations

from dataclasses import replace
import argparse
import json
from math import pi
from pathlib import Path
from typing import Callable

import numpy as np

from .constrained_validation import residual
from .kokuno_complete_curl import KokunoCompleteCurlCorrection
from .kokuno_leading_core_series import KokunoLeadingCoreSeriesCandidate
from .kokuno_mean_correction_cycle import AMPLITUDE_GRID, StressShapedMeanSwirlCorrection
from .kokuno_mean_defect import PhaseMeanDefectContract, cylindrical_theta_component
from .kokuno_radial_stress import CompactRadialStressInverse, sample_real_phase_mean_radial_profiles


TASK = "KOKUNO-A3-CORE-MEAN-CORRECTION-CYCLE-004"
SOURCE_READER = "KokunoYumeto corrected 208-page reconstruction"
SOURCE_EDITION = "2026.09.09-consolidated"
SOURCE_DOI = "10.5281/zenodo.22678406"
SOURCE_RADIAL_FORMULA = "C12 / R6-R9 weighted moment-complement radial stress"
SOURCE_CYCLE_CONTEXT = "finite correction cycle: source stress -> signed covariance-amplitude map"
BASE_AGENT5_PR = 233
BASE_AGENT5_HEAD = "4b77580b2bedbbc9b047471cb77b2229549192f2"
OSCILLATORY_AMPLITUDE = 0.125
OSCILLATORY_PHASE = 0.0
TRAIN_SEED = 9_172_921
HELD_OUT_SEED = 9_172_927
DIVERGENCE_SEED = 9_172_933
TRAIN_COUNT = 12
HELD_OUT_COUNT = 18
DIVERGENCE_COUNT = 64
TRAIN_TIME = 0.5
HELD_OUT_TIMES = (0.375, 0.5, 0.625)
PROFILE_ANNULUS = (0.08, 0.36)
CORE_SAMPLE_RADIAL_BOUNDS = (0.09, 0.22)
CORE_SAMPLE_Z_HALF_WIDTH = 0.12
CORE_CORRECTION_AXIAL_HALF_WIDTH = 0.20
DEFAULT_RADIAL_COUNT = 33
DEFAULT_ANGULAR_COUNT = 8
PointVelocity = Callable[[np.ndarray, float], np.ndarray]


def _vector_rms(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(np.sum(values * values, axis=-1))))


def _vector_max(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    return float(np.max(np.linalg.norm(values, axis=-1)))


def _scalar_rms(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(values * values)))


def _zero_velocity(points: np.ndarray, time: float) -> np.ndarray:
    del time
    return np.zeros_like(points, dtype=float)


def _zero_force(points: np.ndarray, time: float) -> np.ndarray:
    del time
    return np.zeros_like(points, dtype=float)


def _sample_core_points(*, seed: int, count: int) -> np.ndarray:
    """Fresh cylindrical samples kept well inside the finite leading-core domain."""
    if count <= 0:
        raise ValueError("count must be positive")
    r0, r1 = CORE_SAMPLE_RADIAL_BOUNDS
    rng = np.random.default_rng(seed)
    radius = np.sqrt(rng.uniform(r0 * r0, r1 * r1, count))
    angle = rng.uniform(0.0, 2.0 * pi, count)
    z = rng.uniform(-CORE_SAMPLE_Z_HALF_WIDTH, CORE_SAMPLE_Z_HALF_WIDTH, count)
    return np.column_stack((radius * np.cos(angle), radius * np.sin(angle), z))


def _compose(base_velocity: PointVelocity, correction: StressShapedMeanSwirlCorrection) -> PointVelocity:
    def velocity(points: np.ndarray, time: float) -> np.ndarray:
        points = np.asarray(points, dtype=float)
        base = np.asarray(base_velocity(points, time), dtype=float)
        delta = np.asarray(correction.at_points(points, time), dtype=float)
        if base.shape != points.shape or delta.shape != points.shape:
            raise ValueError("base and mean correction must return finite (N,3) values")
        out = base + delta
        if not np.isfinite(out).all():
            raise RuntimeError("corrected leading-core velocity became non-finite")
        return out
    return velocity


def _increment_metrics(values: np.ndarray, points: np.ndarray) -> dict[str, float]:
    theta = cylindrical_theta_component(values, points)
    return {
        "rms": _vector_rms(values),
        "max": _vector_max(values),
        "theta_rms": _scalar_rms(theta),
        "theta_max": float(np.max(np.abs(theta))),
    }


def _pressure_at_points(leading: KokunoLeadingCoreSeriesCandidate) -> Callable[[np.ndarray, float], np.ndarray]:
    def pressure(points: np.ndarray, time: float) -> np.ndarray:
        points = np.asarray(points, dtype=float)
        return np.asarray(
            leading.pressure(points[:, 0], points[:, 1], points[:, 2], time),
            dtype=float,
        )
    return pressure


def _fixed_pressure_gradient(
    leading: KokunoLeadingCoreSeriesCandidate,
    points: np.ndarray,
    time: float,
    contract: PhaseMeanDefectContract,
) -> np.ndarray:
    """Evaluate only grad(p_leading) with the same independent FD operator.

    The phase-mean defect increment itself does not need pressure because the
    unchanged pressure cancels.  Adding this separately to the zero-pressure
    phase-mean operator gives a pressure-inclusive raw core diagnostic without
    re-running or fitting any pressure model.
    """
    result = residual(
        _zero_velocity,
        _pressure_at_points(leading),
        _zero_force,
        points,
        float(time),
        nu=contract.nu,
        step=contract.step,
        time_bounds=contract.time_bounds,
    )
    return np.asarray(result["momentum"], dtype=float)


def _assert_stencil_safe(
    leading: KokunoLeadingCoreSeriesCandidate,
    points: np.ndarray,
    times: tuple[float, ...],
    step: float,
) -> None:
    """Fail closed unless all spatial/time probes used by CR006 remain in core."""
    points = np.asarray(points, dtype=float)
    for time in times:
        for dt in (-step, 0.0, step):
            leading.coordinates(points[:, 0], points[:, 1], points[:, 2], time + dt)
        for axis in range(3):
            for multiple in (-2.0, -1.0, 1.0, 2.0):
                shifted = points.copy()
                shifted[:, axis] += multiple * step
                leading.coordinates(shifted[:, 0], shifted[:, 1], shifted[:, 2], time)


def _independent_correction_divergence(
    correction: StressShapedMeanSwirlCorrection,
    points: np.ndarray,
    *,
    step: float = 1.0e-5,
) -> np.ndarray:
    points = np.asarray(points, dtype=float)
    divergence = np.zeros(len(points), dtype=float)
    for axis in range(3):
        delta = np.eye(3)[axis] * step
        plus = correction.at_points(points + delta, TRAIN_TIME)
        minus = correction.at_points(points - delta, TRAIN_TIME)
        divergence += ((plus - minus) / (2.0 * step))[:, axis]
    return divergence


def _decision(
    *,
    training_ratio: float,
    held_out_ratio: float,
    theta_ratio: float,
    pressure_operator_ratio: float,
    correction_divergence_max: float,
    amplitude: float,
    stress_scale: float,
) -> tuple[bool, dict[str, bool]]:
    values = np.asarray(
        [training_ratio, held_out_ratio, theta_ratio, pressure_operator_ratio,
         correction_divergence_max, amplitude, stress_scale],
        dtype=float,
    )
    if not np.isfinite(values).all():
        raise ValueError("cycle decision inputs must be finite")
    checks = {
        "held_in_total_mean_defect_improves": bool(training_ratio < 1.0),
        "held_out_total_mean_defect_improves": bool(held_out_ratio < 1.0),
        "held_out_theta_mean_defect_does_not_worsen": bool(theta_ratio <= 1.0),
        "held_out_pressure_inclusive_raw_operator_nonworsening_0p1pct": bool(
            pressure_operator_ratio <= 1.001
        ),
        "independent_correction_divergence_max_le_1e-5": bool(
            correction_divergence_max <= 1.0e-5
        ),
        "correction_nontrivial": bool(abs(amplitude) > 0.0 and stress_scale > 1.0e-14),
    }
    return bool(all(checks.values())), checks


def _select_amplitude(
    leading: KokunoLeadingCoreSeriesCandidate,
    oscillatory: KokunoCompleteCurlCorrection,
    inverse: CompactRadialStressInverse,
    contract: PhaseMeanDefectContract,
    points: np.ndarray,
) -> tuple[float, dict[str, float], list[dict]]:
    baseline = contract.evaluate(leading.at_points, oscillatory, points, TRAIN_TIME)
    before = _increment_metrics(baseline.mean_defect_increment, points)
    rows: list[dict] = []
    for amplitude in AMPLITUDE_GRID:
        correction = StressShapedMeanSwirlCorrection(
            inverse=inverse,
            amplitude=float(amplitude),
            axial_half_width=CORE_CORRECTION_AXIAL_HALF_WIDTH,
        )
        corrected = _compose(leading.at_points, correction)
        sample = contract.evaluate(corrected, oscillatory, points, TRAIN_TIME)
        increment_relative_original_base = sample.mean_momentum - baseline.base_momentum
        metrics = _increment_metrics(increment_relative_original_base, points)
        rows.append({
            "amplitude": float(amplitude),
            **metrics,
            "ratio_to_uncorrected_mean_defect_rms": metrics["rms"] / max(
                before["rms"], np.finfo(float).tiny
            ),
            "theta_ratio_to_uncorrected": metrics["theta_rms"] / max(
                before["theta_rms"], np.finfo(float).tiny
            ),
            "composite_phase_divergence_max": float(np.max(np.abs(sample.phase_divergence))),
        })
    selected = min(rows, key=lambda row: row["rms"])
    return float(selected["amplitude"]), before, rows


def generate_core_mean_correction_cycle_report(
    *,
    output: str | Path = "artifacts/kokuno_agent3/core_mean_correction_cycle_report.json",
    radial_count: int = DEFAULT_RADIAL_COUNT,
    angular_count: int = DEFAULT_ANGULAR_COUNT,
    train_count: int = TRAIN_COUNT,
    held_out_count: int = HELD_OUT_COUNT,
) -> dict:
    """Run one bounded finite correction-cycle screen on the actual leading core."""
    if train_count <= 0 or held_out_count <= 0:
        raise ValueError("train_count and held_out_count must be positive")

    leading = KokunoLeadingCoreSeriesCandidate()
    oscillatory = KokunoCompleteCurlCorrection(
        amplitude=OSCILLATORY_AMPLITUDE,
        phase=OSCILLATORY_PHASE,
    )
    contract = PhaseMeanDefectContract(phase_count=8)

    r_inner, r_outer = PROFILE_ANNULUS
    theta_profile, axial_profile, projection = sample_real_phase_mean_radial_profiles(
        leading.at_points,
        oscillatory,
        contract,
        time=TRAIN_TIME,
        z=0.0,
        r_inner=r_inner,
        r_outer=r_outer,
        radial_count=radial_count,
        angular_count=angular_count,
    )
    inverse = CompactRadialStressInverse(theta_profile)

    train_points = _sample_core_points(seed=TRAIN_SEED, count=train_count)
    held_out_points = _sample_core_points(seed=HELD_OUT_SEED, count=held_out_count)
    _assert_stencil_safe(leading, train_points, (TRAIN_TIME,), contract.step)
    _assert_stencil_safe(leading, held_out_points, HELD_OUT_TIMES, contract.step)

    selected_amplitude, training_before, training_rows = _select_amplitude(
        leading, oscillatory, inverse, contract, train_points
    )
    selected_row = next(
        row for row in training_rows if row["amplitude"] == selected_amplitude
    )
    selected = StressShapedMeanSwirlCorrection(
        inverse=inverse,
        amplitude=selected_amplitude,
        axial_half_width=CORE_CORRECTION_AXIAL_HALF_WIDTH,
    )
    corrected_base = _compose(leading.at_points, selected)

    rows: list[dict] = []
    all_before: list[np.ndarray] = []
    all_after: list[np.ndarray] = []
    all_theta_before: list[np.ndarray] = []
    all_theta_after: list[np.ndarray] = []
    all_pressure_before: list[np.ndarray] = []
    all_pressure_after: list[np.ndarray] = []

    for time in HELD_OUT_TIMES:
        baseline = contract.evaluate(leading.at_points, oscillatory, held_out_points, time)
        after_sample = contract.evaluate(corrected_base, oscillatory, held_out_points, time)
        before_increment = baseline.mean_defect_increment
        after_increment = after_sample.mean_momentum - baseline.base_momentum
        before_metrics = _increment_metrics(before_increment, held_out_points)
        after_metrics = _increment_metrics(after_increment, held_out_points)

        grad_p = _fixed_pressure_gradient(leading, held_out_points, time, contract)
        pressure_before = baseline.mean_momentum + grad_p
        pressure_after = after_sample.mean_momentum + grad_p
        pressure_before_rms = _vector_rms(pressure_before)
        pressure_after_rms = _vector_rms(pressure_after)
        correction_values = selected.at_points(held_out_points, time)

        rows.append({
            "time": float(time),
            "mean_defect_before": before_metrics,
            "mean_defect_after": after_metrics,
            "mean_defect_rms_ratio": after_metrics["rms"] / max(
                before_metrics["rms"], np.finfo(float).tiny
            ),
            "theta_mean_defect_rms_ratio": after_metrics["theta_rms"] / max(
                before_metrics["theta_rms"], np.finfo(float).tiny
            ),
            "pressure_inclusive_raw_phase_mean_operator_rms_before": pressure_before_rms,
            "pressure_inclusive_raw_phase_mean_operator_rms_after": pressure_after_rms,
            "pressure_inclusive_raw_phase_mean_operator_rms_ratio": pressure_after_rms / max(
                pressure_before_rms, np.finfo(float).tiny
            ),
            "selected_correction_rms": _vector_rms(correction_values),
            "selected_correction_max": _vector_max(correction_values),
            "composite_phase_divergence_max_before": float(
                np.max(np.abs(baseline.phase_divergence))
            ),
            "composite_phase_divergence_max_after": float(
                np.max(np.abs(after_sample.phase_divergence))
            ),
        })
        all_before.append(before_increment)
        all_after.append(after_increment)
        all_theta_before.append(cylindrical_theta_component(before_increment, held_out_points))
        all_theta_after.append(cylindrical_theta_component(after_increment, held_out_points))
        all_pressure_before.append(pressure_before)
        all_pressure_after.append(pressure_after)

    before_all = np.concatenate(all_before, axis=0)
    after_all = np.concatenate(all_after, axis=0)
    theta_before_all = np.concatenate(all_theta_before, axis=0)
    theta_after_all = np.concatenate(all_theta_after, axis=0)
    pressure_before_all = np.concatenate(all_pressure_before, axis=0)
    pressure_after_all = np.concatenate(all_pressure_after, axis=0)

    held_out_before_rms = _vector_rms(before_all)
    held_out_after_rms = _vector_rms(after_all)
    held_out_ratio = held_out_after_rms / max(held_out_before_rms, np.finfo(float).tiny)
    theta_before_rms = _scalar_rms(theta_before_all)
    theta_after_rms = _scalar_rms(theta_after_all)
    theta_ratio = theta_after_rms / max(theta_before_rms, np.finfo(float).tiny)
    pressure_before_rms = _vector_rms(pressure_before_all)
    pressure_after_rms = _vector_rms(pressure_after_all)
    pressure_ratio = pressure_after_rms / max(pressure_before_rms, np.finfo(float).tiny)
    training_ratio = selected_row["ratio_to_uncorrected_mean_defect_rms"]

    divergence_points = _sample_core_points(seed=DIVERGENCE_SEED, count=DIVERGENCE_COUNT)
    correction_divergence = _independent_correction_divergence(selected, divergence_points)
    correction_divergence_rms = _scalar_rms(correction_divergence)
    correction_divergence_max = float(np.max(np.abs(correction_divergence)))

    accepted, checks = _decision(
        training_ratio=training_ratio,
        held_out_ratio=held_out_ratio,
        theta_ratio=theta_ratio,
        pressure_operator_ratio=pressure_ratio,
        correction_divergence_max=correction_divergence_max,
        amplitude=selected_amplitude,
        stress_scale=selected.stress_scale,
    )

    report = {
        "task": TASK,
        "base": {
            "agent5_pr": BASE_AGENT5_PR,
            "agent5_exact_head": BASE_AGENT5_HEAD,
            "leading_family": "KokunoLeadingCoreSeriesCandidate",
            "leading_sha256": leading.sha256,
            "core_only": True,
            "global_outer_heat_matching_complete": False,
        },
        "source": {
            "reader": SOURCE_READER,
            "edition": SOURCE_EDITION,
            "doi": SOURCE_DOI,
            "radial_formula": SOURCE_RADIAL_FORMULA,
            "cycle_context": SOURCE_CYCLE_CONTEXT,
            "source_exact_stress_to_velocity_lift": False,
            "source_cycle_gain_diagnostic_only": {
                "wave_residual_gain_exponent": 0.39999,
                "tangential_mean_gain_exponent": 0.17,
                "defect_gain_exponent": 0.89996,
                "used_as_numeric_contraction_gate": False,
            },
        },
        "oscillatory": {
            "metadata": oscillatory.metadata(),
            "amplitude_fixed_from_integrated_agent5_checkpoint": OSCILLATORY_AMPLITUDE,
            "phase_fixed_after_agent2_phase_screen": OSCILLATORY_PHASE,
            "phase_or_amplitude_optimized_here": False,
        },
        "real_defect_contract": {
            "nu": contract.nu,
            "derivative_step": contract.step,
            "phase_count": contract.phase_count,
            "time_bounds": list(contract.time_bounds),
            "pressure_for_increment": "unchanged Agent-1 pressure cancels in before/after mean-defect increment",
            "pressure_for_raw_operator": "Agent-1 public pressure gradient added back with the same independent FD operator",
            "force": "zero-force core diagnostic only; no force fitted or selected",
            "surrogate_defect_used": False,
        },
        "stress_shape_input": {
            "projection": projection,
            "theta_e2_metrics": inverse.metrics(),
            "axial_e1_raw_gate_capture_ratio": axial_profile.gate_capture_rms_ratio,
            "profile_annulus": list(PROFILE_ANNULUS),
            "profile_time": TRAIN_TIME,
            "profile_z": 0.0,
            "radial_count": radial_count,
            "angular_count": angular_count,
        },
        "engineering_mean_velocity_lift": selected.metadata(),
        "selection": {
            "objective": "held-in total phase-mean defect relative to the original leading-core operator",
            "amplitude_grid": list(AMPLITUDE_GRID),
            "amplitude_grid_widened_from_agent3_pr232": False,
            "zero_amplitude_is_baseline_not_selectable_success": True,
            "train_seed": TRAIN_SEED,
            "train_count": train_count,
            "train_time": TRAIN_TIME,
            "training_before": training_before,
            "training_rows": training_rows,
            "selected_amplitude": selected_amplitude,
            "selected_training_ratio": training_ratio,
            "selected_amplitude_at_grid_boundary": bool(
                abs(selected_amplitude) == max(abs(value) for value in AMPLITUDE_GRID)
            ),
        },
        "held_out": {
            "seed": HELD_OUT_SEED,
            "count": held_out_count,
            "times": list(HELD_OUT_TIMES),
            "points_reused_from_training": False,
            "rows": rows,
            "aggregate_mean_defect_rms_before": held_out_before_rms,
            "aggregate_mean_defect_rms_after": held_out_after_rms,
            "aggregate_mean_defect_rms_ratio": held_out_ratio,
            "aggregate_theta_mean_defect_rms_before": theta_before_rms,
            "aggregate_theta_mean_defect_rms_after": theta_after_rms,
            "aggregate_theta_mean_defect_rms_ratio": theta_ratio,
            "aggregate_pressure_inclusive_raw_phase_mean_operator_rms_before": pressure_before_rms,
            "aggregate_pressure_inclusive_raw_phase_mean_operator_rms_after": pressure_after_rms,
            "aggregate_pressure_inclusive_raw_phase_mean_operator_rms_ratio": pressure_ratio,
        },
        "correction_checks": {
            "divergence_seed": DIVERGENCE_SEED,
            "divergence_points": DIVERGENCE_COUNT,
            "divergence_fd_step": 1.0e-5,
            "correction_divergence_rms": correction_divergence_rms,
            "correction_divergence_max": correction_divergence_max,
            "correction_nontrivial": bool(
                abs(selected_amplitude) > 0.0 and selected.stress_scale > 1.0e-14
            ),
            "radial_support": [inverse.r_inner, inverse.r_outer],
            "axial_support": [-selected.axial_half_width, selected.axial_half_width],
        },
        "cycle_decision": {
            "accepted_for_next_cycle": accepted,
            "checks": checks,
            "new_theta_guard_motivated_by_pr232_failure": True,
            "if_rejected": "retain the failed core correction receipt and do not widen the amplitude grid in this run",
            "measured_total_mean_defect_contraction": float(1.0 - held_out_ratio),
            "measured_theta_mean_defect_contraction": float(1.0 - theta_ratio),
        },
        "truth_boundary": {
            "uses_actual_agent1_leading_core": True,
            "temporary_capped_bipolar_bridge_used": False,
            "consumes_real_phase_mean_defect": True,
            "surrogate_defect_used": False,
            "finite_correction_cycle_one_step_screen_run": True,
            "stress_to_velocity_lift_is_autonomous_engineering_adapter": True,
            "global_leading_profile_complete": False,
            "compatible_final_forcing_attached": False,
            "formal_full_domain_normalized_pde_gate_assessed": False,
            "normalized_ns_residual_le_1e-3_claimed": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
    }

    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run one bounded Agent-3 mean correction cycle on the actual Kokuno leading core"
    )
    parser.add_argument(
        "--output",
        default="artifacts/kokuno_agent3/core_mean_correction_cycle_report.json",
    )
    parser.add_argument("--radial-count", type=int, default=DEFAULT_RADIAL_COUNT)
    parser.add_argument("--angular-count", type=int, default=DEFAULT_ANGULAR_COUNT)
    parser.add_argument("--train-count", type=int, default=TRAIN_COUNT)
    parser.add_argument("--held-out-count", type=int, default=HELD_OUT_COUNT)
    args = parser.parse_args()
    report = generate_core_mean_correction_cycle_report(
        output=args.output,
        radial_count=args.radial_count,
        angular_count=args.angular_count,
        train_count=args.train_count,
        held_out_count=args.held_out_count,
    )
    print(json.dumps({
        "leading_sha256": report["base"]["leading_sha256"],
        "selected_amplitude": report["selection"]["selected_amplitude"],
        "held_out_mean_defect_ratio": report["held_out"]["aggregate_mean_defect_rms_ratio"],
        "held_out_theta_ratio": report["held_out"]["aggregate_theta_mean_defect_rms_ratio"],
        "pressure_inclusive_raw_operator_ratio": report["held_out"]["aggregate_pressure_inclusive_raw_phase_mean_operator_rms_ratio"],
        "correction_divergence_max": report["correction_checks"]["correction_divergence_max"],
        "accepted_for_next_cycle": report["cycle_decision"]["accepted_for_next_cycle"],
    }, indent=2))


if __name__ == "__main__":
    main()
