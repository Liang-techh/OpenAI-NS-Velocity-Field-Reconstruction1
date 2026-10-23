"""Independent three-stage Kokuno core-composite finite-difference audit.

This validator treats the current public leading velocity/pressure, oscillatory
velocity, and Agent-3 mean correction as black-box value providers.  The
momentum/divergence operator is Agent-4's separate fourth-order Cartesian FD
implementation from ``kokuno_independent_leading_core``; candidate-side
derivative helpers, training tensors, and training losses are not consumed.

The upstream Agent-3 correction screened in PR #241 is rejected, not promoted.
This module replays its fixed amplitude and construction only so the same
leading-only -> leading+oscillatory -> after-correction stages can be compared
on a fresh held-out sample under one independent operator.  The field remains
core-only, with zero forcing used only as a raw diagnostic, so the registered
full-domain 1e-3 PDE gate remains formally unassessed.
"""
from __future__ import annotations

import argparse
import json
from math import pi
from pathlib import Path
from typing import Any, Callable

import numpy as np

from .kokuno_complete_curl import KokunoCompleteCurlCorrection
from .kokuno_independent_leading_core import (
    NU,
    REGISTERED_DIVERGENCE_THRESHOLD,
    REGISTERED_PDE_THRESHOLD,
    REGISTERED_RESIDUAL_SCALE,
    evaluate_fd4,
)
from .kokuno_leading_core_series import KokunoLeadingCoreSeriesCandidate
from .kokuno_mean_correction_cycle import StressShapedMeanSwirlCorrection
from .kokuno_mean_defect import PhaseMeanDefectContract
from .kokuno_radial_stress import (
    CompactRadialStressInverse,
    sample_real_phase_mean_radial_profiles,
)

SCHEMA = "kokuno-independent-core-composite-v1"
TASK = "KOKUNO-VAL-CORE-COMPOSITE-CORRECTION-004"
BASE_AGENT3_PR = 241
BASE_AGENT3_HEAD = "23b4860f2d1ea514aed4a1cce71e9c06b7b596a3"
VALIDATION_SEED = 9_172_941
DEFAULT_TIMES = (0.375, 0.5, 0.625)
DEFAULT_STEPS = (0.02, 0.01, 0.005)
ANNULUS_RADIAL_BOUNDS = (0.10, 0.22)
ANNULUS_Z_HALF_WIDTH = 0.10
AXIS_NEAR_RADIAL_BOUNDS = (2.0e-4, 1.0e-2)
AXIS_NEAR_Z_HALF_WIDTH = 0.10
DEFAULT_ANNULUS_COUNT = 48
DEFAULT_AXIS_NEAR_COUNT = 24
OSCILLATORY_AMPLITUDE = 0.125
OSCILLATORY_PHASE = 0.0
UPSTREAM_REJECTED_MEAN_AMPLITUDE = 0.005
PROFILE_TIME = 0.5
PROFILE_ANNULUS = (0.08, 0.36)
PROFILE_RADIAL_COUNT = 33
PROFILE_ANGULAR_COUNT = 8
MEAN_CORRECTION_AXIAL_HALF_WIDTH = 0.20

VelocityCallable = Callable[[Any, Any, Any, Any], np.ndarray]


def _sample_validation_points(
    *,
    seed: int = VALIDATION_SEED,
    annulus_count: int = DEFAULT_ANNULUS_COUNT,
    axis_near_count: int = DEFAULT_AXIS_NEAR_COUNT,
) -> dict[str, np.ndarray]:
    if annulus_count < 8 or axis_near_count < 4:
        raise ValueError("held-out strata are too small")
    rng = np.random.default_rng(int(seed))

    r0, r1 = ANNULUS_RADIAL_BOUNDS
    radius = np.sqrt(rng.uniform(r0 * r0, r1 * r1, annulus_count))
    angle = rng.uniform(0.0, 2.0 * pi, annulus_count)
    annulus = np.column_stack(
        (
            radius * np.cos(angle),
            radius * np.sin(angle),
            rng.uniform(-ANNULUS_Z_HALF_WIDTH, ANNULUS_Z_HALF_WIDTH, annulus_count),
        )
    )

    a0, a1 = AXIS_NEAR_RADIAL_BOUNDS
    radius = rng.uniform(a0, a1, axis_near_count)
    angle = rng.uniform(0.0, 2.0 * pi, axis_near_count)
    axis_near = np.column_stack(
        (
            radius * np.cos(angle),
            radius * np.sin(angle),
            rng.uniform(-AXIS_NEAR_Z_HALF_WIDTH, AXIS_NEAR_Z_HALF_WIDTH, axis_near_count),
        )
    )
    return {"annulus": annulus, "axis_near": axis_near}


def _point_correction_as_array(
    correction: StressShapedMeanSwirlCorrection,
    x: Any,
    y: Any,
    z: Any,
    time: Any,
) -> np.ndarray:
    x, y, z, time = np.broadcast_arrays(
        np.asarray(x, dtype=float),
        np.asarray(y, dtype=float),
        np.asarray(z, dtype=float),
        np.asarray(time, dtype=float),
    )
    points = np.column_stack((x.ravel(), y.ravel(), z.ravel()))
    flat_time = time.ravel()
    if not np.all(flat_time == flat_time[0]):
        raise ValueError("mean-correction adapter requires a common evaluation time")
    values = correction.at_points(points, float(flat_time[0]))
    return np.asarray(values, dtype=float).reshape(x.shape + (3,))


def _compose_velocity(
    leading: KokunoLeadingCoreSeriesCandidate,
    oscillatory: KokunoCompleteCurlCorrection | None = None,
    mean_correction: StressShapedMeanSwirlCorrection | None = None,
) -> VelocityCallable:
    def velocity(x: Any, y: Any, z: Any, time: Any) -> np.ndarray:
        value = np.asarray(leading.velocity(x, y, z, time), dtype=float)
        if oscillatory is not None:
            value = value + np.asarray(oscillatory.velocity(x, y, z, time), dtype=float)
        if mean_correction is not None:
            value = value + _point_correction_as_array(mean_correction, x, y, z, time)
        if not np.all(np.isfinite(value)):
            raise RuntimeError("composite velocity became non-finite")
        return value

    return velocity


def _build_upstream_rejected_correction() -> tuple[
    KokunoLeadingCoreSeriesCandidate,
    KokunoCompleteCurlCorrection,
    CompactRadialStressInverse,
    StressShapedMeanSwirlCorrection,
    dict[str, float],
]:
    """Replay only the fixed PR-241 correction shape; do not rerun selection."""
    leading = KokunoLeadingCoreSeriesCandidate()
    oscillatory = KokunoCompleteCurlCorrection(
        amplitude=OSCILLATORY_AMPLITUDE,
        phase=OSCILLATORY_PHASE,
    )
    contract = PhaseMeanDefectContract(phase_count=8)
    theta_profile, _, projection = sample_real_phase_mean_radial_profiles(
        leading.at_points,
        oscillatory,
        contract,
        time=PROFILE_TIME,
        z=0.0,
        r_inner=PROFILE_ANNULUS[0],
        r_outer=PROFILE_ANNULUS[1],
        radial_count=PROFILE_RADIAL_COUNT,
        angular_count=PROFILE_ANGULAR_COUNT,
    )
    inverse = CompactRadialStressInverse(theta_profile)
    correction = StressShapedMeanSwirlCorrection(
        inverse=inverse,
        amplitude=UPSTREAM_REJECTED_MEAN_AMPLITUDE,
        axial_half_width=MEAN_CORRECTION_AXIAL_HALF_WIDTH,
    )
    projection_summary = {
        "full_mean_defect_increment_rms": float(projection["full_mean_defect_increment_rms"]),
        "ring_theta_raw_rms": float(projection["ring_theta_raw_rms"]),
        "theta_gate_capture_rms_ratio": float(projection["theta_gate_capture_rms_ratio"]),
    }
    return leading, oscillatory, inverse, correction, projection_summary


def _cylindrical_component_rms(residual: np.ndarray, points: np.ndarray) -> dict[str, float]:
    residual = np.asarray(residual, dtype=float)
    points = np.asarray(points, dtype=float)
    radius = np.hypot(points[:, 0], points[:, 1])
    if np.any(radius <= 0.0):
        raise ValueError("cylindrical component audit requires non-axis points")
    radial = (residual[:, 0] * points[:, 0] + residual[:, 1] * points[:, 1]) / radius
    theta = (-residual[:, 0] * points[:, 1] + residual[:, 1] * points[:, 0]) / radius
    axial = residual[:, 2]
    return {
        "radial_rms": float(np.sqrt(np.mean(radial * radial))),
        "theta_rms": float(np.sqrt(np.mean(theta * theta))),
        "axial_rms": float(np.sqrt(np.mean(axial * axial))),
    }


def _metrics(
    values: dict[str, np.ndarray],
    points: np.ndarray,
    *,
    volume: float | None,
) -> dict[str, float | dict[str, float]]:
    residual = np.asarray(values["residual"], dtype=float)
    divergence = np.asarray(values["divergence"], dtype=float)
    velocity = np.asarray(values["velocity"], dtype=float)
    residual_norm = np.linalg.norm(residual, axis=1)
    speed = np.linalg.norm(velocity, axis=1)
    residual_rms = float(np.sqrt(np.mean(residual_norm * residual_norm)))
    result: dict[str, float | dict[str, float]] = {
        "residual_max": float(np.max(residual_norm)),
        "residual_rms": residual_rms,
        "normalized_residual_max": float(np.max(residual_norm) / REGISTERED_RESIDUAL_SCALE),
        "normalized_residual_rms": float(residual_rms / REGISTERED_RESIDUAL_SCALE),
        "divergence_max": float(np.max(np.abs(divergence))),
        "divergence_rms": float(np.sqrt(np.mean(divergence * divergence))),
        "velocity_rms": float(np.sqrt(np.mean(speed * speed))),
        "velocity_max": float(np.max(speed)),
        "cylindrical_residual_component_rms": _cylindrical_component_rms(residual, points),
    }
    if volume is not None:
        result["volume_l2_residual_estimate"] = float(
            np.sqrt(float(volume) * np.mean(residual_norm * residual_norm))
        )
        result["volume_l2_divergence_estimate"] = float(
            np.sqrt(float(volume) * np.mean(divergence * divergence))
        )
    return result


def _zero_pressure(x: Any, y: Any, z: Any, time: Any) -> np.ndarray:
    x, y, z, time = np.broadcast_arrays(
        np.asarray(x, dtype=float),
        np.asarray(y, dtype=float),
        np.asarray(z, dtype=float),
        np.asarray(time, dtype=float),
    )
    return np.zeros_like(x, dtype=float)


def _mutation_calibration(
    velocity: VelocityCallable,
    pressure: Callable[[Any, Any, Any, Any], np.ndarray],
    points: np.ndarray,
    *,
    time: float = 0.5,
    step: float = 0.005,
    slope: float = 0.02,
) -> dict[str, float]:
    baseline = evaluate_fd4(velocity, pressure, points, time, step)

    def mutated(x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        values = np.array(velocity(x, y, z, t), copy=True)
        values[..., 0] += slope * np.asarray(x, dtype=float)
        return values

    changed = evaluate_fd4(mutated, pressure, points, time, step)
    delta = changed["divergence"] - baseline["divergence"]
    return {
        "injected_ux_linear_slope": float(slope),
        "mean_divergence_shift": float(np.mean(delta)),
        "max_abs_error_from_expected_shift": float(np.max(np.abs(delta - slope))),
    }


def _stage_rows(
    *,
    stages: dict[str, VelocityCallable],
    pressure: Callable[[Any, Any, Any, Any], np.ndarray],
    annulus: np.ndarray,
    axis_near: np.ndarray,
    times: tuple[float, ...],
    steps: tuple[float, ...],
) -> list[dict[str, Any]]:
    all_points = np.vstack((annulus, axis_near))
    n_annulus = len(annulus)
    r0, r1 = ANNULUS_RADIAL_BOUNDS
    annulus_volume = pi * (r1 * r1 - r0 * r0) * (2.0 * ANNULUS_Z_HALF_WIDTH)
    rows: list[dict[str, Any]] = []
    for stage_name, velocity in stages.items():
        for time in times:
            for step in steps:
                values = evaluate_fd4(
                    velocity,
                    pressure,
                    all_points,
                    float(time),
                    float(step),
                    nu=NU,
                )
                annulus_values = {
                    key: np.asarray(value)[:n_annulus] for key, value in values.items()
                }
                axis_values = {
                    key: np.asarray(value)[n_annulus:] for key, value in values.items()
                }
                rows.append(
                    {
                        "stage": stage_name,
                        "time": float(time),
                        "step": float(step),
                        "annulus": _metrics(annulus_values, annulus, volume=annulus_volume),
                        "axis_near": _metrics(axis_values, axis_near, volume=None),
                        "all_points": _metrics(values, all_points, volume=None),
                    }
                )
    return rows


def _finest_summary(
    rows: list[dict[str, Any]],
    *,
    stages: tuple[str, ...],
    times: tuple[float, ...],
    finest_step: float,
) -> dict[str, Any]:
    summary: dict[str, Any] = {}
    for stage in stages:
        selected = [
            row for row in rows if row["stage"] == stage and row["step"] == finest_step
        ]
        if len(selected) != len(times):
            raise RuntimeError("incomplete finest-step stage rows")
        rms = np.asarray([row["all_points"]["residual_rms"] for row in selected], dtype=float)
        maxima = np.asarray([row["all_points"]["residual_max"] for row in selected], dtype=float)
        divergence = np.asarray([row["all_points"]["divergence_max"] for row in selected], dtype=float)
        axial = np.asarray(
            [row["all_points"]["cylindrical_residual_component_rms"]["axial_rms"] for row in selected],
            dtype=float,
        )
        radial = np.asarray(
            [row["all_points"]["cylindrical_residual_component_rms"]["radial_rms"] for row in selected],
            dtype=float,
        )
        theta = np.asarray(
            [row["all_points"]["cylindrical_residual_component_rms"]["theta_rms"] for row in selected],
            dtype=float,
        )
        summary[stage] = {
            "mean_residual_rms_over_times": float(np.mean(rms)),
            "worst_residual_rms": float(np.max(rms)),
            "worst_residual_max": float(np.max(maxima)),
            "worst_normalized_residual_max": float(np.max(maxima) / REGISTERED_RESIDUAL_SCALE),
            "worst_divergence_max": float(np.max(divergence)),
            "worst_component_rms": {
                "radial": float(np.max(radial)),
                "theta": float(np.max(theta)),
                "axial": float(np.max(axial)),
            },
        }

    lead = summary["leading_only"]["mean_residual_rms_over_times"]
    osc = summary["leading_plus_oscillatory"]["mean_residual_rms_over_times"]
    corrected = summary["after_rejected_mean_correction"]["mean_residual_rms_over_times"]
    summary["stage_ratios"] = {
        "leading_plus_oscillatory_over_leading": float(osc / lead),
        "after_correction_over_leading_plus_oscillatory": float(corrected / osc),
        "after_correction_over_leading": float(corrected / lead),
    }
    return summary


def _resolution_sensitivity(
    rows: list[dict[str, Any]],
    *,
    stages: tuple[str, ...],
    times: tuple[float, ...],
    medium_step: float,
    finest_step: float,
) -> dict[str, dict[str, float]]:
    report: dict[str, dict[str, float]] = {}
    for stage in stages:
        medium = [row for row in rows if row["stage"] == stage and row["step"] == medium_step]
        finest = [row for row in rows if row["stage"] == stage and row["step"] == finest_step]
        if len(medium) != len(times) or len(finest) != len(times):
            raise RuntimeError("incomplete resolution ladder")
        medium_rms = float(np.mean([row["all_points"]["residual_rms"] for row in medium]))
        finest_rms = float(np.mean([row["all_points"]["residual_rms"] for row in finest]))
        medium_max = float(max(row["all_points"]["residual_max"] for row in medium))
        finest_max = float(max(row["all_points"]["residual_max"] for row in finest))
        report[stage] = {
            "medium_step": float(medium_step),
            "finest_step": float(finest_step),
            "mean_rms_relative_change": float(
                abs(finest_rms - medium_rms) / max(abs(finest_rms), np.finfo(float).tiny)
            ),
            "worst_max_relative_change": float(
                abs(finest_max - medium_max) / max(abs(finest_max), np.finfo(float).tiny)
            ),
        }
    return report


def _parameter_perturbation(
    *,
    leading: KokunoLeadingCoreSeriesCandidate,
    oscillatory: KokunoCompleteCurlCorrection,
    inverse: CompactRadialStressInverse,
    pressure: Callable[[Any, Any, Any, Any], np.ndarray],
    points: np.ndarray,
    time: float,
    step: float,
) -> dict[str, Any]:
    rows = []
    for amplitude in (0.0045, 0.005, 0.0055):
        correction = StressShapedMeanSwirlCorrection(
            inverse=inverse,
            amplitude=float(amplitude),
            axial_half_width=MEAN_CORRECTION_AXIAL_HALF_WIDTH,
        )
        velocity = _compose_velocity(leading, oscillatory, correction)
        values = evaluate_fd4(velocity, pressure, points, float(time), float(step), nu=NU)
        rms = float(
            np.sqrt(np.mean(np.sum(np.asarray(values["residual"], dtype=float) ** 2, axis=1)))
        )
        rows.append({"amplitude": float(amplitude), "residual_rms": rms})
    center = next(row["residual_rms"] for row in rows if row["amplitude"] == 0.005)
    for row in rows:
        row["ratio_to_upstream_amplitude"] = float(row["residual_rms"] / center)
    return {
        "kind": "fixed_plus_minus_10pct_mean_correction_amplitude_probe",
        "selection_performed": False,
        "time": float(time),
        "step": float(step),
        "rows": rows,
    }


def generate_core_composite_independent_report(
    *,
    output: str | Path = "artifacts/kokuno_agent4/core_composite_independent_report.json",
    annulus_count: int = DEFAULT_ANNULUS_COUNT,
    axis_near_count: int = DEFAULT_AXIS_NEAR_COUNT,
    seed: int = VALIDATION_SEED,
    times: tuple[float, ...] = DEFAULT_TIMES,
    steps: tuple[float, ...] = DEFAULT_STEPS,
) -> dict[str, Any]:
    if len(steps) < 3:
        raise ValueError("independent report requires at least three resolutions")
    if len(times) == 0:
        raise ValueError("at least one held-out time is required")
    if min(steps) <= 0.0 or not np.all(np.isfinite(steps)):
        raise ValueError("steps must be finite and positive")

    leading, oscillatory, inverse, correction, projection = _build_upstream_rejected_correction()
    samples = _sample_validation_points(
        seed=seed,
        annulus_count=annulus_count,
        axis_near_count=axis_near_count,
    )
    annulus = samples["annulus"]
    axis_near = samples["axis_near"]
    all_points = np.vstack((annulus, axis_near))

    leading_only = _compose_velocity(leading)
    leading_plus_oscillatory = _compose_velocity(leading, oscillatory)
    corrected = _compose_velocity(leading, oscillatory, correction)
    stages = {
        "leading_only": leading_only,
        "leading_plus_oscillatory": leading_plus_oscillatory,
        "after_rejected_mean_correction": corrected,
    }

    largest = float(max(steps))
    for time in times:
        leading.velocity(all_points[:, 0], all_points[:, 1], all_points[:, 2], float(time))
        for axis in range(3):
            for multiple in (-2.0, 2.0):
                shifted = np.array(all_points, copy=True)
                shifted[:, axis] += multiple * largest
                leading.velocity(shifted[:, 0], shifted[:, 1], shifted[:, 2], float(time))
        for multiple in (-2.0, 2.0):
            leading.velocity(
                all_points[:, 0],
                all_points[:, 1],
                all_points[:, 2],
                float(time) + multiple * largest,
            )

    rows = _stage_rows(
        stages=stages,
        pressure=leading.pressure,
        annulus=annulus,
        axis_near=axis_near,
        times=times,
        steps=steps,
    )
    finest = float(min(steps))
    sorted_steps = sorted(float(value) for value in steps)
    medium = sorted_steps[1]
    finest_summary = _finest_summary(
        rows,
        stages=tuple(stages),
        times=times,
        finest_step=finest,
    )
    resolution = _resolution_sensitivity(
        rows,
        stages=tuple(stages),
        times=times,
        medium_step=medium,
        finest_step=finest,
    )

    correction_values = correction.at_points(annulus, float(times[0]))
    oscillatory_values = oscillatory.velocity(
        annulus[:, 0], annulus[:, 1], annulus[:, 2], float(times[0])
    )
    correction_divergence = evaluate_fd4(
        lambda x, y, z, t: _point_correction_as_array(correction, x, y, z, t),
        _zero_pressure,
        annulus,
        float(times[0]),
        finest,
        nu=0.0,
    )["divergence"]
    outside_points = np.asarray(
        [
            [0.40, 0.0, 0.0],
            [0.10, 0.0, 0.25],
            [-0.40, 0.0, 0.0],
            [0.12, 0.0, -0.25],
        ],
        dtype=float,
    )
    outside_correction = correction.at_points(outside_points, float(times[0]))

    parameter_perturbation = _parameter_perturbation(
        leading=leading,
        oscillatory=oscillatory,
        inverse=inverse,
        pressure=leading.pressure,
        points=annulus,
        time=0.5,
        step=finest,
    )
    mutation_points = annulus[: min(12, len(annulus))]
    mutation = _mutation_calibration(
        corrected,
        leading.pressure,
        mutation_points,
        time=0.5,
        step=finest,
    )

    worst_local_residual = max(
        value["worst_normalized_residual_max"]
        for key, value in finest_summary.items()
        if key != "stage_ratios"
    )
    worst_local_divergence = max(
        value["worst_divergence_max"]
        for key, value in finest_summary.items()
        if key != "stage_ratios"
    )

    report: dict[str, Any] = {
        "schema": SCHEMA,
        "task_id": TASK,
        "base": {
            "agent3_pr": BASE_AGENT3_PR,
            "agent3_exact_head": BASE_AGENT3_HEAD,
            "leading_sha256": leading.sha256,
            "core_only": True,
            "global_outer_heat_matching_complete": False,
        },
        "operator": {
            "kind": "independent_fourth_order_centered_cartesian_finite_difference",
            "implementation": "kokuno_independent_leading_core.evaluate_fd4",
            "candidate_derivative_helpers_used": False,
            "training_loss_used": False,
            "training_tensors_used": False,
            "nu": NU,
            "forcing": "zero_only_for_raw_core_diagnostic_not_a_fitted_or_final_force",
            "registered_residual_scale": REGISTERED_RESIDUAL_SCALE,
            "registered_pde_threshold": REGISTERED_PDE_THRESHOLD,
            "registered_divergence_threshold": REGISTERED_DIVERGENCE_THRESHOLD,
            "steps": [float(value) for value in steps],
            "times": [float(value) for value in times],
        },
        "held_out": {
            "seed": int(seed),
            "annulus_count": int(annulus_count),
            "axis_near_count": int(axis_near_count),
            "points_reused_from_agent3_training_or_holdout": False,
            "annulus_radial_bounds": list(ANNULUS_RADIAL_BOUNDS),
            "annulus_z_half_width": ANNULUS_Z_HALF_WIDTH,
            "axis_near_radial_bounds": list(AXIS_NEAR_RADIAL_BOUNDS),
            "axis_near_z_half_width": AXIS_NEAR_Z_HALF_WIDTH,
        },
        "upstream_stage_replay": {
            "oscillatory_amplitude": OSCILLATORY_AMPLITUDE,
            "oscillatory_phase": OSCILLATORY_PHASE,
            "mean_correction_amplitude": UPSTREAM_REJECTED_MEAN_AMPLITUDE,
            "mean_correction_upstream_status": "rejected_by_agent3_pr241_not_promoted",
            "mean_correction_radial_support": list(PROFILE_ANNULUS),
            "mean_correction_axial_half_width": MEAN_CORRECTION_AXIAL_HALF_WIDTH,
            "mean_correction_stress_scale": correction.stress_scale,
            "stress_shape_projection_receipt": projection,
            "selection_rerun_here": False,
        },
        "rows": rows,
        "finest_step_summary": finest_summary,
        "resolution_sensitivity": resolution,
        "support_and_nontriviality": {
            "oscillatory_rms_on_annulus": float(
                np.sqrt(np.mean(np.sum(oscillatory_values**2, axis=1)))
            ),
            "mean_correction_rms_on_annulus": float(
                np.sqrt(np.mean(np.sum(correction_values**2, axis=1)))
            ),
            "mean_correction_max_on_annulus": float(
                np.max(np.linalg.norm(correction_values, axis=1))
            ),
            "mean_correction_fd4_divergence_max": float(
                np.max(np.abs(correction_divergence))
            ),
            "mean_correction_outside_declared_support_max": float(
                np.max(np.abs(outside_correction))
            ),
            "global_leading_support_validated": False,
        },
        "parameter_perturbation": parameter_perturbation,
        "mutation_calibration": mutation,
        "gate": {
            "formal_full_domain_gate_assessed": False,
            "formal_full_domain_gate_passed": False,
            "local_core_reference_threshold": REGISTERED_PDE_THRESHOLD,
            "local_worst_normalized_residual_max": float(worst_local_residual),
            "local_reference_threshold_met": bool(
                worst_local_residual <= REGISTERED_PDE_THRESHOLD
            ),
            "local_worst_divergence_max": float(worst_local_divergence),
            "local_divergence_reference_met": bool(
                worst_local_divergence <= REGISTERED_DIVERGENCE_THRESHOLD
            ),
            "reason_formal_gate_unassessed": (
                "core-only leading field; no completed core-to-heat/exterior matching and no compatible final fixed/restricted forcing"
            ),
        },
        "truth_boundary": {
            "complete_kokuno_composite_velocity": False,
            "after_correction_cycle_independently_preflighted": True,
            "global_support_validated": False,
            "formal_normalized_pde_gate_assessed": False,
            "normalized_ns_residual_le_1e-3_claimed": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
    }

    output_path = Path(output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        default="artifacts/kokuno_agent4/core_composite_independent_report.json",
    )
    parser.add_argument("--annulus-count", type=int, default=DEFAULT_ANNULUS_COUNT)
    parser.add_argument("--axis-near-count", type=int, default=DEFAULT_AXIS_NEAR_COUNT)
    args = parser.parse_args()
    report = generate_core_composite_independent_report(
        output=args.output,
        annulus_count=args.annulus_count,
        axis_near_count=args.axis_near_count,
    )
    summary = report["finest_step_summary"]
    print(
        json.dumps(
            {
                "leading_mean_rms": summary["leading_only"]["mean_residual_rms_over_times"],
                "oscillatory_mean_rms": summary["leading_plus_oscillatory"]["mean_residual_rms_over_times"],
                "after_correction_mean_rms": summary["after_rejected_mean_correction"]["mean_residual_rms_over_times"],
                "stage_ratios": summary["stage_ratios"],
                "local_worst_normalized_max": report["gate"]["local_worst_normalized_residual_max"],
                "formal_gate_assessed": report["gate"]["formal_full_domain_gate_assessed"],
                "pde_validated": report["truth_boundary"]["pde_validated"],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
