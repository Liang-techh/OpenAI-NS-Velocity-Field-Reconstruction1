"""Independent black-box audit of Agent-3's signed complete-curl correction.

This module deliberately avoids the CR006 residual operator used by Agent 3.
It instantiates the production correction once, then treats leading, oscillatory,
and signed-correction velocities as public black boxes and reconstructs
Navier--Stokes momentum/divergence with a separate second-order Cartesian FD
operator on fresh off-grid samples and three derivative resolutions.
"""
from __future__ import annotations

import argparse
import json
from math import pi
from pathlib import Path
from typing import Any, Callable

import numpy as np

from .kokuno_signed_amplitude_curl_cycle import _build_actual_signed_correction


TASK = "KOKUNO-A4-INDEPENDENT-SIGNED-CURL-CYCLE-006"
NU = 0.01
RESIDUAL_SCALE = 1.0
RESIDUAL_REFERENCE = 1.0e-3
DIVERGENCE_REFERENCE = 1.0e-5
SEED = 9_173_041
CORRECTION_SEED = 9_173_047
ANNULUS_COUNT = 32
AXIS_COUNT = 16
TIMES = (0.375, 0.5, 0.625)
STEPS = (0.004, 0.002, 0.001)
ANNULUS_R = (0.10, 0.20)
ANNULUS_Z_HALF = 0.12
AXIS_R = (0.002, 0.008)
AXIS_Z_ABS = (0.04, 0.16)
CORRECTION_SCALE_PROBES = (0.9, 1.0, 1.1)
PointVelocity = Callable[[np.ndarray, float], np.ndarray]
PointPressure = Callable[[np.ndarray, float], np.ndarray]


def _velocity_checked(fn: PointVelocity, points: np.ndarray, time: float) -> np.ndarray:
    values = np.asarray(fn(points, float(time)), dtype=float)
    if values.shape != points.shape or not np.isfinite(values).all():
        raise RuntimeError("velocity callable must return finite (N,3)")
    return values


def _pressure_checked(fn: PointPressure, points: np.ndarray, time: float) -> np.ndarray:
    values = np.asarray(fn(points, float(time)), dtype=float)
    if values.shape != (len(points),) or not np.isfinite(values).all():
        raise RuntimeError("pressure callable must return finite (N,)")
    return values


def evaluate_fd2(
    velocity: PointVelocity,
    pressure: PointPressure,
    points: np.ndarray,
    time: float,
    *,
    nu: float = NU,
    step: float,
) -> dict[str, np.ndarray]:
    """Independent second-order Cartesian momentum/divergence operator."""
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3 or not np.isfinite(points).all():
        raise ValueError("points must be a finite (N,3) array")
    if not np.isfinite(time) or not np.isfinite(nu) or nu < 0.0:
        raise ValueError("time and nonnegative viscosity must be finite")
    if not np.isfinite(step) or step <= 0.0:
        raise ValueError("step must be positive and finite")

    u0 = _velocity_checked(velocity, points, time)
    ut = (
        _velocity_checked(velocity, points, time + step)
        - _velocity_checked(velocity, points, time - step)
    ) / (2.0 * step)

    gradient = np.empty((len(points), 3, 3), dtype=float)
    laplacian = np.zeros_like(u0)
    pressure_gradient = np.empty_like(u0)
    eye = np.eye(3)
    for axis in range(3):
        plus_points = points + eye[axis] * step
        minus_points = points - eye[axis] * step
        u_plus = _velocity_checked(velocity, plus_points, time)
        u_minus = _velocity_checked(velocity, minus_points, time)
        gradient[:, :, axis] = (u_plus - u_minus) / (2.0 * step)
        laplacian += (u_plus - 2.0 * u0 + u_minus) / (step * step)
        p_plus = _pressure_checked(pressure, plus_points, time)
        p_minus = _pressure_checked(pressure, minus_points, time)
        pressure_gradient[:, axis] = (p_plus - p_minus) / (2.0 * step)

    convection = np.einsum("nj,nij->ni", u0, gradient)
    residual = ut + convection + pressure_gradient - nu * laplacian
    divergence = np.einsum("nii->n", gradient)
    return {
        "velocity": u0,
        "time_derivative": ut,
        "gradient": gradient,
        "laplacian": laplacian,
        "pressure_gradient": pressure_gradient,
        "convection": convection,
        "residual": residual,
        "divergence": divergence,
    }


def divergence_fd2(velocity: PointVelocity, points: np.ndarray, time: float, step: float) -> np.ndarray:
    points = np.asarray(points, dtype=float)
    divergence = np.zeros(len(points), dtype=float)
    eye = np.eye(3)
    for axis in range(3):
        plus = _velocity_checked(velocity, points + eye[axis] * step, time)
        minus = _velocity_checked(velocity, points - eye[axis] * step, time)
        divergence += ((plus - minus) / (2.0 * step))[:, axis]
    return divergence


def _vector_rms(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(np.sum(values * values, axis=-1))))


def _vector_max(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    return float(np.max(np.linalg.norm(values, axis=-1)))


def _scalar_rms(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(values * values)))


def _cylindrical_components(values: np.ndarray, points: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    points = np.asarray(points, dtype=float)
    radius = np.hypot(points[:, 0], points[:, 1])
    if np.any(radius <= 0.0):
        raise ValueError("cylindrical diagnostic requires r>0")
    radial = (values[:, 0] * points[:, 0] + values[:, 1] * points[:, 1]) / radius
    theta = (-values[:, 0] * points[:, 1] + values[:, 1] * points[:, 0]) / radius
    return np.column_stack((radial, theta, values[:, 2]))


def _sample_annulus(seed: int, count: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    r0, r1 = ANNULUS_R
    radius = np.sqrt(rng.uniform(r0 * r0, r1 * r1, count))
    angle = rng.uniform(0.0, 2.0 * pi, count)
    z = rng.uniform(-ANNULUS_Z_HALF, ANNULUS_Z_HALF, count)
    return np.column_stack((radius * np.cos(angle), radius * np.sin(angle), z))


def _sample_axis_near(seed: int, count: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    r0, r1 = AXIS_R
    radius = np.sqrt(rng.uniform(r0 * r0, r1 * r1, count))
    angle = rng.uniform(0.0, 2.0 * pi, count)
    zabs = rng.uniform(AXIS_Z_ABS[0], AXIS_Z_ABS[1], count)
    sign = rng.choice(np.asarray([-1.0, 1.0]), size=count)
    return np.column_stack((radius * np.cos(angle), radius * np.sin(angle), sign * zabs))


def _stratum_volume(name: str) -> float:
    if name == "annulus":
        r0, r1 = ANNULUS_R
        return float(pi * (r1 * r1 - r0 * r0) * (2.0 * ANNULUS_Z_HALF))
    if name == "axis_near":
        r0, r1 = AXIS_R
        z_length = 2.0 * (AXIS_Z_ABS[1] - AXIS_Z_ABS[0])
        return float(pi * (r1 * r1 - r0 * r0) * z_length)
    raise KeyError(name)


def _stage_metrics(result: dict[str, np.ndarray], points: np.ndarray, volume: float) -> dict[str, Any]:
    residual = result["residual"]
    divergence = result["divergence"]
    cylindrical = _cylindrical_components(residual, points)
    residual_rms = _vector_rms(residual)
    return {
        "sample_rms": residual_rms,
        "sample_max": _vector_max(residual),
        "normalized_sample_rms": residual_rms / RESIDUAL_SCALE,
        "normalized_sample_max": _vector_max(residual) / RESIDUAL_SCALE,
        "stratum_volume_l2": float(np.sqrt(volume * np.mean(np.sum(residual * residual, axis=1)))),
        "divergence_rms": _scalar_rms(divergence),
        "divergence_max": float(np.max(np.abs(divergence))),
        "radial_rms": _scalar_rms(cylindrical[:, 0]),
        "theta_rms": _scalar_rms(cylindrical[:, 1]),
        "axial_rms": _scalar_rms(cylindrical[:, 2]),
        "velocity_rms": _vector_rms(result["velocity"]),
    }


def _aggregate(rows: list[dict[str, Any]], stage: str, step: float) -> dict[str, float]:
    selected = [r for r in rows if r["stage"] == stage and r["step"] == step]
    if not selected:
        raise RuntimeError("missing aggregate rows")
    sample_weight = sum(int(r["count"]) for r in selected)
    rms_sq = sum(int(r["count"]) * float(r["metrics"]["sample_rms"]) ** 2 for r in selected) / sample_weight
    div_rms_sq = sum(int(r["count"]) * float(r["metrics"]["divergence_rms"]) ** 2 for r in selected) / sample_weight
    return {
        "sample_rms": float(np.sqrt(rms_sq)),
        "sample_max": float(max(r["metrics"]["sample_max"] for r in selected)),
        "normalized_sample_rms": float(np.sqrt(rms_sq) / RESIDUAL_SCALE),
        "normalized_sample_max": float(max(r["metrics"]["sample_max"] for r in selected) / RESIDUAL_SCALE),
        "stratified_time_sum_volume_l2": float(np.sqrt(sum(float(r["metrics"]["stratum_volume_l2"]) ** 2 for r in selected))),
        "divergence_rms": float(np.sqrt(div_rms_sq)),
        "divergence_max": float(max(r["metrics"]["divergence_max"] for r in selected)),
    }


def _mutation_calibration(step: float) -> dict[str, float]:
    points = _sample_annulus(SEED + 91, 8)

    def zero_velocity(x: np.ndarray, time: float) -> np.ndarray:
        del time
        return np.zeros_like(x)

    def slope_velocity(x: np.ndarray, time: float) -> np.ndarray:
        del time
        out = np.zeros_like(x)
        out[:, 0] = 0.02 * x[:, 0]
        return out

    def zero_pressure(x: np.ndarray, time: float) -> np.ndarray:
        del time
        return np.zeros(len(x))

    def slope_pressure(x: np.ndarray, time: float) -> np.ndarray:
        del time
        return 0.02 * x[:, 0]

    base = evaluate_fd2(zero_velocity, zero_pressure, points, 0.5, nu=0.0, step=step)
    vel_mut = evaluate_fd2(slope_velocity, zero_pressure, points, 0.5, nu=0.0, step=step)
    p_mut = evaluate_fd2(zero_velocity, slope_pressure, points, 0.5, nu=0.0, step=step)
    divergence_shift = vel_mut["divergence"] - base["divergence"]
    pressure_x_shift = p_mut["residual"][:, 0] - base["residual"][:, 0]
    return {
        "expected_shift": 0.02,
        "mean_divergence_shift": float(np.mean(divergence_shift)),
        "max_divergence_shift_error": float(np.max(np.abs(divergence_shift - 0.02))),
        "mean_pressure_x_momentum_shift": float(np.mean(pressure_x_shift)),
        "max_pressure_x_momentum_shift_error": float(np.max(np.abs(pressure_x_shift - 0.02))),
    }


def generate_report(
    *,
    output: str | Path,
    annulus_count: int = ANNULUS_COUNT,
    axis_count: int = AXIS_COUNT,
) -> dict[str, Any]:
    if annulus_count <= 0 or axis_count <= 0:
        raise ValueError("sample counts must be positive")

    leading, primary, _contract, sampled, compact, signed, adapter_metrics = _build_actual_signed_correction(
        radial_count=33,
        angular_count=8,
        phase_count=8,
    )

    def pressure(points: np.ndarray, time: float) -> np.ndarray:
        points = np.asarray(points, dtype=float)
        return np.asarray(leading.pressure(points[:, 0], points[:, 1], points[:, 2], time), dtype=float)

    def leading_velocity(points: np.ndarray, time: float) -> np.ndarray:
        return np.asarray(leading.at_points(points, time), dtype=float)

    def oscillatory_velocity(points: np.ndarray, time: float) -> np.ndarray:
        return leading_velocity(points, time) + np.asarray(primary.at_points(points, time), dtype=float)

    def corrected_velocity(points: np.ndarray, time: float) -> np.ndarray:
        return oscillatory_velocity(points, time) + np.asarray(signed.at_points(points, time), dtype=float)

    stages: dict[str, PointVelocity] = {
        "leading_only": leading_velocity,
        "leading_plus_oscillatory": oscillatory_velocity,
        "after_signed_correction": corrected_velocity,
    }
    strata = {
        "annulus": _sample_annulus(SEED, annulus_count),
        "axis_near": _sample_axis_near(SEED + 1, axis_count),
    }

    rows: list[dict[str, Any]] = []
    for step in STEPS:
        for stage_name, velocity in stages.items():
            for stratum_name, points in strata.items():
                volume = _stratum_volume(stratum_name)
                for time in TIMES:
                    result = evaluate_fd2(velocity, pressure, points, time, nu=NU, step=step)
                    rows.append({
                        "step": float(step),
                        "stage": stage_name,
                        "stratum": stratum_name,
                        "time": float(time),
                        "count": int(len(points)),
                        "metrics": _stage_metrics(result, points, volume),
                    })

    aggregates: dict[str, dict[str, dict[str, float]]] = {}
    for step in STEPS:
        aggregates[str(step)] = {stage: _aggregate(rows, stage, step) for stage in stages}

    finest = STEPS[-1]
    medium = STEPS[-2]
    finest_agg = aggregates[str(finest)]
    medium_agg = aggregates[str(medium)]
    leading_rms = finest_agg["leading_only"]["sample_rms"]
    osc_rms = finest_agg["leading_plus_oscillatory"]["sample_rms"]
    corrected_rms = finest_agg["after_signed_correction"]["sample_rms"]

    correction_points = _sample_annulus(CORRECTION_SEED, 64)
    correction_divergence = []
    for step in STEPS:
        div = divergence_fd2(signed.at_points, correction_points, 0.5, step)
        correction_divergence.append({
            "step": float(step),
            "rms": _scalar_rms(div),
            "max": float(np.max(np.abs(div))),
        })

    outside = np.asarray([
        [0.04, 0.0, 0.0],
        [0.0, -0.04, 0.05],
        [0.50, 0.0, 0.0],
        [-0.45, 0.10, -0.05],
    ])
    outside_values = signed.at_points(outside, 0.5)

    scale_sensitivity = []
    sensitivity_points = strata["annulus"]
    for scale in CORRECTION_SCALE_PROBES:
        def scaled(points: np.ndarray, time: float, scale: float = scale) -> np.ndarray:
            return oscillatory_velocity(points, time) + scale * np.asarray(signed.at_points(points, time), dtype=float)
        result = evaluate_fd2(scaled, pressure, sensitivity_points, 0.5, nu=NU, step=finest)
        scale_sensitivity.append({
            "scale": float(scale),
            "sample_rms": _vector_rms(result["residual"]),
            "sample_max": _vector_max(result["residual"]),
            "divergence_max": float(np.max(np.abs(result["divergence"]))),
            "used_for_selection": False,
        })

    mutation = _mutation_calibration(finest)
    medium_to_fine = {}
    for stage in stages:
        m = medium_agg[stage]
        f = finest_agg[stage]
        medium_to_fine[stage] = {
            "sample_rms_relative_change": abs(f["sample_rms"] - m["sample_rms"]) / max(abs(f["sample_rms"]), np.finfo(float).tiny),
            "sample_max_relative_change": abs(f["sample_max"] - m["sample_max"]) / max(abs(f["sample_max"]), np.finfo(float).tiny),
            "divergence_max_relative_change": abs(f["divergence_max"] - m["divergence_max"]) / max(abs(f["divergence_max"]), np.finfo(float).tiny),
        }

    report: dict[str, Any] = {
        "task": TASK,
        "validator": {
            "operator": "independent centered Cartesian FD2",
            "candidate_derivative_helpers_used": False,
            "agent3_cr006_residual_operator_used": False,
            "training_loss_or_tensor_used": False,
            "force": "zero raw diagnostic only; no fitted/free forcing",
            "nu": NU,
            "steps": list(STEPS),
            "times": list(TIMES),
            "seed": SEED,
            "normalization_scale": RESIDUAL_SCALE,
            "fixed_residual_reference": RESIDUAL_REFERENCE,
            "fixed_divergence_reference": DIVERGENCE_REFERENCE,
        },
        "candidate": {
            "leading_sha256": leading.sha256,
            "primary_amplitude": float(primary.amplitude),
            "primary_phase": float(primary.phase),
            "signed_profile_active_nodes": int(np.count_nonzero(sampled.active_target_mask)),
            "signed_profile_max_abs_delta_amplitude": float(np.max(np.abs(sampled.delta_amplitude))),
            "compact_fit_relative_rms": float(compact.sampled_fit_relative_rms),
            "compact_fit_design_condition": float(compact.design_condition),
            "agent3_adapter_metrics_bound_not_reused_for_residual": {
                "sampled_fit_relative_rms": float(adapter_metrics["sampled_fit_relative_rms"]),
            },
        },
        "sampling": {
            "annulus_count": annulus_count,
            "axis_count": axis_count,
            "annulus_r": list(ANNULUS_R),
            "annulus_z_half": ANNULUS_Z_HALF,
            "axis_r": list(AXIS_R),
            "axis_abs_z": list(AXIS_Z_ABS),
            "off_grid": True,
            "fresh_vs_agent3_seeds": True,
        },
        "rows": rows,
        "aggregates": aggregates,
        "finest_stage_ratios": {
            "oscillatory_over_leading_rms": osc_rms / max(leading_rms, np.finfo(float).tiny),
            "corrected_over_oscillatory_rms": corrected_rms / max(osc_rms, np.finfo(float).tiny),
            "corrected_over_leading_rms": corrected_rms / max(leading_rms, np.finfo(float).tiny),
        },
        "medium_to_fine_sensitivity": medium_to_fine,
        "signed_correction_only_divergence": correction_divergence,
        "signed_correction_support": {
            "outside_probe_max_velocity": _vector_max(outside_values),
            "outside_probe_exact_zero": bool(np.array_equal(outside_values, np.zeros_like(outside_values))),
            "annulus_nontrivial_rms": _vector_rms(signed.at_points(correction_points, 0.5)),
        },
        "correction_scale_sensitivity": scale_sensitivity,
        "mutation_calibration": mutation,
        "local_reference_checks": {
            "finest_corrected_normalized_sample_rms_le_1e_minus_3": bool(finest_agg["after_signed_correction"]["normalized_sample_rms"] <= RESIDUAL_REFERENCE),
            "finest_corrected_normalized_sample_max_le_1e_minus_3": bool(finest_agg["after_signed_correction"]["normalized_sample_max"] <= RESIDUAL_REFERENCE),
            "finest_corrected_divergence_rms_le_1e_minus_5": bool(finest_agg["after_signed_correction"]["divergence_rms"] <= DIVERGENCE_REFERENCE),
            "finest_corrected_divergence_max_le_1e_minus_5": bool(finest_agg["after_signed_correction"]["divergence_max"] <= DIVERGENCE_REFERENCE),
            "these_are_not_formal_full_domain_acceptance": True,
        },
        "truth_boundary": {
            "formal_full_domain_gate_assessed": False,
            "global_matched_leading_available": False,
            "compatible_final_fixed_or_restricted_forcing_available": False,
            "normalized_ns_residual_le_1e_minus_3_claimed": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
    }
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="Independently FD2-audit Agent-3 signed complete-curl correction")
    parser.add_argument("--output", default="artifacts/kokuno_agent4/independent_signed_curl_cycle_report.json")
    parser.add_argument("--annulus-count", type=int, default=ANNULUS_COUNT)
    parser.add_argument("--axis-count", type=int, default=AXIS_COUNT)
    args = parser.parse_args()
    report = generate_report(output=args.output, annulus_count=args.annulus_count, axis_count=args.axis_count)
    finest = report["aggregates"][str(STEPS[-1])]
    print(json.dumps({
        "task": report["task"],
        "finest": finest,
        "stage_ratios": report["finest_stage_ratios"],
        "signed_correction_divergence": report["signed_correction_only_divergence"][-1],
        "local_reference_checks": report["local_reference_checks"],
        "pde_validated": report["truth_boundary"]["pde_validated"],
    }, indent=2))


if __name__ == "__main__":
    main()
