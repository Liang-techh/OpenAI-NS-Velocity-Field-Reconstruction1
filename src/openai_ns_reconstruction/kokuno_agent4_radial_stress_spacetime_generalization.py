"""Independent Agent-4 spacetime generalization audit for the Kokuno radial-stress seam.

This is deliberately narrower than the final Navier--Stokes gate.  It consumes
only the already-admitted public oscillatory velocity and public time derivative,
rebuilds the oscillatory self-defect with a fresh Cartesian FD8 operator, and
checks the compact radial moment-complement stress equation at three off-grid
spacetime states and three radial resolutions.

It does not call Agent-3's defect operator, stress constructor, radial derivative,
training tensors, or any pressure/forcing fit.  A scientific FAIL is retained in
the emitted receipt rather than changing thresholds after evaluation.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Callable

import numpy as np
from scipy.integrate import cumulative_simpson, simpson

from .kokuno_public_oscillatory_time_derivative import velocity_osc_dt
from .kokuno_public_z_pullback_velocity import default_field, velocity_osc

TASK = "KOKUNO-A4-RADIAL-STRESS-SPACETIME-GENERALIZATION-044"
SCHEMA = "kokuno-a4-radial-stress-spacetime-generalization-v1"
PARENT_AGENT3_PR = 598
PARENT_AGENT3_HEAD = "a97f0882ab0b1fcc0d4b9575b5366fb44b35c233"
ADMITTED_AGENT2_PR = 561
ADMITTED_AGENT2_HEAD = "732800ce4990464b49c8aa32d0dff4580f6684d4"
PRIOR_AGENT4_PR = 595
PRIOR_AGENT4_HEAD = "19296acad4f84ba05d3b095bb8c01bf5a2c95891"

RADIAL_COUNTS = (33, 65, 129)
ANGULAR_COUNT = 24
SPACETIME_STATES = (
    {"t": 0.37, "z": -0.31},
    {"t": 0.50, "z": 0.08},
    {"t": 0.63, "z": 0.31},
)
CARTESIAN_FD8_STEP = 0.0025
NU = 0.01

# Frozen before exact-head CI evaluation.  These are local correction-operator
# guards, not replacements for the project 1e-3 momentum / 1e-5 divergence gates.
FINEST_RELATIVE_RMS_MAX = 2.0e-2
FINEST_RELATIVE_MAX_MAX = 5.0e-2
MIN_REFINEMENT_RATIO = 2.0
MOMENT_COMPLEMENT_RELATIVE_MAX = 1.0e-10
EDGE_RELATIVE_MAX = 1.0e-8
SIGN_FLIP_MUTATION_MIN = 5.0e-1
NONTRIVIAL_SELF_DEFECT_RMS_MIN = 1.0e-8

VelocityEvaluator = Callable[[Any, Any, Any, Any], np.ndarray]


def _scalar_rms(values: np.ndarray) -> float:
    arr = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(arr * arr)))


def _vector_rms(values: np.ndarray) -> float:
    arr = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(np.sum(arr * arr, axis=-1))))


def _fd8_axis_samples(
    evaluator: VelocityEvaluator,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
    *,
    axis: int,
    step: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Centered eighth-order first/second Cartesian derivatives."""
    if axis not in (0, 1, 2):
        raise ValueError("axis must be 0, 1, or 2")
    if not math.isfinite(step) or step <= 0.0:
        raise ValueError("step must be positive and finite")
    coordinates = [np.asarray(x, dtype=float), np.asarray(y, dtype=float), np.asarray(z, dtype=float)]
    samples: dict[int, np.ndarray] = {}
    for offset in (-4, -3, -2, -1, 1, 2, 3, 4):
        shifted = [v.copy() for v in coordinates]
        shifted[axis] = shifted[axis] + offset * step
        samples[offset] = np.asarray(evaluator(shifted[0], shifted[1], shifted[2], t), dtype=float)
    f0 = np.asarray(evaluator(x, y, z, t), dtype=float)
    first = (
        samples[-4] / 280.0
        - 4.0 * samples[-3] / 105.0
        + samples[-2] / 5.0
        - 4.0 * samples[-1] / 5.0
        + 4.0 * samples[1] / 5.0
        - samples[2] / 5.0
        + 4.0 * samples[3] / 105.0
        - samples[4] / 280.0
    ) / step
    second = (
        -samples[-4] / 560.0
        + 8.0 * samples[-3] / 315.0
        - samples[-2] / 5.0
        + 8.0 * samples[-1] / 5.0
        - 205.0 * f0 / 72.0
        + 8.0 * samples[1] / 5.0
        - samples[2] / 5.0
        + 8.0 * samples[3] / 315.0
        - samples[4] / 560.0
    ) / (step * step)
    return first, second


def _self_defect_fd8(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
    *,
    step: float = CARTESIAN_FD8_STEP,
    nu: float = NU,
) -> dict[str, np.ndarray]:
    u = np.asarray(velocity_osc(x, y, z, t), dtype=float)
    u_t = np.asarray(velocity_osc_dt(x, y, z, t), dtype=float)
    first_by_axis: list[np.ndarray] = []
    second_by_axis: list[np.ndarray] = []
    for axis in range(3):
        first, second = _fd8_axis_samples(velocity_osc, x, y, z, t, axis=axis, step=step)
        first_by_axis.append(first)
        second_by_axis.append(second)
    gradient = np.stack(first_by_axis, axis=-2)
    advection = np.einsum("...j,...ji->...i", u, gradient)
    laplacian = second_by_axis[0] + second_by_axis[1] + second_by_axis[2]
    residual = u_t + advection - float(nu) * laplacian
    return {"velocity": u, "velocity_dt": u_t, "raw_self_residual": residual}


def _cartesian_to_cylindrical(vectors: np.ndarray, theta: np.ndarray) -> np.ndarray:
    values = np.asarray(vectors, dtype=float)
    angle = np.asarray(theta, dtype=float)
    c = np.cos(angle)
    s = np.sin(angle)
    radial = c * values[..., 0] + s * values[..., 1]
    azimuthal = -s * values[..., 0] + c * values[..., 1]
    axial = values[..., 2]
    return np.stack((radial, azimuthal, axial), axis=-1)


def _ring_mean_sources(radial_count: int, *, time: float, z_value: float) -> dict[str, Any]:
    field = default_field()
    radial_min = float(field.radial_center - field.radial_halfwidth)
    radial_max = float(field.radial_center + field.radial_halfwidth)
    radii = np.linspace(radial_min, radial_max, int(radial_count))
    sample_radii = radii[1:-1]
    boundary_margin = min(float(sample_radii[0] - radial_min), float(radial_max - sample_radii[-1]))
    if 4.0 * CARTESIAN_FD8_STEP >= boundary_margin:
        raise ValueError("frozen Cartesian FD8 stencil does not fit inside radial support")

    angles = 2.0 * math.pi * (np.arange(ANGULAR_COUNT, dtype=float) + 0.5) / ANGULAR_COUNT
    rr, tt = np.meshgrid(sample_radii, angles, indexing="ij")
    x = (rr * np.cos(tt)).reshape(-1)
    y = (rr * np.sin(tt)).reshape(-1)
    zz = np.full_like(x, float(z_value))
    times = np.full_like(x, float(time))
    operator = _self_defect_fd8(x, y, zz, times)
    cylindrical = _cartesian_to_cylindrical(operator["raw_self_residual"], tt.reshape(-1))
    interior = cylindrical.reshape(radial_count - 2, ANGULAR_COUNT, 3)
    ring_mean = np.mean(interior, axis=1)
    full = np.zeros((radial_count, 3), dtype=float)
    full[1:-1] = ring_mean
    return {
        "radii": radii,
        "theta_source": full[:, 1],
        "axial_source": full[:, 2],
        "boundary_margin": boundary_margin,
        "velocity_rms": _vector_rms(operator["velocity"]),
        "self_defect_rms": _vector_rms(operator["raw_self_residual"]),
    }


def _compact_stress(
    radii: np.ndarray,
    source: np.ndarray,
    *,
    exponent: int,
    center: float,
    halfwidth: float,
) -> dict[str, np.ndarray | float]:
    r = np.asarray(radii, dtype=float)
    f = np.asarray(source, dtype=float)
    s = (r - float(center)) / float(halfwidth)
    raw_bump = np.zeros_like(r)
    mask = np.abs(s) < 1.0
    raw_bump[mask] = np.cos(0.5 * math.pi * s[mask]) ** 8
    bump_norm = float(simpson((r**exponent) * raw_bump, x=r))
    if not math.isfinite(bump_norm) or bump_norm <= 0.0:
        raise RuntimeError("compact bump normalization failed")
    bump = raw_bump / bump_norm
    moment = float(simpson((r**exponent) * f, x=r))
    complement = f - bump * moment
    weighted = (r**exponent) * complement
    cumulative = np.asarray(cumulative_simpson(weighted, x=r, initial=0.0), dtype=float)
    stress = -cumulative / (r**exponent)
    return {"stress": stress, "bump": bump, "moment": moment, "complement": complement}


def _fd8_first_uniform(values: np.ndarray, coordinates: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    y = np.asarray(values, dtype=float)
    x = np.asarray(coordinates, dtype=float)
    if y.ndim != 1 or x.ndim != 1 or y.shape != x.shape or y.size < 11:
        raise ValueError("FD8 derivative requires matching arrays with at least eleven nodes")
    steps = np.diff(x)
    h = float(steps[0])
    if h <= 0.0 or not np.allclose(steps, h, rtol=0.0, atol=256.0 * math.ulp(max(1.0, abs(h)))):
        raise ValueError("radial grid must be uniformly spaced")
    derivative = (
        y[:-8] / 280.0
        - 4.0 * y[1:-7] / 105.0
        + y[2:-6] / 5.0
        - 4.0 * y[3:-5] / 5.0
        + 4.0 * y[5:-3] / 5.0
        - y[6:-2] / 5.0
        + 4.0 * y[7:-1] / 105.0
        - y[8:] / 280.0
    ) / h
    return x[4:-4], derivative


def _audit_channel(
    radii: np.ndarray,
    source: np.ndarray,
    *,
    exponent: int,
    center: float,
    halfwidth: float,
) -> dict[str, float | int]:
    built = _compact_stress(radii, source, exponent=exponent, center=center, halfwidth=halfwidth)
    stress = np.asarray(built["stress"], dtype=float)
    bump = np.asarray(built["bump"], dtype=float)
    moment = float(built["moment"])
    complement = np.asarray(built["complement"], dtype=float)

    r_mid, d_stress = _fd8_first_uniform(stress, radii)
    sigma_mid = stress[4:-4]
    lhs = d_stress + exponent * sigma_mid / r_mid
    rhs_full = -source + bump * moment
    rhs = rhs_full[4:-4]
    error = lhs - rhs
    rhs_rms = _scalar_rms(rhs)
    rhs_max = float(np.max(np.abs(rhs)))
    relative_rms = _scalar_rms(error) / max(rhs_rms, 1.0e-300)
    relative_max = float(np.max(np.abs(error))) / max(rhs_max, 1.0e-300)

    complement_moment = float(simpson((radii**exponent) * complement, x=radii))
    source_weight_norm = float(simpson(np.abs((radii**exponent) * source), x=radii))
    moment_relative = abs(complement_moment) / max(source_weight_norm, 1.0e-300)
    edge_relative = max(abs(float(stress[0])), abs(float(stress[-1]))) / max(_scalar_rms(stress), 1.0)

    mutated_lhs = -d_stress - exponent * sigma_mid / r_mid
    mutation_relative_rms = _scalar_rms(mutated_lhs - rhs) / max(rhs_rms, 1.0e-300)
    return {
        "exponent": exponent,
        "weighted_moment": moment,
        "operator_relative_rms": relative_rms,
        "operator_relative_max": relative_max,
        "moment_complement_relative": moment_relative,
        "edge_relative": edge_relative,
        "sign_flip_mutation_relative_rms": mutation_relative_rms,
        "rhs_rms": rhs_rms,
        "stress_rms": _scalar_rms(stress),
        "fd8_interior_count": int(r_mid.size),
    }


def _channel_convergence(levels: list[dict[str, Any]], key: str) -> dict[str, Any]:
    rms = [float(level[key]["operator_relative_rms"]) for level in levels]
    ratios = [rms[i] / max(rms[i + 1], 1.0e-300) for i in range(len(rms) - 1)]
    fine = levels[-1][key]
    return {
        "relative_rms_by_radial_count": rms,
        "refinement_ratios": ratios,
        "finest_relative_max": float(fine["operator_relative_max"]),
        "finest_moment_complement_relative": float(fine["moment_complement_relative"]),
        "finest_edge_relative": float(fine["edge_relative"]),
        "finest_sign_flip_mutation_relative_rms": float(fine["sign_flip_mutation_relative_rms"]),
        "finest_weighted_moment": float(fine["weighted_moment"]),
    }


def audit_spacetime_generalization() -> dict[str, Any]:
    field = default_field()
    state_reports: list[dict[str, Any]] = []
    failed_guards: list[str] = []

    for state_index, state in enumerate(SPACETIME_STATES):
        levels: list[dict[str, Any]] = []
        for radial_count in RADIAL_COUNTS:
            actual = _ring_mean_sources(radial_count, time=float(state["t"]), z_value=float(state["z"]))
            radii = np.asarray(actual["radii"], dtype=float)
            levels.append(
                {
                    "radial_count": radial_count,
                    "radial_spacing": float(radii[1] - radii[0]),
                    "velocity_rms": float(actual["velocity_rms"]),
                    "self_defect_rms": float(actual["self_defect_rms"]),
                    "theta_e2": _audit_channel(
                        radii,
                        np.asarray(actual["theta_source"], dtype=float),
                        exponent=2,
                        center=float(field.radial_center),
                        halfwidth=float(field.radial_halfwidth),
                    ),
                    "axial_e1": _audit_channel(
                        radii,
                        np.asarray(actual["axial_source"], dtype=float),
                        exponent=1,
                        center=float(field.radial_center),
                        halfwidth=float(field.radial_halfwidth),
                    ),
                }
            )

        convergence = {
            "theta_e2": _channel_convergence(levels, "theta_e2"),
            "axial_e1": _channel_convergence(levels, "axial_e1"),
        }
        state_failed: list[str] = []
        if float(levels[-1]["self_defect_rms"]) < NONTRIVIAL_SELF_DEFECT_RMS_MIN:
            state_failed.append("nontrivial_self_defect_rms")
        for channel_name, c in convergence.items():
            if c["relative_rms_by_radial_count"][-1] > FINEST_RELATIVE_RMS_MAX:
                state_failed.append(f"{channel_name}:finest_relative_rms")
            if c["finest_relative_max"] > FINEST_RELATIVE_MAX_MAX:
                state_failed.append(f"{channel_name}:finest_relative_max")
            if min(c["refinement_ratios"]) < MIN_REFINEMENT_RATIO:
                state_failed.append(f"{channel_name}:minimum_refinement_ratio")
            if c["finest_moment_complement_relative"] > MOMENT_COMPLEMENT_RELATIVE_MAX:
                state_failed.append(f"{channel_name}:moment_complement")
            if c["finest_edge_relative"] > EDGE_RELATIVE_MAX:
                state_failed.append(f"{channel_name}:edge_relative")
            if c["finest_sign_flip_mutation_relative_rms"] < SIGN_FLIP_MUTATION_MIN:
                state_failed.append(f"{channel_name}:sign_flip_mutation")
        failed_guards.extend([f"state{state_index}:{name}" for name in state_failed])
        state_reports.append(
            {
                "state_index": state_index,
                "t": float(state["t"]),
                "z": float(state["z"]),
                "levels": levels,
                "convergence": convergence,
                "failed_guards": state_failed,
                "state_passed": not state_failed,
            }
        )

    return {
        "schema": SCHEMA,
        "task": TASK,
        "provenance": {
            "parent_agent3_pr": PARENT_AGENT3_PR,
            "parent_agent3_head": PARENT_AGENT3_HEAD,
            "admitted_agent2_pr": ADMITTED_AGENT2_PR,
            "admitted_agent2_head": ADMITTED_AGENT2_HEAD,
            "prior_agent4_pr": PRIOR_AGENT4_PR,
            "prior_agent4_head": PRIOR_AGENT4_HEAD,
        },
        "protocol": {
            "radial_counts": list(RADIAL_COUNTS),
            "angular_count": ANGULAR_COUNT,
            "spacetime_states": [dict(s) for s in SPACETIME_STATES],
            "cartesian_derivative": "centered FD8, independently implemented",
            "cartesian_fd8_step": CARTESIAN_FD8_STEP,
            "radial_derivative": "separate centered FD8",
            "nu": NU,
            "guards": {
                "finest_relative_rms_max": FINEST_RELATIVE_RMS_MAX,
                "finest_relative_max_max": FINEST_RELATIVE_MAX_MAX,
                "minimum_refinement_ratio": MIN_REFINEMENT_RATIO,
                "moment_complement_relative_max": MOMENT_COMPLEMENT_RELATIVE_MAX,
                "edge_relative_max": EDGE_RELATIVE_MAX,
                "sign_flip_mutation_min": SIGN_FLIP_MUTATION_MIN,
                "nontrivial_self_defect_rms_min": NONTRIVIAL_SELF_DEFECT_RMS_MIN,
            },
        },
        "states": state_reports,
        "failed_guards": failed_guards,
        "spacetime_radial_stress_generalization_passed": not failed_guards,
        "truth_boundary": {
            "public_black_box_velocity_consumed": True,
            "public_black_box_time_derivative_consumed": True,
            "agent3_self_defect_operator_used": False,
            "agent3_compact_stress_constructor_used": False,
            "agent3_radial_derivative_used": False,
            "pressure_fitted": False,
            "forcing_fitted": False,
            "full_same_cycle_composite_requested_stress_materialized": False,
            "public_velocity_correction_materialized": False,
            "heldout_ns_momentum_residual_assessed": False,
            "formal_full_domain_pde_gate_assessed": False,
            "pde_validated": False,
            "final_normalized_momentum_gate": 1.0e-3,
            "final_normalized_divergence_gate": 1.0e-5,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = audit_spacetime_generalization()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
