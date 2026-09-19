"""Independent Agent-4 audit of the real oscillatory compact radial-stress seam.

This validator consumes only the admitted public oscillatory velocity and its
independently preflighted public time derivative.  It deliberately does *not*
call Agent-3's Cartesian FD4 self-defect operator, compact-stress constructor,
radial derivative, or training/construction tensors.

For the oscillatory component only it independently rebuilds

    R_osc = u_t + (u . grad) u - nu Delta u

with a Cartesian centered FD6 spatial operator, takes physical-angle means of
the theta/axial channels, reconstructs the compact moment-complement stress with
SciPy Simpson quadrature, and checks

    (d_r + e/r) sigma_e = -F + b_e M_e

using a separate centered FD6 radial derivative on three radial resolutions.
This is a correction-operator preflight only.  It is not a full-candidate NS
residual and cannot promote ``pde_validated``.
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

TASK = "KOKUNO-A4-INDEPENDENT-RADIAL-STRESS-AUDIT-043"
SCHEMA = "kokuno-a4-independent-radial-stress-audit-v1"
PARENT_AGENT3_PR = 589
PARENT_AGENT3_HEAD = "2271c52412922c8255cede8aee9a01e68ba55a55"
ADMITTED_AGENT2_PR = 561
ADMITTED_AGENT2_HEAD = "732800ce4990464b49c8aa32d0dff4580f6684d4"
INDEPENDENT_AGENT4_VELOCITY_PR = 563
INDEPENDENT_AGENT4_VELOCITY_HEAD = "9b0f86012c53fa8e32a19f766dbc150931870425"
INDEPENDENT_AGENT4_DT_PR = 580
INDEPENDENT_AGENT4_DT_HEAD = "6f7edb4d66dcb162a8510a40e4227f48f3c57e28"

RADIAL_COUNTS = (25, 49, 97)
ANGULAR_COUNT = 16
TIME = 0.50
Z = 0.08
CARTESIAN_FD6_STEP = 0.003
NU = 0.01

# Frozen local operator guards.  These are deliberately much looser than the
# final PDE gates because they measure only one compact radial inverse seam.
# They are fixed before exact-head CI evaluation and never replace the project
# momentum/divergence thresholds.
FINEST_RELATIVE_RMS_MAX = 2.0e-2
FINEST_RELATIVE_MAX_MAX = 5.0e-2
MIN_REFINEMENT_RATIO = 2.0
MOMENT_COMPLEMENT_RELATIVE_MAX = 1.0e-10
EDGE_RELATIVE_MAX = 1.0e-8
SIGN_FLIP_MUTATION_MIN = 5.0e-1

VelocityEvaluator = Callable[[Any, Any, Any, Any], np.ndarray]


def _vector_rms(values: np.ndarray) -> float:
    arr = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(np.sum(arr * arr, axis=-1))))


def _scalar_rms(values: np.ndarray) -> float:
    arr = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(arr * arr)))


def _fd6_axis_samples(
    evaluator: VelocityEvaluator,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
    *,
    axis: int,
    step: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Return centered FD6 first and second Cartesian derivatives of velocity."""
    if axis not in (0, 1, 2):
        raise ValueError("axis must be 0, 1 or 2")
    if not math.isfinite(step) or step <= 0.0:
        raise ValueError("step must be positive and finite")
    coordinates = [np.asarray(x, dtype=float), np.asarray(y, dtype=float), np.asarray(z, dtype=float)]
    samples: dict[int, np.ndarray] = {}
    for offset in (-3, -2, -1, 1, 2, 3):
        shifted = [v.copy() for v in coordinates]
        shifted[axis] = shifted[axis] + offset * step
        samples[offset] = np.asarray(evaluator(shifted[0], shifted[1], shifted[2], t), dtype=float)
    f0 = np.asarray(evaluator(x, y, z, t), dtype=float)
    first = (
        -samples[-3]
        + 9.0 * samples[-2]
        - 45.0 * samples[-1]
        + 45.0 * samples[1]
        - 9.0 * samples[2]
        + samples[3]
    ) / (60.0 * step)
    second = (
        2.0 * samples[-3]
        - 27.0 * samples[-2]
        + 270.0 * samples[-1]
        - 490.0 * f0
        + 270.0 * samples[1]
        - 27.0 * samples[2]
        + 2.0 * samples[3]
    ) / (180.0 * step * step)
    return first, second


def _independent_self_defect_fd6(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
    *,
    step: float = CARTESIAN_FD6_STEP,
    nu: float = NU,
) -> dict[str, np.ndarray]:
    """Rebuild the public oscillatory self-defect without Agent-3 operators."""
    u = np.asarray(velocity_osc(x, y, z, t), dtype=float)
    u_t = np.asarray(velocity_osc_dt(x, y, z, t), dtype=float)
    first_by_axis: list[np.ndarray] = []
    second_by_axis: list[np.ndarray] = []
    for axis in range(3):
        first, second = _fd6_axis_samples(
            velocity_osc, x, y, z, t, axis=axis, step=step
        )
        first_by_axis.append(first)
        second_by_axis.append(second)
    gradient = np.stack(first_by_axis, axis=-2)  # [..., derivative-axis, component]
    advection = np.einsum("...j,...ji->...i", u, gradient)
    laplacian = second_by_axis[0] + second_by_axis[1] + second_by_axis[2]
    raw = u_t + advection - float(nu) * laplacian
    return {
        "velocity": u,
        "velocity_dt": u_t,
        "advection": advection,
        "laplacian": laplacian,
        "raw_self_residual": raw,
    }


def _cartesian_to_cylindrical(vectors: np.ndarray, theta: np.ndarray) -> np.ndarray:
    values = np.asarray(vectors, dtype=float)
    angle = np.asarray(theta, dtype=float)
    c = np.cos(angle)
    s = np.sin(angle)
    radial = c * values[..., 0] + s * values[..., 1]
    azimuthal = -s * values[..., 0] + c * values[..., 1]
    axial = values[..., 2]
    return np.stack((radial, azimuthal, axial), axis=-1)


def _independent_ring_mean_sources(radial_count: int) -> dict[str, Any]:
    if int(radial_count) != radial_count or radial_count < 17:
        raise ValueError("radial_count must be an integer >=17")
    radial_count = int(radial_count)
    field = default_field()
    radial_min = float(field.radial_center - field.radial_halfwidth)
    radial_max = float(field.radial_center + field.radial_halfwidth)
    radii = np.linspace(radial_min, radial_max, radial_count)
    sample_radii = radii[1:-1]
    boundary_margin = min(
        float(sample_radii[0] - radial_min),
        float(radial_max - sample_radii[-1]),
    )
    if 3.0 * CARTESIAN_FD6_STEP >= boundary_margin:
        raise ValueError("frozen Cartesian FD6 stencil does not fit inside radial support")

    angles = 2.0 * math.pi * (np.arange(ANGULAR_COUNT, dtype=float) + 0.5) / ANGULAR_COUNT
    rr, tt = np.meshgrid(sample_radii, angles, indexing="ij")
    x = (rr * np.cos(tt)).reshape(-1)
    y = (rr * np.sin(tt)).reshape(-1)
    zz = np.full_like(x, Z)
    times = np.full_like(x, TIME)
    operator = _independent_self_defect_fd6(x, y, zz, times)
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


def _independent_compact_stress(
    radii: np.ndarray,
    source: np.ndarray,
    *,
    exponent: int,
    center: float,
    halfwidth: float,
) -> dict[str, np.ndarray | float]:
    """Build the compact moment-complement stress with Simpson quadrature."""
    r = np.asarray(radii, dtype=float)
    f = np.asarray(source, dtype=float)
    if r.ndim != 1 or f.shape != r.shape or r.size < 9:
        raise ValueError("radii/source must be matching one-dimensional arrays")
    s = (r - float(center)) / float(halfwidth)
    raw_bump = np.zeros_like(r)
    mask = np.abs(s) < 1.0
    raw_bump[mask] = np.cos(0.5 * math.pi * s[mask]) ** 8
    bump_norm = float(simpson((r**exponent) * raw_bump, x=r))
    if not math.isfinite(bump_norm) or bump_norm <= 0.0:
        raise RuntimeError("independent compact bump normalization failed")
    bump = raw_bump / bump_norm
    moment = float(simpson((r**exponent) * f, x=r))
    complement = f - bump * moment
    weighted = (r**exponent) * complement
    cumulative = np.asarray(cumulative_simpson(weighted, x=r, initial=0.0), dtype=float)
    stress = -cumulative / (r**exponent)
    return {
        "stress": stress,
        "bump": bump,
        "moment": moment,
        "complement": complement,
    }


def _fd6_first_uniform(values: np.ndarray, coordinates: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    y = np.asarray(values, dtype=float)
    x = np.asarray(coordinates, dtype=float)
    if y.ndim != 1 or x.ndim != 1 or y.shape != x.shape or y.size < 9:
        raise ValueError("FD6 derivative requires matching arrays with at least nine nodes")
    steps = np.diff(x)
    h = float(steps[0])
    if h <= 0.0 or not np.allclose(steps, h, rtol=0.0, atol=128.0 * math.ulp(max(1.0, abs(h)))):
        raise ValueError("radial grid must be uniformly spaced")
    derivative = (
        -y[:-6]
        + 9.0 * y[1:-5]
        - 45.0 * y[2:-4]
        + 45.0 * y[4:-2]
        - 9.0 * y[5:-1]
        + y[6:]
    ) / (60.0 * h)
    return x[3:-3], derivative


def _audit_channel(
    radii: np.ndarray,
    source: np.ndarray,
    *,
    exponent: int,
    center: float,
    halfwidth: float,
) -> dict[str, float | int]:
    built = _independent_compact_stress(
        radii,
        source,
        exponent=exponent,
        center=center,
        halfwidth=halfwidth,
    )
    stress = np.asarray(built["stress"], dtype=float)
    bump = np.asarray(built["bump"], dtype=float)
    moment = float(built["moment"])
    complement = np.asarray(built["complement"], dtype=float)

    r_mid, d_stress = _fd6_first_uniform(stress, radii)
    sigma_mid = stress[3:-3]
    lhs = d_stress + exponent * sigma_mid / r_mid
    rhs_full = -source + bump * moment
    rhs = rhs_full[3:-3]
    error = lhs - rhs
    rhs_rms = _scalar_rms(rhs)
    rhs_max = float(np.max(np.abs(rhs)))
    rel_rms = _scalar_rms(error) / max(rhs_rms, 1.0e-300)
    rel_max = float(np.max(np.abs(error))) / max(rhs_max, 1.0e-300)

    independent_moment_complement = float(simpson((radii**exponent) * complement, x=radii))
    source_weight_norm = float(simpson(np.abs((radii**exponent) * source), x=radii))
    moment_relative = abs(independent_moment_complement) / max(source_weight_norm, 1.0e-300)
    stress_scale = max(_scalar_rms(stress), 1.0)
    edge_relative = max(abs(float(stress[0])), abs(float(stress[-1]))) / stress_scale

    mutated_lhs = -d_stress - exponent * sigma_mid / r_mid
    mutation_error = mutated_lhs - rhs
    mutation_relative_rms = _scalar_rms(mutation_error) / max(rhs_rms, 1.0e-300)

    return {
        "exponent": exponent,
        "weighted_moment": moment,
        "moment_complement_relative": moment_relative,
        "stress_rms": _scalar_rms(stress),
        "stress_max_abs": float(np.max(np.abs(stress))),
        "inner_edge": float(stress[0]),
        "outer_edge": float(stress[-1]),
        "edge_relative": edge_relative,
        "operator_relative_rms": rel_rms,
        "operator_relative_max": rel_max,
        "operator_error_rms": _scalar_rms(error),
        "operator_error_max_abs": float(np.max(np.abs(error))),
        "rhs_rms": rhs_rms,
        "rhs_max_abs": rhs_max,
        "sign_flip_mutation_relative_rms": mutation_relative_rms,
        "fd6_interior_count": int(r_mid.size),
    }


def audit_independent_radial_stress_closure() -> dict[str, Any]:
    field = default_field()
    levels: list[dict[str, Any]] = []
    for radial_count in RADIAL_COUNTS:
        actual = _independent_ring_mean_sources(radial_count)
        radii = np.asarray(actual["radii"], dtype=float)
        theta = _audit_channel(
            radii,
            np.asarray(actual["theta_source"], dtype=float),
            exponent=2,
            center=float(field.radial_center),
            halfwidth=float(field.radial_halfwidth),
        )
        axial = _audit_channel(
            radii,
            np.asarray(actual["axial_source"], dtype=float),
            exponent=1,
            center=float(field.radial_center),
            halfwidth=float(field.radial_halfwidth),
        )
        levels.append(
            {
                "radial_count": radial_count,
                "radial_spacing": float(radii[1] - radii[0]),
                "cartesian_fd6_step": CARTESIAN_FD6_STEP,
                "boundary_margin": float(actual["boundary_margin"]),
                "velocity_rms": float(actual["velocity_rms"]),
                "self_defect_rms": float(actual["self_defect_rms"]),
                "theta_e2": theta,
                "axial_e1": axial,
            }
        )

    convergence: dict[str, Any] = {}
    failed_guards: list[str] = []
    for channel in ("theta_e2", "axial_e1"):
        errors = [float(level[channel]["operator_relative_rms"]) for level in levels]
        ratios = [errors[i] / max(errors[i + 1], 1.0e-300) for i in range(len(errors) - 1)]
        finest = levels[-1][channel]
        convergence[channel] = {
            "relative_rms_by_radial_count": errors,
            "refinement_ratios": ratios,
            "finest_relative_rms": errors[-1],
            "finest_relative_max": float(finest["operator_relative_max"]),
            "finest_moment_complement_relative": float(finest["moment_complement_relative"]),
            "finest_edge_relative": float(finest["edge_relative"]),
            "finest_sign_flip_mutation_relative_rms": float(finest["sign_flip_mutation_relative_rms"]),
        }
        if errors[-1] > FINEST_RELATIVE_RMS_MAX:
            failed_guards.append(f"{channel}:finest_relative_rms")
        if float(finest["operator_relative_max"]) > FINEST_RELATIVE_MAX_MAX:
            failed_guards.append(f"{channel}:finest_relative_max")
        if min(ratios) < MIN_REFINEMENT_RATIO:
            failed_guards.append(f"{channel}:minimum_refinement_ratio")
        if float(finest["moment_complement_relative"]) > MOMENT_COMPLEMENT_RELATIVE_MAX:
            failed_guards.append(f"{channel}:moment_complement_relative")
        if float(finest["edge_relative"]) > EDGE_RELATIVE_MAX:
            failed_guards.append(f"{channel}:edge_relative")
        if float(finest["sign_flip_mutation_relative_rms"]) < SIGN_FLIP_MUTATION_MIN:
            failed_guards.append(f"{channel}:sign_flip_mutation")

    passed = not failed_guards
    return {
        "task": TASK,
        "schema": SCHEMA,
        "provenance": {
            "parent_agent3_pr": PARENT_AGENT3_PR,
            "parent_agent3_head": PARENT_AGENT3_HEAD,
            "admitted_agent2_pr": ADMITTED_AGENT2_PR,
            "admitted_agent2_head": ADMITTED_AGENT2_HEAD,
            "independent_agent4_velocity_pr": INDEPENDENT_AGENT4_VELOCITY_PR,
            "independent_agent4_velocity_head": INDEPENDENT_AGENT4_VELOCITY_HEAD,
            "independent_agent4_dt_pr": INDEPENDENT_AGENT4_DT_PR,
            "independent_agent4_dt_head": INDEPENDENT_AGENT4_DT_HEAD,
            "velocity_provider": "openai_ns_reconstruction.kokuno_public_z_pullback_velocity:velocity_osc",
            "time_derivative_provider": "openai_ns_reconstruction.kokuno_public_oscillatory_time_derivative:velocity_osc_dt",
            "spatial_operator": "new centered Cartesian FD6; no Agent-3 FD4 helper",
            "stress_constructor": "new Simpson/cumulative-Simpson moment-complement implementation",
            "radial_derivative": "new centered uniform-grid FD6",
        },
        "frozen_state": {
            "radial_counts": list(RADIAL_COUNTS),
            "angular_count": ANGULAR_COUNT,
            "time": TIME,
            "z": Z,
            "cartesian_fd6_step": CARTESIAN_FD6_STEP,
            "nu": NU,
        },
        "frozen_guards": {
            "finest_relative_rms_max": FINEST_RELATIVE_RMS_MAX,
            "finest_relative_max_max": FINEST_RELATIVE_MAX_MAX,
            "minimum_refinement_ratio": MIN_REFINEMENT_RATIO,
            "moment_complement_relative_max": MOMENT_COMPLEMENT_RELATIVE_MAX,
            "edge_relative_max": EDGE_RELATIVE_MAX,
            "sign_flip_mutation_relative_rms_min": SIGN_FLIP_MUTATION_MIN,
            "project_momentum_max_l2_gate": 1.0e-3,
            "project_divergence_max_l2_gate": 1.0e-5,
        },
        "levels": levels,
        "convergence": convergence,
        "failed_guards": failed_guards,
        "independent_radial_stress_preflight_passed": passed,
        "truth_boundary": {
            "public_black_box_velocity_consumed": True,
            "public_black_box_time_derivative_consumed": True,
            "agent3_fd4_self_defect_helper_used": False,
            "agent3_compact_stress_constructor_used": False,
            "agent3_radial_derivative_used": False,
            "independent_real_oscillatory_self_defect_consumed": True,
            "independent_compact_radial_stress_reconstruction_audited": True,
            "full_same_cycle_composite_requested_stress_materialized": False,
            "agent1_leading_cross_terms_included": False,
            "matched_pressure_included": False,
            "restricted_forcing_included": False,
            "candidate_finite_head_mean_debt_materialized": False,
            "signed_mean_inverse_input_ready": False,
            "public_velocity_correction_materialized": False,
            "finite_correction_cycle_run": False,
            "heldout_ns_momentum_residual_assessed": False,
            "formal_full_domain_pde_gate_assessed": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
    }


def write_receipt(output: str | Path) -> dict[str, Any]:
    receipt = audit_independent_radial_stress_closure()
    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        default="artifacts/kokuno_agent4_independent_radial_stress_audit_v1.json",
    )
    args = parser.parse_args()
    print(json.dumps(write_receipt(args.output), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
