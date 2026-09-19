"""Independent closure audit for Agent-3 compact radial stress reconstruction.

This module does not create a correction velocity.  It reuses the already
admitted Agent-2 public oscillatory field and Agent-3's actual self-defect
operator, then independently checks the compact radial inverse used by the
mean-correction lane:

    M_e F = integral r^e F dr,
    P_e F = F - b_e M_e F,
    sigma_e(r) = -r^(-e) integral r^e P_e dr,

which should satisfy

    (d_r + e/r) sigma_e = -F + b_e M_e.

The audit uses a direct fourth-order radial derivative that is separate from
the cumulative trapezoid used to construct ``sigma_e``.  It runs on a frozen
radial-resolution ladder and records moment closure, compact endpoint closure,
and direct operator residuals for the real oscillatory theta/e=2 and axial/e=1
mean-defect channels.

This remains component-level evidence only.  Agent-1 leading/cross terms,
matched pressure, restricted forcing, signed covariance inversion and finite
correction-cycle execution remain downstream and fail closed here.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_actual_oscillatory_mean_stress import (
    ADMITTED_AGENT2_HEAD,
    ADMITTED_AGENT2_PR,
    INDEPENDENT_AGENT4_HEAD,
    INDEPENDENT_AGENT4_PR,
    _cartesian_to_cylindrical,
    _compact_radial_stress,
    _fd4_spatial_operator,
)
from .kokuno_public_z_pullback_velocity import default_field

TASK = "KOKUNO-A3-RADIAL-STRESS-CLOSURE-041"
SCHEMA = "kokuno-a3-radial-stress-closure-v1"
PARENT_AGENT3_MEAN_DEBT_PR = 586
PARENT_AGENT3_MEAN_DEBT_HEAD = "ca4da20c90ccd7a5046f79fb2e7965ec05b6086f"

RADIAL_COUNTS = (25, 49, 97)
ANGULAR_COUNT = 32
TIME = 0.50
Z = 0.08
FD4_STEP = 0.005
NU = 0.01

# Preregistered diagnostic guards.  They are intentionally local radial-inverse
# checks and do not replace the final project momentum/divergence thresholds.
FINEST_RELATIVE_RMS_MAX = 5.0e-2
FINEST_RELATIVE_MAX_MAX = 1.0e-1
MIN_REFINEMENT_RATIO = 1.5
MOMENT_RELATIVE_MAX = 1.0e-10
EDGE_RELATIVE_MAX = 1.0e-10


def _trapz(values: np.ndarray, coordinates: np.ndarray) -> float:
    y = np.asarray(values, dtype=float)
    x = np.asarray(coordinates, dtype=float)
    if y.ndim != 1 or x.ndim != 1 or y.shape != x.shape or y.size < 2:
        raise ValueError("trapz requires matching one-dimensional arrays")
    return float(np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(x)))


def _scalar_rms(values: np.ndarray) -> float:
    v = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(v * v)))


def _independent_bump(
    radii: np.ndarray,
    *,
    exponent: int,
    center: float,
    halfwidth: float,
) -> np.ndarray:
    """Rebuild the repository numerical cos^8 bump without calling A3's helper."""
    r = np.asarray(radii, dtype=float)
    s = (r - float(center)) / float(halfwidth)
    raw = np.zeros_like(r)
    mask = np.abs(s) < 1.0
    raw[mask] = np.cos(0.5 * math.pi * s[mask]) ** 8
    weighted = (r**exponent) * raw
    norm = _trapz(weighted, r)
    if not math.isfinite(norm) or norm <= 0.0:
        raise RuntimeError("independent compact bump normalization failed")
    return raw / norm


def _fd4_first_uniform(values: np.ndarray, coordinates: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    y = np.asarray(values, dtype=float)
    x = np.asarray(coordinates, dtype=float)
    if y.ndim != 1 or x.ndim != 1 or y.shape != x.shape or y.size < 9:
        raise ValueError("FD4 derivative requires matching arrays with at least nine nodes")
    steps = np.diff(x)
    h = float(steps[0])
    if h <= 0.0 or not np.allclose(steps, h, rtol=0.0, atol=64.0 * math.ulp(max(1.0, abs(h)))):
        raise ValueError("radial grid must be uniformly spaced")
    derivative = (
        y[:-4]
        - 8.0 * y[1:-3]
        + 8.0 * y[3:-1]
        - y[4:]
    ) / (12.0 * h)
    return x[2:-2], derivative


def _actual_ring_mean_sources(radial_count: int) -> dict[str, Any]:
    """Recompute the real oscillatory ring-mean theta/axial defect profiles."""
    if int(radial_count) != radial_count or radial_count < 17:
        raise ValueError("radial_count must be an integer >=17")
    radial_count = int(radial_count)
    field = default_field()
    radial_min = field.radial_center - field.radial_halfwidth
    radial_max = field.radial_center + field.radial_halfwidth
    radii = np.linspace(radial_min, radial_max, radial_count)
    sample_radii = radii[1:-1]
    boundary_margin = min(
        float(sample_radii[0] - radial_min),
        float(radial_max - sample_radii[-1]),
    )
    if 2.0 * FD4_STEP >= boundary_margin:
        raise ValueError(
            "frozen Cartesian FD4 stencil no longer fits inside radial support at this resolution"
        )

    angles = 2.0 * math.pi * (np.arange(ANGULAR_COUNT, dtype=float) + 0.5) / ANGULAR_COUNT
    rr, tt = np.meshgrid(sample_radii, angles, indexing="ij")
    x = (rr * np.cos(tt)).reshape(-1)
    y = (rr * np.sin(tt)).reshape(-1)
    zz = np.full_like(x, Z)
    times = np.full_like(x, TIME)
    operator = _fd4_spatial_operator(x, y, zz, times, step=FD4_STEP, nu=NU)
    cyl = _cartesian_to_cylindrical(operator["raw_self_residual"], tt.reshape(-1))
    interior = cyl.reshape(radial_count - 2, ANGULAR_COUNT, 3)
    ring_mean = np.mean(interior, axis=1)
    full = np.zeros((radial_count, 3), dtype=float)
    full[1:-1] = ring_mean
    return {
        "radii": radii,
        "theta_source": full[:, 1],
        "axial_source": full[:, 2],
        "field": field,
        "boundary_margin": boundary_margin,
    }


def _audit_channel(
    radii: np.ndarray,
    source: np.ndarray,
    *,
    exponent: int,
    center: float,
    halfwidth: float,
) -> dict[str, Any]:
    built = _compact_radial_stress(
        radii,
        source,
        exponent=exponent,
        bump_center=center,
        bump_halfwidth=halfwidth,
    )
    stress = np.asarray(built["stress"], dtype=float)
    bump = _independent_bump(
        radii,
        exponent=exponent,
        center=center,
        halfwidth=halfwidth,
    )
    moment = _trapz((radii**exponent) * source, radii)
    complement = source - bump * moment
    independent_moment = _trapz((radii**exponent) * complement, radii)

    r_mid, d_stress = _fd4_first_uniform(stress, radii)
    sigma_mid = stress[2:-2]
    lhs = d_stress + exponent * sigma_mid / r_mid
    rhs_full = -source + bump * moment
    rhs = rhs_full[2:-2]
    error = lhs - rhs
    rhs_rms = _scalar_rms(rhs)
    rhs_max = float(np.max(np.abs(rhs)))
    rel_rms = _scalar_rms(error) / max(rhs_rms, 1.0e-300)
    rel_max = float(np.max(np.abs(error))) / max(rhs_max, 1.0e-300)
    source_weight_norm = _trapz(np.abs((radii**exponent) * source), radii)
    moment_relative = abs(independent_moment) / max(source_weight_norm, 1.0e-300)
    stress_scale = max(_scalar_rms(stress), 1.0)
    edge_relative = max(abs(float(stress[0])), abs(float(stress[-1]))) / stress_scale

    return {
        "exponent": exponent,
        "weighted_moment": moment,
        "independent_moment_complement": independent_moment,
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
        "fd4_interior_count": int(r_mid.size),
    }


def audit_actual_radial_stress_closure() -> dict[str, Any]:
    levels: list[dict[str, Any]] = []
    for radial_count in RADIAL_COUNTS:
        actual = _actual_ring_mean_sources(radial_count)
        r = np.asarray(actual["radii"], dtype=float)
        field = actual["field"]
        theta = _audit_channel(
            r,
            np.asarray(actual["theta_source"], dtype=float),
            exponent=2,
            center=field.radial_center,
            halfwidth=field.radial_halfwidth,
        )
        axial = _audit_channel(
            r,
            np.asarray(actual["axial_source"], dtype=float),
            exponent=1,
            center=field.radial_center,
            halfwidth=field.radial_halfwidth,
        )
        levels.append(
            {
                "radial_count": radial_count,
                "radial_spacing": float(r[1] - r[0]),
                "cartesian_fd4_step": FD4_STEP,
                "boundary_margin": float(actual["boundary_margin"]),
                "theta_e2": theta,
                "axial_e1": axial,
            }
        )

    channel_names = ("theta_e2", "axial_e1")
    convergence: dict[str, Any] = {}
    failed_guards: list[str] = []
    for channel in channel_names:
        errors = [float(level[channel]["operator_relative_rms"]) for level in levels]
        ratios = [errors[i] / max(errors[i + 1], 1.0e-300) for i in range(len(errors) - 1)]
        convergence[channel] = {
            "relative_rms_by_radial_count": errors,
            "refinement_ratios": ratios,
            "finest_relative_rms": errors[-1],
            "finest_relative_max": float(levels[-1][channel]["operator_relative_max"]),
            "finest_moment_complement_relative": float(
                levels[-1][channel]["moment_complement_relative"]
            ),
            "finest_edge_relative": float(levels[-1][channel]["edge_relative"]),
        }
        if errors[-1] > FINEST_RELATIVE_RMS_MAX:
            failed_guards.append(f"{channel}:finest_relative_rms")
        if levels[-1][channel]["operator_relative_max"] > FINEST_RELATIVE_MAX_MAX:
            failed_guards.append(f"{channel}:finest_relative_max")
        if min(ratios) < MIN_REFINEMENT_RATIO:
            failed_guards.append(f"{channel}:minimum_refinement_ratio")
        if levels[-1][channel]["moment_complement_relative"] > MOMENT_RELATIVE_MAX:
            failed_guards.append(f"{channel}:moment_complement_relative")
        if levels[-1][channel]["edge_relative"] > EDGE_RELATIVE_MAX:
            failed_guards.append(f"{channel}:edge_relative")

    passed = not failed_guards
    return {
        "task": TASK,
        "schema": SCHEMA,
        "provenance": {
            "parent_agent3_mean_debt_pr": PARENT_AGENT3_MEAN_DEBT_PR,
            "parent_agent3_mean_debt_head": PARENT_AGENT3_MEAN_DEBT_HEAD,
            "admitted_agent2_pr": ADMITTED_AGENT2_PR,
            "admitted_agent2_head": ADMITTED_AGENT2_HEAD,
            "independent_agent4_pr": INDEPENDENT_AGENT4_PR,
            "independent_agent4_head": INDEPENDENT_AGENT4_HEAD,
            "velocity_provider": "openai_ns_reconstruction.kokuno_public_z_pullback_velocity:velocity_osc",
            "actual_defect_operator": "openai_ns_reconstruction.kokuno_actual_oscillatory_mean_stress:_fd4_spatial_operator",
            "stress_constructor": "openai_ns_reconstruction.kokuno_actual_oscillatory_mean_stress:_compact_radial_stress",
            "independent_radial_derivative": "new centered uniform-grid FD4, not used by stress construction",
        },
        "frozen_state": {
            "radial_counts": list(RADIAL_COUNTS),
            "angular_count": ANGULAR_COUNT,
            "time": TIME,
            "z": Z,
            "cartesian_fd4_step": FD4_STEP,
            "nu": NU,
        },
        "frozen_guards": {
            "finest_relative_rms_max": FINEST_RELATIVE_RMS_MAX,
            "finest_relative_max_max": FINEST_RELATIVE_MAX_MAX,
            "minimum_refinement_ratio": MIN_REFINEMENT_RATIO,
            "moment_complement_relative_max": MOMENT_RELATIVE_MAX,
            "edge_relative_max": EDGE_RELATIVE_MAX,
        },
        "levels": levels,
        "convergence": convergence,
        "failed_guards": failed_guards,
        "radial_stress_closure_preflight_passed": passed,
        "truth_boundary": {
            "real_oscillatory_self_defect_component_consumed": True,
            "compact_radial_stress_reconstruction_audited": True,
            "radial_stress_closure_preflight_passed": passed,
            "full_same_cycle_composite_requested_stress_materialized": False,
            "agent1_leading_cross_terms_included": False,
            "matched_pressure_included": False,
            "restricted_forcing_included": False,
            "candidate_finite_head_mean_debt_materialized": False,
            "signed_mean_inverse_input_ready": False,
            "public_velocity_correction_materialized": False,
            "finite_correction_cycle_rerun_allowed": False,
            "finite_correction_cycle_run": False,
            "heldout_ns_residual_assessed": False,
            "residual_reduction_claimed": False,
            "pde_validated": False,
            "paper_exact": False,
            "blowup_proved": False,
        },
    }


def write_receipt(output: str | Path) -> dict[str, Any]:
    receipt = audit_actual_radial_stress_closure()
    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        default="artifacts/kokuno_agent3_radial_stress_closure_v1.json",
    )
    args = parser.parse_args()
    print(json.dumps(write_receipt(args.output), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
