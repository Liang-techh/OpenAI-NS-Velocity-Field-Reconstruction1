"""Independent Agent-4 black-box audit of Agent-2's carrier-resolved provider.

This module preserves the scientific protocol used by Agent 4 on PR #543 while
rebinding only the public candidate under test to Agent 2 PR #551.  The primary
oracle consumes only ``velocity_osc(x,y,z,t)`` from the carrier-resolved public
module.  Cartesian spatial derivatives are reconstructed independently with a
centered fourth-order stencil on held-out off-grid points; no Agent-2 phase,
complete-curl, mass, training-loss, pressure, forcing, or derivative helper is
used as a validation oracle.

This remains an oscillatory-component preflight, not a full Navier--Stokes
validation.  The final project gates stay fixed at normalized momentum max/L2
<= 1e-3 and divergence max/L2 <= 1e-5 and are not assessed here because the
global matched leading velocity/pressure and materialized correction are not yet
available together.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Callable

import numpy as np

from .kokuno_public_carrier_resolved_velocity import (
    KokunoCarrierResolvedCandidateOscillatoryVelocity as Candidate,
    velocity_osc,
)

SCHEMA = "kokuno-agent4-carrier-resolved-independent-audit-v1"
AUDITED_AGENT2_HEAD = "92852046e6e1ec0a5889b53f989df3ea5d79cb3b"
FROZEN_PROTOCOL_ORIGIN_AGENT4_HEAD = "9498209376fd9dc5abf228248967958759e9f64a"
SEED = 9173241
FD4_STEPS = (0.02, 0.01, 0.005)
PHASE_RESOLUTIONS = (32, 64, 128)

# Frozen from the prior Agent-4 black-box audit before evaluating this candidate.
GUARDS = {
    "finest_relative_divergence_rms": 2.0e-5,
    "finest_relative_divergence_max": 1.0e-4,
    "minimum_divergence_refinement_ratio": 3.0,
    "minimum_covariance_rank_ratio": 2.0e-2,
    "maximum_covariance_resolution_drift": 1.0e-6,
    "axis_near_absolute_max": 1.0e-14,
    "radial_exterior_absolute_max": 1.0e-14,
    "project_axial_exterior_absolute_max": 1.0e-12,
    "minimum_nontrivial_velocity_rms": 1.0e-8,
    "minimum_parameter_perturbation_relative_change": 1.0e-4,
    "minimum_divergence_mutation_relative_rms": 5.0e-2,
    "maximum_duplicated_covariance_rank_ratio": 1.0e-8,
}


def _rms(values: np.ndarray) -> float:
    arr = np.asarray(values)
    return float(np.sqrt(np.mean(np.abs(arr) ** 2)))


def _vector_rms(values: np.ndarray) -> float:
    arr = np.asarray(values)
    return float(np.sqrt(np.mean(np.sum(np.abs(arr) ** 2, axis=-1))))


def _relative_change(actual: np.ndarray, reference: np.ndarray) -> float:
    a = np.asarray(actual)
    r = np.asarray(reference)
    return float(np.linalg.norm((a - r).ravel()) / max(np.linalg.norm(r.ravel()), 1.0e-300))


def _fd4_jacobian(
    points: np.ndarray,
    times: np.ndarray,
    field_fn: Callable[[Any, Any, Any, Any], np.ndarray],
    step: float,
) -> np.ndarray:
    """Return du_i/dx_j using an independent centered Cartesian FD4 stencil."""
    points = np.asarray(points, dtype=float)
    times = np.asarray(times, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError("points must have shape (n,3)")
    if times.shape != (points.shape[0],):
        raise ValueError("times must have shape (n,)")
    if not math.isfinite(step) or step <= 0.0:
        raise ValueError("step must be positive and finite")

    jac = np.empty((points.shape[0], 3, 3), dtype=float)
    for axis in range(3):
        offset = np.zeros(3, dtype=float)
        offset[axis] = step
        p2 = points + 2.0 * offset
        p1 = points + offset
        m1 = points - offset
        m2 = points - 2.0 * offset
        f_p2 = np.asarray(field_fn(p2[:, 0], p2[:, 1], p2[:, 2], times), dtype=float)
        f_p1 = np.asarray(field_fn(p1[:, 0], p1[:, 1], p1[:, 2], times), dtype=float)
        f_m1 = np.asarray(field_fn(m1[:, 0], m1[:, 1], m1[:, 2], times), dtype=float)
        f_m2 = np.asarray(field_fn(m2[:, 0], m2[:, 1], m2[:, 2], times), dtype=float)
        jac[:, :, axis] = (-f_p2 + 8.0 * f_p1 - 8.0 * f_m1 + f_m2) / (12.0 * step)
    return jac


def _divergence_metrics(jac: np.ndarray) -> dict[str, float]:
    jac = np.asarray(jac, dtype=float)
    div = np.trace(jac, axis1=1, axis2=2)
    point_grad = np.sqrt(np.sum(jac * jac, axis=(1, 2)))
    grad_rms = max(_rms(point_grad), 1.0e-300)
    grad_max = max(float(np.max(point_grad, initial=0.0)), 1.0e-300)
    return {
        "absolute_rms": _rms(div),
        "absolute_max": float(np.max(np.abs(div), initial=0.0)),
        "relative_rms": _rms(div) / grad_rms,
        "relative_max": float(np.max(np.abs(div), initial=0.0)) / grad_max,
        "gradient_rms": grad_rms,
        "gradient_max": grad_max,
    }


def _heldout_points() -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(SEED)
    n = 24
    radius = rng.uniform(0.42, 1.08, size=n)
    theta = rng.uniform(-math.pi, math.pi, size=n)
    z = rng.uniform(-0.72, 0.72, size=n)
    t = rng.uniform(0.31, 0.69, size=n)
    points = np.stack((radius * np.cos(theta), radius * np.sin(theta), z), axis=-1)
    return points, t


def _cylindrical_covariance(
    field: Candidate,
    *,
    radius: float,
    z: float,
    time: float,
    resolution: int,
) -> np.ndarray:
    theta = 2.0 * math.pi * (np.arange(resolution, dtype=float) + 0.37) / float(resolution)
    x = radius * np.cos(theta)
    y = radius * np.sin(theta)
    vel = np.asarray(
        field.velocity_osc(x, y, np.full_like(theta, z), np.full_like(theta, time)),
        dtype=float,
    )
    ur = np.cos(theta) * vel[:, 0] + np.sin(theta) * vel[:, 1]
    uth = -np.sin(theta) * vel[:, 0] + np.cos(theta) * vel[:, 1]
    uz = vel[:, 2]
    return np.asarray((np.mean(ur * uth), np.mean(ur * uz)), dtype=float)


def _covariance_jacobian(resolution: int, delta: float = 1.0e-4) -> np.ndarray:
    base = dict(
        h=0.005,
        radial_center=0.75,
        radial_halfwidth=0.60,
        time_min=0.25,
        time_max=0.75,
        normal_target_ratio=0.40,
        cross_target_ratio=0.05,
        time_modulation=0.08,
        mode_imaginary_ratio=0.17,
        pulse_samples=129,
        transverse_samples=257,
    )
    columns: list[np.ndarray] = []
    for key in ("normal_target_ratio", "cross_target_ratio"):
        plus = dict(base)
        minus = dict(base)
        plus[key] += delta
        minus[key] -= delta
        c_plus = _cylindrical_covariance(
            Candidate(**plus), radius=0.73, z=0.17, time=0.43, resolution=resolution
        )
        c_minus = _cylindrical_covariance(
            Candidate(**minus), radius=0.73, z=0.17, time=0.43, resolution=resolution
        )
        columns.append((c_plus - c_minus) / (2.0 * delta))
    return np.stack(columns, axis=1)


def _support_metrics() -> dict[str, float]:
    axis_r = np.asarray((0.0, 0.02, 0.08, 0.14), dtype=float)
    axis_theta = np.asarray((0.0, 0.7, -1.4, 2.1), dtype=float)
    axis_vel = np.asarray(
        velocity_osc(
            axis_r * np.cos(axis_theta),
            axis_r * np.sin(axis_theta),
            np.asarray((0.0, 0.2, -0.3, 0.6)),
            np.asarray((0.31, 0.41, 0.57, 0.68)),
        ),
        dtype=float,
    )

    outer_r = np.asarray((1.36, 1.50, 1.80, 1.95), dtype=float)
    outer_theta = np.asarray((0.2, -0.9, 1.7, -2.4), dtype=float)
    outer_vel = np.asarray(
        velocity_osc(
            outer_r * np.cos(outer_theta),
            outer_r * np.sin(outer_theta),
            np.asarray((0.1, -0.4, 0.7, -0.8)),
            np.asarray((0.33, 0.48, 0.59, 0.67)),
        ),
        dtype=float,
    )

    axial_r = np.asarray((0.48, 0.66, 0.84, 1.02, 0.57, 0.93, 0.74, 1.10), dtype=float)
    axial_theta = np.asarray((0.1, 0.7, -1.2, 2.0, -2.4, 1.4, -0.5, 2.6), dtype=float)
    axial_z = np.asarray((2.05, 2.25, 2.45, 2.70, -2.05, -2.25, -2.45, -2.70), dtype=float)
    axial_t = np.asarray((0.31, 0.37, 0.43, 0.49, 0.53, 0.59, 0.63, 0.69), dtype=float)
    axial_vel = np.asarray(
        velocity_osc(axial_r * np.cos(axial_theta), axial_r * np.sin(axial_theta), axial_z, axial_t),
        dtype=float,
    )
    return {
        "axis_near_absolute_max": float(np.max(np.abs(axis_vel), initial=0.0)),
        "radial_exterior_absolute_max": float(np.max(np.abs(outer_vel), initial=0.0)),
        "project_axial_exterior_absolute_max": float(np.max(np.abs(axial_vel), initial=0.0)),
        "project_axial_exterior_velocity_rms": _vector_rms(axial_vel),
    }


def run_audit() -> dict[str, Any]:
    points, times = _heldout_points()
    base_vel = np.asarray(velocity_osc(points[:, 0], points[:, 1], points[:, 2], times), dtype=float)
    nontrivial_rms = _vector_rms(base_vel)

    div_rows: list[dict[str, float]] = []
    jac_finest: np.ndarray | None = None
    for step in FD4_STEPS:
        jac = _fd4_jacobian(points, times, velocity_osc, step)
        div_rows.append({"step": float(step), **_divergence_metrics(jac)})
        jac_finest = jac
    assert jac_finest is not None
    refinement = [
        div_rows[i]["relative_rms"] / max(div_rows[i + 1]["relative_rms"], 1.0e-300)
        for i in range(len(div_rows) - 1)
    ]

    lam = max(float(div_rows[-1]["gradient_rms"]), 1.0) * 0.25
    mutated_div = np.trace(jac_finest, axis1=1, axis2=2) + lam
    mutation_relative_rms = _rms(mutated_div) / max(float(div_rows[-1]["gradient_rms"]), 1.0e-300)

    covariance_rows: list[dict[str, Any]] = []
    jacobians: list[np.ndarray] = []
    for resolution in PHASE_RESOLUTIONS:
        cov_jac = _covariance_jacobian(int(resolution))
        singular = np.linalg.svd(cov_jac, compute_uv=False)
        ratio = float(singular[-1] / max(singular[0], 1.0e-300))
        covariance_rows.append(
            {
                "phase_resolution": int(resolution),
                "singular_values": [float(v) for v in singular],
                "rank_ratio": ratio,
                "jacobian": cov_jac.tolist(),
            }
        )
        jacobians.append(cov_jac)
    cov_ref = jacobians[-1]
    covariance_drift = [
        float(np.linalg.norm(j - cov_ref) / max(np.linalg.norm(cov_ref), 1.0e-300))
        for j in jacobians[:-1]
    ]
    duplicated = np.stack((cov_ref[:, 0], cov_ref[:, 0]), axis=1)
    dup_singular = np.linalg.svd(duplicated, compute_uv=False)
    duplicated_rank_ratio = float(dup_singular[-1] / max(dup_singular[0], 1.0e-300))

    support = _support_metrics()
    perturb_points = points[:12]
    perturb_times = times[:12]
    low = Candidate(h=0.0045).velocity_osc(
        perturb_points[:, 0], perturb_points[:, 1], perturb_points[:, 2], perturb_times
    )
    high = Candidate(h=0.0055).velocity_osc(
        perturb_points[:, 0], perturb_points[:, 1], perturb_points[:, 2], perturb_times
    )
    ref = base_vel[:12]
    parameter_changes = {
        "h_minus_10pct_relative_change": _relative_change(low, ref),
        "h_plus_10pct_relative_change": _relative_change(high, ref),
    }

    finest = div_rows[-1]
    local_divergence_passed = bool(
        finest["relative_rms"] <= GUARDS["finest_relative_divergence_rms"]
        and finest["relative_max"] <= GUARDS["finest_relative_divergence_max"]
        and min(refinement) >= GUARDS["minimum_divergence_refinement_ratio"]
    )
    covariance_passed = bool(
        min(row["rank_ratio"] for row in covariance_rows) >= GUARDS["minimum_covariance_rank_ratio"]
        and max(covariance_drift) <= GUARDS["maximum_covariance_resolution_drift"]
        and duplicated_rank_ratio <= GUARDS["maximum_duplicated_covariance_rank_ratio"]
    )
    radial_axis_support_passed = bool(
        support["axis_near_absolute_max"] <= GUARDS["axis_near_absolute_max"]
        and support["radial_exterior_absolute_max"] <= GUARDS["radial_exterior_absolute_max"]
    )
    project_axial_support_passed = bool(
        support["project_axial_exterior_absolute_max"] <= GUARDS["project_axial_exterior_absolute_max"]
    )
    nontriviality_passed = bool(nontrivial_rms >= GUARDS["minimum_nontrivial_velocity_rms"])
    parameter_perturbation_passed = bool(
        min(parameter_changes.values()) >= GUARDS["minimum_parameter_perturbation_relative_change"]
    )
    mutation_passed = bool(mutation_relative_rms >= GUARDS["minimum_divergence_mutation_relative_rms"])
    local_curl_covariance_preflight_passed = bool(
        local_divergence_passed
        and covariance_passed
        and radial_axis_support_passed
        and nontriviality_passed
        and parameter_perturbation_passed
        and mutation_passed
    )
    project_support_preflight_passed = bool(radial_axis_support_passed and project_axial_support_passed)
    public_oscillatory_preflight_passed = bool(
        local_curl_covariance_preflight_passed and project_support_preflight_passed
    )

    return {
        "schema": SCHEMA,
        "audited_agent2_head": AUDITED_AGENT2_HEAD,
        "frozen_protocol_origin_agent4_head": FROZEN_PROTOCOL_ORIGIN_AGENT4_HEAD,
        "seed": SEED,
        "protocol": {
            "primary_interface": "carrier-resolved module-level velocity_osc(x,y,z,t) only",
            "heldout_offgrid_points": int(points.shape[0]),
            "fd4_steps": list(FD4_STEPS),
            "phase_resolutions": list(PHASE_RESOLUTIONS),
            "scientific_guards_changed_from_agent4_543": False,
            "pressure_or_forcing_fit": False,
            "candidate_internal_derivative_oracle_used": False,
            "training_tensor_or_loss_used": False,
            "project_registered_support": "r<2 and |z|<2, smooth zero extension",
        },
        "guards_frozen_before_actions": dict(GUARDS),
        "divergence": {
            "rows": div_rows,
            "relative_rms_refinement_ratios": refinement,
            "local_divergence_passed": local_divergence_passed,
        },
        "divergence_mutation": {
            "mutation": "external u_x += lambda*x",
            "lambda": lam,
            "relative_rms": mutation_relative_rms,
            "caught": mutation_passed,
        },
        "phase_mean_covariance": {
            "definition": "C=(mean(u_r*u_theta), mean(u_r*u_z)); columns are public normal/cross target finite differences",
            "rows": covariance_rows,
            "relative_drift_to_finest": covariance_drift,
            "duplicated_first_column_rank_ratio": duplicated_rank_ratio,
            "covariance_rank_preflight_passed": covariance_passed,
        },
        "support": {
            **support,
            "radial_axis_support_passed": radial_axis_support_passed,
            "project_axial_support_passed": project_axial_support_passed,
        },
        "nontriviality": {"heldout_velocity_rms": nontrivial_rms, "passed": nontriviality_passed},
        "parameter_perturbation": {**parameter_changes, "passed": parameter_perturbation_passed},
        "local_curl_covariance_preflight_passed": local_curl_covariance_preflight_passed,
        "project_support_preflight_passed": project_support_preflight_passed,
        "public_oscillatory_preflight_passed": public_oscillatory_preflight_passed,
        "truth_boundary": {
            "oscillatory_component_only": True,
            "global_leading_velocity_pressure_available": False,
            "materialized_correction_available": False,
            "formal_full_domain_pde_gate_assessed": False,
            "heldout_ns_residual_assessed": False,
            "pde_validated": False,
            "formal_momentum_normalized_max_l2_gate": 1.0e-3,
            "formal_divergence_max_l2_gate": 1.0e-5,
        },
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    report = run_audit()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
