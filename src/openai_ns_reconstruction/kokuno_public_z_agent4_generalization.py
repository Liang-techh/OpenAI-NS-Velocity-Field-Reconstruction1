"""Independent Agent-4 cross-seed/generalization audit of the accepted public-z field.

This is deliberately separate from the frozen FD4 protocol in Agent 4 PR #563.
It consumes only the public ``velocity_osc(x,y,z,t)`` evaluator and reconstructs
Cartesian derivatives with a fresh sixth-order centered stencil on a disjoint,
stratified held-out cloud.  The cloud includes ordinary interior points plus
radial- and axial-collar interior points while keeping every derivative stencil
strictly inside the declared compact support.

This is an oscillatory-component robustness check, not a full Navier--Stokes
validation.  It does not read Agent-2 phase/curl/mass internals, training data,
pressure, forcing, or derivative helpers, and it does not change the final
project gates (normalized momentum max/L2 <= 1e-3; divergence max/L2 <= 1e-5).
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Callable

import numpy as np

from .kokuno_public_z_pullback_velocity import (
    KokunoPublicZPullbackCandidateOscillatoryVelocity as Candidate,
    velocity_osc,
)

TASK = "KOKUNO-A4-PUBLIC-Z-GENERALIZATION-041"
SCHEMA = "kokuno-agent4-public-z-generalization-v1"
AUDITED_AGENT2_PR = 561
AUDITED_AGENT2_HEAD = "732800ce4990464b49c8aa32d0dff4580f6684d4"
PARENT_AGENT4_PR = 563
PARENT_AGENT4_HEAD = "9b0f86012c53fa8e32a19f766dbc150931870425"
SEED = 9173251
FD6_STEPS = (0.016, 0.008, 0.004)
POINTS_PER_STRATUM = 16

# Frozen before this audit is evaluated.  These are component-robustness guards,
# not replacements for the formal full-candidate PDE gates.
GUARDS = {
    "finest_relative_divergence_rms": 2.0e-5,
    "finest_relative_divergence_max": 1.0e-4,
    "minimum_divergence_refinement_ratio": 4.0,
    "maximum_stratum_finest_relative_divergence_rms": 5.0e-5,
    "support_exterior_absolute_max": 1.0e-12,
    "minimum_nontrivial_velocity_rms": 1.0e-8,
    "minimum_parameter_perturbation_relative_change": 1.0e-4,
    "minimum_divergence_mutation_relative_rms": 5.0e-2,
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


def _fd6_jacobian(
    points: np.ndarray,
    times: np.ndarray,
    field_fn: Callable[[Any, Any, Any, Any], np.ndarray],
    step: float,
) -> np.ndarray:
    """Return du_i/dx_j with a centered sixth-order finite difference."""
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
        m3 = points - 3.0 * offset
        m2 = points - 2.0 * offset
        m1 = points - offset
        p1 = points + offset
        p2 = points + 2.0 * offset
        p3 = points + 3.0 * offset
        f_m3 = np.asarray(field_fn(m3[:, 0], m3[:, 1], m3[:, 2], times), dtype=float)
        f_m2 = np.asarray(field_fn(m2[:, 0], m2[:, 1], m2[:, 2], times), dtype=float)
        f_m1 = np.asarray(field_fn(m1[:, 0], m1[:, 1], m1[:, 2], times), dtype=float)
        f_p1 = np.asarray(field_fn(p1[:, 0], p1[:, 1], p1[:, 2], times), dtype=float)
        f_p2 = np.asarray(field_fn(p2[:, 0], p2[:, 1], p2[:, 2], times), dtype=float)
        f_p3 = np.asarray(field_fn(p3[:, 0], p3[:, 1], p3[:, 2], times), dtype=float)
        jac[:, :, axis] = (
            -f_m3 + 9.0 * f_m2 - 45.0 * f_m1 + 45.0 * f_p1 - 9.0 * f_p2 + f_p3
        ) / (60.0 * step)
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


def _heldout_stratified_points() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Fresh points disjoint in seed and geometry from the PR #563 cloud."""
    rng = np.random.default_rng(SEED)
    n = POINTS_PER_STRATUM

    # Ordinary support interior.
    r_core = rng.uniform(0.38, 1.12, size=n)
    th_core = rng.uniform(-math.pi, math.pi, size=n)
    z_core = rng.uniform(-0.90, 0.90, size=n)
    t_core = rng.uniform(0.30, 0.70, size=n)

    # Radial collar interior.  With max step .016, all +/-3h probes remain in
    # the public annulus 0.15 < R < 1.35.
    r_radial = np.concatenate(
        (rng.uniform(0.22, 0.30, size=n // 2), rng.uniform(1.20, 1.28, size=n - n // 2))
    )
    th_radial = rng.uniform(-math.pi, math.pi, size=n)
    z_radial = rng.uniform(-0.72, 0.72, size=n)
    t_radial = rng.uniform(0.30, 0.70, size=n)

    # Axial collar interior.  +/-3h probes stay below |z|=2.
    r_axial = rng.uniform(0.42, 1.08, size=n)
    th_axial = rng.uniform(-math.pi, math.pi, size=n)
    zmag = rng.uniform(1.62, 1.86, size=n)
    signs = np.where(np.arange(n) % 2 == 0, 1.0, -1.0)
    z_axial = signs * zmag
    t_axial = rng.uniform(0.30, 0.70, size=n)

    def cart(r: np.ndarray, th: np.ndarray, z: np.ndarray) -> np.ndarray:
        return np.stack((r * np.cos(th), r * np.sin(th), z), axis=-1)

    points = np.concatenate(
        (cart(r_core, th_core, z_core), cart(r_radial, th_radial, z_radial), cart(r_axial, th_axial, z_axial)),
        axis=0,
    )
    times = np.concatenate((t_core, t_radial, t_axial))
    labels = np.asarray(["core"] * n + ["radial_collar"] * n + ["axial_collar"] * n)
    return points, times, labels


def _support_exterior_metrics() -> dict[str, float]:
    r = np.asarray((0.00, 0.10, 1.40, 1.70, 0.54, 0.82, 0.66, 1.02), dtype=float)
    th = np.asarray((0.0, 0.7, -0.9, 1.8, -2.2, 0.4, 1.2, -1.5), dtype=float)
    z = np.asarray((0.0, 0.3, -0.2, 0.5, 2.04, -2.12, 2.30, -2.45), dtype=float)
    t = np.asarray((0.31, 0.36, 0.42, 0.48, 0.53, 0.59, 0.64, 0.69), dtype=float)
    vel = np.asarray(velocity_osc(r * np.cos(th), r * np.sin(th), z, t), dtype=float)
    return {
        "absolute_max": float(np.max(np.abs(vel), initial=0.0)),
        "velocity_rms": _vector_rms(vel),
    }


def run_audit() -> dict[str, Any]:
    points, times, labels = _heldout_stratified_points()
    base_vel = np.asarray(velocity_osc(points[:, 0], points[:, 1], points[:, 2], times), dtype=float)
    nontrivial_rms = _vector_rms(base_vel)

    rows: list[dict[str, Any]] = []
    finest_jac: np.ndarray | None = None
    for step in FD6_STEPS:
        jac = _fd6_jacobian(points, times, velocity_osc, step)
        overall = _divergence_metrics(jac)
        strata = {
            name: _divergence_metrics(jac[labels == name])
            for name in ("core", "radial_collar", "axial_collar")
        }
        rows.append({"step": float(step), **overall, "strata": strata})
        finest_jac = jac
    assert finest_jac is not None

    refinement = [
        rows[i]["relative_rms"] / max(rows[i + 1]["relative_rms"], 1.0e-300)
        for i in range(len(rows) - 1)
    ]

    finest = rows[-1]
    finest_stratum_rms = {
        name: float(finest["strata"][name]["relative_rms"])
        for name in ("core", "radial_collar", "axial_collar")
    }

    # External divergence mutation u_x -> u_x + lambda*x.  Its exact added
    # divergence is lambda, so detector sensitivity is checked without touching
    # the candidate implementation.
    lam = max(float(finest["gradient_rms"]), 1.0) * 0.25
    mutated_div = np.trace(finest_jac, axis1=1, axis2=2) + lam
    mutation_relative_rms = _rms(mutated_div) / max(float(finest["gradient_rms"]), 1.0e-300)

    # Public-parameter sensitivity on a subset of the fresh cloud.
    subset = slice(0, 18)
    ref = base_vel[subset]
    low = np.asarray(
        Candidate(h=0.0045).velocity_osc(
            points[subset, 0], points[subset, 1], points[subset, 2], times[subset]
        ),
        dtype=float,
    )
    high = np.asarray(
        Candidate(h=0.0055).velocity_osc(
            points[subset, 0], points[subset, 1], points[subset, 2], times[subset]
        ),
        dtype=float,
    )
    parameter_changes = {
        "h_minus_10pct_relative_change": _relative_change(low, ref),
        "h_plus_10pct_relative_change": _relative_change(high, ref),
    }

    support = _support_exterior_metrics()
    passed = bool(
        float(finest["relative_rms"]) <= GUARDS["finest_relative_divergence_rms"]
        and float(finest["relative_max"]) <= GUARDS["finest_relative_divergence_max"]
        and min(refinement) >= GUARDS["minimum_divergence_refinement_ratio"]
        and max(finest_stratum_rms.values())
        <= GUARDS["maximum_stratum_finest_relative_divergence_rms"]
        and support["absolute_max"] <= GUARDS["support_exterior_absolute_max"]
        and nontrivial_rms >= GUARDS["minimum_nontrivial_velocity_rms"]
        and min(parameter_changes.values())
        >= GUARDS["minimum_parameter_perturbation_relative_change"]
        and mutation_relative_rms >= GUARDS["minimum_divergence_mutation_relative_rms"]
    )

    failed: list[str] = []
    if float(finest["relative_rms"]) > GUARDS["finest_relative_divergence_rms"]:
        failed.append("finest_relative_divergence_rms")
    if float(finest["relative_max"]) > GUARDS["finest_relative_divergence_max"]:
        failed.append("finest_relative_divergence_max")
    if min(refinement) < GUARDS["minimum_divergence_refinement_ratio"]:
        failed.append("minimum_divergence_refinement_ratio")
    if max(finest_stratum_rms.values()) > GUARDS["maximum_stratum_finest_relative_divergence_rms"]:
        failed.append("maximum_stratum_finest_relative_divergence_rms")
    if support["absolute_max"] > GUARDS["support_exterior_absolute_max"]:
        failed.append("support_exterior_absolute_max")
    if nontrivial_rms < GUARDS["minimum_nontrivial_velocity_rms"]:
        failed.append("minimum_nontrivial_velocity_rms")
    if min(parameter_changes.values()) < GUARDS["minimum_parameter_perturbation_relative_change"]:
        failed.append("minimum_parameter_perturbation_relative_change")
    if mutation_relative_rms < GUARDS["minimum_divergence_mutation_relative_rms"]:
        failed.append("minimum_divergence_mutation_relative_rms")

    return {
        "task": TASK,
        "schema": SCHEMA,
        "audited_agent2_pr": AUDITED_AGENT2_PR,
        "audited_agent2_head": AUDITED_AGENT2_HEAD,
        "parent_agent4_pr": PARENT_AGENT4_PR,
        "parent_agent4_head": PARENT_AGENT4_HEAD,
        "seed": SEED,
        "point_count": int(points.shape[0]),
        "points_per_stratum": POINTS_PER_STRATUM,
        "fd6_steps": list(FD6_STEPS),
        "guards_frozen_before_actions": GUARDS,
        "divergence": {
            "rows": rows,
            "refinement_ratios": refinement,
            "finest_stratum_relative_rms": finest_stratum_rms,
            "mutation_relative_rms": mutation_relative_rms,
        },
        "support_exterior": support,
        "heldout_velocity_rms": nontrivial_rms,
        "parameter_perturbations": parameter_changes,
        "failed_guards": failed,
        "public_z_generalization_passed": passed,
        "truth_boundary": {
            "oscillatory_component_only": True,
            "different_derivative_implementation_from_agent2": True,
            "different_derivative_order_from_parent_agent4_fd4": True,
            "training_loss_used": False,
            "candidate_internal_derivative_helper_used": False,
            "pressure_or_forcing_fitted": False,
            "formal_full_domain_pde_gate_assessed": False,
            "heldout_ns_residual_assessed": False,
            "pde_validated": False,
            "formal_momentum_normalized_max_l2_gate": 1.0e-3,
            "formal_divergence_max_l2_gate": 1.0e-5,
        },
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    report = run_audit()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
