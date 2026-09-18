"""Independent Agent-4 audit for the segmented Kokuno leading public contract.

The construction stack now exposes one fail-closed ``velocity(x,y,z,t)`` router for
its reconstructed reference and repaired-I2 stages.  This module treats that router
as a black box for a fresh Cartesian FD2 audit on the reference stage and separately
checks whether the audited I2 heat repair is actually observable through the public
float64 velocity contract.

This is deliberately not a full-domain Navier--Stokes verdict.  The segmented
assembly has unreconstructed radial gaps, no matched I2 pressure, no complete
oscillatory/correction composite and no admissible global forcing contract.  The
fixed project gates are retained, but ``pde_validated`` cannot be promoted here.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_segmented_leading_assembly import KokunoSegmentedLeadingAssembly


SCHEMA = "kokuno-agent4-segmented-leading-independent-audit-v1"
TASK_ID = "KOKUNO-A4-SEGMENTED-PUBLIC-CONTRACT-AUDIT-016"
BASE_PR = 345
BASE_HEAD = "57df8375f134277ebea509847f6b319aa84acfd7"
SEED = 9173071
NU = 0.01
FD_STEPS = (0.004, 0.002, 0.001)
RESIDUAL_GATE = 1.0e-3
DIVERGENCE_GATE = 1.0e-5
MUTATION = 0.02


def _sample_reference_points(
    assembly: KokunoSegmentedLeadingAssembly,
) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Fresh off-grid and axis-near points safely inside the reference stage."""

    rng = np.random.default_rng(SEED)
    D = 0.5 - assembly.h
    times_grid = np.asarray((0.375, 0.5, 0.625), dtype=float)
    regions: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for name, count, x_bounds, eta_bound in (
        ("off_grid", 18, (0.055, 0.285), 0.46),
        ("axis_near", 15, (2.0e-5, 2.5e-3), 0.38),
    ):
        times = rng.choice(times_grid, size=count, replace=True)
        eta = rng.uniform(-eta_bound, eta_bound, size=count)
        q = (1.0 - times) / (1.0 - eta * eta)
        z = np.power(q, D) * eta
        X = rng.uniform(x_bounds[0], x_bounds[1], size=count)
        theta = rng.uniform(0.17, 2.91, size=count)
        radius = np.sqrt(2.0 * q * X)
        points = np.column_stack(
            (radius * np.cos(theta), radius * np.sin(theta), z)
        )
        # Leave a wide margin to the reference/I2 router boundary for every FD stencil.
        routed = assembly._coordinates.evaluate(
            points[:, 0], points[:, 1], points[:, 2], times
        )["X"]
        if not np.all(np.asarray(routed) < 0.32):
            raise RuntimeError("reference held-out sampler left the preregistered interior")
        regions[name] = (points, times)
    return regions


def _velocity(
    assembly: KokunoSegmentedLeadingAssembly,
    points: np.ndarray,
    times: np.ndarray,
    velocity_mutation: float = 0.0,
) -> np.ndarray:
    values = np.asarray(assembly.at_points(points, times), dtype=float).reshape(-1, 3)
    if velocity_mutation:
        values = values.copy()
        values[:, 0] += velocity_mutation * points[:, 0]
    return values


def _pressure(
    assembly: KokunoSegmentedLeadingAssembly,
    points: np.ndarray,
    times: np.ndarray,
    pressure_mutation: float = 0.0,
) -> np.ndarray:
    values = np.asarray(
        assembly.pressure(points[:, 0], points[:, 1], points[:, 2], times),
        dtype=float,
    ).reshape(-1)
    if pressure_mutation:
        values = values + pressure_mutation * points[:, 0]
    return values


def _fd2_reference_operator(
    assembly: KokunoSegmentedLeadingAssembly,
    points: np.ndarray,
    times: np.ndarray,
    step: float,
    *,
    velocity_mutation: float = 0.0,
    pressure_mutation: float = 0.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Centered Cartesian FD2 raw zero-force NS operator on the reference stage."""

    h = float(step)
    velocity0 = _velocity(assembly, points, times, velocity_mutation)
    velocity_t_plus = _velocity(assembly, points, times + h, velocity_mutation)
    velocity_t_minus = _velocity(assembly, points, times - h, velocity_mutation)
    velocity_t = (velocity_t_plus - velocity_t_minus) / (2.0 * h)

    grad_velocity = np.empty((points.shape[0], 3, 3), dtype=float)
    laplacian = np.zeros_like(velocity0)
    grad_pressure = np.empty((points.shape[0], 3), dtype=float)
    pressure0 = _pressure(assembly, points, times, pressure_mutation)

    for axis in range(3):
        plus = points.copy()
        minus = points.copy()
        plus[:, axis] += h
        minus[:, axis] -= h
        velocity_plus = _velocity(assembly, plus, times, velocity_mutation)
        velocity_minus = _velocity(assembly, minus, times, velocity_mutation)
        grad_velocity[:, :, axis] = (velocity_plus - velocity_minus) / (2.0 * h)
        laplacian += (velocity_plus - 2.0 * velocity0 + velocity_minus) / (h * h)

        pressure_plus = _pressure(assembly, plus, times, pressure_mutation)
        pressure_minus = _pressure(assembly, minus, times, pressure_mutation)
        grad_pressure[:, axis] = (pressure_plus - pressure_minus) / (2.0 * h)

    convection = np.einsum("nij,nj->ni", grad_velocity, velocity0)
    residual = velocity_t + convection + grad_pressure - NU * laplacian
    divergence = np.trace(grad_velocity, axis1=1, axis2=2)
    # pressure0 is intentionally evaluated above: it exercises the same public pressure
    # contract at every held-out point even though only its gradient enters the operator.
    if not np.all(np.isfinite(pressure0)):
        raise RuntimeError("non-finite public pressure on reference held-out points")
    return residual, divergence


def _metrics(residual: np.ndarray, divergence: np.ndarray) -> dict[str, Any]:
    vector_norm = np.linalg.norm(residual, axis=1)
    component_rms = np.sqrt(np.mean(residual * residual, axis=0))
    rms = float(np.sqrt(np.mean(vector_norm * vector_norm)))
    return {
        "sample_count": int(residual.shape[0]),
        "max_vector_residual": float(np.max(vector_norm)),
        "sample_l2_vector_residual": rms,
        "normalized_residual": rms,
        "normalization": "registered dimensionless residual scale U_ref^2/L_ref = 1",
        "component_rms": [float(value) for value in component_rms],
        "divergence_max_abs": float(np.max(np.abs(divergence))),
        "divergence_sample_l2": float(np.sqrt(np.mean(divergence * divergence))),
    }


def _relative_change(new: float, old: float) -> float:
    return abs(new - old) / max(abs(new), abs(old), 1.0e-300)


def _solve_q_bisection(z: float, t: float, h: float) -> float:
    """Independent monotone bisection for q-z^2 q^(2h)=1-t."""

    tau = 1.0 - float(t)
    D = 0.5 - float(h)
    if tau <= 0.0:
        raise ValueError("audit probes require t<1")
    q_star = abs(float(z)) ** (1.0 / D)

    def equation(q: float) -> float:
        return q - float(z) * float(z) * q ** (2.0 * h) - tau

    low = q_star
    high = max(q_star + tau + 1.0, 1.0)
    while equation(high) <= 0.0:
        high *= 2.0
        if not math.isfinite(high):
            raise RuntimeError("independent q bracket overflowed")
    for _ in range(180):
        middle = 0.5 * (low + high)
        if equation(middle) > 0.0:
            high = middle
        else:
            low = middle
    return 0.5 * (low + high)


def _i2_probe_points(
    assembly: KokunoSegmentedLeadingAssembly,
) -> tuple[np.ndarray, np.ndarray]:
    outer = assembly.repaired_i2.outer_schedule
    X_star = math.exp(float(outer.log_X_star))
    D = 0.5 - assembly.h
    eta_values = np.asarray((-0.41, -0.13, 0.19, 0.47), dtype=float)
    times = np.asarray((0.375, 0.5, 0.625, 0.4375), dtype=float)
    x_star = np.asarray((0.91, 1.03, 1.24, 1.52), dtype=float)
    theta = np.asarray((0.31, 0.79, 1.37, 2.11), dtype=float)
    q = (1.0 - times) / (1.0 - eta_values * eta_values)
    z = np.power(q, D) * eta_values
    X = X_star * x_star
    radius = np.sqrt(2.0 * q * X)
    points = np.column_stack(
        (radius * np.cos(theta), radius * np.sin(theta), z)
    )
    return points, times


def _independent_i2_base_velocity(
    assembly: KokunoSegmentedLeadingAssembly,
    points: np.ndarray,
    times: np.ndarray,
) -> np.ndarray:
    """Reconstruct the source I2 pure-swirl formula without construction derivatives."""

    outer = assembly.repaired_i2.outer_schedule
    lam = float(outer.lambda_outer)
    h = float(outer.h)
    D = 0.5 - h
    X_star = math.exp(float(outer.log_X_star))
    e_star = math.exp(float(outer.log_e_star))
    result = np.empty((points.shape[0], 3), dtype=float)
    for index, (point, t) in enumerate(zip(points, times)):
        x, y, z = map(float, point)
        q = _solve_q_bisection(z, float(t), h)
        eta = z / q**D
        X = (x * x + y * y) / (2.0 * q)
        x_star = X / X_star
        f_eta = 1.0 / (1.0 + eta * eta)
        E = e_star * f_eta * x_star ** (-0.5 - lam)
        F = E / math.sqrt(2.0 * X)
        factor = q ** (-1.0 - h) * F
        result[index] = (-y * factor, x * factor, 0.0)
    return result


def _i2_public_contract_audit(
    assembly: KokunoSegmentedLeadingAssembly,
) -> dict[str, Any]:
    points, times = _i2_probe_points(assembly)
    public_velocity = np.asarray(assembly.at_points(points, times), dtype=float)

    # This exact comparison is intentionally a *public-contract visibility* check,
    # not an independent PDE oracle.  The uncorrected source patch is an existing
    # public component; equality means the high-precision repair is rounded away by
    # the candidate-facing float64 router at these held-out probes.
    uncorrected_public = np.asarray(
        assembly.repaired_i2.outer_schedule.patch_velocity(
            points[:, 0], points[:, 1], points[:, 2], times
        ),
        dtype=float,
    )
    exact_equal = np.all(public_velocity == uncorrected_public, axis=1)

    independent_base = _independent_i2_base_velocity(assembly, points, times)
    base_norm = np.linalg.norm(independent_base, axis=1)
    independent_relative = np.linalg.norm(public_velocity - independent_base, axis=1) / np.maximum(
        base_norm, 1.0e-300
    )
    speeds = np.linalg.norm(public_velocity, axis=1)
    repair_delta = np.linalg.norm(public_velocity - uncorrected_public, axis=1)
    repair_relative = repair_delta / np.maximum(speeds, 1.0e-300)
    exact_count = int(np.count_nonzero(exact_equal))
    visible = bool(not np.all(exact_equal))
    return {
        "sample_count": int(points.shape[0]),
        "fresh_probe_seed": SEED,
        "public_float64_equals_uncorrected_base_all_probes": bool(np.all(exact_equal)),
        "public_float64_equals_uncorrected_base_fraction": float(np.mean(exact_equal)),
        "router_vs_uncorrected_base_max_abs": float(
            np.max(np.abs(public_velocity - uncorrected_public))
        ),
        "router_vs_uncorrected_base_max_relative_to_public_speed": float(
            np.max(repair_relative)
        ),
        "public_vs_independent_source_base_max_relative": float(
            np.max(independent_relative)
        ),
        "public_speed_min": float(np.min(speeds)),
        "public_speed_max": float(np.max(speeds)),
        "public_float64_heat_repair_visible": visible,
        "module_level_after_heat_repair_attribution_available": False,
        "interpretation": (
            "The high-precision heat repair is independently certified in PR #339. "
            f"At these probes {exact_count}/{points.shape[0]} float64 velocities round "
            "exactly to the uncorrected I2 base and the remaining visibility is only at "
            "float64 rounding scale, so the public float64 router is not a robust contract "
            "for module-level heat-repair attribution."
        ),
    }


def run_audit() -> dict[str, Any]:
    assembly = KokunoSegmentedLeadingAssembly()
    regions = _sample_reference_points(assembly)
    reference: dict[str, Any] = {}
    for region_name, (points, times) in regions.items():
        ladder: list[dict[str, Any]] = []
        for step in FD_STEPS:
            residual, divergence = _fd2_reference_operator(
                assembly, points, times, step
            )
            entry = {"step": step, **_metrics(residual, divergence)}
            ladder.append(entry)
        finest = ladder[-1]
        previous = ladder[-2]
        reference[region_name] = {
            "resolution_ladder": ladder,
            "fine_to_previous_relative_change": {
                "sample_l2_vector_residual": _relative_change(
                    finest["sample_l2_vector_residual"],
                    previous["sample_l2_vector_residual"],
                ),
                "max_vector_residual": _relative_change(
                    finest["max_vector_residual"], previous["max_vector_residual"]
                ),
            },
            "finest_local_residual_gate_met": bool(
                finest["max_vector_residual"] <= RESIDUAL_GATE
                and finest["sample_l2_vector_residual"] <= RESIDUAL_GATE
            ),
            "finest_local_divergence_gate_met": bool(
                finest["divergence_max_abs"] <= DIVERGENCE_GATE
                and finest["divergence_sample_l2"] <= DIVERGENCE_GATE
            ),
        }

    # Calibrate both independent derivative channels on fresh off-grid points.
    points, times = regions["off_grid"]
    base_residual, base_divergence = _fd2_reference_operator(
        assembly, points, times, FD_STEPS[-1]
    )
    _, mutated_divergence = _fd2_reference_operator(
        assembly,
        points,
        times,
        FD_STEPS[-1],
        velocity_mutation=MUTATION,
    )
    pressure_residual, _ = _fd2_reference_operator(
        assembly,
        points,
        times,
        FD_STEPS[-1],
        pressure_mutation=MUTATION,
    )
    divergence_shift = mutated_divergence - base_divergence
    pressure_shift = pressure_residual - base_residual
    mutation = {
        "declared_mutation": MUTATION,
        "mean_divergence_shift": float(np.mean(divergence_shift)),
        "max_divergence_shift_error": float(
            np.max(np.abs(divergence_shift - MUTATION))
        ),
        "mean_pressure_gradient_x_shift": float(np.mean(pressure_shift[:, 0])),
        "max_pressure_gradient_shift_error": float(
            np.max(
                np.abs(
                    pressure_shift
                    - np.asarray((MUTATION, 0.0, 0.0), dtype=float)[None, :]
                )
            )
        ),
    }

    i2 = _i2_public_contract_audit(assembly)
    independent_base_guard = i2["public_vs_independent_source_base_max_relative"] <= 5.0e-10
    mutation_guard = (
        mutation["max_divergence_shift_error"] <= 1.0e-8
        and mutation["max_pressure_gradient_shift_error"] <= 1.0e-8
    )

    return {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "base": {"pr": BASE_PR, "head": BASE_HEAD},
        "seed": SEED,
        "fixed_project_contract": {
            "nu": NU,
            "normalized_momentum_max_gate": RESIDUAL_GATE,
            "normalized_momentum_l2_gate": RESIDUAL_GATE,
            "divergence_max_gate": DIVERGENCE_GATE,
            "divergence_l2_gate": DIVERGENCE_GATE,
            "thresholds_changed": False,
        },
        "independent_operator": {
            "reference_velocity_input": "KokunoSegmentedLeadingAssembly.at_points only",
            "reference_pressure_input": "KokunoSegmentedLeadingAssembly.pressure only",
            "forcing": "zero for this local raw reference-stage diagnostic; no fitted/free forcing",
            "spatial_and_time_derivatives": "independent centered Cartesian FD2",
            "steps": list(FD_STEPS),
            "training_or_construction_derivatives_read": False,
            "training_loss_read": False,
        },
        "reference_stage": reference,
        "mutation_calibration": mutation,
        "i2_public_contract": i2,
        "preregistered_local_guards": {
            "public_vs_independent_I2_source_base_relative_le_5e-10": bool(
                independent_base_guard
            ),
            "mutation_shift_error_le_1e-8": bool(mutation_guard),
            "all_local_implementation_guards_passed": bool(
                independent_base_guard and mutation_guard
            ),
        },
        "routing": {
            "agent1": (
                "The segmented router is independently consumable on its reference stage. "
                "On I2 the high-precision repair changes only a subset of held-out float64 "
                "probes and only at rounding scale, so preserve a precision-aware/rescaled "
                "public contract or explicitly mark module-level after-repair attribution "
                "unavailable. Continue reconstructing the missing global radial stages and "
                "matched pressure."
            ),
            "agent2": (
                "No change: actual-source public oscillatory velocity and an independent "
                "second covariance direction are still required."
            ),
            "agent3": (
                "No correction-cycle promotion: consume a genuinely independent Agent-2 "
                "column only after it passes the frozen spacetime/radial preflight."
            ),
        },
        "truth_boundary": {
            "segmented_public_contract_independently_audited": True,
            "i2_high_precision_heat_repair_prior_independent_audit_pr": 339,
            "i2_float64_public_repair_observable_at_probes": bool(
                i2["public_float64_heat_repair_visible"]
            ),
            "i2_float64_public_repair_robustly_attributable": False,
            "global_leading_profile_reconstructed": False,
            "global_pressure_available": False,
            "complete_leading_oscillatory_correction_composite_available": False,
            "formal_full_domain_pde_gate_assessed": False,
            "normalized_ns_residual_le_1e-3_claimed": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
        },
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        default="artifacts/kokuno_agent4/segmented_leading_independent_audit.json",
    )
    args = parser.parse_args(argv)
    report = run_audit()
    path = Path(args.output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
