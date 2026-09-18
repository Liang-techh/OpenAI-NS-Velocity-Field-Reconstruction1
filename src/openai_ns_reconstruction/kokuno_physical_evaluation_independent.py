"""Independent black-box audit of Agent 2's auxiliary-torus pullback.

The candidate-side implementation exposes the physical pullback Y(r,t) and
its post-evaluation chain-rule derivatives.  This audit deliberately avoids the
Agent-2 FD4 regression.  It treats ``torus_phase`` as a black-box value map,
composes a fresh periodic scalar field on the torus, and reconstructs physical
r/t/z derivatives with an external centered FD6 operator at three frozen step
sizes.

This is a local structural/interface preflight.  It is not a Navier--Stokes
momentum residual, does not bind the missing positive-order/background source
family, and cannot set ``pde_validated=true``.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_physical_evaluation import KokunoPhysicalEvaluationMap

SCHEMA = "kokuno-agent4-physical-evaluation-independent-audit-v1"
TASK_ID = "KOKUNO-A4-PHYSICAL-EVALUATION-INDEPENDENT-AUDIT-023"
BASE_PR = 410
BASE_HEAD = "482c86c9e8cf10dc87e2c857c97112e74c28809f"
SEED = 9173141
STEPS = (0.008, 0.004, 0.002)
RANDOM_POINT_COUNT = 18
AXIS_NEAR_RADII = (0.031, 0.039, 0.052, 0.071, 0.095)
SEAM_RADII = (0.17, 0.29, 0.47, 0.73)
SEAM_TARGETS = (1.0e-4, 1.0 - 1.0e-4, 2.0e-4, 1.0 - 2.0e-4)

# Frozen before execution.  These are local derivative-contract guards, not
# replacements for the formal project PDE gates.
FINE_RELATIVE_RMS_GUARD = 2.0e-7
MIN_REFINEMENT_RATIO_GUARD = 6.0
OMIT_TORUS_MUTATION_RELATIVE_RMS_FLOOR = 5.0e-3
FORMAL_MOMENTUM_GATE = 1.0e-3
FORMAL_DIVERGENCE_GATE = 1.0e-5

CASES = (
    {"v_r": (0.31, -0.27), "v_t": (0.23, 0.41), "d_r": 1.37, "h": 0.004},
    {"v_r": (-0.82, 0.47), "v_t": (0.61, -0.73), "d_r": 0.72, "h": 0.006},
    {"v_r": (0.44, 0.91), "v_t": (-0.57, 0.36), "d_r": 2.15, "h": 0.0025},
)


def _relative_rms(error: Any, reference: Any) -> float:
    err = np.asarray(error, dtype=float)
    ref = np.asarray(reference, dtype=float)
    denom = max(float(np.sqrt(np.mean(ref * ref))), 1.0e-14)
    return float(np.sqrt(np.mean(err * err)) / denom)


def _max_relative_to_rms_scale(error: Any, reference: Any) -> float:
    err = np.asarray(error, dtype=float)
    ref = np.asarray(reference, dtype=float)
    denom = max(float(np.sqrt(np.mean(ref * ref))), 1.0e-14)
    return float(np.max(np.abs(err)) / denom)


def _fd6(fun, x: float, step: float) -> float:
    """Centered sixth-order first derivative, independent of Agent-2 FD4."""
    h = float(step)
    return float(
        (
            -fun(x - 3.0 * h)
            + 9.0 * fun(x - 2.0 * h)
            - 45.0 * fun(x - h)
            + 45.0 * fun(x + h)
            - 9.0 * fun(x + 2.0 * h)
            + fun(x + 3.0 * h)
        )
        / (60.0 * h)
    )


def _lifted_field_and_partials(
    op: KokunoPhysicalEvaluationMap,
    r: Any,
    z: Any,
    t: Any,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Fresh periodic torus field plus partials holding Y fixed.

    Integer torus frequencies make this field smooth across the mod-1 seams,
    which lets the independent finite-difference stencil cross those seams.
    """
    rr, zz, tt = np.broadcast_arrays(
        np.asarray(r, dtype=float), np.asarray(z, dtype=float), np.asarray(t, dtype=float)
    )
    y = op.torus_phase(rr, tt)
    twopi = 2.0 * np.pi
    a1 = twopi * (y[..., 0] + 2.0 * y[..., 1] + 0.17)
    a2 = twopi * (2.0 * y[..., 0] - y[..., 1] - 0.09)
    a3 = twopi * (3.0 * y[..., 0] + y[..., 1] + 0.03)
    value = (
        0.13 * rr * rr
        + 0.07 * zz**3
        - 0.11 * tt * tt
        + 0.20 * rr * tt
        + np.sin(a1)
        + 0.37 * np.cos(a2)
        + 0.11 * np.sin(a3)
    )
    partial_r = 0.26 * rr + 0.20 * tt
    partial_t = -0.22 * tt + 0.20 * rr
    partial_z = 0.21 * zz * zz
    grad_y0 = twopi * np.cos(a1) - 0.74 * twopi * np.sin(a2) + 0.33 * twopi * np.cos(a3)
    grad_y1 = 2.0 * twopi * np.cos(a1) + 0.37 * twopi * np.sin(a2) + 0.11 * twopi * np.cos(a3)
    grad_y = np.stack((grad_y0, grad_y1), axis=-1)
    return value, partial_r, partial_t, partial_z, grad_y


def _lifted_value(op: KokunoPhysicalEvaluationMap, r: float, z: float, t: float) -> float:
    return float(_lifted_field_and_partials(op, r, z, t)[0])


def _build_points(
    rng: np.random.Generator,
    op: KokunoPhysicalEvaluationMap,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    r_random = rng.uniform(0.08, 1.25, RANDOM_POINT_COUNT)
    z_random = rng.uniform(-0.65, 0.65, RANDOM_POINT_COUNT)
    t_random = rng.uniform(-0.45, 0.45, RANDOM_POINT_COUNT)

    r_axis = np.asarray(AXIS_NEAR_RADII, dtype=float)
    z_axis = rng.uniform(-0.4, 0.4, r_axis.size)
    t_axis = rng.uniform(-0.35, 0.35, r_axis.size)

    r_seam = np.asarray(SEAM_RADII, dtype=float)
    z_seam = rng.uniform(-0.5, 0.5, r_seam.size)
    targets = np.asarray(SEAM_TARGETS, dtype=float)
    vt0 = float(op.v_t[0])
    if abs(vt0) < 1.0e-12:
        raise ValueError("audit cases require nonzero first time-torus direction")
    t_seam = []
    for rr, target in zip(r_seam, targets):
        base = rr**op.d_r * float(op.v_r[0])
        integer_shift = round(base - target)
        t_seam.append((target + integer_shift - base) / vt0)

    r = np.concatenate((r_random, r_axis, r_seam))
    z = np.concatenate((z_random, z_axis, z_seam))
    t = np.concatenate((t_random, t_axis, np.asarray(t_seam, dtype=float)))
    seam_mask = np.asarray(
        [False] * RANDOM_POINT_COUNT + [False] * r_axis.size + [True] * r_seam.size,
        dtype=bool,
    )
    if float(np.min(r)) <= 3.0 * max(STEPS):
        raise AssertionError("FD6 radial stencil would cross r=0")
    return r, z, t, seam_mask


def _case_report(
    rng: np.random.Generator,
    case_index: int,
    params: dict[str, Any],
) -> tuple[dict[str, Any], dict[str, Any]]:
    op = KokunoPhysicalEvaluationMap(**params)
    r, z, t, seam_mask = _build_points(rng, op)
    q = rng.uniform(0.22, 1.37, r.size)
    _, partial_r, partial_t, partial_z, grad_y = _lifted_field_and_partials(op, r, z, t)
    public = op.differentiate(
        r=r,
        q=q,
        partial_r=partial_r,
        partial_t=partial_t,
        partial_z=partial_z,
        grad_y=grad_y,
    )

    ladder: list[dict[str, Any]] = []
    for step in STEPS:
        dr_fd = np.asarray(
            [
                _fd6(lambda xx: _lifted_value(op, xx, zz, tt), rr, step)
                for rr, zz, tt in zip(r, z, t)
            ],
            dtype=float,
        )
        dt_fd = np.asarray(
            [
                _fd6(lambda xx: _lifted_value(op, rr, zz, xx), tt, step)
                for rr, zz, tt in zip(r, z, t)
            ],
            dtype=float,
        )
        dz_fd = np.asarray(
            [
                _fd6(lambda xx: _lifted_value(op, rr, xx, tt), zz, step)
                for rr, zz, tt in zip(r, z, t)
            ],
            dtype=float,
        )
        comparisons = {
            "mathsf_r": (public["mathsf_r"], dr_fd),
            "mathsf_t": (public["mathsf_t"], dt_fd),
            "partial_z": (public["partial_z_after_evaluation"], dz_fd),
            "D_r": (public["D_r"], np.sqrt(q) * dr_fd),
            "D_z": (public["D_z"], np.sqrt(q) * dz_fd),
            "mathsf_t_star": (public["mathsf_t_star"], q ** (1.0 + op.h) * dt_fd),
        }
        relative_rms = {
            key: _relative_rms(observed - expected, expected)
            for key, (observed, expected) in comparisons.items()
        }
        max_relative = {
            key: _max_relative_to_rms_scale(observed - expected, expected)
            for key, (observed, expected) in comparisons.items()
        }
        primary_relative_rms = max(
            relative_rms["mathsf_r"], relative_rms["mathsf_t"], relative_rms["partial_z"]
        )
        ladder.append(
            {
                "step": float(step),
                "relative_rms": relative_rms,
                "max_relative_to_rms_scale": max_relative,
                "primary_relative_rms": float(primary_relative_rms),
            }
        )

    primary_errors = [row["primary_relative_rms"] for row in ladder]
    refinement = [
        primary_errors[i] / max(primary_errors[i + 1], np.finfo(float).tiny)
        for i in range(len(primary_errors) - 1)
    ]

    finest_step = STEPS[-1]
    dr_finest = np.asarray(
        [
            _fd6(lambda xx: _lifted_value(op, xx, zz, tt), rr, finest_step)
            for rr, zz, tt in zip(r, z, t)
        ],
        dtype=float,
    )
    dt_finest = np.asarray(
        [
            _fd6(lambda xx: _lifted_value(op, rr, zz, xx), tt, finest_step)
            for rr, zz, tt in zip(r, z, t)
        ],
        dtype=float,
    )
    mutated = op.differentiate(
        r=r,
        q=q,
        partial_r=partial_r,
        partial_t=partial_t,
        partial_z=partial_z,
        grad_y=np.zeros_like(grad_y),
    )
    mutation = {
        "case": case_index,
        "omit_torus_grad_r_relative_rms": _relative_rms(mutated["mathsf_r"] - dr_finest, dr_finest),
        "omit_torus_grad_t_relative_rms": _relative_rms(mutated["mathsf_t"] - dt_finest, dt_finest),
    }

    phase = op.torus_phase(r, t)
    observed_seam = np.any((phase < 5.0e-4) | (phase > 1.0 - 5.0e-4), axis=1)
    report = {
        "case": case_index,
        "parameters": params,
        "point_count": int(r.size),
        "random_off_grid_count": RANDOM_POINT_COUNT,
        "axis_near_count": len(AXIS_NEAR_RADII),
        "designed_torus_seam_count": int(np.sum(seam_mask)),
        "observed_torus_seam_count": int(np.sum(observed_seam)),
        "minimum_radius": float(np.min(r)),
        "maximum_abs_time": float(np.max(np.abs(t))),
        "ladder": ladder,
        "primary_relative_rms_ladder": [float(x) for x in primary_errors],
        "refinement_ratios": [float(x) for x in refinement],
    }
    return report, mutation


def generate_report() -> dict[str, Any]:
    rng = np.random.default_rng(SEED)
    case_reports: list[dict[str, Any]] = []
    mutation_reports: list[dict[str, Any]] = []
    for index, params in enumerate(CASES):
        case_report, mutation_report = _case_report(rng, index, dict(params))
        case_reports.append(case_report)
        mutation_reports.append(mutation_report)

    worst_finest = max(row["primary_relative_rms_ladder"][-1] for row in case_reports)
    worst_refinement = min(min(row["refinement_ratios"]) for row in case_reports)
    weakest_mutation = min(
        min(row["omit_torus_grad_r_relative_rms"], row["omit_torus_grad_t_relative_rms"])
        for row in mutation_reports
    )
    all_seams_hit = all(row["observed_torus_seam_count"] >= len(SEAM_RADII) for row in case_reports)

    checks = {
        "finest_relative_rms": worst_finest <= FINE_RELATIVE_RMS_GUARD,
        "three_level_refinement": worst_refinement >= MIN_REFINEMENT_RATIO_GUARD,
        "omitted_torus_terms_detected": weakest_mutation >= OMIT_TORUS_MUTATION_RELATIVE_RMS_FLOOR,
        "designed_mod1_seams_exercised": all_seams_hit,
        "axis_near_stencil_stays_positive": all(
            row["minimum_radius"] > 3.0 * max(STEPS) for row in case_reports
        ),
    }
    passed = bool(all(checks.values()))
    return {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "base_pr": BASE_PR,
        "base_head": BASE_HEAD,
        "seed": SEED,
        "independent_operator": "periodic black-box torus value composition + centered FD6",
        "sample_contract": {
            "cases": len(CASES),
            "random_off_grid_points_per_case": RANDOM_POINT_COUNT,
            "axis_near_radii": list(AXIS_NEAR_RADII),
            "designed_mod1_seam_radii": list(SEAM_RADII),
            "resolution_ladder": list(STEPS),
            "uses_agent2_fd4": False,
            "uses_training_tensor_or_loss": False,
            "uses_pressure_or_forcing_fit": False,
        },
        "frozen_guards": {
            "finest_relative_rms": FINE_RELATIVE_RMS_GUARD,
            "minimum_refinement_ratio": MIN_REFINEMENT_RATIO_GUARD,
            "omit_torus_mutation_relative_rms_floor": OMIT_TORUS_MUTATION_RELATIVE_RMS_FLOOR,
            "formal_momentum_gate_unchanged": FORMAL_MOMENTUM_GATE,
            "formal_divergence_gate_unchanged": FORMAL_DIVERGENCE_GATE,
        },
        "cases": case_reports,
        "mutation": mutation_reports,
        "summary": {
            "worst_finest_primary_relative_rms": float(worst_finest),
            "worst_refinement_ratio": float(worst_refinement),
            "weakest_omit_torus_mutation_relative_rms": float(weakest_mutation),
            "checks": checks,
            "local_structural_preflight_passed": passed,
            "formal_full_domain_pde_gate_assessed": False,
            "pde_validated": False,
        },
        "truth_boundary": {
            "auxiliary_torus_pullback_independently_preflighted": passed,
            "source_torus_vectors_recovered": False,
            "source_group_constants_recovered": False,
            "actual_positive_order_background_bound": False,
            "actual_auxiliary_torus_mode_family_bound": False,
            "public_xyz_t_velocity_correction_materialized": False,
            "formal_full_domain_pde_gate_assessed": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    report = generate_report()
    text = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    else:
        print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
