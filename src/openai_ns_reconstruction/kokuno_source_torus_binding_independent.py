"""Independent audit of Agent 2's displayed-source auxiliary-torus binding.

This module treats ``KokunoPhysicalEvaluationMap.from_corrected_source_constants``
as a black-box value map.  It independently reconstructs the corrected-reader
constants, composes a fresh periodic torus field, and differentiates only the
black-box values with centered FD6 on three frozen resolutions.

This is a local structural audit, not a Navier--Stokes residual calculation.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_physical_evaluation import KokunoPhysicalEvaluationMap

SCHEMA = "kokuno-agent4-source-torus-binding-independent-audit-v1"
TASK_ID = "KOKUNO-A4-SOURCE-TORUS-BINDING-INDEPENDENT-AUDIT-024"
BASE_PR = 420
BASE_HEAD = "1c77b07a0e8286b8ffa487917f48d54f6017c7e0"
SEED = 9173151
H_CASES = (0.0015, 0.0045, 0.0085)
STEPS = (0.006, 0.003, 0.0015)
RANDOM_POINT_COUNT = 18
AXIS_NEAR_RADII = (0.025, 0.033, 0.047, 0.065, 0.091)
SEAM_RADII = (0.16, 0.28, 0.46, 0.74)
SEAM_TARGETS = (1.0e-4, 1.0 - 1.0e-4, 2.0e-4, 1.0 - 2.0e-4)

# Frozen before execution; local guards only.
PHASE_EMBEDDING_MAX_GUARD = 2.0e-13
SOURCE_BINDING_RELATIVE_GUARD = 2.0e-14
FINE_RELATIVE_RMS_GUARD = 1.0e-7
MIN_REFINEMENT_RATIO_GUARD = 20.0
SIGN_MUTATION_RELATIVE_RMS_FLOOR = 0.5
EXPONENT_MUTATION_RELATIVE_RMS_FLOOR = 1.0e-2
FORMAL_MOMENTUM_GATE = 1.0e-3
FORMAL_DIVERGENCE_GATE = 1.0e-5


def _reference_constants(h: float) -> dict[str, Any]:
    """Independently reconstruct the displayed corrected-reader constants."""
    root2 = math.sqrt(2.0)
    b_g = root2 - 1.0
    lambda_g = 4.0 - root2
    t_g = 4.0 + root2
    rho_g = math.log(lambda_g) / math.log(t_g)
    d_r = 2.0 * ((1.0 + float(h)) * rho_g - float(h) * 1.0e-5)
    return {
        "J_g": np.asarray(((3.0, 1.0), (1.0, 5.0)), dtype=float),
        "b_g": b_g,
        "v_r": np.asarray((1.0, -b_g), dtype=float),
        "v_t": np.asarray((b_g, 1.0), dtype=float),
        "Lambda_g": lambda_g,
        "T_g": t_g,
        "kappa_s": 1.0e-5,
        "rho_g": rho_g,
        "d_r": d_r,
    }


def _relative_rms(error: Any, reference: Any) -> float:
    err = np.asarray(error, dtype=float)
    ref = np.asarray(reference, dtype=float)
    scale = max(float(np.sqrt(np.mean(ref * ref))), 1.0e-14)
    return float(np.sqrt(np.mean(err * err)) / scale)


def _fd6(fun, x: float, step: float) -> float:
    h = float(step)
    return float(
        (-fun(x - 3*h) + 9*fun(x - 2*h) - 45*fun(x - h)
         + 45*fun(x + h) - 9*fun(x + 2*h) + fun(x + 3*h)) / (60*h)
    )


def _value(op: KokunoPhysicalEvaluationMap, r: float, z: float, t: float) -> float:
    y = op.torus_phase(r, t)
    twopi = 2.0 * np.pi
    a1 = twopi * (y[..., 0] + 2.0*y[..., 1] + 0.13)
    a2 = twopi * (3.0*y[..., 0] - y[..., 1] - 0.07)
    a3 = twopi * (2.0*y[..., 0] + 3.0*y[..., 1] + 0.19)
    return float(
        0.09*r*r + 0.04*z**3 - 0.08*t*t + 0.17*r*t
        + np.sin(a1) + 0.29*np.cos(a2) + 0.07*np.sin(a3)
    )


def _analytic_reference(
    op: KokunoPhysicalEvaluationMap,
    r: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
    *,
    v_r: np.ndarray,
    v_t: np.ndarray,
    d_r: float,
) -> tuple[np.ndarray, np.ndarray]:
    y = op.torus_phase(r, t)
    twopi = 2.0 * np.pi
    a1 = twopi * (y[..., 0] + 2.0*y[..., 1] + 0.13)
    a2 = twopi * (3.0*y[..., 0] - y[..., 1] - 0.07)
    a3 = twopi * (2.0*y[..., 0] + 3.0*y[..., 1] + 0.19)
    gy0 = twopi*np.cos(a1) - 0.87*twopi*np.sin(a2) + 0.14*twopi*np.cos(a3)
    gy1 = 2.0*twopi*np.cos(a1) + 0.29*twopi*np.sin(a2) + 0.21*twopi*np.cos(a3)
    gy = np.stack((gy0, gy1), axis=-1)
    partial_r = 0.18*r + 0.17*t
    partial_t = -0.16*t + 0.17*r
    radial_weight = d_r * r ** (d_r - 1.0)
    dr = partial_r + radial_weight * np.einsum("...i,i->...", gy, v_r)
    dt = partial_t + np.einsum("...i,i->...", gy, v_t)
    return dr, dt


def _build_points(rng: np.random.Generator, ref: dict[str, Any]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    r_random = rng.uniform(0.08, 1.25, RANDOM_POINT_COUNT)
    z_random = rng.uniform(-0.65, 0.65, RANDOM_POINT_COUNT)
    t_random = rng.uniform(-0.45, 0.45, RANDOM_POINT_COUNT)

    r_axis = np.asarray(AXIS_NEAR_RADII, dtype=float)
    z_axis = rng.uniform(-0.4, 0.4, r_axis.size)
    t_axis = rng.uniform(-0.35, 0.35, r_axis.size)

    r_seam = np.asarray(SEAM_RADII, dtype=float)
    z_seam = rng.uniform(-0.5, 0.5, r_seam.size)
    t_seam = []
    for rr, target in zip(r_seam, SEAM_TARGETS):
        base = rr ** ref["d_r"] * ref["v_r"][0]
        shift = round(base - target)
        t_seam.append((target + shift - base) / ref["v_t"][0])

    r = np.concatenate((r_random, r_axis, r_seam))
    z = np.concatenate((z_random, z_axis, z_seam))
    t = np.concatenate((t_random, t_axis, np.asarray(t_seam, dtype=float)))
    if float(np.min(r)) <= 3.0 * max(STEPS):
        raise AssertionError("FD6 stencil would cross r=0")
    return r, z, t


def _phase_embedding_error(op: KokunoPhysicalEvaluationMap, r: np.ndarray, t: np.ndarray, ref: dict[str, Any]) -> float:
    observed = op.torus_phase(r, t)
    raw = r[..., None] ** ref["d_r"] * ref["v_r"] + t[..., None] * ref["v_t"]
    obs_embed = np.concatenate((np.cos(2*np.pi*observed), np.sin(2*np.pi*observed)), axis=-1)
    ref_embed = np.concatenate((np.cos(2*np.pi*raw), np.sin(2*np.pi*raw)), axis=-1)
    return float(np.max(np.abs(obs_embed - ref_embed)))


def _case_report(rng: np.random.Generator, h: float) -> dict[str, Any]:
    ref = _reference_constants(h)
    op = KokunoPhysicalEvaluationMap.from_corrected_source_constants(h=h)
    r, z, t = _build_points(rng, ref)

    source_binding_relative = max(
        _relative_rms(np.asarray(op.v_r) - ref["v_r"], ref["v_r"]),
        _relative_rms(np.asarray(op.v_t) - ref["v_t"], ref["v_t"]),
        abs(op.d_r - ref["d_r"]) / abs(ref["d_r"]),
    )
    phase_error = _phase_embedding_error(op, r, t, ref)
    dr_ref, dt_ref = _analytic_reference(op, r, z, t, v_r=ref["v_r"], v_t=ref["v_t"], d_r=ref["d_r"])

    ladder = []
    for step in STEPS:
        dr_fd = np.asarray([_fd6(lambda xx: _value(op, xx, zz, tt), rr, step) for rr, zz, tt in zip(r, z, t)])
        dt_fd = np.asarray([_fd6(lambda xx: _value(op, rr, zz, xx), tt, step) for rr, zz, tt in zip(r, z, t)])
        dr_error = _relative_rms(dr_fd - dr_ref, dr_ref)
        dt_error = _relative_rms(dt_fd - dt_ref, dt_ref)
        ladder.append({"step": float(step), "dr_relative_rms": dr_error, "dt_relative_rms": dt_error,
                       "primary_relative_rms": max(dr_error, dt_error)})

    errors = [row["primary_relative_rms"] for row in ladder]
    refinement = [errors[i] / max(errors[i+1], np.finfo(float).tiny) for i in range(2)]
    finest = STEPS[-1]
    dr_fd = np.asarray([_fd6(lambda xx: _value(op, xx, zz, tt), rr, finest) for rr, zz, tt in zip(r, z, t)])

    wrong_sign_vr = np.asarray((1.0, ref["b_g"]), dtype=float)
    dr_wrong_sign, _ = _analytic_reference(op, r, z, t, v_r=wrong_sign_vr, v_t=ref["v_t"], d_r=ref["d_r"])
    dr_wrong_exp, _ = _analytic_reference(op, r, z, t, v_r=ref["v_r"], v_t=ref["v_t"], d_r=ref["d_r"] + 0.02)

    phase = op.torus_phase(r, t)
    seam_hits = int(np.sum(np.any((phase < 5.0e-4) | (phase > 1.0 - 5.0e-4), axis=1)))
    return {
        "h": h,
        "reference_d_r": ref["d_r"],
        "public_d_r": op.d_r,
        "binding": op.binding,
        "point_count": int(r.size),
        "axis_near_count": len(AXIS_NEAR_RADII),
        "designed_seam_count": len(SEAM_RADII),
        "observed_seam_count": seam_hits,
        "minimum_radius": float(np.min(r)),
        "source_binding_relative_error": source_binding_relative,
        "phase_embedding_max_abs_error": phase_error,
        "ladder": ladder,
        "refinement_ratios": [float(x) for x in refinement],
        "wrong_sign_vr_mutation_relative_rms": _relative_rms(dr_wrong_sign - dr_fd, dr_fd),
        "d_r_plus_0p02_mutation_relative_rms": _relative_rms(dr_wrong_exp - dr_fd, dr_fd),
    }


def generate_report() -> dict[str, Any]:
    rng = np.random.default_rng(SEED)
    cases = [_case_report(rng, h) for h in H_CASES]
    worst_phase = max(row["phase_embedding_max_abs_error"] for row in cases)
    worst_binding = max(row["source_binding_relative_error"] for row in cases)
    worst_finest = max(row["ladder"][-1]["primary_relative_rms"] for row in cases)
    worst_refinement = min(min(row["refinement_ratios"]) for row in cases)
    weakest_sign_mutation = min(row["wrong_sign_vr_mutation_relative_rms"] for row in cases)
    weakest_exp_mutation = min(row["d_r_plus_0p02_mutation_relative_rms"] for row in cases)
    checks = {
        "source_binding_matches_independent_formula": worst_binding <= SOURCE_BINDING_RELATIVE_GUARD,
        "black_box_phase_matches_independent_formula": worst_phase <= PHASE_EMBEDDING_MAX_GUARD,
        "finest_fd6_derivative_error": worst_finest <= FINE_RELATIVE_RMS_GUARD,
        "three_level_refinement": worst_refinement >= MIN_REFINEMENT_RATIO_GUARD,
        "wrong_sign_source_vector_detected": weakest_sign_mutation >= SIGN_MUTATION_RELATIVE_RMS_FLOOR,
        "radial_exponent_mutation_detected": weakest_exp_mutation >= EXPONENT_MUTATION_RELATIVE_RMS_FLOOR,
        "designed_mod1_seams_exercised": all(row["observed_seam_count"] >= len(SEAM_RADII) for row in cases),
        "axis_near_stencil_stays_positive": all(row["minimum_radius"] > 3.0*max(STEPS) for row in cases),
        "source_binding_tag_present": all(row["binding"] == "corrected_source_displayed_constants" for row in cases),
    }
    passed = bool(all(checks.values()))
    return {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "base_pr": BASE_PR,
        "base_head": BASE_HEAD,
        "seed": SEED,
        "independent_operator": "black-box torus_phase values + independent exact source constants + centered FD6",
        "sample_contract": {"h_cases": list(H_CASES), "resolution_ladder": list(STEPS),
                            "random_off_grid_points_per_case": RANDOM_POINT_COUNT,
                            "axis_near_radii": list(AXIS_NEAR_RADII), "designed_mod1_seam_radii": list(SEAM_RADII),
                            "uses_agent2_differentiate": False, "uses_agent2_fd4": False,
                            "uses_training_tensor_or_loss": False, "uses_pressure_or_forcing_fit": False},
        "frozen_guards": {"phase_embedding_max_abs": PHASE_EMBEDDING_MAX_GUARD,
                          "source_binding_relative": SOURCE_BINDING_RELATIVE_GUARD,
                          "finest_relative_rms": FINE_RELATIVE_RMS_GUARD,
                          "minimum_refinement_ratio": MIN_REFINEMENT_RATIO_GUARD,
                          "wrong_sign_mutation_relative_rms_floor": SIGN_MUTATION_RELATIVE_RMS_FLOOR,
                          "exponent_mutation_relative_rms_floor": EXPONENT_MUTATION_RELATIVE_RMS_FLOOR},
        "cases": cases,
        "summary": {"worst_phase_embedding_max_abs_error": worst_phase,
                    "worst_source_binding_relative_error": worst_binding,
                    "worst_finest_primary_relative_rms": worst_finest,
                    "worst_refinement_ratio": worst_refinement,
                    "weakest_wrong_sign_mutation_relative_rms": weakest_sign_mutation,
                    "weakest_exponent_mutation_relative_rms": weakest_exp_mutation},
        "checks": checks,
        "local_structural_preflight_passed": passed,
        "formal_gates": {"normalized_momentum_max_and_L2": FORMAL_MOMENTUM_GATE,
                         "divergence_max_and_L2": FORMAL_DIVERGENCE_GATE,
                         "formal_full_domain_pde_gate_assessed": False,
                         "pde_validated": False},
        "truth_boundary": {"displayed_source_torus_binding_independently_audited": passed,
                           "actual_positive_order_background_bound": False,
                           "actual_auxiliary_torus_mode_family_bound": False,
                           "public_xyz_t_velocity_correction_materialized": False,
                           "full_leading_oscillatory_correction_composite_available": False,
                           "formal_full_domain_pde_gate_assessed": False,
                           "pde_validated": False},
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = generate_report()
    text = json.dumps(report, indent=2, sort_keys=True)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + "\n", encoding="utf-8")
    print(text)
    if not report["local_structural_preflight_passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
