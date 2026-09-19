"""Audit time persistence of the exact #652 axial-stretch advantage over #587.

Issue #673 freezes this one-increment audit before any new time-persistence
result is evaluated.  No basis, parameter, window, time law, pressure, forcing,
normalization or candidate identity is changed.

The public qualitative source is the same OpenAI 2026-09-08 page used by PR
#667.  This module asks only whether the already-selected enstrophy-aspect proxy
remains directionally favorable through the active time interval instead of
being a t=.50 coincidence.  It does not establish pixel correspondence, PDE
validity, paper exactness or exact OpenAI-field identity.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

import agent7_st052m_combined_witness_grid_robustness as grid

TASK_ID = "CR003-ST052M-PUBLIC-SOURCE-TEMPORAL-PERSISTENCE-095"
PREREG_ISSUE = 673
PUBLIC_SOURCE_URL = "https://openai.com/index/navier-stokes-solution/"
PUBLIC_SOURCE_DATE = "2026-09-08"
SOURCE_DIRECTIONAL_PR = 667
SOURCE_DIRECTIONAL_HEAD = "dd4046ac3e570e39e928acda2b00e46d590a4359"
SOURCE_WITNESS_PR = 652
SOURCE_WITNESS_HEAD = "b0d98d0d5fdb3ca6a438fdf6ac0dbaaef7781faf"
SOURCE_TEMPORAL_PR = 587
SOURCE_TEMPORAL_HEAD = "0b93819095f6c8576a7571bdc2d2fbef4154944d"
RESOLUTION = 41
TIMES = (0.25, 0.375, 0.50, 0.625, 0.75)
POST_START_TIMES = TIMES[1:]
START_IDENTITY_TOL = 1.0e-12

TRUTH = {
    "canonical_velocity_changed": False,
    "saved_velocity_changed": False,
    "production_candidate_selected": False,
    "new_spatial_basis_added": False,
    "new_temporal_basis_added": False,
    "witness_retuned": False,
    "parameter_grid_scan_performed": False,
    "optimization_performed": False,
    "pressure_or_force_changed": False,
    "held_out_pde_residual_evaluated": False,
    "parent_pde_receipt_transferred": False,
    "public_image_numeric_target_used": False,
    "pixel_similarity_objective_used": False,
    "hidden_openai_parameters_inferred": False,
    "visual_acceptance_threshold_defined": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "source_correspondence_verified": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _all_finite(obj) -> bool:
    if isinstance(obj, dict):
        return all(_all_finite(v) for v in obj.values())
    if isinstance(obj, (list, tuple)):
        return all(_all_finite(v) for v in obj)
    if isinstance(obj, (float, int, np.floating, np.integer)):
        return bool(np.isfinite(float(obj)))
    return True


def time_record(control_fn, linear_fn, child_fn, time: float) -> dict:
    control_abs = grid.morphology_metrics(control_fn, time, RESOLUTION)
    linear_abs = grid.morphology_metrics(linear_fn, time, RESOLUTION)
    child_abs = grid.morphology_metrics(child_fn, time, RESOLUTION)
    linear_rel = grid.compare_metrics(control_abs, linear_abs)
    child_rel = grid.compare_metrics(control_abs, child_abs)
    return {
        "time": float(time),
        "absolute": {
            "control": control_abs,
            "linear": linear_abs,
            "child": child_abs,
        },
        "relative_to_control": {
            "linear": linear_rel,
            "child": child_rel,
        },
        "child_minus_linear": {
            "full_aspect_ratio": float(
                child_rel["full_aspect_ratio"] - linear_rel["full_aspect_ratio"]
            ),
            "smooth_tip_radial_rms": float(
                child_rel["smooth_tip_radial_rms"]
                - linear_rel["smooth_tip_radial_rms"]
            ),
            "central_radial_rms": float(
                child_rel["central_radial_rms"]
                - linear_rel["central_radial_rms"]
            ),
            "full_enstrophy": float(
                child_rel["full_enstrophy"] - linear_rel["full_enstrophy"]
            ),
        },
    }


def temporal_decision(records: dict[str, dict]) -> dict:
    def rec(t: float) -> dict:
        return records[f"{t:.3f}"]

    start = rec(0.25)
    linear_start = float(
        start["relative_to_control"]["linear"]["full_aspect_ratio"]
    )
    child_start = float(
        start["relative_to_control"]["child"]["full_aspect_ratio"]
    )
    start_identity = bool(
        abs(linear_start) <= START_IDENTITY_TOL
        and abs(child_start) <= START_IDENTITY_TOL
    )

    aspect_increment_by_time = {
        f"{t:.3f}": float(rec(t)["child_minus_linear"]["full_aspect_ratio"])
        for t in POST_START_TIMES
    }
    challenger_advantage_all_times = bool(
        all(v > 0.0 for v in aspect_increment_by_time.values())
    )

    linear_sequence = [
        float(rec(t)["relative_to_control"]["linear"]["full_aspect_ratio"])
        for t in POST_START_TIMES
    ]
    child_sequence = [
        float(rec(t)["relative_to_control"]["child"]["full_aspect_ratio"])
        for t in POST_START_TIMES
    ]
    linear_increasing = bool(
        all(b > a for a, b in zip(linear_sequence, linear_sequence[1:]))
    )
    child_increasing = bool(
        all(b > a for a, b in zip(child_sequence, child_sequence[1:]))
    )
    finite = _all_finite(records)

    preferred = bool(
        start_identity
        and challenger_advantage_all_times
        and linear_increasing
        and child_increasing
        and finite
    )
    return {
        "resolution": RESOLUTION,
        "times": list(TIMES),
        "start_identity_tolerance": START_IDENTITY_TOL,
        "linear_start_aspect_relative_to_control": linear_start,
        "child_start_aspect_relative_to_control": child_start,
        "start_identity_guard": start_identity,
        "aspect_increment_child_vs_linear_by_time": aspect_increment_by_time,
        "challenger_aspect_advantage_all_post_start_times": challenger_advantage_all_times,
        "linear_aspect_relative_to_control_sequence": linear_sequence,
        "child_aspect_relative_to_control_sequence": child_sequence,
        "linear_increasing_elongation_guard": linear_increasing,
        "child_increasing_elongation_guard": child_increasing,
        "all_morphology_values_finite": finite,
        "public_source_temporal_axial_stretch_preferred": preferred,
    }


def run(out: Path) -> dict:
    field, raw = grid.witness.base.morph.prior.replay_st052.reconstruct()
    redistribution_scale, _, _ = grid.witness.base.morph.prior.base.child_scale(field, raw)
    energy_solve = grid.witness.base.morph.prior.comp.solve_energy_beta(
        field, raw, redistribution_scale
    )
    if not energy_solve["root_exists"]:
        raise RuntimeError("#559 frozen energy root disappeared")
    beta = float(energy_solve["beta"])
    if abs(beta - grid.witness.base.EXPECTED_BETA) > 2.0e-10:
        raise RuntimeError("#559 beta identity drift")

    control_fn = lambda p, t: grid.witness.base.control_velocity(
        field, raw, p, t, redistribution_scale
    )
    linear_fn = lambda p, t: grid.witness.base.ramp_velocity(
        field, raw, p, t, redistribution_scale, beta
    )
    epsilon, normalization = grid.witness.swirl.derive_epsilon(linear_fn)
    child_fn = lambda p, t: grid.witness.combined_velocity(
        field, raw, p, t, redistribution_scale, beta, epsilon
    )

    if abs(grid.witness.SWIRL_A - 0.065899695471146) > 0.0:
        raise RuntimeError("#652 frozen swirl coordinate drifted")
    if abs(grid.witness.SHOULDER_LAMBDA - (-0.02)) > 1.0e-15:
        raise RuntimeError("#652 frozen shoulder coordinate drifted")

    records = {
        f"{t:.3f}": time_record(control_fn, linear_fn, child_fn, t)
        for t in TIMES
    }
    decision = temporal_decision(records)
    preferred = bool(decision["public_source_temporal_axial_stretch_preferred"])

    routing = (
        "exact #652 preserves and improves the frozen axial-stretch proxy across the sampled time evolution; keep #652 ahead of #587 for the public-source visualization lane and route the exact child to independent delivery/render comparison without retuning or basis growth"
        if preferred
        else "the #652 axial-stretch advantage is not persistent under the frozen time audit; keep #667 as a midpoint directional result only and do not retune or add a basis from this failure"
    )

    report = {
        "task_id": TASK_ID,
        "prereg_issue": PREREG_ISSUE,
        "public_source": {
            "url": PUBLIC_SOURCE_URL,
            "date": PUBLIC_SOURCE_DATE,
            "qualitative_only": True,
            "statement_scope": "increasing elongation / axial stretching only",
        },
        "candidate_pair": {
            "baseline": {"pr": SOURCE_TEMPORAL_PR, "head": SOURCE_TEMPORAL_HEAD},
            "challenger": {"pr": SOURCE_WITNESS_PR, "head": SOURCE_WITNESS_HEAD},
            "source_directional_audit": {
                "pr": SOURCE_DIRECTIONAL_PR,
                "head": SOURCE_DIRECTIONAL_HEAD,
            },
        },
        "frozen_witness": {
            "swirl_a": grid.witness.SWIRL_A,
            "shoulder_lambda": grid.witness.SHOULDER_LAMBDA,
            "compact_swirl_epsilon": float(epsilon),
            "effective_compact_swirl_amplitude": float(grid.witness.SWIRL_A * epsilon),
            "parameter_grid_scan_performed": False,
            "optimization_performed": False,
        },
        "normalization": normalization,
        "morphology_by_time": records,
        "temporal_axial_stretch_decision": decision,
        "public_source_temporal_axial_stretch_preferred": preferred,
        "source_652_target_free_pareto_rejection_remains_binding": True,
        "basis_growth_justified_by_this_audit": False,
        "routing": routing,
        **TRUTH,
    }

    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    r = run(args.out)
    print(
        "public_source_temporal_axial_stretch_preferred=",
        r["public_source_temporal_axial_stretch_preferred"],
    )
    print("decision=", r["temporal_axial_stretch_decision"])
    print("routing=", r["routing"])


if __name__ == "__main__":
    main()
