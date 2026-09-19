"""Audit segment-wise persistence of exact #652 inward-spiral advantage over #587.

Issue #682 freezes this one-increment audit before any segment-wise trajectory
result is evaluated.  It reuses the exact 48-seed DOP853 material-path protocol
from #587 and splits the already-frozen 33 output samples into four consecutive
time segments.  No basis, parameter, window, time law, pressure, forcing,
normalization or candidate identity is changed.

The public qualitative source is the same OpenAI 2026-09-08 page used by #667
and #674.  This module asks only whether the already-observed inward-spiral /
turns-fidelity advantage of exact #652 is persistent through the active time
interval rather than concentrated in one part of the trajectory.  It does not
establish pixel correspondence, PDE validity, paper exactness or exact
OpenAI-field identity.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

import agent7_st052m_public_source_temporal_persistence as temporal

# Reuse the exact stacked implementation ancestry rather than copying any field
# formula.  temporal.grid is #659; grid.witness is the exact #652 evaluator;
# witness.base is the exact #587 linear-temporal implementation.
grid = temporal.grid
witness = grid.witness
base = witness.base

TASK_ID = "CR003-ST052M-PUBLIC-SOURCE-INWARD-SPIRAL-PERSISTENCE-096"
PREREG_ISSUE = 682
PUBLIC_SOURCE_URL = "https://openai.com/index/navier-stokes-solution/"
PUBLIC_SOURCE_DATE = "2026-09-08"
SOURCE_DIRECTIONAL_PR = 667
SOURCE_TEMPORAL_PERSISTENCE_PR = 674
SOURCE_TEMPORAL_PERSISTENCE_HEAD = "553aaf86b2e341ec79a9c2664f63c120dcdb3054"
SOURCE_WITNESS_PR = 652
SOURCE_WITNESS_HEAD = "b0d98d0d5fdb3ca6a438fdf6ac0dbaaef7781faf"
SOURCE_TEMPORAL_PR = 587
SOURCE_TEMPORAL_HEAD = "0b93819095f6c8576a7571bdc2d2fbef4154944d"

CHECKPOINT_TIMES = (0.25, 0.375, 0.50, 0.625, 0.75)
CHECKPOINT_INDICES = (0, 8, 16, 24, 32)
SEGMENTS = tuple(zip(CHECKPOINT_INDICES[:-1], CHECKPOINT_INDICES[1:]))
WHOLE_TURNS_INCREMENT_EXPECTED = 0.0007867402276163782
WHOLE_TURNS_INCREMENT_TOL = 2.0e-10
WHOLE_INWARD_COUNT_CHANGE_EXPECTED = 0
CONTROL_TURNS_MIN = 1.0e-14

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


def _integrate_positions(velocity_fn):
    """Integrate the exact #587 frozen 48-path protocol once for one field."""
    seeds, metadata = base._seed_table()
    t0, t1 = base.TIME_INTERVAL
    times = np.linspace(t0, t1, base.OUTPUT_SAMPLES)
    if len(times) != 33:
        raise RuntimeError("#587 output-sample count drifted")
    for index, expected in zip(CHECKPOINT_INDICES, CHECKPOINT_TIMES):
        if abs(float(times[index]) - expected) > 1.0e-15:
            raise RuntimeError("frozen 33-sample checkpoint alignment drifted")

    n = len(seeds)

    def rhs(time, flat):
        pts = np.asarray(flat, dtype=float).reshape(n, 3)
        value = np.asarray(velocity_fn(pts, float(time)), dtype=float)
        if value.shape != pts.shape or not np.isfinite(value).all():
            raise ValueError("invalid velocity during path integration")
        return value.reshape(-1)

    sol = solve_ivp(
        rhs,
        (t0, t1),
        seeds.reshape(-1),
        method=base.SOLVER_METHOD,
        t_eval=times,
        rtol=base.SOLVER_RTOL,
        atol=base.SOLVER_ATOL,
        max_step=base.SOLVER_MAX_STEP,
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    positions = np.asarray(sol.y.T, dtype=float).reshape(len(times), n, 3)
    if not np.isfinite(positions).all():
        raise RuntimeError("nonfinite material-path positions")
    return times, positions, metadata


def _segment_summary(positions: np.ndarray, i0: int, i1: int) -> dict:
    radii = np.hypot(positions[:, :, 0], positions[:, :, 1])
    # Unwrap over the full continuous trajectory first, exactly as preregistered.
    angles = np.unwrap(np.arctan2(positions[:, :, 1], positions[:, :, 0]), axis=0)
    radius_change = radii[i1] - radii[i0]
    turns = np.abs((angles[i1] - angles[i0]) / (2.0 * np.pi))
    return {
        "mean_radius_change": float(np.mean(radius_change)),
        "contraction_magnitude": float(abs(np.mean(radius_change))),
        "inward_path_count": int(np.sum(radius_change < 0.0)),
        "mean_absolute_turns": float(np.mean(turns)),
        "maximum_absolute_turns": float(np.max(turns)),
    }


def path_records(velocity_fn) -> dict:
    times, positions, _ = _integrate_positions(velocity_fn)
    segments = {}
    for i0, i1 in SEGMENTS:
        key = f"{times[i0]:.3f}-{times[i1]:.3f}"
        segments[key] = _segment_summary(positions, i0, i1)
    whole = _segment_summary(positions, CHECKPOINT_INDICES[0], CHECKPOINT_INDICES[-1])
    return {
        "checkpoint_times": [float(times[i]) for i in CHECKPOINT_INDICES],
        "checkpoint_indices": list(CHECKPOINT_INDICES),
        "segments": segments,
        "whole_interval": whole,
    }


def _relative(new: float, old: float) -> float:
    if abs(float(old)) <= CONTROL_TURNS_MIN:
        raise ValueError("control segment turns are too small for frozen relative normalization")
    return float(float(new) / float(old) - 1.0)


def persistence_decision(control: dict, linear: dict, child: dict) -> dict:
    rows = {}
    for key in control["segments"]:
        c = control["segments"][key]
        l = linear["segments"][key]
        h = child["segments"][key]
        control_turns_nonzero = bool(abs(float(c["mean_absolute_turns"])) > CONTROL_TURNS_MIN)
        linear_turns_relative = _relative(l["mean_absolute_turns"], c["mean_absolute_turns"])
        child_turns_relative = _relative(h["mean_absolute_turns"], c["mean_absolute_turns"])
        turns_increment = float(child_turns_relative - linear_turns_relative)
        linear_inward = bool(float(l["mean_radius_change"]) < 0.0)
        child_inward = bool(float(h["mean_radius_change"]) < 0.0)
        inward_count_change = int(h["inward_path_count"] - l["inward_path_count"])
        finite = _all_finite({"control": c, "linear": l, "child": h})
        segment_preferred = bool(
            control_turns_nonzero
            and linear_inward
            and child_inward
            and turns_increment > 0.0
            and inward_count_change >= 0
            and finite
        )
        rows[key] = {
            "control": c,
            "linear": l,
            "child": h,
            "linear_turns_relative_to_control": linear_turns_relative,
            "child_turns_relative_to_control": child_turns_relative,
            "child_minus_linear_turns_fidelity_increment": turns_increment,
            "linear_mean_radial_change_negative_guard": linear_inward,
            "child_mean_radial_change_negative_guard": child_inward,
            "inward_path_count_change_child_minus_linear": inward_count_change,
            "inward_path_count_not_lower_guard": bool(inward_count_change >= 0),
            "control_turns_nonzero_guard": control_turns_nonzero,
            "all_segment_values_finite": finite,
            "segment_preferred": segment_preferred,
        }

    whole_c = control["whole_interval"]
    whole_l = linear["whole_interval"]
    whole_h = child["whole_interval"]
    whole_linear_rel = _relative(whole_l["mean_absolute_turns"], whole_c["mean_absolute_turns"])
    whole_child_rel = _relative(whole_h["mean_absolute_turns"], whole_c["mean_absolute_turns"])
    whole_turns_increment = float(whole_child_rel - whole_linear_rel)
    whole_inward_count_change = int(whole_h["inward_path_count"] - whole_l["inward_path_count"])
    whole_consistency = bool(
        abs(whole_turns_increment - WHOLE_TURNS_INCREMENT_EXPECTED)
        <= WHOLE_TURNS_INCREMENT_TOL
        and whole_inward_count_change == WHOLE_INWARD_COUNT_CHANGE_EXPECTED
    )

    preferred = bool(all(row["segment_preferred"] for row in rows.values()) and whole_consistency)
    return {
        "segments": rows,
        "whole_interval_consistency": {
            "linear_turns_relative_to_control": whole_linear_rel,
            "child_turns_relative_to_control": whole_child_rel,
            "child_minus_linear_turns_fidelity_increment": whole_turns_increment,
            "expected_from_667": WHOLE_TURNS_INCREMENT_EXPECTED,
            "absolute_tolerance": WHOLE_TURNS_INCREMENT_TOL,
            "inward_path_count_change_child_minus_linear": whole_inward_count_change,
            "expected_inward_path_count_change_from_667": WHOLE_INWARD_COUNT_CHANGE_EXPECTED,
            "guard": whole_consistency,
        },
        "all_segments_preferred": bool(all(row["segment_preferred"] for row in rows.values())),
        "public_source_temporal_inward_spiral_preferred": preferred,
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
        "control": path_records(control_fn),
        "linear": path_records(linear_fn),
        "child": path_records(child_fn),
    }
    decision = persistence_decision(records["control"], records["linear"], records["child"])
    preferred = bool(decision["public_source_temporal_inward_spiral_preferred"])

    routing = (
        "exact #652 preserves its frozen inward-spiral / turns-fidelity advantage over #587 in every sampled temporal segment; keep #652 ahead for the public-source visualization lane and stop generic basis growth in Agent 7"
        if preferred
        else "the exact #652 whole-interval inward-spiral advantage is not persistent in every frozen temporal segment; retain the segment-level failure without retuning and do not add a basis from this audit alone"
    )

    report = {
        "task_id": TASK_ID,
        "prereg_issue": PREREG_ISSUE,
        "public_source": {
            "url": PUBLIC_SOURCE_URL,
            "date": PUBLIC_SOURCE_DATE,
            "qualitative_only": True,
            "statement_scope": "inward spiraling of the vortex / displayed trajectories only",
        },
        "candidate_pair": {
            "baseline": {"pr": SOURCE_TEMPORAL_PR, "head": SOURCE_TEMPORAL_HEAD},
            "challenger": {"pr": SOURCE_WITNESS_PR, "head": SOURCE_WITNESS_HEAD},
            "source_directional_audit": {"pr": SOURCE_DIRECTIONAL_PR},
            "source_axial_time_persistence": {
                "pr": SOURCE_TEMPORAL_PERSISTENCE_PR,
                "head": SOURCE_TEMPORAL_PERSISTENCE_HEAD,
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
        "trajectory_protocol": {
            "seed_radii": list(base.SEED_RADII),
            "seed_bands": [[name, z] for name, z in base.SEED_BANDS],
            "seed_angles": base.SEED_ANGLES,
            "path_count": 48,
            "solver": base.SOLVER_METHOD,
            "rtol": base.SOLVER_RTOL,
            "atol": base.SOLVER_ATOL,
            "max_step": base.SOLVER_MAX_STEP,
            "output_samples": base.OUTPUT_SAMPLES,
            "checkpoint_times": list(CHECKPOINT_TIMES),
            "checkpoint_indices": list(CHECKPOINT_INDICES),
            "angles_unwrapped_over_full_path_before_segment_differences": True,
        },
        "normalization": normalization,
        "path_records": records,
        "temporal_inward_spiral_decision": decision,
        "public_source_temporal_inward_spiral_preferred": preferred,
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
    report = run(args.out)
    print(
        "public_source_temporal_inward_spiral_preferred=",
        report["public_source_temporal_inward_spiral_preferred"],
    )
    print("decision=", report["temporal_inward_spiral_decision"])
    print("routing=", report["routing"])


if __name__ == "__main__":
    main()
