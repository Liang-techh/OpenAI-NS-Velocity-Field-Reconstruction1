"""Fresh post-freeze path audit for the frozen Agent-7 PR #701 child.

Preregistered in issue #713 before any held-out trajectory result is evaluated.
The #701 coefficient is still calibrated exactly once on its historical 48-path
development protocol; this module then evaluates a disjoint 64-path protocol
with different radii, azimuths, z-bands, output count, tolerances and max step.
No held-out path or position is fed back into alpha.

This is trajectory/radial directional generalization evidence only.  It does
not promote the candidate, establish visual correspondence, or validate PDEs.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
from scipy.integrate import solve_ivp

import agent7_st052m_tip_odd_poloidal_screen as screen

TASK_ID = "CR003-ST052M-TIP-ODD-POLOIDAL-HOLDOUT-099"
PREREG_ISSUE = 713
SOURCE_PARENT_PR = 701
SOURCE_PARENT_HEAD = "e3ecc8654f2c9a2324265bbd01df596b4e503404"
SELECTION_GOVERNANCE_PR = 712
SOURCE_WITNESS_PR = 652
SOURCE_TEMPORAL_PR = 587
PUBLIC_SOURCE_URL = "https://openai.com/index/navier-stokes-solution/"

TIME_INTERVAL = (0.25, 0.75)
HOLDOUT_RADII = (0.50, 0.80, 1.10, 1.40)
HOLDOUT_BANDS = (("shoulder", 0.75), ("tip", 1.25))
HOLDOUT_ANGLES = 4
ANGLE_OFFSET = np.pi / 8.0
OUTPUT_SAMPLES = 37
CHECKPOINT_TIMES = (0.25, 0.375, 0.50, 0.625, 0.75)
CHECKPOINT_INDICES = (0, 9, 18, 27, 36)
SEGMENTS = tuple(zip(CHECKPOINT_INDICES[:-1], CHECKPOINT_INDICES[1:]))
SOLVER_METHOD = "DOP853"
SOLVER_RTOL = 5.0e-10
SOLVER_ATOL = 5.0e-12
SOLVER_MAX_STEP = 0.008
CONTROL_TURNS_MIN = 1.0e-14

TRUTH = {
    "canonical_velocity_changed": False,
    "saved_velocity_changed": False,
    "production_candidate_selected": False,
    "new_spatial_basis_added": False,
    "new_temporal_basis_added": False,
    "source_701_candidate_retuned": False,
    "heldout_evaluation_feedback_into_alpha": False,
    "parameter_grid_scan_performed": False,
    "optimization_performed": False,
    "pressure_or_force_changed": False,
    "held_out_pde_residual_evaluated": False,
    "parent_pde_receipt_transferred": False,
    "public_image_numeric_target_used": False,
    "pixel_similarity_objective_used": False,
    "visual_acceptance_threshold_defined": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "source_correspondence_verified": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _fresh_seed_table() -> tuple[np.ndarray, list[dict[str, Any]]]:
    rows: list[list[float]] = []
    metadata: list[dict[str, Any]] = []
    for radius in HOLDOUT_RADII:
        for angle_index in range(HOLDOUT_ANGLES):
            angle = ANGLE_OFFSET + 2.0 * np.pi * angle_index / HOLDOUT_ANGLES
            x = radius * np.cos(angle)
            y = radius * np.sin(angle)
            for band, z_abs in HOLDOUT_BANDS:
                for sign in (-1, 1):
                    rows.append([x, y, sign * z_abs])
                    metadata.append(
                        {
                            "radius": float(radius),
                            "angle_index": int(angle_index),
                            "angle": float(angle),
                            "band": band,
                            "sign": int(sign),
                            "z_abs": float(z_abs),
                        }
                    )
    seeds = np.asarray(rows, dtype=float)
    if seeds.shape != (64, 3):
        raise RuntimeError("fresh 64-path seed construction drifted")
    return seeds, metadata


def _development_seed_overlap_count(heldout: np.ndarray) -> int:
    dev, _ = screen.loc.parent.base._seed_table()
    heldout = np.asarray(heldout, dtype=float)
    dev = np.asarray(dev, dtype=float)
    return int(
        sum(
            np.any(np.all(heldout == dev_seed[None, :], axis=1))
            for dev_seed in dev
        )
    )


def _integrate_fresh(velocity_fn) -> tuple[np.ndarray, np.ndarray, list[dict[str, Any]]]:
    seeds, metadata = _fresh_seed_table()
    if _development_seed_overlap_count(seeds) != 0:
        raise RuntimeError("held-out seed overlaps historical development seed")

    t0, t1 = TIME_INTERVAL
    times = np.linspace(t0, t1, OUTPUT_SAMPLES)
    for index, expected in zip(CHECKPOINT_INDICES, CHECKPOINT_TIMES):
        if abs(float(times[index]) - expected) > 1.0e-15:
            raise RuntimeError("fresh checkpoint alignment drifted")
    n = len(seeds)

    def rhs(time, flat):
        pts = np.asarray(flat, dtype=float).reshape(n, 3)
        value = np.asarray(velocity_fn(pts, float(time)), dtype=float)
        if value.shape != pts.shape or not np.isfinite(value).all():
            raise ValueError("invalid velocity during fresh path integration")
        return value.reshape(-1)

    sol = solve_ivp(
        rhs,
        (t0, t1),
        seeds.reshape(-1),
        method=SOLVER_METHOD,
        t_eval=times,
        rtol=SOLVER_RTOL,
        atol=SOLVER_ATOL,
        max_step=SOLVER_MAX_STEP,
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    positions = np.asarray(sol.y.T, dtype=float).reshape(len(times), n, 3)
    if not np.isfinite(positions).all():
        raise RuntimeError("nonfinite fresh material-path positions")
    return times, positions, metadata


def _band_summary(
    velocity_fn,
    times: np.ndarray,
    positions: np.ndarray,
    metadata: list[dict[str, Any]],
    i0: int,
    i1: int,
    band: str,
) -> dict[str, Any]:
    idx = np.asarray([i for i, meta in enumerate(metadata) if meta["band"] == band], dtype=int)
    if len(idx) != 32:
        raise RuntimeError("fresh band path count drifted")
    radii = np.hypot(positions[:, :, 0], positions[:, :, 1])
    angles = np.unwrap(np.arctan2(positions[:, :, 1], positions[:, :, 0]), axis=0)
    delta_r = np.asarray(radii[i1, idx] - radii[i0, idx], dtype=float)
    start_ur = screen.loc._radial_velocity(
        velocity_fn, np.asarray(positions[i0, idx], dtype=float), float(times[i0])
    )
    turns = np.abs((angles[i1, idx] - angles[i0, idx]) / (2.0 * np.pi))
    return {
        "path_count": int(len(idx)),
        "inward_path_count": int(np.sum(delta_r < 0.0)),
        "mean_delta_r": float(np.mean(delta_r)),
        "min_delta_r": float(np.min(delta_r)),
        "max_delta_r": float(np.max(delta_r)),
        "mean_start_radial_velocity": float(np.mean(start_ur)),
        "mean_absolute_turns": float(np.mean(turns)),
        "maximum_absolute_turns": float(np.max(turns)),
    }


def _all_summary(positions: np.ndarray, i0: int, i1: int) -> dict[str, Any]:
    radii = np.hypot(positions[:, :, 0], positions[:, :, 1])
    angles = np.unwrap(np.arctan2(positions[:, :, 1], positions[:, :, 0]), axis=0)
    delta_r = np.asarray(radii[i1] - radii[i0], dtype=float)
    turns = np.abs((angles[i1] - angles[i0]) / (2.0 * np.pi))
    return {
        "path_count": int(len(delta_r)),
        "inward_path_count": int(np.sum(delta_r < 0.0)),
        "mean_delta_r": float(np.mean(delta_r)),
        "mean_absolute_turns": float(np.mean(turns)),
        "maximum_absolute_turns": float(np.max(turns)),
    }


def _field_records(velocity_fn) -> dict[str, Any]:
    times, positions, metadata = _integrate_fresh(velocity_fn)
    rows: dict[str, Any] = {}
    for i0, i1 in SEGMENTS:
        key = f"{times[i0]:.3f}-{times[i1]:.3f}"
        rows[key] = {
            "time_start": float(times[i0]),
            "time_end": float(times[i1]),
            "all": _all_summary(positions, i0, i1),
            "shoulder": _band_summary(velocity_fn, times, positions, metadata, i0, i1, "shoulder"),
            "tip": _band_summary(velocity_fn, times, positions, metadata, i0, i1, "tip"),
        }
    i0, i1 = CHECKPOINT_INDICES[0], CHECKPOINT_INDICES[-1]
    whole = {
        "time_start": float(times[i0]),
        "time_end": float(times[i1]),
        "all": _all_summary(positions, i0, i1),
        "shoulder": _band_summary(velocity_fn, times, positions, metadata, i0, i1, "shoulder"),
        "tip": _band_summary(velocity_fn, times, positions, metadata, i0, i1, "tip"),
    }
    return {"segments": rows, "whole_interval": whole}


def _relative(new: float, old: float) -> float:
    old = float(old)
    if abs(old) <= CONTROL_TURNS_MIN:
        raise ValueError("fresh control turns too small for relative normalization")
    return float(float(new) / old - 1.0)


def _finite(obj: Any) -> bool:
    if isinstance(obj, dict):
        return all(_finite(v) for v in obj.values())
    if isinstance(obj, (list, tuple)):
        return all(_finite(v) for v in obj)
    if isinstance(obj, (float, int, np.floating, np.integer)):
        return bool(np.isfinite(float(obj)))
    return True


def directional_decision(
    control: dict[str, Any],
    linear: dict[str, Any],
    child: dict[str, Any],
    candidate: dict[str, Any],
) -> dict[str, Any]:
    rows: dict[str, Any] = {}
    for key in control["segments"]:
        c = control["segments"][key]
        l = linear["segments"][key]
        h = child["segments"][key]
        q = candidate["segments"][key]

        control_turns = float(c["all"]["mean_absolute_turns"])
        control_turns_nonzero = bool(abs(control_turns) > CONTROL_TURNS_MIN)
        linear_rel = _relative(l["all"]["mean_absolute_turns"], control_turns)
        candidate_rel = _relative(q["all"]["mean_absolute_turns"], control_turns)
        turns_increment = float(candidate_rel - linear_rel)

        tip_delta_lower = bool(q["tip"]["mean_delta_r"] < h["tip"]["mean_delta_r"])
        tip_start_ur_lower = bool(
            q["tip"]["mean_start_radial_velocity"]
            < h["tip"]["mean_start_radial_velocity"]
        )
        tip_mean_inward = bool(q["tip"]["mean_delta_r"] < 0.0)
        shoulder_mean_inward = bool(q["shoulder"]["mean_delta_r"] < 0.0)
        shoulder_count_not_lower = bool(
            q["shoulder"]["inward_path_count"] >= h["shoulder"]["inward_path_count"]
        )
        finite = _finite({"control": c, "linear": l, "child": h, "candidate": q})
        preferred = bool(
            control_turns_nonzero
            and tip_delta_lower
            and tip_start_ur_lower
            and tip_mean_inward
            and shoulder_mean_inward
            and shoulder_count_not_lower
            and turns_increment > 0.0
            and finite
        )
        rows[key] = {
            "source_652_tip_mean_delta_r": h["tip"]["mean_delta_r"],
            "candidate_tip_mean_delta_r": q["tip"]["mean_delta_r"],
            "tip_mean_delta_r_strictly_lower": tip_delta_lower,
            "source_652_tip_start_mean_radial_velocity": h["tip"]["mean_start_radial_velocity"],
            "candidate_tip_start_mean_radial_velocity": q["tip"]["mean_start_radial_velocity"],
            "tip_start_mean_radial_velocity_strictly_lower": tip_start_ur_lower,
            "candidate_tip_mean_inward": tip_mean_inward,
            "candidate_shoulder_mean_delta_r": q["shoulder"]["mean_delta_r"],
            "candidate_shoulder_mean_inward": shoulder_mean_inward,
            "source_652_shoulder_inward_path_count": h["shoulder"]["inward_path_count"],
            "candidate_shoulder_inward_path_count": q["shoulder"]["inward_path_count"],
            "candidate_shoulder_inward_count_not_lower": shoulder_count_not_lower,
            "source_652_tip_inward_path_count": h["tip"]["inward_path_count"],
            "candidate_tip_inward_path_count": q["tip"]["inward_path_count"],
            "tip_inward_path_count_change": int(
                q["tip"]["inward_path_count"] - h["tip"]["inward_path_count"]
            ),
            "linear_turns_relative_to_control": linear_rel,
            "candidate_turns_relative_to_control": candidate_rel,
            "candidate_minus_linear_turns_fidelity_increment": turns_increment,
            "turns_fidelity_increment_positive": bool(turns_increment > 0.0),
            "control_turns_nonzero": control_turns_nonzero,
            "all_values_finite": finite,
            "segment_preferred": preferred,
        }

    preferred = bool(all(row["segment_preferred"] for row in rows.values()))
    return {
        "segments": rows,
        "fresh_tip_path_directional_preferred": preferred,
        "all_segments_preferred": preferred,
        "tip_count_gain_any_segment": bool(
            any(row["tip_inward_path_count_change"] > 0 for row in rows.values())
        ),
    }


def _build_frozen_candidate():
    control_fn, linear_fn, child_fn, epsilon, normalization = screen._build_fields()
    # Selection stage: exactly #701's historical development paths, never the holdout paths.
    dev_times, dev_positions, dev_metadata = screen._integrate(child_fn)
    alpha, calibration = screen.derive_alpha_from_frozen_midpoint(
        child_fn, dev_times, dev_positions, dev_metadata
    )
    candidate_fn = lambda p, t: screen.corrected_velocity(child_fn, p, t, alpha)
    return control_fn, linear_fn, child_fn, candidate_fn, alpha, calibration, epsilon, normalization


def run(out: Path) -> dict[str, Any]:
    (
        control_fn,
        linear_fn,
        child_fn,
        candidate_fn,
        alpha,
        calibration,
        epsilon,
        normalization,
    ) = _build_frozen_candidate()

    fresh_seeds, _ = _fresh_seed_table()
    overlap = _development_seed_overlap_count(fresh_seeds)
    if overlap != 0:
        raise RuntimeError("fresh evaluation is not seed-disjoint")

    control = _field_records(control_fn)
    linear = _field_records(linear_fn)
    child = _field_records(child_fn)
    candidate = _field_records(candidate_fn)
    decision = directional_decision(control, linear, child, candidate)

    routing = (
        "the frozen #701 tip-local odd-poloidal channel generalizes to a disjoint post-freeze trajectory protocol: carry this exact channel and coefficient to a separate fixed-render/save-load comparison before any further basis growth"
        if decision["fresh_tip_path_directional_preferred"]
        else "the frozen #701 tip-local odd-poloidal channel does not establish fresh post-freeze inward-tip trajectory generalization; do not retune alpha and do not promote this channel from this audit"
    )

    report = {
        "task_id": TASK_ID,
        "prereg_issue": PREREG_ISSUE,
        "source_parent": {"pr": SOURCE_PARENT_PR, "head": SOURCE_PARENT_HEAD},
        "selection_governance_pr": SELECTION_GOVERNANCE_PR,
        "candidate_lineage": {
            "source_visualization_child_pr": SOURCE_WITNESS_PR,
            "linear_temporal_pr": SOURCE_TEMPORAL_PR,
        },
        "public_source": {"url": PUBLIC_SOURCE_URL, "qualitative_only": True},
        "selection_stage": {
            "alpha": alpha,
            "calibration": calibration,
            "development_path_count": 48,
            "development_radii": list(screen.loc.parent.base.SEED_RADII),
            "development_bands": [list(v) for v in screen.loc.parent.base.SEED_BANDS],
            "development_output_samples": int(screen.loc.parent.base.OUTPUT_SAMPLES),
            "heldout_data_used_for_alpha": False,
            "source_701_dedicated_pass_assumed": False,
        },
        "fresh_protocol": {
            "path_count": 64,
            "radii": list(HOLDOUT_RADII),
            "bands": [list(v) for v in HOLDOUT_BANDS],
            "angle_count": HOLDOUT_ANGLES,
            "angle_offset_radians": float(ANGLE_OFFSET),
            "output_samples": OUTPUT_SAMPLES,
            "checkpoint_indices": list(CHECKPOINT_INDICES),
            "checkpoint_times": list(CHECKPOINT_TIMES),
            "solver": SOLVER_METHOD,
            "rtol": SOLVER_RTOL,
            "atol": SOLVER_ATOL,
            "max_step": SOLVER_MAX_STEP,
            "exact_seed_overlap_with_development": overlap,
            "evaluation_feedback_into_alpha": False,
        },
        "source_652_compact_swirl_epsilon": epsilon,
        "source_652_normalization": normalization,
        "fresh_records": {
            "control": control,
            "linear_587": linear,
            "source_652": child,
            "candidate_701": candidate,
        },
        "decision": decision,
        "fresh_tip_path_directional_preferred": decision[
            "fresh_tip_path_directional_preferred"
        ],
        "source_652_target_free_pareto_rejection_remains_binding": True,
        "source_701_original_screen_result_not_required_for_this_audit": True,
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
    print("alpha=", report["selection_stage"]["alpha"])
    print("fresh_preferred=", report["fresh_tip_path_directional_preferred"])
    print("decision=", report["decision"])
    print("routing=", report["routing"])


if __name__ == "__main__":
    main()
