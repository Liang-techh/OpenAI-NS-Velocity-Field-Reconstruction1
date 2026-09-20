"""Replay one preregistered same-dimension p=9,m=4 tip-envelope child.

Issue #739 freezes this increment before observing any nonlinear result.  The
parent is exact #652.  The spatial correction is the same one-dimensional
axisymmetric vector-potential channel screened in #701, but its axial bump is
replaced by #732's deterministic p=9,m=4 envelope.  Because #738 shows that
peak normalization does not preserve the scalar basis coordinate, the old
#701 alpha is *not* copied.  A new alpha is derived exactly once from the same
historical 48-path development protocol, using the new unit response.

This is expression-capacity / visualization-direction development evidence.
It consumes no #714 fresh-path result, changes no pressure or forcing, and does
not perform a held-out Navier--Stokes acceptance test.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np

import agent7_st052m_tip_odd_poloidal_screen as screen
import agent7_st052m_tip_return_flow_relocation as relocation

TASK_ID = "CR003-ST052M-TIP-RESHAPED-NONLINEAR-REPLAY-103"
PREREG_ISSUE = 739
SOURCE_PARENT_PR = 738
SOURCE_PARENT_HEAD = "8072b099b8886e7c8716347ef2a60066d5985206"
SOURCE_RESHAPE_PR = 732
SOURCE_OLD_CHILD_PR = 701
SOURCE_VISUAL_CHILD_PR = 652
SOURCE_TEMPORAL_PR = 587
RELATED_FRESH_HOLDOUT_PR = 714

SELECTED_P = 9
SELECTED_M = 4
TARGET_PROBE_RADII = (0.6, 0.9, 1.2)
TARGET_PROBE_ABS_Z = 1.55
TARGET_PROBE_ANGLES = (0.0, 0.5 * np.pi, np.pi, 1.5 * np.pi)
TARGET_PROBE_TIMES = (0.375, 0.50, 0.625, 0.75)
DIVERGENCE_STEP = 1.0e-5
DIVERGENCE_MAX = 1.0e-5
CONDITION_MAX = 25.0

TRUTH = {
    "canonical_velocity_changed": False,
    "saved_velocity_changed": False,
    "production_candidate_selected": False,
    "basis_dimension_changed": False,
    "same_dimension_envelope_replayed": True,
    "second_poloidal_basis_added": False,
    "new_temporal_basis_added": False,
    "historical_development_paths_used_for_coefficient": True,
    "fresh_714_path_data_used": False,
    "fresh_714_path_data_used_for_retuning": False,
    "parameter_grid_scan_performed": False,
    "continuous_envelope_fit_performed": False,
    "optimization_performed": False,
    "post_result_damping_or_retuning_performed": False,
    "pressure_or_force_changed": False,
    "held_out_pde_residual_evaluated": False,
    "parent_pde_receipt_transferred": False,
    "public_image_numeric_target_used": False,
    "pixel_similarity_objective_used": False,
    "targeted_visualization_direction_fingerprint_used": True,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "source_correspondence_verified": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _assert_source_lock() -> None:
    p, _, z_star = relocation.select_inner_exponent()
    if p != SELECTED_P or relocation.OUTER_EXPONENT != SELECTED_M:
        raise RuntimeError("#732 selected-envelope identity drifted")
    if z_star < relocation.TARGET_SIGN_CHANGE_ABS_Z:
        raise RuntimeError("#732 sign-location gate drifted")
    if screen.PREREG_ISSUE != 700 or screen.SOURCE_WITNESS_PR != SOURCE_VISUAL_CHILD_PR:
        raise RuntimeError("#701/#652 source lineage drifted")
    if RELATED_FRESH_HOLDOUT_PR != 714:
        raise RuntimeError("fresh holdout identity drifted")


def reshaped_unit_time_correction(points: np.ndarray, time: float) -> np.ndarray:
    pts = np.asarray(points, dtype=float)
    return screen.activation(float(time)) * relocation.reshaped_tip_correction(pts, SELECTED_P)


def reshaped_velocity(base_fn, points: np.ndarray, time: float, alpha: float) -> np.ndarray:
    pts = np.asarray(points, dtype=float)
    return np.asarray(base_fn(pts, float(time)), dtype=float) + float(alpha) * reshaped_unit_time_correction(
        pts, float(time)
    )


def derive_reshaped_alpha(
    child_fn,
    child_times: np.ndarray,
    child_positions: np.ndarray,
    metadata: list[dict[str, Any]],
) -> tuple[float, dict[str, Any]]:
    """Apply issue #739's one-shot development-coordinate rule."""
    if abs(float(child_times[screen.MID_INDEX]) - screen.MID_TIME) > 1.0e-15:
        raise RuntimeError("historical midpoint index drifted")
    tip_idx = np.asarray(
        [i for i, meta in enumerate(metadata) if str(meta["band"]) == "tip"], dtype=int
    )
    if len(tip_idx) != 24:
        raise RuntimeError("historical tip path count drifted")
    pts = np.asarray(child_positions[screen.MID_INDEX, tip_idx], dtype=float)
    base_ur = screen.loc._radial_velocity(child_fn, pts, screen.MID_TIME)
    unit_ur = screen.loc._radial_velocity(reshaped_unit_time_correction, pts, screen.MID_TIME)
    mean_base = float(np.mean(base_ur))
    mean_unit = float(np.mean(unit_ur))
    if not (np.isfinite(mean_base) and np.isfinite(mean_unit)):
        raise RuntimeError("nonfinite reshaped midpoint coordinate calibration")
    if abs(mean_unit) <= np.finfo(float).tiny:
        raise RuntimeError("reshaped unit response is degenerate on frozen development tips")
    alpha = float(-mean_base / mean_unit)
    if not np.isfinite(alpha):
        raise RuntimeError("nonfinite reshaped coefficient")
    return alpha, {
        "time": float(screen.MID_TIME),
        "tip_path_count": int(len(tip_idx)),
        "mean_source_652_radial_velocity": mean_base,
        "mean_reshaped_unit_radial_velocity": mean_unit,
        "derived_alpha_94": alpha,
        "fixed_position_mean_after_linear_cancellation": float(np.mean(base_ur + alpha * unit_ur)),
        "historical_development_protocol_used": True,
        "fresh_714_data_used": False,
        "coefficient_grid_scan_performed": False,
        "optimizer_used": False,
        "post_result_damping_allowed": False,
    }


def _target_probe_points() -> np.ndarray:
    rows = []
    for r in TARGET_PROBE_RADII:
        for z in (-TARGET_PROBE_ABS_Z, TARGET_PROBE_ABS_Z):
            for angle in TARGET_PROBE_ANGLES:
                rows.append([r * np.cos(angle), r * np.sin(angle), z])
    pts = np.asarray(rows, dtype=float)
    if pts.shape != (24, 3):
        raise RuntimeError("targeted outer-tip probe count drifted")
    return pts


def targeted_outer_tip_fingerprint(child_fn, old_candidate_fn, new_candidate_fn) -> dict[str, Any]:
    pts = _target_probe_points()
    rows: dict[str, Any] = {}
    all_pass = True
    for time in TARGET_PROBE_TIMES:
        source_ur = screen.loc._radial_velocity(child_fn, pts, time)
        old_ur = screen.loc._radial_velocity(old_candidate_fn, pts, time)
        new_ur = screen.loc._radial_velocity(new_candidate_fn, pts, time)
        source_mean = float(np.mean(source_ur))
        old_mean = float(np.mean(old_ur))
        new_mean = float(np.mean(new_ur))
        gate = bool(new_mean < source_mean and new_mean < old_mean and new_mean < 0.0)
        all_pass = all_pass and gate
        rows[f"{time:.3f}"] = {
            "source_652_mean_radial_velocity": source_mean,
            "old_701_mean_radial_velocity": old_mean,
            "reshaped_mean_radial_velocity": new_mean,
            "reshaped_minus_source_652": float(new_mean - source_mean),
            "reshaped_minus_old_701": float(new_mean - old_mean),
            "source_652_inward_count": int(np.sum(source_ur < 0.0)),
            "old_701_inward_count": int(np.sum(old_ur < 0.0)),
            "reshaped_inward_count": int(np.sum(new_ur < 0.0)),
            "strictly_lower_than_both_and_inward": gate,
        }
    return {
        "probe_count": int(len(pts)),
        "radii": list(TARGET_PROBE_RADII),
        "abs_z": float(TARGET_PROBE_ABS_Z),
        "angles": [float(v) for v in TARGET_PROBE_ANGLES],
        "times": [float(v) for v in TARGET_PROBE_TIMES],
        "by_time": rows,
        "all_times_strictly_lower_than_both_and_inward": bool(all_pass),
        "scope": "targeted development fingerprint; not blind/global visual validation",
    }


def response_diagnostic(linear_fn, child_fn) -> dict[str, Any]:
    rng = np.random.default_rng(screen.RESPONSE_SEED)
    existing_blocks = []
    new_blocks = []
    for time in screen.RESPONSE_TIMES:
        radius = rng.uniform(0.2, 1.55, size=32)
        angle = rng.uniform(0.0, 2.0 * np.pi, size=32)
        z = rng.uniform(-1.75, 1.75, size=32)
        pts = np.column_stack((radius * np.cos(angle), radius * np.sin(angle), z))
        existing_blocks.append(
            (np.asarray(child_fn(pts, time), float) - np.asarray(linear_fn(pts, time), float)).ravel()
        )
        new_blocks.append(np.asarray(reshaped_unit_time_correction(pts, time), float).ravel())
    a = np.concatenate(existing_blocks)
    b = np.concatenate(new_blocks)
    na, nb = float(np.linalg.norm(a)), float(np.linalg.norm(b))
    if na <= np.finfo(float).tiny or nb <= np.finfo(float).tiny:
        raise RuntimeError("degenerate reshaped response diagnostic")
    matrix = np.column_stack((a / na, b / nb))
    singular = np.linalg.svd(matrix, compute_uv=False)
    condition = float(singular[0] / singular[-1])
    return {
        "parameter_count_increment_relative_to_652": 1,
        "basis_dimension_same_as_701": True,
        "rank": int(np.linalg.matrix_rank(matrix, tol=1.0e-10)),
        "singular_values": singular.tolist(),
        "condition_number": condition,
        "column_cosine": float(np.dot(matrix[:, 0], matrix[:, 1])),
        "existing_response_norm": na,
        "reshaped_unit_response_norm": nb,
        "passes_rank_condition_guard": bool(
            np.linalg.matrix_rank(matrix, tol=1.0e-10) == 2 and condition <= CONDITION_MAX
        ),
    }


def correction_rms_fraction(child_fn, alpha: float) -> dict[str, float]:
    axis, velocity = screen.vp.base.morph.prior.sample_velocity_grid(
        child_fn, screen.MID_TIME, screen.GRID_RESOLUTION
    )
    pts = screen._grid_points(axis)
    perturb = float(alpha) * reshaped_unit_time_correction(pts, screen.MID_TIME)
    base_rms = float(np.sqrt(np.mean(np.sum(np.square(velocity.reshape(-1, 3)), axis=1))))
    perturb_rms = float(np.sqrt(np.mean(np.sum(np.square(perturb), axis=1))))
    if base_rms <= np.finfo(float).tiny:
        raise RuntimeError("degenerate source-652 midpoint RMS")
    return {
        "base_velocity_rms": base_rms,
        "reshaped_correction_rms": perturb_rms,
        "correction_to_base_rms_fraction": float(perturb_rms / base_rms),
    }


def structure_preflight(child_fn, candidate_fn) -> dict[str, Any]:
    rng = np.random.default_rng(screen.STRUCTURE_SEED)
    start_pts = rng.uniform(-1.9, 1.9, size=(96, 3))
    start_identity = float(
        np.max(
            np.abs(
                np.asarray(candidate_fn(start_pts, 0.25), float)
                - np.asarray(child_fn(start_pts, 0.25), float)
            )
        )
    )

    shoulder = np.asarray(
        [[r, 0.0, z] for z in (-0.85, 0.85) for r in (0.6, 0.9, 1.2)], dtype=float
    )
    shoulder_correction = float(
        np.max(np.abs(relocation.reshaped_tip_correction(shoulder, SELECTED_P)))
    )
    outside = np.asarray(
        [
            [1.61, 0.0, 1.25],
            [-1.61, 0.0, -1.25],
            [0.9, 0.0, 0.99],
            [0.9, 0.0, -0.99],
            [0.9, 0.0, 1.81],
            [0.9, 0.0, -1.81],
        ],
        dtype=float,
    )
    outside_correction = float(
        np.max(np.abs(relocation.reshaped_tip_correction(outside, SELECTED_P)))
    )

    random_probes = rng.uniform(-1.45, 1.45, size=(40, 3))
    active_probes = _target_probe_points()
    probes = np.vstack((random_probes, active_probes))
    div_max = 0.0
    for time in TARGET_PROBE_TIMES:
        div = np.zeros(len(probes), dtype=float)
        for axis in range(3):
            shift = np.zeros(3, dtype=float)
            shift[axis] = DIVERGENCE_STEP
            up = np.asarray(candidate_fn(probes + shift, time), float)
            um = np.asarray(candidate_fn(probes - shift, time), float)
            div += (up[:, axis] - um[:, axis]) / (2.0 * DIVERGENCE_STEP)
        div_max = max(div_max, float(np.max(np.abs(div))))

    finite = bool(
        np.isfinite(start_identity)
        and np.isfinite(shoulder_correction)
        and np.isfinite(outside_correction)
        and np.isfinite(div_max)
    )
    passes = bool(
        finite
        and start_identity <= 1.0e-12
        and shoulder_correction == 0.0
        and outside_correction == 0.0
        and div_max <= DIVERGENCE_MAX
    )
    return {
        "start_identity_max_abs": start_identity,
        "shoulder_seed_correction_max_abs": shoulder_correction,
        "outside_correction_max_abs": outside_correction,
        "candidate_cartesian_fd_divergence_max": div_max,
        "divergence_probe_count": int(len(probes)),
        "finite": finite,
        "passes": passes,
    }


def run(out: Path) -> dict[str, Any]:
    _assert_source_lock()
    control_fn, linear_fn, child_fn, source_epsilon, source_normalization = screen._build_fields()

    child_times, child_positions, metadata = screen._integrate(child_fn)
    alpha_94, new_calibration = derive_reshaped_alpha(
        child_fn, child_times, child_positions, metadata
    )
    old_alpha, old_calibration = screen.derive_alpha_from_frozen_midpoint(
        child_fn, child_times, child_positions, metadata
    )

    old_candidate_fn = lambda p, t: screen.corrected_velocity(child_fn, p, t, old_alpha)
    candidate_fn = lambda p, t: reshaped_velocity(child_fn, p, t, alpha_94)

    control_times, control_positions, control_metadata = screen._integrate(control_fn)
    linear_times, linear_positions, linear_metadata = screen._integrate(linear_fn)
    candidate_times, candidate_positions, candidate_metadata = screen._integrate(candidate_fn)
    for name, times, meta in (
        ("control", control_times, control_metadata),
        ("linear", linear_times, linear_metadata),
        ("candidate", candidate_times, candidate_metadata),
    ):
        if not np.array_equal(times, child_times) or meta != metadata:
            raise RuntimeError(f"{name} historical trajectory protocol drifted")

    grouped_child = screen._grouped_records_from_positions(
        child_fn, child_times, child_positions, metadata
    )
    grouped_candidate = screen._grouped_records_from_positions(
        candidate_fn, candidate_times, candidate_positions, metadata
    )
    radial = screen._radial_band_decision(grouped_child, grouped_candidate)

    path_control = screen._path_records_from_positions(control_times, control_positions)
    path_linear = screen._path_records_from_positions(linear_times, linear_positions)
    path_child = screen._path_records_from_positions(child_times, child_positions)
    path_candidate = screen._path_records_from_positions(candidate_times, candidate_positions)
    turns = screen._turns_retention(path_control, path_linear, path_child, path_candidate)

    aspect = screen._aspect_decision(control_fn, linear_fn, candidate_fn)
    response = response_diagnostic(linear_fn, child_fn)
    structure = structure_preflight(child_fn, candidate_fn)
    targeted = targeted_outer_tip_fingerprint(child_fn, old_candidate_fn, candidate_fn)
    rms = correction_rms_fraction(child_fn, alpha_94)

    justified = bool(
        structure["passes"]
        and radial["all_tip_mean_delta_r_strictly_lower"]
        and radial["all_tip_start_mean_radial_velocity_strictly_lower"]
        and radial["shoulder_24_of_24_inward_all_segments"]
        and turns["all_segments_positive"]
        and aspect["all_post_start_positive"]
        and response["passes_rank_condition_guard"]
        and targeted["all_times_strictly_lower_than_both_and_inward"]
    )

    p, q_star, z_star = relocation.select_inner_exponent()
    report = {
        "task_id": TASK_ID,
        "prereg_issue": PREREG_ISSUE,
        "source_parent": {"pr": SOURCE_PARENT_PR, "head": SOURCE_PARENT_HEAD},
        "source_reshape_pr": SOURCE_RESHAPE_PR,
        "source_old_child_pr": SOURCE_OLD_CHILD_PR,
        "related_fresh_holdout_pr": RELATED_FRESH_HOLDOUT_PR,
        "candidate_lineage": {
            "baseline_temporal_pr": SOURCE_TEMPORAL_PR,
            "source_visualization_child_pr": SOURCE_VISUAL_CHILD_PR,
        },
        "frozen_reshaped_basis": {
            "representation": "curl(A), A=(-y*f,x*f,0), f=R(r^2)*z*B_9_4(z^2)",
            "inner_exponent_p": int(p),
            "outer_exponent_m": int(SELECTED_M),
            "physical_q_sign_change": float(q_star),
            "physical_abs_z_sign_change": float(z_star),
            "physical_radial_support_max": float(screen.R_SUPPORT_MAX),
            "physical_axial_support_abs": [float(screen.Z_SUPPORT_INNER), float(screen.Z_SUPPORT_OUTER)],
            "time_activation": "g(t)=2*(t-.25)",
            "basis_dimension_relative_to_701": "unchanged_one_tip_poloidal_channel",
            "coefficient_semantics": "rederived on new unit basis from historical #701 development protocol",
        },
        "source_652_compact_swirl_epsilon": source_epsilon,
        "source_652_normalization": source_normalization,
        "coefficient_calibration": new_calibration,
        "old_701_coefficient_calibration": old_calibration,
        "coefficient_comparison": {
            "old_701_alpha": float(old_alpha),
            "reshaped_alpha_94": float(alpha_94),
            "reshaped_to_old_alpha_ratio": float(alpha_94 / old_alpha) if old_alpha != 0.0 else None,
            "old_alpha_reused_as_same_coordinate": False,
        },
        "midpoint_correction_rms": rms,
        "response_diagnostic": response,
        "structure_preflight": structure,
        "targeted_outer_tip_fingerprint": targeted,
        "historical_radial_band_replay": {
            "source_652": grouped_child,
            "candidate": grouped_candidate,
            "decision": radial,
        },
        "historical_turns_replay": {
            "control": path_control,
            "linear_587": path_linear,
            "source_652": path_child,
            "candidate": path_candidate,
            "decision": turns,
        },
        "axial_stretch_replay": aspect,
        "same_dimension_reshaped_child_visualization_direction_screen_passed": justified,
        "decision": {
            "classification": (
                "same_dimension_reshaped_nonlinear_screen_pass"
                if justified
                else "same_dimension_reshaped_nonlinear_screen_fail"
            ),
            "carry_to_fresh_path_or_fixed_render_comparison": justified,
            "second_poloidal_channel_justified_now": False,
            "retune_alpha_or_p_m_from_this_result_allowed": False,
            "velocity_formula_evaluated_this_increment": True,
            "canonical_or_saved_velocity_promoted": False,
            "closer_visualization_delivery_supported": justified,
            "closer_visualization_delivery_scope": (
                "targeted development fingerprint plus frozen morphology/path guards only; not global visual correspondence"
            ),
        },
        **TRUTH,
    }
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=Path("outputs/agent7-st052m-tip-reshaped-nonlinear-replay/report.json"))
    args = parser.parse_args()
    report = run(args.out)
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
