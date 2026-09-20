"""Replay one preregistered outer-axial-reservoir tip child.

Issue #774 freezes this increment before observing any nonlinear result.  The
parent is exact visualization-oriented #652.  The spatial correction is the
single divergence-free vector-potential channel preflighted in #767: it stays
radially inward through the visible 1<|z|<1.8 tip lobe and places the required
return flow in the unused 1.8<|z|<2 outer axial reservoir.

The scalar coefficient is derived once from the historical #701 development
protocol.  No #714 fresh-path result is consumed, no coefficient scan or
post-result damping is allowed, and no pressure/forcing/PDE acceptance is
rebuilt here.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np

import agent7_st052m_outer_axial_return_reservoir_preflight as reservoir
import agent7_st052m_tip_odd_poloidal_screen as screen
import agent7_st052m_tip_reshaped_nonlinear_replay as p9

TASK_ID = "CR003-ST052M-OUTER-RESERVOIR-NONLINEAR-REPLAY-110"
PREREG_ISSUE = 774
SOURCE_PARENT_PR = 767
SOURCE_PARENT_HEAD = "0930a1b7eeb357dee3e83dd79a3c4fea293e9314"
SOURCE_VISUAL_CHILD_PR = 652
SOURCE_TEMPORAL_PR = 587
SOURCE_P9_PR = 740
RELATED_FRESH_HOLDOUT_PR = 714

VISIBLE_ABS_Z = (1.55, 1.75)
RESERVOIR_ABS_Z = (1.85, 1.95)
PROBE_RADII = (0.6, 0.9, 1.2)
PROBE_ANGLES = (0.0, 0.5 * np.pi, np.pi, 1.5 * np.pi)
PROBE_TIMES = (0.375, 0.50, 0.625, 0.75)
DIVERGENCE_STEP = 1.0e-5
DIVERGENCE_MAX = 1.0e-5
CONDITION_MAX = 25.0

TRUTH = {
    "canonical_velocity_changed": False,
    "saved_velocity_changed": False,
    "production_candidate_selected": False,
    "basis_dimension_changed": False,
    "single_tip_poloidal_channel_replayed": True,
    "second_poloidal_basis_added": False,
    "new_temporal_basis_added": False,
    "historical_development_paths_used_for_coefficient": True,
    "fresh_714_path_data_used": False,
    "fresh_714_path_data_used_for_retuning": False,
    "parameter_grid_scan_performed": False,
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
    reservoir._assert_source_lock()
    p9._assert_source_lock()
    if reservoir.TASK_ID != "CR003-ST052M-OUTER-AXIAL-RETURN-RESERVOIR-PREFLIGHT-108":
        raise RuntimeError("#767 reservoir task identity drifted")
    if reservoir.Z_VISIBLE_OUTER != 1.8 or reservoir.Z_GLOBAL_OUTER != 2.0:
        raise RuntimeError("#767 axial reservoir support drifted")
    if reservoir.R_SUPPORT_MAX != 1.6:
        raise RuntimeError("#767 radial support drifted")
    if RELATED_FRESH_HOLDOUT_PR != 714:
        raise RuntimeError("fresh holdout identity drifted")


def reservoir_unit_time_correction(points: np.ndarray, time: float) -> np.ndarray:
    pts = np.asarray(points, dtype=float)
    return screen.activation(float(time)) * reservoir.reservoir_correction(pts)


def reservoir_velocity(base_fn, points: np.ndarray, time: float, alpha: float) -> np.ndarray:
    pts = np.asarray(points, dtype=float)
    return np.asarray(base_fn(pts, float(time)), dtype=float) + float(alpha) * reservoir_unit_time_correction(
        pts, float(time)
    )


def derive_reservoir_alpha(
    child_fn,
    child_times: np.ndarray,
    child_positions: np.ndarray,
    metadata: list[dict[str, Any]],
) -> tuple[float, dict[str, Any]]:
    """Apply issue #774's one-shot historical development-coordinate rule."""
    if abs(float(child_times[screen.MID_INDEX]) - screen.MID_TIME) > 1.0e-15:
        raise RuntimeError("historical midpoint index drifted")
    tip_idx = np.asarray(
        [i for i, meta in enumerate(metadata) if str(meta["band"]) == "tip"], dtype=int
    )
    if len(tip_idx) != 24:
        raise RuntimeError("historical tip path count drifted")
    pts = np.asarray(child_positions[screen.MID_INDEX, tip_idx], dtype=float)
    base_ur = screen.loc._radial_velocity(child_fn, pts, screen.MID_TIME)
    unit_ur = screen.loc._radial_velocity(reservoir_unit_time_correction, pts, screen.MID_TIME)
    mean_base = float(np.mean(base_ur))
    mean_unit = float(np.mean(unit_ur))
    if not (np.isfinite(mean_base) and np.isfinite(mean_unit)):
        raise RuntimeError("nonfinite reservoir midpoint coordinate calibration")
    if abs(mean_unit) <= np.finfo(float).tiny:
        raise RuntimeError("reservoir unit response is degenerate on frozen development tips")
    alpha = float(-mean_base / mean_unit)
    if not np.isfinite(alpha):
        raise RuntimeError("nonfinite reservoir coefficient")
    return alpha, {
        "time": float(screen.MID_TIME),
        "tip_path_count": int(len(tip_idx)),
        "mean_source_652_radial_velocity": mean_base,
        "mean_reservoir_unit_radial_velocity": mean_unit,
        "derived_alpha_reservoir": alpha,
        "fixed_position_mean_after_linear_cancellation": float(np.mean(base_ur + alpha * unit_ur)),
        "historical_development_protocol_used": True,
        "fresh_714_data_used": False,
        "coefficient_grid_scan_performed": False,
        "optimizer_used": False,
        "post_result_damping_allowed": False,
    }


def _probe_points(abs_z: float) -> np.ndarray:
    rows = []
    for r in PROBE_RADII:
        for z in (-float(abs_z), float(abs_z)):
            for angle in PROBE_ANGLES:
                rows.append([r * np.cos(angle), r * np.sin(angle), z])
    pts = np.asarray(rows, dtype=float)
    if pts.shape != (24, 3):
        raise RuntimeError("frozen probe count drifted")
    return pts


def visible_tip_fingerprint(child_fn, p9_candidate_fn, candidate_fn) -> dict[str, Any]:
    by_abs_z: dict[str, Any] = {}
    all_pass = True
    for abs_z in VISIBLE_ABS_Z:
        pts = _probe_points(abs_z)
        rows: dict[str, Any] = {}
        band_pass = True
        for time in PROBE_TIMES:
            source_ur = screen.loc._radial_velocity(child_fn, pts, time)
            p9_ur = screen.loc._radial_velocity(p9_candidate_fn, pts, time)
            candidate_ur = screen.loc._radial_velocity(candidate_fn, pts, time)
            source_mean = float(np.mean(source_ur))
            p9_mean = float(np.mean(p9_ur))
            candidate_mean = float(np.mean(candidate_ur))
            gate = bool(candidate_mean < 0.0 and candidate_mean < source_mean and candidate_mean < p9_mean)
            band_pass = band_pass and gate
            rows[f"{time:.3f}"] = {
                "source_652_mean_radial_velocity": source_mean,
                "p9_740_mean_radial_velocity": p9_mean,
                "reservoir_child_mean_radial_velocity": candidate_mean,
                "reservoir_minus_source_652": float(candidate_mean - source_mean),
                "reservoir_minus_p9_740": float(candidate_mean - p9_mean),
                "source_652_inward_count": int(np.sum(source_ur < 0.0)),
                "p9_740_inward_count": int(np.sum(p9_ur < 0.0)),
                "reservoir_child_inward_count": int(np.sum(candidate_ur < 0.0)),
                "strictly_lower_than_both_and_inward": gate,
            }
        all_pass = all_pass and band_pass
        by_abs_z[f"{abs_z:.2f}"] = {
            "probe_count": int(len(pts)),
            "by_time": rows,
            "all_times_strictly_lower_than_both_and_inward": bool(band_pass),
        }
    return {
        "radii": list(PROBE_RADII),
        "abs_z_bands": list(VISIBLE_ABS_Z),
        "angles": [float(v) for v in PROBE_ANGLES],
        "times": [float(v) for v in PROBE_TIMES],
        "by_abs_z": by_abs_z,
        "all_visible_bands_all_times_strictly_lower_than_both_and_inward": bool(all_pass),
        "scope": "targeted development morphology fingerprint; not blind/global visual validation",
    }


def outer_reservoir_diagnostic(child_fn, p9_candidate_fn, candidate_fn) -> dict[str, Any]:
    """Record, but do not gate, total-field radial response in the return reservoir."""
    by_abs_z: dict[str, Any] = {}
    for abs_z in RESERVOIR_ABS_Z:
        pts = _probe_points(abs_z)
        rows: dict[str, Any] = {}
        for time in PROBE_TIMES:
            source = screen.loc._radial_velocity(child_fn, pts, time)
            old = screen.loc._radial_velocity(p9_candidate_fn, pts, time)
            new = screen.loc._radial_velocity(candidate_fn, pts, time)
            rows[f"{time:.3f}"] = {
                "source_652_mean_radial_velocity": float(np.mean(source)),
                "p9_740_mean_radial_velocity": float(np.mean(old)),
                "reservoir_child_mean_radial_velocity": float(np.mean(new)),
                "reservoir_child_inward_count": int(np.sum(new < 0.0)),
            }
        by_abs_z[f"{abs_z:.2f}"] = rows
    return {
        "abs_z_bands": list(RESERVOIR_ABS_Z),
        "by_abs_z": by_abs_z,
        "acceptance_gate_applied": False,
        "reason": "outer return-flow location is a project representation choice; no OpenAI/public numeric threshold exists",
    }


def response_diagnostic(linear_fn, child_fn) -> dict[str, Any]:
    rng = np.random.default_rng(screen.RESPONSE_SEED + 110)
    existing_blocks = []
    new_blocks = []
    for time in screen.RESPONSE_TIMES:
        radius = rng.uniform(0.2, 1.55, size=40)
        angle = rng.uniform(0.0, 2.0 * np.pi, size=40)
        z = rng.uniform(-1.95, 1.95, size=40)
        pts = np.column_stack((radius * np.cos(angle), radius * np.sin(angle), z))
        existing_blocks.append(
            (np.asarray(child_fn(pts, time), float) - np.asarray(linear_fn(pts, time), float)).ravel()
        )
        new_blocks.append(np.asarray(reservoir_unit_time_correction(pts, time), float).ravel())
    a = np.concatenate(existing_blocks)
    b = np.concatenate(new_blocks)
    na, nb = float(np.linalg.norm(a)), float(np.linalg.norm(b))
    if na <= np.finfo(float).tiny or nb <= np.finfo(float).tiny:
        raise RuntimeError("degenerate reservoir response diagnostic")
    matrix = np.column_stack((a / na, b / nb))
    singular = np.linalg.svd(matrix, compute_uv=False)
    rank = int(np.linalg.matrix_rank(matrix, tol=1.0e-10))
    condition = float(singular[0] / singular[-1])
    return {
        "parameter_count_increment_relative_to_652": 1,
        "basis_dimension_same_as_prior_single_tip_channel_screens": True,
        "rank": rank,
        "singular_values": singular.tolist(),
        "condition_number": condition,
        "column_cosine": float(np.dot(matrix[:, 0], matrix[:, 1])),
        "existing_response_norm": na,
        "reservoir_unit_response_norm": nb,
        "passes_rank_condition_guard": bool(rank == 2 and condition <= CONDITION_MAX),
    }


def correction_rms_fraction(child_fn, alpha: float) -> dict[str, float]:
    axis, velocity = screen.vp.base.morph.prior.sample_velocity_grid(
        child_fn, screen.MID_TIME, screen.GRID_RESOLUTION
    )
    pts = screen._grid_points(axis)
    perturb = float(alpha) * reservoir_unit_time_correction(pts, screen.MID_TIME)
    base_rms = float(np.sqrt(np.mean(np.sum(np.square(velocity.reshape(-1, 3)), axis=1))))
    perturb_rms = float(np.sqrt(np.mean(np.sum(np.square(perturb), axis=1))))
    if base_rms <= np.finfo(float).tiny:
        raise RuntimeError("degenerate source-652 midpoint RMS")
    return {
        "base_velocity_rms": base_rms,
        "reservoir_correction_rms": perturb_rms,
        "correction_to_base_rms_fraction": float(perturb_rms / base_rms),
    }


def structure_preflight(child_fn, candidate_fn) -> dict[str, Any]:
    rng = np.random.default_rng(screen.STRUCTURE_SEED + 110)
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
        [[r, 0.0, z] for z in (-0.85, 0.85) for r in PROBE_RADII], dtype=float
    )
    shoulder_correction = float(np.max(np.abs(reservoir.reservoir_correction(shoulder))))
    outside = np.asarray(
        [
            [1.61, 0.0, 1.55],
            [-1.61, 0.0, -1.55],
            [0.9, 0.0, 2.0],
            [0.9, 0.0, -2.0],
            [0.9, 0.0, 2.05],
            [0.9, 0.0, -2.05],
        ],
        dtype=float,
    )
    outside_correction = float(np.max(np.abs(reservoir.reservoir_correction(outside))))

    radius = rng.uniform(0.2, 1.50, size=40)
    angle = rng.uniform(0.0, 2.0 * np.pi, size=40)
    z_centers = rng.choice(np.asarray([-1.92, -1.70, -1.35, -0.75, 0.75, 1.35, 1.70, 1.92]), size=40)
    z = z_centers + rng.uniform(-0.015, 0.015, size=40)
    random_probes = np.column_stack((radius * np.cos(angle), radius * np.sin(angle), z))
    fixed_probes = np.vstack([_probe_points(v) for v in (*VISIBLE_ABS_Z, *RESERVOIR_ABS_Z)])
    probes = np.vstack((random_probes, fixed_probes))
    div_max = 0.0
    for time in PROBE_TIMES:
        div = np.zeros(len(probes), dtype=float)
        for axis_index in range(3):
            shift = np.zeros(3, dtype=float)
            shift[axis_index] = DIVERGENCE_STEP
            up = np.asarray(candidate_fn(probes + shift, time), float)
            um = np.asarray(candidate_fn(probes - shift, time), float)
            div += (up[:, axis_index] - um[:, axis_index]) / (2.0 * DIVERGENCE_STEP)
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
        "shoulder_correction_max_abs": shoulder_correction,
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
    alpha_res, calibration = derive_reservoir_alpha(child_fn, child_times, child_positions, metadata)
    alpha_p9, p9_calibration = p9.derive_reshaped_alpha(child_fn, child_times, child_positions, metadata)

    p9_candidate_fn = lambda pts, t: p9.reshaped_velocity(child_fn, pts, t, alpha_p9)
    candidate_fn = lambda pts, t: reservoir_velocity(child_fn, pts, t, alpha_res)

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

    grouped_child = screen._grouped_records_from_positions(child_fn, child_times, child_positions, metadata)
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
    visible = visible_tip_fingerprint(child_fn, p9_candidate_fn, candidate_fn)
    outer = outer_reservoir_diagnostic(child_fn, p9_candidate_fn, candidate_fn)
    rms = correction_rms_fraction(child_fn, alpha_res)

    justified = bool(
        structure["passes"]
        and radial["all_tip_mean_delta_r_strictly_lower"]
        and radial["all_tip_start_mean_radial_velocity_strictly_lower"]
        and radial["shoulder_24_of_24_inward_all_segments"]
        and turns["all_segments_positive"]
        and aspect["all_post_start_positive"]
        and response["passes_rank_condition_guard"]
        and visible["all_visible_bands_all_times_strictly_lower_than_both_and_inward"]
    )

    report = {
        "task_id": TASK_ID,
        "prereg_issue": PREREG_ISSUE,
        "source_parent": {"pr": SOURCE_PARENT_PR, "head": SOURCE_PARENT_HEAD},
        "source_context": {
            "source_visualization_child_pr": SOURCE_VISUAL_CHILD_PR,
            "baseline_temporal_pr": SOURCE_TEMPORAL_PR,
            "p9_comparator_pr": SOURCE_P9_PR,
            "fresh_holdout_pr": RELATED_FRESH_HOLDOUT_PR,
        },
        "frozen_reservoir_basis": {
            "representation": "curl(A), A=(-y*f,x*f,0), f=R(r^2)*Z(z)",
            "radial_support": "r<1.6",
            "visible_abs_z_interval": [1.0, 1.8],
            "return_reservoir_abs_z_interval": [1.8, 2.0],
            "time_activation": "g(t)=2*(t-.25)",
            "basis_dimension_relative_to_prior_tip_child_screens": "unchanged_one_tip_poloidal_channel",
            "coefficient_semantics": "rederived on reservoir unit basis from historical #701 development protocol",
        },
        "source_652_compact_swirl_epsilon": source_epsilon,
        "source_652_normalization": source_normalization,
        "coefficient_calibration": calibration,
        "p9_740_coefficient_calibration": p9_calibration,
        "coefficient_comparison": {
            "reservoir_alpha": float(alpha_res),
            "p9_740_alpha": float(alpha_p9),
            "reservoir_to_p9_alpha_ratio": float(alpha_res / alpha_p9) if alpha_p9 != 0.0 else None,
            "coefficient_scan_performed": False,
        },
        "midpoint_correction_rms": rms,
        "response_diagnostic": response,
        "structure_preflight": structure,
        "visible_tip_total_field_fingerprint": visible,
        "outer_reservoir_total_field_diagnostic": outer,
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
        "outer_reservoir_child_visualization_direction_screen_passed": justified,
        "decision": {
            "classification": (
                "outer_reservoir_nonlinear_screen_pass"
                if justified
                else "outer_reservoir_nonlinear_screen_fail"
            ),
            "carry_exact_child_to_disjoint_fresh_path_or_fixed_render_comparison": justified,
            "second_poloidal_channel_justified_now": False,
            "retune_alpha_or_reservoir_shape_from_this_result_allowed": False,
            "velocity_formula_evaluated_this_increment": True,
            "canonical_or_saved_velocity_promoted": False,
            "closer_visualization_delivery_supported": justified,
            "closer_visualization_delivery_scope": (
                "targeted development visible-tip fingerprint plus frozen morphology/path guards only; not global visual correspondence"
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
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("outputs/agent7-st052m-outer-reservoir-nonlinear-replay/report.json"),
    )
    args = parser.parse_args()
    report = run(args.out)
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
