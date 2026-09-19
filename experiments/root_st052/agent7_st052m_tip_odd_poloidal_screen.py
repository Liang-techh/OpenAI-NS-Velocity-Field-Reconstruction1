"""Screen one tip-band-local odd-poloidal channel on exact ST052-M #652.

Preregistered in issue #700 before child evaluation and stacked directly on
Agent-7 PR #692.  The new spatial channel is the analytic curl of
A=(-y*f,x*f,0) with f=R(r^2)*z*B(z^2).  Its odd axial factor makes f_z even,
so the preregistered positive orientation points radially inward at both frozen
tip seed signs while remaining exactly zero at the shoulder seed band.

The existing linear activation g(t)=2(t-.25) is reused; no new time basis is
introduced.  One coefficient is determined once by a frozen midpoint local
linear cancellation of exact #652's mean tip-band radial velocity.  No scan,
optimizer, post-result damping, pressure/forcing change, or image/pixel numeric
target is allowed.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np

import agent7_st052m_radial_split_localization as loc
import agent7_st052m_vector_potential_poloidal_screen as vp

TASK_ID = "CR003-ST052M-TIP-ODD-POLOIDAL-SCREEN-098"
PREREG_ISSUE = 700
SOURCE_PARENT_PR = 692
SOURCE_PARENT_HEAD = "e3262d4c03dc2d0806eb46639a6e14b2edb5fcb7"
SOURCE_WITNESS_PR = 652
SOURCE_TEMPORAL_PR = 587
PUBLIC_SOURCE_URL = "https://openai.com/index/navier-stokes-solution/"

GRID_RESOLUTION = 41
MID_TIME = 0.50
MID_INDEX = 16
R2_BUMP_A = -1.0
R2_BUMP_B = 2.56
Z2_BUMP_A = 1.0
Z2_BUMP_B = 3.24
R_SUPPORT_MAX = 1.6
Z_SUPPORT_INNER = 1.0
Z_SUPPORT_OUTER = 1.8
DIVERGENCE_STEP = 1.0e-5
DIVERGENCE_MAX = 1.0e-5
CONDITION_MAX = 25.0
RESPONSE_SEED = 9177001
STRUCTURE_SEED = 9177002
RESPONSE_TIMES = (0.375, 0.50, 0.625, 0.75)
POST_START_TIMES = (0.375, 0.50, 0.625, 0.75)

TRUTH = {
    "canonical_velocity_changed": False,
    "saved_velocity_changed": False,
    "production_candidate_selected": False,
    "screened_extra_tip_odd_poloidal_basis": True,
    "new_spatial_basis_promoted": False,
    "new_temporal_basis_added": False,
    "parameter_grid_scan_performed": False,
    "optimization_performed": False,
    "post_result_damping_or_retuning_performed": False,
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


def activation(time: float) -> float:
    time = float(time)
    if not (0.25 - 1.0e-12 <= time <= 0.75 + 1.0e-12):
        raise ValueError("time outside frozen interval")
    return float(2.0 * (time - 0.25))


def tip_odd_poloidal_correction(points: np.ndarray) -> np.ndarray:
    """Analytic curl(A), A=(-y*f,x*f,0), for the frozen tip-local odd channel."""
    pts = np.asarray(points, dtype=float)
    if pts.ndim != 2 or pts.shape[1] != 3:
        raise ValueError("points must have shape (n,3)")
    x, y, z = pts[:, 0], pts[:, 1], pts[:, 2]
    r2 = x * x + y * y
    z2 = z * z

    radial, radial_s = vp._poly_bump(r2, R2_BUMP_A, R2_BUMP_B)
    axial_bump, axial_q = vp._poly_bump(z2, Z2_BUMP_A, Z2_BUMP_B)

    axial = z * axial_bump
    axial_z = axial_bump + 2.0 * z2 * axial_q
    f = radial * axial
    f_s = radial_s * axial
    f_z = radial * axial_z

    correction = np.column_stack(
        (-x * f_z, -y * f_z, 2.0 * f + 2.0 * r2 * f_s)
    )
    if not np.isfinite(correction).all():
        raise RuntimeError("nonfinite tip odd-poloidal correction")
    return correction


def unit_time_correction(points: np.ndarray, time: float) -> np.ndarray:
    return activation(time) * tip_odd_poloidal_correction(np.asarray(points, dtype=float))


def corrected_velocity(base_fn, points: np.ndarray, time: float, alpha: float) -> np.ndarray:
    pts = np.asarray(points, dtype=float)
    return np.asarray(base_fn(pts, float(time)), dtype=float) + float(alpha) * unit_time_correction(
        pts, float(time)
    )


def _integrate(velocity_fn):
    return loc.parent._integrate_positions(velocity_fn)


def _grouped_records_from_positions(
    velocity_fn, times: np.ndarray, positions: np.ndarray, metadata: list[dict[str, Any]]
) -> dict[str, Any]:
    segments: dict[str, Any] = {}
    for i0, i1 in loc.parent.SEGMENTS:
        key = f"{times[i0]:.3f}-{times[i1]:.3f}"
        segments[key] = loc._interval_record(velocity_fn, times, positions, metadata, i0, i1)
    whole = loc._interval_record(
        velocity_fn,
        times,
        positions,
        metadata,
        loc.parent.CHECKPOINT_INDICES[0],
        loc.parent.CHECKPOINT_INDICES[-1],
    )
    return {"segments": segments, "whole_interval": whole}


def _path_records_from_positions(times: np.ndarray, positions: np.ndarray) -> dict[str, Any]:
    segments: dict[str, Any] = {}
    for i0, i1 in loc.parent.SEGMENTS:
        key = f"{times[i0]:.3f}-{times[i1]:.3f}"
        segments[key] = loc.parent._segment_summary(positions, i0, i1)
    whole = loc.parent._segment_summary(
        positions,
        loc.parent.CHECKPOINT_INDICES[0],
        loc.parent.CHECKPOINT_INDICES[-1],
    )
    return {
        "segments": segments,
        "whole_interval": whole,
        "checkpoint_times": [float(times[i]) for i in loc.parent.CHECKPOINT_INDICES],
    }


def derive_alpha_from_frozen_midpoint(
    child_fn,
    child_times: np.ndarray,
    child_positions: np.ndarray,
    metadata: list[dict[str, Any]],
) -> tuple[float, dict[str, Any]]:
    if abs(float(child_times[MID_INDEX]) - MID_TIME) > 1.0e-15:
        raise RuntimeError("frozen midpoint index drifted")
    tip_idx = np.asarray(
        [i for i, meta in enumerate(metadata) if str(meta["band"]) == "tip"], dtype=int
    )
    if len(tip_idx) != 24:
        raise RuntimeError("frozen tip path count drifted")
    pts = np.asarray(child_positions[MID_INDEX, tip_idx], dtype=float)
    base_ur = loc._radial_velocity(child_fn, pts, MID_TIME)
    unit_ur = loc._radial_velocity(unit_time_correction, pts, MID_TIME)
    mean_base = float(np.mean(base_ur))
    mean_unit = float(np.mean(unit_ur))
    if not (np.isfinite(mean_base) and np.isfinite(mean_unit)):
        raise RuntimeError("nonfinite midpoint radial calibration")
    if mean_base <= 0.0:
        raise RuntimeError("#692 diagnosed outward tip mean disappeared")
    if mean_unit >= 0.0 or abs(mean_unit) <= np.finfo(float).tiny:
        raise RuntimeError("preregistered positive correction orientation is not inward")
    alpha = float(-mean_base / mean_unit)
    if not (np.isfinite(alpha) and alpha > 0.0):
        raise RuntimeError("invalid frozen one-scalar alpha")
    fixed_position_mean_after = float(np.mean(base_ur + alpha * unit_ur))
    return alpha, {
        "time": MID_TIME,
        "tip_path_count": int(len(tip_idx)),
        "mean_base_radial_velocity": mean_base,
        "mean_unit_radial_velocity": mean_unit,
        "derived_alpha": alpha,
        "fixed_position_mean_after_linear_cancellation": fixed_position_mean_after,
        "coefficient_grid_scan_performed": False,
        "optimizer_used": False,
        "post_result_damping_allowed": False,
    }


def _grid_points(axis: np.ndarray) -> np.ndarray:
    xx, yy, zz = np.meshgrid(axis, axis, axis, indexing="ij")
    return np.column_stack((xx.ravel(), yy.ravel(), zz.ravel()))


def correction_rms_fraction(child_fn, alpha: float) -> dict[str, float]:
    axis, velocity = vp.base.morph.prior.sample_velocity_grid(
        child_fn, MID_TIME, GRID_RESOLUTION
    )
    pts = _grid_points(axis)
    perturb = float(alpha) * unit_time_correction(pts, MID_TIME)
    base_rms = float(np.sqrt(np.mean(np.sum(np.square(velocity.reshape(-1, 3)), axis=1))))
    perturb_rms = float(np.sqrt(np.mean(np.sum(np.square(perturb), axis=1))))
    if base_rms <= np.finfo(float).tiny:
        raise RuntimeError("degenerate midpoint child RMS")
    return {
        "base_velocity_rms": base_rms,
        "correction_rms": perturb_rms,
        "correction_to_base_rms_fraction": float(perturb_rms / base_rms),
    }


def response_diagnostic(linear_fn, child_fn) -> dict[str, Any]:
    rng = np.random.default_rng(RESPONSE_SEED)
    existing_blocks = []
    new_blocks = []
    for time in RESPONSE_TIMES:
        radius = rng.uniform(0.2, 1.55, size=32)
        angle = rng.uniform(0.0, 2.0 * np.pi, size=32)
        z = rng.uniform(-1.75, 1.75, size=32)
        pts = np.column_stack((radius * np.cos(angle), radius * np.sin(angle), z))
        existing_blocks.append(
            (np.asarray(child_fn(pts, time), float) - np.asarray(linear_fn(pts, time), float)).ravel()
        )
        new_blocks.append(np.asarray(unit_time_correction(pts, time), float).ravel())
    a = np.concatenate(existing_blocks)
    b = np.concatenate(new_blocks)
    na, nb = float(np.linalg.norm(a)), float(np.linalg.norm(b))
    if na <= np.finfo(float).tiny or nb <= np.finfo(float).tiny:
        raise RuntimeError("degenerate response diagnostic")
    matrix = np.column_stack((a / na, b / nb))
    singular = np.linalg.svd(matrix, compute_uv=False)
    return {
        "parameter_count_increment": 1,
        "rank": int(np.linalg.matrix_rank(matrix, tol=1.0e-10)),
        "singular_values": singular.tolist(),
        "condition_number": float(singular[0] / singular[-1]),
        "column_cosine": float(np.dot(matrix[:, 0], matrix[:, 1])),
        "existing_response_norm": na,
        "new_unit_response_norm": nb,
    }


def structure_preflight(child_fn, candidate_fn) -> dict[str, Any]:
    rng = np.random.default_rng(STRUCTURE_SEED)

    start_pts = rng.uniform(-1.9, 1.9, size=(96, 3))
    start_identity = float(
        np.max(
            np.abs(
                np.asarray(candidate_fn(start_pts, 0.25), float)
                - np.asarray(child_fn(start_pts, 0.25), float)
            )
        )
    )

    shoulder = []
    for z in (-0.85, 0.85):
        for r in (0.6, 0.9, 1.2):
            shoulder.append([r, 0.0, z])
    shoulder = np.asarray(shoulder, dtype=float)
    shoulder_correction = float(np.max(np.abs(tip_odd_poloidal_correction(shoulder))))

    outside = np.asarray(
        [
            [1.61, 0.0, 1.2],
            [-1.61, 0.0, -1.2],
            [0.9, 0.0, 0.99],
            [0.9, 0.0, -0.99],
            [0.9, 0.0, 1.81],
            [0.9, 0.0, -1.81],
        ],
        dtype=float,
    )
    correction_outside = float(np.max(np.abs(tip_odd_poloidal_correction(outside))))

    probes = rng.uniform(-1.45, 1.45, size=(40, 3))
    div_max = 0.0
    for time in (0.375, 0.50, 0.625, 0.75):
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
        and np.isfinite(correction_outside)
        and np.isfinite(div_max)
    )
    return {
        "start_identity_max_abs": start_identity,
        "shoulder_seed_correction_max_abs": shoulder_correction,
        "outside_correction_max_abs": correction_outside,
        "candidate_cartesian_fd_divergence_max": div_max,
        "finite": finite,
        "passes": bool(
            finite
            and start_identity <= 1.0e-12
            and shoulder_correction == 0.0
            and correction_outside == 0.0
            and div_max <= DIVERGENCE_MAX
        ),
    }


def _turns_retention(control_path, linear_path, child_path, candidate_path) -> dict[str, Any]:
    rows: dict[str, Any] = {}
    all_positive = True
    for key in control_path["segments"]:
        c = control_path["segments"][key]
        l = linear_path["segments"][key]
        h = child_path["segments"][key]
        q = candidate_path["segments"][key]
        lrel = loc.parent._relative(l["mean_absolute_turns"], c["mean_absolute_turns"])
        hrel = loc.parent._relative(h["mean_absolute_turns"], c["mean_absolute_turns"])
        qrel = loc.parent._relative(q["mean_absolute_turns"], c["mean_absolute_turns"])
        child_increment = float(hrel - lrel)
        candidate_increment = float(qrel - lrel)
        positive = bool(candidate_increment > 0.0)
        all_positive = all_positive and positive
        rows[key] = {
            "linear_relative_to_control": lrel,
            "source_652_relative_to_control": hrel,
            "candidate_relative_to_control": qrel,
            "source_652_minus_linear_fidelity_increment": child_increment,
            "candidate_minus_linear_fidelity_increment": candidate_increment,
            "candidate_retains_positive_fidelity_increment": positive,
        }
    return {"segments": rows, "all_segments_positive": bool(all_positive)}


def _radial_band_decision(source_grouped, candidate_grouped) -> dict[str, Any]:
    rows: dict[str, Any] = {}
    all_tip_delta_better = True
    all_tip_start_ur_better = True
    all_shoulder_inward = True
    tip_count_changes = []
    for key in source_grouped["segments"]:
        s = source_grouped["segments"][key]["groups"]["band"]
        q = candidate_grouped["segments"][key]["groups"]["band"]
        st, qt = s["tip"], q["tip"]
        ss, qs = s["shoulder"], q["shoulder"]
        delta_better = bool(qt["mean_delta_r"] < st["mean_delta_r"])
        ur_better = bool(qt["mean_start_radial_velocity"] < st["mean_start_radial_velocity"])
        shoulder_guard = bool(qs["inward_path_count"] == 24 and qs["path_count"] == 24)
        count_change = int(qt["inward_path_count"] - st["inward_path_count"])
        all_tip_delta_better = all_tip_delta_better and delta_better
        all_tip_start_ur_better = all_tip_start_ur_better and ur_better
        all_shoulder_inward = all_shoulder_inward and shoulder_guard
        tip_count_changes.append(count_change)
        rows[key] = {
            "source_652_tip": st,
            "candidate_tip": qt,
            "source_652_shoulder": ss,
            "candidate_shoulder": qs,
            "tip_mean_delta_r_strictly_lower": delta_better,
            "tip_start_mean_radial_velocity_strictly_lower": ur_better,
            "shoulder_24_of_24_inward_guard": shoulder_guard,
            "tip_inward_path_count_change": count_change,
        }
    return {
        "segments": rows,
        "all_tip_mean_delta_r_strictly_lower": bool(all_tip_delta_better),
        "all_tip_start_mean_radial_velocity_strictly_lower": bool(all_tip_start_ur_better),
        "shoulder_24_of_24_inward_all_segments": bool(all_shoulder_inward),
        "tip_inward_path_count_changes": tip_count_changes,
        "any_tip_inward_path_count_gain": bool(any(v > 0 for v in tip_count_changes)),
    }


def _aspect_decision(control_fn, linear_fn, candidate_fn) -> dict[str, Any]:
    rows: dict[str, Any] = {}
    all_positive = True
    for time in POST_START_TIMES:
        row = loc.parent.temporal.time_record(control_fn, linear_fn, candidate_fn, time)
        increment = float(row["child_minus_linear"]["full_aspect_ratio"])
        positive = bool(increment > 0.0)
        all_positive = all_positive and positive
        rows[f"{time:.3f}"] = {
            "candidate_minus_linear_full_aspect_ratio": increment,
            "positive": positive,
            "candidate_relative_to_control": row["relative_to_control"]["child"]["full_aspect_ratio"],
            "linear_relative_to_control": row["relative_to_control"]["linear"]["full_aspect_ratio"],
        }
    return {"by_time": rows, "all_post_start_positive": bool(all_positive)}


def _build_fields():
    control_fn, linear_fn, child_fn, epsilon, normalization = loc._build_fields()
    if abs(loc.parent.witness.SWIRL_A - 0.065899695471146) > 0.0:
        raise RuntimeError("#652 swirl coordinate drifted")
    if abs(loc.parent.witness.SHOULDER_LAMBDA - (-0.02)) > 1.0e-15:
        raise RuntimeError("#652 shoulder timing coordinate drifted")
    return control_fn, linear_fn, child_fn, epsilon, normalization


def run(out: Path) -> dict[str, Any]:
    control_fn, linear_fn, child_fn, source_epsilon, source_normalization = _build_fields()

    child_times, child_positions, metadata = _integrate(child_fn)
    alpha, calibration = derive_alpha_from_frozen_midpoint(
        child_fn, child_times, child_positions, metadata
    )
    candidate_fn = lambda p, t: corrected_velocity(child_fn, p, t, alpha)

    control_times, control_positions, control_metadata = _integrate(control_fn)
    linear_times, linear_positions, linear_metadata = _integrate(linear_fn)
    candidate_times, candidate_positions, candidate_metadata = _integrate(candidate_fn)
    for name, times, meta in (
        ("control", control_times, control_metadata),
        ("linear", linear_times, linear_metadata),
        ("candidate", candidate_times, candidate_metadata),
    ):
        if not np.array_equal(times, child_times) or meta != metadata:
            raise RuntimeError(f"{name} trajectory protocol identity drifted")

    grouped_child = _grouped_records_from_positions(
        child_fn, child_times, child_positions, metadata
    )
    grouped_candidate = _grouped_records_from_positions(
        candidate_fn, candidate_times, candidate_positions, metadata
    )
    radial = _radial_band_decision(grouped_child, grouped_candidate)

    path_control = _path_records_from_positions(control_times, control_positions)
    path_linear = _path_records_from_positions(linear_times, linear_positions)
    path_child = _path_records_from_positions(child_times, child_positions)
    path_candidate = _path_records_from_positions(candidate_times, candidate_positions)
    turns = _turns_retention(path_control, path_linear, path_child, path_candidate)

    aspect = _aspect_decision(control_fn, linear_fn, candidate_fn)
    response = response_diagnostic(linear_fn, child_fn)
    structure = structure_preflight(child_fn, candidate_fn)
    rms = correction_rms_fraction(child_fn, alpha)

    response_guard = bool(
        response["rank"] == 2 and response["condition_number"] <= CONDITION_MAX
    )
    justified = bool(
        structure["passes"]
        and radial["all_tip_mean_delta_r_strictly_lower"]
        and radial["all_tip_start_mean_radial_velocity_strictly_lower"]
        and radial["shoulder_24_of_24_inward_all_segments"]
        and turns["all_segments_positive"]
        and aspect["all_post_start_positive"]
        and response_guard
    )

    routing = (
        "the one frozen tip-band odd-poloidal channel supplies the previously missing inward radial direction while preserving shoulder inward motion, #652/#587 winding advantage, axial-stretch advantage and structural guards; this justifies carrying exactly this extra spatial channel to one independent fixed-render/save-load comparison before any further basis growth"
        if justified
        else "the frozen tip-band odd-poloidal channel does not satisfy every preregistered visualization/structure guard; retain exact #652 and do not retune this coefficient or grow the basis from this failed screen"
    )

    report = {
        "task_id": TASK_ID,
        "prereg_issue": PREREG_ISSUE,
        "source_parent": {"pr": SOURCE_PARENT_PR, "head": SOURCE_PARENT_HEAD},
        "candidate_lineage": {
            "baseline_pr": SOURCE_TEMPORAL_PR,
            "source_visualization_child_pr": SOURCE_WITNESS_PR,
        },
        "public_source": {"url": PUBLIC_SOURCE_URL, "qualitative_only": True},
        "frozen_basis": {
            "representation": "curl(A), A=(-y*f,x*f,0), f=R(r^2)*z*B(z^2)",
            "radial_bump_r2_interval": [R2_BUMP_A, R2_BUMP_B],
            "physical_radial_support_max": R_SUPPORT_MAX,
            "axial_bump_z2_interval": [Z2_BUMP_A, Z2_BUMP_B],
            "physical_axial_support_abs": [Z_SUPPORT_INNER, Z_SUPPORT_OUTER],
            "time_activation": "g(t)=2*(t-.25)",
            "new_spatial_channel_count": 1,
            "new_temporal_channel_count": 0,
            "coefficient_rule": "one frozen midpoint mean tip radial-velocity cancellation",
            "parameter_grid_scan_performed": False,
            "optimizer_used": False,
        },
        "source_652_compact_swirl_epsilon": source_epsilon,
        "source_652_normalization": source_normalization,
        "coefficient_calibration": calibration,
        "midpoint_correction_rms": rms,
        "response_diagnostic": response,
        "response_rank_condition_guard": response_guard,
        "radial_band_replay": {
            "source_652": grouped_child,
            "candidate": grouped_candidate,
            "decision": radial,
        },
        "turns_replay": {
            "control": path_control,
            "linear_587": path_linear,
            "source_652": path_child,
            "candidate": path_candidate,
            "decision": turns,
        },
        "axial_stretch_replay": aspect,
        "structure_preflight": structure,
        "tip_odd_poloidal_basis_growth_justified": justified,
        "source_652_target_free_pareto_rejection_remains_binding": True,
        "source_683_temporal_inward_spiral_failure_remains_historical_for_652": True,
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
    print("alpha=", report["coefficient_calibration"]["derived_alpha"])
    print("growth_justified=", report["tip_odd_poloidal_basis_growth_justified"])
    print("radial=", report["radial_band_replay"]["decision"])
    print("turns=", report["turns_replay"]["decision"])
    print("aspect=", report["axial_stretch_replay"])
    print("response=", report["response_diagnostic"])
    print("structure=", report["structure_preflight"])
    print("routing=", report["routing"])


if __name__ == "__main__":
    main()
