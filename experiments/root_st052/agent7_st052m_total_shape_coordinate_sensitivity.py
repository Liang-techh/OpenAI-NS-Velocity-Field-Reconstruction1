"""Audit the recalibrated total-child radial-shape tangent for ST052-M.

Frozen in issue #809 and stacked directly on exact Agent-7 PR #801 head.
This is one expression-capacity / sensitivity increment.  It changes no
candidate velocity and selects no new k.

The outer-reservoir channel keeps

    A=(-y*f,x*f,0),  f=R_k(r^2) Z(z),
    R_k(q) proportional to ((q+k)(64/25-q))^4,

with local k in {0.9,1.0,1.1}.  For every k we restore the physical scalar
coordinate by re-deriving alpha with the exact historical #775 development
rule on the same 24 parent-tip positions at t=.50.  The main question is then
whether the resulting total-family shape tangent is independent of the live
scalar-amplitude tangent after that recalibration.

All radial/axial/vorticity quantities below are autonomous repository
fingerprints.  They are not OpenAI-published numerical targets and no sign of a
sensitivity is called an improvement by itself.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Callable

import numpy as np

import agent7_st052m_outer_axial_return_reservoir_preflight as reservoir
import agent7_st052m_outer_reservoir_nonlinear_replay as replay
import agent7_st052m_radial_core_control_preflight as core_control
import agent7_st052m_total_child_axial_stretch_fingerprint as total_axial

TASK_ID = "CR003-ST052M-TOTAL-SHAPE-COORDINATE-SENSITIVITY-115"
PREREG_ISSUE = 809
SOURCE_PARENT_PR = 801
SOURCE_PARENT_HEAD = "5eb09a8ea6a6a4f79b7deee020690eba1612b5e2"
SOURCE_CHILD_PR = 775
SOURCE_CHILD_HEAD = "93887a59729d22113badf2ae4dac2f7d868e2703"
RELATED_FRESH_HOLDOUT_PR = 714

K_MINUS = 0.90
K_LIVE = 1.00
K_PLUS = 1.10
K_VALUES = (K_MINUS, K_LIVE, K_PLUS)
K_STEP = 0.10
IDENTIFIABILITY_COSINE_MAX_ABS = 0.995
IDENTIFIABILITY_CONDITION_MAX = 25.0
RESPONSE_RADII = (0.35, 0.65, 0.95, 1.25, 1.50)
RESPONSE_ABS_Z = (1.25, 1.55, 1.75, 1.85, 1.95)
RESPONSE_ANGLES = tuple(float(v) for v in np.arange(8) * (0.25 * np.pi))
RESPONSE_TIMES = (0.375, 0.50, 0.625, 0.75)
MORPHOLOGY_GRID_RESOLUTION = 33
MORPHOLOGY_TIME = 0.50
ALPHA_REPLAY_TOL = 2.0e-12

TRUTH = {
    "candidate_velocity_changed": False,
    "canonical_velocity_changed": False,
    "saved_velocity_changed": False,
    "basis_dimension_changed": False,
    "new_spatial_basis_added": False,
    "new_temporal_basis_added": False,
    "candidate_k_selected": False,
    "parameter_scan_performed": False,
    "optimization_performed": False,
    "post_result_retuning_performed": False,
    "historical_development_paths_used_only_for_frozen_alpha_semantics": True,
    "fresh_714_path_data_used": False,
    "fresh_714_path_data_used_for_retuning": False,
    "pressure_or_force_changed": False,
    "held_out_pde_residual_evaluated": False,
    "public_image_numeric_target_used": False,
    "public_numeric_core_radius_target_defined": False,
    "pixel_similarity_objective_used": False,
    "total_field_visualization_fingerprints_recorded": True,
    "closer_visualization_delivery_established": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "source_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _assert_source_lock() -> None:
    core_control._assert_source_lock()
    replay._assert_source_lock()
    total_axial._assert_source_lock()
    if core_control.TASK_ID != "CR003-ST052M-RADIAL-CORE-CONTROL-PREFLIGHT-114":
        raise RuntimeError("#801 radial core-control task identity drifted")
    if core_control.PREREG_ISSUE != 800:
        raise RuntimeError("#801 preregistration identity drifted")
    if replay.TASK_ID != "CR003-ST052M-OUTER-RESERVOIR-NONLINEAR-REPLAY-110":
        raise RuntimeError("#775 child task identity drifted")
    if total_axial.TASK_ID != "CR003-ST052M-TOTAL-CHILD-AXIAL-STRETCH-FINGERPRINT-113":
        raise RuntimeError("#794 total axial task identity drifted")
    if tuple(core_control.K_VALUES) != (0.50, 0.75, 1.00, 1.50, 2.00):
        raise RuntimeError("#801 frozen broad k family drifted")
    if reservoir.Z_VISIBLE_OUTER != 1.8 or reservoir.Z_GLOBAL_OUTER != 2.0:
        raise RuntimeError("outer axial reservoir semantics drifted")
    if RELATED_FRESH_HOLDOUT_PR != 714:
        raise RuntimeError("fresh holdout identity drifted")


def reservoir_correction_k(points: np.ndarray, k: float) -> np.ndarray:
    """Exact #767 curl formula with only #801's radial envelope coordinate changed."""
    pts = np.asarray(points, dtype=float)
    if pts.ndim != 2 or pts.shape[1] != 3:
        raise ValueError("points must have shape (n,3)")
    if float(k) <= 0.0:
        raise ValueError("k must be positive")
    x, y, z = pts[:, 0], pts[:, 1], pts[:, 2]
    r2 = x * x + y * y
    radial, radial_q = core_control.radial_profile_k(r2, float(k))
    axial, axial_z = reservoir.axial_potential_and_derivative(z)
    f = radial * axial
    f_q = radial_q * axial
    f_z = radial * axial_z
    correction = np.column_stack((-x * f_z, -y * f_z, 2.0 * f + 2.0 * r2 * f_q))
    if not np.all(np.isfinite(correction)):
        raise RuntimeError("nonfinite k-shaped reservoir correction")
    return correction


def unit_time_correction_k(points: np.ndarray, time: float, k: float) -> np.ndarray:
    return replay.screen.activation(float(time)) * reservoir_correction_k(points, float(k))


def candidate_velocity_k(
    parent_fn: Callable[[np.ndarray, float], np.ndarray],
    points: np.ndarray,
    time: float,
    k: float,
    alpha: float,
) -> np.ndarray:
    pts = np.asarray(points, dtype=float)
    return np.asarray(parent_fn(pts, float(time)), dtype=float) + float(alpha) * unit_time_correction_k(
        pts, float(time), float(k)
    )


def _historical_parent_state() -> tuple[Callable, np.ndarray, np.ndarray, list[dict[str, Any]]]:
    """Build exact #652 parent and the historical trajectory coordinates used by #775."""
    _control_fn, _linear_fn, parent_fn, _epsilon, _normalization = replay.screen._build_fields()
    parent_times, parent_positions, metadata = replay.screen._integrate(parent_fn)
    return parent_fn, parent_times, parent_positions, metadata


def derive_alpha_k(
    parent_fn: Callable,
    parent_times: np.ndarray,
    parent_positions: np.ndarray,
    metadata: list[dict[str, Any]],
    k: float,
) -> tuple[float, dict[str, Any]]:
    """Apply #775's exact one-shot 24-tip rule to the k-shaped unit response."""
    if abs(float(parent_times[replay.screen.MID_INDEX]) - replay.screen.MID_TIME) > 1.0e-15:
        raise RuntimeError("historical midpoint index drifted")
    tip_idx = np.asarray(
        [i for i, meta in enumerate(metadata) if str(meta["band"]) == "tip"], dtype=int
    )
    if len(tip_idx) != 24:
        raise RuntimeError("historical tip path count drifted")
    pts = np.asarray(parent_positions[replay.screen.MID_INDEX, tip_idx], dtype=float)
    base_ur = replay.screen.loc._radial_velocity(parent_fn, pts, replay.screen.MID_TIME)
    unit_fn = lambda p, t: unit_time_correction_k(p, t, float(k))
    unit_ur = replay.screen.loc._radial_velocity(unit_fn, pts, replay.screen.MID_TIME)
    mean_base = float(np.mean(base_ur))
    mean_unit = float(np.mean(unit_ur))
    if not (np.isfinite(mean_base) and np.isfinite(mean_unit)):
        raise RuntimeError("nonfinite historical alpha coordinate")
    if abs(mean_unit) <= np.finfo(float).tiny:
        raise RuntimeError("degenerate k-shaped unit response on historical tips")
    alpha = float(-mean_base / mean_unit)
    if not np.isfinite(alpha) or abs(alpha) <= np.finfo(float).tiny:
        raise RuntimeError("invalid recalibrated alpha")
    return alpha, {
        "k": float(k),
        "time": float(replay.screen.MID_TIME),
        "tip_path_count": int(len(tip_idx)),
        "mean_source_652_radial_velocity": mean_base,
        "mean_unit_radial_velocity": mean_unit,
        "derived_alpha": alpha,
        "fixed_position_mean_after_linear_cancellation": float(np.mean(base_ur + alpha * unit_ur)),
        "historical_development_protocol_used": True,
        "fresh_714_data_used": False,
        "coefficient_scan_performed": False,
        "optimizer_used": False,
    }


def build_local_family() -> tuple[Callable, dict[float, float], dict[float, dict[str, Any]], dict[float, Callable]]:
    """Construct the three frozen total children with independently restored alpha semantics."""
    _assert_source_lock()
    parent_fn, parent_times, parent_positions, metadata = _historical_parent_state()
    alphas: dict[float, float] = {}
    calibrations: dict[float, dict[str, Any]] = {}
    children: dict[float, Callable] = {}
    for k in K_VALUES:
        alpha, calibration = derive_alpha_k(
            parent_fn, parent_times, parent_positions, metadata, float(k)
        )
        alphas[float(k)] = float(alpha)
        calibrations[float(k)] = calibration
        children[float(k)] = (
            lambda pts, time, kk=float(k), aa=float(alpha): candidate_velocity_k(
                parent_fn, pts, time, kk, aa
            )
        )

    # k=1 must reproduce #775's coordinate semantics, not merely look similar.
    live_alpha, live_calibration = replay.derive_reservoir_alpha(
        parent_fn, parent_times, parent_positions, metadata
    )
    rel = abs(alphas[K_LIVE] - float(live_alpha)) / max(abs(float(live_alpha)), np.finfo(float).tiny)
    if rel > ALPHA_REPLAY_TOL:
        raise RuntimeError(f"k=1 alpha does not replay #775: relative mismatch={rel}")
    calibrations[K_LIVE]["exact_775_alpha"] = float(live_alpha)
    calibrations[K_LIVE]["relative_mismatch_vs_exact_775"] = float(rel)
    calibrations[K_LIVE]["exact_775_calibration"] = live_calibration
    return parent_fn, alphas, calibrations, children


def _response_points() -> np.ndarray:
    rows: list[list[float]] = []
    for radius in RESPONSE_RADII:
        for abs_z in RESPONSE_ABS_Z:
            for z_sign in (-1.0, 1.0):
                for angle in RESPONSE_ANGLES:
                    rows.append(
                        [
                            float(radius * math.cos(angle)),
                            float(radius * math.sin(angle)),
                            float(z_sign * abs_z),
                        ]
                    )
    pts = np.asarray(rows, dtype=float)
    expected = len(RESPONSE_RADII) * len(RESPONSE_ABS_Z) * 2 * len(RESPONSE_ANGLES)
    if pts.shape != (expected, 3):
        raise RuntimeError("deterministic response cloud drifted")
    return pts


def total_family_identifiability(
    parent_fn: Callable,
    alphas: dict[float, float],
    children: dict[float, Callable],
) -> dict[str, Any]:
    """Compare the live amplitude tangent with the recalibrated local shape tangent."""
    pts = _response_points()
    amplitude_blocks: list[np.ndarray] = []
    shape_blocks: list[np.ndarray] = []
    for time in RESPONSE_TIMES:
        parent = np.asarray(parent_fn(pts, time), dtype=float)
        live = np.asarray(children[K_LIVE](pts, time), dtype=float)
        minus = np.asarray(children[K_MINUS](pts, time), dtype=float)
        plus = np.asarray(children[K_PLUS](pts, time), dtype=float)
        amplitude_blocks.append((live - parent).ravel())
        shape_blocks.append(((plus - minus) / (K_PLUS - K_MINUS)).ravel())

    amplitude = np.concatenate(amplitude_blocks)
    shape = np.concatenate(shape_blocks)
    n_amp = float(np.linalg.norm(amplitude))
    n_shape = float(np.linalg.norm(shape))
    if min(n_amp, n_shape) <= np.finfo(float).tiny:
        raise RuntimeError("degenerate total-family tangent")
    matrix = np.column_stack((amplitude / n_amp, shape / n_shape))
    singular = np.linalg.svd(matrix, compute_uv=False)
    rank = int(np.linalg.matrix_rank(matrix, tol=1.0e-10))
    cosine = float(np.dot(matrix[:, 0], matrix[:, 1]))
    condition = float(singular[0] / singular[-1])
    passed = bool(
        rank == 2
        and abs(cosine) < IDENTIFIABILITY_COSINE_MAX_ABS
        and condition <= IDENTIFIABILITY_CONDITION_MAX
    )
    return {
        "response_point_count_per_time": int(len(pts)),
        "times": list(RESPONSE_TIMES),
        "amplitude_tangent_definition": "u_k1 - u_parent = alpha_1*g(t)*C_1",
        "shape_tangent_definition": "(u_k1p1-u_k0p9)/0.20 with each alpha_k independently replayed by #775 rule",
        "amplitude_tangent_norm": n_amp,
        "shape_tangent_norm_per_k": n_shape,
        "normalized_column_rank": rank,
        "normalized_column_cosine": cosine,
        "normalized_two_column_singular_values": singular.tolist(),
        "normalized_two_column_condition": condition,
        "cosine_gate_abs_max": IDENTIFIABILITY_COSINE_MAX_ABS,
        "condition_gate_max": IDENTIFIABILITY_CONDITION_MAX,
        "passes": passed,
        "alpha_live": float(alphas[K_LIVE]),
    }


def visible_radial_fingerprints(children: dict[float, Callable]) -> dict[str, Any]:
    by_k: dict[str, Any] = {}
    compact: dict[float, dict[tuple[float, float], float]] = {}
    for k in K_VALUES:
        rows: dict[str, Any] = {}
        compact[float(k)] = {}
        for abs_z in replay.VISIBLE_ABS_Z:
            pts = replay._probe_points(float(abs_z))
            by_time: dict[str, Any] = {}
            for time in replay.PROBE_TIMES:
                ur = replay.screen.loc._radial_velocity(children[float(k)], pts, float(time))
                mean = float(np.mean(ur))
                compact[float(k)][(float(abs_z), float(time))] = mean
                by_time[f"{time:.3f}"] = {
                    "mean_total_radial_velocity": mean,
                    "inward_count": int(np.sum(ur < 0.0)),
                    "probe_count": int(len(ur)),
                }
            rows[f"|z|={abs_z:.2f}"] = by_time
        by_k[f"k={k:.2f}"] = rows

    sensitivity: dict[str, float] = {}
    for abs_z in replay.VISIBLE_ABS_Z:
        for time in replay.PROBE_TIMES:
            key = f"abs_z={abs_z:.2f},t={time:.3f}"
            sensitivity[key] = float(
                (compact[K_PLUS][(float(abs_z), float(time))] - compact[K_MINUS][(float(abs_z), float(time))])
                / (K_PLUS - K_MINUS)
            )
    return {
        "by_k": by_k,
        "symmetric_d_mean_radial_velocity_dk": sensitivity,
        "all_local_children_inward_on_all_registered_means": bool(
            all(value < 0.0 for values in compact.values() for value in values.values())
        ),
        "acceptance_gate_applied": False,
        "reason": "these are autonomous development fingerprints; #775/#794 upstream scientific gates remain separately unresolved",
    }


def _positive_weighted_effective_radius(profile: np.ndarray) -> float | None:
    positive = np.maximum(np.asarray(profile, dtype=float), 0.0)
    total = float(np.sum(positive))
    if total <= np.finfo(float).tiny:
        return None
    return float(np.sqrt(np.sum(np.square(total_axial.RADII) * positive) / total))


def axial_stretch_fingerprints(children: dict[float, Callable]) -> dict[str, Any]:
    by_k: dict[str, Any] = {}
    effective: dict[float, dict[tuple[float, float], float | None]] = {}
    crossings: dict[float, dict[tuple[float, float], float | None]] = {}
    for k in K_VALUES:
        rows: dict[str, Any] = {}
        effective[float(k)] = {}
        crossings[float(k)] = {}
        for abs_z in total_axial.ABS_Z_BANDS:
            by_time: dict[str, Any] = {}
            for time in total_axial.PROBE_TIMES:
                profile = total_axial.axial_stretch_profile(children[float(k)], float(abs_z), float(time))
                eff = _positive_weighted_effective_radius(profile)
                crossing = total_axial.first_positive_to_nonpositive_crossing(total_axial.RADII, profile)
                effective[float(k)][(float(abs_z), float(time))] = eff
                crossings[float(k)][(float(abs_z), float(time))] = crossing
                by_time[f"{time:.3f}"] = {
                    "positive_weighted_effective_radius": eff,
                    "first_total_positive_to_nonpositive_crossing_radius": crossing,
                    "positive_stretch_radius_fraction": float(np.mean(profile > 0.0)),
                    "mean_total_away_from_midplane_w": float(np.mean(profile)),
                    "profile": profile.tolist(),
                }
            rows[f"|z|={abs_z:.2f}"] = by_time
        by_k[f"k={k:.2f}"] = rows

    eff_sensitivity: dict[str, float | None] = {}
    crossing_sensitivity: dict[str, float | None] = {}
    for abs_z in total_axial.ABS_Z_BANDS:
        for time in total_axial.PROBE_TIMES:
            coord = (float(abs_z), float(time))
            key = f"abs_z={abs_z:.2f},t={time:.3f}"
            lo, hi = effective[K_MINUS][coord], effective[K_PLUS][coord]
            eff_sensitivity[key] = None if lo is None or hi is None else float((hi - lo) / (K_PLUS - K_MINUS))
            clo, chi = crossings[K_MINUS][coord], crossings[K_PLUS][coord]
            crossing_sensitivity[key] = None if clo is None or chi is None else float((chi - clo) / (K_PLUS - K_MINUS))
    return {
        "by_k": by_k,
        "symmetric_d_positive_weighted_effective_radius_dk": eff_sensitivity,
        "symmetric_d_total_crossing_radius_dk": crossing_sensitivity,
        "acceptance_gate_applied": False,
        "public_openai_numeric_core_radius_target": None,
    }


def _morphology_metrics(field_fn: Callable) -> dict[str, float]:
    prior = replay.screen.vp.base.morph.prior
    axis, velocity = prior.sample_velocity_grid(field_fn, MORPHOLOGY_TIME, MORPHOLOGY_GRID_RESOLUTION)
    spacing = float(axis[1] - axis[0])
    _magnitude, omega = prior.vorticity(velocity, spacing)
    metrics = replay.screen.vp.base.morph.enstrophy_moment_metrics(omega, axis)
    required = ("full_aspect_ratio", "smooth_tip_radial_rms")
    out: dict[str, float] = {}
    for key in required:
        if key not in metrics:
            raise RuntimeError(f"repository morphology metric missing required key: {key}")
        value = float(metrics[key])
        if not np.isfinite(value):
            raise RuntimeError(f"nonfinite repository morphology metric: {key}")
        out[key] = value
    return out


def vorticity_morphology_fingerprints(children: dict[float, Callable]) -> dict[str, Any]:
    by_k = {f"k={k:.2f}": _morphology_metrics(children[float(k)]) for k in K_VALUES}
    lo = by_k[f"k={K_MINUS:.2f}"]
    live = by_k[f"k={K_LIVE:.2f}"]
    hi = by_k[f"k={K_PLUS:.2f}"]
    sensitivity = {
        key: float((hi[key] - lo[key]) / (K_PLUS - K_MINUS))
        for key in ("full_aspect_ratio", "smooth_tip_radial_rms")
    }
    return {
        "grid_resolution": MORPHOLOGY_GRID_RESOLUTION,
        "time": MORPHOLOGY_TIME,
        "by_k": by_k,
        "live_k_metrics": live,
        "symmetric_d_metric_dk": sensitivity,
        "acceptance_gate_applied": False,
        "reason": "sensitivity is recorded, not optimized against an invented public numeric target",
    }


def run(out: Path) -> dict[str, Any]:
    parent_fn, alphas, calibrations, children = build_local_family()
    identifiability = total_family_identifiability(parent_fn, alphas, children)
    visible = visible_radial_fingerprints(children)
    axial = axial_stretch_fingerprints(children)
    vorticity = vorticity_morphology_fingerprints(children)

    live_alpha = float(alphas[K_LIVE])
    coefficient_rows = {
        f"k={k:.2f}": {
            "alpha": float(alphas[float(k)]),
            "alpha_ratio_to_live_k1": float(alphas[float(k)] / live_alpha),
            "calibration": calibrations[float(k)],
        }
        for k in K_VALUES
    }

    report: dict[str, Any] = {
        "task_id": TASK_ID,
        "prereg_issue": PREREG_ISSUE,
        "source_parent": {"pr": SOURCE_PARENT_PR, "head": SOURCE_PARENT_HEAD},
        "source_frozen_child": {"pr": SOURCE_CHILD_PR, "head": SOURCE_CHILD_HEAD},
        "fresh_holdout_pr_not_consumed": RELATED_FRESH_HOLDOUT_PR,
        "family": {
            "formula": "R_k(q) proportional to ((q+k)(64/25-q))^4 with #767 axial reservoir",
            "k_values": list(K_VALUES),
            "live_k": K_LIVE,
            "coefficient_semantics": "alpha_k replayed independently by exact #775 historical 24-tip t=.50 rule",
            "candidate_selection": None,
            "public_openai_numeric_target": None,
        },
        "coefficients": coefficient_rows,
        "total_family_identifiability": identifiability,
        "visible_radial_fingerprints": visible,
        "total_axial_stretch_fingerprints": axial,
        "vorticity_morphology_fingerprints": vorticity,
        "decision": {
            "same_dimension_total_shape_tangent_identifiable": bool(identifiability["passes"]),
            "second_poloidal_basis_justified_by_this_increment": bool(not identifiability["passes"]),
            "candidate_k_selected": False,
            "actual_velocity_changed": False,
            "closer_visualization_delivery_established": False,
            "if_identifiable_route": "retain one same-dimension radial-envelope reshape as first response to a later preregistered core-width discrepancy",
            "if_not_identifiable_route": "stop same-dimension radial reshaping and test one genuinely independent poloidal basis direction",
            "broader_visual_claim_requires_frozen_candidate_and_disjoint_comparison": True,
        },
        **TRUTH,
    }
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    Path(out).write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    report = run(args.out)
    compact = {
        "coefficients": {
            k: {"alpha": v["alpha"], "alpha_ratio_to_live_k1": v["alpha_ratio_to_live_k1"]}
            for k, v in report["coefficients"].items()
        },
        "total_family_identifiability": report["total_family_identifiability"],
        "vorticity_sensitivity": report["vorticity_morphology_fingerprints"]["symmetric_d_metric_dk"],
        "decision": report["decision"],
    }
    print(json.dumps(compact, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
