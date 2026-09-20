"""Preflight one same-dimension axial-turnover coordinate for ST052-M.

Preregistered in issue #861 before observing new numerical output and stacked
on exact Agent-7 PR #855 head.  This increment changes no selected candidate.
It asks whether moving the outer-reservoir turnover location, while preserving
the same compact support and restoring the exact #775 scalar-coordinate
semantics, supplies a morphology direction independent of scalar amplitude and
the already-audited radial-envelope shape coordinate.

For

    A=(-y*f,x*f,0),  f=R(r^2) Z_s(z),  C_s=curl(A),

the magnitude of the odd C3 axial potential rises from zero to one over
1<|z|<s and falls back to zero over s<|z|<2.  Therefore C_r is inward before
s and outward after s.  Varying s changes where the visible tip lobe hands off
to the mandatory compact-support return reservoir without adding a second
spatial basis, changing global support, or changing the time law.

All total-field fingerprints are autonomous repository diagnostics.  No
OpenAI-published numerical turnover location, axial length, or taper target is
assumed, and no sensitivity sign is called visual correspondence.
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
import agent7_st052m_total_shape_coordinate_sensitivity as radial_shape

TASK_ID = "CR003-ST052M-AXIAL-TURNOVER-SHAPE-PREFLIGHT-122"
PREREG_ISSUE = 861
SOURCE_PARENT_PR = 855
SOURCE_PARENT_HEAD = "47a59dee5ec5bb0a1bcd61c04796638640bb26ed"
SOURCE_CHILD_PR = 775
SOURCE_CHILD_HEAD = "93887a59729d22113badf2ae4dac2f7d868e2703"
SOURCE_RADIAL_SHAPE_PR = 810
SOURCE_RADIAL_SHAPE_HEAD = "b077eaffabfe66ce0b15ac29601a82bfddeb4cf9"
RELATED_FRESH_HOLDOUT_PR = 714

S_MINUS = 1.75
S_LIVE = 1.80
S_PLUS = 1.85
S_VALUES = (S_MINUS, S_LIVE, S_PLUS)
S_STEP = 0.05
Z_INNER = 1.0
Z_OUTER = 2.0
R_SUPPORT_MAX = 1.6
ALPHA_REPLAY_TOL = 2.0e-12

IDENTIFIABILITY_COSINE_MAX_ABS = 0.995
IDENTIFIABILITY_CONDITION_MAX = 25.0
DIVERGENCE_STEP = 1.0e-5
DIVERGENCE_MAX = 1.0e-5

RESPONSE_RADII = radial_shape.RESPONSE_RADII
RESPONSE_ABS_Z = radial_shape.RESPONSE_ABS_Z
RESPONSE_ANGLES = radial_shape.RESPONSE_ANGLES
RESPONSE_TIMES = radial_shape.RESPONSE_TIMES

FINGERPRINT_RADII = (0.60, 0.90, 1.20)
FINGERPRINT_ABS_Z = (1.55, 1.75, 1.80, 1.85, 1.95)
FINGERPRINT_TIMES = (0.375, 0.50, 0.625, 0.75)
CORE_RADII = (0.35, 0.65, 0.90)
ANNULUS_RADII = (1.20, 1.35, 1.50)
FINGERPRINT_ANGLES = tuple(float(v) for v in np.arange(8) * (0.25 * np.pi))
CROSSING_RADII = (0.60, 0.90, 1.20)
CROSSING_Z = tuple(float(v) for v in np.linspace(1.05, 1.98, 187))

TRUTH = {
    "candidate_velocity_changed": False,
    "canonical_velocity_changed": False,
    "saved_velocity_changed": False,
    "basis_dimension_changed": False,
    "new_spatial_basis_added": False,
    "new_temporal_basis_added": False,
    "candidate_turnover_selected": False,
    "parameter_scan_performed": False,
    "optimization_performed": False,
    "post_result_retuning_performed": False,
    "historical_development_paths_used_only_for_frozen_alpha_semantics": True,
    "fresh_714_path_data_used": False,
    "fresh_714_path_data_used_for_retuning": False,
    "pressure_or_force_changed": False,
    "held_out_pde_residual_evaluated": False,
    "st006_comparison_performed": False,
    "public_image_numeric_target_used": False,
    "public_numeric_turnover_target_defined": False,
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
    replay._assert_source_lock()
    radial_shape._assert_source_lock()
    if reservoir.TASK_ID != "CR003-ST052M-OUTER-AXIAL-RETURN-RESERVOIR-PREFLIGHT-108":
        raise RuntimeError("#767 reservoir task identity drifted")
    if replay.TASK_ID != "CR003-ST052M-OUTER-RESERVOIR-NONLINEAR-REPLAY-110":
        raise RuntimeError("#775 child task identity drifted")
    if radial_shape.TASK_ID != "CR003-ST052M-TOTAL-SHAPE-COORDINATE-SENSITIVITY-115":
        raise RuntimeError("#810 radial-shape task identity drifted")
    if reservoir.Z_INNER != Z_INNER or reservoir.Z_GLOBAL_OUTER != Z_OUTER:
        raise RuntimeError("outer-reservoir global support drifted")
    if reservoir.Z_VISIBLE_OUTER != S_LIVE:
        raise RuntimeError("live #767 turnover is no longer 1.8")
    if reservoir.R_SUPPORT_MAX != R_SUPPORT_MAX:
        raise RuntimeError("radial support drifted")
    if RELATED_FRESH_HOLDOUT_PR != 714:
        raise RuntimeError("fresh holdout identity drifted")


def axial_potential_and_derivative_s(z: np.ndarray, turnover: float) -> tuple[np.ndarray, np.ndarray]:
    """Odd C3 Z_s and even dZ_s/dz for a frozen turnover in (1,2)."""
    s_turn = float(turnover)
    if not (Z_INNER < s_turn < Z_OUTER):
        raise ValueError("turnover must lie strictly between 1 and 2")
    zz = np.asarray(z, dtype=float)
    s = np.abs(zz)
    h = np.zeros_like(s)
    hp = np.zeros_like(s)

    visible = (s > Z_INNER) & (s < s_turn)
    if np.any(visible):
        xi = (s[visible] - Z_INNER) / (s_turn - Z_INNER)
        h[visible] = reservoir.septic_smoothstep(xi)
        hp[visible] = reservoir.septic_smoothstep_prime(xi) / (s_turn - Z_INNER)

    outer = (s > s_turn) & (s < Z_OUTER)
    if np.any(outer):
        xi = (s[outer] - s_turn) / (Z_OUTER - s_turn)
        h[outer] = 1.0 - reservoir.septic_smoothstep(xi)
        hp[outer] = -reservoir.septic_smoothstep_prime(xi) / (Z_OUTER - s_turn)

    join = s == s_turn
    h[join] = 1.0
    hp[join] = 0.0

    return np.sign(zz) * h, hp


def reservoir_correction_s(points: np.ndarray, turnover: float) -> np.ndarray:
    """Exact #767 curl formula with only the axial turnover changed."""
    pts = np.asarray(points, dtype=float)
    if pts.ndim != 2 or pts.shape[1] != 3:
        raise ValueError("points must have shape (n,3)")
    x, y, z = pts[:, 0], pts[:, 1], pts[:, 2]
    q = x * x + y * y
    radial, radial_q = reservoir._radial_profile_and_derivative(q)
    axial, axial_z = axial_potential_and_derivative_s(z, float(turnover))
    f = radial * axial
    f_q = radial_q * axial
    f_z = radial * axial_z
    correction = np.column_stack((-x * f_z, -y * f_z, 2.0 * f + 2.0 * q * f_q))
    if not np.all(np.isfinite(correction)):
        raise RuntimeError("nonfinite axial-turnover correction")
    return correction


def unit_time_correction_s(points: np.ndarray, time: float, turnover: float) -> np.ndarray:
    return replay.screen.activation(float(time)) * reservoir_correction_s(points, float(turnover))


def candidate_velocity_s(
    parent_fn: Callable[[np.ndarray, float], np.ndarray],
    points: np.ndarray,
    time: float,
    turnover: float,
    alpha: float,
) -> np.ndarray:
    pts = np.asarray(points, dtype=float)
    return np.asarray(parent_fn(pts, float(time)), dtype=float) + float(alpha) * unit_time_correction_s(
        pts, float(time), float(turnover)
    )


def derive_alpha_s(
    parent_fn: Callable,
    parent_times: np.ndarray,
    parent_positions: np.ndarray,
    metadata: list[dict[str, Any]],
    turnover: float,
) -> tuple[float, dict[str, Any]]:
    """Apply #775's exact historical 24-tip scalar-coordinate rule to Z_s."""
    if abs(float(parent_times[replay.screen.MID_INDEX]) - replay.screen.MID_TIME) > 1.0e-15:
        raise RuntimeError("historical midpoint index drifted")
    tip_idx = np.asarray([i for i, m in enumerate(metadata) if str(m["band"]) == "tip"], dtype=int)
    if len(tip_idx) != 24:
        raise RuntimeError("historical tip path count drifted")
    pts = np.asarray(parent_positions[replay.screen.MID_INDEX, tip_idx], dtype=float)
    base_ur = replay.screen.loc._radial_velocity(parent_fn, pts, replay.screen.MID_TIME)
    unit_fn = lambda p, t: unit_time_correction_s(p, t, float(turnover))
    unit_ur = replay.screen.loc._radial_velocity(unit_fn, pts, replay.screen.MID_TIME)
    mean_base = float(np.mean(base_ur))
    mean_unit = float(np.mean(unit_ur))
    if not (np.isfinite(mean_base) and np.isfinite(mean_unit)):
        raise RuntimeError("nonfinite axial-turnover calibration")
    if abs(mean_unit) <= np.finfo(float).tiny:
        raise RuntimeError("degenerate axial-turnover unit response on historical tips")
    alpha = float(-mean_base / mean_unit)
    if not np.isfinite(alpha) or abs(alpha) <= np.finfo(float).tiny:
        raise RuntimeError("invalid axial-turnover alpha")
    return alpha, {
        "turnover": float(turnover),
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


def build_local_families() -> tuple[
    Callable,
    dict[float, float],
    dict[float, dict[str, Any]],
    dict[float, Callable],
    dict[float, float],
    dict[float, Callable],
]:
    """Build frozen axial-turnover and radial-shape total families on one parent state."""
    _assert_source_lock()
    parent_fn, parent_times, parent_positions, metadata = radial_shape._historical_parent_state()

    s_alphas: dict[float, float] = {}
    s_calibrations: dict[float, dict[str, Any]] = {}
    s_children: dict[float, Callable] = {}
    for turnover in S_VALUES:
        alpha, calibration = derive_alpha_s(
            parent_fn, parent_times, parent_positions, metadata, float(turnover)
        )
        s_alphas[float(turnover)] = float(alpha)
        s_calibrations[float(turnover)] = calibration
        s_children[float(turnover)] = (
            lambda pts, time, ss=float(turnover), aa=float(alpha): candidate_velocity_s(
                parent_fn, pts, time, ss, aa
            )
        )

    live_alpha, live_calibration = replay.derive_reservoir_alpha(
        parent_fn, parent_times, parent_positions, metadata
    )
    live_rel = abs(s_alphas[S_LIVE] - float(live_alpha)) / max(
        abs(float(live_alpha)), np.finfo(float).tiny
    )
    if live_rel > ALPHA_REPLAY_TOL:
        raise RuntimeError(f"s=1.80 alpha does not replay #775: relative mismatch={live_rel}")
    s_calibrations[S_LIVE]["exact_775_alpha"] = float(live_alpha)
    s_calibrations[S_LIVE]["relative_mismatch_vs_exact_775"] = float(live_rel)
    s_calibrations[S_LIVE]["exact_775_calibration"] = live_calibration

    k_alphas: dict[float, float] = {}
    k_children: dict[float, Callable] = {}
    for k in (radial_shape.K_MINUS, radial_shape.K_LIVE, radial_shape.K_PLUS):
        alpha_k, _ = radial_shape.derive_alpha_k(
            parent_fn, parent_times, parent_positions, metadata, float(k)
        )
        k_alphas[float(k)] = float(alpha_k)
        k_children[float(k)] = (
            lambda pts, time, kk=float(k), aa=float(alpha_k): radial_shape.candidate_velocity_k(
                parent_fn, pts, time, kk, aa
            )
        )

    return parent_fn, s_alphas, s_calibrations, s_children, k_alphas, k_children


def _normalize_column(v: np.ndarray, name: str) -> tuple[np.ndarray, float]:
    vec = np.asarray(v, dtype=float).ravel()
    norm = float(np.linalg.norm(vec))
    if not np.isfinite(norm) or norm <= np.finfo(float).tiny:
        raise RuntimeError(f"degenerate tangent column: {name}")
    return vec / norm, norm


def three_coordinate_identifiability(
    parent_fn: Callable,
    s_children: dict[float, Callable],
    k_children: dict[float, Callable],
) -> dict[str, Any]:
    """Amplitude vs radial-shape vs axial-turnover total-family Jacobian."""
    pts = radial_shape._response_points()
    amp_blocks: list[np.ndarray] = []
    radial_blocks: list[np.ndarray] = []
    axial_blocks: list[np.ndarray] = []
    for time in RESPONSE_TIMES:
        parent = np.asarray(parent_fn(pts, time), dtype=float)
        live = np.asarray(s_children[S_LIVE](pts, time), dtype=float)
        k_minus = np.asarray(k_children[radial_shape.K_MINUS](pts, time), dtype=float)
        k_plus = np.asarray(k_children[radial_shape.K_PLUS](pts, time), dtype=float)
        s_minus = np.asarray(s_children[S_MINUS](pts, time), dtype=float)
        s_plus = np.asarray(s_children[S_PLUS](pts, time), dtype=float)
        amp_blocks.append((live - parent).ravel())
        radial_blocks.append(
            ((k_plus - k_minus) / (radial_shape.K_PLUS - radial_shape.K_MINUS)).ravel()
        )
        axial_blocks.append(((s_plus - s_minus) / (S_PLUS - S_MINUS)).ravel())

    amp, n_amp = _normalize_column(np.concatenate(amp_blocks), "amplitude")
    radial, n_radial = _normalize_column(np.concatenate(radial_blocks), "radial shape")
    axial, n_axial = _normalize_column(np.concatenate(axial_blocks), "axial turnover")
    matrix = np.column_stack((amp, radial, axial))
    singular = np.linalg.svd(matrix, compute_uv=False)
    rank = int(np.linalg.matrix_rank(matrix, tol=1.0e-10))
    cosines = {
        "amplitude_radial_shape": float(np.dot(amp, radial)),
        "amplitude_axial_turnover": float(np.dot(amp, axial)),
        "radial_shape_axial_turnover": float(np.dot(radial, axial)),
    }
    max_abs_cosine = float(max(abs(v) for v in cosines.values()))
    condition = float(singular[0] / singular[-1])
    passes = bool(
        rank == 3
        and max_abs_cosine < IDENTIFIABILITY_COSINE_MAX_ABS
        and condition <= IDENTIFIABILITY_CONDITION_MAX
    )
    return {
        "response_point_count_per_time": int(len(pts)),
        "times": list(RESPONSE_TIMES),
        "column_definitions": {
            "amplitude": "u_s1.80-u_parent",
            "radial_shape": "(u_k1.10-u_k0.90)/0.20 with independent #775 alpha semantics",
            "axial_turnover": "(u_s1.85-u_s1.75)/0.10 with independent #775 alpha semantics",
        },
        "raw_tangent_norms": {
            "amplitude": n_amp,
            "radial_shape_per_k": n_radial,
            "axial_turnover_per_s": n_axial,
        },
        "normalized_column_rank": rank,
        "pairwise_cosines": cosines,
        "max_abs_pairwise_cosine": max_abs_cosine,
        "normalized_singular_values": singular.tolist(),
        "normalized_condition": condition,
        "cosine_gate_abs_max": IDENTIFIABILITY_COSINE_MAX_ABS,
        "condition_gate_max": IDENTIFIABILITY_CONDITION_MAX,
        "passes": passes,
    }


def turnover_sign_witness() -> dict[str, Any]:
    """Freeze the direct geometric meaning of s at |z|=1.80."""
    point = np.asarray([[1.0, 0.0, 1.80], [1.0, 0.0, -1.80]], dtype=float)
    values: dict[str, Any] = {}
    means: dict[float, float] = {}
    for turnover in S_VALUES:
        corr = reservoir_correction_s(point, turnover)
        ur = replay.screen.loc._radial_velocity(
            lambda p, _t, ss=float(turnover): reservoir_correction_s(p, ss), point, 0.50
        )
        mean = float(np.mean(ur))
        means[float(turnover)] = mean
        values[f"s={turnover:.2f}"] = {
            "radial_components": np.asarray(ur, dtype=float).tolist(),
            "mean_radial_component": mean,
            "correction_vectors": corr.tolist(),
        }
    tol = 1.0e-13
    passes = bool(means[S_MINUS] > 0.0 and abs(means[S_LIVE]) <= tol and means[S_PLUS] < 0.0)
    return {
        "probe_abs_z": 1.80,
        "probe_radius": 1.0,
        "by_turnover": values,
        "expected_semantics": "s=1.75 outward return; s=1.80 turnover; s=1.85 still inward",
        "live_zero_tolerance": tol,
        "passes": passes,
    }


def _ring_points(radii: tuple[float, ...], abs_z: float) -> np.ndarray:
    rows: list[list[float]] = []
    for r in radii:
        for z_sign in (-1.0, 1.0):
            for angle in FINGERPRINT_ANGLES:
                rows.append([r * math.cos(angle), r * math.sin(angle), z_sign * float(abs_z)])
    return np.asarray(rows, dtype=float)


def total_radial_fingerprints(children: dict[float, Callable]) -> dict[str, Any]:
    by_s: dict[str, Any] = {}
    for turnover in S_VALUES:
        by_z: dict[str, Any] = {}
        for abs_z in FINGERPRINT_ABS_Z:
            pts = _ring_points(FINGERPRINT_RADII, abs_z)
            by_time: dict[str, Any] = {}
            for time in FINGERPRINT_TIMES:
                ur = np.asarray(
                    replay.screen.loc._radial_velocity(children[float(turnover)], pts, float(time)),
                    dtype=float,
                )
                by_time[f"{time:.3f}"] = {
                    "mean_total_radial_velocity": float(np.mean(ur)),
                    "rms_total_radial_velocity": float(np.sqrt(np.mean(np.square(ur)))),
                    "inward_fraction": float(np.mean(ur < 0.0)),
                    "probe_count": int(len(ur)),
                }
            by_z[f"|z|={abs_z:.2f}"] = by_time
        by_s[f"s={turnover:.2f}"] = by_z
    return {
        "radii": list(FINGERPRINT_RADII),
        "abs_z": list(FINGERPRINT_ABS_Z),
        "angles": list(FINGERPRINT_ANGLES),
        "times": list(FINGERPRINT_TIMES),
        "by_turnover": by_s,
        "acceptance_gate_applied": False,
        "reason": "no public numeric axial-turnover or tip-length target exists",
    }


def total_axial_fingerprints(children: dict[float, Callable]) -> dict[str, Any]:
    by_s: dict[str, Any] = {}
    for turnover in S_VALUES:
        by_z: dict[str, Any] = {}
        for abs_z in FINGERPRINT_ABS_Z:
            core_pts = _ring_points(CORE_RADII, abs_z)
            ann_pts = _ring_points(ANNULUS_RADII, abs_z)
            by_time: dict[str, Any] = {}
            for time in FINGERPRINT_TIMES:
                core_vel = np.asarray(children[float(turnover)](core_pts, float(time)), dtype=float)
                ann_vel = np.asarray(children[float(turnover)](ann_pts, float(time)), dtype=float)
                core_away = np.sign(core_pts[:, 2]) * core_vel[:, 2]
                ann_away = np.sign(ann_pts[:, 2]) * ann_vel[:, 2]
                by_time[f"{time:.3f}"] = {
                    "core_mean_away_from_midplane_uz": float(np.mean(core_away)),
                    "core_positive_fraction": float(np.mean(core_away > 0.0)),
                    "annulus_mean_away_from_midplane_uz": float(np.mean(ann_away)),
                    "annulus_positive_fraction": float(np.mean(ann_away > 0.0)),
                }
            by_z[f"|z|={abs_z:.2f}"] = by_time
        by_s[f"s={turnover:.2f}"] = by_z
    return {
        "core_radii": list(CORE_RADII),
        "annulus_radii": list(ANNULUS_RADII),
        "abs_z": list(FINGERPRINT_ABS_Z),
        "times": list(FINGERPRINT_TIMES),
        "by_turnover": by_s,
        "acceptance_gate_applied": False,
    }


def _mean_ring_ur(child_fn: Callable, radius: float, abs_z: float, time: float) -> float:
    pts = _ring_points((float(radius),), float(abs_z))
    ur = replay.screen.loc._radial_velocity(child_fn, pts, float(time))
    return float(np.mean(np.asarray(ur, dtype=float)))


def midpoint_total_radial_crossings(children: dict[float, Callable]) -> dict[str, Any]:
    """Observation-only total-field radial sign crossings on a deterministic z scan."""
    by_s: dict[str, Any] = {}
    z_grid = np.asarray(CROSSING_Z, dtype=float)
    for turnover in S_VALUES:
        rows: dict[str, Any] = {}
        for radius in CROSSING_RADII:
            values = np.asarray(
                [_mean_ring_ur(children[float(turnover)], radius, z, 0.50) for z in z_grid],
                dtype=float,
            )
            crossing = None
            for i in range(len(z_grid) - 1):
                a, b = float(values[i]), float(values[i + 1])
                if a < 0.0 <= b:
                    if b == a:
                        crossing = float(z_grid[i])
                    else:
                        frac = -a / (b - a)
                        crossing = float(z_grid[i] + frac * (z_grid[i + 1] - z_grid[i]))
                    break
            rows[f"r={radius:.2f}"] = {
                "first_inward_to_noninward_crossing_abs_z": crossing,
                "min_mean_total_ur": float(np.min(values)),
                "max_mean_total_ur": float(np.max(values)),
                "starts_inward": bool(values[0] < 0.0),
            }
        by_s[f"s={turnover:.2f}"] = rows
    return {
        "time": 0.50,
        "z_scan": {"start": float(z_grid[0]), "stop": float(z_grid[-1]), "count": int(len(z_grid))},
        "radii": list(CROSSING_RADII),
        "by_turnover": by_s,
        "acceptance_gate_applied": False,
    }


def structural_guards(parent_fn: Callable, children: dict[float, Callable]) -> dict[str, Any]:
    support_points = np.asarray(
        [
            [1.61, 0.0, 1.50],
            [-1.61, 0.0, -1.50],
            [0.90, 0.0, 0.95],
            [0.90, 0.0, -0.95],
            [0.90, 0.0, 2.00],
            [0.90, 0.0, -2.00],
            [0.90, 0.0, 2.05],
            [0.90, 0.0, -2.05],
        ],
        dtype=float,
    )
    div_pts = []
    for r in (0.35, 0.90, 1.35):
        for abs_z in (1.30, 1.60, 1.92):
            for z_sign in (-1.0, 1.0):
                for angle in FINGERPRINT_ANGLES:
                    div_pts.append([r * math.cos(angle), r * math.sin(angle), z_sign * abs_z])
    div_pts = np.asarray(div_pts, dtype=float)

    rng = np.random.default_rng(9173611)
    start_pts = rng.uniform(-1.9, 1.9, size=(96, 3))
    by_s: dict[str, Any] = {}
    all_pass = True
    for turnover in S_VALUES:
        correction_support_max = float(
            np.max(np.abs(reservoir_correction_s(support_points, float(turnover))))
        )
        start_diff = float(
            np.max(
                np.abs(
                    np.asarray(children[float(turnover)](start_pts, 0.25), dtype=float)
                    - np.asarray(parent_fn(start_pts, 0.25), dtype=float)
                )
            )
        )
        div = np.zeros(len(div_pts), dtype=float)
        for axis in range(3):
            shift = np.zeros(3, dtype=float)
            shift[axis] = DIVERGENCE_STEP
            up = reservoir_correction_s(div_pts + shift, float(turnover))
            um = reservoir_correction_s(div_pts - shift, float(turnover))
            div += (up[:, axis] - um[:, axis]) / (2.0 * DIVERGENCE_STEP)
        div_max = float(np.max(np.abs(div)))
        nontrivial = float(
            np.sqrt(np.mean(np.square(reservoir_correction_s(div_pts, float(turnover)))))
        )
        passed = bool(
            correction_support_max == 0.0
            and start_diff <= 1.0e-12
            and div_max <= DIVERGENCE_MAX
            and nontrivial > 1.0e-8
        )
        all_pass = all_pass and passed
        by_s[f"s={turnover:.2f}"] = {
            "outside_support_correction_max_abs": correction_support_max,
            "start_identity_max_abs": start_diff,
            "unit_correction_fd_divergence_max_abs": div_max,
            "unit_correction_rms_on_divergence_cloud": nontrivial,
            "passes": passed,
        }
    return {
        "divergence_step": DIVERGENCE_STEP,
        "divergence_gate_max": DIVERGENCE_MAX,
        "by_turnover": by_s,
        "all_pass": bool(all_pass),
    }


def exact_live_shape_replay() -> dict[str, Any]:
    """Implementation-level check that s=1.80 reproduces #767 exactly."""
    z = np.asarray([-2.05, -2.0, -1.95, -1.85, -1.80, -1.75, -1.55, -1.0, 0.0,
                    1.0, 1.55, 1.75, 1.80, 1.85, 1.95, 2.0, 2.05], dtype=float)
    a0, a1 = reservoir.axial_potential_and_derivative(z)
    b0, b1 = axial_potential_and_derivative_s(z, S_LIVE)
    axial_max = float(max(np.max(np.abs(a0 - b0)), np.max(np.abs(a1 - b1))))

    pts = []
    for r in (0.0, 0.35, 0.90, 1.35, 1.60, 1.65):
        for zz in z:
            pts.append([r, 0.0, float(zz)])
    pts = np.asarray(pts, dtype=float)
    live = reservoir.reservoir_correction(pts)
    replayed = reservoir_correction_s(pts, S_LIVE)
    correction_max = float(np.max(np.abs(live - replayed)))
    passes = bool(axial_max <= 1.0e-14 and correction_max <= 1.0e-14)
    return {
        "axial_potential_derivative_max_abs_difference": axial_max,
        "unit_correction_max_abs_difference": correction_max,
        "gate": 1.0e-14,
        "passes": passes,
    }


def build_report() -> dict[str, Any]:
    _assert_source_lock()
    live_replay = exact_live_shape_replay()
    parent_fn, s_alphas, s_calibrations, s_children, k_alphas, k_children = build_local_families()
    ident = three_coordinate_identifiability(parent_fn, s_children, k_children)
    sign = turnover_sign_witness()
    structure = structural_guards(parent_fn, s_children)
    radial = total_radial_fingerprints(s_children)
    axial = total_axial_fingerprints(s_children)
    crossings = midpoint_total_radial_crossings(s_children)

    calibration_rows: dict[str, Any] = {}
    live_alpha = float(s_alphas[S_LIVE])
    for turnover in S_VALUES:
        row = dict(s_calibrations[float(turnover)])
        row["alpha_over_live"] = float(s_alphas[float(turnover)] / live_alpha)
        calibration_rows[f"s={turnover:.2f}"] = row

    passed = bool(live_replay["passes"] and ident["passes"] and sign["passes"] and structure["all_pass"])
    return {
        "task_id": TASK_ID,
        "prereg_issue": PREREG_ISSUE,
        "source_parent": {"pr": SOURCE_PARENT_PR, "head": SOURCE_PARENT_HEAD},
        "source_child": {"pr": SOURCE_CHILD_PR, "head": SOURCE_CHILD_HEAD},
        "source_radial_shape": {"pr": SOURCE_RADIAL_SHAPE_PR, "head": SOURCE_RADIAL_SHAPE_HEAD},
        "frozen_family": {
            "turnovers": list(S_VALUES),
            "live_turnover": S_LIVE,
            "axial_support": [Z_INNER, Z_OUTER],
            "radial_support_max": R_SUPPORT_MAX,
            "time_law": "existing g1(t)=2*(t-.25)",
            "basis_dimension_increment": 0,
            "selected_candidate": None,
        },
        "live_shape_exact_replay": live_replay,
        "calibration": calibration_rows,
        "radial_shape_alphas_for_identifiability": {f"k={k:.2f}": float(v) for k, v in k_alphas.items()},
        "turnover_sign_witness": sign,
        "three_coordinate_identifiability": ident,
        "structural_guards": structure,
        "visualization_fingerprints": {
            "total_radial": radial,
            "total_axial": axial,
            "midpoint_total_radial_crossings": crossings,
            "numeric_source_target_used": False,
            "improvement_direction_declared": False,
        },
        "held_out_pde_residual": {
            "evaluated": False,
            "st006_comparison_performed": False,
            "reason": "capacity preflight selects no new child and rebuilds no matched pressure/restricted forcing",
        },
        **TRUTH,
        "axial_turnover_shape_preflight_passed": passed,
        "routing_if_pass": (
            "Treat tip axial extent/return-flow location as a same-dimension shape lever. "
            "If later fixed source-observable/render evidence says the tip is too short, a preregistered "
            "s>1.8 child is the minimal direction; if too long/return flow starts too late, test s<1.8. "
            "Do not add a second poloidal basis solely for axial extent before that child test."
        ),
        "routing_if_fail": (
            "Do not retune this family post hoc. Treat turnover shape as non-identifiable under the frozen "
            "semantics and move to one genuinely independent poloidal spatial channel only if a diagnosed "
            "morphology mismatch requires it."
        ),
        "scope_statement": (
            "PASS would establish only same-dimension local morphology capacity and conditioning. It is not "
            "visual/source correspondence, PDE validation, paper exactness, exact OpenAI-field identity, or blow-up evidence."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args()
    report = build_report()
    payload = json.dumps(report, indent=2, sort_keys=True)
    if args.out is not None:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(payload + "\n", encoding="utf-8")
    print(payload)


if __name__ == "__main__":
    main()
