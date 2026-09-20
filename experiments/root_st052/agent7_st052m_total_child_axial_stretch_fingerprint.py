"""Audit total-child axial-stretch morphology for the frozen outer-reservoir child.

Issue #793 freezes this one-increment Agent-7 diagnostic before observing any
new total-child axial result.  The velocity formula is not changed here.  We
reconstruct exact PR #775's frozen child using its historical-development-only
coefficient rule, then compare total parent and total child axial motion on a
fixed radial/axial/time probe family.

The diagnostic is autonomous project evidence for qualitative axial stretching
and annular compensation.  It is not a public OpenAI numerical target, a pixel
objective, PDE validation, or exact-field identification.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Callable

import numpy as np

import agent7_st052m_outer_reservoir_nonlinear_replay as replay

TASK_ID = "CR003-ST052M-TOTAL-CHILD-AXIAL-STRETCH-FINGERPRINT-113"
PREREG_ISSUE = 793
SOURCE_CHILD_PR = 775
SOURCE_CHILD_HEAD = "93887a59729d22113badf2ae4dac2f7d868e2703"
SOURCE_UNIT_AUDIT_PR = 785
SOURCE_SCOPE_PR = 786
SOURCE_SCOPE_HEAD = "f3c1e5a2a86bcc09c9d1f16e9f865c5d57a94ae3"
RELATED_FRESH_HOLDOUT_PR = 714

ABS_Z_BANDS = (1.55, 1.75)
PROBE_TIMES = (0.375, 0.50, 0.625, 0.75)
RADII = np.linspace(0.20, 1.55, 55)
ANGLES = np.arange(8, dtype=float) * (0.25 * np.pi)
CORE_RADIUS_MAX = 0.90
ANNULUS_RADIUS_MIN = 1.20
ANNULUS_RADIUS_MAX = 1.50
CROSSING_RADIUS_MIN = 0.90
CROSSING_RADIUS_MAX = 1.50

TRUTH = {
    "velocity_formula_changed_this_increment": False,
    "candidate_coefficient_changed_this_increment": False,
    "basis_dimension_changed": False,
    "new_spatial_basis_added": False,
    "new_temporal_basis_added": False,
    "frozen_775_total_child_reconstructed": True,
    "historical_development_paths_used_only_to_rederive_frozen_775_alpha": True,
    "fresh_714_path_data_used": False,
    "fresh_714_path_data_used_for_retuning": False,
    "parameter_scan_performed": False,
    "optimization_performed": False,
    "post_result_retuning_performed": False,
    "pressure_or_force_changed": False,
    "held_out_pde_residual_evaluated": False,
    "public_image_numeric_target_used": False,
    "pixel_similarity_objective_used": False,
    "canonical_velocity_changed": False,
    "saved_velocity_changed": False,
    "production_candidate_selected": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "source_correspondence_verified": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _assert_source_lock() -> None:
    replay._assert_source_lock()
    if replay.TASK_ID != "CR003-ST052M-OUTER-RESERVOIR-NONLINEAR-REPLAY-110":
        raise RuntimeError("#775 replay task identity drifted")
    if replay.PREREG_ISSUE != 774:
        raise RuntimeError("#775 preregistration identity drifted")
    if tuple(replay.VISIBLE_ABS_Z) != ABS_Z_BANDS:
        raise RuntimeError("#775 visible-z probe semantics drifted")
    if replay.RELATED_FRESH_HOLDOUT_PR != RELATED_FRESH_HOLDOUT_PR:
        raise RuntimeError("fresh holdout identity drifted")
    if SOURCE_CHILD_HEAD != "93887a59729d22113badf2ae4dac2f7d868e2703":
        raise RuntimeError("frozen #775 exact-head identity drifted")
    if SOURCE_SCOPE_HEAD != "f3c1e5a2a86bcc09c9d1f16e9f865c5d57a94ae3":
        raise RuntimeError("#786 scope-head identity drifted")


def _build_frozen_child() -> tuple[Callable, Callable, float, dict[str, Any]]:
    """Reconstruct #775 exactly without consuming any fresh #714 information."""
    _assert_source_lock()
    _control_fn, _linear_fn, parent_fn, _epsilon, _normalization = replay.screen._build_fields()
    parent_times, parent_positions, metadata = replay.screen._integrate(parent_fn)
    alpha, calibration = replay.derive_reservoir_alpha(
        parent_fn, parent_times, parent_positions, metadata
    )
    child_fn = lambda pts, t: replay.reservoir_velocity(parent_fn, pts, t, alpha)
    return parent_fn, child_fn, float(alpha), calibration


def _ring_cloud(abs_z: float) -> tuple[np.ndarray, np.ndarray]:
    """Return all frozen (radius, +/-z, 8-angle) probes and their z signs."""
    rows: list[list[float]] = []
    signs: list[float] = []
    for radius in RADII:
        for z_sign in (-1.0, 1.0):
            z = z_sign * float(abs_z)
            for angle in ANGLES:
                rows.append(
                    [
                        float(radius * np.cos(angle)),
                        float(radius * np.sin(angle)),
                        z,
                    ]
                )
                signs.append(z_sign)
    points = np.asarray(rows, dtype=float)
    sign_array = np.asarray(signs, dtype=float)
    expected = len(RADII) * 2 * len(ANGLES)
    if points.shape != (expected, 3) or sign_array.shape != (expected,):
        raise RuntimeError("frozen total-child axial probe cloud drifted")
    return points, sign_array


def axial_stretch_profile(field_fn: Callable, abs_z: float, time: float) -> np.ndarray:
    """Average sign(z)*w over both z signs and eight azimuths at each radius."""
    points, z_sign = _ring_cloud(abs_z)
    velocity = np.asarray(field_fn(points, float(time)), dtype=float)
    if velocity.shape != points.shape or not np.all(np.isfinite(velocity)):
        raise RuntimeError("nonfinite or malformed total-field velocity on axial fingerprint cloud")
    stretch = z_sign * velocity[:, 2]
    return stretch.reshape(len(RADII), 2, len(ANGLES)).mean(axis=(1, 2))


def first_positive_to_nonpositive_crossing(
    radii: np.ndarray, profile: np.ndarray
) -> float | None:
    """Locate the first crossing after a positive central prefix, with linear interpolation."""
    r = np.asarray(radii, dtype=float)
    y = np.asarray(profile, dtype=float)
    if r.ndim != 1 or y.shape != r.shape or len(r) < 2:
        raise ValueError("radii/profile must be equal-length one-dimensional arrays")
    if not np.all(np.isfinite(r)) or not np.all(np.isfinite(y)):
        raise ValueError("radii/profile must be finite")
    if not np.all(np.diff(r) > 0.0):
        raise ValueError("radii must be strictly increasing")
    if y[0] <= 0.0:
        return None
    for index in range(1, len(r)):
        if y[index] <= 0.0:
            y0 = float(y[index - 1])
            y1 = float(y[index])
            r0 = float(r[index - 1])
            r1 = float(r[index])
            if y1 == 0.0 or y0 == y1:
                return r1
            weight = y0 / (y0 - y1)
            return float(r0 + weight * (r1 - r0))
    return None


def _contiguous_central_positive_fraction(profile: np.ndarray) -> float:
    y = np.asarray(profile, dtype=float)
    count = 0
    for value in y:
        if value > 0.0:
            count += 1
        else:
            break
    return float(count / len(y))


def _one_fingerprint(
    parent_fn: Callable,
    child_fn: Callable,
    abs_z: float,
    time: float,
) -> dict[str, Any]:
    parent = axial_stretch_profile(parent_fn, abs_z, time)
    child = axial_stretch_profile(child_fn, abs_z, time)
    delta = child - parent

    core_mask = RADII <= CORE_RADIUS_MAX
    annulus_mask = (RADII >= ANNULUS_RADIUS_MIN) & (RADII <= ANNULUS_RADIUS_MAX)
    if not np.any(core_mask) or not np.any(annulus_mask):
        raise RuntimeError("frozen radial masks are empty")

    core_delta_mean = float(np.mean(delta[core_mask]))
    annulus_delta_mean = float(np.mean(delta[annulus_mask]))
    crossing = first_positive_to_nonpositive_crossing(RADII, delta)
    central_positive = bool(np.any(delta[core_mask] > 0.0))
    outer_nonpositive = bool(np.any(delta[annulus_mask] <= 0.0))
    crossing_in_window = bool(
        crossing is not None and CROSSING_RADIUS_MIN < crossing < CROSSING_RADIUS_MAX
    )
    gate = bool(
        core_delta_mean > 0.0
        and annulus_delta_mean < 0.0
        and central_positive
        and outer_nonpositive
        and crossing_in_window
    )

    return {
        "abs_z": float(abs_z),
        "time": float(time),
        "core_delta_mean_away_from_midplane_w": core_delta_mean,
        "annulus_delta_mean_away_from_midplane_w": annulus_delta_mean,
        "parent_core_mean_away_from_midplane_w": float(np.mean(parent[core_mask])),
        "child_core_mean_away_from_midplane_w": float(np.mean(child[core_mask])),
        "parent_annulus_mean_away_from_midplane_w": float(np.mean(parent[annulus_mask])),
        "child_annulus_mean_away_from_midplane_w": float(np.mean(child[annulus_mask])),
        "first_positive_to_nonpositive_delta_crossing_radius": crossing,
        "contiguous_central_positive_radius_fraction": _contiguous_central_positive_fraction(delta),
        "parent_total_positive_stretch_radius_fraction": float(np.mean(parent > 0.0)),
        "child_total_positive_stretch_radius_fraction": float(np.mean(child > 0.0)),
        "delta_has_central_positive_sample": central_positive,
        "delta_has_outer_nonpositive_sample": outer_nonpositive,
        "crossing_inside_preregistered_window": crossing_in_window,
        "passes_directional_gate": gate,
        "radii": RADII.tolist(),
        "parent_axial_stretch_profile": parent.tolist(),
        "child_axial_stretch_profile": child.tolist(),
        "child_minus_parent_axial_stretch_profile": delta.tolist(),
    }


def run(out: Path) -> dict[str, Any]:
    parent_fn, child_fn, alpha, calibration = _build_frozen_child()

    by_abs_z: dict[str, Any] = {}
    all_pass = True
    crossing_values: list[float] = []
    for abs_z in ABS_Z_BANDS:
        rows: dict[str, Any] = {}
        band_pass = True
        for time in PROBE_TIMES:
            row = _one_fingerprint(parent_fn, child_fn, abs_z, time)
            rows[f"{time:.3f}"] = row
            band_pass = band_pass and bool(row["passes_directional_gate"])
            crossing = row["first_positive_to_nonpositive_delta_crossing_radius"]
            if crossing is not None:
                crossing_values.append(float(crossing))
        by_abs_z[f"{abs_z:.2f}"] = {
            "by_time": rows,
            "all_times_pass_directional_gate": bool(band_pass),
        }
        all_pass = all_pass and band_pass

    report = {
        "task_id": TASK_ID,
        "prereg_issue": PREREG_ISSUE,
        "source_child": {"pr": SOURCE_CHILD_PR, "head": SOURCE_CHILD_HEAD},
        "source_unit_capacity_audit_pr": SOURCE_UNIT_AUDIT_PR,
        "source_scope_governance": {"pr": SOURCE_SCOPE_PR, "head": SOURCE_SCOPE_HEAD},
        "fresh_holdout_pr_not_consumed": RELATED_FRESH_HOLDOUT_PR,
        "frozen_child": {
            "formula": "u_#652 + alpha_res * 2*(t-.25) * C_res",
            "alpha_reservoir": alpha,
            "coefficient_calibration": calibration,
            "coefficient_rederived_by_exact_775_rule_without_retuning": True,
        },
        "probe_contract": {
            "abs_z_bands": list(ABS_Z_BANDS),
            "times": list(PROBE_TIMES),
            "radii": RADII.tolist(),
            "azimuth_count": int(len(ANGLES)),
            "z_sign_count": 2,
            "core_radius_max": CORE_RADIUS_MAX,
            "annulus_radius_interval": [ANNULUS_RADIUS_MIN, ANNULUS_RADIUS_MAX],
            "crossing_acceptance_interval": [CROSSING_RADIUS_MIN, CROSSING_RADIUS_MAX],
            "stretch_definition": "sign(z) * total_velocity_w",
            "comparison": "frozen #775 child minus exact #652 parent",
            "public_openai_numeric_target": None,
        },
        "by_abs_z": by_abs_z,
        "crossing_summary": {
            "count": int(len(crossing_values)),
            "minimum": min(crossing_values) if crossing_values else None,
            "maximum": max(crossing_values) if crossing_values else None,
        },
        "total_child_axial_stretch_fingerprint_passed": bool(all_pass),
        "decision": {
            "classification": (
                "frozen_775_total_child_axial_stretch_direction_pass"
                if all_pass
                else "frozen_775_total_child_axial_stretch_direction_fail"
            ),
            "closer_visualization_delivery_supported": bool(all_pass),
            "closer_visualization_delivery_scope": (
                "autonomous total-field axial-stretch direction relative to #652 on preregistered visible-tip probes only; not global/public visual correspondence"
            ),
            "velocity_formula_changed_this_increment": False,
            "second_poloidal_basis_justified_now": False,
            "if_fail_route": (
                "diagnose parent cancellation / same-dimension radial-envelope reshape before adding a second poloidal basis"
            ),
            "disjoint_fixed_render_or_fresh_path_still_required_for_broader_visual_claim": True,
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
        default=Path("outputs/agent7-st052m-total-child-axial-stretch-fingerprint/report.json"),
    )
    args = parser.parse_args()
    report = run(args.out)
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
