"""Replay one frozen toroidal morphology child on the #775 outer-reservoir field.

Preregistered in issue #845 before execution.  This is one minimal Agent-7
basis-growth increment: reconstruct the exact #775 outer-reservoir child, add
the compact pure-toroidal direction preflighted in #837, and freeze one scalar
coefficient by a 10% projection-strengthening rule on a deterministic visible-
tip cloud.

The increment is deliberately visualization/morphology-facing.  The toroidal
channel has cylindrical response (0, u_theta, 0), so this experiment asks
whether one actual +1-parameter child can strengthen swirl/circulation while
leaving the parent's radial-inward and axial-stretch components unchanged on
frozen probes.  It does not rebuild matched pressure/forcing and therefore does
not manufacture an incomparable held-out NS residual.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Callable

import numpy as np

import agent7_st052m_compact_toroidal_swirl_preflight as toroidal
import agent7_st052m_outer_reservoir_nonlinear_replay as reservoir_parent

TASK_ID = "CR003-ST052M-TOROIDAL-CHILD-MORPHOLOGY-REPLAY-120"
PREREG_ISSUE = 845
SOURCE_PARENT_PR = 837
SOURCE_PARENT_HEAD = "4c821c302cf89a7cc7a41730ce7029bdfbe6a7ef"
SOURCE_RESERVOIR_CHILD_PR = 775
SOURCE_RESERVOIR_CHILD_HEAD = "93887a59729d22113badf2ae4dac2f7d868e2703"
RELATED_FRESH_HOLDOUT_PR = 714

CALIBRATION_TIME = 0.50
PROBE_RADII = (0.55, 0.90, 1.20)
PROBE_ABS_Z = (1.55, 1.75)
PROBE_ANGLES = tuple(float(k) * math.pi / 4.0 for k in range(8))
PROBE_TIMES = (0.375, 0.50, 0.625, 0.75)
PROJECTION_STRENGTH = 0.10
ZERO_GATE = 1.0e-12
SWIRL_INCREASE_ABS_MIN = 1.0e-12
DIVERGENCE_STEP = 1.0e-5
DIVERGENCE_MAX = 1.0e-5
NONTRIVIAL_RELATIVE_CORRECTION_MIN = 1.0e-10

TRUTH = {
    "candidate_velocity_formula_evaluated": True,
    "canonical_velocity_changed": False,
    "saved_velocity_changed": False,
    "production_candidate_selected": False,
    "basis_dimension_increment": 1,
    "second_poloidal_basis_added": False,
    "toroidal_basis_added_to_experimental_child": True,
    "new_temporal_basis_added": False,
    "coefficient_scan_performed": False,
    "optimization_performed": False,
    "post_result_retuning_performed": False,
    "fresh_714_path_data_used": False,
    "pressure_or_force_changed": False,
    "held_out_pde_residual_evaluated": False,
    "parent_pde_receipt_transferred": False,
    "public_image_numeric_target_used": False,
    "source_numeric_swirl_target_used": False,
    "pixel_similarity_objective_used": False,
    "renderer_or_camera_used": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "source_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}

VelocityFn = Callable[[np.ndarray, float], np.ndarray]


def _assert_source_lock() -> None:
    toroidal._assert_source_lock()
    reservoir_parent._assert_source_lock()
    if toroidal.TASK_ID != "CR003-ST052M-COMPACT-TOROIDAL-SWIRL-PREFLIGHT-119":
        raise RuntimeError("#837 toroidal parent identity drifted")
    if reservoir_parent.TASK_ID != "CR003-ST052M-OUTER-RESERVOIR-NONLINEAR-REPLAY-110":
        raise RuntimeError("#775 reservoir child identity drifted")
    if reservoir_parent.RELATED_FRESH_HOLDOUT_PR != RELATED_FRESH_HOLDOUT_PR:
        raise RuntimeError("fresh holdout identity drifted")
    for time in (0.25, *PROBE_TIMES):
        if abs(toroidal.activation_g1(time) - reservoir_parent.screen.activation(time)) > 1.0e-15:
            raise RuntimeError("shared g1 time activation drifted")


def _probe_points(abs_z_values: tuple[float, ...] = PROBE_ABS_Z) -> np.ndarray:
    rows: list[list[float]] = []
    for radius in PROBE_RADII:
        for abs_z in abs_z_values:
            for sign in (-1.0, 1.0):
                for angle in PROBE_ANGLES:
                    rows.append(
                        [
                            float(radius * math.cos(angle)),
                            float(radius * math.sin(angle)),
                            float(sign * abs_z),
                        ]
                    )
    points = np.asarray(rows, dtype=float)
    expected = len(PROBE_RADII) * len(abs_z_values) * 2 * len(PROBE_ANGLES)
    if points.shape != (expected, 3):
        raise RuntimeError("frozen visible-tip probe cloud drifted")
    return points


def _build_reservoir_parent() -> tuple[VelocityFn, float, dict[str, Any]]:
    """Reconstruct exact #775 field and its frozen development-only alpha rule."""
    _assert_source_lock()
    _, _, source_652, _, _ = reservoir_parent.screen._build_fields()
    times, positions, metadata = reservoir_parent.screen._integrate(source_652)
    alpha, calibration = reservoir_parent.derive_reservoir_alpha(source_652, times, positions, metadata)

    def parent_fn(points: np.ndarray, time: float) -> np.ndarray:
        return np.asarray(
            reservoir_parent.reservoir_velocity(source_652, np.asarray(points, float), float(time), alpha),
            dtype=float,
        )

    return parent_fn, float(alpha), calibration


def toroidal_time_correction(points: np.ndarray, time: float) -> np.ndarray:
    return float(toroidal.activation_g1(float(time))) * toroidal.toroidal_unit_velocity(
        np.asarray(points, dtype=float)
    )


def derive_beta(parent_fn: VelocityFn) -> tuple[float, dict[str, Any]]:
    """Issue #845 projection-strengthening coefficient, with no scan/optimizer."""
    points = _probe_points()
    parent = toroidal.cylindrical_components(points, parent_fn(points, CALIBRATION_TIME))
    unit = toroidal.cylindrical_components(points, toroidal_time_correction(points, CALIBRATION_TIME))
    p = np.asarray(parent[:, 1], dtype=float)
    c = np.asarray(unit[:, 1], dtype=float)
    denom = float(np.dot(c, c))
    parent_norm = float(np.linalg.norm(p))
    unit_norm = float(np.linalg.norm(c))
    if not np.isfinite(denom) or denom <= np.finfo(float).tiny:
        raise RuntimeError("degenerate toroidal calibration response")
    if not np.isfinite(parent_norm) or parent_norm <= np.finfo(float).tiny:
        raise RuntimeError("degenerate parent swirl on frozen calibration cloud")
    beta_projection = float(np.dot(p, c) / denom)
    beta = float(PROJECTION_STRENGTH * beta_projection)
    correction_norm = float(abs(beta) * unit_norm)
    relative = float(correction_norm / parent_norm)
    if not np.isfinite(beta) or not np.isfinite(relative):
        raise RuntimeError("nonfinite toroidal coefficient calibration")
    nontrivial = bool(relative > NONTRIVIAL_RELATIVE_CORRECTION_MIN)
    cosine = float(np.dot(p, c) / (parent_norm * unit_norm)) if unit_norm > 0.0 else 0.0
    return beta, {
        "time": CALIBRATION_TIME,
        "probe_count": int(len(points)),
        "projection_strength": PROJECTION_STRENGTH,
        "parent_u_theta_l2": parent_norm,
        "unit_time_enveloped_u_theta_l2": unit_norm,
        "parent_unit_cosine": cosine,
        "beta_projection": beta_projection,
        "beta_theta": beta,
        "correction_to_parent_u_theta_l2_ratio": relative,
        "nontrivial_relative_correction_min": NONTRIVIAL_RELATIVE_CORRECTION_MIN,
        "nontrivial": nontrivial,
        "coefficient_scan_performed": False,
        "optimizer_used": False,
        "source_numeric_target_used": False,
        "fresh_714_data_used": False,
    }


def child_velocity(parent_fn: VelocityFn, points: np.ndarray, time: float, beta: float) -> np.ndarray:
    points = np.asarray(points, dtype=float)
    return np.asarray(parent_fn(points, float(time)), dtype=float) + float(beta) * toroidal_time_correction(
        points, float(time)
    )


def _ring_radius_span(points: np.ndarray, cylindrical: np.ndarray) -> float:
    radii = np.hypot(points[:, 0], points[:, 1])
    values: list[float] = []
    for radius in PROBE_RADII:
        mask = np.isclose(radii, radius, rtol=0.0, atol=2.0e-14)
        values.append(float(np.mean(np.abs(cylindrical[mask, 1]))))
    return float(np.ptp(np.asarray(values, dtype=float)))


def morphology_fingerprint(parent_fn: VelocityFn, beta: float) -> dict[str, Any]:
    points = _probe_points()
    radii = np.hypot(points[:, 0], points[:, 1])
    rows: dict[str, Any] = {}
    all_radial_axial_invariant = True
    finite = True
    for time in PROBE_TIMES:
        parent_cart = np.asarray(parent_fn(points, time), dtype=float)
        child_cart = np.asarray(child_velocity(parent_fn, points, time, beta), dtype=float)
        parent = toroidal.cylindrical_components(points, parent_cart)
        child = toroidal.cylindrical_components(points, child_cart)
        radial_change = float(np.max(np.abs(child[:, 0] - parent[:, 0])))
        axial_change = float(np.max(np.abs(child[:, 2] - parent[:, 2])))
        invariance = bool(radial_change <= ZERO_GATE and axial_change <= ZERO_GATE)
        all_radial_axial_invariant = all_radial_axial_invariant and invariance
        finite = finite and bool(np.all(np.isfinite(parent)) and np.all(np.isfinite(child)))

        ptheta = parent[:, 1]
        ctheta = child[:, 1]
        parent_rms = float(np.sqrt(np.mean(ptheta * ptheta)))
        child_rms = float(np.sqrt(np.mean(ctheta * ctheta)))
        parent_mean_abs = float(np.mean(np.abs(ptheta)))
        child_mean_abs = float(np.mean(np.abs(ctheta)))
        parent_angular_rms = float(np.sqrt(np.mean(np.square(ptheta / radii))))
        child_angular_rms = float(np.sqrt(np.mean(np.square(ctheta / radii))))
        parent_ratio = np.zeros_like(ptheta)
        child_ratio = np.zeros_like(ctheta)
        pvalid = np.abs(ptheta) > 1.0e-12
        cvalid = np.abs(ctheta) > 1.0e-12
        parent_ratio[pvalid] = -parent[pvalid, 0] / np.abs(ptheta[pvalid])
        child_ratio[cvalid] = -child[cvalid, 0] / np.abs(ctheta[cvalid])
        rows[f"{time:.3f}"] = {
            "parent_u_theta_rms": parent_rms,
            "child_u_theta_rms": child_rms,
            "u_theta_rms_delta": float(child_rms - parent_rms),
            "parent_mean_abs_u_theta": parent_mean_abs,
            "child_mean_abs_u_theta": child_mean_abs,
            "mean_abs_u_theta_delta": float(child_mean_abs - parent_mean_abs),
            "parent_angular_rate_rms": parent_angular_rms,
            "child_angular_rate_rms": child_angular_rms,
            "angular_rate_rms_delta": float(child_angular_rms - parent_angular_rms),
            "parent_circulation_speed_span_across_radii": _ring_radius_span(points, parent),
            "child_circulation_speed_span_across_radii": _ring_radius_span(points, child),
            "parent_mean_u_r": float(np.mean(parent[:, 0])),
            "child_mean_u_r": float(np.mean(child[:, 0])),
            "parent_inward_fraction": float(np.mean(parent[:, 0] < 0.0)),
            "child_inward_fraction": float(np.mean(child[:, 0] < 0.0)),
            "parent_mean_away_from_midplane_u_z": float(
                np.mean(np.sign(points[:, 2]) * parent[:, 2])
            ),
            "child_mean_away_from_midplane_u_z": float(
                np.mean(np.sign(points[:, 2]) * child[:, 2])
            ),
            "parent_mean_inward_to_circulation_ratio": (
                float(np.mean(parent_ratio[pvalid])) if np.any(pvalid) else 0.0
            ),
            "child_mean_inward_to_circulation_ratio": (
                float(np.mean(child_ratio[cvalid])) if np.any(cvalid) else 0.0
            ),
            "radial_component_change_max_abs": radial_change,
            "axial_component_change_max_abs": axial_change,
            "radial_axial_invariant": invariance,
        }

    midpoint = rows[f"{CALIBRATION_TIME:.3f}"]
    midpoint_swirl_increased = bool(
        midpoint["child_u_theta_rms"]
        > midpoint["parent_u_theta_rms"] + SWIRL_INCREASE_ABS_MIN
    )
    return {
        "probe_radii": list(PROBE_RADII),
        "probe_abs_z": list(PROBE_ABS_Z),
        "probe_azimuth_count": len(PROBE_ANGLES),
        "times": list(PROBE_TIMES),
        "by_time": rows,
        "midpoint_swirl_rms_strictly_increased": midpoint_swirl_increased,
        "midpoint_swirl_increase_abs_min": SWIRL_INCREASE_ABS_MIN,
        "radial_axial_invariant_all_times": bool(all_radial_axial_invariant),
        "finite": bool(finite),
        "source_numeric_target_applied": False,
        "scope": "autonomous velocity-space morphology fingerprint, not source correspondence",
    }


def structure_check(parent_fn: VelocityFn, beta: float) -> dict[str, Any]:
    rng = np.random.default_rng(9178451)
    start_points = rng.uniform(-1.8, 1.8, size=(96, 3))
    start_identity = float(
        np.max(
            np.abs(
                child_velocity(parent_fn, start_points, 0.25, beta)
                - np.asarray(parent_fn(start_points, 0.25), dtype=float)
            )
        )
    )

    outside = np.asarray(
        [
            [1.61, 0.0, 1.55],
            [-1.61, 0.0, -1.55],
            [0.8, 0.0, 0.99],
            [0.8, 0.0, -0.99],
            [0.8, 0.0, 2.0],
            [0.8, 0.0, -2.0],
            [0.8, 0.0, 2.05],
            [0.8, 0.0, -2.05],
        ],
        dtype=float,
    )
    outside_max = float(np.max(np.abs(toroidal.toroidal_unit_velocity(outside))))

    radius = rng.uniform(0.25, 1.45, size=64)
    angle = rng.uniform(0.0, 2.0 * np.pi, size=64)
    z_abs = rng.choice(np.asarray([1.15, 1.45, 1.65, 1.75, 1.88]), size=64)
    z = z_abs * rng.choice(np.asarray([-1.0, 1.0]), size=64)
    points = np.column_stack((radius * np.cos(angle), radius * np.sin(angle), z))
    divergence_max = 0.0
    for time in PROBE_TIMES:
        divergence = np.zeros(len(points), dtype=float)
        for axis in range(3):
            shift = np.zeros(3, dtype=float)
            shift[axis] = DIVERGENCE_STEP
            up = float(beta) * toroidal_time_correction(points + shift, time)
            down = float(beta) * toroidal_time_correction(points - shift, time)
            divergence += (up[:, axis] - down[:, axis]) / (2.0 * DIVERGENCE_STEP)
        divergence_max = max(divergence_max, float(np.max(np.abs(divergence))))

    capacity = toroidal.response_identifiability()
    symbolic = toroidal.symbolic_toroidal_divergence()
    passes = bool(
        start_identity <= ZERO_GATE
        and outside_max <= ZERO_GATE
        and divergence_max <= DIVERGENCE_MAX
        and symbolic["exact_identity"]
        and capacity["passes"]
    )
    return {
        "start_identity_max_abs": start_identity,
        "outside_unit_toroidal_max_abs": outside_max,
        "toroidal_correction_cartesian_fd_divergence_max": divergence_max,
        "divergence_gate": DIVERGENCE_MAX,
        "symbolic_toroidal_divergence_exact": bool(symbolic["exact_identity"]),
        "poloidal_toroidal_capacity": capacity,
        "passes": passes,
    }


def run(out: Path) -> dict[str, Any]:
    _assert_source_lock()
    parent_fn, alpha_reservoir, parent_calibration = _build_reservoir_parent()
    beta, beta_calibration = derive_beta(parent_fn)
    morphology = morphology_fingerprint(parent_fn, beta)
    structure = structure_check(parent_fn, beta)

    passed = bool(
        beta_calibration["nontrivial"]
        and morphology["finite"]
        and morphology["midpoint_swirl_rms_strictly_increased"]
        and morphology["radial_axial_invariant_all_times"]
        and structure["passes"]
    )

    report: dict[str, Any] = {
        "task_id": TASK_ID,
        "prereg_issue": PREREG_ISSUE,
        "source_parent": {"pr": SOURCE_PARENT_PR, "head": SOURCE_PARENT_HEAD},
        "source_reservoir_child": {
            "pr": SOURCE_RESERVOIR_CHILD_PR,
            "head": SOURCE_RESERVOIR_CHILD_HEAD,
            "alpha_reservoir": alpha_reservoir,
            "calibration": parent_calibration,
        },
        "frozen_toroidal_basis": {
            "cartesian": "(-y*H(r^2,z), x*H(r^2,z), 0)",
            "cylindrical": "u_r=0, u_theta=r*H, u_z=0",
            "radial_support": "r<1.6",
            "axial_support": "1<|z|<2",
            "time_activation": "g1(t)=2*(t-.25)",
            "parameter_count_increment": 1,
        },
        "coefficient_calibration": beta_calibration,
        "morphology_fingerprint": morphology,
        "structure_check": structure,
        "toroidal_child_morphology_screen_passed": passed,
        "held_out_pde_residual": {
            "evaluated": False,
            "reason": (
                "#775 has no newly matched pressure/forcing for this +1-parameter child; "
                "an unmatched residual would not be a same-protocol PDE-improvement result"
            ),
            "st006_comparison_performed": False,
        },
        "decision": {
            "classification": "toroidal_child_morphology_pass" if passed else "toroidal_child_morphology_fail",
            "actual_velocity_child_evaluated": True,
            "closer_visualization_delivery_supported": passed,
            "closer_visualization_delivery_scope": (
                "stronger midpoint swirl/circulation along the frozen toroidal direction while "
                "radial/axial probe components remain unchanged; not global/source correspondence"
            ),
            "second_poloidal_basis_justified_by_this_increment": False,
            "retune_beta_from_this_result_allowed": False,
            "pde_improvement_supported": False,
        },
        **TRUTH,
    }
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("outputs/agent7-st052m-toroidal-child-morphology-replay/report.json"),
    )
    args = parser.parse_args()
    print(json.dumps(run(args.out), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
