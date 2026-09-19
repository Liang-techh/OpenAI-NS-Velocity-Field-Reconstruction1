"""Resolution/time-evolution audit for the exact Agent-7 #587 linear ramp.

Preregistered in issue #593 before evaluation. This round changes no velocity
basis or coefficient. It freezes the exact #587 child and asks whether its
threshold-free vorticity morphology evolves with the fixed linear activation in
a directionally monotone, resolution-stable way on 25^3/33^3/41^3 grids.

This is target-free expression-capacity evidence only, not PDE validation,
production selection, or OpenAI-field identification.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

import agent7_st052m_linear_temporal_activation as prior

TASK_ID = "CR003-ST052M-LINEAR-RAMP-CONVERGENCE-086"
PREREG_ISSUE = 593
SOURCE_TEMPORAL_PR = 587
SOURCE_TEMPORAL_HEAD = "0b93819095f6c8576a7571bdc2d2fbef4154944d"
SOURCE_STATIC_PR = 559
SOURCE_STATIC_HEAD = "39b106ad8cb8df2064cabead3a12682089575e74"
EXPECTED_BETA = 0.08837490297155456
GRID_LEVELS = (25, 33, 41)
REFERENCE_TIMES = (0.25, 0.375, 0.50, 0.625, 0.75)

CRITERIA = {
    "identity_relative_abs_max": 1.0e-12,
    "active_aspect_change_min": 0.0,
    "active_tip_radial_change_max": 0.0,
    "central_radial_abs_change_max": 0.005,
    "tip_enstrophy_fraction_retention_min": 0.95,
    "fine_effect_drift_abs_max": 0.0075,
    "finest_adjacent_effect_backslide_max": 0.001,
    "late_static_identity_abs_max": 1.0e-11,
}

TRUTH = {
    "basis_changed": False,
    "velocity_coefficient_changed": False,
    "canonical_velocity_changed": False,
    "saved_velocity_changed": False,
    "production_candidate_selected": False,
    "pressure_or_force_changed": False,
    "held_out_pde_residual_evaluated": False,
    "parent_pde_receipt_transferred": False,
    "material_paths_integrated_this_round": False,
    "public_image_numeric_target_used": False,
    "data_dependent_morphology_threshold_used": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "source_correspondence_verified": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _active_point_pass(relative: dict) -> bool:
    c = CRITERIA
    return bool(
        relative["full_aspect_ratio"] > c["active_aspect_change_min"]
        and relative["smooth_tip_radial_rms"] < c["active_tip_radial_change_max"]
        and abs(relative["central_radial_rms"]) <= c["central_radial_abs_change_max"]
        and relative["tip_enstrophy_fraction_retention"] >= c["tip_enstrophy_fraction_retention_min"]
    )


def _metric_row(control_velocity: np.ndarray, ramp_velocity: np.ndarray, axis: np.ndarray) -> dict:
    spacing = float(axis[1] - axis[0])
    _, omega_control = prior.morph.prior.vorticity(control_velocity, spacing)
    _, omega_ramp = prior.morph.prior.vorticity(ramp_velocity, spacing)
    control_metrics = prior.morph.enstrophy_moment_metrics(omega_control, axis)
    ramp_metrics = prior.morph.enstrophy_moment_metrics(omega_ramp, axis)
    relative = prior.morph.compare_metrics(control_metrics, ramp_metrics)
    return {"control": control_metrics, "ramp": ramp_metrics, "relative": relative}


def run(out: Path) -> dict:
    field, raw = prior.morph.prior.replay_st052.reconstruct()
    redistribution_scale, _, _ = prior.morph.prior.base.child_scale(field, raw)
    energy_solve = prior.morph.prior.comp.solve_energy_beta(field, raw, redistribution_scale)
    if not energy_solve["root_exists"]:
        raise RuntimeError("#559 frozen energy root disappeared")
    beta = float(energy_solve["beta"])
    if abs(beta - EXPECTED_BETA) > 2.0e-10:
        raise RuntimeError("#559 beta identity drift")

    control_fn = lambda p, t: prior.control_velocity(field, raw, p, t, redistribution_scale)
    ramp_fn = lambda p, t: prior.ramp_velocity(field, raw, p, t, redistribution_scale, beta)
    static_fn = lambda p, t: prior.static_velocity(field, raw, p, t, redistribution_scale, beta)

    levels: dict[str, list[dict]] = {}
    point_rules: dict[str, list[bool]] = {}
    identity_rules: dict[str, bool] = {}
    for resolution in GRID_LEVELS:
        rows = []
        rules = []
        for time in REFERENCE_TIMES:
            axis, u_control = prior.morph.prior.sample_velocity_grid(control_fn, time, resolution)
            _, u_ramp = prior.morph.prior.sample_velocity_grid(ramp_fn, time, resolution)
            metrics = _metric_row(u_control, u_ramp, axis)
            relative = metrics["relative"]
            rows.append(
                {
                    "time": float(time),
                    "activation": prior.activation(time),
                    "metrics": metrics,
                    "active_point_pass": True if time == REFERENCE_TIMES[0] else _active_point_pass(relative),
                }
            )
            if time != REFERENCE_TIMES[0]:
                rules.append(_active_point_pass(relative))
        first = rows[0]["metrics"]["relative"]
        identity_rules[str(resolution)] = bool(
            abs(first["full_aspect_ratio"]) <= CRITERIA["identity_relative_abs_max"]
            and abs(first["smooth_tip_radial_rms"]) <= CRITERIA["identity_relative_abs_max"]
        )
        levels[str(resolution)] = rows
        point_rules[str(resolution)] = rules

    fine_drifts = []
    fine_rules = []
    for row33, row41 in zip(levels["33"][1:], levels["41"][1:]):
        if row33["time"] != row41["time"]:
            raise RuntimeError("fine-grid time mismatch")
        r33 = row33["metrics"]["relative"]
        r41 = row41["metrics"]["relative"]
        aspect = float(abs(r41["full_aspect_ratio"] - r33["full_aspect_ratio"]))
        tip = float(abs(r41["smooth_tip_radial_rms"] - r33["smooth_tip_radial_rms"]))
        passed = bool(
            aspect <= CRITERIA["fine_effect_drift_abs_max"]
            and tip <= CRITERIA["fine_effect_drift_abs_max"]
        )
        fine_rules.append(passed)
        fine_drifts.append(
            {
                "time": row41["time"],
                "aspect_effect_abs_drift": aspect,
                "tip_radial_effect_abs_drift": tip,
                "passed": passed,
            }
        )

    finest = levels["41"]
    finest_steps = []
    monotone_rules = []
    for prev, cur in zip(finest, finest[1:]):
        rp = prev["metrics"]["relative"]
        rc = cur["metrics"]["relative"]
        aspect_step = float(rc["full_aspect_ratio"] - rp["full_aspect_ratio"])
        tip_thinning_prev = float(-rp["smooth_tip_radial_rms"])
        tip_thinning_cur = float(-rc["smooth_tip_radial_rms"])
        tip_step = float(tip_thinning_cur - tip_thinning_prev)
        passed = bool(
            aspect_step >= -CRITERIA["finest_adjacent_effect_backslide_max"]
            and tip_step >= -CRITERIA["finest_adjacent_effect_backslide_max"]
        )
        monotone_rules.append(passed)
        finest_steps.append(
            {
                "from_time": prev["time"],
                "to_time": cur["time"],
                "aspect_gain_step": aspect_step,
                "tip_thinning_magnitude_step": tip_step,
                "passed": passed,
            }
        )

    axis41, u_ramp_late = prior.morph.prior.sample_velocity_grid(ramp_fn, 0.75, 41)
    _, u_static_late = prior.morph.prior.sample_velocity_grid(static_fn, 0.75, 41)
    late_static_identity = float(np.max(np.abs(u_ramp_late - u_static_late)))

    clean = bool(
        all(identity_rules.values())
        and all(all(v) for v in point_rules.values())
        and all(fine_rules)
        and all(monotone_rules)
        and late_static_identity <= CRITERIA["late_static_identity_abs_max"]
    )

    report = {
        "task_id": TASK_ID,
        "prereg_issue": PREREG_ISSUE,
        "source_temporal_pr": SOURCE_TEMPORAL_PR,
        "source_temporal_head": SOURCE_TEMPORAL_HEAD,
        "source_static_pr": SOURCE_STATIC_PR,
        "source_static_head": SOURCE_STATIC_HEAD,
        "frozen_transform": {
            "redistribution_gain": prior.morph.prior.base.GAIN,
            "taper_tau_static_endpoint": prior.TAPER_TAU,
            "shoulder_beta_static_endpoint": beta,
            "activation": "g(t)=2*(t-0.25)",
            "post_transform_common_scale": 1.0,
            "parameter_scan_performed": False,
        },
        "grid_levels": list(GRID_LEVELS),
        "reference_times": list(REFERENCE_TIMES),
        "criteria": CRITERIA,
        "levels": levels,
        "identity_rules": identity_rules,
        "point_rules": point_rules,
        "fine_drifts": fine_drifts,
        "finest_adjacent_effect_steps": finest_steps,
        "late_static_identity_abs_max": late_static_identity,
        "clean_linear_ramp_morphology_convergence": clean,
        **TRUTH,
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    report = run(args.out)
    print("clean_linear_ramp_morphology_convergence=", report["clean_linear_ramp_morphology_convergence"])
    print("late_static_identity_abs_max=", report["late_static_identity_abs_max"])
    for level in GRID_LEVELS:
        for row in report["levels"][str(level)]:
            print("morphology=", level, row["time"], row["metrics"]["relative"])
    print("fine_drifts=", report["fine_drifts"])
    print("finest_steps=", report["finest_adjacent_effect_steps"])


if __name__ == "__main__":
    main()
