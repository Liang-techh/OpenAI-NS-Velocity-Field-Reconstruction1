"""Threshold-free morphology-convergence screen for the exact #559 child.

Preregistered in issue #575 before evaluation. This round changes no velocity
basis or coefficient. It freezes the exact compensated ST052-M child already
replayed by Agent 7 #568 and asks whether its morphology signal survives a
fixed, data-independent enstrophy-moment diagnostic on 25^3/33^3/41^3 grids.

This is expression-capacity / visualization-morphology evidence only. It is not
an OpenAI-image fit, PDE validation, production selection, or exact-field claim.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

import agent7_st052m_compensated_visual_fingerprint as prior

TASK_ID = "CR003-ST052M-THRESHOLD-FREE-MORPHOLOGY-084"
PREREG_ISSUE = 575
SOURCE_PRIOR_PR = 568
SOURCE_PRIOR_HEAD = "85f2f0940d80757313651d2ab50c6ef89dc22914"
SOURCE_COMPENSATED_PR = 559
SOURCE_COMPENSATED_HEAD = "39b106ad8cb8df2064cabead3a12682089575e74"
SOURCE_PATH_PR = 565
SOURCE_PATH_HEAD = "196bcc5445039a2b3533128daa5c20cbb151cd58"
EXPECTED_BETA = 0.08837490297155456
GRID_LEVELS = (25, 33, 41)
REFERENCE_TIMES = (0.25, 0.50, 0.75)
TIP_WEIGHT_WINDOW = (0.50, 0.80)
CENTRAL_BAND_MAX = 0.35

CRITERIA = {
    "tip_radial_change_max": -0.005,
    "full_aspect_change_min": 0.005,
    "central_radial_abs_change_max": 0.005,
    "tip_enstrophy_fraction_retention_min": 0.95,
    "fine_effect_drift_abs_max": 0.0075,
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


def tensor_trapezoid_weights(resolution: int) -> np.ndarray:
    """Return product trapezoid weights; the common h^3 factor is omitted."""
    if resolution < 3:
        raise ValueError("resolution must be at least 3")
    w = np.ones(int(resolution), dtype=float)
    w[[0, -1]] = 0.5
    return w[:, None, None] * w[None, :, None] * w[None, None, :]


def smooth_tip_weight(normalized_abs_z: np.ndarray) -> np.ndarray:
    """Fixed C-infinity-style compact bump on the preregistered physical tip band."""
    s = np.asarray(normalized_abs_z, dtype=float)
    a, b = TIP_WEIGHT_WINDOW
    out = np.zeros_like(s)
    mask = (s > a) & (s < b)
    xi = (s[mask] - a) / (b - a)
    out[mask] = np.exp(4.0 - 1.0 / (xi * (1.0 - xi)))
    return out


def _weighted_rms(values: np.ndarray, weight: np.ndarray) -> float:
    total = float(np.sum(weight))
    if total <= np.finfo(float).tiny:
        raise ValueError("empty weighted region")
    return float(np.sqrt(np.sum(weight * np.square(values)) / total))


def enstrophy_moment_metrics(magnitude: np.ndarray, axis: np.ndarray) -> dict:
    """Data-independent enstrophy moments on a fixed Cartesian grid."""
    magnitude = np.asarray(magnitude, dtype=float)
    axis = np.asarray(axis, dtype=float)
    n = axis.size
    if magnitude.shape != (n, n, n):
        raise ValueError("magnitude must have shape (n,n,n) matching axis")
    if n < 3 or not np.isfinite(magnitude).all():
        raise ValueError("invalid vorticity magnitude")

    xx, yy, zz = np.meshgrid(axis, axis, axis, indexing="ij")
    radius = np.hypot(xx, yy)
    abs_z = np.abs(zz)
    normalized_z = abs_z / 2.0
    trap = tensor_trapezoid_weights(n)
    e = np.square(magnitude) * trap
    total = float(np.sum(e))
    if total <= np.finfo(float).tiny:
        raise ValueError("zero full-grid enstrophy")

    axial_rms = _weighted_rms(abs_z, e)
    radial_rms = _weighted_rms(radius, e)
    if radial_rms <= np.finfo(float).tiny:
        raise ValueError("zero radial moment")

    tip_shape = smooth_tip_weight(normalized_z)
    tip_weight = e * tip_shape
    tip_total = float(np.sum(tip_weight))
    if tip_total <= np.finfo(float).tiny:
        raise ValueError("empty smooth tip weight")

    central_weight = e * (normalized_z <= CENTRAL_BAND_MAX)
    central_total = float(np.sum(central_weight))
    if central_total <= np.finfo(float).tiny:
        raise ValueError("empty central band")

    spacing = float(axis[1] - axis[0])
    return {
        "full_axial_rms": axial_rms,
        "full_radial_rms": radial_rms,
        "full_aspect_ratio": float(axial_rms / radial_rms),
        "smooth_tip_radial_rms": _weighted_rms(radius, tip_weight),
        "smooth_tip_axial_rms": _weighted_rms(abs_z, tip_weight),
        "smooth_tip_enstrophy_fraction": float(tip_total / total),
        "central_radial_rms": _weighted_rms(radius, central_weight),
        "full_enstrophy_trapezoid": float(total * spacing**3),
        "smooth_tip_enstrophy_trapezoid": float(tip_total * spacing**3),
    }


def _relative(child: float, control: float) -> float:
    if abs(control) <= np.finfo(float).tiny:
        raise ValueError("zero control in relative comparison")
    return float(child / control - 1.0)


def compare_metrics(control: dict, candidate: dict) -> dict:
    return {
        "full_axial_rms": _relative(candidate["full_axial_rms"], control["full_axial_rms"]),
        "full_radial_rms": _relative(candidate["full_radial_rms"], control["full_radial_rms"]),
        "full_aspect_ratio": _relative(candidate["full_aspect_ratio"], control["full_aspect_ratio"]),
        "smooth_tip_radial_rms": _relative(
            candidate["smooth_tip_radial_rms"], control["smooth_tip_radial_rms"]
        ),
        "smooth_tip_axial_rms": _relative(
            candidate["smooth_tip_axial_rms"], control["smooth_tip_axial_rms"]
        ),
        "central_radial_rms": _relative(candidate["central_radial_rms"], control["central_radial_rms"]),
        "full_enstrophy": _relative(
            candidate["full_enstrophy_trapezoid"], control["full_enstrophy_trapezoid"]
        ),
        "tip_enstrophy_fraction_retention": float(
            candidate["smooth_tip_enstrophy_fraction"] / control["smooth_tip_enstrophy_fraction"]
        ),
    }


def _point_rule(relative: dict) -> bool:
    c = CRITERIA
    return bool(
        relative["smooth_tip_radial_rms"] <= c["tip_radial_change_max"]
        and relative["full_aspect_ratio"] >= c["full_aspect_change_min"]
        and abs(relative["central_radial_rms"]) <= c["central_radial_abs_change_max"]
        and relative["tip_enstrophy_fraction_retention"] >= c["tip_enstrophy_fraction_retention_min"]
    )


def run(out: Path) -> dict:
    field, raw = prior.replay_st052.reconstruct()
    redistribution_scale, _, _ = prior.base.child_scale(field, raw)
    energy_solve = prior.comp.solve_energy_beta(field, raw, redistribution_scale)
    if not energy_solve["root_exists"]:
        raise RuntimeError("#559 frozen energy root disappeared")
    beta = float(energy_solve["beta"])
    if abs(beta - EXPECTED_BETA) > 2.0e-10:
        raise RuntimeError("#559 beta identity drift")

    control = lambda p, t: prior.taper.control_velocity(field, raw, p, t, redistribution_scale)
    compensated = lambda p, t: prior.comp.compensated_velocity(
        field,
        raw,
        p,
        t,
        redistribution_scale=redistribution_scale,
        beta=beta,
    )

    levels: dict[str, list[dict]] = {}
    point_rules: dict[str, list[bool]] = {}
    for resolution in GRID_LEVELS:
        rows = []
        rules = []
        for time in REFERENCE_TIMES:
            axis, u_control = prior.sample_velocity_grid(control, time, resolution)
            _, u_comp = prior.sample_velocity_grid(compensated, time, resolution)
            spacing = float(axis[1] - axis[0])
            metrics = {}
            for label, current in (("control", u_control), ("compensated", u_comp)):
                _, omega_mag = prior.vorticity(current, spacing)
                metrics[label] = enstrophy_moment_metrics(omega_mag, axis)
            relative = compare_metrics(metrics["control"], metrics["compensated"])
            rules.append(_point_rule(relative))
            rows.append({"time": float(time), "metrics": metrics, "relative_compensated_vs_control": relative})
        levels[str(resolution)] = rows
        point_rules[str(resolution)] = rules

    fine_drift_rules = []
    fine_drifts = []
    rows33 = levels["33"]
    rows41 = levels["41"]
    for row33, row41 in zip(rows33, rows41):
        if row33["time"] != row41["time"]:
            raise RuntimeError("fine-grid time mismatch")
        r33 = row33["relative_compensated_vs_control"]
        r41 = row41["relative_compensated_vs_control"]
        tip_drift = float(abs(r41["smooth_tip_radial_rms"] - r33["smooth_tip_radial_rms"]))
        aspect_drift = float(abs(r41["full_aspect_ratio"] - r33["full_aspect_ratio"]))
        passed = bool(
            tip_drift <= CRITERIA["fine_effect_drift_abs_max"]
            and aspect_drift <= CRITERIA["fine_effect_drift_abs_max"]
        )
        fine_drift_rules.append(passed)
        fine_drifts.append(
            {
                "time": row33["time"],
                "tip_radial_relative_effect_abs_drift": tip_drift,
                "full_aspect_relative_effect_abs_drift": aspect_drift,
                "passed": passed,
            }
        )

    clean = bool(all(all(v) for v in point_rules.values()) and all(fine_drift_rules))
    report = {
        "task_id": TASK_ID,
        "prereg_issue": PREREG_ISSUE,
        "source_prior_pr": SOURCE_PRIOR_PR,
        "source_prior_head": SOURCE_PRIOR_HEAD,
        "source_compensated_pr": SOURCE_COMPENSATED_PR,
        "source_compensated_head": SOURCE_COMPENSATED_HEAD,
        "source_path_pr": SOURCE_PATH_PR,
        "source_path_head": SOURCE_PATH_HEAD,
        "frozen_transform": {
            "redistribution_alpha": prior.base.SOURCE_ALPHA,
            "redistribution_gain": prior.base.GAIN,
            "redistribution_inner_window": list(prior.base.INNER_WINDOW),
            "redistribution_outer_window": list(prior.base.OUTER_WINDOW),
            "taper_tau": prior.taper.TAPER_TAU,
            "taper_window_abs_z_over_2": list(prior.taper.WINDOW),
            "shoulder_window_abs_z_over_2": list(prior.comp.SHOULDER_WINDOW),
            "beta": beta,
            "post_transform_common_scale": 1.0,
        },
        "protocol": {
            "times": list(REFERENCE_TIMES),
            "grid_levels": list(GRID_LEVELS),
            "box": [-2.0, 2.0],
            "cartesian_curl_order": 2,
            "integration_rule": "tensor_product_trapezoid_common_h3_cancels_in_ratios",
            "data_dependent_threshold": None,
            "tip_weight_window_abs_z_over_2": list(TIP_WEIGHT_WINDOW),
            "tip_weight_formula": "exp(4-1/(xi*(1-xi))) inside fixed window; zero outside",
            "central_band_abs_z_over_2_max": CENTRAL_BAND_MAX,
            "parameter_scan_performed": False,
            "basis_growth_performed": False,
            "image_or_openai_numeric_target_used": False,
        },
        "energy_replay": {
            "beta": beta,
            "relative_energy_error": energy_solve["relative_energy_error"],
            "compensated_common_scale": 1.0,
        },
        "criteria": CRITERIA,
        "levels": levels,
        "point_rules": point_rules,
        "fine_33_to_41_drifts": fine_drifts,
        "fine_drift_rules": fine_drift_rules,
        "clean_threshold_free_morphology_convergence": clean,
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
    print("beta=", report["frozen_transform"]["beta"])
    print("clean_threshold_free_morphology_convergence=", report["clean_threshold_free_morphology_convergence"])
    for resolution in GRID_LEVELS:
        for row in report["levels"][str(resolution)]:
            r = row["relative_compensated_vs_control"]
            print(
                "N=", resolution,
                "t=", row["time"],
                "axial=", r["full_axial_rms"],
                "radial=", r["full_radial_rms"],
                "aspect=", r["full_aspect_ratio"],
                "tip_radial=", r["smooth_tip_radial_rms"],
                "central=", r["central_radial_rms"],
                "tip_fraction_retention=", r["tip_enstrophy_fraction_retention"],
            )
    for row in report["fine_33_to_41_drifts"]:
        print("fine_drift=", row)


if __name__ == "__main__":
    main()
