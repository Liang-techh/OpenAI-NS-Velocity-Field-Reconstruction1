"""Fixed-grid vorticity-morphology preservation screen for the #559 child.

Preregistered in issue #567 before evaluation.  This round freezes the exact
ST052-M + kappa=.05 redistribution + tau=.05 localized tip Piola transform +
energy-only shoulder swirl compensation from Agent-7 PR #559.  It asks only
whether that compensated child preserves the target-free fixed-grid vorticity
morphology previously seen for the localized taper.  No coefficient/window
scan, image-derived numeric target, pressure/forcing refit, or PDE promotion is
allowed.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

import agent7_st052m_frozen_redistribution as base
import agent7_st052m_local_swirl_energy as comp
import agent7_st052m_localized_taper as taper
import replay_st052

TASK_ID = "CR003-ST052M-COMPENSATED-VISUAL-FINGERPRINT-083"
PREREG_ISSUE = 567
SOURCE_COMPENSATED_PR = 559
SOURCE_COMPENSATED_HEAD = "39b106ad8cb8df2064cabead3a12682089575e74"
SOURCE_RENDER_PR = 556
SOURCE_RENDER_HEAD = "619e3d2ff8c0231c9aced4ce20216c03d78cdef8"
SOURCE_PATH_PR = 565
SOURCE_PATH_HEAD = "196bcc5445039a2b3533128daa5c20cbb151cd58"
EXPECTED_BETA = 0.08837490297155456
GRID_LEVELS = (25, 33)
REFERENCE_TIMES = (0.25, 0.50, 0.75)
VORTICITY_QUANTILE = 0.985
TIP_BAND = (0.50, 0.80)
CENTRAL_BAND_MAX = 0.35

CRITERIA = {
    "fine_top_axial_gain_min": 0.04,
    "fine_top_radial_change_max": -0.015,
    "fine_tip_radial_change_max": -0.007,
    "fine_central_radial_abs_change_max": 0.005,
    "fine_axial_gain_retention_min": 0.70,
    "cross_resolution_sign_agreement_required": True,
}

TRUTH = {
    "canonical_velocity_changed": False,
    "saved_velocity_changed": False,
    "production_candidate_selected": False,
    "production_taper_selected": False,
    "production_compensation_selected": False,
    "pressure_or_force_changed": False,
    "held_out_pde_residual_evaluated": False,
    "parent_pde_receipt_transferred": False,
    "material_paths_integrated_this_round": False,
    "public_image_numeric_target_used": False,
    "visual_acceptance_threshold_defined": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "source_correspondence_verified": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def sample_velocity_grid(velocity, time: float, resolution: int):
    if resolution < 9 or resolution % 2 == 0:
        raise ValueError("resolution must be odd and at least 9")
    axis = np.linspace(-2.0, 2.0, int(resolution), dtype=float)
    xx, yy, zz = np.meshgrid(axis, axis, axis, indexing="ij")
    pts = np.column_stack((xx.ravel(), yy.ravel(), zz.ravel()))
    values = np.asarray(velocity(pts, float(time)), dtype=float)
    if values.shape != pts.shape or not np.isfinite(values).all():
        raise ValueError("velocity must return finite shape-(n,3) values")
    return axis, values.reshape(resolution, resolution, resolution, 3)


def vorticity(field: np.ndarray, spacing: float):
    field = np.asarray(field, dtype=float)
    if field.ndim != 4 or field.shape[-1] != 3:
        raise ValueError("field must have shape (n,n,n,3)")
    u, v, w = (field[..., i] for i in range(3))
    _du_dx, du_dy, du_dz = np.gradient(u, spacing, spacing, spacing, edge_order=2)
    dv_dx, _dv_dy, dv_dz = np.gradient(v, spacing, spacing, spacing, edge_order=2)
    dw_dx, dw_dy, _dw_dz = np.gradient(w, spacing, spacing, spacing, edge_order=2)
    omega = np.stack((dw_dy - dv_dz, du_dz - dw_dx, dv_dx - du_dy), axis=-1)
    return omega, np.linalg.norm(omega, axis=-1)


def _weighted_radial_rms(magnitude, radius, mask):
    weight = np.square(magnitude) * mask
    total = float(np.sum(weight))
    if total <= np.finfo(float).tiny:
        raise ValueError("empty vorticity band")
    return float(np.sqrt(np.sum(weight * np.square(radius)) / total))


def vorticity_metrics(magnitude: np.ndarray, axis: np.ndarray):
    magnitude = np.asarray(magnitude, dtype=float)
    threshold = float(np.quantile(magnitude, VORTICITY_QUANTILE))
    selected = magnitude >= threshold
    xx, yy, zz = np.meshgrid(axis, axis, axis, indexing="ij")
    radius = np.hypot(xx, yy)
    abs_z = np.abs(zz)
    normalized_z = abs_z / 2.0
    tip = (normalized_z >= TIP_BAND[0]) & (normalized_z <= TIP_BAND[1])
    central = normalized_z <= CENTRAL_BAND_MAX
    return {
        "threshold": threshold,
        "selected_points": int(np.count_nonzero(selected)),
        "selected_axial_rms": float(np.sqrt(np.mean(np.square(abs_z[selected])))),
        "selected_radial_rms": float(np.sqrt(np.mean(np.square(radius[selected])))),
        "selected_max_abs_z": float(np.max(abs_z[selected])),
        "tip_enstrophy_weighted_radial_rms": _weighted_radial_rms(magnitude, radius, tip),
        "central_enstrophy_weighted_radial_rms": _weighted_radial_rms(magnitude, radius, central),
        "full_vorticity_rms": float(np.sqrt(np.mean(np.square(magnitude)))),
    }


def _relative(child: float, control: float) -> float:
    if abs(control) <= np.finfo(float).tiny:
        raise ValueError("zero control in relative comparison")
    return float(child / control - 1.0)


def compare_metrics(control: dict, candidate: dict):
    return {
        "top_axial_rms": _relative(candidate["selected_axial_rms"], control["selected_axial_rms"]),
        "top_radial_rms": _relative(candidate["selected_radial_rms"], control["selected_radial_rms"]),
        "tip_radial_rms": _relative(
            candidate["tip_enstrophy_weighted_radial_rms"],
            control["tip_enstrophy_weighted_radial_rms"],
        ),
        "central_radial_rms": _relative(
            candidate["central_enstrophy_weighted_radial_rms"],
            control["central_enstrophy_weighted_radial_rms"],
        ),
        "full_vorticity_rms": _relative(candidate["full_vorticity_rms"], control["full_vorticity_rms"]),
        "top_max_abs_z_delta": float(candidate["selected_max_abs_z"] - control["selected_max_abs_z"]),
    }


def _same_required_sign(coarse: dict, fine: dict) -> bool:
    return bool(
        coarse["top_axial_rms"] > 0.0
        and fine["top_axial_rms"] > 0.0
        and coarse["top_radial_rms"] < 0.0
        and fine["top_radial_rms"] < 0.0
        and coarse["tip_radial_rms"] < 0.0
        and fine["tip_radial_rms"] < 0.0
    )


def run(out: Path):
    field, raw = replay_st052.reconstruct()
    redistribution_scale, _, _ = base.child_scale(field, raw)
    energy_solve = comp.solve_energy_beta(field, raw, redistribution_scale)
    if not energy_solve["root_exists"]:
        raise RuntimeError("#559 frozen energy root disappeared")
    beta = float(energy_solve["beta"])
    if abs(beta - EXPECTED_BETA) > 2.0e-10:
        raise RuntimeError("#559 beta identity drift")
    taper_scale, control_energy, raw_taper_energy = taper.taper_energy_scale(
        field, raw, redistribution_scale
    )

    control = lambda p, t: taper.control_velocity(field, raw, p, t, redistribution_scale)
    normalized_taper = lambda p, t: taper.taper_velocity(
        field,
        raw,
        p,
        t,
        redistribution_scale=redistribution_scale,
        tau=taper.TAPER_TAU,
        taper_scale=taper_scale,
    )
    compensated = lambda p, t: comp.compensated_velocity(
        field,
        raw,
        p,
        t,
        redistribution_scale=redistribution_scale,
        beta=beta,
    )

    levels = {}
    for resolution in GRID_LEVELS:
        rows = []
        for time in REFERENCE_TIMES:
            axis, u_control = sample_velocity_grid(control, time, resolution)
            _, u_taper = sample_velocity_grid(normalized_taper, time, resolution)
            _, u_comp = sample_velocity_grid(compensated, time, resolution)
            spacing = float(axis[1] - axis[0])
            metrics = {}
            for label, current in (
                ("control", u_control),
                ("normalized_taper", u_taper),
                ("compensated", u_comp),
            ):
                _, omega_mag = vorticity(current, spacing)
                metrics[label] = vorticity_metrics(omega_mag, axis)
            rel_taper = compare_metrics(metrics["control"], metrics["normalized_taper"])
            rel_comp = compare_metrics(metrics["control"], metrics["compensated"])
            taper_gain = rel_taper["top_axial_rms"]
            retention = None if taper_gain <= 0.0 else float(rel_comp["top_axial_rms"] / taper_gain)
            rows.append(
                {
                    "time": float(time),
                    "metrics": metrics,
                    "relative_normalized_taper_vs_control": rel_taper,
                    "relative_compensated_vs_control": rel_comp,
                    "compensated_axial_gain_retention_of_taper": retention,
                }
            )
        levels[str(resolution)] = rows

    fine = levels["33"]
    coarse = levels["25"]
    c = CRITERIA
    fine_rules = []
    sign_rules = []
    for coarse_row, fine_row in zip(coarse, fine):
        if coarse_row["time"] != fine_row["time"]:
            raise RuntimeError("grid-level time mismatch")
        r = fine_row["relative_compensated_vs_control"]
        retention = fine_row["compensated_axial_gain_retention_of_taper"]
        fine_rules.append(
            bool(
                r["top_axial_rms"] >= c["fine_top_axial_gain_min"]
                and r["top_radial_rms"] <= c["fine_top_radial_change_max"]
                and r["tip_radial_rms"] <= c["fine_tip_radial_change_max"]
                and abs(r["central_radial_rms"]) <= c["fine_central_radial_abs_change_max"]
                and retention is not None
                and retention >= c["fine_axial_gain_retention_min"]
            )
        )
        sign_rules.append(
            _same_required_sign(
                coarse_row["relative_compensated_vs_control"],
                fine_row["relative_compensated_vs_control"],
            )
        )

    clean = bool(all(fine_rules) and all(sign_rules))
    report = {
        "task_id": TASK_ID,
        "prereg_issue": PREREG_ISSUE,
        "source_compensated_pr": SOURCE_COMPENSATED_PR,
        "source_compensated_head": SOURCE_COMPENSATED_HEAD,
        "source_render_pr": SOURCE_RENDER_PR,
        "source_render_head": SOURCE_RENDER_HEAD,
        "source_path_pr": SOURCE_PATH_PR,
        "source_path_head": SOURCE_PATH_HEAD,
        "frozen_transform": {
            "redistribution_alpha": base.SOURCE_ALPHA,
            "redistribution_gain": base.GAIN,
            "redistribution_inner_window": list(base.INNER_WINDOW),
            "redistribution_outer_window": list(base.OUTER_WINDOW),
            "taper_tau": taper.TAPER_TAU,
            "taper_window_abs_z_over_2": list(taper.WINDOW),
            "shoulder_window_abs_z_over_2": list(comp.SHOULDER_WINDOW),
            "beta": beta,
            "post_transform_common_scale": 1.0,
        },
        "protocol": {
            "times": list(REFERENCE_TIMES),
            "grid_levels": list(GRID_LEVELS),
            "box": [-2.0, 2.0],
            "cartesian_curl_order": 2,
            "vorticity_quantile": VORTICITY_QUANTILE,
            "tip_band_abs_z_over_2": list(TIP_BAND),
            "central_band_abs_z_over_2_max": CENTRAL_BAND_MAX,
            "parameter_scan_performed": False,
            "image_or_openai_numeric_target_used": False,
            "visual_acceptance_threshold_used": False,
        },
        "energy_replay": {
            "beta": beta,
            "relative_energy_error": energy_solve["relative_energy_error"],
            "control_energy": control_energy,
            "raw_taper_energy": raw_taper_energy,
            "normalized_taper_scale": taper_scale,
            "compensated_common_scale": 1.0,
        },
        "criteria": CRITERIA,
        "levels": levels,
        "fine_rules": fine_rules,
        "cross_resolution_sign_rules": sign_rules,
        "clean_compensated_visual_fingerprint_preservation": clean,
        **TRUTH,
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    report = run(args.out)
    print("beta=", report["frozen_transform"]["beta"])
    print("clean_compensated_visual_fingerprint_preservation=", report["clean_compensated_visual_fingerprint_preservation"])
    for resolution in GRID_LEVELS:
        for row in report["levels"][str(resolution)]:
            rel = row["relative_compensated_vs_control"]
            print(
                "N=", resolution,
                "t=", row["time"],
                "axial=", rel["top_axial_rms"],
                "radial=", rel["top_radial_rms"],
                "tip=", rel["tip_radial_rms"],
                "central=", rel["central_radial_rms"],
                "retention=", row["compensated_axial_gain_retention_of_taper"],
            )


if __name__ == "__main__":
    main()
