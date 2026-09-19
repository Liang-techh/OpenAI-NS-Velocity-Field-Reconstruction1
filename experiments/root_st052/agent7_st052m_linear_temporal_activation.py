"""Screen one fixed linear temporal basis on the exact ST052-M local-swirl child.

Preregistered in issue #585 before evaluation.  The static endpoint is exactly
Agent-7 PR #559: ST052-M + frozen kappa=.05 redistribution + localized tip
Piola tau=.05 + shoulder swirl beta=0.08837490297155456, with no common scale.
This round adds exactly one fixed time profile g(t)=2*(t-.25): both tau and beta
are multiplied by g.  No coefficient/window/ramp scan is allowed.

The question is narrow: can a single time-basis direction preserve the exact
late-frame morphology while reducing the cumulative off-midplane material-path
tradeoff measured in Agent-9 PR #574?  This is target-free expression-capacity
evidence only, not PDE validation or OpenAI-field identification.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp

import agent7_st052m_threshold_free_morphology as morph

TASK_ID = "CR003-ST052M-LINEAR-TEMPORAL-ACTIVATION-085"
PREREG_ISSUE = 585
SOURCE_AGENT7_PR = 559
SOURCE_AGENT7_HEAD = "39b106ad8cb8df2064cabead3a12682089575e74"
SOURCE_MORPH_PR = 576
SOURCE_MORPH_HEAD = "0ddf5f6b321ba78e27aa3ca3e073ef593e55869a"
SOURCE_PATH_PR = 574
SOURCE_RENDER_PR = 583

TAPER_TAU = 0.05
EXPECTED_BETA = 0.08837490297155456
REFERENCE_TIMES = (0.25, 0.50, 0.75)
GRID_RESOLUTION = 41
ENERGY_ORDER = 64
DIVERGENCE_STEP = 1.0e-5

TIME_INTERVAL = (0.25, 0.75)
SEED_RADII = (0.6, 0.9, 1.2)
SEED_BANDS = (("shoulder", 0.85), ("tip", 1.15))
SEED_ANGLES = 4
OUTPUT_SAMPLES = 33
SOLVER_METHOD = "DOP853"
SOLVER_RTOL = 1.0e-9
SOLVER_ATOL = 1.0e-11
SOLVER_MAX_STEP = 0.01

CRITERIA = {
    "reference_energy_relative_error_max": 1.0e-10,
    "late_static_identity_abs_max": 1.0e-11,
    "aggregate_mean_turns_relative_min": -0.015,
    "aggregate_pair_axial_change_relative_min": -0.05,
    "tip_pair_axial_change_relative_min": -0.25,
    "late_effect_retention_min": 0.999,
    "mid_effect_retention_min": 0.30,
    "support_max_abs": 1.0e-12,
    "divergence_fd_max": 1.0e-5,
    "response_rank_min": 2,
    "response_condition_max": 8.0,
}

TRUTH = {
    "canonical_velocity_changed": False,
    "saved_velocity_changed": False,
    "production_candidate_selected": False,
    "pressure_or_force_changed": False,
    "held_out_pde_residual_evaluated": False,
    "parent_pde_receipt_transferred": False,
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


def activation(time: float) -> float:
    time = float(time)
    if not (TIME_INTERVAL[0] - 1.0e-12 <= time <= TIME_INTERVAL[1] + 1.0e-12):
        raise ValueError("time outside frozen interval")
    return float(2.0 * (time - 0.25))


def ramp_velocity(field, raw, points, time, redistribution_scale, beta):
    g = activation(time)
    u = morph.prior.taper.taper_velocity(
        field,
        raw,
        np.asarray(points, dtype=float),
        float(time),
        redistribution_scale=redistribution_scale,
        tau=TAPER_TAU * g,
        taper_scale=1.0,
    )
    return morph.prior.comp.apply_swirl_compensation(points, u, beta * g)


def static_velocity(field, raw, points, time, redistribution_scale, beta):
    return morph.prior.comp.compensated_velocity(
        field,
        raw,
        np.asarray(points, dtype=float),
        float(time),
        redistribution_scale=redistribution_scale,
        beta=beta,
    )


def control_velocity(field, raw, points, time, redistribution_scale):
    return morph.prior.taper.control_velocity(
        field, raw, np.asarray(points, dtype=float), float(time), redistribution_scale
    )


def axisymmetric_energy_at_time(velocity_fn, time: float, order: int = ENERGY_ORDER) -> float:
    xr, wr = leggauss(int(order))
    xz, wz = leggauss(int(order))
    radius = xr + 1.0
    z = 2.0 * xz
    rr, zz = np.meshgrid(radius, z, indexing="ij")
    pts = np.column_stack((rr.ravel(), np.zeros(rr.size), zz.ravel()))
    u = np.asarray(velocity_fn(pts, float(time)), float).reshape(rr.shape + (3,))
    weight = wr[:, None] * (2.0 * wz[None, :]) * np.pi * rr
    return float(np.sum(weight * np.sum(u * u, axis=-1)))


def _seed_table():
    rows = []
    metadata = []
    for radius in SEED_RADII:
        for angle_index in range(SEED_ANGLES):
            angle = 2.0 * np.pi * angle_index / SEED_ANGLES
            x = radius * np.cos(angle)
            y = radius * np.sin(angle)
            for band, z_abs in SEED_BANDS:
                for sign in (-1, 1):
                    rows.append([x, y, sign * z_abs])
                    metadata.append(
                        {
                            "radius": float(radius),
                            "angle_index": int(angle_index),
                            "band": band,
                            "sign": int(sign),
                            "z_abs": float(z_abs),
                        }
                    )
    seeds = np.asarray(rows, dtype=float)
    if seeds.shape != (48, 3):
        raise RuntimeError("frozen off-midplane seed construction drift")
    return seeds, metadata


def measure_paths(velocity_fn) -> dict:
    seeds, metadata = _seed_table()
    t0, t1 = TIME_INTERVAL
    times = np.linspace(t0, t1, OUTPUT_SAMPLES)
    n = len(seeds)

    def rhs(time, flat):
        pts = np.asarray(flat, dtype=float).reshape(n, 3)
        value = np.asarray(velocity_fn(pts, float(time)), dtype=float)
        if value.shape != pts.shape or not np.isfinite(value).all():
            raise ValueError("invalid velocity during path integration")
        return value.reshape(-1)

    sol = solve_ivp(
        rhs,
        (t0, t1),
        seeds.reshape(-1),
        method=SOLVER_METHOD,
        t_eval=times,
        rtol=SOLVER_RTOL,
        atol=SOLVER_ATOL,
        max_step=SOLVER_MAX_STEP,
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    pos = np.asarray(sol.y.T, dtype=float).reshape(len(times), n, 3)
    radii = np.hypot(pos[:, :, 0], pos[:, :, 1])
    angles = np.unwrap(np.arctan2(pos[:, :, 1], pos[:, :, 0]), axis=0)
    radius_change = radii[-1] - radii[0]
    turns = np.abs((angles[-1] - angles[0]) / (2.0 * np.pi))

    lookup = {}
    for i, meta in enumerate(metadata):
        key = (meta["radius"], meta["angle_index"], meta["band"])
        lookup.setdefault(key, {})[meta["sign"]] = i
    pair_rows = []
    for (radius, angle_index, band), signs in lookup.items():
        neg, posi = signs[-1], signs[1]
        sep = pos[:, posi, 2] - pos[:, neg, 2]
        pair_rows.append(
            {
                "radius": float(radius),
                "angle_index": int(angle_index),
                "band": band,
                "change": float(sep[-1] - sep[0]),
            }
        )

    def summarize(indices, pairs):
        changes = np.asarray([row["change"] for row in pairs], dtype=float)
        return {
            "path_count": len(indices),
            "pair_count": len(pairs),
            "mean_radius_change": float(np.mean(radius_change[indices])),
            "contraction_magnitude": float(abs(np.mean(radius_change[indices]))),
            "inward_path_count": int(np.sum(radius_change[indices] < 0.0)),
            "mean_absolute_turns": float(np.mean(turns[indices])),
            "maximum_absolute_turns": float(np.max(turns[indices])),
            "mean_pair_axial_separation_change": float(np.mean(changes)),
        }

    all_indices = list(range(n))
    result = summarize(all_indices, pair_rows)
    result["by_seed_band"] = {}
    for band, _ in SEED_BANDS:
        indices = [i for i, meta in enumerate(metadata) if meta["band"] == band]
        pairs = [row for row in pair_rows if row["band"] == band]
        result["by_seed_band"][band] = summarize(indices, pairs)
    return result


def _relative(new: float, old: float) -> float:
    if abs(old) <= np.finfo(float).tiny:
        raise ValueError("zero control metric")
    return float(new / old - 1.0)


def compare_paths(candidate: dict, control: dict) -> dict:
    return {
        "mean_absolute_turns_relative": _relative(
            candidate["mean_absolute_turns"], control["mean_absolute_turns"]
        ),
        "contraction_magnitude_relative": _relative(
            candidate["contraction_magnitude"], control["contraction_magnitude"]
        ),
        "mean_pair_axial_separation_change_relative": _relative(
            candidate["mean_pair_axial_separation_change"],
            control["mean_pair_axial_separation_change"],
        ),
        "inward_path_count_delta": int(
            candidate["inward_path_count"] - control["inward_path_count"]
        ),
    }


def morphology_rows(control_fn, static_fn, ramp_fn) -> list[dict]:
    rows = []
    for time in REFERENCE_TIMES:
        axis, uc = morph.prior.sample_velocity_grid(control_fn, time, GRID_RESOLUTION)
        _, us = morph.prior.sample_velocity_grid(static_fn, time, GRID_RESOLUTION)
        _, ur = morph.prior.sample_velocity_grid(ramp_fn, time, GRID_RESOLUTION)
        spacing = float(axis[1] - axis[0])
        metrics = {}
        for label, vel in (("control", uc), ("static", us), ("ramp", ur)):
            _, omag = morph.prior.vorticity(vel, spacing)
            metrics[label] = morph.enstrophy_moment_metrics(omag, axis)
        static_rel = morph.compare_metrics(metrics["control"], metrics["static"])
        ramp_rel = morph.compare_metrics(metrics["control"], metrics["ramp"])
        aspect_retention = (
            float(ramp_rel["full_aspect_ratio"] / static_rel["full_aspect_ratio"])
            if abs(static_rel["full_aspect_ratio"]) > 1.0e-14
            else None
        )
        tip_retention = (
            float(ramp_rel["smooth_tip_radial_rms"] / static_rel["smooth_tip_radial_rms"])
            if abs(static_rel["smooth_tip_radial_rms"]) > 1.0e-14
            else None
        )
        rows.append(
            {
                "time": float(time),
                "activation": activation(time),
                "static_relative": static_rel,
                "ramp_relative": ramp_rel,
                "aspect_effect_retention": aspect_retention,
                "tip_radial_effect_retention": tip_retention,
                "ramp_static_velocity_max_abs": float(np.max(np.abs(ur - us))),
            }
        )
    return rows


def response_diagnostic(control_fn, static_fn, ramp_fn, seed: int = 9175851) -> dict:
    rng = np.random.default_rng(seed)
    blocks_static = []
    blocks_ramp = []
    for time in REFERENCE_TIMES:
        radius = rng.uniform(0.2, 1.45, size=12)
        angle = rng.uniform(0.0, 2.0 * np.pi, size=12)
        z = rng.uniform(-1.45, 1.45, size=12)
        pts = np.column_stack((radius * np.cos(angle), radius * np.sin(angle), z))
        u0 = control_fn(pts, time)
        blocks_static.append((static_fn(pts, time) - u0).ravel())
        blocks_ramp.append((ramp_fn(pts, time) - u0).ravel())
    a = np.concatenate(blocks_static)
    b = np.concatenate(blocks_ramp)
    matrix = np.column_stack((a / np.linalg.norm(a), b / np.linalg.norm(b)))
    singular = np.linalg.svd(matrix, compute_uv=False)
    rank = int(np.linalg.matrix_rank(matrix, tol=1.0e-10))
    cosine = float(np.dot(matrix[:, 0], matrix[:, 1]))
    return {
        "rank": rank,
        "singular_values": singular.tolist(),
        "condition_number": float(singular[0] / singular[-1]),
        "column_cosine": cosine,
    }


def structure_preflight(ramp_fn, seed: int = 9175852) -> dict:
    outside = np.array(
        [[2.05, 0, 0], [-2.05, 0, 0], [0, 2.05, 0], [0, -2.05, 0],
         [0, 0, 2.05], [0, 0, -2.05], [1.8, 0, 2.05]],
        dtype=float,
    )
    support_max = max(float(np.max(np.abs(ramp_fn(outside, t)))) for t in REFERENCE_TIMES)
    rng = np.random.default_rng(seed)
    pts = rng.uniform(-1.35, 1.35, size=(36, 3))
    div_max = 0.0
    for time in (0.50, 0.75):
        div = np.zeros(len(pts))
        for axis in range(3):
            shift = np.zeros(3)
            shift[axis] = DIVERGENCE_STEP
            up = ramp_fn(pts + shift, time)
            um = ramp_fn(pts - shift, time)
            div += (up[:, axis] - um[:, axis]) / (2.0 * DIVERGENCE_STEP)
        div_max = max(div_max, float(np.max(np.abs(div))))
    return {"support_max_abs": support_max, "divergence_fd_max": div_max}


def run(out: Path) -> dict:
    field, raw = morph.prior.replay_st052.reconstruct()
    redistribution_scale, _, _ = morph.prior.base.child_scale(field, raw)
    energy_solve = morph.prior.comp.solve_energy_beta(field, raw, redistribution_scale)
    if not energy_solve["root_exists"]:
        raise RuntimeError("#559 frozen energy root disappeared")
    beta = float(energy_solve["beta"])
    if abs(beta - EXPECTED_BETA) > 2.0e-10:
        raise RuntimeError("#559 beta identity drift")

    control_fn = lambda p, t: control_velocity(field, raw, p, t, redistribution_scale)
    static_fn = lambda p, t: static_velocity(field, raw, p, t, redistribution_scale, beta)
    ramp_fn = lambda p, t: ramp_velocity(field, raw, p, t, redistribution_scale, beta)

    energies = {}
    for time in REFERENCE_TIMES:
        energies[str(time)] = {
            "control": axisymmetric_energy_at_time(control_fn, time),
            "static": axisymmetric_energy_at_time(static_fn, time),
            "ramp": axisymmetric_energy_at_time(ramp_fn, time),
        }
    reference_energy_relative_error = float(
        energies["0.25"]["ramp"] / energies["0.25"]["control"] - 1.0
    )

    paths_control = measure_paths(control_fn)
    paths_static = measure_paths(static_fn)
    paths_ramp = measure_paths(ramp_fn)
    ramp_vs_control = compare_paths(paths_ramp, paths_control)
    static_vs_control = compare_paths(paths_static, paths_control)
    ramp_band = {
        band: compare_paths(paths_ramp["by_seed_band"][band], paths_control["by_seed_band"][band])
        for band, _ in SEED_BANDS
    }
    static_band = {
        band: compare_paths(paths_static["by_seed_band"][band], paths_control["by_seed_band"][band])
        for band, _ in SEED_BANDS
    }

    mrows = morphology_rows(control_fn, static_fn, ramp_fn)
    row_mid = next(row for row in mrows if row["time"] == 0.50)
    row_late = next(row for row in mrows if row["time"] == 0.75)
    response = response_diagnostic(control_fn, static_fn, ramp_fn)
    structure = structure_preflight(ramp_fn)

    clean = bool(
        abs(reference_energy_relative_error) <= CRITERIA["reference_energy_relative_error_max"]
        and row_late["ramp_static_velocity_max_abs"] <= CRITERIA["late_static_identity_abs_max"]
        and ramp_vs_control["mean_absolute_turns_relative"] >= CRITERIA["aggregate_mean_turns_relative_min"]
        and ramp_vs_control["mean_pair_axial_separation_change_relative"]
        >= CRITERIA["aggregate_pair_axial_change_relative_min"]
        and ramp_band["tip"]["mean_pair_axial_separation_change_relative"]
        >= CRITERIA["tip_pair_axial_change_relative_min"]
        and row_late["aspect_effect_retention"] >= CRITERIA["late_effect_retention_min"]
        and row_late["tip_radial_effect_retention"] >= CRITERIA["late_effect_retention_min"]
        and row_mid["ramp_relative"]["full_aspect_ratio"] > 0.0
        and row_mid["ramp_relative"]["smooth_tip_radial_rms"] < 0.0
        and row_mid["aspect_effect_retention"] >= CRITERIA["mid_effect_retention_min"]
        and row_mid["tip_radial_effect_retention"] >= CRITERIA["mid_effect_retention_min"]
        and structure["support_max_abs"] <= CRITERIA["support_max_abs"]
        and structure["divergence_fd_max"] <= CRITERIA["divergence_fd_max"]
        and response["rank"] >= CRITERIA["response_rank_min"]
        and response["condition_number"] <= CRITERIA["response_condition_max"]
    )

    report = {
        "task_id": TASK_ID,
        "prereg_issue": PREREG_ISSUE,
        "source_agent7_pr": SOURCE_AGENT7_PR,
        "source_agent7_head": SOURCE_AGENT7_HEAD,
        "source_morph_pr": SOURCE_MORPH_PR,
        "source_morph_head": SOURCE_MORPH_HEAD,
        "source_path_pr": SOURCE_PATH_PR,
        "source_render_pr": SOURCE_RENDER_PR,
        "frozen_transform": {
            "redistribution_gain": morph.prior.base.GAIN,
            "taper_tau_static_endpoint": TAPER_TAU,
            "shoulder_beta_static_endpoint": beta,
            "post_transform_common_scale": 1.0,
            "activation": "g(t)=2*(t-0.25)",
            "tau_of_t": "0.05*g(t)",
            "beta_of_t": "0.08837490297155456*g(t)",
            "parameter_scan_performed": False,
        },
        "criteria": CRITERIA,
        "energies": energies,
        "reference_energy_relative_error": reference_energy_relative_error,
        "offmidplane_path_protocol": {
            "time_interval": list(TIME_INTERVAL),
            "seed_radii": list(SEED_RADII),
            "seed_bands_abs_z": {band: z for band, z in SEED_BANDS},
            "seed_angles": SEED_ANGLES,
            "path_count": 48,
            "output_samples": OUTPUT_SAMPLES,
            "solver": SOLVER_METHOD,
            "rtol": SOLVER_RTOL,
            "atol": SOLVER_ATOL,
            "max_step": SOLVER_MAX_STEP,
        },
        "paths": {
            "control": paths_control,
            "static": paths_static,
            "ramp": paths_ramp,
            "static_vs_control": static_vs_control,
            "ramp_vs_control": ramp_vs_control,
            "static_vs_control_by_seed_band": static_band,
            "ramp_vs_control_by_seed_band": ramp_band,
        },
        "morphology_41": mrows,
        "response_diagnostic": response,
        "structure_preflight": structure,
        "clean_linear_temporal_activation_capacity": clean,
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
    print("clean_linear_temporal_activation_capacity=", report["clean_linear_temporal_activation_capacity"])
    print("reference_energy_relative_error=", report["reference_energy_relative_error"])
    print("static_vs_control=", report["paths"]["static_vs_control"])
    print("ramp_vs_control=", report["paths"]["ramp_vs_control"])
    print("ramp_tip=", report["paths"]["ramp_vs_control_by_seed_band"]["tip"])
    for row in report["morphology_41"]:
        print("morphology=", row)
    print("response=", report["response_diagnostic"])
    print("structure=", report["structure_preflight"])


if __name__ == "__main__":
    main()
