"""Target-free fixed-seed 3-D render for the ST052-M localized-taper tradeoff.

This Agent-9 diagnostic compares the exact ST052-M + frozen kappa=.05 swirl
redistribution control with Agent-7's exact localized radial-Piola taper child.
It is visualization evidence only: no OpenAI image is fitted, no visual pass
threshold is defined, and no PDE/pressure/forcing receipt is transferred.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib
import json
from pathlib import Path
import sys
from typing import Callable

import numpy as np
from scipy.interpolate import RegularGridInterpolator

TASK_ID = "CR-A9-061"
BASE_MAIN_SHA = "f0193d66c9d92948b4820ebcb70263673995b324"
ST052_PR = 508
ST052_HEAD = "b3b8bfdbe1077f9ec967d158602951997d81e17d"
REDISTRIBUTION_PR = 528
REDISTRIBUTION_HEAD = "779ffca71066e2864496d37de55a7aafc45d6f57"
PATH_CONTROL_PR = 533
TAPER_PR = 544
TAPER_HEAD = "093c7171cd61c6bd439afa30b2da69598a02d182"
TAPER_PATH_PR = 546
ENERGY_NEUTRAL_TAPER_PR = 549
ENERGY_NEUTRAL_TAPER_HEAD = "62e8c170427d5d830d7f897ba31768e0fc4ce56a"
EXPECTED_REDISTRIBUTION_SCALE = 1.0032534663681094
EXPECTED_TAPER_SCALE = 0.9929528816318394
REFERENCE_TIMES = (0.25, 0.50, 0.75)
SEED_RADII = (0.6, 0.9, 1.2)
SEED_Z = (-0.3, 0.3)
SEED_AZIMUTHS = 8
GRID_RESOLUTION = 33
STREAM_ARCLENGTH = 2.25
STREAM_OUTPUTS_PER_SIDE = 91
SPEED_FLOOR = 1e-8
VORTICITY_QUANTILE = 0.985
TIP_BAND = (0.50, 0.80)  # |z|/2, inherited observational band from PR #544
CENTRAL_BAND_MAX = 0.35
SCIPY_SOURCE_COMMIT = "b12c772edbc1fe0d3db9481cdfcc2e311569cb25"
MATPLOTLIB_RELEASE = "3.10.6"
MATPLOTLIB_SOURCE_COMMIT = "5cd38c3edcdf0792d0e6aded280a9b7a7de6146f"

Velocity = Callable[[np.ndarray, float], np.ndarray]

TRUTH = {
    "canonical_velocity_changed": False,
    "saved_velocity_changed": False,
    "production_candidate_selected": False,
    "production_taper_selected": False,
    "pressure_or_force_changed": False,
    "held_out_pde_residual_evaluated": False,
    "parent_pde_receipt_transferred": False,
    "public_image_numeric_target_used": False,
    "visual_pass_threshold_defined": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "source_correspondence_verified": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def seed_points() -> np.ndarray:
    points = []
    for radius in SEED_RADII:
        for z in SEED_Z:
            for index in range(SEED_AZIMUTHS):
                theta = 2.0 * np.pi * index / SEED_AZIMUTHS
                points.append((radius * np.cos(theta), radius * np.sin(theta), z))
    return np.asarray(points, dtype=np.float64)


def sample_velocity_grid(
    velocity: Velocity,
    *,
    times: tuple[float, ...] = REFERENCE_TIMES,
    resolution: int = GRID_RESOLUTION,
) -> tuple[np.ndarray, np.ndarray]:
    if resolution < 9 or resolution % 2 == 0:
        raise ValueError("resolution must be odd and at least 9")
    axis = np.linspace(-2.0, 2.0, int(resolution), dtype=np.float64)
    xx, yy, zz = np.meshgrid(axis, axis, axis, indexing="ij")
    points = np.column_stack((xx.ravel(), yy.ravel(), zz.ravel()))
    fields = []
    for time in times:
        values = np.asarray(velocity(points, float(time)), dtype=np.float64)
        if values.shape != points.shape or not np.isfinite(values).all():
            raise ValueError("velocity must return finite shape-(n,3) values")
        fields.append(values.reshape(resolution, resolution, resolution, 3))
    result = np.asarray(fields, dtype=np.float64)
    if not np.any(result):
        raise ValueError("zero sampled field is not renderable")
    return axis, result


def grid_interpolator(axis: np.ndarray, field: np.ndarray) -> Callable[[np.ndarray], np.ndarray]:
    axis = np.asarray(axis, dtype=np.float64)
    field = np.asarray(field, dtype=np.float64)
    if field.shape != (len(axis), len(axis), len(axis), 3):
        raise ValueError("field shape must be (n,n,n,3)")
    components = [
        RegularGridInterpolator(
            (axis, axis, axis),
            field[..., component],
            method="linear",
            bounds_error=False,
            fill_value=0.0,
        )
        for component in range(3)
    ]

    def evaluate(points: np.ndarray) -> np.ndarray:
        pts = np.asarray(points, dtype=np.float64)
        one = pts.ndim == 1
        pts = np.atleast_2d(pts)
        out = np.column_stack([interpolator(pts) for interpolator in components])
        return out[0] if one else out

    return evaluate


def vorticity(field: np.ndarray, spacing: float) -> tuple[np.ndarray, np.ndarray]:
    field = np.asarray(field, dtype=np.float64)
    if field.ndim != 4 or field.shape[-1] != 3:
        raise ValueError("field must have shape (n,n,n,3)")
    if not np.isfinite(field).all() or not (spacing > 0.0):
        raise ValueError("invalid field or spacing")
    u, v, w = (field[..., i] for i in range(3))
    _du_dx, du_dy, du_dz = np.gradient(u, spacing, spacing, spacing, edge_order=2)
    dv_dx, _dv_dy, dv_dz = np.gradient(v, spacing, spacing, spacing, edge_order=2)
    dw_dx, dw_dy, _dw_dz = np.gradient(w, spacing, spacing, spacing, edge_order=2)
    omega = np.stack((dw_dy - dv_dz, du_dz - dw_dx, dv_dx - du_dy), axis=-1)
    return omega, np.linalg.norm(omega, axis=-1)


def _normalized_direction(
    evaluate: Callable[[np.ndarray], np.ndarray], points: np.ndarray, sign: float
) -> np.ndarray:
    values = np.asarray(evaluate(points), dtype=np.float64)
    speeds = np.linalg.norm(values, axis=1)
    out = np.zeros_like(values)
    active = speeds >= SPEED_FLOOR
    out[active] = float(sign) * values[active] / speeds[active, None]
    return out


def _rk4_side(
    evaluate: Callable[[np.ndarray], np.ndarray], seeds: np.ndarray, sign: float
) -> np.ndarray:
    seeds = np.asarray(seeds, dtype=np.float64)
    step = STREAM_ARCLENGTH / float(STREAM_OUTPUTS_PER_SIDE - 1)
    states = seeds.copy()
    history = [states.copy()]
    active = np.max(np.abs(states), axis=1) < 1.995
    for _ in range(STREAM_OUTPUTS_PER_SIDE - 1):
        k1 = _normalized_direction(evaluate, states, sign)
        k2 = _normalized_direction(evaluate, states + 0.5 * step * k1, sign)
        k3 = _normalized_direction(evaluate, states + 0.5 * step * k2, sign)
        k4 = _normalized_direction(evaluate, states + step * k3, sign)
        proposal = states + (step / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        proposal_active = np.max(np.abs(proposal), axis=1) < 1.995
        keep = active & proposal_active
        states[keep] = proposal[keep]
        active = keep
        history.append(states.copy())
    return np.stack(history, axis=1)


def integrate_streamlines(
    evaluate: Callable[[np.ndarray], np.ndarray], seeds: np.ndarray
) -> np.ndarray:
    backward = _rk4_side(evaluate, seeds, -1.0)
    forward = _rk4_side(evaluate, seeds, +1.0)
    return np.concatenate((backward[:, :0:-1, :], forward), axis=1)


def streamline_metrics(lines: np.ndarray) -> dict:
    turns, axial_spans, radial_spans, lengths = [], [], [], []
    for line in np.asarray(lines, dtype=np.float64):
        radius = np.hypot(line[:, 0], line[:, 1])
        angle = np.unwrap(np.arctan2(line[:, 1], line[:, 0]))
        turns.append(float(np.sum(np.abs(np.diff(angle))) / (2.0 * np.pi)))
        axial_spans.append(float(np.ptp(line[:, 2])))
        radial_spans.append(float(np.ptp(radius)))
        lengths.append(float(np.sum(np.linalg.norm(np.diff(line, axis=0), axis=1))))
    return {
        "line_count": int(len(lines)),
        "mean_absolute_turns": float(np.mean(turns)),
        "max_absolute_turns": float(np.max(turns)),
        "mean_axial_span": float(np.mean(axial_spans)),
        "mean_radial_span": float(np.mean(radial_spans)),
        "mean_rendered_length": float(np.mean(lengths)),
    }


def _weighted_radial_rms(
    magnitude: np.ndarray, radius: np.ndarray, mask: np.ndarray
) -> float:
    weight = np.square(magnitude) * mask
    total = float(np.sum(weight))
    if total <= np.finfo(float).tiny:
        raise ValueError("empty vorticity band")
    return float(np.sqrt(np.sum(weight * radius * radius) / total))


def vorticity_metrics(magnitude: np.ndarray, axis: np.ndarray) -> dict:
    threshold = float(np.quantile(magnitude, VORTICITY_QUANTILE))
    selected = magnitude >= threshold
    xx, yy, zz = np.meshgrid(axis, axis, axis, indexing="ij")
    radius = np.hypot(xx, yy)
    abs_z = np.abs(zz)
    normalized_z = abs_z / 2.0
    tip = (normalized_z >= TIP_BAND[0]) & (normalized_z <= TIP_BAND[1])
    central = normalized_z <= CENTRAL_BAND_MAX
    return {
        "quantile": VORTICITY_QUANTILE,
        "threshold": threshold,
        "selected_points": int(np.count_nonzero(selected)),
        "selected_axial_rms": float(np.sqrt(np.mean(np.square(abs_z[selected])))),
        "selected_radial_rms": float(np.sqrt(np.mean(np.square(radius[selected])))),
        "selected_max_abs_z": float(np.max(abs_z[selected])),
        "tip_enstrophy_weighted_radial_rms": _weighted_radial_rms(
            magnitude, radius, tip
        ),
        "central_enstrophy_weighted_radial_rms": _weighted_radial_rms(
            magnitude, radius, central
        ),
        "full_vorticity_rms": float(np.sqrt(np.mean(np.square(magnitude)))),
    }


def _relative(child: float, control: float) -> float:
    if abs(control) <= np.finfo(float).tiny:
        raise ValueError("relative comparison has zero control")
    return float(child / control - 1.0)


def _field_digest(axis: np.ndarray, field: np.ndarray) -> str:
    digest = hashlib.sha256()
    for array in (axis, field):
        value = np.ascontiguousarray(np.asarray(array, dtype="<f8"))
        digest.update(np.asarray(value.shape, dtype="<i8").tobytes())
        digest.update(value.tobytes())
    return digest.hexdigest()


def analyze_pair(axis: np.ndarray, control_fields: np.ndarray, child_fields: np.ndarray) -> dict:
    if control_fields.shape != child_fields.shape:
        raise ValueError("control/child grid shapes differ")
    expected = (len(REFERENCE_TIMES), len(axis), len(axis), len(axis), 3)
    if control_fields.shape != expected:
        raise ValueError("unexpected grid contract")
    spacing = float(axis[1] - axis[0])
    seeds = seed_points()
    rows = []
    for index, time in enumerate(REFERENCE_TIMES):
        row = {"time": float(time)}
        metrics = {}
        for label, fields in (("control", control_fields), ("taper", child_fields)):
            current = fields[index]
            _, omega_mag = vorticity(current, spacing)
            lines = integrate_streamlines(grid_interpolator(axis, current), seeds)
            metrics[label] = {
                "streamlines": streamline_metrics(lines),
                "vorticity": vorticity_metrics(omega_mag, axis),
            }
        row.update(metrics)
        cs = metrics["control"]["streamlines"]
        ts = metrics["taper"]["streamlines"]
        cv = metrics["control"]["vorticity"]
        tv = metrics["taper"]["vorticity"]
        row["relative"] = {
            "mean_absolute_turns": _relative(ts["mean_absolute_turns"], cs["mean_absolute_turns"]),
            "max_absolute_turns": _relative(ts["max_absolute_turns"], cs["max_absolute_turns"]),
            "mean_axial_span": _relative(ts["mean_axial_span"], cs["mean_axial_span"]),
            "mean_radial_span": _relative(ts["mean_radial_span"], cs["mean_radial_span"]),
            "top_vorticity_axial_rms": _relative(tv["selected_axial_rms"], cv["selected_axial_rms"]),
            "top_vorticity_radial_rms": _relative(tv["selected_radial_rms"], cv["selected_radial_rms"]),
            "tip_vorticity_radial_rms": _relative(
                tv["tip_enstrophy_weighted_radial_rms"],
                cv["tip_enstrophy_weighted_radial_rms"],
            ),
            "central_vorticity_radial_rms": _relative(
                tv["central_enstrophy_weighted_radial_rms"],
                cv["central_enstrophy_weighted_radial_rms"],
            ),
            "full_vorticity_rms": _relative(tv["full_vorticity_rms"], cv["full_vorticity_rms"]),
        }
        rows.append(row)
    return {
        "task_id": TASK_ID,
        "contract": {
            "times": list(REFERENCE_TIMES),
            "seed_radii": list(SEED_RADII),
            "seed_z": list(SEED_Z),
            "seed_azimuths": SEED_AZIMUTHS,
            "streamline_count": int(len(seeds)),
            "grid_resolution": int(len(axis)),
            "stream_arclength_each_side": STREAM_ARCLENGTH,
            "stream_outputs_per_side": STREAM_OUTPUTS_PER_SIDE,
            "vorticity_quantile": VORTICITY_QUANTILE,
            "tip_band_abs_z_over_2": list(TIP_BAND),
            "central_band_abs_z_over_2_max": CENTRAL_BAND_MAX,
            "image_or_openai_numeric_target_used": False,
            "visual_acceptance_threshold_used": False,
        },
        "identities": {
            "base_main_sha": BASE_MAIN_SHA,
            "st052_pr": ST052_PR,
            "st052_head": ST052_HEAD,
            "redistribution_pr": REDISTRIBUTION_PR,
            "redistribution_head": REDISTRIBUTION_HEAD,
            "path_control_pr": PATH_CONTROL_PR,
            "taper_pr": TAPER_PR,
            "taper_head": TAPER_HEAD,
            "taper_path_pr": TAPER_PATH_PR,
            "energy_neutral_taper_pr": ENERGY_NEUTRAL_TAPER_PR,
            "energy_neutral_taper_head": ENERGY_NEUTRAL_TAPER_HEAD,
            "scipy_source_commit": SCIPY_SOURCE_COMMIT,
            "matplotlib_release": MATPLOTLIB_RELEASE,
            "matplotlib_source_commit": MATPLOTLIB_SOURCE_COMMIT,
        },
        "control_grid_sha256": _field_digest(axis, control_fields),
        "taper_grid_sha256": _field_digest(axis, child_fields),
        "per_time": rows,
        "truth": dict(TRUTH),
    }


def render_pair(
    axis: np.ndarray,
    control_fields: np.ndarray,
    child_fields: np.ndarray,
    output_dir: str | Path,
) -> list[str]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    seeds = seed_points()
    spacing = float(axis[1] - axis[0])
    xx, yy, zz = np.meshgrid(axis, axis, axis, indexing="ij")
    paths = []
    for index, time in enumerate(REFERENCE_TIMES):
        fig = plt.figure(figsize=(12, 6), constrained_layout=True)
        for panel, (label, fields) in enumerate(
            (("ST052-M + kappa=.05", control_fields), ("+ localized taper tau=.05", child_fields)),
            start=1,
        ):
            ax = fig.add_subplot(1, 2, panel, projection="3d")
            current = fields[index]
            _, omega_mag = vorticity(current, spacing)
            threshold = float(np.quantile(omega_mag, VORTICITY_QUANTILE))
            mask = omega_mag >= threshold
            ax.scatter(
                xx[mask], yy[mask], zz[mask], c=omega_mag[mask],
                s=3, alpha=0.18, cmap="magma", linewidths=0,
            )
            lines = integrate_streamlines(grid_interpolator(axis, current), seeds)
            for line in lines:
                ax.plot(line[:, 0], line[:, 1], line[:, 2], linewidth=0.7, alpha=0.78)
            ax.set_xlim(-2, 2)
            ax.set_ylim(-2, 2)
            ax.set_zlim(-2, 2)
            ax.set_box_aspect((1, 1, 1))
            ax.view_init(elev=18, azim=-58)
            ax.set_title(f"{label}   t={time:.2f}")
            ax.set_xlabel("x")
            ax.set_ylabel("y")
            ax.set_zlabel("z")
        target = output / f"st052m_taper_streamline_vorticity_t{time:.2f}.png"
        fig.savefig(target, dpi=170)
        plt.close(fig)
        paths.append(str(target))
    return paths


def save_report(report: dict, path: str | Path) -> str:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    body = json.dumps(report, sort_keys=True, indent=2, allow_nan=False) + "\n"
    target.write_text(body)
    return hashlib.sha256(body.encode("utf-8")).hexdigest()


def load_exact_source(source_root: str | Path) -> tuple[Velocity, Velocity, dict]:
    root = Path(source_root).resolve()
    source_dir = root / "experiments" / "root_st052"
    if not source_dir.is_dir():
        raise ValueError("source root lacks experiments/root_st052")
    sys.path.insert(0, str(source_dir))
    try:
        taper = importlib.import_module("agent7_st052m_localized_taper")
        base = importlib.import_module("agent7_st052m_frozen_redistribution")
        replay = importlib.import_module("replay_st052")
    finally:
        sys.path.pop(0)

    if taper.PARENT_HEAD != REDISTRIBUTION_HEAD:
        raise ValueError("localized-taper parent identity drift")
    if taper.TAPER_TAU != 0.05 or tuple(taper.WINDOW) != (0.50, 0.82):
        raise ValueError("localized-taper transform identity drift")
    if abs(float(base.GAIN) - 0.05) > 1e-15:
        raise ValueError("redistribution gain drift")
    if abs(float(base.SOURCE_ALPHA) - 2.520520814687742) > 1e-14:
        raise ValueError("redistribution alpha drift")

    field, raw = replay.reconstruct()
    redistribution_scale, _, _ = base.child_scale(field, raw)
    taper_scale, control_energy, raw_taper_energy = taper.taper_energy_scale(
        field, raw, redistribution_scale
    )
    if abs(redistribution_scale - EXPECTED_REDISTRIBUTION_SCALE) > 5e-10:
        raise ValueError("redistribution normalization drift")
    if abs(taper_scale - EXPECTED_TAPER_SCALE) > 5e-10:
        raise ValueError("localized-taper normalization drift")

    def control(points: np.ndarray, time: float) -> np.ndarray:
        return taper.control_velocity(field, raw, points, time, redistribution_scale)

    def child(points: np.ndarray, time: float) -> np.ndarray:
        return taper.taper_velocity(
            field,
            raw,
            points,
            time,
            redistribution_scale=redistribution_scale,
            taper_scale=taper_scale,
        )

    source = {
        "source_checkout_expected_head": TAPER_HEAD,
        "agent7_task_id": taper.TASK_ID,
        "parent_head": taper.PARENT_HEAD,
        "taper_tau": float(taper.TAPER_TAU),
        "taper_window_abs_z_over_2": list(taper.WINDOW),
        "redistribution_gain": float(base.GAIN),
        "redistribution_alpha": float(base.SOURCE_ALPHA),
        "redistribution_scale": float(redistribution_scale),
        "taper_scale": float(taper_scale),
        "control_reference_energy": float(control_energy),
        "raw_taper_reference_energy": float(raw_taper_energy),
    }
    return control, child, source


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args(argv)

    control, child, source = load_exact_source(args.source_root)
    control_axis, control_fields = sample_velocity_grid(control)
    child_axis, child_fields = sample_velocity_grid(child)
    if not np.array_equal(control_axis, child_axis):
        raise ValueError("sampled axes differ")

    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output / "control_grid.npz", axis=control_axis, velocity=control_fields)
    np.savez_compressed(output / "localized_taper_grid.npz", axis=child_axis, velocity=child_fields)

    report = analyze_pair(control_axis, control_fields, child_fields)
    report["source_replay"] = source
    report["render_files"] = [
        Path(path).name
        for path in render_pair(control_axis, control_fields, child_fields, output)
    ]
    report_sha = save_report(report, output / "report.json")
    print(
        json.dumps(
            {
                "report_sha256": report_sha,
                "source_replay": source,
                "relative": [row["relative"] for row in report["per_time"]],
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
