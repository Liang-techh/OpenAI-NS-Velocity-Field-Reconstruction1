"""Deterministic fixed-seed 3-D streamline/vorticity comparison for ST051-B.

This module consumes the immutable CR-A9-053 child grid, exactly inverts its
frozen swirl redistribution to recover the sampled ST051-B parent, and performs
only target-free visualization diagnostics. It does not define a visual
acceptance threshold and does not alter or validate the PDE.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Callable

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.io import netcdf_file

TASK_ID = "CR-A9-054"
BASE_MAIN_SHA = "f0193d66c9d92948b4820ebcb70263673995b324"
GRID_SOURCE_PR = 476
GRID_SOURCE_HEAD = "d068e0d4c2cdce23a08aec791af8d2d59e459e65"
GRID_SOURCE_RUN = 35389278384
GRID_SOURCE_ARTIFACT_ID = 10565335507
GRID_SOURCE_SHA256 = "da9e3aa9696f12bfea8e9a3be86b024a493a2242a4c4cd2644e454172a0b89cc"
GRID_SOURCE_RECEIPT_SHA256 = "5738d8156be5fbc5d5d85fa198cbb35d97ac1e3a0810e5e51f29c0b7caeae275"
PARENT_PR = 460
PARENT_HEAD = "4b784f1b8457af2ead49295631d834d4e882000b"
TRANSFER_PR = 469
TRANSFER_HEAD = "95f94188270f02511c8bdf660822f77530c9983a"
PATH_REFERENCE_PR = 471
HEADROOM_PR = 480
HEADROOM_HEAD = "187157d377b14bae56415648ac100342f2b8bbb0"
SOURCE_ALPHA = 2.520520814687742
GAIN = 0.025
INNER_WINDOW = (0.30, 1.05)
OUTER_WINDOW = (0.95, 1.85)
EXPECTED_REFERENCE_SCALE = 1.0014791925672812
REFERENCE_TIMES = (0.25, 0.50, 0.75)
SEED_RADII = (0.6, 0.9, 1.2)
SEED_Z = (-0.3, 0.3)
SEED_AZIMUTHS = 8
GRID_RESOLUTION = 33
STREAM_ARCLENGTH = 2.25
STREAM_OUTPUTS_PER_SIDE = 91
SPEED_FLOOR = 1e-8
VORTICITY_QUANTILE = 0.985
SCIPY_SOURCE_COMMIT = "b12c772edbc1fe0d3db9481cdfcc2e311569cb25"
MATPLOTLIB_RELEASE = "3.10.6"
MATPLOTLIB_SOURCE_COMMIT = "5cd38c3edcdf0792d0e6aded280a9b7a7de6146f"

Velocity = Callable[[np.ndarray, float], np.ndarray]


def compact_bump(radius: np.ndarray, lo: float, hi: float) -> np.ndarray:
    radius = np.asarray(radius, dtype=np.float64)
    mid = 0.5 * (lo + hi)
    half = 0.5 * (hi - lo)
    coordinate = (radius - mid) / half
    result = np.zeros_like(radius)
    mask = np.abs(coordinate) < 1.0
    inside = coordinate[mask]
    result[mask] = np.exp(1.0 - 1.0 / (1.0 - inside * inside))
    return result


def redistribution_profile(radius: np.ndarray) -> np.ndarray:
    return compact_bump(radius, *INNER_WINDOW) - SOURCE_ALPHA * compact_bump(radius, *OUTER_WINDOW)


def load_child_grid_netcdf(path: str | Path) -> tuple[np.ndarray, np.ndarray]:
    """Load the exact CR-A9-053 NetCDF layout into `(axis, velocity)` arrays."""
    with netcdf_file(str(path), mode="r", mmap=False) as nc:
        times = np.asarray(nc.variables["time"][:], dtype=np.float64).copy()
        x = np.asarray(nc.variables["x"][:], dtype=np.float64).copy()
        y = np.asarray(nc.variables["y"][:], dtype=np.float64).copy()
        z = np.asarray(nc.variables["z"][:], dtype=np.float64).copy()
        components = [
            np.asarray(nc.variables[name][:], dtype=np.float64).copy()
            for name in ("u", "v", "w")
        ]
    if not np.array_equal(times, np.asarray(REFERENCE_TIMES, dtype=np.float64)):
        raise ValueError("source artifact time contract drift")
    if not (np.array_equal(x, y) and np.array_equal(x, z)):
        raise ValueError("source artifact Cartesian axes differ")
    if len(x) != GRID_RESOLUTION:
        raise ValueError("source artifact resolution drift")
    velocity = np.stack(components, axis=-1)
    expected_shape = (len(REFERENCE_TIMES), len(x), len(x), len(x), 3)
    if velocity.shape != expected_shape or not np.isfinite(velocity).all() or not np.any(velocity):
        raise ValueError("source artifact velocity grid is invalid")
    return x, velocity


def _transform_grid(
    axis: np.ndarray,
    fields: np.ndarray,
    *,
    inverse: bool,
    gain: float = GAIN,
    scale: float = EXPECTED_REFERENCE_SCALE,
) -> np.ndarray:
    axis = np.asarray(axis, dtype=np.float64)
    fields = np.asarray(fields, dtype=np.float64)
    expected = (len(REFERENCE_TIMES), len(axis), len(axis), len(axis), 3)
    if fields.shape != expected or not np.isfinite(fields).all():
        raise ValueError("unexpected velocity-grid shape")
    if not (np.isfinite(scale) and scale > 0.0):
        raise ValueError("scale must be positive and finite")
    xx, yy, _zz = np.meshgrid(axis, axis, axis, indexing="ij")
    radius = np.hypot(xx, yy)
    factor = 1.0 + float(gain) * redistribution_profile(radius)
    if not np.all(np.isfinite(factor)) or np.min(factor) <= 0.0:
        raise ValueError("swirl redistribution factor must stay positive")

    out = fields.copy()
    radial_mask = radius > 1e-14
    rx = np.zeros_like(radius)
    ry = np.zeros_like(radius)
    rx[radial_mask] = xx[radial_mask] / radius[radial_mask]
    ry[radial_mask] = yy[radial_mask] / radius[radial_mask]

    for index in range(len(REFERENCE_TIMES)):
        ux = fields[index, ..., 0]
        uy = fields[index, ..., 1]
        uz = fields[index, ..., 2]
        ur = rx * ux + ry * uy
        ut = -ry * ux + rx * uy
        if inverse:
            ur = ur / scale
            ut = ut / (scale * factor)
            uz_new = uz / scale
        else:
            ur = scale * ur
            ut = scale * factor * ut
            uz_new = scale * uz
        ux_new = rx * ur - ry * ut
        uy_new = ry * ur + rx * ut
        # The cylindrical decomposition is undefined at r=0; the frozen
        # transform changes no horizontal direction there and applies only scale.
        ux_new[~radial_mask] = ux[~radial_mask] / scale if inverse else scale * ux[~radial_mask]
        uy_new[~radial_mask] = uy[~radial_mask] / scale if inverse else scale * uy[~radial_mask]
        out[index, ..., 0] = ux_new
        out[index, ..., 1] = uy_new
        out[index, ..., 2] = uz_new
    return out


def recover_parent_grid(axis: np.ndarray, child_fields: np.ndarray) -> np.ndarray:
    """Exactly invert the frozen `.025` representation on sampled grid values."""
    return _transform_grid(axis, child_fields, inverse=True)


def apply_frozen_child_transform(axis: np.ndarray, parent_fields: np.ndarray) -> np.ndarray:
    """Reapply the frozen transform, used as a fail-closed round-trip check."""
    return _transform_grid(axis, parent_fields, inverse=False)


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
    field = np.asarray(fields, dtype=np.float64)
    if not np.any(field):
        raise ValueError("zero sampled field is not renderable")
    return axis, field


def grid_interpolator(axis: np.ndarray, field: np.ndarray) -> Callable[[np.ndarray], np.ndarray]:
    axis = np.asarray(axis, dtype=np.float64)
    field = np.asarray(field, dtype=np.float64)
    if field.shape != (len(axis), len(axis), len(axis), 3):
        raise ValueError("field shape must be (n,n,n,3)")
    components = [
        RegularGridInterpolator(
            (axis, axis, axis), field[..., component], method="linear",
            bounds_error=False, fill_value=0.0,
        )
        for component in range(3)
    ]

    def evaluate(points: np.ndarray) -> np.ndarray:
        pts = np.asarray(points, dtype=np.float64)
        one = pts.ndim == 1
        pts = np.atleast_2d(pts)
        result = np.column_stack([interpolator(pts) for interpolator in components])
        return result[0] if one else result

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
    magnitude = np.linalg.norm(omega, axis=-1)
    return omega, magnitude


def _normalized_direction(evaluate: Callable[[np.ndarray], np.ndarray], points: np.ndarray, sign: float) -> np.ndarray:
    values = np.asarray(evaluate(points), dtype=np.float64)
    speeds = np.linalg.norm(values, axis=1)
    out = np.zeros_like(values)
    active = speeds >= SPEED_FLOOR
    out[active] = float(sign) * values[active] / speeds[active, None]
    return out


def _rk4_side(evaluate: Callable[[np.ndarray], np.ndarray], seeds: np.ndarray, sign: float) -> np.ndarray:
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


def integrate_streamlines(evaluate: Callable[[np.ndarray], np.ndarray], seeds: np.ndarray) -> np.ndarray:
    backward = _rk4_side(evaluate, seeds, -1.0)
    forward = _rk4_side(evaluate, seeds, +1.0)
    return np.concatenate((backward[:, :0:-1, :], forward), axis=1)


def streamline_metrics(lines: np.ndarray) -> dict:
    turns, axial_spans, radial_spans, lengths = [], [], [], []
    for line in lines:
        radius = np.hypot(line[:, 0], line[:, 1])
        angle = np.unwrap(np.arctan2(line[:, 1], line[:, 0]))
        turns.append(float(np.sum(np.abs(np.diff(angle))) / (2.0 * np.pi)))
        axial_spans.append(float(np.ptp(line[:, 2])))
        radial_spans.append(float(np.ptp(radius)))
        lengths.append(float(np.sum(np.linalg.norm(np.diff(line, axis=0), axis=1))))
    return {
        "line_count": len(lines),
        "mean_absolute_turns": float(np.mean(turns)),
        "max_absolute_turns": float(np.max(turns)),
        "mean_axial_span": float(np.mean(axial_spans)),
        "mean_radial_span": float(np.mean(radial_spans)),
        "mean_rendered_length": float(np.mean(lengths)),
    }


def _vorticity_metrics(magnitude: np.ndarray, axis: np.ndarray) -> dict:
    threshold = float(np.quantile(magnitude, VORTICITY_QUANTILE))
    mask = magnitude >= threshold
    xx, yy, zz = np.meshgrid(axis, axis, axis, indexing="ij")
    radius = np.hypot(xx[mask], yy[mask])
    axial = np.abs(zz[mask])
    return {
        "quantile": VORTICITY_QUANTILE,
        "threshold": threshold,
        "selected_points": int(np.count_nonzero(mask)),
        "selected_axial_rms": float(np.sqrt(np.mean(axial * axial))),
        "selected_radial_rms": float(np.sqrt(np.mean(radius * radius))),
        "selected_max_abs_z": float(np.max(axial)),
    }


def _field_digest(axis: np.ndarray, field: np.ndarray) -> str:
    h = hashlib.sha256()
    for array in (axis, field):
        value = np.ascontiguousarray(np.asarray(array, dtype="<f8"))
        h.update(np.asarray(value.shape, dtype="<i8").tobytes())
        h.update(value.tobytes())
    return h.hexdigest()


def analyze_pair(parent_axis: np.ndarray, parent_fields: np.ndarray, child_axis: np.ndarray, child_fields: np.ndarray) -> dict:
    if not np.array_equal(parent_axis, child_axis):
        raise ValueError("parent/child grid axes differ")
    if parent_fields.shape != child_fields.shape:
        raise ValueError("parent/child grid shapes differ")
    if parent_fields.shape[0] != len(REFERENCE_TIMES):
        raise ValueError("unexpected time count")
    spacing = float(parent_axis[1] - parent_axis[0])
    seeds = seed_points()
    per_time = []
    for index, time in enumerate(REFERENCE_TIMES):
        row = {"time": time}
        for label, fields in (("parent", parent_fields), ("child", child_fields)):
            current = fields[index]
            _, omega_mag = vorticity(current, spacing)
            lines = integrate_streamlines(grid_interpolator(parent_axis, current), seeds)
            row[label] = {
                "streamlines": streamline_metrics(lines),
                "vorticity": _vorticity_metrics(omega_mag, parent_axis),
            }
        per_time.append(row)
    return {
        "task_id": TASK_ID,
        "contract": {
            "times": list(REFERENCE_TIMES),
            "seed_radii": list(SEED_RADII),
            "seed_z": list(SEED_Z),
            "seed_azimuths": SEED_AZIMUTHS,
            "streamline_count": int(len(seeds)),
            "grid_resolution": int(len(parent_axis)),
            "stream_arclength_each_side": STREAM_ARCLENGTH,
            "stream_outputs_per_side": STREAM_OUTPUTS_PER_SIDE,
            "vorticity_quantile": VORTICITY_QUANTILE,
            "image_or_openai_numeric_target_used": False,
        },
        "identities": {
            "base_main_sha": BASE_MAIN_SHA,
            "grid_source_pr": GRID_SOURCE_PR,
            "grid_source_head": GRID_SOURCE_HEAD,
            "grid_source_run": GRID_SOURCE_RUN,
            "grid_source_artifact_id": GRID_SOURCE_ARTIFACT_ID,
            "grid_source_sha256": GRID_SOURCE_SHA256,
            "grid_source_receipt_sha256": GRID_SOURCE_RECEIPT_SHA256,
            "parent_pr": PARENT_PR,
            "parent_head": PARENT_HEAD,
            "transfer_pr": TRANSFER_PR,
            "transfer_head": TRANSFER_HEAD,
            "path_reference_pr": PATH_REFERENCE_PR,
            "headroom_pr": HEADROOM_PR,
            "headroom_head": HEADROOM_HEAD,
            "source_alpha": SOURCE_ALPHA,
            "gain": GAIN,
            "common_reference_energy_scale": EXPECTED_REFERENCE_SCALE,
            "scipy_source_commit": SCIPY_SOURCE_COMMIT,
            "matplotlib_release": MATPLOTLIB_RELEASE,
            "matplotlib_source_commit": MATPLOTLIB_SOURCE_COMMIT,
        },
        "parent_grid_sha256": _field_digest(parent_axis, parent_fields),
        "child_grid_sha256": _field_digest(child_axis, child_fields),
        "per_time": per_time,
        "truth": {
            "production_candidate_selected": False,
            "visualization_ready": False,
            "visual_correspondence_verified": False,
            "pde_validated": False,
            "source_correspondence_verified": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
    }


def render_pair(axis: np.ndarray, parent_fields: np.ndarray, child_fields: np.ndarray, output_dir: str | Path) -> list[str]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    seeds = seed_points()
    spacing = float(axis[1] - axis[0])
    paths = []
    xx, yy, zz = np.meshgrid(axis, axis, axis, indexing="ij")
    for index, time in enumerate(REFERENCE_TIMES):
        fig = plt.figure(figsize=(12, 6), constrained_layout=True)
        for panel, (label, fields) in enumerate(
            (("ST051-B parent", parent_fields), ("frozen redistribution .025", child_fields)), start=1
        ):
            ax = fig.add_subplot(1, 2, panel, projection="3d")
            current = fields[index]
            _, omega_mag = vorticity(current, spacing)
            threshold = float(np.quantile(omega_mag, VORTICITY_QUANTILE))
            mask = omega_mag >= threshold
            ax.scatter(xx[mask], yy[mask], zz[mask], c=omega_mag[mask], s=3, alpha=0.18, cmap="magma", linewidths=0)
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
        target = output / f"st051b_streamline_vorticity_t{time:.2f}.png"
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


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--child-netcdf", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args(argv)
    axis, child_fields = load_child_grid_netcdf(args.child_netcdf)
    parent_fields = recover_parent_grid(axis, child_fields)
    reapplied = apply_frozen_child_transform(axis, parent_fields)
    roundtrip_max = float(np.max(np.abs(reapplied - child_fields)))
    if roundtrip_max > 5e-13:
        raise ValueError(f"frozen transform round-trip drift: {roundtrip_max}")
    report = analyze_pair(axis, parent_fields, axis, child_fields)
    report["inverse_transform_roundtrip_max_abs"] = roundtrip_max
    report["render_files"] = [Path(path).name for path in render_pair(axis, parent_fields, child_fields, args.output_dir)]
    report_sha = save_report(report, Path(args.output_dir) / "report.json")
    print(json.dumps({"report_sha256": report_sha, "roundtrip_max_abs": roundtrip_max}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
