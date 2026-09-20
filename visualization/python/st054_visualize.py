#!/usr/bin/env python3
"""Render the published continuous ST054 velocity snapshot in Python.

This is a delivery/visualization adapter around the checksum-bound ST054 field
already published on ``main``.  It does not fit a new velocity, compare against
an OpenAI image, or evaluate Navier--Stokes acceptance.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import numpy as np

from openai_ns_reconstruction import available_st054_models, load_st054

SCHEMA = "st054_python_visualization_smoke_v1"
SEED_RADII = (0.6, 0.9, 1.2)
SEED_Z = (-0.3, 0.3)
SEED_AZIMUTHS = 8
STREAMLINE_STEP = 0.025
STREAMLINE_STEPS_EACH_DIRECTION = 90
STREAMLINE_SPEED_FLOOR = 1.0e-12
VORTICITY_GRID_SIZE = 65
VORTICITY_HALF_WIDTH = 1.6

# Screened external/public implementations.  Only Matplotlib public APIs are
# directly used here; no third-party source code is copied into this repository.
EXTERNAL_METHODS = (
    {
        "repo": "matplotlib/matplotlib",
        "screened_commit": "a0acad0712a8571509ce37b864e651255ff9a4fe",
        "license": "Matplotlib License Agreement",
        "classification": "direct_migration_public_api_only",
        "scope": "headless 2-D/3-D rendering and Line3DCollection",
    },
    {
        "repo": "scipy/scipy",
        "screened_commit": "2c2998c8d9d5fc63d3e6e3d330db46ee476379b2",
        "license": "BSD-3-Clause",
        "classification": "screened_not_adopted_for_streamline_integrator",
        "scope": "solve_ivp/DOP853 considered; vectorized frozen-step RK4 retained for deterministic batch streamline geometry",
    },
    {
        "repo": "Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1",
        "screened_commit": "ee80de11775bd45157d114000378b6ee3422ff46",
        "license": "repository-native",
        "classification": "suitable_reimplementation_internal_protocol",
        "scope": "fixed seed rings and normalized-velocity bidirectional streamline protocol from CR-A9-064",
    },
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fixed_seeds() -> np.ndarray:
    """Return the frozen 48 support-interior streamline seeds."""
    rows: list[tuple[float, float, float]] = []
    for radius in SEED_RADII:
        for z in SEED_Z:
            for index in range(SEED_AZIMUTHS):
                phi = 2.0 * np.pi * index / SEED_AZIMUTHS
                rows.append((radius * np.cos(phi), radius * np.sin(phi), z))
    seeds = np.asarray(rows, dtype=float)
    if seeds.shape != (48, 3):
        raise RuntimeError("Frozen ST054 seed contract drift")
    return seeds


def _inside_support(field: Any, points: np.ndarray) -> np.ndarray:
    radius = np.hypot(points[:, 0], points[:, 1])
    return (
        np.isfinite(points).all(axis=1)
        & (radius < float(field.support_radius))
        & (np.abs(points[:, 2]) < float(field.support_z))
    )


def _unit_velocity(field: Any, points: np.ndarray, time: float, direction: float) -> tuple[np.ndarray, np.ndarray]:
    velocity = np.asarray(
        field.velocity(points[:, 0], points[:, 1], points[:, 2], float(time)),
        dtype=float,
    )
    if velocity.shape != points.shape:
        raise RuntimeError("Malformed ST054 velocity shape during streamline integration")
    speed = np.linalg.norm(velocity, axis=1)
    valid = _inside_support(field, points) & np.isfinite(speed) & (speed > STREAMLINE_SPEED_FLOOR)
    tangent = np.zeros_like(velocity)
    tangent[valid] = float(direction) * velocity[valid] / speed[valid, None]
    return tangent, valid


def _trace_direction(field: Any, seeds: np.ndarray, time: float, direction: float) -> np.ndarray:
    """Vectorized fixed-arclength RK4 on the instantaneous normalized field."""
    positions = np.asarray(seeds, dtype=float).copy()
    active = _inside_support(field, positions)
    states = np.full(
        (STREAMLINE_STEPS_EACH_DIRECTION + 1, positions.shape[0], 3),
        np.nan,
        dtype=float,
    )
    states[0, active] = positions[active]

    h = STREAMLINE_STEP
    for step in range(STREAMLINE_STEPS_EACH_DIRECTION):
        k1, ok1 = _unit_velocity(field, positions, time, direction)
        k2, ok2 = _unit_velocity(field, positions + 0.5 * h * k1, time, direction)
        k3, ok3 = _unit_velocity(field, positions + 0.5 * h * k2, time, direction)
        k4, ok4 = _unit_velocity(field, positions + h * k3, time, direction)
        next_positions = positions + h * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
        next_active = active & ok1 & ok2 & ok3 & ok4 & _inside_support(field, next_positions)
        positions[next_active] = next_positions[next_active]
        active = next_active
        states[step + 1, active] = positions[active]
        if not active.any():
            break
    return states


def trace_streamlines(field: Any, time: float) -> list[np.ndarray]:
    seeds = fixed_seeds()
    backward = _trace_direction(field, seeds, time, -1.0)
    forward = _trace_direction(field, seeds, time, +1.0)
    paths: list[np.ndarray] = []
    for index in range(seeds.shape[0]):
        back = backward[:, index]
        back = back[np.isfinite(back).all(axis=1)]
        front = forward[:, index]
        front = front[np.isfinite(front).all(axis=1)]
        if back.shape[0] == 0 or front.shape[0] == 0:
            continue
        path = np.vstack((back[::-1], front[1:]))
        if path.shape[0] >= 2:
            paths.append(path)
    if not paths:
        raise RuntimeError("No nonempty ST054 streamlines were produced")
    return paths


def vorticity_midplane(field: Any, time: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Centered finite-difference omega_z on z=0 for visualization only."""
    axis = np.linspace(-VORTICITY_HALF_WIDTH, VORTICITY_HALF_WIDTH, VORTICITY_GRID_SIZE)
    x, y = np.meshgrid(axis, axis, indexing="xy")
    z = np.zeros_like(x)
    velocity = np.asarray(field.velocity(x, y, z, float(time)), dtype=float)
    if velocity.shape != x.shape + (3,) or not np.isfinite(velocity).all():
        raise RuntimeError("Malformed/non-finite ST054 mid-plane velocity")
    dv_dx = np.gradient(velocity[..., 1], axis, axis=1, edge_order=2)
    du_dy = np.gradient(velocity[..., 0], axis, axis=0, edge_order=2)
    omega_z = dv_dx - du_dy
    if not np.isfinite(omega_z).all():
        raise RuntimeError("Non-finite ST054 mid-plane vorticity")
    return x, y, omega_z


def _streamline_segments(field: Any, paths: list[np.ndarray], time: float) -> tuple[np.ndarray, np.ndarray]:
    segments = np.concatenate(
        [np.stack((path[:-1], path[1:]), axis=1) for path in paths if path.shape[0] >= 2],
        axis=0,
    )
    midpoints = segments.mean(axis=1)
    velocity = np.asarray(
        field.velocity(midpoints[:, 0], midpoints[:, 1], midpoints[:, 2], float(time)),
        dtype=float,
    )
    speed = np.linalg.norm(velocity, axis=1)
    if not (segments.shape[0] > 0 and np.isfinite(speed).all() and np.max(speed) > 0.0):
        raise RuntimeError("Degenerate ST054 streamline rendering payload")
    return segments, speed


def render(model_id: str, time: float, output: Path, receipt_path: Path | None = None) -> dict[str, Any]:
    field = load_st054(model_id)
    if not (field.tmin <= float(time) <= field.tmax):
        raise ValueError(f"time must lie in [{field.tmin}, {field.tmax}]")
    if (field.nu, field.support_radius, field.support_z) != (0.01, 2.0, 2.0):
        raise RuntimeError("Published ST054 physical metadata drift")

    paths = trace_streamlines(field, float(time))
    segments, speed = _streamline_segments(field, paths, float(time))
    x, y, omega_z = vorticity_midplane(field, float(time))

    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(12.0, 5.4), constrained_layout=True)

    ax3d = fig.add_subplot(1, 2, 1, projection="3d")
    speed_min = float(np.min(speed))
    speed_max = float(np.max(speed))
    if speed_max <= speed_min:
        speed_max = speed_min + np.finfo(float).eps
    line_collection = Line3DCollection(
        segments,
        cmap="viridis",
        norm=Normalize(vmin=speed_min, vmax=speed_max),
        linewidth=0.75,
        alpha=0.92,
    )
    line_collection.set_array(speed)
    ax3d.add_collection3d(line_collection)
    limit = 1.9
    ax3d.set(xlim=(-limit, limit), ylim=(-limit, limit), zlim=(-limit, limit))
    ax3d.set_xlabel("x")
    ax3d.set_ylabel("y")
    ax3d.set_zlabel("z")
    ax3d.set_box_aspect((1.0, 1.0, 1.15))
    ax3d.view_init(elev=18.0, azim=-58.0)
    ax3d.set_title(f"{model_id} instantaneous streamlines, t={float(time):.3f}")
    colorbar = fig.colorbar(line_collection, ax=ax3d, shrink=0.72, pad=0.08)
    colorbar.set_label("speed")

    ax2d = fig.add_subplot(1, 2, 2)
    vort_limit = float(np.max(np.abs(omega_z)))
    if vort_limit == 0.0:
        vort_limit = 1.0
    mesh = ax2d.pcolormesh(
        x,
        y,
        omega_z,
        shading="auto",
        cmap="coolwarm",
        vmin=-vort_limit,
        vmax=vort_limit,
    )
    ax2d.set_aspect("equal", adjustable="box")
    ax2d.set_xlabel("x")
    ax2d.set_ylabel("y")
    ax2d.set_title(r"mid-plane $\omega_z=\partial_x v-\partial_y u$ (visual diagnostic)")
    vort_colorbar = fig.colorbar(mesh, ax=ax2d, shrink=0.86)
    vort_colorbar.set_label(r"$\omega_z$")

    fig.suptitle("Published ST054 snapshot — software visualization only; not PDE or OpenAI-image acceptance")
    fig.savefig(output, dpi=150)
    plt.close(fig)

    if not output.is_file() or output.stat().st_size <= 0:
        raise RuntimeError("Python ST054 visualization did not produce a nonempty PNG")

    receipt: dict[str, Any] = {
        "schema": SCHEMA,
        "model_id": model_id,
        "time": float(time),
        "source_commit": field.source_commit,
        "source_mat_sha256": field.mat_sha256,
        "nu": float(field.nu),
        "support_radius": float(field.support_radius),
        "support_z": float(field.support_z),
        "time_interval": [float(field.tmin), float(field.tmax)],
        "streamline_protocol": {
            "seed_count": int(fixed_seeds().shape[0]),
            "radii": list(SEED_RADII),
            "z": list(SEED_Z),
            "azimuth_count": SEED_AZIMUTHS,
            "integration": "bidirectional_vectorized_normalized_velocity_RK4",
            "step": STREAMLINE_STEP,
            "steps_each_direction": STREAMLINE_STEPS_EACH_DIRECTION,
            "speed_floor": STREAMLINE_SPEED_FLOOR,
            "camera_elev": 18.0,
            "camera_azim": -58.0,
        },
        "vorticity_protocol": {
            "plane_z": 0.0,
            "grid_size": VORTICITY_GRID_SIZE,
            "half_width": VORTICITY_HALF_WIDTH,
            "operator": "numpy_gradient_centered_second_order_interior",
        },
        "measurements": {
            "nonempty_streamlines": len(paths),
            "line_segments": int(segments.shape[0]),
            "speed_min": float(np.min(speed)),
            "speed_max": float(np.max(speed)),
            "omega_z_min": float(np.min(omega_z)),
            "omega_z_max": float(np.max(omega_z)),
            "omega_z_rms": float(np.sqrt(np.mean(omega_z * omega_z))),
        },
        "png": {
            "filename": output.name,
            "bytes": int(output.stat().st_size),
            "sha256": _sha256(output),
        },
        "external_method_screen": list(EXTERNAL_METHODS),
        "python_visualization_smoke_ready": True,
        "visualization_ready": False,
        "visual_correspondence_verified": False,
        "source_correspondence_verified": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
        "engineering_scope": "continuous_saved_field_python_render_smoke_only_not_scientific_acceptance",
    }
    if receipt_path is not None:
        receipt_path = Path(receipt_path)
        receipt_path.parent.mkdir(parents=True, exist_ok=True)
        receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", choices=available_st054_models(), default="ST054-Q2")
    parser.add_argument("--time", type=float, default=0.5)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--receipt", type=Path)
    args = parser.parse_args()
    receipt = render(args.model, args.time, args.output, args.receipt)
    print(json.dumps(receipt, sort_keys=True))


if __name__ == "__main__":
    main()
