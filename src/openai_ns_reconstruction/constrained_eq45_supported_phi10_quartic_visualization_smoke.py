"""Minimal truth-bounded visualization smoke for the quartic supported Eq45 field.

This module deliberately does not score or fit a public image.  It consumes only
an already materialized
``Eq45SupportedPhi10QuarticDerivativeBalancedTemporalCandidate`` through its
public Cartesian velocity interface and writes fixed-frame meridional slices that
can be inspected in Python or loaded as raw arrays by downstream tools.

The rendered quantity is speed on the physical ``y=0`` plane with projected
``(u,w)`` streamlines.  The azimuthal component ``v`` is preserved in the saved
NPZ data even though it is not part of the projected streamline direction.
Generating this smoke is evidence that a candidate is renderable/exportable; it
is not evidence of OpenAI visual correspondence, PDE validity, paper exactness,
or blow-up.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable

import numpy as np

from .constrained_eq45_supported_collar_vorticity import governed_supported_seed
from .constrained_eq45_supported_phi10_quartic_temporal_mode import (
    Eq45SupportedPhi10QuarticDerivativeBalancedTemporalCandidate,
)


SCHEMA = "eq45_supported_phi10_quartic_visualization_smoke_v1"
DEFAULT_TIMES = (0.25, 0.50, 0.75)
DEFAULT_EXTENT = 2.0
DEFAULT_GRID_SIZE = 81

_TRUTH_BOUNDARY = {
    "velocity_changed": False,
    "public_image_used": False,
    "camera_or_registration_fitted": False,
    "visual_pass_threshold_defined": False,
    "visualization_smoke_generated": True,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def governed_quartic_candidate(
) -> Eq45SupportedPhi10QuarticDerivativeBalancedTemporalCandidate:
    """Return the frozen quartic visualization/research candidate."""
    return Eq45SupportedPhi10QuarticDerivativeBalancedTemporalCandidate(
        base=governed_supported_seed()
    )


def _validated_grid_size(grid_size: int) -> int:
    if isinstance(grid_size, bool) or not isinstance(grid_size, (int, np.integer)):
        raise ValueError("grid_size must be an integer")
    value = int(grid_size)
    if value < 9 or value % 2 == 0:
        raise ValueError("grid_size must be an odd integer >= 9")
    return value


def _validated_extent(extent: float) -> float:
    value = float(extent)
    if not np.isfinite(value) or value <= 0.0:
        raise ValueError("extent must be finite and positive")
    return value


def _validated_times(
    candidate: Eq45SupportedPhi10QuarticDerivativeBalancedTemporalCandidate,
    times: Iterable[float],
) -> tuple[float, ...]:
    values = tuple(float(value) for value in times)
    if not values:
        raise ValueError("times must contain at least one value")
    if not all(np.isfinite(value) for value in values):
        raise ValueError("times must be finite")
    if len(set(values)) != len(values):
        raise ValueError("times must not contain duplicates")
    if any(value < candidate.time_start or value > candidate.time_end for value in values):
        raise ValueError(
            f"times must lie in [{candidate.time_start}, {candidate.time_end}]"
        )
    return values


def meridional_velocity_slice(
    candidate: Eq45SupportedPhi10QuarticDerivativeBalancedTemporalCandidate,
    time: float,
    *,
    grid_size: int = DEFAULT_GRID_SIZE,
    extent: float = DEFAULT_EXTENT,
) -> dict[str, np.ndarray | float | str]:
    """Sample the public quartic ``[u,v,w]`` on the fixed physical ``y=0`` plane."""
    if not isinstance(candidate, Eq45SupportedPhi10QuarticDerivativeBalancedTemporalCandidate):
        raise TypeError("candidate must be the quartic supported Eq45 temporal candidate")
    size = _validated_grid_size(grid_size)
    radius = _validated_extent(extent)
    scalar_time = _validated_times(candidate, (time,))[0]

    x = np.linspace(-radius, radius, size, dtype=float)
    z = np.linspace(-radius, radius, size, dtype=float)
    xx, zz = np.meshgrid(x, z, indexing="xy")
    yy = np.zeros_like(xx)
    velocity = candidate.velocity_xyz(xx, yy, zz, scalar_time)
    if velocity.shape != (size, size, 3) or not np.all(np.isfinite(velocity)):
        raise RuntimeError("quartic public velocity returned an invalid meridional grid")

    u = np.asarray(velocity[..., 0], dtype=float)
    v = np.asarray(velocity[..., 1], dtype=float)
    w = np.asarray(velocity[..., 2], dtype=float)
    speed = np.sqrt(u * u + v * v + w * w)
    poloidal_speed = np.sqrt(u * u + w * w)
    return {
        "candidate_sha256": candidate.sha256,
        "time": scalar_time,
        "x": x,
        "z": z,
        "u": u,
        "v": v,
        "w": w,
        "speed": speed,
        "poloidal_speed": poloidal_speed,
    }


def _time_tag(time: float) -> str:
    return f"{float(time):.6f}".replace("-", "m").replace(".", "p")


def _write_slice_png(path: Path, data: dict[str, object]) -> None:
    # Import lazily so raw NPZ export still works in non-visual runtime installs.
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    x = np.asarray(data["x"], dtype=float)
    z = np.asarray(data["z"], dtype=float)
    speed = np.asarray(data["speed"], dtype=float)
    u = np.asarray(data["u"], dtype=float)
    w = np.asarray(data["w"], dtype=float)

    figure, axis = plt.subplots(figsize=(7.0, 7.0))
    mesh = axis.pcolormesh(x, z, speed, shading="auto")
    axis.streamplot(x, z, u, w, density=1.35, linewidth=0.65, arrowsize=0.65)
    figure.colorbar(mesh, ax=axis, label="|u|")
    axis.set_xlabel("x")
    axis.set_ylabel("z")
    axis.set_aspect("equal", adjustable="box")
    axis.set_xlim(float(x[0]), float(x[-1]))
    axis.set_ylim(float(z[0]), float(z[-1]))
    axis.set_title(f"Quartic supported Eq45 meridional slice, t={float(data['time']):.4f}")
    figure.tight_layout()
    figure.savefig(path, dpi=150)
    plt.close(figure)


def write_quartic_visualization_smoke(
    output_dir: str | Path,
    *,
    candidate: Eq45SupportedPhi10QuarticDerivativeBalancedTemporalCandidate | None = None,
    times: Iterable[float] = DEFAULT_TIMES,
    grid_size: int = DEFAULT_GRID_SIZE,
    extent: float = DEFAULT_EXTENT,
    write_png: bool = True,
) -> dict[str, object]:
    """Write deterministic NPZ slices and optional PNG render smoke plus a manifest."""
    field = governed_quartic_candidate() if candidate is None else candidate
    if not isinstance(field, Eq45SupportedPhi10QuarticDerivativeBalancedTemporalCandidate):
        raise TypeError("candidate must be the quartic supported Eq45 temporal candidate")
    size = _validated_grid_size(grid_size)
    radius = _validated_extent(extent)
    sample_times = _validated_times(field, times)

    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    if not destination.is_dir():
        raise ValueError("output_dir must resolve to a directory")

    entries: list[dict[str, object]] = []
    for time in sample_times:
        data = meridional_velocity_slice(
            field,
            time,
            grid_size=size,
            extent=radius,
        )
        tag = _time_tag(time)
        npz_name = f"quartic_eq45_meridional_t_{tag}.npz"
        npz_path = destination / npz_name
        np.savez_compressed(
            npz_path,
            candidate_sha256=np.asarray(field.sha256),
            time=np.asarray(time, dtype=float),
            x=np.asarray(data["x"], dtype=float),
            z=np.asarray(data["z"], dtype=float),
            u=np.asarray(data["u"], dtype=float),
            v=np.asarray(data["v"], dtype=float),
            w=np.asarray(data["w"], dtype=float),
            speed=np.asarray(data["speed"], dtype=float),
            poloidal_speed=np.asarray(data["poloidal_speed"], dtype=float),
        )
        entry: dict[str, object] = {
            "time": time,
            "phi10_coefficient": float(field.coefficient_at(time)),
            "npz": npz_name,
            "png": None,
            "max_speed": float(np.max(np.asarray(data["speed"], dtype=float))),
            "rms_speed": float(
                np.sqrt(np.mean(np.asarray(data["speed"], dtype=float) ** 2))
            ),
        }
        if write_png:
            png_name = f"quartic_eq45_meridional_t_{tag}.png"
            _write_slice_png(destination / png_name, data)
            entry["png"] = png_name
        entries.append(entry)

    manifest = {
        "schema": SCHEMA,
        "candidate_sha256": field.sha256,
        "base_supported_sha256": field.base_sha256,
        "candidate_class": type(field).__name__,
        "observable": "fixed_physical_y0_meridional_velocity_slice",
        "render": "speed_background_plus_projected_u_w_streamlines",
        "component_order": ["u", "v", "w"],
        "grid": {
            "grid_size": size,
            "extent": radius,
            "x_bounds": [-radius, radius],
            "y": 0.0,
            "z_bounds": [-radius, radius],
        },
        "times": list(sample_times),
        "files": entries,
        "truth_boundary": dict(_TRUTH_BOUNDARY),
    }
    manifest_path = destination / "quartic_eq45_visualization_smoke_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return manifest


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Render a fixed-frame smoke visualization of the quartic supported Eq45 field."
    )
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--candidate-json")
    parser.add_argument("--grid-size", type=int, default=DEFAULT_GRID_SIZE)
    parser.add_argument("--extent", type=float, default=DEFAULT_EXTENT)
    parser.add_argument("--times", type=float, nargs="+", default=list(DEFAULT_TIMES))
    parser.add_argument("--no-png", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    candidate = (
        Eq45SupportedPhi10QuarticDerivativeBalancedTemporalCandidate.load_json(
            args.candidate_json
        )
        if args.candidate_json
        else governed_quartic_candidate()
    )
    manifest = write_quartic_visualization_smoke(
        args.output_dir,
        candidate=candidate,
        times=args.times,
        grid_size=args.grid_size,
        extent=args.extent,
        write_png=not args.no_png,
    )
    print(json.dumps(manifest, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
