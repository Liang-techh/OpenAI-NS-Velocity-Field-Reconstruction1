"""Truth-bounded Python slice visualization smoke for callable velocity fields.

This module is intentionally candidate-agnostic. It is designed to consume the
final ``velocity(points, time) -> [..., 3]`` interface once the Eq. (4.5)
profile/evaluator stack is composed. It does not infer hidden profile
parameters, fit a field, evaluate a Navier--Stokes residual, or decide visual
correspondence.

Two orthogonal slices are sampled directly from the callable at every requested
time:

* meridional plane ``y = 0`` with vector components ``(u, w)``;
* equatorial plane ``z = 0`` with vector components ``(u, v)``.

The renderer uses one speed scale across all panels and times so apparent growth
is not created by independent per-panel color normalization. A numerically
inactive field fails closed, preventing a zero field from producing a
"successful" visualization smoke artifact.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Sequence

import numpy as np


VelocityCallable = Callable[[np.ndarray, float], np.ndarray]


@dataclass(frozen=True)
class VelocitySliceSequence:
    """Direct velocity samples used by the visualization smoke."""

    times: np.ndarray
    horizontal: np.ndarray
    axial: np.ndarray
    meridional_velocity: np.ndarray
    equatorial_velocity: np.ndarray
    speed_max_by_time: np.ndarray
    global_speed_max: float

    @property
    def metadata(self) -> dict[str, object]:
        return {
            "claim_scope": "visualization_smoke_only",
            "slice_planes": ["y=0", "z=0"],
            "times": self.times.tolist(),
            "grid_size": int(self.horizontal.size),
            "global_speed_max": float(self.global_speed_max),
            "speed_max_by_time": self.speed_max_by_time.tolist(),
            "visual_correspondence_verified": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        }


def _evaluate_velocity(velocity: VelocityCallable, points: np.ndarray, time: float) -> np.ndarray:
    values = np.asarray(velocity(points, float(time)), dtype=float)
    if values.shape != points.shape:
        raise ValueError(
            "velocity(points,time) must return an array with the same (...,3) shape as points"
        )
    if not np.all(np.isfinite(values)):
        raise ValueError("velocity output must be finite")
    return values


def sample_velocity_slice_sequence(
    velocity: VelocityCallable,
    *,
    times: Sequence[float],
    horizontal_extent: float = 2.0,
    axial_extent: float = 2.0,
    grid_size: int = 65,
    min_activity: float = 1e-12,
) -> VelocitySliceSequence:
    """Sample orthogonal slices directly from a public-style velocity callable.

    ``horizontal_extent`` is used for both x and y on the equatorial plane.
    ``axial_extent`` is used for z on the meridional plane. ``times`` must be
    strictly increasing so a rendered sequence has an unambiguous temporal
    ordering.
    """

    if not callable(velocity):
        raise ValueError("velocity must be callable")
    time_array = np.asarray(tuple(times), dtype=float)
    if time_array.ndim != 1 or time_array.size < 2:
        raise ValueError("times must contain at least two values")
    if not np.all(np.isfinite(time_array)) or np.any(np.diff(time_array) <= 0.0):
        raise ValueError("times must be finite and strictly increasing")
    for name, value in (
        ("horizontal_extent", horizontal_extent),
        ("axial_extent", axial_extent),
        ("min_activity", min_activity),
    ):
        if not np.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be positive and finite")
    if not isinstance(grid_size, int) or grid_size < 9 or grid_size % 2 == 0:
        raise ValueError("grid_size must be an odd integer >= 9")

    horizontal = np.linspace(-horizontal_extent, horizontal_extent, grid_size)
    axial = np.linspace(-axial_extent, axial_extent, grid_size)

    mx, mz = np.meshgrid(horizontal, axial, indexing="xy")
    meridional_points = np.stack(
        (mx, np.zeros_like(mx), mz), axis=-1
    ).reshape(-1, 3)

    ex, ey = np.meshgrid(horizontal, horizontal, indexing="xy")
    equatorial_points = np.stack(
        (ex, ey, np.zeros_like(ex)), axis=-1
    ).reshape(-1, 3)

    meridional = []
    equatorial = []
    speed_max_by_time = []
    for time in time_array:
        mer = _evaluate_velocity(velocity, meridional_points, float(time)).reshape(
            grid_size, grid_size, 3
        )
        eq = _evaluate_velocity(velocity, equatorial_points, float(time)).reshape(
            grid_size, grid_size, 3
        )
        meridional.append(mer)
        equatorial.append(eq)
        speed_max_by_time.append(
            max(
                float(np.max(np.linalg.norm(mer, axis=-1))),
                float(np.max(np.linalg.norm(eq, axis=-1))),
            )
        )

    meridional_velocity = np.stack(meridional, axis=0)
    equatorial_velocity = np.stack(equatorial, axis=0)
    speed_max = np.asarray(speed_max_by_time, dtype=float)
    global_speed_max = float(np.max(speed_max))
    if not np.isfinite(global_speed_max) or global_speed_max <= min_activity:
        raise ValueError(
            "sampled velocity is numerically inactive; visualization smoke refuses a zero-field path"
        )

    for array in (
        time_array,
        horizontal,
        axial,
        meridional_velocity,
        equatorial_velocity,
        speed_max,
    ):
        array.setflags(write=False)

    return VelocitySliceSequence(
        times=time_array,
        horizontal=horizontal,
        axial=axial,
        meridional_velocity=meridional_velocity,
        equatorial_velocity=equatorial_velocity,
        speed_max_by_time=speed_max,
        global_speed_max=global_speed_max,
    )


def render_velocity_slice_sequence(
    velocity: VelocityCallable,
    output_path: str | Path,
    *,
    times: Sequence[float],
    horizontal_extent: float = 2.0,
    axial_extent: float = 2.0,
    grid_size: int = 65,
    quiver_stride: int = 5,
    min_activity: float = 1e-12,
    dpi: int = 140,
) -> dict[str, object]:
    """Render a deterministic multi-time orthogonal-slice PNG.

    The top row shows the meridional ``y=0`` plane; the bottom row shows the
    equatorial ``z=0`` plane. All panels use the same speed normalization.
    The returned metadata is visualization-only evidence and keeps all stronger
    claims explicitly false.
    """

    if not isinstance(quiver_stride, int) or quiver_stride <= 0:
        raise ValueError("quiver_stride must be a positive integer")
    if not isinstance(dpi, int) or dpi <= 0:
        raise ValueError("dpi must be a positive integer")

    sequence = sample_velocity_slice_sequence(
        velocity,
        times=times,
        horizontal_extent=horizontal_extent,
        axial_extent=axial_extent,
        grid_size=grid_size,
        min_activity=min_activity,
    )

    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    n_times = sequence.times.size
    figure = Figure(figsize=(4.0 * n_times, 7.0), constrained_layout=True)
    FigureCanvasAgg(figure)
    axes = figure.subplots(2, n_times, squeeze=False)

    horizontal = sequence.horizontal
    axial = sequence.axial
    mx, mz = np.meshgrid(horizontal, axial, indexing="xy")
    ex, ey = np.meshgrid(horizontal, horizontal, indexing="xy")
    stride = slice(None, None, quiver_stride)

    image = None
    for column, time in enumerate(sequence.times):
        mer = sequence.meridional_velocity[column]
        eq = sequence.equatorial_velocity[column]
        mer_speed = np.linalg.norm(mer, axis=-1)
        eq_speed = np.linalg.norm(eq, axis=-1)

        ax = axes[0, column]
        image = ax.pcolormesh(
            horizontal,
            axial,
            mer_speed,
            shading="auto",
            vmin=0.0,
            vmax=sequence.global_speed_max,
        )
        ax.quiver(
            mx[stride, stride],
            mz[stride, stride],
            mer[stride, stride, 0],
            mer[stride, stride, 2],
        )
        ax.set_title(f"y=0, t={time:.6g}")
        ax.set_xlabel("x")
        ax.set_ylabel("z")
        ax.set_aspect("equal")

        ax = axes[1, column]
        ax.pcolormesh(
            horizontal,
            horizontal,
            eq_speed,
            shading="auto",
            vmin=0.0,
            vmax=sequence.global_speed_max,
        )
        ax.quiver(
            ex[stride, stride],
            ey[stride, stride],
            eq[stride, stride, 0],
            eq[stride, stride, 1],
        )
        ax.set_title(f"z=0, t={time:.6g}")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_aspect("equal")

    if image is not None:
        figure.colorbar(image, ax=axes.ravel().tolist(), label="|velocity|")

    target = Path(output_path)
    if target.suffix.lower() != ".png":
        raise ValueError("output_path must end in .png")
    target.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(target, dpi=dpi, format="png")
    figure.clear()

    if not target.exists() or target.stat().st_size == 0:
        raise RuntimeError("visualization PNG was not written")

    metadata = dict(sequence.metadata)
    metadata.update(
        {
            "output_path": str(target),
            "renderer": "matplotlib_agg",
            "shared_speed_normalization": True,
            "quiver_stride": quiver_stride,
            "visualization_smoke_rendered": True,
        }
    )
    return metadata
