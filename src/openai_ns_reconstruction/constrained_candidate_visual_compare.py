"""Truth-bounded shared-frame comparison of multiple velocity candidates.

This module samples already-frozen public-style velocity callables on exactly the
same meridional image frame and times.  It supplies one global speed scale across
all candidate/time panels so per-panel autoscaling cannot make a weaker field
look equally strong.  Optional rendering is a projected ``(u, w)`` diagnostic,
not a true 3-D streamline calculation and not visual-correspondence evidence.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Sequence

import numpy as np
from numpy.typing import NDArray

VelocityCallable = Callable[[NDArray[np.float64], float], NDArray[np.float64]]


@dataclass(frozen=True)
class SharedFrameCandidateComparison:
    """Samples and fixed-frame normalization for multiple frozen candidates."""

    candidate_labels: tuple[str, ...]
    x: NDArray[np.float64]
    z: NDArray[np.float64]
    times: NDArray[np.float64]
    velocity: NDArray[np.float64]
    speed: NDArray[np.float64]
    poloidal_speed: NDArray[np.float64]
    swirl_speed: NDArray[np.float64]
    panel_max_speed: NDArray[np.float64]
    shared_speed_vmax: float
    y_plane: float
    frame_provenance: str
    candidate_provenance: tuple[str, ...]
    normalization_scope: str = "single_global_speed_scale_across_all_candidates_and_times"
    projected_streamline_semantics: str = "meridional_(u,w)_projection_only_not_true_3d_streamline"
    per_panel_autoscaling: bool = False
    camera_fitted: bool = False
    registration_fitted: bool = False
    velocity_changed: bool = False
    visualization_ready: bool = False
    visual_correspondence_verified: bool = False
    pde_validated: bool = False
    paper_exact: bool = False
    openai_field_identified: bool = False
    blowup_proved: bool = False

    @property
    def normalized_speed(self) -> NDArray[np.float64]:
        out = self.speed / self.shared_speed_vmax
        out.setflags(write=False)
        return out


def _strict_axis(values: Sequence[float], name: str, *, minimum_size: int) -> NDArray[np.float64]:
    array = np.asarray(values, dtype=np.float64)
    if array.ndim != 1 or array.size < minimum_size:
        raise ValueError(f"{name} must be a 1-D array with at least {minimum_size} values")
    if not bool(np.all(np.isfinite(array))):
        raise ValueError(f"{name} must contain only finite values")
    if not bool(np.all(np.diff(array) > 0.0)):
        raise ValueError(f"{name} must be strictly increasing")
    return array


def _provenance(value: str, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a non-empty provenance string")
    return value.strip()


def _readonly(array: NDArray[np.float64]) -> NDArray[np.float64]:
    out = np.asarray(array, dtype=np.float64)
    out.setflags(write=False)
    return out


def _candidate_entries(
    candidates: Sequence[tuple[str, VelocityCallable, str]],
) -> tuple[tuple[str, VelocityCallable, str], ...]:
    if not isinstance(candidates, Sequence) or isinstance(candidates, (str, bytes)):
        raise ValueError("candidates must be a sequence of (label, velocity, provenance) triples")
    entries = tuple(candidates)
    if len(entries) < 2:
        raise ValueError("at least two candidates are required for a comparison")

    labels: list[str] = []
    cleaned: list[tuple[str, VelocityCallable, str]] = []
    for index, entry in enumerate(entries):
        if not isinstance(entry, tuple) or len(entry) != 3:
            raise ValueError("each candidate must be a (label, velocity, provenance) triple")
        label, velocity, provenance = entry
        if not isinstance(label, str) or not label.strip():
            raise ValueError(f"candidate {index} label must be a non-empty string")
        label = label.strip()
        if label in labels:
            raise ValueError("candidate labels must be unique")
        if not callable(velocity):
            raise ValueError(f"candidate {label!r} velocity must be callable")
        source = _provenance(provenance, f"candidate {label!r} provenance")
        labels.append(label)
        cleaned.append((label, velocity, source))
    return tuple(cleaned)


def sample_shared_frame_candidates(
    candidates: Sequence[tuple[str, VelocityCallable, str]],
    *,
    x: Sequence[float],
    z: Sequence[float],
    times: Sequence[float],
    y_plane: float = 0.0,
    frame_provenance: str,
    min_panel_activity: float = 1e-12,
) -> SharedFrameCandidateComparison:
    """Sample frozen candidates on one identical meridional ``y=constant`` frame.

    The caller supplies every candidate, axis, time, the plane location, and
    provenance.  The function performs no registration, camera fitting,
    thresholding, hidden-time alignment, coefficient fitting, or candidate
    selection.  Every candidate/time panel is sampled on the exact same points.

    ``min_panel_activity`` is only a fail-closed guard against a zero or
    numerically inactive panel; it is not a visual acceptance threshold.
    """

    entries = _candidate_entries(candidates)
    x_axis = _strict_axis(x, "x", minimum_size=3)
    z_axis = _strict_axis(z, "z", minimum_size=3)
    time_axis = _strict_axis(times, "times", minimum_size=1)
    frame_source = _provenance(frame_provenance, "frame_provenance")

    y_value = float(y_plane)
    activity_floor = float(min_panel_activity)
    if not np.isfinite(y_value):
        raise ValueError("y_plane must be finite")
    if not np.isfinite(activity_floor) or activity_floor <= 0.0:
        raise ValueError("min_panel_activity must be finite and strictly positive")

    xx, zz = np.meshgrid(x_axis, z_axis, indexing="xy")
    points = np.column_stack(
        [xx.ravel(), np.full(xx.size, y_value, dtype=np.float64), zz.ravel()]
    )

    n_candidates = len(entries)
    n_times = time_axis.size
    velocity_grid = np.empty(
        (n_candidates, n_times, z_axis.size, x_axis.size, 3), dtype=np.float64
    )

    panel_max = np.empty((n_candidates, n_times), dtype=np.float64)
    for candidate_index, (label, velocity_fn, _) in enumerate(entries):
        for time_index, time_value in enumerate(time_axis):
            values = np.asarray(velocity_fn(points.copy(), float(time_value)), dtype=np.float64)
            if values.shape != (points.shape[0], 3):
                raise ValueError(
                    f"candidate {label!r} returned shape {values.shape}; expected {(points.shape[0], 3)}"
                )
            if not bool(np.all(np.isfinite(values))):
                raise ValueError(f"candidate {label!r} returned nonfinite velocity values")
            reshaped = values.reshape(z_axis.size, x_axis.size, 3)
            speed = np.linalg.norm(reshaped, axis=-1)
            maximum = float(np.max(speed))
            if maximum <= activity_floor:
                raise ValueError(
                    f"candidate {label!r} is inactive on the declared frame at t={float(time_value):.17g}"
                )
            velocity_grid[candidate_index, time_index] = reshaped
            panel_max[candidate_index, time_index] = maximum

    speed_grid = np.linalg.norm(velocity_grid, axis=-1)
    poloidal_grid = np.hypot(velocity_grid[..., 0], velocity_grid[..., 2])
    swirl_grid = np.abs(velocity_grid[..., 1])
    shared_vmax = float(np.max(speed_grid))
    if not np.isfinite(shared_vmax) or shared_vmax <= activity_floor:
        raise ValueError("shared speed scale is not finite and active")

    return SharedFrameCandidateComparison(
        candidate_labels=tuple(entry[0] for entry in entries),
        x=_readonly(x_axis.copy()),
        z=_readonly(z_axis.copy()),
        times=_readonly(time_axis.copy()),
        velocity=_readonly(velocity_grid),
        speed=_readonly(speed_grid),
        poloidal_speed=_readonly(poloidal_grid),
        swirl_speed=_readonly(swirl_grid),
        panel_max_speed=_readonly(panel_max),
        shared_speed_vmax=shared_vmax,
        y_plane=y_value,
        frame_provenance=frame_source,
        candidate_provenance=tuple(entry[2] for entry in entries),
    )


def render_shared_frame_candidate_comparison(
    comparison: SharedFrameCandidateComparison,
    output_path: str | Path,
    *,
    streamline_density: float = 1.0,
    linewidth: float = 0.55,
    cmap: str = "viridis",
    dpi: int = 140,
) -> Path:
    """Render a side-by-side speed comparison with one shared color normalization.

    The overlay uses only the projected meridional components ``(u, w)``.  These
    lines must not be described as true 3-D streamlines when the swirl component
    is nonzero.  Rendering performs no candidate alignment or parameter fitting.
    """

    if not isinstance(comparison, SharedFrameCandidateComparison):
        raise ValueError("comparison must be a SharedFrameCandidateComparison")
    density = float(streamline_density)
    line_width = float(linewidth)
    if not np.isfinite(density) or density <= 0.0:
        raise ValueError("streamline_density must be finite and positive")
    if not np.isfinite(line_width) or line_width <= 0.0:
        raise ValueError("linewidth must be finite and positive")
    if isinstance(dpi, bool) or not isinstance(dpi, (int, np.integer)) or int(dpi) < 50:
        raise ValueError("dpi must be an integer >= 50")
    if not isinstance(cmap, str) or not cmap.strip():
        raise ValueError("cmap must be a non-empty Matplotlib colormap name")

    path = Path(output_path)
    if path.suffix.lower() != ".png":
        raise ValueError("output_path must end in .png")
    path.parent.mkdir(parents=True, exist_ok=True)

    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize

    norm = Normalize(vmin=0.0, vmax=comparison.shared_speed_vmax, clip=True)
    rows = len(comparison.candidate_labels)
    cols = comparison.times.size
    figure, axes = plt.subplots(
        rows,
        cols,
        figsize=(4.0 * cols, 3.4 * rows),
        squeeze=False,
        sharex=True,
        sharey=True,
        constrained_layout=True,
    )

    image = None
    for candidate_index, label in enumerate(comparison.candidate_labels):
        for time_index, time_value in enumerate(comparison.times):
            axis = axes[candidate_index, time_index]
            speed = comparison.speed[candidate_index, time_index]
            velocity = comparison.velocity[candidate_index, time_index]
            image = axis.pcolormesh(
                comparison.x,
                comparison.z,
                speed,
                shading="auto",
                cmap=cmap,
                norm=norm,
            )
            axis.streamplot(
                comparison.x,
                comparison.z,
                velocity[..., 0],
                velocity[..., 2],
                density=density,
                linewidth=line_width,
                color="white",
                arrowsize=0.6,
            )
            axis.set_title(f"{label}  t={float(time_value):.4g}")
            axis.set_aspect("equal", adjustable="box")
            axis.set_xlabel("x")
            axis.set_ylabel("z")

    if image is None:
        raise ValueError("comparison contained no renderable panels")
    figure.colorbar(image, ax=axes.ravel().tolist(), label="|u| (shared scale)")
    figure.suptitle(
        "Fixed meridional frame; white lines are projected (u,w), not true 3-D streamlines"
    )
    figure.savefig(path, dpi=int(dpi))
    plt.close(figure)
    if not path.exists() or path.stat().st_size <= 0:
        raise RuntimeError("Matplotlib did not create a nonempty PNG")
    return path
