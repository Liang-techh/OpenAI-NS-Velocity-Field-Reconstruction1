"""Resolution-aware finite-window spectrum audit for a frozen public velocity field.

This diagnostic is intentionally visualization-side only.  It samples a field
through ``at_points(points, time)`` on a finite Cartesian box, checks Parseval,
reports radial Fourier-shell energy scales, and measures opposite-face seam
mismatch.  The seam metric is important because an FFT on a non-periodic finite
window can manufacture high-frequency content even when the interior field is
smooth.

No spectral quantity in this module is a Navier--Stokes acceptance criterion or
blow-up test.  In particular, the current Eq. (4.5) seed has not yet received the
physical-space support connection, so its periodic-box spectrum must retain the
seam warning rather than being interpreted as an R^3 spectrum.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable

import numpy as np

from .constrained_eq45_candidate import Eq45VelocityCandidate


SCHEMA = "eq45_public_spectral_audit_v1"
DEFAULT_BOX_HALF_WIDTH = 2.0
DEFAULT_GRID_SIZES = (16, 24, 32)
DEFAULT_TIMES = (0.25, 0.5, 0.75)
DEFAULT_FIXED_SHELL_CUTOFF = 6


def _velocity(field, points, time):
    if not hasattr(field, "at_points") or not callable(field.at_points):
        raise TypeError("field must expose callable at_points(points, time)")
    values = np.asarray(field.at_points(points, time), dtype=float)
    if values.shape != np.asarray(points).shape:
        raise ValueError("at_points must return the same (...,3) shape as points")
    if not np.all(np.isfinite(values)):
        raise ValueError("velocity output must be finite")
    return values


def _validated_levels(grid_sizes: Iterable[int]) -> tuple[int, ...]:
    levels = tuple(int(n) for n in grid_sizes)
    if len(levels) < 3 or any(n < 8 for n in levels):
        raise ValueError("at least three grid sizes >= 8 are required")
    if any(b <= a for a, b in zip(levels, levels[1:])):
        raise ValueError("grid sizes must be strictly increasing")
    return levels


def _validated_times(times: Iterable[float]) -> tuple[float, ...]:
    out = tuple(float(t) for t in times)
    if not out or not np.all(np.isfinite(out)):
        raise ValueError("times must be a nonempty finite sequence")
    if any(b <= a for a, b in zip(out, out[1:])):
        raise ValueError("times must be strictly increasing")
    return out


def _quantile_shell(cumulative: np.ndarray, fraction: float) -> int:
    return int(np.searchsorted(cumulative, fraction, side="left"))


def diagnose_spectrum_level(
    field,
    *,
    time: float,
    grid_size: int,
    box_half_width: float = DEFAULT_BOX_HALF_WIDTH,
    fixed_shell_cutoff: int = DEFAULT_FIXED_SHELL_CUTOFF,
) -> dict:
    """Audit one finite-window FFT level using only the public velocity path."""
    n = int(grid_size)
    L = float(box_half_width)
    cutoff = int(fixed_shell_cutoff)
    if n < 8:
        raise ValueError("grid_size must be >= 8")
    if not np.isfinite(L) or L <= 0.0:
        raise ValueError("box_half_width must be positive and finite")
    if cutoff <= 0 or cutoff >= n // 2:
        raise ValueError("fixed_shell_cutoff must lie between 1 and grid_size/2")
    if not np.isfinite(time):
        raise ValueError("time must be finite")

    axis = -L + (2.0 * L) * np.arange(n, dtype=float) / n
    x, y, z = np.meshgrid(axis, axis, axis, indexing="ij")
    points = np.stack((x, y, z), axis=-1)
    velocity = _velocity(field, points, float(time))

    volume = (2.0 * L) ** 3
    velocity_rms = float(np.sqrt(np.mean(np.sum(velocity * velocity, axis=-1))))
    if not np.isfinite(velocity_rms) or velocity_rms <= 1e-14:
        raise ValueError("sampled velocity is numerically inactive")
    direct_energy = float(0.5 * volume * np.mean(np.sum(velocity * velocity, axis=-1)))

    transformed = np.fft.fftn(velocity, axes=(0, 1, 2)) / n**3
    spectral_power = 0.5 * volume * np.sum(np.abs(transformed) ** 2, axis=-1)
    spectral_energy = float(np.sum(spectral_power))
    parseval_relative_error = abs(spectral_energy - direct_energy) / direct_energy

    modes = np.fft.fftfreq(n) * n
    kx, ky, kz = np.meshgrid(modes, modes, modes, indexing="ij")
    shell_index = np.floor(np.sqrt(kx * kx + ky * ky + kz * kz) + 1e-12).astype(int)
    shell_energy = np.bincount(shell_index.ravel(), weights=spectral_power.ravel())
    cumulative = np.cumsum(shell_energy) / spectral_energy
    dominant_shell = int(1 + np.argmax(shell_energy[1:])) if len(shell_energy) > 1 else 0

    fixed_tail = float(np.sum(spectral_power[shell_index >= cutoff]) / spectral_energy)
    near_nyquist_mask = np.maximum.reduce((np.abs(kx), np.abs(ky), np.abs(kz))) >= 0.8 * (n / 2.0)
    near_nyquist = float(np.sum(spectral_power[near_nyquist_mask]) / spectral_energy)

    squared_jump_sum = 0.0
    seam_count = 0
    seam_max = 0.0
    a, b = np.meshgrid(axis, axis, indexing="ij")
    for normal_axis in range(3):
        minus = np.zeros(a.shape + (3,), dtype=float)
        plus = np.zeros_like(minus)
        tangential = [index for index in range(3) if index != normal_axis]
        minus[..., normal_axis] = -L
        plus[..., normal_axis] = L
        minus[..., tangential[0]] = a
        plus[..., tangential[0]] = a
        minus[..., tangential[1]] = b
        plus[..., tangential[1]] = b
        jump = _velocity(field, plus, float(time)) - _velocity(field, minus, float(time))
        jump_norm = np.linalg.norm(jump, axis=-1)
        squared_jump_sum += float(np.sum(jump_norm * jump_norm))
        seam_count += int(jump_norm.size)
        seam_max = max(seam_max, float(np.max(jump_norm)))

    seam_rms = float(np.sqrt(squared_jump_sum / seam_count))
    fundamental_wavenumber = np.pi / L
    return {
        "grid_size": n,
        "time": float(time),
        "box_half_width": L,
        "direct_energy": direct_energy,
        "spectral_energy": spectral_energy,
        "parseval_relative_error": float(parseval_relative_error),
        "velocity_rms": velocity_rms,
        "dominant_shell": dominant_shell,
        "shell_q50": _quantile_shell(cumulative, 0.50),
        "shell_q90": _quantile_shell(cumulative, 0.90),
        "shell_q99": _quantile_shell(cumulative, 0.99),
        "physical_wavenumber_q50": float(_quantile_shell(cumulative, 0.50) * fundamental_wavenumber),
        "physical_wavenumber_q90": float(_quantile_shell(cumulative, 0.90) * fundamental_wavenumber),
        "fixed_shell_cutoff": cutoff,
        "fixed_shell_tail_fraction": fixed_tail,
        "near_nyquist_energy_fraction": near_nyquist,
        "opposite_face_seam_rms": seam_rms,
        "opposite_face_seam_max": seam_max,
        "seam_rms_over_velocity_rms": float(seam_rms / velocity_rms),
    }


def audit_spectral_resolution(
    field,
    *,
    times: Iterable[float] = DEFAULT_TIMES,
    grid_sizes: Iterable[int] = DEFAULT_GRID_SIZES,
    box_half_width: float = DEFAULT_BOX_HALF_WIDTH,
    fixed_shell_cutoff: int = DEFAULT_FIXED_SHELL_CUTOFF,
) -> dict:
    """Run three-or-more FFT resolutions at every declared visualization time."""
    levels = _validated_levels(grid_sizes)
    sample_times = _validated_times(times)
    if int(fixed_shell_cutoff) >= levels[0] // 2:
        raise ValueError("fixed_shell_cutoff must be resolved on the coarsest grid")

    rows = [
        diagnose_spectrum_level(
            field,
            time=time,
            grid_size=n,
            box_half_width=box_half_width,
            fixed_shell_cutoff=fixed_shell_cutoff,
        )
        for n in levels
        for time in sample_times
    ]

    summaries = []
    for time in sample_times:
        time_rows = [row for row in rows if row["time"] == time]
        finest = time_rows[-1]
        summaries.append(
            {
                "time": time,
                "finest_grid_size": finest["grid_size"],
                "finest_energy": finest["direct_energy"],
                "finest_dominant_shell": finest["dominant_shell"],
                "finest_shell_q50": finest["shell_q50"],
                "finest_shell_q90": finest["shell_q90"],
                "finest_shell_q99": finest["shell_q99"],
                "finest_fixed_shell_tail_fraction": finest["fixed_shell_tail_fraction"],
                "finest_near_nyquist_energy_fraction": finest["near_nyquist_energy_fraction"],
                "finest_seam_rms_over_velocity_rms": finest["seam_rms_over_velocity_rms"],
                "max_energy_relative_delta_to_finest": float(
                    max(abs(row["direct_energy"] - finest["direct_energy"]) / finest["direct_energy"] for row in time_rows[:-1])
                ),
                "max_shell_q90_delta_to_finest": int(
                    max(abs(row["shell_q90"] - finest["shell_q90"]) for row in time_rows[:-1])
                ),
            }
        )

    return {
        "schema": SCHEMA,
        "claim_scope": "finite_window_visualization_spectrum_only",
        "box_half_width": float(box_half_width),
        "grid_sizes": list(levels),
        "times": list(sample_times),
        "fixed_shell_cutoff": int(fixed_shell_cutoff),
        "fixed_shell_cutoff_classification": "autonomous_diagnostic_choice",
        "periodic_interpretation_requires_small_seam": True,
        "rows": rows,
        "time_summaries": summaries,
        "truth_boundary": {
            "physical_support_validated": False,
            "visualization_ready": False,
            "visual_correspondence_verified": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
    }


def audit_checked_eq45_seed(
    candidate_path: str | Path = "artifacts/constrained/eq45_velocity_candidate_seed.json",
) -> dict:
    """Load the frozen Eq. (4.5) seed artifact and audit its public velocity path."""
    field = Eq45VelocityCandidate.load_json(candidate_path)
    report = audit_spectral_resolution(field)
    report["candidate_path"] = str(candidate_path)
    report["candidate_sha256"] = field.sha256
    return report


def main() -> None:
    print(json.dumps(audit_checked_eq45_seed(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
