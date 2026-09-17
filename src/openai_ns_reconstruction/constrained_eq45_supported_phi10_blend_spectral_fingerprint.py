"""Resolution-aware spectrum fingerprint for the supported Phi10 blend family.

This module is visualization-side only.  It samples the public
``at_points(...)->[u,v,w]`` interface on a finite Cartesian box, checks Parseval,
measures radial Fourier-shell content and the opposite-face seam, and compares
the materialized compact--quartic blend at one fixed early off-keyframe.

The field is compactly supported by its parent representation, so the seam check
is retained as an explicit guard rather than silently treating a finite-window
FFT as a clean periodic/R^3 spectrum.  None of these quantities is a PDE gate or
blow-up/singularity test, and no blend weight is selected here.
"""
from __future__ import annotations

import json
from collections.abc import Iterable

import numpy as np

from .constrained_eq45_supported_candidate import Eq45SupportedVelocityCandidate
from .constrained_eq45_supported_collar_vorticity import governed_supported_seed
from .constrained_eq45_supported_phi10_compact_quartic_blend_temporal_mode import (
    Eq45SupportedPhi10CompactQuarticBlendTemporalCandidate,
)

SCHEMA = "eq45_supported_phi10_blend_spectral_fingerprint_v1"
DEFAULT_TIME = 0.3125
DEFAULT_BLEND_WEIGHTS = (0.0, 0.5, 1.0)
DEFAULT_GRID_SIZES = (24, 32, 40)
DEFAULT_BOX_HALF_WIDTH = 2.0
DEFAULT_FIXED_SHELL_CUTOFF = 6


def _validated_levels(grid_sizes: Iterable[int]) -> tuple[int, ...]:
    levels = tuple(int(value) for value in grid_sizes)
    if len(levels) < 3 or any(value < 16 for value in levels):
        raise ValueError("at least three grid sizes >=16 are required")
    if any(b <= a for a, b in zip(levels, levels[1:])):
        raise ValueError("grid sizes must be strictly increasing")
    return levels


def _validated_weights(blend_weights: Iterable[float]) -> tuple[float, ...]:
    weights = tuple(float(value) for value in blend_weights)
    if len(weights) < 3 or not np.all(np.isfinite(weights)):
        raise ValueError("at least three finite blend weights are required")
    if weights[0] != 0.0 or weights[-1] != 1.0:
        raise ValueError("weights must include quartic=0 and compact=1 endpoints")
    if any(not (0.0 <= value <= 1.0) for value in weights):
        raise ValueError("blend weights must lie in [0,1]")
    if any(b <= a for a, b in zip(weights, weights[1:])):
        raise ValueError("blend weights must be strictly increasing")
    return weights


def _velocity(field, points: np.ndarray, time: float) -> np.ndarray:
    if not hasattr(field, "at_points") or not callable(field.at_points):
        raise TypeError("field must expose at_points(points,time)")
    values = np.asarray(field.at_points(points, time), dtype=float)
    if values.shape != points.shape:
        raise ValueError("at_points must return the same (...,3) shape as points")
    if not np.all(np.isfinite(values)):
        raise ValueError("velocity output must be finite")
    return values


def _quantile_shell(cumulative: np.ndarray, fraction: float) -> int:
    return int(np.searchsorted(cumulative, fraction, side="left"))


def diagnose_blend_spectrum_level(
    field,
    *,
    time: float,
    grid_size: int,
    box_half_width: float = DEFAULT_BOX_HALF_WIDTH,
    fixed_shell_cutoff: int = DEFAULT_FIXED_SHELL_CUTOFF,
) -> dict[str, float | int]:
    """Measure one finite-box velocity spectrum through the public API only."""
    n = int(grid_size)
    half_width = float(box_half_width)
    cutoff = int(fixed_shell_cutoff)
    scalar_time = float(time)
    if n < 16:
        raise ValueError("grid_size must be >=16")
    if not np.isfinite(half_width) or half_width <= 0.0:
        raise ValueError("box_half_width must be positive and finite")
    if not np.isfinite(scalar_time):
        raise ValueError("time must be finite")
    if cutoff <= 0 or cutoff >= n // 2:
        raise ValueError("fixed_shell_cutoff must be resolved below Nyquist")

    # Endpoint-excluding samples are the natural DFT grid on a width-2L box.
    axis = -half_width + (2.0 * half_width) * np.arange(n, dtype=float) / n
    x, y, z = np.meshgrid(axis, axis, axis, indexing="ij")
    points = np.stack((x, y, z), axis=-1)
    velocity = _velocity(field, points, scalar_time)
    speed2 = np.sum(velocity * velocity, axis=-1)
    velocity_rms = float(np.sqrt(np.mean(speed2)))
    if not np.isfinite(velocity_rms) or velocity_rms <= 1e-14:
        raise ValueError("sampled velocity is numerically inactive")

    volume = (2.0 * half_width) ** 3
    direct_energy = float(0.5 * volume * np.mean(speed2))
    transformed = np.fft.fftn(velocity, axes=(0, 1, 2)) / n**3
    power = 0.5 * volume * np.sum(np.abs(transformed) ** 2, axis=-1)
    spectral_energy = float(np.sum(power))
    parseval_relative_error = abs(spectral_energy - direct_energy) / direct_energy

    modes = np.fft.fftfreq(n) * n
    kx, ky, kz = np.meshgrid(modes, modes, modes, indexing="ij")
    radius = np.sqrt(kx * kx + ky * ky + kz * kz)
    shell = np.floor(radius + 1e-12).astype(int)
    shell_energy = np.bincount(shell.ravel(), weights=power.ravel())
    cumulative = np.cumsum(shell_energy) / spectral_energy
    dominant_shell = int(1 + np.argmax(shell_energy[1:])) if len(shell_energy) > 1 else 0
    shell_rms = float(np.sqrt(np.sum((shell.astype(float) ** 2) * power) / spectral_energy))
    fixed_tail = float(np.sum(power[shell >= cutoff]) / spectral_energy)
    near_nyquist_mask = np.maximum.reduce((np.abs(kx), np.abs(ky), np.abs(kz))) >= 0.8 * (n / 2.0)
    near_nyquist = float(np.sum(power[near_nyquist_mask]) / spectral_energy)

    # Keep an explicit finite-window seam guard even though this child is
    # designed to vanish on the physical support boundary.
    tangential = np.linspace(-half_width, half_width, n, dtype=float)
    a, b = np.meshgrid(tangential, tangential, indexing="ij")
    jump_square_sum = 0.0
    seam_count = 0
    seam_max = 0.0
    for normal_axis in range(3):
        minus = np.zeros(a.shape + (3,), dtype=float)
        plus = np.zeros_like(minus)
        free_axes = [index for index in range(3) if index != normal_axis]
        minus[..., normal_axis] = -half_width
        plus[..., normal_axis] = half_width
        minus[..., free_axes[0]] = a
        plus[..., free_axes[0]] = a
        minus[..., free_axes[1]] = b
        plus[..., free_axes[1]] = b
        jump = _velocity(field, plus, scalar_time) - _velocity(field, minus, scalar_time)
        jump_norm = np.linalg.norm(jump, axis=-1)
        jump_square_sum += float(np.sum(jump_norm * jump_norm))
        seam_count += int(jump_norm.size)
        seam_max = max(seam_max, float(np.max(jump_norm)))
    seam_rms = float(np.sqrt(jump_square_sum / seam_count))

    fundamental = np.pi / half_width
    q50 = _quantile_shell(cumulative, 0.50)
    q90 = _quantile_shell(cumulative, 0.90)
    q99 = _quantile_shell(cumulative, 0.99)
    return {
        "grid_size": n,
        "direct_energy": direct_energy,
        "spectral_energy": spectral_energy,
        "parseval_relative_error": float(parseval_relative_error),
        "velocity_rms": velocity_rms,
        "dominant_shell": dominant_shell,
        "shell_q50": q50,
        "shell_q90": q90,
        "shell_q99": q99,
        "shell_rms": shell_rms,
        "physical_wavenumber_q90": float(q90 * fundamental),
        "physical_wavenumber_rms": float(shell_rms * fundamental),
        "fixed_shell_cutoff": cutoff,
        "fixed_shell_tail_fraction": fixed_tail,
        "near_nyquist_energy_fraction": near_nyquist,
        "opposite_face_seam_rms": seam_rms,
        "opposite_face_seam_max": seam_max,
        "seam_rms_over_velocity_rms": float(seam_rms / velocity_rms),
    }


def _relative_change(new: float, old: float) -> float:
    denominator = max(abs(float(old)), np.finfo(float).tiny)
    return float(abs(float(new) - float(old)) / denominator)


def audit_supported_phi10_blend_spectral_fingerprint(
    base: Eq45SupportedVelocityCandidate | None = None,
    *,
    time: float = DEFAULT_TIME,
    blend_weights: Iterable[float] = DEFAULT_BLEND_WEIGHTS,
    grid_sizes: Iterable[int] = DEFAULT_GRID_SIZES,
    box_half_width: float = DEFAULT_BOX_HALF_WIDTH,
    fixed_shell_cutoff: int = DEFAULT_FIXED_SHELL_CUTOFF,
) -> dict[str, object]:
    """Audit quartic/mid-blend/compact spectral content at three resolutions."""
    child = governed_supported_seed() if base is None else base
    if not isinstance(child, Eq45SupportedVelocityCandidate):
        raise TypeError("base must be Eq45SupportedVelocityCandidate")
    scalar_time = float(time)
    if not np.isfinite(scalar_time) or not (child.time_start <= scalar_time <= child.time_end):
        raise ValueError("time must lie inside the supported delivery interval")
    levels = _validated_levels(grid_sizes)
    weights = _validated_weights(blend_weights)
    if int(fixed_shell_cutoff) >= levels[0] // 2:
        raise ValueError("fixed_shell_cutoff must be resolved on the coarsest grid")

    rows: list[dict[str, object]] = []
    for weight in weights:
        field = Eq45SupportedPhi10CompactQuarticBlendTemporalCandidate(
            base=child, blend_weight=weight
        )
        for n in levels:
            row = diagnose_blend_spectrum_level(
                field,
                time=scalar_time,
                grid_size=n,
                box_half_width=box_half_width,
                fixed_shell_cutoff=fixed_shell_cutoff,
            )
            row["blend_weight"] = weight
            row["candidate_sha256"] = field.sha256
            rows.append(row)

    summaries: list[dict[str, object]] = []
    for weight in weights:
        weight_rows = [row for row in rows if row["blend_weight"] == weight]
        medium, fine = weight_rows[-2], weight_rows[-1]
        summaries.append(
            {
                "blend_weight": weight,
                "candidate_sha256": fine["candidate_sha256"],
                "finest_grid_size": fine["grid_size"],
                "finest_direct_energy": fine["direct_energy"],
                "finest_dominant_shell": fine["dominant_shell"],
                "finest_shell_q50": fine["shell_q50"],
                "finest_shell_q90": fine["shell_q90"],
                "finest_shell_q99": fine["shell_q99"],
                "finest_shell_rms": fine["shell_rms"],
                "finest_fixed_shell_tail_fraction": fine["fixed_shell_tail_fraction"],
                "finest_near_nyquist_energy_fraction": fine["near_nyquist_energy_fraction"],
                "finest_seam_rms_over_velocity_rms": fine["seam_rms_over_velocity_rms"],
                "medium_to_fine_energy_relative_change": _relative_change(fine["direct_energy"], medium["direct_energy"]),
                "medium_to_fine_shell_rms_relative_change": _relative_change(fine["shell_rms"], medium["shell_rms"]),
                "medium_to_fine_fixed_tail_relative_change": _relative_change(fine["fixed_shell_tail_fraction"], medium["fixed_shell_tail_fraction"]),
                "medium_to_fine_near_nyquist_relative_change": _relative_change(fine["near_nyquist_energy_fraction"], medium["near_nyquist_energy_fraction"]),
            }
        )

    quartic, midpoint, compact = summaries[0], summaries[len(summaries) // 2], summaries[-1]
    return {
        "schema": SCHEMA,
        "task_id": "CR008-EQ45-SUPPORTED-PHI10-BLEND-SPECTRAL-FINGERPRINT-023",
        "claim_scope": "supported_public_velocity_finite_box_spectrum_only",
        "base_supported_sha256": child.sha256,
        "time": scalar_time,
        "blend_weights": list(weights),
        "grid_sizes": list(levels),
        "box_half_width": float(box_half_width),
        "fixed_shell_cutoff": int(fixed_shell_cutoff),
        "fixed_shell_cutoff_classification": "autonomous_diagnostic_choice",
        "rows": rows,
        "summaries": summaries,
        "finest_endpoint_comparison": {
            "compact_over_quartic_shell_rms": float(compact["finest_shell_rms"] / quartic["finest_shell_rms"]),
            "compact_over_quartic_fixed_tail": float(compact["finest_fixed_shell_tail_fraction"] / quartic["finest_fixed_shell_tail_fraction"]),
            "compact_over_quartic_energy": float(compact["finest_direct_energy"] / quartic["finest_direct_energy"]),
            "midpoint_shell_rms_between_endpoints": bool(
                min(quartic["finest_shell_rms"], compact["finest_shell_rms"]) <= midpoint["finest_shell_rms"] <= max(quartic["finest_shell_rms"], compact["finest_shell_rms"])
            ),
            "midpoint_fixed_tail_between_endpoints": bool(
                min(quartic["finest_fixed_shell_tail_fraction"], compact["finest_fixed_shell_tail_fraction"]) <= midpoint["finest_fixed_shell_tail_fraction"] <= max(quartic["finest_fixed_shell_tail_fraction"], compact["finest_fixed_shell_tail_fraction"])
            ),
        },
        "interpretation": {
            "periodic_interpretation_requires_small_seam": True,
            "finite_window_spectral_shift_is_not_blowup_evidence": True,
            "blend_weight_selected": False,
        },
        "truth_boundary": {
            "canonical_velocity_changed": False,
            "physical_support_validated": False,
            "visualization_candidate_only": True,
            "visualization_ready": False,
            "visual_correspondence_verified": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
    }


def main() -> None:
    print(json.dumps(audit_supported_phi10_blend_spectral_fingerprint(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
