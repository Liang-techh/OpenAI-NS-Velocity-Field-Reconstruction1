"""FD6 spatial Laplacian seam for the frozen A2 signed-amplitude correction.

This Agent-2 module consumes an already-materialized
``SignedRadialAmplitudeCorrection`` and exposes Cartesian pure second
derivatives plus ``Delta(delta u)``.  The correction remains vector-potential
first and divergence-free by the upstream complete-curl construction.  No
Agent-3 defect/mean/radial inverse, pressure, forcing, target, gain, residual,
or scientific threshold enters the public evaluator.

The corrected Kokuno reconstruction supplies structural provenance for the
localized complete-curl construction.  The compact spline lift, public-z
pullback, autonomous time modulation, and finite-difference operator used here
are repository realizations rather than paper-exact hidden data.
"""
from __future__ import annotations

import argparse
import inspect
import json
import math
from pathlib import Path
from typing import Any, Callable

import numpy as np

from .kokuno_public_signed_amplitude_complete_curl import (
    SignedRadialAmplitudeCorrection,
    correction_velocity,
    profile_from_delta_a,
)

TASK = "KOKUNO-A2-SIGNED-AMPLITUDE-SPATIAL-LAPLACIAN-052"
SCHEMA = "kokuno-a2-signed-amplitude-spatial-laplacian-v1"
PARENT_AGENT2_PR = 742
PARENT_AGENT2_HEAD = "e4b5be5f768df6e5cdd09773487ec9888e3891b7"
SOURCE_CORRECTED_READER_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_CORRECTED_READER_DATE = "2026-09-09"
FD6_STEPS = (0.004, 0.002, 0.001)

_FD6_D2_OFFSETS = (-3, -2, -1, 0, 1, 2, 3)
_FD6_D2_COEFFICIENTS = (2.0, -27.0, 270.0, -490.0, 270.0, -27.0, 2.0)
_FD6_D2_DENOMINATOR = 180.0

VectorProvider = Callable[[Any, Any, Any, Any], np.ndarray]


def _broadcast_xyzt(x: Any, y: Any, z: Any, t: Any) -> tuple[np.ndarray, ...]:
    arrays = np.broadcast_arrays(
        np.asarray(x, dtype=float),
        np.asarray(y, dtype=float),
        np.asarray(z, dtype=float),
        np.asarray(t, dtype=float),
    )
    if not all(np.all(np.isfinite(array)) for array in arrays):
        raise ValueError("x, y, z and t must be finite and broadcastable")
    return tuple(arrays)


def _validate_step(step: float) -> float:
    h = float(step)
    if not math.isfinite(h) or not (1.0e-5 <= h <= 5.0e-2):
        raise ValueError("spatial_step must be finite and lie in [1e-5, 5e-2]")
    return h


def _provider_value(
    provider: VectorProvider,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
) -> np.ndarray:
    value = np.asarray(provider(x, y, z, t), dtype=float)
    expected = x.shape + (3,)
    if value.shape != expected:
        raise RuntimeError(
            f"vector provider returned shape {value.shape}, expected {expected}"
        )
    if not np.all(np.isfinite(value)):
        raise RuntimeError("vector provider returned non-finite values")
    return value


def _fd6_second_derivatives(
    provider: VectorProvider,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    *,
    spatial_step: float,
) -> np.ndarray:
    """Return pure second derivatives as ``[..., component, axis]``."""
    h = _validate_step(spatial_step)
    xb, yb, zb, tb = _broadcast_xyzt(x, y, z, t)
    base = (xb, yb, zb)
    derivatives: list[np.ndarray] = []
    for axis in range(3):
        accum = np.zeros(xb.shape + (3,), dtype=float)
        for offset, coefficient in zip(
            _FD6_D2_OFFSETS, _FD6_D2_COEFFICIENTS, strict=True
        ):
            shifted = [array for array in base]
            shifted[axis] = shifted[axis] + offset * h
            accum += coefficient * _provider_value(
                provider, shifted[0], shifted[1], shifted[2], tb
            )
        derivatives.append(accum / (_FD6_D2_DENOMINATOR * h * h))
    return np.stack(derivatives, axis=-1)


def evaluate_correction_velocity_laplacian_fd6(
    correction: SignedRadialAmplitudeCorrection,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    *,
    spatial_step: float = 0.002,
) -> dict[str, np.ndarray | float | str]:
    """Evaluate ``delta u``, its pure second derivatives, and its Laplacian."""
    h = _validate_step(spatial_step)
    xb, yb, zb, tb = _broadcast_xyzt(x, y, z, t)

    def provider(xx: Any, yy: Any, zz: Any, tt: Any) -> np.ndarray:
        return correction_velocity(correction, xx, yy, zz, tt)

    velocity = _provider_value(provider, xb, yb, zb, tb)
    second = _fd6_second_derivatives(
        provider, xb, yb, zb, tb, spatial_step=h
    )
    laplacian = np.sum(second, axis=-1)
    return {
        "velocity": velocity,
        "second_derivatives": second,
        "laplacian": laplacian,
        "spatial_step": h,
        "spatial_operator": "centered_cartesian_fd6_second_derivative",
    }


def _vector_rms(values: np.ndarray) -> float:
    arr = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(np.sum(arr * arr, axis=-1))))


def _verification_profile() -> SignedRadialAmplitudeCorrection:
    radii = np.linspace(0.31, 1.19, 21)
    s = (radii - radii[0]) / (radii[-1] - radii[0])
    bump = np.sin(math.pi * s) ** 8
    delta_a = np.stack(
        (
            0.0095 * bump * (1.0 + 0.09 * np.cos(2.0 * math.pi * s)),
            -0.0068 * bump * (1.0 - 0.07 * np.sin(2.0 * math.pi * s)),
        ),
        axis=-1,
    )
    delta_a[[0, -1], :] = 0.0
    return profile_from_delta_a(
        radii,
        delta_a,
        reference_time=0.50,
        producer_kind="analytic-regression-not-agent3-real-candidate",
        provenance=(
            "fresh compact signed radial mechanics profile for correction-Laplacian "
            "stability only; no Agent-3 defect, stress, inverse, or residual input"
        ),
        source_mean_amplitude_differential_certified=False,
    )


def _verification_cloud() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    # Mid-cell radii keep the finest three-point reach away from spline knots.
    r0 = 0.31
    dr = (1.19 - 0.31) / 20.0
    cell_ids = (2, 5, 8, 11, 14, 17)
    radii = np.asarray([r0 + (index + 0.5) * dr for index in cell_ids], dtype=float)
    z_block = np.asarray((-0.58, -0.36, -0.12, 0.14, 0.37, 0.59), dtype=float)
    golden = math.pi * (3.0 - math.sqrt(5.0))
    rows: list[tuple[float, float, float, float]] = []
    count = 0
    for time in (0.39, 0.50, 0.61):
        for radius, z_value in zip(radii, z_block, strict=True):
            angle = 0.193 + (count + 1) * golden
            rows.append(
                (
                    radius * math.cos(angle),
                    radius * math.sin(angle),
                    float(z_value),
                    time,
                )
            )
            count += 1
    array = np.asarray(rows, dtype=float)
    return array[:, 0], array[:, 1], array[:, 2], array[:, 3]


def _manufactured_provider(x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
    xb, yb, zb, tb = _broadcast_xyzt(x, y, z, t)
    del tb
    return np.stack(
        (
            xb**6 + 2.0 * yb**4 + 3.0 * zb**2,
            0.5 * xb**4 - yb**6 + zb**4,
            2.0 * xb**2 + 3.0 * yb**2 + 4.0 * zb**6,
        ),
        axis=-1,
    )


def _manufactured_laplacian(x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
    xb, yb, zb, tb = _broadcast_xyzt(x, y, z, t)
    del tb
    return np.stack(
        (
            30.0 * xb**4 + 24.0 * yb**2 + 6.0,
            6.0 * xb**2 - 30.0 * yb**4 + 12.0 * zb**2,
            10.0 + 120.0 * zb**4,
        ),
        axis=-1,
    )


def verification_receipt() -> dict[str, Any]:
    correction = _verification_profile()
    x, y, z, t = _verification_cloud()
    levels = [
        evaluate_correction_velocity_laplacian_fd6(
            correction, x, y, z, t, spatial_step=step
        )
        for step in FD6_STEPS
    ]
    laplacians = [np.asarray(level["laplacian"], dtype=float) for level in levels]
    differences = [
        _vector_rms(laplacians[0] - laplacians[1]),
        _vector_rms(laplacians[1] - laplacians[2]),
    ]
    refinement = differences[0] / max(differences[1], np.finfo(float).tiny)
    finest_rms = _vector_rms(laplacians[-1])
    finest_max = float(np.max(np.linalg.norm(laplacians[-1], axis=-1)))

    mx = np.asarray((-0.43, 0.27, 0.51, -0.35), dtype=float)
    my = np.asarray((0.21, -0.39, 0.31, 0.46), dtype=float)
    mz = np.asarray((0.34, -0.23, -0.48, 0.19), dtype=float)
    mt = np.asarray((0.39, 0.45, 0.55, 0.61), dtype=float)
    manufactured = np.sum(
        _fd6_second_derivatives(
            _manufactured_provider, mx, my, mz, mt, spatial_step=0.011
        ),
        axis=-1,
    )
    manufactured_exact = _manufactured_laplacian(mx, my, mz, mt)
    manufactured_abs_max = float(np.max(np.abs(manufactured - manufactured_exact)))
    manufactured_scale = max(float(np.max(np.abs(manufactured_exact))), 1.0)
    manufactured_rel_max = manufactured_abs_max / manufactured_scale

    exterior = evaluate_correction_velocity_laplacian_fd6(
        correction,
        np.asarray((0.0, 0.20, 1.35, 0.74), dtype=float),
        np.zeros(4, dtype=float),
        np.asarray((0.0, 0.0, 0.0, 2.05), dtype=float),
        np.full(4, 0.50, dtype=float),
        spatial_step=FD6_STEPS[0],
    )
    exterior_abs_max = max(
        float(np.max(np.abs(np.asarray(exterior[key], dtype=float))))
        for key in ("velocity", "second_derivatives", "laplacian")
    )

    guards = {
        "laplacian_nontrivial_rms": finest_rms >= 1.0e-8,
        # The correction profile is only C2 across spline knots.  The cloud is
        # mid-cell, but this guard asks only for clear multi-resolution
        # stabilization rather than claiming global sixth-order regularity.
        "successive_difference_refinement": refinement >= 4.0,
        "manufactured_relative_max_error": manufactured_rel_max <= 5.0e-8,
        "support_exterior": exterior_abs_max <= 1.0e-12,
    }

    return {
        "task": TASK,
        "schema": SCHEMA,
        "parent_agent2_pr": PARENT_AGENT2_PR,
        "parent_agent2_head": PARENT_AGENT2_HEAD,
        "source_corrected_reader_commit": SOURCE_CORRECTED_READER_COMMIT,
        "source_corrected_reader_date": SOURCE_CORRECTED_READER_DATE,
        "sample_count": int(t.size),
        "verification_times": sorted({float(value) for value in t}),
        "fd6_steps": list(FD6_STEPS),
        "laplacian_rms": [_vector_rms(value) for value in laplacians],
        "laplacian_sampled_max": [
            float(np.max(np.linalg.norm(value, axis=-1))) for value in laplacians
        ],
        "successive_difference_rms": differences,
        "successive_difference_refinement_ratio": refinement,
        "finest_laplacian_rms": finest_rms,
        "finest_laplacian_sampled_max": finest_max,
        "manufactured_polynomial": {
            "step": 0.011,
            "degree_bound": 6,
            "absolute_max_error": manufactured_abs_max,
            "relative_max_error": manufactured_rel_max,
        },
        "support_exterior_absolute_max": exterior_abs_max,
        "guards": guards,
        "failed_guards": [name for name, passed in guards.items() if not passed],
        "provenance": {
            "source_structure": (
                "inherits signed localized complete-curl structure from the pinned "
                "corrected Kokuno reconstruction"
            ),
            "repository_realization": (
                "public-z pullback, compact quintic radial lift, autonomous signed "
                "time modulation, and Cartesian FD6 second-derivative diagnostic"
            ),
            "paper_exact_hidden_data_claimed": False,
            "agent3_lane_reimplemented": False,
        },
        "truth_boundary": truth_boundary(),
    }


def truth_boundary() -> dict[str, Any]:
    signature = inspect.signature(evaluate_correction_velocity_laplacian_fd6)
    forbidden = {
        "residual",
        "defect",
        "stress",
        "pressure",
        "forcing",
        "target",
        "gain",
        "normalized_score",
        "scientific_threshold",
        "delta_y",
    }
    return {
        "correction_velocity_laplacian_executable": True,
        "viscous_term_ready_for_future_correction_diagnostics": True,
        "vector_potential_first_complete_curl_reused": True,
        "forbidden_public_parameters_absent": forbidden.isdisjoint(signature.parameters),
        "candidate_amplitude_phase_support_retuned": False,
        "agent3_mean_radial_chain_reimplemented": False,
        "source_agent2_complete_curl_certified_for_full_candidate": False,
        "independent_agent4_correction_vector_potential_audit_required": True,
        "real_agent3_delta_a_bound": False,
        "real_full_candidate_correction_cycle_run": False,
        "heldout_ns_momentum_residual_assessed": False,
        "residual_reduction_claimed": False,
        "paper_exact": False,
        "pde_validated": False,
        "openai_field_identified": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    receipt = verification_receipt()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2, sort_keys=True))
    if receipt["failed_guards"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
