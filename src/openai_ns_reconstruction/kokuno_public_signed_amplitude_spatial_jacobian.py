"""FD4 spatial-Jacobian seam for the frozen A2 signed-amplitude correction.

This Agent-2 module consumes an already-materialized
``SignedRadialAmplitudeCorrection`` and exposes the Cartesian first spatial
derivatives of ``delta u``.  The upstream correction remains vector-potential
first and is produced by the localized complete-curl construction; this module
does not fit Cartesian velocity components and does not consume Agent-3
residual, defect, mean, stress, inverse, damping, pressure, forcing, target, or
gain data.

The corrected Kokuno reconstruction supplies structural provenance for the
localized signed-amplitude complete-curl construction.  The public-z pullback,
compact radial lift, autonomous temporal lift, and Cartesian FD4 derivative
operator used by the executable repository are repository realizations rather
than paper-exact hidden data.
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

TASK = "KOKUNO-A2-SIGNED-AMPLITUDE-SPATIAL-JACOBIAN-054"
SCHEMA = "kokuno-a2-signed-amplitude-spatial-jacobian-v1"
PARENT_AGENT2_PR = 761
PARENT_AGENT2_HEAD = "d36dd9ecde8c0c69e06c93b412bf3f9e091391fe"
SOURCE_CORRECTED_READER_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_CORRECTED_READER_DATE = "2026-09-09"
FD4_STEPS = (0.004, 0.002, 0.001)

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


def _fd4_first_derivatives(
    provider: VectorProvider,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    *,
    spatial_step: float,
) -> np.ndarray:
    """Return ``[..., component, axis] = d_axis provider_component``."""
    h = _validate_step(spatial_step)
    xb, yb, zb, tb = _broadcast_xyzt(x, y, z, t)
    base = (xb, yb, zb)
    derivatives: list[np.ndarray] = []
    for axis in range(3):
        values: dict[int, np.ndarray] = {}
        for offset in (-2, -1, 1, 2):
            shifted = [array for array in base]
            shifted[axis] = shifted[axis] + offset * h
            values[offset] = _provider_value(
                provider, shifted[0], shifted[1], shifted[2], tb
            )
        derivative = (
            values[-2] - 8.0 * values[-1] + 8.0 * values[1] - values[2]
        ) / (12.0 * h)
        derivatives.append(derivative)
    return np.stack(derivatives, axis=-1)


def _vorticity_from_jacobian(jacobian: np.ndarray) -> np.ndarray:
    j = np.asarray(jacobian, dtype=float)
    if j.shape[-2:] != (3, 3):
        raise ValueError("jacobian must end in (component, axis) = (3, 3)")
    return np.stack(
        (
            j[..., 2, 1] - j[..., 1, 2],
            j[..., 0, 2] - j[..., 2, 0],
            j[..., 1, 0] - j[..., 0, 1],
        ),
        axis=-1,
    )


def evaluate_correction_velocity_jacobian_fd4(
    correction: SignedRadialAmplitudeCorrection,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    *,
    spatial_step: float = 0.002,
) -> dict[str, np.ndarray | float | str]:
    """Evaluate ``delta u``, Cartesian Jacobian, divergence, and vorticity.

    The Jacobian convention is ``[..., component, axis]`` with axis order
    ``(x, y, z)``.  This makes the tensor directly usable in future NS
    convective cross terms without changing the frozen complete-curl field.
    """
    h = _validate_step(spatial_step)
    xb, yb, zb, tb = _broadcast_xyzt(x, y, z, t)

    def provider(xx: Any, yy: Any, zz: Any, tt: Any) -> np.ndarray:
        return correction_velocity(correction, xx, yy, zz, tt)

    velocity = _provider_value(provider, xb, yb, zb, tb)
    jacobian = _fd4_first_derivatives(
        provider, xb, yb, zb, tb, spatial_step=h
    )
    divergence = np.trace(jacobian, axis1=-2, axis2=-1)
    vorticity = _vorticity_from_jacobian(jacobian)
    return {
        "velocity": velocity,
        "jacobian": jacobian,
        "divergence": divergence,
        "vorticity": vorticity,
        "spatial_step": h,
        "spatial_operator": "centered_cartesian_fd4_first_derivative",
        "jacobian_convention": "component_axis",
    }


def _tensor_rms(values: np.ndarray) -> float:
    arr = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(np.sum(arr * arr, axis=(-2, -1)))))


def _vector_rms(values: np.ndarray) -> float:
    arr = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(np.sum(arr * arr, axis=-1))))


def _scalar_rms(values: np.ndarray) -> float:
    arr = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(arr * arr)))


def _verification_profile() -> SignedRadialAmplitudeCorrection:
    radii = np.linspace(0.31, 1.19, 25)
    s = (radii - radii[0]) / (radii[-1] - radii[0])
    bump = np.sin(math.pi * s) ** 8
    delta_a = np.stack(
        (
            0.0089 * bump * (1.0 + 0.08 * np.cos(2.0 * math.pi * s)),
            -0.0063 * bump * (1.0 - 0.06 * np.sin(2.0 * math.pi * s)),
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
            "fresh compact signed radial mechanics profile for correction-Jacobian "
            "stability only; no Agent-3 defect, stress, inverse, damping, or residual input"
        ),
        source_mean_amplitude_differential_certified=False,
    )


def _verification_cloud() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    # Mid-cell radii keep the full +/-2h FD4 reach away from radial spline knots.
    r0 = 0.31
    dr = (1.19 - 0.31) / 24.0
    cell_ids = (2, 6, 10, 13, 17, 21)
    radii = np.asarray([r0 + (index + 0.5) * dr for index in cell_ids], dtype=float)
    z_block = np.asarray((-0.57, -0.35, -0.11, 0.13, 0.36, 0.58), dtype=float)
    golden = math.pi * (3.0 - math.sqrt(5.0))
    rows: list[tuple[float, float, float, float]] = []
    count = 0
    for time in (0.41, 0.50, 0.59):
        for radius, z_value in zip(radii, z_block, strict=True):
            angle = 0.237 + (count + 1) * golden
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
            xb**4 + 2.0 * xb * yb**2 - 0.7 * yb * zb**2 + 0.3 * zb,
            -0.4 * xb**2 * yb + yb**4 + 1.1 * xb * zb - 0.2 * zb**3,
            0.6 * xb * yb * zb + 0.5 * yb**3 + zb**4 - 0.9 * xb,
        ),
        axis=-1,
    )


def _manufactured_jacobian(x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
    xb, yb, zb, tb = _broadcast_xyzt(x, y, z, t)
    del tb
    row0 = np.stack(
        (4.0 * xb**3 + 2.0 * yb**2, 4.0 * xb * yb - 0.7 * zb**2, -1.4 * yb * zb + 0.3),
        axis=-1,
    )
    row1 = np.stack(
        (-0.8 * xb * yb + 1.1 * zb, -0.4 * xb**2 + 4.0 * yb**3, 1.1 * xb - 0.6 * zb**2),
        axis=-1,
    )
    row2 = np.stack(
        (0.6 * yb * zb - 0.9, 0.6 * xb * zb + 1.5 * yb**2, 0.6 * xb * yb + 4.0 * zb**3),
        axis=-1,
    )
    return np.stack((row0, row1, row2), axis=-2)


def verification_receipt() -> dict[str, Any]:
    correction = _verification_profile()
    x, y, z, t = _verification_cloud()
    levels = [
        evaluate_correction_velocity_jacobian_fd4(
            correction, x, y, z, t, spatial_step=step
        )
        for step in FD4_STEPS
    ]
    jacobians = [np.asarray(level["jacobian"], dtype=float) for level in levels]
    differences = [
        _tensor_rms(jacobians[0] - jacobians[1]),
        _tensor_rms(jacobians[1] - jacobians[2]),
    ]
    refinement = differences[0] / max(differences[1], np.finfo(float).tiny)
    finest_jacobian_rms = _tensor_rms(jacobians[-1])
    finest_divergence = np.asarray(levels[-1]["divergence"], dtype=float)
    finest_divergence_rms = _scalar_rms(finest_divergence)
    finest_divergence_max = float(np.max(np.abs(finest_divergence)))
    relative_divergence_rms = finest_divergence_rms / max(
        finest_jacobian_rms, np.finfo(float).tiny
    )
    vorticity_rms = [_vector_rms(np.asarray(level["vorticity"], dtype=float)) for level in levels]

    mx = np.asarray((-0.43, 0.27, 0.51, -0.35), dtype=float)
    my = np.asarray((0.21, -0.39, 0.31, 0.46), dtype=float)
    mz = np.asarray((0.34, -0.23, -0.48, 0.19), dtype=float)
    mt = np.asarray((0.41, 0.47, 0.53, 0.59), dtype=float)
    manufactured = _fd4_first_derivatives(
        _manufactured_provider, mx, my, mz, mt, spatial_step=0.011
    )
    manufactured_exact = _manufactured_jacobian(mx, my, mz, mt)
    manufactured_abs_max = float(np.max(np.abs(manufactured - manufactured_exact)))
    manufactured_scale = max(float(np.max(np.abs(manufactured_exact))), 1.0)
    manufactured_rel_max = manufactured_abs_max / manufactured_scale

    exterior = evaluate_correction_velocity_jacobian_fd4(
        correction,
        np.asarray((0.0, 0.18, 1.36, 0.72), dtype=float),
        np.zeros(4, dtype=float),
        np.asarray((0.0, 0.0, 0.0, 2.05), dtype=float),
        np.full(4, 0.50, dtype=float),
        spatial_step=FD4_STEPS[0],
    )
    exterior_abs_max = max(
        float(np.max(np.abs(np.asarray(exterior[key], dtype=float))))
        for key in ("velocity", "jacobian", "divergence", "vorticity")
    )

    guards = {
        "jacobian_nontrivial_rms": finest_jacobian_rms >= 1.0e-8,
        # FD4 is formally fourth order inside a smooth spline cell.  The public
        # correction is only C2 across cell boundaries, so we preregister a
        # conservative stabilization ratio instead of claiming global order 4.
        "successive_difference_refinement": refinement >= 4.0,
        "divergence_identity_relative_rms": relative_divergence_rms <= 5.0e-3,
        "manufactured_relative_max_error": manufactured_rel_max <= 5.0e-9,
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
        "fd4_steps": list(FD4_STEPS),
        "jacobian_rms": [_tensor_rms(value) for value in jacobians],
        "successive_difference_rms": differences,
        "successive_difference_refinement_ratio": refinement,
        "finest_jacobian_rms": finest_jacobian_rms,
        "finest_divergence_rms": finest_divergence_rms,
        "finest_divergence_sampled_max": finest_divergence_max,
        "finest_relative_divergence_rms": relative_divergence_rms,
        "vorticity_rms": vorticity_rms,
        "manufactured_polynomial": {
            "step": 0.011,
            "degree_bound": 4,
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
                "time modulation, and Cartesian FD4 first-spatial-derivative diagnostic"
            ),
            "paper_exact_hidden_data_claimed": False,
            "agent3_lane_reimplemented": False,
        },
        "truth_boundary": truth_boundary(),
    }


def truth_boundary() -> dict[str, Any]:
    signature = inspect.signature(evaluate_correction_velocity_jacobian_fd4)
    forbidden = {
        "residual",
        "defect",
        "mean",
        "stress",
        "inverse",
        "damping",
        "pressure",
        "forcing",
        "target",
        "gain",
        "normalized_score",
        "scientific_threshold",
        "delta_y",
        "delta_a",
    }
    return {
        "correction_velocity_spatial_jacobian_executable": True,
        "convective_cross_term_derivatives_ready_for_future_diagnostics": True,
        "divergence_and_vorticity_derived_from_same_jacobian": True,
        "vector_potential_first_complete_curl_reused": True,
        "forbidden_public_parameters_absent": forbidden.isdisjoint(signature.parameters),
        "candidate_amplitude_phase_support_retuned": False,
        "agent3_mean_radial_chain_reimplemented": False,
        "source_agent2_complete_curl_certified_for_full_candidate": False,
        "independent_agent4_correction_vector_potential_audit_required": True,
        "real_agent3_delta_a_bound": False,
        "full_composite_velocity_available": False,
        "real_full_candidate_correction_cycle_run": False,
        "heldout_ns_momentum_residual_assessed": False,
        "residual_reduction_claimed": False,
        "paper_exact": False,
        "pde_validated": False,
        "openai_field_identified": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    receipt = verification_receipt()
    payload = json.dumps(receipt, indent=2, sort_keys=True)
    if args.output is None:
        print(payload)
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload + "\n", encoding="utf-8")
    if receipt["failed_guards"]:
        raise SystemExit(f"failed guards: {receipt['failed_guards']}")


if __name__ == "__main__":
    main()
