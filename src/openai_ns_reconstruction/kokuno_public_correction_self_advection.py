"""Nonlinear self-advection seam for the frozen A2 complete-curl correction.

This module adds no correction degree of freedom.  It consumes an already
materialized ``SignedRadialAmplitudeCorrection`` and evaluates

    N_corr = (delta_u . grad) delta_u.

The correction velocity remains vector-potential first and is generated only by
the existing signed-amplitude complete-curl kernel.  The Cartesian Jacobian is
the existing Agent-2 FD4 diagnostic.  This quadratic term is a reusable local
piece of a future composite momentum decomposition; it is not a complete
Navier--Stokes residual and it is not projected into Agent 3's mean/radial lane
here.

Provenance boundary
-------------------
The corrected Kokuno 2026-09-09 reader is structural provenance for the
upstream localized signed-amplitude complete-curl construction.  The compact
radial lift, public-z pullback, autonomous temporal lift, Cartesian FD4
Jacobian, independent directional verifier, and quadratic composition below
are repository realizations, not paper-exact hidden data.
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
from .kokuno_public_signed_amplitude_spatial_jacobian import (
    evaluate_correction_velocity_jacobian_fd4,
)

TASK = "KOKUNO-A2-CORRECTION-SELF-ADVECTION-057"
SCHEMA = "kokuno-a2-correction-self-advection-v1"
PARENT_AGENT2_PR = 788
PARENT_AGENT2_HEAD = "c1e45c12d1842d8e401a4b54ee926576aacfb7f5"
SOURCE_CORRECTED_READER_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_CORRECTED_READER_DATE = "2026-09-09"
CORRECTION_FD4_STEP = 0.001
DIRECTIONAL_FD4_STEPS = (0.004, 0.002, 0.001)

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


def _validate_step(step: float, name: str) -> float:
    value = float(step)
    if not math.isfinite(value) or not (1.0e-5 <= value <= 5.0e-2):
        raise ValueError(f"{name} must be finite and lie in [1e-5, 5e-2]")
    return value


def _vector_value(
    provider: VectorProvider,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
) -> np.ndarray:
    value = np.asarray(provider(x, y, z, t), dtype=float)
    expected = x.shape + (3,)
    if value.shape != expected:
        raise RuntimeError(f"vector provider returned {value.shape}, expected {expected}")
    if not np.all(np.isfinite(value)):
        raise RuntimeError("vector provider returned non-finite values")
    return value


def _self_advection_from_jacobian(
    velocity: np.ndarray,
    jacobian: np.ndarray,
) -> np.ndarray:
    """Return ``(u.grad)u`` for ``J[component, axis]=d_axis u_component``."""
    u = np.asarray(velocity, dtype=float)
    j = np.asarray(jacobian, dtype=float)
    if u.shape[-1:] != (3,):
        raise ValueError("velocity must have shape (...,3)")
    expected = u.shape[:-1] + (3, 3)
    if j.shape != expected:
        raise ValueError("jacobian must have shape velocity.shape[:-1] + (3,3)")
    if not np.all(np.isfinite(u)) or not np.all(np.isfinite(j)):
        raise ValueError("velocity and jacobian must be finite")
    return np.einsum("...ij,...j->...i", j, u)


def evaluate_correction_self_advection_fd4(
    correction: SignedRadialAmplitudeCorrection,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    *,
    spatial_step: float = CORRECTION_FD4_STEP,
) -> dict[str, np.ndarray | float | str]:
    """Evaluate ``delta_u`` and ``(delta_u.grad)delta_u``.

    The caller must provide an already-materialized typed correction witness.
    No residual, defect, mean, stress, pressure, forcing, damping target, gain,
    raw ``delta_y`` or raw ``delta_a`` is accepted.
    """
    h = _validate_step(spatial_step, "spatial_step")
    result = evaluate_correction_velocity_jacobian_fd4(
        correction, x, y, z, t, spatial_step=h
    )
    velocity = np.asarray(result["velocity"], dtype=float)
    jacobian = np.asarray(result["jacobian"], dtype=float)
    self_advection = _self_advection_from_jacobian(velocity, jacobian)
    if not np.all(np.isfinite(self_advection)):
        raise RuntimeError("correction self-advection is non-finite")
    return {
        "correction_velocity": velocity,
        "correction_jacobian": jacobian,
        "correction_divergence": np.asarray(result["divergence"], dtype=float),
        "correction_vorticity": np.asarray(result["vorticity"], dtype=float),
        "correction_self_advection": self_advection,
        "spatial_step": h,
        "spatial_operator": "centered_cartesian_fd4_first_derivative",
        "jacobian_convention": "component_axis",
    }


def _directional_fd4_self_advection(
    provider: VectorProvider,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    *,
    directional_step: float,
) -> np.ndarray:
    """Independent frozen-direction FD4 approximation of ``(u.grad)u``."""
    h = _validate_step(directional_step, "directional_step")
    xb, yb, zb, tb = _broadcast_xyzt(x, y, z, t)
    velocity = _vector_value(provider, xb, yb, zb, tb)
    speed = np.linalg.norm(velocity, axis=-1)
    direction = np.zeros_like(velocity)
    nonzero = speed > 1.0e-14
    direction[nonzero] = velocity[nonzero] / speed[nonzero, None]

    values: dict[int, np.ndarray] = {}
    for offset in (-2, -1, 1, 2):
        displacement = offset * h * direction
        values[offset] = _vector_value(
            provider,
            xb + displacement[..., 0],
            yb + displacement[..., 1],
            zb + displacement[..., 2],
            tb,
        )
    derivative = (
        values[-2] - 8.0 * values[-1] + 8.0 * values[1] - values[2]
    ) / (12.0 * h)
    return speed[..., None] * derivative


def _vector_rms(value: np.ndarray) -> float:
    array = np.asarray(value, dtype=float)
    return float(np.sqrt(np.mean(np.sum(array * array, axis=-1))))


def _relative_rms(reference: np.ndarray, value: np.ndarray) -> float:
    return _vector_rms(np.asarray(value) - np.asarray(reference)) / max(
        _vector_rms(reference), np.finfo(float).tiny
    )


def _relative_sampled_max(reference: np.ndarray, value: np.ndarray) -> float:
    ref = np.asarray(reference, dtype=float)
    val = np.asarray(value, dtype=float)
    error = np.linalg.norm(val - ref, axis=-1)
    scale = max(float(np.max(np.linalg.norm(ref, axis=-1))), np.finfo(float).tiny)
    return float(np.max(error) / scale)


def _verification_profile() -> SignedRadialAmplitudeCorrection:
    radii = np.linspace(0.30, 1.20, 29)
    s = (radii - radii[0]) / (radii[-1] - radii[0])
    bump = np.sin(math.pi * s) ** 8
    delta_a = np.stack(
        (
            0.0081 * bump * (1.0 + 0.06 * np.cos(2.0 * math.pi * s)),
            -0.0059 * bump * (1.0 - 0.04 * np.sin(2.0 * math.pi * s)),
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
            "fresh compact signed radial mechanics profile for A2 correction "
            "self-advection verification only; no Agent-3 defect, mean, stress, "
            "inverse, damping, or residual input"
        ),
        source_mean_amplitude_differential_certified=False,
    )


def _verification_cloud() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    # Cell-midpoint radii keep the full +/-2h directional stencil away from the
    # compact radial spline knots at the largest preregistered step.
    r0 = 0.30
    dr = (1.20 - 0.30) / 28.0
    cell_ids = (3, 7, 11, 16, 20, 24)
    radii = np.asarray([r0 + (index + 0.5) * dr for index in cell_ids], dtype=float)
    z_block = np.asarray((-0.56, -0.34, -0.10, 0.14, 0.37, 0.55), dtype=float)
    golden = math.pi * (3.0 - math.sqrt(5.0))
    rows: list[tuple[float, float, float, float]] = []
    count = 0
    for time in (0.44, 0.50, 0.56):
        for radius, z_value in zip(radii, z_block, strict=True):
            angle = 0.271 + (count + 1) * golden
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


def _manufactured_algebra_check() -> dict[str, float]:
    velocity = np.asarray(
        ((1.0, -2.0, 0.5), (-0.4, 1.2, 2.1), (0.3, -0.7, 1.4)),
        dtype=float,
    )
    jacobian = np.asarray(
        (
            ((1.0, 2.0, 3.0), (0.5, -1.0, 0.0), (2.0, 0.25, -0.5)),
            ((-1.0, 0.3, 1.2), (2.2, -0.7, 0.4), (0.1, 1.5, -2.0)),
            ((0.8, -0.2, 1.1), (-0.6, 1.7, 0.3), (1.4, 0.5, -0.9)),
        ),
        dtype=float,
    )
    actual = _self_advection_from_jacobian(velocity, jacobian)
    expected = np.stack([jacobian[i] @ velocity[i] for i in range(3)], axis=0)
    return {
        "self_advection_absolute_max_error": float(np.max(np.abs(actual - expected)))
    }


def build_receipt() -> dict[str, Any]:
    correction = _verification_profile()
    x, y, z, t = _verification_cloud()
    production = evaluate_correction_self_advection_fd4(
        correction, x, y, z, t, spatial_step=CORRECTION_FD4_STEP
    )
    reference = np.asarray(production["correction_self_advection"], dtype=float)

    def provider(xx: Any, yy: Any, zz: Any, tt: Any) -> np.ndarray:
        return correction_velocity(correction, xx, yy, zz, tt)

    directional = [
        _directional_fd4_self_advection(
            provider, x, y, z, t, directional_step=step
        )
        for step in DIRECTIONAL_FD4_STEPS
    ]
    successive = [
        _vector_rms(directional[0] - directional[1]),
        _vector_rms(directional[1] - directional[2]),
    ]
    refinement = successive[0] / max(successive[1], np.finfo(float).tiny)
    relative_rms = _relative_rms(reference, directional[-1])
    relative_max = _relative_sampled_max(reference, directional[-1])

    exterior = evaluate_correction_self_advection_fd4(
        correction,
        np.asarray((2.20, -2.24, 0.0, 1.72), dtype=float),
        np.asarray((0.0, 0.0, 2.22, 1.72), dtype=float),
        np.asarray((0.0, 0.0, 0.0, 2.20), dtype=float),
        np.full(4, 0.50, dtype=float),
        spatial_step=CORRECTION_FD4_STEP,
    )
    exterior_abs_max = max(
        float(np.max(np.abs(np.asarray(exterior[key], dtype=float))))
        for key in (
            "correction_velocity",
            "correction_jacobian",
            "correction_divergence",
            "correction_vorticity",
            "correction_self_advection",
        )
    )

    manufactured = _manufactured_algebra_check()
    signature = inspect.signature(evaluate_correction_self_advection_fd4)
    forbidden = {
        "residual", "defect", "mean", "stress", "inverse", "pressure",
        "forcing", "target", "gain", "alpha", "damping", "delta_y",
        "delta_a", "normalized_score", "scientific_threshold",
    }
    failed: list[str] = []
    if _vector_rms(reference) < 1.0e-12:
        failed.append("correction_self_advection_nontrivial")
    if refinement < 8.0:
        failed.append("directional_fd4_refinement")
    if relative_rms > 5.0e-3:
        failed.append("finest_directional_relative_rms")
    if relative_max > 1.0e-2:
        failed.append("finest_directional_relative_sampled_max")
    if max(manufactured.values()) > 1.0e-14:
        failed.append("manufactured_algebra")
    if exterior_abs_max > 1.0e-12:
        failed.append("support_exterior")
    if forbidden.intersection(signature.parameters):
        failed.append("forbidden_public_inputs")

    return {
        "schema": SCHEMA,
        "task": TASK,
        "parent_agent2_pr": PARENT_AGENT2_PR,
        "parent_agent2_head": PARENT_AGENT2_HEAD,
        "source_provenance": {
            "corrected_reader_commit": SOURCE_CORRECTED_READER_COMMIT,
            "corrected_reader_date": SOURCE_CORRECTED_READER_DATE,
            "structural_scope": "localized signed-amplitude complete-curl correction",
        },
        "sample_count": int(t.size),
        "sample_times": [0.44, 0.50, 0.56],
        "correction_fd4_step": CORRECTION_FD4_STEP,
        "directional_fd4_steps": list(DIRECTIONAL_FD4_STEPS),
        "correction_self_advection_rms": _vector_rms(reference),
        "directional_successive_difference_rms": successive,
        "directional_refinement_ratio": refinement,
        "directional_finest_relative_rms": relative_rms,
        "directional_finest_relative_sampled_max": relative_max,
        "manufactured_algebra": manufactured,
        "support_exterior_absolute_max": exterior_abs_max,
        "failed_guards": failed,
        "truth_boundary": {
            "correction_velocity_candidate_changed": False,
            "amplitude_phase_carrier_support_retuned": False,
            "vector_potential_complete_curl_reused": True,
            "correction_self_advection_executable": True,
            "correction_self_advection_is_ns_residual": False,
            "agent3_mean_radial_chain_reimplemented": False,
            "agent3_mean_defect_contract_emitted": False,
            "real_agent3_delta_a_bound": False,
            "real_global_leading_field_bound": False,
            "full_composite_velocity_available": False,
            "heldout_ns_momentum_residual_assessed": False,
            "residual_reduction_claimed": False,
            "paper_exact": False,
            "pde_validated": False,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    receipt = build_receipt()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    if receipt["failed_guards"]:
        raise SystemExit(
            "correction self-advection guards failed: "
            + ", ".join(receipt["failed_guards"])
        )


if __name__ == "__main__":
    main()
