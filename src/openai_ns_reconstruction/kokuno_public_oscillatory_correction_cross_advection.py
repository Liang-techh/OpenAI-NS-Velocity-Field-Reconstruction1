"""Mixed nonlinear advection seam for frozen A2 oscillation and correction.

This module adds no new oscillatory or correction degree of freedom.  It combines
only already-materialized Agent-2 fields:

    N_cross = (u_osc . grad) delta_u + (delta_u . grad) u_osc.

The oscillatory velocity remains the frozen public complete-curl field.  The
correction remains vector-potential first and is realized through the existing
signed-amplitude complete-curl kernel.  The mixed term is useful later when a
real full candidate can be assembled, but it is not a complete Navier--Stokes
residual and it is not averaged into the Agent-3 mean/radial lane here.

Provenance boundary
-------------------
The corrected Kokuno 2026-09-09 reader is structural provenance for the
upstream localized oscillatory / complete-curl construction.  The public-z
pullback, compact correction spline, autonomous time lift, Cartesian FD6/FD4
Jacobians, and the mixed quadratic composition below are repository
realizations, not paper-exact hidden data.
"""
from __future__ import annotations

import argparse
import inspect
import json
import math
from pathlib import Path
from typing import Any, Callable

import numpy as np

from .kokuno_public_oscillatory_vorticity_diagnostic import (
    evaluate_vorticity_osc_fd6,
)
from .kokuno_public_signed_amplitude_complete_curl import (
    SignedRadialAmplitudeCorrection,
    correction_velocity,
    profile_from_delta_a,
)
from .kokuno_public_signed_amplitude_spatial_jacobian import (
    evaluate_correction_velocity_jacobian_fd4,
)
from .kokuno_public_z_pullback_velocity import velocity_osc

TASK = "KOKUNO-A2-OSCILLATORY-CORRECTION-CROSS-ADVECTION-056"
SCHEMA = "kokuno-a2-oscillatory-correction-cross-advection-v1"
PARENT_AGENT2_PR = 779
PARENT_AGENT2_HEAD = "a20878fb781fe75698ab22cfbe5c36a78f33ddc5"
SOURCE_CORRECTED_READER_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_CORRECTED_READER_DATE = "2026-09-09"
OSCILLATORY_FD6_STEP = 0.001
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
    if not all(np.all(np.isfinite(a)) for a in arrays):
        raise ValueError("x, y, z and t must be finite and broadcastable")
    return tuple(arrays)


def _validate_step(step: float, name: str) -> float:
    h = float(step)
    if not math.isfinite(h) or not (1.0e-5 <= h <= 5.0e-2):
        raise ValueError(f"{name} must be finite and lie in [1e-5, 5e-2]")
    return h


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


def _cross_advection_from_jacobians(
    oscillatory_velocity: np.ndarray,
    oscillatory_jacobian: np.ndarray,
    correction_velocity_value: np.ndarray,
    correction_jacobian: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return both mixed terms and their sum for J[component, axis]."""
    uo = np.asarray(oscillatory_velocity, dtype=float)
    jo = np.asarray(oscillatory_jacobian, dtype=float)
    uc = np.asarray(correction_velocity_value, dtype=float)
    jc = np.asarray(correction_jacobian, dtype=float)
    if uo.shape[-1:] != (3,) or uc.shape != uo.shape:
        raise ValueError("oscillatory and correction velocities must share shape (...,3)")
    expected_j = uo.shape[:-1] + (3, 3)
    if jo.shape != expected_j or jc.shape != expected_j:
        raise ValueError("Jacobians must have shape velocity.shape[:-1] + (3,3)")
    if not all(np.all(np.isfinite(a)) for a in (uo, jo, uc, jc)):
        raise ValueError("velocities and Jacobians must be finite")

    oscillation_advects_correction = np.einsum("...ij,...j->...i", jc, uo)
    correction_advects_oscillation = np.einsum("...ij,...j->...i", jo, uc)
    cross = oscillation_advects_correction + correction_advects_oscillation
    return oscillation_advects_correction, correction_advects_oscillation, cross


def evaluate_oscillatory_correction_cross_advection(
    correction: SignedRadialAmplitudeCorrection,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    *,
    oscillatory_spatial_step: float = OSCILLATORY_FD6_STEP,
    correction_spatial_step: float = CORRECTION_FD4_STEP,
) -> dict[str, np.ndarray | float | str]:
    """Evaluate the raw A2 oscillatory/correction mixed advection term.

    No residual, mean, pressure, forcing, damping, target, gain, ``delta_y`` or
    raw ``delta_a`` is accepted.  The caller must supply an already-materialized
    typed correction witness.
    """
    ho = _validate_step(oscillatory_spatial_step, "oscillatory_spatial_step")
    hc = _validate_step(correction_spatial_step, "correction_spatial_step")
    osc = evaluate_vorticity_osc_fd6(x, y, z, t, spatial_step=ho)
    corr = evaluate_correction_velocity_jacobian_fd4(
        correction, x, y, z, t, spatial_step=hc
    )

    uo = np.asarray(osc["velocity"], dtype=float)
    jo = np.asarray(osc["velocity_gradient_fd6"], dtype=float)
    uc = np.asarray(corr["velocity"], dtype=float)
    jc = np.asarray(corr["jacobian"], dtype=float)
    term_oc, term_co, cross = _cross_advection_from_jacobians(uo, jo, uc, jc)
    if not np.all(np.isfinite(cross)):
        raise RuntimeError("oscillatory/correction cross advection is non-finite")

    return {
        "oscillatory_velocity": uo,
        "correction_velocity": uc,
        "oscillatory_jacobian": jo,
        "correction_jacobian": jc,
        "oscillation_advects_correction": term_oc,
        "correction_advects_oscillation": term_co,
        "cross_advection": cross,
        "oscillatory_spatial_step": ho,
        "correction_spatial_step": hc,
        "oscillatory_operator": "centered_cartesian_fd6_first_derivative",
        "correction_operator": "centered_cartesian_fd4_first_derivative",
        "jacobian_convention": "component_axis",
    }


def _directional_fd4_advection(
    provider: VectorProvider,
    advecting_velocity: np.ndarray,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    *,
    directional_step: float,
) -> np.ndarray:
    """Independent frozen-direction FD4 approximation of ``(a.grad)provider``."""
    h = _validate_step(directional_step, "directional_step")
    xb, yb, zb, tb = _broadcast_xyzt(x, y, z, t)
    a = np.asarray(advecting_velocity, dtype=float)
    if a.shape != xb.shape + (3,):
        raise ValueError("advecting_velocity must have shape broadcast(x,y,z,t)+(3,)")
    if not np.all(np.isfinite(a)):
        raise ValueError("advecting_velocity must be finite")

    speed = np.linalg.norm(a, axis=-1)
    direction = np.zeros_like(a)
    nonzero = speed > 1.0e-14
    direction[nonzero] = a[nonzero] / speed[nonzero, None]

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


def _independent_directional_cross(
    correction: SignedRadialAmplitudeCorrection,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    *,
    directional_step: float,
) -> np.ndarray:
    xb, yb, zb, tb = _broadcast_xyzt(x, y, z, t)
    uo = _vector_value(velocity_osc, xb, yb, zb, tb)

    def corr_provider(xx: Any, yy: Any, zz: Any, tt: Any) -> np.ndarray:
        return correction_velocity(correction, xx, yy, zz, tt)

    uc = _vector_value(corr_provider, xb, yb, zb, tb)
    osc_advects_corr = _directional_fd4_advection(
        corr_provider, uo, xb, yb, zb, tb, directional_step=directional_step
    )
    corr_advects_osc = _directional_fd4_advection(
        velocity_osc, uc, xb, yb, zb, tb, directional_step=directional_step
    )
    return osc_advects_corr + corr_advects_osc


def _vector_rms(value: np.ndarray) -> float:
    a = np.asarray(value, dtype=float)
    return float(np.sqrt(np.mean(np.sum(a * a, axis=-1))))


def _relative_rms(reference: np.ndarray, value: np.ndarray) -> float:
    return _vector_rms(np.asarray(value) - np.asarray(reference)) / max(
        _vector_rms(reference), np.finfo(float).tiny
    )


def _relative_sampled_max(reference: np.ndarray, value: np.ndarray) -> float:
    ref = np.asarray(reference, dtype=float)
    val = np.asarray(value, dtype=float)
    diff = np.linalg.norm(val - ref, axis=-1)
    scale = max(float(np.max(np.linalg.norm(ref, axis=-1))), np.finfo(float).tiny)
    return float(np.max(diff) / scale)


def _verification_profile() -> SignedRadialAmplitudeCorrection:
    radii = np.linspace(0.32, 1.18, 27)
    s = (radii - radii[0]) / (radii[-1] - radii[0])
    bump = np.sin(math.pi * s) ** 8
    delta_a = np.stack(
        (
            0.0077 * bump * (1.0 + 0.07 * np.cos(2.0 * math.pi * s)),
            -0.0055 * bump * (1.0 - 0.05 * np.sin(2.0 * math.pi * s)),
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
            "fresh compact signed radial mechanics profile for A2 mixed-advection "
            "verification only; no Agent-3 defect/mean/stress/inverse/residual input"
        ),
        source_mean_amplitude_differential_certified=False,
    )


def _verification_cloud() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    r0 = 0.32
    dr = (1.18 - 0.32) / 26.0
    cell_ids = (2, 6, 10, 14, 18, 22)
    radii = np.asarray([r0 + (i + 0.5) * dr for i in cell_ids], dtype=float)
    z_block = np.asarray((-0.55, -0.33, -0.09, 0.15, 0.38, 0.56), dtype=float)
    golden = math.pi * (3.0 - math.sqrt(5.0))
    rows: list[tuple[float, float, float, float]] = []
    count = 0
    for time in (0.43, 0.50, 0.57):
        for radius, z_value in zip(radii, z_block, strict=True):
            angle = 0.311 + (count + 1) * golden
            rows.append(
                (
                    radius * math.cos(angle),
                    radius * math.sin(angle),
                    float(z_value),
                    time,
                )
            )
            count += 1
    a = np.asarray(rows, dtype=float)
    return a[:, 0], a[:, 1], a[:, 2], a[:, 3]


def _manufactured_algebra_check() -> dict[str, float]:
    uo = np.asarray(((1.0, -2.0, 0.5), (-0.4, 1.2, 2.1)), dtype=float)
    uc = np.asarray(((0.3, 0.8, -1.1), (1.4, -0.6, 0.2)), dtype=float)
    jo = np.asarray(
        (
            ((1.0, 2.0, 3.0), (0.5, -1.0, 0.0), (2.0, 0.25, -0.5)),
            ((-1.0, 0.3, 1.2), (2.2, -0.7, 0.4), (0.1, 1.5, -2.0)),
        ),
        dtype=float,
    )
    jc = np.asarray(
        (
            ((0.2, -0.4, 1.1), (1.3, 0.6, -0.8), (-0.5, 0.9, 0.7)),
            ((0.8, 1.1, -0.2), (-0.3, 0.4, 1.7), (1.6, -1.2, 0.5)),
        ),
        dtype=float,
    )
    a, b, cross = _cross_advection_from_jacobians(uo, jo, uc, jc)
    expected_a = np.stack([jc[i] @ uo[i] for i in range(2)], axis=0)
    expected_b = np.stack([jo[i] @ uc[i] for i in range(2)], axis=0)
    expected = expected_a + expected_b
    return {
        "oscillation_advects_correction_absolute_max_error": float(
            np.max(np.abs(a - expected_a))
        ),
        "correction_advects_oscillation_absolute_max_error": float(
            np.max(np.abs(b - expected_b))
        ),
        "cross_absolute_max_error": float(np.max(np.abs(cross - expected))),
    }


def build_receipt() -> dict[str, Any]:
    correction = _verification_profile()
    x, y, z, t = _verification_cloud()
    production = evaluate_oscillatory_correction_cross_advection(
        correction,
        x,
        y,
        z,
        t,
        oscillatory_spatial_step=OSCILLATORY_FD6_STEP,
        correction_spatial_step=CORRECTION_FD4_STEP,
    )
    reference = np.asarray(production["cross_advection"], dtype=float)
    directional = [
        _independent_directional_cross(
            correction, x, y, z, t, directional_step=step
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

    exterior_x = np.asarray((2.20, -2.24, 0.0, 1.72), dtype=float)
    exterior_y = np.asarray((0.0, 0.0, 2.22, 1.72), dtype=float)
    exterior_z = np.asarray((0.0, 0.0, 0.0, 2.20), dtype=float)
    exterior_t = np.full(4, 0.50, dtype=float)
    exterior = evaluate_oscillatory_correction_cross_advection(
        correction,
        exterior_x,
        exterior_y,
        exterior_z,
        exterior_t,
        oscillatory_spatial_step=OSCILLATORY_FD6_STEP,
        correction_spatial_step=CORRECTION_FD4_STEP,
    )
    exterior_abs_max = max(
        float(np.max(np.abs(np.asarray(exterior[key], dtype=float))))
        for key in (
            "oscillatory_velocity",
            "correction_velocity",
            "oscillation_advects_correction",
            "correction_advects_oscillation",
            "cross_advection",
        )
    )

    manufactured = _manufactured_algebra_check()
    signature = inspect.signature(evaluate_oscillatory_correction_cross_advection)
    forbidden = {
        "residual", "defect", "mean", "stress", "pressure", "forcing",
        "target", "gain", "alpha", "damping", "delta_y", "delta_a",
        "normalized_score", "scientific_threshold",
    }
    failed: list[str] = []
    if _vector_rms(reference) < 1.0e-10:
        failed.append("cross_advection_nontrivial")
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
            "structural_scope": "localized oscillatory and signed complete-curl fields",
        },
        "sample_count": int(t.size),
        "sample_times": [0.43, 0.50, 0.57],
        "oscillatory_fd6_step": OSCILLATORY_FD6_STEP,
        "correction_fd4_step": CORRECTION_FD4_STEP,
        "directional_fd4_steps": list(DIRECTIONAL_FD4_STEPS),
        "cross_advection_rms": _vector_rms(reference),
        "directional_successive_difference_rms": successive,
        "directional_refinement_ratio": refinement,
        "directional_finest_relative_rms": relative_rms,
        "directional_finest_relative_sampled_max": relative_max,
        "manufactured_algebra": manufactured,
        "support_exterior_absolute_max": exterior_abs_max,
        "failed_guards": failed,
        "truth_boundary": {
            "oscillatory_velocity_candidate_changed": False,
            "correction_velocity_candidate_changed": False,
            "amplitude_phase_carrier_support_retuned": False,
            "complete_curl_correction_reused": True,
            "cross_advection_executable": True,
            "cross_advection_is_ns_residual": False,
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
        raise SystemExit("mixed-advection guards failed: " + ", ".join(receipt["failed_guards"]))


if __name__ == "__main__":
    main()
