"""Independent black-box audit for the materialized strict-inner Kokuno candidate.

Agent 4 consumes only a candidate-facing ``velocity(x,y,z,t)`` callable for
scientific numerical checks.  Spatial derivatives are reconstructed with a
Richardson-extrapolated centered two-point operator, deliberately distinct from
the Agent-1 analytic/FD6 paths and the Agent-2 FD4/FD8 transport diagnostics.

The scoped artifact is still only PA.10 inner contraction-center + frozen
complete-curl oscillation.  Pressure, restricted forcing, outer/global joining,
and the Agent-3 correction are absent, so this module cannot form or certify a
complete Navier--Stokes momentum residual.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Callable

import numpy as np

TASK = "KOKUNO-A4-STRICT-INNER-CANDIDATE-ARTIFACT-INDEPENDENT-AUDIT-076"
SCHEMA = "kokuno-a4-strict-inner-candidate-artifact-independent-audit-v1"
AGENT2_CANDIDATE_PR = 857
AGENT2_CANDIDATE_HEAD = "a586b7afe4bb47dbfd4b25177586c79d362afff9"
AGENT2_CANDIDATE_BLOB = "82e732af14695c39afa3b631c00b7f20960976da"
AGENT1_REFERENCE_PR = 848
AGENT1_REFERENCE_HEAD = "300919075c07396e4e466f553be78eaae0707df2"
AGENT1_REFERENCE_BLOB = "5076ffa69748c2d7e68289a155b67288eab6dab4"

RICHARDSON_STEPS = (1.6e-3, 8.0e-4, 4.0e-4)
DIVERGENCE_SAMPLED_MAX_GATE = 1.0e-5
DIVERGENCE_SAMPLED_RMS_GATE = 1.0e-5
RESOLUTION_STABILITY_FACTOR = 1.25
RESOLUTION_FLOOR = 2.0e-8
NONTRIVIAL_VELOCITY_RMS_FLOOR = 1.0e-8
FINAL_MOMENTUM_NORMALIZED_GATE = 1.0e-3
FINAL_DIVERGENCE_GATE = 1.0e-5

VelocityProvider = Callable[[Any, Any, Any, Any], np.ndarray]


@dataclass(frozen=True)
class RichardsonDivergenceLevel:
    """One independent black-box spatial-derivative resolution."""

    step: float
    velocity: np.ndarray
    jacobian: np.ndarray
    divergence: np.ndarray


def _broadcast_xyzt(x: Any, y: Any, z: Any, t: Any) -> tuple[np.ndarray, ...]:
    arrays = np.broadcast_arrays(
        np.asarray(x, dtype=float), np.asarray(y, dtype=float),
        np.asarray(z, dtype=float), np.asarray(t, dtype=float),
    )
    if not all(np.all(np.isfinite(a)) for a in arrays):
        raise ValueError("x, y, z and t must be finite and broadcastable")
    return tuple(arrays)


def _validate_step(step: float) -> float:
    value = float(step)
    if not math.isfinite(value) or not 1.0e-6 <= value <= 5.0e-2:
        raise ValueError("step must be finite and lie in [1e-6,5e-2]")
    return value


def _velocity_value(provider: VelocityProvider, x: np.ndarray, y: np.ndarray,
                    z: np.ndarray, t: np.ndarray) -> np.ndarray:
    value = np.asarray(provider(x, y, z, t), dtype=float)
    expected = x.shape + (3,)
    if value.shape != expected:
        raise ValueError(f"velocity provider must return shape {expected}, got {value.shape}")
    if not np.all(np.isfinite(value)):
        raise ValueError("velocity provider returned non-finite values")
    return value


def _centered_two_point(provider: VelocityProvider, x: np.ndarray, y: np.ndarray,
                        z: np.ndarray, t: np.ndarray, *, step: float,
                        axis: int) -> np.ndarray:
    h = _validate_step(step)
    if axis not in (0, 1, 2):
        raise ValueError("spatial axis must be 0, 1 or 2")
    coords = [x, y, z]
    plus = [a.copy() for a in coords]
    minus = [a.copy() for a in coords]
    plus[axis] = plus[axis] + h
    minus[axis] = minus[axis] - h
    up = _velocity_value(provider, plus[0], plus[1], plus[2], t)
    um = _velocity_value(provider, minus[0], minus[1], minus[2], t)
    return (up - um) / (2.0 * h)


def richardson_first_derivative(provider: VelocityProvider, x: np.ndarray,
                                y: np.ndarray, z: np.ndarray, t: np.ndarray,
                                *, step: float, axis: int) -> np.ndarray:
    """Fourth-order derivative from two independent centered two-point levels."""
    h = _validate_step(step)
    coarse = _centered_two_point(provider, x, y, z, t, step=h, axis=axis)
    half = _centered_two_point(provider, x, y, z, t, step=0.5 * h, axis=axis)
    return (4.0 * half - coarse) / 3.0


def richardson_divergence_level(velocity: VelocityProvider, x: Any, y: Any,
                                z: Any, t: Any, *, step: float) -> RichardsonDivergenceLevel:
    """Reconstruct Cartesian Jacobian/divergence from public velocity values only."""
    h = _validate_step(step)
    xb, yb, zb, tb = _broadcast_xyzt(x, y, z, t)
    base = _velocity_value(velocity, xb, yb, zb, tb)
    jac = np.empty(xb.shape + (3, 3), dtype=float)
    for axis in range(3):
        jac[..., :, axis] = richardson_first_derivative(
            velocity, xb, yb, zb, tb, step=h, axis=axis
        )
    divergence = np.trace(jac, axis1=-2, axis2=-1)
    if not np.all(np.isfinite(jac)) or not np.all(np.isfinite(divergence)):
        raise RuntimeError("independent Richardson spatial audit became non-finite")
    return RichardsonDivergenceLevel(h, base, jac, divergence)


def scalar_rms(value: Any) -> float:
    array = np.asarray(value, dtype=float)
    return float(np.sqrt(np.mean(array * array)))


def vector_rms(value: Any) -> float:
    array = np.asarray(value, dtype=float)
    if array.shape[-1:] != (3,):
        raise ValueError("vector_rms expects a final Cartesian component axis of length 3")
    return float(np.sqrt(np.mean(np.sum(array * array, axis=-1))))


def divergence_metrics(level: RichardsonDivergenceLevel) -> dict[str, float]:
    return {
        "step": float(level.step),
        "sampled_max": float(np.max(np.abs(level.divergence))),
        "sampled_rms": scalar_rms(level.divergence),
    }


def frozen_gate_failures(levels: tuple[RichardsonDivergenceLevel,
                                       RichardsonDivergenceLevel,
                                       RichardsonDivergenceLevel]) -> list[str]:
    """Apply preregistered scoped gates without changing the final PDE threshold."""
    if tuple(level.step for level in levels) != RICHARDSON_STEPS:
        raise ValueError("Richardson ladder drifted from the frozen protocol")
    failures: list[str] = []
    fine = levels[-1]
    fine_max = float(np.max(np.abs(fine.divergence)))
    fine_rms = scalar_rms(fine.divergence)
    if fine_max > DIVERGENCE_SAMPLED_MAX_GATE:
        failures.append("divergence_sampled_max")
    if fine_rms > DIVERGENCE_SAMPLED_RMS_GATE:
        failures.append("divergence_sampled_rms")
    if vector_rms(fine.velocity) < NONTRIVIAL_VELOCITY_RMS_FLOOR:
        failures.append("velocity_nontriviality")

    metrics = [divergence_metrics(level) for level in levels]
    for name in ("sampled_max", "sampled_rms"):
        coarse, medium, finest = (m[name] for m in metrics)
        medium_limit = max(RESOLUTION_STABILITY_FACTOR * coarse, RESOLUTION_FLOOR)
        fine_limit = max(RESOLUTION_STABILITY_FACTOR * medium, RESOLUTION_FLOOR)
        if medium > medium_limit:
            failures.append(f"divergence_{name}_coarse_to_medium_stability")
        if finest > fine_limit:
            failures.append(f"divergence_{name}_medium_to_fine_stability")
    return failures


def public_contract() -> dict[str, Any]:
    return {
        "schema": SCHEMA,
        "task": TASK,
        "agent2_candidate_pr": AGENT2_CANDIDATE_PR,
        "agent2_candidate_head": AGENT2_CANDIDATE_HEAD,
        "agent2_candidate_blob": AGENT2_CANDIDATE_BLOB,
        "agent1_reference_pr": AGENT1_REFERENCE_PR,
        "agent1_reference_head": AGENT1_REFERENCE_HEAD,
        "agent1_reference_blob": AGENT1_REFERENCE_BLOB,
        "richardson_steps": list(RICHARDSON_STEPS),
        "independent_reference": {
            "candidate_surface_consumed": "velocity(x,y,z,t) only",
            "centered_two_point_base_operator": True,
            "richardson_extrapolated_order": 4,
            "uses_candidate_internal_tensors": False,
            "uses_agent1_analytic_derivatives": False,
            "uses_agent1_fd6_or_laplacian": False,
            "uses_agent2_fd4_or_fd8_transport_verifier": False,
        },
        "frozen_gates": {
            "divergence_sampled_max": DIVERGENCE_SAMPLED_MAX_GATE,
            "divergence_sampled_rms": DIVERGENCE_SAMPLED_RMS_GATE,
            "resolution_stability_factor": RESOLUTION_STABILITY_FACTOR,
            "resolution_floor": RESOLUTION_FLOOR,
            "velocity_rms_floor": NONTRIVIAL_VELOCITY_RMS_FLOOR,
            "final_momentum_normalized_max_l2": FINAL_MOMENTUM_NORMALIZED_GATE,
            "final_divergence_max_l2": FINAL_DIVERGENCE_GATE,
        },
        "norm_scope": {
            "divergence_l2_is_heldout_sampled_rms": True,
            "whole_domain_volume_l2_assessed": False,
            "complete_ns_residual_assessed": False,
            "same_protocol_st006_comparison_valid": False,
        },
        "truth_boundary": {
            "strict_inner_candidate_artifact_independently_consumed": True,
            "strict_inner_divergence_independently_assessed": True,
            "candidate_manifest_integrity_independently_checked": True,
            "global_leading_velocity_materialized": False,
            "outer_join_materialized": False,
            "matched_pressure_included": False,
            "restricted_forcing_included": False,
            "correction_velocity_included": False,
            "complete_kokuno_candidate_assembled": False,
            "leading_only_ns_residual_assessed": False,
            "leading_plus_oscillatory_ns_residual_assessed": False,
            "after_correction_ns_residual_assessed": False,
            "heldout_ns_momentum_residual_assessed": False,
            "formal_full_domain_pde_gate_assessed": False,
            "whole_domain_boundary_support_gate_assessed": False,
            "pde_validated": False,
            "paper_exact": False,
        },
    }
