"""Independent strict-inner leading+oscillatory advection audit for Kokuno Agent 4.

This module is verifier-only.  It never calls Agent-1/A2 Jacobian,
self-advection, cross-advection, or directional-FD verifier paths to form the
numerical reference.  The reference consumes only a public total-velocity
callable and differentiates that callable in Cartesian coordinates with a
centered sixth-order stencil.

The audited production quantity lives on Agent-2 PR #830.  That quantity is
still restricted to the executable PA.10 inner contraction center plus the
frozen oscillatory field; it is not the missing global leading candidate and
is not a complete Navier--Stokes residual.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Callable

import numpy as np

TASK = "KOKUNO-A4-INNER-LEADING-OSCILLATORY-FULL-ADVECTION-INDEPENDENT-AUDIT-073"
SCHEMA = "kokuno-a4-inner-leading-oscillatory-full-advection-independent-audit-v1"
AGENT2_FULL_ADVECTION_PR = 830
AGENT2_FULL_ADVECTION_HEAD = "3479bed94ac6ef5c67547227f81a224ad43fb409"
AGENT2_FULL_ADVECTION_BLOB = "e44ccbe993d08951f48cfb5b093ee02ea5642dcf"
AGENT1_SELF_ADVECTION_PR = 829
AGENT1_SELF_ADVECTION_HEAD = "177190d770e0306e116674c5a2e1cbf977628760"
AGENT1_SELF_ADVECTION_BLOB = "b75f3b08ef3f1853550d8dfdef912980c15287c3"
FD6_STEPS = (1.6e-3, 8.0e-4, 4.0e-4)
FINE_RELATIVE_RMS_GATE = 2.0e-3
FINE_RELATIVE_MAX_GATE = 5.0e-3
REFINEMENT_RATIO_GATE = 8.0
REFINEMENT_FLOOR = 5.0e-10
NONTRIVIAL_ADVECTION_RMS_FLOOR = 1.0e-10
NONTRIVIAL_VELOCITY_RMS_FLOOR = 1.0e-8

VelocityProvider = Callable[[Any, Any, Any, Any], np.ndarray]


@dataclass(frozen=True)
class FD6AdvectionAuditLevel:
    step: float
    jacobian: np.ndarray
    self_advection: np.ndarray


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


def _validate_step(step: float) -> float:
    value = float(step)
    if not math.isfinite(value) or not 1.0e-6 <= value <= 5.0e-2:
        raise ValueError("step must be finite and lie in [1e-6,5e-2]")
    return value


def _velocity_value(
    provider: VelocityProvider,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
) -> np.ndarray:
    value = np.asarray(provider(x, y, z, t), dtype=float)
    expected = x.shape + (3,)
    if value.shape != expected:
        raise ValueError(f"velocity provider must return shape {expected}, got {value.shape}")
    if not np.all(np.isfinite(value)):
        raise ValueError("velocity provider returned non-finite values")
    return value


def cartesian_fd6_jacobian(
    velocity: VelocityProvider,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    *,
    step: float,
) -> np.ndarray:
    """Differentiate a public velocity callable with an independent Cartesian FD6 stencil.

    Returned convention is ``J[..., component, axis] = d_axis u_component``.
    """
    h = _validate_step(step)
    xb, yb, zb, tb = _broadcast_xyzt(x, y, z, t)
    coords = [xb, yb, zb]
    jac = np.empty(xb.shape + (3, 3), dtype=float)
    coefficients = {-3: -1.0, -2: 9.0, -1: -45.0, 1: 45.0, 2: -9.0, 3: 1.0}

    for axis in range(3):
        accum = np.zeros(xb.shape + (3,), dtype=float)
        for offset, coefficient in coefficients.items():
            shifted = [a.copy() for a in coords]
            shifted[axis] = shifted[axis] + offset * h
            accum += coefficient * _velocity_value(
                velocity, shifted[0], shifted[1], shifted[2], tb
            )
        jac[..., :, axis] = accum / (60.0 * h)

    if not np.all(np.isfinite(jac)):
        raise RuntimeError("FD6 Jacobian became non-finite")
    return jac


def cartesian_fd6_self_advection(
    velocity: VelocityProvider,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    *,
    step: float,
) -> FD6AdvectionAuditLevel:
    xb, yb, zb, tb = _broadcast_xyzt(x, y, z, t)
    base = _velocity_value(velocity, xb, yb, zb, tb)
    jac = cartesian_fd6_jacobian(velocity, xb, yb, zb, tb, step=step)
    adv = np.einsum("...j,...ij->...i", base, jac)
    if not np.all(np.isfinite(adv)):
        raise RuntimeError("FD6 self-advection became non-finite")
    return FD6AdvectionAuditLevel(step=float(step), jacobian=jac, self_advection=adv)


def vector_rms(value: Any) -> float:
    array = np.asarray(value, dtype=float)
    if array.shape[-1] != 3:
        raise ValueError("vector_rms expects a final Cartesian component axis of size 3")
    return float(np.sqrt(np.mean(np.sum(array * array, axis=-1))))


def relative_rms(reference: Any, value: Any) -> float:
    ref = np.asarray(reference, dtype=float)
    val = np.asarray(value, dtype=float)
    return vector_rms(val - ref) / max(vector_rms(ref), np.finfo(float).tiny)


def relative_sampled_max(reference: Any, value: Any) -> float:
    ref = np.asarray(reference, dtype=float)
    val = np.asarray(value, dtype=float)
    numerator = float(np.max(np.linalg.norm(val - ref, axis=-1)))
    denominator = max(float(np.max(np.linalg.norm(ref, axis=-1))), np.finfo(float).tiny)
    return numerator / denominator


def successive_refinement_ratio(coarse: Any, medium: Any, fine: Any) -> tuple[float, float, float]:
    coarse_medium = vector_rms(np.asarray(coarse) - np.asarray(medium))
    medium_fine = vector_rms(np.asarray(medium) - np.asarray(fine))
    ratio = coarse_medium / max(medium_fine, np.finfo(float).tiny)
    return ratio, coarse_medium, medium_fine


def frozen_gate_failures(
    production: Any,
    levels: tuple[FD6AdvectionAuditLevel, FD6AdvectionAuditLevel, FD6AdvectionAuditLevel],
    *,
    total_velocity: Any,
) -> list[str]:
    """Apply the pre-registered scoped A4 gates without altering thresholds."""
    if tuple(level.step for level in levels) != FD6_STEPS:
        raise ValueError("FD6 ladder drifted from the frozen protocol")
    reference = np.asarray(production, dtype=float)
    fine = levels[-1].self_advection
    ratio, _, medium_fine = successive_refinement_ratio(
        levels[0].self_advection, levels[1].self_advection, fine
    )
    failures: list[str] = []
    if relative_rms(reference, fine) > FINE_RELATIVE_RMS_GATE:
        failures.append("fine_relative_rms")
    if relative_sampled_max(reference, fine) > FINE_RELATIVE_MAX_GATE:
        failures.append("fine_relative_sampled_max")
    if ratio < REFINEMENT_RATIO_GATE and medium_fine > REFINEMENT_FLOOR:
        failures.append("fd6_refinement")
    if vector_rms(reference) < NONTRIVIAL_ADVECTION_RMS_FLOOR:
        failures.append("advection_nontriviality")
    if vector_rms(total_velocity) < NONTRIVIAL_VELOCITY_RMS_FLOOR:
        failures.append("velocity_nontriviality")
    return failures


def public_contract() -> dict[str, Any]:
    return {
        "schema": SCHEMA,
        "task": TASK,
        "agent2_full_advection_pr": AGENT2_FULL_ADVECTION_PR,
        "agent2_full_advection_head": AGENT2_FULL_ADVECTION_HEAD,
        "agent2_full_advection_blob": AGENT2_FULL_ADVECTION_BLOB,
        "agent1_self_advection_pr": AGENT1_SELF_ADVECTION_PR,
        "agent1_self_advection_head": AGENT1_SELF_ADVECTION_HEAD,
        "agent1_self_advection_blob": AGENT1_SELF_ADVECTION_BLOB,
        "fd6_steps": list(FD6_STEPS),
        "frozen_gates": {
            "fine_relative_rms": FINE_RELATIVE_RMS_GATE,
            "fine_relative_sampled_max": FINE_RELATIVE_MAX_GATE,
            "refinement_ratio": REFINEMENT_RATIO_GATE,
            "refinement_floor": REFINEMENT_FLOOR,
            "advection_rms_floor": NONTRIVIAL_ADVECTION_RMS_FLOOR,
            "velocity_rms_floor": NONTRIVIAL_VELOCITY_RMS_FLOOR,
            "final_momentum_normalized_max_l2": 1.0e-3,
            "final_divergence_max_l2": 1.0e-5,
        },
        "independent_reference": {
            "uses_public_total_velocity_only": True,
            "cartesian_fd6": True,
            "uses_agent1_self_advection": False,
            "uses_a2_cross_advection": False,
            "uses_a2_oscillatory_self_advection": False,
            "uses_a2_directional_fd4_verifier": False,
        },
        "truth_boundary": {
            "strict_inner_leading_plus_oscillatory_advection_assessed": True,
            "global_leading_velocity_materialized": False,
            "outer_join_materialized": False,
            "matched_pressure_included": False,
            "restricted_forcing_included": False,
            "viscous_laplacian_included": False,
            "correction_velocity_included": False,
            "complete_kokuno_candidate_assembled": False,
            "leading_only_ns_residual_assessed": False,
            "leading_plus_oscillatory_ns_residual_assessed": False,
            "after_correction_ns_residual_assessed": False,
            "heldout_ns_momentum_residual_assessed": False,
            "formal_full_domain_pde_gate_assessed": False,
            "pde_validated": False,
            "paper_exact": False,
        },
    }
