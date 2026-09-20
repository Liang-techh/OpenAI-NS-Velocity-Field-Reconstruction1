"""Independent black-box transport audit for the strict-inner Kokuno artifact.

This Agent-4 validator consumes only a save/reloaded public
``velocity(x,y,z,t)`` callable from the checksum-bound Agent-2 strict-inner
candidate.  It reconstructs the pressure/forcing-free transport precursor

    T = u_t + (u . grad)u - nu * Laplacian(u)

with Richardson-extrapolated centered two-point derivatives.  The operator is
implementation-distinct from Agent-1 analytic derivatives, Agent-2 production
FD6/FD4 derivative paths, and the earlier A4 FD8 transport audit.

This is deliberately *not* a complete Navier-Stokes residual.  No pressure,
forcing, global/outer join, or Agent-3 correction velocity is invented here.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Callable

import numpy as np

TASK = "KOKUNO-A4-STRICT-INNER-TRANSPORT-ARTIFACT-INDEPENDENT-AUDIT-080"
SCHEMA = "kokuno-a4-strict-inner-transport-artifact-independent-audit-v1"
AGENT2_MORPHOLOGY_PR = 887
AGENT2_MORPHOLOGY_HEAD = "946b4d13f63ae6761f02e9e66d4c3c7ed485ffa6"
AGENT2_VORTICITY_PR = 881
AGENT2_VORTICITY_HEAD = "ded8c88dd03b784368cf1a30601faf76a3682c64"
AGENT1_REFERENCE_PR = 848
AGENT1_REFERENCE_HEAD = "300919075c07396e4e466f553be78eaae0707df2"

HELDOUT_SEED = 9173561
SPACE_STEPS = (4.0e-3, 2.0e-3, 1.0e-3)
TIME_STEPS = (3.2e-3, 1.6e-3, 8.0e-4)
NU = 0.01
MEDIUM_TO_FINE_RELATIVE_CHANGE_GATE = 5.0e-2
REFINEMENT_RATIO_GATE = 1.5
REFINEMENT_FLOOR = 5.0e-8
NONTRIVIAL_VELOCITY_RMS_FLOOR = 1.0e-8
FINAL_MOMENTUM_NORMALIZED_GATE = 1.0e-3
FINAL_DIVERGENCE_GATE = 1.0e-5

VelocityProvider = Callable[[Any, Any, Any, Any], np.ndarray]


@dataclass(frozen=True)
class TransportLevel:
    space_step: float
    time_step: float
    velocity: np.ndarray
    velocity_dt: np.ndarray
    jacobian: np.ndarray
    advection: np.ndarray
    laplacian: np.ndarray
    transport: np.ndarray
    divergence: np.ndarray
    normalized_rms: float
    normalized_sampled_max: float


def _validate_step(value: float, name: str) -> float:
    step = float(value)
    if not math.isfinite(step) or not 1.0e-6 <= step <= 5.0e-2:
        raise ValueError(f"{name} must be finite and lie in [1e-6,5e-2]")
    return step


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


def vector_rms(value: Any) -> float:
    array = np.asarray(value, dtype=float)
    if array.shape[-1:] != (3,):
        raise ValueError("vector_rms expects final axis length 3")
    return float(np.sqrt(np.mean(np.sum(array * array, axis=-1))))


def sampled_vector_max(value: Any) -> float:
    array = np.asarray(value, dtype=float)
    if array.shape[-1:] != (3,):
        raise ValueError("sampled_vector_max expects final axis length 3")
    return float(np.max(np.linalg.norm(array, axis=-1)))


def _shifted_velocity(
    velocity: VelocityProvider,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
    *,
    axis: int,
    offset: float,
) -> np.ndarray:
    if axis not in (0, 1, 2, 3):
        raise ValueError("axis must be 0,1,2,3 where 3 denotes time")
    coords = [np.array(x, copy=True), np.array(y, copy=True), np.array(z, copy=True), np.array(t, copy=True)]
    coords[axis] = coords[axis] + float(offset)
    return _velocity_value(velocity, coords[0], coords[1], coords[2], coords[3])


def _centered_first(
    velocity: VelocityProvider,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
    *,
    axis: int,
    step: float,
) -> np.ndarray:
    plus = _shifted_velocity(velocity, x, y, z, t, axis=axis, offset=step)
    minus = _shifted_velocity(velocity, x, y, z, t, axis=axis, offset=-step)
    return (plus - minus) / (2.0 * step)


def _centered_second(
    velocity: VelocityProvider,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
    *,
    axis: int,
    step: float,
) -> np.ndarray:
    if axis not in (0, 1, 2):
        raise ValueError("second derivative axis must be spatial")
    plus = _shifted_velocity(velocity, x, y, z, t, axis=axis, offset=step)
    center = _velocity_value(velocity, x, y, z, t)
    minus = _shifted_velocity(velocity, x, y, z, t, axis=axis, offset=-step)
    return (plus - 2.0 * center + minus) / (step * step)


def richardson_first(
    velocity: VelocityProvider,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    *,
    axis: int,
    step: float,
) -> np.ndarray:
    h = _validate_step(step, "first-derivative step")
    xb, yb, zb, tb = _broadcast_xyzt(x, y, z, t)
    coarse = _centered_first(velocity, xb, yb, zb, tb, axis=axis, step=h)
    fine = _centered_first(velocity, xb, yb, zb, tb, axis=axis, step=0.5 * h)
    out = (4.0 * fine - coarse) / 3.0
    if not np.all(np.isfinite(out)):
        raise RuntimeError("Richardson first derivative became non-finite")
    return out


def richardson_second(
    velocity: VelocityProvider,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    *,
    axis: int,
    step: float,
) -> np.ndarray:
    h = _validate_step(step, "second-derivative step")
    xb, yb, zb, tb = _broadcast_xyzt(x, y, z, t)
    coarse = _centered_second(velocity, xb, yb, zb, tb, axis=axis, step=h)
    fine = _centered_second(velocity, xb, yb, zb, tb, axis=axis, step=0.5 * h)
    out = (4.0 * fine - coarse) / 3.0
    if not np.all(np.isfinite(out)):
        raise RuntimeError("Richardson second derivative became non-finite")
    return out


def transport_level(
    velocity: VelocityProvider,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    *,
    space_step: float,
    time_step: float,
    nu: float = NU,
) -> TransportLevel:
    hs = _validate_step(space_step, "space_step")
    ht = _validate_step(time_step, "time_step")
    viscosity = float(nu)
    if not math.isfinite(viscosity) or viscosity < 0.0:
        raise ValueError("nu must be finite and nonnegative")
    xb, yb, zb, tb = _broadcast_xyzt(x, y, z, t)
    u = _velocity_value(velocity, xb, yb, zb, tb)
    ut = richardson_first(velocity, xb, yb, zb, tb, axis=3, step=ht)
    spatial = [richardson_first(velocity, xb, yb, zb, tb, axis=j, step=hs) for j in range(3)]
    jacobian = np.stack(spatial, axis=-1)
    advection = np.einsum("...j,...ij->...i", u, jacobian)
    laplacian = sum(
        (richardson_second(velocity, xb, yb, zb, tb, axis=j, step=hs) for j in range(3)),
        start=np.zeros_like(u),
    )
    transport = ut + advection - viscosity * laplacian
    divergence = np.trace(jacobian, axis1=-2, axis2=-1)
    rms_scale = vector_rms(ut) + vector_rms(advection) + viscosity * vector_rms(laplacian)
    max_scale = sampled_vector_max(ut) + sampled_vector_max(advection) + viscosity * sampled_vector_max(laplacian)
    return TransportLevel(
        space_step=hs,
        time_step=ht,
        velocity=u,
        velocity_dt=ut,
        jacobian=jacobian,
        advection=advection,
        laplacian=laplacian,
        transport=transport,
        divergence=divergence,
        normalized_rms=vector_rms(transport) / max(rms_scale, 1.0e-30),
        normalized_sampled_max=sampled_vector_max(transport) / max(max_scale, 1.0e-30),
    )


def three_resolution_levels(
    velocity: VelocityProvider,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
) -> tuple[TransportLevel, TransportLevel, TransportLevel]:
    return tuple(
        transport_level(velocity, x, y, z, t, space_step=hs, time_step=ht)
        for hs, ht in zip(SPACE_STEPS, TIME_STEPS)
    )  # type: ignore[return-value]


def resolution_metrics(
    levels: tuple[TransportLevel, TransportLevel, TransportLevel],
) -> tuple[float, float, float, float]:
    coarse, medium, fine = levels
    c2m = vector_rms(coarse.transport - medium.transport)
    m2f = vector_rms(medium.transport - fine.transport)
    relative_m2f = m2f / max(vector_rms(fine.transport), 1.0e-30)
    ratio = c2m / max(m2f, 1.0e-300)
    return c2m, m2f, relative_m2f, ratio


def gate_failures(levels: tuple[TransportLevel, TransportLevel, TransportLevel]) -> list[str]:
    if tuple(level.space_step for level in levels) != SPACE_STEPS:
        raise ValueError("space-step ladder drifted from frozen protocol")
    if tuple(level.time_step for level in levels) != TIME_STEPS:
        raise ValueError("time-step ladder drifted from frozen protocol")
    failures: list[str] = []
    fine = levels[-1]
    if vector_rms(fine.velocity) < NONTRIVIAL_VELOCITY_RMS_FLOOR:
        failures.append("velocity_nontriviality")
    div_max = float(np.max(np.abs(fine.divergence)))
    div_rms = float(np.sqrt(np.mean(fine.divergence * fine.divergence)))
    if div_max > FINAL_DIVERGENCE_GATE:
        failures.append("fine_divergence_sampled_max")
    if div_rms > FINAL_DIVERGENCE_GATE:
        failures.append("fine_divergence_sampled_rms")
    _, m2f, relative_m2f, ratio = resolution_metrics(levels)
    if relative_m2f > MEDIUM_TO_FINE_RELATIVE_CHANGE_GATE:
        failures.append("transport_medium_to_fine_relative_change")
    if m2f > REFINEMENT_FLOOR and ratio < REFINEMENT_RATIO_GATE:
        failures.append("transport_refinement_ratio")
    return failures


def public_contract() -> dict[str, Any]:
    return {
        "schema": SCHEMA,
        "task": TASK,
        "agent2_base_pr": AGENT2_MORPHOLOGY_PR,
        "agent2_base_head": AGENT2_MORPHOLOGY_HEAD,
        "agent2_candidate_pr": AGENT2_VORTICITY_PR,
        "agent2_candidate_head": AGENT2_VORTICITY_HEAD,
        "agent1_reference_pr": AGENT1_REFERENCE_PR,
        "agent1_reference_head": AGENT1_REFERENCE_HEAD,
        "heldout_seed": HELDOUT_SEED,
        "space_steps": list(SPACE_STEPS),
        "time_steps": list(TIME_STEPS),
        "nu": NU,
        "independent_reference": {
            "artifact_must_be_saved_and_reloaded": True,
            "numerical_source": "reloaded public velocity(x,y,z,t) only",
            "operator": "Richardson-extrapolated centered two-point first/second derivatives",
            "agent1_analytic_derivatives_used": False,
            "agent2_production_derivatives_used": False,
            "agent2_engineering_verifiers_used": False,
            "candidate_component_tensors_used": False,
            "pressure_fit_used": False,
            "forcing_fit_used": False,
        },
        "frozen_gates": {
            "medium_to_fine_relative_change": MEDIUM_TO_FINE_RELATIVE_CHANGE_GATE,
            "refinement_ratio": REFINEMENT_RATIO_GATE,
            "refinement_floor": REFINEMENT_FLOOR,
            "velocity_rms_floor": NONTRIVIAL_VELOCITY_RMS_FLOOR,
            "final_momentum_normalized": FINAL_MOMENTUM_NORMALIZED_GATE,
            "final_divergence": FINAL_DIVERGENCE_GATE,
        },
        "truth_boundary": {
            "strict_inner_transport_precursor_assessed": True,
            "transport_precursor_is_complete_ns_residual": False,
            "matched_pressure_included": False,
            "restricted_forcing_included": False,
            "global_leading_velocity_materialized": False,
            "outer_join_materialized": False,
            "agent3_correction_velocity_included": False,
            "complete_kokuno_candidate_assembled": False,
            "leading_only_ns_residual_assessed": False,
            "leading_plus_oscillatory_ns_residual_assessed": False,
            "after_correction_ns_residual_assessed": False,
            "heldout_ns_momentum_residual_assessed": False,
            "formal_full_domain_pde_gate_assessed": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
    }
