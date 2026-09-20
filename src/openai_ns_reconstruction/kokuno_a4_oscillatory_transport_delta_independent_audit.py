"""Independent before/after strict-inner oscillatory transport audit for Kokuno Agent 4.

This verifier consumes public velocity values only. It reconstructs

    T_inner = d_t u_inner + (u_inner . grad)u_inner - nu Delta u_inner
    T_total = d_t u_total + (u_total . grad)u_total - nu Delta u_total
    delta T_osc = T_total - T_inner

with a centered eighth-order Cartesian/time operator, where
``u_total = u_inner + u_osc`` and repository viscosity is frozen at ``nu=0.01``.
It never calls Agent-1 production derivative methods, Agent-2 transport
subterms, or Agent-2's FD4 verifier. Pressure, forcing, outer/global joining,
and the correction field are intentionally absent, so this is not a complete
Navier--Stokes residual audit.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Callable

import numpy as np

TASK = "KOKUNO-A4-OSCILLATORY-TRANSPORT-DELTA-INDEPENDENT-AUDIT-075"
SCHEMA = "kokuno-a4-oscillatory-transport-delta-independent-audit-v1"
AGENT2_DELTA_PR = 849
AGENT2_DELTA_HEAD = "616e61abd2420a7c7aff3590a08393360039cd1a"
AGENT2_DELTA_BLOB = "0f4d1614945688a6f5d0d6f4192979cabbbe9da8"
AGENT2_TRANSPORT_PR = 840
AGENT2_TRANSPORT_HEAD = "016e3c152d7d11e4fedcd79b819df337f8942983"
AGENT2_TRANSPORT_BLOB = "d6082ebd346ba8150321bf3d62f104c107a70300"
AGENT1_LAPLACIAN_PR = 839
AGENT1_LAPLACIAN_HEAD = "e96ee90144976a992b62d76d4361a94eb16bd91e"
AGENT1_LAPLACIAN_BLOB = "74b0e185e8e0bbb08695d493e0a091c305a03006"
VISCOSITY = 0.01
FD8_STEPS = (2.0e-3, 1.0e-3, 5.0e-4)
INNER_FINE_RELATIVE_RMS_GATE = 1.5e-2
INNER_FINE_RELATIVE_MAX_GATE = 4.0e-2
DELTA_FINE_RELATIVE_RMS_GATE = 5.0e-2
DELTA_FINE_RELATIVE_MAX_GATE = 1.0e-1
INNER_TOTAL_REFINEMENT_RATIO_GATE = 8.0
DELTA_REFINEMENT_RATIO_GATE = 6.0
REFINEMENT_FLOOR = 2.0e-8
DIVERGENCE_SAMPLED_MAX_GATE = 1.0e-5
DIVERGENCE_SAMPLED_RMS_GATE = 1.0e-5
NONTRIVIAL_TRANSPORT_RMS_FLOOR = 1.0e-10
NONTRIVIAL_DELTA_RMS_FLOOR = 1.0e-10
NONTRIVIAL_VELOCITY_RMS_FLOOR = 1.0e-8
FINAL_MOMENTUM_NORMALIZED_GATE = 1.0e-3
FINAL_DIVERGENCE_GATE = 1.0e-5

VelocityProvider = Callable[[Any, Any, Any, Any], np.ndarray]

_FIRST_WEIGHTS = {
    -4: 1.0 / 280.0,
    -3: -4.0 / 105.0,
    -2: 1.0 / 5.0,
    -1: -4.0 / 5.0,
    1: 4.0 / 5.0,
    2: -1.0 / 5.0,
    3: 4.0 / 105.0,
    4: -1.0 / 280.0,
}
_SECOND_WEIGHTS = {
    -4: -1.0 / 560.0,
    -3: 8.0 / 315.0,
    -2: -1.0 / 5.0,
    -1: 8.0 / 5.0,
    0: -205.0 / 72.0,
    1: 8.0 / 5.0,
    2: -1.0 / 5.0,
    3: 8.0 / 315.0,
    4: -1.0 / 560.0,
}


@dataclass(frozen=True)
class FD8TransportLevel:
    step: float
    velocity: np.ndarray
    velocity_dt: np.ndarray
    jacobian: np.ndarray
    divergence: np.ndarray
    self_advection: np.ndarray
    laplacian: np.ndarray
    transport: np.ndarray


@dataclass(frozen=True)
class BeforeAfterFD8Level:
    step: float
    inner: FD8TransportLevel
    total: FD8TransportLevel
    delta_transport: np.ndarray


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


def _first_derivative(provider: VelocityProvider, x: np.ndarray, y: np.ndarray,
                      z: np.ndarray, t: np.ndarray, *, step: float,
                      axis: int) -> np.ndarray:
    h = _validate_step(step)
    if axis not in (0, 1, 2, 3):
        raise ValueError("axis must be 0,1,2,3 for x,y,z,t")
    base = [x, y, z, t]
    accum = np.zeros(x.shape + (3,), dtype=float)
    for offset, coefficient in _FIRST_WEIGHTS.items():
        shifted = [a.copy() for a in base]
        shifted[axis] = shifted[axis] + offset * h
        accum += coefficient * _velocity_value(
            provider, shifted[0], shifted[1], shifted[2], shifted[3]
        )
    return accum / h


def _second_spatial_derivative(provider: VelocityProvider, x: np.ndarray,
                               y: np.ndarray, z: np.ndarray, t: np.ndarray,
                               *, step: float, axis: int,
                               base_value: np.ndarray) -> np.ndarray:
    h = _validate_step(step)
    if axis not in (0, 1, 2):
        raise ValueError("spatial axis must be 0,1,2")
    coords = [x, y, z]
    accum = _SECOND_WEIGHTS[0] * base_value
    for offset, coefficient in _SECOND_WEIGHTS.items():
        if offset == 0:
            continue
        shifted = [a.copy() for a in coords]
        shifted[axis] = shifted[axis] + offset * h
        accum += coefficient * _velocity_value(
            provider, shifted[0], shifted[1], shifted[2], t
        )
    return accum / (h * h)


def cartesian_fd8_transport(velocity: VelocityProvider, x: Any, y: Any,
                            z: Any, t: Any, *, step: float) -> FD8TransportLevel:
    """Reconstruct ``u_t + (u.grad)u - nu Delta u`` from velocity values only."""
    h = _validate_step(step)
    xb, yb, zb, tb = _broadcast_xyzt(x, y, z, t)
    base = _velocity_value(velocity, xb, yb, zb, tb)
    velocity_dt = _first_derivative(velocity, xb, yb, zb, tb, step=h, axis=3)
    jacobian = np.empty(xb.shape + (3, 3), dtype=float)
    laplacian = np.zeros_like(base)
    for axis in range(3):
        jacobian[..., :, axis] = _first_derivative(
            velocity, xb, yb, zb, tb, step=h, axis=axis
        )
        laplacian += _second_spatial_derivative(
            velocity, xb, yb, zb, tb, step=h, axis=axis, base_value=base
        )
    divergence = np.trace(jacobian, axis1=-2, axis2=-1)
    self_advection = np.einsum("...j,...ij->...i", base, jacobian)
    transport = velocity_dt + self_advection - VISCOSITY * laplacian
    for name, value in (("velocity_dt", velocity_dt), ("jacobian", jacobian),
                        ("divergence", divergence),
                        ("self_advection", self_advection),
                        ("laplacian", laplacian), ("transport", transport)):
        if not np.all(np.isfinite(value)):
            raise RuntimeError(f"FD8 {name} became non-finite")
    return FD8TransportLevel(h, base, velocity_dt, jacobian, divergence,
                             self_advection, laplacian, transport)


def before_after_fd8_transport(inner_velocity: VelocityProvider,
                               total_velocity: VelocityProvider,
                               x: Any, y: Any, z: Any, t: Any, *,
                               step: float) -> BeforeAfterFD8Level:
    """Independently reconstruct strict-inner transport before and after oscillation."""
    inner = cartesian_fd8_transport(inner_velocity, x, y, z, t, step=step)
    total = cartesian_fd8_transport(total_velocity, x, y, z, t, step=step)
    if inner.step != total.step:
        raise RuntimeError("inner/total FD8 steps drifted")
    delta = total.transport - inner.transport
    if not np.all(np.isfinite(delta)):
        raise RuntimeError("FD8 oscillatory transport delta became non-finite")
    return BeforeAfterFD8Level(inner.step, inner, total, delta)


def vector_rms(value: Any) -> float:
    array = np.asarray(value, dtype=float)
    if array.shape[-1] != 3:
        raise ValueError("vector_rms expects a final Cartesian component axis of size 3")
    return float(np.sqrt(np.mean(np.sum(array * array, axis=-1))))


def scalar_rms(value: Any) -> float:
    array = np.asarray(value, dtype=float)
    return float(np.sqrt(np.mean(array * array)))


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


def successive_refinement_ratio(coarse: Any, medium: Any,
                                fine: Any) -> tuple[float, float, float]:
    coarse_medium = vector_rms(np.asarray(coarse) - np.asarray(medium))
    medium_fine = vector_rms(np.asarray(medium) - np.asarray(fine))
    ratio = coarse_medium / max(medium_fine, np.finfo(float).tiny)
    return ratio, coarse_medium, medium_fine


def transports_at_viscosity(level: BeforeAfterFD8Level,
                            viscosity: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    value = float(viscosity)
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError("viscosity must be positive and finite")
    inner = level.inner.velocity_dt + level.inner.self_advection - value * level.inner.laplacian
    total = level.total.velocity_dt + level.total.self_advection - value * level.total.laplacian
    return inner, total, total - inner


def frozen_gate_failures(production_inner_transport: Any,
                         production_delta_transport: Any,
                         levels: tuple[BeforeAfterFD8Level, BeforeAfterFD8Level,
                                       BeforeAfterFD8Level], *,
                         inner_velocity: Any, total_velocity: Any) -> list[str]:
    """Apply preregistered scoped gates without changing the final PDE threshold."""
    if tuple(level.step for level in levels) != FD8_STEPS:
        raise ValueError("FD8 ladder drifted from the frozen protocol")
    production_inner = np.asarray(production_inner_transport, dtype=float)
    production_delta = np.asarray(production_delta_transport, dtype=float)
    fine = levels[-1]
    if production_inner.shape != fine.inner.transport.shape:
        raise ValueError("production inner transport shape mismatch")
    if production_delta.shape != fine.delta_transport.shape:
        raise ValueError("production delta transport shape mismatch")

    failures: list[str] = []
    if relative_rms(production_inner, fine.inner.transport) > INNER_FINE_RELATIVE_RMS_GATE:
        failures.append("inner_fine_relative_rms")
    if relative_sampled_max(production_inner, fine.inner.transport) > INNER_FINE_RELATIVE_MAX_GATE:
        failures.append("inner_fine_relative_sampled_max")
    if relative_rms(production_delta, fine.delta_transport) > DELTA_FINE_RELATIVE_RMS_GATE:
        failures.append("delta_fine_relative_rms")
    if relative_sampled_max(production_delta, fine.delta_transport) > DELTA_FINE_RELATIVE_MAX_GATE:
        failures.append("delta_fine_relative_sampled_max")

    for name, values, gate in (
        ("inner", tuple(level.inner.transport for level in levels), INNER_TOTAL_REFINEMENT_RATIO_GATE),
        ("total", tuple(level.total.transport for level in levels), INNER_TOTAL_REFINEMENT_RATIO_GATE),
        ("delta", tuple(level.delta_transport for level in levels), DELTA_REFINEMENT_RATIO_GATE),
    ):
        ratio, _, medium_fine = successive_refinement_ratio(*values)
        if ratio < gate and medium_fine > REFINEMENT_FLOOR:
            failures.append(f"{name}_fd8_refinement")

    for name, divergence in (("inner", fine.inner.divergence),
                             ("total", fine.total.divergence)):
        if float(np.max(np.abs(divergence))) > DIVERGENCE_SAMPLED_MAX_GATE:
            failures.append(f"{name}_divergence_sampled_max")
        if scalar_rms(divergence) > DIVERGENCE_SAMPLED_RMS_GATE:
            failures.append(f"{name}_divergence_sampled_rms")

    if vector_rms(production_inner) < NONTRIVIAL_TRANSPORT_RMS_FLOOR:
        failures.append("inner_transport_nontriviality")
    if vector_rms(production_delta) < NONTRIVIAL_DELTA_RMS_FLOOR:
        failures.append("delta_transport_nontriviality")
    if vector_rms(inner_velocity) < NONTRIVIAL_VELOCITY_RMS_FLOOR:
        failures.append("inner_velocity_nontriviality")
    if vector_rms(total_velocity) < NONTRIVIAL_VELOCITY_RMS_FLOOR:
        failures.append("total_velocity_nontriviality")
    return failures


def public_contract() -> dict[str, Any]:
    return {
        "schema": SCHEMA, "task": TASK,
        "agent2_delta_pr": AGENT2_DELTA_PR,
        "agent2_delta_head": AGENT2_DELTA_HEAD,
        "agent2_delta_blob": AGENT2_DELTA_BLOB,
        "agent2_transport_pr": AGENT2_TRANSPORT_PR,
        "agent2_transport_head": AGENT2_TRANSPORT_HEAD,
        "agent2_transport_blob": AGENT2_TRANSPORT_BLOB,
        "agent1_laplacian_pr": AGENT1_LAPLACIAN_PR,
        "agent1_laplacian_head": AGENT1_LAPLACIAN_HEAD,
        "agent1_laplacian_blob": AGENT1_LAPLACIAN_BLOB,
        "viscosity": VISCOSITY, "fd8_steps": list(FD8_STEPS),
        "frozen_gates": {
            "inner_fine_relative_rms": INNER_FINE_RELATIVE_RMS_GATE,
            "inner_fine_relative_sampled_max": INNER_FINE_RELATIVE_MAX_GATE,
            "delta_fine_relative_rms": DELTA_FINE_RELATIVE_RMS_GATE,
            "delta_fine_relative_sampled_max": DELTA_FINE_RELATIVE_MAX_GATE,
            "inner_total_refinement_ratio": INNER_TOTAL_REFINEMENT_RATIO_GATE,
            "delta_refinement_ratio": DELTA_REFINEMENT_RATIO_GATE,
            "refinement_floor": REFINEMENT_FLOOR,
            "divergence_sampled_max": DIVERGENCE_SAMPLED_MAX_GATE,
            "divergence_sampled_rms": DIVERGENCE_SAMPLED_RMS_GATE,
            "inner_transport_rms_floor": NONTRIVIAL_TRANSPORT_RMS_FLOOR,
            "delta_transport_rms_floor": NONTRIVIAL_DELTA_RMS_FLOOR,
            "velocity_rms_floor": NONTRIVIAL_VELOCITY_RMS_FLOOR,
            "final_momentum_normalized_max_l2": FINAL_MOMENTUM_NORMALIZED_GATE,
            "final_divergence_max_l2": FINAL_DIVERGENCE_GATE,
        },
        "independent_reference": {
            "uses_public_inner_velocity_only": True,
            "uses_public_total_velocity_only": True,
            "centered_fd8_time": True,
            "cartesian_fd8_jacobian": True,
            "cartesian_fd8_laplacian": True,
            "cartesian_contraction_advection": True,
            "uses_agent1_production_derivatives": False,
            "uses_agent2_production_subterms": False,
            "uses_agent2_fd4_verifier": False,
        },
        "norm_scope": {
            "before_after_metrics_are_transport_consistency": True,
            "oscillation_rms_ratio_is_observation_only": True,
            "divergence_l2_is_heldout_sampled_rms": True,
            "whole_domain_volume_l2_assessed": False,
            "complete_ns_residual_assessed": False,
            "same_protocol_st006_comparison_valid": False,
        },
        "truth_boundary": {
            "strict_inner_before_after_transport_independently_assessed": True,
            "oscillatory_transport_delta_independently_assessed": True,
            "strict_inner_divergence_independently_assessed": True,
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
            "boundary_support_gate_assessed": False,
            "pde_validated": False,
            "paper_exact": False,
        },
    }
