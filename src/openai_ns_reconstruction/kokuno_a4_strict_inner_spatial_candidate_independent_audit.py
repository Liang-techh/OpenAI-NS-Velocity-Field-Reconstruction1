"""Independent black-box spatial-Jacobian audit for the strict-inner Kokuno artifact.

The object under audit is Agent-2 PR #874's checksum-bound strict-inner spatial
candidate.  Agent 4 saves/reloads that artifact, treats public
``velocity_jacobian(x,y,z,t)`` as the production quantity, and reconstructs an
implementation-distinct reference from **only** reloaded public
``velocity(x,y,z,t)`` values.

The independent derivative is Richardson extrapolation of centered two-point
Cartesian derivatives.  It intentionally does not reuse Agent-1 analytic
Jacobians, Agent-2's fixed oscillatory FD6 Jacobian, Agent-2 #874's FD4
verifier, candidate component tensors, training tensors/losses, pressure fits,
or forcing fits.

This remains a scoped strict-inner audit.  It is not a complete Navier--Stokes
residual and cannot promote ``pde_validated``.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Callable

import numpy as np

TASK = "KOKUNO-A4-STRICT-INNER-SPATIAL-CANDIDATE-INDEPENDENT-AUDIT-078"
SCHEMA = "kokuno-a4-strict-inner-spatial-candidate-independent-audit-v1"
AGENT2_SPATIAL_CANDIDATE_PR = 874
AGENT2_SPATIAL_CANDIDATE_HEAD = "9cb3869b9cfd2d5b8dcb1da222df74d80e12d0c2"
AGENT2_SPATIAL_CANDIDATE_BLOB = "276c943cd1e769c36d1615366266bf3b2c633040"
AGENT2_DIFFERENTIABLE_PARENT_PR = 866
AGENT2_DIFFERENTIABLE_PARENT_HEAD = "68f84128f07b6743d368a9ab7a051e441f5a59b2"
AGENT1_REFERENCE_PR = 848
AGENT1_REFERENCE_HEAD = "300919075c07396e4e466f553be78eaae0707df2"
AGENT1_REFERENCE_BLOB = "5076ffa69748c2d7e68289a155b67288eab6dab4"

HELDOUT_SEED = 9173541
RICHARDSON_SPATIAL_STEPS = (2.4e-3, 1.2e-3, 6.0e-4)
FINE_RELATIVE_RMS_GATE = 3.0e-3
FINE_RELATIVE_SAMPLED_MAX_GATE = 8.0e-3
REFINEMENT_RATIO_GATE = 6.0
REFINEMENT_FLOOR = 5.0e-10
DIVERGENCE_SAMPLED_MAX_GATE = 1.0e-5
DIVERGENCE_SAMPLED_RMS_GATE = 1.0e-5
NONTRIVIAL_JACOBIAN_RMS_FLOOR = 1.0e-10
NONTRIVIAL_VELOCITY_RMS_FLOOR = 1.0e-8
FINAL_MOMENTUM_NORMALIZED_GATE = 1.0e-3
FINAL_DIVERGENCE_GATE = 1.0e-5

VelocityProvider = Callable[[Any, Any, Any, Any], np.ndarray]


@dataclass(frozen=True)
class RichardsonJacobianLevel:
    """One independent Cartesian-Jacobian resolution."""

    step: float
    jacobian: np.ndarray


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
        raise ValueError("spatial step must be finite and lie in [1e-6,5e-2]")
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


def centered_two_point_axis_derivative(
    velocity: VelocityProvider,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    *,
    axis: int,
    step: float,
) -> np.ndarray:
    """Second-order derivative of public velocity along one Cartesian axis."""
    if axis not in (0, 1, 2):
        raise ValueError("axis must be 0, 1 or 2")
    h = _validate_step(step)
    xb, yb, zb, tb = _broadcast_xyzt(x, y, z, t)
    base = [xb, yb, zb]
    plus = [np.array(v, copy=True) for v in base]
    minus = [np.array(v, copy=True) for v in base]
    plus[axis] = plus[axis] + h
    minus[axis] = minus[axis] - h
    vp = _velocity_value(velocity, plus[0], plus[1], plus[2], tb)
    vm = _velocity_value(velocity, minus[0], minus[1], minus[2], tb)
    return (vp - vm) / (2.0 * h)


def richardson_axis_derivative(
    velocity: VelocityProvider,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    *,
    axis: int,
    step: float,
) -> np.ndarray:
    """Fourth-order Richardson derivative from centered h and h/2 levels."""
    h = _validate_step(step)
    coarse = centered_two_point_axis_derivative(
        velocity, x, y, z, t, axis=axis, step=h
    )
    half = centered_two_point_axis_derivative(
        velocity, x, y, z, t, axis=axis, step=0.5 * h
    )
    out = (4.0 * half - coarse) / 3.0
    if not np.all(np.isfinite(out)):
        raise RuntimeError("independent Richardson spatial derivative became non-finite")
    return out


def richardson_jacobian(
    velocity: VelocityProvider,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    *,
    step: float,
) -> np.ndarray:
    """Return J[...,component,axis] from public velocity values only."""
    derivatives = [
        richardson_axis_derivative(velocity, x, y, z, t, axis=axis, step=step)
        for axis in range(3)
    ]
    out = np.stack(derivatives, axis=-1)
    if out.shape[-2:] != (3, 3):
        raise RuntimeError("independent Jacobian has an unexpected shape")
    return out


def richardson_jacobian_level(
    velocity: VelocityProvider,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    *,
    step: float,
) -> RichardsonJacobianLevel:
    h = _validate_step(step)
    return RichardsonJacobianLevel(
        step=h,
        jacobian=richardson_jacobian(velocity, x, y, z, t, step=h),
    )


def tensor_rms(value: Any) -> float:
    array = np.asarray(value, dtype=float)
    if array.shape[-2:] != (3, 3):
        raise ValueError("tensor_rms expects final axes (3,3)")
    return float(np.sqrt(np.mean(np.sum(array * array, axis=(-2, -1)))))


def vector_rms(value: Any) -> float:
    array = np.asarray(value, dtype=float)
    if array.shape[-1:] != (3,):
        raise ValueError("vector_rms expects a final Cartesian component axis of length 3")
    return float(np.sqrt(np.mean(np.sum(array * array, axis=-1))))


def relative_rms_error(reference: Any, observed: Any) -> float:
    ref = np.asarray(reference, dtype=float)
    obs = np.asarray(observed, dtype=float)
    if ref.shape != obs.shape or ref.shape[-2:] != (3, 3):
        raise ValueError("reference and observed must share a (...,3,3) shape")
    return tensor_rms(obs - ref) / max(tensor_rms(ref), 1.0e-30)


def relative_sampled_max_error(reference: Any, observed: Any) -> float:
    ref = np.asarray(reference, dtype=float)
    obs = np.asarray(observed, dtype=float)
    if ref.shape != obs.shape or ref.shape[-2:] != (3, 3):
        raise ValueError("reference and observed must share a (...,3,3) shape")
    diff_norm = np.sqrt(np.sum((obs - ref) ** 2, axis=(-2, -1)))
    ref_norm = np.sqrt(np.sum(ref * ref, axis=(-2, -1)))
    return float(np.max(diff_norm)) / max(float(np.max(ref_norm)), 1.0e-30)


def divergence_from_jacobian(jacobian: Any) -> np.ndarray:
    value = np.asarray(jacobian, dtype=float)
    if value.shape[-2:] != (3, 3):
        raise ValueError("jacobian must have final axes (3,3)")
    return np.trace(value, axis1=-2, axis2=-1)


def vorticity_from_jacobian(jacobian: Any) -> np.ndarray:
    value = np.asarray(jacobian, dtype=float)
    if value.shape[-2:] != (3, 3):
        raise ValueError("jacobian must have final axes (3,3)")
    return np.stack(
        (
            value[..., 2, 1] - value[..., 1, 2],
            value[..., 0, 2] - value[..., 2, 0],
            value[..., 1, 0] - value[..., 0, 1],
        ),
        axis=-1,
    )


def sampled_scalar_rms(value: Any) -> float:
    array = np.asarray(value, dtype=float)
    return float(np.sqrt(np.mean(array * array)))


def component_axis_error_metrics(reference: Any, observed: Any) -> list[dict[str, Any]]:
    ref = np.asarray(reference, dtype=float)
    obs = np.asarray(observed, dtype=float)
    if ref.shape != obs.shape or ref.shape[-2:] != (3, 3):
        raise ValueError("reference and observed must share a (...,3,3) shape")
    component_labels = ("u", "v", "w")
    axis_labels = ("x", "y", "z")
    out: list[dict[str, Any]] = []
    for component in range(3):
        for axis in range(3):
            delta = obs[..., component, axis] - ref[..., component, axis]
            base = ref[..., component, axis]
            absolute_rms = float(np.sqrt(np.mean(delta * delta)))
            base_rms = float(np.sqrt(np.mean(base * base)))
            out.append(
                {
                    "component": component_labels[component],
                    "axis": axis_labels[axis],
                    "absolute_rms": absolute_rms,
                    "relative_rms": absolute_rms / max(base_rms, 1.0e-30),
                    "absolute_sampled_max": float(np.max(np.abs(delta))),
                }
            )
    return out


def level_error_metrics(
    production_jacobian: Any,
    level: RichardsonJacobianLevel,
) -> dict[str, float]:
    production = np.asarray(production_jacobian, dtype=float)
    return {
        "step": float(level.step),
        "relative_rms": relative_rms_error(production, level.jacobian),
        "relative_sampled_max": relative_sampled_max_error(production, level.jacobian),
        "absolute_rms": tensor_rms(level.jacobian - production),
    }


def refinement_ratio(
    coarse: RichardsonJacobianLevel,
    medium: RichardsonJacobianLevel,
    fine: RichardsonJacobianLevel,
) -> tuple[float, float, float]:
    if not coarse.step > medium.step > fine.step:
        raise ValueError("Richardson levels must be ordered coarse > medium > fine")
    coarse_to_medium = tensor_rms(coarse.jacobian - medium.jacobian)
    medium_to_fine = tensor_rms(medium.jacobian - fine.jacobian)
    ratio = coarse_to_medium / max(medium_to_fine, 1.0e-300)
    return coarse_to_medium, medium_to_fine, ratio


def divergence_gate_failures(jacobian: Any) -> list[str]:
    divergence = divergence_from_jacobian(jacobian)
    sampled_max = float(np.max(np.abs(divergence)))
    sampled_rms = sampled_scalar_rms(divergence)
    failures: list[str] = []
    if sampled_max > DIVERGENCE_SAMPLED_MAX_GATE:
        failures.append("divergence_sampled_max")
    if sampled_rms > DIVERGENCE_SAMPLED_RMS_GATE:
        failures.append("divergence_sampled_rms")
    return failures


def frozen_gate_failures(
    production_jacobian: Any,
    levels: tuple[
        RichardsonJacobianLevel,
        RichardsonJacobianLevel,
        RichardsonJacobianLevel,
    ],
) -> list[str]:
    """Apply scoped frozen gates; these are not the final NS momentum gates."""
    if tuple(level.step for level in levels) != RICHARDSON_SPATIAL_STEPS:
        raise ValueError("Richardson spatial ladder drifted from the frozen protocol")
    production = np.asarray(production_jacobian, dtype=float)
    if production.shape != levels[-1].jacobian.shape or production.shape[-2:] != (3, 3):
        raise ValueError("production velocity_jacobian has the wrong shape")
    if not np.all(np.isfinite(production)):
        raise ValueError("production velocity_jacobian is non-finite")

    failures: list[str] = []
    fine = levels[-1]
    if relative_rms_error(production, fine.jacobian) > FINE_RELATIVE_RMS_GATE:
        failures.append("fine_relative_rms")
    if relative_sampled_max_error(production, fine.jacobian) > FINE_RELATIVE_SAMPLED_MAX_GATE:
        failures.append("fine_relative_sampled_max")
    if tensor_rms(production) < NONTRIVIAL_JACOBIAN_RMS_FLOOR:
        failures.append("jacobian_nontriviality")

    _, medium_to_fine, ratio = refinement_ratio(*levels)
    if medium_to_fine > REFINEMENT_FLOOR and ratio < REFINEMENT_RATIO_GATE:
        failures.append("richardson_refinement_ratio")
    failures.extend(divergence_gate_failures(fine.jacobian))
    return failures


def public_contract() -> dict[str, Any]:
    return {
        "schema": SCHEMA,
        "task": TASK,
        "agent2_spatial_candidate_pr": AGENT2_SPATIAL_CANDIDATE_PR,
        "agent2_spatial_candidate_head": AGENT2_SPATIAL_CANDIDATE_HEAD,
        "agent2_spatial_candidate_blob": AGENT2_SPATIAL_CANDIDATE_BLOB,
        "agent2_differentiable_parent_pr": AGENT2_DIFFERENTIABLE_PARENT_PR,
        "agent2_differentiable_parent_head": AGENT2_DIFFERENTIABLE_PARENT_HEAD,
        "agent1_reference_pr": AGENT1_REFERENCE_PR,
        "agent1_reference_head": AGENT1_REFERENCE_HEAD,
        "agent1_reference_blob": AGENT1_REFERENCE_BLOB,
        "heldout_seed": HELDOUT_SEED,
        "richardson_spatial_steps": list(RICHARDSON_SPATIAL_STEPS),
        "independent_reference": {
            "artifact_must_be_saved_and_reloaded": True,
            "candidate_surfaces_consumed": [
                "velocity(x,y,z,t)",
                "velocity_jacobian(x,y,z,t)",
            ],
            "numerical_reference_uses": "reloaded public velocity(x,y,z,t) only",
            "centered_two_point_base_operator": True,
            "richardson_extrapolated_order": 4,
            "uses_agent1_analytic_jacobian_as_reference": False,
            "uses_agent2_oscillatory_fd6_as_reference": False,
            "uses_agent2_fd4_verifier_as_reference": False,
            "uses_candidate_component_tensors": False,
        },
        "frozen_gates": {
            "fine_relative_rms": FINE_RELATIVE_RMS_GATE,
            "fine_relative_sampled_max": FINE_RELATIVE_SAMPLED_MAX_GATE,
            "refinement_ratio": REFINEMENT_RATIO_GATE,
            "refinement_floor": REFINEMENT_FLOOR,
            "jacobian_rms_floor": NONTRIVIAL_JACOBIAN_RMS_FLOOR,
            "velocity_rms_floor": NONTRIVIAL_VELOCITY_RMS_FLOOR,
            "independent_divergence_sampled_max": DIVERGENCE_SAMPLED_MAX_GATE,
            "independent_divergence_sampled_rms": DIVERGENCE_SAMPLED_RMS_GATE,
            "final_momentum_normalized_max_l2": FINAL_MOMENTUM_NORMALIZED_GATE,
            "final_divergence_max_l2": FINAL_DIVERGENCE_GATE,
        },
        "norm_scope": {
            "jacobian_consistency_is_pde_residual": False,
            "divergence_rms_is_heldout_sampled_rms": True,
            "whole_domain_volume_l2_assessed": False,
            "complete_ns_residual_assessed": False,
            "same_protocol_st006_comparison_valid": False,
        },
        "truth_boundary": {
            "strict_inner_spatial_candidate_artifact_independently_consumed": True,
            "strict_inner_velocity_jacobian_independently_assessed": True,
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
