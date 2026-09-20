"""Independent black-box audit for the strict-inner Kokuno vorticity artifact.

The production quantity is Agent-2 PR #881's save/reload-bound public
``vorticity(x,y,z,t)`` surface.  The numerical reference uses only the reloaded
public ``velocity(x,y,z,t)`` callable and a centered Cartesian FD6 curl.

This module intentionally does not use Agent-1 analytic spatial derivatives,
Agent-2's production Jacobian/curl assembly, Agent-2's FD4 verifier, candidate
component tensors, pressure fitting, forcing fitting, or training internals.
It is a scoped strict-inner consistency audit and cannot promote PDE validity.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Callable

import numpy as np

TASK = "KOKUNO-A4-STRICT-INNER-VORTICITY-ARTIFACT-INDEPENDENT-AUDIT-079"
SCHEMA = "kokuno-a4-strict-inner-vorticity-artifact-independent-audit-v1"
AGENT2_VORTICITY_PR = 881
AGENT2_VORTICITY_HEAD = "ded8c88dd03b784368cf1a30601faf76a3682c64"
AGENT2_VORTICITY_BLOB = "024069383d9390e54e101529332bfd2707040659"
AGENT2_SPATIAL_PARENT_PR = 874
AGENT2_SPATIAL_PARENT_HEAD = "9cb3869b9cfd2d5b8dcb1da222df74d80e12d0c2"
AGENT1_REFERENCE_PR = 848
AGENT1_REFERENCE_HEAD = "300919075c07396e4e466f553be78eaae0707df2"
AGENT1_REFERENCE_BLOB = "5076ffa69748c2d7e68289a155b67288eab6dab4"

HELDOUT_SEED = 9173551
FD6_SPATIAL_STEPS = (3.0e-3, 1.5e-3, 7.5e-4)
FINE_RELATIVE_RMS_GATE = 5.0e-3
FINE_RELATIVE_SAMPLED_MAX_GATE = 1.0e-2
REFINEMENT_RATIO_GATE = 12.0
REFINEMENT_FLOOR = 2.0e-10
NONTRIVIAL_VORTICITY_RMS_FLOOR = 1.0e-10
NONTRIVIAL_VELOCITY_RMS_FLOOR = 1.0e-8
FINAL_MOMENTUM_NORMALIZED_GATE = 1.0e-3
FINAL_DIVERGENCE_GATE = 1.0e-5

VelocityProvider = Callable[[Any, Any, Any, Any], np.ndarray]


@dataclass(frozen=True)
class FD6CurlLevel:
    step: float
    vorticity: np.ndarray


def _validate_step(step: float) -> float:
    value = float(step)
    if not math.isfinite(value) or not 1.0e-6 <= value <= 5.0e-2:
        raise ValueError("spatial step must be finite and lie in [1e-6,5e-2]")
    return value


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


def fd6_axis_derivative(
    velocity: VelocityProvider,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    *,
    axis: int,
    step: float,
) -> np.ndarray:
    """Sixth-order centered derivative of public velocity along one axis."""
    if axis not in (0, 1, 2):
        raise ValueError("axis must be 0, 1 or 2")
    h = _validate_step(step)
    xb, yb, zb, tb = _broadcast_xyzt(x, y, z, t)
    base = [xb, yb, zb]
    coefficients = {
        -3: -1.0,
        -2: 9.0,
        -1: -45.0,
        1: 45.0,
        2: -9.0,
        3: 1.0,
    }
    total: np.ndarray | None = None
    for offset, coefficient in coefficients.items():
        shifted = [np.array(v, copy=True) for v in base]
        shifted[axis] = shifted[axis] + offset * h
        sample = _velocity_value(velocity, shifted[0], shifted[1], shifted[2], tb)
        contribution = coefficient * sample
        total = contribution if total is None else total + contribution
    assert total is not None
    out = total / (60.0 * h)
    if not np.all(np.isfinite(out)):
        raise RuntimeError("independent FD6 derivative became non-finite")
    return out


def fd6_jacobian(
    velocity: VelocityProvider,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    *,
    step: float,
) -> np.ndarray:
    derivatives = [
        fd6_axis_derivative(velocity, x, y, z, t, axis=axis, step=step)
        for axis in range(3)
    ]
    out = np.stack(derivatives, axis=-1)
    if out.shape[-2:] != (3, 3):
        raise RuntimeError("independent FD6 Jacobian has unexpected shape")
    return out


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


def fd6_vorticity(
    velocity: VelocityProvider,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    *,
    step: float,
) -> np.ndarray:
    return vorticity_from_jacobian(fd6_jacobian(velocity, x, y, z, t, step=step))


def fd6_level(
    velocity: VelocityProvider,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    *,
    step: float,
) -> FD6CurlLevel:
    h = _validate_step(step)
    return FD6CurlLevel(step=h, vorticity=fd6_vorticity(velocity, x, y, z, t, step=h))


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


def relative_rms_error(reference: Any, observed: Any) -> float:
    ref = np.asarray(reference, dtype=float)
    obs = np.asarray(observed, dtype=float)
    if ref.shape != obs.shape or ref.shape[-1:] != (3,):
        raise ValueError("reference and observed must share a (...,3) shape")
    return vector_rms(obs - ref) / max(vector_rms(ref), 1.0e-30)


def relative_sampled_max_error(reference: Any, observed: Any) -> float:
    ref = np.asarray(reference, dtype=float)
    obs = np.asarray(observed, dtype=float)
    if ref.shape != obs.shape or ref.shape[-1:] != (3,):
        raise ValueError("reference and observed must share a (...,3) shape")
    return sampled_vector_max(obs - ref) / max(sampled_vector_max(ref), 1.0e-30)


def refinement_metrics(
    coarse: FD6CurlLevel,
    medium: FD6CurlLevel,
    fine: FD6CurlLevel,
) -> tuple[float, float, float]:
    if not coarse.step > medium.step > fine.step:
        raise ValueError("FD6 levels must be ordered coarse > medium > fine")
    coarse_to_medium = vector_rms(coarse.vorticity - medium.vorticity)
    medium_to_fine = vector_rms(medium.vorticity - fine.vorticity)
    ratio = coarse_to_medium / max(medium_to_fine, 1.0e-300)
    return coarse_to_medium, medium_to_fine, ratio


def component_error_metrics(reference: Any, observed: Any) -> list[dict[str, float | str]]:
    ref = np.asarray(reference, dtype=float)
    obs = np.asarray(observed, dtype=float)
    if ref.shape != obs.shape or ref.shape[-1:] != (3,):
        raise ValueError("reference and observed must share a (...,3) shape")
    labels = ("omega_x", "omega_y", "omega_z")
    metrics: list[dict[str, float | str]] = []
    for component, label in enumerate(labels):
        delta = obs[..., component] - ref[..., component]
        base = ref[..., component]
        absolute_rms = float(np.sqrt(np.mean(delta * delta)))
        base_rms = float(np.sqrt(np.mean(base * base)))
        metrics.append(
            {
                "component": label,
                "absolute_rms": absolute_rms,
                "relative_rms": absolute_rms / max(base_rms, 1.0e-30),
                "absolute_sampled_max": float(np.max(np.abs(delta))),
            }
        )
    return metrics


def frozen_gate_failures(
    production_vorticity: Any,
    levels: tuple[FD6CurlLevel, FD6CurlLevel, FD6CurlLevel],
    *,
    velocity_rms: float,
) -> list[str]:
    if tuple(level.step for level in levels) != FD6_SPATIAL_STEPS:
        raise ValueError("FD6 spatial ladder drifted from frozen protocol")
    production = np.asarray(production_vorticity, dtype=float)
    if production.shape != levels[-1].vorticity.shape or production.shape[-1:] != (3,):
        raise ValueError("production vorticity has wrong shape")
    if not np.all(np.isfinite(production)):
        raise ValueError("production vorticity is non-finite")

    failures: list[str] = []
    fine = levels[-1].vorticity
    if relative_rms_error(fine, production) > FINE_RELATIVE_RMS_GATE:
        failures.append("fine_relative_rms")
    if relative_sampled_max_error(fine, production) > FINE_RELATIVE_SAMPLED_MAX_GATE:
        failures.append("fine_relative_sampled_max")
    if vector_rms(production) < NONTRIVIAL_VORTICITY_RMS_FLOOR:
        failures.append("vorticity_nontriviality")
    if float(velocity_rms) < NONTRIVIAL_VELOCITY_RMS_FLOOR:
        failures.append("velocity_nontriviality")

    _, medium_to_fine, ratio = refinement_metrics(*levels)
    if medium_to_fine > REFINEMENT_FLOOR and ratio < REFINEMENT_RATIO_GATE:
        failures.append("fd6_refinement_ratio")
    return failures


def public_contract() -> dict[str, Any]:
    return {
        "schema": SCHEMA,
        "task": TASK,
        "agent2_vorticity_pr": AGENT2_VORTICITY_PR,
        "agent2_vorticity_head": AGENT2_VORTICITY_HEAD,
        "agent2_vorticity_workflow_blob": AGENT2_VORTICITY_BLOB,
        "agent2_spatial_parent_pr": AGENT2_SPATIAL_PARENT_PR,
        "agent2_spatial_parent_head": AGENT2_SPATIAL_PARENT_HEAD,
        "agent1_reference_pr": AGENT1_REFERENCE_PR,
        "agent1_reference_head": AGENT1_REFERENCE_HEAD,
        "agent1_reference_blob": AGENT1_REFERENCE_BLOB,
        "heldout_seed": HELDOUT_SEED,
        "fd6_spatial_steps": list(FD6_SPATIAL_STEPS),
        "independent_reference": {
            "artifact_must_be_saved_and_reloaded": True,
            "production_surface": "reloaded public vorticity(x,y,z,t)",
            "numerical_reference_uses": "reloaded public velocity(x,y,z,t) only",
            "operator": "centered Cartesian FD6 curl",
            "agent1_analytic_spatial_derivatives_used": False,
            "agent2_production_jacobian_used": False,
            "agent2_fd4_verifier_used": False,
            "candidate_component_tensors_used": False,
        },
        "frozen_gates": {
            "fine_relative_rms": FINE_RELATIVE_RMS_GATE,
            "fine_relative_sampled_max": FINE_RELATIVE_SAMPLED_MAX_GATE,
            "refinement_ratio": REFINEMENT_RATIO_GATE,
            "refinement_floor": REFINEMENT_FLOOR,
            "vorticity_rms_floor": NONTRIVIAL_VORTICITY_RMS_FLOOR,
            "velocity_rms_floor": NONTRIVIAL_VELOCITY_RMS_FLOOR,
            "final_momentum_normalized": FINAL_MOMENTUM_NORMALIZED_GATE,
            "final_divergence": FINAL_DIVERGENCE_GATE,
        },
        "truth_boundary": {
            "strict_inner_vorticity_artifact_auditable": True,
            "whole_domain_vorticity_morphology_verified": False,
            "global_leading_velocity_materialized": False,
            "outer_join_materialized": False,
            "matched_pressure_included": False,
            "restricted_forcing_included": False,
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
