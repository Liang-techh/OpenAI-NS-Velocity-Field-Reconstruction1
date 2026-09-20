"""Independent black-box time-derivative audit for the strict-inner Kokuno artifact.

The object under audit is Agent-2 PR #866's first-class differentiable strict-inner
candidate.  Agent 4 saves/reloads that artifact and treats its public
``velocity(x,y,z,t)`` and ``velocity_dt(x,y,z,t)`` surfaces as the only candidate-facing
interfaces needed here.  The independent numerical reference differentiates only
``velocity`` at fixed Cartesian coordinates using Richardson extrapolation of a
centered two-point time derivative.  It therefore does not reuse Agent-1 analytic
time-chain rules or Agent-2's centered-FD4 verification path.

This remains a scoped strict-inner audit.  There is no outer/global leading join,
matched pressure, preregistered restricted forcing, Agent-3 correction velocity,
whole-domain support validation, or complete Navier--Stokes momentum residual.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Callable

import numpy as np

TASK = "KOKUNO-A4-STRICT-INNER-DIFFERENTIABLE-CANDIDATE-INDEPENDENT-AUDIT-077"
SCHEMA = "kokuno-a4-strict-inner-differentiable-candidate-independent-audit-v1"
AGENT2_DIFFERENTIABLE_CANDIDATE_PR = 866
AGENT2_DIFFERENTIABLE_CANDIDATE_HEAD = "68f84128f07b6743d368a9ab7a051e441f5a59b2"
AGENT2_DIFFERENTIABLE_CANDIDATE_BLOB = "8c09048a1d15dff321642e21ee699cce7d19831a"
AGENT2_BASE_CANDIDATE_PR = 857
AGENT2_BASE_CANDIDATE_HEAD = "a586b7afe4bb47dbfd4b25177586c79d362afff9"
AGENT1_REFERENCE_PR = 848
AGENT1_REFERENCE_HEAD = "300919075c07396e4e466f553be78eaae0707df2"
AGENT1_REFERENCE_BLOB = "5076ffa69748c2d7e68289a155b67288eab6dab4"

HELDOUT_SEED = 9173531
RICHARDSON_TIME_STEPS = (3.2e-3, 1.6e-3, 8.0e-4)
FINE_RELATIVE_RMS_GATE = 2.0e-3
FINE_RELATIVE_SAMPLED_MAX_GATE = 5.0e-3
REFINEMENT_RATIO_GATE = 6.0
REFINEMENT_FLOOR = 2.0e-10
NONTRIVIAL_VELOCITY_DT_RMS_FLOOR = 1.0e-10
FINAL_MOMENTUM_NORMALIZED_GATE = 1.0e-3
FINAL_DIVERGENCE_GATE = 1.0e-5

VelocityProvider = Callable[[Any, Any, Any, Any], np.ndarray]


@dataclass(frozen=True)
class RichardsonTimeDerivativeLevel:
    """One independent fixed-Cartesian time-derivative resolution."""

    step: float
    velocity_dt: np.ndarray


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
        raise ValueError("time step must be finite and lie in [1e-6,5e-2]")
    return value


def _vector_value(
    provider: VelocityProvider,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
    *,
    name: str,
) -> np.ndarray:
    value = np.asarray(provider(x, y, z, t), dtype=float)
    expected = x.shape + (3,)
    if value.shape != expected:
        raise ValueError(f"{name} provider must return shape {expected}, got {value.shape}")
    if not np.all(np.isfinite(value)):
        raise ValueError(f"{name} provider returned non-finite values")
    return value


def centered_two_point_time_derivative(
    velocity: VelocityProvider,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    *,
    step: float,
) -> np.ndarray:
    """Second-order centered derivative of public velocity at fixed Cartesian x,y,z."""
    h = _validate_step(step)
    xb, yb, zb, tb = _broadcast_xyzt(x, y, z, t)
    up = _vector_value(velocity, xb, yb, zb, tb + h, name="velocity")
    um = _vector_value(velocity, xb, yb, zb, tb - h, name="velocity")
    return (up - um) / (2.0 * h)


def richardson_time_derivative(
    velocity: VelocityProvider,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    *,
    step: float,
) -> np.ndarray:
    """Fourth-order Richardson derivative from two centered two-point levels."""
    h = _validate_step(step)
    coarse = centered_two_point_time_derivative(velocity, x, y, z, t, step=h)
    half = centered_two_point_time_derivative(velocity, x, y, z, t, step=0.5 * h)
    out = (4.0 * half - coarse) / 3.0
    if not np.all(np.isfinite(out)):
        raise RuntimeError("independent Richardson time derivative became non-finite")
    return out


def richardson_time_level(
    velocity: VelocityProvider,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    *,
    step: float,
) -> RichardsonTimeDerivativeLevel:
    h = _validate_step(step)
    return RichardsonTimeDerivativeLevel(
        step=h,
        velocity_dt=richardson_time_derivative(velocity, x, y, z, t, step=h),
    )


def vector_rms(value: Any) -> float:
    array = np.asarray(value, dtype=float)
    if array.shape[-1:] != (3,):
        raise ValueError("vector_rms expects a final Cartesian component axis of length 3")
    return float(np.sqrt(np.mean(np.sum(array * array, axis=-1))))


def relative_rms_error(reference: Any, observed: Any) -> float:
    ref = np.asarray(reference, dtype=float)
    obs = np.asarray(observed, dtype=float)
    if ref.shape != obs.shape or ref.shape[-1:] != (3,):
        raise ValueError("reference and observed must share a Cartesian-vector shape")
    denom = max(vector_rms(ref), 1.0e-30)
    return vector_rms(obs - ref) / denom


def relative_sampled_max_error(reference: Any, observed: Any) -> float:
    ref = np.asarray(reference, dtype=float)
    obs = np.asarray(observed, dtype=float)
    if ref.shape != obs.shape or ref.shape[-1:] != (3,):
        raise ValueError("reference and observed must share a Cartesian-vector shape")
    diff_norm = np.linalg.norm(obs - ref, axis=-1)
    ref_norm = np.linalg.norm(ref, axis=-1)
    denom = max(float(np.max(ref_norm)), 1.0e-30)
    return float(np.max(diff_norm)) / denom


def component_error_metrics(reference: Any, observed: Any) -> list[dict[str, float]]:
    ref = np.asarray(reference, dtype=float)
    obs = np.asarray(observed, dtype=float)
    if ref.shape != obs.shape or ref.shape[-1:] != (3,):
        raise ValueError("reference and observed must share a Cartesian-vector shape")
    labels = ("u_t", "v_t", "w_t")
    out: list[dict[str, float]] = []
    for component, label in enumerate(labels):
        delta = obs[..., component] - ref[..., component]
        ref_component = ref[..., component]
        rms = float(np.sqrt(np.mean(delta * delta)))
        ref_rms = float(np.sqrt(np.mean(ref_component * ref_component)))
        out.append(
            {
                "component": label,
                "absolute_rms": rms,
                "relative_rms": rms / max(ref_rms, 1.0e-30),
                "absolute_sampled_max": float(np.max(np.abs(delta))),
            }
        )
    return out


def level_error_metrics(
    production_velocity_dt: Any,
    level: RichardsonTimeDerivativeLevel,
) -> dict[str, float]:
    production = np.asarray(production_velocity_dt, dtype=float)
    return {
        "step": float(level.step),
        "relative_rms": relative_rms_error(production, level.velocity_dt),
        "relative_sampled_max": relative_sampled_max_error(production, level.velocity_dt),
        "absolute_rms": vector_rms(level.velocity_dt - production),
    }


def refinement_ratio(
    coarse: RichardsonTimeDerivativeLevel,
    medium: RichardsonTimeDerivativeLevel,
    fine: RichardsonTimeDerivativeLevel,
) -> tuple[float, float, float]:
    if not coarse.step > medium.step > fine.step:
        raise ValueError("Richardson levels must be ordered coarse > medium > fine")
    coarse_to_medium = vector_rms(coarse.velocity_dt - medium.velocity_dt)
    medium_to_fine = vector_rms(medium.velocity_dt - fine.velocity_dt)
    ratio = coarse_to_medium / max(medium_to_fine, 1.0e-300)
    return coarse_to_medium, medium_to_fine, ratio


def frozen_gate_failures(
    production_velocity_dt: Any,
    levels: tuple[
        RichardsonTimeDerivativeLevel,
        RichardsonTimeDerivativeLevel,
        RichardsonTimeDerivativeLevel,
    ],
) -> list[str]:
    """Apply frozen scoped consistency gates; these are not the final NS residual gates."""
    if tuple(level.step for level in levels) != RICHARDSON_TIME_STEPS:
        raise ValueError("Richardson time ladder drifted from the frozen protocol")
    production = np.asarray(production_velocity_dt, dtype=float)
    if production.shape != levels[-1].velocity_dt.shape or production.shape[-1:] != (3,):
        raise ValueError("production velocity_dt has the wrong shape")
    if not np.all(np.isfinite(production)):
        raise ValueError("production velocity_dt is non-finite")

    failures: list[str] = []
    fine = levels[-1]
    if relative_rms_error(production, fine.velocity_dt) > FINE_RELATIVE_RMS_GATE:
        failures.append("fine_relative_rms")
    if (
        relative_sampled_max_error(production, fine.velocity_dt)
        > FINE_RELATIVE_SAMPLED_MAX_GATE
    ):
        failures.append("fine_relative_sampled_max")
    if vector_rms(production) < NONTRIVIAL_VELOCITY_DT_RMS_FLOOR:
        failures.append("velocity_dt_nontriviality")

    _, medium_to_fine, ratio = refinement_ratio(*levels)
    if medium_to_fine > REFINEMENT_FLOOR and ratio < REFINEMENT_RATIO_GATE:
        failures.append("richardson_refinement_ratio")
    return failures


def public_contract() -> dict[str, Any]:
    return {
        "schema": SCHEMA,
        "task": TASK,
        "agent2_differentiable_candidate_pr": AGENT2_DIFFERENTIABLE_CANDIDATE_PR,
        "agent2_differentiable_candidate_head": AGENT2_DIFFERENTIABLE_CANDIDATE_HEAD,
        "agent2_differentiable_candidate_blob": AGENT2_DIFFERENTIABLE_CANDIDATE_BLOB,
        "agent2_base_candidate_pr": AGENT2_BASE_CANDIDATE_PR,
        "agent2_base_candidate_head": AGENT2_BASE_CANDIDATE_HEAD,
        "agent1_reference_pr": AGENT1_REFERENCE_PR,
        "agent1_reference_head": AGENT1_REFERENCE_HEAD,
        "agent1_reference_blob": AGENT1_REFERENCE_BLOB,
        "heldout_seed": HELDOUT_SEED,
        "richardson_time_steps": list(RICHARDSON_TIME_STEPS),
        "independent_reference": {
            "artifact_must_be_saved_and_reloaded": True,
            "candidate_surfaces_consumed": ["velocity(x,y,z,t)", "velocity_dt(x,y,z,t)"],
            "numerical_reference_uses": "reloaded public velocity(x,y,z,t) only",
            "centered_two_point_base_operator": True,
            "richardson_extrapolated_order": 4,
            "uses_agent1_analytic_time_derivative_as_reference": False,
            "uses_agent2_fd4_verifier_as_reference": False,
            "uses_candidate_component_tensors": False,
        },
        "frozen_gates": {
            "fine_relative_rms": FINE_RELATIVE_RMS_GATE,
            "fine_relative_sampled_max": FINE_RELATIVE_SAMPLED_MAX_GATE,
            "refinement_ratio": REFINEMENT_RATIO_GATE,
            "refinement_floor": REFINEMENT_FLOOR,
            "velocity_dt_rms_floor": NONTRIVIAL_VELOCITY_DT_RMS_FLOOR,
            "final_momentum_normalized_max_l2": FINAL_MOMENTUM_NORMALIZED_GATE,
            "final_divergence_max_l2": FINAL_DIVERGENCE_GATE,
        },
        "norm_scope": {
            "velocity_dt_consistency_is_pde_residual": False,
            "whole_domain_volume_l2_assessed": False,
            "complete_ns_residual_assessed": False,
            "same_protocol_st006_comparison_valid": False,
        },
        "truth_boundary": {
            "strict_inner_differentiable_candidate_artifact_independently_consumed": True,
            "strict_inner_velocity_dt_independently_assessed": True,
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
