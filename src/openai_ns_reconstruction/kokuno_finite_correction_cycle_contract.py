"""Finite correction-cycle audit contract for Kokuno Agent 3.

This module advances the typed same-cycle momentum-defect seam into a finite-step
transition contract.  It checks the project iteration

    u^(k+1) = u^k + delta_u_k

(and the corresponding time-derivative identity) on raw providers, then
recomputes the Navier--Stokes momentum defect before and after the step through
``kokuno_same_cycle_defect_contract`` on disjoint held-in / held-out point sets.

No residual, defect, stress, correction target, gain, or scientific threshold is
accepted from the caller.  The finite-step gain guard is deliberately simple:
both held-in and held-out raw defect RMS must strictly decrease and their sampled
maxima must not increase (up to one floating-point ULP).  This is a repository
engineering rejection guard for a divergent finite step, *not* Kokuno's formal
correction-cycle gain theorem and not the final normalized 1e-3 PDE gate.

The current deterministic receipt is an analytic mechanics regression only.  It
does not claim that the missing Agent-1 global leading/pressure/forcing handoff
exists or that a real Kokuno composite correction cycle has run.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Callable, Sequence

import numpy as np

from .kokuno_same_cycle_defect_contract import (
    CycleIdentity,
    RestrictedForcingProvider,
    ScalarFieldProvider,
    VectorFieldProvider,
    evaluate_disjoint_heldin_heldout,
)

VectorEvaluator = Callable[[float, float, float, float], Sequence[float]]

# Numerical identity tolerance only.  This is not a PDE/gain threshold.
UPDATE_CLOSURE_ATOL = 2.0e-11


@dataclass(frozen=True)
class CorrectionFieldProvider:
    """One velocity correction linking two consecutive cycle identities."""

    from_identity: CycleIdentity
    to_identity: CycleIdentity
    source_ref: str
    velocity_evaluator: VectorEvaluator
    velocity_dt_evaluator: VectorEvaluator

    def __post_init__(self) -> None:
        if not self.source_ref.strip():
            raise ValueError("source_ref must be non-empty")
        if self.to_identity.cycle_id != self.from_identity.cycle_id:
            raise ValueError("correction identities must share one cycle_id")
        if self.to_identity.cycle_index != self.from_identity.cycle_index + 1:
            raise ValueError("correction must advance exactly one cycle index")
        if self.to_identity.state_token == self.from_identity.state_token:
            raise ValueError("correction must change the cycle state token")


@dataclass(frozen=True)
class FiniteCorrectionCycleReport:
    from_identity: CycleIdentity
    to_identity: CycleIdentity
    correction_source_ref: str
    update_velocity_closure_max: float
    update_velocity_dt_closure_max: float
    correction_vector_rms: float
    correction_vector_max: float
    correction_dt_vector_rms: float
    correction_divergence_rms: float
    correction_divergence_max: float
    held_in_before_rms: float
    held_in_after_rms: float
    held_in_before_max: float
    held_in_after_max: float
    held_out_before_rms: float
    held_out_after_rms: float
    held_out_before_max: float
    held_out_after_max: float
    held_in_rms_ratio: float
    held_in_max_ratio: float
    held_out_rms_ratio: float
    held_out_max_ratio: float
    held_in_divergence_before_rms: float
    held_in_divergence_after_rms: float
    held_in_divergence_before_max: float
    held_in_divergence_after_max: float
    held_out_divergence_before_rms: float
    held_out_divergence_after_rms: float
    held_out_divergence_before_max: float
    held_out_divergence_after_max: float
    correction_nontrivial: bool
    finite_step_gain_guard_passed: bool

    def to_receipt(self) -> dict[str, object]:
        payload = asdict(self)
        payload["from_identity"] = asdict(self.from_identity)
        payload["to_identity"] = asdict(self.to_identity)
        return payload


def _as_points(
    points: Sequence[Sequence[float]] | np.ndarray, *, label: str
) -> np.ndarray:
    arr = np.asarray(points, dtype=float)
    if arr.ndim != 2 or arr.shape[1] != 4 or arr.shape[0] == 0:
        raise ValueError(f"{label} must have shape (N, 4) with N >= 1")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{label} must contain only finite values")
    rows = [tuple(float(v) for v in row) for row in arr]
    if len(set(rows)) != len(rows):
        raise ValueError(f"{label} contains duplicate points")
    return arr


def _eval_vector(
    evaluator: VectorEvaluator, point: np.ndarray, *, source_ref: str
) -> np.ndarray:
    value = np.asarray(evaluator(*map(float, point)), dtype=float)
    if value.shape != (3,) or not np.all(np.isfinite(value)):
        raise ValueError(f"{source_ref} must return a finite 3-vector")
    return value


def _bundle_identity(
    velocity: VectorFieldProvider,
    velocity_dt: VectorFieldProvider,
    pressure: ScalarFieldProvider,
    forcing: RestrictedForcingProvider,
    *,
    label: str,
) -> CycleIdentity:
    identity = velocity.identity
    for name, other in (
        ("velocity_dt", velocity_dt.identity),
        ("pressure", pressure.identity),
        ("restricted_forcing", forcing.identity),
    ):
        if other != identity:
            raise ValueError(f"{label} {name} does not share one cycle identity")
    return identity


def _assert_transition_identity(
    before_velocity: VectorFieldProvider,
    before_velocity_dt: VectorFieldProvider,
    before_pressure: ScalarFieldProvider,
    before_forcing: RestrictedForcingProvider,
    after_velocity: VectorFieldProvider,
    after_velocity_dt: VectorFieldProvider,
    after_pressure: ScalarFieldProvider,
    after_forcing: RestrictedForcingProvider,
    correction: CorrectionFieldProvider,
) -> tuple[CycleIdentity, CycleIdentity]:
    before = _bundle_identity(
        before_velocity,
        before_velocity_dt,
        before_pressure,
        before_forcing,
        label="before",
    )
    after = _bundle_identity(
        after_velocity,
        after_velocity_dt,
        after_pressure,
        after_forcing,
        label="after",
    )
    if before.cycle_id != after.cycle_id:
        raise ValueError("before/after states must share one cycle_id")
    if after.cycle_index != before.cycle_index + 1:
        raise ValueError("finite correction transition must advance one cycle index")
    if after.state_token == before.state_token:
        raise ValueError("finite correction transition must change state_token")
    if correction.from_identity != before or correction.to_identity != after:
        raise ValueError("correction provider does not bind the exact transition identities")
    if before_forcing.source_ref != after_forcing.source_ref:
        raise ValueError("restricted forcing source_ref changed across correction step")
    if before_forcing.restriction_receipt != after_forcing.restriction_receipt:
        raise ValueError("restricted forcing receipt changed across correction step")
    return before, after


def _shift(point: np.ndarray, axis: int, amount: float) -> np.ndarray:
    shifted = point.copy()
    shifted[axis] += amount
    return shifted


def _fd4_component_derivative(
    provider: VectorFieldProvider,
    point: np.ndarray,
    axis: int,
    component: int,
    h: float,
) -> float:
    def component_at(p: np.ndarray) -> float:
        return float(
            _eval_vector(provider.evaluator, p, source_ref=provider.source_ref)[component]
        )

    return (
        -component_at(_shift(point, axis, 2.0 * h))
        + 8.0 * component_at(_shift(point, axis, h))
        - 8.0 * component_at(_shift(point, axis, -h))
        + component_at(_shift(point, axis, -2.0 * h))
    ) / (12.0 * h)


def _divergence_values(
    provider: VectorFieldProvider, points: np.ndarray, h: float
) -> np.ndarray:
    values = []
    for point in points:
        values.append(
            sum(
                _fd4_component_derivative(provider, point, axis, axis, h)
                for axis in range(3)
            )
        )
    return np.asarray(values, dtype=float)


def _divergence_stats(
    provider: VectorFieldProvider, points: np.ndarray, h: float
) -> tuple[float, float]:
    values = _divergence_values(provider, points, h)
    return float(np.sqrt(np.mean(values * values))), float(np.max(np.abs(values)))


def _correction_divergence_stats(
    correction: CorrectionFieldProvider, points: np.ndarray, h: float
) -> tuple[float, float]:
    provider = VectorFieldProvider(
        correction.to_identity,
        correction.source_ref + ":delta_u",
        correction.velocity_evaluator,
    )
    return _divergence_stats(provider, points, h)


def _ratio(after: float, before: float) -> float:
    if before == 0.0:
        return 1.0 if after == 0.0 else float("inf")
    return float(after / before)


def _not_increased(after: float, before: float) -> bool:
    return bool(after <= np.nextafter(float(before), np.inf))


def _update_closure(
    before_velocity: VectorFieldProvider,
    before_velocity_dt: VectorFieldProvider,
    after_velocity: VectorFieldProvider,
    after_velocity_dt: VectorFieldProvider,
    correction: CorrectionFieldProvider,
    points: np.ndarray,
) -> tuple[float, float, float, float]:
    velocity_errors: list[float] = []
    velocity_dt_errors: list[float] = []
    correction_norms: list[float] = []
    correction_dt_norms: list[float] = []
    for point in points:
        before_u = _eval_vector(
            before_velocity.evaluator, point, source_ref=before_velocity.source_ref
        )
        after_u = _eval_vector(
            after_velocity.evaluator, point, source_ref=after_velocity.source_ref
        )
        delta_u = _eval_vector(
            correction.velocity_evaluator, point, source_ref=correction.source_ref
        )
        before_ut = _eval_vector(
            before_velocity_dt.evaluator,
            point,
            source_ref=before_velocity_dt.source_ref,
        )
        after_ut = _eval_vector(
            after_velocity_dt.evaluator,
            point,
            source_ref=after_velocity_dt.source_ref,
        )
        delta_ut = _eval_vector(
            correction.velocity_dt_evaluator,
            point,
            source_ref=correction.source_ref + ":dt",
        )
        velocity_errors.append(float(np.max(np.abs(after_u - before_u - delta_u))))
        velocity_dt_errors.append(
            float(np.max(np.abs(after_ut - before_ut - delta_ut)))
        )
        correction_norms.append(float(np.linalg.norm(delta_u)))
        correction_dt_norms.append(float(np.linalg.norm(delta_ut)))
    correction_norms_array = np.asarray(correction_norms, dtype=float)
    correction_dt_norms_array = np.asarray(correction_dt_norms, dtype=float)
    return (
        max(velocity_errors),
        max(velocity_dt_errors),
        float(np.sqrt(np.mean(correction_norms_array * correction_norms_array))),
        float(np.sqrt(np.mean(correction_dt_norms_array * correction_dt_norms_array))),
    )


def audit_finite_correction_transition(
    before_velocity: VectorFieldProvider,
    before_velocity_dt: VectorFieldProvider,
    before_pressure: ScalarFieldProvider,
    before_restricted_forcing: RestrictedForcingProvider,
    after_velocity: VectorFieldProvider,
    after_velocity_dt: VectorFieldProvider,
    after_pressure: ScalarFieldProvider,
    after_restricted_forcing: RestrictedForcingProvider,
    correction: CorrectionFieldProvider,
    held_in_points: Sequence[Sequence[float]] | np.ndarray,
    held_out_points: Sequence[Sequence[float]] | np.ndarray,
    update_check_points: Sequence[Sequence[float]] | np.ndarray,
    *,
    viscosity: float,
    spatial_step: float,
) -> FiniteCorrectionCycleReport:
    """Recompute before/after raw defect and finite-step gain from field providers."""

    before_identity, after_identity = _assert_transition_identity(
        before_velocity,
        before_velocity_dt,
        before_pressure,
        before_restricted_forcing,
        after_velocity,
        after_velocity_dt,
        after_pressure,
        after_restricted_forcing,
        correction,
    )
    held_in = _as_points(held_in_points, label="held_in_points")
    held_out = _as_points(held_out_points, label="held_out_points")
    update_points = _as_points(update_check_points, label="update_check_points")
    held_in_rows = {tuple(float(v) for v in row) for row in held_in}
    held_out_rows = {tuple(float(v) for v in row) for row in held_out}
    if held_in_rows & held_out_rows:
        raise ValueError("held-in and held-out point sets must be disjoint")
    if not np.isfinite(spatial_step) or spatial_step <= 0.0:
        raise ValueError("spatial_step must be finite and positive")

    velocity_closure, velocity_dt_closure, correction_rms, correction_dt_rms = (
        _update_closure(
            before_velocity,
            before_velocity_dt,
            after_velocity,
            after_velocity_dt,
            correction,
            update_points,
        )
    )
    if velocity_closure > UPDATE_CLOSURE_ATOL:
        raise ValueError("u^(k+1)=u^k+delta_u_k closure failed")
    if velocity_dt_closure > UPDATE_CLOSURE_ATOL:
        raise ValueError("time-derivative correction closure failed")
    correction_values = np.vstack(
        [
            _eval_vector(
                correction.velocity_evaluator,
                point,
                source_ref=correction.source_ref,
            )
            for point in update_points
        ]
    )
    correction_point_norms = np.linalg.norm(correction_values, axis=1)
    correction_max = float(np.max(correction_point_norms))
    correction_nontrivial = bool(correction_rms > 0.0 and correction_max > 0.0)

    before = evaluate_disjoint_heldin_heldout(
        before_velocity,
        before_velocity_dt,
        before_pressure,
        before_restricted_forcing,
        held_in,
        held_out,
        viscosity=viscosity,
        spatial_step=spatial_step,
    )
    after = evaluate_disjoint_heldin_heldout(
        after_velocity,
        after_velocity_dt,
        after_pressure,
        after_restricted_forcing,
        held_in,
        held_out,
        viscosity=viscosity,
        spatial_step=spatial_step,
    )

    held_in_before_rms = before["held_in"].vector_rms
    held_in_after_rms = after["held_in"].vector_rms
    held_out_before_rms = before["held_out"].vector_rms
    held_out_after_rms = after["held_out"].vector_rms
    held_in_before_max = before["held_in"].vector_max
    held_in_after_max = after["held_in"].vector_max
    held_out_before_max = before["held_out"].vector_max
    held_out_after_max = after["held_out"].vector_max

    finite_step_gain_guard_passed = bool(
        correction_nontrivial
        and held_in_after_rms < held_in_before_rms
        and held_out_after_rms < held_out_before_rms
        and _not_increased(held_in_after_max, held_in_before_max)
        and _not_increased(held_out_after_max, held_out_before_max)
    )

    held_in_divergence_before = _divergence_stats(
        before_velocity, held_in, spatial_step
    )
    held_in_divergence_after = _divergence_stats(
        after_velocity, held_in, spatial_step
    )
    held_out_divergence_before = _divergence_stats(
        before_velocity, held_out, spatial_step
    )
    held_out_divergence_after = _divergence_stats(
        after_velocity, held_out, spatial_step
    )
    correction_divergence = _correction_divergence_stats(
        correction, update_points, spatial_step
    )

    return FiniteCorrectionCycleReport(
        from_identity=before_identity,
        to_identity=after_identity,
        correction_source_ref=correction.source_ref,
        update_velocity_closure_max=velocity_closure,
        update_velocity_dt_closure_max=velocity_dt_closure,
        correction_vector_rms=correction_rms,
        correction_vector_max=correction_max,
        correction_dt_vector_rms=correction_dt_rms,
        correction_divergence_rms=correction_divergence[0],
        correction_divergence_max=correction_divergence[1],
        held_in_before_rms=held_in_before_rms,
        held_in_after_rms=held_in_after_rms,
        held_in_before_max=held_in_before_max,
        held_in_after_max=held_in_after_max,
        held_out_before_rms=held_out_before_rms,
        held_out_after_rms=held_out_after_rms,
        held_out_before_max=held_out_before_max,
        held_out_after_max=held_out_after_max,
        held_in_rms_ratio=_ratio(held_in_after_rms, held_in_before_rms),
        held_in_max_ratio=_ratio(held_in_after_max, held_in_before_max),
        held_out_rms_ratio=_ratio(held_out_after_rms, held_out_before_rms),
        held_out_max_ratio=_ratio(held_out_after_max, held_out_before_max),
        held_in_divergence_before_rms=held_in_divergence_before[0],
        held_in_divergence_after_rms=held_in_divergence_after[0],
        held_in_divergence_before_max=held_in_divergence_before[1],
        held_in_divergence_after_max=held_in_divergence_after[1],
        held_out_divergence_before_rms=held_out_divergence_before[0],
        held_out_divergence_after_rms=held_out_divergence_after[0],
        held_out_divergence_before_max=held_out_divergence_before[1],
        held_out_divergence_after_max=held_out_divergence_after[1],
        correction_nontrivial=correction_nontrivial,
        finite_step_gain_guard_passed=finite_step_gain_guard_passed,
    )


def admit_finite_correction_transition(*args: object, **kwargs: object) -> FiniteCorrectionCycleReport:
    """Audit and reject a finite step unless the actual-defect gain guard passes."""

    report = audit_finite_correction_transition(*args, **kwargs)
    if not report.finite_step_gain_guard_passed:
        raise ValueError("actual-defect finite-step gain guard rejected correction")
    return report


def truth_boundary() -> dict[str, object]:
    return {
        "finite_iteration_formula": "u^(k+1)=u^k+delta_u_k",
        "finite_correction_cycle_mechanics_executable": True,
        "before_after_actual_defect_recomputed_from_raw_providers": True,
        "heldin_heldout_defects_recomputed_separately": True,
        "correction_velocity_update_identity_checked": True,
        "correction_time_derivative_update_identity_checked": True,
        "correction_norm_recorded": True,
        "correction_divergence_recorded": True,
        "correction_nontriviality_recorded": True,
        "finite_step_actual_defect_gain_guard_executable": True,
        "divergent_correction_rejected_by_admission_api": True,
        "caller_supplied_residual_allowed": False,
        "caller_supplied_defect_allowed": False,
        "caller_supplied_stress_allowed": False,
        "caller_supplied_target_allowed": False,
        "caller_supplied_gain_allowed": False,
        "caller_supplied_scientific_threshold_allowed": False,
        "gain_guard_scope": "repository engineering raw-defect non-expansion; not formal Kokuno theorem gain",
        "formal_kokuno_correction_cycle_gain_bound_claimed": False,
        "restricted_forcing_semantics_independently_validated_here": False,
        "real_full_candidate_defect_consumed": False,
        "full_same_cycle_composite_requested_stress_materialized": False,
        "candidate_finite_head_mean_debt_materialized": False,
        "signed_mean_inverse_input_ready": False,
        "public_velocity_correction_materialized": False,
        "real_candidate_finite_correction_cycle_run": False,
        "heldout_normalized_ns_residual_assessed": False,
        "residual_reduction_claimed_for_real_candidate": False,
        "pde_validated": False,
        "final_normalized_momentum_gate": 1.0e-3,
        "final_normalized_divergence_gate": 1.0e-5,
    }


def _analytic_cycle(scale_after: float = 0.5) -> tuple[object, ...]:
    before_id = CycleIdentity("analytic-finite-cycle-v1", 0, "before")
    after_id = CycleIdentity("analytic-finite-cycle-v1", 1, "after")
    before_velocity = VectorFieldProvider(
        before_id, "analytic:u=t*y e_x", lambda x, y, z, t: (t * y, 0.0, 0.0)
    )
    before_velocity_dt = VectorFieldProvider(
        before_id, "analytic:u_t=y e_x", lambda x, y, z, t: (y, 0.0, 0.0)
    )
    before_pressure = ScalarFieldProvider(
        before_id, "analytic:p=0", lambda x, y, z, t: 0.0
    )
    before_forcing = RestrictedForcingProvider(
        before_id,
        "analytic:fixed-zero-force",
        "analytic fixed zero forcing; independent of residual",
        lambda x, y, z, t: (0.0, 0.0, 0.0),
    )
    after_velocity = VectorFieldProvider(
        after_id,
        f"analytic:u={scale_after}*t*y e_x",
        lambda x, y, z, t: (scale_after * t * y, 0.0, 0.0),
    )
    after_velocity_dt = VectorFieldProvider(
        after_id,
        f"analytic:u_t={scale_after}*y e_x",
        lambda x, y, z, t: (scale_after * y, 0.0, 0.0),
    )
    after_pressure = ScalarFieldProvider(
        after_id, "analytic:p=0", lambda x, y, z, t: 0.0
    )
    after_forcing = RestrictedForcingProvider(
        after_id,
        "analytic:fixed-zero-force",
        "analytic fixed zero forcing; independent of residual",
        lambda x, y, z, t: (0.0, 0.0, 0.0),
    )
    delta_scale = scale_after - 1.0
    correction = CorrectionFieldProvider(
        before_id,
        after_id,
        f"analytic:delta_scale={delta_scale}",
        lambda x, y, z, t: (delta_scale * t * y, 0.0, 0.0),
        lambda x, y, z, t: (delta_scale * y, 0.0, 0.0),
    )
    return (
        before_velocity,
        before_velocity_dt,
        before_pressure,
        before_forcing,
        after_velocity,
        after_velocity_dt,
        after_pressure,
        after_forcing,
        correction,
    )


def deterministic_receipt() -> dict[str, object]:
    providers = _analytic_cycle(0.5)
    held_in = [
        (-0.31, 0.27, 0.14, 0.37),
        (0.22, -0.19, -0.33, 0.63),
        (0.17, 0.41, -0.12, 0.46),
    ]
    held_out = [
        (0.11, 0.36, -0.21, 0.41),
        (-0.28, -0.17, 0.29, 0.59),
        (0.33, -0.44, 0.18, 0.52),
    ]
    update_points = [
        (0.09, 0.24, -0.16, 0.35),
        (-0.21, -0.32, 0.27, 0.55),
        (0.25, 0.38, 0.11, 0.67),
    ]
    report = admit_finite_correction_transition(
        *providers,
        held_in,
        held_out,
        update_points,
        viscosity=0.01,
        spatial_step=1.0e-3,
    )
    return {
        "schema": "kokuno-a3-finite-correction-cycle-contract-v1",
        "task": "KOKUNO-A3-FINITE-CORRECTION-CYCLE-CONTRACT-047",
        "provenance": {
            "parent_agent3_pr": 627,
            "parent_agent3_head": "0504056a6a5c9fb61193c68c1227b83dbd0a63ea",
            "older_joint_gain_reference_pr": 449,
            "older_joint_gain_reference_head": "383eeb622e40bc4ac46f02a369bd205ca003adc8",
            "purpose": "actual raw-defect before/after cycle mechanics while Agent-1 full composite handoff remains absent",
        },
        "analytic_mechanics_regression": report.to_receipt(),
        "truth_boundary": truth_boundary(),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    receipt = deterministic_receipt()
    text = json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(text, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    main()
