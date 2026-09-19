"""Held-in-only damping search for one Kokuno Agent-3 finite correction step.

The correction direction is fixed before the search.  A small repository-frozen
finite grid of positive damping factors is evaluated using *only* recomputed
held-in Navier--Stokes momentum defect.  The factor with the smallest eligible
held-in RMS is selected.  Only then is the selected transition evaluated on the
untouched held-out cloud through the existing finite-correction transition
contract.

This is a repository engineering line-search/rejection mechanism.  It is not
Kokuno's formal correction-cycle gain theorem, does not tune on held-out data,
and does not alter the final normalized 1e-3 momentum or 1e-5 divergence gates.
Callers cannot inject a residual, defect, stress, target, gain, score, or
scientific threshold.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Callable, Sequence

import numpy as np

from .kokuno_finite_correction_cycle_contract import (
    CorrectionFieldProvider,
    FiniteCorrectionCycleReport,
    audit_finite_correction_transition,
)
from .kokuno_finite_correction_cycle_ledger import CycleStateProviders
from .kokuno_same_cycle_defect_contract import (
    CycleIdentity,
    RestrictedForcingProvider,
    ScalarFieldProvider,
    VectorFieldProvider,
    evaluate_same_cycle_momentum_defect,
)

VectorEvaluator = Callable[[float, float, float, float], Sequence[float]]

# Frozen repository engineering grid.  This is not a source/paper constant.
DAMPING_GRID: tuple[float, ...] = (1.0, 0.5, 0.25, 0.125)


@dataclass(frozen=True)
class CorrectionDirectionProvider:
    """A fixed correction direction, prior to choosing a damping factor."""

    from_identity: CycleIdentity
    source_ref: str
    velocity_evaluator: VectorEvaluator
    velocity_dt_evaluator: VectorEvaluator

    def __post_init__(self) -> None:
        if not self.source_ref.strip():
            raise ValueError("source_ref must be non-empty")


@dataclass(frozen=True)
class HeldInDampingTrial:
    alpha: float
    held_in_before_rms: float
    held_in_after_rms: float
    held_in_before_max: float
    held_in_after_max: float
    held_in_rms_ratio: float
    held_in_max_ratio: float
    held_in_eligible: bool


@dataclass(frozen=True)
class DampedCorrectionSearchReport:
    cycle_id: str
    from_cycle_index: int
    damping_grid: tuple[float, ...]
    trials: tuple[HeldInDampingTrial, ...]
    selected_alpha: float
    selected_trial_index: int
    heldout_used_for_alpha_selection: bool
    heldout_trial_evaluation_count: int
    selected_transition: FiniteCorrectionCycleReport
    selected_step_admitted: bool

    def to_receipt(self) -> dict[str, object]:
        payload = asdict(self)
        payload["trials"] = [asdict(trial) for trial in self.trials]
        payload["selected_transition"] = self.selected_transition.to_receipt()
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


def _ratio(after: float, before: float) -> float:
    if before == 0.0:
        return 1.0 if after == 0.0 else float("inf")
    return float(after / before)


def _validate_grid(damping_grid: Sequence[float]) -> tuple[float, ...]:
    grid = tuple(float(alpha) for alpha in damping_grid)
    if not grid:
        raise ValueError("damping grid must be non-empty")
    if any((not np.isfinite(alpha)) or alpha <= 0.0 or alpha > 1.0 for alpha in grid):
        raise ValueError("damping factors must be finite and lie in (0, 1]")
    if len(set(grid)) != len(grid):
        raise ValueError("damping grid must not contain duplicates")
    if grid != DAMPING_GRID:
        raise ValueError("scientific callers cannot replace the frozen damping grid")
    return grid


def _vector(evaluator: VectorEvaluator, args: tuple[float, float, float, float]) -> np.ndarray:
    value = np.asarray(evaluator(*args), dtype=float)
    if value.shape != (3,) or not np.all(np.isfinite(value)):
        raise ValueError("correction direction must return a finite 3-vector")
    return value


def _trial_identity(before: CycleIdentity, alpha: float) -> CycleIdentity:
    token = f"{before.state_token}|damped-alpha={alpha:.6g}"
    return CycleIdentity(before.cycle_id, before.cycle_index + 1, token)


def _build_trial(
    before: CycleStateProviders,
    direction: CorrectionDirectionProvider,
    alpha: float,
) -> tuple[CycleStateProviders, CorrectionFieldProvider]:
    if direction.from_identity != before.identity:
        raise ValueError("correction direction does not bind the exact before-state identity")
    after_id = _trial_identity(before.identity, alpha)

    def scaled_delta(x: float, y: float, z: float, t: float) -> tuple[float, float, float]:
        args = (x, y, z, t)
        value = alpha * _vector(direction.velocity_evaluator, args)
        return tuple(float(v) for v in value)

    def scaled_delta_dt(x: float, y: float, z: float, t: float) -> tuple[float, float, float]:
        args = (x, y, z, t)
        value = alpha * _vector(direction.velocity_dt_evaluator, args)
        return tuple(float(v) for v in value)

    def after_velocity(x: float, y: float, z: float, t: float) -> tuple[float, float, float]:
        args = (x, y, z, t)
        base = np.asarray(before.velocity.evaluator(*args), dtype=float)
        value = base + alpha * _vector(direction.velocity_evaluator, args)
        if base.shape != (3,) or not np.all(np.isfinite(base)) or not np.all(np.isfinite(value)):
            raise ValueError("before/trial velocity must be a finite 3-vector")
        return tuple(float(v) for v in value)

    def after_velocity_dt(x: float, y: float, z: float, t: float) -> tuple[float, float, float]:
        args = (x, y, z, t)
        base = np.asarray(before.velocity_dt.evaluator(*args), dtype=float)
        value = base + alpha * _vector(direction.velocity_dt_evaluator, args)
        if base.shape != (3,) or not np.all(np.isfinite(base)) or not np.all(np.isfinite(value)):
            raise ValueError("before/trial velocity_dt must be a finite 3-vector")
        return tuple(float(v) for v in value)

    state = CycleStateProviders(
        VectorFieldProvider(after_id, before.velocity.source_ref + f":damped:{alpha:g}", after_velocity),
        VectorFieldProvider(after_id, before.velocity_dt.source_ref + f":damped:{alpha:g}", after_velocity_dt),
        ScalarFieldProvider(after_id, before.pressure.source_ref, before.pressure.evaluator),
        RestrictedForcingProvider(
            after_id,
            before.restricted_forcing.source_ref,
            before.restricted_forcing.restriction_receipt,
            before.restricted_forcing.evaluator,
        ),
    )
    correction = CorrectionFieldProvider(
        before.identity,
        after_id,
        direction.source_ref + f":alpha={alpha:g}",
        scaled_delta,
        scaled_delta_dt,
    )
    return state, correction


def audit_damped_correction_search(
    before: CycleStateProviders,
    direction: CorrectionDirectionProvider,
    held_in_points: Sequence[Sequence[float]] | np.ndarray,
    held_out_points: Sequence[Sequence[float]] | np.ndarray,
    update_check_points: Sequence[Sequence[float]] | np.ndarray,
    *,
    viscosity: float,
    spatial_step: float,
    damping_grid: Sequence[float] = DAMPING_GRID,
) -> DampedCorrectionSearchReport:
    """Select alpha with held-in actual defect only; assess held-out once afterward."""

    grid = _validate_grid(damping_grid)
    if direction.from_identity != before.identity:
        raise ValueError("direction/from-state identity mismatch")
    held_in = _as_points(held_in_points, label="held_in_points")
    held_out = _as_points(held_out_points, label="held_out_points")
    update = _as_points(update_check_points, label="update_check_points")
    held_in_rows = {tuple(float(v) for v in row) for row in held_in}
    held_out_rows = {tuple(float(v) for v in row) for row in held_out}
    if held_in_rows & held_out_rows:
        raise ValueError("held-in and held-out point sets must be disjoint")

    before_defect = evaluate_same_cycle_momentum_defect(
        before.velocity,
        before.velocity_dt,
        before.pressure,
        before.restricted_forcing,
        held_in,
        viscosity=viscosity,
        spatial_step=spatial_step,
    )

    trials: list[HeldInDampingTrial] = []
    trial_states: list[CycleStateProviders] = []
    trial_corrections: list[CorrectionFieldProvider] = []
    for alpha in grid:
        state, correction = _build_trial(before, direction, alpha)
        after_defect = evaluate_same_cycle_momentum_defect(
            state.velocity,
            state.velocity_dt,
            state.pressure,
            state.restricted_forcing,
            held_in,
            viscosity=viscosity,
            spatial_step=spatial_step,
        )
        rms_ratio = _ratio(after_defect.vector_rms, before_defect.vector_rms)
        max_ratio = _ratio(after_defect.vector_max, before_defect.vector_max)
        eligible = bool(
            after_defect.vector_rms < before_defect.vector_rms
            and after_defect.vector_max <= np.nextafter(before_defect.vector_max, np.inf)
        )
        trials.append(
            HeldInDampingTrial(
                alpha=alpha,
                held_in_before_rms=before_defect.vector_rms,
                held_in_after_rms=after_defect.vector_rms,
                held_in_before_max=before_defect.vector_max,
                held_in_after_max=after_defect.vector_max,
                held_in_rms_ratio=rms_ratio,
                held_in_max_ratio=max_ratio,
                held_in_eligible=eligible,
            )
        )
        trial_states.append(state)
        trial_corrections.append(correction)

    eligible_indices = [index for index, trial in enumerate(trials) if trial.held_in_eligible]
    if not eligible_indices:
        raise ValueError("no frozen-grid damping factor improves held-in actual defect")
    selected_index = min(
        eligible_indices,
        key=lambda index: (trials[index].held_in_after_rms, -trials[index].alpha),
    )

    # Held-out data first enters here, after alpha is irreversibly selected from held-in only.
    selected_transition = audit_finite_correction_transition(
        before.velocity,
        before.velocity_dt,
        before.pressure,
        before.restricted_forcing,
        trial_states[selected_index].velocity,
        trial_states[selected_index].velocity_dt,
        trial_states[selected_index].pressure,
        trial_states[selected_index].restricted_forcing,
        trial_corrections[selected_index],
        held_in,
        held_out,
        update,
        viscosity=viscosity,
        spatial_step=spatial_step,
    )
    return DampedCorrectionSearchReport(
        cycle_id=before.identity.cycle_id,
        from_cycle_index=before.identity.cycle_index,
        damping_grid=grid,
        trials=tuple(trials),
        selected_alpha=trials[selected_index].alpha,
        selected_trial_index=selected_index,
        heldout_used_for_alpha_selection=False,
        heldout_trial_evaluation_count=1,
        selected_transition=selected_transition,
        selected_step_admitted=bool(selected_transition.finite_step_gain_guard_passed),
    )


def admit_damped_correction_search(*args: object, **kwargs: object) -> DampedCorrectionSearchReport:
    """Fail closed if the held-in-selected step does not independently pass held-out."""

    report = audit_damped_correction_search(*args, **kwargs)
    if not report.selected_step_admitted:
        raise ValueError("held-in-selected damped correction failed held-out actual-defect admission")
    return report


def truth_boundary() -> dict[str, object]:
    return {
        "heldin_only_frozen_grid_damping_search_executable": True,
        "damping_grid": list(DAMPING_GRID),
        "damping_grid_classification": "repository engineering; not a Kokuno source constant",
        "alpha_selection_uses_heldout": False,
        "heldout_evaluated_only_after_alpha_selection": True,
        "heldout_failure_causes_rejection_without_retuning": True,
        "actual_defect_recomputed_from_raw_providers": True,
        "caller_supplied_residual_allowed": False,
        "caller_supplied_defect_allowed": False,
        "caller_supplied_stress_allowed": False,
        "caller_supplied_target_allowed": False,
        "caller_supplied_gain_allowed": False,
        "caller_supplied_normalized_score_allowed": False,
        "caller_supplied_scientific_threshold_allowed": False,
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


def _analytic_before() -> CycleStateProviders:
    identity = CycleIdentity("analytic-heldin-damping-v1", 0, "scale-1")
    return CycleStateProviders(
        VectorFieldProvider(identity, "analytic:u=t*y e_x", lambda x, y, z, t: (t * y, 0.0, 0.0)),
        VectorFieldProvider(identity, "analytic:u_t=y e_x", lambda x, y, z, t: (y, 0.0, 0.0)),
        ScalarFieldProvider(identity, "analytic:p=0", lambda x, y, z, t: 0.0),
        RestrictedForcingProvider(
            identity,
            "analytic:fixed-zero-forcing",
            "regression-zero-forcing; fixed before residual evaluation",
            lambda x, y, z, t: (0.0, 0.0, 0.0),
        ),
    )


def _analytic_direction() -> CorrectionDirectionProvider:
    before = _analytic_before()
    return CorrectionDirectionProvider(
        before.identity,
        "analytic:oversized-delta-scale=-2.4",
        lambda x, y, z, t: (-2.4 * t * y, 0.0, 0.0),
        lambda x, y, z, t: (-2.4 * y, 0.0, 0.0),
    )


def _regression_points() -> tuple[list[tuple[float, float, float, float]], list[tuple[float, float, float, float]], list[tuple[float, float, float, float]]]:
    held_in = [
        (-0.31, 0.27, 0.14, 0.37),
        (0.22, -0.19, -0.33, 0.63),
        (0.17, 0.41, -0.12, 0.46),
    ]
    held_out = [
        (0.11, 0.36, -0.21, 0.41),
        (-0.28, -0.17, 0.29, 0.59),
        (0.33, 0.23, 0.18, 0.54),
    ]
    update = [
        (-0.13, 0.29, 0.08, 0.35),
        (0.26, -0.24, -0.15, 0.51),
        (0.19, 0.31, 0.27, 0.67),
    ]
    return held_in, held_out, update


def deterministic_receipt() -> dict[str, object]:
    held_in, held_out, update = _regression_points()
    report = admit_damped_correction_search(
        _analytic_before(),
        _analytic_direction(),
        held_in,
        held_out,
        update,
        viscosity=0.01,
        spatial_step=1.0e-3,
    )
    return {
        "schema": "kokuno-a3-heldin-damped-correction-search-v1",
        "task": "KOKUNO-A3-HELDIN-DAMPED-SEARCH-049",
        "provenance": {
            "parent_agent3_pr": 645,
            "parent_agent3_head": "ce93d1681f6ad879b9ca2d33e595ca479cacbf9f",
            "agent1_latest_interface_pr_at_claim": 644,
            "agent2_latest_interface_pr_at_claim": 643,
        },
        "analytic_mechanics_regression": {
            "before_scale": 1.0,
            "undamped_direction_scale": -2.4,
            "expected_selected_alpha": 0.5,
            "expected_selected_residual_ratio": 0.2,
            "report": report.to_receipt(),
        },
        "truth_boundary": truth_boundary(),
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    receipt = deterministic_receipt()
    text = json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(text, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
