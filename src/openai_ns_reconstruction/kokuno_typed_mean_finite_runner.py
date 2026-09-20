"""Finite typed mean-correction runner for Kokuno Agent 3.

This module closes the iteration seam between the already-selected damped next
state from Agent-3 #753 and the source-gain bookkeeping from #762.  It does not
introduce another residual evaluator, radial inverse, complete-curl formula, or
damping rule.

For each finite stage it:
1. receives the exact currently selected raw ``CycleStateProviders``;
2. obtains cycle-local typed mean-correction witnesses from a factory that sees
   only that selected state;
3. runs the existing typed actual-defect -> mean/radial -> signed correction ->
   Agent-2 complete-curl -> held-in damping -> one-shot held-out admission path;
4. feeds the exact admitted after-state directly into the next stage; and
5. attaches the corrected-reader source exponent bookkeeping for that stage.

The runner deliberately aggregates the transition receipts already produced by
the selected-step machinery instead of re-running the held-out audit through the
multi-step ledger.  Therefore held-out remains evaluated exactly once per
selected step and never participates in alpha selection.

This is finite correction-cycle orchestration only.  It is not a real Kokuno
full-candidate run, final normalized PDE validation, a paper-exact theorem, or a
blow-up proof.
"""
from __future__ import annotations

import argparse
import inspect
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Callable, Sequence

import numpy as np

from .kokuno_correction_cycle_gain_diagnostic import (
    SOURCE_READER_DATE,
    SOURCE_READER_HEAD,
    SOURCE_READER_REPO,
    SourceCorrectionStageGain,
    source_correction_stage_gain,
)
from .kokuno_finite_correction_cycle_ledger import CycleStateProviders
from .kokuno_typed_mean_amplitude_differential import MeanReferenceBaseAmplitudes
from .kokuno_typed_mean_inverse_rhs import MeanInverseReferenceScale
from .kokuno_typed_mean_radial_stress import RadialMeanStressGeometry
from .kokuno_typed_mean_reference_inverse import MeanReferenceInverseOperator
from .kokuno_typed_mean_selected_state import (
    TypedSelectedDampedStep,
    materialize_selected_typed_mean_damped_step,
)

TASK = "KOKUNO-A3-TYPED-MEAN-FINITE-RUNNER-063"
SCHEMA = "kokuno-a3-typed-mean-finite-runner-v1"
PARENT_AGENT3_PR = 762
PARENT_AGENT3_HEAD = "aee40eb5c82cf13d5ffb59fae8eba0d6b8da7490"

FINAL_NORMALIZED_MOMENTUM_GATE = 1.0e-3
FINAL_NORMALIZED_DIVERGENCE_GATE = 1.0e-5


@dataclass(frozen=True)
class TypedMeanStageInputs:
    """Cycle-local source/provenance witnesses for one typed correction stage."""

    geometry: RadialMeanStressGeometry
    reference_scale: MeanInverseReferenceScale
    reference_operator: MeanReferenceInverseOperator
    base_amplitudes: MeanReferenceBaseAmplitudes
    agent2_backend: Any


StageInputsFactory = Callable[[CycleStateProviders], TypedMeanStageInputs]


@dataclass(frozen=True)
class TypedMeanFiniteCycleRun:
    """Finite sequence of already-admitted typed mean-correction steps."""

    cycle_id: str
    step_count: int
    initial_identity: object
    final_state: CycleStateProviders
    steps: tuple[TypedSelectedDampedStep, ...]
    source_stages: tuple[SourceCorrectionStageGain, ...]
    selected_alphas: tuple[float, ...]
    held_in_rms_contraction_factors: tuple[float, ...]
    held_out_rms_contraction_factors: tuple[float, ...]
    held_in_max_contraction_factors: tuple[float, ...]
    held_out_max_contraction_factors: tuple[float, ...]
    cumulative_held_in_rms_ratio: float
    cumulative_held_out_rms_ratio: float
    cumulative_held_in_max_ratio: float
    cumulative_held_out_max_ratio: float
    cumulative_correction_rms_budget: float
    maximum_correction_rms: float
    maximum_correction_divergence_max: float
    heldout_trial_evaluation_count_total: int
    heldout_used_for_alpha_selection: bool
    every_step_nontrivial: bool
    every_step_actual_defect_gain_guard_passed: bool
    source_gain_bounds_passed: bool
    admitted: bool

    def to_receipt(self) -> dict[str, object]:
        return {
            "schema": SCHEMA,
            "task": TASK,
            "cycle_id": self.cycle_id,
            "step_count": self.step_count,
            "initial_identity": asdict(self.initial_identity),
            "final_identity": asdict(self.final_state.identity),
            "steps": [step.to_receipt() for step in self.steps],
            "source_stages": [stage.to_receipt() for stage in self.source_stages],
            "selected_alphas": list(self.selected_alphas),
            "held_in_rms_contraction_factors": list(
                self.held_in_rms_contraction_factors
            ),
            "held_out_rms_contraction_factors": list(
                self.held_out_rms_contraction_factors
            ),
            "held_in_max_contraction_factors": list(
                self.held_in_max_contraction_factors
            ),
            "held_out_max_contraction_factors": list(
                self.held_out_max_contraction_factors
            ),
            "cumulative_held_in_rms_ratio": self.cumulative_held_in_rms_ratio,
            "cumulative_held_out_rms_ratio": self.cumulative_held_out_rms_ratio,
            "cumulative_held_in_max_ratio": self.cumulative_held_in_max_ratio,
            "cumulative_held_out_max_ratio": self.cumulative_held_out_max_ratio,
            "cumulative_correction_rms_budget": self.cumulative_correction_rms_budget,
            "maximum_correction_rms": self.maximum_correction_rms,
            "maximum_correction_divergence_max": self.maximum_correction_divergence_max,
            "heldout_trial_evaluation_count_total": (
                self.heldout_trial_evaluation_count_total
            ),
            "heldout_used_for_alpha_selection": self.heldout_used_for_alpha_selection,
            "every_step_nontrivial": self.every_step_nontrivial,
            "every_step_actual_defect_gain_guard_passed": (
                self.every_step_actual_defect_gain_guard_passed
            ),
            "source_gain_bounds_passed": self.source_gain_bounds_passed,
            "admitted": self.admitted,
            "truth_boundary": truth_boundary(),
        }


def _ratio(final: float, initial: float) -> float:
    if initial == 0.0:
        return 1.0 if final == 0.0 else float("inf")
    return float(final / initial)


def _validate_step_count(step_count: int) -> int:
    if isinstance(step_count, bool) or not isinstance(step_count, int):
        raise TypeError("step_count must be an integer")
    if step_count < 1:
        raise ValueError("step_count must be at least one")
    return step_count


def _transition(step: TypedSelectedDampedStep) -> object:
    search = step.bridge_report.damped_search
    if search.heldout_used_for_alpha_selection:
        raise ValueError("held-out data contaminated damping selection")
    if search.heldout_trial_evaluation_count != 1:
        raise ValueError("each selected step must evaluate held-out exactly once")
    if not search.selected_step_admitted:
        raise ValueError("typed selected step was not admitted")
    transition = search.selected_transition
    if transition.from_identity != step.before_identity:
        raise ValueError("selected transition from-identity drifted")
    if transition.to_identity != step.after_identity:
        raise ValueError("selected transition to-identity drifted")
    if not transition.correction_nontrivial:
        raise ValueError("selected transition correction is trivial")
    if not transition.finite_step_gain_guard_passed:
        raise ValueError("selected transition failed actual-defect gain guard")
    return transition


def materialize_typed_mean_finite_cycle(
    initial_state: CycleStateProviders,
    stage_inputs_factory: StageInputsFactory,
    step_count: int,
    held_in_points: Sequence[Sequence[float]] | np.ndarray,
    held_out_points: Sequence[Sequence[float]] | np.ndarray,
    update_check_points: Sequence[Sequence[float]] | np.ndarray,
    *,
    viscosity: float,
    spatial_step: float,
) -> TypedMeanFiniteCycleRun:
    """Iterate the existing typed correction path over exact selected states.

    The factory is invoked only with the exact selected raw state for the
    current stage.  It may bind cycle-local source-certified epsilon/H_ref/base
    amplitudes and the external Agent-2 backend, but it never receives held-out
    values or a residual score from this runner.

    Held-out is evaluated inside the inherited selected-step path exactly once
    after held-in damping selection.  This runner consumes those transition
    reports and never re-audits them.
    """

    if not isinstance(initial_state, CycleStateProviders):
        raise TypeError("initial_state must be CycleStateProviders")
    if not callable(stage_inputs_factory):
        raise TypeError("stage_inputs_factory must be callable")
    count = _validate_step_count(step_count)
    if initial_state.identity.cycle_index != 0:
        raise ValueError("source-structured finite cycle must begin at cycle index 0")

    current = initial_state
    steps: list[TypedSelectedDampedStep] = []
    source_stages: list[SourceCorrectionStageGain] = []

    for expected_stage in range(count):
        if current.identity.cycle_index != expected_stage:
            raise ValueError("selected state cycle index drifted before next stage")
        inputs = stage_inputs_factory(current)
        if not isinstance(inputs, TypedMeanStageInputs):
            raise TypeError("stage_inputs_factory must return TypedMeanStageInputs")

        step = materialize_selected_typed_mean_damped_step(
            current,
            geometry=inputs.geometry,
            reference_scale=inputs.reference_scale,
            reference_operator=inputs.reference_operator,
            base_amplitudes=inputs.base_amplitudes,
            agent2_backend=inputs.agent2_backend,
            held_in_points=held_in_points,
            held_out_points=held_out_points,
            update_check_points=update_check_points,
            viscosity=viscosity,
            spatial_step=spatial_step,
        )
        if step.before_identity != current.identity:
            raise ValueError("typed selected step does not bind the exact current state")
        if step.after_identity.cycle_id != current.identity.cycle_id:
            raise ValueError("typed selected step changed cycle_id")
        if step.after_identity.cycle_index != current.identity.cycle_index + 1:
            raise ValueError("typed selected step must advance exactly one cycle index")
        if step.correction.from_identity != current.identity:
            raise ValueError("selected correction does not bind exact current state")
        if step.correction.to_identity != step.after_identity:
            raise ValueError("selected correction does not bind exact next state")

        _transition(step)
        source = source_correction_stage_gain(current.identity.cycle_index)
        if not source.source_bound_passed:
            raise ValueError("corrected-source stage gain arithmetic failed")

        steps.append(step)
        source_stages.append(source)
        current = step.after_state

    transitions = tuple(_transition(step) for step in steps)
    searches = tuple(step.bridge_report.damped_search for step in steps)

    first = transitions[0]
    last = transitions[-1]
    held_in_rms_factors = tuple(float(item.held_in_rms_ratio) for item in transitions)
    held_out_rms_factors = tuple(float(item.held_out_rms_ratio) for item in transitions)
    held_in_max_factors = tuple(float(item.held_in_max_ratio) for item in transitions)
    held_out_max_factors = tuple(float(item.held_out_max_ratio) for item in transitions)
    heldout_count = sum(int(search.heldout_trial_evaluation_count) for search in searches)
    heldout_used = any(bool(search.heldout_used_for_alpha_selection) for search in searches)
    nontrivial = all(bool(item.correction_nontrivial) for item in transitions)
    gain_passed = all(bool(item.finite_step_gain_guard_passed) for item in transitions)
    source_passed = all(stage.source_bound_passed for stage in source_stages)

    admitted = bool(
        not heldout_used
        and heldout_count == count
        and nontrivial
        and gain_passed
        and source_passed
        and last.held_in_after_rms < first.held_in_before_rms
        and last.held_out_after_rms < first.held_out_before_rms
    )
    if not admitted:
        raise ValueError("finite typed mean correction cycle failed aggregate admission")

    return TypedMeanFiniteCycleRun(
        cycle_id=initial_state.identity.cycle_id,
        step_count=count,
        initial_identity=initial_state.identity,
        final_state=current,
        steps=tuple(steps),
        source_stages=tuple(source_stages),
        selected_alphas=tuple(float(search.selected_alpha) for search in searches),
        held_in_rms_contraction_factors=held_in_rms_factors,
        held_out_rms_contraction_factors=held_out_rms_factors,
        held_in_max_contraction_factors=held_in_max_factors,
        held_out_max_contraction_factors=held_out_max_factors,
        cumulative_held_in_rms_ratio=_ratio(
            last.held_in_after_rms, first.held_in_before_rms
        ),
        cumulative_held_out_rms_ratio=_ratio(
            last.held_out_after_rms, first.held_out_before_rms
        ),
        cumulative_held_in_max_ratio=_ratio(
            last.held_in_after_max, first.held_in_before_max
        ),
        cumulative_held_out_max_ratio=_ratio(
            last.held_out_after_max, first.held_out_before_max
        ),
        cumulative_correction_rms_budget=float(
            sum(float(item.correction_vector_rms) for item in transitions)
        ),
        maximum_correction_rms=max(
            float(item.correction_vector_rms) for item in transitions
        ),
        maximum_correction_divergence_max=max(
            float(item.correction_divergence_max) for item in transitions
        ),
        heldout_trial_evaluation_count_total=heldout_count,
        heldout_used_for_alpha_selection=heldout_used,
        every_step_nontrivial=nontrivial,
        every_step_actual_defect_gain_guard_passed=gain_passed,
        source_gain_bounds_passed=source_passed,
        admitted=True,
    )


def truth_boundary() -> dict[str, object]:
    signature = inspect.signature(materialize_typed_mean_finite_cycle)
    forbidden = {
        "residual",
        "defect",
        "mean_source",
        "stress",
        "requested_stress",
        "debt",
        "delta_c",
        "inverse_rhs",
        "rhs",
        "delta_y",
        "delta_a",
        "target",
        "gain",
        "normalized_score",
        "scientific_threshold",
        "alpha",
        "damping_grid",
    }
    return {
        "parent_agent3_pr": PARENT_AGENT3_PR,
        "parent_agent3_head": PARENT_AGENT3_HEAD,
        "source_reader_repository": SOURCE_READER_REPO,
        "source_reader_head": SOURCE_READER_HEAD,
        "source_reader_date": SOURCE_READER_DATE,
        "typed_mean_finite_runner_executable": True,
        "selected_after_state_fed_directly_to_next_stage": True,
        "actual_defect_recomputed_inside_each_typed_step": True,
        "source_gain_recomputed_from_stage_index": True,
        "runner_reaudits_heldout_transitions": False,
        "heldout_evaluated_once_per_selected_step": True,
        "heldout_used_for_alpha_selection": False,
        "stage_factory_receives_only_selected_state": True,
        "caller_supplied_residual_allowed": False,
        "caller_supplied_defect_allowed": False,
        "caller_supplied_stress_allowed": False,
        "caller_supplied_target_allowed": False,
        "caller_supplied_gain_allowed": False,
        "caller_supplied_normalized_score_allowed": False,
        "caller_supplied_scientific_threshold_allowed": False,
        "forbidden_public_parameters_absent": forbidden.isdisjoint(signature.parameters),
        "formal_kokuno_correction_cycle_theorem_claimed": False,
        "source_reader_is_final_independent_validation": False,
        "real_full_candidate_bound": False,
        "real_candidate_finite_correction_cycle_run": False,
        "heldout_normalized_ns_residual_assessed": False,
        "residual_reduction_claimed_for_real_candidate": False,
        "pde_validated": False,
        "final_normalized_momentum_gate": FINAL_NORMALIZED_MOMENTUM_GATE,
        "final_normalized_divergence_gate": FINAL_NORMALIZED_DIVERGENCE_GATE,
    }


def deterministic_receipt() -> dict[str, object]:
    """Emit only frozen orchestration/source facts, not a fake candidate run."""

    return {
        "schema": SCHEMA,
        "task": TASK,
        "provenance": {
            "parent_agent3_pr": PARENT_AGENT3_PR,
            "parent_agent3_head": PARENT_AGENT3_HEAD,
            "source_reader_repository": SOURCE_READER_REPO,
            "source_reader_head": SOURCE_READER_HEAD,
            "source_reader_date": SOURCE_READER_DATE,
        },
        "source_stage_examples": [
            source_correction_stage_gain(0).to_receipt(),
            source_correction_stage_gain(1).to_receipt(),
        ],
        "runner_contract": {
            "selected_state_feedback": "u^(k+1) raw providers -> stage k+1 typed machinery",
            "heldout_reaudit_by_runner": False,
            "real_candidate_run": False,
        },
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
