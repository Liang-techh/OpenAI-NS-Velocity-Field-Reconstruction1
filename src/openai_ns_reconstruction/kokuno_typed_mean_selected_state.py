"""Expose the already-admitted damped Agent-3 step as the next raw cycle state.

PR #743 closes the typed correction path through held-in-only damping and a
single held-out actual-defect admission, but returns the audit report rather
than the concrete ``CycleStateProviders`` needed to recompute the next
correction from ``u^(k+1)``.  This module closes only that propagation seam.

The selected alpha is *not* supplied by the caller.  We first run #743's
existing admission path, then replay the exact internal damped-trial constructor
at the already-selected alpha to recover the selected after-state and correction
provider.  No held-out data are evaluated a second time and no new damping,
residual, mean, radial-stress, or complete-curl rule is introduced.

This is correction-cycle mechanics only.  It does not provide the missing real
full Kokuno candidate, source-certified epsilon/H_ref/base amplitudes, or a
scientific residual reduction.
"""
from __future__ import annotations

import inspect
from dataclasses import asdict, dataclass
from typing import Any, Sequence

import numpy as np

from . import kokuno_finite_correction_damping as damping
from .kokuno_finite_correction_cycle_contract import CorrectionFieldProvider
from .kokuno_finite_correction_cycle_ledger import CycleStateProviders
from .kokuno_same_cycle_defect_contract import CycleIdentity
from .kokuno_typed_mean_amplitude_differential import MeanReferenceBaseAmplitudes
from .kokuno_typed_mean_damped_cycle_bridge import (
    TypedMeanDampedCycleBridgeReport,
    materialize_and_admit_typed_mean_damped_step,
)
from .kokuno_typed_mean_inverse_rhs import MeanInverseReferenceScale
from .kokuno_typed_mean_radial_stress import RadialMeanStressGeometry
from .kokuno_typed_mean_reference_inverse import MeanReferenceInverseOperator

TASK = "KOKUNO-A3-TYPED-SELECTED-STATE-PROPAGATION-061"
SCHEMA = "kokuno-a3-typed-selected-state-propagation-v1"
PARENT_AGENT3_PR = 743
PARENT_AGENT3_HEAD = "d751a08c8f72a6e29221df03ecea587d0f749f61"
FINAL_NORMALIZED_MOMENTUM_GATE = 1.0e-3
FINAL_NORMALIZED_DIVERGENCE_GATE = 1.0e-5


@dataclass(frozen=True)
class TypedSelectedDampedStep:
    """One admitted typed step plus concrete providers for ``u^(k+1)``."""

    before_identity: CycleIdentity
    after_state: CycleStateProviders
    correction: CorrectionFieldProvider
    bridge_report: TypedMeanDampedCycleBridgeReport

    @property
    def after_identity(self) -> CycleIdentity:
        return self.after_state.identity

    def to_receipt(self) -> dict[str, object]:
        # Raw providers contain callables and intentionally are not serialized.
        return {
            "schema": SCHEMA,
            "task": TASK,
            "before_identity": asdict(self.before_identity),
            "after_identity": asdict(self.after_identity),
            "correction_source_ref": self.correction.source_ref,
            "selected_alpha": self.bridge_report.damped_search.selected_alpha,
            "heldout_trial_evaluation_count": (
                self.bridge_report.damped_search.heldout_trial_evaluation_count
            ),
            "bridge_report": self.bridge_report.to_receipt(),
            "truth_boundary": truth_boundary(),
        }


def _direction_from_bridge(
    before: CycleStateProviders,
    report: TypedMeanDampedCycleBridgeReport,
) -> damping.CorrectionDirectionProvider:
    handoff = report.handoff.correction_provider
    if report.identity != before.identity:
        raise ValueError("typed bridge report does not bind the exact before-state identity")
    if handoff.from_identity != before.identity:
        raise ValueError("typed bridge correction does not bind the exact before-state identity")
    return damping.CorrectionDirectionProvider(
        from_identity=before.identity,
        source_ref=handoff.source_ref + ":undamped-direction",
        velocity_evaluator=handoff.velocity_evaluator,
        velocity_dt_evaluator=handoff.velocity_dt_evaluator,
    )


def _replay_selected_trial(
    before: CycleStateProviders,
    report: TypedMeanDampedCycleBridgeReport,
) -> tuple[CycleStateProviders, CorrectionFieldProvider]:
    search = report.damped_search
    if not search.selected_step_admitted:
        raise ValueError("typed bridge did not admit its selected damped step")
    if search.heldout_used_for_alpha_selection:
        raise ValueError("selected alpha was contaminated by held-out data")
    if search.heldout_trial_evaluation_count != 1:
        raise ValueError("typed bridge must evaluate held-out exactly once")
    if search.selected_alpha not in damping.DAMPING_GRID:
        raise ValueError("selected alpha is not on the frozen damping grid")

    direction = _direction_from_bridge(before, report)
    # Reuse the exact constructor already exercised inside #743/#650.  This
    # materializes providers only; it does not reevaluate held-in or held-out.
    after_state, correction = damping._build_trial(  # noqa: SLF001 - intentional intra-package reuse
        before,
        direction,
        search.selected_alpha,
    )
    selected = search.selected_transition
    if selected.from_identity != before.identity:
        raise ValueError("selected transition does not start from the exact before state")
    if after_state.identity != selected.to_identity:
        raise ValueError("replayed selected after-state identity drifted from audited transition")
    if correction.from_identity != selected.from_identity:
        raise ValueError("replayed correction from-identity drifted from audited transition")
    if correction.to_identity != selected.to_identity:
        raise ValueError("replayed correction to-identity drifted from audited transition")
    if correction.source_ref != selected.correction_source_ref:
        raise ValueError("replayed correction provenance drifted from audited transition")
    if after_state.pressure.source_ref != before.pressure.source_ref:
        raise ValueError("pressure source changed while materializing selected correction state")
    if after_state.restricted_forcing.source_ref != before.restricted_forcing.source_ref:
        raise ValueError("restricted forcing source changed while materializing selected correction state")
    if (
        after_state.restricted_forcing.restriction_receipt
        != before.restricted_forcing.restriction_receipt
    ):
        raise ValueError("restricted forcing receipt changed while materializing selected state")
    return after_state, correction


def materialize_selected_typed_mean_damped_step(
    before: CycleStateProviders,
    geometry: RadialMeanStressGeometry,
    reference_scale: MeanInverseReferenceScale,
    reference_operator: MeanReferenceInverseOperator,
    base_amplitudes: MeanReferenceBaseAmplitudes,
    agent2_backend: Any,
    held_in_points: Sequence[Sequence[float]] | np.ndarray,
    held_out_points: Sequence[Sequence[float]] | np.ndarray,
    update_check_points: Sequence[Sequence[float]] | np.ndarray,
    *,
    viscosity: float,
    spatial_step: float,
) -> TypedSelectedDampedStep:
    """Run the existing typed guarded step and expose its selected next state.

    The public API deliberately has no alpha/damping-grid argument.  Selection
    remains owned by #650/#743 on held-in actual defect only; held-out is used
    exactly once by the inherited admission path.  This function only recovers
    the concrete providers corresponding to that already-admitted transition.
    """

    if not isinstance(before, CycleStateProviders):
        raise TypeError("before must be CycleStateProviders")
    report = materialize_and_admit_typed_mean_damped_step(
        before,
        geometry,
        reference_scale,
        reference_operator,
        base_amplitudes,
        agent2_backend,
        held_in_points,
        held_out_points,
        update_check_points,
        viscosity=viscosity,
        spatial_step=spatial_step,
    )
    after_state, correction = _replay_selected_trial(before, report)
    return TypedSelectedDampedStep(
        before_identity=before.identity,
        after_state=after_state,
        correction=correction,
        bridge_report=report,
    )


def truth_boundary() -> dict[str, object]:
    signature = inspect.signature(materialize_selected_typed_mean_damped_step)
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
        "selected_damped_after_state_materialized": True,
        "selected_correction_provider_materialized": True,
        "state_propagation_reuses_existing_damped_trial_constructor": True,
        "selected_alpha_caller_supplied": False,
        "damping_grid_caller_supplied": False,
        "heldout_used_for_alpha_selection": False,
        "heldout_evaluated_again_during_state_materialization": False,
        "selected_transition_identity_replayed_fail_closed": True,
        "restricted_forcing_identity_preserved": True,
        "caller_supplied_residual_allowed": False,
        "caller_supplied_defect_allowed": False,
        "caller_supplied_stress_allowed": False,
        "caller_supplied_target_allowed": False,
        "caller_supplied_gain_allowed": False,
        "caller_supplied_normalized_score_allowed": False,
        "caller_supplied_scientific_threshold_allowed": False,
        "forbidden_public_parameters_absent": forbidden.isdisjoint(signature.parameters),
        "formal_kokuno_correction_cycle_gain_bound_claimed": False,
        "real_full_candidate_bound": False,
        "real_candidate_finite_correction_cycle_run": False,
        "heldout_normalized_ns_residual_assessed": False,
        "residual_reduction_claimed_for_real_candidate": False,
        "pde_validated": False,
        "final_normalized_momentum_gate": FINAL_NORMALIZED_MOMENTUM_GATE,
        "final_normalized_divergence_gate": FINAL_NORMALIZED_DIVERGENCE_GATE,
    }
