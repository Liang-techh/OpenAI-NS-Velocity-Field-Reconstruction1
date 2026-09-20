"""Bridge the typed Agent-3 correction chain into the guarded finite-step search.

This module closes one plumbing seam only.  It does not implement a new mean
operator, radial inverse, or complete curl.  Starting from one raw same-cycle
state it

    raw fields -> actual defect -> mean/radial stress -> DeltaC/epsilon
      -> H_ref^{-1} -> delta_y -> delta_a
      -> external Agent-2 complete-curl backend -> (delta_u, delta_u_t)
      -> held-in-only frozen damping search -> one held-out finite-step audit.

The caller cannot provide a residual, defect, stress, target, gain, normalized
score, damping grid, or scientific threshold.  The existing Agent-3 modules
recompute those quantities from field providers.  The Agent-2 backend is bound
through the exact seven-field adapter seam introduced in PR #735; no Agent-2
curl formula is copied here.

This bridge is executable mechanics, not evidence that a real full Kokuno
candidate exists or that the final normalized Navier--Stokes gate is satisfied.
"""
from __future__ import annotations

import inspect
from dataclasses import asdict, dataclass
from typing import Any, Sequence

import numpy as np

from .kokuno_agent2_complete_curl_binding import bind_agent2_complete_curl_backend
from .kokuno_finite_correction_cycle_ledger import CycleStateProviders
from .kokuno_finite_correction_damping import (
    CorrectionDirectionProvider,
    DampedCorrectionSearchReport,
    admit_damped_correction_search,
)
from .kokuno_same_cycle_defect_contract import CycleIdentity
from .kokuno_typed_mean_amplitude_differential import MeanReferenceBaseAmplitudes
from .kokuno_typed_mean_inverse_rhs import MeanInverseReferenceScale
from .kokuno_typed_mean_radial_stress import RadialMeanStressGeometry
from .kokuno_typed_mean_reference_inverse import MeanReferenceInverseOperator
from .kokuno_typed_mean_velocity_correction_handoff import (
    TypedMeanVelocityCorrectionHandoff,
    materialize_typed_mean_velocity_correction_handoff,
)

TASK = "KOKUNO-A3-TYPED-DAMPED-CYCLE-BRIDGE-060"
SCHEMA = "kokuno-a3-typed-damped-cycle-bridge-v1"
PARENT_AGENT3_PR = 735
PARENT_AGENT3_HEAD = "4692da3ee204604f1f9942739fc568782a5ee97d"
AGENT2_BACKEND_PR = 734
AGENT2_BACKEND_HEAD = "c535eeec3267869f630c1327506b3683d0417969"
FINAL_NORMALIZED_MOMENTUM_GATE = 1.0e-3
FINAL_NORMALIZED_DIVERGENCE_GATE = 1.0e-5


@dataclass(frozen=True)
class TypedMeanDampedCycleBridgeReport:
    """One typed correction direction plus its guarded finite-step audit."""

    identity: CycleIdentity
    provisional_to_identity: CycleIdentity
    handoff: TypedMeanVelocityCorrectionHandoff
    damped_search: DampedCorrectionSearchReport

    def to_receipt(self) -> dict[str, object]:
        return {
            "schema": SCHEMA,
            "task": TASK,
            "identity": asdict(self.identity),
            "provisional_to_identity": asdict(self.provisional_to_identity),
            "source_agent2_complete_curl_certified": (
                self.handoff.adapter.source_agent2_complete_curl_certified
            ),
            "handoff": self.handoff.to_receipt(),
            "damped_search": self.damped_search.to_receipt(),
            "truth_boundary": truth_boundary(),
        }


def _provisional_to_identity(identity: CycleIdentity) -> CycleIdentity:
    return CycleIdentity(
        identity.cycle_id,
        identity.cycle_index + 1,
        identity.state_token + "|typed-undamped-direction",
    )


def materialize_and_admit_typed_mean_damped_step(
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
) -> TypedMeanDampedCycleBridgeReport:
    """Materialize one correction direction and run the existing guarded step.

    The complete-curl backend is external.  The typed mean/radial chain is
    recomputed internally from ``before``.  Damping selection is delegated to
    the existing frozen-grid held-in-only search; held-out is evaluated exactly
    once after alpha selection and cannot trigger retuning.
    """

    if not isinstance(before, CycleStateProviders):
        raise TypeError("before must be CycleStateProviders")
    identity = before.identity
    adapter = bind_agent2_complete_curl_backend(agent2_backend, identity)
    provisional_to = _provisional_to_identity(identity)
    handoff = materialize_typed_mean_velocity_correction_handoff(
        before.velocity,
        before.velocity_dt,
        before.pressure,
        before.restricted_forcing,
        geometry,
        reference_scale,
        reference_operator,
        base_amplitudes,
        adapter,
        provisional_to,
        viscosity=viscosity,
        spatial_step=spatial_step,
    )
    if handoff.identity != identity:
        raise ValueError("typed correction handoff changed the before-state identity")
    if handoff.correction_provider.from_identity != identity:
        raise ValueError("typed correction provider does not bind the before-state identity")

    # The provisional to-identity exists only to satisfy the typed handoff's
    # one-step correction contract.  The damping selector owns the scientific
    # trial identities because alpha is selected later from held-in data only.
    direction = CorrectionDirectionProvider(
        from_identity=identity,
        source_ref=handoff.correction_provider.source_ref + ":undamped-direction",
        velocity_evaluator=handoff.correction_provider.velocity_evaluator,
        velocity_dt_evaluator=handoff.correction_provider.velocity_dt_evaluator,
    )
    damped = admit_damped_correction_search(
        before,
        direction,
        held_in_points,
        held_out_points,
        update_check_points,
        viscosity=viscosity,
        spatial_step=spatial_step,
    )
    return TypedMeanDampedCycleBridgeReport(
        identity=identity,
        provisional_to_identity=provisional_to,
        handoff=handoff,
        damped_search=damped,
    )


def truth_boundary() -> dict[str, object]:
    signature = inspect.signature(materialize_and_admit_typed_mean_damped_step)
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
        "damping_grid",
    }
    return {
        "parent_agent3_pr": PARENT_AGENT3_PR,
        "parent_agent3_head": PARENT_AGENT3_HEAD,
        "agent2_backend_pr": AGENT2_BACKEND_PR,
        "agent2_backend_head": AGENT2_BACKEND_HEAD,
        "raw_same_cycle_state_required": True,
        "typed_mean_radial_chain_recomputed_internally": True,
        "agent2_complete_curl_bound_without_reimplementation": True,
        "typed_delta_a_not_caller_supplied": True,
        "heldin_only_frozen_grid_damping_reused": True,
        "heldout_used_for_alpha_selection": False,
        "heldout_evaluated_once_after_selection": True,
        "caller_supplied_residual_allowed": False,
        "caller_supplied_defect_allowed": False,
        "caller_supplied_stress_allowed": False,
        "caller_supplied_target_allowed": False,
        "caller_supplied_gain_allowed": False,
        "caller_supplied_normalized_score_allowed": False,
        "caller_supplied_scientific_threshold_allowed": False,
        "caller_supplied_damping_grid_allowed": False,
        "forbidden_public_parameters_absent": forbidden.isdisjoint(signature.parameters),
        "source_certification_copied_from_agent2_backend_not_self_promoted": True,
        "formal_kokuno_correction_cycle_gain_bound_claimed": False,
        "real_full_candidate_bound": False,
        "real_candidate_finite_correction_cycle_run": False,
        "heldout_normalized_ns_residual_assessed": False,
        "residual_reduction_claimed_for_real_candidate": False,
        "pde_validated": False,
        "final_normalized_momentum_gate": FINAL_NORMALIZED_MOMENTUM_GATE,
        "final_normalized_divergence_gate": FINAL_NORMALIZED_DIVERGENCE_GATE,
    }
