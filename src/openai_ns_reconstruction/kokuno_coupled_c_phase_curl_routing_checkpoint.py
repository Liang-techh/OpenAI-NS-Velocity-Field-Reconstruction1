"""Agent-5 routing checkpoint for the coupled-C and mass-bound phase-curl seams.

This is one minimal integration increment on top of Agent-5 v35.  It records
two newer sibling deliveries without rewriting their mathematical lanes:

* Agent 1 now has a freshly repropagated, provenance-labelled coupled-C
  Appendix-B realization that is no longer rejected by the frozen PA.10
  pointwise necessary screen.  That is only a negative obstruction removed;
  PA.16, source B_0/T_sh, the global join and global pressure remain unproved.
* Agent 2 now composes the candidate numerical signed masses into the audited
  phase/complete-curl stack while requiring explicit slow derivatives of
  h_sigma.  It returns Q-scaled by-sign/by-beta/total Cartesian arrays for
  caller-supplied candidate data, but still does not expose a public
  velocity_osc(x,y,z,t) provider or recover hidden source background/modes.

The increment freezes typed handoff contracts for the four Kokuno lanes.  It
is not a candidate artifact and does not assess the held-out Navier--Stokes
PDE gate.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from .kokuno_numeric_prerequisites_routing_checkpoint import (
    AUTONOMOUS_REALIZATION_POLICY,
    FORMAL_GATES,
    SCHEMA as PARENT_SCHEMA,
    TASK_ID as PARENT_TASK_ID,
    build_checkpoint as build_parent_checkpoint,
    checkpoint_sha256,
    validate_checkpoint as validate_parent_checkpoint,
)

SCHEMA = "kokuno-agent5-coupled-c-phase-curl-routing-checkpoint-v36"
TASK_ID = "KOKUNO-A5-COUPLED-C-PHASE-CURL-ROUTING-036"

AGENT1_COUPLED_C_RECEIPT = {
    "pr": 521,
    "head": "c1cef3ac2b3e4d58ab4b721b31a1d2de14abba65",
    "evidence_class": "repository_autonomous_source_compatible",
    "dedicated_run": 35409636108,
    "dedicated_status": "success",
    "latest_standard_run_observed": 35409925427,
    "latest_standard_status_observed": "in_progress",
    "artifact": "kokuno-agent1-coupled-c-appendix-b-v1",
    "artifact_id": 10573590676,
    "artifact_digest": "sha256:8cac9e357d34b4ab15eba8ebc3cb6524cd6255cbe4a176e0f9ff8f905ec150fb",
    "log_C_factor": 32.0,
    "upstream_C_repropagated_before_core": True,
    "posthoc_C_multiplier_applied": False,
    "native_xyz_t_velocity_interface_executable": True,
    "boundary_values_interface_executable": True,
    "sha_json_roundtrip_executable": True,
    "pointwise_PA10_necessary_screen_passed": True,
    "previous_pointwise_PA10_obstruction_removed": True,
    "nontrivial_G_i_guard_passed": True,
    "passing_pointwise_screen_is_source_T_sh_certificate": False,
    "source_complex_C_bound_verified": False,
    "source_all_PA11_C_bounds_verified": False,
    "source_B0_analytic_bound_proved": False,
    "source_T_sh_lower_bound_verified": False,
    "actual_source_incoming_five_moment_discrepancy_bound": False,
    "selected_pa16_handoff_allowed": False,
    "inner_to_outer_join_completed": False,
    "global_pressure_matched": False,
    "global_leading_profile_reconstructed": False,
    "formal_full_domain_pde_gate_assessed": False,
    "pde_validated": False,
    "consumed_in_executable_ancestry": False,
}

AGENT2_PHASE_CURL_RECEIPT = {
    "pr": 522,
    "head": "296ed78b79cade56547aa9b751735207685c80bc",
    "evidence_class": "candidate_derived_composition",
    "dedicated_run": 35409788333,
    "standard_run": 35409788239,
    "dedicated_status": "success",
    "standard_status": "success",
    "focused_passed": 20,
    "focused_seconds": 0.40,
    "candidate_signed_mass_consumed_directly": True,
    "manual_h_sigma_reentry_removed": True,
    "explicit_D_r_h_sigma_required": True,
    "explicit_D_z_h_sigma_required": True,
    "frozen_local_zero_derivative_mode_requires_exact_zero": True,
    "beta_labels_revalidated": True,
    "Q_epsilon_Ls_schedule_revalidated": True,
    "source_numeric_or_paper_exact_laundering_rejected": True,
    "q_scaled_by_sign_cartesian_arrays_executable": True,
    "q_scaled_by_beta_cartesian_arrays_executable": True,
    "q_scaled_total_cartesian_array_executable": True,
    "public_xyz_t_velocity_osc_provider_ready": False,
    "actual_source_h_sigma_pulse_integrals_bound": False,
    "actual_positive_order_background_bound": False,
    "actual_auxiliary_torus_mode_family_bound": False,
    "actual_source_partition_labels_bound": False,
    "formal_full_domain_pde_gate_assessed": False,
    "pde_validated": False,
    "consumed_in_executable_ancestry": False,
}

PREVIOUS_AGENT5_RECEIPT = {
    "pr": 516,
    "head": "02d3c2b5663434101d9517a82d5cabe0fb667df6",
    "schema": "kokuno-agent5-numeric-prerequisites-routing-checkpoint-v35",
    "dedicated_run": 35407618350,
    "standard_run": 35407618244,
    "dedicated_status": "success",
    "standard_status": "success",
    "artifact": "kokuno-agent5-numeric-prerequisites-routing-checkpoint-v35",
    "artifact_id": 10572867184,
    "artifact_digest": "sha256:c070ae7a411fd2918d8a6964a24b3068f678a64f778e42e14f1312db8c686132",
    "checkpoint_sha256": "54721306b7563d6b1c32f8a8473e738e7271848d0e92e887000e009081aa30b5",
    "consumed_in_executable_ancestry": True,
}

TYPED_HANDOFF_CONTRACTS = {
    "leading": {
        "provider": "KokunoCoupledCAppendixBNormalization",
        "required_public_outputs": [
            "boundary_values(eta)",
            "velocity(x,y,z,t)->[...,3]",
            "sha-bound JSON save/load",
        ],
        "current_semantics": "pointwise_PA10_not_excluded_only",
        "next_required_evidence": [
            "source-compatible incoming five-moment bound",
            "source B0/T_sh certificate or separately justified autonomous replacement",
            "PA.16 inner-to-outer join",
            "I3/I4 and global matched pressure",
        ],
        "ready": False,
    },
    "oscillatory": {
        "provider": "KokunoCandidateMassBoundPhaseCurlFamily",
        "required_inputs": [
            "provenance-labelled candidate signed mass",
            "D_r h_sigma",
            "D_z h_sigma",
            "beta labels and Q/epsilon/L_s schedule",
            "background values/derivatives",
            "signed mode prototypes",
            "slow partition/covariance target values/derivatives",
        ],
        "current_outputs": [
            "Q-scaled Cartesian velocity by sign",
            "Q-scaled Cartesian velocity by beta",
            "Q-scaled Cartesian velocity total",
        ],
        "missing_public_output": "velocity_osc(x,y,z,t)",
        "next_required_evidence": [
            "bind one serialized provenance-labelled numerical background/mode/partition realization",
            "public xyz,t provider with coefficient units",
            "Agent-4 black-box complete-curl/divergence and phase-mean covariance audit",
        ],
        "ready": False,
    },
    "correction": {
        "provider": "Agent-3 signed mean/radial correction lane",
        "required_inputs": [
            "same-cycle theta/axial requestedStress",
            "candidate-specific finite-head mean factor/debt",
            "independently audited physical oscillatory covariance interface",
        ],
        "required_guards": [
            "H_ref y = DeltaC/epsilon",
            "bounded inverse",
            "unit/budget",
            "spacetime",
            "radial +d_z sigma_1",
            "quadratic damping",
            "joint gain",
        ],
        "ready": False,
    },
    "independent_validator": {
        "provider": "Agent-4 independent validator lane",
        "next_black_box_target": "public velocity_osc(x,y,z,t)",
        "must_not_reuse": [
            "candidate complete-curl implementation as oracle",
            "optimization/training residual",
            "free forcing fitted to residual",
        ],
        "final_gate": "held-out normalized NS momentum max/L2 <= 1e-3 and divergence max/L2 <= 1e-5",
        "ready_for_final_gate": False,
    },
}


def _copy(value: Any) -> Any:
    return json.loads(json.dumps(value))


def _set_pipeline_status(payload: dict[str, Any], stage: str, status: str) -> None:
    for item in payload["pipeline"]:
        if item["stage"] == stage:
            item["ready"] = False
            item["status"] = status
            return
    raise ValueError(f"missing pipeline stage: {stage}")


def build_checkpoint() -> dict[str, Any]:
    payload = _copy(build_parent_checkpoint())
    payload.pop("checkpoint_sha256", None)
    payload["schema"] = SCHEMA
    payload["task_id"] = TASK_ID

    upstream = payload["upstream"]
    upstream["agent1_coupled_c_sibling"] = _copy(AGENT1_COUPLED_C_RECEIPT)
    upstream["agent2_mass_bound_phase_curl_sibling"] = _copy(AGENT2_PHASE_CURL_RECEIPT)
    upstream["previous_agent5_v35_ancestry"] = _copy(PREVIOUS_AGENT5_RECEIPT)

    states = payload["states"]
    # The old selected shared-C path remains historically obstructed in the
    # parent receipt.  The fresh, separately propagated coupled-C realization
    # is a new path and is only known to pass the pointwise necessary screen.
    states["fresh_coupled_C_pa10_pointwise_not_excluded"] = True
    states["fresh_coupled_C_pa16_handoff_allowed"] = False
    states["candidate_mass_bound_phase_curl_composition_ready"] = True
    states["candidate_h_sigma_spatial_derivatives_explicit"] = True
    states["public_provenance_labelled_xyz_t_oscillatory_velocity_ready"] = False

    payload["typed_handoff_contracts"] = _copy(TYPED_HANDOFF_CONTRACTS)

    _set_pipeline_status(
        payload,
        "leading_candidate",
        "Agent 1 #521 has repropagated a fresh factor-32 coupled-C Appendix-B realization and removed the old pointwise PA.10 obstruction. Treat this only as not excluded by the necessary screen: PA.16 remains forbidden until incoming moments, B0/T_sh, join and global matched-pressure evidence are supplied.",
    )
    _set_pipeline_status(
        payload,
        "oscillatory_augmentation",
        "Agent 2 #522 now consumes the candidate signed h_sigma mass directly inside the phase-bound complete-curl stack and requires explicit D_r/D_z h_sigma; Q-scaled by-sign/by-beta/total Cartesian arrays are executable on caller-supplied candidate data. The remaining integration seam is a serialized provenance-labelled background/mode/partition realization plus a public velocity_osc(x,y,z,t) provider, followed by Agent-4 black-box audit.",
    )
    _set_pipeline_status(
        payload,
        "mean_radial_corrections",
        "Agent 3 #514 remains blocked on same-cycle physical requestedStress and candidate-specific finite-head mean debt. Do not enter the signed inverse or finite cycle merely because the autonomous factor is executable; consume the independently audited public oscillatory covariance interface first.",
    )
    _set_pipeline_status(
        payload,
        "independent_validation",
        "Agent 4 #515 has independently closed the signed-mass quadrature seam. The next high-value target is the first public fully numerical velocity_osc(x,y,z,t): independently test complete-curl/divergence and phase-mean covariance rank before any correction or PDE promotion.",
    )

    payload["routing"] = {
        "new_fact": (
            "The previous leading pointwise PA.10 obstruction has been removed on a freshly repropagated autonomous coupled-C path, and the numerical signed h_sigma mass is now wired into the phase/complete-curl implementation with explicit spatial derivatives. Neither fact yet supplies a global leading profile nor a public xyz,t oscillatory provider."
        ),
        "shortest_next_closure": [
            "Agent 1: continue only the fresh #521 coupled-C path through incoming five-moment control and a valid B0/T_sh route; PA.16/I3/I4/global pressure remain fail-closed despite the pointwise necessary-screen pass.",
            "Agent 2: serialize one provenance-labelled numerical background + signed-mode + partition realization on top of #522, preserve explicit D_r/D_z h_sigma and schedule checks, and expose a deterministic public velocity_osc(x,y,z,t) with by-sign/by-beta/total access and explicit units.",
            "Agent 4: consume only that public velocity_osc interface and independently audit complete-curl/divergence plus phase-mean covariance rank under frozen guards; do not infer PASS from #522 composition equivalence.",
            "Agent 3: compute same-cycle theta/axial requestedStress from the audited candidate state, form only candidate-specific finite-head mean debt, then retain the existing DeltaC/epsilon inverse, budget, spacetime, radial, quadratic and gain guards before any finite correction cycle.",
            "Agent 5: once leading + audited public oscillatory + guarded correction coexist in one executable ancestry, instantiate the deterministic composite candidate artifact with velocity/pressure/restricted-forcing APIs, save/load and Python/MATLAB smoke, then freeze for Agent-4 held-out NS validation and same-protocol ST006 comparison.",
        ],
    }

    payload["checkpoint_sha256"] = checkpoint_sha256(payload)
    return payload


def validate_checkpoint(payload: dict[str, Any]) -> None:
    if payload.get("schema") != SCHEMA or payload.get("task_id") != TASK_ID:
        raise ValueError("checkpoint schema/task mismatch")
    if payload.get("checkpoint_sha256") != checkpoint_sha256(payload):
        raise ValueError("checkpoint sha256 mismatch")
    if payload.get("formal_gates") != FORMAL_GATES:
        raise ValueError("formal PDE gates changed")
    if payload.get("autonomous_realization_policy") != AUTONOMOUS_REALIZATION_POLICY:
        raise ValueError("autonomous realization policy changed")

    parent_view = _copy(payload)
    parent_view["schema"] = PARENT_SCHEMA
    parent_view["task_id"] = PARENT_TASK_ID
    parent_view["checkpoint_sha256"] = checkpoint_sha256(parent_view)
    validate_parent_checkpoint(parent_view)

    states = payload["states"]
    if states.get("fresh_coupled_C_pa10_pointwise_not_excluded") is not True:
        raise ValueError("fresh coupled-C PA.10 pass was removed")
    if states.get("fresh_coupled_C_pa16_handoff_allowed") is not False:
        raise ValueError("fresh coupled-C path was promoted into PA.16")
    if states.get("candidate_mass_bound_phase_curl_composition_ready") is not True:
        raise ValueError("mass-bound phase-curl composition was removed")
    if states.get("candidate_h_sigma_spatial_derivatives_explicit") is not True:
        raise ValueError("explicit h_sigma derivative contract was removed")

    for name in (
        "leading_ready",
        "public_provenance_labelled_xyz_t_oscillatory_velocity_ready",
        "oscillatory_ready",
        "genuinely_independent_second_covariance_column_ready",
        "same_cycle_requested_stress_materialized",
        "candidate_numeric_finite_head_mean_debt_materialized",
        "correction_ready",
        "finite_correction_cycle_run",
        "candidate_artifact_instantiated",
        "velocity_export_ready",
        "formal_full_domain_pde_gate_assessed",
        "pde_validated",
    ):
        if states.get(name) is not False:
            raise ValueError(f"fail-closed state promoted: {name}")

    a1 = payload["upstream"].get("agent1_coupled_c_sibling")
    if a1 != AGENT1_COUPLED_C_RECEIPT:
        raise ValueError("Agent 1 coupled-C receipt changed")
    if a1["dedicated_status"] != "success":
        raise ValueError("Agent 1 dedicated receipt is not green")
    if a1["latest_standard_status_observed"] != "in_progress":
        raise ValueError("Agent 1 pending standard receipt was rewritten")
    if not a1["upstream_C_repropagated_before_core"] or a1["posthoc_C_multiplier_applied"]:
        raise ValueError("Agent 1 coupled-C propagation semantics changed")
    if not a1["pointwise_PA10_necessary_screen_passed"]:
        raise ValueError("fresh Agent 1 pointwise PA.10 necessary screen no longer passes")
    if a1["passing_pointwise_screen_is_source_T_sh_certificate"]:
        raise ValueError("pointwise PA.10 pass was laundered into a source T_sh certificate")
    if a1["selected_pa16_handoff_allowed"]:
        raise ValueError("Agent 1 fresh path was prematurely allowed into PA.16")
    if a1["global_pressure_matched"] or a1["global_leading_profile_reconstructed"]:
        raise ValueError("Agent 1 local pass was promoted to global leading readiness")

    a2 = payload["upstream"].get("agent2_mass_bound_phase_curl_sibling")
    if a2 != AGENT2_PHASE_CURL_RECEIPT:
        raise ValueError("Agent 2 phase-curl receipt changed")
    if a2["dedicated_status"] != "success" or a2["standard_status"] != "success":
        raise ValueError("Agent 2 exact-head CI is not green")
    if not a2["candidate_signed_mass_consumed_directly"] or not a2["manual_h_sigma_reentry_removed"]:
        raise ValueError("Agent 2 candidate-mass composition semantics changed")
    if not a2["explicit_D_r_h_sigma_required"] or not a2["explicit_D_z_h_sigma_required"]:
        raise ValueError("Agent 2 h_sigma derivative requirement was weakened")
    if not a2["frozen_local_zero_derivative_mode_requires_exact_zero"]:
        raise ValueError("Agent 2 frozen-zero derivative guard was weakened")
    if a2["public_xyz_t_velocity_osc_provider_ready"]:
        raise ValueError("Agent 2 array composition was promoted to public xyz,t provider")
    if a2["actual_positive_order_background_bound"] or a2["actual_auxiliary_torus_mode_family_bound"]:
        raise ValueError("Agent 2 candidate data were promoted to recovered source inputs")
    if a2["formal_full_domain_pde_gate_assessed"] or a2["pde_validated"]:
        raise ValueError("Agent 2 composition seam was promoted to PDE validation")

    if payload.get("typed_handoff_contracts") != TYPED_HANDOFF_CONTRACTS:
        raise ValueError("typed Kokuno handoff contracts changed")
    if payload["typed_handoff_contracts"]["leading"]["ready"]:
        raise ValueError("leading typed handoff prematurely ready")
    if payload["typed_handoff_contracts"]["oscillatory"]["ready"]:
        raise ValueError("oscillatory typed handoff prematurely ready")
    if payload["typed_handoff_contracts"]["correction"]["ready"]:
        raise ValueError("correction typed handoff prematurely ready")
    if payload["typed_handoff_contracts"]["independent_validator"]["ready_for_final_gate"]:
        raise ValueError("validator prematurely ready for final PDE gate")

    previous = payload["upstream"].get("previous_agent5_v35_ancestry")
    if previous != PREVIOUS_AGENT5_RECEIPT:
        raise ValueError("previous Agent 5 v35 receipt changed")
    if payload["baseline_vs_kokuno"]["kokuno_current_comparable_full_domain_receipt"] is not None:
        raise ValueError("no comparable full-domain Kokuno receipt exists yet")


def write_checkpoint(path: str | Path) -> dict[str, Any]:
    payload = build_checkpoint()
    validate_checkpoint(payload)
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = write_checkpoint(args.output)
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
