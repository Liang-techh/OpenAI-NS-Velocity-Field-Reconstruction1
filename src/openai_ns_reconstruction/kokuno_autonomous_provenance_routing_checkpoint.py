"""Agent-5 routing checkpoint for source-compatible autonomous realization.

The active project does not require recovery of unpublished Kokuno coefficients,
but it does require strict provenance.  This checkpoint makes that distinction
machine-readable: a repository-autonomous, source-compatible numerical choice
may enter the *candidate* pipeline when labelled as such, while it may never be
promoted to recovered source data or used to discharge an opaque formal-source
numeric witness.

The immediate motivation is the current two-sided frontier:

* Agent 2 has a CI-green autonomous signed-rectangle geometry witness that
  preserves the displayed source structure but is not the hidden source tuple.
* Agent 3 pins the formal finite-head mean-factor dependency and proves that the
  selected Mathlib bump / actual cycle-state stress are not numerically exported.

Treating either fact as a permanent requirement for hidden-data recovery would
stall the user-requested Kokuno-derived candidate.  Treating autonomous values
as source-exact would be a provenance error.  This module encodes the shortest
middle route without changing any PDE gate.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

SCHEMA = "kokuno-agent5-autonomous-provenance-routing-checkpoint-v33"
TASK_ID = "KOKUNO-A5-AUTONOMOUS-PROVENANCE-ROUTING-033"

SOURCE_PROVENANCE = {
    "repository": "KokunoYumeto/yang-mills-interacting-workbench",
    "commit": "143f6773feb424ad9ed3a8d116653200f20346b7",
    "document": "navier-stokes/navier_stokes_workbench.tex",
    "corrected_reader_date": "2026-09-09",
    "zenodo_record": "22678406",
    "paper_exact": False,
}

EVIDENCE_CLASSES = (
    "source_displayed",
    "formal_structure",
    "repository_autonomous_source_compatible",
    "candidate_derived",
    "independent_validation",
)

AUTONOMOUS_REALIZATION_POLICY = {
    "candidate_may_use_repository_autonomous_source_compatible_numerics": True,
    "autonomous_values_may_be_called_recovered_source_values": False,
    "autonomous_values_may_discharge_opaque_formal_numeric_witness": False,
    "candidate_specific_numeric_mean_factor_may_be_defined": True,
    "candidate_specific_mean_factor_must_be_labelled_autonomous": True,
    "formal_actual_missingWeight_claim_must_remain_false_without_machine_linked_source_bump": True,
    "candidate_state_requestedStress_must_be_computed_from_same_candidate_cycle_state": True,
    "independent_validator_must_not_reuse_candidate_construction_or_training_loss": True,
    "forcing_policy": "fixed/restricted only; residual-defined free forcing forbidden",
}

AGENT1_RECEIPT = {
    "pr": 500,
    "head": "c360b7a130e007f160bb2680415cc8f020db2e04",
    "standard_run": 35403080310,
    "standard_status_at_freeze": "in_progress",
    "same_scale_T_sh_screen_implemented": True,
    "source_B0_analytic_bound_proved": False,
    "source_T_sh_lower_bound_verified": False,
    "inner_to_outer_join_completed": False,
    "global_pressure_matched": False,
    "global_leading_profile_reconstructed": False,
    "consumed_in_executable_ancestry": False,
}

AGENT1_LAST_FULLY_GREEN_RECEIPT = {
    "pr": 491,
    "head": "f44182fe310cc87a010f8ae96a7c66d173a3229d",
    "standard_run": 35397685224,
    "standard_status": "success",
    "selected_source_rescaled_prefix_moments_executable": True,
    "PA16_input_tuple_formable": True,
}

AGENT2_RECEIPT = {
    "pr": 501,
    "head": "58b95eafb3429685f5d7e83c22fa7b0a8f9eeba8",
    "dedicated_run": 35402140447,
    "standard_run": 35402140319,
    "dedicated_status": "success",
    "standard_status": "success",
    "evidence_class": "repository_autonomous_source_compatible",
    "autonomous_signed_rectangle_geometry_ready": True,
    "source_signed_rectangle_structure_preserved": True,
    "source_rectangle_separation_guards_certified": True,
    "source_band_Q_epsilon_Ls_schedule_executable": True,
    "source_rectangle_centers_recovered": False,
    "source_rectangle_radius_r0_recovered": False,
    "source_actual_rectangle_labels_instantiated": False,
    "source_actual_partition_labels_instantiated": False,
    "actual_source_h_sigma_pulse_integrals_bound": False,
    "actual_positive_order_background_bound": False,
    "actual_auxiliary_torus_mode_family_bound": False,
    "public_xyz_t_velocity_correction_materialized": False,
    "consumed_in_executable_ancestry": False,
}

AGENT3_RECEIPT = {
    "pr": 503,
    "head": "d4ce470abe655378746e8e3131e2891280476619",
    "dedicated_run": 35402618032,
    "standard_run": 35402618004,
    "dedicated_status": "success",
    "standard_status": "success",
    "evidence_class": "formal_structure",
    "source_algebra_pinned": True,
    "physical_scale_algebra_executable_given_coordinate_q": True,
    "numeric_bump_profile_export_present": False,
    "missing_weight_values_materialized": False,
    "requested_stress_actual_state_values_materialized": False,
    "finite_head_mean_debt_materialized": False,
    "arbitrary_python_bump_may_be_promoted_to_actual": False,
    "real_candidate_defect_consumed": False,
    "signed_mean_inverse_input_ready": False,
    "finite_correction_cycle_rerun_allowed": False,
    "consumed_in_executable_ancestry": True,
}

AGENT4_RECEIPT = {
    "pr": 494,
    "head": "31f7510545488818b3bcf979c5cac3cb5d826a49",
    "dedicated_run": 35398848549,
    "standard_run": 35398848498,
    "dedicated_status": "success",
    "standard_status": "success",
    "evidence_class": "independent_validation",
    "displayed_source_phase_frame_independently_audited": True,
    "local_structural_preflight_passed": True,
    "value_relative_max": 6.661338147750939e-16,
    "fd_finest_relative_rms": 4.3321017055747775e-12,
    "wrong_axial_normalization_relative_rms": 2.0965112461866924,
    "wrong_tilt_mutation_relative_rms": 1.5780867798741534,
    "formal_full_domain_pde_gate_assessed": False,
    "pde_validated": False,
    "consumed_in_executable_ancestry": False,
}

PREVIOUS_AGENT5_RECEIPT = {
    "pr": 495,
    "head": "cdf586f8b173fd99a515ccba57d8038d860c6ea5",
    "schema": "kokuno-agent5-source-phase-mean-routing-checkpoint-v32",
    "dedicated_run": 35399233925,
    "standard_run": 35399233974,
    "dedicated_status": "success",
    "standard_status": "success",
    "consumed_in_executable_ancestry": False,
}

ST006_BASELINE = {
    "candidate": "ST006",
    "held_out_seed": 9172801,
    "points": 4096,
    "times": 6,
    "finest_spatial_step": 0.005,
    "momentum_sampled_max": 0.1082289305112118,
    "momentum_volume_l2": 0.10758432876230622,
    "pde_validated": False,
}

FORMAL_GATES = {
    "held_out_normalized_momentum_max": 1.0e-3,
    "held_out_normalized_momentum_l2": 1.0e-3,
    "held_out_divergence_max": 1.0e-5,
    "held_out_divergence_l2": 1.0e-5,
    "changed_this_round": False,
}

CANDIDATE_ARTIFACT_CONTRACT = {
    "schema_reserved": "kokuno-derived-3d-candidate-v1",
    "instantiated": False,
    "required_component_provenance": [
        "evidence_class",
        "source repository/version when applicable",
        "repository-autonomous parameters when applicable",
        "truth-boundary flags",
        "component payload sha256",
    ],
    "required_public_api": [
        "velocity(x,y,z,t)->[u,v,w]",
        "pressure(x,y,z,t)",
        "restricted_forcing(x,y,z,t)",
        "deterministic save/load",
    ],
    "required_exports_after_instantiation": ["Python smoke", "MATLAB smoke"],
}


def _copy(value: Any) -> Any:
    return json.loads(json.dumps(value))


def checkpoint_sha256(payload: dict[str, Any]) -> str:
    clean = _copy(payload)
    clean.pop("checkpoint_sha256", None)
    encoded = json.dumps(clean, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode()
    return hashlib.sha256(encoded).hexdigest()


def validate_component_provenance(component: dict[str, Any]) -> None:
    evidence_class = component.get("evidence_class")
    if evidence_class not in EVIDENCE_CLASSES:
        raise ValueError("unknown evidence_class")
    source_exact = bool(component.get("source_exact", False))
    recovered_source = bool(component.get("recovered_source_numeric", False))
    discharges_formal = bool(component.get("discharges_opaque_formal_numeric_witness", False))
    if evidence_class == "repository_autonomous_source_compatible":
        if source_exact or recovered_source:
            raise ValueError("autonomous source-compatible data cannot be promoted to recovered source truth")
        if discharges_formal:
            raise ValueError("autonomous source-compatible data cannot discharge an opaque formal numeric witness")


def build_checkpoint() -> dict[str, Any]:
    states = {
        "leading_ready": False,
        "displayed_source_phase_frame_independently_audited": True,
        "autonomous_signed_rectangle_geometry_ready": True,
        "formal_mean_factor_dependency_pinned": True,
        "source_exact_numeric_realization_required_for_candidate": False,
        "candidate_autonomous_realization_route_open": True,
        "candidate_autonomous_mean_profile_allowed": True,
        "candidate_autonomous_mean_profile_selected": False,
        "candidate_state_requestedStress_materialized": False,
        "candidate_numeric_finite_head_mean_debt_materialized": False,
        "actual_source_numeric_bump_profile_recovered": False,
        "actual_source_missing_weight_materialized": False,
        "actual_source_requested_stress_materialized": False,
        "actual_source_finite_head_mean_debt_materialized": False,
        "actual_positive_order_background_bound": False,
        "actual_source_h_sigma_pulse_integrals_bound": False,
        "actual_auxiliary_torus_mode_family_bound": False,
        "public_provenance_labelled_xyz_t_oscillatory_velocity_ready": False,
        "oscillatory_ready": False,
        "genuinely_independent_second_covariance_column_ready": False,
        "signed_mean_inverse_input_ready": False,
        "real_candidate_defect_consumed": False,
        "correction_ready": False,
        "finite_correction_cycle_run": False,
        "candidate_artifact_instantiated": False,
        "velocity_export_ready": False,
        "formal_full_domain_pde_gate_assessed": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
    }
    payload = {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "source_provenance": _copy(SOURCE_PROVENANCE),
        "evidence_classes": list(EVIDENCE_CLASSES),
        "autonomous_realization_policy": _copy(AUTONOMOUS_REALIZATION_POLICY),
        "upstream": {
            "agent1_latest_sibling": _copy(AGENT1_RECEIPT),
            "agent1_last_fully_green_sibling": _copy(AGENT1_LAST_FULLY_GREEN_RECEIPT),
            "agent2_sibling": _copy(AGENT2_RECEIPT),
            "agent3_ancestry": _copy(AGENT3_RECEIPT),
            "agent4_sibling": _copy(AGENT4_RECEIPT),
            "previous_agent5_sibling": _copy(PREVIOUS_AGENT5_RECEIPT),
        },
        "states": states,
        "pipeline": [
            {
                "stage": "source_profile_ingest",
                "ready": False,
                "status": "displayed source structure plus autonomous signed rectangle witness are available, but the numerical positive-order/background and h_sigma pulse realization are not yet bound into one provenance-labelled family",
            },
            {
                "stage": "leading_candidate",
                "ready": False,
                "status": "Agent 1 latest same-scale T_sh screen exists but exact-head CI is still pending at this freeze; global join and matched pressure remain absent",
            },
            {
                "stage": "oscillatory_augmentation",
                "ready": False,
                "status": "source-exact hidden rectangle recovery is no longer treated as a candidate prerequisite; next bind explicitly autonomous/source-compatible rectangle, pulse, background and mode numerics into the audited complete-curl API",
            },
            {
                "stage": "mean_radial_corrections",
                "ready": False,
                "status": "formal actual missingWeight remains opaque; next candidate route may select an autonomous numeric bump/profile only if it is labelled candidate-specific and uses requestedStress from the same candidate cycle state",
            },
            {
                "stage": "finite_correction_cycle",
                "ready": False,
                "status": "blocked until candidate-specific numeric mean debt and independently audited physical covariance family coexist and all inherited inverse/budget/spacetime/radial/quadratic guards pass",
            },
            {
                "stage": "candidate_artifact",
                "ready": False,
                "status": "reserved typed contract only; leading, oscillatory and correction components do not yet coexist in one executable ancestry",
            },
            {
                "stage": "independent_validation",
                "ready": False,
                "status": "Agent 4 local phase/frame audit is green; formal held-out NS gate waits for the frozen global candidate",
            },
            {
                "stage": "report_and_export",
                "ready": False,
                "status": "Python/MATLAB smoke remains blocked until deterministic candidate artifact instantiation",
            },
        ],
        "routing": {
            "new_fact": "Hidden source numerics and candidate readiness are now separate axes: autonomous source-compatible values may advance the labelled candidate, but may not be called recovered source truth or discharge the formal theorem's opaque numeric witness.",
            "shortest_next_closure": [
                "Agent 2: consume #501 autonomous signed rectangle geometry together with explicitly provenance-labelled numerical positive-order/background, h_sigma pulses and auxiliary modes; emit Q-scaled by-sign/by-beta/total velocity_osc(x,y,z,t) and keep all source-recovery flags false.",
                "Agent 3: for the candidate route, select and serialize one repository-autonomous source-compatible numeric bump/profile and compute requestedStress from the same candidate cycle state; materialize candidate-specific finite-head mean debt, while keeping the formal theorem's actual missingWeight/materialization flags false.",
                "Agent 4: independently black-box audit the first fully numerical provenance-labelled oscillatory family and its phase-mean covariance response; do not inherit manufactured/source-compatible construction tests as independent validation.",
                "Agent 1: finish the exact-head same-scale T_sh result, then PA.16 join, I3/I4 and global matched pressure; a failed same-scale screen must remain a failed screen rather than be bypassed by relabelling.",
                "Agent 5: once leading + oscillatory + correction coexist, instantiate the deterministic candidate artifact, save/load, Python/MATLAB smoke, then freeze it for Agent 4 held-out normalized momentum/divergence validation and same-protocol ST006 comparison.",
            ],
        },
        "candidate_artifact_contract": _copy(CANDIDATE_ARTIFACT_CONTRACT),
        "baseline_vs_kokuno": {
            "st006": _copy(ST006_BASELINE),
            "kokuno_current_comparable_full_domain_receipt": None,
            "comparison_status": "not_comparable_yet",
        },
        "formal_gates": _copy(FORMAL_GATES),
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
    policy = payload.get("autonomous_realization_policy", {})
    if policy != AUTONOMOUS_REALIZATION_POLICY:
        raise ValueError("autonomous realization policy changed")
    states = payload.get("states", {})
    required_false = (
        "leading_ready",
        "candidate_autonomous_mean_profile_selected",
        "candidate_state_requestedStress_materialized",
        "candidate_numeric_finite_head_mean_debt_materialized",
        "actual_source_numeric_bump_profile_recovered",
        "actual_source_missing_weight_materialized",
        "actual_source_requested_stress_materialized",
        "actual_source_finite_head_mean_debt_materialized",
        "actual_positive_order_background_bound",
        "actual_source_h_sigma_pulse_integrals_bound",
        "actual_auxiliary_torus_mode_family_bound",
        "public_provenance_labelled_xyz_t_oscillatory_velocity_ready",
        "oscillatory_ready",
        "genuinely_independent_second_covariance_column_ready",
        "signed_mean_inverse_input_ready",
        "real_candidate_defect_consumed",
        "correction_ready",
        "finite_correction_cycle_run",
        "candidate_artifact_instantiated",
        "velocity_export_ready",
        "formal_full_domain_pde_gate_assessed",
        "pde_validated",
    )
    if any(states.get(name) is not False for name in required_false):
        raise ValueError("fail-closed state promoted")
    required_true = (
        "displayed_source_phase_frame_independently_audited",
        "autonomous_signed_rectangle_geometry_ready",
        "formal_mean_factor_dependency_pinned",
        "candidate_autonomous_realization_route_open",
        "candidate_autonomous_mean_profile_allowed",
    )
    if any(states.get(name) is not True for name in required_true):
        raise ValueError("established routing fact was removed")
    if states.get("source_exact_numeric_realization_required_for_candidate") is not False:
        raise ValueError("hidden source recovery was reintroduced as a candidate prerequisite")
    a2 = payload["upstream"]["agent2_sibling"]
    validate_component_provenance(
        {
            "evidence_class": a2["evidence_class"],
            "source_exact": a2["source_rectangle_centers_recovered"],
            "recovered_source_numeric": a2["source_rectangle_radius_r0_recovered"],
            "discharges_opaque_formal_numeric_witness": False,
        }
    )
    a3 = payload["upstream"]["agent3_ancestry"]
    if a3["numeric_bump_profile_export_present"] or a3["finite_head_mean_debt_materialized"]:
        raise ValueError("opaque formal mean factor was promoted numerically")
    if a3["arbitrary_python_bump_may_be_promoted_to_actual"]:
        raise ValueError("autonomous bump was laundered into formal source truth")
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
