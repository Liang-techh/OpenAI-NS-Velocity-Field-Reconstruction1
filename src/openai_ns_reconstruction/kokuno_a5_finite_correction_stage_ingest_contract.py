"""Fail-closed A5 registration of the Agent-3 finite-correction stage contract.

Agent 3 PR #888 defines the typed one-step interface required for a future real
Kokuno correction cycle. It recomputes a complete Navier--Stokes residual on
disjoint held-in/held-out partitions before and after a correction and refuses
incomplete provenance, held-out leakage, protocol drift, residual-as-forcing,
trivial/divergent corrections, or residual growth.

This module registers that downstream seam only. It does not copy Agent-3
mathematics, run a surrogate correction cycle, or invent Agent-4 scientific
evidence. The current candidate remains strict-inner and incomplete, therefore
``correction_ready`` and ``pde_validated`` remain false.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from .kokuno_a5_strict_inner_vorticity_artifact_ingest_contract import (
    deterministic_strict_inner_vorticity_artifact_ingest_contract,
    validate_strict_inner_vorticity_artifact_ingest_contract,
)

SCHEMA = "kokuno-a5-finite-correction-stage-ingest-contract-v1"
TASK = "KOKUNO-A5-FINITE-CORRECTION-STAGE-INGEST-081"
PARENT_A5_PR = 884
PARENT_A5_HEAD = "1d947d4c71812b59b6729befb2be9c46914852cb"
PARENT_A5_SOURCE_BLOB = "c5633e309990c4627c3a2fa858c011c5cb308cb6"

AGENT3_PR = 888
AGENT3_HEAD = "3fd33d3ca20d3ebed7dad872d3264036bdd9cb51"
AGENT3_SOURCE_BLOB = "4fc931fe7703b8b8e05efb41af957b54c9f7a4f9"
AGENT3_WORKFLOW_BLOB = "07da1793048fdfe4369d7327a8dc1d00427819e4"
AGENT3_DEDICATED_RUN = 35530175857
AGENT3_TESTS_RUN = 35530175820
AGENT3_MODULE = "openai_ns_reconstruction.kokuno_finite_correction_stage"
AGENT3_PUBLIC_API = "run_finite_correction_stage(backend,candidate,correction,held_in,held_out)"
AGENT3_PARENT_PR = 882
AGENT3_PARENT_HEAD = "29b385398d0c7636ed8d4049821a534d80728b33"
AGENT3_EXACT_HEAD_CI_CONCLUSION: str | None = None

# Freshness context only: #889 independently audits the saved/reloaded
# strict-inner pressure/forcing-free transport precursor. It is not a finite
# correction-stage audit and not final PDE validation.
LATEST_AGENT4_CONTEXT_PR = 889
LATEST_AGENT4_CONTEXT_HEAD = "574fd5f7d67e5fbf905ef91a483a7e9e8124698d"
LATEST_AGENT4_DEDICATED_RUN = 35530344739
LATEST_AGENT4_TESTS_RUN = 35530344662
AGENT4_DEDICATED_FINITE_CORRECTION_STAGE_AUDIT_PRESENT = False
AGENT4_FINITE_CORRECTION_STAGE_AUDIT_CONCLUSION: str | None = None

FINAL_NORMALIZED_MOMENTUM_GATE = 1.0e-3
FINAL_NORMALIZED_DIVERGENCE_GATE = 1.0e-5

AGENT3_STAGE_PROTOCOL = {
    "public_api": AGENT3_PUBLIC_API,
    "candidate_admission_requires": [
        "complete_ns_defect",
        "corrected_global_leading_join_complete",
        "matched_pressure_gradient_included",
        "restricted_forcing_included",
        "restricted_forcing_preregistered",
        "residual_as_forcing_shortcut_used=false",
    ],
    "correction_provenance_requires": [
        "exact_source_candidate_id_and_sha256",
        "derived_from_complete_ns_defect",
        "fit_partition_id_equals_held_in_partition_id",
        "fit_sample_ids_subset_of_held_in",
        "fit_sample_ids_disjoint_from_held_out",
        "residual_as_forcing_shortcut_used=false",
    ],
    "held_in_held_out_partition_ids_must_differ": True,
    "held_in_held_out_sample_ids_must_be_disjoint": True,
    "backend_recomputes_complete_ns_residual_before_and_after": True,
    "same_residual_protocol_sha_required_before_and_after": True,
    "residual_sample_count_must_match_partition": True,
    "incomplete_ns_residual_rejected": True,
    "reported_momentum_metrics": [
        "normalized_momentum_sample_max",
        "normalized_momentum_grid_l2",
        "normalized_momentum_volume_l2",
    ],
    "reported_divergence_metrics": [
        "normalized_divergence_sample_max",
        "normalized_divergence_grid_l2",
        "normalized_divergence_volume_l2",
    ],
    "reported_correction_metrics": [
        "l2_norm", "max_norm", "divergence_max", "nontriviality"
    ],
    "observed_contraction_factors_recorded_for_held_in_and_held_out": True,
    "trivial_correction_rejected": True,
    "zero_correction_l2_rejected": True,
    "correction_divergence_gate": FINAL_NORMALIZED_DIVERGENCE_GATE,
    "held_in_momentum_growth_rejected": True,
    "held_out_momentum_growth_rejected": True,
    "held_in_divergence_growth_rejected": True,
    "held_out_divergence_growth_rejected": True,
    "caller_supplied_residual_or_gain_allowed": False,
    "caller_supplied_pressure_or_forcing_allowed": False,
    "caller_supplied_scientific_threshold_allowed": False,
    "analytic_kokuno_gain_bindings_complete": False,
    "analytic_kokuno_gain_available": False,
    "observed_gain_is_analytic_kokuno_bound": False,
    "mechanics_only_receipt_is_candidate_residual_evidence": False,
    "mechanics_only_receipt_can_set_pde_validated": False,
    "post_observation_retuning_allowed": False,
}

CURRENT_AGENT3_SCOPED_SOURCE_ADMISSION = {
    "source_kind": "strict-inner pressure/forcing-free transport radial stress",
    "complete_ns_defect": False,
    "pressure_gradient_included": False,
    "restricted_forcing_included": False,
    "scoped_transport_stress_authorized_as_correction_target": False,
    "admitted_to_real_finite_correction_cycle": False,
    "surrogate_defect_substitution_allowed": False,
    "residual_as_forcing_shortcut_allowed": False,
    "real_candidate_finite_correction_cycle_run": False,
    "heldout_normalized_ns_residual_assessed": False,
    "residual_reduction_claimed": False,
    "same_protocol_comparable_to_st006": False,
}


def _digest(payload: Mapping[str, Any]) -> str:
    body = dict(payload)
    body.pop("contract_sha256", None)
    encoded = json.dumps(
        body, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def deterministic_finite_correction_stage_ingest_contract(
    *, exact_head: str | None = None
) -> dict[str, Any]:
    parent = deterministic_strict_inner_vorticity_artifact_ingest_contract(
        exact_head=PARENT_A5_HEAD
    )
    validate_strict_inner_vorticity_artifact_ingest_contract(parent)

    stage_contract_ingest_admitted = bool(
        parent["ingest_status"]["strict_inner_vorticity_artifact_ingest_admitted"]
        and AGENT3_EXACT_HEAD_CI_CONCLUSION == "success"
    )

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "task": TASK,
        "exact_head": exact_head,
        "parent_a5": {
            "pr": PARENT_A5_PR,
            "head": PARENT_A5_HEAD,
            "source_blob_sha": PARENT_A5_SOURCE_BLOB,
        },
        "parent_contract_sha256": parent["contract_sha256"],
        "pipeline_position": {
            "registered_scoped_stage": "finite_correction_stage_contract",
            "upstream_strict_inner_candidate_artifact_registered": True,
            "finite_correction_stage_interface_registered": True,
            "real_finite_correction_cycle_reached": False,
            "full_pipeline_candidate_artifact_stage_reached": False,
            "independent_final_pde_validation_stage_reached": False,
            "reason": (
                "current Kokuno object lacks corrected/global leading join, complete NS defect, "
                "matched pressure, preregistered restricted forcing and authorized correction velocity"
            ),
        },
        "agent3_binding": {
            "pr": AGENT3_PR,
            "head": AGENT3_HEAD,
            "source_blob_sha": AGENT3_SOURCE_BLOB,
            "workflow_blob_sha": AGENT3_WORKFLOW_BLOB,
            "dedicated_run": AGENT3_DEDICATED_RUN,
            "tests_run": AGENT3_TESTS_RUN,
            "module": AGENT3_MODULE,
            "public_api": AGENT3_PUBLIC_API,
            "parent_pr": AGENT3_PARENT_PR,
            "parent_head": AGENT3_PARENT_HEAD,
            "protocol": dict(AGENT3_STAGE_PROTOCOL),
        },
        "agent4_context": {
            "latest_relevant_candidate_audit_pr": LATEST_AGENT4_CONTEXT_PR,
            "latest_relevant_candidate_audit_head": LATEST_AGENT4_CONTEXT_HEAD,
            "latest_relevant_candidate_audit_dedicated_run": LATEST_AGENT4_DEDICATED_RUN,
            "latest_relevant_candidate_audit_tests_run": LATEST_AGENT4_TESTS_RUN,
            "dedicated_finite_correction_stage_audit_present": (
                AGENT4_DEDICATED_FINITE_CORRECTION_STAGE_AUDIT_PRESENT
            ),
            "finite_correction_stage_audit_conclusion": (
                AGENT4_FINITE_CORRECTION_STAGE_AUDIT_CONCLUSION
            ),
            "latest_a4_scope": (
                "saved/reloaded strict-inner pressure/forcing-free transport precursor audit"
            ),
            "latest_a4_scope_is_finite_correction_stage_audit": False,
            "latest_a4_scope_is_final_pde_validation": False,
        },
        "evidence": {
            "parent_vorticity_artifact_ingest_admitted": bool(
                parent["ingest_status"]["strict_inner_vorticity_artifact_ingest_admitted"]
            ),
            "agent3_exact_head_ci_conclusion": AGENT3_EXACT_HEAD_CI_CONCLUSION,
            "agent3_typed_stage_contract_present": True,
            "agent3_mechanics_receipt_is_candidate_residual_evidence": False,
            "agent3_real_candidate_finite_correction_cycle_run": False,
            "agent4_dedicated_finite_correction_stage_audit_present": False,
            "agent4_finite_correction_stage_audit_conclusion": None,
            "finite_correction_stage_contract_ingest_admitted": stage_contract_ingest_admitted,
            "real_correction_scientifically_admitted": False,
        },
        "finite_correction_stage_handoff": {
            "registered": True,
            "status": "registered_unresolved",
            "construction_authority": "Agent3#888",
            "public_stage_api": AGENT3_PUBLIC_API,
            "candidate_and_correction_provenance_typed": True,
            "disjoint_held_in_held_out_enforced": True,
            "held_out_fit_leakage_rejected": True,
            "same_residual_protocol_before_after_enforced": True,
            "complete_ns_residual_required": True,
            "matched_pressure_required": True,
            "restricted_forcing_required": True,
            "restricted_forcing_preregistered_required": True,
            "residual_as_forcing_shortcut_forbidden": True,
            "correction_divergence_gate": FINAL_NORMALIZED_DIVERGENCE_GATE,
            "momentum_growth_rejected_on_held_in_and_held_out": True,
            "divergence_growth_rejected_on_held_in_and_held_out": True,
            "mechanics_only_stage_receipt_is_scientific_candidate_evidence": False,
            "analytic_kokuno_gain_available": False,
            "dedicated_agent4_stage_audit_available": False,
            "usable_for_real_cycle_now": False,
            "usable_as_st006_comparison_now": False,
            "usable_for_final_independent_pde_validation_now": False,
        },
        "current_source_admission": dict(CURRENT_AGENT3_SCOPED_SOURCE_ADMISSION),
        "candidate_api_handoff": {
            **dict(parent["candidate_api_handoff"]),
            "pressure": None,
            "forcing": None,
            "complete_candidate_api_ready": False,
        },
        "artifact_status": {
            **dict(parent["artifact_status"]),
            "upstream_agent3_finite_correction_stage_contract_present": True,
            "real_agent3_correction_velocity_materialized": False,
            "real_finite_correction_cycle_receipt_materialized": False,
            "complete_global_candidate_artifact_materialized": False,
        },
        "ingest_status": {
            **dict(parent["ingest_status"]),
            "typed_finite_correction_stage_contract_registered": True,
            "finite_correction_stage_contract_ingest_admitted": stage_contract_ingest_admitted,
            "real_finite_correction_cycle_admitted": False,
            "full_candidate_artifact_ready": False,
        },
        "stage_state": dict(parent["stage_state"]),
        "baseline_vs_kokuno": dict(parent["baseline_vs_kokuno"]),
        "final_project_gates_unchanged": dict(parent["final_project_gates_unchanged"]),
        "truth_boundary": {
            "finite_correction_stage_interface_registered": True,
            "registration_not_real_correction_cycle": True,
            "registration_not_scientific_admission": True,
            "strict_inner_scope_preserved": True,
            "agent3_ci_promoted": False,
            "agent3_mechanics_receipt_laundered_as_candidate_evidence": False,
            "agent4_stage_audit_invented": False,
            "complete_ns_defect_invented": False,
            "matched_pressure_invented": False,
            "restricted_forcing_invented": False,
            "correction_velocity_invented": False,
            "outer_global_join_invented": False,
            "real_candidate_finite_correction_cycle_run": False,
            "heldout_normalized_ns_residual_assessed": False,
            "residual_reduction_claimed": False,
            "same_protocol_st006_comparison_performed": False,
            "analytic_kokuno_gain_invented": False,
            "heldout_sample_leakage_allowed": False,
            "residual_protocol_drift_allowed": False,
            "free_residual_defined_forcing_allowed": False,
            "threshold_relaxed": False,
            "leading_ready": False,
            "correction_ready": False,
            "velocity_export_ready": False,
            "pde_validated": False,
        },
    }
    payload["contract_sha256"] = _digest(payload)
    return payload


def validate_finite_correction_stage_ingest_contract(payload: Mapping[str, Any]) -> None:
    if not isinstance(payload, Mapping):
        raise ValueError("contract must be a mapping")
    exact_head = payload.get("exact_head")
    if exact_head is not None:
        if not isinstance(exact_head, str) or len(exact_head) != 40:
            raise ValueError("exact_head must be a 40-character commit id")
        if any(c not in "0123456789abcdef" for c in exact_head):
            raise ValueError("exact_head must be lowercase hexadecimal")
    expected = deterministic_finite_correction_stage_ingest_contract(exact_head=exact_head)
    if dict(payload) != expected:
        raise ValueError("finite correction stage ingest contract drifted")
    if payload.get("contract_sha256") != _digest(payload):
        raise ValueError("contract sha256 mismatch")


def write_report(path: str | Path, *, exact_head: str | None = None) -> dict[str, Any]:
    payload = deterministic_finite_correction_stage_ingest_contract(exact_head=exact_head)
    validate_finite_correction_stage_ingest_contract(payload)
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return payload


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--exact-head", default=None)
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)
    payload = write_report(args.output, exact_head=args.exact_head)
    print(json.dumps(payload, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
