"""Agent-5 registration seam for Agent-4 canonical quadrature evidence.

Integration glue only.  This module binds the evidence-bearing successor
validator from Agent 4 PR #903 into the Agent-5 fixed-contract full-NS
validation seam from PR #897.  It does not fabricate a 24/48/96 run, a complete
Kokuno candidate, or a scientific admission.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

SCHEMA_NAME = "kokuno-agent5-canonical-quadrature-evidence-ingest-v1"
TASK_ID = "KOKUNO-A5-CANONICAL-QUADRATURE-EVIDENCE-INGEST-083"

PARENT_A5 = {
    "pr": 897,
    "head": "21ff4a86638a1d8e179e25ed25886a68a34b4e6b",
    "branch": "codex/kokuno-a5-fixed-contract-full-ns-validation-ingest-082",
    "source_path": "src/openai_ns_reconstruction/kokuno_a5_fixed_contract_full_ns_validation_ingest_contract.py",
    "source_blob": "fcf2f8945592782b2dd4008007a970d5d4fdf83d",
}

AGENT4_CANONICAL = {
    "pr": 903,
    "head": "065810d94a37ade8d81acd2c7f622c1ff8cc6bdd",
    "source_path": "src/openai_ns_reconstruction/kokuno_a4_canonical_quadrature_evidence.py",
    "source_blob": "4a14a2e2bd6c06280d992802414109444d41e84b",
    "workflow_path": ".github/workflows/kokuno-agent4-canonical-quadrature-evidence.yml",
    "workflow_blob": "d8dc0e4bc971e94f5b886cc57ffc3a29bade58e6",
    "parent_a4_pr": 896,
    "parent_a4_head": "e8a7712515f75bb8a0a5d86b9a6177ab9447f2b1",
    "cr002_audit_pr": 898,
    "cr002_audit_head": "3b34375303c7237eb8a300a6053aab72e7826bf6",
    "observed_ci": {
        "dedicated": {"run_id": 35536767684, "status": "queued", "conclusion": None},
        "repository_tests": {"run_id": 35536767692, "status": "queued", "conclusion": None},
    },
}

QUADRATURE_ORDERS = [24, 48, 96]
DERIVATIVE_STEPS = [0.02, 0.01, 0.005]
VALIDATION_TIMES = [0.25, 0.3125, 0.4375, 0.5625, 0.6875, 0.75]
QUADRATURE_STABILITY_GATE = 5.0e-2
FINAL_GATE = {
    "momentum_sampled_max": 1.0e-3,
    "momentum_volume_l2": 1.0e-3,
    "divergence_sampled_max": 1.0e-5,
    "divergence_volume_l2": 1.0e-5,
}
ST006_BASELINE = {
    "momentum_sampled_max": 0.1082289305112118,
    "momentum_volume_l2": 0.10758432876230622,
}
READINESS = {
    "leading_ready": False,
    "oscillatory_ready": True,
    "correction_ready": False,
    "velocity_export_ready": False,
    "pde_validated": False,
}


def _canonical(payload: Mapping[str, Any]) -> bytes:
    return json.dumps(dict(payload), sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")


def _sha256(payload: Mapping[str, Any]) -> str:
    return hashlib.sha256(_canonical(payload)).hexdigest()


def _without_sha(contract: Mapping[str, Any]) -> dict[str, Any]:
    payload = copy.deepcopy(dict(contract))
    payload.pop("contract_sha256", None)
    return payload


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def build_contract(*, exact_head: str) -> dict[str, Any]:
    contract: dict[str, Any] = {
        "schema_name": SCHEMA_NAME,
        "task_id": TASK_ID,
        "exact_head": str(exact_head),
        "parent_a5": copy.deepcopy(PARENT_A5),
        "agent4_canonical_quadrature_handoff": {
            **copy.deepcopy(AGENT4_CANONICAL),
            "api": "validate_complete_candidate_with_canonical_quadrature",
            "receipt_schema": "canonical-cr001-quadrature-receipt-v1",
            "candidate_sha256_required": True,
            "stage_identity_required": True,
            "physical_contract_sha256_required": True,
            "validator_operator_sha256_required": True,
            "heldout_protocol_sha256_required": True,
            "derivative_steps": copy.deepcopy(DERIVATIVE_STEPS),
            "validation_times": copy.deepcopy(VALIDATION_TIMES),
            "quadrature_orders_per_axis": copy.deepcopy(QUADRATURE_ORDERS),
            "per_time_per_order_momentum_volume_l2_required": True,
            "per_time_per_order_divergence_volume_l2_required": True,
            "per_time_per_order_kinetic_energy_required": True,
            "receipt_sha256_required": True,
            "independent_of_training_and_held_in_selection_required": True,
            "bare_boolean_authorizes_pde_promotion": False,
            "caller_sample_weights_are_canonical_quadrature": False,
            "quadrature_stability_gate_48_to_96": QUADRATURE_STABILITY_GATE,
            "real_canonical_quadrature_run_completed": False,
            "real_receipt_available": False,
            "scientific_admission": False,
        },
        "admission_seam": {
            "candidate_identity_must_match": True,
            "stage_identity_must_match": True,
            "physical_contract_identity_must_match": True,
            "validator_operator_identity_must_match": True,
            "heldout_protocol_identity_must_match": True,
            "derivative_ladder_must_match": True,
            "validation_times_must_match": True,
            "canonical_quadrature_ladder_must_match": True,
            "typed_receipt_required": True,
            "receipt_digest_must_verify": True,
            "sampled_rms_cannot_be_relabelled_canonical_volume_l2": True,
            "monte_carlo_weights_cannot_replace_canonical_quadrature": True,
            "residual_defined_forcing_forbidden": True,
            "exact_head_ci_success_required_before_scientific_admission": True,
            "current_evidence_admitted": False,
            "global_join_present": False,
            "matched_pressure_present": False,
            "restricted_forcing_present": False,
            "restricted_forcing_preregistered": False,
            "real_correction_velocity_present": False,
            "complete_candidate_api_ready": False,
            "current_candidate_eligible_for_full_ns_validation": False,
            "same_protocol_st006_comparison_available_now": False,
        },
        "readiness": copy.deepcopy(READINESS),
        "final_gate": copy.deepcopy(FINAL_GATE),
        "baseline": {"st006": copy.deepcopy(ST006_BASELINE)},
        "truth_boundary": {
            "kokuno_reconstruction_is_paper_exact": False,
            "queued_ci_is_not_pass": True,
            "typed_receipt_is_not_a_real_quadrature_run": True,
            "sampled_rms_is_not_canonical_volume_l2": True,
            "visual_success_is_not_pde_validation": True,
            "residual_defined_forcing_forbidden": True,
            "threshold_relaxation_allowed": False,
        },
    }
    contract["contract_sha256"] = _sha256(contract)
    validate_contract(contract)
    return contract


def validate_contract(contract: Mapping[str, Any]) -> None:
    value = dict(contract)
    _require(value.get("schema_name") == SCHEMA_NAME, "schema drift")
    _require(value.get("task_id") == TASK_ID, "task drift")
    exact_head = str(value.get("exact_head", ""))
    _require(len(exact_head) == 40 and all(c in "0123456789abcdef" for c in exact_head), "exact head must be a lowercase 40-hex SHA")
    _require(value.get("parent_a5") == PARENT_A5, "parent A5 provenance drift")

    a4 = value.get("agent4_canonical_quadrature_handoff", {})
    for key in AGENT4_CANONICAL:
        _require(a4.get(key) == AGENT4_CANONICAL[key], f"Agent-4 provenance drift: {key}")
    _require(a4.get("api") == "validate_complete_candidate_with_canonical_quadrature", "Agent-4 API drift")
    _require(a4.get("receipt_schema") == "canonical-cr001-quadrature-receipt-v1", "receipt schema drift")
    for key in (
        "candidate_sha256_required", "stage_identity_required", "physical_contract_sha256_required",
        "validator_operator_sha256_required", "heldout_protocol_sha256_required",
        "per_time_per_order_momentum_volume_l2_required", "per_time_per_order_divergence_volume_l2_required",
        "per_time_per_order_kinetic_energy_required", "receipt_sha256_required",
        "independent_of_training_and_held_in_selection_required",
    ):
        _require(a4.get(key) is True, f"canonical evidence requirement weakened: {key}")
    _require(a4.get("derivative_steps") == DERIVATIVE_STEPS, "derivative ladder drift")
    _require(a4.get("validation_times") == VALIDATION_TIMES, "validation times drift")
    _require(a4.get("quadrature_orders_per_axis") == QUADRATURE_ORDERS, "quadrature ladder drift")
    _require(a4.get("quadrature_stability_gate_48_to_96") == QUADRATURE_STABILITY_GATE, "quadrature stability gate drift")
    for key in ("bare_boolean_authorizes_pde_promotion", "caller_sample_weights_are_canonical_quadrature", "real_canonical_quadrature_run_completed", "real_receipt_available", "scientific_admission"):
        _require(a4.get(key) is False, f"false admission boundary drift: {key}")

    seam = value.get("admission_seam", {})
    for key in (
        "candidate_identity_must_match", "stage_identity_must_match", "physical_contract_identity_must_match",
        "validator_operator_identity_must_match", "heldout_protocol_identity_must_match",
        "derivative_ladder_must_match", "validation_times_must_match", "canonical_quadrature_ladder_must_match",
        "typed_receipt_required", "receipt_digest_must_verify",
        "sampled_rms_cannot_be_relabelled_canonical_volume_l2",
        "monte_carlo_weights_cannot_replace_canonical_quadrature", "residual_defined_forcing_forbidden",
        "exact_head_ci_success_required_before_scientific_admission",
    ):
        _require(seam.get(key) is True, f"admission firewall weakened: {key}")
    for key in (
        "current_evidence_admitted", "global_join_present", "matched_pressure_present", "restricted_forcing_present",
        "restricted_forcing_preregistered", "real_correction_velocity_present", "complete_candidate_api_ready",
        "current_candidate_eligible_for_full_ns_validation", "same_protocol_st006_comparison_available_now",
    ):
        _require(seam.get(key) is False, f"missing dependency laundered as present: {key}")

    _require(value.get("readiness") == READINESS, "readiness drift")
    _require(value.get("final_gate") == FINAL_GATE, "final 1e-3/1e-5 gate drift")
    _require(value.get("baseline") == {"st006": ST006_BASELINE}, "ST006 baseline drift")
    truth = value.get("truth_boundary", {})
    for key in ("queued_ci_is_not_pass", "typed_receipt_is_not_a_real_quadrature_run", "sampled_rms_is_not_canonical_volume_l2", "visual_success_is_not_pde_validation", "residual_defined_forcing_forbidden"):
        _require(truth.get(key) is True, f"truth boundary weakened: {key}")
    _require(truth.get("kokuno_reconstruction_is_paper_exact") is False, "Kokuno reconstruction cannot be relabeled paper-exact")
    _require(truth.get("threshold_relaxation_allowed") is False, "threshold relaxation forbidden")

    supplied_sha = value.get("contract_sha256")
    _require(isinstance(supplied_sha, str) and len(supplied_sha) == 64, "missing contract SHA")
    _require(supplied_sha == _sha256(_without_sha(value)), "contract SHA mismatch")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--exact-head", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    contract = build_contract(exact_head=args.exact_head)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(contract, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
