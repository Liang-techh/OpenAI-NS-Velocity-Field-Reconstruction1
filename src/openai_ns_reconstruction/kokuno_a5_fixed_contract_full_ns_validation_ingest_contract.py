"""Agent-5 typed ingest contract for fixed-contract correction -> full-NS validation.

This module is integration glue only.  It binds the Agent-3 fixed-physical-
contract finite-correction surface to the Agent-4 implementation-distinct
black-box full-NS validator without inventing the scientific surfaces that are
still missing from the current Kokuno route.

In particular, registration here does *not* mean that a corrected/global
leading velocity, matched pressure, preregistered restricted forcing, real
correction velocity, canonical volume-L2 assessment, or PDE validation exists.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

SCHEMA_NAME = "kokuno-agent5-fixed-contract-full-ns-validation-ingest-v1"
TASK_ID = "KOKUNO-A5-FIXED-CONTRACT-FULL-NS-VALIDATION-INGEST-082"

PARENT_A5 = {
    "pr": 890,
    "head": "e833441279ed1f1df28c83b91b6226257e478e3f",
    "branch": "codex/kokuno-a5-finite-correction-stage-ingest-081",
}
AGENT3 = {
    "pr": 895,
    "head": "d305b64edcdde4f16dd8a48908f53fa23c2b9baf",
    "source_path": "src/openai_ns_reconstruction/kokuno_finite_correction_fixed_physical_contract.py",
    "source_blob": "0cd41dfcc998064c37075b6086a246e93831ebfd",
    "workflow_path": ".github/workflows/kokuno-agent3-finite-correction-physical-contract.yml",
    "observed_ci": {
        "repository_tests": {"run_id": 35533175829, "status": "queued", "conclusion": None},
        "dedicated": {"run_id": 35533175825, "status": "queued", "conclusion": None},
    },
}
AGENT4 = {
    "pr": 896,
    "head": "e8a7712515f75bb8a0a5d86b9a6177ab9447f2b1",
    "source_path": "src/openai_ns_reconstruction/kokuno_a4_blackbox_full_ns_validator.py",
    "source_blob": "a9e46daf60d909a63261dd10821e1a08101dd872",
    "workflow_path": ".github/workflows/kokuno-agent4-blackbox-full-ns-validator.yml",
    "observed_ci": {
        "repository_tests": {"run_id": 35533703074, "status": "queued", "conclusion": None},
        "dedicated": {"run_id": 35533703072, "status": "queued", "conclusion": None},
    },
}

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
QUADRATURE_LADDER = [24, 48, 96]
DERIVATIVE_STEP_LADDER = [2.0e-2, 1.0e-2, 5.0e-3]


def _canonical(payload: Mapping[str, Any]) -> bytes:
    return json.dumps(
        dict(payload), sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")


def _sha256(payload: Mapping[str, Any]) -> str:
    return hashlib.sha256(_canonical(payload)).hexdigest()


def _contract_payload_without_sha(contract: Mapping[str, Any]) -> dict[str, Any]:
    payload = copy.deepcopy(dict(contract))
    payload.pop("contract_sha256", None)
    return payload


def build_contract(*, exact_head: str) -> dict[str, Any]:
    contract: dict[str, Any] = {
        "schema_name": SCHEMA_NAME,
        "task_id": TASK_ID,
        "exact_head": str(exact_head),
        "parent_a5": copy.deepcopy(PARENT_A5),
        "agent3_fixed_physical_contract_handoff": {
            **copy.deepcopy(AGENT3),
            "api": "run_fixed_physical_contract_finite_correction_stage",
            "physical_contract_sha256_required": True,
            "before_after_same_physical_contract_sha256": True,
            "velocity_only_correction": True,
            "pressure_frozen_through_correction": True,
            "forcing_frozen_through_correction": True,
            "candidate_source_identity_required": True,
            "held_in_held_out_disjoint_required": True,
            "residual_defined_forcing_forbidden": True,
            "joint_pressure_or_forcing_updates_require_separate_preregistered_stage": True,
            "scientific_admission": False,
        },
        "agent4_blackbox_full_ns_validator_handoff": {
            **copy.deepcopy(AGENT4),
            "public_candidate_surfaces": [
                "velocity",
                "pressure",
                "forcing",
                "validation_metadata",
            ],
            "construction_derivatives_allowed": False,
            "construction_defect_tensors_allowed": False,
            "training_state_allowed": False,
            "operator": {
                "space": "centered_cartesian_fourth_order_five_point",
                "time": "centered_or_one_sided_fourth_order",
                "viscosity": 0.01,
                "derivative_step_ladder": copy.deepcopy(DERIVATIVE_STEP_LADDER),
            },
            "heldout": {
                "seed": 914027,
                "sample_count": 4096,
                "validation_times": [0.25, 0.3125, 0.4375, 0.5625, 0.6875, 0.75],
                "evaluation_box": [[-2.0, 2.0], [-2.0, 2.0], [-2.0, 2.0]],
                "time_interval": [0.25, 0.75],
            },
            "canonical_quadrature_ladder": copy.deepcopy(QUADRATURE_LADDER),
            "quadrature_ladder_assessed_required_for_pde_validated": True,
            "physical_contract_sha256_required": True,
            "complete_candidate_admission_required": True,
            "scientific_admission": False,
        },
        "integration_seam": {
            "correction_and_validator_physical_contract_sha256_must_match": True,
            "before_after_physical_contract_identity_locked": True,
            "pressure_identity_immutable_through_velocity_only_correction": True,
            "forcing_identity_immutable_through_velocity_only_correction": True,
            "residual_defined_forcing_forbidden": True,
            "global_join_present": False,
            "matched_pressure_present": False,
            "restricted_forcing_present": False,
            "restricted_forcing_preregistered": False,
            "real_correction_velocity_present": False,
            "complete_candidate_api_ready": False,
            "current_candidate_eligible_for_full_ns_validation": False,
            "same_protocol_st006_comparison_available_now": False,
        },
        "readiness": {
            "leading_ready": False,
            "oscillatory_ready": True,
            "correction_ready": False,
            "velocity_export_ready": False,
            "pde_validated": False,
        },
        "final_gate": copy.deepcopy(FINAL_GATE),
        "baseline": {"st006": copy.deepcopy(ST006_BASELINE)},
        "truth_boundary": {
            "kokuno_reconstruction_is_paper_exact": False,
            "queued_ci_is_not_pass": True,
            "mechanics_or_manufactured_checks_are_not_scientific_admission": True,
            "sampled_rms_is_not_canonical_volume_l2": True,
            "monte_carlo_volume_l2_does_not_replace_canonical_quadrature_ladder": True,
            "visual_success_is_not_pde_validation": True,
            "residual_defined_forcing_forbidden": True,
        },
    }
    contract["contract_sha256"] = _sha256(contract)
    validate_contract(contract)
    return contract


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def validate_contract(contract: Mapping[str, Any]) -> None:
    value = dict(contract)
    _require(value.get("schema_name") == SCHEMA_NAME, "schema drift")
    _require(value.get("task_id") == TASK_ID, "task drift")
    exact_head = str(value.get("exact_head", ""))
    _require(len(exact_head) == 40 and all(c in "0123456789abcdef" for c in exact_head), "exact head must be a lowercase 40-hex SHA")
    _require(value.get("parent_a5") == PARENT_A5, "parent A5 provenance drift")

    a3 = value.get("agent3_fixed_physical_contract_handoff", {})
    for key in ("pr", "head", "source_path", "source_blob", "workflow_path", "observed_ci"):
        _require(a3.get(key) == AGENT3[key], f"Agent-3 provenance drift: {key}")
    _require(a3.get("api") == "run_fixed_physical_contract_finite_correction_stage", "Agent-3 API drift")
    for key in (
        "physical_contract_sha256_required",
        "before_after_same_physical_contract_sha256",
        "velocity_only_correction",
        "pressure_frozen_through_correction",
        "forcing_frozen_through_correction",
        "candidate_source_identity_required",
        "held_in_held_out_disjoint_required",
        "residual_defined_forcing_forbidden",
        "joint_pressure_or_forcing_updates_require_separate_preregistered_stage",
    ):
        _require(a3.get(key) is True, f"Agent-3 firewall weakened: {key}")
    _require(a3.get("scientific_admission") is False, "Agent-3 queued mechanics cannot be admitted")

    a4 = value.get("agent4_blackbox_full_ns_validator_handoff", {})
    for key in ("pr", "head", "source_path", "source_blob", "workflow_path", "observed_ci"):
        _require(a4.get(key) == AGENT4[key], f"Agent-4 provenance drift: {key}")
    _require(a4.get("public_candidate_surfaces") == ["velocity", "pressure", "forcing", "validation_metadata"], "A4 public-surface drift")
    _require(a4.get("construction_derivatives_allowed") is False, "construction derivatives cannot enter the independent validator")
    _require(a4.get("construction_defect_tensors_allowed") is False, "construction defects cannot enter the independent validator")
    _require(a4.get("training_state_allowed") is False, "training state cannot enter the independent validator")
    operator = a4.get("operator", {})
    _require(operator.get("space") == "centered_cartesian_fourth_order_five_point", "A4 spatial operator drift")
    _require(operator.get("time") == "centered_or_one_sided_fourth_order", "A4 time operator drift")
    _require(operator.get("viscosity") == 0.01, "viscosity drift")
    _require(operator.get("derivative_step_ladder") == DERIVATIVE_STEP_LADDER, "derivative ladder drift")
    heldout = a4.get("heldout", {})
    _require(heldout.get("seed") == 914027, "held-out seed drift")
    _require(heldout.get("sample_count") == 4096, "held-out count drift")
    _require(heldout.get("validation_times") == [0.25, 0.3125, 0.4375, 0.5625, 0.6875, 0.75], "validation times drift")
    _require(heldout.get("evaluation_box") == [[-2.0, 2.0], [-2.0, 2.0], [-2.0, 2.0]], "evaluation box drift")
    _require(heldout.get("time_interval") == [0.25, 0.75], "time interval drift")
    _require(a4.get("canonical_quadrature_ladder") == QUADRATURE_LADDER, "canonical quadrature ladder drift")
    _require(a4.get("quadrature_ladder_assessed_required_for_pde_validated") is True, "quadrature firewall weakened")
    _require(a4.get("physical_contract_sha256_required") is True, "validator physical-contract binding weakened")
    _require(a4.get("complete_candidate_admission_required") is True, "complete-candidate admission weakened")
    _require(a4.get("scientific_admission") is False, "Agent-4 queued validator cannot be admitted")

    seam = value.get("integration_seam", {})
    for key in (
        "correction_and_validator_physical_contract_sha256_must_match",
        "before_after_physical_contract_identity_locked",
        "pressure_identity_immutable_through_velocity_only_correction",
        "forcing_identity_immutable_through_velocity_only_correction",
        "residual_defined_forcing_forbidden",
    ):
        _require(seam.get(key) is True, f"integration firewall weakened: {key}")
    for key in (
        "global_join_present",
        "matched_pressure_present",
        "restricted_forcing_present",
        "restricted_forcing_preregistered",
        "real_correction_velocity_present",
        "complete_candidate_api_ready",
        "current_candidate_eligible_for_full_ns_validation",
        "same_protocol_st006_comparison_available_now",
    ):
        _require(seam.get(key) is False, f"missing dependency laundered as present: {key}")

    _require(value.get("readiness") == {
        "leading_ready": False,
        "oscillatory_ready": True,
        "correction_ready": False,
        "velocity_export_ready": False,
        "pde_validated": False,
    }, "readiness drift")
    _require(value.get("final_gate") == FINAL_GATE, "final 1e-3/1e-5 gate drift")
    _require(value.get("baseline") == {"st006": ST006_BASELINE}, "ST006 baseline drift")
    truth = value.get("truth_boundary", {})
    for key in (
        "queued_ci_is_not_pass",
        "mechanics_or_manufactured_checks_are_not_scientific_admission",
        "sampled_rms_is_not_canonical_volume_l2",
        "monte_carlo_volume_l2_does_not_replace_canonical_quadrature_ladder",
        "visual_success_is_not_pde_validation",
        "residual_defined_forcing_forbidden",
    ):
        _require(truth.get(key) is True, f"truth boundary weakened: {key}")
    _require(truth.get("kokuno_reconstruction_is_paper_exact") is False, "Kokuno reconstruction cannot be relabeled paper-exact")

    supplied_sha = value.get("contract_sha256")
    _require(isinstance(supplied_sha, str) and len(supplied_sha) == 64, "missing contract SHA")
    _require(supplied_sha == _sha256(_contract_payload_without_sha(value)), "contract SHA mismatch")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--exact-head", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    contract = build_contract(exact_head=args.exact_head)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(contract, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"contract_sha256": contract["contract_sha256"], "output": str(output)}, sort_keys=True))


if __name__ == "__main__":
    main()
