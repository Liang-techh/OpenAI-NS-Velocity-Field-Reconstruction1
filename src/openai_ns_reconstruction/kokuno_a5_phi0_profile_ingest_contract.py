"""Fail-closed A5 registration of Agent-1's PA.10 Phi_0 profile API."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from .kokuno_a5_candidate_artifact_contract import FINAL_PROJECT_GATES, PIPELINE_ORDER, ST006_BASELINE
from .kokuno_a5_leading_profile_ingest_contract import deterministic_ingest_contract, validate_ingest_contract

SCHEMA = "kokuno-a5-phi0-profile-ingest-contract-v1"
TASK = "KOKUNO-A5-PHI0-PROFILE-INGEST-CONTRACT-068"
PARENT_A5_PR = 772
PARENT_A5_HEAD = "bd5e0d11224327f5b1cd674843e17e8f2d45d642"
AGENT1_PHI0_PR = 778
AGENT1_PHI0_HEAD = "1deec8fb36bd9f8058d179ae3e6895f92aa9ecb3"
AGENT1_PHI0_SOURCE_BLOB = "43198ab5b10331525820d62191a4e79c0873a735"
AGENT1_PHI0_DEDICATED_RUN = 35491602168
AGENT1_PHI0_SCHEMA = "kokuno-pa10-phi0-profile-contract-v1"
AGENT1_PHI0_MODULE = "openai_ns_reconstruction.kokuno_pa10_phi0_profile_contract"
AGENT1_PHI0_CLASS = "KokunoPA10Phi0ProfileContract"
AGENT1_PHI0_CI_PASSED = False
AGENT4_PHI0_AUDIT_PRESENT = False
AGENT4_PHI0_AUDIT_PASSED = False
REQUIRED_METHODS = (
    "f0_triplet", "values", "derivatives", "configuration", "from_configuration",
    "save_configuration", "load_configuration", "report", "save_report",
)
REQUIRED_VALUE_FIELDS = (
    "Y", "eta", "chi", "f0_argument", "Phi_0", "u_0",
    "f0_prime", "f0_second", "f0_tail_upper",
)
REQUIRED_DERIVATIVE_FIELDS = (
    "Phi_0_Y", "Phi_0_eta", "Phi_0_YY", "Phi_0_Y_eta",
    "u_0_Y", "u_0_eta", "u_0_Y_eta",
)


def _digest(payload: Mapping[str, Any]) -> str:
    body = dict(payload)
    body.pop("contract_sha256", None)
    raw = json.dumps(body, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    return hashlib.sha256(raw).hexdigest()


def deterministic_phi0_ingest_contract(*, exact_head: str | None = None) -> dict[str, Any]:
    parent = deterministic_ingest_contract(exact_head=PARENT_A5_HEAD)
    validate_ingest_contract(parent)
    if PIPELINE_ORDER[0] != "source_profile_ingest":
        raise RuntimeError("pipeline drift")
    parent_status = parent["ingest_status"]
    admitted = bool(
        parent_status["source_profile_ingest_admitted"]
        and AGENT1_PHI0_CI_PASSED
        and AGENT4_PHI0_AUDIT_PRESENT
        and AGENT4_PHI0_AUDIT_PASSED
    )
    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "task": TASK,
        "exact_head": exact_head,
        "parent_a5": {"pr": PARENT_A5_PR, "head": PARENT_A5_HEAD},
        "pipeline_stage": "source_profile_ingest",
        "pipeline_substage": "source_phi0_contraction_center_ingest",
        "pipeline_order": list(PIPELINE_ORDER),
        "parent_axis_profile_contract_sha256": parent["contract_sha256"],
        "agent1_phi0_binding": {
            "pr": AGENT1_PHI0_PR,
            "head": AGENT1_PHI0_HEAD,
            "source_blob_sha": AGENT1_PHI0_SOURCE_BLOB,
            "dedicated_run": AGENT1_PHI0_DEDICATED_RUN,
            "schema": AGENT1_PHI0_SCHEMA,
            "module": AGENT1_PHI0_MODULE,
            "class": AGENT1_PHI0_CLASS,
            "required_methods": list(REQUIRED_METHODS),
            "required_value_fields": list(REQUIRED_VALUE_FIELDS),
            "required_derivative_fields": list(REQUIRED_DERIVATIVE_FIELDS),
            "coordinate_contract": {
                "native_similarity_coordinates": ["Y", "eta"],
                "source_Y_interval": [0.0, 4.1],
                "repository_X_identification_asserted": False,
                "physical_F_from_Phi_mapping_resolved": False,
                "Phi_0_is_contraction_center_not_final_corrected_Phi": True,
            },
            "numerical_realization": {
                "f0_series_terms": 64,
                "finite_series_is_repository_numerical_realization": True,
            },
        },
        "evidence": {
            "agent1_exact_head_ci_passed": AGENT1_PHI0_CI_PASSED,
            "agent4_phi0_independent_audit_present": AGENT4_PHI0_AUDIT_PRESENT,
            "agent4_phi0_independent_audit_passed": AGENT4_PHI0_AUDIT_PASSED,
        },
        "ingest_status": {
            "typed_axis_profile_interface_registered": parent_status["typed_profile_interface_registered"],
            "axis_profile_ingest_admitted": parent_status["source_profile_ingest_admitted"],
            "typed_phi0_profile_interface_registered": True,
            "phi0_profile_ingest_admitted": admitted,
            "repository_X_to_source_Y_mapping_resolved": False,
            "physical_F_from_Phi_mapping_resolved": False,
            "fixed_point_correction_materialized": False,
            "global_leading_velocity_materialized": False,
            "leading_ready": False,
        },
        "stage_state": dict(parent["stage_state"]),
        "baseline_vs_kokuno": {
            "retained_repository_baseline": dict(ST006_BASELINE),
            "kokuno_same_protocol_full_candidate_residual_available": False,
            "comparison_performed": False,
        },
        "final_project_gates_unchanged": dict(FINAL_PROJECT_GATES),
        "truth_boundary": {
            "interface_only_not_admission": True,
            "queued_or_missing_evidence_promoted": False,
            "repository_X_identified_with_source_Y": False,
            "physical_F_from_Phi_mapping_invented": False,
            "Phi0_promoted_to_final_corrected_Phi": False,
            "source_native_profile_promoted_to_global_velocity": False,
            "free_residual_defined_forcing_allowed": False,
            "threshold_relaxed": False,
            "kokuno_replay_used_as_final_independent_validation": False,
            "leading_ready": False,
            "correction_ready": False,
            "velocity_export_ready": False,
            "pde_validated": False,
        },
    }
    payload["contract_sha256"] = _digest(payload)
    return payload


def validate_phi0_ingest_contract(payload: Mapping[str, Any]) -> None:
    parent = deterministic_ingest_contract(exact_head=PARENT_A5_HEAD)
    validate_ingest_contract(parent)
    if payload.get("schema") != SCHEMA:
        raise ValueError("schema mismatch")
    if payload.get("parent_a5") != {"pr": PARENT_A5_PR, "head": PARENT_A5_HEAD}:
        raise ValueError("parent drift")
    if payload.get("parent_axis_profile_contract_sha256") != parent["contract_sha256"]:
        raise ValueError("parent digest drift")
    binding = payload.get("agent1_phi0_binding")
    if not isinstance(binding, Mapping):
        raise ValueError("missing binding")
    for key, value in {
        "head": AGENT1_PHI0_HEAD,
        "source_blob_sha": AGENT1_PHI0_SOURCE_BLOB,
        "schema": AGENT1_PHI0_SCHEMA,
    }.items():
        if binding.get(key) != value:
            raise ValueError(f"Agent-1 Phi0 {key} drift")
    coordinate = binding["coordinate_contract"]
    if coordinate["repository_X_identification_asserted"] is not False:
        raise ValueError("X/Y mapping invented")
    if coordinate["physical_F_from_Phi_mapping_resolved"] is not False:
        raise ValueError("physical F/Phi mapping invented")
    if payload.get("evidence") != {
        "agent1_exact_head_ci_passed": False,
        "agent4_phi0_independent_audit_present": False,
        "agent4_phi0_independent_audit_passed": False,
    }:
        raise ValueError("evidence laundering")
    status = payload.get("ingest_status")
    if status != {
        "typed_axis_profile_interface_registered": True,
        "axis_profile_ingest_admitted": False,
        "typed_phi0_profile_interface_registered": True,
        "phi0_profile_ingest_admitted": False,
        "repository_X_to_source_Y_mapping_resolved": False,
        "physical_F_from_Phi_mapping_resolved": False,
        "fixed_point_correction_materialized": False,
        "global_leading_velocity_materialized": False,
        "leading_ready": False,
    }:
        raise ValueError("ingest state laundering")
    if payload.get("stage_state") != parent["stage_state"]:
        raise ValueError("stage state promotion")
    if payload.get("final_project_gates_unchanged") != FINAL_PROJECT_GATES:
        raise ValueError("gate drift")
    baseline = payload["baseline_vs_kokuno"]
    if baseline["kokuno_same_protocol_full_candidate_residual_available"] is not False:
        raise ValueError("fabricated residual")
    if baseline["comparison_performed"] is not False:
        raise ValueError("premature ST006 comparison")
    truth = payload["truth_boundary"]
    for key in (
        "queued_or_missing_evidence_promoted",
        "repository_X_identified_with_source_Y",
        "physical_F_from_Phi_mapping_invented",
        "Phi0_promoted_to_final_corrected_Phi",
        "source_native_profile_promoted_to_global_velocity",
        "free_residual_defined_forcing_allowed",
        "threshold_relaxed",
        "kokuno_replay_used_as_final_independent_validation",
        "leading_ready", "correction_ready", "velocity_export_ready", "pde_validated",
    ):
        if truth[key] is not False:
            raise ValueError("truth-boundary promotion")
    if payload.get("contract_sha256") != _digest(payload):
        raise ValueError("contract digest mismatch")


def save_contract(path: str | Path, *, exact_head: str | None = None) -> dict[str, Any]:
    payload = deterministic_phi0_ingest_contract(exact_head=exact_head)
    validate_phi0_ingest_contract(payload)
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return payload


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exact-head")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    payload = save_contract(args.output, exact_head=args.exact_head)
    print(json.dumps({
        "exact_head": payload["exact_head"],
        "phi0_profile_ingest_admitted": payload["ingest_status"]["phi0_profile_ingest_admitted"],
        "leading_ready": payload["stage_state"]["leading_ready"],
        "contract_sha256": payload["contract_sha256"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
