"""Fail-closed A5 registration of Agent-1's source-C-normalized PA.10 center.

This integration seam harvests only the typed interface introduced by Agent 1
PR #796.  It registers the explicit repository choice C=1000 and the requirement
that the same C be used by the physical F map and the outer schedule, while
keeping scientific admission closed until exact-head CI and a dedicated Agent-4
independent source-C audit both exist and pass.

It does not treat the older Agent-4 physical-center mapping audit (#790), which
was built against #787's pre-normalized configuration, as independent evidence
for the new C=1000 complex-domain certificate.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from .kokuno_a5_candidate_artifact_contract import FINAL_PROJECT_GATES, ST006_BASELINE
from .kokuno_a5_physical_center_ingest_contract import (
    deterministic_physical_center_ingest_contract,
    validate_physical_center_ingest_contract,
)

SCHEMA = "kokuno-a5-source-c-normalization-ingest-contract-v1"
TASK = "KOKUNO-A5-SOURCE-C-NORMALIZATION-INGEST-CONTRACT-070"
PARENT_A5_PR = 791
PARENT_A5_HEAD = "566a3fe7d53918c19016aa11d55c3e071dff779d"

AGENT1_SOURCE_C_PR = 796
AGENT1_SOURCE_C_HEAD = "393abd5976d4c621ce37dad8544a1b645fdbe83b"
AGENT1_SOURCE_C_SOURCE_BLOB = "722640b0e6e43a2bba521776398d2f9cc8d2a052"
AGENT1_SOURCE_C_DEDICATED_RUN = 35497206269
AGENT1_SOURCE_C_SCHEMA = "kokuno-pa10-source-c-normalized-physical-center-v1"
AGENT1_SOURCE_C_MODULE = (
    "openai_ns_reconstruction.kokuno_pa10_source_c_normalized_physical_center"
)
AGENT1_SOURCE_C_CLASS = "KokunoPA10SourceCNormalizedPhysicalCenter"
SELECTED_SOURCE_C = 1000.0

# #796 exact-head CI is still queued at this freeze.  These values are fixed
# evidence fields rather than caller inputs.
AGENT1_EXACT_HEAD_CI_CONCLUSION: str | None = None
AGENT1_SELF_CERTIFICATE_EXECUTED_AND_PASSED: bool | None = None

# No dedicated Agent-4 independent audit of #796's C=1000 complex-domain
# certificate exists at this freeze.  #790 remains a mapping audit of #787.
AGENT4_SOURCE_C_AUDIT_PR: int | None = None
AGENT4_SOURCE_C_AUDIT_HEAD: str | None = None
AGENT4_SOURCE_C_AUDIT_DEDICATED_RUN: int | None = None
AGENT4_SOURCE_C_AUDIT_CONCLUSION: str | None = None
AGENT4_SOURCE_C_INDEPENDENTLY_ADMITTED = False

REQUIRED_MEMBERS = (
    "source_C_certificate",
    "values",
    "derivatives",
    "configuration",
    "from_configuration",
    "save_configuration",
    "load_configuration",
    "truth_boundary",
    "report",
    "save_report",
)
SOURCE_FORMULAS = {
    "phi_star": "phi_*=exp(Lambda int_0^eta zeta_*(w)dw)",
    "zeta_star": "zeta_*=-L H_*/(H_*^2+sigma_*^2)",
    "source_C_normalization": "C>=sup_{eta in Omega}|phi_*(eta)|",
    "normalized_g": "g=phi_*/C; sup_Omega|g|<=1",
    "physical_F": "F=(phi_*/C)Phi=g Phi",
    "same_outer_C": "the same source-normalization C enters the outer schedule",
}


def _digest(payload: Mapping[str, Any]) -> str:
    body = dict(payload)
    body.pop("contract_sha256", None)
    raw = json.dumps(
        body, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def deterministic_source_c_normalization_ingest_contract(
    *, exact_head: str | None = None
) -> dict[str, Any]:
    parent = deterministic_physical_center_ingest_contract(exact_head=PARENT_A5_HEAD)
    validate_physical_center_ingest_contract(parent)

    parent_status = parent["ingest_status"]
    independently_admitted = bool(
        parent_status["physical_center_profile_ingest_admitted"]
        and AGENT1_EXACT_HEAD_CI_CONCLUSION == "success"
        and AGENT1_SELF_CERTIFICATE_EXECUTED_AND_PASSED is True
        and AGENT4_SOURCE_C_AUDIT_PR is not None
        and AGENT4_SOURCE_C_AUDIT_HEAD is not None
        and AGENT4_SOURCE_C_AUDIT_DEDICATED_RUN is not None
        and AGENT4_SOURCE_C_AUDIT_CONCLUSION == "success"
        and AGENT4_SOURCE_C_INDEPENDENTLY_ADMITTED
    )

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "task": TASK,
        "exact_head": exact_head,
        "parent_a5": {"pr": PARENT_A5_PR, "head": PARENT_A5_HEAD},
        "pipeline_stage": "source_profile_ingest",
        "pipeline_substage": "source_C_normalized_physical_center_ingest",
        "parent_physical_center_contract_sha256": parent["contract_sha256"],
        "agent1_source_C_binding": {
            "pr": AGENT1_SOURCE_C_PR,
            "head": AGENT1_SOURCE_C_HEAD,
            "source_blob_sha": AGENT1_SOURCE_C_SOURCE_BLOB,
            "dedicated_run": AGENT1_SOURCE_C_DEDICATED_RUN,
            "schema": AGENT1_SOURCE_C_SCHEMA,
            "module": AGENT1_SOURCE_C_MODULE,
            "class": AGENT1_SOURCE_C_CLASS,
            "selected_source_C": SELECTED_SOURCE_C,
            "required_members": list(REQUIRED_MEMBERS),
            "source_formulas": dict(SOURCE_FORMULAS),
            "selection_boundary": {
                "chosen_from_public_source_condition": True,
                "chosen_from_NS_residual": False,
                "chosen_from_CR001_energy": False,
                "chosen_from_ST006": False,
                "same_C_required_in_physical_F_and_outer_schedule": True,
            },
        },
        "independent_audit_requirement": {
            "dedicated_agent4_source_C_audit_required": True,
            "prior_agent4_physical_mapping_audit_pr": 790,
            "prior_mapping_audit_is_source_C_authority": False,
            "source_C_audit_pr": AGENT4_SOURCE_C_AUDIT_PR,
            "source_C_audit_head": AGENT4_SOURCE_C_AUDIT_HEAD,
            "source_C_audit_dedicated_run": AGENT4_SOURCE_C_AUDIT_DEDICATED_RUN,
        },
        "evidence": {
            "agent1_exact_head_ci_conclusion": AGENT1_EXACT_HEAD_CI_CONCLUSION,
            "agent1_self_certificate_executed_and_passed": (
                AGENT1_SELF_CERTIFICATE_EXECUTED_AND_PASSED
            ),
            "agent4_dedicated_source_C_audit_present": False,
            "agent4_source_C_audit_conclusion": AGENT4_SOURCE_C_AUDIT_CONCLUSION,
            "agent4_source_C_independently_admitted": (
                AGENT4_SOURCE_C_INDEPENDENTLY_ADMITTED
            ),
        },
        "normalization_registration": {
            "selected_source_C_registered": SELECTED_SOURCE_C,
            "same_C_physical_F_and_outer_schedule_required": True,
            "agent1_source_C_self_certificate_interface_registered": True,
            "source_complex_C_normalization_independently_admitted": independently_admitted,
            "source_C_used_as_CR001_energy_normalizer": False,
            "CR001_energy_used_as_source_C_proof": False,
            "residual_driven_C_tuning_allowed": False,
        },
        "ingest_status": {
            **dict(parent_status),
            "typed_source_C_normalized_physical_center_interface_registered": True,
            "source_C_normalized_physical_center_ingest_admitted": independently_admitted,
            "source_complex_C_normalization_agent1_self_certificate_passed": False,
            "source_complex_C_normalization_independently_admitted": independently_admitted,
            "global_leading_velocity_materialized": False,
            "matched_global_pressure_materialized": False,
            "restricted_forcing_composite_validated": False,
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
            "interface_registration_not_scientific_admission": True,
            "queued_agent1_ci_promoted": False,
            "agent1_self_certificate_promoted_without_execution": False,
            "missing_agent4_source_C_audit_invented": False,
            "prior_mapping_audit_reused_as_source_C_certificate": False,
            "source_C_used_as_residual_amplitude_knob": False,
            "source_C_used_as_CR001_energy_normalizer": False,
            "contraction_center_promoted_to_final_fixed_point": False,
            "global_cartesian_velocity_invented": False,
            "matched_global_pressure_invented": False,
            "restricted_forcing_composite_invented": False,
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


def validate_source_c_normalization_ingest_contract(payload: Mapping[str, Any]) -> None:
    parent = deterministic_physical_center_ingest_contract(exact_head=PARENT_A5_HEAD)
    validate_physical_center_ingest_contract(parent)

    if payload.get("schema") != SCHEMA:
        raise ValueError("schema mismatch")
    if payload.get("parent_a5") != {"pr": PARENT_A5_PR, "head": PARENT_A5_HEAD}:
        raise ValueError("parent drift")
    if payload.get("parent_physical_center_contract_sha256") != parent["contract_sha256"]:
        raise ValueError("parent digest drift")

    binding = payload.get("agent1_source_C_binding")
    if not isinstance(binding, Mapping):
        raise ValueError("missing Agent-1 source-C binding")
    expected_binding = {
        "head": AGENT1_SOURCE_C_HEAD,
        "source_blob_sha": AGENT1_SOURCE_C_SOURCE_BLOB,
        "schema": AGENT1_SOURCE_C_SCHEMA,
        "module": AGENT1_SOURCE_C_MODULE,
        "class": AGENT1_SOURCE_C_CLASS,
        "selected_source_C": SELECTED_SOURCE_C,
    }
    for key, value in expected_binding.items():
        if binding.get(key) != value:
            raise ValueError(f"Agent-1 source-C {key} drift")
    if binding.get("source_formulas") != SOURCE_FORMULAS:
        raise ValueError("source formula drift")
    if tuple(binding.get("required_members", ())) != REQUIRED_MEMBERS:
        raise ValueError("required member drift")

    audit = payload.get("independent_audit_requirement")
    expected_audit = {
        "dedicated_agent4_source_C_audit_required": True,
        "prior_agent4_physical_mapping_audit_pr": 790,
        "prior_mapping_audit_is_source_C_authority": False,
        "source_C_audit_pr": None,
        "source_C_audit_head": None,
        "source_C_audit_dedicated_run": None,
    }
    if audit != expected_audit:
        raise ValueError("independent audit laundering")

    expected_evidence = {
        "agent1_exact_head_ci_conclusion": None,
        "agent1_self_certificate_executed_and_passed": None,
        "agent4_dedicated_source_C_audit_present": False,
        "agent4_source_C_audit_conclusion": None,
        "agent4_source_C_independently_admitted": False,
    }
    if payload.get("evidence") != expected_evidence:
        raise ValueError("evidence laundering")

    if payload.get("normalization_registration") != {
        "selected_source_C_registered": SELECTED_SOURCE_C,
        "same_C_physical_F_and_outer_schedule_required": True,
        "agent1_source_C_self_certificate_interface_registered": True,
        "source_complex_C_normalization_independently_admitted": False,
        "source_C_used_as_CR001_energy_normalizer": False,
        "CR001_energy_used_as_source_C_proof": False,
        "residual_driven_C_tuning_allowed": False,
    }:
        raise ValueError("normalization state laundering")

    parent_status = parent["ingest_status"]
    expected_status = {
        **dict(parent_status),
        "typed_source_C_normalized_physical_center_interface_registered": True,
        "source_C_normalized_physical_center_ingest_admitted": False,
        "source_complex_C_normalization_agent1_self_certificate_passed": False,
        "source_complex_C_normalization_independently_admitted": False,
        "global_leading_velocity_materialized": False,
        "matched_global_pressure_materialized": False,
        "restricted_forcing_composite_validated": False,
        "leading_ready": False,
    }
    if payload.get("ingest_status") != expected_status:
        raise ValueError("ingest state laundering")
    if payload.get("stage_state") != parent["stage_state"]:
        raise ValueError("stage state promotion")
    if payload.get("final_project_gates_unchanged") != FINAL_PROJECT_GATES:
        raise ValueError("gate drift")

    baseline = payload.get("baseline_vs_kokuno")
    if not isinstance(baseline, Mapping):
        raise ValueError("missing baseline comparison state")
    if baseline.get("retained_repository_baseline") != ST006_BASELINE:
        raise ValueError("ST006 baseline drift")
    if baseline.get("kokuno_same_protocol_full_candidate_residual_available") is not False:
        raise ValueError("fabricated full-candidate residual")
    if baseline.get("comparison_performed") is not False:
        raise ValueError("premature ST006 comparison")

    truth = payload.get("truth_boundary")
    if not isinstance(truth, Mapping):
        raise ValueError("missing truth boundary")
    if truth.get("interface_registration_not_scientific_admission") is not True:
        raise ValueError("registration/admission boundary lost")
    for key, value in truth.items():
        if key == "interface_registration_not_scientific_admission":
            continue
        if value is not False:
            raise ValueError(f"truth-boundary promotion: {key}")

    if payload.get("contract_sha256") != _digest(payload):
        raise ValueError("contract digest mismatch")


def save_source_c_normalization_ingest_contract(
    path: str | Path, *, exact_head: str | None = None
) -> dict[str, Any]:
    payload = deterministic_source_c_normalization_ingest_contract(exact_head=exact_head)
    validate_source_c_normalization_ingest_contract(payload)
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return payload


def _main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exact-head", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    payload = save_source_c_normalization_ingest_contract(
        args.output, exact_head=args.exact_head
    )
    print(json.dumps({
        "schema": payload["schema"],
        "exact_head": payload["exact_head"],
        "contract_sha256": payload["contract_sha256"],
        "source_C_normalized_physical_center_ingest_admitted": payload["ingest_status"]["source_C_normalized_physical_center_ingest_admitted"],
        "pde_validated": payload["stage_state"]["pde_validated"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
