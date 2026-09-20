"""Agent-5 fail-closed ingest contract for the Agent-1 PA.10 profile surface.

This module is integration glue only.  It registers the exact typed surface that
Agent 5 will accept at the ``source_profile_ingest`` stage once the pinned
Agent-1 exact-head contract and the independent Agent-4 audit have both passed.
It deliberately does not copy Agent-1 mathematics and does not promote a
source-native ``(Y, eta)`` profile into a Cartesian/global leading velocity.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Mapping

from .kokuno_a5_candidate_artifact_contract import (
    FINAL_PROJECT_GATES,
    PINNED_KOKUNO_SOURCE,
    PIPELINE_ORDER,
    ST006_BASELINE,
    deterministic_contract,
)

SCHEMA = "kokuno-a5-leading-profile-ingest-contract-v1"
TASK = "KOKUNO-A5-LEADING-PROFILE-INGEST-CONTRACT-067"
PARENT_A5_PR = 764
PARENT_A5_HEAD = "2c91fd48f178fecbd7a151c5767818f805afcc14"

AGENT1_PROFILE_PR = 768
AGENT1_PROFILE_HEAD = "67f3c2670920bb160a23931ce113b0f05cd1b748"
AGENT1_PROFILE_SCHEMA = "kokuno-pa10-leading-axis-profile-contract-v1"
AGENT1_PROFILE_MODULE = (
    "openai_ns_reconstruction.kokuno_pa10_leading_axis_profile_contract"
)
AGENT1_PROFILE_CLASS = "KokunoPA10LeadingAxisProfileContract"
AGENT1_DEDICATED_RUN = 35489353259

AGENT4_PROFILE_AUDIT_PR = 771
AGENT4_PROFILE_AUDIT_HEAD = "8bc24386b71f9361df851de3d1d2de492003f4a4"

# Frozen observations at this integration cut.  A later A5 increment may replace
# them only with immutable exact-head evidence; callers cannot override them.
AGENT1_EXACT_HEAD_CI_PASSED_AT_FREEZE = False
AGENT4_INDEPENDENT_PROFILE_AUDIT_PASSED_AT_FREEZE = False

REQUIRED_PROFILE_METHODS = (
    "values",
    "derivatives",
    "configuration",
    "from_configuration",
    "save_configuration",
    "load_configuration",
    "report",
)
REQUIRED_VALUE_FIELDS = (
    "Y",
    "eta",
    "d",
    "L",
    "U_star",
    "H_star",
    "W_star",
    "Z_star",
    "Pi_0",
    "Pi_0_eta",
    "chi",
    "zeta_star",
    "u_0",
    "u_0_Y",
)
REQUIRED_DERIVATIVE_FIELDS = (
    "d_eta",
    "L_eta",
    "U_star_eta",
    "H_star_eta",
    "W_star_eta",
    "Z_star_eta",
    "Pi_0_eta_eta",
    "chi_eta",
    "zeta_star_eta",
    "u_0_Y",
    "u_0_eta",
    "u_0_Y_eta",
)


@dataclass(frozen=True)
class LeadingProfileEvidencePin:
    """Immutable sibling-evidence identity, not caller-selected science state."""

    agent1_pr: int = AGENT1_PROFILE_PR
    agent1_head: str = AGENT1_PROFILE_HEAD
    agent1_dedicated_run: int = AGENT1_DEDICATED_RUN
    agent1_exact_head_ci_passed: bool = AGENT1_EXACT_HEAD_CI_PASSED_AT_FREEZE
    agent4_pr: int = AGENT4_PROFILE_AUDIT_PR
    agent4_head: str = AGENT4_PROFILE_AUDIT_HEAD
    agent4_independent_audit_passed: bool = (
        AGENT4_INDEPENDENT_PROFILE_AUDIT_PASSED_AT_FREEZE
    )


def _canonical(payload: Mapping[str, Any]) -> bytes:
    return json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")


def _digest(payload: Mapping[str, Any]) -> str:
    body = dict(payload)
    body.pop("contract_sha256", None)
    return hashlib.sha256(_canonical(body)).hexdigest()


def _profile_ingest_status() -> dict[str, bool]:
    evidence = LeadingProfileEvidencePin()
    admitted = bool(
        evidence.agent1_exact_head_ci_passed
        and evidence.agent4_independent_audit_passed
    )
    return {
        "typed_profile_interface_registered": True,
        "agent1_exact_head_profile_contract_passed": bool(
            evidence.agent1_exact_head_ci_passed
        ),
        "agent4_independent_profile_audit_passed": bool(
            evidence.agent4_independent_audit_passed
        ),
        "source_profile_ingest_admitted": admitted,
        # Admission of this source-native profile surface can never, by itself,
        # promote the global leading candidate or Cartesian velocity.
        "repository_X_to_source_Y_mapping_resolved": False,
        "global_leading_velocity_materialized": False,
        "matched_global_pressure_materialized": False,
        "restricted_forcing_semantics_independently_validated": False,
        "leading_ready": False,
    }


def deterministic_ingest_contract(*, exact_head: str | None = None) -> dict[str, Any]:
    parent = deterministic_contract(exact_head=PARENT_A5_HEAD)
    parent_state = parent["stage_state"]
    if parent_state != {
        "leading_ready": False,
        "oscillatory_ready": True,
        "correction_ready": False,
        "velocity_export_ready": False,
        "pde_validated": False,
    }:
        raise RuntimeError("parent A5 stage state drifted; fail closed")

    if PIPELINE_ORDER[0] != "source_profile_ingest":
        raise RuntimeError("candidate pipeline no longer starts at source_profile_ingest")

    status = _profile_ingest_status()
    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "task": TASK,
        "exact_head": exact_head,
        "parent_a5": {"pr": PARENT_A5_PR, "head": PARENT_A5_HEAD},
        "pipeline_stage": "source_profile_ingest",
        "pipeline_order": list(PIPELINE_ORDER),
        "source_provenance": dict(PINNED_KOKUNO_SOURCE),
        "agent1_profile_binding": {
            "pr": AGENT1_PROFILE_PR,
            "head": AGENT1_PROFILE_HEAD,
            "schema": AGENT1_PROFILE_SCHEMA,
            "module": AGENT1_PROFILE_MODULE,
            "class": AGENT1_PROFILE_CLASS,
            "required_methods": list(REQUIRED_PROFILE_METHODS),
            "required_value_fields": list(REQUIRED_VALUE_FIELDS),
            "required_derivative_fields": list(REQUIRED_DERIVATIVE_FIELDS),
            "coordinate_contract": {
                "native_similarity_coordinates": ["Y", "eta"],
                "source_Y_interval": [0.0, 4.1],
                "repository_X_identification_asserted": False,
                "u_0_is_PA10_contraction_center_not_final_velocity": True,
            },
        },
        "independent_evidence_pin": asdict(LeadingProfileEvidencePin()),
        "ingest_status": status,
        "stage_state": dict(parent_state),
        "baseline_vs_kokuno": {
            "retained_repository_baseline": dict(ST006_BASELINE),
            "kokuno_same_protocol_full_candidate_residual_available": False,
            "comparison_performed": False,
        },
        "final_project_gates_unchanged": dict(FINAL_PROJECT_GATES),
        "truth_boundary": {
            "interface_contract_only_not_profile_admission": True,
            "agent1_queued_ci_laundered_to_pass": False,
            "agent4_queued_audit_laundered_to_pass": False,
            "source_native_profile_promoted_to_global_velocity": False,
            "repository_X_identified_with_source_Y": False,
            "surrogate_promoted_to_candidate": False,
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


def validate_ingest_contract(payload: Mapping[str, Any]) -> None:
    """Fail closed on state laundering or source/profile API drift."""
    if payload.get("schema") != SCHEMA:
        raise ValueError("schema mismatch")
    if payload.get("pipeline_stage") != "source_profile_ingest":
        raise ValueError("pipeline-stage mismatch")
    binding = payload.get("agent1_profile_binding")
    if not isinstance(binding, Mapping):
        raise ValueError("missing Agent-1 profile binding")
    if binding.get("head") != AGENT1_PROFILE_HEAD:
        raise ValueError("Agent-1 exact-head drift")
    if binding.get("schema") != AGENT1_PROFILE_SCHEMA:
        raise ValueError("Agent-1 schema drift")
    coordinate = binding.get("coordinate_contract")
    if not isinstance(coordinate, Mapping):
        raise ValueError("missing coordinate contract")
    if coordinate.get("repository_X_identification_asserted") is not False:
        raise ValueError("repository X/source Y mapping may not be invented")

    evidence = payload.get("independent_evidence_pin")
    if not isinstance(evidence, Mapping):
        raise ValueError("missing evidence pin")
    if evidence != asdict(LeadingProfileEvidencePin()):
        raise ValueError("sibling evidence may not be caller-mutated")

    status = payload.get("ingest_status")
    if status != _profile_ingest_status():
        raise ValueError("profile-ingest state laundering detected")
    if status["source_profile_ingest_admitted"] is not False:
        raise ValueError("queued sibling evidence cannot admit source-profile ingest")
    if status["leading_ready"] is not False:
        raise ValueError("profile ingest alone cannot make leading_ready true")

    stage = payload.get("stage_state")
    if stage != {
        "leading_ready": False,
        "oscillatory_ready": True,
        "correction_ready": False,
        "velocity_export_ready": False,
        "pde_validated": False,
    }:
        raise ValueError("A5 stage state promotion is not authorized")

    gates = payload.get("final_project_gates_unchanged")
    if gates != FINAL_PROJECT_GATES:
        raise ValueError("scientific gate drift")
    baseline = payload.get("baseline_vs_kokuno")
    if not isinstance(baseline, Mapping):
        raise ValueError("missing baseline boundary")
    if baseline.get("kokuno_same_protocol_full_candidate_residual_available") is not False:
        raise ValueError("full Kokuno residual was fabricated")
    if baseline.get("comparison_performed") is not False:
        raise ValueError("ST006 comparison is unavailable before a full candidate")

    truth = payload.get("truth_boundary")
    if not isinstance(truth, Mapping):
        raise ValueError("missing truth boundary")
    forbidden_true = (
        "agent1_queued_ci_laundered_to_pass",
        "agent4_queued_audit_laundered_to_pass",
        "source_native_profile_promoted_to_global_velocity",
        "repository_X_identified_with_source_Y",
        "surrogate_promoted_to_candidate",
        "free_residual_defined_forcing_allowed",
        "threshold_relaxed",
        "kokuno_replay_used_as_final_independent_validation",
        "leading_ready",
        "correction_ready",
        "velocity_export_ready",
        "pde_validated",
    )
    if any(truth.get(key) is not False for key in forbidden_true):
        raise ValueError("truth-boundary promotion detected")
    if payload.get("contract_sha256") != _digest(payload):
        raise ValueError("contract digest mismatch")


def save_contract(path: str | Path, *, exact_head: str | None = None) -> dict[str, Any]:
    payload = deterministic_ingest_contract(exact_head=exact_head)
    validate_ingest_contract(payload)
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
        "schema": payload["schema"],
        "exact_head": payload["exact_head"],
        "source_profile_ingest_admitted": payload["ingest_status"]["source_profile_ingest_admitted"],
        "leading_ready": payload["stage_state"]["leading_ready"],
        "contract_sha256": payload["contract_sha256"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
