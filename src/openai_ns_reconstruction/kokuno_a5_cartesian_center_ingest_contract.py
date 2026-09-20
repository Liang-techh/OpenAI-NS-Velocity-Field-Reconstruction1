"""Fail-closed A5 registration of the PA.10 inner Cartesian center velocity.

This seam harvests Agent-1 PR #803's first candidate-facing Cartesian
``velocity(x,y,z,t)`` for the source-C-normalized PA.10 contraction center and
pins Agent-4 PR #806 as the corresponding independent audit surface.

Registration is deliberately narrower than scientific admission.  The field is
inner-only: no final fixed-point correction, outer/global join, matched pressure,
restricted forcing composite, complete Kokuno candidate, or held-out momentum
validation exists here.  Consequently this module cannot promote
``leading_ready``, ``velocity_export_ready`` or ``pde_validated``.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from .kokuno_a5_candidate_artifact_contract import FINAL_PROJECT_GATES, ST006_BASELINE
from .kokuno_a5_source_c_normalization_ingest_contract import (
    deterministic_source_c_normalization_ingest_contract,
    validate_source_c_normalization_ingest_contract,
)

SCHEMA = "kokuno-a5-cartesian-center-ingest-contract-v1"
TASK = "KOKUNO-A5-CARTESIAN-CENTER-INGEST-CONTRACT-071"
PARENT_A5_PR = 799
PARENT_A5_HEAD = "a4fe0ab22b4f8ef499751453963502706f9e2d66"

AGENT1_CARTESIAN_CENTER_PR = 803
AGENT1_CARTESIAN_CENTER_HEAD = "0666fff748d7c2659774b0aacba70042d045aab5"
AGENT1_CARTESIAN_CENTER_SOURCE_BLOB = "f404eff536b91248f76095c111d788482a19f6a3"
AGENT1_CARTESIAN_CENTER_DEDICATED_RUN = 35499772525
AGENT1_CARTESIAN_CENTER_SCHEMA = "kokuno-pa10-cartesian-center-velocity-v1"
AGENT1_CARTESIAN_CENTER_MODULE = (
    "openai_ns_reconstruction.kokuno_pa10_cartesian_center_velocity"
)
AGENT1_CARTESIAN_CENTER_CLASS = "KokunoPA10CartesianCenterVelocity"

AGENT4_CARTESIAN_CENTER_PR = 806
AGENT4_CARTESIAN_CENTER_HEAD = "1f2c121483dac630a51d95823df59303a56e78da"
AGENT4_CARTESIAN_CENTER_SOURCE_BLOB = "ef913488dee576ffe79dd359a55d881e62127bb6"
AGENT4_CARTESIAN_CENTER_DEDICATED_RUN = 35500500957
AGENT4_CARTESIAN_CENTER_SCHEMA = (
    "kokuno-a4-pa10-cartesian-center-independent-audit-v1"
)

# Both exact-head workflows are queued at this freeze.  These are evidence
# constants, not caller-controlled switches.
AGENT1_EXACT_HEAD_CI_CONCLUSION: str | None = None
AGENT4_EXACT_HEAD_CI_CONCLUSION: str | None = None
AGENT4_INNER_CARTESIAN_CENTER_INDEPENDENTLY_AUDITED = False

REQUIRED_AGENT1_MEMBERS = (
    "similarity_coordinates",
    "cartesian_from_similarity",
    "values",
    "velocity",
    "field_configuration",
    "field_sha256",
    "evidence_configuration",
    "save_configuration",
    "from_configuration",
    "load_configuration",
    "truth_boundary",
    "report",
    "save_report",
)

SOURCE_FORMULAS = {
    "coordinates": (
        "tau=1-t; z=q^D eta; tau=q(1-eta^2); "
        "X=(x^2+y^2)/(2q); A=1/2+h; D=1/2-h"
    ),
    "q_inverse": (
        "q-z^2 q^(2h)-tau=0; q_*=|z|^(1/D); "
        "g'(q)>=1-2h for q>=q_*"
    ),
    "cartesian_velocity": (
        "u1=(v0/(2q))x-q^(-A-1/2)Fy; "
        "u2=(v0/(2q))y+q^(-A-1/2)Fx; u3=q^(-A)U"
    ),
    "profile_input": "F=F_0(X,eta); U=U_0(X,eta); v0=v_0(X,eta)",
}


def _digest(payload: Mapping[str, Any]) -> str:
    body = dict(payload)
    body.pop("contract_sha256", None)
    raw = json.dumps(
        body, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def deterministic_cartesian_center_ingest_contract(
    *, exact_head: str | None = None
) -> dict[str, Any]:
    parent = deterministic_source_c_normalization_ingest_contract(
        exact_head=PARENT_A5_HEAD
    )
    validate_source_c_normalization_ingest_contract(parent)

    parent_status = parent["ingest_status"]
    independently_admitted = bool(
        parent_status["source_C_normalized_physical_center_ingest_admitted"]
        and parent_status["source_complex_C_normalization_independently_admitted"]
        and AGENT1_EXACT_HEAD_CI_CONCLUSION == "success"
        and AGENT4_EXACT_HEAD_CI_CONCLUSION == "success"
        and AGENT4_INNER_CARTESIAN_CENTER_INDEPENDENTLY_AUDITED
    )

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "task": TASK,
        "exact_head": exact_head,
        "parent_a5": {"pr": PARENT_A5_PR, "head": PARENT_A5_HEAD},
        "pipeline_stage": "source_profile_ingest",
        "pipeline_substage": "inner_cartesian_center_velocity_ingest",
        "parent_source_C_contract_sha256": parent["contract_sha256"],
        "agent1_cartesian_center_binding": {
            "pr": AGENT1_CARTESIAN_CENTER_PR,
            "head": AGENT1_CARTESIAN_CENTER_HEAD,
            "source_blob_sha": AGENT1_CARTESIAN_CENTER_SOURCE_BLOB,
            "dedicated_run": AGENT1_CARTESIAN_CENTER_DEDICATED_RUN,
            "schema": AGENT1_CARTESIAN_CENTER_SCHEMA,
            "module": AGENT1_CARTESIAN_CENTER_MODULE,
            "class": AGENT1_CARTESIAN_CENTER_CLASS,
            "required_members": list(REQUIRED_AGENT1_MEMBERS),
            "source_formulas": dict(SOURCE_FORMULAS),
            "public_velocity_signature": "velocity(x,y,z,t)->[...,3]",
            "registered_time_interval": [0.25, 0.75],
            "domain_boundary": {
                "inner_source_X_only": True,
                "outside_inner_X_raises": True,
                "outer_global_continuation_materialized": False,
            },
            "identity_boundary": {
                "field_semantic_sha256_exposed": True,
                "field_identity_separate_from_evidence_receipt": True,
                "source_C_certificate_execution_knobs_excluded_from_field_identity": True,
            },
        },
        "agent4_cartesian_center_binding": {
            "pr": AGENT4_CARTESIAN_CENTER_PR,
            "head": AGENT4_CARTESIAN_CENTER_HEAD,
            "source_blob_sha": AGENT4_CARTESIAN_CENTER_SOURCE_BLOB,
            "dedicated_run": AGENT4_CARTESIAN_CENTER_DEDICATED_RUN,
            "schema": AGENT4_CARTESIAN_CENTER_SCHEMA,
            "independent_q_solver": "long-double safeguarded Newton + bracket contraction",
            "independent_velocity_reconstruction": True,
            "public_velocity_only_fd4_divergence": True,
            "fd4_spatial_steps": [0.004, 0.002, 0.001],
            "mapping_sample_count": 2048,
            "divergence_sample_count": 48,
        },
        "evidence": {
            "agent1_exact_head_ci_conclusion": AGENT1_EXACT_HEAD_CI_CONCLUSION,
            "agent4_exact_head_ci_conclusion": AGENT4_EXACT_HEAD_CI_CONCLUSION,
            "agent4_inner_cartesian_center_independently_audited": (
                AGENT4_INNER_CARTESIAN_CENTER_INDEPENDENTLY_AUDITED
            ),
            "source_C_independent_admission_inherited": bool(
                parent_status["source_complex_C_normalization_independently_admitted"]
            ),
        },
        "candidate_api_handoff": {
            "velocity": "Agent1#803.KokunoPA10CartesianCenterVelocity.velocity",
            "velocity_dt": None,
            "pressure": None,
            "forcing": None,
            "complete_candidate_api_ready": False,
            "reason": (
                "inner Cartesian center exists, but global leading join, matched pressure, "
                "restricted forcing, and time-derivative composite are not materialized"
            ),
        },
        "ingest_status": {
            **dict(parent_status),
            "typed_inner_cartesian_center_velocity_interface_registered": True,
            "typed_inner_cartesian_center_independent_audit_registered": True,
            "inner_cartesian_center_velocity_ingest_admitted": independently_admitted,
            "inner_cartesian_center_field_identity_registered": True,
            "outer_global_leading_join_materialized": False,
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
            "inner_cartesian_center_is_final_fixed_point": False,
            "source_C_independent_admission_invented": False,
            "queued_agent1_ci_promoted": False,
            "queued_agent4_audit_promoted": False,
            "inner_only_velocity_promoted_to_global_leading": False,
            "outer_global_join_invented": False,
            "velocity_dt_invented": False,
            "matched_global_pressure_invented": False,
            "restricted_forcing_composite_invented": False,
            "complete_candidate_artifact_invented": False,
            "leading_only_ns_residual_invented": False,
            "full_kokuno_residual_invented": False,
            "premature_ST006_comparison_performed": False,
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


def validate_cartesian_center_ingest_contract(payload: Mapping[str, Any]) -> None:
    parent = deterministic_source_c_normalization_ingest_contract(
        exact_head=PARENT_A5_HEAD
    )
    validate_source_c_normalization_ingest_contract(parent)

    if payload.get("schema") != SCHEMA:
        raise ValueError("schema mismatch")
    if payload.get("parent_a5") != {"pr": PARENT_A5_PR, "head": PARENT_A5_HEAD}:
        raise ValueError("parent drift")
    if payload.get("parent_source_C_contract_sha256") != parent["contract_sha256"]:
        raise ValueError("parent digest drift")

    a1 = payload.get("agent1_cartesian_center_binding")
    if not isinstance(a1, Mapping):
        raise ValueError("missing Agent-1 Cartesian-center binding")
    expected_a1 = {
        "pr": AGENT1_CARTESIAN_CENTER_PR,
        "head": AGENT1_CARTESIAN_CENTER_HEAD,
        "source_blob_sha": AGENT1_CARTESIAN_CENTER_SOURCE_BLOB,
        "dedicated_run": AGENT1_CARTESIAN_CENTER_DEDICATED_RUN,
        "schema": AGENT1_CARTESIAN_CENTER_SCHEMA,
        "module": AGENT1_CARTESIAN_CENTER_MODULE,
        "class": AGENT1_CARTESIAN_CENTER_CLASS,
    }
    for key, value in expected_a1.items():
        if a1.get(key) != value:
            raise ValueError(f"Agent-1 Cartesian-center {key} drift")
    if tuple(a1.get("required_members", ())) != REQUIRED_AGENT1_MEMBERS:
        raise ValueError("Agent-1 required member drift")
    if a1.get("source_formulas") != SOURCE_FORMULAS:
        raise ValueError("source formula drift")
    if a1.get("public_velocity_signature") != "velocity(x,y,z,t)->[...,3]":
        raise ValueError("velocity signature drift")
    if a1.get("registered_time_interval") != [0.25, 0.75]:
        raise ValueError("time interval drift")

    a4 = payload.get("agent4_cartesian_center_binding")
    if not isinstance(a4, Mapping):
        raise ValueError("missing Agent-4 Cartesian-center binding")
    expected_a4 = {
        "pr": AGENT4_CARTESIAN_CENTER_PR,
        "head": AGENT4_CARTESIAN_CENTER_HEAD,
        "source_blob_sha": AGENT4_CARTESIAN_CENTER_SOURCE_BLOB,
        "dedicated_run": AGENT4_CARTESIAN_CENTER_DEDICATED_RUN,
        "schema": AGENT4_CARTESIAN_CENTER_SCHEMA,
        "fd4_spatial_steps": [0.004, 0.002, 0.001],
        "mapping_sample_count": 2048,
        "divergence_sample_count": 48,
    }
    for key, value in expected_a4.items():
        if a4.get(key) != value:
            raise ValueError(f"Agent-4 Cartesian-center {key} drift")

    if payload.get("evidence") != {
        "agent1_exact_head_ci_conclusion": None,
        "agent4_exact_head_ci_conclusion": None,
        "agent4_inner_cartesian_center_independently_audited": False,
        "source_C_independent_admission_inherited": False,
    }:
        raise ValueError("evidence laundering")

    api = payload.get("candidate_api_handoff")
    if not isinstance(api, Mapping):
        raise ValueError("missing candidate API handoff")
    if api.get("velocity") != "Agent1#803.KokunoPA10CartesianCenterVelocity.velocity":
        raise ValueError("velocity handoff drift")
    for key in ("velocity_dt", "pressure", "forcing"):
        if api.get(key) is not None:
            raise ValueError(f"invented candidate API component: {key}")
    if api.get("complete_candidate_api_ready") is not False:
        raise ValueError("premature complete candidate API")

    parent_status = parent["ingest_status"]
    expected_status = {
        **dict(parent_status),
        "typed_inner_cartesian_center_velocity_interface_registered": True,
        "typed_inner_cartesian_center_independent_audit_registered": True,
        "inner_cartesian_center_velocity_ingest_admitted": False,
        "inner_cartesian_center_field_identity_registered": True,
        "outer_global_leading_join_materialized": False,
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


def save_cartesian_center_ingest_contract(
    path: str | Path, *, exact_head: str | None = None
) -> dict[str, Any]:
    payload = deterministic_cartesian_center_ingest_contract(exact_head=exact_head)
    validate_cartesian_center_ingest_contract(payload)
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return payload


def _main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exact-head", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    payload = save_cartesian_center_ingest_contract(
        args.output, exact_head=args.exact_head
    )
    print(
        json.dumps(
            {
                "schema": payload["schema"],
                "exact_head": payload["exact_head"],
                "contract_sha256": payload["contract_sha256"],
                "inner_cartesian_center_velocity_ingest_admitted": payload[
                    "ingest_status"
                ]["inner_cartesian_center_velocity_ingest_admitted"],
                "leading_ready": payload["stage_state"]["leading_ready"],
                "pde_validated": payload["stage_state"]["pde_validated"],
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
