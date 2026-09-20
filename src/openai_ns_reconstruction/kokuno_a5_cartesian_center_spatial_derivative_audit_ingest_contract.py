"""Fail-closed A5 registration of the Agent-4 spatial-derivative audit.

This layer is intentionally glue only.  Agent 1 #819 owns the analytic inner
PA.10 Cartesian Jacobian/divergence/vorticity surface; Agent 4 #823 owns an
implementation-distinct public-velocity-only FD6 audit of that surface.  A5
binds the exact heads/protocol into the pipeline without copying either lane's
mathematics and without promoting queued evidence to scientific admission.

The Agent-4 divergence ``L2`` quantity here is a sampled RMS on the frozen
inner held-out points.  It is not the canonical CR001 volume-weighted spatial
L2 norm and it is not a final full-candidate PDE validation result.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from .kokuno_a5_cartesian_center_spatial_derivative_ingest_contract import (
    deterministic_cartesian_center_spatial_derivative_ingest_contract,
    validate_cartesian_center_spatial_derivative_ingest_contract,
)

SCHEMA = "kokuno-a5-cartesian-center-spatial-derivative-audit-ingest-contract-v1"
TASK = "KOKUNO-A5-CARTESIAN-CENTER-SPATIAL-DERIVATIVE-AUDIT-INGEST-074"
PARENT_A5_PR = 824
PARENT_A5_HEAD = "089e7589f3f0a9ba50f3cc764bd20c348f926654"

AGENT4_AUDIT_PR = 823
AGENT4_AUDIT_HEAD = "56f4cf63e79fd08fe005ea3dce061f9a2f8ce50e"
AGENT4_AUDIT_SOURCE_BLOB = "9e6a4c95ef0f52d952c3278d028a201d85efa578"
AGENT4_AUDIT_DEDICATED_RUN = 35506076292
AGENT4_AUDIT_TESTS_RUN = 35506076340
AGENT4_AUDIT_SCHEMA = (
    "kokuno-a4-pa10-cartesian-center-spatial-derivatives-independent-audit-v1"
)
AGENT4_AUDIT_MODULE = (
    "openai_ns_reconstruction."
    "kokuno_a4_pa10_cartesian_center_spatial_derivatives_independent_audit"
)

# Frozen facts at this integration cut.  A dedicated A4 surface now exists,
# but its exact-head Actions are queued and no scientific receipt is admitted.
AGENT4_EXACT_HEAD_CI_CONCLUSION: str | None = None
AGENT4_INDEPENDENT_AUDIT_CONCLUSION: str | None = None

AGENT4_PROTOCOL = {
    "independent_reference": "public-velocity-only-centered-FD6",
    "seed": 9173481,
    "random_sample_count": 384,
    "axis_and_axis_near_probe_count": 5,
    "space_steps": [1.2e-3, 6.0e-4, 3.0e-4],
    "jacobian_relative_rms_gate": 2.0e-6,
    "jacobian_relative_max_gate": 2.0e-5,
    "vorticity_relative_rms_gate": 3.0e-6,
    "vorticity_relative_max_gate": 3.0e-5,
    "refinement_ratio_gate": 8.0,
    "refinement_floor": 2.0e-9,
    "divergence_sampled_max_gate": 1.0e-5,
    "divergence_sampled_rms_gate": 1.0e-5,
    "jacobian_rms_nontriviality_floor": 1.0e-8,
    "vorticity_rms_nontriviality_floor": 1.0e-8,
    "public_closure_gate": 1.0e-13,
    "negative_controls": [
        "0.99x-public-jacobian",
        "public-vorticity-sign-flip",
        "+1e-3-to-du_dx",
    ],
    "post_observation_retuning_allowed": False,
}

DIVERGENCE_NORM_SCOPE = {
    "sampled_max": "max(abs(divergence)) on frozen inner held-out samples",
    "sampled_l2_rms": "sqrt(mean(divergence**2)) on frozen inner held-out samples",
    "canonical_cr001_volume_weighted_spatial_l2_established": False,
    "whole_domain_global_candidate_divergence_established": False,
    "same_metric_as_final_cr001_divergence_l2": False,
}


def _digest(payload: Mapping[str, Any]) -> str:
    body = dict(payload)
    body.pop("contract_sha256", None)
    raw = json.dumps(
        body, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def deterministic_cartesian_center_spatial_derivative_audit_ingest_contract(
    *, exact_head: str | None = None
) -> dict[str, Any]:
    parent = deterministic_cartesian_center_spatial_derivative_ingest_contract(
        exact_head=PARENT_A5_HEAD
    )
    validate_cartesian_center_spatial_derivative_ingest_contract(parent)
    parent_status = parent["ingest_status"]

    admitted = bool(
        parent_status["inner_cartesian_center_spatial_derivative_ingest_admitted"]
        and AGENT4_EXACT_HEAD_CI_CONCLUSION == "success"
        and AGENT4_INDEPENDENT_AUDIT_CONCLUSION == "pass"
    )

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "task": TASK,
        "exact_head": exact_head,
        "parent_a5": {"pr": PARENT_A5_PR, "head": PARENT_A5_HEAD},
        "pipeline_stage": "source_profile_ingest",
        "pipeline_substage": (
            "inner_cartesian_center_spatial_derivative_independent_audit_ingest"
        ),
        "parent_spatial_derivative_contract_sha256": parent["contract_sha256"],
        "agent4_spatial_derivative_audit_binding": {
            "pr": AGENT4_AUDIT_PR,
            "head": AGENT4_AUDIT_HEAD,
            "source_blob_sha": AGENT4_AUDIT_SOURCE_BLOB,
            "dedicated_run": AGENT4_AUDIT_DEDICATED_RUN,
            "tests_run": AGENT4_AUDIT_TESTS_RUN,
            "schema": AGENT4_AUDIT_SCHEMA,
            "module": AGENT4_AUDIT_MODULE,
            "audited_agent1_pr": 819,
            "audited_agent1_head": (
                "cb54e6e1a9061cacf78e446dc66bc992f64a0f8d"
            ),
            "protocol": dict(AGENT4_PROTOCOL),
            "divergence_norm_scope": dict(DIVERGENCE_NORM_SCOPE),
        },
        "evidence": {
            "parent_inner_cartesian_center_spatial_derivative_ingest_admitted": bool(
                parent_status[
                    "inner_cartesian_center_spatial_derivative_ingest_admitted"
                ]
            ),
            "agent4_dedicated_spatial_derivative_audit_present": True,
            "agent4_exact_head_ci_conclusion": AGENT4_EXACT_HEAD_CI_CONCLUSION,
            "agent4_independent_audit_conclusion": AGENT4_INDEPENDENT_AUDIT_CONCLUSION,
            "agent4_scientific_receipt_admitted": False,
        },
        "candidate_api_handoff": dict(parent["candidate_api_handoff"]),
        "differential_operator_handoff": {
            **dict(parent["differential_operator_handoff"]),
            "independent_audit_available": True,
            "independent_audit_status": "registered_unresolved",
            "independent_audit_authority": "Agent4#823",
            "usable_for_inner_spatial_derivative_admission_if_passes": True,
            "usable_for_final_independent_pde_validation": False,
            "reason": (
                "A4 #823 audits only the inner PA.10 spatial derivative surface; "
                "its exact-head CI is unresolved, its divergence L2 is sampled RMS, "
                "and no global velocity/pressure/restricted-forcing candidate exists"
            ),
        },
        "ingest_status": {
            **dict(parent_status),
            "typed_inner_cartesian_center_spatial_derivative_independent_audit_registered": True,
            "inner_cartesian_center_spatial_derivative_ingest_admitted": admitted,
            "outer_global_leading_join_materialized": False,
            "global_leading_velocity_materialized": False,
            "matched_global_pressure_materialized": False,
            "restricted_forcing_composite_validated": False,
            "leading_ready": False,
        },
        "stage_state": dict(parent["stage_state"]),
        "baseline_vs_kokuno": dict(parent["baseline_vs_kokuno"]),
        "final_project_gates_unchanged": dict(parent["final_project_gates_unchanged"]),
        "truth_boundary": {
            "a4_audit_surface_registered": True,
            "interface_registration_not_scientific_admission": True,
            "agent4_exact_head_ci_promoted": False,
            "agent4_scientific_receipt_invented": False,
            "inner_spatial_derivatives_independently_admitted": False,
            "sampled_divergence_rms_laundered_as_cr001_volume_l2": False,
            "a4_inner_audit_used_as_final_pde_validator": False,
            "source_center_promoted_to_final_fixed_point": False,
            "outer_global_join_invented": False,
            "global_leading_velocity_invented": False,
            "matched_global_pressure_invented": False,
            "restricted_forcing_composite_invented": False,
            "complete_candidate_api_promoted": False,
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


def validate_cartesian_center_spatial_derivative_audit_ingest_contract(
    payload: Mapping[str, Any]
) -> None:
    parent = deterministic_cartesian_center_spatial_derivative_ingest_contract(
        exact_head=PARENT_A5_HEAD
    )
    validate_cartesian_center_spatial_derivative_ingest_contract(parent)

    if payload.get("schema") != SCHEMA:
        raise ValueError("schema mismatch")
    if payload.get("parent_a5") != {"pr": PARENT_A5_PR, "head": PARENT_A5_HEAD}:
        raise ValueError("parent drift")
    if payload.get("parent_spatial_derivative_contract_sha256") != parent[
        "contract_sha256"
    ]:
        raise ValueError("parent digest drift")

    a4 = payload.get("agent4_spatial_derivative_audit_binding")
    if not isinstance(a4, Mapping):
        raise ValueError("missing Agent-4 audit binding")
    expected = {
        "pr": AGENT4_AUDIT_PR,
        "head": AGENT4_AUDIT_HEAD,
        "source_blob_sha": AGENT4_AUDIT_SOURCE_BLOB,
        "dedicated_run": AGENT4_AUDIT_DEDICATED_RUN,
        "tests_run": AGENT4_AUDIT_TESTS_RUN,
        "schema": AGENT4_AUDIT_SCHEMA,
        "module": AGENT4_AUDIT_MODULE,
        "audited_agent1_pr": 819,
        "audited_agent1_head": "cb54e6e1a9061cacf78e446dc66bc992f64a0f8d",
    }
    for key, value in expected.items():
        if a4.get(key) != value:
            raise ValueError(f"Agent-4 {key} drift")
    if a4.get("protocol") != AGENT4_PROTOCOL:
        raise ValueError("Agent-4 protocol drift")
    if a4.get("divergence_norm_scope") != DIVERGENCE_NORM_SCOPE:
        raise ValueError("Agent-4 divergence norm scope drift")

    expected_evidence = {
        "parent_inner_cartesian_center_spatial_derivative_ingest_admitted": False,
        "agent4_dedicated_spatial_derivative_audit_present": True,
        "agent4_exact_head_ci_conclusion": None,
        "agent4_independent_audit_conclusion": None,
        "agent4_scientific_receipt_admitted": False,
    }
    if payload.get("evidence") != expected_evidence:
        raise ValueError("evidence laundering")

    if payload.get("candidate_api_handoff") != parent["candidate_api_handoff"]:
        raise ValueError("candidate API drift")

    ops = payload.get("differential_operator_handoff")
    if not isinstance(ops, Mapping):
        raise ValueError("missing differential operator handoff")
    for name in ("velocity_jacobian", "divergence", "vorticity"):
        if ops.get(name) != parent["differential_operator_handoff"][name]:
            raise ValueError(f"{name} handoff drift")
    required_ops = {
        "independent_audit_available": True,
        "independent_audit_status": "registered_unresolved",
        "independent_audit_authority": "Agent4#823",
        "usable_for_inner_spatial_derivative_admission_if_passes": True,
        "usable_for_final_independent_pde_validation": False,
    }
    for key, value in required_ops.items():
        if ops.get(key) != value:
            raise ValueError(f"audit handoff {key} drift")

    status = payload.get("ingest_status")
    if not isinstance(status, Mapping):
        raise ValueError("missing ingest status")
    expected_status = dict(parent["ingest_status"])
    expected_status.update(
        {
            "typed_inner_cartesian_center_spatial_derivative_independent_audit_registered": True,
            "inner_cartesian_center_spatial_derivative_ingest_admitted": False,
            "outer_global_leading_join_materialized": False,
            "global_leading_velocity_materialized": False,
            "matched_global_pressure_materialized": False,
            "restricted_forcing_composite_validated": False,
            "leading_ready": False,
        }
    )
    if dict(status) != expected_status:
        raise ValueError("ingest status promotion/drift")

    if payload.get("stage_state") != parent["stage_state"]:
        raise ValueError("stage state drift")
    if payload.get("baseline_vs_kokuno") != parent["baseline_vs_kokuno"]:
        raise ValueError("baseline comparison drift")
    if payload.get("final_project_gates_unchanged") != parent[
        "final_project_gates_unchanged"
    ]:
        raise ValueError("project gate drift")

    truth = payload.get("truth_boundary")
    if not isinstance(truth, Mapping):
        raise ValueError("missing truth boundary")
    if truth.get("a4_audit_surface_registered") is not True:
        raise ValueError("A4 audit registration lost")
    if truth.get("interface_registration_not_scientific_admission") is not True:
        raise ValueError("registration/admission boundary lost")
    for key, value in truth.items():
        if key in {
            "a4_audit_surface_registered",
            "interface_registration_not_scientific_admission",
        }:
            continue
        if value is not False:
            raise ValueError(f"truth promotion: {key}")

    if payload.get("contract_sha256") != _digest(payload):
        raise ValueError("contract digest mismatch")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--exact-head", default=None)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)

    payload = deterministic_cartesian_center_spatial_derivative_audit_ingest_contract(
        exact_head=args.exact_head
    )
    validate_cartesian_center_spatial_derivative_audit_ingest_contract(payload)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"schema": SCHEMA, "contract_sha256": payload["contract_sha256"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
