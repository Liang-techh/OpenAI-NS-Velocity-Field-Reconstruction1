"""Fail-closed A5 registration of the first strict-inner Kokuno candidate artifact.

Agent 2 #857 materializes a checksum-bound callable/manifest for the currently
legal strict-inner composition

    u_candidate = u_inner_PA10_center + u_osc_complete_curl.

Agent 4 #859 independently saves/reloads that manifest and audits only the
reloaded public ``velocity(x,y,z,t)`` surface with a Richardson-extrapolated
Cartesian derivative path.  This module records that cross-lane handoff without
copying sibling mathematics and without promoting the scoped object into the
full candidate-artifact stage.

The object is still inner-only.  It has no corrected/global leading join,
matched pressure, preregistered restricted forcing, Agent-3 correction velocity,
whole-domain support proof, or complete Navier--Stokes momentum residual.
Consequently ``velocity_export_ready`` and ``pde_validated`` remain false.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from .kokuno_a5_oscillatory_transport_delta_ingest_contract import (
    deterministic_oscillatory_transport_delta_ingest_contract,
    validate_oscillatory_transport_delta_ingest_contract,
)

SCHEMA = "kokuno-a5-strict-inner-candidate-artifact-ingest-contract-v1"
TASK = "KOKUNO-A5-STRICT-INNER-CANDIDATE-ARTIFACT-INGEST-077"
PARENT_A5_PR = 852
PARENT_A5_HEAD = "0e6e3e81791050ffd7fa6770ae1a9f52817daaed"

AGENT2_CANDIDATE_PR = 857
AGENT2_CANDIDATE_HEAD = "a586b7afe4bb47dbfd4b25177586c79d362afff9"
AGENT2_CANDIDATE_SOURCE_BLOB = "82e732af14695c39afa3b631c00b7f20960976da"
AGENT2_CANDIDATE_DEDICATED_RUN = 35516967454
AGENT2_CANDIDATE_TESTS_RUN = 35516967416
AGENT2_CANDIDATE_SCHEMA = "kokuno-a2-strict-inner-leading-oscillatory-candidate-v1"
AGENT2_CANDIDATE_MODULE = (
    "openai_ns_reconstruction.kokuno_strict_inner_leading_oscillatory_candidate"
)
AGENT2_CANDIDATE_CLASS = "KokunoStrictInnerLeadingOscillatoryCandidate"
AGENT2_CANDIDATE_VELOCITY = "velocity"
AGENT2_CANDIDATE_SAVE = "save_manifest"
AGENT2_CANDIDATE_LOAD = "load_manifest"

AGENT4_ARTIFACT_AUDIT_PR = 859
AGENT4_ARTIFACT_AUDIT_HEAD = "5b047f681eda2fe830f0d9663946e5c2c07e5cdc"
AGENT4_ARTIFACT_AUDIT_SOURCE_BLOB = "d65d4931f2734de32e101f21d48c7b8d2381cbde"
AGENT4_ARTIFACT_AUDIT_DEDICATED_RUN = 35517764181
AGENT4_ARTIFACT_AUDIT_TESTS_RUN = 35517764175
AGENT4_ARTIFACT_AUDIT_SCHEMA = (
    "kokuno-a4-strict-inner-candidate-artifact-independent-audit-v1"
)
AGENT4_ARTIFACT_AUDIT_MODULE = (
    "openai_ns_reconstruction.kokuno_a4_strict_inner_candidate_artifact_independent_audit"
)

# Exact-head Actions remain unresolved at this integration cut.  Queued is not PASS.
AGENT2_CANDIDATE_EXACT_HEAD_CI_CONCLUSION: str | None = None
AGENT4_ARTIFACT_AUDIT_EXACT_HEAD_CI_CONCLUSION: str | None = None
AGENT4_ARTIFACT_AUDIT_CONCLUSION: str | None = None

AGENT2_ARTIFACT_PROTOCOL = {
    "scope": "strict-inner PA10 contraction-center plus frozen complete-curl oscillation",
    "composition": "velocity=u_inner_strict+u_osc_complete_curl",
    "public_surface": "velocity(x,y,z,t)->[...,3]",
    "manifest_save_method": AGENT2_CANDIDATE_SAVE,
    "manifest_load_method": AGENT2_CANDIDATE_LOAD,
    "candidate_semantic_sha256_bound": True,
    "manifest_sha256_bound": True,
    "agent1_inner_domain_failure_propagates": True,
    "zero_extension_or_outer_taper_invented": False,
    "frozen_probe_count": 12,
    "batch_scalar_consistency_gate": 1.0e-11,
    "axis_or_outside_oscillatory_support_zero_gate": 1.0e-14,
    "velocity_accepts_residual_pressure_forcing_or_correction_knobs": False,
    "post_observation_retuning_allowed": False,
}

AGENT4_ARTIFACT_AUDIT_PROTOCOL = {
    "independent_reference": "manifest-reloaded-public-velocity-only-Richardson-centered-two-point",
    "seed": 9173521,
    "random_sample_count": 96,
    "axis_and_axis_near_probe_count": 6,
    "nominal_richardson_steps": [1.6e-3, 8.0e-4, 4.0e-4],
    "richardson_extrapolated_order": 4,
    "divergence_sampled_max_gate": 1.0e-5,
    "divergence_sampled_rms_gate": 1.0e-5,
    "resolution_stability_factor": 1.25,
    "resolution_floor": 2.0e-8,
    "velocity_rms_nontriviality_floor": 1.0e-8,
    "negative_controls": [
        "public-velocity-plus-(1e-3*x,0,0)-divergence-injection",
        "candidate-semantic-sha-mutation-with-refreshed-manifest-sha",
        "momentum-gate-1e-3-to-1.001e-3-with-refreshed-manifest-sha",
    ],
    "uses_candidate_internal_derivatives": False,
    "uses_agent1_analytic_or_fd6_derivatives": False,
    "uses_agent2_fd4_or_fd8_transport_verifier": False,
    "post_observation_retuning_allowed": False,
}

NORM_SCOPE_FIREWALL = {
    "agent4_divergence_rms_is_heldout_sampled_rms": True,
    "agent4_divergence_rms_is_whole_domain_volume_weighted_L2": False,
    "whole_domain_boundary_support_assessed": False,
    "complete_ns_momentum_residual_assessed": False,
    "formal_cr001_pde_gate_assessed": False,
    "same_protocol_ST006_comparison_legal": False,
}


def _digest(payload: Mapping[str, Any]) -> str:
    body = dict(payload)
    body.pop("contract_sha256", None)
    raw = json.dumps(
        body, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def deterministic_strict_inner_candidate_artifact_ingest_contract(
    *, exact_head: str | None = None
) -> dict[str, Any]:
    parent = deterministic_oscillatory_transport_delta_ingest_contract(
        exact_head=PARENT_A5_HEAD
    )
    validate_oscillatory_transport_delta_ingest_contract(parent)

    # Registration is deliberately distinct from admission.  In particular,
    # a future A4 scoped divergence PASS is still not a full PDE receipt.
    admitted = bool(
        parent["ingest_status"]["oscillatory_transport_delta_ingest_admitted"]
        and AGENT2_CANDIDATE_EXACT_HEAD_CI_CONCLUSION == "success"
        and AGENT4_ARTIFACT_AUDIT_EXACT_HEAD_CI_CONCLUSION == "success"
        and AGENT4_ARTIFACT_AUDIT_CONCLUSION == "pass"
    )

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "task": TASK,
        "exact_head": exact_head,
        "parent_a5": {"pr": PARENT_A5_PR, "head": PARENT_A5_HEAD},
        "parent_contract_sha256": parent["contract_sha256"],
        "pipeline_position": {
            "registered_scoped_stage": "strict_inner_candidate_artifact",
            "full_pipeline_candidate_artifact_stage_reached": False,
            "reason": (
                "the materialized artifact is inner PA10 center + oscillation only; "
                "global leading, pressure, forcing and correction remain absent"
            ),
        },
        "agent2_candidate_binding": {
            "pr": AGENT2_CANDIDATE_PR,
            "head": AGENT2_CANDIDATE_HEAD,
            "source_blob_sha": AGENT2_CANDIDATE_SOURCE_BLOB,
            "dedicated_run": AGENT2_CANDIDATE_DEDICATED_RUN,
            "tests_run": AGENT2_CANDIDATE_TESTS_RUN,
            "schema": AGENT2_CANDIDATE_SCHEMA,
            "module": AGENT2_CANDIDATE_MODULE,
            "class": AGENT2_CANDIDATE_CLASS,
            "protocol": dict(AGENT2_ARTIFACT_PROTOCOL),
        },
        "agent4_artifact_audit_binding": {
            "pr": AGENT4_ARTIFACT_AUDIT_PR,
            "head": AGENT4_ARTIFACT_AUDIT_HEAD,
            "source_blob_sha": AGENT4_ARTIFACT_AUDIT_SOURCE_BLOB,
            "dedicated_run": AGENT4_ARTIFACT_AUDIT_DEDICATED_RUN,
            "tests_run": AGENT4_ARTIFACT_AUDIT_TESTS_RUN,
            "schema": AGENT4_ARTIFACT_AUDIT_SCHEMA,
            "module": AGENT4_ARTIFACT_AUDIT_MODULE,
            "audited_agent2_pr": AGENT2_CANDIDATE_PR,
            "audited_agent2_head": AGENT2_CANDIDATE_HEAD,
            "protocol": dict(AGENT4_ARTIFACT_AUDIT_PROTOCOL),
            "norm_scope_firewall": dict(NORM_SCOPE_FIREWALL),
        },
        "evidence": {
            "parent_oscillatory_transport_delta_ingest_admitted": bool(
                parent["ingest_status"]["oscillatory_transport_delta_ingest_admitted"]
            ),
            "agent2_candidate_exact_head_ci_conclusion": AGENT2_CANDIDATE_EXACT_HEAD_CI_CONCLUSION,
            "agent4_dedicated_artifact_audit_present": True,
            "agent4_artifact_audit_exact_head_ci_conclusion": AGENT4_ARTIFACT_AUDIT_EXACT_HEAD_CI_CONCLUSION,
            "agent4_artifact_audit_conclusion": AGENT4_ARTIFACT_AUDIT_CONCLUSION,
            "agent4_artifact_scientific_receipt_admitted": False,
            "strict_inner_candidate_artifact_scientifically_admitted": False,
        },
        "strict_inner_candidate_artifact_handoff": {
            "registered": True,
            "status": "registered_unresolved",
            "construction_authority": "Agent2#857",
            "independent_audit_available": True,
            "independent_audit_authority": "Agent4#859",
            "public_velocity_binding": (
                "Agent2#857.KokunoStrictInnerLeadingOscillatoryCandidate.velocity"
            ),
            "save_manifest_binding": (
                "Agent2#857.KokunoStrictInnerLeadingOscillatoryCandidate.save_manifest"
            ),
            "load_manifest_binding": (
                "Agent2#857.KokunoStrictInnerLeadingOscillatoryCandidate.load_manifest"
            ),
            "semantic_candidate_sha256_registered": True,
            "manifest_sha256_registered": True,
            "python_callable_surface_registered": True,
            "matlab_export_materialized": False,
            "velocity_dt_in_artifact_surface": False,
            "pressure_in_artifact_surface": False,
            "restricted_forcing_in_artifact_surface": False,
            "agent3_correction_velocity_in_artifact": False,
            "outer_global_join_in_artifact": False,
            "complete_candidate_artifact": False,
            "eligible_for_velocity_export_ready": False,
            "usable_as_complete_ns_defect_input": False,
            "usable_for_final_independent_pde_validation": False,
        },
        "candidate_api_handoff": dict(parent["candidate_api_handoff"]),
        "differential_operator_handoff": dict(parent["differential_operator_handoff"]),
        "transport_operator_handoff": dict(parent["transport_operator_handoff"]),
        "oscillatory_transport_delta_handoff": dict(parent["oscillatory_transport_delta_handoff"]),
        "artifact_status": {
            "upstream_strict_inner_artifact_implementation_present": True,
            "upstream_manifest_save_load_surface_present": True,
            "upstream_independent_artifact_audit_present": True,
            "a5_artifact_schema_binding_registered": True,
            "a5_materialized_candidate_payload": False,
            "strict_inner_candidate_artifact_ingest_admitted": admitted,
            "complete_global_candidate_artifact_materialized": False,
            "candidate_parameter_digest_for_full_pipeline": None,
            "candidate_validation_receipt_for_full_pipeline": None,
        },
        "staged_comparison": {
            **dict(parent["staged_comparison"]),
            "strict_inner_velocity_artifact_available_upstream": True,
            "strict_inner_artifact_independent_divergence_audit_registered": True,
            "strict_inner_artifact_full_ns_residual_available": False,
            "strict_inner_artifact_same_protocol_ST006_comparison_legal": False,
        },
        "ingest_status": {
            **dict(parent["ingest_status"]),
            "typed_strict_inner_candidate_artifact_registered": True,
            "typed_strict_inner_candidate_artifact_independent_audit_registered": True,
            "strict_inner_candidate_artifact_ingest_admitted": admitted,
            "full_candidate_artifact_ready": False,
        },
        "stage_state": dict(parent["stage_state"]),
        "baseline_vs_kokuno": dict(parent["baseline_vs_kokuno"]),
        "final_project_gates_unchanged": dict(parent["final_project_gates_unchanged"]),
        "truth_boundary": {
            "strict_inner_candidate_artifact_surface_registered": True,
            "agent4_artifact_audit_surface_registered": True,
            "manifest_identity_boundary_registered": True,
            "registration_not_scientific_admission": True,
            "strict_inner_scope_preserved": True,
            "agent2_candidate_ci_promoted": False,
            "agent4_artifact_audit_ci_promoted": False,
            "agent4_artifact_receipt_invented": False,
            "strict_inner_artifact_independently_admitted": False,
            "sampled_rms_laundered_as_volume_l2": False,
            "scoped_divergence_laundered_as_pde_validation": False,
            "full_candidate_artifact_invented": False,
            "full_ns_residual_invented": False,
            "premature_ST006_comparison_performed": False,
            "matched_pressure_invented": False,
            "restricted_forcing_invented": False,
            "correction_velocity_invented": False,
            "outer_global_join_invented": False,
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


def validate_strict_inner_candidate_artifact_ingest_contract(
    payload: Mapping[str, Any]
) -> None:
    if not isinstance(payload, Mapping):
        raise ValueError("contract must be a mapping")
    exact_head = payload.get("exact_head")
    if exact_head is not None:
        if not isinstance(exact_head, str) or len(exact_head) != 40:
            raise ValueError("exact_head must be a 40-character commit id")
        if any(c not in "0123456789abcdef" for c in exact_head):
            raise ValueError("exact_head must be lowercase hexadecimal")
    expected = deterministic_strict_inner_candidate_artifact_ingest_contract(
        exact_head=exact_head
    )
    if dict(payload) != expected:
        raise ValueError("strict-inner candidate artifact ingest contract drifted")
    if payload.get("contract_sha256") != _digest(payload):
        raise ValueError("contract sha256 mismatch")


def write_report(path: str | Path, *, exact_head: str | None = None) -> dict[str, Any]:
    payload = deterministic_strict_inner_candidate_artifact_ingest_contract(
        exact_head=exact_head
    )
    validate_strict_inner_candidate_artifact_ingest_contract(payload)
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return payload


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--exact-head", default=None)
    parser.add_argument(
        "--output",
        default=(
            "artifacts/kokuno_agent5/strict_inner_candidate_artifact_ingest_v1/report.json"
        ),
    )
    args = parser.parse_args(argv)
    write_report(args.output, exact_head=args.exact_head)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
