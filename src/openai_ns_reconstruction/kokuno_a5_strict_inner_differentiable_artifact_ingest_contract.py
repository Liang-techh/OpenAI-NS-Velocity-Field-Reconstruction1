"""Fail-closed A5 registration of the strict-inner differentiable candidate artifact.

Agent 2 #866 extends the existing strict-inner ``u_inner + u_osc`` candidate
artifact with a checksum-bound fixed-Cartesian ``velocity_dt`` surface while
preserving the #857 velocity semantic identity. Agent 4 #869 independently
saves/reloads that differentiable artifact and reconstructs the time derivative
from only the reloaded public ``velocity`` using a Richardson-extrapolated
centered two-point operator.

This module records that typed cross-lane handoff. It does not copy sibling
mathematics and it does not promote the scoped inner object into a global
candidate, a complete Navier--Stokes residual, or PDE validation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from .kokuno_a5_strict_inner_candidate_artifact_ingest_contract import (
    deterministic_strict_inner_candidate_artifact_ingest_contract,
    validate_strict_inner_candidate_artifact_ingest_contract,
)

SCHEMA = "kokuno-a5-strict-inner-differentiable-artifact-ingest-contract-v1"
TASK = "KOKUNO-A5-STRICT-INNER-DIFFERENTIABLE-ARTIFACT-INGEST-078"
PARENT_A5_PR = 860
PARENT_A5_HEAD = "d295670e51802abc345801bcdfef0c99bbf91e98"

AGENT2_PR = 866
AGENT2_HEAD = "68f84128f07b6743d368a9ab7a051e441f5a59b2"
AGENT2_SOURCE_BLOB = "8c09048a1d15dff321642e21ee699cce7d19831a"
AGENT2_DEDICATED_RUN = 35520204098
AGENT2_TESTS_RUN = 35520204018
AGENT2_SCHEMA = "kokuno-a2-strict-inner-differentiable-candidate-v1"
AGENT2_MODULE = (
    "openai_ns_reconstruction."
    "kokuno_strict_inner_leading_oscillatory_differentiable_candidate"
)
AGENT2_CLASS = "KokunoStrictInnerLeadingOscillatoryDifferentiableCandidate"

AGENT4_PR = 869
AGENT4_HEAD = "ef10076fee8f5f4c5bb589ca7d2cd90ae0c889a4"
AGENT4_SOURCE_BLOB = "ff32f0cd641ff59f1d4e2fa7535c100087ce00f6"
AGENT4_DEDICATED_RUN = 35520902624
AGENT4_TESTS_RUN = 35520902790
AGENT4_SCHEMA = "kokuno-a4-strict-inner-differentiable-candidate-independent-audit-v1"
AGENT4_MODULE = (
    "openai_ns_reconstruction."
    "kokuno_a4_strict_inner_differentiable_candidate_independent_audit"
)

# Exact-head Actions are unresolved at this integration cut. Queued is not PASS.
AGENT2_EXACT_HEAD_CI_CONCLUSION: str | None = None
AGENT4_EXACT_HEAD_CI_CONCLUSION: str | None = None
AGENT4_INDEPENDENT_AUDIT_CONCLUSION: str | None = None

AGENT2_PROTOCOL = {
    "scope": "strict-inner PA10 contraction-center plus frozen complete-curl oscillation",
    "velocity_surface": "velocity(x,y,z,t)->[...,3]",
    "velocity_dt_surface": "velocity_dt(x,y,z,t)->[...,3]",
    "base_velocity_candidate_identity_preserved": True,
    "velocity_candidate_sha256_preserved": True,
    "differentiable_sha256_bound": True,
    "manifest_sha256_bound": True,
    "save_manifest_present": True,
    "load_manifest_present": True,
    "production_velocity_dt_uses_finite_difference": False,
    "agent1_inner_domain_failure_propagates": True,
    "outer_taper_or_zero_extension_invented": False,
    "fd4_time_steps": [4.0e-3, 2.0e-3, 1.0e-3],
    "fd4_refinement_ratio_gate": 8.0,
    "fine_relative_rms_gate": 2.0e-3,
    "fine_relative_sampled_max_gate": 5.0e-3,
    "velocity_dt_rms_floor": 1.0e-10,
    "post_observation_retuning_allowed": False,
}

AGENT4_PROTOCOL = {
    "independent_reference": (
        "saved-reloaded-public-velocity-only-Richardson-centered-two-point-time-derivative"
    ),
    "seed": 9173531,
    "random_sample_count": 128,
    "axis_and_axis_near_probe_count": 6,
    "nominal_richardson_time_steps": [3.2e-3, 1.6e-3, 8.0e-4],
    "richardson_extrapolated_order": 4,
    "fine_relative_rms_gate": 2.0e-3,
    "fine_relative_sampled_max_gate": 5.0e-3,
    "refinement_ratio_gate": 6.0,
    "refinement_floor": 2.0e-10,
    "velocity_dt_rms_floor": 1.0e-10,
    "negative_controls": [
        "production-velocity_dt-times-0.99",
        "production-axial-velocity_dt-sign-flip",
        "differentiable-sha-mutation-with-refreshed-manifest-sha",
        "momentum-gate-1e-3-to-1.001e-3-with-refreshed-manifest-sha",
    ],
    "uses_agent1_analytic_time_derivative_as_reference": False,
    "uses_agent2_fd4_verifier_as_reference": False,
    "uses_candidate_component_tensors": False,
    "post_observation_retuning_allowed": False,
}

NORM_SCOPE_FIREWALL = {
    "velocity_dt_consistency_is_complete_ns_residual": False,
    "velocity_dt_consistency_is_cr001_momentum_norm": False,
    "whole_domain_volume_weighted_l2_assessed": False,
    "whole_domain_boundary_support_assessed": False,
    "same_protocol_st006_comparison_legal": False,
    "formal_full_domain_pde_gate_assessed": False,
}


def _digest(payload: Mapping[str, Any]) -> str:
    body = dict(payload)
    body.pop("contract_sha256", None)
    raw = json.dumps(
        body, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def deterministic_strict_inner_differentiable_artifact_ingest_contract(
    *, exact_head: str | None = None
) -> dict[str, Any]:
    parent = deterministic_strict_inner_candidate_artifact_ingest_contract(
        exact_head=PARENT_A5_HEAD
    )
    validate_strict_inner_candidate_artifact_ingest_contract(parent)

    admitted = bool(
        parent["ingest_status"]["strict_inner_candidate_artifact_ingest_admitted"]
        and AGENT2_EXACT_HEAD_CI_CONCLUSION == "success"
        and AGENT4_EXACT_HEAD_CI_CONCLUSION == "success"
        and AGENT4_INDEPENDENT_AUDIT_CONCLUSION == "pass"
    )

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "task": TASK,
        "exact_head": exact_head,
        "parent_a5": {"pr": PARENT_A5_PR, "head": PARENT_A5_HEAD},
        "parent_contract_sha256": parent["contract_sha256"],
        "pipeline_position": {
            "registered_scoped_stage": "strict_inner_differentiable_candidate_artifact",
            "velocity_surface_registered": True,
            "velocity_dt_surface_registered": True,
            "full_pipeline_candidate_artifact_stage_reached": False,
            "reason": (
                "artifact remains strict-inner and lacks corrected/global leading join, "
                "matched pressure, restricted forcing and Agent-3 correction velocity"
            ),
        },
        "agent2_binding": {
            "pr": AGENT2_PR,
            "head": AGENT2_HEAD,
            "source_blob_sha": AGENT2_SOURCE_BLOB,
            "dedicated_run": AGENT2_DEDICATED_RUN,
            "tests_run": AGENT2_TESTS_RUN,
            "schema": AGENT2_SCHEMA,
            "module": AGENT2_MODULE,
            "class": AGENT2_CLASS,
            "protocol": dict(AGENT2_PROTOCOL),
        },
        "agent4_binding": {
            "pr": AGENT4_PR,
            "head": AGENT4_HEAD,
            "source_blob_sha": AGENT4_SOURCE_BLOB,
            "dedicated_run": AGENT4_DEDICATED_RUN,
            "tests_run": AGENT4_TESTS_RUN,
            "schema": AGENT4_SCHEMA,
            "module": AGENT4_MODULE,
            "audited_agent2_pr": AGENT2_PR,
            "audited_agent2_head": AGENT2_HEAD,
            "protocol": dict(AGENT4_PROTOCOL),
            "norm_scope_firewall": dict(NORM_SCOPE_FIREWALL),
        },
        "evidence": {
            "parent_strict_inner_artifact_ingest_admitted": bool(
                parent["ingest_status"]["strict_inner_candidate_artifact_ingest_admitted"]
            ),
            "agent2_exact_head_ci_conclusion": AGENT2_EXACT_HEAD_CI_CONCLUSION,
            "agent4_dedicated_differentiable_artifact_audit_present": True,
            "agent4_exact_head_ci_conclusion": AGENT4_EXACT_HEAD_CI_CONCLUSION,
            "agent4_independent_audit_conclusion": AGENT4_INDEPENDENT_AUDIT_CONCLUSION,
            "agent4_scientific_receipt_admitted": False,
            "strict_inner_differentiable_artifact_scientifically_admitted": False,
        },
        "differentiable_artifact_handoff": {
            "registered": True,
            "status": "registered_unresolved",
            "construction_authority": "Agent2#866",
            "independent_audit_available": True,
            "independent_audit_authority": "Agent4#869",
            "public_velocity_binding": (
                "Agent2#866.KokunoStrictInnerLeadingOscillatoryDifferentiableCandidate.velocity"
            ),
            "public_velocity_dt_binding": (
                "Agent2#866.KokunoStrictInnerLeadingOscillatoryDifferentiableCandidate.velocity_dt"
            ),
            "save_manifest_binding": (
                "Agent2#866.KokunoStrictInnerLeadingOscillatoryDifferentiableCandidate.save_manifest"
            ),
            "load_manifest_binding": (
                "Agent2#866.KokunoStrictInnerLeadingOscillatoryDifferentiableCandidate.load_manifest"
            ),
            "base_velocity_candidate_sha256_preserved": True,
            "differentiable_sha256_registered": True,
            "manifest_sha256_registered": True,
            "python_velocity_and_velocity_dt_callable_surface_registered": True,
            "matlab_export_materialized": False,
            "pressure_in_artifact_surface": False,
            "restricted_forcing_in_artifact_surface": False,
            "agent3_correction_velocity_in_artifact": False,
            "outer_global_join_in_artifact": False,
            "complete_candidate_artifact": False,
            "eligible_for_velocity_export_ready": False,
            "usable_as_complete_ns_defect_input": False,
            "usable_for_final_independent_pde_validation": False,
        },
        "candidate_api_handoff": {
            **dict(parent["candidate_api_handoff"]),
            "strict_inner_artifact_velocity": "Agent2#866.velocity",
            "strict_inner_artifact_velocity_dt": "Agent2#866.velocity_dt",
            "pressure": None,
            "forcing": None,
            "complete_candidate_api_ready": False,
        },
        "artifact_status": {
            **dict(parent["artifact_status"]),
            "upstream_strict_inner_differentiable_artifact_implementation_present": True,
            "upstream_velocity_dt_surface_present": True,
            "upstream_differentiable_manifest_save_load_surface_present": True,
            "upstream_independent_velocity_dt_artifact_audit_present": True,
            "strict_inner_differentiable_artifact_ingest_admitted": admitted,
            "complete_global_candidate_artifact_materialized": False,
        },
        "ingest_status": {
            **dict(parent["ingest_status"]),
            "typed_strict_inner_differentiable_artifact_registered": True,
            "typed_strict_inner_differentiable_artifact_independent_audit_registered": True,
            "strict_inner_differentiable_artifact_ingest_admitted": admitted,
            "full_candidate_artifact_ready": False,
        },
        "stage_state": dict(parent["stage_state"]),
        "baseline_vs_kokuno": dict(parent["baseline_vs_kokuno"]),
        "final_project_gates_unchanged": dict(parent["final_project_gates_unchanged"]),
        "truth_boundary": {
            "strict_inner_differentiable_artifact_surface_registered": True,
            "agent4_differentiable_artifact_audit_surface_registered": True,
            "registration_not_scientific_admission": True,
            "strict_inner_scope_preserved": True,
            "agent2_ci_promoted": False,
            "agent4_ci_promoted": False,
            "agent4_receipt_invented": False,
            "strict_inner_differentiable_artifact_independently_admitted": False,
            "velocity_dt_consistency_laundered_as_pde_residual": False,
            "whole_domain_volume_l2_invented": False,
            "full_candidate_artifact_invented": False,
            "full_ns_residual_invented": False,
            "premature_st006_comparison_performed": False,
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


def validate_strict_inner_differentiable_artifact_ingest_contract(
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
    expected = deterministic_strict_inner_differentiable_artifact_ingest_contract(
        exact_head=exact_head
    )
    if dict(payload) != expected:
        raise ValueError("strict-inner differentiable artifact ingest contract drifted")
    if payload.get("contract_sha256") != _digest(payload):
        raise ValueError("contract sha256 mismatch")


def write_report(path: str | Path, *, exact_head: str | None = None) -> dict[str, Any]:
    payload = deterministic_strict_inner_differentiable_artifact_ingest_contract(
        exact_head=exact_head
    )
    validate_strict_inner_differentiable_artifact_ingest_contract(payload)
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
