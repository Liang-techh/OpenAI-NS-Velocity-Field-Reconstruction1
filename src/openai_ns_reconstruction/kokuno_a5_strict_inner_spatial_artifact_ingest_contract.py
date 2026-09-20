"""Fail-closed A5 registration of the strict-inner spatial candidate artifact.

Agent 2 #874 extends the checksum-bound strict-inner ``u_inner + u_osc``
artifact from public ``velocity`` + ``velocity_dt`` to a first-class public
``velocity_jacobian`` surface while preserving the existing velocity and time-
derivative semantic identities.  The new spatial identity binds the analytic
Agent-1 inner Jacobian, the frozen Agent-2 oscillatory FD6 realization, and the
additive composition law.

This module records that typed handoff only.  At this integration cut there is
no dedicated Agent-4 independent audit for #874.  Agent-2's construction-side
FD4 comparison therefore cannot be promoted into independent scientific
admission.  The artifact remains strict-inner and still lacks the corrected /
global leading join, Agent-3 correction velocity, matched pressure and
preregistered restricted forcing required for a complete Navier--Stokes
candidate.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from .kokuno_a5_strict_inner_differentiable_artifact_ingest_contract import (
    deterministic_strict_inner_differentiable_artifact_ingest_contract,
    validate_strict_inner_differentiable_artifact_ingest_contract,
)

SCHEMA = "kokuno-a5-strict-inner-spatial-artifact-ingest-contract-v1"
TASK = "KOKUNO-A5-STRICT-INNER-SPATIAL-ARTIFACT-INGEST-079"
PARENT_A5_PR = 870
PARENT_A5_HEAD = "b155d5a170ad56bb5a2dde9d5f4e0f43d1e2b474"

AGENT2_PR = 874
AGENT2_HEAD = "9cb3869b9cfd2d5b8dcb1da222df74d80e12d0c2"
AGENT2_SOURCE_BLOB = "276c943cd1e769c36d1615366266bf3b2c633040"
AGENT2_DEDICATED_RUN = 35523243070
AGENT2_TESTS_RUN = 35523243046
AGENT2_SCHEMA = "kokuno-a2-strict-inner-spatial-candidate-v1"
AGENT2_MODULE = (
    "openai_ns_reconstruction."
    "kokuno_strict_inner_leading_oscillatory_spatial_candidate"
)
AGENT2_CLASS = "KokunoStrictInnerLeadingOscillatorySpatialCandidate"
AGENT2_PARENT_PR = 866
AGENT2_PARENT_HEAD = "68f84128f07b6743d368a9ab7a051e441f5a59b2"

# Exact-head Actions were queued at the freeze. Queued is not PASS.
AGENT2_EXACT_HEAD_CI_CONCLUSION: str | None = None

# No Agent-4 PR currently targets #874's saved/reloaded spatial artifact.
# A2's own public-velocity FD4 comparison is an engineering verifier, not the
# independent validator authority required for scientific admission.
AGENT4_DEDICATED_SPATIAL_ARTIFACT_AUDIT_PRESENT = False
AGENT4_EXACT_HEAD_CI_CONCLUSION: str | None = None
AGENT4_INDEPENDENT_AUDIT_CONCLUSION: str | None = None

AGENT2_PROTOCOL = {
    "scope": "strict-inner PA10 contraction-center plus frozen complete-curl oscillation",
    "velocity_surface": "velocity(x,y,z,t)->[...,3]",
    "velocity_dt_surface": "velocity_dt(x,y,z,t)->[...,3]",
    "velocity_jacobian_surface": "velocity_jacobian(x,y,z,t)->[...,3,3]",
    "velocity_candidate_sha256_preserved": True,
    "differentiable_sha256_preserved": True,
    "spatial_sha256_bound": True,
    "manifest_sha256_bound": True,
    "save_manifest_present": True,
    "load_manifest_present": True,
    "inner_jacobian_semantics": "analytic Agent1 PA10 strict-inner Cartesian Jacobian",
    "oscillatory_jacobian_operator": "centered_cartesian_fd6",
    "oscillatory_jacobian_fixed_step": 1.0e-3,
    "caller_spatial_step_knob_present": False,
    "jacobian_composition": "J_total=J_inner_strict+J_osc_complete_curl",
    "construction_side_fd4_sample_count": 12,
    "construction_side_fd4_source_X_range": [0.166, 0.301],
    "construction_side_fd4_steps": [4.0e-3, 2.0e-3, 1.0e-3],
    "construction_side_fd4_refinement_ratio_gate": 6.0,
    "construction_side_fd4_refinement_floor": 2.0e-10,
    "construction_side_fd4_fine_relative_rms_gate": 5.0e-3,
    "construction_side_fd4_fine_relative_sampled_max_gate": 1.0e-2,
    "composite_jacobian_rms_floor": 1.0e-10,
    "oscillatory_jacobian_rms_floor": 1.0e-10,
    "jacobian_additive_closure_gate": 1.0e-13,
    "batch_scalar_jacobian_replay_gate": 1.0e-10,
    "manifest_rebound_velocity_exact": True,
    "manifest_rebound_velocity_dt_exact": True,
    "manifest_rebound_velocity_jacobian_exact": True,
    "agent1_inner_domain_failure_propagates": True,
    "outer_taper_or_zero_extension_invented": False,
    "post_observation_retuning_allowed": False,
}

NORM_SCOPE_FIREWALL = {
    "agent2_fd4_is_independent_agent4_validation": False,
    "jacobian_consistency_is_complete_ns_residual": False,
    "jacobian_consistency_is_cr001_momentum_norm": False,
    "derived_divergence_diagnostic_is_final_divergence_gate": False,
    "derived_vorticity_diagnostic_is_pde_acceptance_metric": False,
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


def deterministic_strict_inner_spatial_artifact_ingest_contract(
    *, exact_head: str | None = None
) -> dict[str, Any]:
    parent = deterministic_strict_inner_differentiable_artifact_ingest_contract(
        exact_head=PARENT_A5_HEAD
    )
    validate_strict_inner_differentiable_artifact_ingest_contract(parent)

    admitted = bool(
        parent["ingest_status"]["strict_inner_differentiable_artifact_ingest_admitted"]
        and AGENT2_EXACT_HEAD_CI_CONCLUSION == "success"
        and AGENT4_DEDICATED_SPATIAL_ARTIFACT_AUDIT_PRESENT
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
            "registered_scoped_stage": "strict_inner_spatial_candidate_artifact",
            "velocity_surface_registered": True,
            "velocity_dt_surface_registered": True,
            "velocity_jacobian_surface_registered": True,
            "full_pipeline_candidate_artifact_stage_reached": False,
            "reason": (
                "artifact remains strict-inner, has no dedicated Agent-4 spatial-artifact audit, "
                "and lacks corrected/global leading join, matched pressure, restricted forcing "
                "and Agent-3 correction velocity"
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
            "parent_pr": AGENT2_PARENT_PR,
            "parent_head": AGENT2_PARENT_HEAD,
            "protocol": dict(AGENT2_PROTOCOL),
            "norm_scope_firewall": dict(NORM_SCOPE_FIREWALL),
        },
        "agent4_binding": {
            "dedicated_spatial_artifact_audit_present": (
                AGENT4_DEDICATED_SPATIAL_ARTIFACT_AUDIT_PRESENT
            ),
            "authority": None,
            "audited_agent2_pr": None,
            "audited_agent2_head": None,
            "exact_head_ci_conclusion": AGENT4_EXACT_HEAD_CI_CONCLUSION,
            "independent_audit_conclusion": AGENT4_INDEPENDENT_AUDIT_CONCLUSION,
            "reason": (
                "latest matching independent artifact audit is not yet present; "
                "Agent4#869 audits the #866 velocity_dt seam, not Agent2#874 velocity_jacobian"
            ),
        },
        "evidence": {
            "parent_differentiable_artifact_ingest_admitted": bool(
                parent["ingest_status"]["strict_inner_differentiable_artifact_ingest_admitted"]
            ),
            "agent2_exact_head_ci_conclusion": AGENT2_EXACT_HEAD_CI_CONCLUSION,
            "agent2_construction_side_fd4_verifier_present": True,
            "agent4_dedicated_spatial_artifact_audit_present": (
                AGENT4_DEDICATED_SPATIAL_ARTIFACT_AUDIT_PRESENT
            ),
            "agent4_exact_head_ci_conclusion": AGENT4_EXACT_HEAD_CI_CONCLUSION,
            "agent4_independent_audit_conclusion": AGENT4_INDEPENDENT_AUDIT_CONCLUSION,
            "agent4_scientific_receipt_admitted": False,
            "strict_inner_spatial_artifact_scientifically_admitted": False,
        },
        "spatial_artifact_handoff": {
            "registered": True,
            "status": "registered_waiting_for_independent_a4_audit",
            "construction_authority": "Agent2#874",
            "independent_audit_available": False,
            "independent_audit_authority": None,
            "public_velocity_binding": (
                "Agent2#874.KokunoStrictInnerLeadingOscillatorySpatialCandidate.velocity"
            ),
            "public_velocity_dt_binding": (
                "Agent2#874.KokunoStrictInnerLeadingOscillatorySpatialCandidate.velocity_dt"
            ),
            "public_velocity_jacobian_binding": (
                "Agent2#874.KokunoStrictInnerLeadingOscillatorySpatialCandidate.velocity_jacobian"
            ),
            "save_manifest_binding": (
                "Agent2#874.KokunoStrictInnerLeadingOscillatorySpatialCandidate.save_manifest"
            ),
            "load_manifest_binding": (
                "Agent2#874.KokunoStrictInnerLeadingOscillatorySpatialCandidate.load_manifest"
            ),
            "base_velocity_candidate_sha256_preserved": True,
            "differentiable_sha256_preserved": True,
            "spatial_sha256_registered": True,
            "manifest_sha256_registered": True,
            "python_velocity_velocity_dt_and_jacobian_callable_surface_registered": True,
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
            "strict_inner_artifact_velocity": "Agent2#874.velocity",
            "strict_inner_artifact_velocity_dt": "Agent2#874.velocity_dt",
            "strict_inner_artifact_velocity_jacobian": "Agent2#874.velocity_jacobian",
            "pressure": None,
            "forcing": None,
            "complete_candidate_api_ready": False,
        },
        "artifact_status": {
            **dict(parent["artifact_status"]),
            "upstream_strict_inner_spatial_artifact_implementation_present": True,
            "upstream_velocity_jacobian_surface_present": True,
            "upstream_spatial_manifest_save_load_surface_present": True,
            "upstream_independent_spatial_artifact_audit_present": False,
            "strict_inner_spatial_artifact_ingest_admitted": admitted,
            "complete_global_candidate_artifact_materialized": False,
        },
        "ingest_status": {
            **dict(parent["ingest_status"]),
            "typed_strict_inner_spatial_artifact_registered": True,
            "typed_strict_inner_spatial_artifact_independent_audit_registered": False,
            "strict_inner_spatial_artifact_ingest_admitted": admitted,
            "full_candidate_artifact_ready": False,
        },
        "stage_state": dict(parent["stage_state"]),
        "baseline_vs_kokuno": dict(parent["baseline_vs_kokuno"]),
        "final_project_gates_unchanged": dict(parent["final_project_gates_unchanged"]),
        "truth_boundary": {
            "strict_inner_spatial_artifact_surface_registered": True,
            "agent2_construction_side_fd4_verifier_registered": True,
            "registration_not_scientific_admission": True,
            "strict_inner_scope_preserved": True,
            "agent4_spatial_artifact_audit_surface_registered": False,
            "agent2_ci_promoted": False,
            "agent2_fd4_laundered_as_independent_a4_validation": False,
            "agent4_receipt_invented": False,
            "strict_inner_spatial_artifact_independently_admitted": False,
            "jacobian_consistency_laundered_as_pde_residual": False,
            "derived_divergence_laundered_as_final_divergence_gate": False,
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


def validate_strict_inner_spatial_artifact_ingest_contract(
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
    expected = deterministic_strict_inner_spatial_artifact_ingest_contract(
        exact_head=exact_head
    )
    if dict(payload) != expected:
        raise ValueError("strict-inner spatial artifact ingest contract drifted")
    if payload.get("contract_sha256") != _digest(payload):
        raise ValueError("contract sha256 mismatch")


def write_report(path: str | Path, *, exact_head: str | None = None) -> dict[str, Any]:
    payload = deterministic_strict_inner_spatial_artifact_ingest_contract(
        exact_head=exact_head
    )
    validate_strict_inner_spatial_artifact_ingest_contract(payload)
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
