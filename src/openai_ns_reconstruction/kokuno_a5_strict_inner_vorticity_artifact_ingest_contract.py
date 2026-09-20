"""Fail-closed A5 registration of the strict-inner vorticity candidate artifact.

Agent 2 PR #881 extends the checksum-bound strict-inner ``u_inner + u_osc``
artifact from public ``velocity`` / ``velocity_dt`` / ``velocity_jacobian`` to a
first-class public ``vorticity`` surface.  Production vorticity is the curl of
the already-frozen #874 spatial Jacobian, so this increment introduces no new
production derivative realization.

Agent 4 PR #883 independently saves/reloads that artifact and reconstructs curl
from only the reloaded public ``velocity`` using a centered Cartesian FD6 path,
distinct from Agent-1 analytic derivatives, Agent-2 production curl assembly,
and Agent-2's construction-side FD4 verifier.

This module registers that construction/audit handoff without copying sibling
mathematics.  Exact-head Actions are unresolved at the freeze, so registration
is not scientific admission.  The artifact remains strict-inner and lacks the
corrected/global leading join, Agent-3 correction velocity, matched pressure,
and preregistered restricted forcing required for a complete Navier--Stokes
candidate.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from .kokuno_a5_strict_inner_spatial_artifact_ingest_contract import (
    deterministic_strict_inner_spatial_artifact_ingest_contract,
    validate_strict_inner_spatial_artifact_ingest_contract,
)

SCHEMA = "kokuno-a5-strict-inner-vorticity-artifact-ingest-contract-v1"
TASK = "KOKUNO-A5-STRICT-INNER-VORTICITY-ARTIFACT-INGEST-080"
PARENT_A5_PR = 877
PARENT_A5_HEAD = "550a7039390979f2705bb75da98b43cc7259af02"

AGENT2_PR = 881
AGENT2_HEAD = "ded8c88dd03b784368cf1a30601faf76a3682c64"
AGENT2_SOURCE_BLOB = "9dd3b7400fea8f2d4462cda86561103e29a1d053"
AGENT2_WORKFLOW_BLOB = "024069383d9390e54e101529332bfd2707040659"
AGENT2_DEDICATED_RUN = 35526134983
AGENT2_TESTS_RUN = 35526134949
AGENT2_SCHEMA = "kokuno-a2-strict-inner-vorticity-candidate-v1"
AGENT2_MODULE = (
    "openai_ns_reconstruction."
    "kokuno_strict_inner_leading_oscillatory_vorticity_candidate"
)
AGENT2_CLASS = "KokunoStrictInnerLeadingOscillatoryVorticityCandidate"
AGENT2_PARENT_PR = 874
AGENT2_PARENT_HEAD = "9cb3869b9cfd2d5b8dcb1da222df74d80e12d0c2"

AGENT4_PR = 883
AGENT4_HEAD = "772acaced47beb8d6d7cc4745833cbde06f10940"
AGENT4_SOURCE_BLOB = "c0121785a487771711f6940426025d25731b4b3a"
AGENT4_WORKFLOW_BLOB = "d47c069b81e6052fc112a3e2f3a1a3b32e9558a3"
AGENT4_DEDICATED_RUN = 35526934729
AGENT4_TESTS_RUN = 35526934711
AGENT4_SCHEMA = "kokuno-a4-strict-inner-vorticity-artifact-independent-audit-v1"
AGENT4_MODULE = (
    "openai_ns_reconstruction."
    "kokuno_a4_strict_inner_vorticity_artifact_independent_audit"
)

# Exact-head Actions were queued at this freeze. Queued is not PASS.
AGENT2_EXACT_HEAD_CI_CONCLUSION: str | None = None
AGENT4_EXACT_HEAD_CI_CONCLUSION: str | None = None
AGENT4_INDEPENDENT_AUDIT_CONCLUSION: str | None = None

AGENT2_PROTOCOL = {
    "scope": "strict-inner PA10 contraction-center plus frozen complete-curl oscillation",
    "velocity_surface": "velocity(x,y,z,t)->[...,3]",
    "velocity_dt_surface": "velocity_dt(x,y,z,t)->[...,3]",
    "velocity_jacobian_surface": "velocity_jacobian(x,y,z,t)->[...,3,3]",
    "vorticity_surface": "vorticity(x,y,z,t)->[...,3]",
    "vorticity_definition": (
        "curl from inherited spatial Jacobian with J[component,axis]=partial_axis velocity_component"
    ),
    "new_production_derivative_realization_introduced": False,
    "velocity_candidate_sha256_preserved": True,
    "differentiable_sha256_preserved": True,
    "spatial_sha256_preserved": True,
    "vorticity_sha256_bound": True,
    "manifest_sha256_bound": True,
    "save_manifest_present": True,
    "load_manifest_present": True,
    "construction_side_reference": "public-composite-velocity-only-centered-fd4-curl",
    "construction_side_fd4_sample_count": 12,
    "construction_side_fd4_source_X_range": [0.168, 0.298],
    "construction_side_fd4_steps": [4.0e-3, 2.0e-3, 1.0e-3],
    "construction_side_fd4_refinement_ratio_gate": 6.0,
    "construction_side_fd4_refinement_floor": 2.0e-10,
    "construction_side_fd4_fine_relative_rms_gate": 5.0e-3,
    "construction_side_fd4_fine_relative_sampled_max_gate": 1.0e-2,
    "production_vorticity_rms_floor": 1.0e-10,
    "batch_scalar_vorticity_replay_gate": 1.0e-10,
    "manifest_rebound_vorticity_exact": True,
    "three_resolution_morphology_verified": False,
    "whole_domain_vorticity_morphology_verified": False,
    "agent1_inner_domain_failure_propagates": True,
    "outer_taper_or_zero_extension_invented": False,
    "post_observation_retuning_allowed": False,
}

AGENT2_NORM_SCOPE_FIREWALL = {
    "agent2_fd4_is_independent_agent4_validation": False,
    "vorticity_consistency_is_complete_ns_residual": False,
    "vorticity_consistency_is_cr001_momentum_norm": False,
    "derived_divergence_diagnostic_is_final_divergence_gate": False,
    "vorticity_consistency_is_whole_domain_morphology_validation": False,
    "whole_domain_volume_weighted_l2_assessed": False,
    "whole_domain_boundary_support_assessed": False,
    "same_protocol_st006_comparison_legal": False,
    "formal_full_domain_pde_gate_assessed": False,
}

AGENT4_PROTOCOL = {
    "independent_reference": "saved-reloaded-public-velocity-only-centered-cartesian-fd6-curl",
    "seed": 9173551,
    "random_offgrid_sample_count": 128,
    "exact_axis_probe_count": 3,
    "axis_near_probe_count": 3,
    "source_X_range": [0.18, 0.27],
    "source_eta_range": [-0.22, 0.22],
    "time_range": [0.47, 0.53],
    "fd6_spatial_steps": [3.0e-3, 1.5e-3, 7.5e-4],
    "fine_relative_rms_gate": 5.0e-3,
    "fine_relative_sampled_max_gate": 1.0e-2,
    "refinement_ratio_gate": 12.0,
    "refinement_floor": 2.0e-10,
    "production_vorticity_rms_floor": 1.0e-10,
    "public_velocity_rms_floor": 1.0e-8,
    "independent_divergence_sampled_max_gate": 1.0e-5,
    "independent_divergence_sampled_rms_gate": 1.0e-5,
    "negative_controls": [
        "production-vorticity-times-0.90",
        "largest-rms-vorticity-component-sign-flip",
        "vorticity-sha-mutation-with-refreshed-manifest-sha",
        "momentum-gate-1e-3-to-1.001e-3-with-refreshed-manifest-sha",
    ],
    "time_plus_1e_minus_4_response_is_observation_only": True,
    "uses_agent1_analytic_spatial_derivatives_as_reference": False,
    "uses_agent2_production_jacobian_as_reference": False,
    "uses_agent2_production_curl_as_reference": False,
    "uses_agent2_fd4_verifier_as_reference": False,
    "uses_candidate_component_tensors_as_reference": False,
    "post_observation_retuning_allowed": False,
}

AGENT4_NORM_SCOPE_FIREWALL = {
    "vorticity_consistency_is_complete_ns_residual": False,
    "vorticity_consistency_is_cr001_momentum_norm": False,
    "vorticity_consistency_is_whole_domain_morphology_validation": False,
    "divergence_rms_is_heldout_sampled_rms": True,
    "divergence_rms_is_whole_domain_volume_weighted_l2": False,
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


def deterministic_strict_inner_vorticity_artifact_ingest_contract(
    *, exact_head: str | None = None
) -> dict[str, Any]:
    parent = deterministic_strict_inner_spatial_artifact_ingest_contract(
        exact_head=PARENT_A5_HEAD
    )
    validate_strict_inner_spatial_artifact_ingest_contract(parent)

    admitted = bool(
        parent["ingest_status"]["strict_inner_spatial_artifact_ingest_admitted"]
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
            "registered_scoped_stage": "strict_inner_vorticity_candidate_artifact",
            "velocity_surface_registered": True,
            "velocity_dt_surface_registered": True,
            "velocity_jacobian_surface_registered": True,
            "vorticity_surface_registered": True,
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
            "workflow_blob_sha": AGENT2_WORKFLOW_BLOB,
            "dedicated_run": AGENT2_DEDICATED_RUN,
            "tests_run": AGENT2_TESTS_RUN,
            "schema": AGENT2_SCHEMA,
            "module": AGENT2_MODULE,
            "class": AGENT2_CLASS,
            "parent_pr": AGENT2_PARENT_PR,
            "parent_head": AGENT2_PARENT_HEAD,
            "protocol": dict(AGENT2_PROTOCOL),
            "norm_scope_firewall": dict(AGENT2_NORM_SCOPE_FIREWALL),
        },
        "agent4_binding": {
            "pr": AGENT4_PR,
            "head": AGENT4_HEAD,
            "source_blob_sha": AGENT4_SOURCE_BLOB,
            "workflow_blob_sha": AGENT4_WORKFLOW_BLOB,
            "dedicated_run": AGENT4_DEDICATED_RUN,
            "tests_run": AGENT4_TESTS_RUN,
            "schema": AGENT4_SCHEMA,
            "module": AGENT4_MODULE,
            "dedicated_vorticity_artifact_audit_present": True,
            "authority": "Agent4#883",
            "audited_agent2_pr": AGENT2_PR,
            "audited_agent2_head": AGENT2_HEAD,
            "protocol": dict(AGENT4_PROTOCOL),
            "norm_scope_firewall": dict(AGENT4_NORM_SCOPE_FIREWALL),
        },
        "evidence": {
            "parent_spatial_artifact_ingest_admitted": bool(
                parent["ingest_status"]["strict_inner_spatial_artifact_ingest_admitted"]
            ),
            "agent2_exact_head_ci_conclusion": AGENT2_EXACT_HEAD_CI_CONCLUSION,
            "agent2_construction_side_fd4_verifier_present": True,
            "agent4_dedicated_vorticity_artifact_audit_present": True,
            "agent4_exact_head_ci_conclusion": AGENT4_EXACT_HEAD_CI_CONCLUSION,
            "agent4_independent_audit_conclusion": AGENT4_INDEPENDENT_AUDIT_CONCLUSION,
            "agent4_scientific_receipt_admitted": False,
            "strict_inner_vorticity_artifact_scientifically_admitted": False,
        },
        "vorticity_artifact_handoff": {
            "registered": True,
            "status": "registered_unresolved",
            "construction_authority": "Agent2#881",
            "independent_audit_available": True,
            "independent_audit_authority": "Agent4#883",
            "public_velocity_binding": (
                "Agent2#881.KokunoStrictInnerLeadingOscillatoryVorticityCandidate.velocity"
            ),
            "public_velocity_dt_binding": (
                "Agent2#881.KokunoStrictInnerLeadingOscillatoryVorticityCandidate.velocity_dt"
            ),
            "public_velocity_jacobian_binding": (
                "Agent2#881.KokunoStrictInnerLeadingOscillatoryVorticityCandidate.velocity_jacobian"
            ),
            "public_vorticity_binding": (
                "Agent2#881.KokunoStrictInnerLeadingOscillatoryVorticityCandidate.vorticity"
            ),
            "save_manifest_binding": (
                "Agent2#881.KokunoStrictInnerLeadingOscillatoryVorticityCandidate.save_manifest"
            ),
            "load_manifest_binding": (
                "Agent2#881.KokunoStrictInnerLeadingOscillatoryVorticityCandidate.load_manifest"
            ),
            "base_velocity_candidate_sha256_preserved": True,
            "differentiable_sha256_preserved": True,
            "spatial_sha256_preserved": True,
            "vorticity_sha256_registered": True,
            "manifest_sha256_registered": True,
            "python_velocity_velocity_dt_jacobian_vorticity_callable_surface_registered": True,
            "matlab_export_materialized": False,
            "pressure_in_artifact_surface": False,
            "restricted_forcing_in_artifact_surface": False,
            "agent3_correction_velocity_in_artifact": False,
            "outer_global_join_in_artifact": False,
            "complete_candidate_artifact": False,
            "eligible_for_velocity_export_ready": False,
            "usable_as_complete_ns_defect_input": False,
            "usable_for_final_independent_pde_validation": False,
            "usable_as_whole_domain_vorticity_morphology_evidence": False,
        },
        "candidate_api_handoff": {
            **dict(parent["candidate_api_handoff"]),
            "strict_inner_artifact_velocity": "Agent2#881.velocity",
            "strict_inner_artifact_velocity_dt": "Agent2#881.velocity_dt",
            "strict_inner_artifact_velocity_jacobian": "Agent2#881.velocity_jacobian",
            "strict_inner_artifact_vorticity": "Agent2#881.vorticity",
            "pressure": None,
            "forcing": None,
            "complete_candidate_api_ready": False,
        },
        "artifact_status": {
            **dict(parent["artifact_status"]),
            "upstream_strict_inner_vorticity_artifact_implementation_present": True,
            "upstream_vorticity_surface_present": True,
            "upstream_vorticity_manifest_save_load_surface_present": True,
            "upstream_independent_vorticity_artifact_audit_present": True,
            "strict_inner_vorticity_artifact_ingest_admitted": admitted,
            "complete_global_candidate_artifact_materialized": False,
        },
        "ingest_status": {
            **dict(parent["ingest_status"]),
            "typed_strict_inner_vorticity_artifact_registered": True,
            "typed_strict_inner_vorticity_artifact_independent_audit_registered": True,
            "strict_inner_vorticity_artifact_ingest_admitted": admitted,
            "full_candidate_artifact_ready": False,
        },
        "stage_state": dict(parent["stage_state"]),
        "baseline_vs_kokuno": dict(parent["baseline_vs_kokuno"]),
        "final_project_gates_unchanged": dict(parent["final_project_gates_unchanged"]),
        "truth_boundary": {
            "strict_inner_vorticity_artifact_surface_registered": True,
            "agent2_construction_side_fd4_verifier_registered": True,
            "agent4_vorticity_artifact_audit_surface_registered": True,
            "registration_not_scientific_admission": True,
            "strict_inner_scope_preserved": True,
            "three_resolution_morphology_verified": False,
            "whole_domain_vorticity_morphology_verified": False,
            "agent2_ci_promoted": False,
            "agent4_ci_promoted": False,
            "agent2_fd4_laundered_as_independent_a4_validation": False,
            "agent4_receipt_invented": False,
            "strict_inner_vorticity_artifact_independently_admitted": False,
            "vorticity_consistency_laundered_as_pde_residual": False,
            "vorticity_consistency_laundered_as_whole_domain_morphology": False,
            "sampled_divergence_laundered_as_whole_domain_l2": False,
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


def validate_strict_inner_vorticity_artifact_ingest_contract(
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
    expected = deterministic_strict_inner_vorticity_artifact_ingest_contract(
        exact_head=exact_head
    )
    if dict(payload) != expected:
        raise ValueError("strict-inner vorticity artifact ingest contract drifted")
    if payload.get("contract_sha256") != _digest(payload):
        raise ValueError("contract sha256 mismatch")


def write_report(path: str | Path, *, exact_head: str | None = None) -> dict[str, Any]:
    payload = deterministic_strict_inner_vorticity_artifact_ingest_contract(
        exact_head=exact_head
    )
    validate_strict_inner_vorticity_artifact_ingest_contract(payload)
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
