from __future__ import annotations

import copy

import pytest

from openai_ns_reconstruction.kokuno_a5_strict_inner_spatial_artifact_ingest_contract import (
    deterministic_strict_inner_spatial_artifact_ingest_contract,
)
from openai_ns_reconstruction.kokuno_a5_strict_inner_vorticity_artifact_ingest_contract import (
    AGENT2_HEAD,
    AGENT2_SOURCE_BLOB,
    AGENT2_WORKFLOW_BLOB,
    AGENT4_HEAD,
    AGENT4_SOURCE_BLOB,
    AGENT4_WORKFLOW_BLOB,
    deterministic_strict_inner_vorticity_artifact_ingest_contract,
    validate_strict_inner_vorticity_artifact_ingest_contract,
)


EXACT_HEAD = "0" * 40
PARENT_HEAD = "550a7039390979f2705bb75da98b43cc7259af02"


def _payload():
    return deterministic_strict_inner_vorticity_artifact_ingest_contract(
        exact_head=EXACT_HEAD
    )


def test_vorticity_artifact_and_a4_audit_are_registered_but_unresolved() -> None:
    payload = _payload()
    validate_strict_inner_vorticity_artifact_ingest_contract(payload)

    a2 = payload["agent2_binding"]
    assert a2["head"] == AGENT2_HEAD
    assert a2["source_blob_sha"] == AGENT2_SOURCE_BLOB
    assert a2["workflow_blob_sha"] == AGENT2_WORKFLOW_BLOB
    assert a2["schema"] == "kokuno-a2-strict-inner-vorticity-candidate-v1"
    assert a2["class"] == "KokunoStrictInnerLeadingOscillatoryVorticityCandidate"
    assert a2["parent_pr"] == 874

    a4 = payload["agent4_binding"]
    assert a4["head"] == AGENT4_HEAD
    assert a4["source_blob_sha"] == AGENT4_SOURCE_BLOB
    assert a4["workflow_blob_sha"] == AGENT4_WORKFLOW_BLOB
    assert a4["schema"] == "kokuno-a4-strict-inner-vorticity-artifact-independent-audit-v1"
    assert a4["dedicated_vorticity_artifact_audit_present"] is True
    assert a4["authority"] == "Agent4#883"
    assert a4["audited_agent2_pr"] == 881
    assert a4["audited_agent2_head"] == AGENT2_HEAD

    evidence = payload["evidence"]
    assert evidence["parent_spatial_artifact_ingest_admitted"] is False
    assert evidence["agent2_exact_head_ci_conclusion"] is None
    assert evidence["agent2_construction_side_fd4_verifier_present"] is True
    assert evidence["agent4_dedicated_vorticity_artifact_audit_present"] is True
    assert evidence["agent4_exact_head_ci_conclusion"] is None
    assert evidence["agent4_independent_audit_conclusion"] is None
    assert evidence["agent4_scientific_receipt_admitted"] is False
    assert evidence["strict_inner_vorticity_artifact_scientifically_admitted"] is False

    handoff = payload["vorticity_artifact_handoff"]
    assert handoff["registered"] is True
    assert handoff["status"] == "registered_unresolved"
    assert handoff["construction_authority"] == "Agent2#881"
    assert handoff["independent_audit_authority"] == "Agent4#883"
    assert handoff["public_velocity_binding"].endswith(".velocity")
    assert handoff["public_velocity_dt_binding"].endswith(".velocity_dt")
    assert handoff["public_velocity_jacobian_binding"].endswith(".velocity_jacobian")
    assert handoff["public_vorticity_binding"].endswith(".vorticity")
    assert handoff["vorticity_sha256_registered"] is True
    assert handoff["manifest_sha256_registered"] is True
    assert handoff["pressure_in_artifact_surface"] is False
    assert handoff["restricted_forcing_in_artifact_surface"] is False
    assert handoff["agent3_correction_velocity_in_artifact"] is False
    assert handoff["outer_global_join_in_artifact"] is False
    assert handoff["complete_candidate_artifact"] is False
    assert handoff["eligible_for_velocity_export_ready"] is False
    assert handoff["usable_as_whole_domain_vorticity_morphology_evidence"] is False

    assert payload["ingest_status"]["typed_strict_inner_vorticity_artifact_registered"] is True
    assert payload["ingest_status"]["typed_strict_inner_vorticity_artifact_independent_audit_registered"] is True
    assert payload["ingest_status"]["strict_inner_vorticity_artifact_ingest_admitted"] is False
    assert payload["stage_state"] == {
        "leading_ready": False,
        "oscillatory_ready": True,
        "correction_ready": False,
        "velocity_export_ready": False,
        "pde_validated": False,
    }


def test_vorticity_construction_and_independent_audit_keep_scope_separate() -> None:
    payload = _payload()

    a2 = payload["agent2_binding"]["protocol"]
    assert a2["vorticity_surface"] == "vorticity(x,y,z,t)->[...,3]"
    assert a2["new_production_derivative_realization_introduced"] is False
    assert a2["construction_side_reference"] == "public-composite-velocity-only-centered-fd4-curl"
    assert a2["construction_side_fd4_steps"] == [4e-3, 2e-3, 1e-3]
    assert a2["construction_side_fd4_refinement_ratio_gate"] == 6.0
    assert a2["construction_side_fd4_fine_relative_rms_gate"] == 5e-3
    assert a2["construction_side_fd4_fine_relative_sampled_max_gate"] == 1e-2
    assert a2["three_resolution_morphology_verified"] is False
    assert a2["whole_domain_vorticity_morphology_verified"] is False

    a2_firewall = payload["agent2_binding"]["norm_scope_firewall"]
    assert a2_firewall["agent2_fd4_is_independent_agent4_validation"] is False
    assert a2_firewall["vorticity_consistency_is_complete_ns_residual"] is False
    assert a2_firewall["vorticity_consistency_is_cr001_momentum_norm"] is False
    assert a2_firewall["vorticity_consistency_is_whole_domain_morphology_validation"] is False
    assert a2_firewall["whole_domain_volume_weighted_l2_assessed"] is False

    a4 = payload["agent4_binding"]["protocol"]
    assert a4["seed"] == 9173551
    assert a4["random_offgrid_sample_count"] == 128
    assert a4["exact_axis_probe_count"] == 3
    assert a4["axis_near_probe_count"] == 3
    assert a4["fd6_spatial_steps"] == [3e-3, 1.5e-3, 7.5e-4]
    assert a4["fine_relative_rms_gate"] == 5e-3
    assert a4["fine_relative_sampled_max_gate"] == 1e-2
    assert a4["refinement_ratio_gate"] == 12.0
    assert a4["independent_divergence_sampled_max_gate"] == 1e-5
    assert a4["independent_divergence_sampled_rms_gate"] == 1e-5
    assert a4["uses_agent1_analytic_spatial_derivatives_as_reference"] is False
    assert a4["uses_agent2_production_jacobian_as_reference"] is False
    assert a4["uses_agent2_production_curl_as_reference"] is False
    assert a4["uses_agent2_fd4_verifier_as_reference"] is False

    a4_firewall = payload["agent4_binding"]["norm_scope_firewall"]
    assert a4_firewall["vorticity_consistency_is_complete_ns_residual"] is False
    assert a4_firewall["vorticity_consistency_is_whole_domain_morphology_validation"] is False
    assert a4_firewall["divergence_rms_is_heldout_sampled_rms"] is True
    assert a4_firewall["divergence_rms_is_whole_domain_volume_weighted_l2"] is False
    assert a4_firewall["same_protocol_st006_comparison_legal"] is False

    api = payload["candidate_api_handoff"]
    assert api["strict_inner_artifact_velocity"] == "Agent2#881.velocity"
    assert api["strict_inner_artifact_velocity_dt"] == "Agent2#881.velocity_dt"
    assert api["strict_inner_artifact_velocity_jacobian"] == "Agent2#881.velocity_jacobian"
    assert api["strict_inner_artifact_vorticity"] == "Agent2#881.vorticity"
    assert api["pressure"] is None
    assert api["forcing"] is None
    assert api["complete_candidate_api_ready"] is False

    parent = deterministic_strict_inner_spatial_artifact_ingest_contract(
        exact_head=PARENT_HEAD
    )
    assert payload["final_project_gates_unchanged"] == parent["final_project_gates_unchanged"]
    assert payload["baseline_vs_kokuno"] == parent["baseline_vs_kokuno"]
    truth = payload["truth_boundary"]
    assert truth["three_resolution_morphology_verified"] is False
    assert truth["whole_domain_vorticity_morphology_verified"] is False
    assert truth["vorticity_consistency_laundered_as_pde_residual"] is False
    assert truth["vorticity_consistency_laundered_as_whole_domain_morphology"] is False
    assert truth["sampled_divergence_laundered_as_whole_domain_l2"] is False
    assert truth["free_residual_defined_forcing_allowed"] is False
    assert truth["threshold_relaxed"] is False


@pytest.mark.parametrize(
    "mutate",
    [
        lambda p: p["evidence"].__setitem__("agent2_exact_head_ci_conclusion", "success"),
        lambda p: p["evidence"].__setitem__("agent4_exact_head_ci_conclusion", "success"),
        lambda p: p["evidence"].__setitem__("agent4_independent_audit_conclusion", "pass"),
        lambda p: p["agent2_binding"].__setitem__("head", "1" * 40),
        lambda p: p["agent2_binding"].__setitem__("source_blob_sha", "2" * 40),
        lambda p: p["agent2_binding"].__setitem__("workflow_blob_sha", "3" * 40),
        lambda p: p["agent4_binding"].__setitem__("head", "4" * 40),
        lambda p: p["agent4_binding"].__setitem__("source_blob_sha", "5" * 40),
        lambda p: p["agent4_binding"].__setitem__("workflow_blob_sha", "6" * 40),
        lambda p: p["agent2_binding"]["protocol"].__setitem__("three_resolution_morphology_verified", True),
        lambda p: p["agent2_binding"]["norm_scope_firewall"].__setitem__("agent2_fd4_is_independent_agent4_validation", True),
        lambda p: p["agent4_binding"]["protocol"].__setitem__("refinement_ratio_gate", 6.0),
        lambda p: p["agent4_binding"]["norm_scope_firewall"].__setitem__("divergence_rms_is_whole_domain_volume_weighted_l2", True),
        lambda p: p["vorticity_artifact_handoff"].__setitem__("pressure_in_artifact_surface", True),
        lambda p: p["vorticity_artifact_handoff"].__setitem__("restricted_forcing_in_artifact_surface", True),
        lambda p: p["vorticity_artifact_handoff"].__setitem__("agent3_correction_velocity_in_artifact", True),
        lambda p: p["vorticity_artifact_handoff"].__setitem__("outer_global_join_in_artifact", True),
        lambda p: p["vorticity_artifact_handoff"].__setitem__("complete_candidate_artifact", True),
        lambda p: p["vorticity_artifact_handoff"].__setitem__("eligible_for_velocity_export_ready", True),
        lambda p: p["vorticity_artifact_handoff"].__setitem__("usable_as_whole_domain_vorticity_morphology_evidence", True),
        lambda p: p["stage_state"].__setitem__("leading_ready", True),
        lambda p: p["stage_state"].__setitem__("correction_ready", True),
        lambda p: p["stage_state"].__setitem__("velocity_export_ready", True),
        lambda p: p["stage_state"].__setitem__("pde_validated", True),
        lambda p: p["truth_boundary"].__setitem__("whole_domain_vorticity_morphology_verified", True),
        lambda p: p["truth_boundary"].__setitem__("free_residual_defined_forcing_allowed", True),
        lambda p: p["truth_boundary"].__setitem__("threshold_relaxed", True),
    ],
)
def test_truth_boundary_mutations_fail_closed(mutate) -> None:
    payload = copy.deepcopy(_payload())
    mutate(payload)
    with pytest.raises(ValueError, match="drifted"):
        validate_strict_inner_vorticity_artifact_ingest_contract(payload)


def test_invalid_exact_head_fails_closed() -> None:
    payload = deterministic_strict_inner_vorticity_artifact_ingest_contract(
        exact_head=None
    )
    payload["exact_head"] = "not-a-commit"
    with pytest.raises(ValueError, match="40-character"):
        validate_strict_inner_vorticity_artifact_ingest_contract(payload)
