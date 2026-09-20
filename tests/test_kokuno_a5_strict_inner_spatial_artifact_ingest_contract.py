from __future__ import annotations

import copy

import pytest

from openai_ns_reconstruction.kokuno_a5_strict_inner_differentiable_artifact_ingest_contract import (
    deterministic_strict_inner_differentiable_artifact_ingest_contract,
)
from openai_ns_reconstruction.kokuno_a5_strict_inner_spatial_artifact_ingest_contract import (
    AGENT2_HEAD,
    AGENT2_SOURCE_BLOB,
    deterministic_strict_inner_spatial_artifact_ingest_contract,
    validate_strict_inner_spatial_artifact_ingest_contract,
)


EXACT_HEAD = "0" * 40
PARENT_HEAD = "b155d5a170ad56bb5a2dde9d5f4e0f43d1e2b474"


def _payload():
    return deterministic_strict_inner_spatial_artifact_ingest_contract(
        exact_head=EXACT_HEAD
    )


def test_spatial_artifact_is_registered_but_waits_for_agent4() -> None:
    payload = _payload()
    validate_strict_inner_spatial_artifact_ingest_contract(payload)

    assert payload["agent2_binding"]["head"] == AGENT2_HEAD
    assert payload["agent2_binding"]["source_blob_sha"] == AGENT2_SOURCE_BLOB
    assert payload["agent2_binding"]["schema"] == "kokuno-a2-strict-inner-spatial-candidate-v1"
    assert payload["agent2_binding"]["class"] == "KokunoStrictInnerLeadingOscillatorySpatialCandidate"

    a4 = payload["agent4_binding"]
    assert a4["dedicated_spatial_artifact_audit_present"] is False
    assert a4["authority"] is None
    assert a4["audited_agent2_pr"] is None
    assert a4["exact_head_ci_conclusion"] is None
    assert a4["independent_audit_conclusion"] is None

    evidence = payload["evidence"]
    assert evidence["agent2_exact_head_ci_conclusion"] is None
    assert evidence["agent2_construction_side_fd4_verifier_present"] is True
    assert evidence["agent4_dedicated_spatial_artifact_audit_present"] is False
    assert evidence["agent4_scientific_receipt_admitted"] is False
    assert evidence["strict_inner_spatial_artifact_scientifically_admitted"] is False

    handoff = payload["spatial_artifact_handoff"]
    assert handoff["registered"] is True
    assert handoff["status"] == "registered_waiting_for_independent_a4_audit"
    assert handoff["construction_authority"] == "Agent2#874"
    assert handoff["independent_audit_available"] is False
    assert handoff["independent_audit_authority"] is None
    assert handoff["public_velocity_binding"].endswith(".velocity")
    assert handoff["public_velocity_dt_binding"].endswith(".velocity_dt")
    assert handoff["public_velocity_jacobian_binding"].endswith(".velocity_jacobian")
    assert handoff["base_velocity_candidate_sha256_preserved"] is True
    assert handoff["differentiable_sha256_preserved"] is True
    assert handoff["spatial_sha256_registered"] is True
    assert handoff["manifest_sha256_registered"] is True
    assert handoff["pressure_in_artifact_surface"] is False
    assert handoff["restricted_forcing_in_artifact_surface"] is False
    assert handoff["agent3_correction_velocity_in_artifact"] is False
    assert handoff["outer_global_join_in_artifact"] is False
    assert handoff["complete_candidate_artifact"] is False
    assert handoff["eligible_for_velocity_export_ready"] is False
    assert handoff["usable_for_final_independent_pde_validation"] is False

    assert payload["ingest_status"]["strict_inner_spatial_artifact_ingest_admitted"] is False
    assert payload["stage_state"] == {
        "leading_ready": False,
        "oscillatory_ready": True,
        "correction_ready": False,
        "velocity_export_ready": False,
        "pde_validated": False,
    }


def test_agent2_spatial_engineering_checks_are_not_laundered() -> None:
    payload = _payload()
    protocol = payload["agent2_binding"]["protocol"]
    assert protocol["oscillatory_jacobian_operator"] == "centered_cartesian_fd6"
    assert protocol["oscillatory_jacobian_fixed_step"] == 1.0e-3
    assert protocol["caller_spatial_step_knob_present"] is False
    assert protocol["construction_side_fd4_sample_count"] == 12
    assert protocol["construction_side_fd4_steps"] == [4e-3, 2e-3, 1e-3]
    assert protocol["construction_side_fd4_refinement_ratio_gate"] == 6.0
    assert protocol["construction_side_fd4_fine_relative_rms_gate"] == 5e-3
    assert protocol["construction_side_fd4_fine_relative_sampled_max_gate"] == 1e-2

    firewall = payload["agent2_binding"]["norm_scope_firewall"]
    assert firewall["agent2_fd4_is_independent_agent4_validation"] is False
    assert firewall["jacobian_consistency_is_complete_ns_residual"] is False
    assert firewall["jacobian_consistency_is_cr001_momentum_norm"] is False
    assert firewall["derived_divergence_diagnostic_is_final_divergence_gate"] is False
    assert firewall["whole_domain_volume_weighted_l2_assessed"] is False
    assert firewall["same_protocol_st006_comparison_legal"] is False
    assert firewall["formal_full_domain_pde_gate_assessed"] is False

    api = payload["candidate_api_handoff"]
    assert api["strict_inner_artifact_velocity"] == "Agent2#874.velocity"
    assert api["strict_inner_artifact_velocity_dt"] == "Agent2#874.velocity_dt"
    assert api["strict_inner_artifact_velocity_jacobian"] == "Agent2#874.velocity_jacobian"
    assert api["pressure"] is None
    assert api["forcing"] is None
    assert api["complete_candidate_api_ready"] is False

    parent = deterministic_strict_inner_differentiable_artifact_ingest_contract(
        exact_head=PARENT_HEAD
    )
    assert payload["final_project_gates_unchanged"] == parent["final_project_gates_unchanged"]
    assert payload["truth_boundary"]["agent2_fd4_laundered_as_independent_a4_validation"] is False
    assert payload["truth_boundary"]["jacobian_consistency_laundered_as_pde_residual"] is False
    assert payload["truth_boundary"]["free_residual_defined_forcing_allowed"] is False
    assert payload["truth_boundary"]["threshold_relaxed"] is False


@pytest.mark.parametrize(
    "mutate",
    [
        lambda p: p["evidence"].__setitem__("agent2_exact_head_ci_conclusion", "success"),
        lambda p: p["evidence"].__setitem__("agent4_dedicated_spatial_artifact_audit_present", True),
        lambda p: p["agent4_binding"].__setitem__("authority", "Agent4#fake"),
        lambda p: p["agent2_binding"].__setitem__("head", "1" * 40),
        lambda p: p["agent2_binding"].__setitem__("source_blob_sha", "2" * 40),
        lambda p: p["agent2_binding"]["protocol"].__setitem__("oscillatory_jacobian_fixed_step", 2e-3),
        lambda p: p["agent2_binding"]["norm_scope_firewall"].__setitem__("agent2_fd4_is_independent_agent4_validation", True),
        lambda p: p["spatial_artifact_handoff"].__setitem__("independent_audit_available", True),
        lambda p: p["spatial_artifact_handoff"].__setitem__("pressure_in_artifact_surface", True),
        lambda p: p["spatial_artifact_handoff"].__setitem__("restricted_forcing_in_artifact_surface", True),
        lambda p: p["spatial_artifact_handoff"].__setitem__("agent3_correction_velocity_in_artifact", True),
        lambda p: p["spatial_artifact_handoff"].__setitem__("outer_global_join_in_artifact", True),
        lambda p: p["spatial_artifact_handoff"].__setitem__("complete_candidate_artifact", True),
        lambda p: p["spatial_artifact_handoff"].__setitem__("eligible_for_velocity_export_ready", True),
        lambda p: p["stage_state"].__setitem__("leading_ready", True),
        lambda p: p["stage_state"].__setitem__("correction_ready", True),
        lambda p: p["stage_state"].__setitem__("velocity_export_ready", True),
        lambda p: p["stage_state"].__setitem__("pde_validated", True),
        lambda p: p["truth_boundary"].__setitem__("free_residual_defined_forcing_allowed", True),
        lambda p: p["truth_boundary"].__setitem__("threshold_relaxed", True),
    ],
)
def test_truth_boundary_mutations_fail_closed(mutate) -> None:
    payload = copy.deepcopy(_payload())
    mutate(payload)
    with pytest.raises(ValueError, match="drifted"):
        validate_strict_inner_spatial_artifact_ingest_contract(payload)


def test_invalid_exact_head_fails_closed() -> None:
    payload = deterministic_strict_inner_spatial_artifact_ingest_contract(
        exact_head=None
    )
    payload["exact_head"] = "not-a-commit"
    with pytest.raises(ValueError, match="40-character"):
        validate_strict_inner_spatial_artifact_ingest_contract(payload)
