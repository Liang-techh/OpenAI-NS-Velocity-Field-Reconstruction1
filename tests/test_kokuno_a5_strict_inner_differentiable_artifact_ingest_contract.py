from __future__ import annotations

import copy

import pytest

from openai_ns_reconstruction.kokuno_a5_strict_inner_candidate_artifact_ingest_contract import (
    deterministic_strict_inner_candidate_artifact_ingest_contract,
)
from openai_ns_reconstruction.kokuno_a5_strict_inner_differentiable_artifact_ingest_contract import (
    AGENT2_HEAD,
    AGENT2_SOURCE_BLOB,
    AGENT4_HEAD,
    AGENT4_SOURCE_BLOB,
    deterministic_strict_inner_differentiable_artifact_ingest_contract,
    validate_strict_inner_differentiable_artifact_ingest_contract,
)


EXACT_HEAD = "0" * 40
PARENT_HEAD = "d295670e51802abc345801bcdfef0c99bbf91e98"


def _payload():
    return deterministic_strict_inner_differentiable_artifact_ingest_contract(
        exact_head=EXACT_HEAD
    )


def test_registered_differentiable_artifact_remains_fail_closed() -> None:
    payload = _payload()
    validate_strict_inner_differentiable_artifact_ingest_contract(payload)

    assert payload["agent2_binding"]["head"] == AGENT2_HEAD
    assert payload["agent2_binding"]["source_blob_sha"] == AGENT2_SOURCE_BLOB
    assert payload["agent4_binding"]["head"] == AGENT4_HEAD
    assert payload["agent4_binding"]["source_blob_sha"] == AGENT4_SOURCE_BLOB

    handoff = payload["differentiable_artifact_handoff"]
    assert handoff["registered"] is True
    assert handoff["public_velocity_binding"].endswith(".velocity")
    assert handoff["public_velocity_dt_binding"].endswith(".velocity_dt")
    assert handoff["base_velocity_candidate_sha256_preserved"] is True
    assert handoff["differentiable_sha256_registered"] is True
    assert handoff["pressure_in_artifact_surface"] is False
    assert handoff["restricted_forcing_in_artifact_surface"] is False
    assert handoff["agent3_correction_velocity_in_artifact"] is False
    assert handoff["outer_global_join_in_artifact"] is False
    assert handoff["complete_candidate_artifact"] is False
    assert handoff["eligible_for_velocity_export_ready"] is False
    assert handoff["usable_for_final_independent_pde_validation"] is False

    evidence = payload["evidence"]
    assert evidence["agent2_exact_head_ci_conclusion"] is None
    assert evidence["agent4_exact_head_ci_conclusion"] is None
    assert evidence["agent4_independent_audit_conclusion"] is None
    assert evidence["agent4_scientific_receipt_admitted"] is False

    assert payload["ingest_status"]["strict_inner_differentiable_artifact_ingest_admitted"] is False
    assert payload["stage_state"]["leading_ready"] is False
    assert payload["stage_state"]["correction_ready"] is False
    assert payload["stage_state"]["velocity_export_ready"] is False
    assert payload["stage_state"]["pde_validated"] is False


def test_velocity_dt_audit_is_not_laundered_into_pde_validation() -> None:
    payload = _payload()
    firewall = payload["agent4_binding"]["norm_scope_firewall"]
    assert firewall["velocity_dt_consistency_is_complete_ns_residual"] is False
    assert firewall["velocity_dt_consistency_is_cr001_momentum_norm"] is False
    assert firewall["whole_domain_volume_weighted_l2_assessed"] is False
    assert firewall["same_protocol_st006_comparison_legal"] is False

    api = payload["candidate_api_handoff"]
    assert api["strict_inner_artifact_velocity"] == "Agent2#866.velocity"
    assert api["strict_inner_artifact_velocity_dt"] == "Agent2#866.velocity_dt"
    assert api["pressure"] is None
    assert api["forcing"] is None
    assert api["complete_candidate_api_ready"] is False

    parent = deterministic_strict_inner_candidate_artifact_ingest_contract(
        exact_head=PARENT_HEAD
    )
    assert (
        payload["final_project_gates_unchanged"]
        == parent["final_project_gates_unchanged"]
    )
    assert payload["truth_boundary"]["free_residual_defined_forcing_allowed"] is False
    assert payload["truth_boundary"]["threshold_relaxed"] is False


@pytest.mark.parametrize(
    "mutate",
    [
        lambda p: p["evidence"].__setitem__("agent2_exact_head_ci_conclusion", "success"),
        lambda p: p["evidence"].__setitem__("agent4_exact_head_ci_conclusion", "success"),
        lambda p: p["evidence"].__setitem__("agent4_independent_audit_conclusion", "pass"),
        lambda p: p["agent2_binding"].__setitem__("head", "1" * 40),
        lambda p: p["agent4_binding"].__setitem__("head", "2" * 40),
        lambda p: p["differentiable_artifact_handoff"].__setitem__("pressure_in_artifact_surface", True),
        lambda p: p["differentiable_artifact_handoff"].__setitem__("restricted_forcing_in_artifact_surface", True),
        lambda p: p["differentiable_artifact_handoff"].__setitem__("outer_global_join_in_artifact", True),
        lambda p: p["differentiable_artifact_handoff"].__setitem__("complete_candidate_artifact", True),
        lambda p: p["differentiable_artifact_handoff"].__setitem__("eligible_for_velocity_export_ready", True),
        lambda p: p["agent4_binding"]["norm_scope_firewall"].__setitem__("velocity_dt_consistency_is_complete_ns_residual", True),
        lambda p: p["agent4_binding"]["norm_scope_firewall"].__setitem__("whole_domain_volume_weighted_l2_assessed", True),
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
        validate_strict_inner_differentiable_artifact_ingest_contract(payload)


def test_invalid_exact_head_fails_closed() -> None:
    payload = deterministic_strict_inner_differentiable_artifact_ingest_contract(
        exact_head=None
    )
    payload["exact_head"] = "not-a-commit"
    with pytest.raises(ValueError, match="40-character"):
        validate_strict_inner_differentiable_artifact_ingest_contract(payload)
