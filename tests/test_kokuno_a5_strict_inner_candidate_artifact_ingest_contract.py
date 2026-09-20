from __future__ import annotations

from copy import deepcopy
import json

import pytest

from openai_ns_reconstruction.kokuno_a5_strict_inner_candidate_artifact_ingest_contract import (
    AGENT2_ARTIFACT_PROTOCOL,
    AGENT2_CANDIDATE_HEAD,
    AGENT4_ARTIFACT_AUDIT_HEAD,
    AGENT4_ARTIFACT_AUDIT_PROTOCOL,
    NORM_SCOPE_FIREWALL,
    PARENT_A5_HEAD,
    SCHEMA,
    deterministic_strict_inner_candidate_artifact_ingest_contract,
    main,
    validate_strict_inner_candidate_artifact_ingest_contract,
)


def _payload() -> dict:
    return deterministic_strict_inner_candidate_artifact_ingest_contract(
        exact_head="f" * 40
    )


def _reject(mutator) -> None:
    payload = _payload()
    mutator(payload)
    with pytest.raises(ValueError):
        validate_strict_inner_candidate_artifact_ingest_contract(payload)


def test_registers_first_scoped_candidate_artifact_without_full_candidate_promotion() -> None:
    payload = _payload()
    validate_strict_inner_candidate_artifact_ingest_contract(payload)

    assert payload["schema"] == SCHEMA
    assert payload["parent_a5"] == {"pr": 852, "head": PARENT_A5_HEAD}
    assert payload["pipeline_position"]["registered_scoped_stage"] == "strict_inner_candidate_artifact"
    assert payload["pipeline_position"]["full_pipeline_candidate_artifact_stage_reached"] is False

    a2 = payload["agent2_candidate_binding"]
    assert a2["pr"] == 857
    assert a2["head"] == AGENT2_CANDIDATE_HEAD
    assert a2["source_blob_sha"] == "82e732af14695c39afa3b631c00b7f20960976da"
    assert a2["schema"] == "kokuno-a2-strict-inner-leading-oscillatory-candidate-v1"
    assert a2["class"] == "KokunoStrictInnerLeadingOscillatoryCandidate"
    assert a2["protocol"] == AGENT2_ARTIFACT_PROTOCOL

    a4 = payload["agent4_artifact_audit_binding"]
    assert a4["pr"] == 859
    assert a4["head"] == AGENT4_ARTIFACT_AUDIT_HEAD
    assert a4["source_blob_sha"] == "d65d4931f2734de32e101f21d48c7b8d2381cbde"
    assert a4["schema"] == "kokuno-a4-strict-inner-candidate-artifact-independent-audit-v1"
    assert a4["protocol"] == AGENT4_ARTIFACT_AUDIT_PROTOCOL
    assert a4["norm_scope_firewall"] == NORM_SCOPE_FIREWALL

    assert payload["evidence"] == {
        "parent_oscillatory_transport_delta_ingest_admitted": False,
        "agent2_candidate_exact_head_ci_conclusion": None,
        "agent4_dedicated_artifact_audit_present": True,
        "agent4_artifact_audit_exact_head_ci_conclusion": None,
        "agent4_artifact_audit_conclusion": None,
        "agent4_artifact_scientific_receipt_admitted": False,
        "strict_inner_candidate_artifact_scientifically_admitted": False,
    }

    handoff = payload["strict_inner_candidate_artifact_handoff"]
    assert handoff["registered"] is True
    assert handoff["status"] == "registered_unresolved"
    assert handoff["construction_authority"] == "Agent2#857"
    assert handoff["independent_audit_authority"] == "Agent4#859"
    assert handoff["semantic_candidate_sha256_registered"] is True
    assert handoff["manifest_sha256_registered"] is True
    assert handoff["python_callable_surface_registered"] is True
    assert handoff["velocity_dt_in_artifact_surface"] is False
    assert handoff["pressure_in_artifact_surface"] is False
    assert handoff["restricted_forcing_in_artifact_surface"] is False
    assert handoff["agent3_correction_velocity_in_artifact"] is False
    assert handoff["outer_global_join_in_artifact"] is False
    assert handoff["complete_candidate_artifact"] is False
    assert handoff["eligible_for_velocity_export_ready"] is False
    assert handoff["usable_as_complete_ns_defect_input"] is False
    assert handoff["usable_for_final_independent_pde_validation"] is False

    artifact = payload["artifact_status"]
    assert artifact["upstream_strict_inner_artifact_implementation_present"] is True
    assert artifact["upstream_manifest_save_load_surface_present"] is True
    assert artifact["upstream_independent_artifact_audit_present"] is True
    assert artifact["a5_artifact_schema_binding_registered"] is True
    assert artifact["a5_materialized_candidate_payload"] is False
    assert artifact["strict_inner_candidate_artifact_ingest_admitted"] is False
    assert artifact["complete_global_candidate_artifact_materialized"] is False
    assert artifact["candidate_parameter_digest_for_full_pipeline"] is None
    assert artifact["candidate_validation_receipt_for_full_pipeline"] is None

    api = payload["candidate_api_handoff"]
    assert api["pressure"] is None
    assert api["forcing"] is None
    assert api["complete_candidate_api_ready"] is False

    status = payload["ingest_status"]
    assert status["typed_strict_inner_candidate_artifact_registered"] is True
    assert status["typed_strict_inner_candidate_artifact_independent_audit_registered"] is True
    assert status["strict_inner_candidate_artifact_ingest_admitted"] is False
    assert status["full_candidate_artifact_ready"] is False

    assert payload["stage_state"] == {
        "leading_ready": False,
        "oscillatory_ready": True,
        "correction_ready": False,
        "velocity_export_ready": False,
        "pde_validated": False,
    }
    assert payload["baseline_vs_kokuno"]["comparison_performed"] is False
    assert payload["final_project_gates_unchanged"] == {
        "normalized_momentum_max": 1e-3,
        "normalized_momentum_L2": 1e-3,
        "divergence_max": 1e-5,
        "divergence_L2": 1e-5,
    }


def test_rejects_queued_evidence_laundering() -> None:
    _reject(lambda p: p["evidence"].__setitem__(
        "agent2_candidate_exact_head_ci_conclusion", "success"
    ))
    _reject(lambda p: p["evidence"].__setitem__(
        "agent4_artifact_audit_exact_head_ci_conclusion", "success"
    ))
    _reject(lambda p: p["evidence"].__setitem__("agent4_artifact_audit_conclusion", "pass"))
    _reject(lambda p: p["evidence"].__setitem__(
        "agent4_artifact_scientific_receipt_admitted", True
    ))
    _reject(lambda p: p["artifact_status"].__setitem__(
        "strict_inner_candidate_artifact_ingest_admitted", True
    ))


def test_rejects_candidate_or_audit_identity_protocol_drift() -> None:
    _reject(lambda p: p["agent2_candidate_binding"].__setitem__("head", "0" * 40))
    _reject(lambda p: p["agent2_candidate_binding"].__setitem__("source_blob_sha", "1" * 40))
    _reject(lambda p: p["agent4_artifact_audit_binding"].__setitem__("head", "2" * 40))
    _reject(lambda p: p["agent4_artifact_audit_binding"].__setitem__("source_blob_sha", "3" * 40))
    _reject(lambda p: p["agent2_candidate_binding"]["protocol"].__setitem__(
        "batch_scalar_consistency_gate", 1e-8
    ))
    _reject(lambda p: p["agent4_artifact_audit_binding"]["protocol"].__setitem__(
        "nominal_richardson_steps", [0.002, 0.001, 0.0005]
    ))


def test_rejects_norm_scope_or_pde_laundering() -> None:
    _reject(lambda p: p["agent4_artifact_audit_binding"]["norm_scope_firewall"].__setitem__(
        "agent4_divergence_rms_is_whole_domain_volume_weighted_L2", True
    ))
    _reject(lambda p: p["agent4_artifact_audit_binding"]["norm_scope_firewall"].__setitem__(
        "formal_cr001_pde_gate_assessed", True
    ))
    _reject(lambda p: p["agent4_artifact_audit_binding"]["norm_scope_firewall"].__setitem__(
        "same_protocol_ST006_comparison_legal", True
    ))
    _reject(lambda p: p["strict_inner_candidate_artifact_handoff"].__setitem__(
        "usable_for_final_independent_pde_validation", True
    ))
    _reject(lambda p: p["truth_boundary"].__setitem__(
        "scoped_divergence_laundered_as_pde_validation", True
    ))


def test_rejects_full_candidate_or_export_promotion() -> None:
    for key in (
        "velocity_dt_in_artifact_surface",
        "pressure_in_artifact_surface",
        "restricted_forcing_in_artifact_surface",
        "agent3_correction_velocity_in_artifact",
        "outer_global_join_in_artifact",
        "complete_candidate_artifact",
        "eligible_for_velocity_export_ready",
        "usable_as_complete_ns_defect_input",
    ):
        _reject(lambda p, key=key: p["strict_inner_candidate_artifact_handoff"].__setitem__(key, True))
    _reject(lambda p: p["artifact_status"].__setitem__("a5_materialized_candidate_payload", True))
    _reject(lambda p: p["artifact_status"].__setitem__("complete_global_candidate_artifact_materialized", True))
    _reject(lambda p: p["ingest_status"].__setitem__("full_candidate_artifact_ready", True))
    _reject(lambda p: p["stage_state"].__setitem__("leading_ready", True))
    _reject(lambda p: p["stage_state"].__setitem__("correction_ready", True))
    _reject(lambda p: p["stage_state"].__setitem__("velocity_export_ready", True))
    _reject(lambda p: p["stage_state"].__setitem__("pde_validated", True))


def test_rejects_baseline_threshold_forcing_or_truth_escape_hatch() -> None:
    _reject(lambda p: p["staged_comparison"].__setitem__(
        "strict_inner_artifact_same_protocol_ST006_comparison_legal", True
    ))
    _reject(lambda p: p["baseline_vs_kokuno"].__setitem__("comparison_performed", True))
    _reject(lambda p: p["final_project_gates_unchanged"].__setitem__(
        "normalized_momentum_max", 1e-2
    ))
    _reject(lambda p: p["truth_boundary"].__setitem__(
        "free_residual_defined_forcing_allowed", True
    ))
    _reject(lambda p: p["truth_boundary"].__setitem__("full_ns_residual_invented", True))
    _reject(lambda p: p["truth_boundary"].__setitem__(
        "premature_ST006_comparison_performed", True
    ))


def test_contract_digest_detects_mutation() -> None:
    payload = _payload()
    payload["contract_sha256"] = "0" * 64
    with pytest.raises(ValueError):
        validate_strict_inner_candidate_artifact_ingest_contract(payload)


def test_cli_writes_valid_exact_head_report(tmp_path) -> None:
    output = tmp_path / "report.json"
    exact_head = "a" * 40
    assert main(["--exact-head", exact_head, "--output", str(output)]) == 0
    payload = json.loads(output.read_text())
    assert payload["exact_head"] == exact_head
    validate_strict_inner_candidate_artifact_ingest_contract(payload)


def test_copy_mutation_does_not_modify_fresh_contract() -> None:
    first = _payload()
    second = deepcopy(first)
    second["strict_inner_candidate_artifact_handoff"]["matched_pressure_invented"] = True
    assert "matched_pressure_invented" not in first["strict_inner_candidate_artifact_handoff"]
