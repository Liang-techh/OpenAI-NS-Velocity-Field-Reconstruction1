from __future__ import annotations

from copy import deepcopy
import json

import pytest

from openai_ns_reconstruction.kokuno_a5_oscillatory_transport_delta_ingest_contract import (
    AGENT2_DELTA_HEAD,
    AGENT2_DELTA_PROTOCOL,
    AGENT4_DELTA_AUDIT_HEAD,
    AGENT4_DELTA_PROTOCOL,
    NORM_SCOPE_FIREWALL,
    PARENT_A5_HEAD,
    SCHEMA,
    deterministic_oscillatory_transport_delta_ingest_contract,
    main,
    validate_oscillatory_transport_delta_ingest_contract,
)


def _payload() -> dict:
    return deterministic_oscillatory_transport_delta_ingest_contract(
        exact_head="f" * 40
    )


def _reject(mutator) -> None:
    payload = _payload()
    mutator(payload)
    with pytest.raises(ValueError):
        validate_oscillatory_transport_delta_ingest_contract(payload)


def test_contract_registers_before_after_delta_without_pde_promotion() -> None:
    payload = _payload()
    validate_oscillatory_transport_delta_ingest_contract(payload)

    assert payload["schema"] == SCHEMA
    assert payload["parent_a5"] == {"pr": 843, "head": PARENT_A5_HEAD}
    assert payload["agent2_delta_binding"]["head"] == AGENT2_DELTA_HEAD
    assert payload["agent2_delta_binding"]["protocol"] == AGENT2_DELTA_PROTOCOL
    assert payload["agent4_delta_audit_binding"]["head"] == AGENT4_DELTA_AUDIT_HEAD
    assert payload["agent4_delta_audit_binding"]["protocol"] == AGENT4_DELTA_PROTOCOL
    assert payload["agent4_delta_audit_binding"]["norm_scope_firewall"] == NORM_SCOPE_FIREWALL

    assert payload["evidence"] == {
        "parent_strict_inner_transport_ingest_admitted": False,
        "agent2_delta_exact_head_ci_conclusion": None,
        "agent4_dedicated_delta_audit_present": True,
        "agent4_delta_exact_head_ci_conclusion": None,
        "agent4_delta_audit_conclusion": None,
        "agent4_delta_scientific_receipt_admitted": False,
        "oscillatory_transport_delta_scientifically_admitted": False,
        "staged_transport_ratio_numeric_value": None,
    }

    handoff = payload["oscillatory_transport_delta_handoff"]
    assert handoff["registered"] is True
    assert handoff["status"] == "registered_unresolved"
    assert handoff["construction_authority"] == "Agent2#849"
    assert handoff["independent_audit_available"] is True
    assert handoff["independent_audit_authority"] == "Agent4#851"
    assert handoff["quantity"] == "T(u_inner+u_osc)-T(u_inner)"
    assert handoff["ratio_R_definition"] == "RMS(T_inner+osc)/RMS(T_inner)"
    assert handoff["ratio_R_observation_only"] is True
    assert handoff["ratio_R_available_now"] is False
    assert handoff["complete_ns_momentum_residual"] is False
    assert handoff["usable_as_full_candidate_defect_now"] is False
    assert handoff["usable_for_final_independent_pde_validation"] is False

    staged = payload["staged_comparison"]
    assert staged["leading_only_transport_surface_available"] is True
    assert staged["leading_plus_oscillatory_transport_surface_available"] is True
    assert staged["oscillatory_transport_delta_surface_available"] is True
    assert staged["independent_before_after_delta_audit_registered"] is True
    assert staged["numeric_before_after_ratio_available"] is False
    assert staged["leading_only_full_ns_residual_available"] is False
    assert staged["leading_plus_oscillatory_full_ns_residual_available"] is False
    assert staged["after_correction_full_ns_residual_available"] is False
    assert staged["same_protocol_ST006_comparison_legal"] is False

    assert payload["candidate_api_handoff"]["pressure"] is None
    assert payload["candidate_api_handoff"]["forcing"] is None
    assert payload["candidate_api_handoff"]["complete_candidate_api_ready"] is False

    status = payload["ingest_status"]
    assert status["typed_oscillatory_transport_delta_registered"] is True
    assert status["typed_oscillatory_transport_delta_independent_audit_registered"] is True
    assert status["oscillatory_transport_delta_ingest_admitted"] is False
    assert status["staged_transport_attribution_ready"] is False

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


def test_rejects_queued_evidence_or_ratio_laundering() -> None:
    _reject(lambda p: p["evidence"].__setitem__(
        "agent2_delta_exact_head_ci_conclusion", "success"
    ))
    _reject(lambda p: p["evidence"].__setitem__(
        "agent4_delta_exact_head_ci_conclusion", "success"
    ))
    _reject(lambda p: p["evidence"].__setitem__("agent4_delta_audit_conclusion", "pass"))
    _reject(lambda p: p["evidence"].__setitem__(
        "agent4_delta_scientific_receipt_admitted", True
    ))
    _reject(lambda p: p["evidence"].__setitem__("staged_transport_ratio_numeric_value", 0.9))
    _reject(lambda p: p["staged_comparison"].__setitem__(
        "numeric_before_after_ratio_available", True
    ))


def test_rejects_upstream_identity_or_protocol_drift() -> None:
    _reject(lambda p: p["agent2_delta_binding"].__setitem__("head", "0" * 40))
    _reject(lambda p: p["agent2_delta_binding"].__setitem__("source_blob_sha", "1" * 40))
    _reject(lambda p: p["agent4_delta_audit_binding"].__setitem__("head", "2" * 40))
    _reject(lambda p: p["agent4_delta_audit_binding"].__setitem__("source_blob_sha", "3" * 40))
    _reject(lambda p: p["agent2_delta_binding"]["protocol"].__setitem__("viscosity", 0.02))
    _reject(lambda p: p["agent4_delta_audit_binding"]["protocol"].__setitem__(
        "fd8_steps", [0.004, 0.002, 0.001]
    ))


def test_rejects_norm_scope_or_residual_laundering() -> None:
    _reject(lambda p: p["agent4_delta_audit_binding"]["norm_scope_firewall"].__setitem__(
        "transport_delta_is_complete_ns_momentum_residual", True
    ))
    _reject(lambda p: p["agent4_delta_audit_binding"]["norm_scope_firewall"].__setitem__(
        "transport_ratio_is_cr001_acceptance_metric", True
    ))
    _reject(lambda p: p["agent4_delta_audit_binding"]["norm_scope_firewall"].__setitem__(
        "agent4_divergence_rms_is_canonical_volume_weighted_L2", True
    ))
    _reject(lambda p: p["oscillatory_transport_delta_handoff"].__setitem__(
        "complete_ns_momentum_residual", True
    ))
    _reject(lambda p: p["oscillatory_transport_delta_handoff"].__setitem__(
        "usable_as_full_candidate_defect_now", True
    ))
    _reject(lambda p: p["staged_comparison"].__setitem__(
        "leading_plus_oscillatory_full_ns_residual_available", True
    ))
    _reject(lambda p: p["staged_comparison"].__setitem__(
        "same_protocol_ST006_comparison_legal", True
    ))


def test_rejects_missing_terms_or_stage_promotion() -> None:
    for key in (
        "pressure_included",
        "restricted_forcing_included",
        "correction_velocity_included",
        "outer_global_join_included",
    ):
        _reject(lambda p, key=key: p["oscillatory_transport_delta_handoff"].__setitem__(key, True))
    _reject(lambda p: p["ingest_status"].__setitem__(
        "oscillatory_transport_delta_ingest_admitted", True
    ))
    _reject(lambda p: p["ingest_status"].__setitem__("staged_transport_attribution_ready", True))
    _reject(lambda p: p["stage_state"].__setitem__("leading_ready", True))
    _reject(lambda p: p["stage_state"].__setitem__("correction_ready", True))
    _reject(lambda p: p["stage_state"].__setitem__("velocity_export_ready", True))
    _reject(lambda p: p["stage_state"].__setitem__("pde_validated", True))


def test_rejects_threshold_forcing_or_truth_escape_hatch() -> None:
    _reject(lambda p: p["baseline_vs_kokuno"].__setitem__("comparison_performed", True))
    _reject(lambda p: p["final_project_gates_unchanged"].__setitem__(
        "normalized_momentum_max", 1e-2
    ))
    _reject(lambda p: p["truth_boundary"].__setitem__(
        "free_residual_defined_forcing_allowed", True
    ))
    _reject(lambda p: p["truth_boundary"].__setitem__(
        "transport_ratio_laundered_as_cr001_acceptance", True
    ))
    _reject(lambda p: p["truth_boundary"].__setitem__("full_ns_residual_invented", True))


def test_contract_digest_detects_mutation() -> None:
    payload = _payload()
    payload["contract_sha256"] = "0" * 64
    with pytest.raises(ValueError):
        validate_oscillatory_transport_delta_ingest_contract(payload)


def test_cli_writes_valid_exact_head_report(tmp_path) -> None:
    output = tmp_path / "report.json"
    exact_head = "a" * 40
    assert main(["--exact-head", exact_head, "--output", str(output)]) == 0
    payload = json.loads(output.read_text())
    assert payload["exact_head"] == exact_head
    validate_oscillatory_transport_delta_ingest_contract(payload)


def test_copy_mutation_does_not_modify_fresh_contract() -> None:
    first = _payload()
    second = deepcopy(first)
    second["oscillatory_transport_delta_handoff"]["pressure_included"] = True
    assert first["oscillatory_transport_delta_handoff"]["pressure_included"] is False
