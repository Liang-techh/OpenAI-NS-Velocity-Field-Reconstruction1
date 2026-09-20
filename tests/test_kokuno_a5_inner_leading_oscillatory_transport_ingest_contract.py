from __future__ import annotations

from copy import deepcopy
import json

import pytest

from openai_ns_reconstruction.kokuno_a5_inner_leading_oscillatory_transport_ingest_contract import (
    AGENT1_LAPLACIAN_HEAD,
    AGENT2_TRANSPORT_HEAD,
    PARENT_A5_HEAD,
    SCHEMA,
    TRANSPORT_PROTOCOL,
    deterministic_inner_leading_oscillatory_transport_ingest_contract,
    main,
    validate_inner_leading_oscillatory_transport_ingest_contract,
)


def _payload() -> dict:
    return deterministic_inner_leading_oscillatory_transport_ingest_contract(
        exact_head="f" * 40
    )


def _reject(mutator) -> None:
    payload = _payload()
    mutator(payload)
    with pytest.raises(ValueError):
        validate_inner_leading_oscillatory_transport_ingest_contract(payload)


def test_contract_registers_transport_without_promoting_residual_or_readiness() -> None:
    payload = _payload()
    validate_inner_leading_oscillatory_transport_ingest_contract(payload)

    assert payload["schema"] == SCHEMA
    assert payload["parent_a5"] == {"pr": 833, "head": PARENT_A5_HEAD}
    assert payload["agent1_laplacian_binding"]["head"] == AGENT1_LAPLACIAN_HEAD
    assert payload["agent2_transport_binding"]["head"] == AGENT2_TRANSPORT_HEAD
    assert payload["agent2_transport_binding"]["protocol"] == TRANSPORT_PROTOCOL

    evidence = payload["evidence"]
    assert evidence == {
        "parent_inner_spatial_derivative_ingest_admitted": False,
        "agent1_laplacian_exact_head_ci_conclusion": None,
        "agent2_transport_exact_head_ci_conclusion": None,
        "agent4_dedicated_transport_audit_present": False,
        "agent4_transport_audit_conclusion": None,
        "strict_inner_transport_scientifically_admitted": False,
    }

    handoff = payload["transport_operator_handoff"]
    assert handoff["registered"] is True
    assert handoff["status"] == "registered_unresolved"
    assert handoff["authority"] == "Agent2#840"
    assert handoff["includes_velocity_dt"] is True
    assert handoff["includes_full_strict_inner_advection"] is True
    assert handoff["includes_base_viscosity"] is True
    assert handoff["includes_pressure_gradient"] is False
    assert handoff["includes_restricted_forcing"] is False
    assert handoff["includes_correction_velocity"] is False
    assert handoff["complete_ns_momentum_residual"] is False
    assert handoff["usable_as_actual_defect_input_now"] is False
    assert handoff["usable_for_final_independent_pde_validation"] is False

    assert payload["candidate_api_handoff"]["pressure"] is None
    assert payload["candidate_api_handoff"]["forcing"] is None
    assert payload["candidate_api_handoff"]["complete_candidate_api_ready"] is False

    status = payload["ingest_status"]
    assert status["typed_strict_inner_leading_oscillatory_transport_registered"] is True
    assert status["strict_inner_leading_oscillatory_transport_ingest_admitted"] is False
    assert status["strict_inner_transport_ready"] is False
    assert status["leading_ready"] is False

    assert payload["stage_state"] == {
        "leading_ready": False,
        "oscillatory_ready": True,
        "correction_ready": False,
        "velocity_export_ready": False,
        "pde_validated": False,
    }
    assert payload["baseline_vs_kokuno"][
        "kokuno_same_protocol_full_candidate_residual_available"
    ] is False
    assert payload["baseline_vs_kokuno"]["comparison_performed"] is False
    assert payload["final_project_gates_unchanged"] == {
        "normalized_momentum_max": 1e-3,
        "normalized_momentum_L2": 1e-3,
        "divergence_max": 1e-5,
        "divergence_L2": 1e-5,
    }


def test_rejects_queued_evidence_laundering() -> None:
    _reject(lambda p: p["evidence"].__setitem__(
        "agent2_transport_exact_head_ci_conclusion", "success"
    ))
    _reject(lambda p: p["evidence"].__setitem__(
        "agent4_dedicated_transport_audit_present", True
    ))
    _reject(lambda p: p["evidence"].__setitem__(
        "strict_inner_transport_scientifically_admitted", True
    ))


def test_rejects_upstream_identity_or_protocol_drift() -> None:
    _reject(lambda p: p["agent2_transport_binding"].__setitem__("head", "0" * 40))
    _reject(lambda p: p["agent2_transport_binding"].__setitem__(
        "source_blob_sha", "1" * 40
    ))
    _reject(lambda p: p["agent1_laplacian_binding"].__setitem__("head", "2" * 40))
    _reject(lambda p: p["agent2_transport_binding"]["protocol"].__setitem__(
        "viscosity", 0.02
    ))
    _reject(lambda p: p["agent2_transport_binding"]["protocol"].__setitem__(
        "independent_fd4_steps", [0.008, 0.004, 0.002]
    ))


def test_rejects_transport_being_laundered_into_complete_ns_defect() -> None:
    _reject(lambda p: p["transport_operator_handoff"].__setitem__(
        "includes_pressure_gradient", True
    ))
    _reject(lambda p: p["transport_operator_handoff"].__setitem__(
        "includes_restricted_forcing", True
    ))
    _reject(lambda p: p["transport_operator_handoff"].__setitem__(
        "complete_ns_momentum_residual", True
    ))
    _reject(lambda p: p["transport_operator_handoff"].__setitem__(
        "usable_as_actual_defect_input_now", True
    ))
    _reject(lambda p: p["staged_residual_availability"].__setitem__(
        "leading_plus_oscillatory_full_ns_residual_available", True
    ))


def test_rejects_candidate_api_or_stage_promotion() -> None:
    _reject(lambda p: p["candidate_api_handoff"].__setitem__("pressure", "invented"))
    _reject(lambda p: p["candidate_api_handoff"].__setitem__("forcing", "invented"))
    _reject(lambda p: p["ingest_status"].__setitem__(
        "strict_inner_leading_oscillatory_transport_ingest_admitted", True
    ))
    _reject(lambda p: p["stage_state"].__setitem__("leading_ready", True))
    _reject(lambda p: p["stage_state"].__setitem__("correction_ready", True))
    _reject(lambda p: p["stage_state"].__setitem__("velocity_export_ready", True))
    _reject(lambda p: p["stage_state"].__setitem__("pde_validated", True))


def test_rejects_ST006_comparison_threshold_or_forcing_escape_hatch() -> None:
    _reject(lambda p: p["baseline_vs_kokuno"].__setitem__("comparison_performed", True))
    _reject(lambda p: p["final_project_gates_unchanged"].__setitem__(
        "normalized_momentum_max", 1e-2
    ))
    _reject(lambda p: p["truth_boundary"].__setitem__(
        "free_residual_defined_forcing_allowed", True
    ))
    _reject(lambda p: p["truth_boundary"].__setitem__(
        "transport_relabelled_as_complete_ns_residual", True
    ))


def test_contract_digest_detects_metadata_mutation() -> None:
    payload = _payload()
    payload["contract_sha256"] = "0" * 64
    with pytest.raises(ValueError):
        validate_inner_leading_oscillatory_transport_ingest_contract(payload)


def test_cli_writes_valid_exact_head_report(tmp_path) -> None:
    output = tmp_path / "report.json"
    exact_head = "a" * 40
    assert main(["--exact-head", exact_head, "--output", str(output)]) == 0
    payload = json.loads(output.read_text())
    assert payload["exact_head"] == exact_head
    validate_inner_leading_oscillatory_transport_ingest_contract(payload)


def test_copy_mutation_does_not_modify_fresh_contract() -> None:
    first = _payload()
    second = deepcopy(first)
    second["transport_operator_handoff"]["includes_pressure_gradient"] = True
    assert first["transport_operator_handoff"]["includes_pressure_gradient"] is False
