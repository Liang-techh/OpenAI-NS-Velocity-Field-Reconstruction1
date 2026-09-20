from __future__ import annotations

import copy

import pytest

from openai_ns_reconstruction.kokuno_a5_cartesian_center_ingest_contract import (
    AGENT1_CARTESIAN_CENTER_HEAD,
    AGENT1_CARTESIAN_CENTER_SOURCE_BLOB,
    AGENT4_CARTESIAN_CENTER_HEAD,
    AGENT4_CARTESIAN_CENTER_SOURCE_BLOB,
    SCHEMA,
    deterministic_cartesian_center_ingest_contract,
    validate_cartesian_center_ingest_contract,
)


def _payload() -> dict:
    return deterministic_cartesian_center_ingest_contract(exact_head="deadbeef")


def test_contract_is_fail_closed_and_valid() -> None:
    payload = _payload()
    validate_cartesian_center_ingest_contract(payload)

    assert payload["schema"] == SCHEMA
    assert payload["agent1_cartesian_center_binding"]["head"] == AGENT1_CARTESIAN_CENTER_HEAD
    assert (
        payload["agent1_cartesian_center_binding"]["source_blob_sha"]
        == AGENT1_CARTESIAN_CENTER_SOURCE_BLOB
    )
    assert payload["agent4_cartesian_center_binding"]["head"] == AGENT4_CARTESIAN_CENTER_HEAD
    assert (
        payload["agent4_cartesian_center_binding"]["source_blob_sha"]
        == AGENT4_CARTESIAN_CENTER_SOURCE_BLOB
    )

    evidence = payload["evidence"]
    assert evidence["agent1_exact_head_ci_conclusion"] is None
    assert evidence["agent4_exact_head_ci_conclusion"] is None
    assert evidence["agent4_inner_cartesian_center_independently_audited"] is False
    assert evidence["source_C_independent_admission_inherited"] is False

    api = payload["candidate_api_handoff"]
    assert api["velocity"] == "Agent1#803.KokunoPA10CartesianCenterVelocity.velocity"
    assert api["velocity_dt"] is None
    assert api["pressure"] is None
    assert api["forcing"] is None
    assert api["complete_candidate_api_ready"] is False

    status = payload["ingest_status"]
    assert status["typed_inner_cartesian_center_velocity_interface_registered"] is True
    assert status["typed_inner_cartesian_center_independent_audit_registered"] is True
    assert status["inner_cartesian_center_velocity_ingest_admitted"] is False
    assert status["inner_cartesian_center_field_identity_registered"] is True
    assert status["outer_global_leading_join_materialized"] is False
    assert status["global_leading_velocity_materialized"] is False
    assert status["matched_global_pressure_materialized"] is False
    assert status["restricted_forcing_composite_validated"] is False

    assert payload["stage_state"] == {
        "leading_ready": False,
        "oscillatory_ready": True,
        "correction_ready": False,
        "velocity_export_ready": False,
        "pde_validated": False,
    }
    assert payload["baseline_vs_kokuno"]["comparison_performed"] is False
    assert (
        payload["baseline_vs_kokuno"][
            "kokuno_same_protocol_full_candidate_residual_available"
        ]
        is False
    )


def test_rejects_queued_evidence_laundering() -> None:
    for field, value in (
        ("agent1_exact_head_ci_conclusion", "success"),
        ("agent4_exact_head_ci_conclusion", "success"),
        ("agent4_inner_cartesian_center_independently_audited", True),
        ("source_C_independent_admission_inherited", True),
    ):
        payload = _payload()
        payload["evidence"][field] = value
        with pytest.raises(ValueError):
            validate_cartesian_center_ingest_contract(payload)


def test_rejects_candidate_api_invention() -> None:
    for field, value in (
        ("velocity_dt", "invented.velocity_dt"),
        ("pressure", "invented.pressure"),
        ("forcing", "residual_defined_force"),
        ("complete_candidate_api_ready", True),
    ):
        payload = _payload()
        payload["candidate_api_handoff"][field] = value
        with pytest.raises(ValueError):
            validate_cartesian_center_ingest_contract(payload)


def test_rejects_scientific_state_promotion() -> None:
    for field in (
        "inner_cartesian_center_velocity_ingest_admitted",
        "outer_global_leading_join_materialized",
        "global_leading_velocity_materialized",
        "matched_global_pressure_materialized",
        "restricted_forcing_composite_validated",
        "leading_ready",
    ):
        payload = _payload()
        payload["ingest_status"][field] = True
        with pytest.raises(ValueError):
            validate_cartesian_center_ingest_contract(payload)

    for field in ("leading_ready", "correction_ready", "velocity_export_ready", "pde_validated"):
        payload = _payload()
        payload["stage_state"][field] = True
        with pytest.raises(ValueError):
            validate_cartesian_center_ingest_contract(payload)


def test_rejects_baseline_or_gate_laundering() -> None:
    payload = _payload()
    payload["baseline_vs_kokuno"]["comparison_performed"] = True
    with pytest.raises(ValueError):
        validate_cartesian_center_ingest_contract(payload)

    payload = _payload()
    payload["baseline_vs_kokuno"]["kokuno_same_protocol_full_candidate_residual_available"] = True
    with pytest.raises(ValueError):
        validate_cartesian_center_ingest_contract(payload)

    payload = _payload()
    payload["final_project_gates_unchanged"]["normalized_momentum_max"] = 1.1e-3
    with pytest.raises(ValueError):
        validate_cartesian_center_ingest_contract(payload)


def test_rejects_pin_truth_or_digest_drift() -> None:
    payload = _payload()
    payload["agent1_cartesian_center_binding"]["source_blob_sha"] = "0" * 40
    with pytest.raises(ValueError):
        validate_cartesian_center_ingest_contract(payload)

    payload = _payload()
    payload["agent4_cartesian_center_binding"]["source_blob_sha"] = "0" * 40
    with pytest.raises(ValueError):
        validate_cartesian_center_ingest_contract(payload)

    payload = _payload()
    payload["truth_boundary"]["velocity_export_ready"] = True
    with pytest.raises(ValueError):
        validate_cartesian_center_ingest_contract(payload)

    payload = _payload()
    payload["contract_sha256"] = "0" * 64
    with pytest.raises(ValueError):
        validate_cartesian_center_ingest_contract(payload)


def test_digest_detects_other_mutation() -> None:
    payload = _payload()
    mutated = copy.deepcopy(payload)
    mutated["agent1_cartesian_center_binding"]["registered_time_interval"] = [0.25, 0.74]
    with pytest.raises(ValueError):
        validate_cartesian_center_ingest_contract(mutated)
