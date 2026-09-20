from __future__ import annotations

import copy

import pytest

from openai_ns_reconstruction.kokuno_a5_cartesian_center_velocity_dt_ingest_contract import (
    AGENT1_VELOCITY_DT_HEAD,
    AGENT1_VELOCITY_DT_SOURCE_BLOB,
    AGENT4_VELOCITY_DT_HEAD,
    AGENT4_VELOCITY_DT_SOURCE_BLOB,
    SCHEMA,
    deterministic_cartesian_center_velocity_dt_ingest_contract,
    validate_cartesian_center_velocity_dt_ingest_contract,
)


def _payload() -> dict:
    return deterministic_cartesian_center_velocity_dt_ingest_contract(
        exact_head="deadbeef"
    )


def test_contract_registers_velocity_dt_but_remains_fail_closed() -> None:
    payload = _payload()
    validate_cartesian_center_velocity_dt_ingest_contract(payload)

    assert payload["schema"] == SCHEMA
    assert payload["agent1_velocity_dt_binding"]["head"] == AGENT1_VELOCITY_DT_HEAD
    assert (
        payload["agent1_velocity_dt_binding"]["source_blob_sha"]
        == AGENT1_VELOCITY_DT_SOURCE_BLOB
    )
    assert payload["agent4_velocity_dt_binding"]["head"] == AGENT4_VELOCITY_DT_HEAD
    assert (
        payload["agent4_velocity_dt_binding"]["source_blob_sha"]
        == AGENT4_VELOCITY_DT_SOURCE_BLOB
    )

    evidence = payload["evidence"]
    assert evidence["parent_inner_cartesian_center_velocity_ingest_admitted"] is False
    assert evidence["agent1_exact_head_ci_conclusion"] is None
    assert evidence["agent4_exact_head_ci_conclusion"] is None
    assert (
        evidence[
            "agent4_inner_cartesian_center_velocity_dt_independently_audited"
        ]
        is False
    )

    api = payload["candidate_api_handoff"]
    assert api["velocity"] == "Agent1#803.KokunoPA10CartesianCenterVelocity.velocity"
    assert (
        api["velocity_dt"]
        == "Agent1#811.KokunoPA10CartesianCenterVelocityTimeDerivative.velocity_dt"
    )
    assert api["pressure"] is None
    assert api["forcing"] is None
    assert api["complete_candidate_api_ready"] is False

    status = payload["ingest_status"]
    assert status["typed_inner_cartesian_center_velocity_dt_interface_registered"] is True
    assert (
        status["typed_inner_cartesian_center_velocity_dt_independent_audit_registered"]
        is True
    )
    assert status["inner_cartesian_center_velocity_dt_ingest_admitted"] is False
    assert status["inner_cartesian_center_derivative_identity_registered"] is True
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


def test_rejects_parent_or_sibling_evidence_laundering() -> None:
    for field, value in (
        ("parent_inner_cartesian_center_velocity_ingest_admitted", True),
        ("agent1_exact_head_ci_conclusion", "success"),
        ("agent4_exact_head_ci_conclusion", "success"),
        ("agent4_inner_cartesian_center_velocity_dt_independently_audited", True),
    ):
        payload = _payload()
        payload["evidence"][field] = value
        with pytest.raises(ValueError):
            validate_cartesian_center_velocity_dt_ingest_contract(payload)


def test_rejects_pressure_forcing_or_complete_api_invention() -> None:
    for field, value in (
        ("pressure", "invented.pressure"),
        ("forcing", "residual_defined_force"),
        ("complete_candidate_api_ready", True),
    ):
        payload = _payload()
        payload["candidate_api_handoff"][field] = value
        with pytest.raises(ValueError):
            validate_cartesian_center_velocity_dt_ingest_contract(payload)


def test_rejects_velocity_or_velocity_dt_binding_drift() -> None:
    payload = _payload()
    payload["candidate_api_handoff"]["velocity"] = "invented.velocity"
    with pytest.raises(ValueError):
        validate_cartesian_center_velocity_dt_ingest_contract(payload)

    payload = _payload()
    payload["candidate_api_handoff"]["velocity_dt"] = "invented.velocity_dt"
    with pytest.raises(ValueError):
        validate_cartesian_center_velocity_dt_ingest_contract(payload)


def test_rejects_scientific_state_promotion() -> None:
    for field in (
        "inner_cartesian_center_velocity_dt_ingest_admitted",
        "outer_global_leading_join_materialized",
        "global_leading_velocity_materialized",
        "matched_global_pressure_materialized",
        "restricted_forcing_composite_validated",
        "leading_ready",
    ):
        payload = _payload()
        payload["ingest_status"][field] = True
        with pytest.raises(ValueError):
            validate_cartesian_center_velocity_dt_ingest_contract(payload)

    for field in (
        "leading_ready",
        "correction_ready",
        "velocity_export_ready",
        "pde_validated",
    ):
        payload = _payload()
        payload["stage_state"][field] = True
        with pytest.raises(ValueError):
            validate_cartesian_center_velocity_dt_ingest_contract(payload)


def test_rejects_a4_protocol_or_pin_drift() -> None:
    payload = _payload()
    payload["agent1_velocity_dt_binding"]["source_blob_sha"] = "0" * 40
    with pytest.raises(ValueError):
        validate_cartesian_center_velocity_dt_ingest_contract(payload)

    payload = _payload()
    payload["agent4_velocity_dt_binding"]["head"] = "0" * 40
    with pytest.raises(ValueError):
        validate_cartesian_center_velocity_dt_ingest_contract(payload)

    payload = _payload()
    payload["agent4_velocity_dt_binding"]["fd4_time_steps"] = [4e-4, 2e-4, 2e-4]
    with pytest.raises(ValueError):
        validate_cartesian_center_velocity_dt_ingest_contract(payload)

    payload = _payload()
    payload["agent4_velocity_dt_binding"]["mutations"] = ["velocity_dt_x0.99"]
    with pytest.raises(ValueError):
        validate_cartesian_center_velocity_dt_ingest_contract(payload)


def test_rejects_baseline_gate_truth_or_digest_laundering() -> None:
    payload = _payload()
    payload["baseline_vs_kokuno"]["comparison_performed"] = True
    with pytest.raises(ValueError):
        validate_cartesian_center_velocity_dt_ingest_contract(payload)

    payload = _payload()
    payload["final_project_gates_unchanged"]["normalized_momentum_max"] = 1.1e-3
    with pytest.raises(ValueError):
        validate_cartesian_center_velocity_dt_ingest_contract(payload)

    payload = _payload()
    payload["truth_boundary"]["pde_validated"] = True
    with pytest.raises(ValueError):
        validate_cartesian_center_velocity_dt_ingest_contract(payload)

    payload = _payload()
    payload["contract_sha256"] = "0" * 64
    with pytest.raises(ValueError):
        validate_cartesian_center_velocity_dt_ingest_contract(payload)


def test_digest_detects_other_mutation() -> None:
    payload = _payload()
    mutated = copy.deepcopy(payload)
    mutated["agent1_velocity_dt_binding"]["registered_time_interval"] = [0.25, 0.74]
    with pytest.raises(ValueError):
        validate_cartesian_center_velocity_dt_ingest_contract(mutated)
