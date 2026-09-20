from __future__ import annotations

import copy

import pytest

from openai_ns_reconstruction.kokuno_a5_cartesian_center_spatial_derivative_ingest_contract import (
    AGENT1_SPATIAL_DERIVATIVE_HEAD,
    AGENT1_SPATIAL_DERIVATIVE_SOURCE_BLOB,
    SCHEMA,
    deterministic_cartesian_center_spatial_derivative_ingest_contract,
    validate_cartesian_center_spatial_derivative_ingest_contract,
)


def _payload() -> dict:
    return deterministic_cartesian_center_spatial_derivative_ingest_contract(
        exact_head="deadbeef"
    )


def test_registers_spatial_derivative_surface_but_remains_fail_closed() -> None:
    payload = _payload()
    validate_cartesian_center_spatial_derivative_ingest_contract(payload)

    assert payload["schema"] == SCHEMA
    a1 = payload["agent1_spatial_derivative_binding"]
    assert a1["head"] == AGENT1_SPATIAL_DERIVATIVE_HEAD
    assert a1["source_blob_sha"] == AGENT1_SPATIAL_DERIVATIVE_SOURCE_BLOB
    assert a1["identity_boundary"] == {
        "underlying_field_sha256_preserved": True,
        "temporal_derivative_sha256_preserved": True,
        "spatial_derivative_sha256_separate": True,
        "production_finite_difference_used": False,
    }

    evidence = payload["evidence"]
    assert evidence["parent_inner_cartesian_center_velocity_dt_ingest_admitted"] is False
    assert evidence["agent1_exact_head_ci_conclusion"] is None
    assert evidence["agent4_dedicated_spatial_derivative_audit_present"] is False
    assert (
        evidence[
            "agent4_inner_cartesian_center_spatial_derivatives_independently_audited"
        ]
        is False
    )
    assert payload["agent4_spatial_derivative_audit_binding"] is None

    parent_api = payload["candidate_api_handoff"]
    assert parent_api["velocity"] == "Agent1#803.KokunoPA10CartesianCenterVelocity.velocity"
    assert (
        parent_api["velocity_dt"]
        == "Agent1#811.KokunoPA10CartesianCenterVelocityTimeDerivative.velocity_dt"
    )
    assert parent_api["pressure"] is None
    assert parent_api["forcing"] is None
    assert parent_api["complete_candidate_api_ready"] is False

    ops = payload["differential_operator_handoff"]
    assert ops["velocity_jacobian"].endswith(".velocity_jacobian")
    assert ops["divergence"].endswith(".divergence")
    assert ops["vorticity"].endswith(".vorticity")
    assert ops["independent_audit_available"] is False
    assert ops["usable_for_final_independent_pde_validation"] is False

    status = payload["ingest_status"]
    assert status["typed_inner_cartesian_center_spatial_derivative_interface_registered"] is True
    assert (
        status[
            "typed_inner_cartesian_center_spatial_derivative_independent_audit_registered"
        ]
        is False
    )
    assert status["inner_cartesian_center_spatial_derivative_ingest_admitted"] is False
    assert status["leading_ready"] is False

    assert payload["stage_state"] == {
        "leading_ready": False,
        "oscillatory_ready": True,
        "correction_ready": False,
        "velocity_export_ready": False,
        "pde_validated": False,
    }
    assert payload["baseline_vs_kokuno"]["comparison_performed"] is False


def test_rejects_parent_or_agent1_evidence_laundering() -> None:
    for field, value in (
        ("parent_inner_cartesian_center_velocity_dt_ingest_admitted", True),
        ("agent1_exact_head_ci_conclusion", "success"),
        ("agent4_dedicated_spatial_derivative_audit_present", True),
        ("agent4_inner_cartesian_center_spatial_derivatives_independently_audited", True),
    ):
        payload = _payload()
        payload["evidence"][field] = value
        with pytest.raises(ValueError):
            validate_cartesian_center_spatial_derivative_ingest_contract(payload)


def test_rejects_invented_agent4_audit_or_validator_promotion() -> None:
    payload = _payload()
    payload["agent4_spatial_derivative_audit_binding"] = {"pr": 999}
    with pytest.raises(ValueError):
        validate_cartesian_center_spatial_derivative_ingest_contract(payload)

    payload = _payload()
    payload["differential_operator_handoff"]["independent_audit_available"] = True
    with pytest.raises(ValueError):
        validate_cartesian_center_spatial_derivative_ingest_contract(payload)

    payload = _payload()
    payload["differential_operator_handoff"]["usable_for_final_independent_pde_validation"] = True
    with pytest.raises(ValueError):
        validate_cartesian_center_spatial_derivative_ingest_contract(payload)


def test_rejects_spatial_operator_or_pin_drift() -> None:
    payload = _payload()
    payload["agent1_spatial_derivative_binding"]["source_blob_sha"] = "0" * 40
    with pytest.raises(ValueError):
        validate_cartesian_center_spatial_derivative_ingest_contract(payload)

    payload = _payload()
    payload["agent1_spatial_derivative_binding"]["engineering_protocol"][
        "velocity_jacobian_vs_fd4_fine_relative_max_gate"
    ] = 6e-5
    with pytest.raises(ValueError):
        validate_cartesian_center_spatial_derivative_ingest_contract(payload)

    payload = _payload()
    payload["differential_operator_handoff"]["divergence"] = "invented.divergence"
    with pytest.raises(ValueError):
        validate_cartesian_center_spatial_derivative_ingest_contract(payload)


def test_rejects_candidate_api_or_scientific_state_promotion() -> None:
    for field, value in (
        ("pressure", "invented.pressure"),
        ("forcing", "residual_defined_force"),
        ("complete_candidate_api_ready", True),
    ):
        payload = _payload()
        payload["candidate_api_handoff"][field] = value
        with pytest.raises(ValueError):
            validate_cartesian_center_spatial_derivative_ingest_contract(payload)

    for field in (
        "inner_cartesian_center_spatial_derivative_ingest_admitted",
        "outer_global_leading_join_materialized",
        "global_leading_velocity_materialized",
        "matched_global_pressure_materialized",
        "restricted_forcing_composite_validated",
        "leading_ready",
    ):
        payload = _payload()
        payload["ingest_status"][field] = True
        with pytest.raises(ValueError):
            validate_cartesian_center_spatial_derivative_ingest_contract(payload)

    for field in (
        "leading_ready",
        "correction_ready",
        "velocity_export_ready",
        "pde_validated",
    ):
        payload = _payload()
        payload["stage_state"][field] = True
        with pytest.raises(ValueError):
            validate_cartesian_center_spatial_derivative_ingest_contract(payload)


def test_rejects_baseline_gate_truth_or_digest_laundering() -> None:
    payload = _payload()
    payload["baseline_vs_kokuno"]["comparison_performed"] = True
    with pytest.raises(ValueError):
        validate_cartesian_center_spatial_derivative_ingest_contract(payload)

    payload = _payload()
    payload["final_project_gates_unchanged"]["normalized_momentum_max"] = 1.1e-3
    with pytest.raises(ValueError):
        validate_cartesian_center_spatial_derivative_ingest_contract(payload)

    payload = _payload()
    payload["truth_boundary"]["pde_validated"] = True
    with pytest.raises(ValueError):
        validate_cartesian_center_spatial_derivative_ingest_contract(payload)

    payload = _payload()
    payload["contract_sha256"] = "0" * 64
    with pytest.raises(ValueError):
        validate_cartesian_center_spatial_derivative_ingest_contract(payload)


def test_digest_detects_other_mutation() -> None:
    payload = _payload()
    mutated = copy.deepcopy(payload)
    mutated["agent1_spatial_derivative_binding"]["registered_time_interval"] = [0.25, 0.74]
    with pytest.raises(ValueError):
        validate_cartesian_center_spatial_derivative_ingest_contract(mutated)
