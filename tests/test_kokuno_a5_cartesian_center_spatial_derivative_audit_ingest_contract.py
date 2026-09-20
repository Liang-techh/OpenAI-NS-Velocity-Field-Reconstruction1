from __future__ import annotations

import copy

import pytest

from openai_ns_reconstruction.kokuno_a5_cartesian_center_spatial_derivative_audit_ingest_contract import (
    AGENT4_AUDIT_HEAD,
    AGENT4_AUDIT_SOURCE_BLOB,
    AGENT4_PROTOCOL,
    DIVERGENCE_NORM_SCOPE,
    SCHEMA,
    deterministic_cartesian_center_spatial_derivative_audit_ingest_contract,
    validate_cartesian_center_spatial_derivative_audit_ingest_contract,
)


def _payload() -> dict:
    return deterministic_cartesian_center_spatial_derivative_audit_ingest_contract(
        exact_head="deadbeef"
    )


def test_registers_real_agent4_audit_without_promoting_it() -> None:
    payload = _payload()
    validate_cartesian_center_spatial_derivative_audit_ingest_contract(payload)

    assert payload["schema"] == SCHEMA
    a4 = payload["agent4_spatial_derivative_audit_binding"]
    assert a4["pr"] == 823
    assert a4["head"] == AGENT4_AUDIT_HEAD
    assert a4["source_blob_sha"] == AGENT4_AUDIT_SOURCE_BLOB
    assert a4["protocol"] == AGENT4_PROTOCOL
    assert a4["divergence_norm_scope"] == DIVERGENCE_NORM_SCOPE

    evidence = payload["evidence"]
    assert evidence["agent4_dedicated_spatial_derivative_audit_present"] is True
    assert evidence["agent4_exact_head_ci_conclusion"] is None
    assert evidence["agent4_independent_audit_conclusion"] is None
    assert evidence["agent4_scientific_receipt_admitted"] is False

    ops = payload["differential_operator_handoff"]
    assert ops["independent_audit_available"] is True
    assert ops["independent_audit_status"] == "registered_unresolved"
    assert ops["independent_audit_authority"] == "Agent4#823"
    assert ops["usable_for_inner_spatial_derivative_admission_if_passes"] is True
    assert ops["usable_for_final_independent_pde_validation"] is False

    status = payload["ingest_status"]
    assert (
        status[
            "typed_inner_cartesian_center_spatial_derivative_independent_audit_registered"
        ]
        is True
    )
    assert status["inner_cartesian_center_spatial_derivative_ingest_admitted"] is False
    assert status["leading_ready"] is False

    api = payload["candidate_api_handoff"]
    assert api["pressure"] is None
    assert api["forcing"] is None
    assert api["complete_candidate_api_ready"] is False
    assert payload["stage_state"] == {
        "leading_ready": False,
        "oscillatory_ready": True,
        "correction_ready": False,
        "velocity_export_ready": False,
        "pde_validated": False,
    }
    assert payload["baseline_vs_kokuno"]["comparison_performed"] is False


def test_rejects_agent4_ci_or_scientific_admission_laundering() -> None:
    mutations = (
        ("parent_inner_cartesian_center_spatial_derivative_ingest_admitted", True),
        ("agent4_exact_head_ci_conclusion", "success"),
        ("agent4_independent_audit_conclusion", "pass"),
        ("agent4_scientific_receipt_admitted", True),
    )
    for field, value in mutations:
        payload = _payload()
        payload["evidence"][field] = value
        with pytest.raises(ValueError):
            validate_cartesian_center_spatial_derivative_audit_ingest_contract(payload)


def test_rejects_agent4_pin_or_protocol_drift() -> None:
    payload = _payload()
    payload["agent4_spatial_derivative_audit_binding"]["head"] = "0" * 40
    with pytest.raises(ValueError):
        validate_cartesian_center_spatial_derivative_audit_ingest_contract(payload)

    payload = _payload()
    payload["agent4_spatial_derivative_audit_binding"]["source_blob_sha"] = "0" * 40
    with pytest.raises(ValueError):
        validate_cartesian_center_spatial_derivative_audit_ingest_contract(payload)

    payload = _payload()
    payload["agent4_spatial_derivative_audit_binding"]["protocol"][
        "jacobian_relative_max_gate"
    ] = 3e-5
    with pytest.raises(ValueError):
        validate_cartesian_center_spatial_derivative_audit_ingest_contract(payload)


def test_rejects_divergence_norm_scope_laundering() -> None:
    payload = _payload()
    scope = payload["agent4_spatial_derivative_audit_binding"]["divergence_norm_scope"]
    scope["canonical_cr001_volume_weighted_spatial_l2_established"] = True
    with pytest.raises(ValueError):
        validate_cartesian_center_spatial_derivative_audit_ingest_contract(payload)

    payload = _payload()
    scope = payload["agent4_spatial_derivative_audit_binding"]["divergence_norm_scope"]
    scope["same_metric_as_final_cr001_divergence_l2"] = True
    with pytest.raises(ValueError):
        validate_cartesian_center_spatial_derivative_audit_ingest_contract(payload)


def test_rejects_inner_audit_promotion_to_final_validator() -> None:
    payload = _payload()
    payload["differential_operator_handoff"][
        "usable_for_final_independent_pde_validation"
    ] = True
    with pytest.raises(ValueError):
        validate_cartesian_center_spatial_derivative_audit_ingest_contract(payload)

    payload = _payload()
    payload["differential_operator_handoff"]["independent_audit_status"] = "pass"
    with pytest.raises(ValueError):
        validate_cartesian_center_spatial_derivative_audit_ingest_contract(payload)

    payload = _payload()
    payload["ingest_status"][
        "inner_cartesian_center_spatial_derivative_ingest_admitted"
    ] = True
    with pytest.raises(ValueError):
        validate_cartesian_center_spatial_derivative_audit_ingest_contract(payload)


def test_rejects_candidate_stage_baseline_or_gate_promotion() -> None:
    payload = _payload()
    payload["candidate_api_handoff"]["pressure"] = "invented.pressure"
    with pytest.raises(ValueError):
        validate_cartesian_center_spatial_derivative_audit_ingest_contract(payload)

    payload = _payload()
    payload["candidate_api_handoff"]["forcing"] = "residual_defined_force"
    with pytest.raises(ValueError):
        validate_cartesian_center_spatial_derivative_audit_ingest_contract(payload)

    payload = _payload()
    payload["stage_state"]["pde_validated"] = True
    with pytest.raises(ValueError):
        validate_cartesian_center_spatial_derivative_audit_ingest_contract(payload)

    payload = _payload()
    payload["baseline_vs_kokuno"]["comparison_performed"] = True
    with pytest.raises(ValueError):
        validate_cartesian_center_spatial_derivative_audit_ingest_contract(payload)

    payload = _payload()
    payload["final_project_gates_unchanged"]["normalized_momentum_max"] = 1.1e-3
    with pytest.raises(ValueError):
        validate_cartesian_center_spatial_derivative_audit_ingest_contract(payload)


def test_rejects_truth_or_digest_laundering() -> None:
    payload = _payload()
    payload["truth_boundary"]["pde_validated"] = True
    with pytest.raises(ValueError):
        validate_cartesian_center_spatial_derivative_audit_ingest_contract(payload)

    payload = _payload()
    payload["truth_boundary"][
        "sampled_divergence_rms_laundered_as_cr001_volume_l2"
    ] = True
    with pytest.raises(ValueError):
        validate_cartesian_center_spatial_derivative_audit_ingest_contract(payload)

    payload = _payload()
    payload["contract_sha256"] = "0" * 64
    with pytest.raises(ValueError):
        validate_cartesian_center_spatial_derivative_audit_ingest_contract(payload)


def test_digest_detects_other_mutation() -> None:
    payload = _payload()
    mutated = copy.deepcopy(payload)
    mutated["agent4_spatial_derivative_audit_binding"]["tests_run"] += 1
    with pytest.raises(ValueError):
        validate_cartesian_center_spatial_derivative_audit_ingest_contract(mutated)
