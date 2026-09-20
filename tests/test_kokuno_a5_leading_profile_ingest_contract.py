from __future__ import annotations

import copy

import pytest

from openai_ns_reconstruction.kokuno_a5_leading_profile_ingest_contract import (
    AGENT1_PROFILE_HEAD,
    AGENT1_PROFILE_SCHEMA,
    AGENT4_PROFILE_AUDIT_HEAD,
    REQUIRED_DERIVATIVE_FIELDS,
    REQUIRED_PROFILE_METHODS,
    REQUIRED_VALUE_FIELDS,
    deterministic_ingest_contract,
    validate_ingest_contract,
)


def test_leading_profile_ingest_contract_is_registered_but_not_admitted() -> None:
    payload = deterministic_ingest_contract(exact_head="deadbeef")
    validate_ingest_contract(payload)

    assert payload["pipeline_stage"] == "source_profile_ingest"
    assert payload["pipeline_order"][0] == "source_profile_ingest"
    binding = payload["agent1_profile_binding"]
    assert binding["head"] == AGENT1_PROFILE_HEAD
    assert binding["schema"] == AGENT1_PROFILE_SCHEMA
    assert tuple(binding["required_methods"]) == REQUIRED_PROFILE_METHODS
    assert tuple(binding["required_value_fields"]) == REQUIRED_VALUE_FIELDS
    assert tuple(binding["required_derivative_fields"]) == REQUIRED_DERIVATIVE_FIELDS
    assert binding["coordinate_contract"]["native_similarity_coordinates"] == ["Y", "eta"]
    assert binding["coordinate_contract"]["source_Y_interval"] == [0.0, 4.1]
    assert binding["coordinate_contract"]["repository_X_identification_asserted"] is False

    evidence = payload["independent_evidence_pin"]
    assert evidence["agent1_head"] == AGENT1_PROFILE_HEAD
    assert evidence["agent4_head"] == AGENT4_PROFILE_AUDIT_HEAD
    assert evidence["agent1_exact_head_ci_passed"] is False
    assert evidence["agent4_independent_audit_passed"] is False

    status = payload["ingest_status"]
    assert status["typed_profile_interface_registered"] is True
    assert status["source_profile_ingest_admitted"] is False
    assert status["repository_X_to_source_Y_mapping_resolved"] is False
    assert status["global_leading_velocity_materialized"] is False
    assert status["leading_ready"] is False

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


@pytest.mark.parametrize(
    ("section", "key", "value"),
    [
        ("ingest_status", "source_profile_ingest_admitted", True),
        ("ingest_status", "leading_ready", True),
        ("stage_state", "leading_ready", True),
        ("stage_state", "velocity_export_ready", True),
        ("stage_state", "pde_validated", True),
        ("truth_boundary", "agent1_queued_ci_laundered_to_pass", True),
        ("truth_boundary", "agent4_queued_audit_laundered_to_pass", True),
        ("truth_boundary", "repository_X_identified_with_source_Y", True),
        ("truth_boundary", "free_residual_defined_forcing_allowed", True),
        ("truth_boundary", "threshold_relaxed", True),
    ],
)
def test_state_laundering_fails_closed(section: str, key: str, value: object) -> None:
    payload = deterministic_ingest_contract()
    mutated = copy.deepcopy(payload)
    mutated[section][key] = value
    with pytest.raises(ValueError):
        validate_ingest_contract(mutated)


def test_sibling_evidence_cannot_be_mutated_by_caller() -> None:
    payload = deterministic_ingest_contract()
    mutated = copy.deepcopy(payload)
    mutated["independent_evidence_pin"]["agent1_exact_head_ci_passed"] = True
    with pytest.raises(ValueError, match="evidence"):
        validate_ingest_contract(mutated)

    mutated = copy.deepcopy(payload)
    mutated["independent_evidence_pin"]["agent4_independent_audit_passed"] = True
    with pytest.raises(ValueError, match="evidence"):
        validate_ingest_contract(mutated)


def test_coordinate_mapping_and_profile_head_drift_fail_closed() -> None:
    payload = deterministic_ingest_contract()
    mutated = copy.deepcopy(payload)
    mutated["agent1_profile_binding"]["coordinate_contract"][
        "repository_X_identification_asserted"
    ] = True
    with pytest.raises(ValueError, match="mapping"):
        validate_ingest_contract(mutated)

    mutated = copy.deepcopy(payload)
    mutated["agent1_profile_binding"]["head"] = "0" * 40
    with pytest.raises(ValueError, match="exact-head"):
        validate_ingest_contract(mutated)


def test_full_candidate_residual_or_baseline_comparison_cannot_be_fabricated() -> None:
    payload = deterministic_ingest_contract()
    mutated = copy.deepcopy(payload)
    mutated["baseline_vs_kokuno"][
        "kokuno_same_protocol_full_candidate_residual_available"
    ] = True
    with pytest.raises(ValueError, match="fabricated"):
        validate_ingest_contract(mutated)

    mutated = copy.deepcopy(payload)
    mutated["baseline_vs_kokuno"]["comparison_performed"] = True
    with pytest.raises(ValueError, match="ST006"):
        validate_ingest_contract(mutated)


def test_contract_digest_detects_metadata_tampering() -> None:
    payload = deterministic_ingest_contract()
    mutated = copy.deepcopy(payload)
    mutated["agent1_profile_binding"]["required_methods"].append("velocity")
    with pytest.raises(ValueError, match="digest"):
        validate_ingest_contract(mutated)
