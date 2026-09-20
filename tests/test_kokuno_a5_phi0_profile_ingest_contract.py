from __future__ import annotations

import copy

import pytest

from openai_ns_reconstruction.kokuno_a5_phi0_profile_ingest_contract import (
    AGENT1_PHI0_HEAD,
    AGENT1_PHI0_SCHEMA,
    AGENT1_PHI0_SOURCE_BLOB,
    REQUIRED_DERIVATIVE_FIELDS,
    REQUIRED_METHODS,
    REQUIRED_VALUE_FIELDS,
    deterministic_phi0_ingest_contract,
    validate_phi0_ingest_contract,
)


def test_phi0_profile_ingest_is_registered_but_fail_closed() -> None:
    payload = deterministic_phi0_ingest_contract(exact_head="deadbeef")
    validate_phi0_ingest_contract(payload)

    binding = payload["agent1_phi0_binding"]
    assert binding["head"] == AGENT1_PHI0_HEAD
    assert binding["source_blob_sha"] == AGENT1_PHI0_SOURCE_BLOB
    assert binding["schema"] == AGENT1_PHI0_SCHEMA
    assert tuple(binding["required_methods"]) == REQUIRED_METHODS
    assert tuple(binding["required_value_fields"]) == REQUIRED_VALUE_FIELDS
    assert tuple(binding["required_derivative_fields"]) == REQUIRED_DERIVATIVE_FIELDS
    assert binding["coordinate_contract"]["native_similarity_coordinates"] == ["Y", "eta"]
    assert binding["coordinate_contract"]["source_Y_interval"] == [0.0, 4.1]
    assert binding["coordinate_contract"]["repository_X_identification_asserted"] is False
    assert binding["coordinate_contract"]["physical_F_from_Phi_mapping_resolved"] is False
    assert binding["numerical_realization"]["f0_series_terms"] == 64

    assert payload["evidence"] == {
        "agent1_exact_head_ci_passed": False,
        "agent4_phi0_independent_audit_present": False,
        "agent4_phi0_independent_audit_passed": False,
    }
    assert payload["ingest_status"] == {
        "typed_axis_profile_interface_registered": True,
        "axis_profile_ingest_admitted": False,
        "typed_phi0_profile_interface_registered": True,
        "phi0_profile_ingest_admitted": False,
        "repository_X_to_source_Y_mapping_resolved": False,
        "physical_F_from_Phi_mapping_resolved": False,
        "fixed_point_correction_materialized": False,
        "global_leading_velocity_materialized": False,
        "leading_ready": False,
    }
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
        ("ingest_status", "phi0_profile_ingest_admitted", True),
        ("ingest_status", "leading_ready", True),
        ("stage_state", "leading_ready", True),
        ("stage_state", "velocity_export_ready", True),
        ("stage_state", "pde_validated", True),
        ("truth_boundary", "queued_or_missing_evidence_promoted", True),
        ("truth_boundary", "repository_X_identified_with_source_Y", True),
        ("truth_boundary", "physical_F_from_Phi_mapping_invented", True),
        ("truth_boundary", "Phi0_promoted_to_final_corrected_Phi", True),
        ("truth_boundary", "free_residual_defined_forcing_allowed", True),
        ("truth_boundary", "threshold_relaxed", True),
    ],
)
def test_state_laundering_fails_closed(section: str, key: str, value: object) -> None:
    payload = deterministic_phi0_ingest_contract()
    mutated = copy.deepcopy(payload)
    mutated[section][key] = value
    with pytest.raises(ValueError):
        validate_phi0_ingest_contract(mutated)


def test_sibling_evidence_and_coordinate_mapping_cannot_be_invented() -> None:
    payload = deterministic_phi0_ingest_contract()

    mutated = copy.deepcopy(payload)
    mutated["evidence"]["agent1_exact_head_ci_passed"] = True
    with pytest.raises(ValueError, match="evidence"):
        validate_phi0_ingest_contract(mutated)

    mutated = copy.deepcopy(payload)
    mutated["evidence"]["agent4_phi0_independent_audit_present"] = True
    with pytest.raises(ValueError, match="evidence"):
        validate_phi0_ingest_contract(mutated)

    mutated = copy.deepcopy(payload)
    mutated["agent1_phi0_binding"]["coordinate_contract"][
        "repository_X_identification_asserted"
    ] = True
    with pytest.raises(ValueError, match="X/Y"):
        validate_phi0_ingest_contract(mutated)

    mutated = copy.deepcopy(payload)
    mutated["agent1_phi0_binding"]["coordinate_contract"][
        "physical_F_from_Phi_mapping_resolved"
    ] = True
    with pytest.raises(ValueError, match="physical F/Phi"):
        validate_phi0_ingest_contract(mutated)


def test_full_candidate_residual_and_baseline_comparison_cannot_be_fabricated() -> None:
    payload = deterministic_phi0_ingest_contract()

    mutated = copy.deepcopy(payload)
    mutated["baseline_vs_kokuno"][
        "kokuno_same_protocol_full_candidate_residual_available"
    ] = True
    with pytest.raises(ValueError, match="fabricated"):
        validate_phi0_ingest_contract(mutated)

    mutated = copy.deepcopy(payload)
    mutated["baseline_vs_kokuno"]["comparison_performed"] = True
    with pytest.raises(ValueError, match="ST006"):
        validate_phi0_ingest_contract(mutated)


def test_contract_digest_detects_metadata_tampering() -> None:
    payload = deterministic_phi0_ingest_contract()
    mutated = copy.deepcopy(payload)
    mutated["agent1_phi0_binding"]["required_methods"].append("velocity")
    with pytest.raises(ValueError, match="digest"):
        validate_phi0_ingest_contract(mutated)
