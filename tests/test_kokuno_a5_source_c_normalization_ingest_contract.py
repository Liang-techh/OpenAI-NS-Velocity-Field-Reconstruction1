from __future__ import annotations

import copy

import pytest

from openai_ns_reconstruction.kokuno_a5_source_c_normalization_ingest_contract import (
    AGENT1_SOURCE_C_HEAD,
    AGENT1_SOURCE_C_SCHEMA,
    AGENT1_SOURCE_C_SOURCE_BLOB,
    REQUIRED_MEMBERS,
    SELECTED_SOURCE_C,
    SOURCE_FORMULAS,
    deterministic_source_c_normalization_ingest_contract,
    validate_source_c_normalization_ingest_contract,
)


def test_source_c_interface_registered_but_not_admitted() -> None:
    payload = deterministic_source_c_normalization_ingest_contract(exact_head="deadbeef")
    validate_source_c_normalization_ingest_contract(payload)

    binding = payload["agent1_source_C_binding"]
    assert binding["head"] == AGENT1_SOURCE_C_HEAD
    assert binding["source_blob_sha"] == AGENT1_SOURCE_C_SOURCE_BLOB
    assert binding["schema"] == AGENT1_SOURCE_C_SCHEMA
    assert binding["selected_source_C"] == SELECTED_SOURCE_C == 1000.0
    assert tuple(binding["required_members"]) == REQUIRED_MEMBERS
    assert binding["source_formulas"] == SOURCE_FORMULAS
    assert binding["selection_boundary"] == {
        "chosen_from_public_source_condition": True,
        "chosen_from_NS_residual": False,
        "chosen_from_CR001_energy": False,
        "chosen_from_ST006": False,
        "same_C_required_in_physical_F_and_outer_schedule": True,
    }

    assert payload["independent_audit_requirement"] == {
        "dedicated_agent4_source_C_audit_required": True,
        "prior_agent4_physical_mapping_audit_pr": 790,
        "prior_mapping_audit_is_source_C_authority": False,
        "source_C_audit_pr": None,
        "source_C_audit_head": None,
        "source_C_audit_dedicated_run": None,
    }
    assert payload["evidence"] == {
        "agent1_exact_head_ci_conclusion": None,
        "agent1_self_certificate_executed_and_passed": None,
        "agent4_dedicated_source_C_audit_present": False,
        "agent4_source_C_audit_conclusion": None,
        "agent4_source_C_independently_admitted": False,
    }

    normalization = payload["normalization_registration"]
    assert normalization["selected_source_C_registered"] == 1000.0
    assert normalization["agent1_source_C_self_certificate_interface_registered"] is True
    assert normalization["source_complex_C_normalization_independently_admitted"] is False
    assert normalization["residual_driven_C_tuning_allowed"] is False

    status = payload["ingest_status"]
    assert status["typed_source_C_normalized_physical_center_interface_registered"] is True
    assert status["source_C_normalized_physical_center_ingest_admitted"] is False
    assert status["source_complex_C_normalization_agent1_self_certificate_passed"] is False
    assert status["source_complex_C_normalization_independently_admitted"] is False
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
        ("evidence", "agent1_exact_head_ci_conclusion", "success"),
        ("evidence", "agent1_self_certificate_executed_and_passed", True),
        ("evidence", "agent4_dedicated_source_C_audit_present", True),
        ("evidence", "agent4_source_C_audit_conclusion", "success"),
        ("evidence", "agent4_source_C_independently_admitted", True),
        ("ingest_status", "source_C_normalized_physical_center_ingest_admitted", True),
        ("ingest_status", "source_complex_C_normalization_agent1_self_certificate_passed", True),
        ("ingest_status", "source_complex_C_normalization_independently_admitted", True),
        ("ingest_status", "leading_ready", True),
        ("stage_state", "leading_ready", True),
        ("stage_state", "velocity_export_ready", True),
        ("stage_state", "pde_validated", True),
        ("truth_boundary", "queued_agent1_ci_promoted", True),
        ("truth_boundary", "missing_agent4_source_C_audit_invented", True),
        ("truth_boundary", "prior_mapping_audit_reused_as_source_C_certificate", True),
        ("truth_boundary", "source_C_used_as_residual_amplitude_knob", True),
        ("truth_boundary", "free_residual_defined_forcing_allowed", True),
        ("truth_boundary", "threshold_relaxed", True),
    ],
)
def test_evidence_or_state_laundering_fails_closed(
    section: str, key: str, value: object
) -> None:
    payload = deterministic_source_c_normalization_ingest_contract()
    mutated = copy.deepcopy(payload)
    mutated[section][key] = value
    with pytest.raises(ValueError):
        validate_source_c_normalization_ingest_contract(mutated)


def test_old_mapping_audit_cannot_be_reused_as_source_c_authority() -> None:
    payload = deterministic_source_c_normalization_ingest_contract()
    mutated = copy.deepcopy(payload)
    mutated["independent_audit_requirement"]["prior_mapping_audit_is_source_C_authority"] = True
    with pytest.raises(ValueError, match="audit"):
        validate_source_c_normalization_ingest_contract(mutated)


def test_selected_c_cannot_be_repurposed_as_residual_or_energy_normalizer() -> None:
    payload = deterministic_source_c_normalization_ingest_contract()

    mutated = copy.deepcopy(payload)
    mutated["normalization_registration"]["residual_driven_C_tuning_allowed"] = True
    with pytest.raises(ValueError, match="normalization"):
        validate_source_c_normalization_ingest_contract(mutated)

    mutated = copy.deepcopy(payload)
    mutated["normalization_registration"]["source_C_used_as_CR001_energy_normalizer"] = True
    with pytest.raises(ValueError, match="normalization"):
        validate_source_c_normalization_ingest_contract(mutated)


def test_full_candidate_residual_and_st006_comparison_cannot_be_fabricated() -> None:
    payload = deterministic_source_c_normalization_ingest_contract()

    mutated = copy.deepcopy(payload)
    mutated["baseline_vs_kokuno"]["kokuno_same_protocol_full_candidate_residual_available"] = True
    with pytest.raises(ValueError, match="fabricated"):
        validate_source_c_normalization_ingest_contract(mutated)

    mutated = copy.deepcopy(payload)
    mutated["baseline_vs_kokuno"]["comparison_performed"] = True
    with pytest.raises(ValueError, match="ST006"):
        validate_source_c_normalization_ingest_contract(mutated)


def test_contract_digest_detects_metadata_tampering() -> None:
    payload = deterministic_source_c_normalization_ingest_contract()
    mutated = copy.deepcopy(payload)
    mutated["agent1_source_C_binding"]["required_members"].append("velocity")
    with pytest.raises(ValueError, match="required member|digest"):
        validate_source_c_normalization_ingest_contract(mutated)
