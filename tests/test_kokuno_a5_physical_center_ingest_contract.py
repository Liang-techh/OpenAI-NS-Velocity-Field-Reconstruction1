from __future__ import annotations

import copy

import pytest

from openai_ns_reconstruction.kokuno_a5_physical_center_ingest_contract import (
    AGENT1_PHYSICAL_HEAD,
    AGENT1_PHYSICAL_SCHEMA,
    AGENT1_PHYSICAL_SOURCE_BLOB,
    AGENT4_PHYSICAL_AUDIT_HEAD,
    AGENT4_PHYSICAL_AUDIT_SCHEMA,
    AGENT4_PHYSICAL_AUDIT_SOURCE_BLOB,
    REQUIRED_DERIVATIVE_FIELDS,
    REQUIRED_METHODS,
    REQUIRED_VALUE_FIELDS,
    SOURCE_FORMULAS,
    deterministic_physical_center_ingest_contract,
    validate_physical_center_ingest_contract,
)


def test_physical_center_interface_is_registered_but_not_admitted() -> None:
    payload = deterministic_physical_center_ingest_contract(exact_head="deadbeef")
    validate_physical_center_ingest_contract(payload)

    binding = payload["agent1_physical_center_binding"]
    assert binding["head"] == AGENT1_PHYSICAL_HEAD
    assert binding["source_blob_sha"] == AGENT1_PHYSICAL_SOURCE_BLOB
    assert binding["schema"] == AGENT1_PHYSICAL_SCHEMA
    assert tuple(binding["required_methods"]) == REQUIRED_METHODS
    assert tuple(binding["required_value_fields"]) == REQUIRED_VALUE_FIELDS
    assert tuple(binding["required_derivative_fields"]) == REQUIRED_DERIVATIVE_FIELDS
    assert binding["source_formulas"] == SOURCE_FORMULAS
    assert binding["scope"] == {
        "contraction_center_only": True,
        "final_corrected_fixed_point": False,
        "cartesian_spacetime_velocity": False,
        "matched_global_pressure": False,
    }

    audit = payload["agent4_independent_audit_binding"]
    assert audit["head"] == AGENT4_PHYSICAL_AUDIT_HEAD
    assert audit["source_blob_sha"] == AGENT4_PHYSICAL_AUDIT_SOURCE_BLOB
    assert audit["schema"] == AGENT4_PHYSICAL_AUDIT_SCHEMA

    assert payload["evidence"] == {
        "agent1_exact_head_ci_conclusion": None,
        "agent4_audit_present": True,
        "agent4_exact_head_ci_conclusion": None,
        "agent4_mapping_formula_independently_consistent": None,
        "agent4_configured_C_real_axis_obstruction_detected": None,
        "agent4_physical_center_profile_independently_admitted": False,
        "source_complex_C_normalization_certified": False,
    }
    assert payload["mapping_registration"] == {
        "source_Y_equals_Lambda_X_formula_registered": True,
        "physical_F_from_Phi_formula_registered": True,
        "physical_U_from_u_formula_registered": True,
        "repository_X_equals_source_Y_asserted": False,
        "mapping_independently_admitted": False,
    }

    status = payload["ingest_status"]
    assert status["typed_physical_center_profile_interface_registered"] is True
    assert status["physical_center_profile_ingest_admitted"] is False
    assert status["source_Y_equals_Lambda_X_formula_registered"] is True
    assert status["physical_F_from_Phi_formula_registered"] is True
    assert status["physical_U_from_u_formula_registered"] is True
    assert status["source_complex_C_normalization_certified"] is False
    assert status["global_leading_velocity_materialized"] is False
    assert status["matched_global_pressure_materialized"] is False
    assert status["restricted_forcing_composite_validated"] is False
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
        ("ingest_status", "physical_center_profile_ingest_admitted", True),
        ("ingest_status", "source_complex_C_normalization_certified", True),
        ("ingest_status", "global_leading_velocity_materialized", True),
        ("ingest_status", "matched_global_pressure_materialized", True),
        ("ingest_status", "restricted_forcing_composite_validated", True),
        ("ingest_status", "leading_ready", True),
        ("stage_state", "leading_ready", True),
        ("stage_state", "velocity_export_ready", True),
        ("stage_state", "pde_validated", True),
        ("truth_boundary", "queued_sibling_ci_promoted", True),
        ("truth_boundary", "agent4_unknown_result_invented", True),
        ("truth_boundary", "repository_X_equals_source_Y_asserted", True),
        ("truth_boundary", "source_complex_C_normalization_claimed", True),
        ("truth_boundary", "free_residual_defined_forcing_allowed", True),
        ("truth_boundary", "threshold_relaxed", True),
    ],
)
def test_state_laundering_fails_closed(section: str, key: str, value: object) -> None:
    payload = deterministic_physical_center_ingest_contract()
    mutated = copy.deepcopy(payload)
    mutated[section][key] = value
    with pytest.raises(ValueError):
        validate_physical_center_ingest_contract(mutated)


def test_queued_sibling_evidence_cannot_be_invented() -> None:
    payload = deterministic_physical_center_ingest_contract()

    for key, value in (
        ("agent1_exact_head_ci_conclusion", "success"),
        ("agent4_exact_head_ci_conclusion", "success"),
        ("agent4_mapping_formula_independently_consistent", True),
        ("agent4_configured_C_real_axis_obstruction_detected", False),
        ("agent4_physical_center_profile_independently_admitted", True),
        ("source_complex_C_normalization_certified", True),
    ):
        mutated = copy.deepcopy(payload)
        mutated["evidence"][key] = value
        with pytest.raises(ValueError, match="evidence"):
            validate_physical_center_ingest_contract(mutated)


def test_formula_registration_is_not_x_equals_y_or_scientific_admission() -> None:
    payload = deterministic_physical_center_ingest_contract()

    mutated = copy.deepcopy(payload)
    mutated["mapping_registration"]["repository_X_equals_source_Y_asserted"] = True
    with pytest.raises(ValueError, match="mapping"):
        validate_physical_center_ingest_contract(mutated)

    mutated = copy.deepcopy(payload)
    mutated["mapping_registration"]["mapping_independently_admitted"] = True
    with pytest.raises(ValueError, match="mapping"):
        validate_physical_center_ingest_contract(mutated)

    mutated = copy.deepcopy(payload)
    mutated["agent1_physical_center_binding"]["source_formulas"][
        "source_radial_mapping"
    ] = "Y=X"
    with pytest.raises(ValueError, match="formula"):
        validate_physical_center_ingest_contract(mutated)


def test_contraction_center_cannot_be_promoted_to_final_global_candidate() -> None:
    payload = deterministic_physical_center_ingest_contract()

    for key in (
        "contraction_center_promoted_to_final_fixed_point",
        "source_native_center_promoted_to_cartesian_velocity",
        "matched_global_pressure_invented",
        "restricted_forcing_composite_invented",
    ):
        mutated = copy.deepcopy(payload)
        mutated["truth_boundary"][key] = True
        with pytest.raises(ValueError, match="truth-boundary"):
            validate_physical_center_ingest_contract(mutated)


def test_full_candidate_residual_and_baseline_comparison_cannot_be_fabricated() -> None:
    payload = deterministic_physical_center_ingest_contract()

    mutated = copy.deepcopy(payload)
    mutated["baseline_vs_kokuno"][
        "kokuno_same_protocol_full_candidate_residual_available"
    ] = True
    with pytest.raises(ValueError, match="fabricated"):
        validate_physical_center_ingest_contract(mutated)

    mutated = copy.deepcopy(payload)
    mutated["baseline_vs_kokuno"]["comparison_performed"] = True
    with pytest.raises(ValueError, match="ST006"):
        validate_physical_center_ingest_contract(mutated)


def test_contract_digest_detects_metadata_tampering() -> None:
    payload = deterministic_physical_center_ingest_contract()
    mutated = copy.deepcopy(payload)
    mutated["agent1_physical_center_binding"]["required_methods"].append("velocity")
    with pytest.raises(ValueError, match="digest"):
        validate_physical_center_ingest_contract(mutated)
