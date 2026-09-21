from __future__ import annotations

import copy

import pytest

from openai_ns_reconstruction.audit_kokuno_source_sensitivity_upstream_input_scope import (
    audit_contract,
    load_contract,
    mechanics_only_upstream_forcing_jet_witness,
)


def test_contract_baseline_is_clean_without_live_replay() -> None:
    contract = load_contract()
    assert audit_contract(contract, check_live_files=False) == []


def test_contract_live_files_are_consistent() -> None:
    contract = load_contract()
    assert audit_contract(contract, check_live_files=True) == []


def test_upstream_forcing_jet_witness_is_nonvacuous() -> None:
    receipt = mechanics_only_upstream_forcing_jet_witness()
    assert receipt["shared_generated_state"] > 0.0
    assert receipt["generated_dr_state_a"] != receipt["generated_dr_state_b"]
    assert receipt[
        "generated_sensitivity_remains_conditional_on_upstream_forcing_jet"
    ] is True


@pytest.mark.parametrize(
    ("section", "key", "bad_value"),
    [
        ("parent_interface_truth", "caller_supplies_mode_state_t_m", True),
        ("parent_interface_truth", "caller_supplies_amplitude_dr_t_m", True),
        ("parent_interface_truth", "caller_supplies_amplitude_dz_t_m", True),
        ("parent_interface_truth", "caller_supplies_corrected_background", False),
        ("parent_interface_truth", "caller_supplies_mode_forcing_value", False),
        ("parent_interface_truth", "caller_supplies_full_normalized_dr_forcing", False),
        ("parent_interface_truth", "caller_supplies_full_normalized_dz_forcing", False),
        (
            "parent_interface_truth",
            "deterministic_corrected_background_provider_materialized",
            True,
        ),
        (
            "parent_interface_truth",
            "deterministic_source_forcing_provider_materialized",
            True,
        ),
        ("parent_interface_truth", "source_amplitude_mode_provider_complete", True),
        ("parent_interface_truth", "self_contained_velocity_xyzt_provider", True),
        ("parent_interface_truth", "complete_ns_residual_assessed", True),
        ("parent_interface_truth", "pde_validated", True),
        ("parent_interface_truth", "paper_exact", True),
        ("parent_interface_truth", "openai_field_identified", True),
        (
            "input_retirement",
            "generated_directional_sensitivities_are_conditional_on_upstream_forcing_jets",
            False,
        ),
        (
            "input_retirement",
            "generated_directional_sensitivities_are_conditional_on_corrected_background",
            False,
        ),
        (
            "representation_boundary",
            "generated_amplitude_directional_jet_is_not_upstream_input_free",
            False,
        ),
        (
            "representation_boundary",
            "internal_sensitivity_integration_is_not_deterministic_source_forcing_materialization",
            False,
        ),
        (
            "representation_boundary",
            "complete_curl_consumption_is_not_unified_cartesian_velocity_delivery",
            False,
        ),
        (
            "canonical_eq45_state_must_remain_independent",
            "velocity_export_ready",
            False,
        ),
        (
            "canonical_eq45_state_must_remain_independent",
            "pde_validated",
            True,
        ),
        ("canonical_cr001", "momentum_max", 0.01),
        ("canonical_cr001", "divergence_max", 1e-4),
    ],
)
def test_forbidden_mutations_fail_closed(
    section: str, key: str, bad_value: object
) -> None:
    contract = load_contract()
    mutated = copy.deepcopy(contract)
    mutated[section][key] = bad_value
    assert audit_contract(mutated, check_live_files=False)


@pytest.mark.parametrize(
    "key",
    [
        "generated_Dr_t_m_and_Dz_t_m_implies_no_output_determining_source_inputs_remain",
        "generated_amplitude_directional_jet_implies_source_amplitude_mode_provider_complete",
        "sensitivity_solver_implies_deterministic_source_forcing_provider",
        "complete_curl_consumes_generated_jet_implies_self_contained_velocity_xyzt_provider",
        "scoped_source_sensitivity_mechanics_implies_visual_correspondence",
        "scoped_source_sensitivity_mechanics_implies_complete_ns_residual",
        "scoped_source_sensitivity_mechanics_implies_pde_validated",
        "scoped_source_sensitivity_mechanics_implies_paper_exact",
        "scoped_source_sensitivity_mechanics_implies_openai_field_identified",
        "source_route_incompleteness_implies_eq45_velocity_export_not_ready",
    ],
)
def test_forbidden_promotion_flip_is_rejected(key: str) -> None:
    contract = load_contract()
    mutated = copy.deepcopy(contract)
    mutated["forbidden_promotions"][key] = True
    assert audit_contract(mutated, check_live_files=False)


def test_remaining_upstream_inputs_cannot_be_erased() -> None:
    contract = load_contract()
    mutated = copy.deepcopy(contract)
    mutated["input_retirement"]["still_caller_supplied_output_determining_inputs"] = []
    assert audit_contract(mutated, check_live_files=False)


def test_four_way_provenance_cannot_be_collapsed() -> None:
    contract = load_contract()
    mutated = copy.deepcopy(contract)
    mutated["four_way_provenance"]["public_source_fact"] = []
    assert audit_contract(mutated, check_live_files=False)
