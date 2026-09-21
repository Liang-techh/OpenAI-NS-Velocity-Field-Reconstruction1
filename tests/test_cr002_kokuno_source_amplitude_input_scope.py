from __future__ import annotations

import copy

import pytest

from openai_ns_reconstruction.audit_kokuno_source_amplitude_input_scope import (
    audit_contract,
    load_contract,
    mechanics_only_amplitude_dependence_witness,
)


def test_live_contract_is_clean():
    contract = load_contract()
    assert audit_contract(contract) == []


@pytest.mark.parametrize(
    ("section", "key", "value"),
    [
        ("parent_interface_truth", "caller_supplies_background_phase_derivative_jet", False),
        ("parent_interface_truth", "caller_supplies_amplitude_directional_jet", False),
        ("parent_interface_truth", "source_amplitude_ode_materialized", True),
        ("parent_interface_truth", "self_contained_velocity_xyzt_provider", True),
        ("parent_interface_truth", "source_to_runtime_parameter_map_complete", True),
        ("parent_interface_truth", "pde_validated", True),
        ("forbidden_promotions", "no_caller_coefficient_jets_implies_no_caller_source_jets", True),
        ("forbidden_promotions", "directional_transversality_check_implies_source_amplitude_ode_solved", True),
        ("forbidden_promotions", "manufactured_transverse_amplitude_implies_source_amplitude_realization", True),
        ("forbidden_promotions", "analytic_phase_vector_jet_implies_background_realization_complete", True),
        ("forbidden_promotions", "annular_source_formula_implies_axis_safe_global_cartesian_runtime", True),
        ("forbidden_promotions", "source_formula_execution_implies_velocity_export_ready", True),
        ("forbidden_promotions", "source_formula_execution_implies_visual_correspondence_verified", True),
        ("forbidden_promotions", "source_formula_execution_implies_paper_exact", True),
        ("forbidden_promotions", "source_formula_execution_implies_openai_field_identified", True),
    ],
)
def test_forbidden_scope_promotions_fail_closed(section, key, value):
    contract = load_contract()
    mutated = copy.deepcopy(contract)
    mutated[section][key] = value
    assert audit_contract(mutated, check_live_files=False)


@pytest.mark.parametrize(
    ("key", "value"),
    [
        ("nu", 0.02),
        ("forcing_mode", "free_residual_defined_force"),
        ("reference_energy_abs_tolerance", 0.1),
        ("momentum_max", 0.01),
        ("momentum_L2", 0.01),
        ("divergence_max", 0.001),
        ("divergence_L2", 0.001),
        ("residual_defined_free_force_forbidden", False),
        ("candidate_collapse_forbidden", False),
        ("post_hoc_threshold_relaxation_forbidden", False),
    ],
)
def test_cr001_drift_fails_closed(key, value):
    contract = load_contract()
    mutated = copy.deepcopy(contract)
    mutated["canonical_cr001"][key] = value
    assert audit_contract(mutated, check_live_files=False)


@pytest.mark.parametrize(
    ("key", "value"),
    [
        ("velocity_export_ready", False),
        ("visualization_ready", True),
        ("visual_correspondence_verified", True),
        ("pde_validated", True),
        ("paper_exact", True),
        ("openai_field_identified", True),
    ],
)
def test_eq45_state_laundering_fails_closed(key, value):
    contract = load_contract()
    mutated = copy.deepcopy(contract)
    mutated["canonical_eq45_state_must_remain_independent"][key] = value
    assert audit_contract(mutated, check_live_files=False)


def test_remaining_background_input_field_set_is_locked():
    contract = load_contract()
    mutated = copy.deepcopy(contract)
    mutated["remaining_caller_inputs"]["background_phase_derivative_jet"]["fields"].pop()
    assert audit_contract(mutated, check_live_files=False)


def test_remaining_amplitude_input_field_set_is_locked():
    contract = load_contract()
    mutated = copy.deepcopy(contract)
    mutated["remaining_caller_inputs"]["amplitude_directional_jet"]["fields"] = ["t_m"]
    assert audit_contract(mutated, check_live_files=False)


def test_mechanics_only_amplitude_dependence_witness_is_nonvacuous():
    witness = mechanics_only_amplitude_dependence_witness()
    assert witness["outputs_differ"] is True
    assert witness["coefficient_a"] != witness["coefficient_b"]
    assert "not Kokuno/OpenAI candidate data" in witness["classification"]
