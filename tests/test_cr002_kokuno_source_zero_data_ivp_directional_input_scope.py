from __future__ import annotations

import copy

import pytest

from openai_ns_reconstruction.audit_kokuno_source_zero_data_ivp_directional_input_scope import (
    assert_clean,
    audit_contract,
    load_contract,
    mechanics_only_partial_input_retirement_witness,
)


def _mutated(contract: dict, path: tuple[str, ...], value):
    candidate = copy.deepcopy(contract)
    node = candidate
    for key in path[:-1]:
        node = node[key]
    node[path[-1]] = value
    return candidate


def test_registered_contract_is_clean_without_live_io():
    contract = load_contract()
    assert audit_contract(contract, check_live_files=False) == []


def test_registered_contract_is_clean_against_exact_live_parent_and_cr001():
    assert_clean()


def test_mechanics_witness_is_nonvacuous_and_autonomous_only():
    receipt = mechanics_only_partial_input_retirement_witness()
    assert receipt["caller_forcing_remains_output_determining"] is True
    assert receipt["same_point_state_does_not_determine_directional_jet"] is True
    assert (
        receipt["state_for_constant_forcing_1"]
        != receipt["state_for_constant_forcing_2"]
    )
    assert receipt["directional_jet_a"] != receipt["directional_jet_b"]
    assert receipt["classification"].startswith("autonomous mechanics only")


@pytest.mark.parametrize(
    ("path", "value"),
    [
        (("exact_parent", "pr"), 1025),
        (("exact_parent", "head"), "0" * 40),
        (("exact_parent", "module_git_blob"), "0" * 40),
        (
            ("parent_interface_truth", "source_zero_data_ivp_numerical_solver_materialized"),
            False,
        ),
        (("parent_interface_truth", "caller_supplies_mode_state_t_m"), True),
        (("parent_interface_truth", "caller_supplies_mode_forcing_f_m"), False),
        (("parent_interface_truth", "caller_supplies_corrected_background"), False),
        (
            ("parent_interface_truth", "source_amplitude_directional_jet_materialized"),
            True,
        ),
        (
            ("parent_interface_truth", "source_amplitude_mode_provider_complete"),
            True,
        ),
        (
            ("parent_interface_truth", "self_contained_velocity_xyzt_provider"),
            True,
        ),
        (
            ("parent_interface_truth", "source_exact_duhamel_frame_B_materialized"),
            True,
        ),
        (
            ("parent_interface_truth", "source_exact_duhamel_propagator_Vm_materialized"),
            True,
        ),
        (("parent_interface_truth", "complete_ns_residual_assessed"), True),
        (("parent_interface_truth", "pde_validated"), True),
        (("input_retirement", "retired_caller_inputs"), []),
        (
            ("input_retirement", "still_caller_supplied_output_determining_inputs"),
            ["corrected_background_jet"],
        ),
        (
            ("input_retirement", "still_unmaterialized_amplitude_directional_outputs"),
            [],
        ),
        (
            ("input_retirement", "retiring_t_m_does_not_retire_directional_jets"),
            False,
        ),
        (
            ("representation_boundary", "caller_free_t_m_is_not_caller_free_source_realization"),
            False,
        ),
        (
            (
                "representation_boundary",
                "source_ivp_provider_is_not_unified_cartesian_velocity_provider",
            ),
            False,
        ),
        (
            (
                "forbidden_promotions",
                "zero_data_ivp_materialized_implies_amplitude_directional_jet_materialized",
            ),
            True,
        ),
        (
            (
                "forbidden_promotions",
                "caller_no_longer_supplies_t_m_implies_no_output_determining_source_inputs_remain",
            ),
            True,
        ),
        (
            (
                "forbidden_promotions",
                "fixed_rk4_ladder_implies_exact_B_frame_or_Vm_propagator",
            ),
            True,
        ),
        (
            (
                "forbidden_promotions",
                "generated_t_m_implies_source_amplitude_mode_provider_complete",
            ),
            True,
        ),
        (
            (
                "forbidden_promotions",
                "generated_t_m_implies_self_contained_velocity_xyzt_provider",
            ),
            True,
        ),
        (
            ("forbidden_promotions", "source_ivp_mechanics_implies_pde_validated"),
            True,
        ),
        (
            ("forbidden_promotions", "source_ivp_mechanics_implies_paper_exact"),
            True,
        ),
        (
            (
                "forbidden_promotions",
                "source_ivp_mechanics_implies_openai_field_identified",
            ),
            True,
        ),
        (("canonical_cr001", "nu"), 0.02),
        (("canonical_cr001", "forcing_mode"), "free_residual"),
        (("canonical_cr001", "momentum_max"), 0.01),
        (("canonical_cr001", "divergence_L2"), 1e-4),
        (
            ("canonical_cr001", "residual_defined_free_force_forbidden"),
            False,
        ),
        (("canonical_cr001", "candidate_collapse_forbidden"), False),
        (
            ("canonical_cr001", "post_hoc_threshold_relaxation_forbidden"),
            False,
        ),
        (
            ("canonical_eq45_state_must_remain_independent", "velocity_export_ready"),
            False,
        ),
        (
            (
                "canonical_eq45_state_must_remain_independent",
                "visual_correspondence_verified",
            ),
            True,
        ),
        (
            ("canonical_eq45_state_must_remain_independent", "pde_validated"),
            True,
        ),
        (
            ("canonical_eq45_state_must_remain_independent", "paper_exact"),
            True,
        ),
    ],
)
def test_fail_closed_mutations_are_rejected(path, value):
    contract = load_contract()
    candidate = _mutated(contract, path, value)
    assert audit_contract(candidate, check_live_files=False)


def test_four_way_provenance_cannot_be_collapsed_or_erased():
    contract = load_contract()

    collapsed = copy.deepcopy(contract)
    collapsed["four_way_provenance"]["public_source_fact"] = []
    assert audit_contract(collapsed, check_live_files=False)

    relabeled = copy.deepcopy(contract)
    relabeled["four_way_provenance"]["source_fact"] = relabeled[
        "four_way_provenance"
    ].pop("public_source_fact")
    assert audit_contract(relabeled, check_live_files=False)
