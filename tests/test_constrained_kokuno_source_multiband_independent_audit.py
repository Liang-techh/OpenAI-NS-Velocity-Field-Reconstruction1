from __future__ import annotations

import pytest

from openai_ns_reconstruction.kokuno_source_multiband_independent_audit import (
    BASE_HEAD,
    BASE_PR,
    CASES,
    FD4_STEPS,
    FROZEN_GUARDS,
    SEED,
    build_report,
)


@pytest.fixture(scope="module")
def report():
    return build_report()


def test_frozen_451_contract_is_unchanged(report):
    assert SEED == 9173171
    assert FD4_STEPS == (0.032, 0.016, 0.008)
    assert FROZEN_GUARDS == {
        "source_schedule_relative_error_max": 5.0e-14,
        "finest_cartesian_curl_relative_rms_max": 2.0e-5,
        "curl_refinement_ratio_min": 4.0,
        "finest_normalized_divergence_rms_max": 2.0e-5,
        "finest_normalized_divergence_point_max": 8.0e-5,
        "divergence_refinement_ratio_min": 3.0,
        "shared_first_band_schedule_mutation_relative_rms_min": 0.10,
        "label_schedule_swap_mutation_relative_rms_min": 0.10,
    }
    assert BASE_PR == 461
    assert BASE_HEAD == "aafb67af44c40794701ed0d66dbf0a08c1e45fb3"
    assert report["dependency"]["previous_agent4_pr"] == 451


def test_independent_numerical_and_permutation_guards_all_pass(report):
    assert report["source_schedule_relative_error_max"] <= FROZEN_GUARDS[
        "source_schedule_relative_error_max"
    ]
    assert report["curl_resolution_ladder"][-1]["curl_relative_rms"] <= FROZEN_GUARDS[
        "finest_cartesian_curl_relative_rms_max"
    ]
    assert report["divergence_resolution_ladder"][-1][
        "normalized_divergence_rms"
    ] <= FROZEN_GUARDS["finest_normalized_divergence_rms_max"]
    assert report["divergence_resolution_ladder"][-1][
        "normalized_divergence_point_max"
    ] <= FROZEN_GUARDS["finest_normalized_divergence_point_max"]
    assert all(report["checks"].values())
    assert report["mutation"]["consistent_public_label_permutation_failures"] == 0
    assert report["mutation"]["consistent_public_label_permutation_relative_max_when_evaluable"] == 0.0
    assert report["mutation"]["consistent_public_label_permutation_bitwise_all"] is True
    assert report["local_structural_preflight_passed"] is True


def test_resolution_and_mutation_controls_remain_discriminating(report):
    assert min(report["curl_refinement_ratios"]) >= FROZEN_GUARDS[
        "curl_refinement_ratio_min"
    ]
    assert min(report["divergence_refinement_ratios"]) >= FROZEN_GUARDS[
        "divergence_refinement_ratio_min"
    ]
    assert report["mutation"]["shared_first_band_schedule_relative_rms"] >= FROZEN_GUARDS[
        "shared_first_band_schedule_mutation_relative_rms_min"
    ]
    assert report["mutation"]["label_schedule_swap_relative_rms"] >= FROZEN_GUARDS[
        "label_schedule_swap_mutation_relative_rms_min"
    ]
    transition = next(case for case in report["parameter_cases"] if case["name"] == "k_transition")
    assert [item["k"] for item in transition["schedules"]] == [2, 3, 3]
    assert tuple(CASES[-1]["ells"]) == (221, 223, 225)


def test_truth_boundary_remains_fail_closed(report):
    truth = report["truth_boundary"]
    assert all(value is False for value in truth.values())
    assert report["formal_project_gates"]["normalized_momentum_max_L2"] == 1.0e-3
    assert report["formal_project_gates"]["divergence_max_L2"] == 1.0e-5
    assert report["formal_project_gates"]["formal_full_domain_pde_gate_assessed"] is False
    assert report["independent_operator"]["training_tensor_or_loss_read"] is False
    assert report["independent_operator"]["agent2_complete_curl_helper_used_by_reference"] is False
    assert report["independent_operator"]["free_forcing_used"] is False
