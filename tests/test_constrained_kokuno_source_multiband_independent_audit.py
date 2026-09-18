from __future__ import annotations

import pytest

from openai_ns_reconstruction.kokuno_source_multiband_independent_audit import (
    CASES,
    FROZEN_GUARDS,
    build_report,
)


@pytest.fixture(scope="module")
def report():
    return build_report()


def test_fresh_multiband_cartesian_preflight_passes_frozen_guards(report):
    assert report["local_structural_preflight_passed"] is True
    assert all(report["checks"].values())
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


def test_resolution_and_mutation_controls_are_discriminating(report):
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
