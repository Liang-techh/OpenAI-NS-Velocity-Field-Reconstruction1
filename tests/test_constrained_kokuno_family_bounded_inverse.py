from __future__ import annotations

import numpy as np

from openai_ns_reconstruction.kokuno_family_bounded_inverse import (
    evaluate_family_bounded_inverse,
    generate_actual_core_report,
)
from openai_ns_reconstruction.kokuno_missing_covariance_column_target import (
    build_missing_covariance_column_target,
)


def _receipt(target=(0.5, 0.5)):
    radii = np.linspace(0.1, 0.5, 5)
    target_stress = np.tile(np.asarray(target, dtype=float), (len(radii), 1))
    current = np.tile(np.array([1.0, 0.0]), (len(radii), 1))
    return build_missing_covariance_column_target(
        radii,
        target_stress,
        current,
        target_floor_fraction=1.0e-12,
        response_floor_fraction=1.0e-12,
    )


def test_single_additional_column_matches_prior_budget_semantics():
    receipt = _receipt()
    additional = np.tile(np.array([[0.0, 1.0]]), (len(receipt.radii), 1, 1))
    report = evaluate_family_bounded_inverse(
        receipt,
        additional,
        coefficient_labels=["transverse"],
        current_coefficient_budget=0.6,
        additional_family_l1_budget=0.6,
        budget_unit_matches_jacobian=True,
        coefficient_unit="test_amplitude",
    )
    assert report["all_required_nodes_rank2"] is True
    assert report["algebraic_stress_fit_within_tolerance"] is True
    assert report["current_update_within_budget_on_all_required"] is True
    assert report["additional_family_l1_within_budget_on_all_required"] is True
    assert report["family_bounded_inverse_preflight_passed"] is True
    assert report["max_additional_coefficient_l1_on_required"] == 0.5


def test_duplicate_family_column_is_rejected_as_rank_deficient():
    receipt = _receipt()
    additional = np.tile(np.array([[2.0, 0.0]]), (len(receipt.radii), 1, 1))
    report = evaluate_family_bounded_inverse(
        receipt,
        additional,
        coefficient_labels=["duplicate"],
        current_coefficient_budget=2.0,
        additional_family_l1_budget=2.0,
        budget_unit_matches_jacobian=True,
        coefficient_unit="test_amplitude",
    )
    assert report["rank2_required_nodes"] == 0
    assert report["family_bounded_inverse_preflight_passed"] is False


def test_l1_guard_blocks_budget_laundering_by_duplicate_label_split():
    receipt = _receipt(target=(0.5, 1.0))
    # Minimum-l2 least squares splits the transverse coefficient into 0.5 + 0.5.
    # Each label is individually below 0.6, but the aggregate correction is 1.0.
    additional = np.tile(
        np.array([[0.0, 1.0], [0.0, 1.0]]),
        (len(receipt.radii), 1, 1),
    )
    report = evaluate_family_bounded_inverse(
        receipt,
        additional,
        coefficient_labels=["b0", "b1"],
        current_coefficient_budget=0.6,
        additional_family_l1_budget=0.6,
        budget_unit_matches_jacobian=True,
        coefficient_unit="test_amplitude",
    )
    assert report["all_required_nodes_rank2"] is True
    assert report["max_additional_coefficient_linf_on_required"] < 0.6
    assert report["max_additional_coefficient_l1_on_required"] > 0.99
    assert report["additional_family_l1_within_budget_on_all_required"] is False
    assert report["family_bounded_inverse_preflight_passed"] is False


def test_unit_mismatch_keeps_an_algebraically_good_family_fail_closed():
    receipt = _receipt()
    additional = np.tile(np.array([[0.0, 1.0]]), (len(receipt.radii), 1, 1))
    report = evaluate_family_bounded_inverse(
        receipt,
        additional,
        coefficient_labels=["transverse"],
        current_coefficient_budget=0.6,
        additional_family_l1_budget=0.6,
        budget_unit_matches_jacobian=False,
        coefficient_unit="unknown_source_units",
    )
    assert report["all_required_nodes_rank2"] is True
    assert report["algebraic_stress_fit_within_tolerance"] is True
    assert report["budget_unit_matches_jacobian"] is False
    assert report["family_bounded_inverse_preflight_passed"] is False


def test_real_existing_family_remains_rank_degenerate_and_cycle_blocked(tmp_path):
    report = generate_actual_core_report(
        output=tmp_path / "family_inverse.json",
        target_output=tmp_path / "target.json",
    )
    calibration = report["existing_one_label_family_calibration"]
    inverse = calibration["family_inverse"]
    assert report["inputs"]["surrogate_defect_used"] is False
    assert calibration["relative_error_to_existing_current_response"] < 2.0e-12
    assert inverse["rank2_required_nodes"] == 0
    assert inverse["family_bounded_inverse_preflight_passed"] is False
    assert report["routing"]["actual_source_multilabel_family_consumed"] is False
    assert report["routing"]["source_coefficient_unit_mapping_available"] is False
    assert report["routing"]["finite_correction_cycle_rerun_allowed"] is False
    assert report["truth_boundary"]["finite_correction_cycle_run"] is False
    assert report["truth_boundary"]["residual_reduction_claimed"] is False
    assert report["truth_boundary"]["pde_validated"] is False
