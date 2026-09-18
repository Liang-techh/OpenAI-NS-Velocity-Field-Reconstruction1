import numpy as np

from openai_ns_reconstruction.kokuno_missing_covariance_column_target import (
    build_missing_covariance_column_target,
)
from openai_ns_reconstruction.kokuno_second_column_bounded_inverse import (
    current_signed_coefficient_budget,
    evaluate_bounded_second_column,
    required_second_column_envelope,
)


def _synthetic_receipt():
    radii = np.linspace(0.08, 0.36, 9)
    scale = 1.0 + 0.1 * np.arange(9, dtype=float)
    current = np.column_stack((scale, np.zeros(9)))
    alpha = np.linspace(0.01, 0.03, 9)
    missing_y = np.linspace(0.15, 0.23, 9)
    target = alpha[:, None] * current + np.column_stack((np.zeros(9), missing_y))
    return build_missing_covariance_column_target(radii, target, current)


def test_required_envelope_is_exact_second_coefficient_floor():
    receipt = _synthetic_receipt()
    budget = current_signed_coefficient_budget(receipt)
    np.testing.assert_allclose(budget, 0.03, rtol=0.0, atol=1.0e-15)

    envelope = required_second_column_envelope(receipt, coefficient_budget=budget)
    needed = receipt.active_target_mask
    expected = np.linalg.norm(receipt.missing_response[needed], axis=1) / budget
    np.testing.assert_allclose(
        envelope["required_transverse_response"], expected, rtol=1.0e-14, atol=1.0e-14
    )
    expected_det = np.linalg.norm(receipt.current_response[needed], axis=1) * expected
    np.testing.assert_allclose(
        envelope["required_determinant_abs_min"], expected_det, rtol=1.0e-14, atol=1.0e-14
    )


def test_raw_missing_direction_spans_target_but_violates_small_signed_budget():
    receipt = _synthetic_receipt()
    budget = current_signed_coefficient_budget(receipt)
    audit = evaluate_bounded_second_column(
        receipt, receipt.missing_response, coefficient_budget=budget
    )
    assert audit["all_required_nodes_rank2"]
    assert audit["algebraic_stress_fit_within_tolerance"]
    np.testing.assert_allclose(
        audit["max_abs_second_signed_coefficient_on_required"], 1.0, rtol=1.0e-13, atol=1.0e-13
    )
    assert not audit["second_signed_update_within_budget_on_all_required"]
    assert not audit["bounded_inverse_preflight_passed"]


def test_scaled_missing_direction_hits_budget_and_passes_algebraic_preflight():
    receipt = _synthetic_receipt()
    budget = current_signed_coefficient_budget(receipt)
    second = receipt.missing_response / budget
    audit = evaluate_bounded_second_column(receipt, second, coefficient_budget=budget)
    assert audit["all_required_nodes_rank2"]
    assert audit["both_signed_updates_within_budget_on_all_required"]
    assert audit["algebraic_stress_fit_within_tolerance"]
    assert audit["bounded_inverse_preflight_passed"]
    np.testing.assert_allclose(
        audit["max_abs_second_signed_coefficient_on_required"], budget, rtol=1.0e-13, atol=1.0e-13
    )


def test_duplicate_column_is_rejected_even_with_large_budget():
    receipt = _synthetic_receipt()
    audit = evaluate_bounded_second_column(
        receipt, receipt.current_response, coefficient_budget=10.0
    )
    assert audit["rank2_required_nodes"] == 0
    assert audit["rank_deficient_required_nodes"] == audit["nodes_requiring_second_direction"]
    assert not audit["bounded_inverse_preflight_passed"]


def test_bad_budgets_and_shapes_fail_closed():
    receipt = _synthetic_receipt()
    with np.testing.assert_raises(ValueError):
        required_second_column_envelope(receipt, coefficient_budget=0.0)
    with np.testing.assert_raises(ValueError):
        evaluate_bounded_second_column(
            receipt, np.ones((8, 2)), coefficient_budget=0.03
        )
    with np.testing.assert_raises(ValueError):
        evaluate_bounded_second_column(
            receipt,
            receipt.missing_response,
            coefficient_budget=0.03,
            relative_rank_tolerance=0.1,
        )
