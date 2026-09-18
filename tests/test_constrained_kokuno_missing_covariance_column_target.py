import numpy as np

from openai_ns_reconstruction.kokuno_missing_covariance_column_target import (
    build_missing_covariance_column_target,
    evaluate_second_column_candidate,
)


def _synthetic_contract():
    radii = np.linspace(0.08, 0.36,  nine := 9)
    current = np.column_stack(
        (
            1.0 + 0.1 * np.arange(nine, dtype=float),
            0.25 + 0.02 * np.arange(nine, dtype=float),
        )
    )
    perpendicular = np.column_stack((-current[:, 1], current[:, 0]))
    perpendicular /= np.linalg.norm(perpendicular, axis=1)[:, None]
    missing = (0.20 + 0.01 * np.arange(nine, dtype=float))[:, None] * perpendicular
    target = 0.45 * current + missing
    return build_missing_covariance_column_target(radii, target, current)


def test_missing_target_is_pointwise_orthogonal_and_exact_abstract_second_column():
    receipt = _synthetic_contract()
    active = receipt.safe_active_mask
    dot = np.sum(
        receipt.current_response[active] * receipt.missing_response[active], axis=1
    )
    np.testing.assert_allclose(dot, 0.0, rtol=0.0, atol=2.0e-16)
    np.testing.assert_allclose(
        receipt.current_projection + receipt.missing_response,
        receipt.target_stress,
        rtol=1.0e-14,
        atol=1.0e-14,
    )

    audit = evaluate_second_column_candidate(receipt, receipt.missing_response)
    assert audit["all_active_nodes_rank2"]
    assert audit["rank_deficient_active_nodes"] == 0
    assert audit["two_column_relative_residual_rms"] < 1.0e-13
    assert audit["incremental_relative_residual_reduction"] > 0.999999999999


def test_duplicate_current_column_stays_rank_one_and_cannot_remove_missing_target():
    receipt = _synthetic_contract()
    audit = evaluate_second_column_candidate(receipt, receipt.current_response)
    assert audit["rank2_active_nodes"] == 0
    assert audit["rank_deficient_active_nodes"] == audit["active_target_nodes"]
    np.testing.assert_allclose(
        audit["two_column_relative_residual_rms"],
        audit["rank1_relative_residual_rms"],
        rtol=1.0e-13,
        atol=1.0e-13,
    )
    assert audit["incremental_relative_residual_reduction"] < 1.0e-12


def test_small_current_response_is_left_unresolved_instead_of_divided_through():
    radii = np.linspace(0.08, 0.36, 9)
    target = np.column_stack((np.linspace(1.0, 1.8, 9), np.linspace(0.4, 0.8, 9)))
    current = np.ones((9, 2), dtype=float)
    current[4] = 1.0e-9
    receipt = build_missing_covariance_column_target(
        radii,
        target,
        current,
        target_floor_fraction=1.0e-3,
        response_floor_fraction=1.0e-3,
    )
    assert receipt.active_target_mask[4]
    assert not receipt.safe_response_mask[4]
    assert receipt.unresolved_active_mask[4]
    np.testing.assert_allclose(receipt.current_projection[4], 0.0)
    np.testing.assert_allclose(receipt.missing_response[4], target[4])
    assert receipt.current_coefficient[4] == 0.0


def test_contract_rejects_zero_or_malformed_response_inputs():
    radii = np.linspace(0.08, 0.36, 9)
    target = np.ones((9, 2), dtype=float)
    with np.testing.assert_raises(ValueError):
        build_missing_covariance_column_target(radii, target, np.zeros((9, 2)))
    with np.testing.assert_raises(ValueError):
        build_missing_covariance_column_target(radii, target[:, 0], np.ones((9, 2)))

    receipt = _synthetic_contract()
    with np.testing.assert_raises(ValueError):
        evaluate_second_column_candidate(receipt, np.ones((8, 2)))
    with np.testing.assert_raises(ValueError):
        evaluate_second_column_candidate(
            receipt,
            receipt.missing_response,
            relative_rank_tolerance=0.1,
        )
