from __future__ import annotations

import math

import pytest

from openai_ns_reconstruction.kokuno_autonomous_mean_factor import (
    KokunoAutonomousMeanFactor,
    autonomous_dyadic_masks,
    autonomous_missing_weight,
    autonomous_profile,
)


def test_autonomous_profile_is_peak_normalized_and_compact() -> None:
    assert autonomous_profile(0.0) == pytest.approx(1.0)
    assert autonomous_profile(0.5) > 0.0
    assert autonomous_profile(-0.5) == pytest.approx(autonomous_profile(0.5))
    assert autonomous_profile(1.0) == 0.0
    assert autonomous_profile(-1.0) == 0.0
    assert autonomous_profile(2.0) == 0.0


def test_autonomous_masks_form_squared_partition() -> None:
    for q in (2.0**-7, 1.3 * 2.0**-5, 1.7 * 2.0**3):
        masks = autonomous_dyadic_masks(q)
        assert 1 <= len(masks) <= 2
        assert math.fsum(value * value for value in masks.values()) == pytest.approx(
            1.0, abs=2.0e-15
        )


def test_head_weight_changes_only_when_active_translate_enters_head() -> None:
    q = 1.3 * 2.0**-5
    masks = autonomous_dyadic_masks(q)
    assert tuple(sorted(masks)) == (4, 5)
    assert autonomous_missing_weight(4, q) == 0.0
    assert autonomous_missing_weight(5, q) == pytest.approx(
        masks[4] * masks[4], abs=2.0e-15
    )
    assert autonomous_missing_weight(6, q) == pytest.approx(1.0, abs=2.0e-15)


def test_receipt_keeps_formal_and_autonomous_factors_separate() -> None:
    report = KokunoAutonomousMeanFactor(
        prepared_n=5,
        band=5,
        coordinate_q=1.3,
    ).evaluate()

    assert report["physical_q"] == pytest.approx(0.040625)
    assert report["active_mask_indices"] == (4, 5)
    assert report["autonomous_missing_weight"] == pytest.approx(
        0.28410624176179694, abs=2.0e-15
    )
    assert report["squared_partition_closure_error"] <= 2.0e-15
    assert report["profile_is_repository_autonomous"] is True
    assert report["formal_theorem_machine_bump_identity_claimed"] is False
    assert report["formal_missing_weight_equality_claimed"] is False
    assert report["theorem_missing_weight_replaced"] is False
    assert report["theorem_missing_weight_materialized"] is False
    assert report["requested_stress_actual_state_values_materialized"] is False
    assert report["finite_head_mean_debt_materialized"] is False
    assert report["surrogate_defect_used"] is False
    assert report["real_candidate_defect_consumed"] is False
    assert report["signed_mean_inverse_input_ready"] is False
    assert report["finite_correction_cycle_rerun_allowed"] is False
    assert report["residual_reduction_claimed"] is False
    assert report["pde_validated"] is False


def test_invalid_scales_fail_closed() -> None:
    with pytest.raises(ValueError):
        autonomous_dyadic_masks(0.0)
    with pytest.raises(ValueError):
        autonomous_dyadic_masks(float("nan"))
    with pytest.raises(ValueError):
        autonomous_profile(float("inf"))
    with pytest.raises(ValueError):
        KokunoAutonomousMeanFactor(prepared_n=5, band=5, coordinate_q=-1.0).evaluate()
