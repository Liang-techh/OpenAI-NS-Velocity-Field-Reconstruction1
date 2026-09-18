import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_source_signed_covariance_bridge import (
    TARGET_UNITS,
    KokunoSourceSignedCovarianceBridge,
    structural_receipt,
)


def test_epsilon_scaled_reference_inverse_reconstructs_physical_target():
    bridge = KokunoSourceSignedCovarianceBridge()
    eps = np.array([4.0e-2, 1.0e-2, 2.5e-3])
    out = bridge.screen_target(
        A_c=2.0,
        u_star=3.0,
        h_plus=1.0,
        h_minus=1.0,
        epsilon=eps,
        target_N=-8.0e-2,
        target_K=2.4e-2,
        target_kind="structural_calibration",
        target_provenance="unit_test_predeclared_structural_target",
        direction_gap_eta=0.2,
    )

    expected_target = np.broadcast_to(np.array([-8.0e-2, 2.4e-2]), (3, 2))
    expected_y = np.array([[0.4, 0.6], [1.6, 2.4], [6.4, 9.6]])
    np.testing.assert_allclose(out["physical_reconstructed_target"], expected_target, rtol=0.0, atol=2e-15)
    np.testing.assert_allclose(out["squared_amplitudes"], expected_y, rtol=0.0, atol=3e-15)
    assert out["max_abs_physical_reconstruction_error"] <= 2e-15
    assert out["min_cone_margin_fraction"] == pytest.approx(0.8)
    assert out["max_reference_column_condition_number"] == pytest.approx(1.5)
    assert out["reference_covariance_rank_two"] is True
    assert out["reference_inverse_preflight_passed"] is True


def test_missing_epsilon_division_would_be_wrong_physical_units():
    bridge = KokunoSourceSignedCovarianceBridge()
    out = bridge.screen_target(
        A_c=2.0,
        u_star=3.0,
        h_plus=1.0,
        h_minus=1.0,
        epsilon=0.04,
        target_N=-0.08,
        target_K=0.024,
        target_kind="structural_calibration",
        target_provenance="epsilon_unit_regression",
    )
    np.testing.assert_allclose(out["source_reference_target"], [-2.0, 0.6], rtol=0.0, atol=2e-15)
    np.testing.assert_allclose(out["physical_reconstructed_target"], [-0.08, 0.024], rtol=0.0, atol=2e-15)
    assert out["epsilon_scaling_identity"] == "physical_target = epsilon * H_ref * y"


def test_surrogate_and_ambiguous_target_metadata_fail_closed():
    bridge = KokunoSourceSignedCovarianceBridge()
    base = dict(
        A_c=2.0,
        u_star=3.0,
        h_plus=1.0,
        h_minus=1.0,
        epsilon=0.04,
        target_N=-0.08,
        target_K=0.024,
        target_kind="real_candidate_phase_mean_defect",
        target_provenance="candidate_sha256:0123456789abcdef",
    )
    with pytest.raises(ValueError, match="surrogate defects"):
        bridge.screen_target(**base, surrogate_defect_used=True)
    with pytest.raises(ValueError, match="target_provenance"):
        bridge.screen_target(**{**base, "target_provenance": ""})
    with pytest.raises(ValueError, match="target_units"):
        bridge.screen_target(**base, target_units="source_reference_T0")
    with pytest.raises(ValueError, match="target_kind"):
        bridge.screen_target(**{**base, "target_kind": "validation_reused_as_training"})


def test_positive_cone_still_rejects_inadmissible_physical_target():
    bridge = KokunoSourceSignedCovarianceBridge()
    with pytest.raises(ValueError, match="positive two-sign source covariance cone"):
        bridge.screen_target(
            A_c=2.0,
            u_star=3.0,
            h_plus=1.0,
            h_minus=1.0,
            epsilon=0.04,
            target_N=-0.08,
            target_K=0.30,
            target_kind="structural_calibration",
            target_provenance="inadmissible_cone_regression",
        )


def test_real_target_label_does_not_promote_physical_correction_readiness():
    bridge = KokunoSourceSignedCovarianceBridge()
    out = bridge.screen_target(
        A_c=2.0,
        u_star=3.0,
        h_plus=1.0,
        h_minus=1.0,
        epsilon=0.04,
        target_N=-0.08,
        target_K=0.024,
        target_kind="real_candidate_phase_mean_defect",
        target_provenance="candidate_sha256:0123456789abcdef; independent_operator=true",
    )
    assert out["real_candidate_defect_consumed"] is True
    assert out["actual_source_h_sigma_pulse_integrals_bound"] is False
    assert out["actual_signed_auxiliary_rectangles_bound"] is False
    assert out["physical_complete_curl_signed_family_bound"] is False
    assert out["public_velocity_correction_materialized"] is False
    assert out["finite_correction_cycle_rerun_allowed"] is False
    assert out["heldout_ns_residual_assessed"] is False
    assert out["residual_reduction_claimed"] is False
    assert out["pde_validated"] is False


def test_structural_receipt_is_explicitly_not_real_defect_evidence():
    receipt = structural_receipt()
    assert receipt["epsilon"] == pytest.approx([0.04, 0.01, 0.0025])
    assert receipt["min_cone_margin_fraction"] == pytest.approx(0.8)
    assert receipt["max_reference_column_condition_number"] == pytest.approx(1.5)
    assert receipt["max_abs_physical_reconstruction_error"] <= 2e-15
    assert receipt["reference_covariance_rank_two"] is True
    assert receipt["reference_inverse_preflight_passed"] is True
    assert receipt["real_candidate_defect_consumed"] is False
    assert receipt["physical_complete_curl_signed_family_bound"] is False
    assert receipt["finite_correction_cycle_rerun_allowed"] is False
    assert receipt["residual_reduction_claimed"] is False
    assert receipt["pde_validated"] is False


def test_expected_target_units_are_frozen():
    assert TARGET_UNITS == "physical_phase_mean_covariance_correction"
