import math

import pytest

from openai_ns_reconstruction.constrained_eq45_bipolar_axial_cap_joint_renorm import (
    audit_axial_cap_joint_renormalization,
)


def test_axial_cap_joint_renorm_restores_energy_and_preserves_tip_capacity():
    report = audit_axial_cap_joint_renormalization(
        quadrature_orders=(48, 96), morphology_grid_size=33
    )

    assert report["basis_growth"] == {
        "new_basis_shapes_added": 0,
        "new_temporal_degrees_added": 0,
        "new_scalar_fit_parameters_added": 0,
        "coefficient_selected": False,
    }
    assert report["truth_boundary"]["candidate_artifact_changed"] is False
    assert report["truth_boundary"]["held_out_pde_residual_evaluated"] is False
    assert report["truth_boundary"]["visual_correspondence_verified"] is False
    assert report["truth_boundary"]["pde_validated"] is False

    target = report["energy_model"]["target_energy"]
    negative = report["negative_trial"]
    positive = report["positive_trial"]
    assert negative["coefficient_before_common_scale"] == -0.25
    assert positive["coefficient_before_common_scale"] == 0.25

    for row in (negative, positive):
        assert math.isclose(row["jointly_renormalized_reference_energy"], target, rel_tol=0.0, abs_tol=2e-12)
        assert 0.0 < row["common_velocity_scale"] < 1.0
        assert row["profile_bound_screen"]["simple_common_scale_representation_preflight_passed"] is True
        assert row["profile_bound_screen"]["maximum_abs_stored_profile_coefficient_after_scale"] < 4.0
        assert abs(row["profile_bound_screen"]["effective_axial_cap_coefficient_after_common_scale"]) < 0.25
        assert row["core_sign_screen"]["all_times_satisfied"] is True
        assert row["validation_time_energy"]["all_times_satisfied"] is True
        assert row["gross_tip_reach_extended_on_this_grid"] is True
        delta = row["morphology_delta_from_baseline"]
        assert delta["axial_rms_over_Zp_delta"] > 0.0
        assert delta["axial_q90_over_Zp_delta"] > 0.0
        assert delta["axial_q99_over_Zp_delta"] > 0.0

    assert report["both_signs_simple_representation_preflight_passed"] is True
    assert report["both_signs_gross_tip_reach_extended_on_this_grid"] is True


def test_axial_cap_joint_renorm_scales_match_refined_energy_quadratic():
    report = audit_axial_cap_joint_renormalization(
        quadrature_orders=(96, 192), morphology_grid_size=25
    )
    negative = report["negative_trial"]
    positive = report["positive_trial"]

    assert negative["raw_reference_energy"] == pytest.approx(1.013801194448826, rel=2e-10, abs=2e-12)
    assert positive["raw_reference_energy"] == pytest.approx(1.013666732563905, rel=2e-10, abs=2e-12)
    assert negative["common_velocity_scale"] == pytest.approx(0.993170018451513, rel=2e-10, abs=2e-12)
    assert positive["common_velocity_scale"] == pytest.approx(0.993235887775290, rel=2e-10, abs=2e-12)
    assert report["energy_model"]["maximum_relative_quadratic_coefficient_change"] < 2e-4


def test_axial_cap_joint_renorm_rejects_invalid_inputs():
    with pytest.raises(ValueError, match="trial_magnitude"):
        audit_axial_cap_joint_renormalization(trial_magnitude=0.0)
    with pytest.raises(ValueError, match="trial_magnitude"):
        audit_axial_cap_joint_renormalization(trial_magnitude=float("nan"))
    with pytest.raises(ValueError, match="morphology_grid_size"):
        audit_axial_cap_joint_renormalization(morphology_grid_size=24)
