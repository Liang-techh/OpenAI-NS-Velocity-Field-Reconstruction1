import math

import pytest

from openai_ns_reconstruction.constrained_eq45_bipolar_compact_poloidal_joint_renorm import (
    audit_compact_poloidal_joint_renormalization,
)


def test_joint_renormalization_is_energy_exact_but_sign_asymmetric_at_profile_bound():
    report = audit_compact_poloidal_joint_renormalization(
        quadrature_orders=(48, 96), morphology_grid_size=25
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

    negative = report["negative_trial"]
    positive = report["positive_trial"]
    assert negative["coefficient_before_common_scale"] == -0.25
    assert positive["coefficient_before_common_scale"] == 0.25

    target = report["energy_model"]["target_energy"]
    assert math.isclose(negative["jointly_renormalized_reference_energy"], target, rel_tol=0.0, abs_tol=2e-12)
    assert math.isclose(positive["jointly_renormalized_reference_energy"], target, rel_tol=0.0, abs_tol=2e-12)
    assert negative["common_velocity_scale"] > 1.0
    assert positive["common_velocity_scale"] < 1.0

    # The current field already has a stored profile coefficient at the inherited
    # limit, so the scale-up sign fails the simplest common coefficient scaling,
    # while the scale-down sign stays inside that same unchanged limit.
    assert negative["profile_bound_screen"]["simple_common_scale_representation_preflight_passed"] is False
    assert positive["profile_bound_screen"]["simple_common_scale_representation_preflight_passed"] is True
    assert negative["profile_bound_screen"]["maximum_abs_stored_profile_coefficient_after_scale"] > 4.0
    assert positive["profile_bound_screen"]["maximum_abs_stored_profile_coefficient_after_scale"] < 4.0

    assert negative["core_sign_screen"]["all_times_satisfied"] is True
    assert positive["core_sign_screen"]["all_times_satisfied"] is True
    assert negative["validation_time_energy"]["all_times_satisfied"] is True
    assert positive["validation_time_energy"]["all_times_satisfied"] is True

    # This is a bulk-shape knob, not the cap/tip mode: the previously observed
    # finite-amplitude center-compact trial must not be laundered into a tip claim.
    for row in (negative, positive):
        delta = row["morphology_delta_from_baseline"]
        assert abs(delta["axial_rms_over_Zp_delta"]) > 0.0
        assert delta["axial_q90_over_Zp_delta"] == pytest.approx(0.0, abs=1e-15)
        assert delta["axial_q99_over_Zp_delta"] == pytest.approx(0.0, abs=1e-15)


def test_joint_renormalization_rejects_invalid_trial_magnitude():
    with pytest.raises(ValueError, match="trial_magnitude"):
        audit_compact_poloidal_joint_renormalization(trial_magnitude=0.0)
    with pytest.raises(ValueError, match="trial_magnitude"):
        audit_compact_poloidal_joint_renormalization(trial_magnitude=float("nan"))
