import pytest

from openai_ns_reconstruction.constrained_eq45_supported_phi10_blend_axial_taper_morphology import (
    TASK_ID,
    audit_axial_taper_vorticity_morphology,
)


def _summary_by_q(report):
    return {
        float(row["axial_plateau_q"]): row
        for row in report["finest_summaries"]
    }


def test_axial_taper_changes_only_subthreshold_axial_collar_vorticity_at_early_time():
    report = audit_axial_taper_vorticity_morphology()
    summaries = _summary_by_q(report)

    assert report["task_id"] == TASK_ID
    assert report["resolutions"] == [49, 65, 81]
    assert report["axial_plateau_q_values"] == [0.64, 0.7225, 0.81]
    assert report["time"] == pytest.approx(0.3125)
    assert report["blend_weight"] == pytest.approx(0.5)
    assert report["early_delta"] == pytest.approx(-1.4)

    assert summaries[0.64]["axial_identity_half_height"] == pytest.approx(1.6)
    assert summaries[0.7225]["axial_identity_half_height"] == pytest.approx(1.7)
    assert summaries[0.81]["axial_identity_half_height"] == pytest.approx(1.8)

    # The 25%-of-core envelope is unchanged even at the widest axial plateau.
    # The smooth bulk geometry likewise moves by much less than 0.01%.
    extended = summaries[0.81]
    assert extended["baseline_radial_q99_ratio"] == pytest.approx(1.0, abs=1.0e-12)
    assert extended["baseline_axial_q99_ratio"] == pytest.approx(1.0, abs=1.0e-12)
    assert extended["baseline_aspect_q99_ratio"] == pytest.approx(1.0, abs=1.0e-12)
    assert abs(extended["baseline_vorticity2_weighted_radial_rms_ratio"] - 1.0) < 1.0e-4
    assert abs(extended["baseline_vorticity2_weighted_axial_rms_ratio"] - 1.0) < 1.0e-4
    assert abs(extended["baseline_vorticity2_weighted_aspect_ratio"] - 1.0) < 1.0e-4
    assert extended["superlevel_axial_collar_voxel_fraction"] == 0.0
    assert extended["superlevel_whole_collar_voxel_fraction"] == 0.0

    # The only material response is below the core threshold in the fixed
    # |z|>1.6 collar.  It is not used to select q because this tiny outer signal
    # is much less resolution-stable than the smooth core metrics.
    middle = summaries[0.7225]
    assert middle["baseline_axial_outer_vorticity2_fraction_ratio"] > 1.10
    assert extended["baseline_axial_outer_vorticity2_fraction_ratio"] > 1.30
    assert extended["baseline_whole_grid_collar_vorticity2_fraction_ratio"] < 1.0001
    assert extended["axial_outer_vorticity2_fraction_mid_to_fine_relative_change"] > 0.05

    for row in summaries.values():
        assert row["radial_q99_mid_to_fine_relative_change"] < 0.01
        assert row["axial_q99_mid_to_fine_relative_change"] == pytest.approx(0.0, abs=1.0e-15)
        assert row["vorticity2_weighted_radial_rms_mid_to_fine_relative_change"] < 0.01
        assert row["vorticity2_weighted_axial_rms_mid_to_fine_relative_change"] < 0.01
        assert row["vorticity2_weighted_aspect_mid_to_fine_relative_change"] < 0.01

    truth = report["truth_boundary"]
    assert truth["visualization_candidate_only"] is True
    for key in (
        "velocity_changed_by_diagnostic",
        "axial_taper_value_selected",
        "canonical_velocity_changed",
        "force_or_pressure_fitted",
        "pde_objective_used_to_choose_axial_taper",
        "public_image_fitted",
        "physical_support_validated",
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        assert truth[key] is False


def test_axial_taper_morphology_audit_rejects_non_governed_grids():
    with pytest.raises(ValueError, match="begin at the governed"):
        audit_axial_taper_vorticity_morphology(
            resolutions=(9, 11, 13),
            axial_plateau_q_values=(0.70, 0.81),
        )
    with pytest.raises(ValueError, match="exceeds the autonomous"):
        audit_axial_taper_vorticity_morphology(
            resolutions=(9, 11, 13),
            axial_plateau_q_values=(0.64, 0.82),
        )
