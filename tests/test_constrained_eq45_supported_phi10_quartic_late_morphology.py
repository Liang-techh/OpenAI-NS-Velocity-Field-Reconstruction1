import numpy as np

from openai_ns_reconstruction.constrained_eq45_supported_phi10_quartic_late_morphology import (
    DEFAULT_RESOLUTIONS,
    DEFAULT_TIME,
    audit_quartic_late_offkeyframe_morphology,
)


def test_quartic_late_offkeyframe_morphology_is_nearly_static_and_resolved():
    report = audit_quartic_late_offkeyframe_morphology()

    assert report["schema"] == "eq45_supported_phi10_quartic_late_offkeyframe_morphology_v1"
    assert report["resolutions"] == list(DEFAULT_RESOLUTIONS)
    assert report["time"] == DEFAULT_TIME
    np.testing.assert_allclose(report["tau"], 0.75, rtol=0.0, atol=1.0e-15)
    np.testing.assert_allclose(report["coefficients"]["static"], -0.3, rtol=0.0, atol=1.0e-15)
    np.testing.assert_allclose(report["coefficients"]["cubic"], -0.321875, rtol=0.0, atol=1.0e-15)
    np.testing.assert_allclose(
        report["coefficients"]["quartic"],
        -0.2986481741573034,
        rtol=0.0,
        atol=1.0e-15,
    )
    np.testing.assert_allclose(
        report["quartic_nullspace_coefficient"],
        -0.2831460674157303,
        rtol=0.0,
        atol=1.0e-15,
    )

    fine = report["rows"][-1]
    summary = report["finest_summary"]

    # The derivative-balanced quartic coefficient excursion is only about 6.18%
    # of the cubic late-half excursion.  The actual smooth 3-D vorticity geometry
    # follows the same localization without assuming an exactly linear response.
    for key in (
        "vorticity2_weighted_radial_rms",
        "vorticity2_weighted_axial_rms",
        "vorticity2_weighted_aspect",
    ):
        excursion = summary[f"quartic_cubic_{key}_excursion_fraction"]
        assert 0.05 < excursion < 0.08

    # On 81^3 the quartic field is substantially closer to the static supported
    # morphology than cubic: all smooth geometry changes are below 0.02%.
    assert abs(summary["quartic_baseline_vorticity2_weighted_radial_rms_ratio"] - 1.0) < 2.0e-4
    assert abs(summary["quartic_baseline_vorticity2_weighted_axial_rms_ratio"] - 1.0) < 2.0e-4
    assert abs(summary["quartic_baseline_vorticity2_weighted_aspect_ratio"] - 1.0) < 2.0e-4
    assert abs(summary["quartic_baseline_radial_q99_ratio"] - 1.0) < 1.0e-12
    assert abs(summary["quartic_baseline_axial_q99_ratio"] - 1.0) < 1.0e-12

    # The late quartic off-keyframe does not create a support-collar isosurface.
    # Its already-tiny collar-vorticity budget stays within 0.5% of static.
    assert fine["quartic"]["superlevel_collar_voxel_fraction"] == 0.0
    assert fine["quartic"]["whole_grid_collar_vorticity2_fraction"] < 3.0e-5
    assert abs(summary["quartic_baseline_whole_grid_collar_vorticity2_fraction_ratio"] - 1.0) < 0.005
    assert 0.05 < summary["quartic_cubic_whole_grid_collar_vorticity2_fraction_excursion_fraction"] < 0.08

    # Smooth metrics remain stable from 65^3 -> 81^3.  Radial q99 receives a
    # wider 2% gate because thresholded Cartesian quantiles move in grid-sized
    # jumps.  Axial q99 is deliberately auxiliary: it moves one z-plane while
    # the smooth axial RMS stays below 1% change.
    assert summary["quartic_radial_q99_mid_to_fine_relative_change"] < 0.02
    assert summary["quartic_vorticity2_weighted_radial_rms_mid_to_fine_relative_change"] < 0.01
    assert summary["quartic_vorticity2_weighted_axial_rms_mid_to_fine_relative_change"] < 0.01
    assert summary["quartic_vorticity2_weighted_aspect_mid_to_fine_relative_change"] < 0.01
    assert summary["quartic_whole_grid_collar_vorticity2_fraction_mid_to_fine_relative_change"] < 0.05

    truth = report["truth_boundary"]
    assert truth["visualization_candidate_only"] is True
    for key in (
        "canonical_velocity_changed",
        "diagnostic_velocity_changed",
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
