import numpy as np

from openai_ns_reconstruction.constrained_eq45_supported_phi10_cubic_late_morphology import (
    DEFAULT_RESOLUTIONS,
    DEFAULT_TIME,
    audit_cubic_late_offkeyframe_morphology,
)


def test_cubic_late_offkeyframe_morphology_is_small_and_resolved():
    report = audit_cubic_late_offkeyframe_morphology()

    assert report["schema"] == "eq45_supported_phi10_cubic_late_offkeyframe_morphology_v1"
    assert report["resolutions"] == list(DEFAULT_RESOLUTIONS)
    assert report["time"] == DEFAULT_TIME
    np.testing.assert_allclose(report["tau"], 0.75, rtol=0.0, atol=1.0e-15)
    np.testing.assert_allclose(report["coefficients"]["static"], -0.3, rtol=0.0, atol=1.0e-15)
    np.testing.assert_allclose(report["coefficients"]["quadratic"], -0.16875, rtol=0.0, atol=1.0e-15)
    np.testing.assert_allclose(report["coefficients"]["cubic"], -0.321875, rtol=0.0, atol=1.0e-15)

    fine = report["rows"][-1]
    summary = report["finest_summary"]

    # The cubic coefficient excursion from static is one sixth the magnitude of
    # the quadratic late-half excursion (and has opposite sign).  The smooth
    # whole-domain morphology response follows that localization closely without
    # assuming exact linearity in the velocity-to-vorticity map.
    for key in (
        "vorticity2_weighted_radial_rms",
        "vorticity2_weighted_axial_rms",
        "vorticity2_weighted_aspect",
    ):
        excursion = summary[f"cubic_quadratic_{key}_excursion_fraction"]
        assert 0.12 < excursion < 0.23

    # At the finest grid the cubic field is visually almost static by smooth
    # vorticity-weighted geometry, while the quadratic comparison moves farther.
    assert abs(summary["cubic_baseline_vorticity2_weighted_radial_rms_ratio"] - 1.0) < 0.003
    assert abs(summary["cubic_baseline_vorticity2_weighted_axial_rms_ratio"] - 1.0) < 0.003
    assert abs(summary["cubic_baseline_vorticity2_weighted_aspect_ratio"] - 1.0) < 0.003
    assert abs(summary["cubic_baseline_radial_q99_ratio"] - 1.0) < 0.02
    assert abs(summary["cubic_baseline_axial_q99_ratio"] - 1.0) < 0.02

    # The late cubic off-keyframe does not create a support-collar isosurface.
    # Its tiny collar-vorticity budget is lower than static here, while the
    # quadratic late-half field increases the same budget.
    assert fine["cubic"]["superlevel_collar_voxel_fraction"] == 0.0
    assert fine["cubic"]["whole_grid_collar_vorticity2_fraction"] < 3.0e-5
    assert fine["cubic"]["whole_grid_collar_vorticity2_fraction"] < fine["baseline"]["whole_grid_collar_vorticity2_fraction"]
    assert fine["quadratic"]["whole_grid_collar_vorticity2_fraction"] > fine["baseline"]["whole_grid_collar_vorticity2_fraction"]
    assert summary["cubic_quadratic_whole_grid_collar_vorticity2_fraction_excursion_fraction"] < 0.25

    # Smooth metrics remain stable from 65^3 -> 81^3.  Radial q99 is allowed a
    # slightly wider 2% gate because thresholded Cartesian quantiles move in
    # grid-sized jumps.  Axial q99 is intentionally not a convergence gate: it
    # moves by one z-plane even while the smooth axial RMS stays below 1% change.
    assert summary["cubic_radial_q99_mid_to_fine_relative_change"] < 0.02
    assert summary["cubic_vorticity2_weighted_radial_rms_mid_to_fine_relative_change"] < 0.01
    assert summary["cubic_vorticity2_weighted_axial_rms_mid_to_fine_relative_change"] < 0.01
    assert summary["cubic_vorticity2_weighted_aspect_mid_to_fine_relative_change"] < 0.01
    assert summary["cubic_whole_grid_collar_vorticity2_fraction_mid_to_fine_relative_change"] < 0.05

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
