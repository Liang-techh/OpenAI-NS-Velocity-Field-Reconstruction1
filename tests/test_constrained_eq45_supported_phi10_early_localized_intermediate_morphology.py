import numpy as np

from openai_ns_reconstruction.constrained_eq45_supported_phi10_early_localized_intermediate_morphology import (
    DEFAULT_RESOLUTIONS,
    DEFAULT_TIME,
    audit_early_localized_intermediate_morphology,
)


def test_early_localized_intermediate_morphology_is_small_and_resolved():
    report = audit_early_localized_intermediate_morphology()

    assert report["schema"] == "eq45_supported_phi10_early_localized_intermediate_morphology_v1"
    assert report["resolutions"] == list(DEFAULT_RESOLUTIONS)
    assert report["time"] == DEFAULT_TIME
    np.testing.assert_allclose(report["tau"], 0.5, rtol=0.0, atol=1.0e-15)
    np.testing.assert_allclose(report["coefficients"]["static"], -0.3, rtol=0.0, atol=1.0e-15)
    np.testing.assert_allclose(report["coefficients"]["quadratic"], -0.125, rtol=0.0, atol=1.0e-15)
    np.testing.assert_allclose(report["coefficients"]["affine"], 0.4, rtol=0.0, atol=1.0e-15)

    summary = report["finest_summary"]
    fine = report["rows"][-1]

    # The quadratic schedule keeps the intermediate whole-domain morphology much
    # closer to the static field than the older affine schedule.  Smooth
    # vorticity-weighted metrics retain about one quarter of the affine change,
    # matching the one-quarter coefficient excursion without assuming linearity.
    for key in (
        "vorticity2_weighted_radial_rms",
        "vorticity2_weighted_axial_rms",
        "vorticity2_weighted_aspect",
    ):
        collateral = summary[f"quadratic_affine_{key}_collateral_fraction"]
        assert 0.20 < collateral < 0.32

    assert summary["quadratic_affine_radial_q99_collateral_fraction"] < 0.55
    assert summary["quadratic_affine_radial_outer_vorticity2_fraction_collateral_fraction"] < 0.35
    assert summary["quadratic_affine_whole_grid_collar_vorticity2_fraction_collateral_fraction"] < 0.35

    # At the finest grid the quadratic field moves the smooth bulk geometry only
    # modestly from static while the affine comparison changes it materially more.
    assert abs(summary["quadratic_baseline_vorticity2_weighted_radial_rms_ratio"] - 1.0) < 0.02
    assert abs(summary["quadratic_baseline_vorticity2_weighted_axial_rms_ratio"] - 1.0) < 0.01
    assert abs(summary["quadratic_baseline_vorticity2_weighted_aspect_ratio"] - 1.0) < 0.02
    assert summary["quadratic_baseline_radial_q99_ratio"] < 1.04
    assert summary["affine_baseline_radial_q99_ratio"] > summary["quadratic_baseline_radial_q99_ratio"]

    # No 25%-core superlevel point enters the support collar at this time.  The
    # tiny vorticity-squared collar fraction rises versus static, but is still
    # much smaller than in the affine field.
    assert fine["quadratic"]["superlevel_collar_voxel_fraction"] == 0.0
    assert fine["quadratic"]["whole_grid_collar_vorticity2_fraction"] < 3.0e-4
    assert fine["quadratic"]["whole_grid_collar_vorticity2_fraction"] < fine["affine"]["whole_grid_collar_vorticity2_fraction"]

    # Smooth metrics and radial q99 are stable from 65^3 -> 81^3.  Do not impose
    # a convergence gate on axial q99: its 0.05 Cartesian z-plane quantization
    # makes that thresholded statistic jump by one plane even while axial RMS is
    # stable below 1%.
    assert summary["quadratic_radial_q99_mid_to_fine_relative_change"] < 0.01
    assert summary["quadratic_vorticity2_weighted_radial_rms_mid_to_fine_relative_change"] < 0.01
    assert summary["quadratic_vorticity2_weighted_axial_rms_mid_to_fine_relative_change"] < 0.01
    assert summary["quadratic_vorticity2_weighted_aspect_mid_to_fine_relative_change"] < 0.01
    assert summary["quadratic_whole_grid_collar_vorticity2_fraction_mid_to_fine_relative_change"] < 0.05

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
