import numpy as np

from openai_ns_reconstruction.constrained_eq45_supported_phi10_temporal_envelope import (
    DEFAULT_RESOLUTIONS,
    DEFAULT_SLOPE,
    DEFAULT_TIMES,
    audit_supported_phi10_temporal_envelope,
)


def test_supported_phi10_temporal_envelope_is_resolved_and_removes_early_radial_branch():
    report = audit_supported_phi10_temporal_envelope()

    assert report["schema"] == "eq45_supported_phi10_temporal_envelope_audit_v1"
    assert report["slope"] == DEFAULT_SLOPE
    assert report["resolutions"] == list(DEFAULT_RESOLUTIONS)
    assert report["times"] == list(DEFAULT_TIMES)
    assert len(report["rows"]) == len(DEFAULT_RESOLUTIONS) * len(DEFAULT_TIMES)

    summaries = {row["time"]: row for row in report["time_summaries"]}
    early = summaries[0.25]
    middle = summaries[0.50]
    late = summaries[0.75]

    # The selected temporal direction removes the resolved early support-induced
    # radial branch rather than only moving weight inside the same envelope.
    assert early["trial_baseline_radial_outer_vorticity2_ratio"] < 0.40
    assert early["trial_baseline_collar_vorticity2_ratio"] < 0.40
    assert early["trial_baseline_radial_q99_ratio"] < 0.70
    assert abs(early["trial_baseline_axial_q99_ratio"] - 1.0) < 0.02
    assert early["trial_baseline_aspect_q99_ratio"] > 1.45
    assert early["trial_baseline_radial_rms_ratio"] < 0.90
    assert early["trial_baseline_weighted_aspect_ratio"] > 1.08

    # The morphology change is resolved at the two finest grids.
    assert early["trial_radial_q99_mid_to_fine_relative_change"] < 0.02
    assert early["trial_axial_q99_mid_to_fine_relative_change"] < 0.02
    assert early["trial_radial_rms_mid_to_fine_relative_change"] < 0.01
    assert early["trial_axial_rms_mid_to_fine_relative_change"] < 0.01

    # The affine wrapper is exactly the static supported candidate at midpoint.
    for key, value in middle.items():
        if key.startswith("trial_baseline_"):
            np.testing.assert_allclose(value, 1.0, rtol=0.0, atol=1.0e-14)

    # The late superlevel envelope is changed, so do not claim morphology
    # identity.  The outer-support vorticity remains absolutely negligible,
    # which distinguishes the central-envelope change from a new support collar.
    fine_late = next(
        row
        for row in report["rows"]
        if row["resolution"] == DEFAULT_RESOLUTIONS[-1] and row["time"] == 0.75
    )
    assert fine_late["trial"]["radial_outer_vorticity2_fraction"] < 1.0e-4
    assert fine_late["trial"]["whole_grid_collar_vorticity2_fraction"] < 1.0e-4
    assert fine_late["trial"]["superlevel_collar_voxel_fraction"] == 0.0
    assert late["trial_radial_q99_mid_to_fine_relative_change"] < 0.02
    assert late["trial_radial_rms_mid_to_fine_relative_change"] < 0.01
    assert late["trial_axial_rms_mid_to_fine_relative_change"] < 0.02

    truth = report["truth_boundary"]
    assert truth["trial_velocity_changed"] is True
    assert truth["production_slope_promoted"] is False
    for key in (
        "canonical_velocity_changed",
        "openai_image_fitted",
        "physical_support_validated",
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        assert truth[key] is False
