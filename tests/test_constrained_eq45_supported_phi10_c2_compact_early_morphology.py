import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_supported_phi10_c2_compact_early_morphology import (
    audit_c2_compact_early_morphology,
)


def test_c2_compact_early_morphology_is_resolved_and_more_radially_compact_than_quartic():
    report = audit_c2_compact_early_morphology()
    assert report["time"] == pytest.approx(0.3125)
    assert report["tau"] == pytest.approx(-0.75)
    assert report["compact_return_time"] == pytest.approx(0.4190158676569885)
    assert report["coefficients"] == pytest.approx(
        {
            "static": -0.3,
            "quartic": -0.9494908707865168,
            "compact": -1.3266671113047983,
        }
    )

    fine = report["rows"][-1]
    baseline = fine["baseline"]
    quartic = fine["quartic"]
    compact = fine["compact"]
    summary = report["finest_summary"]

    assert compact["radial_q99"] < quartic["radial_q99"] < baseline["radial_q99"]
    assert (
        compact["vorticity2_weighted_radial_rms"]
        < quartic["vorticity2_weighted_radial_rms"]
        < baseline["vorticity2_weighted_radial_rms"]
    )
    assert compact["vorticity2_weighted_aspect"] > quartic["vorticity2_weighted_aspect"] > baseline[
        "vorticity2_weighted_aspect"
    ]
    assert (
        compact["whole_grid_collar_vorticity2_fraction"]
        < quartic["whole_grid_collar_vorticity2_fraction"]
        < baseline["whole_grid_collar_vorticity2_fraction"]
    )

    # The apparent aspect gain is radial tightening, not extra axial reach.
    assert compact["axial_q99"] == pytest.approx(quartic["axial_q99"])
    assert quartic["axial_q99"] == pytest.approx(baseline["axial_q99"])
    assert compact["vorticity2_weighted_axial_rms"] < quartic["vorticity2_weighted_axial_rms"]

    # The compact window must buy a material morphology change beyond quartic.
    assert summary["compact_quartic_radial_q99_ratio"] < 0.98
    assert summary["compact_quartic_vorticity2_weighted_radial_rms_ratio"] < 0.98
    assert summary["compact_quartic_vorticity2_weighted_aspect_ratio"] > 1.01
    assert summary["compact_quartic_whole_grid_collar_vorticity2_fraction_ratio"] < 0.80
    assert summary["compact_vs_quartic_radial_q99_excursion_fraction"] > 1.5

    for name in ("baseline", "quartic", "compact"):
        assert fine[name]["superlevel_collar_voxel_fraction"] == 0.0

    # Smooth morphology measures are stable on the 65^3 -> 81^3 refinement.
    assert summary["compact_radial_q99_mid_to_fine_relative_change"] < 0.005
    assert summary["compact_vorticity2_weighted_radial_rms_mid_to_fine_relative_change"] < 0.005
    assert summary["compact_vorticity2_weighted_axial_rms_mid_to_fine_relative_change"] < 0.01
    assert summary["compact_vorticity2_weighted_aspect_mid_to_fine_relative_change"] < 0.01
    assert summary["compact_whole_grid_collar_vorticity2_fraction_mid_to_fine_relative_change"] < 0.05


def test_c2_compact_early_morphology_fails_closed_outside_active_window():
    with pytest.raises(ValueError, match="active compact C2 transition"):
        audit_c2_compact_early_morphology(time=0.4375)
    with pytest.raises(ValueError, match="strictly increasing"):
        audit_c2_compact_early_morphology(resolutions=(49, 49, 81))
    with pytest.raises(ValueError, match="superlevel_fraction"):
        audit_c2_compact_early_morphology(superlevel_fraction=1.0)


def test_c2_compact_early_morphology_truth_boundary_stays_false():
    report = audit_c2_compact_early_morphology(resolutions=(17, 21, 25))
    truth = report["truth_boundary"]
    assert truth["visualization_candidate_only"] is True
    for key in (
        "canonical_velocity_changed",
        "diagnostic_velocity_changed",
        "compact_temporal_shape_promoted",
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
    assert np.isfinite(report["finest_summary"]["compact_baseline_vorticity2_weighted_radial_rms_ratio"])
