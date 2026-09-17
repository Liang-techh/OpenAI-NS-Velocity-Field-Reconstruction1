import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_supported_phi10_compact_quartic_blend_capacity import (
    audit_supported_phi10_compact_quartic_blend_capacity,
    blended_phi10_delta,
)


def test_compact_quartic_blend_gives_smooth_radial_morphology_control_without_spatial_growth():
    report = audit_supported_phi10_compact_quartic_blend_capacity()
    assert report["time"] == pytest.approx(0.3125)
    assert report["tau"] == pytest.approx(-0.75)
    assert report["representation"]["spatial_basis_increment"] == 0
    assert report["representation"]["temporal_scalar_degree_increment_if_materialized"] == 1
    assert report["representation"]["blend_weight_selected"] is False

    rows = report["rows"]
    coefficients = [row["phi10_coefficient"] for row in rows]
    assert coefficients[0] == pytest.approx(-0.9494908707865168)
    assert coefficients[-1] == pytest.approx(-1.3266671113047983)
    assert coefficients[2] == pytest.approx(0.5 * (coefficients[0] + coefficients[-1]))

    radial = [row["metrics"]["vorticity2_weighted_radial_rms"] for row in rows]
    aspect = [row["metrics"]["vorticity2_weighted_aspect"] for row in rows]
    collar = [row["metrics"]["whole_grid_collar_vorticity2_fraction"] for row in rows]
    axial = [row["metrics"]["vorticity2_weighted_axial_rms"] for row in rows]

    assert all(b < a for a, b in zip(radial, radial[1:]))
    assert all(b > a for a, b in zip(aspect, aspect[1:]))
    assert all(b < a for a, b in zip(collar, collar[1:]))
    assert axial[-1] < axial[0]

    # The visual fingerprint should move progressively rather than behaving like
    # another hidden discrete branch.  We deliberately keep this loose because
    # morphology metrics are nonlinear functions of the public velocity field.
    for key in (
        "vorticity2_weighted_radial_rms",
        "vorticity2_weighted_aspect",
        "whole_grid_collar_vorticity2_fraction",
    ):
        progress = [row[f"{key}_endpoint_progress"] for row in rows]
        assert progress[0] == pytest.approx(0.0, abs=1e-12)
        assert progress[-1] == pytest.approx(1.0, abs=1e-12)
        assert all(b > a for a, b in zip(progress, progress[1:]))
        assert report["max_endpoint_progress_minus_weight"][key] < 0.35

    # Phi enters the public velocity representation affinely along this one-mode
    # path, so the local velocity response is exceptionally well conditioned.
    for row in rows:
        assert row["public_velocity_affine_linearity_rms"] < 1.0e-12
        assert row["public_velocity_affine_linearity_max_abs"] < 1.0e-11


def test_compact_quartic_blend_preserves_anchor_and_peak_slope_cap():
    weights = (0.0, 0.25, 0.5, 0.75, 1.0)
    for weight in weights:
        assert blended_phi10_delta(-1.0, blend_weight=weight) == pytest.approx(-1.4)
        assert blended_phi10_delta(0.0, blend_weight=weight) == pytest.approx(0.0, abs=1e-14)
        assert blended_phi10_delta(0.5, blend_weight=weight) == pytest.approx(0.0, abs=1e-14)
        assert blended_phi10_delta(1.0, blend_weight=weight) == pytest.approx(0.0, abs=1e-14)

    report = audit_supported_phi10_compact_quartic_blend_capacity(resolution=25)
    derivative = report["temporal_derivative"]
    cap = max(
        derivative["quartic_peak_abs_d_delta_dtau"],
        derivative["compact_peak_abs_d_delta_dtau"],
    )
    for row in derivative["rows"]:
        assert row["peak_abs_d_delta_dtau"] <= cap * (1.0 + 2.0e-4)
        assert np.isfinite(row["integral_d_delta_dtau_squared"])
        assert row["integral_d_delta_dtau_squared"] > 0.0


def test_compact_quartic_blend_fails_closed_and_keeps_truth_boundary():
    with pytest.raises(ValueError, match=r"\[0,1\]"):
        blended_phi10_delta(-0.5, blend_weight=1.1)
    with pytest.raises(ValueError, match="strictly increasing"):
        audit_supported_phi10_compact_quartic_blend_capacity(
            blend_weights=(0.0, 0.5, 0.5, 1.0), resolution=25
        )
    with pytest.raises(ValueError, match="endpoints"):
        audit_supported_phi10_compact_quartic_blend_capacity(
            blend_weights=(0.1, 0.5, 1.0), resolution=25
        )
    with pytest.raises(ValueError, match="early transition"):
        audit_supported_phi10_compact_quartic_blend_capacity(time=0.625, resolution=25)

    report = audit_supported_phi10_compact_quartic_blend_capacity(resolution=25)
    truth = report["truth_boundary"]
    assert truth["screen_temporal_blend_degree_added"] is True
    assert truth["visualization_candidate_only"] is True
    for key in (
        "canonical_velocity_changed",
        "diagnostic_velocity_changed",
        "new_spatial_basis_added",
        "blend_weight_fitted",
        "force_or_pressure_fitted",
        "pde_objective_used_to_choose_blend",
        "public_image_fitted",
        "production_temporal_shape_promoted",
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        assert truth[key] is False
