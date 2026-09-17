import pytest

from openai_ns_reconstruction.constrained_eq45_supported_phi10_temporal_replay import (
    DEFAULT_SLOPES,
    replay_supported_phi10_temporal_slopes,
)


def test_phi10_temporal_replay_selects_bounded_nonlinear_trial():
    report = replay_supported_phi10_temporal_slopes()
    selected = report["selected_trial"]
    assert selected is not None

    assert report["baseline_early_radial_outer_fraction"] == pytest.approx(
        0.057728804031023005, rel=2.0e-6, abs=1.0e-10
    )
    assert report["baseline_late_radial_outer_fraction"] == pytest.approx(
        1.977917180848768e-06, rel=2.0e-6, abs=1.0e-12
    )
    assert report["late_guard_limit"] == pytest.approx(1.0e-4, rel=0.0, abs=0.0)

    assert selected["slope"] == pytest.approx(1.4, rel=0.0, abs=0.0)
    assert selected["coefficient_start"] == pytest.approx(-1.7, rel=0.0, abs=1.0e-14)
    assert selected["coefficient_end"] == pytest.approx(1.1, rel=0.0, abs=1.0e-14)
    assert selected["early_radial_outer_fraction"] == pytest.approx(
        0.021671584584770098, rel=2.0e-6, abs=1.0e-10
    )
    assert selected["early_relative_change"] == pytest.approx(
        -0.624596681872641, rel=2.0e-6, abs=1.0e-9
    )
    assert selected["late_radial_outer_fraction"] == pytest.approx(
        6.677067957832705e-06, rel=2.0e-6, abs=1.0e-12
    )
    assert selected["midpoint_exact"] is True
    assert selected["feasible"] is True

    # The selected nonlinear trial removes most of the resolved early radial
    # branch while staying comfortably below the preregistered late guard.
    assert selected["early_radial_outer_fraction"] < 0.40 * report[
        "baseline_early_radial_outer_fraction"
    ]
    assert selected["late_radial_outer_fraction"] < report["late_guard_limit"]
    assert all(row["feasible"] for row in report["trials"])

    truth = report["truth_boundary"]
    assert truth["canonical_velocity_changed"] is False
    assert truth["trial_velocity_changed"] is True
    assert truth["production_slope_promoted"] is False
    assert truth["new_spatial_basis_added"] is False
    assert truth["openai_image_fitted"] is False
    for key in (
        "physical_support_validated",
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        assert truth[key] is False


def test_phi10_replay_rejects_invalid_slope_grid():
    with pytest.raises(ValueError, match="positive finite"):
        replay_supported_phi10_temporal_slopes(slopes=(0.8, 0.0))
    with pytest.raises(ValueError, match="unique"):
        replay_supported_phi10_temporal_slopes(slopes=(0.8, 0.8))


def test_default_slope_grid_stays_inside_screened_range():
    assert DEFAULT_SLOPES == (0.8, 1.1, 1.4, 1.7, 2.0)
