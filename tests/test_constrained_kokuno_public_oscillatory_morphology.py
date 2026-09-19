import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_public_oscillatory_morphology import (
    FROZEN_FD4_STEPS,
    build_morphology_receipt,
    fd4_vorticity,
    frozen_offgrid_points,
)


def test_fd4_vorticity_calibrates_against_independent_polynomial_field():
    def manufactured(x, y, z, t):
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        z = np.asarray(z, dtype=float)
        t = np.asarray(t, dtype=float)
        assert x.shape == y.shape == z.shape == t.shape
        return np.column_stack((y * y + 0.0 * t, z * z, x * x))

    points = frozen_offgrid_points()
    omega = fd4_vorticity(
        manufactured,
        points["x"],
        points["y"],
        points["z"],
        points["t"],
        step=0.013,
    )
    expected = np.column_stack(
        (-2.0 * points["z"], -2.0 * points["x"], -2.0 * points["y"])
    )
    assert np.allclose(omega, expected, rtol=2.0e-13, atol=2.0e-13)


@pytest.fixture(scope="module")
def morphology_receipt():
    return build_morphology_receipt()


def test_frozen_cloud_is_offgrid_and_has_fd4_support_margin():
    points = frozen_offgrid_points()
    radius = np.hypot(points["x"], points["y"])
    assert points["x"].shape == (60,)
    assert np.unique(np.column_stack(tuple(points.values())), axis=0).shape[0] == 60
    assert float(np.min(radius)) >= 0.42 - 1.0e-14
    assert float(np.max(radius)) <= 1.12 + 1.0e-14
    assert float(np.max(np.abs(points["z"]))) <= 1.25
    # Largest centered FD4 shift is 2*0.02=.04.  These checks keep every
    # Cartesian stencil comfortably away from the public radial/axial cutoffs.
    assert float(np.min(radius)) - 0.04 > 0.15
    assert float(np.max(radius)) + 0.04 < 1.35
    assert float(np.max(np.abs(points["z"]))) + 0.04 < 2.0
    assert set(np.unique(points["t"])) == {0.31, 0.50, 0.69}


def test_public_receipt_reports_three_resolution_target_free_morphology(morphology_receipt):
    receipt = morphology_receipt
    assert receipt["schema"] == "kokuno-a2-public-oscillatory-morphology-v1"
    assert receipt["candidate"]["agent2_head"] == (
        "732800ce4990464b49c8aa32d0dff4580f6684d4"
    )
    assert receipt["protocol"]["fd4_steps"] == list(FROZEN_FD4_STEPS)
    assert receipt["protocol"]["sample_count"] == 60
    assert receipt["protocol"]["target_image_or_visual_score_used"] is False
    assert set(receipt["results_by_step"]) == {"0.02", "0.01", "0.005"}

    for metrics in receipt["results_by_step"].values():
        scalars = (
            metrics["velocity_rms"],
            metrics["velocity_speed_max"],
            metrics["vorticity_rms"],
            metrics["vorticity_abs_max"],
            metrics["axial_enstrophy_fraction"],
            metrics["radial_enstrophy_fraction"],
            metrics["mean_abs_vorticity_axis_alignment"],
        )
        assert np.all(np.isfinite(scalars))
        assert metrics["velocity_rms"] > 0.0
        assert metrics["vorticity_rms"] > 0.0
        assert metrics["vorticity_abs_max"] > 0.0
        assert len(metrics["vorticity_component_rms_xyz"]) == 3
        assert np.all(np.asarray(metrics["vorticity_component_rms_xyz"]) > 0.0)
        assert metrics["axial_enstrophy_fraction"] > 0.0
        assert metrics["radial_enstrophy_fraction"] > 0.0
        assert np.isclose(
            metrics["axial_enstrophy_fraction"]
            + metrics["radial_enstrophy_fraction"],
            1.0,
            rtol=0.0,
            atol=2.0e-15,
        )
        assert 0.0 <= metrics["mean_abs_vorticity_axis_alignment"] <= 1.0

    stability = receipt["offgrid_resolution_stability"]
    assert stability["coarse_to_medium_relative_vorticity_rms_change"] > 0.0
    assert stability["medium_to_fine_relative_vorticity_rms_change"] > 0.0
    assert np.all(np.isfinite(tuple(stability.values())))


def test_receipt_preserves_scientific_truth_boundary(morphology_receipt):
    receipt = morphology_receipt
    a4 = receipt["independent_agent4_preflight_reference"]
    assert a4["agent4_head"] == "9b0f86012c53fa8e32a19f766dbc150931870425"
    assert a4["workflow"] == 35422203621
    assert a4["public_oscillatory_preflight_passed"] is True
    assert "external independent receipt" in a4["note"]

    truth = receipt["truth_boundary"]
    assert truth["velocity_candidate_changed"] is False
    assert truth["source_formula_changed"] is False
    assert truth["agent4_scientific_guards_changed"] is False
    assert truth["three_resolution_vorticity_morphology_assessed"] is True
    assert truth["offgrid_vorticity_stability_assessed"] is True
    assert truth["morphology_acceptance_threshold_preregistered"] is False
    assert truth["visual_correspondence_assessed"] is False
    assert truth["heldout_ns_momentum_residual_assessed"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False
