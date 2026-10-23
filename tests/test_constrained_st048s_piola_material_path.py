import copy

import numpy as np
import pytest

from openai_ns_reconstruction import constrained_st048s_piola_material_path as mod


def test_piola_coordinate_map_identity_and_orientation():
    z = np.linspace(-2.2, 2.2, 401)
    mapped0, jac0 = mod._warp_z_and_jacobian(z, 0.0)
    assert np.array_equal(mapped0, z)
    assert np.array_equal(jac0, np.ones_like(z))

    mapped, jac = mod._warp_z_and_jacobian(z, mod.PIOLA_BETA)
    assert np.all(jac > 0.0)
    assert np.all(np.diff(mapped) > 0.0)
    assert mapped[0] == pytest.approx(z[0])
    assert mapped[-1] == pytest.approx(z[-1])
    assert mod._warp_z_and_jacobian(np.array([0.0]), mod.PIOLA_BETA)[1][0] == pytest.approx(0.925)


def test_piola_velocity_matches_definition_and_beta_zero_identity():
    def parent(points, time):
        points = np.asarray(points, dtype=float)
        return np.column_stack((points[:, 0] + time, 2.0 * points[:, 1], -points[:, 2]))

    points = np.array([[0.2, -0.4, 0.3], [0.6, 0.1, -0.7]], dtype=float)
    zero = mod.PiolaWarpVelocity(parent, beta=0.0, scale=1.0)
    np.testing.assert_allclose(zero(points, 0.5), parent(points, 0.5), rtol=0.0, atol=0.0)

    warped = mod.PiolaWarpVelocity(parent, beta=mod.PIOLA_BETA, scale=1.25)
    mapped_z, jac = mod._warp_z_and_jacobian(points[:, 2], mod.PIOLA_BETA)
    mapped = points.copy()
    mapped[:, 2] = mapped_z
    expected = 1.25 * parent(mapped, 0.5)
    expected[:, :2] *= jac[:, None]
    np.testing.assert_allclose(warped(points, 0.5), expected, rtol=2e-15, atol=2e-15)


def test_material_path_contract_detects_inward_rotation_and_axial_growth():
    # Linear analytic flow stays well inside the registered box for the fixed seeds.
    a = -0.10
    omega = 0.35
    c = 0.20

    def velocity(points, time):
        points = np.asarray(points, dtype=float)
        x, y, z = points.T
        return np.column_stack((a * x - omega * y, omega * x + a * y, c * z))

    report = mod.measure_material_paths(velocity)
    assert report["path_count"] == 48
    assert report["paired_material_line_count"] == 24
    assert report["inward_path_count"] == 48
    assert report["outward_path_count"] == 0
    assert report["pair_axial_separation_growth_count"] == 24
    assert report["pair_axial_separation_shrink_count"] == 0
    assert report["mean_absolute_turns"] == pytest.approx(omega * 0.5 / (2.0 * np.pi), rel=2e-9)
    assert report["mean_pair_axial_separation_ratio"] == pytest.approx(np.exp(c * 0.5), rel=2e-9)


def test_fail_closed_inputs_and_truth_boundary():
    with pytest.raises(ValueError):
        mod._warp_z_and_jacobian(np.array([0.0]), 0.21)
    with pytest.raises(ValueError):
        mod.PiolaWarpVelocity(lambda p, t: p, beta=0.1, scale=0.0)
    with pytest.raises(TypeError):
        mod.measure_material_paths(None)

    truth = copy.deepcopy(mod.TRUTH_BOUNDARY)
    for key in (
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "source_correspondence_verified",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        assert truth[key] is False
    assert truth["visual_acceptance_threshold_defined"] is False
    assert truth["comparison_is_descriptive_not_acceptance"] is True
    assert truth["pressure_or_force_transferred_from_parent"] is False


def test_comparison_is_descriptive_and_has_no_aggregate_score():
    summary = {
        "mean_radius_change": -0.09,
        "mean_absolute_turns": 0.03,
        "maximum_absolute_turns": 0.05,
        "mean_pair_axial_separation_change": 0.08,
        "mean_pair_axial_separation_ratio": 1.13,
        "inward_path_count": 48,
        "pair_axial_separation_growth_count": 17,
        "pair_axial_separation_shrink_count": 7,
    }
    reference = {
        "mean_radius_change": -0.08,
        "mean_absolute_turns": 0.02,
        "maximum_absolute_turns": 0.04,
        "mean_pair_separation_change": 0.07,
        "mean_pair_separation_ratio": 1.12,
        "inward_path_count": 48,
        "pair_growth_count": 16,
        "pair_shrink_count": 8,
    }
    result = mod._comparison(summary, reference)
    assert "score" not in result
    assert "pass" not in result
    assert result["mean_absolute_turns_ratio"] == pytest.approx(1.5)
    assert result["pair_growth_count_delta"] == 1
