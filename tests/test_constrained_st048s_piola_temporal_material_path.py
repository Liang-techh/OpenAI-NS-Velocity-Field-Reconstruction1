import copy

import numpy as np
import pytest

from openai_ns_reconstruction import constrained_st048s_piola_temporal_material_path as mod


def test_temporal_beta_schedule_and_orientation():
    assert mod.beta_at_time(0.25, 0.025) == pytest.approx(0.10)
    assert mod.beta_at_time(0.50, 0.025) == pytest.approx(0.075)
    assert mod.beta_at_time(0.75, 0.025) == pytest.approx(0.10)

    z = np.linspace(-2.2, 2.2, 401)
    for time in (0.25, 0.50, 0.75):
        mapped, jac = mod._warp_z_and_jacobian(z, mod.beta_at_time(time, 0.025))
        assert np.all(jac > 0.0)
        assert np.all(np.diff(mapped) > 0.0)
        assert mapped[0] == pytest.approx(z[0])
        assert mapped[-1] == pytest.approx(z[-1])


def test_gamma_zero_replays_fixed_beta_definition():
    def parent(points, time):
        points = np.asarray(points, dtype=float)
        return np.column_stack((points[:, 0] + time, 2.0 * points[:, 1], -points[:, 2]))

    points = np.array([[0.2, -0.4, 0.3], [0.6, 0.1, -0.7]], dtype=float)
    warped = mod.TemporalPiolaWarpVelocity(parent, gamma=0.0, scale=1.25)
    beta = mod.beta_at_time(0.4, 0.0)
    mapped_z, jac = mod._warp_z_and_jacobian(points[:, 2], beta)
    mapped = points.copy()
    mapped[:, 2] = mapped_z
    expected = 1.25 * parent(mapped, 0.4)
    expected[:, :2] *= jac[:, None]
    np.testing.assert_allclose(warped(points, 0.4), expected, rtol=2e-15, atol=2e-15)


def test_material_path_contract_detects_inward_rotation_and_axial_growth():
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
        mod.beta_at_time(0.20, 0.025)
    with pytest.raises(ValueError):
        mod.beta_at_time(0.50, 0.13)
    with pytest.raises(ValueError):
        mod.TemporalPiolaWarpVelocity(lambda p, t: p, gamma=0.025, scale=0.0)
    with pytest.raises(TypeError):
        mod.measure_material_paths(None)

    truth = copy.deepcopy(mod.TRUTH_BOUNDARY)
    for key in (
        "production_candidate_selected",
        "production_gamma_selected",
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
    assert truth["held_out_pde_residual_transferred_from_parent"] is False


def test_comparison_is_descriptive_and_has_no_aggregate_score():
    summary = {
        "mean_radius_change": -0.09,
        "mean_absolute_turns": 0.03,
        "maximum_absolute_turns": 0.05,
        "mean_pair_axial_separation_change": 0.08,
        "mean_pair_axial_separation_ratio": 1.13,
        "pair_axial_separation_growth_count": 17,
        "pair_axial_separation_shrink_count": 7,
    }
    reference = {
        "mean_radius_change": -0.08,
        "mean_absolute_turns": 0.02,
        "maximum_absolute_turns": 0.04,
        "mean_pair_separation_change": 0.07,
        "mean_pair_separation_ratio": 1.12,
        "pair_growth_count": 16,
        "pair_shrink_count": 8,
    }
    result = mod._comparison(summary, reference)
    assert "score" not in result
    assert "pass" not in result
    assert result["mean_absolute_turns_ratio"] == pytest.approx(1.5)
    assert result["pair_growth_count_delta"] == 1
