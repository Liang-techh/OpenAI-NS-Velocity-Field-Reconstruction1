import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_supported_collar_vorticity import (
    governed_supported_seed,
)
from openai_ns_reconstruction.constrained_eq45_supported_phi10_c2_compact_temporal_capacity import (
    audit_supported_phi10_c2_compact_temporal_capacity,
    compact_phi10_delta,
    compact_phi10_delta_d2tau,
    compact_phi10_delta_dtau,
    compact_phi10_snapshot,
    compact_window_parameters,
)


def test_slope_capped_c2_window_has_checked_capacity_tradeoff():
    parameters = compact_window_parameters()

    assert parameters["quartic_nullspace_coefficient"] == pytest.approx(
        -0.2831460674157303, rel=0.0, abs=1.0e-14
    )
    assert parameters["quartic_peak_abs_d_delta_dtau"] == pytest.approx(
        3.8827715355805243, rel=1.0e-12
    )
    assert parameters["quartic_peak_abs_d_delta_dtau_tau"] == pytest.approx(-1.0)
    assert parameters["window_tau_length"] == pytest.approx(
        0.6760634706279541, rel=1.0e-12
    )
    assert parameters["return_tau"] == pytest.approx(
        -0.32393652937204587, rel=1.0e-12
    )
    assert parameters["compact_to_quartic_peak_slope_ratio"] == pytest.approx(
        1.0, rel=1.0e-12
    )
    assert parameters["compact_to_quartic_slope_energy_ratio"] == pytest.approx(
        1.3118013881205237, rel=1.0e-11
    )
    assert parameters["compact_to_quartic_peak_curvature_ratio"] == pytest.approx(
        2.243968273244653, rel=1.0e-11
    )

    return_tau = parameters["return_tau"]
    np.testing.assert_allclose(compact_phi10_delta(-1.0), -1.4, rtol=0.0, atol=1.0e-14)
    np.testing.assert_allclose(compact_phi10_delta(return_tau), 0.0, rtol=0.0, atol=1.0e-14)
    np.testing.assert_allclose(compact_phi10_delta(0.0), 0.0, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(compact_phi10_delta(1.0), 0.0, rtol=0.0, atol=0.0)

    for tau in (-1.0, return_tau, 0.0, 0.5, 1.0):
        np.testing.assert_allclose(compact_phi10_delta_dtau(tau), 0.0, rtol=0.0, atol=1.0e-13)
        np.testing.assert_allclose(compact_phi10_delta_d2tau(tau), 0.0, rtol=0.0, atol=1.0e-12)


def test_compact_window_public_velocity_preserves_early_and_eliminates_post_return_collateral():
    report = audit_supported_phi10_c2_compact_temporal_capacity()
    parameters = report["parameters"]
    endpoint = report["endpoint_checks"]

    assert parameters["return_time"] == pytest.approx(
        0.4190158676569885, rel=1.0e-12
    )
    assert endpoint["early_public_velocity_rms_difference_from_quartic"] < 1.0e-13
    assert endpoint["early_public_velocity_max_abs_difference_from_quartic"] < 1.0e-12
    assert endpoint["compact_d_delta_dtau_at_start"] == pytest.approx(0.0, abs=1.0e-13)
    assert endpoint["compact_d2_delta_dtau2_at_start"] == pytest.approx(0.0, abs=1.0e-12)
    assert endpoint["compact_d_delta_dtau_at_return"] == pytest.approx(0.0, abs=1.0e-13)
    assert endpoint["compact_d2_delta_dtau2_at_return"] == pytest.approx(0.0, abs=1.0e-12)

    rows = {row["time"]: row for row in report["off_keyframe_public_velocity_checks"]}

    # The slope cap does not make the entire early transition uniformly smaller:
    # the compact C2 window holds the early state longer before returning sharply.
    assert 1.45 < rows[0.3125]["compact_to_quartic_public_velocity_delta_rms_ratio"] < 1.75

    # By the representative mid-early off-keyframe it is already less intrusive
    # than the quartic, while still remaining a nonzero public velocity field.
    assert 0.55 < rows[0.375]["compact_to_quartic_public_velocity_delta_rms_ratio"] < 0.75
    assert rows[0.375]["compact_public_velocity_delta_rms"] > 0.0

    # The compact family has returned before .4375 and is then exactly the static
    # supported public velocity, not merely close in coefficient space.
    for time in (0.4375, 0.5625, 0.6875):
        assert rows[time]["after_compact_return"] is True
        assert rows[time]["compact_public_velocity_delta_rms"] < 1.0e-14
        assert rows[time]["compact_public_velocity_max_abs_from_static"] < 1.0e-13


def test_compact_snapshot_changes_only_existing_phi10_and_truth_boundary_stays_closed():
    base = governed_supported_seed()
    trial = compact_phi10_snapshot(base, 0.375)
    mode = (1, 0)
    index = base.parent.profile_basis.mode_indices.index(mode)

    base_phi = np.asarray(base.parent.profile_basis.phi_coefficients)
    trial_phi = np.asarray(trial.parent.profile_basis.phi_coefficients)
    changed = np.flatnonzero(base_phi != trial_phi)
    np.testing.assert_array_equal(changed, np.array([index]))
    np.testing.assert_allclose(
        np.delete(trial_phi, index), np.delete(base_phi, index), rtol=0.0, atol=0.0
    )
    assert trial.taper == base.taper

    report = audit_supported_phi10_c2_compact_temporal_capacity()
    truth = report["truth_boundary"]
    assert truth["canonical_velocity_changed"] is False
    assert truth["screen_trial_velocity_changed"] is True
    assert truth["new_spatial_basis_added"] is False
    assert truth["new_fitted_temporal_parameter_added"] is False
    assert truth["force_or_pressure_fitted"] is False
    assert truth["pde_objective_used_to_choose_window"] is False
    assert truth["public_image_fitted"] is False
    assert truth["production_temporal_shape_promoted"] is False
    assert truth["visualization_ready"] is False
    assert truth["visual_correspondence_verified"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False

    with pytest.raises(ValueError):
        compact_window_parameters(early_delta=0.0)
    with pytest.raises(ValueError):
        compact_phi10_delta(-1.01)
    with pytest.raises(ValueError):
        compact_phi10_delta_dtau(np.nan)
    with pytest.raises(ValueError):
        compact_phi10_snapshot(base, 0.375, early_delta=-10.0)
    with pytest.raises(ValueError):
        audit_supported_phi10_c2_compact_temporal_capacity(probes=np.zeros((2, 3)))
