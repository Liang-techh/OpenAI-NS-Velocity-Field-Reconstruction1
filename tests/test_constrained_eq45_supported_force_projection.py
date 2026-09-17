import numpy as np

from openai_ns_reconstruction.constrained_eq45_supported_force_projection import (
    FORCE_BOUNDS,
    _curl_force,
    _fit_probe_set,
    audit_supported_eq45_restricted_force,
)
from openai_ns_reconstruction.constrained_eq45_supported_vorticity import _probe_set
from openai_ns_reconstruction.constrained_force import RestrictedForce


def test_restricted_force_curl_is_linear_in_registered_two_parameter_family():
    points = np.asarray(
        [[0.31, -0.19, 0.27], [1.71, 0.22, -0.33], [0.42, -0.14, 1.71]],
        dtype=float,
    )
    times = np.asarray([0.39, 0.51, 0.63], dtype=float)
    step = 0.005
    a = 1.7
    c = 2.3
    combined = _curl_force(RestrictedForce(a=a, c=c), points, times, step)
    separate = (
        a * _curl_force(RestrictedForce(a=1.0, c=0.0), points, times, step)
        + c * _curl_force(RestrictedForce(a=0.0, c=1.0), points, times, step)
    )
    np.testing.assert_allclose(combined, separate, rtol=1e-10, atol=1e-10)


def test_supported_force_fit_and_holdout_probe_tuples_are_disjoint():
    fit_points, fit_times, fit_labels = _fit_probe_set()
    hold_points, hold_times, hold_labels = _probe_set()
    assert fit_points.shape == hold_points.shape == (16, 3)
    assert fit_times.shape == hold_times.shape == (16,)
    assert set(fit_labels) == set(hold_labels) == {
        "plateau",
        "radial_collar",
        "axial_collar",
        "corner_collar",
    }
    for point, time in zip(fit_points, fit_times):
        same_point = np.all(hold_points == point, axis=1)
        same_time = hold_times == time
        assert not np.any(same_point & same_time)


def test_supported_eq45_restricted_force_projection_is_bounded_frozen_and_truthful():
    report = audit_supported_eq45_restricted_force()

    assert report["schema"] == "eq45_supported_restricted_force_projection_v1"
    assert report["supported_child_sha256"] == (
        "2fdff812c22131d56eac1b7d6e455187ff3207500c6385aa9abf282a8e3d7b1d"
    )
    assert report["parent_sha256"] == (
        "48f1845fcfd71ec95495c98bb7aac3fca4653e748a8856a5421234e9601525a7"
    )
    assert report["force_family"] == "preregistered_restricted_two_parameter_family"
    assert report["force_parameters_fitted"] == ["a", "c"]
    assert report["force_bounds"] == list(FORCE_BOUNDS)
    assert report["fit_and_holdout_separate"] is True
    assert report["holdout_force_refit"] is False
    assert report["untapered_force_transferred"] is False
    assert report["pressure_fitted"] is False
    assert report["new_force_direction_added"] is False
    assert report["residual_defined_force_allowed"] is False

    fit = report["fit"]
    assert FORCE_BOUNDS[0] <= fit["a"] <= FORCE_BOUNDS[1]
    assert FORCE_BOUNDS[0] <= fit["c"] <= FORCE_BOUNDS[1]
    assert fit["design_rank"] == 2
    assert np.isfinite(fit["design_condition"])
    assert fit["design_condition"] < 1.0e6
    assert fit["probe_count"] == 16
    assert fit["after"]["rms"] < fit["before"]["rms"]
    assert fit["rms_reduction_fraction"] > 0.0

    rows = report["holdout_rows"]
    np.testing.assert_allclose(
        [row["spatial_step"] for row in rows],
        [0.02, 0.01, 0.005],
        rtol=0.0,
        atol=0.0,
    )
    assert report["holdout_probe_count"] == 16
    # Replays the exact zero-force supported-child baseline from Agent 3 #129.
    np.testing.assert_allclose(
        rows[-1]["before"]["rms"],
        3.58728819303,
        rtol=2e-9,
        atol=2e-9,
    )
    for row in rows:
        assert np.isfinite(row["after"]["rms"])
        assert np.isfinite(row["after"]["max"])
        # Same-family projection must not be allowed to hide a catastrophic holdout regression.
        assert row["after"]["rms"] < 1.5 * row["before"]["rms"]
        assert set(row["by_region"]) == {
            "plateau",
            "radial_collar",
            "axial_collar",
            "corner_collar",
        }

    assert report["formal_pde_gate_assessed"] is False
    assert report["velocity_changed"] is False
    assert report["support_transform_changed"] is False
    assert report["visualization_ready"] is False
    assert report["visual_correspondence_verified"] is False
    assert report["pde_validated"] is False
    assert report["paper_exact"] is False
    assert report["openai_field_identified"] is False
    assert report["blowup_proved"] is False
