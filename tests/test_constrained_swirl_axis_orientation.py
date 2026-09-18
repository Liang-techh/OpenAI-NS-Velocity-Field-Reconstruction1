import math
import numpy as np
import pytest

from openai_ns_reconstruction.constrained_swirl_axis_orientation import diagnose_swirl_axis_orientation

SHA = "a" * 64


def solid_rotation(axis=(0.0, 0.0, 1.0), omega=2.0, center=(0.0, 0.0, 0.0), scale=1.0):
    a = np.asarray(axis, dtype=float)
    a = a / np.linalg.norm(a)
    c = np.asarray(center, dtype=float)

    def velocity(points, time):
        p = np.asarray(points, dtype=float) - c
        return scale * omega * np.cross(np.broadcast_to(a, p.shape), p)

    return velocity


def test_known_vertical_rotation_axis_is_recovered():
    result = diagnose_swirl_axis_orientation(
        solid_rotation(axis=(0, 0, 1), omega=2.5),
        time=0.5,
        grid_size=9,
        candidate_sha256=SHA,
        provenance="analytic rigid rotation",
    )
    assert result["swirl_axis_detected"] is True
    assert result["swirl_active_fraction"] == pytest.approx(1.0)
    assert result["lambda_ci_rms"] == pytest.approx(2.5, rel=1e-12, abs=1e-12)
    assert result["weighted_mean_abs_axis_alignment"] == pytest.approx(1.0, abs=1e-12)
    assert result["weighted_mean_tilt_degrees"] == pytest.approx(0.0, abs=1e-10)
    assert result["principal_axis_abs_alignment"] == pytest.approx(1.0, abs=1e-12)
    assert result["orientation_concentration"] == pytest.approx(1.0, abs=1e-12)
    assert result["divergence_rms_same_operator"] == pytest.approx(0.0, abs=1e-12)
    assert all(value is False for value in result["truth_states"].values())


def test_oblique_axis_and_velocity_scaling_leave_orientation_unchanged():
    axis = np.array([1.0, -2.0, 3.0])
    axis /= np.linalg.norm(axis)
    kwargs = dict(
        time=0.625,
        bounds=((-0.8, 0.8), (-0.8, 0.8), (-0.8, 0.8)),
        grid_size=9,
        target_axis=axis,
        candidate_sha256=SHA,
        provenance="oblique rigid rotation",
    )
    a = diagnose_swirl_axis_orientation(solid_rotation(axis=axis, omega=1.7, scale=1.0), **kwargs)
    b = diagnose_swirl_axis_orientation(solid_rotation(axis=axis, omega=1.7, scale=3.25), **kwargs)
    assert a["principal_axis_abs_alignment"] == pytest.approx(1.0, abs=2e-12)
    assert a["weighted_mean_abs_axis_alignment"] == pytest.approx(1.0, abs=2e-12)
    assert b["principal_axis_abs_alignment"] == pytest.approx(a["principal_axis_abs_alignment"], abs=2e-12)
    assert b["weighted_mean_tilt_degrees"] == pytest.approx(a["weighted_mean_tilt_degrees"], abs=1e-9)
    assert b["lambda_ci_rms"] == pytest.approx(3.25 * a["lambda_ci_rms"], rel=1e-12)


def test_simple_shear_has_no_false_swirl_axis():
    def shear(points, time):
        p = np.asarray(points, dtype=float)
        out = np.zeros_like(p)
        out[:, 0] = 2.0 * p[:, 1]
        return out

    result = diagnose_swirl_axis_orientation(
        shear,
        time=0.5,
        grid_size=9,
        candidate_sha256=SHA,
        provenance="analytic simple shear",
    )
    assert result["velocity_rms"] > 0.0
    assert result["lambda_ci_max"] == pytest.approx(0.0, abs=1e-14)
    assert result["swirl_active_points"] == 0
    assert result["swirl_axis_detected"] is False
    assert result["principal_unsigned_axis"] is None
    assert result["weighted_mean_abs_axis_alignment"] is None


def test_wrong_declared_axis_reports_geometric_mismatch_without_fitting_it_away():
    result = diagnose_swirl_axis_orientation(
        solid_rotation(axis=(1, 0, 0), omega=1.0),
        time=0.5,
        grid_size=7,
        target_axis=(0, 0, 1),
        candidate_sha256=SHA,
        provenance="horizontal rigid rotation versus declared z axis",
    )
    assert result["swirl_axis_detected"] is True
    assert result["weighted_mean_abs_axis_alignment"] == pytest.approx(0.0, abs=1e-12)
    assert result["weighted_mean_tilt_degrees"] == pytest.approx(90.0, abs=1e-10)
    assert result["principal_axis_abs_alignment"] == pytest.approx(0.0, abs=1e-12)


def test_fail_closed_inputs_and_velocity_outputs():
    good = solid_rotation()
    base = dict(time=0.5, grid_size=7, candidate_sha256=SHA, provenance="guard test")
    with pytest.raises(ValueError):
        diagnose_swirl_axis_orientation(good, **{**base, "time": 0.2})
    with pytest.raises(ValueError):
        diagnose_swirl_axis_orientation(good, **{**base, "grid_size": 8})
    with pytest.raises(ValueError):
        diagnose_swirl_axis_orientation(good, **{**base, "bounds": ((-2.1, 1), (-1, 1), (-1, 1))})
    with pytest.raises(ValueError):
        diagnose_swirl_axis_orientation(good, **{**base, "target_axis": (0, 0, 0)})
    with pytest.raises(ValueError):
        diagnose_swirl_axis_orientation(good, **{**base, "candidate_sha256": "bad"})
    with pytest.raises(ValueError):
        diagnose_swirl_axis_orientation(good, **{**base, "provenance": " "})

    def zero(points, time):
        return np.zeros_like(points, dtype=float)

    with pytest.raises(ValueError):
        diagnose_swirl_axis_orientation(zero, **base)

    def bad_shape(points, time):
        return np.zeros((len(points), 2))

    with pytest.raises(ValueError):
        diagnose_swirl_axis_orientation(bad_shape, **base)

    def nonfinite(points, time):
        out = good(points, time)
        out[0, 0] = math.nan
        return out

    with pytest.raises(ValueError):
        diagnose_swirl_axis_orientation(nonfinite, **base)
