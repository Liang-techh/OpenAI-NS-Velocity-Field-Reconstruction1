import numpy as np
import pytest

from openai_ns_reconstruction.constrained_streamline_resolution import (
    audit_streamline_grid_resolution,
)


COMMON = {
    "seed_xyz": (0.4, 0.0, -0.2),
    "time": 0.5,
    "field_identity": "candidate-sha256:abc",
    "integration_contract": "same RK4 direction-field contract / same seed / same bounds",
    "resolution_provenance": "frozen Cartesian grid points per axis",
    "provenance": "CR-A9-028 regression",
}


def _line(n, y=0.0, reverse=False):
    x = np.linspace(0.0, 2.0, n)
    if reverse:
        x = x[::-1]
    return np.column_stack((x, np.full(n, y), np.zeros(n)))


def test_same_geometry_at_different_sampling_is_zero_discrepancy():
    report = audit_streamline_grid_resolution(
        {24: _line(4), 48: _line(7), 96: _line(13)},
        sample_count=65,
        **COMMON,
    )
    assert report.finest_resolution == 96
    np.testing.assert_allclose(report.arclengths, 2.0, atol=1e-14)
    np.testing.assert_allclose(report.symmetric_hausdorff, 0.0, atol=1e-14)
    np.testing.assert_allclose(report.pointwise_rms, 0.0, atol=1e-14)
    assert report.metadata["visualization_ready"] is False
    assert report.metadata["pde_validated"] is False
    assert report.resolutions.flags.writeable is False
    assert report.resampled_points.flags.writeable is False


def test_known_lateral_offsets_remain_visible_in_shared_physical_frame():
    report = audit_streamline_grid_resolution(
        {32: _line(5, 0.2), 48: _line(9, 0.1), 64: _line(17, 0.0)},
        sample_count=65,
        **COMMON,
    )
    np.testing.assert_allclose(report.symmetric_hausdorff, [0.2, 0.1, 0.0], atol=1e-14)
    np.testing.assert_allclose(report.pointwise_rms, [0.2, 0.1, 0.0], atol=1e-14)
    np.testing.assert_allclose(report.start_error, [0.2, 0.1, 0.0], atol=1e-14)
    np.testing.assert_allclose(report.end_error, [0.2, 0.1, 0.0], atol=1e-14)
    np.testing.assert_allclose(
        report.normalized_hausdorff_by_finest_arclength,
        [0.1, 0.05, 0.0],
        atol=1e-14,
    )
    assert report.metadata["automatic_alignment"] is False


def test_orientation_is_not_silently_reversed():
    report = audit_streamline_grid_resolution(
        {32: _line(7, reverse=True), 48: _line(9), 64: _line(17)},
        sample_count=65,
        **COMMON,
    )
    assert report.symmetric_hausdorff[0] < 1e-14
    assert report.pointwise_rms[0] > 1.0
    assert report.start_error[0] == pytest.approx(2.0)
    assert report.end_error[0] == pytest.approx(2.0)
    assert report.metadata["automatic_reversal"] is False


def test_helix_discrepancy_decreases_with_grid_resolution():
    def helix(n):
        s = np.linspace(0.0, 2.0 * np.pi, n)
        return np.column_stack((np.cos(s), np.sin(s), 0.2 * s))

    report = audit_streamline_grid_resolution(
        {16: helix(9), 32: helix(25), 64: helix(129)},
        sample_count=129,
        **COMMON,
    )
    assert report.symmetric_hausdorff[1] < report.symmetric_hausdorff[0]
    assert report.pointwise_rms[1] < report.pointwise_rms[0]
    assert report.arclength_relative_error_to_finest[1] < report.arclength_relative_error_to_finest[0]
    assert report.symmetric_hausdorff[-1] == 0.0
    assert report.pointwise_rms[-1] == 0.0


def test_fail_closed_malformed_or_unprovenanced_inputs():
    valid = {16: _line(4), 32: _line(5), 64: _line(7)}

    with pytest.raises(ValueError):
        audit_streamline_grid_resolution({16: _line(4), 32: _line(5)}, sample_count=65, **COMMON)
    with pytest.raises(ValueError):
        audit_streamline_grid_resolution(
            {16: _line(4), 32: _line(5), 64: _line(7)},
            sample_count=64,
            **COMMON,
        )

    duplicate = _line(5)
    duplicate[2] = duplicate[1]
    with pytest.raises(ValueError):
        audit_streamline_grid_resolution(
            {16: duplicate, 32: _line(5), 64: _line(7)}, sample_count=65, **COMMON
        )

    nonfinite = _line(5)
    nonfinite[2, 1] = np.nan
    with pytest.raises(ValueError):
        audit_streamline_grid_resolution(
            {16: nonfinite, 32: _line(5), 64: _line(7)}, sample_count=65, **COMMON
        )

    bad_resolution = {16.0: _line(4), 32: _line(5), 64: _line(7)}
    with pytest.raises(ValueError):
        audit_streamline_grid_resolution(bad_resolution, sample_count=65, **COMMON)

    for field, value in [
        ("field_identity", ""),
        ("integration_contract", ""),
        ("resolution_provenance", ""),
        ("provenance", ""),
    ]:
        kwargs = dict(COMMON)
        kwargs[field] = value
        with pytest.raises(ValueError):
            audit_streamline_grid_resolution(valid, sample_count=65, **kwargs)

    kwargs = dict(COMMON)
    kwargs["seed_xyz"] = (0.0, np.nan, 0.0)
    with pytest.raises(ValueError):
        audit_streamline_grid_resolution(valid, sample_count=65, **kwargs)

    kwargs = dict(COMMON)
    kwargs["time"] = np.inf
    with pytest.raises(ValueError):
        audit_streamline_grid_resolution(valid, sample_count=65, **kwargs)
