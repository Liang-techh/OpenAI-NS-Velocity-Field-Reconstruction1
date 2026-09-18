from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_st051b_redistribution_visualization_grid import (
    GAIN,
    INNER_WINDOW,
    OUTER_WINDOW,
    SOURCE_ALPHA,
    _transform_velocity_values,
    compact_bump,
    export_velocity_netcdf,
    load_visualization_netcdf,
    redistribution_profile,
)


def test_frozen_profile_has_intended_inner_mid_outer_signs():
    values = redistribution_profile(np.array([0.6, 0.9, 1.2]))
    assert values[0] > 0.0
    assert values[1] > 0.0
    assert values[2] < 0.0
    assert SOURCE_ALPHA == pytest.approx(2.520520814687742, rel=0, abs=0)
    assert GAIN == pytest.approx(0.025, rel=0, abs=0)
    assert INNER_WINDOW == (0.30, 1.05)
    assert OUTER_WINDOW == (0.95, 1.85)


def test_compact_bump_is_zero_outside_window_and_one_at_midpoint():
    lo, hi = 0.30, 1.05
    mid = 0.5 * (lo + hi)
    values = compact_bump(np.array([lo, mid, hi, lo - 0.1, hi + 0.1]), lo, hi)
    assert np.array_equal(values[[0, 2, 3, 4]], np.zeros(4))
    assert values[1] == pytest.approx(1.0, rel=0, abs=1e-15)


def test_transform_changes_only_azimuthal_component_before_common_scale():
    points = np.array([[0.6, 0.0, 0.2], [0.0, 0.9, -0.2], [0.0, 0.0, 0.1]])
    values = np.array([[2.0, 3.0, 4.0], [-3.0, 2.0, 5.0], [0.0, 0.0, 6.0]])
    scale = 1.001
    out = _transform_velocity_values(points, values, gain=GAIN, scale=scale)

    # At x>0, x is radial and y is azimuthal.
    factor_r06 = 1.0 + GAIN * float(redistribution_profile(np.array([0.6]))[0])
    assert out[0, 0] == pytest.approx(scale * 2.0)
    assert out[0, 1] == pytest.approx(scale * 3.0 * factor_r06)
    assert out[0, 2] == pytest.approx(scale * 4.0)

    # On +y, y is radial and -x is azimuthal; radial and axial stay common-scaled.
    factor_r09 = 1.0 + GAIN * float(redistribution_profile(np.array([0.9]))[0])
    assert out[1, 1] == pytest.approx(scale * 2.0)
    assert out[1, 0] == pytest.approx(scale * -3.0 * factor_r09)
    assert out[1, 2] == pytest.approx(scale * 5.0)

    # Axis is not given an arbitrary cylindrical direction.
    assert np.array_equal(out[2], scale * values[2])


def test_netcdf_roundtrip_is_nonzero_and_checksum_guarded(tmp_path: Path):
    def fake_velocity(points: np.ndarray, time: float) -> np.ndarray:
        points = np.asarray(points, float)
        return np.column_stack(
            (
                (1.0 + time) * points[:, 1],
                -(1.0 + time) * points[:, 0],
                0.25 * points[:, 2],
            )
        )

    path = tmp_path / "grid.nc"
    result = export_velocity_netcdf(fake_velocity, path, resolution=5)
    loaded = load_visualization_netcdf(path)
    assert result["resolution"] == 5
    assert loaded["u"].shape == (3, 5, 5, 5)
    assert loaded["grid_sha256"] == result["grid_sha256"]
    assert max(result["speed_max_by_time"]) > 0.0
    assert np.max(np.abs(loaded["speed"] - np.sqrt(loaded["u"] ** 2 + loaded["v"] ** 2 + loaded["w"] ** 2))) < 1e-12


def test_grid_contract_rejects_even_or_zero_velocity(tmp_path: Path):
    def zero(points: np.ndarray, time: float) -> np.ndarray:
        return np.zeros_like(points)

    with pytest.raises(ValueError, match="odd"):
        export_velocity_netcdf(zero, tmp_path / "even.nc", resolution=6)
    with pytest.raises(ValueError, match="exact-zero"):
        export_velocity_netcdf(zero, tmp_path / "zero.nc", resolution=5)
