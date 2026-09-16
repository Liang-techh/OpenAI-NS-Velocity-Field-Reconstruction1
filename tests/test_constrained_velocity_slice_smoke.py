from pathlib import Path

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_velocity_slice_smoke import (
    render_velocity_slice_sequence,
    sample_velocity_slice_sequence,
)


def _inward_helical_field(points, time):
    points = np.asarray(points, dtype=float)
    x = points[..., 0]
    y = points[..., 1]
    z = points[..., 2]
    a = 0.20 + 0.10 * float(time)
    swirl = 0.90 + 0.20 * float(time)
    # Divergence-free linear calibration: u_r=-a*r, u_theta=swirl*r, w=2*a*z.
    return np.stack(
        (
            -a * x - swirl * y,
            -a * y + swirl * x,
            2.0 * a * z,
        ),
        axis=-1,
    )


def test_sample_sequence_uses_direct_velocity_and_shared_scale():
    result = sample_velocity_slice_sequence(
        _inward_helical_field,
        times=(0.25, 0.5, 0.75),
        horizontal_extent=1.5,
        axial_extent=2.0,
        grid_size=21,
    )

    assert result.meridional_velocity.shape == (3, 21, 21, 3)
    assert result.equatorial_velocity.shape == (3, 21, 21, 3)
    assert np.all(np.diff(result.speed_max_by_time) > 0.0)
    assert result.global_speed_max == pytest.approx(result.speed_max_by_time[-1])
    assert not result.meridional_velocity.flags.writeable
    assert not result.equatorial_velocity.flags.writeable

    metadata = result.metadata
    assert metadata["claim_scope"] == "visualization_smoke_only"
    assert metadata["visual_correspondence_verified"] is False
    assert metadata["pde_validated"] is False
    assert metadata["paper_exact"] is False


def test_render_sequence_writes_nonempty_png(tmp_path: Path):
    output = tmp_path / "eq45_slice_smoke.png"
    metadata = render_velocity_slice_sequence(
        _inward_helical_field,
        output,
        times=(0.25, 0.5, 0.75),
        horizontal_extent=1.0,
        axial_extent=1.5,
        grid_size=17,
        quiver_stride=4,
        dpi=90,
    )

    assert output.exists()
    assert output.stat().st_size > 1000
    assert metadata["visualization_smoke_rendered"] is True
    assert metadata["shared_speed_normalization"] is True
    assert metadata["visual_correspondence_verified"] is False
    assert metadata["pde_validated"] is False


def test_fail_closed_on_zero_field_bad_output_and_bad_controls(tmp_path: Path):
    def zero_field(points, time):
        return np.zeros_like(points, dtype=float)

    def wrong_shape(points, time):
        return np.zeros(points.shape[:-1] + (2,), dtype=float)

    with pytest.raises(ValueError, match="inactive"):
        sample_velocity_slice_sequence(zero_field, times=(0.25, 0.75), grid_size=17)
    with pytest.raises(ValueError, match="same"):
        sample_velocity_slice_sequence(wrong_shape, times=(0.25, 0.75), grid_size=17)
    with pytest.raises(ValueError, match="strictly increasing"):
        sample_velocity_slice_sequence(_inward_helical_field, times=(0.5, 0.5), grid_size=17)
    with pytest.raises(ValueError, match="odd"):
        sample_velocity_slice_sequence(_inward_helical_field, times=(0.25, 0.75), grid_size=16)
    with pytest.raises(ValueError, match=".png"):
        render_velocity_slice_sequence(
            _inward_helical_field,
            tmp_path / "bad.svg",
            times=(0.25, 0.75),
            grid_size=17,
        )
