import numpy as np
import pytest

from openai_ns_reconstruction.constrained_parallel_vector_core_trace import (
    diagnose_parallel_vector_core_trace,
)


def _helical_rotation(points, time, *, omega=2.0, x0=0.375, y0=-0.1875, axial=0.8):
    del time
    p = np.asarray(points, dtype=float)
    x = p[:, 0] - x0
    y = p[:, 1] - y0
    return np.column_stack((-omega * y, omega * x, np.full(p.shape[0], axial)))


def test_offset_helical_rotation_recovers_axial_probe_trace():
    out = diagnose_parallel_vector_core_trace(
        _helical_rotation,
        0.5,
        box_half_width=1.5,
        grid_size=17,
        provenance="analytic offset helical rotation",
    )
    pts = np.asarray(out.trace_points)
    assert out.trace_point_count == 15
    assert out.valid_slice_fraction == pytest.approx(1.0)
    assert np.allclose(pts[:, 0], 0.375, atol=1e-14)
    assert np.allclose(pts[:, 1], -0.1875, atol=1e-14)
    assert out.trace_swirl_strength_median == pytest.approx(2.0, abs=1e-12)
    assert out.trace_parallel_residual_max < 1e-12
    assert out.transverse_step_rms < 1e-14
    assert out.divergence_rms < 1e-12
    assert not out.visualization_ready
    assert not out.pde_validated
    assert not out.openai_field_identified


def test_simple_shear_has_no_false_swirl_trace():
    def shear(points, time):
        del time
        p = np.asarray(points, dtype=float)
        return np.column_stack(
            (1.7 * p[:, 1], np.zeros(p.shape[0]), np.ones(p.shape[0]))
        )

    out = diagnose_parallel_vector_core_trace(
        shear,
        0.5,
        grid_size=15,
        provenance="analytic simple shear",
    )
    assert out.swirling_strength_max == pytest.approx(0.0)
    assert out.trace_point_count == 0
    assert out.trace_points == ()
    assert out.transverse_radius_rms is None
    assert out.trace_parallel_residual_median is None


def test_velocity_amplitude_scaling_preserves_trace_geometry_and_normalized_residual():
    base = diagnose_parallel_vector_core_trace(
        _helical_rotation,
        0.5,
        grid_size=17,
        provenance="base",
    )

    def scaled(points, time):
        return 3.25 * _helical_rotation(points, time)

    scaled_out = diagnose_parallel_vector_core_trace(
        scaled,
        0.5,
        grid_size=17,
        provenance="scaled",
    )

    assert scaled_out.trace_points == base.trace_points
    assert np.allclose(
        scaled_out.trace_parallel_residual,
        base.trace_parallel_residual,
        atol=1e-13,
    )
    assert scaled_out.swirling_strength_max == pytest.approx(
        3.25 * base.swirling_strength_max,
        rel=1e-12,
    )


@pytest.mark.parametrize(
    "velocity,time,box_half_width,grid_size,provenance",
    [
        (lambda p, t: np.zeros_like(p), 0.5, 1.5, 17, "zero"),
        (lambda p, t: np.zeros((len(p), 2)), 0.5, 1.5, 17, "bad shape"),
        (lambda p, t: np.full_like(p, np.nan), 0.5, 1.5, 17, "nan"),
        (_helical_rotation, 0.2, 1.5, 17, "bad time"),
        (_helical_rotation, 0.5, 2.1, 17, "bad box"),
        (_helical_rotation, 0.5, 1.5, 16, "even"),
        (_helical_rotation, 0.5, 1.5, 17, ""),
    ],
)
def test_fail_closed_inputs(velocity, time, box_half_width, grid_size, provenance):
    with pytest.raises(ValueError):
        diagnose_parallel_vector_core_trace(
            velocity,
            time,
            box_half_width=box_half_width,
            grid_size=grid_size,
            provenance=provenance,
        )
