from __future__ import annotations

import inspect

import numpy as np
import pytest

import openai_ns_reconstruction.kokuno_oscillatory_batch_axis_safety as batchmod


def test_batch_masks_axis_and_support_before_interior_evaluator(monkeypatch):
    field = batchmod.default_field()
    r0 = float(field.radial_inner)
    r1 = float(field.radial_outer)
    z0 = float(field.axial_lower)
    z1 = float(field.axial_upper)
    calls = []

    def guarded_velocity(x, y, z, t):
        x, y, z, t = np.broadcast_arrays(
            np.asarray(x, dtype=float),
            np.asarray(y, dtype=float),
            np.asarray(z, dtype=float),
            np.asarray(t, dtype=float),
        )
        radius = np.hypot(x, y)
        assert np.all(radius > r0)
        assert np.all(radius < r1)
        assert np.all(z > z0)
        assert np.all(z < z1)
        calls.append(int(x.size))
        return np.stack((x + 0.1 * t, y - 0.2 * t, z + 0.3 * t), axis=-1)

    monkeypatch.setattr(batchmod, "velocity_osc", guarded_velocity)
    points = np.asarray(
        [
            [(0.0, 0.0, 0.0), (0.42, 0.18, -0.4), (r0, 0.0, 0.1), (0.65, -0.2, 0.5)],
            [(0.5 * r0, 0.0, 0.0), (-0.55, 0.3, -0.7), (r1, 0.0, 0.0), (0.73, 0.22, 0.9)],
        ],
        dtype=float,
    )
    times = np.asarray(((0.31, 0.35, 0.39, 0.43), (0.47, 0.51, 0.55, 0.59)))
    out = batchmod.velocity_osc_batch(points, times)
    assert out.shape == points.shape

    radius = np.hypot(points[..., 0], points[..., 1])
    inside = (
        (radius > r0)
        & (radius < r1)
        & (points[..., 2] > z0)
        & (points[..., 2] < z1)
    )
    assert calls == [int(np.count_nonzero(inside))]
    assert np.array_equal(out[~inside], np.zeros_like(out[~inside]))
    expected_inside = np.stack(
        (
            points[..., 0] + 0.1 * times,
            points[..., 1] - 0.2 * times,
            points[..., 2] + 0.3 * times,
        ),
        axis=-1,
    )
    assert np.array_equal(out[inside], expected_inside[inside])


def test_actual_batch_matches_scalar_public_velocity_on_interior_points():
    points = np.asarray(
        [
            (0.42, 0.19, -1.10),
            (0.63, -0.27, -0.45),
            (-0.51, 0.38, 0.15),
            (-0.74, -0.22, 0.65),
            (0.91, 0.31, 1.10),
            (-1.02, 0.24, -0.80),
        ],
        dtype=float,
    )
    times = np.asarray((0.31, 0.37, 0.43, 0.49, 0.55, 0.61))
    batch = batchmod.velocity_osc_batch(points, times)
    scalar = np.stack(
        [
            np.asarray(batchmod.velocity_osc(x, y, z, t), dtype=float).reshape(3)
            for (x, y, z), t in zip(points, times)
        ]
    )
    assert batch.shape == (points.shape[0], 3)
    assert np.all(np.isfinite(batch))
    assert np.allclose(batch, scalar, rtol=0.0, atol=batchmod.SCALAR_REPLAY_ATOL)


def test_actual_axis_tiny_radius_and_support_boundaries_are_exact_zero():
    field = batchmod.default_field()
    r0 = float(field.radial_inner)
    r1 = float(field.radial_outer)
    z0 = float(field.axial_lower)
    z1 = float(field.axial_upper)
    tiny = np.nextafter(0.0, 1.0)
    points = np.asarray(
        [
            (0.0, 0.0, 0.0),
            (tiny, -tiny, 0.25),
            (0.5 * r0, 0.0, -0.4),
            (r0, 0.0, 0.0),
            (r1, 0.0, 0.0),
            (0.5 * (r0 + r1), 0.0, z0),
            (0.5 * (r0 + r1), 0.0, z1),
            (1.01 * r1, 0.0, 0.3),
        ],
        dtype=float,
    )
    out = batchmod.velocity_osc_batch(points, 0.47)
    assert out.shape == points.shape
    assert np.all(np.isfinite(out))
    assert np.array_equal(out, np.zeros_like(out))


def test_batch_supports_higher_rank_shapes_and_scalar_time():
    points = np.asarray(
        [
            [(0.42, 0.19, -0.6), (0.63, -0.27, -0.2), (-0.51, 0.38, 0.2)],
            [(-0.74, -0.22, 0.6), (0.91, 0.31, 1.0), (-1.02, 0.24, -1.0)],
        ],
        dtype=float,
    )
    out = batchmod.velocity_osc_batch(points, 0.53)
    assert out.shape == (2, 3, 3)
    assert np.all(np.isfinite(out))


def test_invalid_shape_nonfinite_and_time_broadcast_fail_closed():
    with pytest.raises(ValueError, match="shape"):
        batchmod.velocity_osc_batch(np.zeros((4, 2)), 0.4)
    with pytest.raises(ValueError, match="finite"):
        batchmod.velocity_osc_batch(np.asarray(((0.4, 0.2, np.nan),)), 0.4)
    with pytest.raises(ValueError, match="finite"):
        batchmod.velocity_osc_batch(np.asarray(((0.4, 0.2, 0.1),)), np.nan)
    with pytest.raises(ValueError, match="broadcast"):
        batchmod.velocity_osc_batch(np.zeros((2, 3)), np.zeros(3))


def test_receipt_and_public_contract_preserve_scientific_boundary():
    receipt = batchmod.materialize_batch_axis_safety_receipt()
    diagnostic = receipt["diagnostic"]
    assert diagnostic["scalar_replay_max_abs"] <= batchmod.SCALAR_REPLAY_ATOL
    assert diagnostic["masked_speed_max"] == 0.0
    assert diagnostic["interior_vector_rms"] >= batchmod.INTERIOR_SIGNAL_FLOOR
    assert diagnostic["batch_shape"] == [8, 3]
    assert receipt["semantic_sha256"] == batchmod.batch_axis_safety_sha256()

    contract = batchmod.public_contract()
    assert contract["forbidden_inputs_present"] == []
    assert contract["vectorized_batch_api"] is True
    assert contract["shape_preserving"] is True
    assert contract["time_broadcast_supported"] is True
    assert contract["axis_safe_by_strict_support_mask"] is True
    assert contract["interior_uses_existing_public_velocity"] is True
    assert contract["complete_curl_reimplemented"] is False
    assert contract["mean_projection_performed"] is False
    assert contract["radial_inverse_performed"] is False
    assert contract["correction_velocity_constructed"] is False
    assert contract["complete_ns_residual"] is False
    assert contract["paper_exact"] is False
    assert contract["pde_validated"] is False

    params = inspect.signature(batchmod.velocity_osc_batch).parameters
    assert list(params) == ["points", "time"]
