from pathlib import Path

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_candidate_visual_compare import (
    render_shared_frame_candidate_comparison,
    sample_shared_frame_candidates,
)


def _linear_field(scale):
    def velocity(points, time):
        points = np.asarray(points, dtype=float)
        a = float(scale) * (1.0 + 0.2 * (float(time) - 0.25))
        x = points[:, 0]
        y = points[:, 1]
        z = points[:, 2]
        return np.column_stack((-a * x, 0.5 * a * x + 0.1 * a * y, a * z))

    return velocity


def _sample(candidates=None, **kwargs):
    if candidates is None:
        candidates = [
            ("base", _linear_field(1.0), "manufactured-base-v1"),
            ("strong", _linear_field(2.0), "manufactured-strong-v1"),
        ]
    defaults = dict(
        x=np.linspace(-1.0, 1.0, 9),
        z=np.linspace(-1.5, 1.5, 11),
        times=[0.25, 0.5],
        frame_provenance="fixed-y0-frame-v1",
    )
    defaults.update(kwargs)
    return sample_shared_frame_candidates(candidates, **defaults)


def test_shared_scale_preserves_candidate_amplitude_difference():
    result = _sample()
    assert result.velocity.shape == (2, 2, 11, 9, 3)
    assert result.speed.shape == (2, 2, 11, 9)
    np.testing.assert_allclose(result.panel_max_speed[1], 2.0 * result.panel_max_speed[0])
    assert result.shared_speed_vmax == pytest.approx(np.max(result.panel_max_speed[1]))
    assert np.max(result.normalized_speed[0]) == pytest.approx(0.5)
    assert np.max(result.normalized_speed[1]) == pytest.approx(1.0)
    assert result.per_panel_autoscaling is False
    assert result.normalization_scope == "single_global_speed_scale_across_all_candidates_and_times"


def test_identical_points_times_and_declared_y_plane_are_used_for_every_candidate():
    calls = {"a": [], "b": []}

    def recorder(name, multiplier):
        def velocity(points, time):
            calls[name].append((points.copy(), float(time)))
            return multiplier * np.column_stack((points[:, 0] + 2.0, points[:, 1] + 1.0, points[:, 2] + 3.0))

        return velocity

    result = _sample(
        candidates=[
            ("a", recorder("a", 1.0), "candidate-a"),
            ("b", recorder("b", 1.5), "candidate-b"),
        ],
        y_plane=0.375,
        times=[0.25, 0.375, 0.75],
    )
    assert result.y_plane == pytest.approx(0.375)
    assert len(calls["a"]) == len(calls["b"]) == 3
    for (points_a, time_a), (points_b, time_b) in zip(calls["a"], calls["b"]):
        np.testing.assert_array_equal(points_a, points_b)
        np.testing.assert_allclose(points_a[:, 1], 0.375)
        assert time_a == time_b


def test_swirl_is_retained_but_projected_streamline_semantics_are_explicit():
    result = _sample()
    assert np.max(result.swirl_speed) > 0.0
    np.testing.assert_allclose(
        result.poloidal_speed,
        np.hypot(result.velocity[..., 0], result.velocity[..., 2]),
    )
    assert result.projected_streamline_semantics == "meridional_(u,w)_projection_only_not_true_3d_streamline"
    assert result.visualization_ready is False
    assert result.visual_correspondence_verified is False
    assert result.pde_validated is False
    assert result.openai_field_identified is False


def test_arrays_are_read_only_and_candidate_order_is_not_rewritten():
    result = _sample(
        candidates=[
            ("zeta", _linear_field(1.0), "zeta-v1"),
            ("alpha", _linear_field(1.2), "alpha-v1"),
        ]
    )
    assert result.candidate_labels == ("zeta", "alpha")
    assert result.candidate_provenance == ("zeta-v1", "alpha-v1")
    for array in (result.x, result.z, result.times, result.velocity, result.speed, result.panel_max_speed):
        assert not array.flags.writeable
    with pytest.raises(ValueError):
        result.speed[0, 0, 0, 0] = 0.0


def test_render_uses_one_fixed_comparison_payload(tmp_path: Path):
    pytest.importorskip("matplotlib")
    result = _sample(x=np.linspace(-1.0, 1.0, 17), z=np.linspace(-1.5, 1.5, 19))
    output = render_shared_frame_candidate_comparison(
        result,
        tmp_path / "comparison.png",
        streamline_density=0.7,
        linewidth=0.4,
        dpi=80,
    )
    assert output.exists()
    assert output.stat().st_size > 1000


@pytest.mark.parametrize(
    "candidates,kwargs,match",
    [
        ([("only", _linear_field(1.0), "only-v1")], {}, "at least two"),
        (
            [("dup", _linear_field(1.0), "a"), ("dup", _linear_field(2.0), "b")],
            {},
            "unique",
        ),
        (
            [("a", lambda p, t: np.zeros((len(p), 3)), "a"), ("b", _linear_field(1.0), "b")],
            {},
            "inactive",
        ),
        (
            [("a", lambda p, t: np.zeros((len(p), 2)), "a"), ("b", _linear_field(1.0), "b")],
            {},
            "returned shape",
        ),
        (
            [("a", lambda p, t: np.full((len(p), 3), np.nan), "a"), ("b", _linear_field(1.0), "b")],
            {},
            "nonfinite",
        ),
        (
            [("a", _linear_field(1.0), "a"), ("b", _linear_field(2.0), "b")],
            {"times": [0.5, 0.25]},
            "strictly increasing",
        ),
        (
            [("a", _linear_field(1.0), "a"), ("b", _linear_field(2.0), "b")],
            {"x": [-1.0, 0.0, 0.0, 1.0]},
            "strictly increasing",
        ),
        (
            [("a", _linear_field(1.0), ""), ("b", _linear_field(2.0), "b")],
            {},
            "provenance",
        ),
    ],
)
def test_fail_closed_inputs(candidates, kwargs, match):
    with pytest.raises(ValueError, match=match):
        _sample(candidates=candidates, **kwargs)


def test_render_rejects_non_png_and_bad_parameters(tmp_path: Path):
    result = _sample()
    with pytest.raises(ValueError, match=".png"):
        render_shared_frame_candidate_comparison(result, tmp_path / "bad.jpg")
    with pytest.raises(ValueError, match="positive"):
        render_shared_frame_candidate_comparison(result, tmp_path / "bad.png", streamline_density=0.0)
