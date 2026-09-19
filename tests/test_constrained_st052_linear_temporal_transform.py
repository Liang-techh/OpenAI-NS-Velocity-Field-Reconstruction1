from __future__ import annotations

import json

import numpy as np
import pytest

from openai_ns_reconstruction.st052_linear_temporal_transform import (
    DEFAULT_TEMPORAL_SPEC,
    STATIC_TRANSFORM_SPEC_SHA256,
    St052LinearTemporalAdapter,
    St052LinearTemporalTransformSpec,
    activation,
    effective_static_spec,
    temporal_transformed_velocity_points,
)
from openai_ns_reconstruction.st052_local_swirl_transform import (
    DEFAULT_SPEC as STATIC_SPEC,
    redistributed_control,
    transformed_velocity_points,
)


def _base_velocity(points: np.ndarray, time: float) -> np.ndarray:
    points = np.asarray(points, dtype=float)
    x, y, z = points.T
    scale = 1.0 + 0.1 * float(time)
    return np.column_stack((-scale * y, scale * x, 0.2 * z))


def _points() -> np.ndarray:
    return np.array(
        [
            [0.65, 0.10, 0.85],
            [0.90, -0.20, -0.85],
            [1.15, 0.15, 1.15],
            [0.72, -0.33, -1.15],
            [0.40, 0.20, 0.30],
        ],
        dtype=float,
    )


def test_static_transform_binding_is_exact() -> None:
    assert STATIC_SPEC.sha256() == STATIC_TRANSFORM_SPEC_SHA256
    DEFAULT_TEMPORAL_SPEC.validate_static_binding()
    assert DEFAULT_TEMPORAL_SPEC.source_temporal_pr == 587
    assert (
        DEFAULT_TEMPORAL_SPEC.source_temporal_head
        == "0b93819095f6c8576a7571bdc2d2fbef4154944d"
    )


def test_activation_and_endpoint_representation_identities() -> None:
    pts = _points()
    assert activation(0.25) == 0.0
    assert activation(0.50) == 0.5
    assert activation(0.75) == 1.0

    early = temporal_transformed_velocity_points(_base_velocity, pts, 0.25)
    control = redistributed_control(_base_velocity, pts, 0.25, STATIC_SPEC)
    np.testing.assert_array_equal(early, control)

    late = temporal_transformed_velocity_points(_base_velocity, pts, 0.75)
    static = transformed_velocity_points(_base_velocity, pts, 0.75, STATIC_SPEC)
    np.testing.assert_array_equal(late, static)

    middle = temporal_transformed_velocity_points(_base_velocity, pts, 0.50)
    middle_control = redistributed_control(_base_velocity, pts, 0.50, STATIC_SPEC)
    middle_static = transformed_velocity_points(_base_velocity, pts, 0.50, STATIC_SPEC)
    assert np.max(np.abs(middle - middle_control)) > 1.0e-8
    assert np.max(np.abs(middle - middle_static)) > 1.0e-8


def test_effective_static_spec_scales_only_temporal_pair() -> None:
    spec = effective_static_spec(0.50)
    assert spec.taper_tau == pytest.approx(0.025)
    assert spec.shoulder_beta == pytest.approx(0.5 * STATIC_SPEC.shoulder_beta)
    assert spec.redistribution_gain == STATIC_SPEC.redistribution_gain
    assert spec.redistribution_alpha == STATIC_SPEC.redistribution_alpha
    assert spec.redistribution_scale == STATIC_SPEC.redistribution_scale
    assert spec.post_transform_common_scale == 1.0


def test_temporal_spec_round_trip_and_checksum_failure(tmp_path) -> None:
    path = tmp_path / "temporal_spec.json"
    DEFAULT_TEMPORAL_SPEC.save(path)
    loaded = St052LinearTemporalTransformSpec.load(path)
    assert loaded == DEFAULT_TEMPORAL_SPEC
    assert loaded.sha256() == DEFAULT_TEMPORAL_SPEC.sha256()

    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["spec"]["activation_slope"] = 1.5
    path.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(ValueError, match="checksum"):
        St052LinearTemporalTransformSpec.load(path)


def test_temporal_spec_rejects_static_identity_laundering() -> None:
    bad = St052LinearTemporalTransformSpec(static_child_pr=560)
    with pytest.raises(ValueError, match="child PR"):
        bad.validate_static_binding()


def test_adapter_broadcasts_and_rejects_parent_metadata_drift() -> None:
    adapter = St052LinearTemporalAdapter(
        _base_velocity,
        parent_candidate_id="ST052-M",
        parent_source_head="b3b8bfdbe1077f9ec967d158602951997d81e17d",
    )
    out = adapter.velocity(
        np.array([0.65, 1.15]),
        np.array([0.10, 0.15]),
        np.array([0.85, 1.15]),
        np.array([0.25, 0.75]),
    )
    assert out.shape == (2, 3)
    assert np.isfinite(out).all()

    with pytest.raises(ValueError, match="parent candidate"):
        St052LinearTemporalAdapter(
            _base_velocity,
            parent_candidate_id="not-ST052-M",
            parent_source_head="b3b8bfdbe1077f9ec967d158602951997d81e17d",
        )
    with pytest.raises(ValueError, match="parent source"):
        St052LinearTemporalAdapter(
            _base_velocity,
            parent_candidate_id="ST052-M",
            parent_source_head="deadbeef",
        )


def test_invalid_time_and_nonfinite_inputs_fail_closed() -> None:
    pts = _points()
    with pytest.raises(ValueError, match="outside"):
        temporal_transformed_velocity_points(_base_velocity, pts, 0.20)

    adapter = St052LinearTemporalAdapter(
        _base_velocity,
        parent_candidate_id="ST052-M",
        parent_source_head="b3b8bfdbe1077f9ec967d158602951997d81e17d",
    )
    with pytest.raises(ValueError, match="finite"):
        adapter.velocity(np.nan, 0.0, 0.0, 0.5)
