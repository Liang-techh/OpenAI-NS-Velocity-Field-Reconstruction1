import json

import numpy as np
import pytest

from openai_ns_reconstruction.st052_local_swirl_transform import (
    DEFAULT_SPEC,
    St052LocalSwirlAdapter,
    St052LocalSwirlTransformSpec,
    redistributed_control,
    transformed_velocity_points,
)


def base_velocity(points, time):
    p = np.asarray(points, dtype=float)
    x, y, z = p.T
    inside = (np.hypot(x, y) < 2.0) & (np.abs(z) < 2.0)
    out = np.column_stack(
        (
            -0.18 * x + 0.03 * z,
            0.31 * x + 0.07 * y,
            0.16 * z - 0.02 * x,
        )
    )
    return np.where(inside[:, None], (1.0 + 0.1 * (float(time) - 0.25)) * out, 0.0)


def adapter():
    return St052LocalSwirlAdapter(
        base_velocity,
        parent_candidate_id="ST052-M",
        parent_source_head=DEFAULT_SPEC.parent_source_head,
    )


def test_transform_spec_roundtrip_and_checksum(tmp_path):
    path = tmp_path / "transform.json"
    DEFAULT_SPEC.save(path)
    loaded = St052LocalSwirlTransformSpec.load(path)
    assert loaded == DEFAULT_SPEC
    assert loaded.sha256() == DEFAULT_SPEC.sha256()

    obj = json.loads(path.read_text())
    obj["spec"]["shoulder_beta"] += 1.0e-4
    path.write_text(json.dumps(obj))
    with pytest.raises(ValueError, match="checksum"):
        St052LocalSwirlTransformSpec.load(path)


def test_adapter_rejects_parent_identity_drift():
    with pytest.raises(ValueError, match="candidate identity"):
        St052LocalSwirlAdapter(
            base_velocity,
            parent_candidate_id="not-ST052-M",
            parent_source_head=DEFAULT_SPEC.parent_source_head,
        )
    with pytest.raises(ValueError, match="source head"):
        St052LocalSwirlAdapter(
            base_velocity,
            parent_candidate_id="ST052-M",
            parent_source_head="0" * 40,
        )


def test_central_band_reduces_to_frozen_redistribution():
    points = np.array([[0.6, 0.0, -0.3], [0.9, 0.0, 0.0], [1.2, 0.0, 0.3]])
    expected = redistributed_control(base_velocity, points, 0.5)
    actual = transformed_velocity_points(base_velocity, points, 0.5)
    np.testing.assert_allclose(actual, expected, rtol=0.0, atol=1.0e-14)


def test_transform_is_active_in_shoulder_and_tip_regions():
    points = np.array([[0.9, 0.0, 0.85], [0.9, 0.0, 1.15]])
    control = redistributed_control(base_velocity, points, 0.5)
    child = transformed_velocity_points(base_velocity, points, 0.5)
    assert np.linalg.norm(child[0] - control[0]) > 1.0e-6
    assert np.linalg.norm(child[1] - control[1]) > 1.0e-6


def test_velocity_xyz_broadcast_and_time_guard():
    field = adapter()
    x = np.array([0.6, 0.9, 1.2])[:, None]
    z = np.array([-0.3, 0.85])[None, :]
    out = field.velocity(x, 0.0, z, 0.5)
    assert out.shape == (3, 2, 3)
    for i in range(3):
        for j in range(2):
            expected = field.points(np.array([[x[i, 0], 0.0, z[0, j]]]), 0.5)[0]
            np.testing.assert_allclose(out[i, j], expected, rtol=0.0, atol=1.0e-14)
    with pytest.raises(ValueError, match="time outside"):
        field.velocity(0.6, 0.0, 0.0, 0.9)


def test_axis_and_exterior_are_finite():
    field = adapter()
    points = np.array([[0.0, 0.0, 0.85], [2.1, 0.0, 0.0], [0.0, 0.0, 2.1]])
    out = field.points(points, 0.5)
    assert np.isfinite(out).all()
    np.testing.assert_allclose(out[1:], 0.0, rtol=0.0, atol=0.0)
