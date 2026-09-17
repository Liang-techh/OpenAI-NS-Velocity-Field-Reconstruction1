from pathlib import Path

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_candidate import Eq45VelocityCandidate
from openai_ns_reconstruction.constrained_eq45_supported_candidate import (
    Eq45SupportedVelocityCandidate,
)
from openai_ns_reconstruction.eq45_delivery import default_field as untapered_default_field
from openai_ns_reconstruction.eq45_supported_delivery import (
    Eq45SupportedDeliveryField,
    components,
    default_field,
    u,
    v,
    velocity,
    w,
)


ROOT = Path(__file__).resolve().parents[1]
SEED = ROOT / "artifacts" / "constrained" / "eq45_velocity_candidate_seed.json"


def _expected_child() -> Eq45SupportedVelocityCandidate:
    return Eq45SupportedVelocityCandidate(parent=Eq45VelocityCandidate.load_json(SEED))


def test_supported_named_entrypoint_matches_child_and_preserves_distinct_identity():
    expected_child = _expected_child()
    field = default_field()

    assert field.sha256 == expected_child.sha256
    assert field.sha256 != untapered_default_field().sha256

    expected = expected_child.velocity_xyz(0.1, -0.2, 0.15, 0.5)
    actual = velocity(0.1, -0.2, 0.15, 0.5)
    np.testing.assert_allclose(actual, expected, rtol=0.0, atol=0.0)
    assert actual.shape == (3,)

    comp = components(0.1, -0.2, 0.15, 0.5)
    assert all(isinstance(value, float) for value in comp)
    np.testing.assert_allclose(comp, expected, rtol=0.0, atol=0.0)
    assert u(0.1, -0.2, 0.15, 0.5) == comp[0]
    assert v(0.1, -0.2, 0.15, 0.5) == comp[1]
    assert w(0.1, -0.2, 0.15, 0.5) == comp[2]


def test_supported_named_entrypoint_broadcast_grid_support_and_roundtrip(tmp_path):
    field = Eq45SupportedDeliveryField(candidate=_expected_child())
    x = np.array([0.0, 0.2])[:, None]
    z = np.array([-0.1, 0.1])[None, :]
    values = field.velocity(x, 0.05, z, 0.5)
    assert values.shape == (2, 2, 3)

    points = np.stack(np.broadcast_arrays(x, np.asarray(0.05), z), axis=-1)
    np.testing.assert_allclose(values, field.at_points(points, 0.5), rtol=0.0, atol=0.0)

    np.testing.assert_allclose(field.velocity(2.0, 0.0, 0.0, 0.5), 0.0, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(field.velocity(0.0, 0.0, 2.0, 0.5), 0.0, rtol=0.0, atol=0.0)

    axes = np.array([-0.2, 0.0, 0.2])
    times = np.array([0.25, 0.5, 0.75])
    grid = field.grid(axes, axes, axes, times)
    assert grid.shape == (3, 3, 3, 3, 3)

    saved = tmp_path / "eq45_supported.json"
    field.save_candidate(saved)
    reloaded = Eq45SupportedDeliveryField.load_candidate(saved)
    assert reloaded.sha256 == field.sha256
    np.testing.assert_allclose(reloaded.grid(axes, axes, axes, times), grid, rtol=0.0, atol=0.0)


def test_supported_named_metadata_keeps_child_truth_boundary_and_fails_closed():
    field = Eq45SupportedDeliveryField(candidate=_expected_child())
    metadata = field.metadata()
    truth = metadata["truth_boundary"]

    assert metadata["family"] == "eq45_supported_velocity_candidate_v1"
    assert metadata["entrypoint"] == "openai_ns_reconstruction.eq45_supported_delivery:velocity"
    assert metadata["candidate_sha256"] == field.sha256
    assert metadata["parent_sha256"] == field.candidate.parent_sha256
    assert metadata["grid_layout"] == ["time", "x", "y", "z", "component"]
    assert truth["velocity_export_ready"] is True
    assert truth["physical_support_connection_implemented"] is True
    assert truth["physical_support_validated"] is False
    assert truth["visualization_ready"] is False
    assert truth["visual_correspondence_verified"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False

    with pytest.raises(ValueError, match="not both"):
        Eq45SupportedDeliveryField(field.candidate, candidate_path=SEED)
    with pytest.raises(TypeError, match="Eq45SupportedVelocityCandidate"):
        Eq45SupportedDeliveryField(object())
    with pytest.raises(ValueError, match="declared interval"):
        field.velocity(0.0, 0.0, 0.0, 0.9)
