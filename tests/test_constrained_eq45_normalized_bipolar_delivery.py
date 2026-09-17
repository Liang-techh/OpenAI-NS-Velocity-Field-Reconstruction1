from pathlib import Path

import numpy as np

from openai_ns_reconstruction.constrained_eq45_supported_candidate import (
    Eq45SupportedVelocityCandidate,
)
from openai_ns_reconstruction.eq45_normalized_bipolar_delivery import (
    ENERGY_NORMALIZATION_SCALE,
    NORMALIZED_BIPOLAR_SHA256,
    SOURCE_BIPOLAR_SHA256,
    Eq45NormalizedBipolarDeliveryField,
    components,
    default_field,
    normalized_bipolar_candidate,
    u,
    v,
    velocity,
    w,
)


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "artifacts" / "bipolar_energy" / "normalized_candidate.json"


def test_normalized_bipolar_named_entrypoint_replays_committed_identity():
    expected = Eq45SupportedVelocityCandidate.load_json(ARTIFACT)
    candidate = normalized_bipolar_candidate()
    field = default_field()

    assert field.sha256 == NORMALIZED_BIPOLAR_SHA256
    assert candidate.sha256 == NORMALIZED_BIPOLAR_SHA256
    assert candidate.to_dict() == expected.to_dict()

    point = (0.1, 0.0, 0.1, 0.5)
    expected_value = expected.velocity_xyz(*point)
    actual = velocity(*point)
    np.testing.assert_allclose(actual, expected_value, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(components(*point), expected_value, rtol=0.0, atol=0.0)
    assert u(*point) == expected_value[0]
    assert v(*point) == expected_value[1]
    assert w(*point) == expected_value[2]


def test_normalized_bipolar_delivery_keeps_opposite_axial_core_and_exports_grid(tmp_path):
    field = Eq45NormalizedBipolarDeliveryField()
    plus = field.velocity(0.1, 0.0, 0.1, 0.5)
    minus = field.velocity(0.1, 0.0, -0.1, 0.5)
    assert plus[2] > 0.0
    assert minus[2] < 0.0

    axes = np.array([-0.2, 0.0, 0.2])
    times = np.array([0.25, 0.5, 0.75])
    grid = field.grid(axes, axes, axes, times)
    assert grid.shape == (3, 3, 3, 3, 3)
    assert np.isfinite(grid).all()
    assert np.max(np.abs(grid)) > 0.0

    saved = tmp_path / "normalized_bipolar.json"
    field.save_candidate(saved)
    reloaded = Eq45NormalizedBipolarDeliveryField.load_candidate(saved)
    assert reloaded.sha256 == field.sha256
    np.testing.assert_allclose(reloaded.grid(axes, axes, axes, times), grid, rtol=0.0, atol=0.0)


def test_normalized_bipolar_metadata_keeps_truth_states_separate():
    field = default_field()
    metadata = field.metadata()
    truth = metadata["truth_boundary"]

    assert metadata["candidate_sha256"] == NORMALIZED_BIPOLAR_SHA256
    assert metadata["source_bipolar_sha256"] == SOURCE_BIPOLAR_SHA256
    assert metadata["energy_normalization_scale"] == ENERGY_NORMALIZATION_SCALE
    assert metadata["candidate_role"] == "experimental_source_aligned_energy_normalized"
    assert metadata["selection_status"] == "not_canonical_not_visualization_selected"
    assert metadata["entrypoint"] == (
        "openai_ns_reconstruction.eq45_normalized_bipolar_delivery:velocity"
    )
    assert truth["velocity_export_ready"] is True
    assert truth["physical_support_connection_implemented"] is True
    assert truth["physical_support_validated"] is False
    assert truth["visualization_ready"] is False
    assert truth["visual_correspondence_verified"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False
