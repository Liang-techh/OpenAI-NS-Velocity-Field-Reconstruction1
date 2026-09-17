import json
from pathlib import Path

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_supported_artifact import (
    CANONICAL_PARENT_SHA256,
    SUPPORTED_CANDIDATE_SHA256,
    default_candidate_path,
    load_checked_supported_candidate,
    write_checked_supported_candidate,
)


def test_supported_artifact_replays_public_velocity_and_grid():
    candidate = load_checked_supported_candidate()

    assert candidate.parent_sha256 == CANONICAL_PARENT_SHA256
    assert candidate.sha256 == SUPPORTED_CANDIDATE_SHA256

    points = np.array(
        [
            [0.3, 0.2, -0.4],
            [1.75, 0.0, 0.3],
            [2.0, 0.0, 0.0],
        ]
    )
    values = candidate.at_points(points, np.array([0.25, 0.5, 0.75]))
    assert values.shape == (3, 3)
    assert np.all(np.isfinite(values))
    np.testing.assert_allclose(values[-1], 0.0, rtol=0.0, atol=0.0)

    axis = np.array([-2.0, 0.0, 2.0])
    grid = candidate.grid(axis, axis, axis, np.array([0.25, 0.75]))
    assert grid.shape == (2, 3, 3, 3, 3)
    assert np.all(np.isfinite(grid))


def test_supported_artifact_regeneration_is_byte_identical(tmp_path):
    regenerated = tmp_path / "supported.json"
    written = write_checked_supported_candidate(regenerated)

    assert written.sha256 == SUPPORTED_CANDIDATE_SHA256
    assert regenerated.read_bytes() == default_candidate_path().read_bytes()


def test_supported_artifact_rejects_valid_but_noncanonical_taper(tmp_path):
    payload = json.loads(default_candidate_path().read_text(encoding="utf-8"))
    payload["physical_taper"]["radial_plateau_q"] = 0.65
    mutated = tmp_path / "mutated.json"
    mutated.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    with pytest.raises(ValueError, match="SHA drifted"):
        load_checked_supported_candidate(mutated)
