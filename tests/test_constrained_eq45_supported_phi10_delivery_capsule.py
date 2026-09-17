import json

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_supported_phi10_delivery_capsule import (
    EXPECTED_BASE_SHA256,
    EXPECTED_SLOPE,
    build_trial_candidate,
    delivery_capsule,
    load_trial_recipe,
    write_trial_bundle,
)
from openai_ns_reconstruction.constrained_eq45_supported_temporal_mode import (
    Eq45SupportedAffineTemporalModeCandidate,
)


def test_phi10_temporal_capsule_keeps_export_visual_and_pde_states_independent():
    capsule = delivery_capsule()
    assert capsule["base_supported_sha256"] == EXPECTED_BASE_SHA256
    assert capsule["temporal_mode"] == {
        "family": "phi",
        "index": [1, 0],
        "slope": EXPECTED_SLOPE,
    }
    assert len(capsule["candidate_sha256"]) == 64
    assert len(capsule["recipe_sha256"]) == 64

    status = capsule["status"]
    assert status["velocity_export_ready"] is True
    assert status["visualization_candidate_only"] is True
    assert status["visualization_ready"] is False
    assert status["visual_correspondence_verified"] is False
    assert status["physical_support_validated"] is False
    assert status["pde_validated"] is False
    assert status["paper_exact"] is False
    assert status["openai_field_identified"] is False
    assert status["canonical_velocity_changed"] is False
    assert status["production_slope_promoted"] is False
    assert capsule["pending"]["whole_domain_morphology_acceptance"] == "PR146_open_unconsumed"
    assert capsule["pde_evidence"]["uniform_generalization"] is False
    assert capsule["pde_evidence"]["formal_pde_gate_assessed"] is False


def test_phi10_temporal_bundle_roundtrips_public_velocity_and_grid(tmp_path):
    paths = write_trial_bundle(tmp_path / "bundle")
    capsule = json.loads(paths["capsule"].read_text(encoding="utf-8"))
    loaded = Eq45SupportedAffineTemporalModeCandidate.load_json(paths["candidate"])
    fresh = build_trial_candidate()

    assert loaded.sha256 == capsule["candidate_sha256"] == fresh.sha256
    assert loaded.base_sha256 == EXPECTED_BASE_SHA256

    points = np.array(
        [
            [0.37, 0.11, -0.22],
            [1.25, 0.00, 0.35],
            [1.72, 0.20, 0.40],
        ],
        dtype=float,
    )
    times = np.array([0.25, 0.50, 0.75], dtype=float)
    np.testing.assert_array_equal(
        loaded.at_points(points, times),
        fresh.at_points(points, times),
    )
    np.testing.assert_array_equal(
        loaded.at_points(points, 0.50),
        loaded.base.at_points(points, 0.50),
    )

    grid = loaded.grid(
        np.array([-0.5, 0.5]),
        np.array([-0.4, 0.4]),
        np.array([-0.3, 0.3]),
        times,
    )
    assert grid.shape == (3, 2, 2, 2, 3)
    assert np.all(np.isfinite(grid))


def test_phi10_temporal_recipe_fails_closed_on_promotion_or_identity_drift(tmp_path):
    recipe = load_trial_recipe()

    promoted = json.loads(json.dumps(recipe))
    promoted["status"]["visualization_ready"] = True
    promoted_path = tmp_path / "promoted.json"
    promoted_path.write_text(json.dumps(promoted), encoding="utf-8")
    with pytest.raises(ValueError, match="truth-boundary"):
        load_trial_recipe(promoted_path)

    changed = json.loads(json.dumps(recipe))
    changed["base_supported_sha256"] = "0" * 64
    changed_path = tmp_path / "changed.json"
    changed_path.write_text(json.dumps(changed), encoding="utf-8")
    with pytest.raises(ValueError, match="base supported identity"):
        load_trial_recipe(changed_path)
