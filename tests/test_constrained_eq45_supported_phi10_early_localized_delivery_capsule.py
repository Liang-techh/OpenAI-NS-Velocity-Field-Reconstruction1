import copy
import json

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_candidate import Eq45VelocityCandidate
from openai_ns_reconstruction.constrained_eq45_supported_candidate import (
    Eq45SupportedVelocityCandidate,
)
from openai_ns_reconstruction.constrained_eq45_supported_phi10_early_localized_delivery_capsule import (
    EXPECTED_BASE_SHA256,
    EXPECTED_CANDIDATE_SHA256,
    EXPECTED_COEFFICIENTS,
    EXPECTED_TIMES,
    build_candidate,
    delivery_capsule,
    load_recipe,
    write_delivery_bundle,
)
from openai_ns_reconstruction.constrained_eq45_supported_phi10_early_localized_temporal_mode import (
    Eq45SupportedPhi10EarlyLocalizedTemporalCandidate,
)


def _probes():
    return np.asarray(
        [
            [0.37, 0.11, -0.22],
            [0.95, -0.35, 0.48],
            [1.45, 0.15, -0.62],
            [1.72, 0.00, 0.30],
        ],
        dtype=float,
    )


def test_early_localized_capsule_binds_exact_public_candidate():
    recipe = load_recipe()
    candidate = build_candidate(recipe)
    capsule = delivery_capsule()

    assert candidate.base_sha256 == EXPECTED_BASE_SHA256
    assert candidate.sha256 == EXPECTED_CANDIDATE_SHA256
    assert capsule["candidate_sha256"] == EXPECTED_CANDIDATE_SHA256
    assert capsule["base_supported_sha256"] == EXPECTED_BASE_SHA256
    assert tuple(capsule["delivery"]["reference_times"]) == EXPECTED_TIMES
    assert tuple(float(candidate.coefficient_at(t)) for t in EXPECTED_TIMES) == pytest.approx(
        EXPECTED_COEFFICIENTS, rel=0.0, abs=1.0e-14
    )

    status = capsule["status"]
    assert status["velocity_export_ready"] is True
    assert status["visualization_candidate_only"] is True
    assert status["visualization_ready"] is False
    assert status["visual_correspondence_verified"] is False
    assert status["pde_validated"] is False
    assert status["paper_exact"] is False
    assert status["openai_field_identified"] is False
    assert status["blowup_proved"] is False

    pde = capsule["pde_evidence"]
    assert pde["fresh_seed_ci_success"] is True
    assert pde["fresh_seed_numeric_summary_bound_here"] is False
    assert pde["formal_pde_gate_assessed"] is False
    assert pde["pde_validated"] is False
    assert pde["early_localized_zero_force_rms"] > pde["static_zero_force_rms"]


def test_early_localized_capsule_replays_endpoint_public_velocity_exactly():
    candidate = build_candidate()
    static = Eq45SupportedVelocityCandidate(parent=Eq45VelocityCandidate.seed())
    points = _probes()

    midpoint = candidate.at_points(points, 0.5)
    late = candidate.at_points(points, 0.75)
    assert np.array_equal(midpoint, static.at_points(points, 0.5))
    assert np.array_equal(late, static.at_points(points, 0.75))

    early = candidate.at_points(points, 0.25)
    assert np.all(np.isfinite(early))
    assert not np.array_equal(early, static.at_points(points, 0.25))


def test_early_localized_delivery_bundle_roundtrips_public_uvws_and_grid(tmp_path):
    original = build_candidate()
    paths = write_delivery_bundle(tmp_path)
    loaded = Eq45SupportedPhi10EarlyLocalizedTemporalCandidate.load_json(paths["candidate"])
    capsule = json.loads(paths["capsule"].read_text(encoding="utf-8"))

    assert loaded.sha256 == original.sha256 == capsule["candidate_sha256"]
    points = _probes()
    times = np.asarray([0.25, 0.5, 0.625, 0.75], dtype=float)
    assert np.array_equal(loaded.at_points(points, times), original.at_points(points, times))

    grid = loaded.grid(
        np.asarray([-0.4, 0.4]),
        np.asarray([-0.3, 0.3]),
        np.asarray([-0.2, 0.2]),
        np.asarray(EXPECTED_TIMES),
    )
    assert grid.shape == (4, 2, 2, 2, 3)
    assert np.all(np.isfinite(grid))


def test_early_localized_recipe_fails_closed_on_identity_or_claim_promotion(tmp_path):
    recipe = load_recipe()

    cases = []
    mutated = copy.deepcopy(recipe)
    mutated["candidate_sha256"] = "0" * 64
    cases.append(mutated)

    mutated = copy.deepcopy(recipe)
    mutated["temporal_mode"]["early_delta"] = -1.3
    cases.append(mutated)

    mutated = copy.deepcopy(recipe)
    mutated["status"]["visualization_ready"] = True
    cases.append(mutated)

    mutated = copy.deepcopy(recipe)
    mutated["status"]["pde_validated"] = True
    cases.append(mutated)

    mutated = copy.deepcopy(recipe)
    mutated["selection_evidence"]["visual_correspondence_established"] = True
    cases.append(mutated)

    mutated = copy.deepcopy(recipe)
    mutated["pde_evidence"]["formal_pde_gate_assessed"] = True
    cases.append(mutated)

    for index, payload in enumerate(cases):
        path = tmp_path / f"mutated-{index}.json"
        path.write_text(json.dumps(payload), encoding="utf-8")
        with pytest.raises(ValueError):
            load_recipe(path)
