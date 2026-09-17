import copy
import json

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_candidate import Eq45VelocityCandidate
from openai_ns_reconstruction.constrained_eq45_supported_candidate import (
    Eq45SupportedVelocityCandidate,
)
from openai_ns_reconstruction.constrained_eq45_supported_phi10_cubic_delivery_capsule import (
    EXPECTED_BASE_SHA256,
    EXPECTED_CANDIDATE_SHA256,
    EXPECTED_COEFFICIENTS,
    EXPECTED_TIMES,
    build_candidate,
    delivery_capsule,
    load_recipe,
    write_delivery_bundle,
)
from openai_ns_reconstruction.constrained_eq45_supported_phi10_cubic_temporal_mode import (
    Eq45SupportedPhi10CubicLocalizedTemporalCandidate,
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


def test_cubic_capsule_binds_exact_public_candidate_and_convergence_evidence():
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

    temporal = capsule["temporal_derivative_evidence"]
    assert temporal["probe_seed"] == 914509
    assert temporal["second_order_steps"] == [0.01, 0.005, 0.0025]
    assert temporal["normalized_rms_errors"][-1] < 1.1e-4
    assert min(temporal["observed_rms_orders"]) > 1.9
    assert temporal["snapshot_identity_not_derivative_identity"] is True
    assert temporal["formal_pde_gate_assessed"] is False
    assert temporal["pde_validated"] is False

    sibling = capsule["sibling_evidence"]
    assert sibling["restricted_force_crosscheck_pr"] == 160
    assert sibling["ancestry_status"] == "open_unconsumed_sibling_evidence"
    assert sibling["numeric_summary_bound_here"] is False

    status = capsule["status"]
    assert status["velocity_export_ready"] is True
    assert status["visualization_candidate_only"] is True
    assert status["visualization_ready"] is False
    assert status["visual_correspondence_verified"] is False
    assert status["pde_validated"] is False
    assert status["paper_exact"] is False
    assert status["openai_field_identified"] is False
    assert status["blowup_proved"] is False


def test_cubic_capsule_replays_static_keyframes_and_nontrivial_off_keyframes():
    candidate = build_candidate()
    static = Eq45SupportedVelocityCandidate(parent=Eq45VelocityCandidate.seed())
    points = _probes()

    for time in (0.5, 0.625, 0.75):
        assert np.array_equal(
            candidate.at_points(points, time), static.at_points(points, time)
        )

    early = candidate.at_points(points, 0.25)
    off_keyframe = candidate.at_points(points, 0.375)
    assert np.all(np.isfinite(early))
    assert np.all(np.isfinite(off_keyframe))
    assert not np.array_equal(early, static.at_points(points, 0.25))
    assert not np.array_equal(off_keyframe, static.at_points(points, 0.375))


def test_cubic_delivery_bundle_roundtrips_public_uvws_and_grid(tmp_path):
    original = build_candidate()
    paths = write_delivery_bundle(tmp_path)
    loaded = Eq45SupportedPhi10CubicLocalizedTemporalCandidate.load_json(paths["candidate"])
    capsule = json.loads(paths["capsule"].read_text(encoding="utf-8"))

    assert loaded.sha256 == original.sha256 == capsule["candidate_sha256"]
    points = np.asarray(
        [
            [0.37, 0.11, -0.22],
            [0.95, -0.35, 0.48],
            [1.45, 0.15, -0.62],
            [1.72, 0.00, 0.30],
            [0.44, -0.52, 1.68],
            [1.70, 0.17, -1.72],
        ],
        dtype=float,
    )
    times = np.asarray(EXPECTED_TIMES, dtype=float)
    assert np.array_equal(loaded.at_points(points, times), original.at_points(points, times))

    grid = loaded.grid(
        np.asarray([-0.4, 0.4]),
        np.asarray([-0.3, 0.3]),
        np.asarray([-0.2, 0.2]),
        np.asarray(EXPECTED_TIMES),
    )
    assert grid.shape == (6, 2, 2, 2, 3)
    assert np.all(np.isfinite(grid))


def test_cubic_recipe_fails_closed_on_identity_evidence_or_claim_promotion(tmp_path):
    recipe = load_recipe()
    cases = []

    mutated = copy.deepcopy(recipe)
    mutated["candidate_sha256"] = "0" * 64
    cases.append(mutated)

    mutated = copy.deepcopy(recipe)
    mutated["temporal_mode"]["early_delta"] = -1.3
    cases.append(mutated)

    mutated = copy.deepcopy(recipe)
    mutated["temporal_derivative_evidence"]["probe_seed"] = 914510
    cases.append(mutated)

    mutated = copy.deepcopy(recipe)
    mutated["temporal_derivative_evidence"]["formal_pde_gate_assessed"] = True
    cases.append(mutated)

    mutated = copy.deepcopy(recipe)
    mutated["sibling_evidence"]["numeric_summary_bound_here"] = True
    cases.append(mutated)

    mutated = copy.deepcopy(recipe)
    mutated["representation_evidence"]["visual_correspondence_established"] = True
    cases.append(mutated)

    mutated = copy.deepcopy(recipe)
    mutated["status"]["visualization_ready"] = True
    cases.append(mutated)

    mutated = copy.deepcopy(recipe)
    mutated["status"]["pde_validated"] = True
    cases.append(mutated)

    for index, payload in enumerate(cases):
        path = tmp_path / f"mutated-{index}.json"
        path.write_text(json.dumps(payload), encoding="utf-8")
        with pytest.raises(ValueError):
            load_recipe(path)
