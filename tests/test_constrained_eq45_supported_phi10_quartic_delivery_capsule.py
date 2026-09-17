import copy
import json

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_candidate import Eq45VelocityCandidate
from openai_ns_reconstruction.constrained_eq45_supported_candidate import (
    Eq45SupportedVelocityCandidate,
)
from openai_ns_reconstruction.constrained_eq45_supported_phi10_quartic_delivery_capsule import (
    EXPECTED_BASE_SHA256,
    EXPECTED_CANDIDATE_SHA256,
    EXPECTED_COEFFICIENTS,
    EXPECTED_NULLSPACE_COEFFICIENT,
    EXPECTED_TIMES,
    build_candidate,
    delivery_capsule,
    load_recipe,
    write_delivery_bundle,
)
from openai_ns_reconstruction.constrained_eq45_supported_phi10_quartic_temporal_mode import (
    Eq45SupportedPhi10QuarticDerivativeBalancedTemporalCandidate,
)


def _probes():
    return np.asarray(
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


def test_quartic_capsule_binds_exact_candidate_and_force_evidence():
    recipe = load_recipe()
    candidate = build_candidate(recipe)
    capsule = delivery_capsule()

    assert candidate.base_sha256 == EXPECTED_BASE_SHA256
    assert candidate.sha256 == EXPECTED_CANDIDATE_SHA256
    assert candidate.nullspace_coefficient == EXPECTED_NULLSPACE_COEFFICIENT
    assert capsule["candidate_sha256"] == EXPECTED_CANDIDATE_SHA256
    assert capsule["base_supported_sha256"] == EXPECTED_BASE_SHA256
    assert tuple(capsule["delivery"]["reference_times"]) == EXPECTED_TIMES
    assert tuple(float(candidate.coefficient_at(t)) for t in EXPECTED_TIMES) == pytest.approx(
        EXPECTED_COEFFICIENTS, rel=0.0, abs=1.0e-14
    )

    force = capsule["restricted_force_evidence"]
    assert force["audit_pr"] == 169
    assert force["spatial_derivative_steps"] == [0.02, 0.01, 0.005]
    assert force["frozen_force_coefficients"]["c"] == pytest.approx(
        0.07770084654400304, rel=0.0, abs=1.0e-15
    )
    assert force["heldout_rms_before"][-1] == pytest.approx(
        3.7016892575139453, rel=0.0, abs=1.0e-14
    )
    assert force["heldout_rms_after"][-1] == pytest.approx(
        3.700414843717181, rel=0.0, abs=1.0e-14
    )
    assert force["formal_pde_gate_assessed"] is False
    assert force["pde_validated"] is False

    siblings = capsule["sibling_evidence"]
    assert siblings["seed_generalization_status"] == "claimed_unconsumed"
    assert siblings["late_morphology_pr"] == 170
    assert siblings["late_morphology_status"] == "open_unconsumed_calibration_pending"
    assert siblings["numeric_summary_bound_here"] is False

    status = capsule["status"]
    assert status["velocity_export_ready"] is True
    assert status["visualization_candidate_only"] is True
    assert status["visualization_ready"] is False
    assert status["visual_correspondence_verified"] is False
    assert status["pde_validated"] is False
    assert status["paper_exact"] is False
    assert status["openai_field_identified"] is False
    assert status["blowup_proved"] is False


def test_quartic_capsule_replays_static_keyframes_and_off_keyframes():
    candidate = build_candidate()
    static = Eq45SupportedVelocityCandidate(parent=Eq45VelocityCandidate.seed())
    points = _probes()

    for time in (0.5, 0.625, 0.75):
        assert np.array_equal(
            candidate.at_points(points, time), static.at_points(points, time)
        )

    for time in (0.25, 0.375, 0.5625, 0.6875):
        values = candidate.at_points(points, time)
        assert np.all(np.isfinite(values))
        assert not np.array_equal(values, static.at_points(points, time))


def test_quartic_delivery_bundle_roundtrips_public_uvws_and_grid(tmp_path):
    original = build_candidate()
    paths = write_delivery_bundle(tmp_path)
    loaded = Eq45SupportedPhi10QuarticDerivativeBalancedTemporalCandidate.load_json(
        paths["candidate"]
    )
    capsule = json.loads(paths["capsule"].read_text(encoding="utf-8"))

    assert loaded.sha256 == original.sha256 == capsule["candidate_sha256"]
    points = _probes()
    times = np.asarray(EXPECTED_TIMES[: points.shape[0]], dtype=float)
    assert np.array_equal(loaded.at_points(points, times), original.at_points(points, times))

    grid = loaded.grid(
        np.asarray([-0.4, 0.4]),
        np.asarray([-0.3, 0.3]),
        np.asarray([-0.2, 0.2]),
        np.asarray(EXPECTED_TIMES),
    )
    assert grid.shape == (7, 2, 2, 2, 3)
    assert np.all(np.isfinite(grid))


def test_quartic_recipe_fails_closed_on_identity_evidence_or_claim_promotion(tmp_path):
    recipe = load_recipe()
    cases = []

    mutated = copy.deepcopy(recipe)
    mutated["candidate_sha256"] = "0" * 64
    cases.append(mutated)

    mutated = copy.deepcopy(recipe)
    mutated["temporal_mode"]["nullspace_coefficient"] = -0.28
    cases.append(mutated)

    mutated = copy.deepcopy(recipe)
    mutated["restricted_force_evidence"]["audit_exact_head"] = "0" * 40
    cases.append(mutated)

    mutated = copy.deepcopy(recipe)
    mutated["restricted_force_evidence"]["formal_pde_gate_assessed"] = True
    cases.append(mutated)

    mutated = copy.deepcopy(recipe)
    mutated["sibling_evidence"]["numeric_summary_bound_here"] = True
    cases.append(mutated)

    mutated = copy.deepcopy(recipe)
    mutated["sibling_evidence"]["late_morphology_status"] = "consumed"
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
