import copy
import json
from pathlib import Path

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_candidate import Eq45VelocityCandidate
from openai_ns_reconstruction.constrained_eq45_optimized_artifact import (
    CANONICAL_SEED_SHA256,
    OPTIMIZED_CANDIDATE_SHA256,
    candidate_from_optimization_receipt,
    load_checked_candidate,
    write_checked_candidate,
)


ROOT = Path(__file__).resolve().parents[1]
OPTIMIZATION = ROOT / "artifacts" / "constrained" / "eq45_profile_force_optimization.json"
SEED = ROOT / "artifacts" / "constrained" / "eq45_velocity_candidate_seed.json"
SAVED = ROOT / "artifacts" / "constrained" / "eq45_profile_force_optimized_candidate.json"


def _optimization_payload():
    return json.loads(OPTIMIZATION.read_text(encoding="utf-8"))


def _write_json(path, payload):
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def test_optimized_artifact_replays_optimizer_receipt_through_public_velocity(tmp_path):
    seed = Eq45VelocityCandidate.load_json(SEED)
    rebuilt = candidate_from_optimization_receipt(OPTIMIZATION, SEED)
    saved = load_checked_candidate(SAVED, OPTIMIZATION, SEED)

    assert seed.sha256 == CANONICAL_SEED_SHA256
    assert rebuilt.sha256 == OPTIMIZED_CANDIDATE_SHA256
    assert saved.sha256 == OPTIMIZED_CANDIDATE_SHA256
    assert saved.to_dict() == rebuilt.to_dict()

    mode = seed.profile_basis.mode_indices.index((0, 2))
    assert rebuilt.profile_basis.phi_coefficients[mode] == pytest.approx(3.8740984117202455)
    assert rebuilt.profile_basis.swirl_coefficients[mode] == pytest.approx(3.9916795387937536)
    for index in range(len(seed.profile_basis.mode_indices)):
        if index == mode:
            continue
        assert rebuilt.profile_basis.phi_coefficients[index] == seed.profile_basis.phi_coefficients[index]
        assert rebuilt.profile_basis.swirl_coefficients[index] == seed.profile_basis.swirl_coefficients[index]

    points = np.array(
        [
            [0.13, -0.07, 0.11],
            [-0.22, 0.18, -0.09],
            [0.31, 0.04, 0.16],
            [-0.08, -0.27, -0.14],
        ],
        dtype=float,
    )
    times = np.array([0.29, 0.41, 0.58, 0.71], dtype=float)
    direct = saved.at_points(points, times)
    xyz = saved.velocity_xyz(points[:, 0], points[:, 1], points[:, 2], times)
    np.testing.assert_allclose(direct, xyz, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(direct, rebuilt.at_points(points, times), rtol=0.0, atol=0.0)
    assert np.all(np.isfinite(direct))
    assert np.linalg.norm(direct) > 0.0

    regenerated_path = tmp_path / "candidate.json"
    regenerated = write_checked_candidate(regenerated_path, OPTIMIZATION, SEED)
    assert regenerated.sha256 == OPTIMIZED_CANDIDATE_SHA256
    assert regenerated_path.read_text(encoding="utf-8") == SAVED.read_text(encoding="utf-8")

    truth = saved.to_dict()["truth_boundary"]
    assert truth["velocity_export_ready"] is True
    assert truth["callable_serializable"] is True
    for key in (
        "physical_support_validated",
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        assert truth[key] is False


def test_optimized_artifact_receipt_mutations_fail_closed(tmp_path):
    mutations = []

    wrong_sha = copy.deepcopy(_optimization_payload())
    wrong_sha["optimized_candidate"]["sha256"] = "0" * 64
    mutations.append(wrong_sha)

    promoted_seed = copy.deepcopy(_optimization_payload())
    promoted_seed["optimized_candidate"]["canonical_seed_promoted"] = True
    mutations.append(promoted_seed)

    promoted_pde = copy.deepcopy(_optimization_payload())
    promoted_pde["status"]["pde_validated"] = True
    mutations.append(promoted_pde)

    wrong_parameter_set = copy.deepcopy(_optimization_payload())
    wrong_parameter_set["contract"]["profile_parameters"] = ["Phi(0,2)", "F(0,4)"]
    mutations.append(wrong_parameter_set)

    out_of_bounds = copy.deepcopy(_optimization_payload())
    out_of_bounds["fit"]["Phi_02"] = 4.5
    mutations.append(out_of_bounds)

    for index, payload in enumerate(mutations):
        path = tmp_path / f"mutation_{index}.json"
        _write_json(path, payload)
        with pytest.raises(ValueError):
            candidate_from_optimization_receipt(path, SEED)


def test_checked_saved_artifact_detects_candidate_drift(tmp_path):
    payload = json.loads(SAVED.read_text(encoding="utf-8"))
    payload["profile_basis"]["phi_coefficients"][0] = 1.01
    drifted = tmp_path / "drifted_candidate.json"
    _write_json(drifted, payload)
    with pytest.raises(ValueError, match="saved optimized candidate"):
        load_checked_candidate(drifted, OPTIMIZATION, SEED)
