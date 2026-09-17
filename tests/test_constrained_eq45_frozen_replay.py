import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.constrained_eq45_frozen_replay import (
    EXPECTED_CANDIDATE_SHA256,
    replay_eq45_candidate,
    save_replay_receipt,
)


ROOT = Path(__file__).resolve().parents[1]
SEED = ROOT / "artifacts/constrained/eq45_velocity_candidate_seed.json"


def test_checked_eq45_seed_replays_deterministically_through_public_velocity(tmp_path):
    first = replay_eq45_candidate(SEED)
    second = replay_eq45_candidate(SEED)

    assert first["candidate"]["sha256"] == EXPECTED_CANDIDATE_SHA256
    assert first["candidate"]["evaluator"] == "Eq45VelocityCandidate.at_points"
    assert first["replay"]["probe_fingerprint_sha256"] == second["replay"]["probe_fingerprint_sha256"]
    assert first["replay"]["grid_fingerprint_sha256"] == second["replay"]["grid_fingerprint_sha256"]
    assert first["replay"]["probe_speed_rms"] > 1e-6
    assert first["replay"]["reference_times"] == [0.25, 0.5, 0.75]

    states = first["states"]
    assert states["velocity_export_ready"] is True
    assert states["visualization_ready"] is False
    assert states["physical_support_validated"] is False
    assert states["pde_validated"] is False
    assert states["visual_correspondence_verified"] is False
    assert states["paper_exact"] is False
    assert states["openai_field_identified"] is False

    output = tmp_path / "receipt.json"
    save_replay_receipt(first, output)
    loaded = json.loads(output.read_text(encoding="utf-8"))
    assert loaded["replay"]["probe_fingerprint_sha256"] == first["replay"]["probe_fingerprint_sha256"]


def test_candidate_parameter_drift_fails_frozen_hash_contract(tmp_path):
    payload = json.loads(SEED.read_text(encoding="utf-8"))
    payload["profile_basis"]["phi_coefficients"][0] += 0.125
    mutated = tmp_path / "mutated.json"
    mutated.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(ValueError, match="SHA256"):
        replay_eq45_candidate(mutated)


def test_receipt_cannot_promote_visual_or_pde_state(tmp_path):
    receipt = replay_eq45_candidate(SEED, probe_count=16, grid_size=3)

    receipt["states"]["pde_validated"] = True
    with pytest.raises(ValueError, match="unsupported scientific state"):
        save_replay_receipt(receipt, tmp_path / "pde.json")

    receipt = replay_eq45_candidate(SEED, probe_count=16, grid_size=3)
    receipt["states"]["visualization_ready"] = True
    with pytest.raises(ValueError, match="unsupported scientific state"):
        save_replay_receipt(receipt, tmp_path / "visual.json")


def test_replay_inputs_fail_closed():
    with pytest.raises(ValueError, match="probe_count"):
        replay_eq45_candidate(SEED, probe_count=7)
    with pytest.raises(ValueError, match="odd integer"):
        replay_eq45_candidate(SEED, grid_size=4)
    with pytest.raises(ValueError, match="hexadecimal"):
        replay_eq45_candidate(SEED, expected_candidate_sha256="z" * 64)
