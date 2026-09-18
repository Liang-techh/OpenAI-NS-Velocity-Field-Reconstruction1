from __future__ import annotations

import json

import pytest

from openai_ns_reconstruction.kokuno_covariance_rank_routing_checkpoint import (
    FORMAL_GATES,
    build_checkpoint,
    checkpoint_sha256,
    validate_checkpoint,
    write_checkpoint,
)


def test_checkpoint_preserves_structural_pass_and_covariance_rank_reject() -> None:
    payload = build_checkpoint()
    validate_checkpoint(payload)

    states = payload["states"]
    assert states["source_multiband_public_contract_structural_preflight_passed"]
    assert states["source_multiband_independent_cartesian_audit_passed"]
    assert states["supplied_multiband_raw_velocity_rank_two"]
    assert not states["supplied_multiband_phase_mean_covariance_rank_two"]
    assert not states["genuinely_independent_second_covariance_column_ready"]
    assert not states["correction_ready"]
    assert not states["pde_validated"]

    a4 = payload["upstream"]["agent4"]
    assert a4["simultaneous_beta_permutation_points"] == 36
    assert a4["simultaneous_beta_permutation_failures"] == 0
    assert a4["simultaneous_beta_permutation_relative_max"] == 0.0
    assert a4["local_structural_preflight_passed"]

    a3 = payload["upstream"]["agent3"]
    assert a3["raw_velocity_rank_ratio"] > 0.5
    assert a3["phase_mean_covariance_rank_two_cells"] == 0
    assert a3["radial_cells"] == 5
    assert not a3["local_supplied_family_covariance_rank_two"]


def test_checkpoint_keeps_formal_gates_and_st006_truth_boundary() -> None:
    payload = build_checkpoint()
    assert payload["formal_gates"] == FORMAL_GATES
    assert payload["formal_gates"]["held_out_normalized_momentum_max"] == 1.0e-3
    assert payload["formal_gates"]["held_out_divergence_max"] == 1.0e-5

    baseline = payload["baseline_vs_kokuno"]["st006"]
    assert baseline["momentum_sampled_max"] == pytest.approx(0.1082289305112118)
    assert baseline["momentum_volume_l2"] == pytest.approx(0.10758432876230622)
    assert not baseline["pde_validated"]
    assert payload["baseline_vs_kokuno"]["kokuno_current_comparable_full_domain_receipt"] is None


def test_checkpoint_fails_closed_if_second_direction_or_pde_state_is_promoted() -> None:
    payload = build_checkpoint()
    payload["states"]["genuinely_independent_second_covariance_column_ready"] = True
    payload["checkpoint_sha256"] = checkpoint_sha256(payload)
    with pytest.raises(ValueError, match="fail-closed state promoted"):
        validate_checkpoint(payload)

    payload = build_checkpoint()
    payload["states"]["pde_validated"] = True
    payload["checkpoint_sha256"] = checkpoint_sha256(payload)
    with pytest.raises(ValueError, match="fail-closed state promoted"):
        validate_checkpoint(payload)


def test_checkpoint_sha_and_roundtrip(tmp_path) -> None:
    path = tmp_path / "checkpoint.json"
    payload = write_checkpoint(path)
    loaded = json.loads(path.read_text(encoding="utf-8"))
    assert loaded == payload
    assert loaded["checkpoint_sha256"] == checkpoint_sha256(loaded)

    loaded["upstream"]["agent3"]["phase_mean_covariance_rank_two_cells"] = 1
    with pytest.raises(ValueError, match="checkpoint sha256 mismatch"):
        validate_checkpoint(loaded)
