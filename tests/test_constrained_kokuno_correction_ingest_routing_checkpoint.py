from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.kokuno_correction_ingest_routing_checkpoint import (
    A4_GENERALIZATION_METRICS,
    FORMAL_GATES,
    SCHEMA,
    _canonical_sha256,
    build_checkpoint,
    validate_checkpoint,
)

AUDIT = Path(
    "artifacts/constrained/kokuno_agent3/agent4_563_public_z_pullback_independent_audit.json"
)


def _checkpoint() -> dict:
    return build_checkpoint(AUDIT.read_bytes())


def _resign(payload: dict) -> dict:
    payload["checkpoint_sha256"] = _canonical_sha256(payload)
    return payload


def test_v41_promotes_only_correction_ingest_permission() -> None:
    checkpoint = _checkpoint()
    assert checkpoint["schema"] == SCHEMA
    states = checkpoint["states"]

    assert states["oscillatory_ready"] is True
    assert states["oscillatory_generalization_independently_passed"] is True
    assert states["correction_receipt_identity_pinned_to_current_pass"] is True
    assert states["correction_ingest_allowed"] is True
    assert states["same_cycle_requested_stress_materialization_allowed"] is True
    assert states["candidate_finite_head_mean_debt_materialization_allowed"] is True

    assert states["leading_ready"] is False
    assert states["same_cycle_requested_stress_materialized"] is False
    assert states["candidate_numeric_finite_head_mean_debt_materialized"] is False
    assert states["real_candidate_defect_consumed"] is False
    assert states["signed_mean_inverse_input_ready"] is False
    assert states["correction_ready"] is False
    assert states["finite_correction_cycle_run"] is False
    assert states["candidate_artifact_instantiated"] is False
    assert states["velocity_export_ready"] is False
    assert states["formal_full_domain_pde_gate_assessed"] is False
    assert states["pde_validated"] is False
    assert checkpoint["formal_gates"] == FORMAL_GATES
    assert checkpoint["baseline_vs_kokuno"]["kokuno_current_comparable_full_domain_receipt"] is None
    validate_checkpoint(checkpoint)


def test_v41_roundtrip_is_deterministic() -> None:
    checkpoint = _checkpoint()
    encoded = json.dumps(checkpoint, sort_keys=True, separators=(",", ":"), allow_nan=False)
    replay = json.loads(encoded)
    assert json.dumps(replay, sort_keys=True, separators=(",", ":"), allow_nan=False) == encoded
    assert replay["checkpoint_sha256"] == _canonical_sha256(replay)
    validate_checkpoint(replay)


def test_v41_rejects_downstream_scientific_overpromotion() -> None:
    for key in (
        "same_cycle_requested_stress_materialized",
        "candidate_numeric_finite_head_mean_debt_materialized",
        "real_candidate_defect_consumed",
        "signed_mean_inverse_input_ready",
        "correction_ready",
        "finite_correction_cycle_run",
        "candidate_artifact_instantiated",
        "formal_full_domain_pde_gate_assessed",
        "pde_validated",
    ):
        mutated = copy.deepcopy(_checkpoint())
        mutated["states"][key] = True
        _resign(mutated)
        with pytest.raises(ValueError):
            validate_checkpoint(mutated)


def test_v41_rejects_identity_or_supplemental_receipt_drift() -> None:
    mutated = copy.deepcopy(_checkpoint())
    mutated["upstream"]["agent3_current_pass_admission"]["head"] = "0" * 40
    _resign(mutated)
    with pytest.raises(ValueError, match="Agent-3 execution identity mismatch"):
        validate_checkpoint(mutated)

    mutated = copy.deepcopy(_checkpoint())
    mutated["upstream"]["agent4_supplemental_generalization"]["metrics"][
        "relative_divergence_rms_by_step"
    ][-1] = A4_GENERALIZATION_METRICS["relative_divergence_rms_by_step"][-1] * 2.0
    _resign(mutated)
    with pytest.raises(ValueError, match="supplemental frozen metrics changed"):
        validate_checkpoint(mutated)


def test_v41_rejects_gate_or_baseline_laundering() -> None:
    mutated = copy.deepcopy(_checkpoint())
    mutated["formal_gates"]["held_out_normalized_momentum_max"] = 2.0e-3
    _resign(mutated)
    with pytest.raises(ValueError, match="formal PDE gates changed"):
        validate_checkpoint(mutated)

    mutated = copy.deepcopy(_checkpoint())
    mutated["baseline_vs_kokuno"]["kokuno_current_comparable_full_domain_receipt"] = {
        "fabricated": True
    }
    _resign(mutated)
    with pytest.raises(ValueError, match="no comparable Kokuno full-domain receipt"):
        validate_checkpoint(mutated)
