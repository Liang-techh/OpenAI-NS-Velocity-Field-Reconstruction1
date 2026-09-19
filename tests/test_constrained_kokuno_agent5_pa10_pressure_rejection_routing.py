from __future__ import annotations

import copy
import json

import pytest

from openai_ns_reconstruction.kokuno_agent4_pa10_pressure_primitive_audit import run_audit
from openai_ns_reconstruction.kokuno_agent5_pa10_pressure_rejection_routing import (
    EXPECTED_FAILED_GUARDS,
    FORMAL_GATES,
    _canonical_sha256,
    build_checkpoint,
    validate_checkpoint,
)


@pytest.fixture(scope="module")
def audit() -> dict:
    return run_audit()


@pytest.fixture(scope="module")
def checkpoint(audit: dict) -> dict:
    result = build_checkpoint(audit)
    validate_checkpoint(result)
    return result


def _resign(checkpoint: dict) -> dict:
    checkpoint["checkpoint_sha256"] = _canonical_sha256(checkpoint)
    return checkpoint


def test_routes_exact_independent_scientific_reject(checkpoint: dict) -> None:
    upstream = checkpoint["upstream"]["agent4_independent_pressure_audit"]
    assert upstream["scientific_verdict"] == "REJECT"
    assert upstream["failed_guards"] == EXPECTED_FAILED_GUARDS
    states = checkpoint["states"]
    assert states["selected_pressure_primitive_public_api_ready"] is True
    assert states["selected_pressure_primitive_independently_audited"] is False
    assert states["selected_pressure_primitive_independent_preflight_passed"] is False
    assert states["agent1_public_pressure_value_path_repair_required"] is True
    assert states["leading_ready"] is False


def test_keeps_informative_and_blocked_channels_separate(checkpoint: dict) -> None:
    audit = checkpoint["upstream"]["agent4_independent_pressure_audit"]
    assert audit["p_Y"]["relative_rms_ladder"][-1] == pytest.approx(2.22085e-14, rel=2e-5)
    assert audit["Phi_eta"]["relative_rms_ladder"][-1] == pytest.approx(3.98927e-6, rel=2e-5)
    assert audit["Phi_eta"]["independent_channel_passed"] is True
    assert audit["Y_Phi_Y"]["reference_rms"] == pytest.approx(2.48457e-27, rel=2e-5)
    assert audit["Y_Phi_Y"]["binary64_identifiable_under_frozen_stationary_cloud"] is False
    assert audit["p_eta"]["relative_rms_ladder"] == pytest.approx(
        [0.999114, 0.903906, 1.226775], rel=2e-5
    )
    assert audit["p_eta"]["reference_rms"] == pytest.approx(1.570226741493718e15, rel=1e-12)
    assert audit["p_eta"]["substantive_blocker"] is True


def test_preserves_final_project_gates_and_downstream_false(checkpoint: dict) -> None:
    assert checkpoint["formal_gates_unchanged"] == FORMAL_GATES
    states = checkpoint["states"]
    assert states["oscillatory_ready"] is True
    assert states["correction_ingest_allowed"] is True
    for key in (
        "same_cycle_requested_stress_materialized",
        "candidate_numeric_finite_head_mean_debt_materialized",
        "real_full_candidate_defect_consumed",
        "signed_mean_inverse_input_ready",
        "public_velocity_correction_materialized",
        "correction_ready",
        "finite_correction_cycle_run",
        "candidate_artifact_instantiated",
        "velocity_export_ready",
        "formal_full_domain_pde_gate_assessed",
        "heldout_ns_residual_assessed",
        "pde_validated",
    ):
        assert states[key] is False


def test_checkpoint_is_json_roundtrip_and_deterministic(audit: dict, checkpoint: dict) -> None:
    roundtrip = json.loads(json.dumps(checkpoint, sort_keys=True, allow_nan=False))
    validate_checkpoint(roundtrip)
    rebuilt = build_checkpoint(audit)
    assert rebuilt["checkpoint_sha256"] == checkpoint["checkpoint_sha256"]
    assert rebuilt == checkpoint


def test_rejects_pressure_pass_laundering(checkpoint: dict) -> None:
    bad = copy.deepcopy(checkpoint)
    bad["states"]["selected_pressure_primitive_independently_audited"] = True
    bad["states"]["selected_pressure_primitive_independent_preflight_passed"] = True
    _resign(bad)
    with pytest.raises(ValueError, match="laundered"):
        validate_checkpoint(bad)


def test_rejects_leading_or_pde_promotion(checkpoint: dict) -> None:
    for key in ("leading_ready", "correction_ready", "pde_validated"):
        bad = copy.deepcopy(checkpoint)
        bad["states"][key] = True
        _resign(bad)
        with pytest.raises(ValueError, match="over-promoted"):
            validate_checkpoint(bad)


def test_rejects_final_threshold_relaxation(checkpoint: dict) -> None:
    bad = copy.deepcopy(checkpoint)
    bad["formal_gates_unchanged"]["held_out_normalized_momentum_l2"] = 2.0e-3
    _resign(bad)
    with pytest.raises(ValueError, match="formal project gates changed"):
        validate_checkpoint(bad)


def test_rejects_changed_agent4_failure_signature(audit: dict) -> None:
    bad = copy.deepcopy(audit)
    bad["failed_guards"] = list(EXPECTED_FAILED_GUARDS[:-1])
    with pytest.raises(ValueError, match="failure signature changed"):
        build_checkpoint(bad)


def test_rejects_agent4_scientific_pass_laundering(audit: dict) -> None:
    bad = copy.deepcopy(audit)
    bad["selected_pressure_primitive_independent_preflight_passed"] = True
    bad["truth_boundary"]["selected_center_pressure_primitive_independently_audited"] = True
    bad["failed_guards"] = []
    with pytest.raises(ValueError, match="scientific REJECT"):
        build_checkpoint(bad)
