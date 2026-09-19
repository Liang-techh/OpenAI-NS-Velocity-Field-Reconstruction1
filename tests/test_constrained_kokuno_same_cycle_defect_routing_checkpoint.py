from __future__ import annotations

import copy
import inspect

import pytest

from openai_ns_reconstruction.kokuno_same_cycle_defect_contract import (
    evaluate_disjoint_heldin_heldout,
    evaluate_same_cycle_momentum_defect,
)
from openai_ns_reconstruction.kokuno_same_cycle_defect_routing_checkpoint import (
    AGENT3_ARTIFACT_DIGEST,
    AGENT3_ARTIFACT_ID,
    AGENT3_DEDICATED_RUN,
    AGENT3_HEAD,
    AGENT3_PR,
    FORMAL_GATES,
    ST006_MOMENTUM_L2,
    ST006_MOMENTUM_MAX,
    _canonical_sha256,
    build_checkpoint,
    validate_checkpoint,
)


def _resign(checkpoint: dict[str, object]) -> dict[str, object]:
    checkpoint["checkpoint_sha256"] = _canonical_sha256(checkpoint)
    return checkpoint


def test_v47_checkpoint_builds_and_validates() -> None:
    checkpoint = build_checkpoint()
    validate_checkpoint(checkpoint)

    states = checkpoint["states"]
    assert states["same_cycle_momentum_defect_contract_ready"] is True
    assert states["actual_momentum_defect_formula_executable"] is True
    assert states["same_cycle_identity_enforced"] is True
    assert states["restricted_forcing_provider_required"] is True
    assert states["heldin_heldout_disjointness_enforced"] is True
    assert states["real_full_candidate_defect_consumed"] is False
    assert states["correction_ready"] is False
    assert states["candidate_artifact_instantiated"] is False
    assert states["pde_validated"] is False


def test_v47_pins_exact_agent3_execution_and_artifact() -> None:
    checkpoint = build_checkpoint()
    a3 = checkpoint["upstream"]["agent3_same_cycle_defect_contract"]
    assert a3["pr"] == AGENT3_PR
    assert a3["head"] == AGENT3_HEAD
    assert a3["dedicated_workflow_run"] == AGENT3_DEDICATED_RUN
    assert a3["dedicated_workflow_conclusion"] == "success"
    assert a3["artifact_id"] == AGENT3_ARTIFACT_ID
    assert a3["artifact_zip_digest"] == AGENT3_ARTIFACT_DIGEST
    assert a3["analytic_expected_residual"] == [2.0, 0.0, 3.0]
    assert a3["analytic_maximum_absolute_error"] <= 2.0e-11


def test_v47_keeps_contract_api_raw_and_noninjectable() -> None:
    defect_names = set(inspect.signature(evaluate_same_cycle_momentum_defect).parameters)
    split_names = set(inspect.signature(evaluate_disjoint_heldin_heldout).parameters)
    forbidden = {
        "residual",
        "defect",
        "stress",
        "target",
        "pressure_gradient",
        "gain",
        "normalized_score",
    }
    assert defect_names.isdisjoint(forbidden)
    assert split_names.isdisjoint(forbidden)
    assert "restricted_forcing" in defect_names
    assert "restricted_forcing" in split_names
    assert {"held_in_points", "held_out_points"}.issubset(split_names)


def test_v47_preserves_baseline_and_final_gates() -> None:
    checkpoint = build_checkpoint()
    assert checkpoint["formal_gates_unchanged"] == FORMAL_GATES
    baseline = checkpoint["baseline_vs_kokuno"]
    assert baseline["st006_same_protocol_baseline"]["momentum_sampled_max"] == ST006_MOMENTUM_MAX
    assert baseline["st006_same_protocol_baseline"]["momentum_volume_l2"] == ST006_MOMENTUM_L2
    assert baseline["kokuno_full_same_protocol_residual_available"] is False
    assert baseline["current_agent3_contract_regression_is_comparable_to_st006"] is False


def test_v47_rejects_premature_pde_or_correction_promotion() -> None:
    for state in (
        "real_full_candidate_defect_consumed",
        "same_cycle_requested_stress_materialized",
        "correction_ready",
        "candidate_artifact_instantiated",
        "formal_full_domain_pde_gate_assessed",
        "pde_validated",
    ):
        mutated = copy.deepcopy(build_checkpoint())
        mutated["states"][state] = True
        _resign(mutated)
        with pytest.raises(ValueError, match="premature downstream promotion"):
            validate_checkpoint(mutated)


def test_v47_rejects_residual_or_forcing_truth_laundering() -> None:
    mutated = copy.deepcopy(build_checkpoint())
    mutated["typed_handoffs"]["same_cycle_defect"]["caller_supplied_residual_allowed"] = True
    _resign(mutated)
    with pytest.raises(ValueError, match="caller residual laundering"):
        validate_checkpoint(mutated)

    mutated = copy.deepcopy(build_checkpoint())
    mutated["states"]["restricted_forcing_semantics_independently_validated"] = True
    _resign(mutated)
    with pytest.raises(ValueError, match="premature downstream promotion"):
        validate_checkpoint(mutated)


def test_v47_rejects_threshold_or_agent3_identity_drift() -> None:
    mutated = copy.deepcopy(build_checkpoint())
    mutated["formal_gates_unchanged"]["held_out_normalized_momentum_max"] = 2.0e-3
    _resign(mutated)
    with pytest.raises(ValueError, match="formal PDE gates changed"):
        validate_checkpoint(mutated)

    mutated = copy.deepcopy(build_checkpoint())
    mutated["upstream"]["agent3_same_cycle_defect_contract"]["artifact_zip_digest"] = "sha256:wrong"
    _resign(mutated)
    with pytest.raises(ValueError, match="wrong Agent-3 artifact digest"):
        validate_checkpoint(mutated)
