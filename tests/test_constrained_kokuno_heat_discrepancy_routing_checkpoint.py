import json

import pytest

from openai_ns_reconstruction.kokuno_heat_discrepancy_routing_checkpoint import (
    STATES,
    _sha,
    build_checkpoint,
    load_checkpoint,
    run_independent_heat_audit,
    validate_checkpoint,
    write_bundle,
)


def _resign(payload):
    unsigned = dict(payload)
    unsigned.pop("checkpoint_sha256", None)
    payload["checkpoint_sha256"] = _sha(unsigned)
    return payload


def test_independent_heat_bug_is_bound_and_not_promoted():
    audit = run_independent_heat_audit()
    rel = audit["max_relative_by_component"]

    assert 2.3e-3 < rel[0] < 2.6e-3
    assert rel[1] < 2.0e-8
    assert rel[2] < 2.0e-6
    assert audit["agent1_heat_discrepancy_eligible_for_repair_target"] is False
    assert audit["used_for_parameter_selection"] is False
    assert audit["used_as_formal_full_domain_pde_gate"] is False

    checkpoint = validate_checkpoint(build_checkpoint(audit))
    assert checkpoint["states"] == STATES
    assert checkpoint["states"]["heat_discrepancy_promotion_blocked"] is True
    assert checkpoint["states"]["leading_ready"] is False
    assert checkpoint["states"]["correction_ready"] is False
    assert checkpoint["states"]["velocity_export_ready"] is False
    assert checkpoint["states"]["pde_validated"] is False
    assert checkpoint["upstream"]["agent2"]["consumed_in_executable_ancestry"] is False
    assert checkpoint["upstream"]["agent3"]["consumed_in_executable_ancestry"] is False


def test_checkpoint_roundtrip_and_resigned_promotions_fail_closed(tmp_path):
    bundle = tmp_path / "bundle"
    checkpoint = write_bundle(bundle)
    loaded = load_checkpoint(bundle / "heat_discrepancy_routing_checkpoint.json")
    assert loaded["checkpoint_sha256"] == checkpoint["checkpoint_sha256"]

    promote_heat = json.loads(json.dumps(loaded))
    promote_heat["routing"]["use_agent1_280_heat_discrepancy_as_repair_target"] = True
    _resign(promote_heat)
    with pytest.raises(ValueError, match="blocked routing step"):
        validate_checkpoint(promote_heat)

    promote_pde = json.loads(json.dumps(loaded))
    promote_pde["states"]["pde_validated"] = True
    _resign(promote_pde)
    with pytest.raises(ValueError, match="scientific state vector"):
        validate_checkpoint(promote_pde)


def test_fixed_gate_and_missing_second_direction_truth_are_preserved():
    checkpoint = build_checkpoint(run_independent_heat_audit())
    assert checkpoint["fixed_gates"] == {
        "held_out_normalized_full_momentum_residual": 1.0e-3,
        "divergence_max": 1.0e-5,
        "changed": False,
    }
    agent3 = checkpoint["upstream"]["agent3"]
    assert agent3["active_nodes"] == 27
    assert agent3["active_nodes_requiring_second_direction"] == 25
    assert agent3["missing_second_column_target_measured"] is True
    assert agent3["new_oscillatory_column_constructed"] is False
