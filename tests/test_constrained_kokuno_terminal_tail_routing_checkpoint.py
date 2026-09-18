import copy

import pytest

from openai_ns_reconstruction.kokuno_terminal_tail_routing_checkpoint import (
    FIXED_GATES,
    STATES,
    _sha,
    load_checkpoint,
    validate_checkpoint,
    write_bundle,
)


def test_terminal_tail_routing_bundle_round_trip(tmp_path):
    checkpoint = write_bundle(tmp_path)
    loaded = load_checkpoint(tmp_path / "terminal_tail_routing_checkpoint.json")
    assert loaded == checkpoint
    assert checkpoint["fixed_gates"] == FIXED_GATES
    assert checkpoint["states"] == STATES

    states = checkpoint["states"]
    assert states["terminal_tail_schedule_bound_to_outer_schedule"] is True
    assert states["terminal_tail_schedule_independently_preflighted"] is True
    assert states["terminal_scales_ready_in_log_space"] is True
    assert states["default_terminal_heat_inputs_float_materializable"] is False
    assert states["physical_three_bump_target_ready"] is False
    assert states["global_leading_profile_reconstructed"] is False
    assert states["second_public_covariance_column_realized"] is False
    assert states["leading_ready"] is False
    assert states["correction_ready"] is False
    assert states["velocity_export_ready"] is False
    assert states["formal_full_domain_pde_gate_assessed"] is False
    assert states["pde_validated"] is False

    typed = checkpoint["typed_schedule_component"]
    assert typed["default_complete_heat_materialization_ready"] is False
    assert typed["full_domain_velocity_api_ready"] is False
    assert typed["pressure_api_ready"] is False
    assert typed["forcing_api_ready"] is False

    scales = checkpoint["default_terminal_scales"]
    assert scales["log_X_tail"] > 709.0
    assert scales["X_tail"] is None
    assert scales["c_inf"] > 0.0
    assert scales["complete_heat_inputs_materializable"] is False
    assert scales["source_hidden_numeric_choices_recovered"] is False

    audit = checkpoint["independent_terminal_tail_audit"]
    assert audit["summary"]["structural_preflight_passed"] is True
    assert audit["default_extreme_public_materialized"] is False
    assert audit["moderate_public_materialized"] is True
    assert audit["formal_full_domain_pde_gate_assessed"] is False
    assert audit["st006_directly_comparable"] is False

    rows = checkpoint["baseline_vs_kokuno"]
    assert rows[0]["name"] == "ST006"
    assert rows[0]["momentum_sampled_max"] == pytest.approx(0.1082289305112118)
    assert rows[1]["directly_comparable_to_ST006"] is False
    assert rows[1]["momentum_sampled_max"] is None


def test_terminal_tail_checkpoint_rejects_resigned_pde_promotion(tmp_path):
    checkpoint = write_bundle(tmp_path)
    promoted = copy.deepcopy(checkpoint)
    promoted["states"]["pde_validated"] = True
    unsigned = dict(promoted)
    unsigned.pop("checkpoint_sha256")
    promoted["checkpoint_sha256"] = _sha(unsigned)
    with pytest.raises(ValueError, match="scientific state vector changed"):
        validate_checkpoint(promoted)


def test_terminal_tail_checkpoint_rejects_resigned_st006_relabel(tmp_path):
    checkpoint = write_bundle(tmp_path)
    promoted = copy.deepcopy(checkpoint)
    promoted["baseline_vs_kokuno"][1]["directly_comparable_to_ST006"] = True
    unsigned = dict(promoted)
    unsigned.pop("checkpoint_sha256")
    promoted["checkpoint_sha256"] = _sha(unsigned)
    with pytest.raises(ValueError, match="ST006-comparable"):
        validate_checkpoint(promoted)
