import copy

import pytest

from openai_ns_reconstruction.kokuno_outer_patch_routing_checkpoint import (
    FIXED_GATES,
    ST006_REFERENCE,
    STATES,
    _sha,
    load_checkpoint,
    validate_checkpoint,
    write_bundle,
)


def test_outer_patch_routing_bundle_round_trip(tmp_path):
    checkpoint = write_bundle(tmp_path)
    loaded = load_checkpoint(tmp_path / "outer_patch_routing_checkpoint.json")

    assert loaded == checkpoint
    assert checkpoint["fixed_gates"] == FIXED_GATES
    assert checkpoint["st006_reference"] == ST006_REFERENCE
    assert checkpoint["states"] == STATES

    states = checkpoint["states"]
    assert states["X_star_e_star_binding_ready"] is True
    assert states["local_I2_patch_velocity_independently_validated"] is True
    assert states["terminal_tail_schedule_bound_to_outer_schedule"] is False
    assert states["physical_three_bump_target_ready"] is False
    assert states["global_leading_profile_reconstructed"] is False
    assert states["leading_ready"] is False
    assert states["correction_ready"] is False
    assert states["velocity_export_ready"] is False
    assert states["formal_full_domain_pde_gate_assessed"] is False
    assert states["pde_validated"] is False

    component = checkpoint["typed_component"]
    assert component["full_domain_candidate"] is False
    assert component["pressure_api_ready"] is False
    assert component["forcing_api_ready"] is False

    scales = checkpoint["derived_default_scales"]
    low, high = scales["I2_log_interval"]
    assert low < scales["log_X_star"] < high
    assert scales["X_star"] > 0.0
    assert scales["e_star"] > 0.0
    assert scales["autonomous_choice"] is True
    assert scales["hidden_source_numbers_recovered"] is False

    rows = checkpoint["baseline_vs_kokuno"]
    assert rows[0]["name"] == "ST006"
    assert rows[0]["momentum_sampled_max"] == pytest.approx(0.1082289305112118)
    assert rows[1]["directly_comparable_to_ST006"] is False
    assert rows[1]["momentum_sampled_max"] is None


def test_checkpoint_rejects_resigned_scientific_promotion(tmp_path):
    checkpoint = write_bundle(tmp_path)
    promoted = copy.deepcopy(checkpoint)
    promoted["states"]["pde_validated"] = True
    unsigned = dict(promoted)
    unsigned.pop("checkpoint_sha256")
    promoted["checkpoint_sha256"] = _sha(unsigned)

    with pytest.raises(ValueError, match="scientific state vector changed"):
        validate_checkpoint(promoted)


def test_checkpoint_rejects_resigned_local_to_global_promotion(tmp_path):
    checkpoint = write_bundle(tmp_path)
    promoted = copy.deepcopy(checkpoint)
    promoted["typed_component"]["full_domain_candidate"] = True
    unsigned = dict(promoted)
    unsigned.pop("checkpoint_sha256")
    promoted["checkpoint_sha256"] = _sha(unsigned)

    with pytest.raises(ValueError, match="local I2 component was promoted"):
        validate_checkpoint(promoted)
