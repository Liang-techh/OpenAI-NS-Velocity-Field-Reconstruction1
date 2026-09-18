import json

import pytest

from openai_ns_reconstruction.kokuno_heat_v2_routing_checkpoint import (
    FIXED_GATES,
    STATES,
    ST006_REFERENCE,
    _sha,
    build_checkpoint,
    load_checkpoint,
    run_bound_heat_v2_audit,
    validate_checkpoint,
    write_bundle,
)


def _resign(payload):
    unsigned = dict(payload)
    unsigned.pop("checkpoint_sha256", None)
    payload["checkpoint_sha256"] = _sha(unsigned)
    return payload


def test_corrected_heat_v2_is_locally_unblocked_but_global_route_stays_closed():
    audit = run_bound_heat_v2_audit()
    assert audit["structural_preflight_passed"] is True
    assert audit["parameter_case_count"] == 3
    assert audit["training_loss_used"] is False
    assert audit["used_for_parameter_selection"] is False
    assert audit["used_as_formal_full_domain_pde_gate"] is False
    assert audit["st006_directly_comparable"] is False

    checkpoint = validate_checkpoint(build_checkpoint(audit))
    assert checkpoint["states"] == STATES
    assert checkpoint["states"]["delta_cp_tail_factor_blocker_resolved"] is True
    assert checkpoint["states"]["corrected_heat_discrepancy_independently_validated"] is True
    assert checkpoint["states"]["physical_three_bump_target_ready"] is False
    assert checkpoint["states"]["global_leading_profile_reconstructed"] is False
    assert checkpoint["states"]["correction_ready"] is False
    assert checkpoint["states"]["velocity_export_ready"] is False
    assert checkpoint["states"]["pde_validated"] is False
    assert checkpoint["routing"]["accept_corrected_heat_v2_as_local_source_moment_formula"] is True
    assert checkpoint["routing"]["use_corrected_heat_v2_as_physical_three_bump_target_now"] is False


def test_bundle_roundtrip_and_resigned_promotions_fail_closed(tmp_path):
    target = tmp_path / "bundle"
    checkpoint = write_bundle(target)
    loaded = load_checkpoint(target / "heat_v2_routing_checkpoint.json")
    assert loaded["checkpoint_sha256"] == checkpoint["checkpoint_sha256"]
    assert (target / "heat_v2_robustness_report.json").is_file()

    promote_target = json.loads(json.dumps(loaded))
    promote_target["routing"]["use_corrected_heat_v2_as_physical_three_bump_target_now"] = True
    _resign(promote_target)
    with pytest.raises(ValueError, match="blocked routing step"):
        validate_checkpoint(promote_target)

    promote_pde = json.loads(json.dumps(loaded))
    promote_pde["states"]["pde_validated"] = True
    _resign(promote_pde)
    with pytest.raises(ValueError, match="scientific state vector"):
        validate_checkpoint(promote_pde)

    relabel_st006 = json.loads(json.dumps(loaded))
    relabel_st006["st006_reference"]["directly_comparable_to_this_heat_moment_audit"] = True
    _resign(relabel_st006)
    with pytest.raises(ValueError, match="ST006 reference"):
        validate_checkpoint(relabel_st006)


def test_fixed_gates_st006_and_sibling_provenance_are_bound():
    checkpoint = build_checkpoint(run_bound_heat_v2_audit())
    assert checkpoint["fixed_gates"] == FIXED_GATES == {
        "held_out_normalized_full_momentum_residual": 1.0e-3,
        "divergence_max": 1.0e-5,
        "changed": False,
    }
    assert checkpoint["st006_reference"] == ST006_REFERENCE
    assert checkpoint["st006_reference"]["momentum_sampled_max"] == pytest.approx(
        0.1082289305112118
    )
    assert checkpoint["st006_reference"]["volume_l2"] == pytest.approx(
        0.10758432876230622
    )
    assert checkpoint["st006_reference"]["pde_validated"] is False

    assert checkpoint["upstream"]["agent1"]["consumed_in_executable_ancestry"] is True
    assert checkpoint["upstream"]["agent4"]["consumed_in_executable_ancestry"] is True
    assert checkpoint["upstream"]["agent2"]["consumed_in_executable_ancestry"] is False
    assert checkpoint["upstream"]["agent3"]["consumed_in_executable_ancestry"] is False
    assert checkpoint["upstream"]["agent3"]["one_column_cycle_accepted"] is False
    assert checkpoint["upstream"]["agent3"]["held_out_total_defect_ratio"] > 1.0
