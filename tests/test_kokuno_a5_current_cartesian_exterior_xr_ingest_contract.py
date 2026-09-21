import copy
import json

import pytest

from openai_ns_reconstruction.kokuno_a5_current_cartesian_exterior_xr_ingest_contract import (
    AGENT1_EXTERIOR_TO_XR,
    AGENT2_IDENTITY_SAVE_LOAD_SIBLING,
    AGENT3_RADIAL_FORCE_SIBLING,
    LATEST_MATCHING_A4,
    PARENT_A5,
    READINESS,
    TRUTH_BOUNDARY,
    _sha256,
    build_contract,
    validate_contract,
    write_contract,
)


HEAD = "1" * 40


def _resign(payload):
    out = copy.deepcopy(payload)
    out.pop("contract_sha256", None)
    out["contract_sha256"] = _sha256(out)
    return out


def test_valid_contract_pins_exact_current_lineage_and_truth_boundary():
    payload = build_contract(HEAD)
    assert validate_contract(payload) == []

    assert payload["parent_a5"] == PARENT_A5
    assert payload["agent1_exterior_to_xr"] == AGENT1_EXTERIOR_TO_XR
    assert payload["agent2_identity_save_load_sibling"] == AGENT2_IDENTITY_SAVE_LOAD_SIBLING
    assert payload["agent3_radial_force_sibling"] == AGENT3_RADIAL_FORCE_SIBLING
    assert payload["latest_matching_a4"] == LATEST_MATCHING_A4
    assert payload["readiness"] == READINESS
    assert payload["truth_boundary"] == TRUTH_BOUNDARY

    truth = payload["truth_boundary"]
    assert truth["current_cartesian_leading_velocity_materialized_through_xr"] is True
    assert truth["velocity_beyond_xh_materialized"] is True
    assert truth["velocity_beyond_xr_materialized"] is False
    assert truth["post_xr_rf40_current_lineage_materialized"] is False
    assert truth["outer_global_leading_velocity_materialized"] is False
    assert truth["current_partial_composite_extended_through_xr"] is False
    assert truth["matched_cartesian_pressure_materialized"] is False
    assert truth["restricted_forcing_materialized"] is False
    assert truth["heldout_normalized_ns_residual_assessed"] is False
    assert payload["readiness"] == {
        "leading_ready": False,
        "oscillatory_ready": True,
        "correction_ready": False,
        "velocity_export_ready": False,
        "pde_validated": False,
    }


def test_stale_checksum_mutation_fails_closed():
    payload = build_contract(HEAD)
    payload["truth_boundary"]["velocity_beyond_xr_materialized"] = True
    errors = validate_contract(payload)
    assert "contract_sha256_mismatch" in errors
    assert "truth_boundary_drift" in errors
    assert "post_xr_scope_promoted" in errors


@pytest.mark.parametrize(
    ("key", "error"),
    [
        ("velocity_beyond_xr_materialized", "post_xr_scope_promoted"),
        ("outer_global_leading_velocity_materialized", "global_leading_promoted"),
        ("current_partial_composite_extended_through_xr", "composite_lineage_laundered_to_xr"),
        ("global_compact_support_completed", "truth_global_compact_support_completed_promoted"),
        ("matched_cartesian_pressure_materialized", "truth_matched_cartesian_pressure_materialized_promoted"),
        ("restricted_forcing_materialized", "truth_restricted_forcing_materialized_promoted"),
        ("complete_ns_defect_materialized", "truth_complete_ns_defect_materialized_promoted"),
        ("real_agent3_ns_correction_velocity_materialized", "truth_real_agent3_ns_correction_velocity_materialized_promoted"),
        ("heldout_normalized_ns_residual_assessed", "truth_heldout_normalized_ns_residual_assessed_promoted"),
        ("pde_validated", "truth_pde_validated_promoted"),
    ],
)
def test_rehashed_scientific_promotions_still_fail_closed(key, error):
    payload = build_contract(HEAD)
    if key == "pde_validated":
        payload["readiness"][key] = True
        payload = _resign(payload)
        errors = validate_contract(payload)
        assert "readiness_pde_validated_promoted" in errors
    else:
        payload["truth_boundary"][key] = True
        payload = _resign(payload)
        errors = validate_contract(payload)
        assert error in errors


def test_rehashed_agent2_partial_composite_cannot_be_laundered_to_xr():
    payload = build_contract(HEAD)
    payload["agent2_identity_save_load_sibling"]["consumed_by_this_increment"] = True
    payload["agent2_identity_save_load_sibling"]["composite_extended_to_agent1_980_X_R"] = True
    payload = _resign(payload)
    errors = validate_contract(payload)
    assert "agent2_identity_save_load_sibling_drift" in errors
    assert "agent2_sibling_silently_consumed" in errors
    assert "agent2_partial_composite_silently_extended" in errors


def test_rehashed_agent3_radial_force_cannot_be_retargeted_to_new_exterior():
    payload = build_contract(HEAD)
    payload["agent3_radial_force_sibling"]["consumed_by_this_increment"] = True
    payload["agent3_radial_force_sibling"]["consumes_agent1_980_exterior"] = True
    payload = _resign(payload)
    errors = validate_contract(payload)
    assert "agent3_radial_force_sibling_drift" in errors
    assert "agent3_sibling_silently_consumed" in errors
    assert "agent3_lineage_silently_retargeted" in errors


def test_missing_agent4_exterior_audit_cannot_be_promoted():
    payload = build_contract(HEAD)
    payload["latest_matching_a4"]["agent1_980_independent_a4_audit_available"] = True
    payload = _resign(payload)
    errors = validate_contract(payload)
    assert "latest_matching_a4_drift" in errors
    assert "missing_a4_audit_promoted" in errors


def test_invalid_exact_head_rejected():
    with pytest.raises(ValueError, match="40-hex"):
        build_contract("not-a-commit")


def test_write_contract_round_trip(tmp_path):
    path = tmp_path / "receipt.json"
    payload = write_contract(path, HEAD)
    loaded = json.loads(path.read_text(encoding="utf-8"))
    assert loaded == payload
    assert validate_contract(loaded) == []
