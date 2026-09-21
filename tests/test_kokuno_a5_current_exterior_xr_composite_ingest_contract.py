from __future__ import annotations

import copy
import json

import pytest

from openai_ns_reconstruction import (
    kokuno_a5_current_exterior_xr_composite_ingest_contract as contract,
)


EXACT_HEAD = "a" * 40


def _rehash(payload):
    body = copy.deepcopy(payload)
    body.pop("contract_sha256", None)
    payload["contract_sha256"] = contract._sha256(body)
    return payload


def test_contract_round_trip_and_narrow_truth_advance(tmp_path):
    payload = contract.build_contract(EXACT_HEAD)
    assert contract.validate_contract(payload) == []

    truth = payload["truth_boundary"]
    assert truth["current_cartesian_leading_velocity_materialized_through_xr"] is True
    assert truth["current_leading_plus_oscillatory_velocity_through_xr_materialized"] is True
    assert truth["identity_preserving_composite_save_load_through_xr_available"] is True
    assert truth["current_nonlinear_m0_mean_attribution_through_xr_materialized"] is True
    assert truth["agent4_989_independent_composite_divergence_audit_registered"] is True
    assert truth["agent4_989_independent_composite_divergence_audit_admitted"] is False

    assert truth["agent1_986_leading_only_post_xr_first_turn_materialized_as_sibling"] is True
    assert truth["current_leading_plus_oscillatory_velocity_post_xr_materialized"] is False
    assert truth["current_nonlinear_mean_post_xr_materialized"] is False
    assert truth["outer_global_leading_velocity_materialized"] is False
    assert truth["matched_cartesian_pressure_materialized"] is False
    assert truth["restricted_forcing_materialized"] is False
    assert truth["complete_ns_defect_materialized"] is False
    assert truth["real_agent3_ns_correction_velocity_materialized"] is False
    assert truth["heldout_normalized_ns_residual_assessed"] is False

    assert payload["readiness"] == {
        "leading_ready": False,
        "oscillatory_ready": True,
        "correction_ready": False,
        "velocity_export_ready": False,
        "pde_validated": False,
    }

    out = tmp_path / "receipt.json"
    written = contract.write_contract(out, EXACT_HEAD)
    rebound = json.loads(out.read_text(encoding="utf-8"))
    assert rebound == written == payload
    assert contract.validate_contract(rebound) == []


def test_exact_upstream_identities_and_lineage_are_frozen():
    payload = contract.build_contract(EXACT_HEAD)
    assert payload["parent_a5"]["head"] == "dae4797fbd270a4f676b113f5ca3f7801d71a47e"
    assert payload["agent2_xr_composite"]["head"] == "4be2c9ee898c44dd1ad2217e90601b161fe81964"
    assert payload["agent3_xr_nonlinear_mean"]["head"] == "fa95dee3709a82326167af1a2cb42b9b5f9a88a5"
    assert payload["agent4_xr_composite_audit"]["head"] == "84b8aab2f5a56dbfc8ef02b6053d19bea910a959"
    assert payload["agent1_post_xr_first_turn_sibling"]["head"] == "23c00b98526e187ff04d432745300a084ec859f2"

    assert payload["agent3_xr_nonlinear_mean"]["consumed_agent2_head"] == payload["agent2_xr_composite"]["head"]
    assert payload["agent4_xr_composite_audit"]["audited_agent2_head"] == payload["agent2_xr_composite"]["head"]
    assert payload["agent4_xr_composite_audit"]["covers_agent1_986_post_XR_first_turn"] is False
    assert payload["agent1_post_xr_first_turn_sibling"]["consumed_by_agent2_987"] is False


def test_fixed_science_and_baseline_are_unchanged():
    payload = contract.build_contract(EXACT_HEAD)
    assert payload["frozen_science"]["viscosity"] == 0.01
    assert payload["frozen_science"]["residual_defined_free_forcing_forbidden"] is True
    assert payload["final_gate"] == {
        "normalized_momentum_sampled_max": 1.0e-3,
        "normalized_momentum_volume_l2": 1.0e-3,
        "divergence_sampled_max": 1.0e-5,
        "divergence_volume_l2": 1.0e-5,
        "canonical_volume_quadrature_ladder": [24, 48, 96],
    }
    assert payload["st006_baseline"]["momentum_sampled_max"] == pytest.approx(
        0.1082289305112118
    )
    assert payload["st006_baseline"]["momentum_volume_l2"] == pytest.approx(
        0.10758432876230622
    )


@pytest.mark.parametrize(
    ("section", "key", "value", "expected_fragment"),
    [
        ("readiness", "pde_validated", True, "pde_validated"),
        ("readiness", "leading_ready", True, "leading_ready"),
        ("readiness", "correction_ready", True, "correction_ready"),
        (
            "truth_boundary",
            "agent4_989_independent_composite_divergence_audit_admitted",
            True,
            "agent4_989_independent_composite_divergence_audit_admitted",
        ),
        (
            "truth_boundary",
            "current_leading_plus_oscillatory_velocity_post_xr_materialized",
            True,
            "current_leading_plus_oscillatory_velocity_post_xr_materialized",
        ),
        (
            "truth_boundary",
            "matched_cartesian_pressure_materialized",
            True,
            "matched_cartesian_pressure_materialized",
        ),
        (
            "truth_boundary",
            "current_xr_nonlinear_mean_authorized_as_correction_target",
            True,
            "current_xr_nonlinear_mean_authorized_as_correction_target",
        ),
        (
            "truth_boundary",
            "heldout_normalized_ns_residual_assessed",
            True,
            "heldout_normalized_ns_residual_assessed",
        ),
    ],
)
def test_checksum_valid_scientific_promotions_fail_closed(
    section, key, value, expected_fragment
):
    payload = contract.build_contract(EXACT_HEAD)
    payload[section][key] = value
    _rehash(payload)
    errors = contract.validate_contract(payload)
    assert errors
    assert any(expected_fragment in error for error in errors)


def test_a4_scope_cannot_be_laundered_into_post_xr_audit():
    payload = contract.build_contract(EXACT_HEAD)
    payload["agent4_xr_composite_audit"]["covers_agent1_986_post_XR_first_turn"] = True
    _rehash(payload)
    errors = contract.validate_contract(payload)
    assert "agent4_scope_laundering" in errors


def test_agent1_986_cannot_be_relabelled_as_consumed_by_agent2_987():
    payload = contract.build_contract(EXACT_HEAD)
    payload["agent1_post_xr_first_turn_sibling"]["consumed_by_agent2_987"] = True
    _rehash(payload)
    errors = contract.validate_contract(payload)
    assert "agent1_986_lineage_laundering" in errors


@pytest.mark.parametrize(
    ("key", "value", "expected"),
    [
        ("normalized_momentum_sampled_max", 1.1e-3, "momentum_gate_changed"),
        ("normalized_momentum_volume_l2", 1.1e-3, "momentum_l2_gate_changed"),
        ("divergence_sampled_max", 1.1e-5, "divergence_gate_changed"),
        ("divergence_volume_l2", 1.1e-5, "divergence_l2_gate_changed"),
        ("canonical_volume_quadrature_ladder", [16, 32, 64], "canonical_quadrature_changed"),
    ],
)
def test_gate_drift_fails_closed(key, value, expected):
    payload = contract.build_contract(EXACT_HEAD)
    payload["final_gate"][key] = value
    _rehash(payload)
    assert expected in contract.validate_contract(payload)


def test_free_forcing_firewall_cannot_be_relaxed():
    payload = contract.build_contract(EXACT_HEAD)
    payload["frozen_science"]["residual_defined_free_forcing_forbidden"] = False
    _rehash(payload)
    assert "free_forcing_firewall_relaxed" in contract.validate_contract(payload)


def test_invalid_exact_head_is_rejected():
    with pytest.raises(ValueError):
        contract.build_contract("not-a-sha")


def test_plain_checksum_mutation_is_detected():
    payload = contract.build_contract(EXACT_HEAD)
    payload["pipeline_position"]["stage"] = "pretend complete PDE"
    errors = contract.validate_contract(payload)
    assert "contract_sha256_mismatch" in errors
    assert "pipeline_position_drift" in errors
