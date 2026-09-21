from __future__ import annotations

import copy

import pytest

from openai_ns_reconstruction import (
    kokuno_a5_current_rf40_first_turn_composite_ingest_contract as mod,
)


EXACT_HEAD = "1" * 40


def _rehash(contract: dict) -> dict:
    payload = copy.deepcopy(contract)
    payload.pop("registration_sha256", None)
    contract["registration_sha256"] = mod._sha256(payload)
    return contract


def test_contract_registers_exact_first_turn_chain_and_stays_fail_closed() -> None:
    contract = mod.build_contract(EXACT_HEAD)
    mod.verify_contract(contract)

    assert contract["parent_a5"]["pr"] == 990
    assert contract["agent2_first_turn_composite"]["pr"] == 992
    assert contract["agent3_first_turn_nonlinear_mean"]["pr"] == 994
    assert contract["agent4_first_turn_audit"]["pr"] == 995
    assert contract["agent1_axial_shutdown_sibling"]["pr"] == 993

    assert (
        contract["agent3_first_turn_nonlinear_mean"]["consumed_agent2_head"]
        == contract["agent2_first_turn_composite"]["head"]
    )
    assert (
        contract["agent4_first_turn_audit"]["audited_agent2_head"]
        == contract["agent2_first_turn_composite"]["head"]
    )

    truth = contract["truth_boundary"]
    assert truth["current_leading_plus_oscillatory_velocity_through_rf40_first_turn_materialized"]
    assert truth["identity_preserving_first_turn_composite_save_load_available"]
    assert truth["current_nonlinear_m0_mean_through_rf40_first_turn_materialized"]
    assert truth["agent4_995_independent_first_turn_composite_divergence_audit_registered"]
    assert not truth["agent4_995_independent_first_turn_composite_divergence_audit_admitted"]
    assert truth["agent1_993_axial_shutdown_leading_only_materialized_as_sibling"]

    assert contract["readiness"] == {
        "leading_ready": False,
        "oscillatory_ready": True,
        "correction_ready": False,
        "velocity_export_ready": False,
        "pde_validated": False,
    }
    assert contract["final_gate"]["normalized_momentum_sampled_max"] == 1.0e-3
    assert contract["final_gate"]["normalized_momentum_volume_l2"] == 1.0e-3
    assert contract["final_gate"]["divergence_sampled_max"] == 1.0e-5
    assert contract["final_gate"]["divergence_volume_l2"] == 1.0e-5
    assert contract["final_gate"]["canonical_volume_quadrature_ladder"] == [24, 48, 96]
    assert contract["frozen_science"]["residual_defined_free_forcing_forbidden"] is True


def test_exact_heads_and_blobs_are_checksum_bound() -> None:
    contract = mod.build_contract(EXACT_HEAD)
    assert contract["agent2_first_turn_composite"]["head"] == (
        "2a2fad307a674ec1eb5a2bc426c0d01c3a0e7377"
    )
    assert contract["agent2_first_turn_composite"]["source_blob"] == (
        "c126dace161e5dae1745185c79152736817994b8"
    )
    assert contract["agent3_first_turn_nonlinear_mean"]["head"] == (
        "6088e3d3055e2792ba20fbd703fcc238ce1b27b5"
    )
    assert contract["agent3_first_turn_nonlinear_mean"]["source_blob"] == (
        "0f0e1f5c31cf33e28fb6f2851c2df7250bc582d1"
    )
    assert contract["agent4_first_turn_audit"]["head"] == (
        "70cd44ed0b5cabbcf5103df8959ff008377d17fa"
    )
    assert contract["agent4_first_turn_audit"]["source_blob"] == (
        "3881a0bb98e6369cf8fde85556c9d7ac0c2db7fc"
    )
    assert contract["agent1_axial_shutdown_sibling"]["head"] == (
        "2ac6460b483efb1c07f2fa65e7fed781a32f2718"
    )
    assert contract["agent1_axial_shutdown_sibling"]["source_blob"] == (
        "057a0514c8a941c3e59922b8158f480434b4441e"
    )


@pytest.mark.parametrize(
    "key",
    [
        "agent4_995_independent_first_turn_composite_divergence_audit_admitted",
        "current_leading_plus_oscillatory_velocity_through_axial_shutdown_materialized",
        "current_nonlinear_mean_through_axial_shutdown_materialized",
        "independent_axial_shutdown_composite_audit_available",
        "complete_post_xr_rf40_current_lineage_materialized",
        "matched_cartesian_pressure_materialized",
        "restricted_forcing_materialized",
        "complete_ns_defect_materialized",
        "current_nonlinear_mean_authorized_as_correction_target",
        "real_agent3_ns_correction_velocity_materialized",
        "complete_candidate_api_ready",
        "heldout_normalized_ns_residual_assessed",
        "same_protocol_st006_comparison_available_now",
        "residual_reduction_claimed",
        "scientific_admission",
        "paper_exact",
        "pde_validated" if False else "blowup_proved",
    ],
)
def test_truth_promotions_fail_closed(key: str) -> None:
    contract = mod.build_contract(EXACT_HEAD)
    contract["truth_boundary"][key] = True
    _rehash(contract)
    with pytest.raises(ValueError):
        mod.verify_contract(contract)


def test_axial_shutdown_sibling_cannot_be_laundered_into_matched_lineage() -> None:
    contract = mod.build_contract(EXACT_HEAD)
    contract["agent1_axial_shutdown_sibling"]["consumed_by_agent2_992"] = True
    _rehash(contract)
    with pytest.raises(ValueError):
        mod.verify_contract(contract)


def test_a4_registration_cannot_be_relabelled_as_scientific_admission() -> None:
    contract = mod.build_contract(EXACT_HEAD)
    contract["agent4_first_turn_audit"]["scientific_admission"] = True
    _rehash(contract)
    with pytest.raises(ValueError):
        mod.verify_contract(contract)


def test_fixed_gates_and_forcing_policy_cannot_drift() -> None:
    contract = mod.build_contract(EXACT_HEAD)
    contract["final_gate"]["normalized_momentum_sampled_max"] = 2.0e-3
    _rehash(contract)
    with pytest.raises(ValueError):
        mod.verify_contract(contract)

    contract = mod.build_contract(EXACT_HEAD)
    contract["frozen_science"]["residual_defined_free_forcing_forbidden"] = False
    _rehash(contract)
    with pytest.raises(ValueError):
        mod.verify_contract(contract)


def test_digest_and_exact_head_are_enforced() -> None:
    contract = mod.build_contract(EXACT_HEAD)
    contract["registration_sha256"] = "0" * 64
    with pytest.raises(ValueError):
        mod.verify_contract(contract)

    with pytest.raises(ValueError):
        mod.build_contract("not-a-git-sha")
