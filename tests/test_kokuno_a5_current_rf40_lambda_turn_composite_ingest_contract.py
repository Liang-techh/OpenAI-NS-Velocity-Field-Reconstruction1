from __future__ import annotations

import copy

import pytest

from openai_ns_reconstruction import (
    kokuno_a5_current_rf40_lambda_turn_composite_ingest_contract as mod,
)


def test_exact_lambda_turn_lineage_is_bound() -> None:
    receipt = mod.materialize_registration_receipt()
    mod.enforce_registration_receipt(receipt)
    reg = receipt["registration"]
    assert reg["parent_a5"]["head"] == "bb3d943461820fe43c0bdb2063cb11b428c4ea1c"
    assert reg["agent2_lambda_turn_composite"]["head"] == (
        "e9697cdfe2249e1b8bc2f9241bd16d456c60418f"
    )
    assert reg["agent2_lambda_turn_composite"]["leading_head"] == (
        "43b295444b1e9558222d757cb551e04385eddcb5"
    )
    assert reg["agent4_lambda_turn_audit"]["audited_agent2_head"] == (
        reg["agent2_lambda_turn_composite"]["head"]
    )


def test_truth_advance_stops_at_lambda_turn_composite_and_scoped_audit() -> None:
    truth = mod.materialize_registration_receipt()["registration"]["truth_boundary"]
    assert truth["current_leading_plus_oscillatory_velocity_through_rf40_lambda_turn_materialized"] is True
    assert truth["identity_preserving_lambda_turn_composite_save_load_available"] is True
    assert truth["agent4_1007_independent_lambda_turn_composite_divergence_audit_registered"] is True
    assert truth["agent4_1007_independent_lambda_turn_composite_divergence_audit_admitted"] is False
    for key in (
        "current_nonlinear_mean_through_axial_shutdown_materialized",
        "current_nonlinear_mean_through_lambda_turn_materialized",
        "current_nonlinear_radial_stress_through_lambda_turn_materialized",
        "current_nonlinear_radial_force_through_lambda_turn_materialized",
        "leading_plus_oscillatory_velocity_through_power_law_materialized",
        "independent_power_law_composite_audit_available",
        "complete_post_xr_rf40_current_lineage_materialized",
        "cone_i1_i2_i3_i4_outer_overlays_completed",
        "outer_global_leading_velocity_materialized",
        "matched_cartesian_pressure_materialized",
        "restricted_forcing_materialized",
        "complete_ns_defect_materialized",
        "real_agent3_ns_correction_velocity_materialized",
        "real_candidate_finite_correction_cycle_run",
        "complete_candidate_api_ready",
        "heldout_normalized_ns_residual_assessed",
        "same_protocol_st006_comparison_available_now",
        "residual_reduction_claimed",
        "scientific_admission",
    ):
        assert truth[key] is False, key


def test_newer_a1_and_lagging_a3_are_explicitly_nonconsumed() -> None:
    reg = mod.materialize_registration_receipt()["registration"]
    a1 = reg["agent1_power_law_sibling"]
    a3 = reg["agent3_first_turn_radial_force_sibling"]
    assert a1["pr"] == 1005
    assert "X_4" in a1["materialized_domain"]
    assert a1["consumed_by_agent2_1004"] is False
    assert a1["covered_by_agent4_1007"] is False
    assert a1["consumed_by_this_increment"] is False
    assert a3["pr"] == 1006
    assert "X_1" in a3["materialized_domain"]
    assert a3["rf40_axial_shutdown_current_lineage_consumed"] is False
    assert a3["rf40_lambda_turn_current_lineage_consumed"] is False
    assert a3["rf40_power_law_current_lineage_consumed"] is False
    assert a3["authorized_as_correction_target"] is False
    assert a3["real_cartesian_correction_velocity"] is False
    assert a3["consumed_by_this_increment"] is False


def test_a4_scope_cannot_be_laundered_into_full_ns_or_canonical_admission() -> None:
    a4 = mod.materialize_registration_receipt()["registration"]["agent4_lambda_turn_audit"]
    assert a4["seed"] == 9173731
    assert a4["fd2_step_ladder"] == [0.02, 0.01, 0.005]
    assert a4["strict_lambda_turn_heldout_points"] == 72
    assert a4["scoped_divergence_gate"] == 1.0e-5
    assert a4["uses_agent2_production_jacobian_or_divergence"] is False
    assert a4["uses_agent1_source_derivative_helpers"] is False
    assert a4["uses_pressure_forcing_or_residual_helper"] is False
    assert a4["canonical_24_48_96_whole_domain_admission"] is False
    assert a4["momentum_or_full_ns_evidence"] is False
    assert a4["scientific_admission"] is False


def test_readiness_baseline_and_final_gates_remain_fixed() -> None:
    reg = mod.materialize_registration_receipt()["registration"]
    assert reg["readiness"] == {
        "leading_ready": False,
        "oscillatory_ready": True,
        "correction_ready": False,
        "velocity_export_ready": False,
        "pde_validated": False,
    }
    assert reg["final_gate"] == {
        "normalized_momentum_sampled_max": 1.0e-3,
        "normalized_momentum_volume_l2": 1.0e-3,
        "divergence_sampled_max": 1.0e-5,
        "divergence_volume_l2": 1.0e-5,
        "canonical_volume_quadrature_ladder": [24, 48, 96],
    }
    assert reg["st006_baseline"]["momentum_sampled_max"] == pytest.approx(
        0.1082289305112118, rel=0.0, abs=0.0
    )
    assert reg["st006_baseline"]["momentum_volume_l2"] == pytest.approx(
        0.10758432876230622, rel=0.0, abs=0.0
    )
    assert reg["frozen_science"]["residual_defined_free_forcing_forbidden"] is True


def test_mutations_fail_closed_even_if_registration_digest_is_recomputed() -> None:
    receipt = mod.materialize_registration_receipt()

    promoted = copy.deepcopy(receipt)
    promoted["registration"]["truth_boundary"]["scientific_admission"] = True
    promoted["registration_sha256"] = mod._sha256(promoted["registration"])
    with pytest.raises(ValueError, match="registration payload drifted"):
        mod.enforce_registration_receipt(promoted)

    retargeted = copy.deepcopy(receipt)
    retargeted["registration"]["agent4_lambda_turn_audit"]["audited_agent2_head"] = "0" * 40
    retargeted["registration_sha256"] = mod._sha256(retargeted["registration"])
    with pytest.raises(ValueError, match="registration payload drifted|A4/A2"):
        mod.enforce_registration_receipt(retargeted)

    laundered_a3 = copy.deepcopy(receipt)
    laundered_a3["registration"]["agent3_first_turn_radial_force_sibling"][
        "rf40_lambda_turn_current_lineage_consumed"
    ] = True
    laundered_a3["registration_sha256"] = mod._sha256(laundered_a3["registration"])
    with pytest.raises(ValueError, match="registration payload drifted"):
        mod.enforce_registration_receipt(laundered_a3)

    weakened = copy.deepcopy(receipt)
    weakened["registration"]["final_gate"]["normalized_momentum_sampled_max"] = 2.0e-3
    weakened["registration_sha256"] = mod._sha256(weakened["registration"])
    with pytest.raises(ValueError, match="registration payload drifted"):
        mod.enforce_registration_receipt(weakened)


def test_checksum_schema_and_task_id_are_enforced() -> None:
    receipt = mod.materialize_registration_receipt()

    bad = copy.deepcopy(receipt)
    bad["registration_sha256"] = "0" * 64
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        mod.enforce_registration_receipt(bad)

    bad = copy.deepcopy(receipt)
    bad["schema"] = "wrong"
    with pytest.raises(ValueError, match="schema drifted"):
        mod.enforce_registration_receipt(bad)

    bad = copy.deepcopy(receipt)
    bad["task_id"] = "wrong"
    with pytest.raises(ValueError, match="task id drifted"):
        mod.enforce_registration_receipt(bad)
