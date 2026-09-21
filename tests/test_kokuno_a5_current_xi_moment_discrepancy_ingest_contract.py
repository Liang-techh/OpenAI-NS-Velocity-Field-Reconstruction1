from __future__ import annotations

import copy
import json

import pytest

from openai_ns_reconstruction.kokuno_a5_current_xi_moment_discrepancy_ingest_contract import (
    AGENT1_XI_PREFIX_MOMENTS,
    AGENT2_SIBLING_STATUS,
    AGENT3_XI_DISCREPANCY_BRIDGE,
    AGENT4_XI_MOMENT_AUDIT,
    FINAL_GATE,
    PARENT_A5,
    READINESS,
    ST006_BASELINE,
    TRUTH_BOUNDARY,
    build_contract,
    validate_contract,
    write_contract,
)

HEAD = "0123456789abcdef0123456789abcdef01234567"


def test_exact_parent_and_agent_pins_are_frozen() -> None:
    contract = build_contract(HEAD)
    assert contract["parent_a5"] == PARENT_A5
    assert contract["agent1_xi_prefix_moments"] == AGENT1_XI_PREFIX_MOMENTS
    assert contract["agent2_sibling_status"] == AGENT2_SIBLING_STATUS
    assert contract["agent3_xi_discrepancy_bridge"] == AGENT3_XI_DISCREPANCY_BRIDGE
    assert contract["agent4_xi_moment_audit"] == AGENT4_XI_MOMENT_AUDIT

    assert PARENT_A5["head"] == "a732424661a615f24b7034f733bf5f281ac9a2d2"
    assert PARENT_A5["source_blob"] == "fc0f463a451fd26d0d7a5e90570a2c1698f725e6"
    assert PARENT_A5["test_blob"] == "8372a6a612839077d1d7e88f50dfe61edeab7c19"
    assert PARENT_A5["workflow_blob"] == "5391fe52e78b280628e7e62c71f7686842b1a3ef"

    assert AGENT1_XI_PREFIX_MOMENTS["head"] == "2797eed8e3b9b3374c676d8f517f05e3f0fc3e6b"
    assert AGENT1_XI_PREFIX_MOMENTS["source_blob"] == "92dce65c9a8793e06497a7e7347c832861032238"
    assert AGENT1_XI_PREFIX_MOMENTS["test_blob"] == "169a6282e1fd0e20225c669ff4473958d6e8e014"
    assert AGENT1_XI_PREFIX_MOMENTS["workflow_blob"] == "58a1b0c07cd7137d943a437d78f8db223d584652"

    assert AGENT3_XI_DISCREPANCY_BRIDGE["head"] == "87c1ef765d6bf3613c7c6c301669ffccf5260117"
    assert AGENT3_XI_DISCREPANCY_BRIDGE["source_blob"] == "223058e2132bb22298c6a642b0f2a55562181e25"
    assert AGENT3_XI_DISCREPANCY_BRIDGE["test_blob"] == "875bd551924a988b6adbc8f8f98180cddd78b04c"
    assert AGENT3_XI_DISCREPANCY_BRIDGE["workflow_blob"] == "23493753b08408826d54464a2646f95a3f8c33aa"

    assert AGENT4_XI_MOMENT_AUDIT["head"] == "d366a8aa3b36bb41825e5f5927548cc8bab58eaf"
    assert AGENT4_XI_MOMENT_AUDIT["source_blob"] == "b503283bfe7f5e36a4f1da82d4fec559fcf3e171"
    assert AGENT4_XI_MOMENT_AUDIT["test_blob"] == "6bdb06c8fbeec5a42b1af5091c46cc7f4837a61d"
    assert AGENT4_XI_MOMENT_AUDIT["workflow_blob"] == "3b3af059ac84e3ae6314a8aba25f9e19ca5d6fd0"


def test_current_xi_source_moment_discrepancy_is_real_but_not_complete_ns_defect() -> None:
    contract = build_contract(HEAD)
    a1 = contract["agent1_xi_prefix_moments"]
    a3 = contract["agent3_xi_discrepancy_bridge"]
    seam = contract["xi_moment_discrepancy_ingest_seam"]
    truth = contract["truth_boundary"]

    assert a1["X_i"] == 110.0
    assert a1["source_normalization_C"] == 1000.0
    assert a1["moment_names"] == ["M", "I", "J", "S", "C_p"]
    assert a1["candidate_side_actual_minus_ideal_discrepancy_materialized"] is True
    assert a1["G_i_ell_i_bound_into_pa16_shaped_input"] is True
    assert a1["pa16_repair_applied"] is False

    assert a3["upstream_agent1_head"] == a1["head"]
    assert a3["public_pa16_row_transform"] == ["M", "J-4*eta*I", "I", "S-8*eta*M", "C_p"]
    assert a3["current_candidate_side_discrepancy_consumed"] is True
    assert a3["pa16_inverse_or_fixed_point_solved"] is False
    assert a3["pa16_repair_applied"] is False
    assert a3["discrepancy_from_complete_ns_defect"] is False
    assert a3["authorized_for_gain_gated_ns_stage"] is False

    assert seam["registration_only"] is True
    assert seam["current_candidate_side_discrepancy_registered"] is True
    assert seam["agent3_row_transform_handoff_registered"] is True
    assert seam["source_moment_discrepancy_is_complete_ns_defect"] is False
    assert seam["source_moment_discrepancy_authorizes_correction_stage"] is False
    assert seam["row_transform_is_pa16_repair"] is False

    assert truth["current_candidate_side_xi_prefix_moment_discrepancy_materialized"] is True
    assert truth["discrepancy_from_complete_ns_defect"] is False
    assert truth["authorized_for_gain_gated_ns_stage"] is False
    assert truth["current_real_ns_defect_five_moment_discrepancy_materialized"] is False
    assert truth["pa16_repair_applied"] is False
    assert truth["five_moment_repair_closed"] is False


def test_matching_a4_audit_protocol_is_independent_but_not_admitted() -> None:
    contract = build_contract(HEAD)
    a4 = contract["agent4_xi_moment_audit"]
    seam = contract["xi_moment_discrepancy_ingest_seam"]

    assert a4["audited_agent1_head"] == AGENT1_XI_PREFIX_MOMENTS["head"]
    assert a4["saved_reloaded_public_profile_values_only"] is True
    assert a4["production_moment_integrator_used_as_reference"] is False
    assert "Simpson" in a4["independent_reference"]
    assert a4["protocol"]["seed"] == 9173631
    assert a4["protocol"]["fresh_off_grid_eta_count"] == 10
    assert a4["protocol"]["simpson_panels_per_public_segment"] == [48, 96, 192]
    assert a4["protocol"]["relative_rms_gate"] == 5.0e-3
    assert a4["protocol"]["relative_sampled_max_gate"] == 2.0e-2
    assert a4["protocol"]["refinement_ratio_gate"] == 6.0
    assert a4["protocol"]["refinement_floor"] == 5.0e-9
    assert a4["protocol"]["nontrivial_discrepancy_norm_floor"] == 1.0e-14
    assert a4["independent_xi_moment_audit_registered"] is True
    assert a4["independent_xi_moment_audit_admitted"] is False
    assert a4["complete_ns_residual_audit"] is False

    assert seam["matching_independent_a4_audit_protocol_registered"] is True
    assert seam["a4_exact_head_pass_required_before_audit_admission"] is True
    assert seam["queued_or_unresolved_ci_is_pass"] is False


def test_agent2_is_bound_as_sibling_without_crossing_lane_ownership() -> None:
    a2 = build_contract(HEAD)["agent2_sibling_status"]
    assert a2["pr"] == 948
    assert a2["head"] == "fc89770a9ee2ed1897e53ca83f57b995e7b1250d"
    assert a2["closes_current_xi_moment_seam"] is False
    assert a2["creates_correction_velocity"] is False
    assert a2["complete_ns_residual"] is False


def test_readiness_science_gates_and_st006_are_unchanged() -> None:
    contract = build_contract(HEAD)
    assert contract["readiness"] == READINESS == {
        "leading_ready": False,
        "oscillatory_ready": True,
        "correction_ready": False,
        "velocity_export_ready": False,
        "pde_validated": False,
    }
    assert contract["final_gate"] == FINAL_GATE
    assert FINAL_GATE["normalized_momentum_sampled_max"] == 1.0e-3
    assert FINAL_GATE["normalized_momentum_volume_l2"] == 1.0e-3
    assert FINAL_GATE["divergence_sampled_max"] == 1.0e-5
    assert FINAL_GATE["divergence_volume_l2"] == 1.0e-5
    assert contract["frozen_science"]["residual_defined_free_forcing_forbidden"] is True
    assert contract["st006_baseline"] == ST006_BASELINE
    assert ST006_BASELINE["momentum_sampled_max"] == pytest.approx(0.1082289305112118)
    assert ST006_BASELINE["momentum_volume_l2"] == pytest.approx(0.10758432876230622)


def test_truth_boundary_keeps_global_and_pde_promotions_fail_closed() -> None:
    truth = build_contract(HEAD)["truth_boundary"]
    assert truth == TRUTH_BOUNDARY
    for key in (
        "agent4_independent_xi_moment_audit_admitted",
        "upstream_ci_admitted_as_pass",
        "source_prepared_appendixA_pressure_materialized",
        "source_prepared_appendixA_full_reference_stress_pair_materialized",
        "source_admitted_global_kappa0_smallness",
        "source_T_sh_lower_bound_verified",
        "pa16_repair_applied",
        "five_moment_repair_closed",
        "discrepancy_from_complete_ns_defect",
        "authorized_for_gain_gated_ns_stage",
        "current_real_ns_defect_five_moment_discrepancy_materialized",
        "actual_inner_to_outer_global_join_materialized",
        "corrected_global_cartesian_leading_velocity_materialized",
        "matched_cartesian_pressure_materialized",
        "cartesian_matched_pressure_gradient_materialized",
        "restricted_forcing_materialized",
        "real_agent3_correction_velocity_materialized",
        "real_candidate_finite_correction_cycle_run",
        "complete_candidate_api_ready",
        "same_protocol_full_ns_residual_available",
        "current_candidate_eligible_for_full_ns_validation",
        "same_protocol_st006_comparison_available_now",
        "scientific_admission",
        "paper_exact",
    ):
        assert truth[key] is False, key


@pytest.mark.parametrize(
    ("path", "value"),
    [
        (("parent_a5", "head"), "f" * 40),
        (("agent1_xi_prefix_moments", "head"), "e" * 40),
        (("agent1_xi_prefix_moments", "source_blob"), "d" * 40),
        (("agent1_xi_prefix_moments", "X_i"), 109.999),
        (("agent1_xi_prefix_moments", "source_normalization_C"), 999.0),
        (("agent1_xi_prefix_moments", "pa16_repair_applied"), True),
        (("agent3_xi_discrepancy_bridge", "head"), "c" * 40),
        (("agent3_xi_discrepancy_bridge", "source_blob"), "b" * 40),
        (("agent3_xi_discrepancy_bridge", "public_pa16_row_transform", 1), "J+4*eta*I"),
        (("agent3_xi_discrepancy_bridge", "discrepancy_from_complete_ns_defect"), True),
        (("agent3_xi_discrepancy_bridge", "authorized_for_gain_gated_ns_stage"), True),
        (("agent4_xi_moment_audit", "head"), "a" * 40),
        (("agent4_xi_moment_audit", "protocol", "seed"), 1),
        (("agent4_xi_moment_audit", "protocol", "simpson_panels_per_public_segment"), [24, 48, 96]),
        (("agent4_xi_moment_audit", "protocol", "relative_rms_gate"), 6.0e-3),
        (("agent4_xi_moment_audit", "independent_xi_moment_audit_admitted"), True),
        (("xi_moment_discrepancy_ingest_seam", "queued_or_unresolved_ci_is_pass"), True),
        (("xi_moment_discrepancy_ingest_seam", "source_moment_discrepancy_is_complete_ns_defect"), True),
        (("xi_moment_discrepancy_ingest_seam", "source_moment_discrepancy_authorizes_correction_stage"), True),
        (("xi_moment_discrepancy_ingest_seam", "row_transform_is_pa16_repair"), True),
        (("final_gate", "normalized_momentum_volume_l2"), 1.001e-3),
        (("final_gate", "divergence_volume_l2"), 1.001e-5),
        (("truth_boundary", "agent4_independent_xi_moment_audit_admitted"), True),
        (("truth_boundary", "pa16_repair_applied"), True),
        (("truth_boundary", "discrepancy_from_complete_ns_defect"), True),
        (("truth_boundary", "actual_inner_to_outer_global_join_materialized"), True),
        (("truth_boundary", "matched_cartesian_pressure_materialized"), True),
        (("truth_boundary", "restricted_forcing_materialized"), True),
        (("truth_boundary", "real_agent3_correction_velocity_materialized"), True),
        (("truth_boundary", "same_protocol_full_ns_residual_available"), True),
        (("readiness", "leading_ready"), True),
        (("readiness", "correction_ready"), True),
        (("readiness", "pde_validated"), True),
    ],
)
def test_fail_closed_against_provenance_gate_and_truth_mutations(path: tuple[object, ...], value: object) -> None:
    mutated = copy.deepcopy(build_contract(HEAD))
    cursor = mutated
    for key in path[:-1]:
        cursor = cursor[key]
    cursor[path[-1]] = value
    with pytest.raises(ValueError):
        validate_contract(mutated)


def test_fail_closed_against_queued_evidence_laundering() -> None:
    for lane in (
        "agent1_xi_prefix_moments",
        "agent2_sibling_status",
        "agent3_xi_discrepancy_bridge",
        "agent4_xi_moment_audit",
    ):
        mutated = copy.deepcopy(build_contract(HEAD))
        evidence = mutated[lane]["observed_ci_at_registration"]["dedicated"]
        evidence["status"] = "completed"
        evidence["conclusion"] = "success"
        with pytest.raises(ValueError):
            validate_contract(mutated)


def test_contract_hash_and_round_trip_are_stable(tmp_path) -> None:
    first = build_contract(HEAD)
    second = build_contract(HEAD)
    assert first == second
    assert first["contract_sha256"] == second["contract_sha256"]
    other = build_contract("89abcdef0123456789abcdef0123456789abcdef")
    assert other["contract_sha256"] != first["contract_sha256"]
    validate_contract(first)

    path = tmp_path / "current_xi_moment_discrepancy_ingest.json"
    written = write_contract(path, exact_head=HEAD)
    assert json.loads(path.read_text(encoding="utf-8")) == written
    validate_contract(written)


def test_invalid_exact_head_is_rejected() -> None:
    with pytest.raises(ValueError):
        build_contract("not-a-sha")
