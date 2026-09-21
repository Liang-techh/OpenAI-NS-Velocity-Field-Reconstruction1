from __future__ import annotations

import copy
import json

import pytest

from openai_ns_reconstruction.kokuno_a5_actual_final_bridge_xi110_ingest_contract import (
    AGENT1_ACTUAL_FINAL_BRIDGE_XI110,
    AGENT4_ACTUAL_FINAL_BRIDGE_AUDIT,
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


def test_exact_parent_agent1_and_agent4_pins_are_frozen() -> None:
    contract = build_contract(HEAD)
    assert contract["parent_a5"] == PARENT_A5
    assert contract["agent1_actual_final_bridge_xi110"] == AGENT1_ACTUAL_FINAL_BRIDGE_XI110
    assert contract["agent4_actual_final_bridge_audit"] == AGENT4_ACTUAL_FINAL_BRIDGE_AUDIT

    assert PARENT_A5["head"] == "44b9349ef07d947d0da7eff053ae2409348374b8"
    assert PARENT_A5["source_blob"] == "531106d200b4163018f0331537b7bf6706042bd0"
    assert PARENT_A5["test_blob"] == "32864416a1270cfacb5cd896a36ecec28a64dda8"
    assert PARENT_A5["workflow_blob"] == "7fee1034269a6dd7cb2c183c2d5e16ec4c3a285d"

    assert AGENT1_ACTUAL_FINAL_BRIDGE_XI110["head"] == "f430f8f97df3bfafce56bf094047c769ecec922d"
    assert AGENT1_ACTUAL_FINAL_BRIDGE_XI110["source_blob"] == "0cbd49ad0726c1573f2e649c96a3b516d8df053e"
    assert AGENT1_ACTUAL_FINAL_BRIDGE_XI110["test_blob"] == "227d132a7ef1ea74323cfeb6c94f0deda490d3c1"
    assert AGENT1_ACTUAL_FINAL_BRIDGE_XI110["workflow_blob"] == "99c03716ac33a029931a67bc90d9511401963f89"

    assert AGENT4_ACTUAL_FINAL_BRIDGE_AUDIT["head"] == "d090572fdc766626c417cb53b8acc8632b2bdb63"
    assert AGENT4_ACTUAL_FINAL_BRIDGE_AUDIT["source_blob"] == "afe00af9274408051bb01ec9924ac287c538fe5f"
    assert AGENT4_ACTUAL_FINAL_BRIDGE_AUDIT["test_blob"] == "0a5de622f718321c74de8733718daf989928a97b"
    assert AGENT4_ACTUAL_FINAL_BRIDGE_AUDIT["workflow_blob"] == "2e96ff193b603cd283ac8e342434a0336a9df931"


def test_actual_bridge_and_matching_a4_protocol_are_registered_without_admission() -> None:
    contract = build_contract(HEAD)
    a1 = contract["agent1_actual_final_bridge_xi110"]
    a4 = contract["agent4_actual_final_bridge_audit"]
    seam = contract["actual_final_bridge_ingest_seam"]

    assert a1["source_coordinate_range"] == [100.0, 110.0]
    assert a1["X_i"] == 110.0
    assert a1["kappa_0"] == 0.1
    assert a1["axial_shutdown_log_width"] == 0.02
    assert a1["angular_settle_log_width"] == 0.02
    assert a1["bridge_quadrature_order"] == 48
    assert a1["quadrature_replay_order"] == 64
    assert a1["candidate_side_actual_final_interpolation_materialized"] is True
    assert a1["G_i_registered"] is True
    assert a1["ell_i_registered"] is True
    assert a1["semantic_identity_and_save_load_registered"] is True
    assert a1["global_cartesian_velocity_materialized"] is False
    assert a1["matched_cartesian_pressure_materialized"] is False
    assert a1["restricted_forcing_materialized"] is False

    assert a4["audited_agent1_head"] == a1["head"]
    assert a4["saved_reloaded_public_values_only"] is True
    assert a4["production_derivatives_used_as_reference"] is False
    assert a4["endpoint_safe_independent_method_registered"] is True
    assert a4["negative_result_receipt_retained_before_gate"] is True
    assert a4["protocol"]["seed"] == 9173621
    assert a4["protocol"]["random_per_region"] == 40
    assert a4["protocol"]["interior_sample_count"] == 130
    assert a4["protocol"]["endpoint_sample_count"] == 7
    assert a4["protocol"]["interior_halfwidths"] == [0.12, 0.06, 0.03]
    assert a4["protocol"]["endpoint_steps"] == [0.08, 0.04, 0.02]
    assert a4["protocol"]["derivative_relative_rms_gate"] == 5.0e-3
    assert a4["protocol"]["derivative_relative_max_gate"] == 2.0e-2
    assert a4["protocol"]["refinement_ratio_gate"] == 8.0
    assert a4["protocol"]["xi_final_slope_abs_gate"] == 2.0e-4
    assert a4["actual_final_bridge_independent_a4_audit_registered"] is True
    assert a4["actual_final_bridge_independent_a4_audit_admitted"] is False
    assert a4["actual_Xi_endpoint_derivative_independently_verified"] is False
    assert a4["complete_ns_residual_audit"] is False

    assert seam["registration_only"] is True
    assert seam["actual_candidate_source_coordinate_bridge_registered"] is True
    assert seam["G_i_ell_i_handoff_registered"] is True
    assert seam["matching_independent_a4_audit_protocol_registered"] is True
    assert seam["a4_exact_head_pass_required_before_audit_admission"] is True
    assert seam["queued_or_unresolved_ci_is_pass"] is False
    assert seam["profile_success_substitutes_for_cartesian_delivery"] is False
    assert seam["residual_defined_forcing_forbidden"] is True


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


def test_truth_boundary_registers_a4_method_but_not_a4_pass_or_global_science() -> None:
    truth = build_contract(HEAD)["truth_boundary"]
    assert truth == TRUTH_BOUNDARY
    assert truth["candidate_side_actual_final_bridge_100_to_Xi_registered"] is True
    assert truth["candidate_side_actual_G_i_at_Xi_registered"] is True
    assert truth["candidate_side_actual_ell_i_at_Xi_registered"] is True
    assert truth["semantic_identity_and_save_load_registered"] is True
    assert truth["actual_final_bridge_independent_a4_audit_registered"] is True
    assert truth["Xi_endpoint_independent_derivative_method_registered"] is True

    for key in (
        "upstream_ci_admitted_as_pass",
        "actual_final_bridge_independent_a4_audit_admitted",
        "actual_Xi_endpoint_derivative_independently_verified",
        "source_coordinate_profiles_are_cartesian_velocity",
        "source_prepared_appendixA_pressure_materialized",
        "source_prepared_appendixA_full_reference_stress_pair_materialized",
        "source_admitted_global_kappa0_smallness",
        "final_bridge_cone_admissibility_independently_certified",
        "actual_upstream_five_moment_discrepancy_materialized",
        "five_moment_repair_closed",
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
        (("agent1_actual_final_bridge_xi110", "head"), "e" * 40),
        (("agent1_actual_final_bridge_xi110", "source_blob"), "d" * 40),
        (("agent1_actual_final_bridge_xi110", "X_i"), 109.999),
        (("agent1_actual_final_bridge_xi110", "kappa_0"), 0.1001),
        (("agent1_actual_final_bridge_xi110", "axial_shutdown_log_width"), 0.0201),
        (("agent1_actual_final_bridge_xi110", "bridge_quadrature_order"), 32),
        (("agent4_actual_final_bridge_audit", "head"), "c" * 40),
        (("agent4_actual_final_bridge_audit", "source_blob"), "b" * 40),
        (("agent4_actual_final_bridge_audit", "protocol", "seed"), 1),
        (("agent4_actual_final_bridge_audit", "protocol", "endpoint_steps"), [0.1, 0.05, 0.025]),
        (("agent4_actual_final_bridge_audit", "protocol", "derivative_relative_rms_gate"), 6.0e-3),
        (("agent4_actual_final_bridge_audit", "actual_final_bridge_independent_a4_audit_admitted"), True),
        (("agent4_actual_final_bridge_audit", "actual_Xi_endpoint_derivative_independently_verified"), True),
        (("actual_final_bridge_ingest_seam", "queued_or_unresolved_ci_is_pass"), True),
        (("actual_final_bridge_ingest_seam", "profile_success_substitutes_for_cartesian_delivery"), True),
        (("final_gate", "normalized_momentum_volume_l2"), 1.001e-3),
        (("final_gate", "divergence_volume_l2"), 1.001e-5),
        (("truth_boundary", "actual_final_bridge_independent_a4_audit_admitted"), True),
        (("truth_boundary", "actual_Xi_endpoint_derivative_independently_verified"), True),
        (("truth_boundary", "actual_inner_to_outer_global_join_materialized"), True),
        (("truth_boundary", "corrected_global_cartesian_leading_velocity_materialized"), True),
        (("truth_boundary", "matched_cartesian_pressure_materialized"), True),
        (("truth_boundary", "restricted_forcing_materialized"), True),
        (("truth_boundary", "real_agent3_correction_velocity_materialized"), True),
        (("truth_boundary", "same_protocol_full_ns_residual_available"), True),
        (("truth_boundary", "scientific_admission"), True),
        (("readiness", "leading_ready"), True),
        (("readiness", "pde_validated"), True),
    ],
)
def test_fail_closed_against_provenance_gate_and_truth_mutations(path: tuple[str, ...], value: object) -> None:
    mutated = copy.deepcopy(build_contract(HEAD))
    cursor = mutated
    for key in path[:-1]:
        cursor = cursor[key]
    cursor[path[-1]] = value
    with pytest.raises(ValueError):
        validate_contract(mutated)


def test_fail_closed_against_queued_evidence_laundering() -> None:
    for lane in ("agent1_actual_final_bridge_xi110", "agent4_actual_final_bridge_audit"):
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

    path = tmp_path / "actual_final_bridge_ingest.json"
    written = write_contract(path, exact_head=HEAD)
    assert json.loads(path.read_text(encoding="utf-8")) == written
    validate_contract(written)


def test_invalid_exact_head_is_rejected() -> None:
    with pytest.raises(ValueError):
        build_contract("not-a-sha")
