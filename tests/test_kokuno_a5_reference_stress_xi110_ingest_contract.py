from __future__ import annotations

import copy

import pytest

from openai_ns_reconstruction.kokuno_a5_reference_stress_xi110_ingest_contract import (
    AGENT1_REFERENCE_STRESS_XI110,
    AGENT4_REFERENCE_STRESS_XI110_AUDIT,
    FINAL_GATE,
    PARENT_A5,
    READINESS,
    ST006_BASELINE,
    TRUTH_BOUNDARY,
    build_contract,
    validate_contract,
)

HEAD = "0123456789abcdef0123456789abcdef01234567"


def test_exact_parent_and_upstream_pins_are_frozen() -> None:
    contract = build_contract(HEAD)
    assert contract["parent_a5"] == PARENT_A5
    assert contract["agent1_reference_stress_xi110"] == AGENT1_REFERENCE_STRESS_XI110
    assert contract["agent4_reference_stress_xi110_audit"] == AGENT4_REFERENCE_STRESS_XI110_AUDIT
    assert PARENT_A5["head"] == "744a4fbc642b861f61da91de713dd06adc357f09"
    assert PARENT_A5["source_blob"] == "502e64e97903853a16025815eabe31bf8b2111d0"
    assert PARENT_A5["test_blob"] == "cb6ab8142651b1879bd1a4d597fb7f9e3f6f9266"
    assert PARENT_A5["workflow_blob"] == "57b90219fcc45be891b6a826bd53a0c908ca2a65"
    assert AGENT1_REFERENCE_STRESS_XI110["head"] == "15791891fbbcfdaa614930791a843983afde2b4d"
    assert AGENT1_REFERENCE_STRESS_XI110["source_blob"] == "d534da9efed25df9d281bd0c78fc18ec9b706c09"
    assert AGENT1_REFERENCE_STRESS_XI110["test_blob"] == "369f9201303375a446a3d1bf04d019f5c0029d29"
    assert AGENT1_REFERENCE_STRESS_XI110["workflow_blob"] == "1144c17afa4d918f213d5477cd51e2c8b8ba4ed5"
    assert AGENT4_REFERENCE_STRESS_XI110_AUDIT["head"] == "5580526b783500be3701fbddabd6d7c7f90a5a61"
    assert AGENT4_REFERENCE_STRESS_XI110_AUDIT["source_blob"] == "43b73888c2ceeeb3c406cfc3c884eb0935855118"
    assert AGENT4_REFERENCE_STRESS_XI110_AUDIT["test_blob"] == "b1be7ada5ccfad758692052398efd63e695b8206"
    assert AGENT4_REFERENCE_STRESS_XI110_AUDIT["workflow_blob"] == "de2285dfe7575285629f1ff66dc7ba72412c7497"


def test_reference_stress_and_independent_audit_are_registered_without_bridge_promotion() -> None:
    contract = build_contract(HEAD)
    a1 = contract["agent1_reference_stress_xi110"]
    a4 = contract["agent4_reference_stress_xi110_audit"]
    seam = contract["reference_stress_xi110_ingest_seam"]
    assert a1["reference_range"] == [100.0, 110.0]
    assert a1["X_i"] == 110.0
    assert a1["quadrature_order"] == 64
    assert a1["Xi_handoff_registered"] is True
    assert a1["source_prepared_appendixA_pressure"] is False
    assert a1["actual_candidate_fixed_kappa_FUE_bridge_100_to_Xi"] is False
    assert a1["actual_final_interpolation_100_to_Xi"] is False
    assert a4["saved_reloaded_public_values_only"] is True
    assert a4["production_derivatives_used_as_reference"] is False
    assert a4["actual_candidate_bridge_audit"] is False
    assert a4["protocol"]["seed"] == 9173611
    assert a4["protocol"]["random_offgrid_count"] == 128
    assert a4["protocol"]["total_sample_count"] == 138
    assert a4["protocol"]["physical_halfwidths"] == [0.24, 0.12, 0.06]
    assert a4["protocol"]["derivative_relative_rms_gate"] == 5.0e-3
    assert a4["protocol"]["derivative_relative_max_gate"] == 2.0e-2
    assert a4["protocol"]["refinement_ratio_gate"] == 20.0
    assert seam["registration_only"] is True
    assert seam["this_is_actual_candidate_fixed_kappa_FUE_bridge"] is False
    assert seam["queued_or_unresolved_ci_is_pass"] is False
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
    assert contract["st006_baseline"] == ST006_BASELINE
    assert ST006_BASELINE["momentum_sampled_max"] == pytest.approx(0.1082289305112118)
    assert ST006_BASELINE["momentum_volume_l2"] == pytest.approx(0.10758432876230622)


def test_truth_boundary_registers_only_reference_prerequisite() -> None:
    truth = build_contract(HEAD)["truth_boundary"]
    assert truth == TRUTH_BOUNDARY
    assert truth["autonomous_reference_stress_to_Xi_registered"] is True
    assert truth["reference_stress_Xi_handoff_registered"] is True
    assert truth["independent_reference_stress_Xi_audit_registered"] is True
    for key in (
        "upstream_ci_admitted_as_pass",
        "source_prepared_appendixA_pressure_materialized",
        "source_prepared_appendixA_full_reference_stress_pair_materialized",
        "actual_final_interpolation_100_to_Xi_materialized",
        "fixed_kappa_FUE_X100_to_Xi_bridge_materialized",
        "actual_G_i_at_Xi_materialized",
        "actual_ell_i_at_Xi_materialized",
        "actual_upstream_five_moment_discrepancy_materialized",
        "five_moment_repair_closed",
        "corrected_global_cartesian_leading_velocity_materialized",
        "matched_cartesian_pressure_materialized",
        "cartesian_matched_pressure_gradient_materialized",
        "restricted_forcing_materialized",
        "real_agent3_correction_velocity_materialized",
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
        (("agent1_reference_stress_xi110", "head"), "f" * 40),
        (("agent1_reference_stress_xi110", "source_blob"), "e" * 40),
        (("agent1_reference_stress_xi110", "X_i"), 109.999),
        (("agent1_reference_stress_xi110", "quadrature_order"), 32),
        (("agent1_reference_stress_xi110", "actual_candidate_fixed_kappa_FUE_bridge_100_to_Xi"), True),
        (("agent4_reference_stress_xi110_audit", "independent_operator"), "source ODE derivative"),
        (("agent4_reference_stress_xi110_audit", "protocol", "seed"), 1),
        (("agent4_reference_stress_xi110_audit", "protocol", "derivative_relative_rms_gate"), 6.0e-3),
        (("final_gate", "normalized_momentum_volume_l2"), 1.001e-3),
        (("final_gate", "divergence_volume_l2"), 1.001e-5),
        (("truth_boundary", "actual_final_interpolation_100_to_Xi_materialized"), True),
        (("truth_boundary", "fixed_kappa_FUE_X100_to_Xi_bridge_materialized"), True),
        (("truth_boundary", "corrected_global_cartesian_leading_velocity_materialized"), True),
        (("truth_boundary", "matched_cartesian_pressure_materialized"), True),
        (("truth_boundary", "restricted_forcing_materialized"), True),
        (("truth_boundary", "real_agent3_correction_velocity_materialized"), True),
        (("truth_boundary", "scientific_admission"), True),
        (("readiness", "pde_validated"), True),
        (("reference_stress_xi110_ingest_seam", "queued_or_unresolved_ci_is_pass"), True),
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
    mutated = copy.deepcopy(build_contract(HEAD))
    evidence = mutated["agent4_reference_stress_xi110_audit"]["observed_ci_at_registration"]["dedicated"]
    evidence["status"] = "completed"
    evidence["conclusion"] = "success"
    with pytest.raises(ValueError):
        validate_contract(mutated)


def test_contract_hash_is_stable_and_mutation_sensitive() -> None:
    first = build_contract(HEAD)
    second = build_contract(HEAD)
    assert first == second
    assert first["contract_sha256"] == second["contract_sha256"]
    other = build_contract("89abcdef0123456789abcdef0123456789abcdef")
    assert other["contract_sha256"] != first["contract_sha256"]
    validate_contract(first)


def test_invalid_exact_head_is_rejected() -> None:
    with pytest.raises(ValueError):
        build_contract("not-a-sha")
