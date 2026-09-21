from __future__ import annotations

import copy

import pytest

from openai_ns_reconstruction.kokuno_a5_fixed_kappa_x100_ingest_contract import (
    AGENT1_FIXED_KAPPA_X100,
    AGENT4_FIXED_KAPPA_X100_AUDIT,
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
    assert contract["agent1_fixed_kappa_x100"] == AGENT1_FIXED_KAPPA_X100
    assert contract["agent4_fixed_kappa_x100_audit"] == AGENT4_FIXED_KAPPA_X100_AUDIT
    assert PARENT_A5["head"] == "8bb1139f2b2d9824fd3bbac693ffa4bf2e71dea2"
    assert AGENT1_FIXED_KAPPA_X100["head"] == "1f3dec711f8ce2052276cd951b2f44d2849579bf"
    assert AGENT1_FIXED_KAPPA_X100["source_blob"] == "818d965d82bd92505bd6d3e805e698e1e0da64d8"
    assert AGENT1_FIXED_KAPPA_X100["test_blob"] == "dcf0e4b02000364eb77abb85751a4f76fecb9da7"
    assert AGENT1_FIXED_KAPPA_X100["workflow_blob"] == "87352d0a12e72a898f0a3b74bf98dccbd13c65b5"
    assert AGENT4_FIXED_KAPPA_X100_AUDIT["head"] == "602a8bb17d78421cddd493e2a6688a7e4c622e1e"
    assert AGENT4_FIXED_KAPPA_X100_AUDIT["source_blob"] == "d6df5ca24c92fe8e2a608dca983719711ca2f85d"
    assert AGENT4_FIXED_KAPPA_X100_AUDIT["test_blob"] == "95710d62c4889f4822dd80eb3b14912a3ca8e667"
    assert AGENT4_FIXED_KAPPA_X100_AUDIT["workflow_blob"] == "9693c6347d4904b63690db83e8b51e630e5b190b"


def test_fixed_kappa_candidate_and_audit_protocol_are_registered_without_promotion() -> None:
    contract = build_contract(HEAD)
    a1 = contract["agent1_fixed_kappa_x100"]
    a4 = contract["agent4_fixed_kappa_x100_audit"]
    seam = contract["fixed_kappa_x100_ingest_seam"]
    assert a1["continuation_range"] == [8.0, 100.0]
    assert a1["selected_kappa0"] == 0.1
    assert a1["continuation_quadrature_order"] == 48
    assert a1["public_values"] == ["F", "U", "E"]
    assert a1["public_radial_derivatives"] == ["F_X", "U_X", "E_X"]
    assert a1["source_prepared_appendixA_pressure"] is False
    assert a1["source_exact_fixed_kappa_continuation"] is False
    assert a4["saved_reloaded_public_values_only"] is True
    assert a4["production_derivatives_used_as_reference"] is False
    assert a4["protocol"]["seed"] == 9173601
    assert a4["protocol"]["relative_x_steps"] == [1.6e-3, 8.0e-4, 4.0e-4]
    assert a4["protocol"]["derivative_relative_rms_gate"] == 5.0e-3
    assert a4["protocol"]["derivative_relative_max_gate"] == 2.0e-2
    assert a4["protocol"]["refinement_ratio_gate"] == 20.0
    assert seam["registration_only"] is True
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


def test_truth_boundary_keeps_all_missing_global_objects_false() -> None:
    contract = build_contract(HEAD)
    truth = contract["truth_boundary"]
    assert truth == TRUTH_BOUNDARY
    for key in (
        "source_prepared_appendixA_pressure_materialized",
        "source_exact_fixed_kappa_continuation_materialized",
        "final_X100_to_X110_interpolation_materialized",
        "Xi110_handoff_materialized",
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
        (("agent1_fixed_kappa_x100", "selected_kappa0"), 0.2),
        (("agent1_fixed_kappa_x100", "head"), "f" * 40),
        (("agent1_fixed_kappa_x100", "source_blob"), "e" * 40),
        (("agent4_fixed_kappa_x100_audit", "independent_operator"), "production analytic derivative"),
        (("agent4_fixed_kappa_x100_audit", "protocol", "seed"), 1),
        (("final_gate", "normalized_momentum_volume_l2"), 1.001e-3),
        (("truth_boundary", "final_X100_to_X110_interpolation_materialized"), True),
        (("truth_boundary", "corrected_global_cartesian_leading_velocity_materialized"), True),
        (("truth_boundary", "matched_cartesian_pressure_materialized"), True),
        (("truth_boundary", "pde_validated"), True),
        (("fixed_kappa_x100_ingest_seam", "queued_or_unresolved_ci_is_pass"), True),
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
    evidence = mutated["agent4_fixed_kappa_x100_audit"]["observed_ci_at_registration"]["dedicated"]
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
