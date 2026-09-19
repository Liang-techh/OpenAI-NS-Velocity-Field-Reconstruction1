from __future__ import annotations

from openai_ns_reconstruction.kokuno_pa10_source_algebraic_slots_independent_audit import (
    run_independent_audit,
)


def test_algebraic_ordinary_slots_independent_preflight_passes() -> None:
    report = run_independent_audit()
    assert report["source_R2_algebraic_ordinary_slots_independent_preflight_passed"] is True
    assert report["failed_guards"] == []

    ratios = report["public_to_independent_ratios"]
    assert min(ratios["multipliers"].values()) >= 1.0
    for slot in ratios["numerator_slots"].values():
        assert slot["norm"] >= 1.0
        assert slot["lipschitz"] >= 1.0
    for slot in ratios["post_J1_slots"].values():
        assert slot["norm"] >= 1.0
        assert slot["lipschitz"] >= 1.0


def test_fresh_offgrid_axis_values_and_mutation() -> None:
    axis = run_independent_audit()["fresh_axis_value_audit"]
    assert axis["fresh_random_offgrid_points"] == 4096
    assert axis["value_path_passed"] is True
    assert max(axis["max_abs_error"].values()) <= 2e-15
    assert axis["value_path_mutation_detected"] is True
    assert axis["U_star_slope_4_to_3p99_mutation_response"] >= 5e-3


def test_preregistered_downward_mutations_are_detected() -> None:
    mutation = run_independent_audit()["mutation"]
    assert all(mutation["multiplier_0p999_detected"].values())
    assert all(mutation["numerator_0p99_detected"].values())
    assert all(mutation["post_J1_0p99_detected"].values())
    assert mutation["axis_value_slope_mutation_detected"] is True


def test_truth_boundary_remains_fail_closed() -> None:
    report = run_independent_audit()
    truth = report["truth_boundary"]
    assert truth["source_R2_algebraic_ordinary_slots_independently_audited"] is True
    assert truth["all_R2_ordinary_slots_source_bound"] is False
    assert truth["all_R2_ordinary_slots_independently_audited"] is False
    assert truth["all_R1_ordinary_slots_source_bound"] is False
    assert truth["source_operator_constant_M_machine_bound"] is False
    assert truth["source_operator_constant_K_machine_bound"] is False
    assert truth["global_pressure_matched"] is False
    assert truth["global_leading_profile_reconstructed"] is False
    assert truth["leading_only_ns_residual_assessed"] is False
    assert truth["leading_plus_oscillatory_ns_residual_assessed"] is False
    assert truth["after_correction_ns_residual_assessed"] is False
    assert truth["heldout_ns_momentum_residual_assessed"] is False
    assert truth["formal_full_domain_pde_gate_assessed"] is False
    assert truth["pde_validated"] is False
    assert report["latest_agent1_derivative_R2_subset_independently_audited_here"] is False


def test_final_project_gates_are_unchanged() -> None:
    gates = run_independent_audit()["final_project_gates_unchanged"]
    assert gates == {
        "normalized_momentum_max": 1e-3,
        "normalized_momentum_L2": 1e-3,
        "divergence_max": 1e-5,
        "divergence_L2": 1e-5,
    }
