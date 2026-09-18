from openai_ns_reconstruction.kokuno_rational_rectangle_independent_audit import (
    FROZEN_GUARDS,
    SEED,
    run_audit,
)


def test_agent4_independent_rational_rectangle_audit_passes_frozen_guards():
    report = run_audit(seed=SEED)
    assert report["local_structural_preflight_passed"] is True
    assert report["formal_full_domain_pde_gate_assessed"] is False
    assert report["pde_validated"] is False
    assert len(report["cases"]) == 3

    summary = report["summary"]
    assert summary["source_delta_max_mismatches"] == 0
    assert summary["exact_minimum_mismatches"] == 0
    assert summary["C_J_relative_error_max"] <= FROZEN_GUARDS["C_J_relative_error_max"]
    assert summary["C_v_relative_error_max"] <= FROZEN_GUARDS["C_v_relative_error_max"]
    assert summary["collision_mutation_misses"] == 0
    assert summary["oversized_radius_mutation_misses"] == 0
    assert summary["truth_boundary_fail_closed"] is True


def test_fresh_cases_are_nontrivial_exact_and_mutation_sensitive():
    report = run_audit(seed=SEED)
    for case in report["cases"]:
        exact = case["independent_d_c_squared"]
        assert exact["numerator"] > 0
        assert exact["denominator"] > 0
        assert case["public_exact_minimum_match"] is True
        assert case["delta_max_public"] == case["delta_max_independent"]
        assert 0.0 < case["injectivity_left_side"] < 1.0
        assert 0.0 < case["separation_left_over_d_c"] < 1.0
        assert case["collision_mutation_detected"] is True
        assert case["oversized_radius_mutation_detected"] is True
        assert all(value is False for value in case["truth_boundary"].values())


def test_audit_is_deterministic_for_frozen_seed():
    first = run_audit(seed=SEED)
    second = run_audit(seed=SEED)
    assert first["sha256"] == second["sha256"]
    assert first["cases"] == second["cases"]
