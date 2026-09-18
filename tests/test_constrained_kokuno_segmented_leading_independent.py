import numpy as np

from openai_ns_reconstruction.kokuno_segmented_leading_independent import (
    DIVERGENCE_GATE,
    FD_STEPS,
    MUTATION,
    RESIDUAL_GATE,
    run_audit,
)


def test_segmented_independent_audit_preserves_fixed_gates_and_calibrates_operator():
    report = run_audit()
    contract = report["fixed_project_contract"]
    assert contract["normalized_momentum_max_gate"] == RESIDUAL_GATE == 1.0e-3
    assert contract["normalized_momentum_l2_gate"] == RESIDUAL_GATE
    assert contract["divergence_max_gate"] == DIVERGENCE_GATE == 1.0e-5
    assert contract["divergence_l2_gate"] == DIVERGENCE_GATE
    assert contract["thresholds_changed"] is False
    assert report["independent_operator"]["steps"] == list(FD_STEPS)
    assert report["independent_operator"]["training_or_construction_derivatives_read"] is False
    assert report["independent_operator"]["training_loss_read"] is False

    mutation = report["mutation_calibration"]
    assert np.isclose(mutation["mean_divergence_shift"], MUTATION, atol=1.0e-10)
    assert mutation["max_divergence_shift_error"] <= 1.0e-8
    assert np.isclose(
        mutation["mean_pressure_gradient_x_shift"], MUTATION, atol=1.0e-10
    )
    assert mutation["max_pressure_gradient_shift_error"] <= 1.0e-8
    assert report["preregistered_local_guards"]["mutation_shift_error_le_1e-8"] is True


def test_reference_stage_uses_three_resolutions_without_promoting_local_result():
    report = run_audit()
    for region_name in ("off_grid", "axis_near"):
        region = report["reference_stage"][region_name]
        ladder = region["resolution_ladder"]
        assert [entry["step"] for entry in ladder] == list(FD_STEPS)
        assert all(entry["sample_count"] > 0 for entry in ladder)
        assert all(np.isfinite(entry["max_vector_residual"]) for entry in ladder)
        assert all(np.isfinite(entry["sample_l2_vector_residual"]) for entry in ladder)
        assert all(np.isfinite(entry["divergence_max_abs"]) for entry in ladder)
        assert all(np.isfinite(entry["divergence_sample_l2"]) for entry in ladder)
        assert region["fine_to_previous_relative_change"]["sample_l2_vector_residual"] >= 0.0
        assert region["fine_to_previous_relative_change"]["max_vector_residual"] >= 0.0

    truth = report["truth_boundary"]
    assert truth["formal_full_domain_pde_gate_assessed"] is False
    assert truth["pde_validated"] is False
    assert truth["global_pressure_available"] is False
    assert truth["complete_leading_oscillatory_correction_composite_available"] is False


def test_i2_public_float64_contract_is_independently_consistent_but_repair_invisible():
    report = run_audit()
    i2 = report["i2_public_contract"]
    assert i2["sample_count"] == 4
    assert i2["public_vs_independent_source_base_max_relative"] <= 5.0e-10
    assert report["preregistered_local_guards"][
        "public_vs_independent_I2_source_base_relative_le_5e-10"
    ] is True

    # Scientific rejection/limitation, not a software failure: the candidate-facing
    # float64 router rounds the independently certified heat correction away.
    assert i2["public_float64_equals_uncorrected_base_all_probes"] is True
    assert i2["public_float64_equals_uncorrected_base_fraction"] == 1.0
    assert i2["router_vs_uncorrected_base_max_abs"] == 0.0
    assert i2["public_float64_heat_repair_visible"] is False
    assert i2["module_level_after_heat_repair_attribution_available"] is False
    assert report["truth_boundary"]["i2_float64_public_repair_observable_at_probes"] is False
    assert report["preregistered_local_guards"]["all_local_implementation_guards_passed"] is True
