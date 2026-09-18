from __future__ import annotations

import json

from openai_ns_reconstruction.kokuno_bounded_family_coordinate_independent import (
    BASE_HEAD,
    COVARIANCE_FINE_RELATIVE_RMS_GUARD,
    COVARIANCE_MIN_REFINEMENT_RATIO_GUARD,
    FORMAL_DIVERGENCE_GATE,
    FORMAL_MOMENTUM_GATE,
    LABEL_MISALIGNMENT_RESPONSE_RELATIVE_FLOOR,
    RESPONSE_STEPS,
    SCHEMA,
    SEED,
    generate_report,
)


def test_independent_report_is_deterministic_and_truthful():
    first = generate_report()
    second = generate_report()
    assert first == second
    assert first["schema"] == SCHEMA
    assert first["base_head"] == BASE_HEAD
    assert first["seed"] == SEED
    assert first["sample_contract"]["response_steps"] == list(RESPONSE_STEPS)
    assert first["sample_contract"]["uses_public_tangent_arrays"] is False
    assert first["sample_contract"]["uses_training_tensor_or_loss"] is False
    assert first["sample_contract"]["uses_pressure_or_forcing_fit"] is False
    assert first["structural_preflight_passed"] == all(first["checks"].values())

    truth = first["truth_boundary"]
    assert truth["actual_source_positive_order_background_bound"] is False
    assert truth["actual_source_xyz_t_oscillatory_velocity_ready"] is False
    assert truth["genuinely_independent_second_covariance_column_ready"] is False
    assert truth["formal_full_domain_pde_gate_assessed"] is False
    assert truth["pde_validated"] is False


def test_frozen_guards_and_mutation_are_reported_without_threshold_drift():
    report = generate_report()
    assert report["guards"]["covariance_fine_relative_rms"] == COVARIANCE_FINE_RELATIVE_RMS_GUARD
    assert report["guards"]["covariance_min_refinement_ratio"] == COVARIANCE_MIN_REFINEMENT_RATIO_GUARD
    assert report["guards"]["label_misalignment_response_relative_floor"] == LABEL_MISALIGNMENT_RESPONSE_RELATIVE_FLOOR
    assert report["guards"]["formal_momentum_gate_unchanged"] == FORMAL_MOMENTUM_GATE
    assert report["guards"]["formal_divergence_gate_unchanged"] == FORMAL_DIVERGENCE_GATE

    for direction in ("common_direction", "band_direction"):
        ladder = report[direction]["forward_response_ladder"]
        assert [item["step"] for item in ladder] == list(RESPONSE_STEPS)
        assert len(report[direction]["refinement_ratios"]) == 2

    observed = report["label_misalignment_band_response_relative_change"]
    assert report["checks"]["label_misalignment_mutation_detected"] == (
        observed >= LABEL_MISALIGNMENT_RESPONSE_RELATIVE_FLOOR
    )


def test_report_json_roundtrip_has_no_nonfinite_values():
    report = generate_report()
    encoded = json.dumps(report, allow_nan=False, sort_keys=True)
    decoded = json.loads(encoded)
    assert decoded["task_id"] == report["task_id"]
    assert decoded["checks"] == report["checks"]
