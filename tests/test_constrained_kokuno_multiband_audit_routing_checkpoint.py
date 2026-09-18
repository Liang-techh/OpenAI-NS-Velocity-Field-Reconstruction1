import copy

import pytest

from openai_ns_reconstruction.kokuno_multiband_audit_routing_checkpoint import (
    build_checkpoint,
    checkpoint_sha256,
    validate_checkpoint,
)


def test_checkpoint_is_deterministic_and_fail_closed():
    first = build_checkpoint()
    second = build_checkpoint()
    assert first == second
    validate_checkpoint(first)
    assert first["checkpoint_sha256"] == checkpoint_sha256(first)
    assert first["states"]["source_scheduled_multiband_supplied_mode_interface_ready"]
    assert not first["states"]["source_multiband_independent_cartesian_audit_assessed"]
    assert not first["states"]["source_multiband_independent_cartesian_audit_passed"]
    assert not first["states"]["oscillatory_ready"]
    assert not first["states"]["correction_ready"]
    assert not first["states"]["candidate_artifact_instantiated"]
    assert not first["states"]["pde_validated"]


def test_agent2_receipt_is_green_but_not_source_bound():
    payload = build_checkpoint()
    a2 = payload["upstream"]["agent2"]
    assert a2["pr"] == 450
    assert a2["head"] == "054d5b380a4ca2fcf5c2c8cf5e5a6532e82a5703"
    assert a2["dedicated_status"] == "success"
    assert a2["standard_status"] == "success"
    assert a2["source_scheduled_multiband_supplied_mode_family_executable"]
    assert a2["requires_two_distinct_ell_bands"]
    assert a2["per_band_Q_epsilon_schedule_bound"]
    assert a2["by_beta_cartesian_columns_exposed"]
    assert a2["by_band_cartesian_columns_exposed"]
    assert not a2["actual_positive_order_background_bound"]
    assert not a2["actual_auxiliary_torus_mode_family_bound"]
    assert not a2["public_source_bound_xyz_t_velocity_ready"]
    assert not a2["genuinely_independent_second_covariance_column_ready"]


def test_agent4_failure_is_contract_mismatch_not_scientific_verdict():
    payload = build_checkpoint()
    a4 = payload["upstream"]["agent4_failed_audit"]
    assert a4["pr"] == 451
    assert a4["dedicated_run"] == 35376636412
    assert a4["standard_run"] == 35376636366
    assert a4["dedicated_status"] == "failure"
    assert a4["standard_status"] == "failure"
    assert a4["focused_tests_passed_before_error"] == 7
    assert a4["focused_tests_errors"] == 3
    assert not a4["scientific_report_generated"]
    assert a4["failure_phase"] == "mutation_control_report_construction"
    assert a4["scientific_preflight_status"] == "not_assessed_due_to_contract_mismatch"
    assert "band aggregation does not reproduce" in a4["failure_exception"]


def test_agent4_frozen_guards_are_retained_verbatim():
    payload = build_checkpoint()
    guards = payload["agent4_frozen_local_guards"]
    assert guards["source_schedule_relative_error_max"] == 5.0e-14
    assert guards["finest_cartesian_curl_relative_rms_max"] == 2.0e-5
    assert guards["curl_refinement_ratio_min"] == 4.0
    assert guards["finest_normalized_divergence_rms_max"] == 2.0e-5
    assert guards["finest_normalized_divergence_point_max"] == 8.0e-5
    assert guards["divergence_refinement_ratio_min"] == 3.0
    assert guards["shared_first_band_schedule_mutation_relative_rms_min"] == 0.10
    assert guards["label_schedule_swap_mutation_relative_rms_min"] == 0.10
    assert not guards["changed_after_failure"]


def test_formal_gates_and_st006_baseline_are_unchanged():
    payload = build_checkpoint()
    gates = payload["formal_gates"]
    baseline = payload["baseline_vs_kokuno"]["st006"]
    assert gates["held_out_normalized_momentum_max"] == 1.0e-3
    assert gates["held_out_normalized_momentum_l2"] == 1.0e-3
    assert gates["held_out_divergence_max"] == 1.0e-5
    assert gates["held_out_divergence_l2"] == 1.0e-5
    assert not gates["changed_this_round"]
    assert baseline["momentum_sampled_max"] == 0.1082289305112118
    assert baseline["momentum_volume_l2"] == 0.10758432876230622
    assert payload["baseline_vs_kokuno"]["kokuno_current_comparable_full_domain_receipt"] is None


def test_validator_rejects_posthoc_guard_or_truth_promotion():
    payload = build_checkpoint()

    bad = copy.deepcopy(payload)
    bad["states"]["source_multiband_independent_cartesian_audit_passed"] = True
    bad["checkpoint_sha256"] = checkpoint_sha256(bad)
    with pytest.raises(ValueError, match="fail-closed"):
        validate_checkpoint(bad)

    bad = copy.deepcopy(payload)
    bad["upstream"]["agent4_failed_audit"]["scientific_report_generated"] = True
    bad["checkpoint_sha256"] = checkpoint_sha256(bad)
    with pytest.raises(ValueError, match="scientific assessment"):
        validate_checkpoint(bad)

    bad = copy.deepcopy(payload)
    bad["agent4_frozen_local_guards"]["label_schedule_swap_mutation_relative_rms_min"] = 0.0
    bad["checkpoint_sha256"] = checkpoint_sha256(bad)
    with pytest.raises(ValueError, match="frozen local guards"):
        validate_checkpoint(bad)

    bad = copy.deepcopy(payload)
    bad["formal_gates"]["held_out_normalized_momentum_max"] = 2.0e-3
    bad["checkpoint_sha256"] = checkpoint_sha256(bad)
    with pytest.raises(ValueError, match="formal PDE gates"):
        validate_checkpoint(bad)
