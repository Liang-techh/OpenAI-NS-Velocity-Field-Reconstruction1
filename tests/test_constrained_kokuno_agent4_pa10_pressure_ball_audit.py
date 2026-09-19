import hashlib
import json

from openai_ns_reconstruction.kokuno_agent4_pa10_pressure_ball_audit import (
    FINAL_PROJECT_GATES,
    run_audit,
)


def test_independent_source_weight_and_public_bound_audit_passes():
    report = run_audit()
    assert report["source_weight_audit"]["failures"] == []
    assert report["source_weight_audit"]["first_global_algebraic_guard"] is True
    assert report["source_weight_audit"]["second_global_factor_guard"] is True
    assert report["eta_factor_underbounds"] == []
    assert report["public_pressure_underbounds"] == []
    assert report["minimum_public_to_independent_upper_ratio"] >= 1.0
    assert report["pressure_ball_operator_algebra_independent_preflight_passed"] is True
    assert report["failed_guards"] == []


def test_one_percent_downward_mutation_is_detected_for_every_checked_bound():
    report = run_audit()
    mutation = report["downward_mutation"]
    assert mutation["detected"] == mutation["total"]
    assert mutation["total"] == report["case_count"] * 8


def test_truth_boundary_and_receipt_hash_are_fail_closed():
    report = run_audit()
    truth = report["truth_boundary"]
    assert truth["pressure_ball_operator_algebra_independently_audited"] is True
    assert truth["diagnostic_inputs_are_source_radius_one_ball_bounds"] is False
    assert truth["source_rho_machine_bound"] is False
    assert truth["source_g_coefficient_norm_machine_bound"] is False
    assert truth["source_pressure_radius_one_ball_norm_machine_bound"] is False
    assert truth["source_R1_R2_machine_bound"] is False
    assert truth["global_pressure_matched"] is False
    assert truth["global_leading_profile_reconstructed"] is False
    assert truth["heldout_ns_momentum_residual_assessed"] is False
    assert truth["formal_full_domain_pde_gate_assessed"] is False
    assert truth["pde_validated"] is False
    assert report["frozen_protocol"]["final_project_gates_unchanged"] == FINAL_PROJECT_GATES

    body = dict(report)
    expected = body.pop("receipt_sha256")
    raw = json.dumps(body, sort_keys=True, separators=(",", ":"), allow_nan=False)
    assert hashlib.sha256(raw.encode("utf-8")).hexdigest() == expected
