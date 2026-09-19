import hashlib
import json

from openai_ns_reconstruction.kokuno_agent4_pa10_mixed_jnu_independent_audit import (
    RANDOM_CASES,
    run_audit,
    save_report,
)


def test_independent_mixed_jnu_audit_passes_all_preregistered_guards():
    report = run_audit()
    assert report["mixed_jnu_independent_preflight_passed"] is True
    assert report["failed_guards"] == []
    assert all(report["guards"].values())


def test_uses_fresh_random_stress_and_different_constant_path():
    report = run_audit()
    assert report["protocol"]["random_cases"] == RANDOM_CASES == 64
    assert report["protocol"]["independent_pi_upper"] == "355/113"
    assert report["protocol"]["agent1_pi_upper_not_reused"] is True
    op = report["operator_constant"]
    assert op["public_product_constant"] > op["independent_tight_product_constant"]
    assert min(op["mixed_factor_ratios_extras_0_to_3"]) >= 1.0
    assert all(op["mixed_factor_downward_mutations_detected"])


def test_random_norm_and_lipschitz_bounds_have_mutation_sensitivity():
    stress = run_audit()["random_ball_stress"]
    assert stress["min_norm_ratio_public_over_independent"] >= 1.0
    assert stress["min_lipschitz_ratio_public_over_independent"] >= 1.0
    assert stress["downward_mutations_detected"] == stress["downward_mutations_total"]
    assert stress["downward_mutations_total"] >= 2 * RANDOM_CASES - 1


def test_conditional_r_shapes_are_checked_but_not_promoted():
    report = run_audit()
    for block in report["conditional_r1_r2_shapes"].values():
        assert min(block.values()) >= 1.0
    truth = report["truth_boundary"]
    assert truth["conditional_mixed_operator_algebra_independently_audited"] is True
    assert truth["source_u_ball_independently_audited"] is False
    assert truth["full_R1_R2_independently_audited"] is False
    assert truth["M_K_independently_audited"] is False
    assert truth["global_pressure_independently_audited"] is False
    assert truth["global_leading_velocity_independently_audited"] is False
    assert truth["heldout_ns_momentum_residual_assessed"] is False
    assert truth["pde_validated"] is False


def test_final_project_gates_are_not_retuned():
    protocol = run_audit()["protocol"]
    assert protocol["final_momentum_gate_unchanged"] == 1.0e-3
    assert protocol["final_divergence_gate_unchanged"] == 1.0e-5


def test_receipt_hash_and_save_are_deterministic(tmp_path):
    first = run_audit()
    second = run_audit()
    assert first == second
    identity = dict(first)
    expected = identity.pop("receipt_sha256")
    raw = json.dumps(identity, sort_keys=True, separators=(",", ":"), allow_nan=False)
    assert hashlib.sha256(raw.encode()).hexdigest() == expected
    path = tmp_path / "mixed_jnu_a4.json"
    saved = save_report(path)
    assert json.loads(path.read_text()) == saved
