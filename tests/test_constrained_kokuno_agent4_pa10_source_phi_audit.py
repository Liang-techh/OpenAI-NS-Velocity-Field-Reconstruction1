import hashlib
import json

from openai_ns_reconstruction.kokuno_agent4_pa10_source_phi_audit import run_audit


def test_independent_source_phi_audit_passes_all_frozen_guards():
    report = run_audit()
    assert report["source_phi_ball_independent_preflight_passed"] is True
    assert report["failed_guards"] == []
    assert all(report["guards"].values())
    assert all(report["mutation"].values())


def test_independent_methods_and_offgrid_complex_stress_are_active():
    report = run_audit()
    protocol = report["frozen_protocol"]
    assert protocol["seed"] == 9173321
    assert protocol["complex_offgrid_samples"] == 20000
    assert protocol["independent_inverse_method"] == "mpmath besseli closed form"
    assert protocol["independent_weight_method"] == "mpmath infinite nsum"
    chi = report["independent_chi"]
    assert 0.0 < chi["direct_offgrid_max_H_variation"] <= chi["variation_bound"]
    assert 0.0 < chi["direct_offgrid_max_abs_chi"] <= chi["public_chi_complex_sup_upper"]


def test_public_bounds_dominate_independent_high_precision_reconstruction():
    report = run_audit()
    coeff = report["independent_coefficient_multiplier"]
    inv = report["independent_inverse_and_ball"]
    assert coeff["public_weight_sum_upper"] >= coeff["high_precision_infinite_weight_sum"]
    assert coeff["public_chi_multiplier_upper"] >= coeff["independent_chi_multiplier"]
    assert inv["public_inverse_upper"] >= inv["bessel_inverse_exact"]
    assert inv["public_Phi0_upper"] >= inv["independent_Phi0_norm"]
    assert inv["public_radius_one_ball_upper"] >= inv["independent_radius_one_ball_norm"]
    assert inv["public_Phi_lipschitz_upper"] == 1.0


def test_truth_boundary_stays_fail_closed_for_downstream_pde_claims():
    report = run_audit()
    truth = report["truth_boundary"]
    assert truth["source_axis_domain_independent_A4_admission_required_separately"] is True
    assert truth["source_Phi_ball_independently_audited"] is True
    assert truth["source_mixed_Y_Phi_Y_radius_one_ball_bound_machine_bound"] is False
    assert truth["source_mixed_Y_u_Y_radius_one_ball_bound_machine_bound"] is False
    assert truth["source_R1_R2_machine_bound"] is False
    assert truth["source_operator_M_K_machine_bound"] is False
    assert truth["global_pressure_matched"] is False
    assert truth["global_leading_profile_reconstructed"] is False
    assert truth["leading_only_ns_residual_assessed"] is False
    assert truth["leading_plus_oscillatory_ns_residual_assessed"] is False
    assert truth["after_correction_ns_residual_assessed"] is False
    assert truth["heldout_ns_momentum_residual_assessed"] is False
    assert truth["formal_full_domain_pde_gate_assessed"] is False
    assert truth["pde_validated"] is False


def test_receipt_hash_is_deterministic_and_project_gates_are_unchanged():
    first = run_audit()
    second = run_audit()
    assert first == second
    identity = dict(first)
    expected = identity.pop("receipt_sha256")
    raw = json.dumps(identity, sort_keys=True, separators=(",", ":"), allow_nan=False)
    assert hashlib.sha256(raw.encode("utf-8")).hexdigest() == expected
    gates = first["frozen_protocol"]["final_project_gates_unchanged"]
    assert gates == {
        "normalized_momentum_max": 1.0e-3,
        "normalized_momentum_L2": 1.0e-3,
        "divergence_max": 1.0e-5,
        "divergence_L2": 1.0e-5,
    }
