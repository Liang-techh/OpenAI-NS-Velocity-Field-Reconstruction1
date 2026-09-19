import math

from openai_ns_reconstruction.kokuno_carrier_resolved_agent4_audit import (
    AUDITED_AGENT2_HEAD,
    FD4_STEPS,
    GUARDS,
    PHASE_RESOLUTIONS,
    SEED,
    run_audit,
)


def test_protocol_identity_and_frozen_guards():
    assert AUDITED_AGENT2_HEAD == "92852046e6e1ec0a5889b53f989df3ea5d79cb3b"
    assert SEED == 9173241
    assert FD4_STEPS == (0.02, 0.01, 0.005)
    assert PHASE_RESOLUTIONS == (32, 64, 128)
    assert GUARDS["finest_relative_divergence_rms"] == 2.0e-5
    assert GUARDS["finest_relative_divergence_max"] == 1.0e-4
    assert GUARDS["minimum_divergence_refinement_ratio"] == 3.0
    assert GUARDS["minimum_covariance_rank_ratio"] == 2.0e-2


def test_black_box_audit_is_finite_and_fail_closed_for_full_pde_claims():
    report = run_audit()
    assert report["audited_agent2_head"] == AUDITED_AGENT2_HEAD
    assert report["protocol"]["scientific_guards_changed_from_agent4_543"] is False
    assert report["protocol"]["pressure_or_forcing_fit"] is False
    assert report["protocol"]["candidate_internal_derivative_oracle_used"] is False
    assert report["protocol"]["training_tensor_or_loss_used"] is False

    rows = report["divergence"]["rows"]
    assert [row["step"] for row in rows] == list(FD4_STEPS)
    for row in rows:
        for key in (
            "absolute_rms",
            "absolute_max",
            "relative_rms",
            "relative_max",
            "gradient_rms",
            "gradient_max",
        ):
            assert math.isfinite(row[key])
            assert row[key] >= 0.0

    cov = report["phase_mean_covariance"]
    assert [row["phase_resolution"] for row in cov["rows"]] == list(PHASE_RESOLUTIONS)
    for row in cov["rows"]:
        assert math.isfinite(row["rank_ratio"])
        assert row["rank_ratio"] >= 0.0

    assert report["divergence_mutation"]["caught"] is True
    assert report["nontriviality"]["passed"] is True
    assert report["parameter_perturbation"]["passed"] is True

    truth = report["truth_boundary"]
    assert truth["formal_full_domain_pde_gate_assessed"] is False
    assert truth["heldout_ns_residual_assessed"] is False
    assert truth["pde_validated"] is False
    assert truth["formal_momentum_normalized_max_l2_gate"] == 1.0e-3
    assert truth["formal_divergence_max_l2_gate"] == 1.0e-5

    expected_local = bool(
        report["divergence"]["local_divergence_passed"]
        and cov["covariance_rank_preflight_passed"]
        and report["support"]["radial_axis_support_passed"]
        and report["nontriviality"]["passed"]
        and report["parameter_perturbation"]["passed"]
        and report["divergence_mutation"]["caught"]
    )
    assert report["local_curl_covariance_preflight_passed"] is expected_local
    assert report["public_oscillatory_preflight_passed"] is bool(
        expected_local and report["project_support_preflight_passed"]
    )
