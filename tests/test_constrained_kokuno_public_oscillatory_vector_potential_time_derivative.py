import numpy as np

from openai_ns_reconstruction.kokuno_public_oscillatory_vector_potential_time_derivative import (
    evaluate_vector_potential_osc_dt,
    vector_potential_osc_dt,
    verification_receipt,
)


def test_vector_potential_dt_public_shapes_and_support() -> None:
    x = np.asarray([0.55, 0.0, 1.8, 0.70])
    y = np.asarray([0.10, 0.0, 0.0, 0.05])
    z = np.asarray([0.20, 0.0, 0.10, 2.1])
    t = np.asarray([0.34, 0.50, 0.66, 0.50])

    evaluated = evaluate_vector_potential_osc_dt(x, y, z, t)
    total = evaluated["vector_potential_dt_cartesian_total"]
    by_beta = evaluated["vector_potential_dt_cartesian_by_beta"]
    by_sign = evaluated["vector_potential_dt_cartesian_by_beta_sign"]

    assert total.shape == (4, 3)
    assert by_beta.shape[:1] == (4,)
    assert by_sign.shape[:1] == (4,)
    np.testing.assert_allclose(total, np.sum(by_beta, axis=-2), rtol=0.0, atol=1e-11)
    np.testing.assert_allclose(by_beta, np.sum(by_sign, axis=-2), rtol=0.0, atol=1e-11)
    assert np.linalg.norm(total[0]) > 0.0
    assert np.all(total[1:] == 0.0)
    assert evaluated["velocity_candidate_changed"] is False
    assert evaluated["vector_potential_candidate_changed"] is False
    assert evaluated["source_formula_changed"] is False
    assert evaluated["paper_exact"] is False
    assert evaluated["pde_validated"] is False


def test_vector_potential_dt_scalar_and_broadcast_interface() -> None:
    scalar = vector_potential_osc_dt(0.62, 0.08, 0.31, 0.47)
    assert scalar.shape == (3,)
    assert np.all(np.isfinite(scalar))
    assert np.linalg.norm(scalar) > 0.0

    x = np.asarray([[0.62], [0.74]])
    y = np.asarray([0.08, -0.11, 0.15])
    out = vector_potential_osc_dt(x, y, 0.31, 0.47)
    assert out.shape == (2, 3, 3)
    assert np.all(np.isfinite(out))


def test_vector_potential_dt_cross_checks_and_truth_boundary() -> None:
    receipt = verification_receipt()
    assert receipt["parent_agent2_head"] == "aa79cb5eb601ad179c79fb8a4f3e21a093d94c31"
    assert receipt["admitted_agent2_head"] == "732800ce4990464b49c8aa32d0dff4580f6684d4"
    assert receipt["public_time_derivative_head"] == "6c8e71c800a17c0e9df1802042f9feebf9dbce04"
    assert receipt["sample_count"] == 18
    assert receipt["failed_guards"] == []
    assert min(receipt["time_rms_refinement_ratios"]) >= 20.0
    assert min(receipt["curl_dt_rms_refinement_ratios"]) >= 20.0
    assert receipt["time_fd6_comparison"][-1]["relative_rms_error"] <= 5e-8
    assert receipt["time_fd6_comparison"][-1]["relative_max_error"] <= 1e-7
    assert receipt["curl_dt_fd6_comparison"][-1]["relative_rms_error"] <= 5e-6
    assert receipt["curl_dt_fd6_comparison"][-1]["relative_max_error"] <= 1e-5
    assert receipt["support_exterior_absolute_max"] == 0.0

    truth = receipt["truth_boundary"]
    assert truth["public_vector_potential_time_derivative_materialized"] is True
    assert truth["agent2_self_check_only"] is True
    assert truth["independent_agent4_vector_potential_audit_required"] is True
    assert truth["velocity_candidate_changed"] is False
    assert truth["vector_potential_candidate_changed"] is False
    assert truth["source_formula_changed"] is False
    assert truth["heldout_ns_momentum_residual_assessed"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
