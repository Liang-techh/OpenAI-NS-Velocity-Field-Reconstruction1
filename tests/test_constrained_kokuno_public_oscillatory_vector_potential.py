import numpy as np

from openai_ns_reconstruction.kokuno_public_oscillatory_vector_potential import (
    FD4_STEPS,
    evaluate_vector_potential_osc,
    vector_potential_osc,
    verification_receipt,
)


def test_public_vector_potential_shapes_aggregation_and_support():
    x = np.asarray((0.62, -0.71, 0.0, 1.8, 0.8))
    y = np.asarray((0.17, 0.24, 0.0, 0.0, 0.0))
    z = np.asarray((0.35, -0.82, 0.0, 0.0, 2.1))
    t = np.asarray((0.33, 0.57, 0.50, 0.50, 0.50))
    out = evaluate_vector_potential_osc(x, y, z, t)

    total = np.asarray(out["vector_potential_cartesian_total"])
    by_beta = np.asarray(out["vector_potential_cartesian_by_beta"])
    by_sign = np.asarray(out["vector_potential_cartesian_by_beta_sign"])
    assert total.shape == (5, 3)
    assert by_beta.shape[:1] == (5,)
    assert by_sign.shape[:1] == (5,)
    assert by_sign.shape[-2:] == (2, 3)
    np.testing.assert_allclose(np.sum(by_sign, axis=-2), by_beta, rtol=0.0, atol=1.0e-12)
    np.testing.assert_allclose(np.sum(by_beta, axis=-2), total, rtol=0.0, atol=1.0e-12)
    assert np.linalg.norm(total[0]) > 0.0
    assert np.linalg.norm(total[1]) > 0.0
    np.testing.assert_array_equal(total[2:], np.zeros((3, 3)))
    np.testing.assert_array_equal(vector_potential_osc(0.0, 0.0, 0.0, 0.5), np.zeros(3))
    assert out["source_localized_potential_formula_reused"] is True
    assert out["velocity_candidate_changed"] is False
    assert out["paper_exact"] is False


def test_public_vector_potential_curl_matches_frozen_velocity_at_three_resolutions():
    receipt = verification_receipt()
    assert receipt["fd4_steps"] == list(FD4_STEPS)
    assert receipt["sample_count"] == 18
    assert receipt["velocity_rms"] > 0.0
    assert receipt["vector_potential_rms"] > 0.0
    assert receipt["support_exterior_absolute_max"] == 0.0

    rows = receipt["curl_comparison"]
    assert len(rows) == 3
    assert rows[-1]["relative_rms_error"] < 1.0e-4
    assert rows[-1]["relative_max_error"] < 2.0e-4
    assert rows[-1]["relative_rms_error"] < rows[0]["relative_rms_error"]
    ratios = receipt["rms_refinement_ratios"]
    assert len(ratios) == 2
    assert ratios[0] > 8.0
    assert ratios[1] > 8.0

    truth = receipt["truth_boundary"]
    assert truth["velocity_candidate_changed"] is False
    assert truth["source_formula_changed"] is False
    assert truth["public_vector_potential_materialized"] is True
    assert truth["independent_agent4_vector_potential_audit_required"] is True
    assert truth["heldout_ns_momentum_residual_assessed"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False


def test_public_vector_potential_rejects_out_of_interval_time():
    for bad in (0.249, 0.751):
        try:
            vector_potential_osc(0.7, 0.0, 0.0, bad)
        except ValueError as exc:
            assert "outside" in str(exc)
        else:
            raise AssertionError("out-of-interval time must fail closed")
