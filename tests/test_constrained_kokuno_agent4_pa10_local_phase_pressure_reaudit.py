import numpy as np

from openai_ns_reconstruction.kokuno_agent4_pa10_local_phase_pressure_reaudit import (
    AGENT1_CANDIDATE_SHA256,
    FINAL_PROJECT_GATES,
    _fd6_eta,
    run_audit,
)


class _PolynomialPublicPrimitive:
    def evaluate(self, y, eta):
        y_arr, eta_arr = np.broadcast_arrays(
            np.asarray(y, dtype=float), np.asarray(eta, dtype=float)
        )
        p = eta_arr**5 + 0.25 * y_arr
        return {"p": p}


def test_independent_fd6_operator_is_calibrated_on_polynomial():
    primitive = _PolynomialPublicPrimitive()
    y = np.asarray([0.7, 1.3, 2.1], dtype=float)
    eta = np.asarray([-0.31, 0.17, 0.43], dtype=float)
    estimate = _fd6_eta(primitive, "p", y, eta, 0.013)
    reference = 5.0 * eta**4
    np.testing.assert_allclose(estimate, reference, rtol=1e-10, atol=1e-12)


def test_reaudit_is_identity_bound_and_truth_boundary_stays_closed():
    report = run_audit()
    assert report["candidate_sha256"] == AGENT1_CANDIDATE_SHA256
    assert report["frozen_protocol"]["final_project_gates_unchanged"] == FINAL_PROJECT_GATES
    assert report["independence_contract"]["candidate_local_phase_quadrature_used_as_oracle"] is False
    assert report["independence_contract"]["candidate_internal_radial_integrals_used_as_oracle"] is False
    assert report["independence_contract"]["training_or_construction_tensor_read"] is False
    assert isinstance(report["repaired_public_p_eta_independent_preflight_passed"], bool)

    truth = report["truth_boundary"]
    assert truth["source_pressure_radius_one_ball_bound_assessed"] is False
    assert truth["global_pressure_matched"] is False
    assert truth["global_leading_profile_reconstructed"] is False
    assert truth["full_same_cycle_requested_stress_materialized"] is False
    assert truth["public_velocity_correction_materialized"] is False
    assert truth["heldout_ns_momentum_residual_assessed"] is False
    assert truth["formal_full_domain_pde_gate_assessed"] is False
    assert truth["pde_validated"] is False
