from __future__ import annotations

import inspect

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_public_signed_amplitude_complete_curl import (
    profile_from_delta_a,
)
from openai_ns_reconstruction.kokuno_public_signed_amplitude_laplacian_identity import (
    FD4_STEPS,
    FD6_COMPARISON_STEP,
    _verification_profile,
    evaluate_correction_laplacian_vector_identity,
    truth_boundary,
    verification_receipt,
)


def _profile(scale: float = 1.0):
    radii = np.linspace(0.33, 1.17, 23)
    s = (radii - radii[0]) / (radii[-1] - radii[0])
    bump = np.sin(np.pi * s) ** 8
    delta = np.stack((scale * 0.008 * bump, -scale * 0.0055 * bump), axis=-1)
    delta[[0, -1], :] = 0.0
    return profile_from_delta_a(
        radii,
        delta,
        reference_time=0.50,
        producer_kind="focused-regression",
        provenance="test-only signed radial correction",
        source_mean_amplitude_differential_certified=False,
    )


def test_public_identity_api_has_no_scientific_escape_hatch() -> None:
    signature = inspect.signature(evaluate_correction_laplacian_vector_identity)
    forbidden = {
        "residual", "defect", "stress", "pressure", "forcing", "target",
        "gain", "normalized_score", "scientific_threshold", "delta_y",
        "delta_a", "alpha", "damping_grid",
    }
    assert forbidden.isdisjoint(signature.parameters)


def test_zero_correction_replays_exact_zero_on_both_laplacian_routes() -> None:
    correction = _profile(scale=0.0)
    result = evaluate_correction_laplacian_vector_identity(
        correction,
        np.asarray((0.46, 0.71, 0.95)),
        np.asarray((0.12, -0.16, 0.14)),
        np.asarray((-0.39, 0.06, 0.47)),
        np.asarray((0.40, 0.50, 0.60)),
        first_derivative_step=FD4_STEPS[-1],
        laplacian_step=FD6_COMPARISON_STEP,
    )
    for key in (
        "velocity", "grad_div_fd4", "curl_curl_fd4",
        "identity_laplacian_fd4", "laplacian_fd6", "identity_minus_fd6",
    ):
        assert np.array_equal(np.asarray(result[key]), np.zeros((3, 3)))


def test_nonzero_identity_is_finite_and_closes_algebraically() -> None:
    correction = _verification_profile()
    result = evaluate_correction_laplacian_vector_identity(
        correction,
        np.asarray((0.45, 0.68, 0.94)),
        np.asarray((0.11, -0.13, 0.17)),
        np.asarray((-0.42, 0.09, 0.51)),
        np.asarray((0.40, 0.50, 0.60)),
        first_derivative_step=FD4_STEPS[-1],
        laplacian_step=FD6_COMPARISON_STEP,
    )
    grad_div = np.asarray(result["grad_div_fd4"])
    curl_curl = np.asarray(result["curl_curl_fd4"])
    identity = np.asarray(result["identity_laplacian_fd4"])
    assert np.all(np.isfinite(identity))
    assert np.array_equal(identity, grad_div - curl_curl)
    assert np.max(np.linalg.norm(identity, axis=-1)) > 1.0e-8


def test_step_and_time_domain_guards_fail_closed() -> None:
    correction = _profile()
    with pytest.raises(ValueError, match="first_derivative_step"):
        evaluate_correction_laplacian_vector_identity(
            correction, 0.6, 0.0, 0.1, 0.5, first_derivative_step=0.0
        )
    with pytest.raises(ValueError, match="laplacian_step"):
        evaluate_correction_laplacian_vector_identity(
            correction, 0.6, 0.0, 0.1, 0.5, laplacian_step=float("inf")
        )
    with pytest.raises(ValueError):
        evaluate_correction_laplacian_vector_identity(
            correction, 0.6, 0.0, 0.1, 0.10
        )


def test_truth_boundary_stays_fail_closed() -> None:
    truth = truth_boundary()
    assert truth["correction_laplacian_vector_identity_executable"] is True
    assert truth["fd6_parent_laplacian_cross_checked_by_independent_route"] is True
    assert truth["vector_potential_first_complete_curl_reused"] is True
    assert truth["candidate_amplitude_phase_support_retuned"] is False
    assert truth["agent3_mean_radial_chain_reimplemented"] is False
    assert truth["source_agent2_complete_curl_certified_for_full_candidate"] is False
    assert truth["independent_agent4_correction_vector_potential_audit_required"] is True
    assert truth["real_agent3_delta_a_bound"] is False
    assert truth["real_full_candidate_correction_cycle_run"] is False
    assert truth["heldout_ns_momentum_residual_assessed"] is False
    assert truth["residual_reduction_claimed"] is False
    assert truth["paper_exact"] is False
    assert truth["pde_validated"] is False


def test_preregistered_receipt_contract() -> None:
    receipt = verification_receipt()
    assert receipt["schema"] == "kokuno-a2-signed-amplitude-laplacian-vector-identity-v1"
    assert receipt["sample_count"] == 18
    assert receipt["fd4_steps"] == list(FD4_STEPS)
    assert receipt["fd6_comparison_step"] == FD6_COMPARISON_STEP
    assert receipt["support_exterior_absolute_max"] <= 1.0e-12
    assert receipt["manufactured_polynomial"]["relative_max_error"] <= 1.0e-8
    assert set(receipt["guards"]) == {
        "identity_laplacian_nontrivial_rms",
        "successive_difference_refinement",
        "identity_vs_fd6_relative_rms",
        "identity_vs_fd6_relative_max",
        "manufactured_relative_max_error",
        "support_exterior",
    }
    assert receipt["truth_boundary"]["heldout_ns_momentum_residual_assessed"] is False
