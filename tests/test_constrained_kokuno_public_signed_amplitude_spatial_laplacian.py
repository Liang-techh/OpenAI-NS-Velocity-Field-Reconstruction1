from __future__ import annotations

import inspect

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_public_signed_amplitude_complete_curl import (
    profile_from_delta_a,
)
from openai_ns_reconstruction.kokuno_public_signed_amplitude_spatial_laplacian import (
    FD6_STEPS,
    _verification_profile,
    evaluate_correction_velocity_laplacian_fd6,
    truth_boundary,
    verification_receipt,
)


def _profile(scale: float = 1.0):
    radii = np.linspace(0.31, 1.19, 21)
    s = (radii - radii[0]) / (radii[-1] - radii[0])
    bump = np.sin(np.pi * s) ** 8
    delta = np.stack((scale * 0.009 * bump, -scale * 0.006 * bump), axis=-1)
    delta[[0, -1], :] = 0.0
    return profile_from_delta_a(
        radii,
        delta,
        reference_time=0.50,
        producer_kind="focused-regression",
        provenance="test-only signed radial correction",
        source_mean_amplitude_differential_certified=False,
    )


def test_public_laplacian_api_has_no_scientific_escape_hatch() -> None:
    signature = inspect.signature(evaluate_correction_velocity_laplacian_fd6)
    forbidden = {
        "residual",
        "defect",
        "stress",
        "pressure",
        "forcing",
        "target",
        "gain",
        "normalized_score",
        "scientific_threshold",
        "delta_y",
    }
    assert forbidden.isdisjoint(signature.parameters)


def test_zero_correction_has_exact_zero_velocity_and_laplacian() -> None:
    correction = _profile(scale=0.0)
    x = np.asarray((0.43, 0.67, 0.91))
    y = np.asarray((0.11, -0.14, 0.18))
    z = np.asarray((-0.41, 0.07, 0.52))
    t = np.asarray((0.39, 0.50, 0.61))
    result = evaluate_correction_velocity_laplacian_fd6(
        correction, x, y, z, t, spatial_step=0.002
    )
    assert result["velocity"].shape == (3, 3)
    assert result["second_derivatives"].shape == (3, 3, 3)
    assert result["laplacian"].shape == (3, 3)
    assert np.array_equal(result["velocity"], np.zeros((3, 3)))
    assert np.array_equal(result["second_derivatives"], np.zeros((3, 3, 3)))
    assert np.array_equal(result["laplacian"], np.zeros((3, 3)))


def test_nonzero_laplacian_is_finite_and_is_trace_of_pure_seconds() -> None:
    correction = _verification_profile()
    x = np.asarray((0.44, 0.69, 0.96))
    y = np.asarray((0.10, -0.15, 0.12))
    z = np.asarray((-0.43, 0.08, 0.49))
    t = np.asarray((0.39, 0.50, 0.61))
    result = evaluate_correction_velocity_laplacian_fd6(
        correction, x, y, z, t, spatial_step=FD6_STEPS[-1]
    )
    second = np.asarray(result["second_derivatives"])
    laplacian = np.asarray(result["laplacian"])
    assert np.all(np.isfinite(second))
    assert np.all(np.isfinite(laplacian))
    assert np.array_equal(laplacian, np.sum(second, axis=-1))
    assert np.max(np.linalg.norm(laplacian, axis=-1)) > 1.0e-8


def test_step_and_time_domain_guards_fail_closed() -> None:
    correction = _profile()
    with pytest.raises(ValueError, match="spatial_step"):
        evaluate_correction_velocity_laplacian_fd6(
            correction, 0.6, 0.0, 0.1, 0.5, spatial_step=0.0
        )
    with pytest.raises(ValueError, match="spatial_step"):
        evaluate_correction_velocity_laplacian_fd6(
            correction, 0.6, 0.0, 0.1, 0.5, spatial_step=float("inf")
        )
    with pytest.raises(ValueError):
        evaluate_correction_velocity_laplacian_fd6(
            correction, 0.6, 0.0, 0.1, 0.10, spatial_step=0.002
        )


def test_truth_boundary_stays_fail_closed() -> None:
    truth = truth_boundary()
    assert truth["correction_velocity_laplacian_executable"] is True
    assert truth["viscous_term_ready_for_future_correction_diagnostics"] is True
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
    assert receipt["schema"] == "kokuno-a2-signed-amplitude-spatial-laplacian-v1"
    assert receipt["sample_count"] == 18
    assert receipt["fd6_steps"] == list(FD6_STEPS)
    assert receipt["support_exterior_absolute_max"] <= 1.0e-12
    assert receipt["manufactured_polynomial"]["relative_max_error"] <= 5.0e-8
    # Scientific refinement/nontriviality guards are intentionally evaluated by
    # the immutable receipt/CI and are not weakened here after observing them.
    assert set(receipt["guards"]) == {
        "laplacian_nontrivial_rms",
        "successive_difference_refinement",
        "manufactured_relative_max_error",
        "support_exterior",
    }
    assert receipt["truth_boundary"]["heldout_ns_momentum_residual_assessed"] is False
