from __future__ import annotations

import inspect

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_public_signed_amplitude_complete_curl import (
    profile_from_delta_a,
)
from openai_ns_reconstruction.kokuno_public_signed_amplitude_spatial_jacobian import (
    FD4_STEPS,
    _fd4_first_derivatives,
    _manufactured_jacobian,
    _manufactured_provider,
    _verification_profile,
    evaluate_correction_velocity_jacobian_fd4,
    truth_boundary,
    verification_receipt,
)


def _profile(scale: float = 1.0):
    radii = np.linspace(0.31, 1.19, 25)
    s = (radii - radii[0]) / (radii[-1] - radii[0])
    bump = np.sin(np.pi * s) ** 8
    delta = np.stack((scale * 0.0085 * bump, -scale * 0.0060 * bump), axis=-1)
    delta[[0, -1], :] = 0.0
    return profile_from_delta_a(
        radii,
        delta,
        reference_time=0.50,
        producer_kind="focused-regression",
        provenance="test-only signed radial correction",
        source_mean_amplitude_differential_certified=False,
    )


def test_public_jacobian_api_has_no_scientific_escape_hatch() -> None:
    signature = inspect.signature(evaluate_correction_velocity_jacobian_fd4)
    forbidden = {
        "residual",
        "defect",
        "mean",
        "stress",
        "inverse",
        "damping",
        "pressure",
        "forcing",
        "target",
        "gain",
        "normalized_score",
        "scientific_threshold",
        "delta_y",
        "delta_a",
    }
    assert forbidden.isdisjoint(signature.parameters)


def test_zero_correction_has_exact_zero_velocity_jacobian_divergence_vorticity() -> None:
    correction = _profile(scale=0.0)
    x = np.asarray((0.43, 0.67, 0.91))
    y = np.asarray((0.11, -0.14, 0.18))
    z = np.asarray((-0.41, 0.07, 0.52))
    t = np.asarray((0.41, 0.50, 0.59))
    result = evaluate_correction_velocity_jacobian_fd4(
        correction, x, y, z, t, spatial_step=0.002
    )
    assert result["velocity"].shape == (3, 3)
    assert result["jacobian"].shape == (3, 3, 3)
    assert result["divergence"].shape == (3,)
    assert result["vorticity"].shape == (3, 3)
    assert np.array_equal(result["velocity"], np.zeros((3, 3)))
    assert np.array_equal(result["jacobian"], np.zeros((3, 3, 3)))
    assert np.array_equal(result["divergence"], np.zeros(3))
    assert np.array_equal(result["vorticity"], np.zeros((3, 3)))


def test_jacobian_convention_trace_and_curl_are_consistent() -> None:
    correction = _verification_profile()
    x = np.asarray((0.44, 0.69, 0.96))
    y = np.asarray((0.10, -0.15, 0.12))
    z = np.asarray((-0.43, 0.08, 0.49))
    t = np.asarray((0.41, 0.50, 0.59))
    result = evaluate_correction_velocity_jacobian_fd4(
        correction, x, y, z, t, spatial_step=FD4_STEPS[-1]
    )
    jac = np.asarray(result["jacobian"])
    divergence = np.asarray(result["divergence"])
    vorticity = np.asarray(result["vorticity"])
    expected_vorticity = np.stack(
        (
            jac[..., 2, 1] - jac[..., 1, 2],
            jac[..., 0, 2] - jac[..., 2, 0],
            jac[..., 1, 0] - jac[..., 0, 1],
        ),
        axis=-1,
    )
    assert np.all(np.isfinite(jac))
    assert np.array_equal(divergence, np.trace(jac, axis1=-2, axis2=-1))
    assert np.array_equal(vorticity, expected_vorticity)
    assert np.max(np.abs(jac)) > 1.0e-8


def test_manufactured_fd4_jacobian_fixes_component_axis_orientation() -> None:
    x = np.asarray((-0.43, 0.27, 0.51, -0.35), dtype=float)
    y = np.asarray((0.21, -0.39, 0.31, 0.46), dtype=float)
    z = np.asarray((0.34, -0.23, -0.48, 0.19), dtype=float)
    t = np.asarray((0.41, 0.47, 0.53, 0.59), dtype=float)
    observed = _fd4_first_derivatives(
        _manufactured_provider, x, y, z, t, spatial_step=0.011
    )
    expected = _manufactured_jacobian(x, y, z, t)
    scale = max(float(np.max(np.abs(expected))), 1.0)
    assert np.max(np.abs(observed - expected)) / scale <= 5.0e-9


def test_step_and_time_domain_guards_fail_closed() -> None:
    correction = _profile()
    with pytest.raises(ValueError, match="spatial_step"):
        evaluate_correction_velocity_jacobian_fd4(
            correction, 0.6, 0.0, 0.1, 0.5, spatial_step=0.0
        )
    with pytest.raises(ValueError, match="spatial_step"):
        evaluate_correction_velocity_jacobian_fd4(
            correction, 0.6, 0.0, 0.1, 0.5, spatial_step=float("inf")
        )
    with pytest.raises(ValueError):
        evaluate_correction_velocity_jacobian_fd4(
            correction, 0.6, 0.0, 0.1, 0.10, spatial_step=0.002
        )


def test_truth_boundary_stays_fail_closed() -> None:
    truth = truth_boundary()
    assert truth["correction_velocity_spatial_jacobian_executable"] is True
    assert truth["convective_cross_term_derivatives_ready_for_future_diagnostics"] is True
    assert truth["divergence_and_vorticity_derived_from_same_jacobian"] is True
    assert truth["vector_potential_first_complete_curl_reused"] is True
    assert truth["candidate_amplitude_phase_support_retuned"] is False
    assert truth["agent3_mean_radial_chain_reimplemented"] is False
    assert truth["source_agent2_complete_curl_certified_for_full_candidate"] is False
    assert truth["independent_agent4_correction_vector_potential_audit_required"] is True
    assert truth["real_agent3_delta_a_bound"] is False
    assert truth["full_composite_velocity_available"] is False
    assert truth["real_full_candidate_correction_cycle_run"] is False
    assert truth["heldout_ns_momentum_residual_assessed"] is False
    assert truth["residual_reduction_claimed"] is False
    assert truth["paper_exact"] is False
    assert truth["pde_validated"] is False


def test_preregistered_receipt_contract() -> None:
    receipt = verification_receipt()
    assert receipt["schema"] == "kokuno-a2-signed-amplitude-spatial-jacobian-v1"
    assert receipt["sample_count"] == 18
    assert receipt["fd4_steps"] == list(FD4_STEPS)
    assert receipt["support_exterior_absolute_max"] <= 1.0e-12
    assert receipt["manufactured_polynomial"]["relative_max_error"] <= 5.0e-9
    assert set(receipt["guards"]) == {
        "jacobian_nontrivial_rms",
        "successive_difference_refinement",
        "divergence_identity_relative_rms",
        "manufactured_relative_max_error",
        "support_exterior",
    }
    assert receipt["truth_boundary"]["heldout_ns_momentum_residual_assessed"] is False
