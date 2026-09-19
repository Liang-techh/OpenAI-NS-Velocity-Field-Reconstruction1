import inspect

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_public_oscillatory_laplacian_identity import (
    _fd4_laplacian_vector_identity,
    _manufactured_laplacian,
    _manufactured_provider,
    build_receipt,
    evaluate_laplacian_vector_identity,
)


def test_fd4_manufactured_vector_identity_laplacian() -> None:
    x = np.asarray([-0.41, 0.23, 0.52, -0.31])
    y = np.asarray([0.19, -0.37, 0.28, 0.44])
    z = np.asarray([0.33, -0.21, -0.46, 0.17])
    t = np.asarray([0.37, 0.45, 0.55, 0.63])
    _, _, lap = _fd4_laplacian_vector_identity(
        _manufactured_provider, x, y, z, t, spatial_step=0.013
    )
    exact = _manufactured_laplacian(x, y, z, t)
    np.testing.assert_allclose(lap, exact, rtol=5e-10, atol=5e-10)


def test_public_identity_batch_shapes_and_algebra() -> None:
    x = np.asarray([0.44, 0.63, 0.85])
    y = np.asarray([0.16, -0.22, 0.30])
    z = np.asarray([-0.39, 0.11, 0.55])
    t = np.asarray([0.37, 0.50, 0.63])
    result = evaluate_laplacian_vector_identity(
        x,
        y,
        z,
        t,
        first_derivative_step=0.009,
        laplacian_step=0.009,
    )
    for key in (
        "velocity",
        "grad_div_fd4",
        "curl_curl_fd4",
        "identity_laplacian_fd4",
        "minus_curl_curl_fd4",
        "laplacian_fd6",
        "identity_minus_fd6",
    ):
        assert np.asarray(result[key]).shape == (3, 3)
        assert np.all(np.isfinite(np.asarray(result[key])))
    np.testing.assert_allclose(
        np.asarray(result["identity_laplacian_fd4"]),
        np.asarray(result["grad_div_fd4"]) - np.asarray(result["curl_curl_fd4"]),
        rtol=0.0,
        atol=0.0,
    )
    np.testing.assert_allclose(
        np.asarray(result["identity_minus_fd6"]),
        np.asarray(result["identity_laplacian_fd4"])
        - np.asarray(result["laplacian_fd6"]),
        rtol=0.0,
        atol=0.0,
    )


def test_registered_support_exterior_nested_stencil_is_zero() -> None:
    result = evaluate_laplacian_vector_identity(
        np.asarray([0.0, 1.80, 0.65]),
        np.asarray([0.0, 0.0, 0.0]),
        np.asarray([0.0, 0.0, 2.35]),
        np.asarray([0.50, 0.50, 0.50]),
        first_derivative_step=0.018,
        laplacian_step=0.018,
    )
    for key in (
        "velocity",
        "grad_div_fd4",
        "curl_curl_fd4",
        "identity_laplacian_fd4",
        "laplacian_fd6",
        "identity_minus_fd6",
    ):
        assert float(np.max(np.abs(np.asarray(result[key])))) <= 1.0e-12


@pytest.mark.parametrize("step", [0.0, -0.01, 1e-7, 0.1, np.nan, np.inf])
def test_first_derivative_step_fails_closed(step: float) -> None:
    with pytest.raises(ValueError):
        evaluate_laplacian_vector_identity(
            0.5, 0.1, 0.0, 0.5, first_derivative_step=step
        )


@pytest.mark.parametrize("step", [0.0, -0.01, 1e-7, 0.1, np.nan, np.inf])
def test_laplacian_step_fails_closed(step: float) -> None:
    with pytest.raises(ValueError):
        evaluate_laplacian_vector_identity(
            0.5, 0.1, 0.0, 0.5, laplacian_step=step
        )


def test_public_signature_has_no_residual_or_correction_inputs() -> None:
    parameters = set(inspect.signature(evaluate_laplacian_vector_identity).parameters)
    forbidden = {
        "residual",
        "defect",
        "stress",
        "pressure",
        "forcing",
        "target",
        "gain",
        "correction",
        "threshold",
    }
    assert parameters.isdisjoint(forbidden)
    assert parameters == {
        "x",
        "y",
        "z",
        "t",
        "first_derivative_step",
        "laplacian_step",
    }


def test_preregistered_receipt_guards_and_truth_boundary() -> None:
    receipt = build_receipt()
    assert receipt["parent_agent2_pr"] == 708
    assert (
        receipt["parent_agent2_head"]
        == "9f8bae37c4d3b2559bee6db050655702e4c058e9"
    )
    assert receipt["velocity_pr"] == 561
    assert receipt["first_derivative_steps"] == [0.018, 0.009, 0.0045]
    assert receipt["comparison_laplacian_step"] == 0.0045
    assert receipt["failed_guards"] == []
    assert receipt["self_diagnostic_passed"] is True
    assert receipt["finest_identity_laplacian_rms"] > 1.0e-8
    assert receipt["support_exterior_absolute_max"] <= 1.0e-12

    truth = receipt["truth_boundary"]
    assert truth["velocity_candidate_changed"] is False
    assert truth["oscillatory_coefficients_retuned"] is False
    assert truth["source_formula_changed"] is False
    assert truth["laplacian_vector_identity_self_diagnostic_only"] is True
    assert truth["viscous_term_formula_or_pressure_fit_changed"] is False
    assert truth["independent_agent4_validation_replaced"] is False
    assert truth["independent_agent4_vector_potential_audit_required"] is True
    assert truth["agent3_momentum_defect_or_correction_duplicated"] is False
    assert truth["heldout_ns_momentum_residual_assessed"] is False
    assert truth["formal_full_domain_pde_gate_assessed"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False
