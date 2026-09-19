import inspect

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_public_oscillatory_spatial_laplacian import (
    _fd6_second_derivatives,
    _manufactured_laplacian,
    _manufactured_provider,
    build_receipt,
    evaluate_velocity_laplacian_fd6,
)


def test_fd6_manufactured_polynomial_laplacian() -> None:
    x = np.asarray([-0.41, 0.23, 0.52, -0.31])
    y = np.asarray([0.19, -0.37, 0.28, 0.44])
    z = np.asarray([0.33, -0.21, -0.46, 0.17])
    t = np.asarray([0.36, 0.44, 0.56, 0.64])
    second = _fd6_second_derivatives(
        _manufactured_provider, x, y, z, t, spatial_step=0.013
    )
    lap = np.sum(second, axis=-1)
    exact = _manufactured_laplacian(x, y, z, t)
    np.testing.assert_allclose(lap, exact, rtol=5e-10, atol=5e-10)


def test_public_laplacian_batch_shapes_and_trace() -> None:
    x = np.asarray([0.42, 0.61, 0.83])
    y = np.asarray([0.17, -0.23, 0.31])
    z = np.asarray([-0.41, 0.08, 0.52])
    t = np.asarray([0.36, 0.50, 0.64])
    result = evaluate_velocity_laplacian_fd6(
        x, y, z, t, spatial_step=0.009
    )
    assert np.asarray(result["velocity"]).shape == (3, 3)
    assert np.asarray(result["second_derivatives"]).shape == (3, 3, 3)
    assert np.asarray(result["laplacian"]).shape == (3, 3)
    np.testing.assert_allclose(
        np.asarray(result["laplacian"]),
        np.sum(np.asarray(result["second_derivatives"]), axis=-1),
        rtol=0.0,
        atol=0.0,
    )
    assert np.all(np.isfinite(np.asarray(result["laplacian"])))


def test_registered_support_exterior_stencil_is_zero() -> None:
    result = evaluate_velocity_laplacian_fd6(
        np.asarray([0.0, 1.75, 0.65]),
        np.asarray([0.0, 0.0, 0.0]),
        np.asarray([0.0, 0.0, 2.30]),
        np.asarray([0.50, 0.50, 0.50]),
        spatial_step=0.018,
    )
    for key in ("velocity", "second_derivatives", "laplacian"):
        assert float(np.max(np.abs(np.asarray(result[key])))) <= 1.0e-12


@pytest.mark.parametrize("step", [0.0, -0.01, 1e-7, 0.1, np.nan, np.inf])
def test_spatial_step_fails_closed(step: float) -> None:
    with pytest.raises(ValueError):
        evaluate_velocity_laplacian_fd6(0.5, 0.1, 0.0, 0.5, spatial_step=step)


def test_public_signature_has_no_residual_or_correction_inputs() -> None:
    parameters = set(inspect.signature(evaluate_velocity_laplacian_fd6).parameters)
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
    assert parameters == {"x", "y", "z", "t", "spatial_step"}


def test_preregistered_receipt_guards_and_truth_boundary() -> None:
    receipt = build_receipt()
    assert receipt["parent_agent2_pr"] == 694
    assert (
        receipt["parent_agent2_head"]
        == "9c4d1feb514fe44ca16793130f7bd1acfda95e26"
    )
    assert receipt["velocity_pr"] == 561
    assert receipt["spatial_steps"] == [0.018, 0.009, 0.0045]
    assert receipt["failed_guards"] == []
    assert receipt["self_diagnostic_passed"] is True
    assert receipt["finest_laplacian_rms"] > 1.0e-8
    assert receipt["support_exterior_absolute_max"] <= 1.0e-12

    truth = receipt["truth_boundary"]
    assert truth["velocity_candidate_changed"] is False
    assert truth["oscillatory_coefficients_retuned"] is False
    assert truth["source_formula_changed"] is False
    assert truth["numerical_spatial_laplacian_diagnostic_only"] is True
    assert truth["viscous_term_ready_for_future_composite_diagnostics"] is True
    assert truth["independent_agent4_validation_replaced"] is False
    assert truth["independent_agent4_vector_potential_audit_required"] is True
    assert truth["agent3_momentum_defect_or_correction_duplicated"] is False
    assert truth["heldout_ns_momentum_residual_assessed"] is False
    assert truth["formal_full_domain_pde_gate_assessed"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False
