from __future__ import annotations

import inspect

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_public_oscillatory_vorticity_diagnostic import (
    _fd6_velocity_gradient,
    build_receipt,
    evaluate_vorticity_osc_fd6,
)


def _polynomial_velocity(x, y, z, t):
    xb, yb, zb, tb = np.broadcast_arrays(
        np.asarray(x, dtype=float),
        np.asarray(y, dtype=float),
        np.asarray(z, dtype=float),
        np.asarray(t, dtype=float),
    )
    return np.stack(
        (
            xb**6 + 2.0 * yb,
            yb**5 + 3.0 * zb,
            zb**4 + tb * xb,
        ),
        axis=-1,
    )


def test_centered_fd6_gradient_is_exact_on_low_degree_polynomial():
    x = np.asarray([0.31, -0.27, 0.44])
    y = np.asarray([-0.22, 0.35, 0.18])
    z = np.asarray([0.17, -0.41, 0.29])
    t = np.asarray([0.34, 0.50, 0.66])
    grad = _fd6_velocity_gradient(
        _polynomial_velocity, x, y, z, t, spatial_step=0.01
    )

    expected = np.zeros_like(grad)
    expected[..., 0, 0] = 6.0 * x**5
    expected[..., 0, 1] = 2.0
    expected[..., 1, 1] = 5.0 * y**4
    expected[..., 1, 2] = 3.0
    expected[..., 2, 0] = t
    expected[..., 2, 2] = 4.0 * z**3
    np.testing.assert_allclose(grad, expected, rtol=2e-10, atol=2e-10)


def test_public_vorticity_diagnostic_is_batchable_finite_and_nontrivial():
    result = evaluate_vorticity_osc_fd6(
        np.asarray([0.51, 0.72, -0.63]),
        np.asarray([0.24, -0.31, 0.43]),
        np.asarray([-0.47, 0.16, 0.81]),
        np.asarray([0.34, 0.50, 0.66]),
        spatial_step=0.0045,
    )
    assert result["velocity"].shape == (3, 3)
    assert result["velocity_gradient_fd6"].shape == (3, 3, 3)
    assert result["vorticity_fd6"].shape == (3, 3)
    assert result["divergence_fd6"].shape == (3,)
    assert all(
        np.all(np.isfinite(np.asarray(result[key], dtype=float)))
        for key in (
            "velocity",
            "velocity_gradient_fd6",
            "vorticity_fd6",
            "divergence_fd6",
        )
    )
    assert np.linalg.norm(np.asarray(result["vorticity_fd6"])) > 0.0


def test_fd6_stencil_remains_zero_far_outside_registered_support():
    result = evaluate_vorticity_osc_fd6(
        np.asarray([0.0, 1.55, 0.60]),
        np.asarray([0.0, 0.0, 0.0]),
        np.asarray([0.0, 0.0, 2.15]),
        np.asarray([0.50, 0.50, 0.50]),
        spatial_step=0.018,
    )
    for key in ("velocity", "velocity_gradient_fd6", "vorticity_fd6", "divergence_fd6"):
        assert np.max(np.abs(np.asarray(result[key], dtype=float))) == 0.0


def test_public_api_has_no_residual_stress_forcing_or_target_injection():
    parameters = set(inspect.signature(evaluate_vorticity_osc_fd6).parameters)
    assert parameters == {"x", "y", "z", "t", "spatial_step"}
    forbidden = {"residual", "defect", "stress", "forcing", "pressure", "target", "gain"}
    assert parameters.isdisjoint(forbidden)


def test_step_guard_is_fail_closed():
    with pytest.raises(ValueError, match="spatial_step"):
        evaluate_vorticity_osc_fd6(0.5, 0.2, 0.1, 0.5, spatial_step=0.0)
    with pytest.raises(ValueError, match="spatial_step"):
        evaluate_vorticity_osc_fd6(0.5, 0.2, 0.1, 0.5, spatial_step=0.08)


def test_receipt_preserves_truth_boundary_and_three_resolution_ladders():
    receipt = build_receipt()
    assert receipt["parent_agent2_head"] == "dfea7e7297d07154adbb59c30e74c81410cf0d58"
    assert receipt["derivative_resolution_ladder"]["steps"] == [0.018, 0.009, 0.0045]
    assert [v["resolution"] for v in receipt["morphology_resolution_ladder"]] == [17, 25, 33]
    assert receipt["self_diagnostic_passed"] == (receipt["failed_guards"] == [])
    truth = receipt["truth_boundary"]
    assert truth["velocity_candidate_changed"] is False
    assert truth["oscillatory_coefficients_retuned"] is False
    assert truth["source_formula_changed"] is False
    assert truth["numerical_fd6_diagnostic_only"] is True
    assert truth["independent_agent4_validation_replaced"] is False
    assert truth["agent3_momentum_defect_or_correction_duplicated"] is False
    assert truth["heldout_ns_momentum_residual_assessed"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
