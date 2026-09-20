from __future__ import annotations

import inspect

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_public_oscillatory_self_advection import (
    _manufactured_algebra_check,
    _self_advection_from_jacobian,
    evaluate_oscillatory_self_advection_fd6,
)
from openai_ns_reconstruction.kokuno_public_oscillatory_vorticity_diagnostic import (
    evaluate_vorticity_osc_fd6,
)


def test_public_api_has_no_residual_or_mean_correction_escape_hatch() -> None:
    signature = inspect.signature(evaluate_oscillatory_self_advection_fd6)
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
        "alpha",
        "normalized_score",
        "scientific_threshold",
        "delta_y",
        "delta_a",
    }
    assert forbidden.isdisjoint(signature.parameters)


def test_self_advection_tensor_orientation_is_component_by_derivative_axis() -> None:
    velocity = np.asarray(
        [[1.0, 2.0, 3.0], [-2.0, 0.5, 4.0]], dtype=float
    )
    jacobian = np.asarray(
        [
            [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]],
            [[0.5, -1.0, 2.0], [3.0, 0.0, -2.0], [1.5, 4.0, -0.5]],
        ],
        dtype=float,
    )
    observed = _self_advection_from_jacobian(velocity, jacobian)
    expected = np.stack(
        [jacobian[i] @ velocity[i] for i in range(velocity.shape[0])], axis=0
    )
    assert np.array_equal(observed, expected)


def test_manufactured_algebra_is_exact_to_roundoff() -> None:
    receipt = _manufactured_algebra_check()
    assert receipt["absolute_max_error"] <= 1.0e-14


def test_public_self_advection_exactly_replays_existing_fd6_gradient() -> None:
    x = np.asarray((0.44, 0.67, 0.91), dtype=float)
    y = np.asarray((0.13, -0.16, 0.19), dtype=float)
    z = np.asarray((-0.47, 0.06, 0.52), dtype=float)
    t = np.asarray((0.42, 0.50, 0.58), dtype=float)
    step = 0.001

    result = evaluate_oscillatory_self_advection_fd6(
        x, y, z, t, spatial_step=step
    )
    parent = evaluate_vorticity_osc_fd6(x, y, z, t, spatial_step=step)
    velocity = np.asarray(parent["velocity"], dtype=float)
    gradient = np.asarray(parent["velocity_gradient_fd6"], dtype=float)
    expected = np.einsum("...ij,...j->...i", gradient, velocity)

    assert result["velocity"].shape == (3, 3)
    assert result["velocity_gradient"].shape == (3, 3, 3)
    assert result["self_advection"].shape == (3, 3)
    assert np.array_equal(result["velocity"], velocity)
    assert np.array_equal(result["velocity_gradient"], gradient)
    assert np.array_equal(result["self_advection"], expected)
    assert np.all(np.isfinite(expected))
    assert np.max(np.abs(expected)) > 1.0e-8


def test_registered_support_exterior_is_exact_zero_for_field_and_quadratic_term() -> None:
    x = np.asarray((2.20, -2.24, 0.0, 0.0), dtype=float)
    y = np.asarray((0.0, 0.0, 2.22, -2.26), dtype=float)
    z = np.asarray((0.0, 0.0, 0.0, 0.0), dtype=float)
    t = np.full(x.shape, 0.50, dtype=float)
    result = evaluate_oscillatory_self_advection_fd6(
        x, y, z, t, spatial_step=0.001
    )
    assert np.array_equal(result["velocity"], np.zeros((4, 3)))
    assert np.array_equal(result["self_advection"], np.zeros((4, 3)))


def test_invalid_step_and_time_fail_closed() -> None:
    with pytest.raises(ValueError, match="spatial_step"):
        evaluate_oscillatory_self_advection_fd6(
            0.6, 0.0, 0.1, 0.5, spatial_step=0.0
        )
    with pytest.raises(ValueError, match="spatial_step"):
        evaluate_oscillatory_self_advection_fd6(
            0.6, 0.0, 0.1, 0.5, spatial_step=float("inf")
        )
    with pytest.raises(ValueError):
        evaluate_oscillatory_self_advection_fd6(
            0.6, 0.0, 0.1, 0.10, spatial_step=0.001
        )
