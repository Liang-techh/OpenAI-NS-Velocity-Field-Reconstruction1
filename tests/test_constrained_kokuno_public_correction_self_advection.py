from __future__ import annotations

import inspect

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_public_correction_self_advection import (
    _manufactured_algebra_check,
    _self_advection_from_jacobian,
    _verification_profile,
    evaluate_correction_self_advection_fd4,
)
from openai_ns_reconstruction.kokuno_public_signed_amplitude_complete_curl import (
    profile_from_delta_a,
)
from openai_ns_reconstruction.kokuno_public_signed_amplitude_spatial_jacobian import (
    evaluate_correction_velocity_jacobian_fd4,
)


def test_public_api_has_no_residual_mean_or_raw_delta_a_escape_hatch() -> None:
    signature = inspect.signature(evaluate_correction_self_advection_fd4)
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
        [[1.0, 2.0, -1.0], [-0.5, 0.7, 1.3]], dtype=float
    )
    jacobian = np.asarray(
        [
            [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]],
            [[0.5, -1.0, 2.0], [3.0, 0.0, -2.0], [1.5, 4.0, -0.5]],
        ],
        dtype=float,
    )
    actual = _self_advection_from_jacobian(velocity, jacobian)
    expected = np.stack([jacobian[i] @ velocity[i] for i in range(2)], axis=0)
    assert np.array_equal(actual, expected)


def test_manufactured_self_advection_algebra_is_exact_to_roundoff() -> None:
    receipt = _manufactured_algebra_check()
    assert max(receipt.values()) <= 1.0e-14


def test_public_self_advection_exactly_replays_existing_correction_jacobian_seam() -> None:
    correction = _verification_profile()
    x = np.asarray((0.47, 0.62, 0.79), dtype=float)
    y = np.asarray((0.13, -0.18, 0.20), dtype=float)
    z = np.asarray((-0.38, 0.06, 0.42), dtype=float)
    t = np.asarray((0.44, 0.50, 0.56), dtype=float)

    result = evaluate_correction_self_advection_fd4(
        correction, x, y, z, t, spatial_step=0.001
    )
    parent = evaluate_correction_velocity_jacobian_fd4(
        correction, x, y, z, t, spatial_step=0.001
    )
    velocity = np.asarray(parent["velocity"], dtype=float)
    jacobian = np.asarray(parent["jacobian"], dtype=float)
    expected = np.einsum("...ij,...j->...i", jacobian, velocity)

    assert np.array_equal(result["correction_velocity"], velocity)
    assert np.array_equal(result["correction_jacobian"], jacobian)
    assert np.array_equal(result["correction_divergence"], parent["divergence"])
    assert np.array_equal(result["correction_vorticity"], parent["vorticity"])
    assert np.array_equal(result["correction_self_advection"], expected)
    assert np.all(np.isfinite(result["correction_self_advection"]))


def test_zero_complete_curl_correction_gives_exact_zero_self_advection() -> None:
    radii = np.linspace(0.30, 1.20, 9)
    correction = profile_from_delta_a(
        radii,
        np.zeros((radii.size, 2), dtype=float),
        reference_time=0.50,
        producer_kind="zero-regression",
        provenance="exact zero signed profile for correction self-advection regression",
        source_mean_amplitude_differential_certified=False,
    )
    result = evaluate_correction_self_advection_fd4(
        correction,
        np.asarray((0.49, 0.72), dtype=float),
        np.asarray((0.12, -0.18), dtype=float),
        np.asarray((-0.22, 0.31), dtype=float),
        np.asarray((0.44, 0.56), dtype=float),
        spatial_step=0.001,
    )
    assert np.array_equal(result["correction_velocity"], np.zeros((2, 3)))
    assert np.array_equal(result["correction_jacobian"], np.zeros((2, 3, 3)))
    assert np.array_equal(result["correction_self_advection"], np.zeros((2, 3)))


def test_registered_support_exterior_is_exact_zero() -> None:
    correction = _verification_profile()
    result = evaluate_correction_self_advection_fd4(
        correction,
        np.asarray((2.20, -2.24, 0.0, 1.72), dtype=float),
        np.asarray((0.0, 0.0, 2.22, 1.72), dtype=float),
        np.asarray((0.0, 0.0, 0.0, 2.20), dtype=float),
        np.full(4, 0.50, dtype=float),
        spatial_step=0.001,
    )
    for key, shape in (
        ("correction_velocity", (4, 3)),
        ("correction_jacobian", (4, 3, 3)),
        ("correction_divergence", (4,)),
        ("correction_vorticity", (4, 3)),
        ("correction_self_advection", (4, 3)),
    ):
        assert np.array_equal(np.asarray(result[key]), np.zeros(shape))


def test_invalid_step_and_time_fail_closed() -> None:
    correction = _verification_profile()
    with pytest.raises(ValueError, match="spatial_step"):
        evaluate_correction_self_advection_fd4(
            correction, 0.6, 0.0, 0.1, 0.5, spatial_step=0.0
        )
    with pytest.raises(ValueError):
        evaluate_correction_self_advection_fd4(
            correction, 0.6, 0.0, 0.1, 0.10, spatial_step=0.001
        )


def test_bad_tensor_shapes_and_nonfinite_values_fail_closed() -> None:
    with pytest.raises(ValueError, match="velocity"):
        _self_advection_from_jacobian(np.zeros((2, 2)), np.zeros((2, 3, 3)))
    with pytest.raises(ValueError, match="jacobian"):
        _self_advection_from_jacobian(np.zeros((2, 3)), np.zeros((2, 3, 2)))
    bad = np.zeros((2, 3))
    bad[0, 0] = np.nan
    with pytest.raises(ValueError, match="finite"):
        _self_advection_from_jacobian(bad, np.zeros((2, 3, 3)))
