from __future__ import annotations

import inspect

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_public_oscillatory_correction_cross_advection import (
    _cross_advection_from_jacobians,
    _manufactured_algebra_check,
    _verification_profile,
    evaluate_oscillatory_correction_cross_advection,
)
from openai_ns_reconstruction.kokuno_public_oscillatory_vorticity_diagnostic import (
    evaluate_vorticity_osc_fd6,
)
from openai_ns_reconstruction.kokuno_public_signed_amplitude_complete_curl import (
    profile_from_delta_a,
)
from openai_ns_reconstruction.kokuno_public_signed_amplitude_spatial_jacobian import (
    evaluate_correction_velocity_jacobian_fd4,
)


def test_public_api_has_no_residual_mean_or_raw_delta_a_escape_hatch() -> None:
    signature = inspect.signature(evaluate_oscillatory_correction_cross_advection)
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


def test_cross_advection_tensor_orientation_is_component_by_derivative_axis() -> None:
    uo = np.asarray([[1.0, 2.0, -1.0], [-0.5, 0.7, 1.3]], dtype=float)
    uc = np.asarray([[0.4, -0.2, 0.9], [1.1, -0.8, 0.3]], dtype=float)
    jo = np.asarray(
        [
            [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]],
            [[0.5, -1.0, 2.0], [3.0, 0.0, -2.0], [1.5, 4.0, -0.5]],
        ],
        dtype=float,
    )
    jc = np.asarray(
        [
            [[-0.2, 0.4, 0.6], [1.1, -0.7, 0.3], [0.8, 0.2, -0.9]],
            [[0.6, 1.2, -0.4], [-0.3, 0.9, 1.5], [2.0, -0.5, 0.7]],
        ],
        dtype=float,
    )
    term_oc, term_co, cross = _cross_advection_from_jacobians(uo, jo, uc, jc)
    expected_oc = np.stack([jc[i] @ uo[i] for i in range(2)], axis=0)
    expected_co = np.stack([jo[i] @ uc[i] for i in range(2)], axis=0)
    assert np.array_equal(term_oc, expected_oc)
    assert np.array_equal(term_co, expected_co)
    assert np.array_equal(cross, expected_oc + expected_co)


def test_manufactured_mixed_algebra_is_exact_to_roundoff() -> None:
    receipt = _manufactured_algebra_check()
    assert max(receipt.values()) <= 1.0e-14


def test_public_cross_term_exactly_replays_existing_a2_jacobian_seams() -> None:
    correction = _verification_profile()
    x = np.asarray((0.46, 0.61, 0.78), dtype=float)
    y = np.asarray((0.14, -0.19, 0.21), dtype=float)
    z = np.asarray((-0.39, 0.07, 0.43), dtype=float)
    t = np.asarray((0.43, 0.50, 0.57), dtype=float)

    result = evaluate_oscillatory_correction_cross_advection(
        correction,
        x,
        y,
        z,
        t,
        oscillatory_spatial_step=0.001,
        correction_spatial_step=0.001,
    )
    osc = evaluate_vorticity_osc_fd6(x, y, z, t, spatial_step=0.001)
    corr = evaluate_correction_velocity_jacobian_fd4(
        correction, x, y, z, t, spatial_step=0.001
    )
    uo = np.asarray(osc["velocity"], dtype=float)
    jo = np.asarray(osc["velocity_gradient_fd6"], dtype=float)
    uc = np.asarray(corr["velocity"], dtype=float)
    jc = np.asarray(corr["jacobian"], dtype=float)
    expected_oc = np.einsum("...ij,...j->...i", jc, uo)
    expected_co = np.einsum("...ij,...j->...i", jo, uc)

    assert np.array_equal(result["oscillatory_velocity"], uo)
    assert np.array_equal(result["correction_velocity"], uc)
    assert np.array_equal(result["oscillatory_jacobian"], jo)
    assert np.array_equal(result["correction_jacobian"], jc)
    assert np.array_equal(result["oscillation_advects_correction"], expected_oc)
    assert np.array_equal(result["correction_advects_oscillation"], expected_co)
    assert np.array_equal(result["cross_advection"], expected_oc + expected_co)
    assert np.all(np.isfinite(result["cross_advection"]))


def test_zero_complete_curl_correction_gives_exact_zero_mixed_term() -> None:
    radii = np.linspace(0.32, 1.18, 9)
    correction = profile_from_delta_a(
        radii,
        np.zeros((radii.size, 2), dtype=float),
        reference_time=0.50,
        producer_kind="zero-regression",
        provenance="exact zero signed profile for mixed-advection regression",
        source_mean_amplitude_differential_certified=False,
    )
    result = evaluate_oscillatory_correction_cross_advection(
        correction,
        np.asarray((0.49, 0.72), dtype=float),
        np.asarray((0.12, -0.18), dtype=float),
        np.asarray((-0.22, 0.31), dtype=float),
        np.asarray((0.43, 0.57), dtype=float),
        oscillatory_spatial_step=0.001,
        correction_spatial_step=0.001,
    )
    assert np.array_equal(result["correction_velocity"], np.zeros((2, 3)))
    assert np.array_equal(result["correction_jacobian"], np.zeros((2, 3, 3)))
    assert np.array_equal(result["oscillation_advects_correction"], np.zeros((2, 3)))
    assert np.array_equal(result["correction_advects_oscillation"], np.zeros((2, 3)))
    assert np.array_equal(result["cross_advection"], np.zeros((2, 3)))


def test_registered_support_exterior_is_exact_zero_for_both_fields_and_cross_term() -> None:
    correction = _verification_profile()
    x = np.asarray((2.20, -2.24, 0.0, 1.72), dtype=float)
    y = np.asarray((0.0, 0.0, 2.22, 1.72), dtype=float)
    z = np.asarray((0.0, 0.0, 0.0, 2.20), dtype=float)
    t = np.full(4, 0.50, dtype=float)
    result = evaluate_oscillatory_correction_cross_advection(
        correction, x, y, z, t,
        oscillatory_spatial_step=0.001,
        correction_spatial_step=0.001,
    )
    for key in (
        "oscillatory_velocity",
        "correction_velocity",
        "oscillation_advects_correction",
        "correction_advects_oscillation",
        "cross_advection",
    ):
        assert np.array_equal(np.asarray(result[key]), np.zeros((4, 3)))


def test_invalid_steps_and_time_fail_closed() -> None:
    correction = _verification_profile()
    with pytest.raises(ValueError, match="oscillatory_spatial_step"):
        evaluate_oscillatory_correction_cross_advection(
            correction, 0.6, 0.0, 0.1, 0.5, oscillatory_spatial_step=0.0
        )
    with pytest.raises(ValueError, match="correction_spatial_step"):
        evaluate_oscillatory_correction_cross_advection(
            correction, 0.6, 0.0, 0.1, 0.5, correction_spatial_step=float("inf")
        )
    with pytest.raises(ValueError):
        evaluate_oscillatory_correction_cross_advection(
            correction,
            0.6,
            0.0,
            0.1,
            0.10,
            oscillatory_spatial_step=0.001,
            correction_spatial_step=0.001,
        )
