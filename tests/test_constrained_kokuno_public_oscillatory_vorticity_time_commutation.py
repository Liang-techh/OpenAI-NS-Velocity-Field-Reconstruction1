from __future__ import annotations

import inspect

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_public_oscillatory_vorticity_time_commutation import (
    _curl_fd4,
    _time_fd6_of_curl_fd4,
    build_receipt,
    evaluate_vorticity_time_commutation,
)


def _manufactured_velocity(x, y, z, t):
    xb, yb, zb, tb = np.broadcast_arrays(
        np.asarray(x, dtype=float),
        np.asarray(y, dtype=float),
        np.asarray(z, dtype=float),
        np.asarray(t, dtype=float),
    )
    return np.stack(
        (
            tb * yb**4 + zb,
            tb * zb**3 + xb,
            tb * xb**2 + yb,
        ),
        axis=-1,
    )


def _manufactured_velocity_dt(x, y, z, t):
    xb, yb, zb, _ = np.broadcast_arrays(
        np.asarray(x, dtype=float),
        np.asarray(y, dtype=float),
        np.asarray(z, dtype=float),
        np.asarray(t, dtype=float),
    )
    return np.stack((yb**4, zb**3, xb**2), axis=-1)


def _manufactured_curl_dt(x, y, z):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    z = np.asarray(z, dtype=float)
    return np.stack((-3.0 * z**2, -2.0 * x, -4.0 * y**3), axis=-1)


def test_fd4_curl_and_time_fd6_commute_on_manufactured_polynomial():
    x = np.asarray([0.31, -0.27, 0.44])
    y = np.asarray([-0.22, 0.35, 0.18])
    z = np.asarray([0.17, -0.41, 0.29])
    t = np.asarray([0.36, 0.50, 0.64])
    expected = _manufactured_curl_dt(x, y, z)

    curl_dt = _curl_fd4(
        _manufactured_velocity_dt, x, y, z, t, spatial_step=0.01
    )
    time_curl = _time_fd6_of_curl_fd4(
        _manufactured_velocity,
        x,
        y,
        z,
        t,
        spatial_step=0.01,
        time_step=0.01,
    )
    np.testing.assert_allclose(curl_dt, expected, rtol=2e-10, atol=2e-10)
    np.testing.assert_allclose(time_curl, expected, rtol=2e-10, atol=2e-10)


def test_public_commutation_diagnostic_is_batchable_finite_and_nontrivial():
    result = evaluate_vorticity_time_commutation(
        np.asarray([0.52, 0.71, -0.62]),
        np.asarray([0.23, -0.33, 0.41]),
        np.asarray([-0.44, 0.19, 0.79]),
        np.asarray([0.35, 0.50, 0.65]),
        spatial_step=0.0045,
        time_step=0.004,
    )
    for key in (
        "vorticity_dt_time_fd6_of_curl_fd4",
        "vorticity_dt_curl_fd4_of_velocity_dt",
        "commutator",
    ):
        value = np.asarray(result[key], dtype=float)
        assert value.shape == (3, 3)
        assert np.all(np.isfinite(value))
    assert np.linalg.norm(
        np.asarray(result["vorticity_dt_curl_fd4_of_velocity_dt"], dtype=float)
    ) > 0.0


def test_public_api_has_no_residual_stress_forcing_or_target_injection():
    parameters = set(inspect.signature(evaluate_vorticity_time_commutation).parameters)
    assert parameters == {"x", "y", "z", "t", "spatial_step", "time_step"}
    forbidden = {
        "residual",
        "defect",
        "stress",
        "forcing",
        "pressure",
        "target",
        "gain",
        "correction",
    }
    assert parameters.isdisjoint(forbidden)


def test_step_and_time_stencil_guards_fail_closed():
    with pytest.raises(ValueError, match="spatial_step"):
        evaluate_vorticity_time_commutation(
            0.5, 0.2, 0.1, 0.5, spatial_step=0.0, time_step=0.004
        )
    with pytest.raises(ValueError, match="time_step"):
        evaluate_vorticity_time_commutation(
            0.5, 0.2, 0.1, 0.5, spatial_step=0.0045, time_step=0.0
        )
    with pytest.raises(ValueError, match="registered candidate interval"):
        evaluate_vorticity_time_commutation(
            0.5, 0.2, 0.1, 0.26, spatial_step=0.0045, time_step=0.01
        )


def test_exterior_stencils_preserve_registered_zero_support():
    result = evaluate_vorticity_time_commutation(
        np.asarray([0.0, 1.60, 0.60]),
        np.asarray([0.0, 0.0, 0.0]),
        np.asarray([0.0, 0.0, 2.16]),
        np.asarray([0.50, 0.50, 0.50]),
        spatial_step=0.0045,
        time_step=0.004,
    )
    for key in (
        "vorticity_dt_time_fd6_of_curl_fd4",
        "vorticity_dt_curl_fd4_of_velocity_dt",
        "commutator",
    ):
        assert np.max(np.abs(np.asarray(result[key], dtype=float))) == 0.0


def test_receipt_preserves_frozen_candidate_and_truth_boundary():
    receipt = build_receipt()
    assert receipt["parent_agent2_pr"] == 669
    assert receipt["parent_agent2_head"] == "f40436f2eddb77db94b07921fefa1e09ab94b134"
    assert receipt["velocity_pr"] == 561
    assert receipt["velocity_dt_pr"] == 579
    assert receipt["time_steps"] == [0.016, 0.008, 0.004]
    assert receipt["spatial_stability_steps"] == [0.018, 0.009, 0.0045]
    assert receipt["self_diagnostic_passed"] == (receipt["failed_guards"] == [])

    truth = receipt["truth_boundary"]
    assert truth["velocity_candidate_changed"] is False
    assert truth["velocity_time_derivative_changed"] is False
    assert truth["oscillatory_coefficients_retuned"] is False
    assert truth["source_formula_changed"] is False
    assert truth["numerical_curl_time_commutation_diagnostic_only"] is True
    assert truth["independent_agent4_validation_replaced"] is False
    assert truth["independent_agent4_vector_potential_audit_required"] is True
    assert truth["agent3_momentum_defect_or_correction_duplicated"] is False
    assert truth["heldout_ns_momentum_residual_assessed"] is False
    assert truth["formal_full_domain_pde_gate_assessed"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
