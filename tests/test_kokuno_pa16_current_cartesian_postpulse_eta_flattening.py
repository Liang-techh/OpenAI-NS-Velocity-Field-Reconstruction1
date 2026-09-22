from __future__ import annotations

import copy
import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_pa16_current_cartesian_postpulse_eta_flattening import (
    T_F_AUTONOMOUS,
    KokunoPA16CurrentCartesianPostPulseEtaFlattening,
    _smooth_step_and_derivative,
)


def _candidate() -> KokunoPA16CurrentCartesianPostPulseEtaFlattening:
    return KokunoPA16CurrentCartesianPostPulseEtaFlattening()


def test_public_flat_step_is_executable_and_has_flat_endpoints() -> None:
    s = np.array([-0.2, 0.0, 0.2, 0.5, 0.8, 1.0, 1.2])
    sigma, derivative = _smooth_step_and_derivative(s)
    np.testing.assert_array_equal(sigma[[0, 1]], np.zeros(2))
    np.testing.assert_array_equal(sigma[[-2, -1]], np.ones(2))
    np.testing.assert_array_equal(derivative[[0, 1, -2, -1]], np.zeros(4))
    assert sigma[2] < sigma[3] < sigma[4]
    assert sigma[3] == pytest.approx(0.5, rel=0.0, abs=2e-15)
    assert derivative[3] == pytest.approx(8.0, rel=0.0, abs=2e-14)

    h = 1.0e-6
    probes = np.array([0.2, 0.37, 0.5, 0.71, 0.8])
    sm, _ = _smooth_step_and_derivative(probes - h)
    sp, _ = _smooth_step_and_derivative(probes + h)
    _, analytic = _smooth_step_and_derivative(probes)
    np.testing.assert_allclose((sp - sm) / (2.0 * h), analytic, rtol=2e-7, atol=2e-10)


def test_exact_parent_replay_through_xi13() -> None:
    c = _candidate()
    eta = np.array([-0.43, 0.0, 0.39])
    for offset in (-120.0, -20.0, -1.0e-6, 0.0):
        log_X = c.log_X_pulse_end + offset
        got = c.similarity_profile_values_logX(np.full(eta.shape, log_X), eta)
        ref = c.parent.similarity_profile_values_logX(np.full(eta.shape, log_X), eta)
        for child_key, parent_key in (
            ("F_current_leading_postpulse_eta_flattening", "F_current_leading_with_pulse_end"),
            ("U_current_leading_postpulse_eta_flattening", "U_current_leading_with_pulse_end"),
            ("E_current_leading_postpulse_eta_flattening", "E_current_leading_with_pulse_end"),
            (
                "M_over_X_current_leading_postpulse_eta_flattening",
                "M_over_X_current_leading_with_pulse_end",
            ),
            (
                "M_eta_over_X_current_leading_postpulse_eta_flattening",
                "M_eta_over_X_current_leading_with_pulse_end",
            ),
            ("v0_current_leading_postpulse_eta_flattening", "v0_current_leading_with_pulse_end"),
        ):
            np.testing.assert_allclose(got[child_key], ref[parent_key], rtol=0.0, atol=0.0)


def test_postpulse_formula_sets_U_zero_but_retains_nontrivial_swirl() -> None:
    c = _candidate()
    eta = np.array([-0.62, -0.15, 0.0, 0.41, 0.88])
    log_X = np.full(eta.shape, c.log_X_pulse_end + 0.37 * c.T_f)
    p = c.similarity_profile_values_logX(log_X, eta)
    np.testing.assert_array_equal(
        p["U_current_leading_postpulse_eta_flattening"], np.zeros(eta.shape)
    )
    assert np.all(p["E_current_leading_postpulse_eta_flattening"] > 0.0)
    assert np.all(p["F_current_leading_postpulse_eta_flattening"] >= 0.0)
    assert np.max(p["E_current_leading_postpulse_eta_flattening"]) > 0.0
    assert np.ptp(p["E_current_leading_postpulse_eta_flattening"]) > 0.0


def test_eta_dependence_flattens_at_Tf_and_Tf_is_explicitly_autonomous() -> None:
    c = _candidate()
    assert c.T_f == T_F_AUTONOMOUS == 100.0
    eta = np.linspace(-1.0, 1.0,  nine := 9)
    end = c._postpulse_profile_logX(
        np.full(eta.shape, c.log_X_flatten_end), eta
    )
    E = np.asarray(end["E"], dtype=float)
    relative_spread = (np.max(E) - np.min(E)) / np.max(E)
    assert relative_spread < 2e-13
    np.testing.assert_allclose(end["sigma"], np.ones_like(eta), rtol=0.0, atol=0.0)
    np.testing.assert_allclose(end["E_eta"], np.zeros_like(eta), rtol=0.0, atol=0.0)

    truth = c.truth_boundary
    assert truth["repository_autonomous_T_f_materialized"] is True
    assert truth["source_exact_T_f_recovered"] is False
    assert truth["source_terminal_hold_after_eta_flattening_materialized"] is False
    assert truth["source_relative_swirl_bumps_materialized"] is False
    assert truth["source_exterior_heat_replacement_materialized"] is False
    assert truth["pde_validated"] is False


def test_public_ell_window_holds_on_deterministic_eta_y_grid() -> None:
    c = _candidate()
    eta = np.linspace(-1.0, 1.0, 41)
    y = np.linspace(0.0, c.T_f, 2001)
    yy, ee = np.meshgrid(y, eta, indexing="ij")
    p = c._postpulse_profile_logX(c.log_X_pulse_end + yy, ee)
    ell = np.asarray(p["ell"], dtype=float)
    assert np.min(ell) >= -c.lambda_value - 0.1 - 2e-14
    assert np.max(ell) <= -c.lambda_value + 2e-14
    assert np.min(ell) < -c.lambda_value - 0.05


def test_current_M_and_M_eta_primitive_identities_hold_with_public_U_zero() -> None:
    c = _candidate()
    eta = np.array([-0.37, 0.19, 0.71])
    h = 2.0e-6
    log_X = c.log_X_pulse_end + 43.0
    pm = c.similarity_profile_values_logX(np.full(eta.shape, log_X - h), eta)
    p0 = c.similarity_profile_values_logX(np.full(eta.shape, log_X), eta)
    pp = c.similarity_profile_values_logX(np.full(eta.shape, log_X + h), eta)
    for key in (
        "M_over_X_current_leading_postpulse_eta_flattening",
        "M_eta_over_X_current_leading_postpulse_eta_flattening",
    ):
        centered = (pp[key] - pm[key]) / (2.0 * h)
        np.testing.assert_allclose(centered, -p0[key], rtol=2e-7, atol=2e-300)


def test_analytic_logX_derivatives_match_independent_centered_difference() -> None:
    c = _candidate()
    eta = np.array([-0.53, -0.17, 0.26, 0.67])
    log_X = c.log_X_pulse_end + 38.0
    h = 1.0e-6
    pm = c.similarity_profile_values_logX(np.full(eta.shape, log_X - h), eta)
    pp = c.similarity_profile_values_logX(np.full(eta.shape, log_X + h), eta)
    d = c.similarity_log_radial_derivatives(np.full(eta.shape, log_X), eta)
    for value_key, deriv_key in (
        (
            "E_current_leading_postpulse_eta_flattening",
            "E_DlogX_current_leading_postpulse_eta_flattening",
        ),
        (
            "F_current_leading_postpulse_eta_flattening",
            "F_DlogX_current_leading_postpulse_eta_flattening",
        ),
        (
            "U_current_leading_postpulse_eta_flattening",
            "U_DlogX_current_leading_postpulse_eta_flattening",
        ),
    ):
        centered = (pp[value_key] - pm[value_key]) / (2.0 * h)
        np.testing.assert_allclose(centered, d[deriv_key], rtol=4e-7, atol=3e-300)
    centered_log_F = (
        pp["log_F_current_leading_postpulse_eta_flattening"]
        - pm["log_F_current_leading_postpulse_eta_flattening"]
    ) / (2.0 * h)
    np.testing.assert_allclose(
        centered_log_F,
        d["log_F_DlogX_current_leading_postpulse_eta_flattening"],
        rtol=3e-7,
        atol=5e-8,
    )


def test_cartesian_velocity_reaches_flattening_endpoint_and_axis_stays_regular() -> None:
    c = _candidate()
    radius = math.exp(c.log_radius_q1_flatten_end)
    v = c.velocity(
        np.array([radius, 0.0]),
        np.array([0.0, 0.0]),
        np.array([0.0, 0.0]),
        np.array([0.0, 0.0]),
    )
    assert v.shape == (2, 3)
    assert np.all(np.isfinite(v))
    np.testing.assert_array_equal(v[1, :2], np.zeros(2))
    assert abs(v[0, 2]) == 0.0
    assert np.linalg.norm(v[0, :2]) > 0.0

    coords = c.similarity_coordinates_logX(radius, 0.0, 0.0, 0.0)
    assert float(coords["log_X"]) == pytest.approx(c.log_X_flatten_end, rel=0.0, abs=4e-12)


def test_configuration_semantic_identity_and_truth_guards(tmp_path) -> None:
    c = _candidate()
    path = tmp_path / "postpulse_eta_flattening.json"
    payload = c.save_configuration(path)
    loaded = KokunoPA16CurrentCartesianPostPulseEtaFlattening.load_configuration(path)
    assert loaded.configuration() == payload
    assert loaded.semantic_sha256 == c.semantic_sha256

    mutated = copy.deepcopy(payload)
    mutated["bound_scope"]["T_f"] = 80.0
    with pytest.raises(ValueError):
        KokunoPA16CurrentCartesianPostPulseEtaFlattening.from_configuration(mutated)

    mutated = copy.deepcopy(payload)
    mutated["bound_scope"]["T_f_role"] = "source_exact"
    with pytest.raises(ValueError):
        KokunoPA16CurrentCartesianPostPulseEtaFlattening.from_configuration(mutated)

    mutated = copy.deepcopy(payload)
    mutated["truth_boundary"]["pde_validated"] = True
    with pytest.raises(ValueError):
        KokunoPA16CurrentCartesianPostPulseEtaFlattening.from_configuration(mutated)


def test_fail_closed_beyond_bounded_eta_flattening_stage() -> None:
    c = _candidate()
    with pytest.raises(ValueError):
        c.similarity_profile_values_logX(c.log_X_flatten_end + 1.0e-4, 0.0)


def test_deterministic_report_preserves_truth_boundary() -> None:
    c = _candidate()
    report = c.flattening_report(radial_samples=1001)
    assert report["T_f"] == 100.0
    assert report["min_ell_on_deterministic_grid"] >= -c.lambda_value - 0.1 - 2e-14
    assert report["max_ell_on_deterministic_grid"] <= -c.lambda_value + 2e-14
    assert report["endpoint_E_relative_eta_spread"] < 2e-13
    assert report["truth_boundary"]["pde_validated"] is False
