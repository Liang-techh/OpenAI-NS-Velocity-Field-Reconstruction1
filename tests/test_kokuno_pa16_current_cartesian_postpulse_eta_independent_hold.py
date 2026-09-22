from __future__ import annotations

import copy
import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_pa16_current_cartesian_postpulse_eta_independent_hold import (
    KokunoPA16CurrentCartesianPostPulseEtaIndependentHold,
)


def _candidate() -> KokunoPA16CurrentCartesianPostPulseEtaIndependentHold:
    return KokunoPA16CurrentCartesianPostPulseEtaIndependentHold()


def test_public_hold_length_is_exact_source_formula() -> None:
    c = _candidate()
    assert c.hold_length == pytest.approx(
        30.0 * math.log(1.0 / c.lambda_value), rel=0.0, abs=0.0
    )
    assert c.lambda_value == pytest.approx(0.05, rel=0.0, abs=1e-15)
    assert c.hold_length == pytest.approx(89.87196820661973, rel=0.0, abs=2e-14)


def test_exact_parent_replay_through_flattening_endpoint() -> None:
    c = _candidate()
    eta = np.array([-0.47, 0.0, 0.52])
    for offset in (-120.0, -20.0, -1.0e-6, 0.0):
        log_X = c.log_X_flatten_end + offset
        got = c.similarity_profile_values_logX(np.full(eta.shape, log_X), eta)
        ref = c.parent.similarity_profile_values_logX(np.full(eta.shape, log_X), eta)
        for child_key, parent_key in (
            (
                "F_current_leading_postpulse_hold",
                "F_current_leading_postpulse_eta_flattening",
            ),
            (
                "U_current_leading_postpulse_hold",
                "U_current_leading_postpulse_eta_flattening",
            ),
            (
                "E_current_leading_postpulse_hold",
                "E_current_leading_postpulse_eta_flattening",
            ),
            (
                "M_over_X_current_leading_postpulse_hold",
                "M_over_X_current_leading_postpulse_eta_flattening",
            ),
            (
                "M_eta_over_X_current_leading_postpulse_hold",
                "M_eta_over_X_current_leading_postpulse_eta_flattening",
            ),
            (
                "v0_current_leading_postpulse_hold",
                "v0_current_leading_postpulse_eta_flattening",
            ),
        ):
            np.testing.assert_allclose(
                got[child_key], ref[parent_key], rtol=0.0, atol=0.0
            )


def test_hold_keeps_U_zero_and_retains_nontrivial_swirl() -> None:
    c = _candidate()
    eta = np.array([-0.62, -0.15, 0.0, 0.41, 0.88])
    log_X = np.full(eta.shape, c.log_X_flatten_end + 0.43 * c.hold_length)
    p = c.similarity_profile_values_logX(log_X, eta)
    np.testing.assert_array_equal(
        p["U_current_leading_postpulse_hold"], np.zeros(eta.shape)
    )
    assert np.all(p["E_current_leading_postpulse_hold"] > 0.0)
    assert np.all(p["F_current_leading_postpulse_hold"] >= 0.0)
    assert np.max(p["E_current_leading_postpulse_hold"]) > 0.0

    truth = c.truth_boundary
    assert truth["source_terminal_hold_after_eta_flattening_materialized"] is True
    assert truth["current_cartesian_postpulse_eta_independent_hold_composed"] is True
    assert truth["source_relative_swirl_bumps_materialized"] is False
    assert truth["source_exterior_heat_replacement_materialized"] is False
    assert truth["pde_validated"] is False


def test_hold_preserves_eta_flattening_up_to_parent_float_roundoff() -> None:
    c = _candidate()
    eta = np.linspace(-1.0, 1.0, 11)
    for frac in (0.0, 0.31, 0.73, 1.0):
        p = c._hold_profile_logX(
            np.full(eta.shape, c.log_X_flatten_end + frac * c.hold_length), eta
        )
        E = np.asarray(p["E"], dtype=float)
        relative_spread = (np.max(E) - np.min(E)) / np.max(E)
        assert relative_spread < 3e-13
        np.testing.assert_array_equal(p["E_eta"], np.zeros_like(eta))


def test_public_hold_has_exact_ell_minus_lambda() -> None:
    c = _candidate()
    eta = np.array([-0.83, -0.24, 0.0, 0.37, 0.91])
    for frac in (0.07, 0.33, 0.68, 0.94):
        d = c.similarity_log_radial_derivatives(
            np.full(eta.shape, c.log_X_flatten_end + frac * c.hold_length), eta
        )
        np.testing.assert_allclose(
            d["ell_current_leading_postpulse_hold"],
            -c.lambda_value,
            rtol=0.0,
            atol=0.0,
        )


def test_current_M_and_M_eta_primitive_identities_hold_with_public_U_zero() -> None:
    c = _candidate()
    eta = np.array([-0.37, 0.19, 0.71])
    h = 2.0e-6
    log_X = c.log_X_flatten_end + 0.47 * c.hold_length
    pm = c.similarity_profile_values_logX(np.full(eta.shape, log_X - h), eta)
    p0 = c.similarity_profile_values_logX(np.full(eta.shape, log_X), eta)
    pp = c.similarity_profile_values_logX(np.full(eta.shape, log_X + h), eta)
    for key in (
        "M_over_X_current_leading_postpulse_hold",
        "M_eta_over_X_current_leading_postpulse_hold",
    ):
        centered = (pp[key] - pm[key]) / (2.0 * h)
        np.testing.assert_allclose(centered, -p0[key], rtol=2e-7, atol=2e-300)


def test_analytic_logX_derivatives_match_independent_centered_difference() -> None:
    c = _candidate()
    eta = np.array([-0.53, -0.17, 0.26, 0.67])
    log_X = c.log_X_flatten_end + 0.39 * c.hold_length
    h = 1.0e-6
    pm = c.similarity_profile_values_logX(np.full(eta.shape, log_X - h), eta)
    pp = c.similarity_profile_values_logX(np.full(eta.shape, log_X + h), eta)
    d = c.similarity_log_radial_derivatives(np.full(eta.shape, log_X), eta)
    for value_key, deriv_key in (
        (
            "E_current_leading_postpulse_hold",
            "E_DlogX_current_leading_postpulse_hold",
        ),
        (
            "F_current_leading_postpulse_hold",
            "F_DlogX_current_leading_postpulse_hold",
        ),
        (
            "U_current_leading_postpulse_hold",
            "U_DlogX_current_leading_postpulse_hold",
        ),
    ):
        centered = (pp[value_key] - pm[value_key]) / (2.0 * h)
        np.testing.assert_allclose(centered, d[deriv_key], rtol=4e-7, atol=3e-300)

    centered_log_F = (
        pp["log_F_current_leading_postpulse_hold"]
        - pm["log_F_current_leading_postpulse_hold"]
    ) / (2.0 * h)
    np.testing.assert_allclose(
        centered_log_F,
        d["log_F_DlogX_current_leading_postpulse_hold"],
        rtol=3e-7,
        atol=5e-8,
    )


def test_cartesian_velocity_reaches_hold_endpoint_and_axis_stays_regular() -> None:
    c = _candidate()
    radius = math.exp(c.log_radius_q1_hold_end)
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
    assert float(coords["log_X"]) == pytest.approx(
        c.log_X_hold_end, rel=0.0, abs=4e-12
    )


def test_configuration_semantic_identity_and_truth_guards(tmp_path) -> None:
    c = _candidate()
    path = tmp_path / "postpulse_eta_independent_hold.json"
    payload = c.save_configuration(path)
    loaded = KokunoPA16CurrentCartesianPostPulseEtaIndependentHold.load_configuration(
        path
    )
    assert loaded.configuration() == payload
    assert loaded.semantic_sha256 == c.semantic_sha256

    mutated = copy.deepcopy(payload)
    mutated["bound_scope"]["hold_length"] *= 0.99
    with pytest.raises(ValueError):
        KokunoPA16CurrentCartesianPostPulseEtaIndependentHold.from_configuration(
            mutated
        )

    mutated = copy.deepcopy(payload)
    mutated["bound_scope"]["hold_length_role"] = "repository_tuned"
    with pytest.raises(ValueError):
        KokunoPA16CurrentCartesianPostPulseEtaIndependentHold.from_configuration(
            mutated
        )

    mutated = copy.deepcopy(payload)
    mutated["bound_scope"]["T_f_role"] = "source_exact"
    with pytest.raises(ValueError):
        KokunoPA16CurrentCartesianPostPulseEtaIndependentHold.from_configuration(
            mutated
        )

    mutated = copy.deepcopy(payload)
    mutated["truth_boundary"]["source_relative_swirl_bumps_materialized"] = True
    with pytest.raises(ValueError):
        KokunoPA16CurrentCartesianPostPulseEtaIndependentHold.from_configuration(
            mutated
        )

    mutated = copy.deepcopy(payload)
    mutated["truth_boundary"]["pde_validated"] = True
    with pytest.raises(ValueError):
        KokunoPA16CurrentCartesianPostPulseEtaIndependentHold.from_configuration(
            mutated
        )


def test_fail_closed_before_relative_swirl_region() -> None:
    c = _candidate()
    with pytest.raises(ValueError):
        c.similarity_profile_values_logX(c.log_X_hold_end + 1.0e-4, 0.0)


def test_deterministic_hold_report_preserves_truth_boundary() -> None:
    c = _candidate()
    report = c.hold_report()
    assert report["hold_length_public_formula"] == pytest.approx(
        30.0 * math.log(1.0 / c.lambda_value), rel=0.0, abs=0.0
    )
    assert report["start_E_relative_eta_spread"] < 3e-13
    assert report["end_E_relative_eta_spread"] < 3e-13
    assert report["ell_on_hold"] == -c.lambda_value
    assert report["max_abs_U_hold_endpoint"] == 0.0
    assert report["truth_boundary"][
        "source_terminal_hold_after_eta_flattening_materialized"
    ] is True
    assert report["truth_boundary"]["source_relative_swirl_bumps_materialized"] is False
    assert report["truth_boundary"]["pde_validated"] is False
