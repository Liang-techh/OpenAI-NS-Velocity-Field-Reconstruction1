from __future__ import annotations

import copy
import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_pa16_current_cartesian_postpulse_eta_independent_hold_prefix import (
    FIRST_RELATIVE_SWIRL_CENTER_FROM_END,
    PRE_BUMP_GUARD_FROM_HOLD_END_AUTONOMOUS,
    RELATIVE_SWIRL_WIDTH_PUBLIC,
    SECOND_RELATIVE_SWIRL_CENTER_FROM_END,
    KokunoPA16CurrentCartesianPostPulseEtaIndependentHoldPrefix,
)


def _candidate() -> KokunoPA16CurrentCartesianPostPulseEtaIndependentHoldPrefix:
    return KokunoPA16CurrentCartesianPostPulseEtaIndependentHoldPrefix()


def test_public_hold_length_and_relative_swirl_geometry_are_recorded() -> None:
    c = _candidate()
    expected = 30.0 * math.log(1.0 / c.lambda_value)
    assert c.hold_length == pytest.approx(expected, rel=0.0, abs=0.0)
    assert c.lambda_value == pytest.approx(0.05, rel=0.0, abs=1e-15)
    assert c.hold_length == pytest.approx(89.87196820661973, rel=0.0, abs=2e-14)
    assert RELATIVE_SWIRL_WIDTH_PUBLIC == 0.3
    assert FIRST_RELATIVE_SWIRL_CENTER_FROM_END == 3.0
    assert SECOND_RELATIVE_SWIRL_CENTER_FROM_END == 1.0
    assert c.first_relative_swirl_center_s == pytest.approx(
        c.hold_length - 3.0, rel=0.0, abs=0.0
    )
    assert c.second_relative_swirl_center_s == pytest.approx(
        c.hold_length - 1.0, rel=0.0, abs=0.0
    )
    assert PRE_BUMP_GUARD_FROM_HOLD_END_AUTONOMOUS == 3.5
    assert c.prefix_length == pytest.approx(
        c.hold_length - 3.5, rel=0.0, abs=0.0
    )
    assert c.prefix_length < c.first_relative_swirl_center_s - RELATIVE_SWIRL_WIDTH_PUBLIC


def test_exact_parent_replay_through_flattening_endpoint() -> None:
    c = _candidate()
    eta = np.array([-0.47, 0.0, 0.52])
    for offset in (-120.0, -20.0, -1.0e-6, 0.0):
        log_X = c.log_X_flatten_end + offset
        got = c.similarity_profile_values_logX(np.full(eta.shape, log_X), eta)
        ref = c.parent.similarity_profile_values_logX(np.full(eta.shape, log_X), eta)
        for child_key, parent_key in (
            (
                "F_current_leading_pre_relative_swirl",
                "F_current_leading_postpulse_eta_flattening",
            ),
            (
                "U_current_leading_pre_relative_swirl",
                "U_current_leading_postpulse_eta_flattening",
            ),
            (
                "E_current_leading_pre_relative_swirl",
                "E_current_leading_postpulse_eta_flattening",
            ),
            (
                "M_over_X_current_leading_pre_relative_swirl",
                "M_over_X_current_leading_postpulse_eta_flattening",
            ),
            (
                "M_eta_over_X_current_leading_pre_relative_swirl",
                "M_eta_over_X_current_leading_postpulse_eta_flattening",
            ),
            (
                "v0_current_leading_pre_relative_swirl",
                "v0_current_leading_postpulse_eta_flattening",
            ),
        ):
            np.testing.assert_allclose(
                got[child_key], ref[parent_key], rtol=0.0, atol=0.0
            )


def test_unedited_baseline_is_executable_over_full_public_hold() -> None:
    c = _candidate()
    eta = np.array([-0.62, -0.15, 0.0, 0.41, 0.88])
    for frac in (0.0, 0.43, 0.91, 1.0):
        log_X = np.full(
            eta.shape, c.log_X_flatten_end + frac * c.hold_length
        )
        p = c.unedited_hold_profile_logX(log_X, eta)
        np.testing.assert_array_equal(p["U_unedited"], np.zeros(eta.shape))
        assert np.all(p["E_unedited"] > 0.0)
        assert np.all(p["F_unedited"] >= 0.0)
        np.testing.assert_allclose(
            p["ell_unedited"], -c.lambda_value, rtol=0.0, atol=0.0
        )

    truth = c.truth_boundary
    assert truth["source_unedited_eta_independent_hold_baseline_materialized"] is True
    assert truth["current_cartesian_pre_relative_swirl_hold_prefix_composed"] is True
    assert truth["source_terminal_hold_after_eta_flattening_materialized"] is False
    assert truth["source_relative_swirl_bumps_materialized"] is False
    assert truth["source_exterior_heat_replacement_materialized"] is False
    assert truth["pde_validated"] is False


def test_unedited_baseline_preserves_eta_flattening_to_float_roundoff() -> None:
    c = _candidate()
    eta = np.linspace(-1.0, 1.0, 11)
    for frac in (0.0, 0.31, 0.73, 1.0):
        p = c.unedited_hold_profile_logX(
            np.full(eta.shape, c.log_X_flatten_end + frac * c.hold_length), eta
        )
        E = np.asarray(p["E_unedited"], dtype=float)
        relative_spread = (np.max(E) - np.min(E)) / np.max(E)
        assert relative_spread < 3e-13
        np.testing.assert_array_equal(
            p["E_eta_unedited"], np.zeros_like(eta)
        )


def test_current_primitives_obey_U_zero_identity_inside_candidate_prefix() -> None:
    c = _candidate()
    eta = np.array([-0.37, 0.19, 0.71])
    h = 2.0e-6
    log_X = c.log_X_flatten_end + 0.47 * c.prefix_length
    pm = c.similarity_profile_values_logX(np.full(eta.shape, log_X - h), eta)
    p0 = c.similarity_profile_values_logX(np.full(eta.shape, log_X), eta)
    pp = c.similarity_profile_values_logX(np.full(eta.shape, log_X + h), eta)
    for key in (
        "M_over_X_current_leading_pre_relative_swirl",
        "M_eta_over_X_current_leading_pre_relative_swirl",
    ):
        centered = (pp[key] - pm[key]) / (2.0 * h)
        np.testing.assert_allclose(centered, -p0[key], rtol=2e-7, atol=2e-300)


def test_analytic_logX_derivatives_match_centered_difference_inside_prefix() -> None:
    c = _candidate()
    eta = np.array([-0.53, -0.17, 0.26, 0.67])
    log_X = c.log_X_flatten_end + 0.39 * c.prefix_length
    h = 1.0e-6
    pm = c.similarity_profile_values_logX(np.full(eta.shape, log_X - h), eta)
    pp = c.similarity_profile_values_logX(np.full(eta.shape, log_X + h), eta)
    d = c.similarity_log_radial_derivatives(np.full(eta.shape, log_X), eta)
    for value_key, deriv_key in (
        (
            "E_current_leading_pre_relative_swirl",
            "E_DlogX_current_leading_pre_relative_swirl",
        ),
        (
            "F_current_leading_pre_relative_swirl",
            "F_DlogX_current_leading_pre_relative_swirl",
        ),
        (
            "U_current_leading_pre_relative_swirl",
            "U_DlogX_current_leading_pre_relative_swirl",
        ),
    ):
        centered = (pp[value_key] - pm[value_key]) / (2.0 * h)
        np.testing.assert_allclose(centered, d[deriv_key], rtol=4e-7, atol=3e-300)

    centered_log_F = (
        pp["log_F_current_leading_pre_relative_swirl"]
        - pm["log_F_current_leading_pre_relative_swirl"]
    ) / (2.0 * h)
    np.testing.assert_allclose(
        centered_log_F,
        d["log_F_DlogX_current_leading_pre_relative_swirl"],
        rtol=3e-7,
        atol=5e-8,
    )
    np.testing.assert_allclose(
        d["ell_current_leading_pre_relative_swirl"],
        -c.lambda_value,
        rtol=0.0,
        atol=0.0,
    )


def test_cartesian_velocity_reaches_safe_pre_bump_endpoint_and_axis_is_regular() -> None:
    c = _candidate()
    radius = math.exp(c.log_radius_q1_prefix_end)
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
        c.log_X_prefix_end, rel=0.0, abs=4e-12
    )


def test_full_unedited_helper_does_not_promote_cartesian_candidate_into_bumps() -> None:
    c = _candidate()
    helper = c.unedited_hold_profile_logX(c.log_X_hold_end, 0.0)
    assert float(helper["U_unedited"]) == 0.0
    assert float(helper["E_unedited"]) > 0.0

    with pytest.raises(ValueError):
        c.similarity_profile_values_logX(c.log_X_prefix_end + 1.0e-4, 0.0)
    radius_beyond = math.exp(
        0.5 * (math.log(2.0) + c.log_X_prefix_end + 1.0e-3)
    )
    with pytest.raises(ValueError):
        c.velocity(radius_beyond, 0.0, 0.0, 0.0)


def test_configuration_semantic_identity_and_truth_guards(tmp_path) -> None:
    c = _candidate()
    path = tmp_path / "postpulse_hold_prefix.json"
    payload = c.save_configuration(path)
    loaded = KokunoPA16CurrentCartesianPostPulseEtaIndependentHoldPrefix.load_configuration(
        path
    )
    assert loaded.configuration() == payload
    assert loaded.semantic_sha256 == c.semantic_sha256

    mutated = copy.deepcopy(payload)
    mutated["bound_scope"]["hold_length"] *= 0.99
    with pytest.raises(ValueError):
        KokunoPA16CurrentCartesianPostPulseEtaIndependentHoldPrefix.from_configuration(
            mutated
        )

    mutated = copy.deepcopy(payload)
    mutated["bound_scope"]["pre_bump_guard_from_hold_end"] = 3.0
    with pytest.raises(ValueError):
        KokunoPA16CurrentCartesianPostPulseEtaIndependentHoldPrefix.from_configuration(
            mutated
        )

    mutated = copy.deepcopy(payload)
    mutated["bound_scope"]["pre_bump_guard_role"] = "source_exact"
    with pytest.raises(ValueError):
        KokunoPA16CurrentCartesianPostPulseEtaIndependentHoldPrefix.from_configuration(
            mutated
        )

    mutated = copy.deepcopy(payload)
    mutated["bound_scope"]["T_f_role"] = "source_exact"
    with pytest.raises(ValueError):
        KokunoPA16CurrentCartesianPostPulseEtaIndependentHoldPrefix.from_configuration(
            mutated
        )

    mutated = copy.deepcopy(payload)
    mutated["truth_boundary"]["source_relative_swirl_bumps_materialized"] = True
    with pytest.raises(ValueError):
        KokunoPA16CurrentCartesianPostPulseEtaIndependentHoldPrefix.from_configuration(
            mutated
        )

    mutated = copy.deepcopy(payload)
    mutated["truth_boundary"]["pde_validated"] = True
    with pytest.raises(ValueError):
        KokunoPA16CurrentCartesianPostPulseEtaIndependentHoldPrefix.from_configuration(
            mutated
        )


def test_deterministic_report_separates_unedited_baseline_from_candidate_scope() -> None:
    c = _candidate()
    report = c.hold_prefix_report()
    assert report["hold_length_public_formula"] == pytest.approx(
        30.0 * math.log(1.0 / c.lambda_value), rel=0.0, abs=0.0
    )
    assert report["relative_swirl_width_public"] == 0.3
    assert report["first_relative_swirl_center_s"] == pytest.approx(
        c.hold_length - 3.0, rel=0.0, abs=0.0
    )
    assert report["second_relative_swirl_center_s"] == pytest.approx(
        c.hold_length - 1.0, rel=0.0, abs=0.0
    )
    assert report["candidate_prefix_length"] == c.prefix_length
    assert report["start_E_relative_eta_spread"] < 3e-13
    assert report["prefix_end_E_relative_eta_spread"] < 3e-13
    assert report["unedited_full_end_E_relative_eta_spread"] < 3e-13
    assert report["ell_unedited"] == -c.lambda_value
    assert report["max_abs_U_prefix_endpoint"] == 0.0
    assert report["truth_boundary"][
        "source_unedited_eta_independent_hold_baseline_materialized"
    ] is True
    assert report["truth_boundary"][
        "source_terminal_hold_after_eta_flattening_materialized"
    ] is False
    assert report["truth_boundary"]["source_relative_swirl_bumps_materialized"] is False
    assert report["truth_boundary"]["pde_validated"] is False
