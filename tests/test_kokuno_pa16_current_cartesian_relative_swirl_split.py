from __future__ import annotations

import copy
import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_pa16_current_cartesian_relative_swirl_split import (
    PARENT_EXACT_HEAD,
    KokunoPA16CurrentCartesianRelativeSwirlSplit,
)


def _candidate() -> KokunoPA16CurrentCartesianRelativeSwirlSplit:
    return KokunoPA16CurrentCartesianRelativeSwirlSplit()


def test_exact_1168_target_and_1154_solver_are_consumed_without_truth_promotion():
    c = _candidate()
    eta = np.array([-0.2, 0.0, 0.2])
    target = c.target.total_target(eta)
    solution = c.compensator.solve(target)
    assert np.max(np.abs(solution.angular_residual)) <= 2e-14
    assert np.max(np.abs(solution.pressure_residual)) <= 2e-14
    assert np.max(np.abs(target)) > 0.0
    truth = c.truth_boundary
    assert truth["relative_swirl_profile_delta_channel_materialized"] is True
    assert truth["relative_swirl_cartesian_delta_channel_materialized"] is True
    assert truth["binary64_total_relative_swirl_sum_is_resolved"] is False
    assert truth["current_cartesian_relative_swirl_composed"] is False
    assert truth["pde_validated"] is False


def test_split_profile_preserves_nonzero_sub_epsilon_edit_instead_of_rounding_it_away():
    c = _candidate()
    eta = np.array([-0.2, 0.0, 0.2])
    s = np.full(eta.shape, c.compensator.y1)
    p = c.split_profile_logX(c.log_X_flatten_end + s, eta)
    edit = np.asarray(p["relative_edit"])
    assert np.all(np.isfinite(edit))
    assert np.max(np.abs(edit)) > 0.0
    assert np.max(np.abs(edit)) < np.finfo(float).eps
    np.testing.assert_array_equal(1.0 + edit, np.ones_like(edit))
    assert np.max(np.abs(p["delta_E"])) > 0.0
    assert np.max(np.abs(p["delta_F"])) > 0.0
    np.testing.assert_array_equal(p["delta_U"], np.zeros_like(edit))
    np.testing.assert_array_equal(p["delta_v0"], np.zeros_like(edit))


def test_split_delta_logX_derivatives_match_centered_difference():
    c = _candidate()
    eta = np.array([-0.17, 0.11])
    log_X = c.log_X_flatten_end + c.compensator.y1 + 0.037
    h = 1.0e-6
    pm = c.split_profile_logX(np.full(eta.shape, log_X - h), eta)
    p0 = c.split_profile_logX(np.full(eta.shape, log_X), eta)
    pp = c.split_profile_logX(np.full(eta.shape, log_X + h), eta)
    for value_key, deriv_key in (
        ("delta_E", "delta_E_DlogX"),
        ("delta_F", "delta_F_DlogX"),
    ):
        centered = (pp[value_key] - pm[value_key]) / (2.0 * h)
        scale = np.maximum(np.abs(centered), np.abs(p0[deriv_key]))
        np.testing.assert_allclose(
            centered / np.maximum(scale, 1e-300),
            p0[deriv_key] / np.maximum(scale, 1e-300),
            rtol=2e-6,
            atol=2e-7,
        )


def test_eta_jet_keeps_the_small_correction_executable():
    c = _candidate()
    eta = 0.13
    log_X = c.log_X_flatten_end + c.compensator.y2 + 0.021
    jet = c.split_profile_eta_jet_logX(log_X, eta)
    assert math.isfinite(float(jet["relative_edit_eta"]))
    assert math.isfinite(float(jet["delta_E_eta"]))
    h = 4.0e-5
    plus = c.split_profile_logX(log_X, eta + h)
    minus = c.split_profile_logX(log_X, eta - h)
    independent = (float(plus["delta_E"]) - float(minus["delta_E"])) / (2.0 * h)
    actual = float(jet["delta_E_eta"])
    scale = max(abs(independent), abs(actual), 1e-300)
    assert abs(actual - independent) / scale < 2e-3


def test_cartesian_split_returns_representable_nonzero_delta_and_regular_axis():
    c = _candidate()
    log_X = c.log_X_flatten_end + c.compensator.y1
    radius = math.exp(0.5 * (math.log(2.0) + log_X))
    base, delta = c.velocity_split(
        np.array([radius, 0.0]),
        np.array([0.0, 0.0]),
        np.array([0.0, 0.0]),
        np.array([0.0, 0.0]),
    )
    assert base.shape == delta.shape == (2, 3)
    assert np.all(np.isfinite(base))
    assert np.all(np.isfinite(delta))
    assert np.max(np.abs(delta[0])) > 0.0
    np.testing.assert_array_equal(delta[1], np.zeros(3))
    np.testing.assert_array_equal(base[1, :2], np.zeros(2))


def test_representation_report_exposes_float64_resolution_blocker():
    c = _candidate()
    report = c.representation_report()
    assert report["max_abs_relative_edit"] > 0.0
    assert report["edit_is_below_float64_relative_epsilon"] is True
    assert report["max_abs_rounded_1_plus_edit_minus_1"] == 0.0
    assert report["max_abs_angular_residual"] <= 2e-14
    assert report["max_abs_pressure_residual"] <= 2e-14
    assert report["truth_boundary"]["current_cartesian_relative_swirl_composed"] is False


def test_configuration_roundtrip_and_provenance_fail_closed(tmp_path):
    c = _candidate()
    assert c.configuration()["parent_exact_head"] == PARENT_EXACT_HEAD
    path = tmp_path / "relative_swirl_split.json"
    payload = c.save_configuration(path)
    loaded = KokunoPA16CurrentCartesianRelativeSwirlSplit.load_configuration(path)
    assert loaded.configuration() == payload
    assert loaded.semantic_sha256 == c.semantic_sha256

    mutated = copy.deepcopy(payload)
    mutated["truth_boundary"]["current_cartesian_relative_swirl_composed"] = True
    with pytest.raises(ValueError, match="configuration/provenance drift"):
        KokunoPA16CurrentCartesianRelativeSwirlSplit.from_configuration(mutated)

    mutated = copy.deepcopy(payload)
    mutated["truth_boundary"]["source_exact_pointwise_bump_shape_recovered"] = True
    with pytest.raises(ValueError, match="configuration/provenance drift"):
        KokunoPA16CurrentCartesianRelativeSwirlSplit.from_configuration(mutated)
