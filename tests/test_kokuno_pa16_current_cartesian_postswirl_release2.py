from __future__ import annotations

import inspect
import json
import math
from decimal import Decimal

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_pa16_current_cartesian_postswirl_release2 import (
    DECIMAL_DIGITS,
    RELEASE2_LENGTH,
    SOURCE_BLOB,
    SOURCE_COMMIT,
    SOURCE_H_UPPER_BOUND,
    KokunoPA16CurrentCartesianPostSwirlRelease2,
)


def test_public_schedule_h_binding_truth_and_scope():
    field = KokunoPA16CurrentCartesianPostSwirlRelease2()
    truth = field.truth_boundary
    assert SOURCE_COMMIT == "143f6773feb424ad9ed3a8d116653200f20346b7"
    assert SOURCE_BLOB == "205a99807302e21a51c5eaf223390c0dfc42bcd0"
    assert DECIMAL_DIGITS == 96
    assert RELEASE2_LENGTH == 1.0
    assert SOURCE_H_UPPER_BOUND == 0.01
    assert field.h_value == pytest.approx(field.A - 0.5, rel=0.0, abs=0.0)
    assert field.h_value == pytest.approx(0.5 - field.D, rel=0.0, abs=2e-15)
    assert field.h_value == pytest.approx(field.parent.h_value, rel=0.0, abs=0.0)
    assert 0.0 < field.h_value < SOURCE_H_UPPER_BOUND

    assert truth["source_l_minus1_hold_materialized"] is True
    assert truth["source_l_minus1_to_minus_h_transition_materialized"] is True
    assert truth["current_cartesian_l_minus1_to_minus_h_composed"] is True
    assert truth["current_h_bound_from_similarity_exponent"] is True
    assert truth["source_exact_h_recovered"] is False
    assert truth["overflow_safe_log_F_materialized"] is True
    assert truth["decimal_output_encoding_materialized"] is True
    assert truth["end_to_end_96digit_arithmetic_materialized"] is False
    for name in (
        "source_terminal_multiplier_materialized",
        "source_exterior_heat_replacement_materialized",
        "outer_global_leading_velocity_materialized",
        "unified_global_cartesian_velocity_export_ready",
        "matched_global_pressure_materialized",
        "restricted_forcing_materialized",
        "heldout_ns_residual_assessed",
        "same_protocol_comparable_to_st006",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        assert truth[name] is False


def test_exact_lminus1_hold_endpoint_replay_and_release2_endpoint_scaling():
    field = KokunoPA16CurrentCartesianPostSwirlRelease2()
    eta = np.asarray([-0.7, -0.25, 0.0, 0.31, 0.76])
    child0 = field.profile_logX(field.log_X_release2_start, eta)
    parent0 = field.parent.profile_logX(field.log_X_release2_start, eta)
    for key in ("E", "F", "log_F", "U", "M_over_X", "M_eta_over_X", "v0"):
        np.testing.assert_array_equal(child0[key], parent0[key])
    np.testing.assert_array_equal(child0["ell"], -np.ones_like(eta))

    child1 = field.profile_logX(field.log_X_release2_end, eta)
    h = field.h_value
    expected_e = math.exp(-1.0 - 0.5 * h)
    expected_f = math.exp(-1.5 - 0.5 * h)
    expected_primitive = math.exp(-1.0)
    np.testing.assert_allclose(
        child1["E"] / child0["E"], expected_e, rtol=8e-14, atol=0.0
    )
    f_ratio_log_chart = np.exp(
        np.asarray(child1["log_F"]) - np.asarray(child0["log_F"])
    )
    np.testing.assert_allclose(f_ratio_log_chart, expected_f, rtol=8e-14, atol=0.0)
    np.testing.assert_allclose(
        child1["M_over_X"],
        child0["M_over_X"] * expected_primitive,
        rtol=8e-14,
        atol=1e-300,
    )
    np.testing.assert_allclose(
        child1["M_eta_over_X"],
        child0["M_eta_over_X"] * expected_primitive,
        rtol=8e-14,
        atol=1e-300,
    )
    np.testing.assert_allclose(child1["sigma_integral"], 0.5, rtol=0.0, atol=2e-14)
    np.testing.assert_allclose(child1["ell"], -h, rtol=0.0, atol=2e-15)
    assert np.all(child1["U"] == 0.0)


def test_public_step_slope_and_analytic_logX_derivatives():
    field = KokunoPA16CurrentCartesianPostSwirlRelease2()
    h = field.h_value
    eta = np.asarray([-0.42, 0.0, 0.53])
    start = field.profile_logX(field.log_X_release2_start, eta)
    mid = field.profile_logX(field.log_X_release2_start + 0.5, eta)
    end = field.profile_logX(field.log_X_release2_end, eta)
    np.testing.assert_array_equal(start["ell"], -np.ones_like(eta))
    np.testing.assert_allclose(mid["ell"], -1.0 + 0.5 * (1.0 - h), atol=2e-15)
    np.testing.assert_allclose(end["ell"], -h, atol=2e-15)

    eps = 2e-6
    for frac in (0.15, 0.5, 0.85):
        x0 = field.log_X_release2_start + frac
        pm = field.profile_logX(x0 - eps, eta)
        p0 = field.profile_logX(x0, eta)
        pp = field.profile_logX(x0 + eps, eta)
        for key, dkey in (
            ("E", "DlogX_E"),
            ("F", "DlogX_F"),
            ("log_F", "DlogX_log_F"),
            ("M_over_X", "DlogX_M_over_X"),
            ("M_eta_over_X", "DlogX_M_eta_over_X"),
            ("v0", "DlogX_v0"),
        ):
            fd = (np.asarray(pp[key]) - np.asarray(pm[key])) / (2.0 * eps)
            np.testing.assert_allclose(np.asarray(p0[dkey]), fd, rtol=5e-7, atol=1e-9)
        assert np.all(p0["DlogX_U"] == 0.0)


def test_public_U_zero_primitive_identity_and_underflow_safe_logF():
    field = KokunoPA16CurrentCartesianPostSwirlRelease2()
    p = field.profile_logX(
        field.log_X_release2_start + 0.63,
        np.asarray([-0.6, 0.0, 0.6]),
    )
    np.testing.assert_array_equal(p["DlogX_M_over_X"], -p["M_over_X"])
    np.testing.assert_array_equal(p["DlogX_M_eta_over_X"], -p["M_eta_over_X"])
    np.testing.assert_allclose(p["DlogX_log_F"], p["ell"] - 1.0, rtol=0.0, atol=0.0)
    assert np.all(p["U"] == 0.0)
    assert np.all(np.isfinite(np.asarray(p["log_F"])))
    assert np.all(np.asarray(p["F"]) >= 0.0)


def test_decimal_cartesian_seam_endpoint_axis_and_stage_guard():
    field = KokunoPA16CurrentCartesianPostSwirlRelease2()
    r0 = math.exp(0.5 * (math.log(2.0) + field.log_X_release2_start))
    parent = field.parent.velocity(r0, 0.0, 0.0, 0.0)
    child = field.velocity(r0, 0.0, 0.0, 0.0)
    assert child.dtype == object
    assert all(isinstance(v, Decimal) for v in child.reshape(-1))
    assert list(child.reshape(-1)) == list(parent.reshape(-1))

    r1 = math.exp(0.5 * (math.log(2.0) + field.log_X_release2_end))
    endpoint = field.velocity(
        r1 / math.sqrt(2.0), r1 / math.sqrt(2.0), 0.0, 0.0
    )
    assert endpoint.shape == (3,)
    assert all(isinstance(v, Decimal) for v in endpoint)
    assert any(v != 0 for v in endpoint[:2])
    assert endpoint[2] == 0

    axis = field.velocity(0.0, 0.0, 0.0, 0.0)
    assert axis.shape == (3,)
    assert all(isinstance(v, Decimal) for v in axis)
    assert axis[0] == 0
    assert axis[1] == 0

    with pytest.raises(ValueError):
        field.profile_logX(field.log_X_release2_start - 1e-3, 0.0)
    with pytest.raises(ValueError):
        field.profile_logX(field.log_X_release2_end + 1e-3, 0.0)


def test_save_load_semantic_identity_report_and_mutation_fail_closed(tmp_path):
    field = KokunoPA16CurrentCartesianPostSwirlRelease2()
    path = tmp_path / "candidate.json"
    payload = field.save_configuration(path)
    loaded = KokunoPA16CurrentCartesianPostSwirlRelease2.load_configuration(path)
    assert loaded.configuration() == payload
    assert loaded.semantic_sha256 == field.semantic_sha256

    report = field.release2_report()
    h = field.h_value
    assert report["parent_child_decimal_seam_exact"] is True
    assert report["h_current"] == pytest.approx(h, rel=0.0, abs=0.0)
    assert report["release2_length"] == 1.0
    assert report["ell_start"] == pytest.approx(-1.0, abs=0.0)
    assert report["ell_end"] == pytest.approx(-h, abs=2e-15)
    assert report["sigma_integral_end"] == pytest.approx(0.5, abs=2e-14)
    assert report["E_ratio_end_min"] == pytest.approx(
        math.exp(-1.0 - 0.5 * h), rel=8e-14
    )
    assert report["E_ratio_end_max"] == pytest.approx(
        math.exp(-1.0 - 0.5 * h), rel=8e-14
    )
    assert report["F_ratio_end_min"] == pytest.approx(
        math.exp(-1.5 - 0.5 * h), rel=8e-14
    )
    assert report["F_ratio_end_max"] == pytest.approx(
        math.exp(-1.5 - 0.5 * h), rel=8e-14
    )
    assert report["end_to_end_96digit_arithmetic_materialized"] is False

    bad = json.loads(json.dumps(payload))
    bad["h_current_from_A_minus_half"] *= 1.01
    with pytest.raises(ValueError):
        KokunoPA16CurrentCartesianPostSwirlRelease2.from_configuration(bad)

    bad2 = json.loads(json.dumps(payload))
    bad2["truth_boundary"]["source_terminal_multiplier_materialized"] = True
    with pytest.raises(ValueError):
        KokunoPA16CurrentCartesianPostSwirlRelease2.from_configuration(bad2)


def test_public_api_has_no_h_or_scientific_tuning_controls():
    bad_words = {
        "h",
        "residual",
        "forcing",
        "pressure",
        "threshold",
        "tolerance",
        "optimizer",
        "viscosity",
        "gain",
        "step",
        "precision",
    }
    field = KokunoPA16CurrentCartesianPostSwirlRelease2()
    for method_name in ("profile_logX", "velocity", "release2_report"):
        sig = inspect.signature(getattr(field, method_name))
        assert not (bad_words & set(sig.parameters))
