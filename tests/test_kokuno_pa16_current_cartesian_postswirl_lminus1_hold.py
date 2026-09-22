from __future__ import annotations

import inspect
import json
import math
from decimal import Decimal

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_pa16_current_cartesian_postswirl_lminus1_hold import (
    DECIMAL_DIGITS,
    SOURCE_BLOB,
    SOURCE_COMMIT,
    SOURCE_H_UPPER_BOUND,
    KokunoPA16CurrentCartesianPostSwirlLMinus1Hold,
)


def test_source_h_binding_truth_and_scope():
    field = KokunoPA16CurrentCartesianPostSwirlLMinus1Hold()
    truth = field.truth_boundary
    assert SOURCE_COMMIT == "143f6773feb424ad9ed3a8d116653200f20346b7"
    assert SOURCE_BLOB == "205a99807302e21a51c5eaf223390c0dfc42bcd0"
    assert DECIMAL_DIGITS == 96
    assert SOURCE_H_UPPER_BOUND == 0.01
    assert field.h_value == pytest.approx(field.A - 0.5, rel=0.0, abs=0.0)
    assert field.h_value == pytest.approx(0.5 - field.D, rel=0.0, abs=2e-15)
    assert 0.0 < field.h_value < SOURCE_H_UPPER_BOUND
    assert field.hold_length == pytest.approx(
        4.0 * math.log(1.0 / field.h_value), rel=0.0, abs=2e-14
    )
    assert truth["source_l_minus1_hold_materialized"] is True
    assert truth["current_cartesian_l_minus1_hold_composed"] is True
    assert truth["current_h_bound_from_similarity_exponent"] is True
    assert truth["source_exact_h_recovered"] is False
    assert truth["decimal_output_encoding_materialized"] is True
    assert truth["end_to_end_96digit_arithmetic_materialized"] is False
    for name in (
        "source_l_minus1_to_minus_h_transition_materialized",
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


def test_exact_release1_endpoint_replay_and_public_hold_endpoint_powers():
    field = KokunoPA16CurrentCartesianPostSwirlLMinus1Hold()
    eta = np.asarray([-0.7, -0.25, 0.0, 0.31, 0.76])
    child0 = field.profile_logX(field.log_X_hold_start, eta)
    parent0 = field.parent.profile_logX(field.log_X_hold_start, eta)
    for key in ("E", "F", "log_F", "U", "M_over_X", "M_eta_over_X", "v0"):
        np.testing.assert_array_equal(child0[key], parent0[key])
    np.testing.assert_array_equal(child0["ell"], -np.ones_like(eta))

    child1 = field.profile_logX(field.log_X_hold_end, eta)
    h = field.h_value
    np.testing.assert_allclose(child1["E"] / child0["E"], h**6, rtol=8e-14, atol=0.0)
    np.testing.assert_allclose(child1["F"] / child0["F"], h**8, rtol=8e-14, atol=0.0)
    np.testing.assert_allclose(
        child1["M_over_X"] / child0["M_over_X"], h**4, rtol=8e-14, atol=0.0
    )
    np.testing.assert_allclose(
        child1["M_eta_over_X"] / child0["M_eta_over_X"], h**4, rtol=8e-14, atol=0.0
    )
    np.testing.assert_array_equal(child1["ell"], -np.ones_like(eta))
    assert np.all(child1["U"] == 0.0)


def test_analytic_logX_derivatives_match_centered_differences():
    field = KokunoPA16CurrentCartesianPostSwirlLMinus1Hold()
    eps = 2e-6
    eta = np.asarray([-0.42, 0.0, 0.53])
    for frac in (0.15, 0.5, 0.85):
        x0 = field.log_X_hold_start + frac * field.hold_length
        pm = field.profile_logX(x0 - eps, eta)
        p0 = field.profile_logX(x0, eta)
        pp = field.profile_logX(x0 + eps, eta)
        for key, dkey in (
            ("E", "DlogX_E"),
            ("F", "DlogX_F"),
            ("M_over_X", "DlogX_M_over_X"),
            ("M_eta_over_X", "DlogX_M_eta_over_X"),
            ("v0", "DlogX_v0"),
        ):
            fd = (np.asarray(pp[key]) - np.asarray(pm[key])) / (2.0 * eps)
            np.testing.assert_allclose(np.asarray(p0[dkey]), fd, rtol=4e-7, atol=1e-13)
        assert np.all(p0["DlogX_U"] == 0.0)


def test_public_U_zero_primitive_identity_and_nontrivial_swirl():
    field = KokunoPA16CurrentCartesianPostSwirlLMinus1Hold()
    p = field.profile_logX(
        field.log_X_hold_start + 0.63 * field.hold_length,
        np.asarray([-0.6, 0.0, 0.6]),
    )
    np.testing.assert_array_equal(p["DlogX_M_over_X"], -p["M_over_X"])
    np.testing.assert_array_equal(p["DlogX_M_eta_over_X"], -p["M_eta_over_X"])
    assert np.all(p["U"] == 0.0)
    assert np.all(np.asarray(p["F"]) > 0.0)


def test_decimal_cartesian_seam_endpoint_axis_and_stage_guard():
    field = KokunoPA16CurrentCartesianPostSwirlLMinus1Hold()
    r0 = math.exp(0.5 * (math.log(2.0) + field.log_X_hold_start))
    parent = field.parent.velocity(r0, 0.0, 0.0, 0.0)
    child = field.velocity(r0, 0.0, 0.0, 0.0)
    assert child.dtype == object
    assert all(isinstance(v, Decimal) for v in child.reshape(-1))
    assert list(child.reshape(-1)) == list(parent.reshape(-1))

    r1 = math.exp(0.5 * (math.log(2.0) + field.log_X_hold_end))
    endpoint = field.velocity(r1 / math.sqrt(2.0), r1 / math.sqrt(2.0), 0.0, 0.0)
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
        field.profile_logX(field.log_X_hold_start - 1e-3, 0.0)
    with pytest.raises(ValueError):
        field.profile_logX(field.log_X_hold_end + 1e-3, 0.0)


def test_save_load_semantic_identity_report_and_mutation_fail_closed(tmp_path):
    field = KokunoPA16CurrentCartesianPostSwirlLMinus1Hold()
    path = tmp_path / "candidate.json"
    payload = field.save_configuration(path)
    loaded = KokunoPA16CurrentCartesianPostSwirlLMinus1Hold.load_configuration(path)
    assert loaded.configuration() == payload
    assert loaded.semantic_sha256 == field.semantic_sha256

    report = field.hold_report()
    h = field.h_value
    assert report["parent_child_decimal_seam_exact"] is True
    assert report["h_current"] == pytest.approx(h, rel=0.0, abs=0.0)
    assert report["hold_length"] == pytest.approx(4.0 * math.log(1.0 / h), rel=0.0, abs=2e-14)
    assert report["E_ratio_end_min"] == pytest.approx(h**6, rel=8e-14)
    assert report["E_ratio_end_max"] == pytest.approx(h**6, rel=8e-14)
    assert report["F_ratio_end_min"] == pytest.approx(h**8, rel=8e-14)
    assert report["F_ratio_end_max"] == pytest.approx(h**8, rel=8e-14)
    assert report["end_to_end_96digit_arithmetic_materialized"] is False

    bad = json.loads(json.dumps(payload))
    bad["h_current_from_A_minus_half"] *= 1.01
    with pytest.raises(ValueError):
        KokunoPA16CurrentCartesianPostSwirlLMinus1Hold.from_configuration(bad)

    bad2 = json.loads(json.dumps(payload))
    bad2["truth_boundary"]["source_l_minus1_to_minus_h_transition_materialized"] = True
    with pytest.raises(ValueError):
        KokunoPA16CurrentCartesianPostSwirlLMinus1Hold.from_configuration(bad2)


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
    field = KokunoPA16CurrentCartesianPostSwirlLMinus1Hold()
    for method_name in ("profile_logX", "velocity", "hold_report"):
        sig = inspect.signature(getattr(field, method_name))
        assert not (bad_words & set(sig.parameters))
