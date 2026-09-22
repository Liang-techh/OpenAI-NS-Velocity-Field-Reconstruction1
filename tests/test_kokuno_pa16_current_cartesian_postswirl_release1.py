from __future__ import annotations

import inspect
import json
import math
from decimal import Decimal

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_pa16_current_cartesian_postswirl_release1 import (
    DECIMAL_DIGITS,
    RELEASE1_LENGTH,
    SOURCE_BLOB,
    SOURCE_COMMIT,
    KokunoPA16CurrentCartesianPostSwirlRelease1,
    _sigma,
    _sigma_integral,
    _sigma_prime,
)


def test_source_truth_boundary_and_stage_geometry():
    field = KokunoPA16CurrentCartesianPostSwirlRelease1()
    truth = field.truth_boundary
    assert SOURCE_COMMIT == "143f6773feb424ad9ed3a8d116653200f20346b7"
    assert SOURCE_BLOB == "205a99807302e21a51c5eaf223390c0dfc42bcd0"
    assert RELEASE1_LENGTH == 1.0
    assert DECIMAL_DIGITS == 96
    assert truth["post_relative_swirl_release1_materialized"] is True
    assert truth["post_relative_swirl_release1_public_velocity_materialized"] is True
    assert truth["source_release1_exact_schedule_materialized"] is True
    for name in (
        "source_l_minus1_hold_materialized",
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


def test_public_sigma_endpoints_symmetry_integral_and_derivative():
    x = np.asarray([0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0])
    s = _sigma(x)
    assert s[0] == 0.0
    assert s[-1] == 1.0
    assert s[3] == pytest.approx(0.5, abs=2e-15)
    np.testing.assert_allclose(s + _sigma(1.0 - x), 1.0, rtol=0.0, atol=3e-15)
    assert float(_sigma_integral(np.asarray(0.0))) == 0.0
    assert float(_sigma_integral(np.asarray(1.0))) == pytest.approx(0.5, abs=2e-14)
    h = 2e-6
    probes = np.asarray([0.2, 0.35, 0.5, 0.68, 0.82])
    fd = (_sigma(probes + h) - _sigma(probes - h)) / (2.0 * h)
    np.testing.assert_allclose(_sigma_prime(probes), fd, rtol=3e-7, atol=2e-9)


def test_exact_relative_swirl_endpoint_seam_and_release_endpoint_ratios():
    field = KokunoPA16CurrentCartesianPostSwirlRelease1()
    eta = np.asarray([-0.6, -0.2, 0.0, 0.3, 0.7])
    child0 = field.profile_logX(field.log_X_rel, eta)
    parent0 = field.split.split_profile_logX(field.log_X_rel, eta)
    assert np.all(np.asarray(parent0["relative_edit"]) == 0.0)
    np.testing.assert_array_equal(child0["E"], parent0["E_base"])
    np.testing.assert_array_equal(child0["F"], parent0["F_base"])
    np.testing.assert_array_equal(child0["M_over_X"], parent0["M_over_X_base"])
    np.testing.assert_array_equal(child0["M_eta_over_X"], parent0["M_eta_over_X_base"])
    np.testing.assert_array_equal(child0["v0"], parent0["v0_base"])
    assert np.all(child0["U"] == 0.0)

    child1 = field.profile_logX(field.log_X_release1_end, eta)
    lam = field.lambda_value
    expected_E = math.exp(-1.0 - 0.5 * lam)
    expected_F = math.exp(-1.5 - 0.5 * lam)
    np.testing.assert_allclose(
        child1["E"] / child0["E"], expected_E, rtol=3e-14, atol=0.0
    )
    np.testing.assert_allclose(
        child1["F"] / child0["F"], expected_F, rtol=3e-14, atol=0.0
    )
    np.testing.assert_allclose(child0["ell"], -lam, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(child1["ell"], -1.0, rtol=0.0, atol=2e-15)


def test_analytic_logX_derivatives_match_centered_differences():
    field = KokunoPA16CurrentCartesianPostSwirlRelease1()
    h = 2e-6
    eta = np.asarray([-0.35, 0.0, 0.42])
    for s in (0.2, 0.5, 0.8):
        x0 = field.log_X_rel + s
        pm = field.profile_logX(x0 - h, eta)
        p0 = field.profile_logX(x0, eta)
        pp = field.profile_logX(x0 + h, eta)
        for key, dkey in (
            ("E", "DlogX_E"),
            ("F", "DlogX_F"),
            ("M_over_X", "DlogX_M_over_X"),
            ("M_eta_over_X", "DlogX_M_eta_over_X"),
            ("v0", "DlogX_v0"),
        ):
            fd = (np.asarray(pp[key]) - np.asarray(pm[key])) / (2.0 * h)
            np.testing.assert_allclose(
                np.asarray(p0[dkey]), fd, rtol=3e-7, atol=1e-13
            )
        assert np.all(p0["U"] == 0.0)
        assert np.all(p0["DlogX_U"] == 0.0)


def test_primitive_transport_is_the_public_U_zero_identity():
    field = KokunoPA16CurrentCartesianPostSwirlRelease1()
    eta = np.asarray([-0.8, -0.1, 0.4, 0.85])
    p = field.profile_logX(field.log_X_rel + 0.63, eta)
    np.testing.assert_array_equal(p["DlogX_M_over_X"], -p["M_over_X"])
    np.testing.assert_array_equal(
        p["DlogX_M_eta_over_X"], -p["M_eta_over_X"]
    )


def test_decimal_cartesian_velocity_replays_parent_seam_and_is_nontrivial():
    field = KokunoPA16CurrentCartesianPostSwirlRelease1()
    r0 = math.exp(0.5 * (math.log(2.0) + field.log_X_rel))
    parent = field.parent.velocity(r0, 0.0, 0.0, 0.0)
    child = field.velocity(r0, 0.0, 0.0, 0.0)
    assert child.dtype == object
    assert all(isinstance(v, Decimal) for v in child.reshape(-1))
    assert list(child.reshape(-1)) == list(parent.reshape(-1))

    r1 = math.exp(
        0.5 * (math.log(2.0) + field.log_X_rel + 0.5)
    )
    value = field.velocity(r1 / math.sqrt(2.0), r1 / math.sqrt(2.0), 0.0, 0.0)
    assert value.shape == (3,)
    assert all(isinstance(v, Decimal) for v in value)
    assert any(v != 0 for v in value[:2])
    assert value[2] == 0


def test_axis_regular_path_and_beyond_stage_fail_closed():
    field = KokunoPA16CurrentCartesianPostSwirlRelease1()
    axis = field.velocity(0.0, 0.0, 0.0, 0.0)
    assert axis.shape == (3,)
    assert all(isinstance(v, Decimal) for v in axis)
    assert axis[0] == 0
    assert axis[1] == 0

    with pytest.raises(ValueError):
        field.profile_logX(field.log_X_rel - 1e-3, 0.0)
    with pytest.raises(ValueError):
        field.profile_logX(field.log_X_release1_end + 1e-3, 0.0)


def test_save_load_semantic_identity_and_report(tmp_path):
    field = KokunoPA16CurrentCartesianPostSwirlRelease1()
    path = tmp_path / "candidate.json"
    payload = field.save_configuration(path)
    loaded = KokunoPA16CurrentCartesianPostSwirlRelease1.load_configuration(path)
    assert loaded.configuration() == payload
    assert loaded.semantic_sha256 == field.semantic_sha256
    report = field.release_report()
    assert report["parent_child_decimal_seam_exact"] is True
    assert report["sigma_integral_end"] == pytest.approx(0.5, abs=2e-14)
    assert report["ell_start"] == pytest.approx(-field.lambda_value, abs=0.0)
    assert report["ell_end"] == pytest.approx(-1.0, abs=2e-15)
    assert report["E_ratio_end_min"] == pytest.approx(math.exp(-1.025), rel=3e-14)
    assert report["E_ratio_end_max"] == pytest.approx(math.exp(-1.025), rel=3e-14)
    assert report["F_ratio_end_min"] == pytest.approx(math.exp(-1.525), rel=3e-14)
    assert report["F_ratio_end_max"] == pytest.approx(math.exp(-1.525), rel=3e-14)

    bad = json.loads(json.dumps(payload))
    bad["truth_boundary"]["source_l_minus1_hold_materialized"] = True
    with pytest.raises(ValueError):
        KokunoPA16CurrentCartesianPostSwirlRelease1.from_configuration(bad)


def test_public_api_has_no_scientific_tuning_controls():
    bad_words = {
        "residual",
        "forcing",
        "pressure",
        "threshold",
        "tolerance",
        "optimizer",
        "viscosity",
        "gain",
        "step",
    }
    field = KokunoPA16CurrentCartesianPostSwirlRelease1()
    for method_name in ("profile_logX", "velocity", "release_report"):
        sig = inspect.signature(getattr(field, method_name))
        assert not (bad_words & set(sig.parameters))
