from __future__ import annotations

import inspect
import json
import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_public_terminal_multiplier_target import (
    C_O_AUTONOMOUS,
    PARENT_EXACT_HEAD,
    QP_QUADRATURE_ORDER,
    SOURCE_COMMIT,
    TERMINAL_LENGTH,
)
from openai_ns_reconstruction.kokuno_public_terminal_multiplier_target import (
    KokunoPublicTerminalMultiplierTarget,
)


def test_terminal_target_source_and_truth_boundary():
    c = KokunoPublicTerminalMultiplierTarget()
    assert PARENT_EXACT_HEAD == "18e75e6e2df43206147db34788c68ee6a7fa0f37"
    assert SOURCE_COMMIT == "143f6773feb424ad9ed3a8d116653200f20346b7"
    assert TERMINAL_LENGTH == 3.0
    assert C_O_AUTONOMOUS == 1.0 / 64.0
    assert QP_QUADRATURE_ORDER == 192
    truth = c.truth_boundary
    assert truth["source_l_minus1_to_minus_h_transition_materialized"] is True
    assert truth["public_terminal_multiplier_target_api_materialized"] is True
    assert truth["repository_autonomous_terminal_c_o_materialized"] is True
    assert truth["source_exact_terminal_c_o_recovered"] is False
    assert truth["current_q_s_release2_endpoint_materialized"] is False
    assert truth["current_l_minus_h_matching_bridge_materialized"] is False
    assert truth["current_cartesian_terminal_multiplier_composed"] is False
    assert truth["source_terminal_multiplier_materialized"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False


def test_terminal_multiplier_endpoints_and_monotonicity():
    c = KokunoPublicTerminalMultiplierTarget()
    p = c.terminal_multiplier(np.array([0.0, 1.0, 2.0, 3.0]))
    rho = c.rho_o
    np.testing.assert_allclose(p["psi_o"], [1.0, 1.0, 0.5, 0.0], rtol=0.0, atol=2e-15)
    np.testing.assert_allclose(p["f_o"], [1.0-rho, 1.0-rho, 1.0-0.5*rho, 1.0], rtol=0.0, atol=2e-15)
    assert np.all(np.asarray(p["f_o_prime"]) >= 0.0)
    assert float(np.asarray(c.terminal_multiplier(0.0)["ell"])) == -c.h_value
    assert float(np.asarray(c.terminal_multiplier(3.0)["ell"])) == -c.h_value


def test_public_slope_inequality_has_wide_margin():
    c = KokunoPublicTerminalMultiplierTarget()
    report = c.source_inequality_report()
    assert report["sampled_inequality_satisfied"] is True
    assert report["certificate_inequality_satisfied"] is True
    assert report["sigma_prime_sampled_max"] <= report["sigma_prime_certificate_bound"]
    assert report["actual_sampled_max_fprime_over_f"] < 0.30 * report["public_upper_limit_h_over_4"]
    assert report["certified_max_fprime_over_f"] < 0.30 * report["public_upper_limit_h_over_4"]


def test_f_and_ell_analytic_derivatives_match_centered_difference():
    c = KokunoPublicTerminalMultiplierTarget()
    y = np.array([1.18, 1.65, 2.0, 2.37, 2.81])
    eps = 2e-6
    p = c.terminal_multiplier(y)
    pp = c.terminal_multiplier(y + eps)
    pm = c.terminal_multiplier(y - eps)

    fd_f = (np.asarray(pp["f_o"]) - np.asarray(pm["f_o"])) / (2.0 * eps)
    np.testing.assert_allclose(fd_f, p["f_o_prime"], rtol=3e-7, atol=3e-12)

    fd_ell = (np.asarray(pp["ell"]) - np.asarray(pm["ell"])) / (2.0 * eps)
    np.testing.assert_allclose(fd_ell, p["D_y_ell"], rtol=2e-5, atol=2e-10)


def test_logH_logE_derivatives_match_public_slopes():
    c = KokunoPublicTerminalMultiplierTarget()
    y = np.array([1.25, 1.75, 2.25, 2.75])
    eps = 1e-6
    ell = np.asarray(c.terminal_multiplier(y)["ell"])

    d_logh = (c.terminal_logH_ratio(y + eps) - c.terminal_logH_ratio(y - eps)) / (2*eps)
    d_loge = (c.terminal_logE_ratio(y + eps) - c.terminal_logE_ratio(y - eps)) / (2*eps)
    np.testing.assert_allclose(d_logh, ell, rtol=5e-8, atol=5e-10)
    np.testing.assert_allclose(d_loge, ell - 0.5, rtol=5e-8, atol=5e-10)


def test_qp_direct_and_reduced_identity_agree_and_are_positive():
    c = KokunoPublicTerminalMultiplierTarget()
    qd = c.q_p_direct()
    qr = c.q_p_reduced()
    assert qd > 0.0
    assert qr > 0.0
    assert abs(qd / qr - 1.0) < 2e-14
    assert math.isclose(c.q_p, qr, rel_tol=0.0, abs_tol=0.0)
    # Frozen current h=.005 + autonomous c_o=1/64 deterministic realization.
    assert math.isclose(qr, 5.746912488789217e-4, rel_tol=2e-12, abs_tol=0.0)


def test_report_records_bridge_sequence_without_claiming_composition():
    c = KokunoPublicTerminalMultiplierTarget()
    r = c.target_report()
    assert r["q_p_reduced"] > 0.0
    assert r["q_p_over_h"] > 0.0
    assert r["q_p_relative_replay_error"] < 2e-14
    assert r["f_o_start"] < 1.0
    assert r["f_o_end"] == 1.0
    truth = r["truth_boundary"]
    assert truth["current_l_minus_h_matching_bridge_materialized"] is False
    assert truth["current_cartesian_terminal_multiplier_composed"] is False
    assert truth["heldout_ns_residual_assessed"] is False


def test_terminal_stage_guards():
    c = KokunoPublicTerminalMultiplierTarget()
    with pytest.raises(ValueError):
        c.terminal_multiplier(-1e-4)
    with pytest.raises(ValueError):
        c.terminal_multiplier(3.0001)
    with pytest.raises(ValueError):
        c.terminal_multiplier(np.nan)


def test_configuration_roundtrip_and_semantic_identity(tmp_path):
    c = KokunoPublicTerminalMultiplierTarget()
    path = tmp_path / "terminal_target.json"
    payload = c.save_configuration(path)
    loaded = KokunoPublicTerminalMultiplierTarget.load_configuration(path)
    assert loaded.configuration() == payload
    assert loaded.semantic_sha256 == c.semantic_sha256

    broken = json.loads(path.read_text())
    broken["c_o_autonomous"] *= 2.0
    with pytest.raises(ValueError):
        KokunoPublicTerminalMultiplierTarget.from_configuration(broken)


def test_public_constructor_exposes_no_scientific_tuning_knobs():
    params = set(inspect.signature(KokunoPublicTerminalMultiplierTarget).parameters)
    forbidden = {
        "c_o", "rho_o", "h", "residual", "forcing", "pressure", "optimizer",
        "threshold", "tolerance", "quadrature_order", "precision", "gain",
    }
    assert not (params & forbidden)
    assert params <= {"parent"}
