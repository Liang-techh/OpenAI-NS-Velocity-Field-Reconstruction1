import copy
import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_pa16_current_cartesian_main_pulse import (
    KokunoPA16CurrentCartesianMainPulse,
)
from openai_ns_reconstruction.kokuno_source_main_axial_pulse_kernel import MAIN_XI_MAX


def _candidate():
    return KokunoPA16CurrentCartesianMainPulse()


def _probe_xi(candidate):
    return min(0.2, 0.25 * candidate.xi_materializable_max)


def test_current_lineage_hits_float_X_guard_before_source_xi_11():
    c = _candidate()
    assert c.lambda_value == pytest.approx(0.05, rel=0.0, abs=0.0)
    assert 8.0 < c.xi_materializable_max < MAIN_XI_MAX
    assert c.log_X_materializable_end < math.log(np.finfo(float).max)
    assert c.X_materializable_end < np.finfo(float).max / 2.0
    with pytest.raises(ValueError):
        c.X_from_xi(MAIN_XI_MAX)


def test_exact_parent_replay_through_pulse_entry():
    c = _candidate()
    eta = np.array([-0.4, 0.0, 0.35])
    X = np.array([0.9 * c.X_p, c.X_p, 0.75 * c.X_p])
    got = c.similarity_profile_values(X, eta)
    ref = c.parent.similarity_profile_values(X, eta)
    pairs = (
        ("F_current_leading_with_main_pulse", "F_current_leading_to_pulse_entry"),
        ("U_current_leading_with_main_pulse", "U_current_leading_to_pulse_entry"),
        ("E_current_leading_with_main_pulse", "E_current_leading_to_pulse_entry"),
        (
            "M_over_X_current_leading_with_main_pulse",
            "M_over_X_current_leading_to_pulse_entry",
        ),
        (
            "M_eta_over_X_current_leading_with_main_pulse",
            "M_eta_over_X_current_leading_to_pulse_entry",
        ),
        ("v0_current_leading_with_main_pulse", "v0_current_leading_to_pulse_entry"),
    )
    for child_key, parent_key in pairs:
        np.testing.assert_array_equal(got[child_key], ref[parent_key])


def test_active_pulse_reuses_public_kernel_F_E_U_and_is_nontrivial():
    c = _candidate()
    xi = _probe_xi(c)
    X = np.full(3, c.X_from_xi(xi))
    eta = np.array([-0.45, 0.0, 0.4])
    got = c.similarity_profile_values(X, eta)
    E_entry, _, _, _ = c._pulse_entry_state(eta)
    ref = c.kernel.profile_values(
        X,
        eta,
        X_p=c.X_p,
        E_entry=E_entry,
        lambda_value=c.lambda_value,
    )
    np.testing.assert_allclose(
        got["F_current_leading_with_main_pulse"],
        ref["F_source_main_pulse_principal"],
        rtol=2e-14,
        atol=0.0,
    )
    np.testing.assert_allclose(
        got["E_current_leading_with_main_pulse"],
        ref["E_source_main_pulse_principal"],
        rtol=2e-14,
        atol=0.0,
    )
    np.testing.assert_allclose(
        got["U_current_leading_with_main_pulse"],
        ref["U_source_main_pulse_principal"],
        rtol=2e-14,
        atol=0.0,
    )
    assert np.all(got["U_current_leading_with_main_pulse"] > 0.0)
    assert np.all(got["delta_M_over_X_main_pulse"] > 0.0)


def test_current_primitive_obeys_dlogX_identity_inside_pulse():
    c = _candidate()
    xi = _probe_xi(c)
    X0 = c.X_from_xi(xi)
    eta = np.array([-0.35, 0.15, 0.42])
    h = 2.0e-5
    Xm = X0 * math.exp(-h)
    Xp = X0 * math.exp(h)
    pm = c.similarity_profile_values(np.full(eta.shape, Xm), eta)
    p0 = c.similarity_profile_values(np.full(eta.shape, X0), eta)
    pp = c.similarity_profile_values(np.full(eta.shape, Xp), eta)
    m_log_derivative = (
        pp["M_over_X_current_leading_with_main_pulse"]
        - pm["M_over_X_current_leading_with_main_pulse"]
    ) / (2.0 * h)
    expected = (
        p0["U_current_leading_with_main_pulse"]
        - p0["M_over_X_current_leading_with_main_pulse"]
    )
    np.testing.assert_allclose(m_log_derivative, expected, rtol=3e-7, atol=2e-18)

    E_entry, E_entry_eta, _, _ = c._pulse_entry_state(eta)
    source_eta = c.kernel.eta_derivatives(
        np.full(eta.shape, X0),
        eta,
        X_p=c.X_p,
        E_entry=E_entry,
        E_entry_eta=E_entry_eta,
        lambda_value=c.lambda_value,
    )
    m_eta_log_derivative = (
        pp["M_eta_over_X_current_leading_with_main_pulse"]
        - pm["M_eta_over_X_current_leading_with_main_pulse"]
    ) / (2.0 * h)
    expected_eta = (
        source_eta["U_eta_source_main_pulse_principal"]
        - p0["M_eta_over_X_current_leading_with_main_pulse"]
    )
    np.testing.assert_allclose(
        m_eta_log_derivative, expected_eta, rtol=4e-7, atol=2e-18
    )


def test_analytic_profile_radial_derivatives_match_centered_logX_difference():
    c = _candidate()
    xi = _probe_xi(c)
    X0 = c.X_from_xi(xi)
    eta = np.array([-0.3, 0.2])
    h = 1.0e-6
    Xm, Xp = X0 * math.exp(-h), X0 * math.exp(h)
    pm = c.similarity_profile_values(np.full(eta.shape, Xm), eta)
    pp = c.similarity_profile_values(np.full(eta.shape, Xp), eta)
    deriv = c.similarity_radial_derivatives(np.full(eta.shape, X0), eta)
    for value_key, deriv_key in (
        ("F_current_leading_with_main_pulse", "F_current_leading_with_main_pulse_X"),
        ("U_current_leading_with_main_pulse", "U_current_leading_with_main_pulse_X"),
        ("E_current_leading_with_main_pulse", "E_current_leading_with_main_pulse_X"),
    ):
        centered_dlog = (pp[value_key] - pm[value_key]) / (2.0 * h)
        analytic_dlog = X0 * deriv[deriv_key]
        np.testing.assert_allclose(centered_dlog, analytic_dlog, rtol=2e-6, atol=1e-300)


def test_vectorized_cartesian_velocity_enters_active_pulse_and_axis_is_regular():
    c = _candidate()
    xi = _probe_xi(c)
    X = c.X_from_xi(xi)
    # z=0,t=0 gives q=1, eta=0, so x=sqrt(2X) realizes the pulse probe.
    radius = math.sqrt(2.0 * X)
    v = c.velocity(
        np.array([radius, 0.0]),
        np.array([0.0, 0.0]),
        np.array([0.0, 0.0]),
        np.array([0.0, 0.0]),
    )
    assert v.shape == (2, 3)
    assert np.all(np.isfinite(v))
    assert v[0, 2] > 0.0
    np.testing.assert_array_equal(v[1, :2], np.zeros(2))


def test_configuration_roundtrip_semantic_identity_and_truth_guards(tmp_path):
    c = _candidate()
    path = tmp_path / "candidate.json"
    payload = c.save_configuration(path)
    loaded = KokunoPA16CurrentCartesianMainPulse.load_configuration(path)
    assert loaded.configuration() == payload
    assert loaded.semantic_sha256 == c.semantic_sha256

    truth = c.truth_boundary
    assert truth["current_cartesian_main_pulse_prefix_materialized"] is True
    assert truth["current_pulse_M_and_M_eta_history_integrated"] is True
    assert truth["source_exact_amplitude_root_materialized"] is False
    assert truth["source_pulse_end_MJ_corrections_materialized"] is False
    assert truth["full_source_xi_11_current_cartesian_materialized"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False

    mutated = copy.deepcopy(payload)
    mutated["bound_scope"]["amplitude_role"] = "source_exact"
    with pytest.raises(ValueError):
        KokunoPA16CurrentCartesianMainPulse.from_configuration(mutated)

    mutated = copy.deepcopy(payload)
    mutated["bound_scope"]["full_source_xi_11_current_cartesian_materialized"] = True
    with pytest.raises(ValueError):
        KokunoPA16CurrentCartesianMainPulse.from_configuration(mutated)


def test_fail_closed_beyond_finite_X_prefix():
    c = _candidate()
    with pytest.raises(ValueError):
        c.similarity_profile_values(
            np.nextafter(c.X_materializable_end, np.inf), 0.0
        )
