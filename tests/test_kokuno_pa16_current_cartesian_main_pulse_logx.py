import copy
import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_pa16_current_cartesian_main_pulse_logx import (
    KokunoPA16CurrentCartesianMainPulseLogX,
)
from openai_ns_reconstruction.kokuno_source_main_axial_pulse_kernel import MAIN_XI_MAX


def _candidate():
    return KokunoPA16CurrentCartesianMainPulseLogX()


def test_logX_path_crosses_float64_X_ceiling_but_physical_radius_stays_finite():
    c = _candidate()
    assert c.lambda_value == pytest.approx(0.05, rel=0.0, abs=0.0)
    assert c.parent.xi_materializable_max < MAIN_XI_MAX
    assert c.log_X_source_end > math.log(np.finfo(float).max)
    assert c.log_radius_q1_source_end < math.log(np.finfo(float).max)
    assert math.isfinite(math.exp(c.log_radius_q1_source_end))
    assert c.truth_boundary["overflow_safe_logX_similarity_materialized"] is True
    assert c.truth_boundary["full_source_xi_11_current_cartesian_materialized"] is True


def test_parent_replay_on_representable_active_pulse_prefix():
    c = _candidate()
    xi = min(0.25, 0.25 * c.parent.xi_materializable_max)
    log_X = c.log_X_from_xi(xi)
    X = math.exp(log_X)
    eta = np.array([-0.35, 0.0, 0.4])
    got = c.similarity_profile_values_logX(np.full(eta.shape, log_X), eta)
    ref = c.parent.similarity_profile_values(np.full(eta.shape, X), eta)
    for child_key, parent_key in (
        ("F_current_leading_with_main_pulse_logX", "F_current_leading_with_main_pulse"),
        ("U_current_leading_with_main_pulse_logX", "U_current_leading_with_main_pulse"),
        ("E_current_leading_with_main_pulse_logX", "E_current_leading_with_main_pulse"),
        (
            "M_over_X_current_leading_with_main_pulse_logX",
            "M_over_X_current_leading_with_main_pulse",
        ),
        (
            "M_eta_over_X_current_leading_with_main_pulse_logX",
            "M_eta_over_X_current_leading_with_main_pulse",
        ),
        ("v0_current_leading_with_main_pulse_logX", "v0_current_leading_with_main_pulse"),
    ):
        np.testing.assert_allclose(got[child_key], ref[parent_key], rtol=8e-14, atol=0.0)


def test_full_public_endpoint_xi_11_is_executable_without_forming_X():
    c = _candidate()
    log_X = c.log_X_from_xi(MAIN_XI_MAX)
    eta = np.array([-0.4, 0.0, 0.35])
    got = c.similarity_profile_values_logX(np.full(eta.shape, log_X), eta)
    assert np.all(np.isfinite(got["F_current_leading_with_main_pulse_logX"]))
    assert np.all(np.isfinite(got["U_current_leading_with_main_pulse_logX"]))
    assert np.all(np.isfinite(got["E_current_leading_with_main_pulse_logX"]))
    assert np.all(np.isfinite(got["M_over_X_current_leading_with_main_pulse_logX"]))
    assert np.all(got["F_current_leading_with_main_pulse_logX"] > 0.0)
    np.testing.assert_allclose(
        got["xi_current_main_pulse_logX"], MAIN_XI_MAX, rtol=0.0, atol=3e-13
    )
    # Public R0 has switched off by xi=11; the primitive history remains nonzero.
    np.testing.assert_allclose(
        got["U_current_leading_with_main_pulse_logX"], 0.0, rtol=0.0, atol=0.0
    )
    assert np.any(np.abs(got["M_over_X_current_leading_with_main_pulse_logX"]) > 0.0)


def test_log_radial_derivatives_match_independent_centered_logX_difference_beyond_X_ceiling():
    c = _candidate()
    xi = 0.5 * (c.parent.xi_materializable_max + MAIN_XI_MAX)
    log_X = c.log_X_from_xi(xi)
    assert log_X > math.log(np.finfo(float).max)
    eta = np.array([-0.25, 0.2])
    h = 1.0e-6
    pm = c.similarity_profile_values_logX(np.full(eta.shape, log_X - h), eta)
    pp = c.similarity_profile_values_logX(np.full(eta.shape, log_X + h), eta)
    d = c.similarity_log_radial_derivatives(np.full(eta.shape, log_X), eta)
    for value_key, deriv_key in (
        ("F_current_leading_with_main_pulse_logX", "F_DlogX_current_leading_with_main_pulse"),
        ("U_current_leading_with_main_pulse_logX", "U_DlogX_current_leading_with_main_pulse"),
        ("E_current_leading_with_main_pulse_logX", "E_DlogX_current_leading_with_main_pulse"),
    ):
        centered = (pp[value_key] - pm[value_key]) / (2.0 * h)
        np.testing.assert_allclose(centered, d[deriv_key], rtol=3e-6, atol=1e-300)


def test_current_primitive_logX_identity_holds_beyond_float64_X():
    c = _candidate()
    xi = 0.5 * (c.parent.xi_materializable_max + 10.0)
    log_X = c.log_X_from_xi(xi)
    assert log_X > math.log(np.finfo(float).max)
    eta = np.array([-0.3, 0.15, 0.4])
    h = 2.0e-5
    pm = c.similarity_profile_values_logX(np.full(eta.shape, log_X - h), eta)
    p0 = c.similarity_profile_values_logX(np.full(eta.shape, log_X), eta)
    pp = c.similarity_profile_values_logX(np.full(eta.shape, log_X + h), eta)
    centered = (
        pp["M_over_X_current_leading_with_main_pulse_logX"]
        - pm["M_over_X_current_leading_with_main_pulse_logX"]
    ) / (2.0 * h)
    expected = (
        p0["U_current_leading_with_main_pulse_logX"]
        - p0["M_over_X_current_leading_with_main_pulse_logX"]
    )
    np.testing.assert_allclose(centered, expected, rtol=4e-7, atol=3e-18)

    # A_principal is eta independent, so U_eta is the entry-shape derivative times
    # the same pulse factor.  This gives an implementation-distinct M_eta identity.
    E_entry, E_entry_eta, _, _ = c._pulse_entry_state(eta)
    ratio = E_entry_eta / E_entry
    U_eta = p0["U_current_leading_with_main_pulse_logX"] * ratio
    centered_eta = (
        pp["M_eta_over_X_current_leading_with_main_pulse_logX"]
        - pm["M_eta_over_X_current_leading_with_main_pulse_logX"]
    ) / (2.0 * h)
    expected_eta = U_eta - p0["M_eta_over_X_current_leading_with_main_pulse_logX"]
    np.testing.assert_allclose(centered_eta, expected_eta, rtol=5e-7, atol=3e-18)


def test_cartesian_velocity_reaches_xi_11_with_finite_radius_and_axis_stays_regular():
    c = _candidate()
    # z=t=0 gives q=1 and eta=0; unlike X, r=sqrt(2X) remains finite at xi=11.
    radius = math.exp(c.log_radius_q1_source_end)
    v = c.velocity(
        np.array([radius, 0.0]),
        np.array([0.0, 0.0]),
        np.array([0.0, 0.0]),
        np.array([0.0, 0.0]),
    )
    assert v.shape == (2, 3)
    assert np.all(np.isfinite(v))
    np.testing.assert_array_equal(v[1, :2], np.zeros(2))
    coords = c.similarity_coordinates_logX(radius, 0.0, 0.0, 0.0)
    xi = c.lambda_value * (float(coords["log_X"]) - c.log_X_p)
    assert xi == pytest.approx(MAIN_XI_MAX, rel=0.0, abs=3e-12)


def test_configuration_semantic_identity_and_provenance_guards(tmp_path):
    c = _candidate()
    path = tmp_path / "logx_candidate.json"
    payload = c.save_configuration(path)
    loaded = KokunoPA16CurrentCartesianMainPulseLogX.load_configuration(path)
    assert loaded.configuration() == payload
    assert loaded.semantic_sha256 == c.semantic_sha256
    truth = c.truth_boundary
    assert truth["full_source_xi_11_current_cartesian_materialized"] is True
    assert truth["source_exact_amplitude_root_materialized"] is False
    assert truth["source_pulse_end_MJ_corrections_materialized"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False

    mutated = copy.deepcopy(payload)
    mutated["bound_scope"]["amplitude_role"] = "source_exact"
    with pytest.raises(ValueError):
        KokunoPA16CurrentCartesianMainPulseLogX.from_configuration(mutated)

    mutated = copy.deepcopy(payload)
    mutated["bound_scope"]["source_pulse_end_MJ_corrections_materialized"] = True
    with pytest.raises(ValueError):
        KokunoPA16CurrentCartesianMainPulseLogX.from_configuration(mutated)


def test_fail_closed_beyond_public_xi_11_endpoint():
    c = _candidate()
    with pytest.raises(ValueError):
        c.similarity_profile_values_logX(c.log_X_from_xi(MAIN_XI_MAX) + 1.0e-4, 0.0)
    with pytest.raises(ValueError):
        c.log_X_from_xi(MAIN_XI_MAX + 1.0e-6)
