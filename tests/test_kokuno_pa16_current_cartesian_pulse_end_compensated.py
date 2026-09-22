from __future__ import annotations

import copy
import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_pa16_current_cartesian_pulse_end_compensated import (
    KokunoPA16CurrentCartesianPulseEndCompensated,
)
from openai_ns_reconstruction.kokuno_public_pulse_end_compensator import (
    MAIN_XI_END,
    PULSE_XI_END,
)


def _candidate() -> KokunoPA16CurrentCartesianPulseEndCompensated:
    return KokunoPA16CurrentCartesianPulseEndCompensated()


def test_exact_parent_replay_through_xi11() -> None:
    c = _candidate()
    eta = np.array([-0.37, 0.0, 0.41])
    for xi in (0.0, 5.0, 10.5, MAIN_XI_END):
        log_X = c.log_X_from_xi(xi)
        got = c.similarity_profile_values_logX(np.full(eta.shape, log_X), eta)
        ref = c.leading.similarity_profile_values_logX(np.full(eta.shape, log_X), eta)
        for child_key, parent_key in (
            ("F_current_leading_with_pulse_end", "F_current_leading_with_main_pulse_logX"),
            ("U_current_leading_with_pulse_end", "U_current_leading_with_main_pulse_logX"),
            ("E_current_leading_with_pulse_end", "E_current_leading_with_main_pulse_logX"),
            (
                "M_over_X_current_leading_with_pulse_end",
                "M_over_X_current_leading_with_main_pulse_logX",
            ),
            (
                "M_eta_over_X_current_leading_with_pulse_end",
                "M_eta_over_X_current_leading_with_main_pulse_logX",
            ),
            ("v0_current_leading_with_pulse_end", "v0_current_leading_with_main_pulse_logX"),
        ):
            np.testing.assert_allclose(got[child_key], ref[parent_key], rtol=0.0, atol=0.0)


def test_end_bumps_are_executable_nontrivial_and_endpoint_closes_axial_pulse() -> None:
    c = _candidate()
    eta = np.array([-0.42, 0.0, 0.33])
    gap = c.similarity_profile_values_logX(
        np.full(eta.shape, c.log_X_from_xi(12.0)), eta
    )
    np.testing.assert_array_equal(
        gap["U_current_leading_with_pulse_end"], np.zeros(eta.shape)
    )

    xi1 = c.lambda_value * c.compensator.y1
    xi2 = c.lambda_value * c.compensator.y2
    p1 = c.similarity_profile_values_logX(
        np.full(eta.shape, c.log_X_from_xi(xi1)), eta
    )
    p2 = c.similarity_profile_values_logX(
        np.full(eta.shape, c.log_X_from_xi(xi2)), eta
    )
    assert np.max(np.abs(p1["U_current_leading_with_pulse_end"])) > 0.0
    assert np.max(np.abs(p2["U_current_leading_with_pulse_end"])) > 0.0

    end = c.similarity_profile_values_logX(
        np.full(eta.shape, c.log_X_from_xi(PULSE_XI_END)), eta
    )
    np.testing.assert_array_equal(
        end["U_current_leading_with_pulse_end"], np.zeros(eta.shape)
    )
    np.testing.assert_allclose(
        end["xi_current_pulse_end"], PULSE_XI_END, rtol=0.0, atol=3e-13
    )
    assert np.all(np.isfinite(end["M_over_X_current_leading_with_pulse_end"]))
    assert np.all(np.isfinite(end["M_eta_over_X_current_leading_with_pulse_end"]))


def test_current_M_primitive_identity_holds_inside_each_end_bump() -> None:
    c = _candidate()
    eta = np.array([-0.31, 0.17, 0.44])
    h = 2.0e-6
    for y in (c.compensator.y1 + 0.031, c.compensator.y2 + 0.029):
        log_X = c.log_X_p + y
        pm = c.similarity_profile_values_logX(np.full(eta.shape, log_X - h), eta)
        p0 = c.similarity_profile_values_logX(np.full(eta.shape, log_X), eta)
        pp = c.similarity_profile_values_logX(np.full(eta.shape, log_X + h), eta)
        centered = (
            pp["M_over_X_current_leading_with_pulse_end"]
            - pm["M_over_X_current_leading_with_pulse_end"]
        ) / (2.0 * h)
        expected = (
            p0["U_current_leading_with_pulse_end"]
            - p0["M_over_X_current_leading_with_pulse_end"]
        )
        np.testing.assert_allclose(centered, expected, rtol=3e-6, atol=3e-18)


def test_current_M_eta_primitive_identity_holds_inside_end_bump() -> None:
    c = _candidate()
    eta = np.array([-0.29, 0.21])
    h = 2.0e-6
    log_X = c.log_X_p + c.compensator.y2 + 0.033
    pm = c.end_eta_jet_logX(np.full(eta.shape, log_X - h), eta)
    p0 = c.end_eta_jet_logX(np.full(eta.shape, log_X), eta)
    pp = c.end_eta_jet_logX(np.full(eta.shape, log_X + h), eta)
    centered = (
        pp["M_eta_over_X_current_leading_with_pulse_end"]
        - pm["M_eta_over_X_current_leading_with_pulse_end"]
    ) / (2.0 * h)
    expected = (
        p0["U_eta_current_leading_with_pulse_end"]
        - p0["M_eta_over_X_current_leading_with_pulse_end"]
    )
    np.testing.assert_allclose(centered, expected, rtol=4e-6, atol=4e-18)


def test_analytic_log_radial_derivatives_match_centered_difference_in_end_bump() -> None:
    c = _candidate()
    eta = np.array([-0.23, 0.19])
    log_X = c.log_X_p + c.compensator.y1 + 0.041
    h = 1.0e-6
    pm = c.similarity_profile_values_logX(np.full(eta.shape, log_X - h), eta)
    pp = c.similarity_profile_values_logX(np.full(eta.shape, log_X + h), eta)
    d = c.similarity_log_radial_derivatives(np.full(eta.shape, log_X), eta)
    for value_key, deriv_key in (
        ("U_current_leading_with_pulse_end", "U_DlogX_current_leading_with_pulse_end"),
        ("E_current_leading_with_pulse_end", "E_DlogX_current_leading_with_pulse_end"),
    ):
        centered = (pp[value_key] - pm[value_key]) / (2.0 * h)
        np.testing.assert_allclose(centered, d[deriv_key], rtol=3e-6, atol=1e-300)
    centered_log_F = (
        pp["log_F_current_leading_with_pulse_end"]
        - pm["log_F_current_leading_with_pulse_end"]
    ) / (2.0 * h)
    np.testing.assert_allclose(
        centered_log_F,
        d["log_F_DlogX_current_leading_with_pulse_end"],
        rtol=2e-7,
        atol=2e-8,
    )


def test_public_two_row_endpoint_closure_is_machine_small_but_not_pde_validation() -> None:
    c = _candidate()
    eta = np.array([-0.83, -0.37, 0.19, 0.73])
    end = c._end_profile_logX(
        np.full(eta.shape, c.log_X_pulse_end, dtype=float), eta
    )
    assert np.max(np.abs(end["public_M_row_scaled"])) < 3e-15
    assert np.max(np.abs(end["public_J_row_scaled"])) < 3e-15
    truth = c.truth_boundary
    assert truth["current_cartesian_end_compensation_composed"] is True
    assert truth["current_pulse_endpoint_xi13_materialized"] is True
    assert truth["source_exact_amplitude_root_materialized"] is False
    assert truth["source_exact_bump_shape_recovered"] is False
    assert truth["heldout_ns_residual_assessed"] is False
    assert truth["pde_validated"] is False


def test_cartesian_velocity_reaches_xi13_and_axis_stays_regular() -> None:
    c = _candidate()
    radius = math.exp(c.log_radius_q1_pulse_end)
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
    assert xi == pytest.approx(PULSE_XI_END, rel=0.0, abs=4e-12)


def test_configuration_semantic_identity_and_provenance_guards(tmp_path) -> None:
    c = _candidate()
    path = tmp_path / "pulse_end_candidate.json"
    payload = c.save_configuration(path)
    loaded = KokunoPA16CurrentCartesianPulseEndCompensated.load_configuration(path)
    assert loaded.configuration() == payload
    assert loaded.semantic_sha256 == c.semantic_sha256

    mutated = copy.deepcopy(payload)
    mutated["bound_scope"]["amplitude_role"] = "source_exact"
    with pytest.raises(ValueError):
        KokunoPA16CurrentCartesianPulseEndCompensated.from_configuration(mutated)

    mutated = copy.deepcopy(payload)
    mutated["bound_scope"]["bump_shape_role"] = "source_exact"
    with pytest.raises(ValueError):
        KokunoPA16CurrentCartesianPulseEndCompensated.from_configuration(mutated)

    mutated = copy.deepcopy(payload)
    mutated["truth_boundary"]["pde_validated"] = True
    with pytest.raises(ValueError):
        KokunoPA16CurrentCartesianPulseEndCompensated.from_configuration(mutated)


def test_fail_closed_beyond_public_xi13_endpoint() -> None:
    c = _candidate()
    with pytest.raises(ValueError):
        c.similarity_profile_values_logX(c.log_X_pulse_end + 1.0e-4, 0.0)
    with pytest.raises(ValueError):
        c.log_X_from_xi(PULSE_XI_END + 1.0e-6)
