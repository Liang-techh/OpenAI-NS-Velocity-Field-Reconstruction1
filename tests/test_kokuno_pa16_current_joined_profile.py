import copy
import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_inner_join_exit import (
    LOG_JOIN_EXIT,
    LOG_REPAIR_START,
)
from openai_ns_reconstruction.kokuno_pa16_current_joined_profile import (
    KokunoPA16CurrentJoinedProfile,
)


@pytest.fixture(scope="module")
def profile():
    return KokunoPA16CurrentJoinedProfile(quadrature_points=64)


def test_prefix_replays_exact_current_parent_and_vectorizes(profile):
    X = np.asarray([0.0, 0.25 * profile.X_i, 0.75 * profile.X_i, profile.X_i])
    eta = np.asarray([-0.75, -0.25, 0.25, 0.75])
    got = profile.profile_values(X, eta)
    parent = profile.moments.profile_values(X, eta)

    np.testing.assert_allclose(
        got["F_current_joined"], parent["F_actual_prefix"], rtol=0, atol=0
    )
    np.testing.assert_allclose(
        got["U_current_joined"], parent["U_actual_prefix"], rtol=0, atol=0
    )
    np.testing.assert_allclose(
        got["E_current_joined"], parent["E_actual_prefix"], rtol=0, atol=0
    )
    np.testing.assert_allclose(
        got["H_current_joined"],
        np.sqrt(2.0 * X) * parent["E_actual_prefix"],
        rtol=0,
        atol=0,
    )
    assert np.all(got["F_current_joined"] > 0.0)
    assert np.all(np.isfinite(got["U_current_joined"]))


def test_geometry_is_fail_closed_or_Xi_value_and_jet_handoff_are_continuous(profile):
    eta = 0.25
    if not profile.route_ready:
        with pytest.raises(ValueError, match="geometry"):
            profile.profile_values(np.nextafter(profile.X_i, math.inf), eta)
        report = profile.Xi_handoff_report(eta)
        assert report["route_ready"] is False
        return

    report = profile.Xi_handoff_report(eta)
    assert report["route_ready"] is True
    np.testing.assert_allclose(
        report["right_F_U_E"], report["left_F_U_E"], rtol=3e-12, atol=2e-13
    )
    np.testing.assert_allclose(
        report["right_FX_UX_EX"], report["left_FX_UX_EX"], rtol=3e-10, atol=2e-13
    )

    # The public wrapper must also be continuous immediately to the right of Xi.
    Xp = np.nextafter(profile.X_i, math.inf)
    right = profile.profile_values(Xp, eta)
    left = profile.profile_values(profile.X_i, eta)
    assert float(right["F_current_joined"]) > 0.0
    assert abs(float(right["U_current_joined"]) - float(left["U_current_joined"])) < 1e-12


def test_join_exit_is_exact_ideal_pair_with_public_radial_slopes(profile):
    if not profile.route_ready:
        pytest.skip("current numerical T_sh geometry is fail-closed")

    for eta in (-0.6, 0.0, 0.6):
        report = profile.ideal_exit_report(eta)
        assert report["route_ready"] is True
        assert report["U_exit"] == pytest.approx(report["U_ideal"], rel=0, abs=2e-13)
        assert report["E_exit"] == pytest.approx(report["E_ideal"], rel=3e-13, abs=0)
        assert report["F_exit"] == pytest.approx(report["F_ideal"], rel=3e-13, abs=0)
        assert report["Xh_U_X"] == pytest.approx(0.0, abs=2e-13)
        assert report["Xh_E_X_over_E"] == pytest.approx(0.1, rel=0, abs=5e-13)
        assert report["Xh_F_X_over_F"] == pytest.approx(-0.4, rel=0, abs=8e-13)


def test_existing_PA16_solve_is_consumed_inside_repair_annulus(profile):
    if not profile.route_ready:
        pytest.skip("current numerical T_sh geometry is fail-closed")

    eta = 0.25
    log_x = np.asarray([-5.8, -5.5, -5.2])
    X = np.exp(profile.log_X_R + log_x)
    got = profile.profile_values(X, eta)
    deriv = profile.radial_derivatives(X, eta)

    assert np.all(np.isfinite(got["U_current_joined"]))
    assert np.all(np.isfinite(got["E_current_joined"]))
    assert np.all(got["F_current_joined"] > 0.0)
    assert np.all(np.isfinite(deriv["F_current_joined_X"]))
    assert np.all(np.isfinite(deriv["U_current_joined_X"]))

    # The repair interval is not silently replaced by the ideal field.  Use all
    # three interior points so an accidental bump node cannot mask the repair.
    f = 1.0 / (1.0 + eta * eta)
    E_ideal = np.exp(profile.moments.outer_schedule.log_P_star) * f * np.exp(0.1 * log_x)
    delta_U = np.max(np.abs(got["U_current_joined"] - 4.0 * eta))
    delta_E_rel = np.max(np.abs(got["E_current_joined"] / E_ideal - 1.0))
    assert max(float(delta_U), float(delta_E_rel)) > 1e-14


def test_physical_X_radial_derivative_replays_independent_centered_difference(profile):
    if not profile.route_ready:
        pytest.skip("current numerical T_sh geometry is fail-closed")

    eta = -0.35
    # This lies in the fixed public axial-restoration interval and outside the
    # PA.16 repair annulus, so the finite-difference replay is implementation-distinct
    # without invoking a second nonlinear repair solve.
    X = math.exp(profile.log_X_R - 7.5)
    h = 2.0e-6 * X
    minus = profile.profile_values(X - h, eta)
    plus = profile.profile_values(X + h, eta)
    deriv = profile.radial_derivatives(X, eta)

    U_fd = (float(plus["U_current_joined"]) - float(minus["U_current_joined"])) / (2.0 * h)
    E_fd = (float(plus["E_current_joined"]) - float(minus["E_current_joined"])) / (2.0 * h)
    F_fd = (float(plus["F_current_joined"]) - float(minus["F_current_joined"])) / (2.0 * h)

    U_scale = max(1.0, abs(float(deriv["U_current_joined_X"]) * X))
    assert abs((U_fd - float(deriv["U_current_joined_X"])) * X) <= 2e-6 * U_scale
    assert E_fd == pytest.approx(float(deriv["E_current_joined_X"]), rel=2e-6, abs=0)
    assert F_fd == pytest.approx(float(deriv["F_current_joined_X"]), rel=2e-6, abs=0)


def test_configuration_roundtrip_and_truth_boundary_remain_scoped(profile, tmp_path):
    path = tmp_path / "joined-profile.json"
    profile.save_configuration(path)
    rebound = KokunoPA16CurrentJoinedProfile.load_configuration(path)
    assert rebound.configuration() == profile.configuration()
    assert rebound.semantic_sha256 == profile.semantic_sha256

    truth = profile.truth_boundary
    assert truth["current_candidate_side_join_Xi_to_Xh_materialized"] is True
    assert truth["current_candidate_side_PA16_coefficients_consumed"] is True
    assert truth["physical_X_F_U_E_profile_executable_through_Xh"] is True
    assert truth["source_B0_analytic_bound_proved"] is False
    assert truth["source_T_sh_lower_bound_verified"] is False
    assert truth["source_global_inner_to_outer_join_admitted"] is False
    assert truth["cartesian_velocity_materialized"] is False
    assert truth["matched_global_pressure_materialized"] is False
    assert truth["restricted_forcing_materialized"] is False
    assert truth["heldout_ns_residual_assessed"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False

    payload = copy.deepcopy(profile.configuration())
    payload["schema"] = "wrong-schema"
    with pytest.raises(ValueError, match="schema"):
        KokunoPA16CurrentJoinedProfile.from_configuration(payload)


def test_domain_and_parameter_guards(profile):
    with pytest.raises(ValueError):
        profile.profile_values(-1.0, 0.0)
    with pytest.raises(ValueError):
        profile.profile_values(np.nextafter(profile.X_h, math.inf), 0.0)
    with pytest.raises(ValueError):
        profile.profile_values(1.0, 1.01)
    with pytest.raises(ValueError):
        profile.radial_derivatives(0.0, 0.0)
    with pytest.raises(ValueError):
        KokunoPA16CurrentJoinedProfile(quadrature_points=16)
    with pytest.raises(ValueError):
        KokunoPA16CurrentJoinedProfile(absolute_tolerance=0.0)


def test_public_repair_interval_constants_are_the_expected_source_annulus():
    assert LOG_REPAIR_START == -6.0
    assert LOG_JOIN_EXIT == -5.0
