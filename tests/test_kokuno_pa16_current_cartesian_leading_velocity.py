import copy
import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_pa16_current_cartesian_leading_velocity import (
    KokunoPA16CurrentCartesianLeadingVelocity,
)


@pytest.fixture(scope="module")
def field():
    return KokunoPA16CurrentCartesianLeadingVelocity(moment_quadrature_order=48)


def _cartesian_probe(field, X, eta, t, theta=0.0):
    X, eta, t, theta = np.broadcast_arrays(
        np.asarray(X, dtype=float),
        np.asarray(eta, dtype=float),
        np.asarray(t, dtype=float),
        np.asarray(theta, dtype=float),
    )
    tau = 1.0 - t
    q = tau / (1.0 - eta * eta)
    z = np.power(q, field.D) * eta
    r = np.sqrt(2.0 * q * X)
    return r * np.cos(theta), r * np.sin(theta), z, t


def test_natural_lane_replays_existing_governed_center_v0(field):
    X = 0.45 * field.joined.moments.X_0
    eta = 0.31
    got = field.similarity_profile_values(X, eta)
    old = field.geometry.source_center.values(X, eta)

    assert float(got["F_current_joined"]) == pytest.approx(
        float(old["F_0"]), rel=3e-12, abs=2e-14
    )
    assert float(got["U_current_joined"]) == pytest.approx(
        float(old["U_0"]), rel=3e-12, abs=2e-14
    )
    assert float(got["M_over_X_current_joined"]) == pytest.approx(
        float(old["M_0_over_X"]), rel=2e-10, abs=2e-12
    )
    assert float(got["v0_current_joined"]) == pytest.approx(
        float(old["v_0"]), rel=2e-7, abs=2e-9
    )


def test_M_primitive_continues_across_Xi_and_dM_dX_replays_U(field):
    eta = -0.27
    X = 2.0 * field.X_i
    h = 2.0e-5 * X

    minus = field.similarity_profile_values(X - h, eta)
    center = field.similarity_profile_values(X, eta)
    plus = field.similarity_profile_values(X + h, eta)

    M_minus = (X - h) * float(minus["M_over_X_current_joined"])
    M_plus = (X + h) * float(plus["M_over_X_current_joined"])
    M_X_fd = (M_plus - M_minus) / (2.0 * h)
    assert M_X_fd == pytest.approx(
        float(center["U_current_joined"]), rel=2e-8, abs=2e-9
    )

    left = field.similarity_profile_values(field.X_i, eta)
    right = field.similarity_profile_values(np.nextafter(field.X_i, math.inf), eta)
    assert float(right["M_over_X_current_joined"]) == pytest.approx(
        float(left["M_over_X_current_joined"]), rel=3e-12, abs=3e-12
    )


def test_velocity_interface_vectorizes_and_reuses_governed_q_inverse(field):
    X = np.asarray([0.0, 0.4 * field.joined.moments.X_0, 2.0 * field.X_i])
    eta = np.asarray([-0.35, 0.15, 0.42])
    t = np.asarray([0.3, 0.5, 0.7])
    theta = np.asarray([0.0, 0.7, 1.4])
    x, y, z, t = _cartesian_probe(field, X, eta, t, theta)

    coords = field.similarity_coordinates(x, y, z, t)
    geometry_coords = field.geometry.similarity_coordinates(x, y, z, t)
    np.testing.assert_array_equal(coords["q"], geometry_coords["q"])
    np.testing.assert_array_equal(coords["X"], geometry_coords["X"])
    np.testing.assert_array_equal(coords["eta"], geometry_coords["eta"])
    np.testing.assert_allclose(coords["X"], X, rtol=4e-14, atol=2e-13)
    np.testing.assert_allclose(coords["eta"], eta, rtol=4e-14, atol=2e-13)

    velocity = field.velocity(x, y, z, t)
    assert velocity.shape == (3, 3)
    assert np.all(np.isfinite(velocity))
    assert np.any(np.abs(velocity) > 0.0)


def test_cartesian_axis_transverse_velocity_is_exact_zero(field):
    eta = np.asarray([-0.6, 0.0, 0.6])
    t = np.asarray([0.3, 0.5, 0.7])
    x, y, z, t = _cartesian_probe(field, np.zeros(3), eta, t)
    velocity = field.velocity(x, y, z, t)

    assert np.all(velocity[:, :2] == 0.0)
    assert np.all(np.isfinite(velocity))
    assert np.any(np.abs(velocity[:, 2]) > 0.0)


def test_beyond_Xh_fails_closed(field):
    X = 1.01 * field.X_h
    x, y, z, t = _cartesian_probe(field, X, 0.2, 0.5, 0.3)
    with pytest.raises(ValueError, match="beyond current X_h"):
        field.velocity(x, y, z, t)


def test_configuration_roundtrip_hash_and_truth_boundary(field, tmp_path):
    path = tmp_path / "current-cartesian-leading.json"
    field.save_configuration(path)
    rebound = KokunoPA16CurrentCartesianLeadingVelocity.load_configuration(path)
    assert rebound.configuration() == field.configuration()
    assert rebound.semantic_sha256 == field.semantic_sha256

    truth = field.truth_boundary
    assert truth["governed_source_coordinate_inverse_reused"] is True
    assert truth["current_joined_incompressibility_primitive_materialized"] is True
    assert truth["current_joined_cartesian_spacetime_leading_velocity_through_Xh_materialized"] is True
    assert truth["eta_derivative_is_repository_numerical_realization"] is True
    assert truth["outer_global_leading_velocity_materialized"] is False
    assert truth["velocity_beyond_Xh_materialized"] is False
    assert truth["unified_global_cartesian_velocity_export_ready"] is False
    assert truth["matched_global_pressure_materialized"] is False
    assert truth["restricted_forcing_materialized"] is False
    assert truth["heldout_ns_residual_assessed"] is False
    assert truth["same_protocol_comparable_to_st006"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False

    payload = copy.deepcopy(field.configuration())
    payload["schema"] = "wrong-schema"
    with pytest.raises(ValueError, match="schema"):
        KokunoPA16CurrentCartesianLeadingVelocity.from_configuration(payload)


def test_parameter_guards(field):
    with pytest.raises(ValueError):
        KokunoPA16CurrentCartesianLeadingVelocity(moment_quadrature_order=8)
    with pytest.raises(ValueError):
        KokunoPA16CurrentCartesianLeadingVelocity(eta_fd_step=0.0)
    with pytest.raises(ValueError):
        field.similarity_profile_values(-1.0, 0.0)
    with pytest.raises(ValueError):
        field.similarity_profile_values(1.0, 1.01)
