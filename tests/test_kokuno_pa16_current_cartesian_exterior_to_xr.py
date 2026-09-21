import copy
import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_pa16_current_cartesian_exterior_to_xr import (
    KokunoPA16CurrentCartesianExteriorToXR,
)


@pytest.fixture(scope="module")
def field():
    return KokunoPA16CurrentCartesianExteriorToXR()


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


def _relative_error(a, b):
    a = float(np.asarray(a))
    b = float(np.asarray(b))
    return abs(a - b) / max(abs(a), abs(b), 1.0e-300)


def test_current_outer_scale_identity_is_reused(field):
    assert math.log(field.X_R) == pytest.approx(
        field.outer_base.log_X_R, rel=0.0, abs=2e-12
    )
    assert field.X_h / field.X_R == pytest.approx(math.exp(-5.0), rel=2e-14)
    assert field.h == pytest.approx(field.outer_base.h, rel=0.0, abs=2e-15)
    assert (
        float(field.current.joined.moments.outer_schedule.C)
        == float(field.outer_base.outer_schedule.C)
    )


def test_inner_domain_replays_exact_current_cartesian_field(field):
    X = 2.0 * field.current.X_i
    eta = -0.23
    t = 0.47
    theta = 0.81
    x, y, z, t = _cartesian_probe(field, X, eta, t, theta)

    got_profile = field.similarity_profile_values(X, eta)
    old_profile = field.current.similarity_profile_values(X, eta)
    np.testing.assert_array_equal(
        got_profile["F_current_exterior"], old_profile["F_current_joined"]
    )
    np.testing.assert_array_equal(
        got_profile["U_current_exterior"], old_profile["U_current_joined"]
    )
    np.testing.assert_array_equal(
        got_profile["M_over_X_current_exterior"], old_profile["M_over_X_current_joined"]
    )
    np.testing.assert_array_equal(
        got_profile["v0_current_exterior"], old_profile["v0_current_joined"]
    )
    np.testing.assert_allclose(
        field.velocity(x, y, z, t), field.current.velocity(x, y, z, t), rtol=0.0, atol=0.0
    )


def test_Xh_value_and_first_profile_jet_match_public_ideal_branch(field):
    eta = 0.31
    current = field.current.similarity_profile_values(field.X_h, eta)
    ideal = field.outer_base.profile_values(field.X_h, eta)
    current_jet = field.current.joined.radial_derivatives(field.X_h, eta)

    assert _relative_error(current["F_current_joined"], ideal["F"]) <= 2e-8
    assert float(current["U_current_joined"]) == pytest.approx(
        float(ideal["U"]), rel=2e-8, abs=2e-10
    )
    assert _relative_error(current["E_current_joined"], ideal["E"]) <= 2e-8

    assert _relative_error(current_jet["F_current_joined_X"], ideal["F_X"]) <= 2e-8
    assert float(current_jet["U_current_joined_X"]) == pytest.approx(
        float(ideal["U_X"]), rel=2e-8, abs=2e-12
    )
    assert _relative_error(current_jet["E_current_joined_X"], ideal["E_X"]) <= 2e-8

    unified = field.similarity_profile_values(field.X_h, eta)
    np.testing.assert_array_equal(
        unified["M_over_X_current_exterior"], current["M_over_X_current_joined"]
    )
    np.testing.assert_array_equal(
        unified["M_eta_over_X_current_exterior"], current["M_eta_over_X_current_joined"]
    )
    np.testing.assert_array_equal(
        unified["v0_current_exterior"], current["v0_current_joined"]
    )


def test_exterior_reuses_ideal_F_U_E_and_carries_actual_primitive_memory(field):
    eta = -0.27
    X = field.X_h * math.exp(2.3)
    got = field.similarity_profile_values(X, eta)
    ideal = field.outer_base.profile_values(X, eta)
    seam = field.current.similarity_profile_values(field.X_h, eta)

    np.testing.assert_array_equal(got["F_current_exterior"], ideal["F"])
    np.testing.assert_array_equal(got["U_current_exterior"], ideal["U"])
    np.testing.assert_array_equal(got["E_current_exterior"], ideal["E"])
    assert str(np.asarray(got["region"]).item()) == "public_exterior_preservation_to_XR"

    r = field.X_h / X
    expected_m = 4.0 * eta + r * (
        float(seam["M_over_X_current_joined"]) - 4.0 * eta
    )
    expected_m_eta = 4.0 + r * (
        float(seam["M_eta_over_X_current_joined"]) - 4.0
    )
    assert float(got["M_over_X_current_exterior"]) == pytest.approx(
        expected_m, rel=2e-13, abs=2e-13
    )
    assert float(got["M_eta_over_X_current_exterior"]) == pytest.approx(
        expected_m_eta, rel=2e-13, abs=2e-13
    )


def test_exterior_carried_primitive_replays_dM_dX_equals_U(field):
    eta = 0.22
    X = field.X_h * math.exp(3.0)
    h = 2.0e-5
    Xm = X * math.exp(-h)
    Xp = X * math.exp(h)
    minus = field.similarity_profile_values(Xm, eta)
    center = field.similarity_profile_values(X, eta)
    plus = field.similarity_profile_values(Xp, eta)

    Mm = Xm * float(minus["M_over_X_current_exterior"])
    Mp = Xp * float(plus["M_over_X_current_exterior"])
    dM = (Mp - Mm) / (Xp - Xm)
    assert dM == pytest.approx(
        float(center["U_current_exterior"]), rel=3e-9, abs=3e-10
    )


def test_outer_radial_derivatives_are_public_reference_derivatives(field):
    X = field.X_h * math.exp(1.7)
    eta = 0.19
    got = field.similarity_radial_derivatives(X, eta)
    ideal = field.outer_base.profile_values(X, eta)
    np.testing.assert_array_equal(got["F_current_exterior_X"], ideal["F_X"])
    np.testing.assert_array_equal(got["U_current_exterior_X"], ideal["U_X"])
    np.testing.assert_array_equal(got["E_current_exterior_X"], ideal["E_X"])


def test_velocity_vectorizes_across_Xh_and_fails_closed_after_XR(field):
    X = np.asarray(
        [
            0.4 * field.current.joined.moments.X_0,
            field.X_h,
            field.X_h * math.exp(1.0),
            field.X_R,
        ]
    )
    eta = np.asarray([-0.35, -0.1, 0.2, 0.45])
    t = np.asarray([0.31, 0.43, 0.57, 0.69])
    theta = np.asarray([0.0, 0.4, 0.8, 1.2])
    x, y, z, t = _cartesian_probe(field, X, eta, t, theta)
    velocity = field.velocity(x, y, z, t)
    assert velocity.shape == (4, 3)
    assert np.all(np.isfinite(velocity))
    assert np.any(np.abs(velocity) > 0.0)

    scalar = np.stack(
        [field.velocity(x[i], y[i], z[i], t[i]) for i in range(X.size)], axis=0
    )
    np.testing.assert_allclose(velocity, scalar, rtol=2e-13, atol=2e-13)

    xb, yb, zb, tb = _cartesian_probe(field, field.X_R * 1.001, 0.2, 0.5, 0.3)
    with pytest.raises(ValueError, match="beyond current X_R"):
        field.velocity(xb, yb, zb, tb)


def test_configuration_roundtrip_hash_and_truth_boundary(field, tmp_path):
    path = tmp_path / "current-cartesian-exterior-to-xr.json"
    field.save_configuration(path)
    rebound = KokunoPA16CurrentCartesianExteriorToXR.load_configuration(path)
    assert rebound.configuration() == field.configuration()
    assert rebound.semantic_sha256 == field.semantic_sha256

    truth = field.truth_boundary
    assert truth["current_cartesian_leading_through_Xh_consumed"] is True
    assert truth["candidate_side_exterior_preservation_Xh_to_XR_materialized"] is True
    assert truth["current_cartesian_spacetime_leading_velocity_through_XR_materialized"] is True
    assert truth["velocity_beyond_Xh_materialized"] is True
    assert truth["current_incompressibility_memory_carried"] is True
    assert truth["source_exact_exterior_preservation_certified"] is False
    assert truth["post_XR_RF40_current_lineage_materialized"] is False
    assert truth["outer_global_leading_velocity_materialized"] is False
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
        KokunoPA16CurrentCartesianExteriorToXR.from_configuration(payload)


def test_similarity_domain_guards(field):
    with pytest.raises(ValueError):
        field.similarity_profile_values(-1.0, 0.0)
    with pytest.raises(ValueError):
        field.similarity_profile_values(field.X_R * 1.0001, 0.0)
    with pytest.raises(ValueError):
        field.similarity_profile_values(field.X_h, 1.01)
    with pytest.raises(ValueError):
        field.similarity_radial_derivatives(0.0, 0.0)
