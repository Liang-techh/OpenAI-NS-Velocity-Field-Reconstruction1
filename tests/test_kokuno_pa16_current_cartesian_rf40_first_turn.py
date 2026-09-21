import copy
import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_pa16_current_cartesian_rf40_first_turn import (
    KokunoPA16CurrentCartesianRF40FirstTurn,
)


@pytest.fixture(scope="module")
def field():
    return KokunoPA16CurrentCartesianRF40FirstTurn()


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


def test_current_XR_scale_and_first_turn_endpoint_are_reused(field):
    assert math.log(field.X_R) == pytest.approx(
        field.outer_base.log_X_R, rel=0.0, abs=2e-12
    )
    assert field.log_X_1 == pytest.approx(
        field.outer_base.log_stage1_end, rel=0.0, abs=2e-12
    )
    assert field.X_1 / field.X_R == pytest.approx(math.e, rel=8e-14)
    assert field.h == pytest.approx(field.outer_base.h, rel=0.0, abs=2e-15)


def test_inherited_domain_replays_exact_current_exterior_field(field):
    X = field.current.X_h * math.exp(2.0)
    eta = -0.21
    t = 0.47
    theta = 0.63
    x, y, z, t = _cartesian_probe(field, X, eta, t, theta)

    got = field.similarity_profile_values(X, eta)
    old = field.current.similarity_profile_values(X, eta)
    np.testing.assert_array_equal(
        got["F_current_rf40_first_turn"], old["F_current_exterior"]
    )
    np.testing.assert_array_equal(
        got["U_current_rf40_first_turn"], old["U_current_exterior"]
    )
    np.testing.assert_array_equal(
        got["M_over_X_current_rf40_first_turn"], old["M_over_X_current_exterior"]
    )
    np.testing.assert_array_equal(
        got["v0_current_rf40_first_turn"], old["v0_current_exterior"]
    )
    np.testing.assert_allclose(
        field.velocity(x, y, z, t), field.current.velocity(x, y, z, t),
        rtol=0.0, atol=0.0,
    )


def test_XR_value_and_first_profile_jet_match_public_RF40_right_branch(field):
    eta = 0.29
    current = field.current.similarity_profile_values(field.X_R, eta)
    public = field.outer_base.profile_values(field.X_R, eta)
    current_jet = field.current.similarity_radial_derivatives(field.X_R, eta)

    assert _relative_error(current["F_current_exterior"], public["F"]) <= 2e-8
    assert float(current["U_current_exterior"]) == pytest.approx(
        float(public["U"]), rel=2e-8, abs=2e-10
    )
    assert _relative_error(current["E_current_exterior"], public["E"]) <= 2e-8

    assert _relative_error(current_jet["F_current_exterior_X"], public["F_X"]) <= 2e-8
    assert float(current_jet["U_current_exterior_X"]) == pytest.approx(
        float(public["U_X"]), rel=2e-8, abs=2e-12
    )
    assert _relative_error(current_jet["E_current_exterior_X"], public["E_X"]) <= 2e-8

    unified = field.similarity_profile_values(field.X_R, eta)
    np.testing.assert_array_equal(
        unified["M_over_X_current_rf40_first_turn"],
        current["M_over_X_current_exterior"],
    )
    np.testing.assert_array_equal(
        unified["M_eta_over_X_current_rf40_first_turn"],
        current["M_eta_over_X_current_exterior"],
    )


def test_RF40_first_turn_reuses_public_F_U_E_and_carries_current_memory(field):
    eta = -0.26
    X = field.X_R * math.exp(0.43)
    got = field.similarity_profile_values(X, eta)
    public = field.outer_base.profile_values(X, eta)
    seam = field.current.similarity_profile_values(field.X_R, eta)

    assert str(np.asarray(public["stage"]).item()) == "first_l_turn"
    assert str(np.asarray(got["region"]).item()) == "public_RF40_first_l_turn"
    np.testing.assert_array_equal(got["F_current_rf40_first_turn"], public["F"])
    np.testing.assert_array_equal(got["U_current_rf40_first_turn"], public["U"])
    np.testing.assert_array_equal(got["E_current_rf40_first_turn"], public["E"])

    ratio = field.X_R / X
    expected_m = 4.0 * eta + ratio * (
        float(seam["M_over_X_current_exterior"]) - 4.0 * eta
    )
    expected_m_eta = 4.0 + ratio * (
        float(seam["M_eta_over_X_current_exterior"]) - 4.0
    )
    assert float(got["M_over_X_current_rf40_first_turn"]) == pytest.approx(
        expected_m, rel=2e-13, abs=2e-13
    )
    assert float(got["M_eta_over_X_current_rf40_first_turn"]) == pytest.approx(
        expected_m_eta, rel=2e-13, abs=2e-13
    )


def test_carried_primitive_independently_replays_dM_dX_equals_U(field):
    eta = 0.24
    X = field.X_R * math.exp(0.51)
    h = 2.0e-5
    Xm = X * math.exp(-h)
    Xp = X * math.exp(h)
    minus = field.similarity_profile_values(Xm, eta)
    center = field.similarity_profile_values(X, eta)
    plus = field.similarity_profile_values(Xp, eta)

    Mm = Xm * float(minus["M_over_X_current_rf40_first_turn"])
    Mp = Xp * float(plus["M_over_X_current_rf40_first_turn"])
    dM = (Mp - Mm) / (Xp - Xm)
    assert dM == pytest.approx(
        float(center["U_current_rf40_first_turn"]), rel=3e-9, abs=3e-10
    )


def test_first_turn_radial_derivatives_are_public_RF40_derivatives(field):
    X = field.X_R * math.exp(0.61)
    eta = 0.18
    got = field.similarity_radial_derivatives(X, eta)
    public = field.outer_base.profile_values(X, eta)
    np.testing.assert_array_equal(got["F_current_rf40_first_turn_X"], public["F_X"])
    np.testing.assert_array_equal(got["U_current_rf40_first_turn_X"], public["U_X"])
    np.testing.assert_array_equal(got["E_current_rf40_first_turn_X"], public["E_X"])

    # A value-only replay would miss a wrong derivative implementation.
    eps = 2.5e-5
    xm = X * math.exp(-eps)
    xp = X * math.exp(eps)
    fm = float(field.similarity_profile_values(xm, eta)["F_current_rf40_first_turn"])
    fp = float(field.similarity_profile_values(xp, eta)["F_current_rf40_first_turn"])
    fd = (fp - fm) / (xp - xm)
    assert fd == pytest.approx(
        float(got["F_current_rf40_first_turn_X"]), rel=2e-7, abs=1e-300
    )


def test_first_turn_endpoint_has_public_slopes(field):
    # Stay a few ulps inside the first-turn endpoint so this regression tests
    # the first-turn formula rather than relying on stage-label tie breaking.
    X = field.X_1 * math.exp(-2.0e-10)
    eta = 0.23
    vals = field.similarity_profile_values(X, eta)
    deriv = field.similarity_radial_derivatives(X, eta)
    F = float(vals["F_current_rf40_first_turn"])
    E = float(vals["E_current_rf40_first_turn"])
    assert X * float(deriv["F_current_rf40_first_turn_X"]) / F == pytest.approx(
        -1.0, abs=5e-8
    )
    assert X * float(deriv["E_current_rf40_first_turn_X"]) / E == pytest.approx(
        -0.5, abs=5e-8
    )
    assert X * float(deriv["U_current_rf40_first_turn_X"]) == pytest.approx(
        0.0, abs=2e-12
    )


def test_velocity_vectorizes_across_XR_and_axis_is_regular(field):
    X = np.asarray(
        [
            0.4 * field.current.current.joined.moments.X_0,
            0.8 * field.X_R,
            field.X_R,
            field.X_R * math.exp(0.37),
        ]
    )
    eta = np.asarray([-0.32, -0.08, 0.19, 0.41])
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

    # Exact Cartesian axis: transverse components vanish structurally.
    axis_velocity = field.velocity(0.0, 0.0, 0.12, 0.5)
    assert np.all(np.isfinite(axis_velocity))
    assert float(axis_velocity[0]) == 0.0
    assert float(axis_velocity[1]) == 0.0


def test_fails_closed_after_first_turn_endpoint(field):
    xb, yb, zb, tb = _cartesian_probe(
        field, field.X_1 * 1.001, 0.2, 0.5, 0.3
    )
    with pytest.raises(ValueError, match="beyond current RF40 first-turn endpoint"):
        field.velocity(xb, yb, zb, tb)
    with pytest.raises(ValueError):
        field.similarity_profile_values(field.X_1 * 1.0001, 0.0)


def test_configuration_roundtrip_hash_and_truth_boundary(field, tmp_path):
    path = tmp_path / "current-cartesian-rf40-first-turn.json"
    field.save_configuration(path)
    rebound = KokunoPA16CurrentCartesianRF40FirstTurn.load_configuration(path)
    assert rebound.configuration() == field.configuration()
    assert rebound.semantic_sha256 == field.semantic_sha256

    truth = field.truth_boundary
    assert truth["current_cartesian_leading_through_XR_consumed"] is True
    assert truth["public_RF40_first_turn_formula_reused"] is True
    assert truth["current_lineage_RF40_first_turn_materialized"] is True
    assert truth["current_incompressibility_memory_carried_through_RF40_first_turn"] is True
    assert truth["full_post_XR_RF40_current_lineage_materialized"] is False
    assert truth["RF40_axial_shutdown_current_lineage_materialized"] is False
    assert truth["current_lineage_cone_I1_I2_I3_I4_outer_overlays_completed"] is False
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
        KokunoPA16CurrentCartesianRF40FirstTurn.from_configuration(payload)


def test_similarity_domain_guards(field):
    with pytest.raises(ValueError):
        field.similarity_profile_values(-1.0, 0.0)
    with pytest.raises(ValueError):
        field.similarity_profile_values(field.X_R, 1.01)
    with pytest.raises(ValueError):
        field.similarity_radial_derivatives(0.0, 0.0)
