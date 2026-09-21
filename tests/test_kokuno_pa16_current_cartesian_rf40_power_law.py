import copy
import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_pa16_current_cartesian_rf40_power_law import (
    KokunoPA16CurrentCartesianRF40PowerLaw,
)


@pytest.fixture(scope="module")
def field():
    return KokunoPA16CurrentCartesianRF40PowerLaw()


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


def test_power_endpoint_matches_public_RF40c_schedule(field):
    assert field.log_X_3 == pytest.approx(
        field.outer_base.log_X_w, rel=0.0, abs=2e-12
    )
    assert field.log_X_4 == pytest.approx(
        field.outer_base.log_X_end, rel=0.0, abs=2e-12
    )
    assert field.log_X_4 - field.log_X_3 == pytest.approx(
        field.T_w, rel=0.0, abs=2e-12
    )
    assert field.X_4 / field.X_3 == pytest.approx(
        math.exp(field.T_w), rel=8e-13
    )


def test_inherited_domain_replays_exact_lambda_turn_candidate(field):
    X = field.X_3 * math.exp(-0.21)
    eta = -0.17
    got = field.similarity_profile_values(X, eta)
    old = field.current.similarity_profile_values(X, eta)
    for new_key, old_key in (
        ("F_current_rf40_power_law", "F_current_rf40_lambda_turn"),
        ("U_current_rf40_power_law", "U_current_rf40_lambda_turn"),
        ("E_current_rf40_power_law", "E_current_rf40_lambda_turn"),
        ("M_over_X_current_rf40_power_law", "M_over_X_current_rf40_lambda_turn"),
        ("M_eta_over_X_current_rf40_power_law", "M_eta_over_X_current_rf40_lambda_turn"),
        ("v0_current_rf40_power_law", "v0_current_rf40_lambda_turn"),
    ):
        np.testing.assert_array_equal(got[new_key], old[old_key])


def test_Xw_value_and_first_profile_jet_match_public_power_entry(field):
    eta = 0.29
    current = field.current.similarity_profile_values(field.X_3, eta)
    public = field.outer_base.profile_values(field.X_3, eta)
    current_jet = field.current.similarity_radial_derivatives(field.X_3, eta)

    assert _relative_error(current["F_current_rf40_lambda_turn"], public["F"]) <= 3e-8
    assert float(current["U_current_rf40_lambda_turn"]) == pytest.approx(
        float(public["U"]), rel=0.0, abs=3e-12
    )
    assert _relative_error(current["E_current_rf40_lambda_turn"], public["E"]) <= 3e-8
    assert _relative_error(
        current_jet["F_current_rf40_lambda_turn_X"], public["F_X"]
    ) <= 3e-8
    assert float(current_jet["U_current_rf40_lambda_turn_X"]) == pytest.approx(
        float(public["U_X"]), rel=0.0, abs=3e-12
    )
    assert _relative_error(
        current_jet["E_current_rf40_lambda_turn_X"], public["E_X"]
    ) <= 3e-8


def test_power_law_reuses_public_profile_and_constant_lambda(field):
    eta = -0.26
    X = field.X_3 * math.exp(0.37 * field.T_w)
    got = field.similarity_profile_values(X, eta)
    public = field.outer_base.profile_values(X, eta)

    assert str(np.asarray(public["stage"]).item()) == "power_law"
    assert str(np.asarray(got["region"]).item()) == "public_RF40_power_law"
    np.testing.assert_array_equal(got["F_current_rf40_power_law"], public["F"])
    np.testing.assert_array_equal(got["U_current_rf40_power_law"], public["U"])
    np.testing.assert_array_equal(got["E_current_rf40_power_law"], public["E"])
    assert float(got["U_current_rf40_power_law"]) == 0.0
    assert float(public["ell"]) == pytest.approx(
        -field.lambda_outer, rel=0.0, abs=2e-12
    )


def test_current_primitive_memory_decays_exactly_as_Xw_over_X(field):
    eta = 0.24
    seam = field.current.similarity_profile_values(field.X_3, eta)
    seam_m = float(seam["M_over_X_current_rf40_lambda_turn"])
    seam_me = float(seam["M_eta_over_X_current_rf40_lambda_turn"])

    for fraction in (0.08, 0.43, 0.81):
        X = field.X_3 * math.exp(fraction * field.T_w)
        got = field.similarity_profile_values(X, eta)
        ratio = field.X_3 / X
        assert float(got["M_over_X_current_rf40_power_law"]) == pytest.approx(
            ratio * seam_m, rel=5e-13, abs=1e-300
        )
        assert float(got["M_eta_over_X_current_rf40_power_law"]) == pytest.approx(
            ratio * seam_me, rel=5e-13, abs=1e-300
        )


def test_carried_physical_M_is_constant_and_log_derivative_is_zero(field):
    eta = 0.24
    X = field.X_3 * math.exp(0.55 * field.T_w)
    eps = 2.0e-5
    Xm = X * math.exp(-eps)
    Xp = X * math.exp(eps)
    minus = field.similarity_profile_values(Xm, eta)
    center = field.similarity_profile_values(X, eta)
    plus = field.similarity_profile_values(Xp, eta)

    Mm = Xm * float(minus["M_over_X_current_rf40_power_law"])
    Mc = X * float(center["M_over_X_current_rf40_power_law"])
    Mp = Xp * float(plus["M_over_X_current_rf40_power_law"])
    assert Mm == pytest.approx(Mc, rel=8e-13)
    assert Mp == pytest.approx(Mc, rel=8e-13)
    log_derivative = (Mp - Mm) / (2.0 * eps * Mc)
    assert log_derivative == pytest.approx(0.0, abs=2e-9)
    assert float(center["U_current_rf40_power_law"]) == 0.0


def test_power_law_radial_derivatives_reuse_public_and_fd_replay(field):
    X = field.X_3 * math.exp(0.34 * field.T_w)
    eta = 0.18
    got = field.similarity_radial_derivatives(X, eta)
    public = field.outer_base.profile_values(X, eta)
    np.testing.assert_array_equal(got["F_current_rf40_power_law_X"], public["F_X"])
    np.testing.assert_array_equal(got["U_current_rf40_power_law_X"], public["U_X"])
    np.testing.assert_array_equal(got["E_current_rf40_power_law_X"], public["E_X"])

    eps = 2.5e-5
    xm = X * math.exp(-eps)
    xp = X * math.exp(eps)
    fm = float(field.similarity_profile_values(xm, eta)["F_current_rf40_power_law"])
    fp = float(field.similarity_profile_values(xp, eta)["F_current_rf40_power_law"])
    fd = (fp - fm) / (xp - xm)
    assert fd == pytest.approx(
        float(got["F_current_rf40_power_law_X"]),
        rel=5e-7,
        abs=1e-300,
    )


def test_power_law_slopes_hold_at_interior_and_endpoint(field):
    eta = 0.23
    for fraction in (0.12, 0.67, 1.0):
        X = (
            field.X_4
            if fraction == 1.0
            else field.X_3 * math.exp(fraction * field.T_w)
        )
        vals = field.similarity_profile_values(X, eta)
        deriv = field.similarity_radial_derivatives(X, eta)
        public = field.outer_base.profile_values(X, eta)
        assert str(np.asarray(public["stage"]).item()) == "power_law"
        F = float(vals["F_current_rf40_power_law"])
        E = float(vals["E_current_rf40_power_law"])
        assert X * float(deriv["F_current_rf40_power_law_X"]) / F == pytest.approx(
            -1.0 - field.lambda_outer, abs=5e-12
        )
        assert X * float(deriv["E_current_rf40_power_law_X"]) / E == pytest.approx(
            -0.5 - field.lambda_outer, abs=5e-12
        )
        assert X * float(deriv["U_current_rf40_power_law_X"]) == pytest.approx(
            0.0, abs=2e-12
        )


def test_velocity_vectorizes_across_Xw_and_axis_is_regular(field):
    X = np.asarray(
        [
            field.X_3 * math.exp(-0.2),
            field.X_3,
            field.X_3 * math.exp(0.08 * field.T_w),
            field.X_3 * math.exp(0.42 * field.T_w),
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
        [field.velocity(x[i], y[i], z[i], t[i]) for i in range(X.size)],
        axis=0,
    )
    np.testing.assert_allclose(velocity, scalar, rtol=2e-12, atol=2e-12)

    axis_velocity = field.velocity(0.0, 0.0, 0.12, 0.5)
    assert np.all(np.isfinite(axis_velocity))
    assert float(axis_velocity[0]) == 0.0
    assert float(axis_velocity[1]) == 0.0


def test_fails_closed_after_power_law_endpoint(field):
    xb, yb, zb, tb = _cartesian_probe(
        field, field.X_4 * 1.001, 0.2, 0.5, 0.3
    )
    with pytest.raises(ValueError, match="beyond current RF40 power-law endpoint"):
        field.velocity(xb, yb, zb, tb)
    with pytest.raises(ValueError):
        field.similarity_profile_values(field.X_4 * 1.0001, 0.0)


def test_configuration_roundtrip_hash_and_truth_boundary(field, tmp_path):
    path = tmp_path / "current-cartesian-rf40-power-law.json"
    field.save_configuration(path)
    rebound = KokunoPA16CurrentCartesianRF40PowerLaw.load_configuration(path)
    assert rebound.configuration() == field.configuration()
    assert rebound.semantic_sha256 == field.semantic_sha256

    truth = field.truth_boundary
    assert truth["current_cartesian_leading_through_RF40_lambda_turn_consumed"] is True
    assert truth["public_RF40_power_law_formula_reused"] is True
    assert truth["current_lineage_RF40_power_law_materialized"] is True
    assert truth["current_RF40_base_backbone_through_power_law_endpoint_materialized"] is True
    assert truth["current_incompressibility_memory_carried_through_RF40_power_law"] is True
    assert truth["full_post_XR_RF40_current_lineage_materialized"] is False
    assert truth["current_lineage_cone_I1_I2_I3_I4_outer_overlays_completed"] is False
    assert truth["source_terminal_tail_schedule_bound"] is False
    assert truth["outer_global_leading_velocity_materialized"] is False
    assert truth["unified_global_cartesian_velocity_export_ready"] is False
    assert truth["heldout_ns_residual_assessed"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False

    payload = copy.deepcopy(field.configuration())
    payload["schema"] = "wrong-schema"
    with pytest.raises(ValueError, match="schema"):
        KokunoPA16CurrentCartesianRF40PowerLaw.from_configuration(payload)


def test_similarity_domain_guards(field):
    with pytest.raises(ValueError):
        field.similarity_profile_values(-1.0, 0.0)
    with pytest.raises(ValueError):
        field.similarity_profile_values(field.X_3, 1.01)
    with pytest.raises(ValueError):
        field.similarity_radial_derivatives(0.0, 0.0)
