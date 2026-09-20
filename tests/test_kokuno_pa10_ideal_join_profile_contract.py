from __future__ import annotations

import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_pa10_ideal_join_profile_contract import (
    DEFAULT_P_STAR,
    KokunoPA10IdealJoinProfileContract,
)


def _relative(actual: np.ndarray, expected: np.ndarray) -> float:
    actual = np.asarray(actual, dtype=float)
    expected = np.asarray(expected, dtype=float)
    scale = np.maximum(np.abs(expected), np.finfo(float).tiny)
    return float(np.max(np.abs(actual - expected) / scale))


def test_source_scale_and_join_are_literal_public_formulas() -> None:
    profile = KokunoPA10IdealJoinProfileContract()
    assert profile.C == 1000.0
    assert profile.P_star == DEFAULT_P_STAR == 1.0
    assert profile.X_R == pytest.approx(110.0 * (profile.C * profile.P_star) ** 10, rel=2e-15)
    assert profile.X_h / profile.X_R == pytest.approx(math.exp(-5.0), rel=2e-15)


def test_vectorized_ideal_pair_and_E_F_identity() -> None:
    profile = KokunoPA10IdealJoinProfileContract()
    X = profile.X_R * np.array([[math.exp(-6.0)], [math.exp(-5.5)], [math.exp(-5.0)]])
    eta = np.array([[-0.8, 0.0, 0.6]])
    out = profile.values(X, eta)
    assert out["U_ideal"].shape == (3, 3)
    assert np.allclose(out["U_ideal"], 4.0 * np.broadcast_to(eta, (3, 3)), rtol=0.0, atol=0.0)
    np.testing.assert_allclose(
        out["E_ideal"], np.sqrt(2.0 * out["X"]) * out["F_ideal"], rtol=2e-15, atol=0.0
    )
    assert float(np.max(np.abs(out["E_ideal"]))) > 0.0


def test_analytic_first_derivatives_against_independent_centered_differences() -> None:
    profile = KokunoPA10IdealJoinProfileContract()
    X = profile.X_R * np.array([math.exp(-5.8), math.exp(-5.3), math.exp(-5.05)])
    eta = np.array([-0.55, 0.15, 0.72])
    deriv = profile.derivatives(X, eta)

    hX = 1.0e-5 * X
    plus = profile.values(X + hX, eta)
    minus = profile.values(X - hX, eta)
    E_X_fd = (plus["E_ideal"] - minus["E_ideal"]) / (2.0 * hX)
    F_X_fd = (plus["F_ideal"] - minus["F_ideal"]) / (2.0 * hX)

    he = 1.0e-6
    plus_e = profile.values(X, eta + he)
    minus_e = profile.values(X, eta - he)
    E_eta_fd = (plus_e["E_ideal"] - minus_e["E_ideal"]) / (2.0 * he)
    F_eta_fd = (plus_e["F_ideal"] - minus_e["F_ideal"]) / (2.0 * he)

    assert _relative(deriv["E_ideal_X"], E_X_fd) < 2e-8
    assert _relative(deriv["F_ideal_X"], F_X_fd) < 2e-8
    assert _relative(deriv["E_ideal_eta"], E_eta_fd) < 2e-8
    assert _relative(deriv["F_ideal_eta"], F_eta_fd) < 2e-8

    values = profile.values(X, eta)
    np.testing.assert_allclose(X * deriv["E_ideal_X"] / values["E_ideal"], 0.1, rtol=2e-15, atol=2e-15)
    np.testing.assert_allclose(X * deriv["F_ideal_X"] / values["F_ideal"], -0.4, rtol=2e-15, atol=2e-15)


def test_five_normalized_moment_targets_at_public_join() -> None:
    profile = KokunoPA10IdealJoinProfileContract()
    eta = np.array([-0.4, 0.0, 0.7])
    X = np.full_like(eta, profile.X_h)
    x = math.exp(-5.0)
    f = 1.0 / (1.0 + eta * eta)
    m = profile.normalized_moments(X, eta)

    np.testing.assert_allclose(m["hat_M"], 4.0 * eta * x, rtol=2e-15, atol=0.0)
    expected_I = (5.0 * math.sqrt(2.0) / 8.0) * profile.P_star * f * x ** (8.0 / 5.0)
    np.testing.assert_allclose(m["hat_I"], expected_I, rtol=2e-15, atol=0.0)
    np.testing.assert_allclose(m["hat_J"], 4.0 * eta * expected_I, rtol=2e-15, atol=0.0)
    np.testing.assert_allclose(
        m["hat_S"],
        16.0 * eta * eta * x - (5.0 / 12.0) * profile.P_star**2 * f**2 * x ** (6.0 / 5.0),
        rtol=2e-15,
        atol=0.0,
    )
    np.testing.assert_allclose(
        m["C_p"], (5.0 / 2.0) * profile.P_star**2 * f**2 * x ** (1.0 / 5.0), rtol=2e-15, atol=0.0
    )


def test_join_gap_is_explicit_not_a_fabricated_interpolation() -> None:
    profile = KokunoPA10IdealJoinProfileContract()
    gap = profile.join_gap_diagnostic()
    assert gap["source_join_is_outside_current_inner_domain"] is True
    assert gap["X_h_over_inner_X_max"] > 1.0
    assert gap["log_X_h_over_inner_X_max"] > 0.0
    assert gap["inner_X_max"] == pytest.approx(profile.inner_X_interval[1], rel=0.0, abs=0.0)


def test_serialization_hash_and_truth_boundaries(tmp_path) -> None:
    profile = KokunoPA10IdealJoinProfileContract()
    path = tmp_path / "ideal_join.json"
    profile.save_configuration(path)
    loaded = KokunoPA10IdealJoinProfileContract.load_configuration(path)
    assert loaded.configuration() == profile.configuration()
    assert loaded.semantic_sha256 == profile.semantic_sha256
    assert KokunoPA10IdealJoinProfileContract(P_star=1.1).semantic_sha256 != profile.semantic_sha256

    truth = profile.truth_boundary
    assert truth["public_ideal_join_target_executable"] is True
    assert truth["public_five_normalized_moment_targets_executable"] is True
    assert truth["P_star_is_repository_autonomous_positive_realization"] is True
    assert truth["source_P_star_numerically_identified"] is False
    assert truth["five_bump_nonlinear_repair_materialized"] is False
    assert truth["source_join_neighborhood_width_known"] is False
    assert truth["outer_global_leading_velocity_materialized"] is False
    assert truth["heldout_ns_residual_assessed"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False


def test_invalid_inputs_fail_closed() -> None:
    with pytest.raises(ValueError):
        KokunoPA10IdealJoinProfileContract(P_star=0.0)
    with pytest.raises(ValueError):
        KokunoPA10IdealJoinProfileContract(P_star=-1.0)
    profile = KokunoPA10IdealJoinProfileContract()
    with pytest.raises(ValueError):
        profile.values(0.0, 0.0)
    with pytest.raises(ValueError):
        profile.values(-1.0, 0.0)
