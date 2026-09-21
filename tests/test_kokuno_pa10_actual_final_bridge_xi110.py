from __future__ import annotations

import copy

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_pa10_actual_final_bridge_xi110 import (
    A_FINAL,
    ANGULAR_SETTLE_LOG_WIDTH,
    AXIAL_SHUTDOWN_LOG_WIDTH,
    BRIDGE_QUADRATURE_ORDER,
    L_FINAL,
    KokunoPA10ActualFinalBridgeToXi110,
)


def test_source_geometry_and_exact_x100_value_derivative_handoff() -> None:
    bridge = KokunoPA10ActualFinalBridgeToXi110()
    assert bridge.X_b_ref == 100.0
    assert bridge.X_i == 110.0
    assert AXIAL_SHUTDOWN_LOG_WIDTH + ANGULAR_SETTLE_LOG_WIDTH < bridge.y_i

    eta = np.asarray([-0.6, -0.2, 0.0, 0.3, 0.6])
    X = np.full_like(eta, bridge.X_b_ref)
    got = bridge.values(X, eta)
    got_d = bridge.radial_derivatives(X, eta)
    parent = bridge.actual_x100.handoff_at_X100(eta)

    np.testing.assert_allclose(
        got["F_final_bridge"], parent["F_fixed_kappa"], rtol=0.0, atol=0.0
    )
    np.testing.assert_allclose(
        got["U_final_bridge"], parent["U_fixed_kappa"], rtol=0.0, atol=0.0
    )
    np.testing.assert_allclose(
        got["E_final_bridge"], parent["E_fixed_kappa"], rtol=0.0, atol=0.0
    )
    np.testing.assert_allclose(
        got_d["F_final_bridge_X"], parent["F_fixed_kappa_X"], rtol=2e-12, atol=2e-12
    )
    np.testing.assert_allclose(
        got_d["U_final_bridge_X"], parent["U_fixed_kappa_X"], rtol=2e-12, atol=2e-12
    )


def test_piecewise_source_schedule_and_final_slopes() -> None:
    bridge = KokunoPA10ActualFinalBridgeToXi110()
    eta = np.asarray([-0.45, 0.0, 0.45])

    start = bridge.schedule(np.full_like(eta, 100.0), eta)
    np.testing.assert_allclose(start["beta"], 1.0, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(start["angular_sigma"], 0.0, rtol=0.0, atol=0.0)

    axial_end = bridge.schedule(np.full_like(eta, bridge.X_axial_end), eta)
    np.testing.assert_allclose(axial_end["beta"], 0.0, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(axial_end["D_X_U"], 0.0, rtol=0.0, atol=0.0)

    angular_end = bridge.schedule(np.full_like(eta, bridge.X_angular_end), eta)
    np.testing.assert_allclose(angular_end["a"], A_FINAL, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(angular_end["l_profile"], L_FINAL, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(angular_end["D_X_U"], 0.0, rtol=0.0, atol=0.0)

    final = bridge.schedule(np.full_like(eta, 110.0), eta)
    np.testing.assert_allclose(final["a"], 0.8, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(final["l_profile"], 0.6, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(final["D_X_U"], 0.0, rtol=0.0, atol=0.0)


def test_radial_derivatives_against_independent_centered_x_difference() -> None:
    bridge = KokunoPA10ActualFinalBridgeToXi110()
    X = np.asarray([101.2, 103.0, 107.0])
    eta = np.asarray([-0.4, 0.05, 0.4])
    h = 2.0e-3
    plus = bridge.values(X + h, eta)
    minus = bridge.values(X - h, eta)
    deriv = bridge.radial_derivatives(X, eta)

    for value_key, derivative_key in (
        ("F_final_bridge", "F_final_bridge_X"),
        ("U_final_bridge", "U_final_bridge_X"),
        ("E_final_bridge", "E_final_bridge_X"),
    ):
        fd = (plus[value_key] - minus[value_key]) / (2.0 * h)
        np.testing.assert_allclose(deriv[derivative_key], fd, rtol=5e-3, atol=3e-7)


def test_quadrature_replay_and_nontrivial_Xi_boundary_profiles() -> None:
    bridge = KokunoPA10ActualFinalBridgeToXi110()
    eta = np.asarray([-0.6, -0.25, 0.0, 0.25, 0.6])
    X = np.full_like(eta, bridge.X_i)
    base = bridge.values(X, eta, order=BRIDGE_QUADRATURE_ORDER)
    replay = bridge.values(X, eta, order=64)
    for key in ("F_final_bridge", "U_final_bridge", "E_final_bridge"):
        np.testing.assert_allclose(base[key], replay[key], rtol=2e-3, atol=2e-8)
        assert np.all(np.isfinite(base[key]))

    handoff = bridge.handoff_at_Xi(eta)
    np.testing.assert_allclose(handoff["G_i"], handoff["U_final_bridge"], rtol=0.0, atol=0.0)
    np.testing.assert_allclose(
        handoff["ell_i"],
        np.log(bridge.C * handoff["E_final_bridge"]),
        rtol=0.0,
        atol=0.0,
    )
    np.testing.assert_allclose(
        bridge.X_i * handoff["F_final_bridge_X"] / handoff["F_final_bridge"],
        -0.5 * A_FINAL,
        rtol=2e-13,
        atol=2e-13,
    )
    np.testing.assert_allclose(handoff["U_final_bridge_X"], 0.0, rtol=0.0, atol=0.0)
    assert float(np.linalg.norm(handoff["G_i"])) > 1e-10
    assert float(np.linalg.norm(handoff["F_final_bridge"])) > 1e-10


def test_vectorization_serialization_guards_and_truth_boundary(tmp_path) -> None:
    bridge = KokunoPA10ActualFinalBridgeToXi110()
    X = np.asarray([[100.0], [104.5], [110.0]])
    eta = np.asarray([[-0.4, 0.0, 0.4]])
    got = bridge.values(X, eta)
    assert got["F_final_bridge"].shape == (3, 3)
    assert got["U_final_bridge"].shape == (3, 3)

    path = tmp_path / "actual_final_bridge_xi110.json"
    payload = bridge.save_configuration(path)
    loaded = KokunoPA10ActualFinalBridgeToXi110.load_configuration(path)
    assert loaded.configuration() == payload
    assert loaded.semantic_sha256 == bridge.semantic_sha256

    bad = copy.deepcopy(payload)
    bad["angular_settle_log_width"] = 0.0201
    with pytest.raises(ValueError):
        KokunoPA10ActualFinalBridgeToXi110.from_configuration(bad)
    with pytest.raises(ValueError):
        bridge.values(99.9, 0.0)
    with pytest.raises(ValueError):
        bridge.values(110.1, 0.0)
    with pytest.raises(ValueError):
        bridge.values(105.0, 1.0)

    truth = bridge.truth_boundary
    assert truth["public_final_bridge_formula_executable"]
    assert truth["candidate_side_actual_final_interpolation_100_to_Xi_materialized"]
    assert truth["candidate_side_actual_G_i_at_Xi_materialized"]
    assert truth["candidate_side_actual_ell_i_at_Xi_materialized"]
    assert not truth["source_prepared_full_reference_stress_pair_materialized"]
    assert not truth["source_admitted_global_kappa0_smallness"]
    assert not truth["final_bridge_cone_admissibility_independently_certified"]
    assert not truth["actual_upstream_five_moment_discrepancy_materialized"]
    assert not truth["outer_global_leading_velocity_materialized"]
    assert not truth["heldout_ns_residual_assessed"]
    assert not truth["pde_validated"]
    assert not truth["paper_exact"]
    assert not truth["openai_field_identified"]
