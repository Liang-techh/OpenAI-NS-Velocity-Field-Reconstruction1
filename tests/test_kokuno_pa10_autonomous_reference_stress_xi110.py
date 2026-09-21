from __future__ import annotations

import copy

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_pa10_autonomous_reference_stress_xi110 import (
    STRESS_QUADRATURE_ORDER,
    KokunoPA10AutonomousReferenceStressToXi110,
)


def test_source_geometry_and_exact_x100_handoff() -> None:
    ext = KokunoPA10AutonomousReferenceStressToXi110()
    assert ext.X_b_ref == 100.0
    assert ext.X_i == 110.0

    eta = np.asarray([-0.6, -0.2, 0.0, 0.35, 0.65])
    X = np.full_like(eta, 100.0)
    got = ext.values(X, eta)
    parent_ns = ext.parent.values(X, eta)
    parent_p1 = ext.parent.p1_reference.values(X, eta)["p1_reference"]
    parent_ref = ext.parent.reference.values(X, eta)
    parent_pi = ext.parent.pressure.pressure(X, eta)

    np.testing.assert_allclose(got["p1_reference"], parent_p1, rtol=0.0, atol=2e-12)
    np.testing.assert_allclose(
        got["N_s_reference_autonomous"],
        parent_ns["N_s_reference_autonomous"],
        rtol=0.0,
        atol=2e-12,
    )
    np.testing.assert_allclose(
        got["n_s_reference_autonomous"],
        parent_ns["n_s_reference_autonomous"],
        rtol=0.0,
        atol=2e-12,
    )
    np.testing.assert_allclose(
        got["p2_reference_autonomous"],
        parent_ns["p2_reference_autonomous"],
        rtol=0.0,
        atol=2e-12,
    )
    np.testing.assert_allclose(got["F_reference"], parent_ref["F_reference"], rtol=0.0, atol=0.0)
    np.testing.assert_allclose(got["U_reference"], parent_ref["U_reference"], rtol=0.0, atol=0.0)
    np.testing.assert_allclose(got["Pi_reference_autonomous"], parent_pi, rtol=0.0, atol=2e-12)


def test_reference_plateau_and_pressure_radial_identity() -> None:
    ext = KokunoPA10AutonomousReferenceStressToXi110()
    eta = np.asarray([-0.55, -0.1, 0.2, 0.6])
    x1 = np.asarray([100.5, 102.0, 106.0, 109.5])
    x2 = np.asarray([101.0, 104.0, 108.0, 110.0])
    v1 = ext.values(x1, eta)
    v2 = ext.values(x2, eta)
    np.testing.assert_allclose(v1["F_reference"], v2["F_reference"], rtol=0.0, atol=0.0)
    np.testing.assert_allclose(v1["U_reference"], v2["U_reference"], rtol=0.0, atol=0.0)

    state = ext.reference_state(x1, eta)
    deriv = ext.radial_derivatives(x1, eta)
    np.testing.assert_allclose(state["ell_reference"], 1.0, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(
        deriv["Pi_reference_autonomous_X"],
        v1["F_reference"] ** 2,
        rtol=1e-14,
        atol=1e-14,
    )
    np.testing.assert_allclose(deriv["F_reference_X"], 0.0, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(deriv["U_reference_X"], 0.0, rtol=0.0, atol=0.0)


def test_reference_stress_odes_against_independent_centered_x_difference() -> None:
    ext = KokunoPA10AutonomousReferenceStressToXi110()
    X = np.asarray([102.5, 105.0, 108.0])
    eta = np.asarray([-0.45, 0.0, 0.4])
    h = 2.0e-3
    plus = ext.values(X + h, eta)
    minus = ext.values(X - h, eta)
    d = ext.radial_derivatives(X, eta)

    p1_fd = (plus["p1_reference"] - minus["p1_reference"]) / (2.0 * h)
    Ns_fd = (
        plus["N_s_reference_autonomous"] - minus["N_s_reference_autonomous"]
    ) / (2.0 * h)
    ns_fd = (
        plus["n_s_reference_autonomous"] - minus["n_s_reference_autonomous"]
    ) / (2.0 * h)
    np.testing.assert_allclose(d["p1_reference_X"], p1_fd, rtol=4e-3, atol=2e-6)
    np.testing.assert_allclose(
        d["N_s_reference_autonomous_X"], Ns_fd, rtol=4e-3, atol=2e-6
    )
    np.testing.assert_allclose(
        d["n_s_reference_autonomous_X"], ns_fd, rtol=4e-3, atol=2e-6
    )


def test_quadrature_replay_and_nontrivial_xi_handoff() -> None:
    ext = KokunoPA10AutonomousReferenceStressToXi110()
    eta = np.asarray([-0.65, -0.3, 0.0, 0.25, 0.6])
    X = np.full_like(eta, 110.0)
    base = ext.values(X, eta, order=STRESS_QUADRATURE_ORDER)
    replay = ext.values(X, eta, order=96)
    for key in (
        "p1_reference",
        "N_s_reference_autonomous",
        "n_s_reference_autonomous",
        "p2_reference_autonomous",
    ):
        np.testing.assert_allclose(base[key], replay[key], rtol=2e-3, atol=2e-7)
        assert np.all(np.isfinite(base[key]))
    assert float(np.linalg.norm(base["p1_reference"])) > 1e-8
    assert float(np.linalg.norm(base["n_s_reference_autonomous"])) > 1e-10
    assert float(np.linalg.norm(base["F_reference"])) > 1e-8


def test_vectorization_serialization_domain_and_truth_boundary(tmp_path) -> None:
    ext = KokunoPA10AutonomousReferenceStressToXi110()
    X = np.asarray([[100.0], [105.0], [110.0]])
    eta = np.asarray([[-0.4, 0.0, 0.45]])
    got = ext.values(X, eta)
    assert got["p1_reference"].shape == (3, 3)
    assert got["n_s_reference_autonomous"].shape == (3, 3)

    path = tmp_path / "reference_stress_xi110.json"
    payload = ext.save_configuration(path)
    loaded = KokunoPA10AutonomousReferenceStressToXi110.load_configuration(path)
    assert loaded.configuration() == payload
    assert loaded.semantic_sha256 == ext.semantic_sha256

    bad = copy.deepcopy(payload)
    bad["source_X_i"] = 109.999
    with pytest.raises(ValueError):
        KokunoPA10AutonomousReferenceStressToXi110.from_configuration(bad)
    with pytest.raises(ValueError):
        ext.values(99.9, 0.0)
    with pytest.raises(ValueError):
        ext.values(110.1, 0.0)
    with pytest.raises(ValueError):
        ext.values(105.0, 1.0)

    truth = ext.truth_boundary
    assert truth["autonomous_pressure_reference_stress_pair_to_Xi_executable"]
    assert truth["candidate_side_reference_Xi_handoff_executable"]
    assert not truth["source_prepared_full_reference_stress_pair_materialized"]
    assert not truth["actual_final_interpolation_100_to_Xi_materialized"]
    assert not truth["outer_global_leading_velocity_materialized"]
    assert not truth["heldout_ns_residual_assessed"]
    assert not truth["pde_validated"]
    assert not truth["paper_exact"]
    assert not truth["openai_field_identified"]
