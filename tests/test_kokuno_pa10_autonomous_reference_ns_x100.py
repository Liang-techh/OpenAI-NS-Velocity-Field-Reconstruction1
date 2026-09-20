from __future__ import annotations

import json

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_pa10_autonomous_reference_ns_x100 import (
    ETA_FD_STEP,
    NS_QUADRATURE_ORDER,
    KokunoPA10AutonomousReferenceNsToX100,
)


def test_stress_free_X0_handoff_replays_public_ns_identity() -> None:
    obj = KokunoPA10AutonomousReferenceNsToX100()
    eta = np.asarray([-0.7, -0.2, 0.0, 0.3, 0.7])
    initial = obj.initial_values(eta)
    source = obj.reference.collar.source_normalized
    X = np.full_like(eta, obj.X_0)
    U_X = source.derivatives(X, eta)["U_0_X"]
    L = 1.0 - 2.0 * obj.h * eta * eta

    np.testing.assert_allclose(initial["n_s_initial"], -2.0 * U_X, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(initial["N_s_initial"], L * initial["n_s_initial"], rtol=0.0, atol=0.0)

    propagated = obj.values(X, eta)
    np.testing.assert_allclose(
        propagated["n_s_reference_autonomous"], initial["n_s_initial"], rtol=2.0e-15, atol=2.0e-15
    )
    np.testing.assert_allclose(
        propagated["N_s_reference_autonomous"], initial["N_s_initial"], rtol=2.0e-15, atol=2.0e-15
    )


def test_axial_source_replays_displayed_formula_term_by_term() -> None:
    obj = KokunoPA10AutonomousReferenceNsToX100()
    X = np.asarray([0.7, 3.0, 17.0, 75.0])
    eta = np.asarray([-0.5, -0.1, 0.25, 0.6])
    terms = obj.source_terms(X, eta)

    U = terms["U_reference"]
    direct = (
        -terms["W_reference"] * X * terms["U_reference_X"]
        - obj.A * (1.0 - 2.0 * eta * U) * U
        - terms["H_c_reference"] * terms["U_reference_eta"]
        - terms["d"] * terms["Pi_reference_autonomous_eta"]
        + 4.0 * obj.A * eta * terms["Pi_reference_autonomous"]
        + 2.0 * eta * X * terms["Pi_reference_autonomous_X"]
    )
    np.testing.assert_allclose(
        terms["S_n_reference_autonomous"], direct, rtol=2.0e-15, atol=2.0e-15
    )
    assert np.all(np.isfinite(direct))
    assert float(np.sqrt(np.mean(direct * direct))) > 1.0e-10


def test_Ns_and_ns_radial_derivatives_replay_source_ode_with_centered_difference() -> None:
    obj = KokunoPA10AutonomousReferenceNsToX100()
    probes = [(1.0, -0.45), (8.0, 0.15), (55.0, 0.5)]
    for X, eta in probes:
        h = 1.0e-3 * max(1.0, X)
        minus = obj.values(X - h, eta)
        plus = obj.values(X + h, eta)
        fd_N = float(
            (plus["N_s_reference_autonomous"] - minus["N_s_reference_autonomous"])
            / (2.0 * h)
        )
        fd_n = float(
            (plus["n_s_reference_autonomous"] - minus["n_s_reference_autonomous"])
            / (2.0 * h)
        )
        production = obj.radial_derivatives(X, eta)
        exact_N = float(production["N_s_reference_autonomous_X"])
        exact_n = float(production["n_s_reference_autonomous_X"])
        assert abs(fd_N - exact_N) <= 2.0e-3 * max(1.0, abs(exact_N))
        assert abs(fd_n - exact_n) <= 2.0e-3 * max(1.0, abs(exact_n))


def test_p2_identity_vectorization_and_X100_handoff_are_nontrivial() -> None:
    obj = KokunoPA10AutonomousReferenceNsToX100()
    eta = np.asarray([-0.65, -0.2, 0.0, 0.25, 0.65])
    handoff = obj.handoff_at_X100(eta)
    E = obj.reference.values(np.full_like(eta, obj.X_b_ref), eta)["E_reference"]
    expected_p2 = obj.X_b_ref * handoff["n_s_reference_autonomous"] / E

    np.testing.assert_allclose(
        handoff["p2_reference_autonomous"], expected_p2, rtol=2.0e-15, atol=2.0e-15
    )
    for key in (
        "p1_reference",
        "n_s_reference_autonomous",
        "p2_reference_autonomous",
        "S_n_reference_autonomous",
        "n_s_reference_autonomous_X",
    ):
        assert handoff[key].shape == eta.shape
        assert np.all(np.isfinite(handoff[key]))
    assert float(np.sqrt(np.mean(handoff["n_s_reference_autonomous"] ** 2))) > 1.0e-10


def test_fixed_quadrature_and_eta_derivative_realizations_are_frozen() -> None:
    obj = KokunoPA10AutonomousReferenceNsToX100()
    cfg = obj.configuration()
    assert cfg["Ns_quadrature_order"] == NS_QUADRATURE_ORDER == 64
    assert cfg["eta_fd_step"] == ETA_FD_STEP == 2.0e-5

    # A higher-order radial replay is a numerical stability diagnostic, not an
    # independent PDE validation.  Keep the gate deliberately engineering-only.
    eta = np.asarray([-0.45, 0.0, 0.45])
    X = np.asarray([7.0, 40.0, 100.0])
    default = obj.values(X, eta)["n_s_reference_autonomous"]
    replay = obj.values(X, eta, order=96)["n_s_reference_autonomous"]
    scale = np.maximum(1.0, np.abs(replay))
    assert float(np.max(np.abs(default - replay) / scale)) <= 2.0e-3


def test_serialization_hash_and_truth_boundary_fail_closed(tmp_path) -> None:
    obj = KokunoPA10AutonomousReferenceNsToX100()
    path = tmp_path / "autonomous_reference_ns.json"
    cfg = obj.save_configuration(path)
    loaded = KokunoPA10AutonomousReferenceNsToX100.load_configuration(path)
    assert loaded.configuration() == cfg
    assert loaded.semantic_sha256 == obj.semantic_sha256
    assert len(obj.semantic_sha256) == 64

    bad = json.loads(path.read_text())
    bad["Ns_quadrature_order"] = 96
    with pytest.raises(ValueError):
        KokunoPA10AutonomousReferenceNsToX100.from_configuration(bad)

    truth = obj.truth_boundary
    assert truth["public_Sn_formula_executable"] is True
    assert truth["public_Ns_ode_executable_to_X100"] is True
    assert truth["autonomous_pressure_reference_ns_to_X100_executable"] is True
    assert truth["autonomous_pressure_reference_stress_pair_p1r_nsr_executable"] is True
    assert truth["source_prepared_appendixA_Pi0_materialized"] is False
    assert truth["source_prepared_reference_nsr_materialized"] is False
    assert truth["source_prepared_full_reference_stress_pair_materialized"] is False
    assert truth["selected_kappa0_global_source_smallness_admitted"] is False
    assert truth["kappa0_continuation_to_X100_materialized"] is False
    assert truth["outer_global_leading_velocity_materialized"] is False
    assert truth["matched_global_pressure_materialized"] is False
    assert truth["heldout_ns_residual_assessed"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False


def test_domain_guards_fail_closed() -> None:
    obj = KokunoPA10AutonomousReferenceNsToX100()
    with pytest.raises(ValueError):
        obj.values(obj.X_0 - 1.0e-8, 0.0)
    with pytest.raises(ValueError):
        obj.values(obj.X_b_ref + 1.0e-5, 0.0)
    with pytest.raises(ValueError):
        obj.values(1.0, -1.0 + ETA_FD_STEP)
    with pytest.raises(ValueError):
        obj.values(1.0, 1.0 - ETA_FD_STEP)
