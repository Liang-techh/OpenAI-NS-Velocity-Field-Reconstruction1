from __future__ import annotations

import json

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_pa10_autonomous_axis_pressure_seed_x100 import (
    FROZEN_P_STAR,
    KokunoPA10AutonomousAxisPressureSeedToX100,
)


def _fd4_eta(
    obj: KokunoPA10AutonomousAxisPressureSeedToX100,
    X: float,
    eta: float,
    h: float,
) -> float:
    fm2 = obj.pressure(X, eta - 2.0 * h)
    fm1 = obj.pressure(X, eta - h)
    fp1 = obj.pressure(X, eta + h)
    fp2 = obj.pressure(X, eta + 2.0 * h)
    return float((fm2 - 8.0 * fm1 + 8.0 * fp1 - fp2) / (12.0 * h))


def test_autonomous_axis_seed_saturates_displayed_envelope_and_sign_condition() -> None:
    obj = KokunoPA10AutonomousAxisPressureSeedToX100()
    eta = np.asarray([-0.8, -0.35, 0.0, 0.35, 0.8])
    check = obj.admissibility(eta)

    np.testing.assert_allclose(check["Pi0_seed"], check["source_envelope"], rtol=0.0, atol=0.0)
    np.testing.assert_allclose(check["envelope_margin"], 0.0, rtol=0.0, atol=0.0)
    assert check["eta_times_Pi0_seed_eta"][2] == 0.0
    assert np.all(check["eta_times_Pi0_seed_eta"][[0, 1, 3, 4]] > 0.0)

    np.testing.assert_allclose(obj.axis_pressure(-eta), obj.axis_pressure(eta), rtol=0.0, atol=0.0)
    np.testing.assert_allclose(
        obj.axis_pressure_eta(-eta), -obj.axis_pressure_eta(eta), rtol=0.0, atol=0.0
    )


def test_axis_seed_analytic_eta_derivative_against_independent_fd4() -> None:
    obj = KokunoPA10AutonomousAxisPressureSeedToX100()
    h = 2.0e-5
    for eta in (-0.7, -0.2, 0.0, 0.3, 0.7):
        fm2 = obj.axis_pressure(eta - 2.0 * h)
        fm1 = obj.axis_pressure(eta - h)
        fp1 = obj.axis_pressure(eta + h)
        fp2 = obj.axis_pressure(eta + 2.0 * h)
        fd = float((fm2 - 8.0 * fm1 + 8.0 * fp1 - fp2) / (12.0 * h))
        exact = float(obj.axis_pressure_eta(eta))
        assert abs(fd - exact) <= 2.0e-10 * max(1.0, abs(exact))


def test_autonomous_reference_pressure_composes_parent_increment_exactly() -> None:
    obj = KokunoPA10AutonomousAxisPressureSeedToX100()
    X = np.asarray([0.0, 0.2, 2.0, 20.0, 100.0])
    eta = np.asarray([-0.6, -0.25, 0.0, 0.3, 0.65])

    cp = obj.pressure_increment.pressure_increment(X, eta)
    pi0 = obj.axis_pressure(eta)
    total = obj.pressure(X, eta)
    np.testing.assert_allclose(total - pi0, cp, rtol=0.0, atol=2.0e-16)
    np.testing.assert_allclose(
        obj.radial_derivative(X, eta),
        obj.pressure_increment.radial_derivative(X, eta),
        rtol=0.0,
        atol=0.0,
    )


def test_total_eta_derivative_replays_with_independent_fd4() -> None:
    obj = KokunoPA10AutonomousAxisPressureSeedToX100()
    probes = [(0.5, -0.5), (5.0, -0.1), (40.0, 0.25), (100.0, 0.55)]
    for X, eta in probes:
        production = float(obj.eta_derivative(X, eta))
        independent = _fd4_eta(obj, X, eta, 8.0e-5)
        assert abs(production - independent) <= 4.0e-5 * max(1.0, abs(production))


def test_vectorization_handoff_and_truth_boundary() -> None:
    obj = KokunoPA10AutonomousAxisPressureSeedToX100()
    eta = np.asarray([-0.6, 0.0, 0.6])
    handoff = obj.handoff_at_X100(eta)
    assert handoff["Pi0_autonomous_seed"].shape == eta.shape
    assert handoff["Pi_reference_autonomous"].shape == eta.shape
    assert handoff["Pi_reference_autonomous_X"].shape == eta.shape
    assert handoff["Pi_reference_autonomous_eta"].shape == eta.shape
    assert np.all(np.isfinite(handoff["Pi_reference_autonomous"]))
    assert np.all(handoff["Pi_reference_autonomous_X"] > 0.0)

    truth = obj.truth_boundary
    assert truth["repository_autonomous_axis_pressure_seed_materialized"] is True
    assert truth["autonomous_reference_pressure_to_X100_executable"] is True
    assert truth["source_prepared_appendixA_Pi0_materialized"] is False
    assert truth["absolute_axis_pressure_Pi0_materialized"] is False
    assert truth["source_large_P_star_regime_admitted"] is False
    assert truth["outer_tail_pressure_compatibility_verified"] is False
    assert truth["matched_global_pressure_materialized"] is False
    assert truth["reference_nsr_materialized"] is False
    assert truth["heldout_ns_residual_assessed"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False


def test_serialization_hash_and_autonomous_P_star_are_frozen(tmp_path) -> None:
    obj = KokunoPA10AutonomousAxisPressureSeedToX100()
    assert obj.P_star == FROZEN_P_STAR == 1.0
    cfg = obj.configuration()
    assert cfg["P_star"] == 1.0
    assert len(obj.semantic_sha256) == 64

    path = tmp_path / "autonomous_axis_pressure.json"
    obj.save_configuration(path)
    loaded = KokunoPA10AutonomousAxisPressureSeedToX100.load_configuration(path)
    assert loaded.configuration() == cfg
    assert loaded.semantic_sha256 == obj.semantic_sha256

    bad = json.loads(path.read_text())
    bad["P_star"] = 2.0
    with pytest.raises(ValueError):
        KokunoPA10AutonomousAxisPressureSeedToX100.from_configuration(bad)
    with pytest.raises(ValueError):
        KokunoPA10AutonomousAxisPressureSeedToX100(P_star=2.0)


def test_domain_and_parent_eta_stencil_fail_closed() -> None:
    obj = KokunoPA10AutonomousAxisPressureSeedToX100()
    with pytest.raises(ValueError):
        obj.axis_pressure(1.0001)
    with pytest.raises(ValueError):
        obj.axis_pressure(-1.0001)
    with pytest.raises(ValueError):
        obj.pressure(-1.0e-6, 0.0)
    with pytest.raises(ValueError):
        obj.pressure(obj.X_b_ref + 1.0e-5, 0.0)

    eta_lo, eta_hi = obj.eta_interval
    with pytest.raises(ValueError):
        obj.eta_derivative(1.0, eta_lo + 1.0e-5)
    with pytest.raises(ValueError):
        obj.eta_derivative(1.0, eta_hi - 1.0e-5)
