from __future__ import annotations

import json

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_pa10_reference_pressure_increment_x100 import (
    ETA_FD_STEP,
    RADIAL_QUADRATURE_ORDER,
    KokunoPA10ReferencePressureIncrementToX100,
)


def _fd4_x(obj: KokunoPA10ReferencePressureIncrementToX100, X: float, eta: float, h: float) -> float:
    fm2 = obj.pressure_increment(X - 2.0 * h, eta)
    fm1 = obj.pressure_increment(X - h, eta)
    fp1 = obj.pressure_increment(X + h, eta)
    fp2 = obj.pressure_increment(X + 2.0 * h, eta)
    return float((fm2 - 8.0 * fm1 + 8.0 * fp1 - fp2) / (12.0 * h))


def test_source_pressure_increment_radial_identity_against_independent_fd4() -> None:
    obj = KokunoPA10ReferencePressureIncrementToX100()
    probes = [(0.02, -0.4), (0.5, 0.0), (2.0, 0.35), (25.0, -0.2), (80.0, 0.6)]
    for X, eta in probes:
        h = min(2.0e-4, X / 8.0, (obj.X_b_ref - X) / 8.0)
        fd = _fd4_x(obj, X, eta, h)
        exact = float(obj.radial_derivative(X, eta))
        scale = max(1.0, abs(exact))
        assert abs(fd - exact) / scale < 3.0e-6


def test_pressure_increment_matches_E_squared_over_2X() -> None:
    obj = KokunoPA10ReferencePressureIncrementToX100()
    X = np.asarray([0.03, 0.4, 3.0, 20.0, 100.0])
    eta = np.asarray([-0.7, -0.2, 0.0, 0.3, 0.7])
    ref = obj.reference.values(X, eta)
    lhs = obj.radial_derivative(X, eta)
    rhs = ref["E_reference"] ** 2 / (2.0 * X)
    np.testing.assert_allclose(lhs, rhs, rtol=2.0e-13, atol=2.0e-15)


def test_zero_axis_increment_and_positive_monotone_radial_increment() -> None:
    obj = KokunoPA10ReferencePressureIncrementToX100()
    for eta in (-0.7, 0.0, 0.7):
        assert float(obj.pressure_increment(0.0, eta)) == 0.0
        values = obj.pressure_increment(np.asarray([0.0, 0.1, 1.0, 10.0, 100.0]), eta)
        assert np.all(np.isfinite(values))
        assert np.all(np.diff(values) > 0.0)


def test_eta_derivative_is_finite_and_replays_with_independent_step() -> None:
    obj = KokunoPA10ReferencePressureIncrementToX100()
    X = np.asarray([0.5, 5.0, 50.0, 100.0])
    eta = np.asarray([-0.5, -0.1, 0.25, 0.55])
    production = obj.eta_derivative(X, eta)

    h = 4.0 * ETA_FD_STEP
    independent = (
        obj.pressure_increment(X, eta - 2.0 * h)
        - 8.0 * obj.pressure_increment(X, eta - h)
        + 8.0 * obj.pressure_increment(X, eta + h)
        - obj.pressure_increment(X, eta + 2.0 * h)
    ) / (12.0 * h)
    assert np.all(np.isfinite(production))
    np.testing.assert_allclose(production, independent, rtol=2.0e-5, atol=2.0e-8)


def test_vectorization_handoff_and_truth_boundary() -> None:
    obj = KokunoPA10ReferencePressureIncrementToX100()
    eta = np.asarray([-0.6, 0.0, 0.6])
    handoff = obj.handoff_at_X100(eta)
    assert handoff["C_p_reference"].shape == eta.shape
    assert handoff["C_p_reference_X"].shape == eta.shape
    assert handoff["C_p_reference_eta"].shape == eta.shape
    assert np.all(handoff["C_p_reference"] > 0.0)

    truth = obj.truth_boundary
    assert truth["reference_pressure_increment_to_X100_executable"] is True
    assert truth["absolute_axis_pressure_Pi0_materialized"] is False
    assert truth["absolute_reference_pressure_materialized"] is False
    assert truth["reference_nsr_materialized"] is False
    assert truth["matched_global_pressure_materialized"] is False
    assert truth["heldout_ns_residual_assessed"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False


def test_serialization_semantic_hash_and_frozen_numerics(tmp_path) -> None:
    obj = KokunoPA10ReferencePressureIncrementToX100()
    cfg = obj.configuration()
    assert cfg["radial_quadrature_order"] == RADIAL_QUADRATURE_ORDER
    assert cfg["eta_fd_step"] == ETA_FD_STEP
    assert len(obj.semantic_sha256) == 64

    path = tmp_path / "pressure_increment.json"
    obj.save_configuration(path)
    loaded = KokunoPA10ReferencePressureIncrementToX100.load_configuration(path)
    assert loaded.configuration() == cfg
    assert loaded.semantic_sha256 == obj.semantic_sha256

    bad = json.loads(path.read_text())
    bad["eta_fd_step"] *= 2.0
    with pytest.raises(ValueError):
        KokunoPA10ReferencePressureIncrementToX100.from_configuration(bad)


def test_domain_and_eta_stencil_fail_closed() -> None:
    obj = KokunoPA10ReferencePressureIncrementToX100()
    with pytest.raises(ValueError):
        obj.pressure_increment(-1.0e-6, 0.0)
    with pytest.raises(ValueError):
        obj.pressure_increment(obj.X_b_ref + 1.0e-5, 0.0)

    eta_lo, eta_hi = obj.eta_interval
    with pytest.raises(ValueError):
        obj.eta_derivative(1.0, eta_lo + ETA_FD_STEP)
    with pytest.raises(ValueError):
        obj.eta_derivative(1.0, eta_hi - ETA_FD_STEP)
