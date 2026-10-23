from __future__ import annotations

import json
import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_reference_continuation import (
    KokunoReferenceContinuationCandidate,
    source_smooth_step,
)


def test_source_step_and_transition_bound() -> None:
    points = np.array([-1.0, 0.0, 0.25, 0.5, 0.75, 1.0, 2.0])
    values = source_smooth_step(points)
    assert values[0] == 0.0
    assert values[1] == 0.0
    assert values[-2] == 1.0
    assert values[-1] == 1.0
    assert values[2] < values[3] < values[4]
    assert values[3] == pytest.approx(0.5, abs=1e-15)

    maximum = 0.5 * math.log(4.1 / 4.0)
    with pytest.raises(ValueError):
        KokunoReferenceContinuationCandidate(log_transition_width=maximum)
    candidate = KokunoReferenceContinuationCandidate(log_transition_width=0.005)
    assert candidate.X2 < candidate.max_source_transition_X


def test_reference_continuation_replays_core_and_freezes_after_transition() -> None:
    candidate = KokunoReferenceContinuationCandidate()
    eta = np.array([-0.6, -0.1, 0.0, 0.35, 0.7])
    X_core = 0.8 * candidate.X1

    np.testing.assert_allclose(
        candidate.F(X_core, eta), candidate.core.series.F(X_core, eta), rtol=0.0, atol=0.0
    )
    np.testing.assert_allclose(
        candidate.U(X_core, eta), candidate.core.series.U(X_core, eta), rtol=0.0, atol=0.0
    )

    F2 = candidate.F(candidate.X2, eta)
    U2 = candidate.U(candidate.X2, eta)
    np.testing.assert_allclose(candidate.F(1.5 * candidate.X2, eta), F2, rtol=2e-13, atol=2e-13)
    np.testing.assert_allclose(candidate.F(5.0 * candidate.X2, eta), F2, rtol=2e-13, atol=2e-13)
    np.testing.assert_allclose(candidate.U(1.5 * candidate.X2, eta), U2, rtol=2e-13, atol=2e-13)
    np.testing.assert_allclose(candidate.U(5.0 * candidate.X2, eta), U2, rtol=2e-13, atol=2e-13)
    np.testing.assert_allclose(candidate.F_X(2.0 * candidate.X2, eta), 0.0, atol=0.0)
    np.testing.assert_allclose(candidate.U_X(2.0 * candidate.X2, eta), 0.0, atol=0.0)


def test_transition_radial_derivatives_follow_source_prescription_and_fd() -> None:
    candidate = KokunoReferenceContinuationCandidate(quadrature_points=32)
    X = candidate.X0 * math.exp(1.5 * candidate.log_transition_width)
    eta = 0.23
    y = math.log(X / candidate.X0)
    gate = 1.0 - float(
        source_smooth_step(
            (y - candidate.log_transition_width) / candidate.log_transition_width
        )
    )

    F = float(candidate.F(X, eta))
    Fnat = float(candidate.core.series.F(X, eta))
    source_F_X = F * gate * float(candidate.core.series.F_radial_derivative(X, eta)) / Fnat
    source_U_X = gate * float(candidate.core.series.U_radial_derivative(X, eta))
    assert float(candidate.F_X(X, eta)) == pytest.approx(source_F_X, rel=2e-12, abs=2e-12)
    assert float(candidate.U_X(X, eta)) == pytest.approx(source_U_X, rel=2e-12, abs=2e-12)

    step = 2.0e-6
    fd_F = (float(candidate.F(X + step, eta)) - float(candidate.F(X - step, eta))) / (2.0 * step)
    fd_U = (float(candidate.U(X + step, eta)) - float(candidate.U(X - step, eta))) / (2.0 * step)
    assert fd_F == pytest.approx(source_F_X, rel=2e-5, abs=2e-6)
    assert fd_U == pytest.approx(source_U_X, rel=2e-5, abs=2e-6)


def test_eta_derivatives_and_pressure_identity_are_independently_resolved() -> None:
    candidate = KokunoReferenceContinuationCandidate(quadrature_points=32)
    X = candidate.X0 * math.exp(1.65 * candidate.log_transition_width)
    eta = 0.19
    d_eta = 2.0e-5
    fd_F_eta = (float(candidate.F(X, eta + d_eta)) - float(candidate.F(X, eta - d_eta))) / (2.0 * d_eta)
    fd_U_eta = (float(candidate.U(X, eta + d_eta)) - float(candidate.U(X, eta - d_eta))) / (2.0 * d_eta)
    assert float(candidate.F_eta(X, eta)) == pytest.approx(fd_F_eta, rel=2e-4, abs=2e-6)
    assert float(candidate.U_eta(X, eta)) == pytest.approx(fd_U_eta, rel=2e-4, abs=2e-6)

    outer_X = 0.75
    dX = 2.0e-5
    fd_Pi = (float(candidate.Pi(outer_X + dX, eta)) - float(candidate.Pi(outer_X - dX, eta))) / (2.0 * dX)
    expected = float(candidate.F(outer_X, eta)) ** 2
    assert fd_Pi == pytest.approx(expected, rel=2e-6, abs=2e-7)


def test_velocity_extends_beyond_old_core_without_zero_collapse() -> None:
    candidate = KokunoReferenceContinuationCandidate(quadrature_points=16)
    core = candidate.core

    # At z=0,t=.5, q=.5 and X=r^2.  This point remains well inside X1.
    point_core = np.array([[0.35, 0.0, 0.0]])
    np.testing.assert_allclose(
        candidate.at_points(point_core, 0.5), core.at_points(point_core, 0.5),
        rtol=2e-10, atol=2e-10,
    )

    # r=.8 gives X=.64, beyond the old core guard Lambda*X<=4.1.
    with pytest.raises(ValueError):
        core.at_points(np.array([[0.8, 0.0, 0.0]]), 0.5)
    velocity = candidate.at_points(np.array([[0.8, 0.0, 0.0]]), 0.5)
    assert velocity.shape == (1, 3)
    assert np.all(np.isfinite(velocity))
    assert np.linalg.norm(velocity[0]) > 1.0e-6

    axis_velocity = candidate.velocity(0.0, 0.0, 0.0, 0.5)
    assert axis_velocity[0] == pytest.approx(0.0, abs=1e-14)
    assert axis_velocity[1] == pytest.approx(0.0, abs=1e-14)
    assert abs(float(axis_velocity[2])) > 1.0e-6


def test_serialization_roundtrip_is_fail_closed(tmp_path) -> None:
    candidate = KokunoReferenceContinuationCandidate(quadrature_points=16)
    path = candidate.save_json(tmp_path / "candidate.json")
    replay = KokunoReferenceContinuationCandidate.load_json(path)
    assert replay.sha256 == candidate.sha256
    probes = np.array([[0.25, 0.1, 0.05], [0.8, 0.0, 0.0]])
    np.testing.assert_allclose(replay.at_points(probes, 0.5), candidate.at_points(probes, 0.5), rtol=0.0, atol=0.0)

    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["truth_boundary"]["global_leading_profile_reconstructed"] = True
    with pytest.raises(ValueError):
        KokunoReferenceContinuationCandidate.from_payload(payload)
