from __future__ import annotations

import json

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_leading_core_series import (
    KokunoLeadingCoreSeriesCandidate,
)


@pytest.fixture(scope="module")
def candidate() -> KokunoLeadingCoreSeriesCandidate:
    # Keep the regression lightweight while exercising the same finite-series API.
    return KokunoLeadingCoreSeriesCandidate(
        sigma=0.5, maxdegree=6, eta_nodes=65, quadrature_points=16
    )


def test_native_coordinates_and_axis_seed_are_bound(candidate):
    c = candidate
    values = c.coordinates(0.12, -0.08, 0.1, 0.5)
    q = float(values["q"])
    eta = float(values["eta"])
    X = float(values["X"])
    assert q > 0.0
    assert abs(q - 0.1**2 * q ** (2.0 * c.h) - 0.5) < 1e-11
    assert abs(eta) < 1.0
    assert c.Lambda * X < 4.1

    axis = c.profile_values(0.0, np.array([-0.4, 0.0, 0.3]))
    expected_u = 4.0 * np.array([-0.4, 0.0, 0.3]) + c.j0
    np.testing.assert_allclose(axis["U"], expected_u, rtol=0.0, atol=2e-12)
    np.testing.assert_allclose(axis["Pi_X"], axis["F"] ** 2, rtol=2e-12, atol=2e-12)

    center = c.profile_values(0.0, 0.0)
    np.testing.assert_allclose(center["F"], 0.5, rtol=0.0, atol=2e-13)
    np.testing.assert_allclose(center["U"], 0.02, rtol=0.0, atol=2e-13)
    np.testing.assert_allclose(center["Pi_X"], 0.25, rtol=0.0, atol=2e-13)
    origin_velocity = c.velocity(0.0, 0.0, 0.0, 0.5)
    np.testing.assert_allclose(
        origin_velocity,
        [0.0, 0.0, 0.02838246712400765],
        rtol=2e-13,
        atol=2e-13,
    )


def test_profile_partials_match_independent_centered_differences(candidate):
    c = candidate
    X = 0.08
    eta = 0.17
    values = c.profile_values(X, eta)
    hX = 2e-6
    he = 2e-6

    f_x_fd = (
        c.profile_values(X + hX, eta)["F"]
        - c.profile_values(X - hX, eta)["F"]
    ) / (2.0 * hX)
    u_x_fd = (
        c.profile_values(X + hX, eta)["U"]
        - c.profile_values(X - hX, eta)["U"]
    ) / (2.0 * hX)
    f_e_fd = (
        c.profile_values(X, eta + he)["F"]
        - c.profile_values(X, eta - he)["F"]
    ) / (2.0 * he)
    u_e_fd = (
        c.profile_values(X, eta + he)["U"]
        - c.profile_values(X, eta - he)["U"]
    ) / (2.0 * he)

    np.testing.assert_allclose(values["F_X"], f_x_fd, rtol=2e-5, atol=2e-7)
    np.testing.assert_allclose(values["U_X"], u_x_fd, rtol=2e-5, atol=2e-7)
    np.testing.assert_allclose(values["F_eta"], f_e_fd, rtol=2e-4, atol=2e-6)
    np.testing.assert_allclose(values["U_eta"], u_e_fd, rtol=2e-4, atol=2e-6)


def test_vectorized_velocity_at_points_and_grid_replay(candidate):
    c = candidate
    points = np.array(
        [[0.0, 0.0, 0.0], [0.08, 0.03, 0.04], [-0.1, 0.05, -0.06]],
        dtype=float,
    )
    times = np.array([0.4, 0.5, 0.6])
    velocity = c.at_points(points, times)
    assert velocity.shape == (3, 3)
    assert np.all(np.isfinite(velocity))
    assert np.linalg.norm(velocity) > 1e-4

    for i in range(3):
        scalar = c.series.velocity(
            *points[i], float(times[i]), quadrature_points=c.quadrature_points
        )
        np.testing.assert_allclose(velocity[i], scalar, rtol=0.0, atol=2e-13)

    grid = c.grid([-0.06, 0.06], [-0.04, 0.04], [-0.05, 0.05], [0.45, 0.55])
    assert grid.shape == (2, 2, 2, 2, 3)
    assert np.all(np.isfinite(grid))


def test_axis_regular_and_nontrivial(candidate):
    c = candidate
    axis = c.velocity(0.0, 0.0, 0.08, 0.5)
    np.testing.assert_allclose(axis[:2], 0.0, rtol=0.0, atol=0.0)
    assert abs(float(axis[2])) > 1e-4

    v1 = c.velocity(1e-5, 0.0, 0.08, 0.5)
    v2 = c.velocity(2e-5, 0.0, 0.08, 0.5)
    transverse1 = np.linalg.norm(v1[:2])
    transverse2 = np.linalg.norm(v2[:2])
    assert transverse1 > 0.0
    assert 1.8 < transverse2 / transverse1 < 2.2


def test_core_domain_fails_closed(candidate):
    c = candidate
    with pytest.raises(ValueError, match=r"Lambda\*X <= 4.1"):
        c.velocity(1.0, 0.0, 0.0, 0.5)
    with pytest.raises(ValueError, match="0 <= t < 1"):
        c.velocity(0.0, 0.0, 0.0, 1.0)
    with pytest.raises(ValueError, match="shape"):
        c.at_points(np.zeros((4, 2)), 0.5)


def test_json_sha_roundtrip_and_truth_boundary(candidate, tmp_path):
    c = candidate
    path = c.save_json(tmp_path / "candidate.json")
    loaded = KokunoLeadingCoreSeriesCandidate.load_json(path)
    assert loaded.sha256 == c.sha256
    probe = np.array([[0.04, -0.03, 0.02], [0.07, 0.02, -0.04]])
    np.testing.assert_allclose(
        loaded.at_points(probe, 0.5), c.at_points(probe, 0.5), rtol=0.0, atol=0.0
    )

    payload = json.loads(path.read_text())
    assert payload["truth_boundary"]["finite_nonlinear_core_series_executable"] is True
    assert payload["truth_boundary"]["global_leading_profile_reconstructed"] is False
    assert payload["truth_boundary"]["pde_validated"] is False
    assert payload["truth_boundary"]["paper_exact"] is False

    payload["truth_boundary"]["pde_validated"] = True
    with pytest.raises(ValueError, match="truth-boundary"):
        KokunoLeadingCoreSeriesCandidate.from_payload(payload)
