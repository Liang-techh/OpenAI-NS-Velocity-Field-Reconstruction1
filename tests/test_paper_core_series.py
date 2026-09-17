import numpy as np
import pytest

from openai_ns_reconstruction.paper_core_reference import PaperCoreReference
from openai_ns_reconstruction.paper_core_series import (
    ChebyshevEtaGrid,
    PaperCoreSeries,
)


def test_chebyshev_eta_differentiation_and_interpolation_are_polynomial_exact():
    grid = ChebyshevEtaGrid(33)
    eta = grid.eta
    values = np.vstack((eta**2, eta**3, 2.0 - eta + 4.0 * eta**4))
    expected = np.vstack((2.0 * eta, 3.0 * eta**2, -1.0 + 16.0 * eta**3))
    np.testing.assert_allclose(grid.differentiate(values), expected, atol=2e-12)

    points = np.array([-0.73, -0.1, 0.0, 0.42, 0.91])
    interpolated = grid.interpolate(values, points)
    expected_points = np.vstack(
        (points**2, points**3, 2.0 - points + 4.0 * points**4)
    )
    np.testing.assert_allclose(interpolated, expected_points, atol=2e-12)


def test_first_radial_rows_and_pressure_datum_are_preserved():
    reference = PaperCoreReference()
    series = PaperCoreSeries(reference, maxdegree=2, eta_nodes=65)
    eta = series.grid.eta

    np.testing.assert_allclose(series.phi[0], 1.0, atol=0.0)
    np.testing.assert_allclose(series.u[0], 4.0 * eta + reference.j, atol=2e-14)
    np.testing.assert_allclose(
        series.pi[0],
        -reference.pressure_scale**2 / (1.0 + eta**2) ** 2,
        atol=2e-14,
    )
    np.testing.assert_allclose(series.pi[1], series.g_values**2, rtol=2e-13)

    points = np.array([-0.61, -0.13, 0.0, 0.27, 0.84])
    np.testing.assert_allclose(
        series.F(0.0, points),
        np.array([reference.amplitude(float(e)) for e in points]),
        rtol=3e-13,
        atol=2e-15,
    )
    np.testing.assert_allclose(series.U(0.0, points), 4.0 * points + reference.j)
    np.testing.assert_allclose(
        series.Pi(0.0, points),
        -reference.pressure_scale**2 / (1.0 + points**2) ** 2,
        atol=2e-13,
    )
    np.testing.assert_allclose(series.dU_deta(0.0, points), 4.0, atol=2e-12)
    np.testing.assert_allclose(
        series.Pi_radial_derivative(0.07, points), series.F(0.07, points) ** 2
    )


def test_finite_profile_adapter_is_vectorized_and_bounded():
    series = PaperCoreSeries(maxdegree=2, eta_nodes=65)
    X = np.array([0.0, 0.01, 0.07])
    eta = np.array([-0.4, 0.0, 0.5])
    profile = series.profile

    np.testing.assert_allclose(profile.U(X, eta), series.U(X, eta))
    np.testing.assert_allclose(profile.F(X, eta), series.F(X, eta))
    np.testing.assert_allclose(profile.Pi(X, eta), series.Pi(X, eta))
    np.testing.assert_allclose(
        [profile.radial_average_U(float(x), float(e)) for x, e in zip(X, eta)],
        series.radial_average_U(X, eta),
    )
    velocity = series.velocity(0.01, -0.02, 0.1, 0.5)
    assert velocity.shape == (3,)
    assert np.all(np.isfinite(velocity))
    assert profile.paper_exact is False
    assert "finite" in profile.name


def test_finite_degree_and_eta_grid_caps_fail_closed():
    with pytest.raises(ValueError):
        PaperCoreSeries(maxdegree=17)
    with pytest.raises(ValueError):
        PaperCoreSeries(eta_nodes=515)
    with pytest.raises(ValueError):
        PaperCoreSeries().F(0.5, 0.0)
