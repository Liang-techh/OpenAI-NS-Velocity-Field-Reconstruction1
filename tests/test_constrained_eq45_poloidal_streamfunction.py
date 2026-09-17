import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_poloidal_streamfunction import (
    eq45_poloidal_divergence_numerator,
    eq45_streamfunction_poloidal_jet,
)


def test_streamfunction_coupling_cancels_profile_divergence_for_random_jets():
    rng = np.random.default_rng(20260916)
    shape = (257,)
    X = rng.uniform(0.0, 4.0, size=shape)
    eta = rng.uniform(-0.95, 0.95, size=shape)
    h = rng.uniform(1e-4, 0.2, size=shape)
    phi, phi_X, phi_eta, phi_XX, phi_Xeta = rng.normal(size=(5,) + shape)

    jet = eq45_streamfunction_poloidal_jet(
        X, eta, h, phi, phi_X, phi_eta, phi_XX, phi_Xeta
    )
    numerator = eq45_poloidal_divergence_numerator(
        X, eta, h, jet.v0, jet.v0_X, jet.U, jet.U_X, jet.U_eta
    )

    assert np.max(np.abs(numerator)) < 5e-13


def test_independent_v0_U_profiles_are_not_generically_incompressible():
    # A constant axial profile U=1 with v0=0 is a valid pointwise Eq. (4.5)
    # velocity channel, but it violates the coupled poloidal divergence identity
    # away from eta=0. This is the obstruction that the scalar streamfunction
    # parameterization removes by construction.
    X = np.array([0.0, 0.5, 2.0])
    eta = np.array([0.2, -0.4, 0.7])
    h = 0.005
    zeros = np.zeros_like(X)
    ones = np.ones_like(X)

    numerator = eq45_poloidal_divergence_numerator(
        X, eta, h, zeros, zeros, ones, zeros, zeros
    )
    expected = -(1.0 + 2.0 * h) * eta

    assert np.allclose(numerator, expected, rtol=0.0, atol=2e-15)
    assert np.min(np.abs(numerator)) > 0.1


def test_axis_regular_map_and_fail_closed_source_domain():
    X = np.array([0.0, 0.0, 1e-12, 0.25])
    eta = np.array([-0.8, 0.0, 0.8, 0.3])
    h = 0.005
    phi = 1.0 + 0.2 * X + 0.1 * eta
    phi_X = np.full_like(X, 0.2)
    phi_eta = np.full_like(X, 0.1)
    phi_XX = np.zeros_like(X)
    phi_Xeta = np.zeros_like(X)

    jet = eq45_streamfunction_poloidal_jet(
        X, eta, h, phi, phi_X, phi_eta, phi_XX, phi_Xeta
    )
    assert np.all(np.isfinite(jet.v0))
    assert np.all(np.isfinite(jet.U))
    assert np.all(np.isfinite(jet.v0_X))

    with pytest.raises(ValueError):
        eq45_streamfunction_poloidal_jet(
            -1e-3, 0.0, h, 1.0, 0.0, 0.0, 0.0, 0.0
        )
    with pytest.raises(ValueError):
        eq45_streamfunction_poloidal_jet(
            0.2, 1.0, h, 1.0, 0.0, 0.0, 0.0, 0.0
        )
    with pytest.raises(ValueError):
        eq45_streamfunction_poloidal_jet(
            0.2, 0.0, 0.5, 1.0, 0.0, 0.0, 0.0, 0.0
        )
