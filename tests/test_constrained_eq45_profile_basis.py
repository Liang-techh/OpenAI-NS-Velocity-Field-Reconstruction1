import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_profile_basis import (
    Eq45CompactProfileBasis,
)


def test_seed_is_bounded_nonzero_vectorized_and_axis_regular():
    basis = Eq45CompactProfileBasis.seed()
    x = np.array([0.0, 0.1, 0.7, 1.5])
    eta = np.array([0.0, 0.25, -0.4, 0.6])
    jets = basis.evaluate(x, eta)

    assert jets.stacked().shape == (4, 4)
    assert np.all(np.isfinite(jets.stacked()))
    assert np.max(np.abs(jets.stacked())) > 0.0

    axis = basis.evaluate(np.zeros(9), np.linspace(-0.8, 0.8, 9))
    assert np.all(np.isfinite(axis.stacked()))
    assert basis.parameter_count == 12
    assert all(abs(v) <= basis.coefficient_limit for v in basis.phi_coefficients)
    assert all(abs(v) <= basis.coefficient_limit for v in basis.swirl_coefficients)


def test_analytic_phi_derivatives_match_independent_centered_difference():
    basis = Eq45CompactProfileBasis.seed()
    x = np.array([0.2, 0.7, 1.4, 2.2])
    eta = np.array([-0.55, -0.2, 0.25, 0.6])
    h = 1e-6

    analytic = basis.evaluate(x, eta)
    phi_x_fd = (
        basis.evaluate(x + h, eta).phi - basis.evaluate(x - h, eta).phi
    ) / (2.0 * h)
    phi_eta_fd = (
        basis.evaluate(x, eta + h).phi - basis.evaluate(x, eta - h).phi
    ) / (2.0 * h)

    assert np.max(np.abs(analytic.phi_x - phi_x_fd)) < 3e-9
    assert np.max(np.abs(analytic.phi_eta - phi_eta_fd)) < 3e-9


def test_compact_support_and_boundary_jets_are_exactly_zero():
    basis = Eq45CompactProfileBasis.seed()
    x = np.array([basis.x_cut, basis.x_cut + 0.2, 0.5, 0.5])
    eta = np.array([0.0, 0.0, basis.eta_cut, -basis.eta_cut])
    jets = basis.evaluate(x, eta).stacked()
    assert np.array_equal(jets, np.zeros_like(jets))

    outside = basis.evaluate(
        np.array([0.25, 1.0, 2.0]),
        np.array([1.2, -1.1, 1.5]),
    ).stacked()
    assert np.array_equal(outside, np.zeros_like(outside))


def test_serialization_round_trip_and_fail_closed_bounds(tmp_path):
    basis = Eq45CompactProfileBasis.seed()
    target = tmp_path / "profile.json"
    basis.save_json(target)
    restored = Eq45CompactProfileBasis.load_json(target)

    assert restored == basis
    probe_x = np.linspace(0.0, 3.5, 17)
    probe_eta = np.linspace(-0.9, 0.9, 17)
    assert np.array_equal(
        restored.evaluate(probe_x, probe_eta).stacked(),
        basis.evaluate(probe_x, probe_eta).stacked(),
    )

    payload = basis.to_dict()
    payload["phi_coefficients"][0] = basis.coefficient_limit + 0.01
    with pytest.raises(ValueError, match="coefficient_limit"):
        Eq45CompactProfileBasis.from_dict(payload)

    payload = basis.to_dict()
    payload["schema"] = "paper_exact_profile"
    with pytest.raises(ValueError, match="schema"):
        Eq45CompactProfileBasis.from_dict(payload)

    with pytest.raises(ValueError, match="nonnegative"):
        basis.evaluate(np.array([-1e-4]), np.array([0.0]))
