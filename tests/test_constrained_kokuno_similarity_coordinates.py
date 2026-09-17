from __future__ import annotations

import copy

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_similarity_coordinates import (
    KokunoNativeSimilarityCoordinates,
)


def test_source_construction_h_bound_is_strict() -> None:
    for h in (0.0, -1.0e-4, 1.0e-2, 2.0e-2, np.nan):
        with pytest.raises(ValueError):
            KokunoNativeSimilarityCoordinates(h=h)


def test_native_q_relation_and_unique_source_domain_on_grid() -> None:
    coordinates = KokunoNativeSimilarityCoordinates()
    z = np.linspace(-2.0, 2.0, 65)[:, None]
    t = np.array([0.25, 0.375, 0.5, 0.625, 0.75])[None, :]
    q = coordinates.solve_q(z, t)
    q_star = coordinates.q_star(z)
    residual = q - z * z * q ** (2.0 * coordinates.h) - (1.0 - t)

    assert q.shape == (65, 5)
    assert np.all(np.isfinite(q))
    assert np.all(q > 0.0)
    assert np.all(q > q_star)
    assert np.max(np.abs(residual)) < 2.0e-12


def test_eta_and_similarity_relations_replay_source_formulas() -> None:
    coordinates = KokunoNativeSimilarityCoordinates(h=0.006)
    z = np.array([-1.4, -0.3, 0.0, 0.8, 1.7])
    t = np.array([0.3, 0.45, 0.6, 0.7, 0.4])
    q = coordinates.solve_q(z, t)
    eta = coordinates.eta(z, t)

    np.testing.assert_allclose(z, q**coordinates.D * eta, rtol=2.0e-13, atol=2.0e-13)
    np.testing.assert_allclose(1.0 - t, q * (1.0 - eta * eta), rtol=2.0e-13, atol=2.0e-13)
    assert np.all(np.abs(eta) < 1.0)

    q_axis = coordinates.solve_q(0.0, t)
    np.testing.assert_array_equal(q_axis, 1.0 - t)
    np.testing.assert_array_equal(coordinates.eta(0.0, t), np.zeros_like(t))


def test_vectorized_evaluate_broadcasts_and_returns_source_jacobian() -> None:
    coordinates = KokunoNativeSimilarityCoordinates()
    x = np.linspace(-0.8, 0.8, 7)[:, None]
    y = 0.25
    z = np.linspace(-1.0, 1.0, 5)[None, :]
    result = coordinates.evaluate(x, y, z, 0.5)

    assert set(result) == {
        "q", "eta", "X", "d", "L", "q_t", "eta_t", "X_t", "q_z", "eta_z", "X_z"
    }
    for values in result.values():
        assert values.shape == (7, 5)
        assert np.all(np.isfinite(values))
    assert np.all(result["d"] > 0.0)
    assert np.all(result["L"] > 0.0)
    assert np.all(result["X"] >= 0.0)


def test_exact_jacobian_matches_independent_centered_differences() -> None:
    coordinates = KokunoNativeSimilarityCoordinates(h=0.004)
    x, y, z, t = 0.43, -0.31, 0.72, 0.46
    base = coordinates.evaluate(x, y, z, t)
    eps = 1.0e-6

    for field, derivative in (("q", "q_t"), ("eta", "eta_t"), ("X", "X_t")):
        plus = coordinates.evaluate(x, y, z, t + eps)[field]
        minus = coordinates.evaluate(x, y, z, t - eps)[field]
        finite_difference = (plus - minus) / (2.0 * eps)
        np.testing.assert_allclose(finite_difference, base[derivative], rtol=2.0e-7, atol=2.0e-9)

    for field, derivative in (("q", "q_z"), ("eta", "eta_z"), ("X", "X_z")):
        plus = coordinates.evaluate(x, y, z + eps, t)[field]
        minus = coordinates.evaluate(x, y, z - eps, t)[field]
        finite_difference = (plus - minus) / (2.0 * eps)
        np.testing.assert_allclose(finite_difference, base[derivative], rtol=2.0e-7, atol=2.0e-9)


def test_solver_is_stable_under_tighter_numerical_tolerances() -> None:
    ordinary = KokunoNativeSimilarityCoordinates(rtol=1.0e-11, atol=1.0e-13)
    tight = KokunoNativeSimilarityCoordinates(rtol=1.0e-14, atol=1.0e-15)
    z = np.linspace(-2.0, 2.0, 81)[:, None]
    t = np.linspace(0.25, 0.75, 9)[None, :]

    q_ordinary = ordinary.solve_q(z, t)
    q_tight = tight.solve_q(z, t)
    np.testing.assert_allclose(q_ordinary, q_tight, rtol=3.0e-12, atol=3.0e-13)
    assert np.max(np.abs(tight.relation_residual(z, t))) < 2.0e-13


def test_native_relation_is_not_silently_relabelled_as_repo_negative_exponent_convention() -> None:
    coordinates = KokunoNativeSimilarityCoordinates(h=0.005)
    z = np.array([0.5, 1.0, 1.5])
    t = 0.4
    q = coordinates.solve_q(z, t)
    native_residual = q - z * z * q ** (2.0 * coordinates.h) - (1.0 - t)
    other_residual = q - z * z * q ** (-2.0 * coordinates.h) - (1.0 - t)

    assert np.max(np.abs(native_residual)) < 1.0e-12
    assert np.max(np.abs(other_residual)) > 1.0e-4
    payload = coordinates.to_payload()
    assert payload["truth_boundary"]["cross_convention_equivalence_claimed"] is False


def test_serialization_round_trip_and_fail_closed_metadata(tmp_path) -> None:
    coordinates = KokunoNativeSimilarityCoordinates(h=0.006, rtol=1.0e-12, atol=1.0e-14)
    path = coordinates.save_json(tmp_path / "coordinates.json")
    replay = KokunoNativeSimilarityCoordinates.load_json(path)

    assert replay == coordinates
    assert replay.sha256 == coordinates.sha256
    np.testing.assert_array_equal(replay.solve_q([0.0, 0.4, 1.1], 0.55), coordinates.solve_q([0.0, 0.4, 1.1], 0.55))

    payload = coordinates.to_payload()
    payload["sha256"] = coordinates.sha256
    mutated = copy.deepcopy(payload)
    mutated["coordinates"]["mapping_status"] = "equivalent"
    with pytest.raises(ValueError):
        KokunoNativeSimilarityCoordinates.from_payload(mutated)

    mutated = copy.deepcopy(payload)
    mutated["solver"]["method"] = "unbounded generic root search"
    with pytest.raises(ValueError):
        KokunoNativeSimilarityCoordinates.from_payload(mutated)


def test_invalid_time_nonfinite_and_solver_settings_fail_closed() -> None:
    coordinates = KokunoNativeSimilarityCoordinates()
    for t in (1.0, 1.1, np.nan):
        with pytest.raises(ValueError):
            coordinates.solve_q(0.5, t)
    with pytest.raises(ValueError):
        coordinates.solve_q(np.inf, 0.5)
    with pytest.raises(ValueError):
        KokunoNativeSimilarityCoordinates(rtol=0.0)
    with pytest.raises(ValueError):
        KokunoNativeSimilarityCoordinates(atol=0.0)
    with pytest.raises(ValueError):
        KokunoNativeSimilarityCoordinates(max_iterations=4)
