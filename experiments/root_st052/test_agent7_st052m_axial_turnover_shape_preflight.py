from __future__ import annotations

import math

import numpy as np
import pytest

import agent7_st052m_axial_turnover_shape_preflight as audit
import agent7_st052m_outer_axial_return_reservoir_preflight as reservoir


def test_live_turnover_exactly_replays_existing_outer_reservoir() -> None:
    receipt = audit.exact_live_shape_replay()
    assert receipt["passes"] is True
    assert receipt["axial_potential_derivative_max_abs_difference"] <= 1.0e-14
    assert receipt["unit_correction_max_abs_difference"] <= 1.0e-14


def test_turnover_sign_witness_has_frozen_direction() -> None:
    witness = audit.turnover_sign_witness()
    assert witness["passes"] is True
    minus = witness["by_turnover"]["s=1.75"]["mean_radial_component"]
    live = witness["by_turnover"]["s=1.80"]["mean_radial_component"]
    plus = witness["by_turnover"]["s=1.85"]["mean_radial_component"]
    assert minus > 0.0
    assert abs(live) <= witness["live_zero_tolerance"]
    assert plus < 0.0


def test_axial_family_keeps_exact_compact_support() -> None:
    points = np.asarray(
        [
            [1.60, 0.0, 1.50],
            [1.61, 0.0, 1.50],
            [0.90, 0.0, 1.00],
            [0.90, 0.0, -1.00],
            [0.90, 0.0, 2.00],
            [0.90, 0.0, -2.00],
            [0.90, 0.0, 2.05],
            [0.90, 0.0, -2.05],
        ],
        dtype=float,
    )
    for turnover in audit.S_VALUES:
        correction = audit.reservoir_correction_s(points, turnover)
        assert np.max(np.abs(correction)) == 0.0


def test_unit_corrections_are_divergence_free_to_frozen_fd_gate() -> None:
    rows = []
    for radius in (0.35, 0.90, 1.35):
        for abs_z in (1.30, 1.60, 1.92):
            for z_sign in (-1.0, 1.0):
                for k in range(8):
                    angle = k * math.pi / 4.0
                    rows.append(
                        [radius * math.cos(angle), radius * math.sin(angle), z_sign * abs_z]
                    )
    points = np.asarray(rows, dtype=float)
    h = audit.DIVERGENCE_STEP
    for turnover in audit.S_VALUES:
        div = np.zeros(len(points), dtype=float)
        for axis in range(3):
            shift = np.zeros(3, dtype=float)
            shift[axis] = h
            up = audit.reservoir_correction_s(points + shift, turnover)
            um = audit.reservoir_correction_s(points - shift, turnover)
            div += (up[:, axis] - um[:, axis]) / (2.0 * h)
        assert float(np.max(np.abs(div))) <= audit.DIVERGENCE_MAX


def test_live_axial_potential_matches_original_on_dense_nonsingular_grid() -> None:
    z = np.linspace(-2.05, 2.05, 4001)
    expected_z, expected_dz = reservoir.axial_potential_and_derivative(z)
    actual_z, actual_dz = audit.axial_potential_and_derivative_s(z, audit.S_LIVE)
    np.testing.assert_allclose(actual_z, expected_z, rtol=0.0, atol=2.0e-14)
    np.testing.assert_allclose(actual_dz, expected_dz, rtol=0.0, atol=2.0e-13)


def test_invalid_turnover_rejected() -> None:
    with pytest.raises(ValueError, match="strictly between"):
        audit.axial_potential_and_derivative_s(np.asarray([1.5]), 1.0)
    with pytest.raises(ValueError, match="strictly between"):
        audit.axial_potential_and_derivative_s(np.asarray([1.5]), 2.0)
