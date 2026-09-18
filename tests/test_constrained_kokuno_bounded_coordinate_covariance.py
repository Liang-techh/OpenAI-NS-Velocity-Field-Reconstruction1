from __future__ import annotations

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_bounded_coordinate_covariance import (
    EXPECTED_PARAMETER_UNIT,
    convert_amplitude_budget_to_fractional,
    covariance_tangent_from_velocity_and_tangents,
    measure_bounded_coordinate_covariance,
)
from openai_ns_reconstruction.kokuno_bounded_family_coefficients import (
    KokunoBoundedFamilyCoefficientCoordinates,
)


def _covariance(v):
    v = np.asarray(v, dtype=float)
    return np.array(
        [np.mean(v[..., 0] * v[..., 1]), np.mean(v[..., 0] * v[..., 2])],
        dtype=float,
    )


def _synthetic_family():
    angles = 2.0 * np.pi * np.arange(16, dtype=float) / 16.0
    w0 = np.stack(
        (
            0.18 + 0.03 * np.cos(angles),
            0.11 * np.sin(angles),
            0.07 + 0.02 * np.cos(2.0 * angles),
        ),
        axis=-1,
    )
    w1 = np.stack(
        (
            -0.04 + 0.02 * np.sin(angles),
            0.09 + 0.03 * np.cos(angles),
            -0.05 * np.sin(2.0 * angles),
        ),
        axis=-1,
    )
    velocity = np.stack((w0, w1), axis=-2)
    labels = [(5, (0, 0, 0)), (6, (1, 0, 0))]
    return velocity, labels


def test_frozen_amplitude_budget_converts_without_numeric_reuse():
    frozen = 0.012513055889617838
    converted = convert_amplitude_budget_to_fractional(frozen, 0.125)
    assert np.isclose(converted, 0.1001044471169427, rtol=0.0, atol=2e-16)
    assert converted != frozen

    with pytest.raises(ValueError, match="nonzero"):
        convert_amplitude_budget_to_fractional(frozen, 0.0)
    with pytest.raises(ValueError, match="positivity contract"):
        convert_amplitude_budget_to_fractional(0.08, 0.125)


def test_coordinate_covariance_tangents_match_independent_centered_difference():
    velocity, labels = _synthetic_family()
    contract = KokunoBoundedFamilyCoefficientCoordinates(max_l1_update=0.2)
    base = contract.evaluate(velocity, labels)
    total = base["velocity_physical_cylindrical_total_base"]
    tangents = base["velocity_parameter_tangent_cylindrical"]
    analytic = covariance_tangent_from_velocity_and_tangents(total, tangents)

    h = 1.0e-6
    for index, name in enumerate(("delta_common", "delta_band")):
        kwargs_plus = {name: h}
        kwargs_minus = {name: -h}
        plus = contract.evaluate(velocity, labels, **kwargs_plus)[
            "velocity_physical_cylindrical_total_modulated"
        ]
        minus = contract.evaluate(velocity, labels, **kwargs_minus)[
            "velocity_physical_cylindrical_total_modulated"
        ]
        fd = (_covariance(plus) - _covariance(minus)) / (2.0 * h)
        np.testing.assert_allclose(analytic[index], fd, rtol=3e-9, atol=3e-11)

    np.testing.assert_allclose(analytic[0], 2.0 * _covariance(total), rtol=0.0, atol=2e-15)


def test_one_band_has_no_fake_band_covariance_direction():
    velocity, _ = _synthetic_family()
    labels = [(5, (0, 0, 0)), (5, (1, 0, 0))]
    contract = KokunoBoundedFamilyCoefficientCoordinates(max_l1_update=0.15)
    out = contract.evaluate(velocity, labels)
    analytic = covariance_tangent_from_velocity_and_tangents(
        out["velocity_physical_cylindrical_total_base"],
        out["velocity_parameter_tangent_cylindrical"],
    )
    assert out["band_contrast_active"] is False
    np.testing.assert_array_equal(analytic[1], np.zeros(2))
    np.testing.assert_allclose(
        analytic[0],
        2.0 * _covariance(out["velocity_physical_cylindrical_total_base"]),
        rtol=0.0,
        atol=2e-15,
    )


def test_measurement_rejects_changed_coordinate_units():
    class BadUnits(KokunoBoundedFamilyCoefficientCoordinates):
        @staticmethod
        def parameter_units():
            return ("bad", "bad")

    def family(points, time, phase):
        points = np.asarray(points, dtype=float)
        r = np.hypot(points[:, 0], points[:, 1])
        theta = np.arctan2(points[:, 1], points[:, 0]) + float(phase)
        w0 = np.stack((0.1 + 0.01 * r, 0.03 * np.cos(theta), 0.02 * np.sin(theta)), axis=-1)
        w1 = np.stack((0.04 * np.sin(theta), 0.02 + 0.01 * r, -0.03 * np.cos(theta)), axis=-1)
        by_beta = np.stack((w0, w1), axis=-2)
        return {
            "beta_labels": ((5, (0, 0, 0)), (6, (1, 0, 0))),
            "velocity_physical_cylindrical_by_beta": by_beta,
            "velocity_physical_cylindrical_total": np.sum(by_beta, axis=-2),
        }

    with pytest.raises(ValueError, match="parameter units"):
        measure_bounded_coordinate_covariance(
            family,
            np.linspace(0.1, 0.3, 9),
            coordinate_contract=BadUnits(max_l1_update=0.1),
            angular_count=8,
            phase_count=4,
        )


def test_declared_repository_unit_is_exact_agent2_unit():
    contract = KokunoBoundedFamilyCoefficientCoordinates(max_l1_update=0.1)
    assert contract.parameter_units() == (EXPECTED_PARAMETER_UNIT, EXPECTED_PARAMETER_UNIT)
