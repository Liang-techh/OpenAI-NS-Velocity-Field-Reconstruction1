from dataclasses import replace

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_complete_curl import KokunoCompleteCurlCorrection
from openai_ns_reconstruction.kokuno_covariance_tangent_preflight import (
    complete_curl_phase_evaluator,
    measure_phase_mean_covariance_tangent,
)
from openai_ns_reconstruction.kokuno_missing_covariance_column_target import (
    build_missing_covariance_column_target,
)
from openai_ns_reconstruction.kokuno_signed_covariance_inverse import (
    PROFILE_TIME,
    ROUTED_OSCILLATORY_AMPLITUDE,
    ROUTED_OSCILLATORY_PHASE,
)
from openai_ns_reconstruction.kokuno_velocity_column_preflight import (
    evaluate_velocity_column_preflight,
    measure_multiplicative_covariance_response,
    scaled_phase_velocity,
)


def _analytic_unit_column(points, time, phase):
    points = np.asarray(points, dtype=float)
    x, y, _ = points.T
    r = np.hypot(x, y)
    c = np.cos(phase)
    er = np.column_stack((x / r, y / r, np.zeros_like(r)))
    et = np.column_stack((-y / r, x / r, np.zeros_like(r)))
    ez = np.column_stack((np.zeros_like(r), np.zeros_like(r), np.ones_like(r)))
    return (
        (r * c)[:, None] * er
        + (2.0 * r * c)[:, None] * et
        + (3.0 * r * c)[:, None] * ez
    )


def test_multiplicative_column_matches_analytic_covariance_derivative():
    radii = np.linspace(0.1, 0.4, 9)
    amplitude = 0.25
    response = measure_multiplicative_covariance_response(
        _analytic_unit_column,
        radii,
        amplitude=amplitude,
        time=0.5,
        z=0.08,
        angular_count=8,
        phase_count=16,
    )
    expected = np.column_stack(
        (2.0 * amplitude * radii**2, 3.0 * amplitude * radii**2)
    )
    assert np.max(np.abs(response - expected)) < 2.0e-14


def test_velocity_column_adapter_matches_existing_complete_curl_amplitude_tangent():
    radii = np.linspace(0.08, 0.36, 17)
    unit = replace(
        KokunoCompleteCurlCorrection(),
        amplitude=1.0,
        phase=ROUTED_OSCILLATORY_PHASE,
    )
    base = replace(unit, amplitude=ROUTED_OSCILLATORY_AMPLITUDE)
    new_response = measure_multiplicative_covariance_response(
        complete_curl_phase_evaluator(unit),
        radii,
        amplitude=ROUTED_OSCILLATORY_AMPLITUDE,
        time=PROFILE_TIME,
        z=0.08,
        angular_count=8,
        phase_count=8,
    )
    existing_response = measure_phase_mean_covariance_tangent(
        complete_curl_phase_evaluator(base),
        complete_curl_phase_evaluator(unit),
        radii,
        time=PROFILE_TIME,
        z=0.08,
        angular_count=8,
        phase_count=8,
    )
    scale = max(float(np.linalg.norm(existing_response)), np.finfo(float).tiny)
    assert np.linalg.norm(new_response - existing_response) / scale < 2.0e-13


def test_duplicate_velocity_column_is_rejected_before_finite_cycle():
    radii = np.linspace(0.1, 0.4, 9)
    current = measure_multiplicative_covariance_response(
        _analytic_unit_column,
        radii,
        amplitude=0.25,
        time=0.5,
        z=0.08,
        angular_count=8,
        phase_count=16,
    )
    target = current + np.column_stack((3.0e-3 * np.ones(9), -2.0e-3 * np.ones(9)))
    receipt = build_missing_covariance_column_target(radii, target, current)
    audit = evaluate_velocity_column_preflight(
        receipt,
        _analytic_unit_column,
        amplitude=0.25,
        time=0.5,
        z=0.08,
        angular_count=8,
        phase_count=16,
    )
    assert audit["finite_cycle_rerun_allowed"] is False
    bounded = audit["bounded_preflight"]["bounded_inverse"]
    assert bounded["all_required_nodes_rank2"] is False
    assert bounded["rank2_required_nodes"] == 0
    assert bounded["bounded_inverse_preflight_passed"] is False


def test_scaled_velocity_and_column_measurement_fail_closed():
    with pytest.raises(ValueError):
        scaled_phase_velocity(_analytic_unit_column, 0.0)
    with pytest.raises(ValueError):
        scaled_phase_velocity(_analytic_unit_column, np.nan)

    radii = np.linspace(0.1, 0.4, 9)

    def bad_shape(points, time, phase):
        return np.ones((len(points), 2))

    with pytest.raises(ValueError):
        measure_multiplicative_covariance_response(
            bad_shape, radii, amplitude=0.25, angular_count=8, phase_count=8
        )

    def bad_finite(points, time, phase):
        out = np.ones_like(points)
        out[0, 0] = np.nan
        return out

    with pytest.raises(ValueError):
        measure_multiplicative_covariance_response(
            bad_finite, radii, amplitude=0.25, angular_count=8, phase_count=8
        )
