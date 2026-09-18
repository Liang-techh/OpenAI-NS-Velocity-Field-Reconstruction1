import numpy as np
import pytest
from dataclasses import replace

from openai_ns_reconstruction.kokuno_complete_curl import KokunoCompleteCurlCorrection
from openai_ns_reconstruction.kokuno_covariance_tangent_preflight import (
    complete_curl_phase_evaluator,
    evaluate_covariance_tangent_preflight,
    measure_phase_mean_covariance_tangent,
)
from openai_ns_reconstruction.kokuno_missing_covariance_column_target import (
    build_missing_covariance_column_target,
)
from openai_ns_reconstruction.kokuno_phase_orbit_covariance_rank import (
    measure_phase_mean_covariance_vector,
)
from openai_ns_reconstruction.kokuno_signed_covariance_inverse import (
    PROFILE_TIME,
    ROUTED_OSCILLATORY_AMPLITUDE,
    ROUTED_OSCILLATORY_PHASE,
)


def test_complete_curl_amplitude_tangent_matches_covariance_derivative():
    radii = np.linspace(0.08, 0.36, 17)
    base = replace(
        KokunoCompleteCurlCorrection(),
        amplitude=ROUTED_OSCILLATORY_AMPLITUDE,
        phase=ROUTED_OSCILLATORY_PHASE,
    )
    tangent = replace(base, amplitude=1.0)
    measured = measure_phase_mean_covariance_tangent(
        complete_curl_phase_evaluator(base),
        complete_curl_phase_evaluator(tangent),
        radii,
        time=PROFILE_TIME,
        z=0.08,
        angular_count=8,
        phase_count=8,
    )
    unit_covariance = measure_phase_mean_covariance_vector(
        base,
        radii,
        phase_offset=0.0,
        time=PROFILE_TIME,
        z=0.08,
        angular_count=8,
        phase_count=8,
    )
    expected = 2.0 * ROUTED_OSCILLATORY_AMPLITUDE * unit_covariance
    scale = max(float(np.linalg.norm(expected)), np.finfo(float).tiny)
    assert np.linalg.norm(measured - expected) / scale < 2.0e-13


def test_existing_amplitude_tangent_is_rejected_as_second_direction():
    radii = np.linspace(0.08, 0.36, 9)
    current = np.column_stack((1.0 + 0.2 * radii, np.zeros(9)))
    target = current + np.column_stack((np.zeros(9), 0.15 + 0.05 * radii))
    receipt = build_missing_covariance_column_target(radii, target, current)
    audit = evaluate_covariance_tangent_preflight(receipt, current)
    assert audit["finite_cycle_rerun_allowed"] is False
    bounded = audit["bounded_inverse"]
    assert bounded["all_required_nodes_rank2"] is False
    assert bounded["rank2_required_nodes"] == 0
    assert bounded["bounded_inverse_preflight_passed"] is False


def test_product_rule_tangent_on_analytic_ring_field():
    radii = np.linspace(0.1, 0.4, 9)

    def base(points, time, phase):
        x, y, _ = points.T
        r = np.hypot(x, y)
        c = np.cos(phase)
        s = np.sin(phase)
        er = np.column_stack((x / r, y / r, np.zeros_like(r)))
        et = np.column_stack((-y / r, x / r, np.zeros_like(r)))
        ez = np.column_stack((np.zeros_like(r), np.zeros_like(r), np.ones_like(r)))
        return (
            (r * c)[:, None] * er
            + (2.0 * r * c)[:, None] * et
            + (3.0 * r * s)[:, None] * ez
        )

    def tangent(points, time, phase):
        x, y, _ = points.T
        r = np.hypot(x, y)
        c = np.cos(phase)
        s = np.sin(phase)
        er = np.column_stack((x / r, y / r, np.zeros_like(r)))
        et = np.column_stack((-y / r, x / r, np.zeros_like(r)))
        ez = np.column_stack((np.zeros_like(r), np.zeros_like(r), np.ones_like(r)))
        return (
            (0.5 * r * c)[:, None] * er
            - (r * c)[:, None] * et
            + (2.0 * r * s)[:, None] * ez
        )

    response = measure_phase_mean_covariance_tangent(
        base,
        tangent,
        radii,
        time=0.5,
        z=0.08,
        angular_count=8,
        phase_count=16,
    )
    assert np.max(np.abs(response)) < 2.0e-14


def test_tangent_measurement_fails_closed_on_bad_output():
    radii = np.linspace(0.1, 0.4, 9)

    def good(points, time, phase):
        return np.ones_like(points)

    def bad_shape(points, time, phase):
        return np.ones((len(points), 2))

    with pytest.raises(ValueError):
        measure_phase_mean_covariance_tangent(
            good, bad_shape, radii, angular_count=8, phase_count=8
        )

    def bad_finite(points, time, phase):
        out = np.ones_like(points)
        out[0, 0] = np.nan
        return out

    with pytest.raises(ValueError):
        measure_phase_mean_covariance_tangent(
            good, bad_finite, radii, angular_count=8, phase_count=8
        )
