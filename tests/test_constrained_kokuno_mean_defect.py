from dataclasses import replace

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_complete_curl import KokunoCompleteCurlCorrection
from openai_ns_reconstruction.kokuno_mean_defect import (
    PhaseMeanDefectContract,
    cylindrical_theta_component,
)


def _points():
    return np.asarray(
        [
            [0.21, 0.16, 0.11],
            [-0.31, 0.19, -0.17],
            [0.38, -0.26, 0.22],
            [-0.24, -0.33, 0.07],
            [0.14, -0.41, -0.28],
            [0.43, 0.24, 0.19],
        ],
        dtype=float,
    )


def _zero_base(points, time):
    del time
    return np.zeros_like(points)


def _constant_base(points, time):
    del time
    value = np.asarray([0.17, -0.11, 0.07])
    return np.broadcast_to(value, np.asarray(points).shape).copy()


def _correction(amplitude=0.2, phase=0.31):
    return KokunoCompleteCurlCorrection(
        amplitude=amplitude,
        wave_vector=(5.0, -3.0, 4.0),
        polarization=(1.0, 2.0, 1.0),
        center=(0.0, 0.0, 0.0),
        half_widths=(0.9, 0.8, 0.7),
        omega=2.3,
        phase=phase,
    )


def test_phase_mean_increment_scales_quadratically_with_complete_curl_amplitude():
    contract = PhaseMeanDefectContract(phase_count=8, step=0.005)
    small = contract.evaluate(_zero_base, _correction(amplitude=0.2), _points(), 0.5)
    large = contract.evaluate(_zero_base, _correction(amplitude=0.4), _points(), 0.5)

    assert np.linalg.norm(small.mean_defect_increment) > 0.0
    assert np.allclose(large.mean_defect_increment, 4.0 * small.mean_defect_increment, rtol=2.0e-9, atol=2.0e-10)


def test_uniform_phase_mean_is_invariant_under_initial_phase_and_cancels_linear_cross_terms():
    contract = PhaseMeanDefectContract(phase_count=8, step=0.005)
    first = contract.evaluate(_constant_base, _correction(phase=-0.7), _points(), 0.5)
    second = contract.evaluate(_constant_base, _correction(phase=1.1), _points(), 0.5)
    zero_base = contract.evaluate(_zero_base, _correction(phase=0.2), _points(), 0.5)

    assert np.allclose(first.mean_defect_increment, second.mean_defect_increment, rtol=3.0e-9, atol=3.0e-10)
    assert np.allclose(first.mean_defect_increment, zero_base.mean_defect_increment, rtol=3.0e-9, atol=3.0e-10)


def test_eight_and_sixteen_phase_projectors_agree_for_quadratic_ns_operator():
    coarse = PhaseMeanDefectContract(phase_count=8).evaluate(_constant_base, _correction(), _points(), 0.5)
    fine = PhaseMeanDefectContract(phase_count=16).evaluate(_constant_base, _correction(), _points(), 0.5)

    assert np.allclose(coarse.mean_defect_increment, fine.mean_defect_increment, rtol=3.0e-9, atol=3.0e-10)


def test_theta_projection_and_typed_contract_fail_closed():
    points = _points()
    vectors = np.column_stack((-points[:, 1], points[:, 0], np.zeros(len(points))))
    theta = cylindrical_theta_component(vectors, points)
    assert np.allclose(theta, np.hypot(points[:, 0], points[:, 1]))

    contract = PhaseMeanDefectContract()
    with pytest.raises(TypeError):
        contract.evaluate(np.zeros((4, 3)), _correction(), points, 0.5)
    with pytest.raises(TypeError):
        contract.evaluate(_zero_base, object(), points, 0.5)
    with pytest.raises(ValueError):
        replace(contract, phase_count=3)
    with pytest.raises(ValueError):
        cylindrical_theta_component(np.zeros((1, 3)), np.asarray([[0.0, 0.0, 0.2]]))
