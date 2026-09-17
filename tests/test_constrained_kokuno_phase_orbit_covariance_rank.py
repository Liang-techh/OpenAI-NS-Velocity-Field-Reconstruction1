from dataclasses import replace
from math import pi

import numpy as np

from openai_ns_reconstruction.kokuno_complete_curl import KokunoCompleteCurlCorrection
from openai_ns_reconstruction.kokuno_phase_orbit_covariance_rank import (
    analyze_response_matrix,
    measure_phase_mean_covariance_vector,
    source_reference_determinant,
    source_reference_response_matrix,
)


def test_source_reference_matrix_has_exact_nonzero_determinant():
    matrix = source_reference_response_matrix(
        A_c=2.0,
        u_star=3.0,
        h_plus=0.4,
        h_minus=0.7,
    )
    expected = source_reference_determinant(
        A_c=2.0,
        u_star=3.0,
        h_plus=0.4,
        h_minus=0.7,
    )
    np.testing.assert_allclose(np.linalg.det(matrix), expected, rtol=1.0e-14, atol=1.0e-14)
    assert np.linalg.matrix_rank(matrix) == 2
    assert expected < 0.0


def test_rank_audit_distinguishes_independent_and_duplicate_columns():
    independent = np.asarray(
        [
            [1.0, 0.0],
            [0.0, 1.0],
            [1.0, 1.0],
            [2.0, -1.0],
        ]
    )
    target = np.asarray([0.3, -0.2, 0.1, 0.8])
    full = analyze_response_matrix(independent, target)
    assert full.rank == 2
    assert np.isfinite(full.condition)
    assert full.relative_target_projection_rms < 1.0e-14

    duplicate = np.column_stack((independent[:, 0], independent[:, 0]))
    rank_one = analyze_response_matrix(duplicate, target)
    assert rank_one.rank == 1
    assert np.isinf(rank_one.condition)
    assert rank_one.relative_target_projection_rms > 0.1


def test_uniform_phase_average_makes_phase_copies_same_covariance_column():
    correction = replace(KokunoCompleteCurlCorrection(), amplitude=0.125, phase=0.0)
    radii = np.linspace(0.10, 0.30, 9)
    first = measure_phase_mean_covariance_vector(
        correction,
        radii,
        phase_offset=0.0,
        time=0.5,
        z=0.08,
        angular_count=8,
        phase_count=8,
    )
    second = measure_phase_mean_covariance_vector(
        correction,
        radii,
        phase_offset=0.5 * pi,
        time=0.5,
        z=0.08,
        angular_count=8,
        phase_count=8,
    )
    assert np.linalg.norm(first) > 1.0e-8
    np.testing.assert_allclose(second, first, rtol=2.0e-13, atol=2.0e-15)

    response0 = 2.0 * correction.amplitude * first[1:-1].reshape(-1)
    response1 = 2.0 * correction.amplitude * second[1:-1].reshape(-1)
    matrix = np.column_stack((response0, response1))
    audit = analyze_response_matrix(matrix, response0)
    assert audit.rank == 1
    assert np.isinf(audit.condition)


def test_contract_rejects_invalid_source_or_rank_inputs():
    with np.testing.assert_raises(ValueError):
        source_reference_response_matrix(A_c=1.0, u_star=0.0, h_plus=1.0, h_minus=1.0)
    with np.testing.assert_raises(ValueError):
        analyze_response_matrix(np.zeros((4, 2)), np.ones(4))
    correction = KokunoCompleteCurlCorrection()
    with np.testing.assert_raises(ValueError):
        measure_phase_mean_covariance_vector(
            correction,
            np.asarray([0.1, 0.2, 0.3]),
            phase_offset=0.0,
            angular_count=7,
            phase_count=8,
        )
