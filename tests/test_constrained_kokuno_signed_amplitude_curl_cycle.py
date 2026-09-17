import numpy as np

from openai_ns_reconstruction.kokuno_complete_curl import KokunoCompleteCurlCorrection
from openai_ns_reconstruction.kokuno_signed_amplitude_curl_cycle import (
    CHEBYSHEV_DEGREE,
    KokunoSignedAmplitudeCurlCorrection,
    fit_compact_chebyshev_signed_amplitude,
)
from openai_ns_reconstruction.kokuno_signed_covariance_inverse import (
    build_rank1_signed_amplitude_profile,
)


def _manufactured_sample():
    radii = np.linspace(0.08, 0.36, 33)
    s = 2.0 * (radii - radii[0]) / (radii[-1] - radii[0]) - 1.0
    gate = np.maximum(1.0 - s * s, 0.0) ** 5
    delta = gate * (0.003 - 0.0012 * s + 0.0007 * (2.0 * s * s - 1.0))
    unit_covariance = 0.30 + 0.04 * np.cos(2.0 * radii)
    derivative = 2.0 * 0.125 * unit_covariance
    target_stress = derivative * delta
    return build_rank1_signed_amplitude_profile(
        radii,
        target_stress,
        unit_covariance,
        base_amplitude=0.125,
    )


def test_compact_chebyshev_adapter_recovers_manufactured_c4_profile():
    sampled = _manufactured_sample()
    compact = fit_compact_chebyshev_signed_amplitude(sampled, degree=CHEBYSHEV_DEGREE)
    metrics = compact.metrics(sampled.radii, sampled.delta_amplitude)
    # The upstream signed-profile contract intentionally zeros target nodes below
    # its predeclared active-stress floor, so the sampled target is not exactly
    # the generating polynomial at the first few collar nodes.  The compact fit
    # should nevertheless recover it far inside the 10% production guard.
    assert metrics["sampled_fit_relative_rms"] < 5.0e-4
    assert metrics["inner_value"] == 0.0
    assert metrics["outer_value"] == 0.0
    assert metrics["inner_radial_derivative"] == 0.0
    assert metrics["outer_radial_derivative"] == 0.0
    outside = compact(np.asarray([0.01, 0.07, 0.37, 0.8]))
    assert np.array_equal(outside, np.zeros_like(outside))


def _curl_from_potential(correction, points, time, step=2.0e-6):
    jac = np.empty((len(points), 3, 3), dtype=float)
    for axis in range(3):
        offset = np.eye(3)[axis] * step
        plus = correction.vector_potential(
            *(points + offset).T,
            np.full(len(points), time),
        )
        minus = correction.vector_potential(
            *(points - offset).T,
            np.full(len(points), time),
        )
        jac[:, :, axis] = (plus - minus) / (2.0 * step)
    return np.column_stack((
        jac[:, 2, 1] - jac[:, 1, 2],
        jac[:, 0, 2] - jac[:, 2, 0],
        jac[:, 1, 0] - jac[:, 0, 1],
    ))


def test_signed_velocity_is_full_curl_of_spatially_scaled_agent2_potential():
    sampled = _manufactured_sample()
    compact = fit_compact_chebyshev_signed_amplitude(sampled)
    primary = KokunoCompleteCurlCorrection(amplitude=0.125, phase=0.17)
    correction = KokunoSignedAmplitudeCurlCorrection(primary=primary, radial_amplitude=compact)
    points = np.asarray([
        [0.14, 0.03, 0.05],
        [0.17, -0.08, -0.04],
        [-0.19, 0.06, 0.02],
        [0.21, 0.09, -0.03],
    ])
    analytic = correction.at_points(points, 0.5)
    numeric = _curl_from_potential(correction, points, 0.5)
    scale = max(float(np.sqrt(np.mean(analytic * analytic))), 1.0e-14)
    assert float(np.sqrt(np.mean((analytic - numeric) ** 2))) / scale < 2.0e-7


def test_signed_complete_curl_is_divergence_free_and_nontrivial():
    sampled = _manufactured_sample()
    compact = fit_compact_chebyshev_signed_amplitude(sampled)
    primary = KokunoCompleteCurlCorrection(amplitude=0.125, phase=0.0)
    correction = KokunoSignedAmplitudeCurlCorrection(primary=primary, radial_amplitude=compact)
    points = np.asarray([
        [0.13, 0.04, 0.01],
        [0.16, -0.06, -0.02],
        [-0.18, 0.07, 0.03],
        [0.20, 0.08, -0.04],
    ])
    step = 2.0e-6
    divergence = np.zeros(len(points), dtype=float)
    for axis in range(3):
        offset = np.eye(3)[axis] * step
        plus = correction.at_points(points + offset, 0.5)
        minus = correction.at_points(points - offset, 0.5)
        divergence += ((plus - minus) / (2.0 * step))[:, axis]
    assert float(np.max(np.abs(divergence))) < 2.0e-7
    assert float(np.max(np.linalg.norm(correction.at_points(points, 0.5), axis=1))) > 1.0e-6


def test_metadata_keeps_source_exact_and_pde_claims_out_of_velocity_adapter():
    sampled = _manufactured_sample()
    compact = fit_compact_chebyshev_signed_amplitude(sampled)
    correction = KokunoSignedAmplitudeCurlCorrection(
        primary=KokunoCompleteCurlCorrection(amplitude=0.125),
        radial_amplitude=compact,
    )
    metadata = correction.metadata()
    assert metadata["agent2_curl_reimplemented"] is False
    assert metadata["source_exact_two_column_map"] is False
    assert metadata["chebyshev_degree"] == CHEBYSHEV_DEGREE
