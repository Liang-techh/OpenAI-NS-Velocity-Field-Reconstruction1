import numpy as np

from openai_ns_reconstruction.kokuno_complete_curl import KokunoCompleteCurlCorrection
from openai_ns_reconstruction.kokuno_signed_covariance_inverse import (
    SOURCE_COMMIT,
    SOURCE_DOI,
    build_rank1_signed_amplitude_profile,
    finite_difference_amplitude_derivative_error,
    measure_unit_radial_theta_covariance,
)


def _compact_target(radii: np.ndarray) -> np.ndarray:
    s = 2.0 * (radii - radii[0]) / (radii[-1] - radii[0]) - 1.0
    base = np.maximum(1.0 - s * s, 0.0)
    return 2.0e-4 * base**5 * (1.0 + 0.15 * s)


def test_rank1_signed_profile_reconstructs_well_conditioned_manufactured_target():
    radii = np.linspace(0.08, 0.36, 33)
    target = _compact_target(radii)
    unit_covariance = 0.25 + 0.05 * np.cos(3.0 * radii)
    profile = build_rank1_signed_amplitude_profile(
        radii,
        target,
        unit_covariance,
        base_amplitude=0.125,
    )
    metrics = profile.metrics()
    assert metrics["rank1_realizable_on_active_nodes"] is True
    assert metrics["unresolved_active_nodes"] == 0
    assert metrics["linear_reconstruction_relative_rms"] < 1.0e-12
    assert metrics["delta_amplitude_edge_max_abs"] == 0.0
    assert metrics["delta_amplitude_max_abs"] > 0.0


def test_rank1_signed_profile_fails_closed_when_response_loses_rank():
    radii = np.linspace(0.08, 0.36, 33)
    target = _compact_target(radii)
    response = np.ones_like(radii)
    middle = len(radii) // 2
    response[middle] = 0.0
    profile = build_rank1_signed_amplitude_profile(
        radii,
        target,
        response,
        base_amplitude=0.125,
    )
    metrics = profile.metrics()
    assert metrics["rank1_realizable_on_active_nodes"] is False
    assert metrics["unresolved_active_nodes"] >= 1
    assert metrics["linear_reconstruction_relative_rms"] > 0.0


def test_public_complete_curl_covariance_derivative_is_quadratic_in_amplitude():
    correction = KokunoCompleteCurlCorrection(amplitude=0.125, phase=0.0)
    radii = np.linspace(0.10, 0.30, 9)
    covariance = measure_unit_radial_theta_covariance(
        correction,
        radii,
        angular_count=8,
        phase_count=8,
    )
    assert covariance.shape == radii.shape
    assert np.isfinite(covariance).all()
    calibration = finite_difference_amplitude_derivative_error(
        correction,
        radii,
        angular_count=8,
        phase_count=8,
    )
    assert calibration["derivative_error_relative_rms"] < 1.0e-9


def test_source_provenance_is_pinned_to_corrected_public_reader():
    assert SOURCE_DOI == "10.5281/zenodo.22678406"
    assert SOURCE_COMMIT == "143f6773feb424ad9ed3a8d116653200f20346b7"
