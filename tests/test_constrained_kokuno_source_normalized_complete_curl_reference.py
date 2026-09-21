from __future__ import annotations

import inspect

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_corrected_oscillation_source_ledger import (
    default_kokuno_corrected_oscillation_source_ledger,
)
from openai_ns_reconstruction.kokuno_source_normalized_complete_curl_reference import (
    PARENT_A2_HEAD,
    PARENT_SOURCE_LEDGER_BLOB,
    NormalizedPotentialJet,
    complete_harmonic_amplitude_from_coefficient_jet,
    harmonic_vector_potential_coefficient,
    normalized_curl_from_potential_jet,
    normalized_divergence_from_vector_jet,
    source_reference_contract,
)


def test_source_reference_is_bound_to_exact_corrected_ledger_parent() -> None:
    contract = source_reference_contract()
    ledger = default_kokuno_corrected_oscillation_source_ledger()

    assert PARENT_A2_HEAD == "bea6517b0f5a827bb2b35f05736e66a54177ed23"
    assert PARENT_SOURCE_LEDGER_BLOB == "241ee6a9210aa2f4c8e3e9678e7b2175731c2b4a"
    assert contract["source_commit"] == ledger.source_commit
    assert contract["source_workbench_blob"] == ledger.source_workbench_blob
    assert contract["source_normalized_curl_formula_executable"] is True
    assert contract["source_complete_harmonic_amplitude_algebra_executable"] is True


def test_normalized_curl_replays_manufactured_auxiliary_independent_field() -> None:
    # Auxiliary-independent specialization: D_r=partial_R, D_z=epsilon partial_Z.
    radius = np.array([0.7, 1.1, 1.8, 2.4])
    theta = np.array([-0.8, -0.1, 0.6, 1.3])
    zed = np.array([-0.4, 0.2, 0.5, -0.3])
    epsilon = 0.17

    a_r = radius**2 * np.sin(theta) + zed
    a_theta = radius * zed * np.cos(theta)
    a_z = radius**2 * zed + np.cos(theta)
    jet = NormalizedPotentialJet(
        a_r=a_r,
        a_theta=a_theta,
        a_z=a_z,
        dtheta_a_r=radius**2 * np.cos(theta),
        dtheta_a_z=-np.sin(theta),
        dr_a_z=2.0 * radius * zed,
        dr_a_theta=zed * np.cos(theta),
        dz_a_r=epsilon * np.ones_like(radius),
        dz_a_theta=epsilon * radius * np.cos(theta),
    )

    actual = normalized_curl_from_potential_jet(radius, jet)
    expected = np.stack(
        (
            -np.sin(theta) / radius - epsilon * radius * np.cos(theta),
            epsilon - 2.0 * radius * zed,
            (2.0 * zed - radius) * np.cos(theta),
        ),
        axis=-1,
    )
    np.testing.assert_allclose(actual, expected, rtol=0.0, atol=2.0e-14)


def test_manufactured_normalized_curl_has_exact_source_divergence_cancellation() -> None:
    radius = np.array([0.62, 0.95, 1.37, 2.05])
    theta = np.array([-0.9, -0.25, 0.45, 1.05])
    epsilon = 0.13

    curl_r = -np.sin(theta) / radius - epsilon * radius * np.cos(theta)
    dr_curl_r = np.sin(theta) / radius**2 - epsilon * np.cos(theta)
    dtheta_curl_theta = np.zeros_like(radius)
    dz_curl_z = 2.0 * epsilon * np.cos(theta)

    divergence = normalized_divergence_from_vector_jet(
        radius,
        curl_r,
        dr_curl_r,
        dtheta_curl_theta,
        dz_curl_z,
    )
    np.testing.assert_allclose(divergence, 0.0, rtol=0.0, atol=2.0e-14)


def test_harmonic_coefficient_matches_corrected_source_cross_product_formula() -> None:
    n_phi = np.array([2.0, -1.0, 0.5])
    t_m = np.array([1.0, 2.0, 0.0], dtype=np.complex128)

    actual = harmonic_vector_potential_coefficient(n_phi, t_m, k=5.0, m=1)
    expected = 1j * np.cross(n_phi, t_m) / (5.0 * np.dot(n_phi, n_phi))
    np.testing.assert_allclose(actual, expected, rtol=0.0, atol=2.0e-15)
    np.testing.assert_allclose(np.dot(n_phi, t_m), 0.0, rtol=0.0, atol=0.0)


def test_complete_harmonic_amplitude_retains_source_longitudinal_remainder() -> None:
    n_phi = np.array([2.0, -1.0, 0.5])
    t_m = np.array([1.0, 2.0, 0.0], dtype=np.complex128)
    dr_c = np.array([0.2 + 0.1j, 0.7 - 0.2j, -0.4 + 0.3j])
    dz_c = np.array([-0.3 + 0.2j, 0.6 + 0.1j, 0.5 - 0.4j])

    result = complete_harmonic_amplitude_from_coefficient_jet(
        2.0,
        n_phi,
        t_m,
        dr_c,
        dz_c,
        k=5.0,
        m=1,
    )
    expected_remainder = np.array(
        [
            -dz_c[1],
            dz_c[0] - dr_c[2],
            dr_c[1] + result.coefficient[1] / 2.0,
        ]
    )
    np.testing.assert_allclose(result.remainder, expected_remainder, rtol=0.0, atol=2.0e-15)
    np.testing.assert_allclose(result.amplitude, t_m + expected_remainder, rtol=0.0, atol=2.0e-15)
    np.testing.assert_allclose(
        result.longitudinal_amplitude,
        result.longitudinal_remainder,
        rtol=0.0,
        atol=2.0e-15,
    )
    assert abs(result.longitudinal_remainder) > 1.0e-3


def test_batch_broadcasting_is_shape_preserving() -> None:
    radius = np.array([[0.8], [1.3]])
    theta = np.array([[0.1, 0.4, 0.9]])
    zeros = np.zeros((2, 3))
    jet = NormalizedPotentialJet(
        a_r=zeros,
        a_theta=np.broadcast_to(radius, (2, 3)),
        a_z=zeros,
        dtheta_a_r=zeros,
        dtheta_a_z=zeros,
        dr_a_z=zeros,
        dr_a_theta=np.ones((2, 3)),
        dz_a_r=zeros,
        dz_a_theta=zeros,
    )
    result = normalized_curl_from_potential_jet(radius, jet)

    assert result.shape == (2, 3, 3)
    np.testing.assert_allclose(result[..., 0], 0.0)
    np.testing.assert_allclose(result[..., 1], 0.0)
    np.testing.assert_allclose(result[..., 2], 2.0)
    assert theta.shape == (1, 3)  # keep a nontrivial broadcast-axis witness in the test.


def test_annular_reference_fails_closed_at_axis_and_nonfinite_radius() -> None:
    jet = NormalizedPotentialJet(
        a_r=0.0,
        a_theta=0.0,
        a_z=0.0,
        dtheta_a_r=0.0,
        dtheta_a_z=0.0,
        dr_a_z=0.0,
        dr_a_theta=0.0,
        dz_a_r=0.0,
        dz_a_theta=0.0,
    )
    with pytest.raises(ValueError, match="radius > 0"):
        normalized_curl_from_potential_jet(0.0, jet)
    with pytest.raises(ValueError, match="finite"):
        normalized_curl_from_potential_jet(np.nan, jet)


def test_harmonic_reference_rejects_nontransverse_or_degenerate_inputs() -> None:
    with pytest.raises(ValueError, match="n_phi dot t_m"):
        harmonic_vector_potential_coefficient(
            [1.0, 0.0, 0.0], [1.0, 0.0, 0.0], k=4.0, m=1
        )
    with pytest.raises(ValueError, match="nonzero"):
        harmonic_vector_potential_coefficient(
            [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], k=4.0, m=1
        )
    with pytest.raises(ValueError, match="positive carrier"):
        harmonic_vector_potential_coefficient(
            [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], k=0.0, m=1
        )
    with pytest.raises(ValueError, match="nonzero integer"):
        harmonic_vector_potential_coefficient(
            [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], k=4.0, m=0
        )


def test_no_source_to_runtime_or_pde_promotion_is_exposed() -> None:
    contract = source_reference_contract()
    for key in (
        "current_runtime_phase_source_exact",
        "current_runtime_support_source_exact",
        "current_runtime_complete_curl_source_equivalence_verified",
        "source_to_runtime_parameter_map_complete",
        "paper_exact",
        "matched_pressure_materialized",
        "restricted_forcing_materialized",
        "complete_ns_residual_assessed",
        "residual_reduction_claimed",
        "pde_validated",
    ):
        assert contract[key] is False

    assert "tolerance" not in inspect.signature(
        harmonic_vector_potential_coefficient
    ).parameters
    assert "tolerance" not in inspect.signature(
        complete_harmonic_amplitude_from_coefficient_jet
    ).parameters
