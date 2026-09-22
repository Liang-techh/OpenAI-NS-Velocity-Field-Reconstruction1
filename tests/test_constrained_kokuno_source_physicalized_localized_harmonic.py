from __future__ import annotations

import inspect
import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_source_physicalized_localized_harmonic import (
    TASK,
    physical_to_source_chart_rzt,
    physicalize_source_localized_harmonic,
    source_chart_to_physical_rzt,
    source_physical_scaling,
    source_physicalized_localized_harmonic_contract,
)
from openai_ns_reconstruction.kokuno_source_zero_data_localized_harmonic import (
    SourceLocalizedZeroDataHarmonic,
)


def _fake_harmonic() -> SourceLocalizedZeroDataHarmonic:
    complex_potential = np.asarray(
        [[1.0 + 2.0j, -0.5j, 0.25 - 0.75j], [0.4j, 0.7 + 0.2j, -1.1j]],
        dtype=np.complex128,
    )
    complex_velocity = np.asarray(
        [[0.5 - 0.25j, 0.2 + 0.7j, -0.3j], [-0.1 + 0.6j, 0.9j, 0.4 + 0.1j]],
        dtype=np.complex128,
    )
    real_cyl = 2.0 * np.real(complex_velocity)
    # Any finite Cartesian basis-rotated witness is enough here; the parent
    # module separately owns and tests the cylindrical->Cartesian rotation.
    real_cart = np.asarray([[0.8, -0.2, 0.0], [0.1, 0.6, 0.8]], dtype=float)
    return SourceLocalizedZeroDataHarmonic(
        amplitude_sensitivity=None,
        phase_jet=None,
        coefficient_jet=None,
        complete_amplitude=None,
        localization_weight=np.ones(2),
        dr_localization_weight=np.zeros(2),
        dz_localization_weight=np.zeros(2),
        localized_coefficient=np.zeros((2, 3), dtype=np.complex128),
        localized_dr_coefficient=np.zeros((2, 3), dtype=np.complex128),
        localized_dz_coefficient=np.zeros((2, 3), dtype=np.complex128),
        complex_vector_potential=complex_potential,
        complex_velocity_cylindrical=complex_velocity,
        real_pair_velocity_cylindrical=real_cyl,
        real_pair_velocity_chart_cartesian=real_cart,
    )


def test_corrected_scaling_formula_and_curl_factor_close() -> None:
    scaling = source_physical_scaling(ell=20, h=0.006)
    assert scaling.A == pytest.approx(0.506, abs=0.0)
    assert scaling.D == pytest.approx(0.494, abs=0.0)
    assert scaling.Q == math.ldexp(1.0, -20)
    assert scaling.epsilon == pytest.approx(scaling.Q**0.006, rel=2e-15)
    assert scaling.vector_potential_scale == pytest.approx(
        scaling.Q ** (0.5 - scaling.A), rel=2e-15
    )
    assert scaling.velocity_scale == pytest.approx(scaling.Q ** (-scaling.A), rel=2e-15)
    assert scaling.vector_potential_scale * scaling.curl_derivative_scale == pytest.approx(
        scaling.velocity_scale, rel=4e-15
    )
    assert scaling.Q ** (-scaling.D) == pytest.approx(
        scaling.curl_derivative_scale * scaling.epsilon, rel=4e-15
    )


def test_manufactured_physical_curl_matches_scaled_normalized_curl() -> None:
    # Independent analytic witness: A_*=(0,R Z,0), axisymmetric.
    # curl_* A_* = (-epsilon R, 0, 2 Z).
    scaling = source_physical_scaling(ell=16, h=0.004)
    R = 0.73
    Z = -0.41
    normalized = np.asarray((-scaling.epsilon * R, 0.0, 2.0 * Z))
    scaled = scaling.velocity_scale * normalized

    # Physical A_theta = Q^(1/2-A) R Z, with
    # R=r/Q^(1/2), Z=z/Q^D. Direct physical cylindrical curl gives:
    direct = np.asarray(
        (
            -scaling.vector_potential_scale * scaling.Q ** (-scaling.D) * R,
            0.0,
            2.0 * scaling.vector_potential_scale * scaling.Q ** (-0.5) * Z,
        )
    )
    np.testing.assert_allclose(direct, scaled, rtol=5e-15, atol=0.0)


def test_chart_physical_roundtrip_is_batch_safe() -> None:
    R = np.asarray((0.4, 0.9, 1.3))
    Z = np.asarray((-0.5, 0.0, 0.7))
    T = np.asarray((0.2, 0.8, 1.1))
    r, z, t = source_chart_to_physical_rzt(R, Z, T, ell=12, h=0.003)
    R2, Z2, T2 = physical_to_source_chart_rzt(r, z, t, ell=12, h=0.003)
    np.testing.assert_allclose(R2, R, rtol=2e-15, atol=0.0)
    np.testing.assert_allclose(Z2, Z, rtol=2e-15, atol=0.0)
    np.testing.assert_allclose(T2, T, rtol=2e-12, atol=2e-12)


def test_physicalize_scales_potential_and_velocity_without_reinterpreting_parent() -> None:
    harmonic = _fake_harmonic()
    result = physicalize_source_localized_harmonic(harmonic, ell=14, h=0.005)
    s = result.scaling
    assert result.source_harmonic is harmonic
    np.testing.assert_allclose(
        result.complex_vector_potential_physical,
        s.vector_potential_scale * harmonic.complex_vector_potential,
        rtol=0.0,
        atol=0.0,
    )
    np.testing.assert_allclose(
        result.complex_velocity_cylindrical_physical,
        s.velocity_scale * harmonic.complex_velocity_cylindrical,
        rtol=0.0,
        atol=0.0,
    )
    np.testing.assert_allclose(
        result.real_pair_velocity_cylindrical_physical,
        s.velocity_scale * harmonic.real_pair_velocity_cylindrical,
        rtol=0.0,
        atol=0.0,
    )
    np.testing.assert_allclose(
        result.real_pair_velocity_cartesian_physical,
        s.velocity_scale * harmonic.real_pair_velocity_chart_cartesian,
        rtol=0.0,
        atol=0.0,
    )


def test_contract_keeps_source_runtime_and_pde_boundaries_fail_closed() -> None:
    contract = source_physicalized_localized_harmonic_contract()
    assert contract["task"] == TASK == "K2-OSC-102"
    assert contract["source_chart_rzt_map_materialized"] is True
    assert contract["source_physical_Q_scaling_applied"] is True
    assert contract["repository_curl_scaling_identity_verified"] is True
    assert contract["physical_real_pair_velocity_materialized"] is True
    assert contract["caller_supplies_mode_forcing_directional_jet"] is True
    assert contract["caller_supplies_corrected_background"] is True
    assert contract["project_domain_source_input_provider_materialized"] is False
    assert contract["global_axis_safe_source_velocity_materialized"] is False
    assert contract["self_contained_velocity_xyzt_provider"] is False
    assert contract["complete_ns_residual_assessed"] is False
    assert contract["residual_reduction_claimed"] is False
    assert contract["paper_exact"] is False
    assert contract["pde_validated"] is False


def test_fail_closed_inputs_and_no_scientific_tuning_surface() -> None:
    for bad_ell in (0, -1, 1.5, True):
        with pytest.raises(ValueError):
            source_physical_scaling(ell=bad_ell, h=0.004)  # type: ignore[arg-type]
    for bad_h in (0.0, 0.01, -1e-3, math.inf, math.nan):
        with pytest.raises(ValueError):
            source_physical_scaling(ell=8, h=bad_h)
    with pytest.raises(ValueError):
        source_physical_scaling(ell=5000, h=0.004)
    with pytest.raises(ValueError):
        source_chart_to_physical_rzt(np.nan, 0.0, 0.0, ell=8, h=0.004)
    with pytest.raises(ValueError):
        physicalize_source_localized_harmonic(object(), ell=8, h=0.004)  # type: ignore[arg-type]

    forbidden = {
        "residual",
        "forcing",
        "pressure",
        "gain",
        "damping",
        "tolerance",
        "rtol",
        "atol",
        "panels",
        "steps",
        "target",
    }
    for fn in (
        source_physical_scaling,
        source_chart_to_physical_rzt,
        physical_to_source_chart_rzt,
        physicalize_source_localized_harmonic,
    ):
        assert forbidden.isdisjoint(inspect.signature(fn).parameters)
