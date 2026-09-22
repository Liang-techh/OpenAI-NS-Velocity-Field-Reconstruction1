from __future__ import annotations

import inspect

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_source_physicalized_localized_harmonic import (
    SourcePhysicalScaling,
)
from openai_ns_reconstruction.kokuno_source_axis_safe_physical_velocity import (
    TASK,
    AxisSafeSourceHarmonicProvider,
    KokunoProviderDrivenAxisSafePhysicalVelocity,
    source_axis_safe_physical_velocity_contract,
)
from openai_ns_reconstruction.kokuno_source_zero_data_localized_harmonic import (
    SourceLocalizedZeroDataHarmonic,
)


def _harmonic_from_chart_cartesian(chart_cartesian: np.ndarray) -> SourceLocalizedZeroDataHarmonic:
    chart_cartesian = np.asarray(chart_cartesian, dtype=float)
    if chart_cartesian.ndim == 1:
        chart_cartesian = chart_cartesian.reshape(1, 3)
    n = chart_cartesian.shape[0]
    complex_velocity = 0.5 * chart_cartesian.astype(np.complex128)
    return SourceLocalizedZeroDataHarmonic(
        amplitude_sensitivity=None,
        phase_jet=None,
        coefficient_jet=None,
        complete_amplitude=None,
        localization_weight=np.ones(n),
        dr_localization_weight=np.zeros(n),
        dz_localization_weight=np.zeros(n),
        localized_coefficient=np.zeros((n, 3), dtype=np.complex128),
        localized_dr_coefficient=np.zeros((n, 3), dtype=np.complex128),
        localized_dz_coefficient=np.zeros((n, 3), dtype=np.complex128),
        complex_vector_potential=np.zeros((n, 3), dtype=np.complex128),
        complex_velocity_cylindrical=complex_velocity,
        real_pair_velocity_cylindrical=chart_cartesian.copy(),
        real_pair_velocity_chart_cartesian=chart_cartesian.copy(),
    )


def _provider(call_log: list[dict[str, np.ndarray]], *, radius: float = 0.25):
    def evaluate(
        R: np.ndarray,
        theta: np.ndarray,
        Z: np.ndarray,
        T: np.ndarray,
        scaling: SourcePhysicalScaling,
    ) -> SourceLocalizedZeroDataHarmonic:
        call_log.append(
            {
                "R": np.asarray(R).copy(),
                "theta": np.asarray(theta).copy(),
                "Z": np.asarray(Z).copy(),
                "T": np.asarray(T).copy(),
                "Q": np.asarray([scaling.Q]),
            }
        )
        chart_cart = np.stack(
            (
                R * np.cos(theta),
                R * np.sin(theta),
                Z + 0.5 * T,
            ),
            axis=-1,
        )
        return _harmonic_from_chart_cartesian(chart_cart)

    return AxisSafeSourceHarmonicProvider(
        provider_id="manufactured-axis-safe-provider",
        provider_semantic_sha256="1" * 64,
        axis_zero_radius_chart=radius,
        axis_zero_certified=True,
        evaluate_chart=evaluate,
    )


def test_axis_core_is_exact_zero_and_provider_is_not_called() -> None:
    calls: list[dict[str, np.ndarray]] = []
    field = KokunoProviderDrivenAxisSafePhysicalVelocity(
        _provider(calls),
        ell=8,
        h=0.004,
    )
    core = field.axis_zero_radius_physical
    x = np.asarray((0.0, 0.5 * core, -0.75 * core))
    y = np.asarray((0.0, 0.0, 0.0))
    velocity = field.velocity(x, y, 0.0, 0.75)
    np.testing.assert_array_equal(velocity, np.zeros((3, 3)))
    assert calls == []


def test_off_axis_batch_maps_to_chart_and_scales_real_pair_velocity() -> None:
    calls: list[dict[str, np.ndarray]] = []
    field = KokunoProviderDrivenAxisSafePhysicalVelocity(
        _provider(calls),
        ell=10,
        h=0.003,
    )
    core = field.axis_zero_radius_physical
    x = np.asarray((0.0, 2.0 * core, -3.0 * core))
    y = np.asarray((0.0, 0.0, 4.0 * core))
    z = np.asarray((0.0, -0.02, 0.03))
    t = np.asarray((0.9, 0.8, 0.7))

    velocity = field.velocity(x, y, z, t)
    assert velocity.shape == (3, 3)
    np.testing.assert_array_equal(velocity[0], np.zeros(3))
    assert len(calls) == 1
    call = calls[0]
    assert call["R"].shape == (2,)
    assert np.all(call["R"] > field.axis_zero_radius_chart)

    chart_cart = np.stack(
        (
            call["R"] * np.cos(call["theta"]),
            call["R"] * np.sin(call["theta"]),
            call["Z"] + 0.5 * call["T"],
        ),
        axis=-1,
    )
    expected = field.scaling.velocity_scale * chart_cart
    np.testing.assert_allclose(velocity[1:], expected, rtol=2e-15, atol=2e-15)


def test_scalar_public_velocity_surface_returns_three_vector() -> None:
    calls: list[dict[str, np.ndarray]] = []
    field = KokunoProviderDrivenAxisSafePhysicalVelocity(
        _provider(calls, radius=0.1),
        ell=6,
        h=0.005,
    )
    x = 2.0 * field.axis_zero_radius_physical
    velocity = field.velocity(x, 0.0, 0.01, 0.85)
    assert velocity.shape == (3,)
    assert np.all(np.isfinite(velocity))
    assert len(calls) == 1
    assert calls[0]["R"].shape == (1,)


def test_provider_batch_shape_mismatch_fails_closed() -> None:
    def bad_evaluate(
        R: np.ndarray,
        theta: np.ndarray,
        Z: np.ndarray,
        T: np.ndarray,
        scaling: SourcePhysicalScaling,
    ) -> SourceLocalizedZeroDataHarmonic:
        del R, theta, Z, T, scaling
        return _harmonic_from_chart_cartesian(np.asarray([[1.0, 2.0, 3.0]]))

    provider = AxisSafeSourceHarmonicProvider(
        provider_id="bad-shape",
        provider_semantic_sha256="2" * 64,
        axis_zero_radius_chart=0.1,
        axis_zero_certified=True,
        evaluate_chart=bad_evaluate,
    )
    field = KokunoProviderDrivenAxisSafePhysicalVelocity(provider, ell=7, h=0.004)
    core = field.axis_zero_radius_physical
    with pytest.raises(ValueError, match="batch shape"):
        field.velocity(np.asarray((2.0 * core, 3.0 * core)), 0.0, 0.0, 0.8)


def test_provider_contract_and_identity_are_reported_without_source_exact_promotion() -> None:
    calls: list[dict[str, np.ndarray]] = []
    provider = _provider(calls)
    field = KokunoProviderDrivenAxisSafePhysicalVelocity(provider, ell=9, h=0.004)
    report = field.report()
    contract = report["contract"]

    assert report["provider_id"] == provider.provider_id
    assert report["provider_semantic_sha256"] == provider.provider_semantic_sha256
    assert contract["task"] == TASK == "K2-OSC-103"
    assert contract["provider_driven_velocity_xyzt_materialized"] is True
    assert contract["batch_evaluation_materialized"] is True
    assert contract["certified_axis_zero_core_enforced"] is True
    assert contract["axis_core_skips_singular_cylindrical_evaluation"] is True
    assert contract["axis_zero_core_radius_is_repository_provider_contract"] is True
    assert contract["source_exact_axis_support_radius_recovered"] is False
    assert contract["project_domain_source_input_provider_materialized"] is False
    assert contract["global_axis_safe_source_velocity_materialized"] is False
    assert contract["self_contained_velocity_xyzt_provider"] is False
    assert contract["complete_ns_residual_assessed"] is False
    assert contract["residual_reduction_claimed"] is False
    assert contract["paper_exact"] is False
    assert contract["pde_validated"] is False


def test_invalid_provider_or_inputs_fail_closed_and_no_scientific_tuning_surface() -> None:
    def evaluate(
        R: np.ndarray,
        theta: np.ndarray,
        Z: np.ndarray,
        T: np.ndarray,
        scaling: SourcePhysicalScaling,
    ) -> SourceLocalizedZeroDataHarmonic:
        del theta, Z, T, scaling
        return _harmonic_from_chart_cartesian(np.zeros((R.size, 3)))

    base = dict(
        provider_id="provider",
        provider_semantic_sha256="3" * 64,
        axis_zero_radius_chart=0.1,
        axis_zero_certified=True,
        evaluate_chart=evaluate,
    )
    for mutation in (
        {"provider_id": ""},
        {"provider_semantic_sha256": "not-a-sha"},
        {"axis_zero_radius_chart": 0.0},
        {"axis_zero_radius_chart": np.inf},
        {"axis_zero_certified": False},
        {"evaluate_chart": None},
    ):
        payload = dict(base)
        payload.update(mutation)
        provider = AxisSafeSourceHarmonicProvider(**payload)  # type: ignore[arg-type]
        with pytest.raises(ValueError):
            KokunoProviderDrivenAxisSafePhysicalVelocity(provider, ell=8, h=0.004)

    field = KokunoProviderDrivenAxisSafePhysicalVelocity(
        AxisSafeSourceHarmonicProvider(**base),
        ell=8,
        h=0.004,
    )
    with pytest.raises(ValueError):
        field.velocity(np.nan, 0.0, 0.0, 0.8)

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
    assert forbidden.isdisjoint(inspect.signature(field.velocity).parameters)
    assert forbidden.isdisjoint(
        inspect.signature(KokunoProviderDrivenAxisSafePhysicalVelocity).parameters
    )


def test_contract_direct_function_is_fail_closed() -> None:
    contract = source_axis_safe_physical_velocity_contract()
    assert contract["source_chart_to_physical_scaling_consumed"] is True
    assert contract["caller_supplies_corrected_background"] is True
    assert contract["caller_supplies_mode_forcing_directional_jet"] is True
