from __future__ import annotations

import copy
import inspect
import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_source_axis_safe_physical_velocity import (
    AxisSafeSourceHarmonicProvider,
)
from openai_ns_reconstruction.kokuno_source_bounded_multiharmonic_velocity import (
    BoundedHarmonicTerm,
    KokunoBoundedMultiHarmonicPhysicalVelocity,
    MAX_L1_AMPLITUDE,
    MAX_TERMS,
    bounded_multiharmonic_contract,
    public_contract,
)
from openai_ns_reconstruction.kokuno_source_physicalized_localized_harmonic import (
    physical_to_source_chart_rzt,
    source_physical_scaling,
)
from openai_ns_reconstruction.kokuno_source_zero_data_localized_harmonic import (
    SourceLocalizedZeroDataHarmonic,
)


def _cart_from_cyl(theta: np.ndarray, cylindrical: np.ndarray) -> np.ndarray:
    theta = np.asarray(theta, dtype=float)
    cylindrical = np.asarray(cylindrical, dtype=float)
    c = np.cos(theta)
    s = np.sin(theta)
    ur, uth, uz = (cylindrical[..., i] for i in range(3))
    return np.stack((ur * c - uth * s, ur * s + uth * c, uz), axis=-1)


def _manufactured_harmonic(
    R: np.ndarray,
    theta: np.ndarray,
    Z: np.ndarray,
    T: np.ndarray,
    *,
    multiplier: float,
) -> SourceLocalizedZeroDataHarmonic:
    R = np.asarray(R, dtype=float)
    theta = np.asarray(theta, dtype=float)
    Z = np.asarray(Z, dtype=float)
    T = np.asarray(T, dtype=float)
    n = R.size
    one = np.ones(n, dtype=float)

    potential = multiplier * np.stack(
        (
            (0.20 + 0.10j) * one,
            R + 1j * Z,
            T + 0.20j * R,
        ),
        axis=-1,
    ).astype(np.complex128)
    velocity = multiplier * np.stack(
        (
            (1.00 + 0.50j) * R,
            (0.30 - 0.20j) * Z,
            (0.10 + 0.40j) * T,
        ),
        axis=-1,
    ).astype(np.complex128)
    real_cyl = 2.0 * np.real(velocity)
    real_cart = _cart_from_cyl(theta, real_cyl)
    zeros = np.zeros((n, 3), dtype=np.complex128)
    return SourceLocalizedZeroDataHarmonic(
        amplitude_sensitivity=None,
        phase_jet=None,
        coefficient_jet=None,
        complete_amplitude=None,
        localization_weight=np.ones(n),
        dr_localization_weight=np.zeros(n),
        dz_localization_weight=np.zeros(n),
        localized_coefficient=zeros.copy(),
        localized_dr_coefficient=zeros.copy(),
        localized_dz_coefficient=zeros.copy(),
        complex_vector_potential=potential,
        complex_velocity_cylindrical=velocity,
        real_pair_velocity_cylindrical=real_cyl,
        real_pair_velocity_chart_cartesian=real_cart,
    )


def _provider(
    *,
    provider_id: str,
    digest_char: str,
    core: float,
    multiplier: float,
    calls: dict[str, int] | None = None,
) -> AxisSafeSourceHarmonicProvider:
    def evaluate(R, theta, Z, T, scaling):
        del scaling
        if calls is not None:
            calls["count"] = calls.get("count", 0) + 1
        return _manufactured_harmonic(R, theta, Z, T, multiplier=multiplier)

    return AxisSafeSourceHarmonicProvider(
        provider_id=provider_id,
        provider_semantic_sha256=digest_char * 64,
        axis_zero_radius_chart=core,
        axis_zero_certified=True,
        evaluate_chart=evaluate,
    )


def _manual_expected(
    field: KokunoBoundedMultiHarmonicPhysicalVelocity,
    x,
    y,
    z,
    t,
    multipliers: dict[str, float],
):
    x, y, z, t = np.broadcast_arrays(
        np.asarray(x, dtype=float),
        np.asarray(y, dtype=float),
        np.asarray(z, dtype=float),
        np.asarray(t, dtype=float),
    )
    theta = np.arctan2(y, x)
    radius = np.hypot(x, y)
    R, Z, T = physical_to_source_chart_rzt(
        radius,
        z,
        t,
        ell=field.scaling.ell,
        h=field.scaling.h,
    )
    scale = source_physical_scaling(ell=field.scaling.ell, h=field.scaling.h)
    potential = np.zeros(R.shape + (3,), dtype=np.complex128)
    velocity = np.zeros(R.shape + (3,), dtype=np.complex128)
    active_count = np.zeros(R.shape, dtype=np.int16)

    for term in field.terms:
        active = R > term.provider.axis_zero_radius_chart
        if not np.any(active):
            continue
        base = _manufactured_harmonic(
            R[active], theta[active], Z[active], T[active],
            multiplier=multipliers[term.provider.provider_id],
        )
        coeff = term.amplitude * np.exp(1j * term.phase_offset)
        potential[active] += coeff * scale.vector_potential_scale * base.complex_vector_potential
        velocity[active] += coeff * scale.velocity_scale * base.complex_velocity_cylindrical
        active_count[active] += 1

    real_cyl = 2.0 * np.real(velocity)
    return potential, velocity, real_cyl, _cart_from_cyl(theta, real_cyl), active_count


def test_two_term_family_matches_independent_complex_sum_and_real_pair_rotation():
    p1 = _provider(provider_id="m1", digest_char="1", core=0.08, multiplier=1.0)
    p2 = _provider(provider_id="m2", digest_char="2", core=0.16, multiplier=0.6)
    field = KokunoBoundedMultiHarmonicPhysicalVelocity(
        (
            BoundedHarmonicTerm(p1, amplitude=0.25, phase_offset=0.40),
            BoundedHarmonicTerm(p2, amplitude=0.50, phase_offset=-0.70),
        ),
        ell=8,
        h=0.004,
    )

    radial = field.scaling.radial_scale
    x = np.asarray([0.0, 0.12 * radial, 0.22 * radial, 0.50 * radial])
    y = np.asarray([0.0, 0.00, 0.08 * radial, -0.14 * radial])
    z = np.asarray([0.0, -0.01, 0.02, 0.03])
    t = np.asarray([0.90, 0.83, 0.76, 0.69])

    got = field.evaluate(x, y, z, t)
    expected = _manual_expected(field, x, y, z, t, {"m1": 1.0, "m2": 0.6})
    assert np.allclose(got.complex_vector_potential_cylindrical_physical, expected[0], rtol=2e-14, atol=2e-14)
    assert np.allclose(got.complex_velocity_cylindrical_physical, expected[1], rtol=2e-14, atol=2e-14)
    assert np.allclose(got.real_pair_velocity_cylindrical_physical, expected[2], rtol=2e-14, atol=2e-14)
    assert np.allclose(got.real_pair_velocity_cartesian_physical, expected[3], rtol=2e-14, atol=2e-14)
    assert np.array_equal(got.active_term_count, expected[4])
    assert np.allclose(field.velocity(x, y, z, t), expected[3], rtol=2e-14, atol=2e-14)


def test_same_complex_coefficient_is_applied_to_potential_and_complete_curl_velocity():
    provider = _provider(provider_id="phase", digest_char="3", core=0.05, multiplier=1.0)
    amplitude = 0.73
    phase = 1.17
    field = KokunoBoundedMultiHarmonicPhysicalVelocity(
        (BoundedHarmonicTerm(provider, amplitude=amplitude, phase_offset=phase),),
        ell=7,
        h=0.003,
    )
    core = field.scaling.radial_scale * provider.axis_zero_radius_chart
    x, y, z, t = 3.0 * core, -1.5 * core, 0.012, 0.81
    out = field.evaluate(x, y, z, t)

    R, Z, T = physical_to_source_chart_rzt(
        np.hypot(x, y), z, t, ell=7, h=0.003
    )
    theta = np.arctan2(y, x)
    source = _manufactured_harmonic(
        np.atleast_1d(R), np.atleast_1d(theta), np.atleast_1d(Z), np.atleast_1d(T), multiplier=1.0
    )
    scaling = field.scaling
    c = amplitude * np.exp(1j * phase)
    expected_A = c * scaling.vector_potential_scale * source.complex_vector_potential[0]
    expected_u = c * scaling.velocity_scale * source.complex_velocity_cylindrical[0]
    assert np.allclose(out.complex_vector_potential_cylindrical_physical, expected_A)
    assert np.allclose(out.complex_velocity_cylindrical_physical, expected_u)


def test_axis_certificates_skip_each_provider_and_return_exact_zero():
    calls1, calls2 = {"count": 0}, {"count": 0}
    p1 = _provider(provider_id="axis1", digest_char="4", core=0.20, multiplier=1.0, calls=calls1)
    p2 = _provider(provider_id="axis2", digest_char="5", core=0.30, multiplier=0.5, calls=calls2)
    field = KokunoBoundedMultiHarmonicPhysicalVelocity(
        (
            BoundedHarmonicTerm(p1, 0.4, 0.0),
            BoundedHarmonicTerm(p2, 0.4, 0.2),
        ),
        ell=8,
        h=0.004,
    )
    radial = field.scaling.radial_scale
    x = np.asarray([0.0, 0.10 * radial, -0.18 * radial])
    got = field.evaluate(x, np.zeros(3), np.zeros(3), np.full(3, 0.8))
    assert np.array_equal(got.real_pair_velocity_cartesian_physical, np.zeros((3, 3)))
    assert np.array_equal(got.active_term_count, np.zeros(3, dtype=np.int16))
    assert calls1["count"] == 0
    assert calls2["count"] == 0


def test_mixed_axis_cores_activate_only_the_eligible_terms():
    calls1, calls2 = {"count": 0}, {"count": 0}
    p1 = _provider(provider_id="small-core", digest_char="6", core=0.10, multiplier=1.0, calls=calls1)
    p2 = _provider(provider_id="large-core", digest_char="7", core=0.30, multiplier=0.5, calls=calls2)
    field = KokunoBoundedMultiHarmonicPhysicalVelocity(
        (
            BoundedHarmonicTerm(p1, 0.4, 0.1),
            BoundedHarmonicTerm(p2, 0.4, -0.1),
        ),
        ell=8,
        h=0.004,
    )
    radial = field.scaling.radial_scale
    x = np.asarray([0.05, 0.20, 0.45]) * radial
    got = field.evaluate(x, np.zeros(3), np.zeros(3), np.full(3, 0.8))
    assert np.array_equal(got.active_term_count, np.asarray([0, 1, 2], dtype=np.int16))
    assert calls1["count"] == 1
    assert calls2["count"] == 1


def test_scalar_public_velocity_returns_three_vector():
    p = _provider(provider_id="scalar", digest_char="8", core=0.05, multiplier=1.0)
    field = KokunoBoundedMultiHarmonicPhysicalVelocity(
        (BoundedHarmonicTerm(p, 0.5, 0.3),), ell=8, h=0.004
    )
    core = field.scaling.radial_scale * p.axis_zero_radius_chart
    value = field.velocity(3.0 * core, core, 0.01, 0.8)
    assert value.shape == (3,)
    assert np.all(np.isfinite(value))


def test_bounds_reject_parameter_explosion_duplicate_identity_and_collapse():
    providers = [
        _provider(provider_id=f"p{i}", digest_char=hex(i + 1)[2:], core=0.05, multiplier=1.0)
        for i in range(MAX_TERMS + 1)
    ]
    with pytest.raises(ValueError, match="between 1"):
        KokunoBoundedMultiHarmonicPhysicalVelocity((), ell=8, h=0.004)
    with pytest.raises(ValueError, match="between 1"):
        KokunoBoundedMultiHarmonicPhysicalVelocity(
            tuple(BoundedHarmonicTerm(p, 0.1, 0.0) for p in providers), ell=8, h=0.004
        )
    with pytest.raises(ValueError, match="\[0,1\]"):
        KokunoBoundedMultiHarmonicPhysicalVelocity(
            (BoundedHarmonicTerm(providers[0], -0.1, 0.0),), ell=8, h=0.004
        )
    with pytest.raises(ValueError, match="\[0,1\]"):
        KokunoBoundedMultiHarmonicPhysicalVelocity(
            (BoundedHarmonicTerm(providers[0], 1.1, 0.0),), ell=8, h=0.004
        )
    with pytest.raises(ValueError, match="phase_offset"):
        KokunoBoundedMultiHarmonicPhysicalVelocity(
            (BoundedHarmonicTerm(providers[0], 0.3, math.pi + 1e-3),), ell=8, h=0.004
        )
    with pytest.raises(ValueError, match="L1 amplitude"):
        KokunoBoundedMultiHarmonicPhysicalVelocity(
            (
                BoundedHarmonicTerm(providers[0], 0.7, 0.0),
                BoundedHarmonicTerm(providers[1], 0.4, 0.0),
            ),
            ell=8,
            h=0.004,
        )
    with pytest.raises(ValueError, match="collapse"):
        KokunoBoundedMultiHarmonicPhysicalVelocity(
            (BoundedHarmonicTerm(providers[0], 0.0, 0.0),), ell=8, h=0.004
        )
    duplicate = _provider(provider_id="p0", digest_char="1", core=0.08, multiplier=0.4)
    with pytest.raises(ValueError, match="duplicate provider identity"):
        KokunoBoundedMultiHarmonicPhysicalVelocity(
            (
                BoundedHarmonicTerm(providers[0], 0.3, 0.0),
                BoundedHarmonicTerm(duplicate, 0.3, 0.2),
            ),
            ell=8,
            h=0.004,
        )


def test_parent_provider_contract_is_still_fail_closed():
    bad_sha = AxisSafeSourceHarmonicProvider(
        provider_id="bad",
        provider_semantic_sha256="not-a-sha",
        axis_zero_radius_chart=0.1,
        axis_zero_certified=True,
        evaluate_chart=lambda *args: None,
    )
    with pytest.raises(ValueError, match="64-hex"):
        KokunoBoundedMultiHarmonicPhysicalVelocity(
            (BoundedHarmonicTerm(bad_sha, 0.2, 0.0),), ell=8, h=0.004
        )
    uncertified = AxisSafeSourceHarmonicProvider(
        provider_id="uncertified",
        provider_semantic_sha256="a" * 64,
        axis_zero_radius_chart=0.1,
        axis_zero_certified=False,
        evaluate_chart=lambda *args: None,
    )
    with pytest.raises(ValueError, match="certify exact zero"):
        KokunoBoundedMultiHarmonicPhysicalVelocity(
            (BoundedHarmonicTerm(uncertified, 0.2, 0.0),), ell=8, h=0.004
        )


def test_configuration_semantic_identity_binds_amplitude_phase_and_provider_identity():
    p1 = _provider(provider_id="identity", digest_char="b", core=0.05, multiplier=1.0)
    a = KokunoBoundedMultiHarmonicPhysicalVelocity(
        (BoundedHarmonicTerm(p1, 0.5, 0.1),), ell=8, h=0.004
    )
    b = KokunoBoundedMultiHarmonicPhysicalVelocity(
        (BoundedHarmonicTerm(p1, 0.5, 0.1),), ell=8, h=0.004
    )
    c = KokunoBoundedMultiHarmonicPhysicalVelocity(
        (BoundedHarmonicTerm(p1, 0.5, 0.2),), ell=8, h=0.004
    )
    p2 = _provider(provider_id="identity-2", digest_char="c", core=0.05, multiplier=1.0)
    d = KokunoBoundedMultiHarmonicPhysicalVelocity(
        (BoundedHarmonicTerm(p2, 0.5, 0.1),), ell=8, h=0.004
    )
    assert a.semantic_sha256 == b.semantic_sha256
    assert a.semantic_sha256 != c.semantic_sha256
    assert a.semantic_sha256 != d.semantic_sha256
    assert len(a.semantic_sha256) == 64


def test_nonfinite_inputs_fail_closed():
    p = _provider(provider_id="finite", digest_char="d", core=0.05, multiplier=1.0)
    field = KokunoBoundedMultiHarmonicPhysicalVelocity(
        (BoundedHarmonicTerm(p, 0.5, 0.0),), ell=8, h=0.004
    )
    with pytest.raises(ValueError, match="x must be finite"):
        field.velocity(np.nan, 0.0, 0.0, 0.8)
    with pytest.raises(ValueError, match="t must be finite"):
        field.velocity(0.1, 0.0, 0.0, np.inf)


def test_contract_separates_source_structure_from_autonomous_parameters_and_pde_truth():
    contract = bounded_multiharmonic_contract()
    assert contract["finite_complete_curl_harmonic_family_materialized"] is True
    assert contract["same_coefficient_applied_to_vector_potential_and_curl_velocity"] is True
    assert contract["real_conjugate_pair_sum_materialized"] is True
    assert contract["linearity_preserves_complete_curl_divergence_free_contract"] is True
    assert contract["bounded_term_count"] is True
    assert contract["bounded_amplitude_l1"] is True
    assert contract["autonomous_amplitude_phase_parameters"] is True
    assert contract["source_exact_amplitude_phase_recovered"] is False
    assert contract["project_domain_source_input_provider_materialized"] is False
    assert contract["self_contained_velocity_xyzt_provider"] is False
    assert contract["complete_ns_residual_assessed"] is False
    assert contract["residual_reduction_claimed"] is False
    assert contract["paper_exact"] is False
    assert contract["pde_validated"] is False

    public = public_contract()
    assert public["max_terms"] == MAX_TERMS
    assert public["max_l1_amplitude"] == MAX_L1_AMPLITUDE
    assert public["forbidden_velocity_inputs_present"] == []
    assert set(public["velocity_inputs"]) == {"self", "x", "y", "z", "t"}


def test_report_contains_no_mutable_scientific_target():
    p = _provider(provider_id="report", digest_char="e", core=0.05, multiplier=1.0)
    field = KokunoBoundedMultiHarmonicPhysicalVelocity(
        (BoundedHarmonicTerm(p, 0.5, 0.0),), ell=8, h=0.004
    )
    report = field.report()
    text = repr(report).lower()
    for forbidden in ("residual_target", "forcing_gain", "pressure_fit", "optimizer_step"):
        assert forbidden not in text
    assert report["term_count"] == 1
    assert report["amplitude_l1"] == 0.5
