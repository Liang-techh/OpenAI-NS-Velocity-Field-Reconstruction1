from __future__ import annotations

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
)
from openai_ns_reconstruction.kokuno_source_physicalized_localized_harmonic import (
    source_chart_to_physical_rzt,
)
from openai_ns_reconstruction.kokuno_source_support_certified_multiharmonic_velocity import (
    KokunoSupportCertifiedBoundedMultiHarmonicPhysicalVelocity,
    SupportCertifiedHarmonicTerm,
    public_contract,
    support_certified_multiharmonic_contract,
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


def _flat_bump_and_derivative(s: np.ndarray, ds: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    s = np.asarray(s, dtype=float)
    ds = np.asarray(ds, dtype=float)
    value = np.zeros_like(s)
    derivative = np.zeros_like(s)
    active = s > 0.0
    value[active] = np.exp(-1.0 / s[active])
    derivative[active] = value[active] * ds[active] / (s[active] ** 2)
    return value, derivative


def _manufactured_complete_curl_harmonic(
    R: np.ndarray,
    theta: np.ndarray,
    Z: np.ndarray,
    T: np.ndarray,
    *,
    core: float,
    rmax: float,
    zmin: float,
    zmax: float,
    multiplier: float,
) -> SourceLocalizedZeroDataHarmonic:
    """C-infinity compact A_z bump and its analytic toroidal curl."""
    R = np.asarray(R, dtype=float)
    theta = np.asarray(theta, dtype=float)
    Z = np.asarray(Z, dtype=float)
    T = np.asarray(T, dtype=float)
    del T
    n = R.size

    s_r = (R - core) * (rmax - R)
    ds_r = (rmax + core) - 2.0 * R
    b_r, db_r = _flat_bump_and_derivative(s_r, ds_r)

    s_z = (Z - zmin) * (zmax - Z)
    ds_z = (zmax + zmin) - 2.0 * Z
    b_z, _ = _flat_bump_and_derivative(s_z, ds_z)

    complex_scale = multiplier * (1.0 + 0.2j)
    a_z = complex_scale * b_r * b_z
    d_r_a_z = complex_scale * db_r * b_z

    potential = np.zeros((n, 3), dtype=np.complex128)
    potential[:, 2] = a_z
    velocity = np.zeros((n, 3), dtype=np.complex128)
    velocity[:, 1] = -d_r_a_z

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


def _provider_and_term(
    *,
    provider_id: str = "m1",
    digest_char: str = "1",
    core: float = 0.10,
    rmax: float = 0.90,
    zmin: float = -0.60,
    zmax: float = 0.70,
    amplitude: float = 0.40,
    phase: float = 0.30,
    multiplier: float = 1.0,
    calls: dict | None = None,
) -> tuple[AxisSafeSourceHarmonicProvider, SupportCertifiedHarmonicTerm]:
    def evaluate(R, theta, Z, T, scaling):
        del scaling
        if calls is not None:
            calls["calls"] = calls.get("calls", 0) + 1
            calls["points"] = calls.get("points", 0) + int(np.asarray(R).size)
        return _manufactured_complete_curl_harmonic(
            R,
            theta,
            Z,
            T,
            core=core,
            rmax=rmax,
            zmin=zmin,
            zmax=zmax,
            multiplier=multiplier,
        )

    provider = AxisSafeSourceHarmonicProvider(
        provider_id=provider_id,
        provider_semantic_sha256=digest_char * 64,
        axis_zero_radius_chart=core,
        axis_zero_certified=True,
        evaluate_chart=evaluate,
    )
    bounded = BoundedHarmonicTerm(
        provider=provider,
        amplitude=amplitude,
        phase_offset=phase,
    )
    support = SupportCertifiedHarmonicTerm(
        harmonic=bounded,
        radial_support_max_chart=rmax,
        z_support_min_chart=zmin,
        z_support_max_chart=zmax,
        spatial_zero_outside_certified=True,
        smooth_zero_extension_certified=True,
    )
    return provider, support


def _physical(field, R, Z, T, theta=0.0):
    R, Z, T, theta = np.broadcast_arrays(
        np.asarray(R, dtype=float),
        np.asarray(Z, dtype=float),
        np.asarray(T, dtype=float),
        np.asarray(theta, dtype=float),
    )
    r, z, t = source_chart_to_physical_rzt(
        R, Z, T, ell=field.scaling.ell, h=field.scaling.h
    )
    return r * np.cos(theta), r * np.sin(theta), z, t


def test_certified_axis_and_outer_spatial_zero_skip_provider_entirely():
    calls: dict[str, int] = {}
    _, term = _provider_and_term(calls=calls)
    field = KokunoSupportCertifiedBoundedMultiHarmonicPhysicalVelocity(
        (term,), ell=8, h=0.004
    )

    R = np.asarray([0.05, 0.91, 0.40, 0.40])
    Z = np.asarray([0.00, 0.00, -0.61, 0.71])
    x, y, z, t = _physical(field, R, Z, np.full(4, 0.5))
    out = field.evaluate(x, y, z, t)

    assert np.array_equal(out.real_pair_velocity_cartesian_physical, np.zeros((4, 3)))
    assert np.array_equal(out.complex_velocity_cylindrical_physical, np.zeros((4, 3)))
    assert np.array_equal(out.active_term_count, np.zeros(4, dtype=np.int16))
    assert calls == {}


def test_boundaries_are_inside_certificate_but_just_outside_is_skipped():
    calls: dict[str, int] = {}
    _, term = _provider_and_term(calls=calls)
    field = KokunoSupportCertifiedBoundedMultiHarmonicPhysicalVelocity(
        (term,), ell=8, h=0.004
    )
    eps = 1.0e-7
    R = np.asarray([0.10, 0.10 + eps, 0.90, 0.90 + eps, 0.40, 0.40])
    Z = np.asarray([0.00, 0.00, 0.00, 0.00, -0.60, 0.70])
    x, y, z, t = _physical(field, R, Z, np.full(R.shape, 0.4))
    out = field.evaluate(x, y, z, t)

    assert np.array_equal(out.active_term_count, np.asarray([0, 1, 1, 0, 1, 1]))
    assert calls["calls"] == 1
    assert calls["points"] == 4
    assert np.array_equal(
        out.real_pair_velocity_cartesian_physical[[0, 2, 3, 4, 5]], np.zeros((5, 3))
    )
    assert np.linalg.norm(out.real_pair_velocity_cartesian_physical[1]) > 0.0


def test_active_interior_matches_parent_bounded_family_without_support_modification():
    calls_parent: dict[str, int] = {}
    _, term = _provider_and_term(calls=calls_parent, amplitude=0.37, phase=-0.41)
    support_field = KokunoSupportCertifiedBoundedMultiHarmonicPhysicalVelocity(
        (term,), ell=7, h=0.003
    )
    parent_field = KokunoBoundedMultiHarmonicPhysicalVelocity(
        (term.harmonic,), ell=7, h=0.003
    )

    R = np.asarray([0.25, 0.40, 0.70])
    Z = np.asarray([-0.20, 0.10, 0.30])
    T = np.asarray([0.3, 0.5, 0.7])
    theta = np.asarray([0.0, 0.7, -1.1])
    x, y, z, t = _physical(support_field, R, Z, T, theta)

    got = support_field.evaluate(x, y, z, t)
    expected = parent_field.evaluate(x, y, z, t)
    assert np.allclose(
        got.complex_vector_potential_cylindrical_physical,
        expected.complex_vector_potential_cylindrical_physical,
        rtol=2e-14,
        atol=2e-14,
    )
    assert np.allclose(
        got.complex_velocity_cylindrical_physical,
        expected.complex_velocity_cylindrical_physical,
        rtol=2e-14,
        atol=2e-14,
    )
    assert np.allclose(
        got.real_pair_velocity_cartesian_physical,
        expected.real_pair_velocity_cartesian_physical,
        rtol=2e-14,
        atol=2e-14,
    )


def test_two_support_envelopes_have_independent_active_term_counts():
    _, t1 = _provider_and_term(
        provider_id="m1",
        digest_char="2",
        core=0.10,
        rmax=0.55,
        zmin=-0.50,
        zmax=0.10,
        amplitude=0.30,
        phase=0.10,
    )
    _, t2 = _provider_and_term(
        provider_id="m2",
        digest_char="3",
        core=0.15,
        rmax=0.95,
        zmin=-0.10,
        zmax=0.60,
        amplitude=0.50,
        phase=-0.20,
        multiplier=0.7,
    )
    field = KokunoSupportCertifiedBoundedMultiHarmonicPhysicalVelocity(
        (t1, t2), ell=8, h=0.004
    )

    R = np.asarray([0.30, 0.30, 0.70, 1.00])
    Z = np.asarray([-0.30, 0.00, 0.20, 0.00])
    x, y, z, t = _physical(field, R, Z, np.full(4, 0.5))
    out = field.evaluate(x, y, z, t)
    assert np.array_equal(out.active_term_count, np.asarray([1, 2, 1, 0]))
    assert np.array_equal(out.real_pair_velocity_cartesian_physical[3], np.zeros(3))
    assert np.linalg.norm(out.real_pair_velocity_cartesian_physical[:3]) > 0.0


def test_support_certificate_is_semantic_identity_bound():
    _, term = _provider_and_term()
    field1 = KokunoSupportCertifiedBoundedMultiHarmonicPhysicalVelocity(
        (term,), ell=8, h=0.004
    )
    term2 = SupportCertifiedHarmonicTerm(
        harmonic=term.harmonic,
        radial_support_max_chart=0.91,
        z_support_min_chart=term.z_support_min_chart,
        z_support_max_chart=term.z_support_max_chart,
        spatial_zero_outside_certified=True,
        smooth_zero_extension_certified=True,
    )
    field2 = KokunoSupportCertifiedBoundedMultiHarmonicPhysicalVelocity(
        (term2,), ell=8, h=0.004
    )
    assert field1.semantic_sha256 != field2.semantic_sha256
    assert (
        field1.configuration()["terms"][0]["support_certificate"]["classification"]
        == "repository_provider_execution_metadata"
    )


@pytest.mark.parametrize(
    "kwargs",
    [
        {"radial_support_max_chart": 0.10},
        {"radial_support_max_chart": math.inf},
        {"z_support_min_chart": 0.70, "z_support_max_chart": 0.60},
        {"z_support_min_chart": -math.inf},
        {"spatial_zero_outside_certified": False},
        {"smooth_zero_extension_certified": False},
    ],
)
def test_invalid_support_certificates_fail_closed(kwargs):
    _, term = _provider_and_term()
    data = {
        "harmonic": term.harmonic,
        "radial_support_max_chart": term.radial_support_max_chart,
        "z_support_min_chart": term.z_support_min_chart,
        "z_support_max_chart": term.z_support_max_chart,
        "spatial_zero_outside_certified": True,
        "smooth_zero_extension_certified": True,
    }
    data.update(kwargs)
    bad = SupportCertifiedHarmonicTerm(**data)
    with pytest.raises(ValueError):
        KokunoSupportCertifiedBoundedMultiHarmonicPhysicalVelocity(
            (bad,), ell=8, h=0.004
        )


def test_parent_bounded_amplitude_and_duplicate_provider_guards_are_reused():
    p, term = _provider_and_term(amplitude=0.8)
    _, second = _provider_and_term(
        provider_id="m2", digest_char="4", amplitude=0.4, phase=0.0
    )
    with pytest.raises(ValueError, match="L1 amplitude"):
        KokunoSupportCertifiedBoundedMultiHarmonicPhysicalVelocity(
            (term, second), ell=8, h=0.004
        )

    duplicate = SupportCertifiedHarmonicTerm(
        harmonic=BoundedHarmonicTerm(p, amplitude=0.1, phase_offset=0.0),
        radial_support_max_chart=0.80,
        z_support_min_chart=-0.50,
        z_support_max_chart=0.50,
        spatial_zero_outside_certified=True,
        smooth_zero_extension_certified=True,
    )
    with pytest.raises(ValueError, match="duplicate provider"):
        KokunoSupportCertifiedBoundedMultiHarmonicPhysicalVelocity(
            (
                SupportCertifiedHarmonicTerm(
                    harmonic=BoundedHarmonicTerm(p, amplitude=0.1, phase_offset=0.2),
                    radial_support_max_chart=0.85,
                    z_support_min_chart=-0.55,
                    z_support_max_chart=0.55,
                    spatial_zero_outside_certified=True,
                    smooth_zero_extension_certified=True,
                ),
                duplicate,
            ),
            ell=8,
            h=0.004,
        )


def test_nonfinite_public_inputs_fail_closed():
    _, term = _provider_and_term()
    field = KokunoSupportCertifiedBoundedMultiHarmonicPhysicalVelocity(
        (term,), ell=8, h=0.004
    )
    with pytest.raises(ValueError, match="finite"):
        field.velocity(np.nan, 0.0, 0.0, 0.5)


def test_contract_keeps_support_provenance_and_pde_boundaries_explicit():
    contract = support_certified_multiharmonic_contract()
    assert contract["provider_certified_outer_spatial_support_materialized"] is True
    assert contract["outer_support_exact_zero_skips_provider_evaluation"] is True
    assert contract["support_applied_as_execution_skip_not_velocity_hard_mask"] is True
    assert contract["smooth_zero_extension_certificate_required"] is True
    assert contract["source_exact_outer_support_geometry_recovered"] is False
    assert contract["source_exact_support_numbers_claimed"] is False
    assert contract["pulse_time_support_certificate_materialized"] is False
    assert contract["self_contained_velocity_xyzt_provider"] is False
    assert contract["complete_ns_residual_assessed"] is False
    assert contract["oscillation_before_after_ns_residual_compared"] is False
    assert contract["residual_reduction_claimed"] is False
    assert contract["paper_exact"] is False
    assert contract["pde_validated"] is False

    public = public_contract()
    assert public["forbidden_velocity_inputs_present"] == []
    params = set(
        inspect.signature(
            KokunoSupportCertifiedBoundedMultiHarmonicPhysicalVelocity.velocity
        ).parameters
    )
    assert params == {"self", "x", "y", "z", "t"}
