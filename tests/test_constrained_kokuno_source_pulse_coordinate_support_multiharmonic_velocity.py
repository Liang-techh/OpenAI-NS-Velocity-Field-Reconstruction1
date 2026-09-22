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
)
from openai_ns_reconstruction.kokuno_source_physicalized_localized_harmonic import (
    source_chart_to_physical_rzt,
)
from openai_ns_reconstruction.kokuno_source_pulse_coordinate_support_multiharmonic_velocity import (
    KokunoPulseCoordinateSupportBoundedMultiHarmonicPhysicalVelocity,
    PulseCoordinateSupportCertifiedHarmonicTerm,
    public_contract,
    pulse_coordinate_support_multiharmonic_contract,
)
from openai_ns_reconstruction.kokuno_source_support_certified_multiharmonic_velocity import (
    KokunoSupportCertifiedBoundedMultiHarmonicPhysicalVelocity,
    SupportCertifiedHarmonicTerm,
)
from openai_ns_reconstruction.kokuno_source_zero_data_localized_harmonic import (
    SourceLocalizedZeroDataHarmonic,
)


def _flat_bump_and_derivative(s: np.ndarray, ds: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    s = np.asarray(s, dtype=float)
    ds = np.asarray(ds, dtype=float)
    value = np.zeros_like(s)
    derivative = np.zeros_like(s)
    active = s > 0.0
    value[active] = np.exp(-1.0 / s[active])
    derivative[active] = value[active] * ds[active] / (s[active] ** 2)
    return value, derivative


def _cart_from_cyl(theta: np.ndarray, cylindrical: np.ndarray) -> np.ndarray:
    c = np.cos(theta)
    s = np.sin(theta)
    ur, uth, uz = (cylindrical[..., i] for i in range(3))
    return np.stack((ur * c - uth * s, ur * s + uth * c, uz), axis=-1)


def _manufactured_harmonic(R, theta, Z, T, *, core, rmax, zmin, zmax, tmin, tmax):
    """Smooth compact A_z(R,Z,T) with analytic toroidal spatial curl."""
    R = np.asarray(R, dtype=float)
    theta = np.asarray(theta, dtype=float)
    Z = np.asarray(Z, dtype=float)
    T = np.asarray(T, dtype=float)
    n = R.size

    sr = (R - core) * (rmax - R)
    dsr = (rmax + core) - 2.0 * R
    br, dbr = _flat_bump_and_derivative(sr, dsr)
    sz = (Z - zmin) * (zmax - Z)
    bz, _ = _flat_bump_and_derivative(sz, (zmax + zmin) - 2.0 * Z)
    st = (T - tmin) * (tmax - T)
    bt, _ = _flat_bump_and_derivative(st, (tmax + tmin) - 2.0 * T)

    scale = 1.0 + 0.2j
    az = scale * br * bz * bt
    dr_az = scale * dbr * bz * bt
    potential = np.zeros((n, 3), dtype=np.complex128)
    potential[:, 2] = az
    velocity = np.zeros((n, 3), dtype=np.complex128)
    velocity[:, 1] = -dr_az
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


def _term(
    *,
    harmonic_id="m1",
    harmonic_sha_char="1",
    pulse_id="pulse-map-1",
    pulse_sha_char="a",
    core=0.10,
    rmax=0.90,
    zmin=-0.60,
    zmax=0.70,
    tmin=0.20,
    tmax=0.80,
    pulse_offset=2.0,
    amplitude=0.40,
    phase=0.30,
    harmonic_calls=None,
    pulse_calls=None,
):
    def evaluate_harmonic(R, theta, Z, T, scaling):
        del scaling
        if harmonic_calls is not None:
            harmonic_calls["calls"] = harmonic_calls.get("calls", 0) + 1
            harmonic_calls["points"] = harmonic_calls.get("points", 0) + int(np.asarray(R).size)
        return _manufactured_harmonic(
            R, theta, Z, T,
            core=core, rmax=rmax, zmin=zmin, zmax=zmax, tmin=tmin, tmax=tmax,
        )

    provider = AxisSafeSourceHarmonicProvider(
        provider_id=harmonic_id,
        provider_semantic_sha256=harmonic_sha_char * 64,
        axis_zero_radius_chart=core,
        axis_zero_certified=True,
        evaluate_chart=evaluate_harmonic,
    )
    bounded = BoundedHarmonicTerm(provider=provider, amplitude=amplitude, phase_offset=phase)
    spatial = SupportCertifiedHarmonicTerm(
        harmonic=bounded,
        radial_support_max_chart=rmax,
        z_support_min_chart=zmin,
        z_support_max_chart=zmax,
        spatial_zero_outside_certified=True,
        smooth_zero_extension_certified=True,
    )

    def evaluate_v(R, theta, Z, T, scaling):
        del R, theta, Z, scaling
        if pulse_calls is not None:
            pulse_calls["calls"] = pulse_calls.get("calls", 0) + 1
            pulse_calls["points"] = pulse_calls.get("points", 0) + int(np.asarray(T).size)
        return np.asarray(T, dtype=float) + pulse_offset

    pulse = PulseCoordinateSupportCertifiedHarmonicTerm(
        spatial=spatial,
        pulse_coordinate_provider_id=pulse_id,
        pulse_coordinate_provider_semantic_sha256=pulse_sha_char * 64,
        pulse_v_min=tmin + pulse_offset,
        pulse_v_max=tmax + pulse_offset,
        pulse_zero_outside_certified=True,
        smooth_pulse_zero_extension_certified=True,
        evaluate_pulse_coordinate=evaluate_v,
    )
    return spatial, pulse


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


def test_pulse_coordinate_exterior_is_exact_zero_and_skips_harmonic_provider():
    harmonic_calls = {}
    pulse_calls = {}
    _, term = _term(harmonic_calls=harmonic_calls, pulse_calls=pulse_calls)
    field = KokunoPulseCoordinateSupportBoundedMultiHarmonicPhysicalVelocity(
        (term,), ell=8, h=0.004
    )
    R = np.full(3, 0.45)
    Z = np.zeros(3)
    T = np.asarray([0.10, 0.50, 0.90])
    x, y, z, t = _physical(field, R, Z, T)
    out = field.evaluate(x, y, z, t)

    assert np.array_equal(out.active_term_count, np.asarray([0, 1, 0]))
    assert np.array_equal(out.real_pair_velocity_cartesian_physical[[0, 2]], np.zeros((2, 3)))
    assert np.linalg.norm(out.real_pair_velocity_cartesian_physical[1]) > 0.0
    assert pulse_calls == {"calls": 1, "points": 3}
    assert harmonic_calls == {"calls": 1, "points": 1}


def test_source_pulse_coordinate_is_not_hardcoded_to_chart_T():
    _, term = _term(pulse_offset=2.0)
    field = KokunoPulseCoordinateSupportBoundedMultiHarmonicPhysicalVelocity(
        (term,), ell=8, h=0.004
    )
    x, y, z, t = _physical(field, [0.45], [0.0], [0.50])
    out = field.evaluate(x, y, z, t)
    assert int(out.active_term_count[0]) == 1
    cfg = field.configuration()["pulse_terms"][0]["pulse_support_certificate"]
    assert cfg["v_min"] == pytest.approx(2.20)
    assert cfg["v_max"] == pytest.approx(2.80)
    assert cfg["source_exact_coordinate_map_claimed"] is False


def test_active_pulse_interior_matches_exact_spatial_parent_algebra():
    spatial, term = _term(amplitude=0.37, phase=-0.41)
    field = KokunoPulseCoordinateSupportBoundedMultiHarmonicPhysicalVelocity(
        (term,), ell=7, h=0.003
    )
    parent = KokunoSupportCertifiedBoundedMultiHarmonicPhysicalVelocity(
        (spatial,), ell=7, h=0.003
    )
    R = np.asarray([0.30, 0.45, 0.70])
    Z = np.asarray([-0.20, 0.10, 0.30])
    T = np.asarray([0.30, 0.50, 0.70])
    theta = np.asarray([0.1, 0.7, -1.1])
    x, y, z, t = _physical(field, R, Z, T, theta)
    got = field.evaluate(x, y, z, t)
    expected = parent.evaluate(x, y, z, t)
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


def test_spatial_exterior_skips_even_pulse_coordinate_provider():
    harmonic_calls = {}
    pulse_calls = {}
    _, term = _term(harmonic_calls=harmonic_calls, pulse_calls=pulse_calls)
    field = KokunoPulseCoordinateSupportBoundedMultiHarmonicPhysicalVelocity(
        (term,), ell=8, h=0.004
    )
    R = np.asarray([0.05, 0.95, 0.45, 0.45])
    Z = np.asarray([0.0, 0.0, -0.70, 0.80])
    x, y, z, t = _physical(field, R, Z, np.full(4, 0.5))
    out = field.evaluate(x, y, z, t)
    assert np.array_equal(out.active_term_count, np.zeros(4, dtype=np.int16))
    assert np.array_equal(out.real_pair_velocity_cartesian_physical, np.zeros((4, 3)))
    assert pulse_calls == {}
    assert harmonic_calls == {}


def test_two_terms_keep_independent_pulse_support_and_active_counts():
    _, p1 = _term(
        harmonic_id="m1", harmonic_sha_char="2", pulse_id="v1", pulse_sha_char="b",
        tmin=0.20, tmax=0.55, pulse_offset=1.0, amplitude=0.30, phase=0.1,
    )
    _, p2 = _term(
        harmonic_id="m2", harmonic_sha_char="3", pulse_id="v2", pulse_sha_char="c",
        tmin=0.45, tmax=0.80, pulse_offset=-1.0, amplitude=0.50, phase=-0.2,
    )
    field = KokunoPulseCoordinateSupportBoundedMultiHarmonicPhysicalVelocity(
        (p1, p2), ell=8, h=0.004
    )
    T = np.asarray([0.30, 0.50, 0.70, 0.90])
    x, y, z, t = _physical(field, np.full(4, 0.45), np.zeros(4), T)
    out = field.evaluate(x, y, z, t)
    assert np.array_equal(out.active_term_count, np.asarray([1, 2, 1, 0]))
    assert np.array_equal(out.real_pair_velocity_cartesian_physical[3], np.zeros(3))


def test_pulse_certificate_and_coordinate_provider_are_semantic_identity_bound():
    _, term = _term()
    field1 = KokunoPulseCoordinateSupportBoundedMultiHarmonicPhysicalVelocity(
        (term,), ell=8, h=0.004
    )
    term2 = PulseCoordinateSupportCertifiedHarmonicTerm(
        spatial=term.spatial,
        pulse_coordinate_provider_id=term.pulse_coordinate_provider_id,
        pulse_coordinate_provider_semantic_sha256="d" * 64,
        pulse_v_min=term.pulse_v_min,
        pulse_v_max=term.pulse_v_max,
        pulse_zero_outside_certified=True,
        smooth_pulse_zero_extension_certified=True,
        evaluate_pulse_coordinate=term.evaluate_pulse_coordinate,
    )
    field2 = KokunoPulseCoordinateSupportBoundedMultiHarmonicPhysicalVelocity(
        (term2,), ell=8, h=0.004
    )
    assert field1.semantic_sha256 != field2.semantic_sha256

    term3 = PulseCoordinateSupportCertifiedHarmonicTerm(
        spatial=term.spatial,
        pulse_coordinate_provider_id=term.pulse_coordinate_provider_id,
        pulse_coordinate_provider_semantic_sha256=term.pulse_coordinate_provider_semantic_sha256,
        pulse_v_min=term.pulse_v_min - 0.01,
        pulse_v_max=term.pulse_v_max,
        pulse_zero_outside_certified=True,
        smooth_pulse_zero_extension_certified=True,
        evaluate_pulse_coordinate=term.evaluate_pulse_coordinate,
    )
    field3 = KokunoPulseCoordinateSupportBoundedMultiHarmonicPhysicalVelocity(
        (term3,), ell=8, h=0.004
    )
    assert field1.semantic_sha256 != field3.semantic_sha256


@pytest.mark.parametrize(
    "changes,match",
    [
        ({"pulse_coordinate_provider_id": ""}, "nonempty"),
        ({"pulse_coordinate_provider_semantic_sha256": "BAD"}, "64-hex"),
        ({"pulse_v_min": 1.0, "pulse_v_max": 1.0}, "v_min < v_max"),
        ({"pulse_v_min": -math.inf}, "v_min < v_max"),
        ({"pulse_zero_outside_certified": False}, "exact zero"),
        ({"smooth_pulse_zero_extension_certified": False}, "smooth zero"),
        ({"evaluate_pulse_coordinate": None}, "callable"),
    ],
)
def test_invalid_pulse_certificates_fail_closed(changes, match):
    _, term = _term()
    data = {
        "spatial": term.spatial,
        "pulse_coordinate_provider_id": term.pulse_coordinate_provider_id,
        "pulse_coordinate_provider_semantic_sha256": term.pulse_coordinate_provider_semantic_sha256,
        "pulse_v_min": term.pulse_v_min,
        "pulse_v_max": term.pulse_v_max,
        "pulse_zero_outside_certified": True,
        "smooth_pulse_zero_extension_certified": True,
        "evaluate_pulse_coordinate": term.evaluate_pulse_coordinate,
    }
    data.update(changes)
    bad = PulseCoordinateSupportCertifiedHarmonicTerm(**data)
    with pytest.raises(ValueError, match=match):
        KokunoPulseCoordinateSupportBoundedMultiHarmonicPhysicalVelocity(
            (bad,), ell=8, h=0.004
        )


def test_bad_pulse_provider_shape_and_nonfinite_values_fail_closed():
    _, term = _term()

    def bad_shape(R, theta, Z, T, scaling):
        del R, theta, Z, T, scaling
        return np.asarray([0.4, 0.5])

    bad = PulseCoordinateSupportCertifiedHarmonicTerm(
        spatial=term.spatial,
        pulse_coordinate_provider_id="bad-shape",
        pulse_coordinate_provider_semantic_sha256="e" * 64,
        pulse_v_min=0.0,
        pulse_v_max=1.0,
        pulse_zero_outside_certified=True,
        smooth_pulse_zero_extension_certified=True,
        evaluate_pulse_coordinate=bad_shape,
    )
    field = KokunoPulseCoordinateSupportBoundedMultiHarmonicPhysicalVelocity(
        (bad,), ell=8, h=0.004
    )
    x, y, z, t = _physical(field, [0.45], [0.0], [0.5])
    with pytest.raises(ValueError, match="batch shape"):
        field.velocity(x, y, z, t)

    def bad_finite(R, theta, Z, T, scaling):
        del theta, Z, T, scaling
        return np.full(np.asarray(R).shape, np.nan)

    bad2 = PulseCoordinateSupportCertifiedHarmonicTerm(
        spatial=term.spatial,
        pulse_coordinate_provider_id="bad-finite",
        pulse_coordinate_provider_semantic_sha256="f" * 64,
        pulse_v_min=0.0,
        pulse_v_max=1.0,
        pulse_zero_outside_certified=True,
        smooth_pulse_zero_extension_certified=True,
        evaluate_pulse_coordinate=bad_finite,
    )
    field2 = KokunoPulseCoordinateSupportBoundedMultiHarmonicPhysicalVelocity(
        (bad2,), ell=8, h=0.004
    )
    with pytest.raises(ValueError, match="finite"):
        field2.velocity(x, y, z, t)


def test_parent_bounded_amplitude_and_duplicate_harmonic_provider_guards_remain_active():
    _, p1 = _term(amplitude=0.8)
    _, p2 = _term(
        harmonic_id="m2", harmonic_sha_char="4", pulse_id="v2", pulse_sha_char="4",
        amplitude=0.4,
    )
    with pytest.raises(ValueError, match="L1 amplitude"):
        KokunoPulseCoordinateSupportBoundedMultiHarmonicPhysicalVelocity(
            (p1, p2), ell=8, h=0.004
        )

    _, duplicate = _term(
        harmonic_id="m1", harmonic_sha_char="1", pulse_id="v-other", pulse_sha_char="5",
        amplitude=0.1,
    )
    with pytest.raises(ValueError, match="duplicate provider"):
        KokunoPulseCoordinateSupportBoundedMultiHarmonicPhysicalVelocity(
            (p1, duplicate), ell=8, h=0.004
        )


def test_contract_keeps_source_and_physical_time_claims_fail_closed():
    contract = pulse_coordinate_support_multiharmonic_contract()
    assert contract["provider_certified_pulse_coordinate_support_materialized"] is True
    assert contract["pulse_coordinate_exact_zero_skips_harmonic_provider_evaluation"] is True
    assert contract["pulse_coordinate_support_applied_as_execution_skip_not_velocity_hard_mask"] is True
    assert contract["source_exact_pulse_coordinate_mapping_recovered"] is False
    assert contract["source_exact_pulse_support_interval_recovered"] is False
    assert contract["source_pulse_coordinate_identified_with_chart_T"] is False
    assert contract["source_pulse_coordinate_identified_with_physical_t"] is False
    assert contract["physical_time_support_certificate_materialized"] is False
    assert contract["self_contained_velocity_xyzt_provider"] is False
    assert contract["complete_ns_residual_assessed"] is False
    assert contract["paper_exact"] is False
    assert contract["pde_validated"] is False

    public = public_contract()
    assert public["forbidden_velocity_inputs_present"] == []
    params = set(
        inspect.signature(
            KokunoPulseCoordinateSupportBoundedMultiHarmonicPhysicalVelocity.velocity
        ).parameters
    )
    assert params == {"self", "x", "y", "z", "t"}
