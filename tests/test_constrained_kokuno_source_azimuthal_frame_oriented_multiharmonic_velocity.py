from __future__ import annotations

from dataclasses import replace
import inspect
import math
import runpy

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_source_azimuthal_frame_oriented_multiharmonic_velocity import (
    KokunoAzimuthalFrameOrientedBoundedMultiHarmonicPhysicalVelocity,
    azimuthal_frame_orientation_contract,
    public_contract,
)
from openai_ns_reconstruction.kokuno_source_pulse_coordinate_support_multiharmonic_velocity import (
    KokunoPulseCoordinateSupportBoundedMultiHarmonicPhysicalVelocity,
)


_PARENT_TEST = "tests/test_constrained_kokuno_source_pulse_coordinate_support_multiharmonic_velocity.py"


def _parent_helpers():
    ns = runpy.run_path(_PARENT_TEST)
    return ns["_term"], ns["_physical"]


def _parent_field(*, harmonic_calls=None, pulse_calls=None):
    term_fn, _ = _parent_helpers()
    _, term = term_fn(harmonic_calls=harmonic_calls, pulse_calls=pulse_calls)
    return KokunoPulseCoordinateSupportBoundedMultiHarmonicPhysicalVelocity(
        (term,), ell=8, h=0.004
    )


def _physical(parent, R, Z, T, theta):
    _, physical_fn = _parent_helpers()
    return physical_fn(parent, R, Z, T, theta)


def _cyl_to_cart(theta, vector):
    theta = np.asarray(theta, dtype=float)
    vector = np.asarray(vector)
    c = np.cos(theta)
    s = np.sin(theta)
    vr = vector[..., 0]
    vt = vector[..., 1]
    vz = vector[..., 2]
    return np.stack((vr * c - vt * s, vr * s + vt * c, vz), axis=-1)


def _rotate_xy(vector, alpha):
    vector = np.asarray(vector)
    c = math.cos(alpha)
    s = math.sin(alpha)
    return np.stack(
        (
            c * vector[..., 0] - s * vector[..., 1],
            s * vector[..., 0] + c * vector[..., 1],
            vector[..., 2],
        ),
        axis=-1,
    )


def _manual_parent_coordinates(x, y, alpha):
    c = math.cos(alpha)
    s = math.sin(alpha)
    return c * x + s * y, -s * x + c * y


def test_zero_orientation_is_exact_parent_field_in_cartesian_basis():
    parent = _parent_field()
    field = KokunoAzimuthalFrameOrientedBoundedMultiHarmonicPhysicalVelocity(
        parent, frame_azimuth=0.0
    )
    R = np.asarray([0.30, 0.45, 0.70])
    Z = np.asarray([-0.20, 0.10, 0.30])
    T = np.asarray([0.30, 0.50, 0.70])
    theta = np.asarray([0.1, 0.7, -1.1])
    x, y, z, t = _physical(parent, R, Z, T, theta)

    got = field.evaluate(x, y, z, t)
    ref = parent.evaluate(x, y, z, t)
    ref_A = _cyl_to_cart(theta, ref.complex_vector_potential_cylindrical_physical)
    ref_u = _cyl_to_cart(theta, ref.complex_velocity_cylindrical_physical)

    assert np.allclose(got.complex_vector_potential_cartesian_physical, ref_A, rtol=0, atol=2e-14)
    assert np.allclose(got.complex_velocity_cartesian_physical, ref_u, rtol=0, atol=2e-14)
    assert np.allclose(
        got.real_pair_velocity_cartesian_physical,
        ref.real_pair_velocity_cartesian_physical,
        rtol=0,
        atol=2e-14,
    )
    assert np.array_equal(got.active_term_count, ref.active_term_count)


def test_nonzero_orientation_matches_independent_rigid_coordinate_and_vector_rotation():
    parent = _parent_field()
    alpha = 0.63
    field = KokunoAzimuthalFrameOrientedBoundedMultiHarmonicPhysicalVelocity(
        parent, frame_azimuth=alpha
    )
    R = np.asarray([0.31, 0.48, 0.68])
    Z = np.asarray([-0.19, 0.08, 0.28])
    T = np.asarray([0.31, 0.52, 0.68])
    theta_lab = np.asarray([-0.4, 0.8, 1.7])
    x, y, z, t = _physical(parent, R, Z, T, theta_lab)

    got = field.evaluate(x, y, z, t)
    xp, yp = _manual_parent_coordinates(x, y, alpha)
    theta_parent = np.arctan2(yp, xp)
    raw = parent.evaluate(xp, yp, z, t)
    manual_A = _rotate_xy(
        _cyl_to_cart(theta_parent, raw.complex_vector_potential_cylindrical_physical),
        alpha,
    )
    manual_u = _rotate_xy(
        _cyl_to_cart(theta_parent, raw.complex_velocity_cylindrical_physical),
        alpha,
    )
    manual_real = _rotate_xy(raw.real_pair_velocity_cartesian_physical, alpha)

    assert np.allclose(got.complex_vector_potential_cartesian_physical, manual_A, rtol=0, atol=3e-14)
    assert np.allclose(got.complex_velocity_cartesian_physical, manual_u, rtol=0, atol=3e-14)
    assert np.allclose(got.real_pair_velocity_cartesian_physical, manual_real, rtol=0, atol=3e-14)
    assert np.allclose(got.real_pair_velocity_cartesian_physical, 2.0 * np.real(manual_u), rtol=0, atol=3e-14)
    assert np.array_equal(got.active_term_count, raw.active_term_count)


def test_orientation_is_applied_before_theta_sensitive_pulse_coordinate_provider():
    term_fn, physical_fn = _parent_helpers()
    _, base = term_fn(pulse_offset=0.0, pulse_sha_char="d")

    def pulse_theta(R, theta, Z, T, scaling):
        del R, Z, T, scaling
        return np.asarray(theta, dtype=float)

    term = replace(
        base,
        pulse_coordinate_provider_id="pulse-theta-map",
        pulse_coordinate_provider_semantic_sha256="e" * 64,
        pulse_v_min=0.40,
        pulse_v_max=0.60,
        evaluate_pulse_coordinate=pulse_theta,
    )
    parent = KokunoPulseCoordinateSupportBoundedMultiHarmonicPhysicalVelocity(
        (term,), ell=8, h=0.004
    )
    x, y, z, t = physical_fn(parent, [0.45], [0.0], [0.50], [0.50])

    unrotated = KokunoAzimuthalFrameOrientedBoundedMultiHarmonicPhysicalVelocity(
        parent, frame_azimuth=0.0
    ).evaluate(x, y, z, t)
    rotated = KokunoAzimuthalFrameOrientedBoundedMultiHarmonicPhysicalVelocity(
        parent, frame_azimuth=0.40
    ).evaluate(x, y, z, t)

    # Lab theta=.50 is active at alpha=0.  In the rotated parent frame theta=.10,
    # so the certified provider says the localized harmonic is exactly zero.
    assert int(unrotated.active_term_count[0]) == 1
    assert np.linalg.norm(unrotated.real_pair_velocity_cartesian_physical[0]) > 0.0
    assert int(rotated.active_term_count[0]) == 0
    assert np.array_equal(rotated.real_pair_velocity_cartesian_physical, np.zeros((1, 3)))


def test_axis_and_outer_spatial_support_remain_exact_zero_without_provider_calls():
    harmonic_calls = {}
    pulse_calls = {}
    parent = _parent_field(harmonic_calls=harmonic_calls, pulse_calls=pulse_calls)
    field = KokunoAzimuthalFrameOrientedBoundedMultiHarmonicPhysicalVelocity(
        parent, frame_azimuth=-1.17
    )
    R = np.asarray([0.05, 0.95, 0.45, 0.45])
    Z = np.asarray([0.0, 0.0, -0.70, 0.80])
    T = np.full(4, 0.5)
    theta = np.asarray([0.2, -1.3, 2.1, -2.4])
    x, y, z, t = _physical(parent, R, Z, T, theta)
    out = field.evaluate(x, y, z, t)

    assert np.array_equal(out.active_term_count, np.zeros(4, dtype=np.int16))
    assert np.array_equal(out.real_pair_velocity_cartesian_physical, np.zeros((4, 3)))
    assert pulse_calls == {}
    assert harmonic_calls == {}


def test_scalar_batch_finite_and_real_pair_closure():
    parent = _parent_field()
    field = KokunoAzimuthalFrameOrientedBoundedMultiHarmonicPhysicalVelocity(
        parent, frame_azimuth=0.27
    )
    x, y, z, t = _physical(parent, [0.45], [0.05], [0.50], [0.8])
    scalar = field.evaluate(float(x[0]), float(y[0]), float(z[0]), float(t[0]))
    batch = field.evaluate(x, y, z, t)

    assert scalar.real_pair_velocity_cartesian_physical.shape == (3,)
    assert batch.real_pair_velocity_cartesian_physical.shape == (1, 3)
    assert np.all(np.isfinite(batch.real_pair_velocity_cartesian_physical))
    assert np.allclose(
        batch.real_pair_velocity_cartesian_physical,
        2.0 * np.real(batch.complex_velocity_cartesian_physical),
        rtol=0,
        atol=2e-14,
    )
    assert np.allclose(scalar.real_pair_velocity_cartesian_physical, batch.real_pair_velocity_cartesian_physical[0])


def test_orientation_bound_parent_type_and_nonfinite_inputs_fail_closed():
    parent = _parent_field()
    for bad in [math.pi + 1e-12, -math.pi - 1e-12, math.inf, math.nan]:
        with pytest.raises(ValueError):
            KokunoAzimuthalFrameOrientedBoundedMultiHarmonicPhysicalVelocity(
                parent, frame_azimuth=bad
            )
    with pytest.raises(ValueError):
        KokunoAzimuthalFrameOrientedBoundedMultiHarmonicPhysicalVelocity(
            object(), frame_azimuth=0.0
        )
    field = KokunoAzimuthalFrameOrientedBoundedMultiHarmonicPhysicalVelocity(
        parent, frame_azimuth=math.pi
    )
    with pytest.raises(ValueError):
        field.velocity(np.nan, 0.0, 0.0, 0.5)


def test_semantic_identity_binds_orientation_and_exact_parent_identity():
    parent = _parent_field()
    f0 = KokunoAzimuthalFrameOrientedBoundedMultiHarmonicPhysicalVelocity(
        parent, frame_azimuth=0.0
    )
    f1 = KokunoAzimuthalFrameOrientedBoundedMultiHarmonicPhysicalVelocity(
        parent, frame_azimuth=0.25
    )
    assert f0.semantic_sha256 != f1.semantic_sha256
    cfg = f1.configuration()
    assert cfg["parent"]["head"] == "0cbf74972af6589f7eb9a2edfe21c957d7512a5c"
    assert cfg["parent"]["semantic_sha256"] == parent.semantic_sha256
    assert cfg["frame_orientation"]["classification"] == "repository_autonomous_candidate_design"
    assert cfg["frame_orientation"]["source_exact_claimed"] is False


def test_truth_boundary_and_public_velocity_tuning_firewall():
    c = azimuthal_frame_orientation_contract()
    p = public_contract()
    assert c["bounded_azimuthal_frame_orientation_materialized"] is True
    assert c["same_rigid_rotation_applied_to_vector_potential_and_complete_curl"] is True
    assert c["rigid_rotation_preserves_divergence_free_contract_by_covariance"] is True
    assert c["pulse_coordinate_provider_evaluated_in_rotated_parent_frame"] is True
    assert c["independent_cartesian_component_fit_used"] is False
    assert c["independent_scale_parameter_materialized"] is False
    assert c["source_exact_frame_vectors_recovered"] is False
    assert c["source_exact_orientation_recovered"] is False
    assert c["self_contained_velocity_xyzt_provider"] is False
    assert c["complete_ns_residual_assessed"] is False
    assert c["paper_exact"] is False
    assert c["pde_validated"] is False
    assert p["forbidden_velocity_inputs_present"] == []
    assert set(inspect.signature(
        KokunoAzimuthalFrameOrientedBoundedMultiHarmonicPhysicalVelocity.velocity
    ).parameters) == {"self", "x", "y", "z", "t"}
