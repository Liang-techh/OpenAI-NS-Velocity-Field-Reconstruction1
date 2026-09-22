from __future__ import annotations

import inspect
import math
import runpy

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_source_azimuthal_frame_oriented_multiharmonic_velocity import (
    KokunoAzimuthalFrameOrientedBoundedMultiHarmonicPhysicalVelocity,
)
from openai_ns_reconstruction.kokuno_source_bounded_scale_oriented_multiharmonic_velocity import (
    KokunoBoundedIsotropicScaleOrientedMultiHarmonicPhysicalVelocity,
    SPATIAL_SCALE_MAX,
    SPATIAL_SCALE_MIN,
    bounded_isotropic_scale_contract,
    public_contract,
)


_ORIENTED_TEST = "tests/test_constrained_kokuno_source_azimuthal_frame_oriented_multiharmonic_velocity.py"


def _helpers():
    ns = runpy.run_path(_ORIENTED_TEST)
    return ns["_parent_field"], ns["_physical"]


def _oriented(*, alpha=0.43, harmonic_calls=None, pulse_calls=None):
    parent_field, _ = _helpers()
    raw = parent_field(harmonic_calls=harmonic_calls, pulse_calls=pulse_calls)
    oriented = KokunoAzimuthalFrameOrientedBoundedMultiHarmonicPhysicalVelocity(
        raw, frame_azimuth=alpha
    )
    return raw, oriented


def _physical(raw, R, Z, T, theta):
    _, physical_fn = _helpers()
    return physical_fn(raw, R, Z, T, theta)


def test_unit_scale_exactly_replays_oriented_parent():
    raw, parent = _oriented(alpha=0.37)
    field = KokunoBoundedIsotropicScaleOrientedMultiHarmonicPhysicalVelocity(
        parent, spatial_scale=1.0
    )
    x, y, z, t = _physical(
        raw,
        np.asarray([0.31, 0.47, 0.67]),
        np.asarray([-0.18, 0.03, 0.24]),
        np.asarray([0.33, 0.51, 0.67]),
        np.asarray([-0.7, 0.4, 1.6]),
    )
    got = field.evaluate(x, y, z, t)
    ref = parent.evaluate(x, y, z, t)

    assert np.array_equal(got.active_term_count, ref.active_term_count)
    assert np.allclose(
        got.complex_vector_potential_cartesian_physical,
        ref.complex_vector_potential_cartesian_physical,
        rtol=0,
        atol=2e-14,
    )
    assert np.allclose(
        got.complex_velocity_cartesian_physical,
        ref.complex_velocity_cartesian_physical,
        rtol=0,
        atol=2e-14,
    )
    assert np.allclose(
        got.real_pair_velocity_cartesian_physical,
        ref.real_pair_velocity_cartesian_physical,
        rtol=0,
        atol=2e-14,
    )


def test_nonunit_scale_matches_independent_coordinate_and_potential_dilation():
    raw, parent = _oriented(alpha=-0.52)
    scale = 1.37
    field = KokunoBoundedIsotropicScaleOrientedMultiHarmonicPhysicalVelocity(
        parent, spatial_scale=scale
    )
    xp, yp, zp, t = _physical(
        raw,
        np.asarray([0.32, 0.46, 0.63]),
        np.asarray([-0.16, 0.07, 0.21]),
        np.asarray([0.34, 0.52, 0.66]),
        np.asarray([-0.4, 0.8, 1.4]),
    )
    ref = parent.evaluate(xp, yp, zp, t)
    got = field.evaluate(scale * xp, scale * yp, scale * zp, t)

    assert np.array_equal(got.active_term_count, ref.active_term_count)
    assert np.allclose(
        got.complex_vector_potential_cartesian_physical,
        scale * ref.complex_vector_potential_cartesian_physical,
        rtol=0,
        atol=5e-14,
    )
    assert np.allclose(
        got.complex_velocity_cartesian_physical,
        ref.complex_velocity_cartesian_physical,
        rtol=0,
        atol=4e-14,
    )
    assert np.allclose(
        got.real_pair_velocity_cartesian_physical,
        ref.real_pair_velocity_cartesian_physical,
        rtol=0,
        atol=4e-14,
    )


def test_scale_is_applied_before_support_and_pulse_provider_evaluation():
    harmonic_calls = {}
    pulse_calls = {}
    raw, parent = _oriented(
        alpha=0.29,
        harmonic_calls=harmonic_calls,
        pulse_calls=pulse_calls,
    )
    scale = 1.5
    field = KokunoBoundedIsotropicScaleOrientedMultiHarmonicPhysicalVelocity(
        parent, spatial_scale=scale
    )

    # These are specified in the parent chart, then dilated into lab coordinates.
    # R=.45 is active; .05 is in the certified axis core; .95 is outside the
    # certified radial support; Z=.80 is outside z support.
    xp, yp, zp, t = _physical(
        raw,
        np.asarray([0.45, 0.05, 0.95, 0.45]),
        np.asarray([0.00, 0.00, 0.00, 0.80]),
        np.asarray([0.50, 0.50, 0.50, 0.50]),
        np.asarray([0.70, -0.20, 1.10, -1.00]),
    )
    out = field.evaluate(scale * xp, scale * yp, scale * zp, t)

    assert int(out.active_term_count[0]) == 1
    assert np.linalg.norm(out.real_pair_velocity_cartesian_physical[0]) > 0.0
    assert np.array_equal(out.active_term_count[1:], np.zeros(3, dtype=np.int16))
    assert np.array_equal(
        out.real_pair_velocity_cartesian_physical[1:], np.zeros((3, 3))
    )
    assert pulse_calls == {"calls": 1, "points": 1}
    assert harmonic_calls == {"calls": 1, "points": 1}


def test_real_pair_closure_scalar_batch_and_time_is_not_scaled():
    raw, parent = _oriented(alpha=0.61)
    scale = 0.83
    field = KokunoBoundedIsotropicScaleOrientedMultiHarmonicPhysicalVelocity(
        parent, spatial_scale=scale
    )
    xp, yp, zp, t = _physical(raw, [0.45], [0.05], [0.50], [0.8])
    x, y, z = scale * xp, scale * yp, scale * zp
    scalar = field.evaluate(float(x[0]), float(y[0]), float(z[0]), float(t[0]))
    batch = field.evaluate(x, y, z, t)
    ref = parent.evaluate(xp, yp, zp, t)

    assert scalar.real_pair_velocity_cartesian_physical.shape == (3,)
    assert batch.real_pair_velocity_cartesian_physical.shape == (1, 3)
    assert np.allclose(
        batch.real_pair_velocity_cartesian_physical,
        2.0 * np.real(batch.complex_velocity_cartesian_physical),
        rtol=0,
        atol=2e-14,
    )
    assert np.allclose(
        batch.real_pair_velocity_cartesian_physical,
        ref.real_pair_velocity_cartesian_physical,
        rtol=0,
        atol=3e-14,
    )
    assert np.allclose(
        scalar.real_pair_velocity_cartesian_physical,
        batch.real_pair_velocity_cartesian_physical[0],
        rtol=0,
        atol=2e-14,
    )


def test_scaled_public_potential_has_matching_complete_curl_by_independent_fd():
    raw, parent = _oriented(alpha=0.23)
    scale = 1.21
    field = KokunoBoundedIsotropicScaleOrientedMultiHarmonicPhysicalVelocity(
        parent, spatial_scale=scale
    )
    xp, yp, zp, t = _physical(raw, [0.46], [0.04], [0.50], [0.73])
    x0, y0, z0, t0 = [float(scale * value[0]) for value in (xp, yp, zp)] + [float(t[0])]

    # Independent Cartesian centered difference of the public complex vector
    # potential.  The step is test-only and not a public/scientific input.
    h = 2.0e-6

    def A(x, y, z):
        return np.asarray(
            field.evaluate(x, y, z, t0).complex_vector_potential_cartesian_physical,
            dtype=np.complex128,
        )

    dA_dx = (A(x0 + h, y0, z0) - A(x0 - h, y0, z0)) / (2.0 * h)
    dA_dy = (A(x0, y0 + h, z0) - A(x0, y0 - h, z0)) / (2.0 * h)
    dA_dz = (A(x0, y0, z0 + h) - A(x0, y0, z0 - h)) / (2.0 * h)
    curl = np.asarray(
        (
            dA_dy[2] - dA_dz[1],
            dA_dz[0] - dA_dx[2],
            dA_dx[1] - dA_dy[0],
        ),
        dtype=np.complex128,
    )
    u = np.asarray(
        field.evaluate(x0, y0, z0, t0).complex_velocity_cartesian_physical,
        dtype=np.complex128,
    )
    denom = max(float(np.linalg.norm(u)), 1e-12)
    assert float(np.linalg.norm(curl - u)) / denom < 2e-3


def test_scale_bounds_parent_type_and_nonfinite_inputs_fail_closed():
    _, parent = _oriented()
    for bad in [
        SPATIAL_SCALE_MIN - 1e-12,
        SPATIAL_SCALE_MAX + 1e-12,
        0.0,
        -1.0,
        math.inf,
        math.nan,
    ]:
        with pytest.raises(ValueError):
            KokunoBoundedIsotropicScaleOrientedMultiHarmonicPhysicalVelocity(
                parent, spatial_scale=bad
            )
    with pytest.raises(ValueError):
        KokunoBoundedIsotropicScaleOrientedMultiHarmonicPhysicalVelocity(
            object(), spatial_scale=1.0
        )
    field = KokunoBoundedIsotropicScaleOrientedMultiHarmonicPhysicalVelocity(
        parent, spatial_scale=SPATIAL_SCALE_MAX
    )
    with pytest.raises(ValueError):
        field.velocity(np.nan, 0.0, 0.0, 0.5)


def test_semantic_identity_binds_scale_and_exact_parent_identity():
    _, parent = _oriented(alpha=-0.31)
    f0 = KokunoBoundedIsotropicScaleOrientedMultiHarmonicPhysicalVelocity(
        parent, spatial_scale=0.9
    )
    f1 = KokunoBoundedIsotropicScaleOrientedMultiHarmonicPhysicalVelocity(
        parent, spatial_scale=1.1
    )
    assert f0.semantic_sha256 != f1.semantic_sha256
    cfg = f1.configuration()
    assert cfg["stack_parent"]["head"] == "3bdf0d6a33a629c4262875b146e1378d550023d0"
    assert cfg["oriented_parent"]["head"] == "59d3e6d64fcdfeaf145cbee67e68b8a25e89489c"
    assert cfg["oriented_parent"]["semantic_sha256"] == parent.semantic_sha256
    assert cfg["spatial_scale"]["classification"] == "repository_autonomous_candidate_design"
    assert cfg["spatial_scale"]["source_exact_claimed"] is False


def test_truth_boundary_and_public_velocity_tuning_firewall():
    contract = bounded_isotropic_scale_contract()
    public = public_contract()
    assert contract["bounded_isotropic_physical_scale_materialized"] is True
    assert contract["independent_scale_parameter_materialized"] is True
    assert contract["curl_contract_preserved_by_constant_dilation_chain_rule"] is True
    assert contract["divergence_free_contract_preserved_by_constant_dilation"] is True
    assert contract["post_velocity_hard_support_mask_used"] is False
    assert contract["source_exact_scale_recovered"] is False
    assert contract["source_exact_scale_claimed"] is False
    assert contract["self_contained_velocity_xyzt_provider"] is False
    assert contract["complete_ns_residual_assessed"] is False
    assert contract["residual_reduction_claimed"] is False
    assert contract["paper_exact"] is False
    assert contract["pde_validated"] is False
    assert public["forbidden_velocity_inputs_present"] == []
    assert set(
        inspect.signature(
            KokunoBoundedIsotropicScaleOrientedMultiHarmonicPhysicalVelocity.velocity
        ).parameters
    ) == {"self", "x", "y", "z", "t"}
