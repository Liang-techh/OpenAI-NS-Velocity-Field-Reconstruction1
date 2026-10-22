import numpy as np
import pytest

from openai_ns_reconstruction.constrained_swirling_strength import (
    audit_swirling_strength_resolution,
    diagnose_swirling_strength,
)


def _rigid_rotation(points, time):
    omega = 2.5 + 0.0 * time
    x, y, z = points.T
    return np.column_stack((-omega * y, omega * x, 0.25 * z))


def _pure_shear(points, time):
    gamma = 3.0 + 0.0 * time
    _, y, _ = points.T
    return np.column_stack((gamma * y, np.zeros_like(y), np.zeros_like(y)))


def test_rigid_rotation_recovers_lambda_ci_and_centered_geometry():
    nine = 9
    result = diagnose_swirling_strength(
        _rigid_rotation,
        0.5,
        box_half_width=1.2,
        grid_size=nine,
        provenance="analytic rigid rotation",
    )
    assert result.grid_size == nine
    assert result.swirling_detected is True
    assert result.active_fraction == pytest.approx(1.0, abs=1e-14)
    assert result.swirling_strength_rms == pytest.approx(2.5, rel=1e-12, abs=1e-12)
    assert result.swirling_strength_max == pytest.approx(2.5, rel=1e-12, abs=1e-12)
    assert result.swirl_weighted_centroid == pytest.approx((0.0, 0.0, 0.0), abs=1e-14)
    assert result.divergence_rms == pytest.approx(0.25, rel=1e-12, abs=1e-12)
    assert result.velocity_changed is False
    assert result.visualization_ready is False
    assert result.pde_validated is False


def test_nonzero_shear_has_vorticity_but_zero_swirling_strength():
    result = diagnose_swirling_strength(
        _pure_shear,
        0.5,
        grid_size=9,
        provenance="analytic simple shear",
    )
    assert result.velocity_rms_speed > 0.0
    assert result.swirling_detected is False
    assert result.active_fraction == 0.0
    assert result.swirling_strength_rms == 0.0
    assert result.swirling_strength_max == 0.0
    assert result.swirl_weighted_centroid is None
    assert result.swirl_weighted_radius_rms is None
    assert result.swirl_weighted_abs_z_mean is None


def test_rotation_plus_real_strain_keeps_imaginary_rotation_rate():
    def field(points, time):
        x, y, z = points.T
        omega = 1.75
        a = 0.4
        return np.column_stack((a * x - omega * y, omega * x + a * y, -2.0 * a * z))

    result = diagnose_swirling_strength(field, 0.4, grid_size=11, provenance="linear spiral strain")
    assert result.swirling_strength_rms == pytest.approx(1.75, rel=1e-12, abs=1e-12)
    assert result.swirling_strength_max == pytest.approx(1.75, rel=1e-12, abs=1e-12)
    assert result.divergence_rms == pytest.approx(0.0, abs=1e-12)


def test_three_grid_audit_uses_same_field_without_acceptance_threshold():
    audit = audit_swirling_strength_resolution(
        _rigid_rotation,
        0.5,
        (7, 9, 13),
        box_half_width=1.1,
        provenance="rigid rotation resolution audit",
    )
    assert audit.grid_sizes == (7, 9, 13)
    assert audit.finest_grid_size == 13
    assert audit.swirling_rms_abs_delta_to_finest == pytest.approx((0.0, 0.0, 0.0), abs=1e-12)
    assert audit.swirling_max_abs_delta_to_finest == pytest.approx((0.0, 0.0, 0.0), abs=1e-12)
    assert audit.active_fraction_abs_delta_to_finest == pytest.approx((0.0, 0.0, 0.0), abs=1e-14)
    assert audit.visualization_ready is False
    assert audit.visual_correspondence_verified is False
    assert audit.pde_validated is False


def test_fail_closed_inputs_and_exact_zero_velocity():
    def zero(points, time):
        return np.zeros_like(points)

    def bad_shape(points, time):
        return np.zeros((len(points), 2))

    def bad_values(points, time):
        out = np.ones_like(points)
        out[0, 0] = np.nan
        return out

    with pytest.raises(ValueError, match="exact-zero"):
        diagnose_swirling_strength(zero, 0.5, grid_size=9, provenance="zero")
    with pytest.raises(ValueError, match="shape"):
        diagnose_swirling_strength(bad_shape, 0.5, grid_size=9, provenance="bad shape")
    with pytest.raises(ValueError, match="finite"):
        diagnose_swirling_strength(bad_values, 0.5, grid_size=9, provenance="bad values")
    with pytest.raises(ValueError, match="odd integer"):
        diagnose_swirling_strength(_rigid_rotation, 0.5, grid_size=8, provenance="even grid")
    with pytest.raises(ValueError, match="provenance"):
        diagnose_swirling_strength(_rigid_rotation, 0.5, grid_size=9, provenance="")
    with pytest.raises(ValueError, match="at least three"):
        audit_swirling_strength_resolution(_rigid_rotation, 0.5, (7, 9), provenance="too short")
    with pytest.raises(ValueError, match="strictly increasing"):
        audit_swirling_strength_resolution(_rigid_rotation, 0.5, (7, 13, 9), provenance="unordered")
