import numpy as np
import pytest

from openai_ns_reconstruction.constrained_vorticity_isosurface import (
    audit_vorticity_isosurface_resolution,
    diagnose_vorticity_isosurface,
)


def solid_rotation(points, time):
    p = np.asarray(points, dtype=float)
    return np.column_stack((-p[:, 1], p[:, 0], np.zeros(p.shape[0])))


def gaussian_vortices(centers):
    centers = tuple(float(c) for c in centers)
    sigma2 = 0.16
    z2 = 1.2

    def velocity(points, time):
        p = np.asarray(points, dtype=float)
        x, y, z = p[:, 0], p[:, 1], p[:, 2]
        u = np.zeros_like(x)
        v = np.zeros_like(x)
        for center in centers:
            g = np.exp(-((x - center) ** 2 + y**2) / sigma2 - z**2 / z2)
            u += -2.0 * y / sigma2 * g
            v += 2.0 * (x - center) / sigma2 * g
        return np.column_stack((u, v, np.zeros_like(x)))

    return velocity


def test_solid_rotation_is_one_connected_superlevel_set():
    report = diagnose_vorticity_isosurface(
        solid_rotation,
        0.5,
        grid_size=17,
        box_half_width=1.0,
        iso_fraction=0.5,
    )
    assert report.omega_peak == pytest.approx(2.0, abs=2e-14)
    assert report.component_count == 1
    assert report.largest_component_fraction == pytest.approx(1.0)
    assert report.active_volume_fraction == pytest.approx(1.0)
    assert report.pde_validated is False
    assert report.paper_exact is False


def test_two_separated_vortex_cores_are_not_merged_by_metric():
    report = diagnose_vorticity_isosurface(
        gaussian_vortices((-0.72, 0.72)),
        0.5,
        grid_size=33,
        box_half_width=1.5,
        iso_fraction=0.72,
    )
    assert report.component_count == 2
    assert 0.45 < report.largest_component_fraction < 0.55
    assert report.largest_axial_span > 0.3
    assert report.largest_radial_extent > 0.5


def test_single_tube_resolution_audit_is_stable():
    audit = audit_vorticity_isosurface_resolution(
        gaussian_vortices((0.0,)),
        0.5,
        grid_sizes=(17, 25, 33),
        box_half_width=1.5,
        iso_fraction=0.6,
    )
    assert audit.component_count_stable
    assert all(level.component_count == 1 for level in audit.levels)
    assert audit.max_largest_component_fraction_delta < 1e-12
    assert audit.max_radial_extent_relative_delta < 0.35
    assert audit.max_axial_span_relative_delta < 0.35
    assert audit.pde_validated is False


def test_fail_closed_inputs():
    with pytest.raises(ValueError):
        diagnose_vorticity_isosurface(solid_rotation, 0.5, grid_size=16)
    with pytest.raises(ValueError):
        audit_vorticity_isosurface_resolution(solid_rotation, 0.5, grid_sizes=(17, 25))
    with pytest.raises(ValueError):
        diagnose_vorticity_isosurface(lambda p, t: np.zeros((len(p), 3)), 0.5)
    with pytest.raises(ValueError):
        diagnose_vorticity_isosurface(lambda p, t: np.zeros((len(p), 2)), 0.5)
