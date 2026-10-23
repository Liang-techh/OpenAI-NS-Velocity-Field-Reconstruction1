import numpy as np
import pytest

from openai_ns_reconstruction.constrained_public_vortex_features import (
    audit_public_vortex_feature_resolution,
    diagnose_public_vortex_features,
)


def _cylindrical_linear(*, radial, swirl, stretch):
    def velocity(points, time):
        points = np.asarray(points, dtype=float)
        x, y, z = points.T
        return np.column_stack(
            (
                radial * x - swirl * y,
                radial * y + swirl * x,
                stretch * z,
            )
        )

    return velocity


def test_public_observables_are_positive_for_inward_spiral_and_axial_stretch():
    velocity = _cylindrical_linear(radial=-0.6, swirl=1.4, stretch=0.8)
    report = diagnose_public_vortex_features(
        velocity,
        0.5,
        radii=(0.4, 0.8, 1.2),
        z_half_levels=(0.3, 0.7),
        angular_samples=32,
    )
    assert report.inward_radial_margin > 0.2
    assert report.circulation_fraction > 0.6
    assert report.inward_spiral_score > 0.2
    assert report.axial_stretch_margin > 0.1
    assert report.active_fraction == pytest.approx(1.0)
    assert report.pde_validated is False
    assert report.paper_exact is False
    assert report.openai_field_identified is False


def test_mutation_flips_inward_and_stretch_signs_without_erasing_swirl():
    desired = diagnose_public_vortex_features(
        _cylindrical_linear(radial=-0.5, swirl=1.0, stretch=0.7),
        0.5,
        radii=(0.5, 1.0),
        z_half_levels=(0.4, 0.8),
    )
    mutated = diagnose_public_vortex_features(
        _cylindrical_linear(radial=0.5, swirl=1.0, stretch=-0.7),
        0.5,
        radii=(0.5, 1.0),
        z_half_levels=(0.4, 0.8),
    )
    assert desired.inward_radial_margin > 0.0
    assert desired.inward_spiral_score > 0.0
    assert desired.axial_stretch_margin > 0.0
    assert mutated.inward_radial_margin < 0.0
    assert mutated.inward_spiral_score < 0.0
    assert mutated.axial_stretch_margin < 0.0
    assert mutated.circulation_fraction > 0.5


def test_axisymmetric_features_are_stable_across_three_angular_resolutions():
    audit = audit_public_vortex_feature_resolution(
        _cylindrical_linear(radial=-0.6, swirl=1.4, stretch=0.8),
        0.5,
        radii=(0.4, 0.8, 1.2),
        z_half_levels=(0.3, 0.7),
        angular_samples=(12, 24, 48),
    )
    assert audit.angular_samples == (12, 24, 48)
    for delta in audit.absolute_deltas_to_finest:
        assert max(delta) < 2e-15


def test_pure_rotation_does_not_fake_inward_or_axial_stretch_and_inputs_fail_closed():
    report = diagnose_public_vortex_features(
        _cylindrical_linear(radial=0.0, swirl=1.0, stretch=0.0),
        0.5,
        radii=(0.5, 1.0),
        z_half_levels=(0.25, 0.75),
    )
    assert report.inward_radial_margin == pytest.approx(0.0, abs=1e-15)
    assert report.inward_spiral_score == pytest.approx(0.0, abs=1e-15)
    assert report.axial_stretch_margin == pytest.approx(0.0, abs=1e-15)
    assert report.circulation_fraction == pytest.approx(1.0, abs=1e-12)

    with pytest.raises(ValueError, match="numerically inactive"):
        diagnose_public_vortex_features(
            lambda points, time: np.zeros_like(points),
            0.5,
            radii=(0.5,),
            z_half_levels=(0.5,),
        )
    with pytest.raises(ValueError, match="shape"):
        diagnose_public_vortex_features(
            lambda points, time: np.zeros((len(points), 2)),
            0.5,
            radii=(0.5,),
            z_half_levels=(0.5,),
        )
    with pytest.raises(ValueError, match="at least three"):
        audit_public_vortex_feature_resolution(
            _cylindrical_linear(radial=-1.0, swirl=1.0, stretch=1.0),
            0.5,
            radii=(0.5,),
            z_half_levels=(0.5,),
            angular_samples=(16, 32),
        )
