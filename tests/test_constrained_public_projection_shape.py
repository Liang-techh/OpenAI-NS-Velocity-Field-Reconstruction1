import numpy as np
import pytest

from openai_ns_reconstruction.constrained_public_projection_shape import (
    compare_public_projection_shapes,
)


def ellipse_points(n=256, a=0.55, b=1.8):
    theta = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
    return np.column_stack((a * np.cos(theta), b * np.sin(theta)))


def test_translation_and_isotropic_scale_are_removed_without_claim_promotion():
    reference = ellipse_points()
    candidate = 3.7 * reference + np.array([12.0, -4.0])
    report = compare_public_projection_shapes(candidate, reference)

    assert report.symmetric_chamfer_rms < 2e-15
    assert report.symmetric_hausdorff < 5e-15
    assert report.log_aspect_error < 1e-15
    assert report.alignment.endswith("scale_only")
    assert report.claim_scope == "public_observable_projection_geometry_only"
    assert report.velocity_changed is False
    assert report.visualization_ready is False
    assert report.visual_correspondence_verified is False
    assert report.pde_validated is False
    assert report.paper_exact is False
    assert report.openai_field_identified is False
    assert report.blowup_proved is False


def test_axial_shape_change_and_orientation_mismatch_remain_visible():
    reference = ellipse_points()
    squeezed = reference.copy()
    squeezed[:, 1] *= 0.55
    squeezed_report = compare_public_projection_shapes(squeezed, reference)

    assert squeezed_report.symmetric_chamfer_rms > 0.04
    assert squeezed_report.log_aspect_error > 0.5

    rotated = reference[:, ::-1].copy()
    rotated[:, 0] *= -1.0
    rotated_report = compare_public_projection_shapes(rotated, reference)

    # The comparator deliberately does not optimize a camera/rotation.
    assert rotated_report.symmetric_chamfer_rms > 0.08
    assert rotated_report.log_aspect_error > 2.0


def test_sparse_outlier_is_exposed_by_hausdorff_and_q90_stays_smaller():
    reference = ellipse_points()
    candidate = np.vstack((reference, np.array([[2.5, 2.5]])))
    report = compare_public_projection_shapes(candidate, reference)

    assert report.symmetric_hausdorff > 0.2
    assert report.candidate_to_reference_q90 < report.symmetric_hausdorff
    assert report.symmetric_chamfer_rms < report.symmetric_hausdorff


def test_fail_closed_inputs():
    reference = ellipse_points(32)

    with pytest.raises(ValueError):
        compare_public_projection_shapes(np.zeros((8, 2)), reference)
    with pytest.raises(ValueError):
        compare_public_projection_shapes(reference[:3], reference)
    with pytest.raises(ValueError):
        compare_public_projection_shapes(np.column_stack((np.zeros(32), np.arange(32))), reference)

    bad = reference.copy()
    bad[0, 0] = np.nan
    with pytest.raises(ValueError):
        compare_public_projection_shapes(bad, reference)

    with pytest.raises(TypeError):
        compare_public_projection_shapes(reference, reference, min_points=8.5)
    with pytest.raises(ValueError):
        compare_public_projection_shapes(reference, reference, min_points=3)
