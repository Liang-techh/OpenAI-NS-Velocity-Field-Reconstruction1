import numpy as np
import pytest

from openai_ns_reconstruction.constrained_public_mask_distance import (
    compare_presegmented_masks,
)


def ring_mask(shape=(41, 41), outer=13.0, inner=5.0):
    yy, xx = np.indices(shape, dtype=float)
    cy = 0.5 * (shape[0] - 1)
    cx = 0.5 * (shape[1] - 1)
    rr = np.hypot(yy - cy, xx - cx)
    return (rr <= outer) & (rr >= inner)


def test_identical_mask_is_exact_without_claim_promotion():
    reference = ring_mask()
    report = compare_presegmented_masks(reference.copy(), reference)

    assert report.intersection_over_union == 1.0
    assert report.dice_coefficient == 1.0
    assert report.symmetric_surface_mean == 0.0
    assert report.symmetric_surface_rms == 0.0
    assert report.symmetric_surface_q90 == 0.0
    assert report.symmetric_hausdorff == 0.0
    assert report.alignment == "none_same_pixel_frame_required"
    assert report.segmentation_performed is False
    assert report.threshold_selected is False
    assert report.camera_fitted is False
    assert report.visualization_ready is False
    assert report.visual_correspondence_verified is False
    assert report.pde_validated is False
    assert report.paper_exact is False
    assert report.openai_field_identified is False
    assert report.blowup_proved is False


def test_same_frame_shift_is_not_registered_away():
    reference = ring_mask()
    candidate = np.zeros_like(reference)
    candidate[2:, :] = reference[:-2, :]
    report = compare_presegmented_masks(candidate, reference)

    assert report.intersection_over_union < 0.8
    assert report.dice_coefficient < 0.9
    assert report.symmetric_surface_mean > 0.2
    assert report.symmetric_hausdorff >= 2.0
    assert report.normalized_hausdorff > 0.0


def test_hole_mismatch_remains_visible_beyond_outer_envelope():
    reference = ring_mask(outer=14.0, inner=7.0)
    yy, xx = np.indices(reference.shape, dtype=float)
    rr = np.hypot(yy - 20.0, xx - 20.0)
    filled_candidate = rr <= 14.0

    report = compare_presegmented_masks(filled_candidate, reference)

    assert report.intersection_over_union < 0.8
    assert report.area_ratio_candidate_to_reference > 1.2
    assert report.symmetric_surface_mean > 0.5
    assert report.symmetric_hausdorff >= 6.0


def test_declared_anisotropic_spacing_changes_only_distance_units():
    reference = ring_mask()
    candidate = np.zeros_like(reference)
    candidate[:, 1:] = reference[:, :-1]

    unit = compare_presegmented_masks(candidate, reference)
    stretched = compare_presegmented_masks(
        candidate, reference, pixel_spacing=(1.0, 2.0)
    )

    assert stretched.intersection_over_union == unit.intersection_over_union
    assert stretched.dice_coefficient == unit.dice_coefficient
    assert stretched.symmetric_surface_mean > unit.symmetric_surface_mean
    assert stretched.pixel_spacing_axial_transverse == (1.0, 2.0)


def test_frame_contact_is_reported_not_interpreted_as_support():
    reference = ring_mask()
    candidate = reference.copy()
    candidate[:, 0] = True
    report = compare_presegmented_masks(candidate, reference)

    assert report.candidate_touches_image_frame is True
    assert report.reference_touches_image_frame is False


@pytest.mark.parametrize(
    "candidate, reference, spacing, error",
    [
        (np.zeros((8, 8)), ring_mask((8, 8), 3, 1), (1, 1), ValueError),
        (np.ones((8, 8)), ring_mask((8, 8), 3, 1), (1, 1), ValueError),
        (np.full((8, 8), 0.5), ring_mask((8, 8), 3, 1), (1, 1), ValueError),
        (np.full((8, 8), np.nan), ring_mask((8, 8), 3, 1), (1, 1), ValueError),
        (ring_mask((9, 9), 3, 1), ring_mask((8, 8), 3, 1), (1, 1), ValueError),
        (ring_mask((8, 8), 3, 1), ring_mask((8, 8), 3, 1), (0, 1), ValueError),
        (ring_mask((8, 8), 3, 1), ring_mask((8, 8), 3, 1), (1,), ValueError),
    ],
)
def test_fail_closed_inputs(candidate, reference, spacing, error):
    with pytest.raises(error):
        compare_presegmented_masks(candidate, reference, pixel_spacing=spacing)
