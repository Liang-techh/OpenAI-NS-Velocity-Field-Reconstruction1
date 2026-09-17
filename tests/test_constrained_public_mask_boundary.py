import numpy as np
import pytest

from openai_ns_reconstruction.constrained_public_mask_boundary import (
    CLAIM_SCOPE,
    extract_public_mask_boundary,
)


def test_rectangle_boundary_and_axis_convention():
    mask = np.zeros((10, 12), dtype=bool)
    mask[2:8, 3:10] = True

    result = extract_public_mask_boundary(mask)

    assert result.image_shape == (10, 12)
    assert result.foreground_pixels == 42
    assert result.boundary_pixels == 22
    assert result.foreground_fraction == pytest.approx(42 / 120)
    assert result.touches_image_frame is False
    assert result.connectivity == 4
    assert result.claim_scope == CLAIM_SCOPE
    assert result.segmentation_performed is False
    assert result.camera_fitted is False
    assert result.visualization_ready is False
    assert result.visual_correspondence_verified is False
    assert result.pde_validated is False
    assert result.paper_exact is False
    assert result.openai_field_identified is False
    assert result.blowup_proved is False
    assert result.points.flags.writeable is False

    points = {tuple(p) for p in result.points}
    # image row 2 becomes axial coordinate 7; row 7 becomes axial coordinate 2
    assert (3.0, 7.0) in points
    assert (9.0, 2.0) in points
    # center foreground pixel is not on the one-pixel boundary
    assert (6.0, 4.0) not in points


def test_inner_hole_is_preserved_as_visible_boundary():
    mask = np.zeros((11, 11), dtype=np.uint8)
    mask[1:10, 1:10] = 1
    mask[4:7, 4:7] = 0

    result = extract_public_mask_boundary(mask, min_boundary_points=8)
    points = {tuple(p) for p in result.points}

    # The outer silhouette and the visible hole boundary are both retained.
    assert (1.0, 9.0) in points
    assert (4.0, 7.0) in points
    assert (5.0, 7.0) in points
    assert (6.0, 7.0) in points
    assert result.touches_image_frame is False


def test_frame_touch_is_reported_not_hidden():
    mask = np.zeros((8, 9), dtype=bool)
    mask[0:5, 2:7] = True

    result = extract_public_mask_boundary(mask)
    assert result.touches_image_frame is True
    assert np.max(result.points[:, 1]) == 7.0


@pytest.mark.parametrize(
    "mask, error",
    [
        (np.zeros((5, 5), dtype=bool), "foreground"),
        (np.ones((5, 5), dtype=bool), "background"),
        (np.full((5, 5), 0.25), "already be segmented"),
        (
            np.array(
                [[0.0, np.nan, 1.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]]
            ),
            "finite",
        ),
        (np.zeros((2, 5), dtype=bool), "at least 3"),
        (np.zeros((3, 3, 1), dtype=bool), "two-dimensional"),
    ],
)
def test_invalid_masks_fail_closed(mask, error):
    with pytest.raises((TypeError, ValueError), match=error):
        extract_public_mask_boundary(mask, min_boundary_points=4)


def test_invalid_boundary_requirement_and_tiny_trace_fail_closed():
    mask = np.zeros((7, 7), dtype=bool)
    mask[3, 3] = True

    with pytest.raises(ValueError, match="fewer than required"):
        extract_public_mask_boundary(mask, min_boundary_points=4)

    with pytest.raises(TypeError, match="integer"):
        extract_public_mask_boundary(mask, min_boundary_points=True)

    with pytest.raises(ValueError, match="at least 4"):
        extract_public_mask_boundary(mask, min_boundary_points=3)
