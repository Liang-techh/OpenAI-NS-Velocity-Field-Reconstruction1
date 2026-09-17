import numpy as np
import pytest

from openai_ns_reconstruction.constrained_mask_width_profile import (
    measure_presegmented_mask_width_profile,
)


def measure(mask, **kwargs):
    return measure_presegmented_mask_width_profile(
        mask,
        mask_provenance="explicit synthetic binary mask",
        frame_provenance="fixed xz orthographic frame",
        **kwargs,
    )


def test_rectangle_has_flat_profile_and_physical_extent():
    mask = np.zeros((12, 14), dtype=bool)
    mask[3:9, 4:10] = True
    result = measure(mask, pixel_spacing=(2.0, 0.5), sample_count=9)
    np.testing.assert_allclose(result.span_width_physical, 3.0)
    np.testing.assert_allclose(result.occupied_width_physical, 3.0)
    np.testing.assert_allclose(result.span_width_fraction, 1.0)
    np.testing.assert_allclose(result.occupancy_fraction, 1.0)
    assert result.axial_extent_physical == 12.0
    assert result.maximum_span_width_physical == 3.0
    assert result.minimum_span_width_physical == 3.0
    assert result.bottom_tip_span_fraction == 1.0
    assert result.midplane_span_fraction == 1.0
    assert result.top_tip_span_fraction == 1.0
    assert result.maximum_width_axial_location == pytest.approx(0.0)


def test_taper_profile_preserves_knots_without_overshoot():
    mask = np.zeros((9, 11), dtype=bool)
    widths = [1, 3, 5, 3, 1]
    for row, width in zip(range(2, 7), widths):
        center = 5
        half = width // 2
        mask[row, center - half:center + half + 1] = True
    result = measure(mask, sample_count=5)
    np.testing.assert_allclose(result.span_width_fraction, [0.2, 0.6, 1.0, 0.6, 0.2])
    assert result.maximum_width_axial_location == pytest.approx(0.0)
    assert result.bottom_tip_span_fraction == pytest.approx(0.2)
    assert result.top_tip_span_fraction == pytest.approx(0.2)


def test_hole_changes_occupancy_not_outer_span():
    mask = np.zeros((9, 11), dtype=bool)
    mask[2:7, 3:8] = True
    mask[4, 4:7] = False
    result = measure(mask, sample_count=5)
    np.testing.assert_allclose(result.span_width_fraction, 1.0)
    np.testing.assert_allclose(result.span_width_physical, 5.0)
    assert result.occupancy_fraction[2] == pytest.approx(0.4)
    assert result.midplane_span_fraction == 1.0


def test_translation_does_not_change_width_descriptor_but_frame_contact_is_separate():
    first = np.zeros((12, 14), dtype=bool)
    first[3:9, 3:9] = True
    second = np.zeros_like(first)
    second[3:9, 5:11] = True
    a = measure(first, sample_count=11)
    b = measure(second, sample_count=11)
    np.testing.assert_allclose(a.span_width_fraction, b.span_width_fraction)
    np.testing.assert_allclose(a.occupancy_fraction, b.occupancy_fraction)
    assert a.axial_extent_physical == b.axial_extent_physical


def test_axial_up_orientation_distinguishes_bottom_and_top_tip():
    mask = np.zeros((9, 11), dtype=bool)
    widths_top_to_bottom = [1, 3, 5, 5, 5]
    for row, width in zip(range(2, 7), widths_top_to_bottom):
        half = width // 2
        mask[row, 5 - half:5 + half + 1] = True
    result = measure(mask, sample_count=5)
    assert result.bottom_tip_span_fraction == pytest.approx(1.0)
    assert result.top_tip_span_fraction == pytest.approx(0.2)
    assert result.normalized_axial[0] == -1.0
    assert result.normalized_axial[-1] == 1.0


def test_anisotropic_spacing_changes_physical_widths_not_fractions():
    mask = np.zeros((9, 11), dtype=bool)
    mask[2:7, 3:8] = True
    a = measure(mask, pixel_spacing=(1.0, 1.0), sample_count=7)
    b = measure(mask, pixel_spacing=(3.0, 0.25), sample_count=7)
    np.testing.assert_allclose(a.span_width_fraction, b.span_width_fraction)
    assert b.axial_extent_physical == pytest.approx(3.0 * a.axial_extent_physical)
    assert b.maximum_span_width_physical == pytest.approx(0.25 * a.maximum_span_width_physical)


def test_frame_contact_is_reported_not_reinterpreted_as_support():
    mask = np.zeros((9, 11), dtype=bool)
    mask[2:7, 0:4] = True
    result = measure(mask)
    assert result.touches_image_frame is True
    assert result.visualization_ready is False
    assert result.visual_correspondence_verified is False
    assert result.pde_validated is False


def test_output_arrays_are_read_only():
    mask = np.zeros((9, 11), dtype=bool)
    mask[2:7, 3:8] = True
    result = measure(mask)
    for array in (
        result.normalized_axial,
        result.span_width_physical,
        result.occupied_width_physical,
        result.span_width_fraction,
        result.occupancy_fraction,
    ):
        assert not array.flags.writeable
        with pytest.raises(ValueError):
            array[0] = 0.0


@pytest.mark.parametrize(
    "mask",
    [
        np.zeros((7, 7), dtype=bool),
        np.ones((7, 7), dtype=bool),
        np.full((7, 7), 0.5),
        np.full((7, 7), np.nan),
        np.ones((2, 7), dtype=bool),
    ],
)
def test_invalid_masks_fail_closed(mask):
    with pytest.raises(ValueError):
        measure(mask)


def test_noncontiguous_active_rows_fail_closed():
    mask = np.zeros((9, 11), dtype=bool)
    mask[2, 3:8] = True
    mask[4, 3:8] = True
    mask[6, 3:8] = True
    with pytest.raises(ValueError, match="contiguous"):
        measure(mask)


def test_too_few_active_rows_fail_closed():
    mask = np.zeros((9, 11), dtype=bool)
    mask[3:5, 3:8] = True
    with pytest.raises(ValueError, match="three active"):
        measure(mask)


def test_transversely_degenerate_mask_fails_closed():
    mask = np.zeros((9, 11), dtype=bool)
    mask[2:7, 5] = True
    with pytest.raises(ValueError, match="transversely degenerate"):
        measure(mask)


@pytest.mark.parametrize("spacing", [(0.0, 1.0), (1.0, -1.0), (np.nan, 1.0), (1.0, np.inf)])
def test_invalid_spacing_fails_closed(spacing):
    mask = np.zeros((9, 11), dtype=bool)
    mask[2:7, 3:8] = True
    with pytest.raises(ValueError):
        measure(mask, pixel_spacing=spacing)


@pytest.mark.parametrize("count", [4, 6, 3, 4099, 7.0, True])
def test_invalid_sample_count_fails_closed(count):
    mask = np.zeros((9, 11), dtype=bool)
    mask[2:7, 3:8] = True
    with pytest.raises(ValueError):
        measure(mask, sample_count=count)


def test_missing_provenance_fails_closed():
    mask = np.zeros((9, 11), dtype=bool)
    mask[2:7, 3:8] = True
    with pytest.raises(ValueError):
        measure_presegmented_mask_width_profile(
            mask,
            mask_provenance=" ",
            frame_provenance="fixed frame",
        )


def test_truth_boundary_flags_remain_false():
    mask = np.zeros((9, 11), dtype=bool)
    mask[2:7, 3:8] = True
    result = measure(mask)
    assert result.segmentation_performed is False
    assert result.threshold_selected is False
    assert result.camera_fitted is False
    assert result.registration_fitted is False
    assert result.velocity_changed is False
    assert result.visualization_ready is False
    assert result.visual_correspondence_verified is False
    assert result.pde_validated is False
    assert result.paper_exact is False
    assert result.openai_field_identified is False
    assert result.blowup_proved is False
