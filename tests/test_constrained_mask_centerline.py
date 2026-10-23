import numpy as np
import pytest

from openai_ns_reconstruction.constrained_mask_centerline import (
    measure_presegmented_mask_centerline,
)


def _measure(mask, **kwargs):
    return measure_presegmented_mask_centerline(
        mask,
        mask_provenance="fixture-mask-v1",
        frame_provenance="fixed-frame-v1",
        **kwargs,
    )


def test_straight_centered_rectangle_has_zero_tilt_and_bending():
    mask = np.zeros((11, 13), dtype=bool)
    mask[2:9, 4:9] = True
    result = _measure(mask, sample_count=9)
    assert result.active_row_count == 7
    assert result.endpoint_axial_separation_physical == pytest.approx(6.0)
    assert result.axial_extent_physical == pytest.approx(7.0)
    assert result.global_transverse_centroid_physical == pytest.approx(6.5)
    assert result.bottom_center_physical == pytest.approx(6.5)
    assert result.midplane_center_physical == pytest.approx(6.5)
    assert result.top_center_physical == pytest.approx(6.5)
    assert result.end_to_end_tilt_slope == pytest.approx(0.0)
    assert result.lateral_drift_physical == pytest.approx(0.0)
    assert result.rms_chord_deviation_physical == pytest.approx(0.0)
    assert result.max_chord_deviation_physical == pytest.approx(0.0)
    assert result.straightness_ratio == pytest.approx(1.0)


def test_fixed_frame_horizontal_translation_is_not_removed():
    base = np.zeros((11, 15), dtype=bool)
    base[2:9, 4:8] = True
    shifted = np.zeros_like(base)
    shifted[2:9, 7:11] = True
    a = _measure(base, sample_count=9)
    b = _measure(shifted, sample_count=9)
    assert b.global_transverse_centroid_physical - a.global_transverse_centroid_physical == pytest.approx(3.0)
    assert b.midplane_center_physical - a.midplane_center_physical == pytest.approx(3.0)
    assert a.rms_chord_deviation_physical == pytest.approx(b.rms_chord_deviation_physical)


def test_linear_tilt_remains_visible_but_has_zero_chord_bending():
    mask = np.zeros((9, 15), dtype=bool)
    for row, center in zip(range(2, 7), [8, 7, 6, 5, 4]):
        mask[row, center - 1 : center + 2] = True
    result = _measure(mask, sample_count=9)
    assert result.top_center_physical - result.bottom_center_physical == pytest.approx(4.0)
    assert result.end_to_end_tilt_slope == pytest.approx(1.0)
    assert result.rms_chord_deviation_physical == pytest.approx(0.0)
    assert result.max_chord_deviation_physical == pytest.approx(0.0)
    assert result.straightness_ratio == pytest.approx(1.0)


def test_bent_centerline_reports_chord_deviation_and_arclength_excess():
    mask = np.zeros((9, 17), dtype=bool)
    top_to_bottom_centers = [5, 6, 8, 6, 5]
    for row, center in zip(range(2, 7), top_to_bottom_centers):
        mask[row, center - 1 : center + 2] = True
    result = _measure(mask, sample_count=9)
    assert result.end_to_end_tilt_slope == pytest.approx(0.0)
    assert result.midplane_center_physical - result.bottom_center_physical == pytest.approx(3.0)
    assert result.max_chord_deviation_physical == pytest.approx(3.0)
    assert result.rms_chord_deviation_physical > 0.0
    assert result.straightness_ratio > 1.0


def test_hole_in_row_uses_foreground_mass_centroid_not_outer_span_midpoint():
    mask = np.zeros((9, 15), dtype=bool)
    mask[2:7, 4:10] = True
    mask[4, 5:8] = False
    result = _measure(mask, sample_count=5)
    assert result.midplane_center_physical == pytest.approx(7.5)
    assert result.global_transverse_centroid_physical != pytest.approx(7.0)


def test_anisotropic_spacing_scales_physical_tilt_and_frame_fraction():
    mask = np.zeros((9, 15), dtype=bool)
    for row, center in zip(range(2, 7), [8, 7, 6, 5, 4]):
        mask[row, center - 1 : center + 2] = True
    result = _measure(mask, pixel_spacing=(2.0, 0.5), sample_count=9)
    assert result.endpoint_axial_separation_physical == pytest.approx(8.0)
    assert result.top_center_physical - result.bottom_center_physical == pytest.approx(2.0)
    assert result.end_to_end_tilt_slope == pytest.approx(0.25)
    assert result.frame_transverse_extent_physical == pytest.approx(7.5)
    assert result.lateral_drift_frame_fraction == pytest.approx(2.0 / 7.5)


def test_axial_up_orientation_is_explicit():
    mask = np.zeros((9, 13), dtype=bool)
    for row, center in zip(range(2, 7), [8, 7, 6, 5, 4]):
        mask[row, center] = True
        mask[row, center + 1] = True
    result = _measure(mask, sample_count=5)
    assert result.normalized_axial[0] == pytest.approx(-1.0)
    assert result.normalized_axial[-1] == pytest.approx(1.0)
    assert result.axial_position_physical[0] < result.axial_position_physical[-1]
    assert result.top_center_physical > result.bottom_center_physical


def test_frame_contact_is_reported_not_cropped():
    mask = np.zeros((9, 11), dtype=bool)
    mask[2:7, 0:3] = True
    result = _measure(mask, sample_count=5)
    assert result.touches_image_frame is True


def test_returned_profiles_are_read_only_and_truth_flags_remain_false():
    mask = np.zeros((9, 11), dtype=bool)
    mask[2:7, 4:7] = True
    result = _measure(mask, sample_count=5)
    for array in (
        result.normalized_axial,
        result.axial_position_physical,
        result.transverse_center_physical,
        result.chord_residual_physical,
    ):
        assert array.flags.writeable is False
    assert result.segmentation_performed is False
    assert result.threshold_selected is False
    assert result.camera_fitted is False
    assert result.registration_fitted is False
    assert result.rotation_fitted is False
    assert result.velocity_changed is False
    assert result.visualization_ready is False
    assert result.visual_correspondence_verified is False
    assert result.pde_validated is False
    assert result.paper_exact is False
    assert result.openai_field_identified is False
    assert result.blowup_proved is False


@pytest.mark.parametrize(
    "mask",
    [
        np.zeros((7, 7), dtype=bool),
        np.ones((7, 7), dtype=bool),
        np.zeros((2, 7), dtype=bool),
        np.zeros((7, 2), dtype=bool),
    ],
)
def test_degenerate_masks_fail_closed(mask):
    with pytest.raises(ValueError):
        _measure(mask)


def test_grayscale_mask_fails_closed():
    mask = np.zeros((7, 7), dtype=float)
    mask[2:5, 2:5] = 0.5
    with pytest.raises(ValueError, match="grayscale"):
        _measure(mask)


def test_nonfinite_mask_fails_closed():
    mask = np.zeros((7, 7), dtype=float)
    mask[2:5, 2:5] = 1.0
    mask[3, 3] = np.nan
    with pytest.raises(ValueError, match="finite"):
        _measure(mask)


def test_only_two_active_rows_fail_closed():
    mask = np.zeros((7, 7), dtype=bool)
    mask[2:4, 2:5] = True
    with pytest.raises(ValueError, match="three active"):
        _measure(mask)


def test_axially_disconnected_foreground_fails_closed():
    mask = np.zeros((9, 9), dtype=bool)
    mask[1:3, 3:6] = True
    mask[5:8, 3:6] = True
    with pytest.raises(ValueError, match="contiguous"):
        _measure(mask)


@pytest.mark.parametrize(
    "spacing",
    [(0.0, 1.0), (1.0, 0.0), (-1.0, 1.0), (1.0, np.nan), (1.0, np.inf)],
)
def test_invalid_spacing_fails_closed(spacing):
    mask = np.zeros((7, 7), dtype=bool)
    mask[1:6, 2:5] = True
    with pytest.raises(ValueError, match="spacing"):
        _measure(mask, pixel_spacing=spacing)


@pytest.mark.parametrize("sample_count", [4, 6, 0, 4099, True, 5.5])
def test_invalid_sample_count_fails_closed(sample_count):
    mask = np.zeros((7, 7), dtype=bool)
    mask[1:6, 2:5] = True
    with pytest.raises(ValueError, match="sample_count"):
        _measure(mask, sample_count=sample_count)


@pytest.mark.parametrize("field", ["mask_provenance", "frame_provenance"])
def test_missing_provenance_fails_closed(field):
    mask = np.zeros((7, 7), dtype=bool)
    mask[1:6, 2:5] = True
    kwargs = {
        "mask_provenance": "fixture-mask-v1",
        "frame_provenance": "fixed-frame-v1",
    }
    kwargs[field] = "   "
    with pytest.raises(ValueError, match=field):
        measure_presegmented_mask_centerline(mask, **kwargs)
