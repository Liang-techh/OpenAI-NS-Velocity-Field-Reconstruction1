import numpy as np
import pytest

from openai_ns_reconstruction.constrained_mask_connectivity import (
    measure_presegmented_mask_connectivity,
)


def measure(mask, **kwargs):
    return measure_presegmented_mask_connectivity(
        mask,
        mask_provenance=kwargs.pop("mask_provenance", "explicit synthetic mask"),
        frame_provenance=kwargs.pop("frame_provenance", "fixed xz raster frame"),
        **kwargs,
    )


def test_solid_rectangle_is_one_unfragmented_component():
    mask = np.zeros((11, 13), dtype=bool)
    mask[3:8, 4:10] = True
    result = measure(mask)
    assert result.component_count == 1
    assert result.component_sizes_pixels.tolist() == [30]
    assert result.foreground_pixels == 30
    assert result.largest_component_fraction == pytest.approx(1.0)
    assert result.detached_foreground_fraction == pytest.approx(0.0)
    assert result.axial_segment_count == 1
    assert result.axial_gap_rows == 0
    assert result.maximum_transverse_runs == 1
    assert result.fragmented_active_row_fraction == pytest.approx(0.0)
    assert not result.touches_image_frame


def test_separated_horizontal_layers_are_counted_without_bridging():
    mask = np.zeros((17, 15), dtype=bool)
    for row in (2, 5, 8, 11, 14):
        mask[row, 4:11] = True
    result = measure(mask)
    assert result.component_count == 5
    assert result.component_sizes_pixels.tolist() == [7, 7, 7, 7, 7]
    assert result.axial_segment_count == 5
    assert result.axial_gap_rows == 8
    assert result.largest_axial_gap_rows == 2
    assert result.largest_component_fraction == pytest.approx(0.2)
    assert result.detached_foreground_fraction == pytest.approx(0.8)
    assert result.maximum_transverse_runs == 1


def test_side_by_side_components_share_one_axial_segment_but_split_rows():
    mask = np.zeros((13, 17), dtype=bool)
    mask[4:9, 2:6] = True
    mask[4:9, 11:15] = True
    result = measure(mask)
    assert result.component_count == 2
    assert result.axial_segment_count == 1
    assert result.axial_gap_rows == 0
    assert result.maximum_transverse_runs == 2
    assert result.fragmented_active_row_fraction == pytest.approx(1.0)
    assert result.largest_component_fraction == pytest.approx(0.5)


def test_hole_stays_one_component_but_creates_multirun_rows():
    mask = np.zeros((15, 15), dtype=bool)
    mask[3:12, 3:12] = True
    mask[6:9, 6:9] = False
    result = measure(mask)
    assert result.component_count == 1
    assert result.axial_segment_count == 1
    assert result.largest_component_fraction == pytest.approx(1.0)
    assert result.maximum_transverse_runs == 2
    assert result.fragmented_active_row_fraction == pytest.approx(3 / 9)


def test_diagonal_only_contact_is_separate_under_fixed_four_neighbor_rule():
    mask = np.zeros((9, 9), dtype=bool)
    mask[2:4, 2:4] = True
    mask[4:6, 4:6] = True
    result = measure(mask)
    assert result.component_count == 2
    assert result.connectivity == "fixed_4_neighbor"
    assert result.component_sizes_pixels.tolist() == [4, 4]


def test_anisotropic_spacing_affects_area_and_axial_gap_length_only_as_declared():
    mask = np.zeros((12, 12), dtype=bool)
    mask[2:4, 3:7] = True
    mask[7:9, 3:7] = True
    result = measure(mask, pixel_spacing=(0.25, 2.0))
    assert result.foreground_pixels == 16
    assert result.foreground_area_physical == pytest.approx(8.0)
    assert result.axial_gap_rows == 3
    assert result.largest_axial_gap_physical == pytest.approx(0.75)


def test_frame_contact_is_reported_not_interpreted_as_support():
    mask = np.zeros((9, 9), dtype=bool)
    mask[0:3, 3:6] = True
    result = measure(mask)
    assert result.touches_image_frame


def test_component_sizes_are_sorted_and_read_only():
    mask = np.zeros((14, 14), dtype=bool)
    mask[2:6, 2:6] = True
    mask[8:10, 8:11] = True
    result = measure(mask)
    assert result.component_sizes_pixels.tolist() == [16, 6]
    assert not result.component_sizes_pixels.flags.writeable
    with pytest.raises(ValueError):
        result.component_sizes_pixels[0] = 99


def test_truth_flags_stay_false():
    mask = np.zeros((9, 9), dtype=bool)
    mask[2:7, 3:6] = True
    result = measure(mask)
    assert result.segmentation_performed is False
    assert result.threshold_selected is False
    assert result.gap_bridging_performed is False
    assert result.component_cleanup_performed is False
    assert result.camera_fitted is False
    assert result.registration_fitted is False
    assert result.velocity_changed is False
    assert result.visualization_ready is False
    assert result.visual_correspondence_verified is False
    assert result.pde_validated is False
    assert result.paper_exact is False
    assert result.openai_field_identified is False
    assert result.blowup_proved is False


@pytest.mark.parametrize(
    "bad_mask",
    [
        np.zeros((5, 5), dtype=bool),
        np.ones((5, 5), dtype=bool),
        np.array([[0.0, 0.5, 1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]]),
        np.array([[0.0, np.nan, 1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]]),
        np.zeros((5,), dtype=bool),
        np.zeros((2, 4), dtype=bool),
    ],
)
def test_bad_masks_fail_closed(bad_mask):
    with pytest.raises(ValueError):
        measure(bad_mask)


@pytest.mark.parametrize(
    "spacing",
    [(0.0, 1.0), (1.0, 0.0), (-1.0, 1.0), (np.inf, 1.0), (1.0, np.nan)],
)
def test_bad_spacing_fails_closed(spacing):
    mask = np.zeros((7, 7), dtype=bool)
    mask[2:5, 2:5] = True
    with pytest.raises(ValueError):
        measure(mask, pixel_spacing=spacing)


def test_bad_spacing_length_fails_closed():
    mask = np.zeros((7, 7), dtype=bool)
    mask[2:5, 2:5] = True
    with pytest.raises(ValueError):
        measure(mask, pixel_spacing=(1.0,))


@pytest.mark.parametrize(
    "kwargs",
    [
        {"mask_provenance": ""},
        {"frame_provenance": "   "},
        {"mask_provenance": None},
        {"frame_provenance": None},
    ],
)
def test_missing_provenance_fails_closed(kwargs):
    mask = np.zeros((7, 7), dtype=bool)
    mask[2:5, 2:5] = True
    with pytest.raises(ValueError):
        measure(mask, **kwargs)
