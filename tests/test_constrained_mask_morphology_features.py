import numpy as np
import pytest

from openai_ns_reconstruction.constrained_mask_morphology_features import (
    measure_presegmented_mask_morphology,
)


KW = dict(
    mask_id="public-frame-1",
    segmentation_provenance="explicit public-observable mask v1",
    frame_provenance="fixed image frame v1",
)


def _rect():
    mask = np.zeros((11, 13), dtype=bool)
    mask[3:8, 4:10] = True
    return mask


def test_rectangle_axis_fixed_features():
    result = measure_presegmented_mask_morphology(_rect(), **KW)
    assert result.foreground_pixels == 30
    assert result.area_fraction == pytest.approx(30 / 143)
    assert result.physical_area == pytest.approx(30.0)
    assert result.transverse_extent == pytest.approx(6.0)
    assert result.axial_extent == pytest.approx(5.0)
    assert result.transverse_centroid == pytest.approx(0.5)
    assert result.axial_centroid == pytest.approx(0.0)
    assert result.bbox_aspect_axial_over_transverse == pytest.approx(5 / 6)
    assert result.hole_pixels == 0
    assert result.hole_area_fraction_of_filled == 0.0
    assert not result.touches_image_frame


def test_anisotropic_spacing_scales_axes_and_area():
    result = measure_presegmented_mask_morphology(
        _rect(), pixel_spacing=(2.0, 0.5), **KW
    )
    assert result.axial_extent == pytest.approx(10.0)
    assert result.transverse_extent == pytest.approx(3.0)
    assert result.physical_area == pytest.approx(30.0)
    assert result.bbox_aspect_axial_over_transverse == pytest.approx(10 / 3)


def test_internal_hole_is_measured_not_filled_away():
    mask = np.zeros((11, 11), dtype=bool)
    mask[2:9, 2:9] = True
    mask[4:7, 4:7] = False
    result = measure_presegmented_mask_morphology(mask, **KW)
    assert result.foreground_pixels == 40
    assert result.hole_pixels == 9
    assert result.hole_area_fraction_of_filled == pytest.approx(9 / 49)


def test_translation_changes_centroid_but_not_shape_features():
    a = np.zeros((13, 15), dtype=bool)
    a[3:7, 4:9] = True
    b = np.zeros_like(a)
    b[5:9, 7:12] = True
    ma = measure_presegmented_mask_morphology(a, **KW)
    mb = measure_presegmented_mask_morphology(b, **{**KW, "mask_id": "shifted"})
    assert ma.transverse_extent == mb.transverse_extent
    assert ma.axial_extent == mb.axial_extent
    assert ma.transverse_rms == pytest.approx(mb.transverse_rms)
    assert ma.axial_rms == pytest.approx(mb.axial_rms)
    assert ma.transverse_centroid != mb.transverse_centroid
    assert ma.axial_centroid != mb.axial_centroid


def test_frame_contact_is_reported_without_becoming_support_claim():
    mask = _rect()
    mask[:, 0] = False
    mask[3:6, 0:2] = True
    result = measure_presegmented_mask_morphology(mask, **KW)
    assert result.touches_image_frame
    assert result.truth_boundary["visualization_ready"] is False
    assert result.truth_boundary["pde_validated"] is False
    assert result.truth_boundary["openai_field_identified"] is False


@pytest.mark.parametrize(
    "mask, override, message",
    [
        (np.zeros((5, 5), dtype=bool), {}, "foreground"),
        (np.ones((5, 5), dtype=bool), {}, "entire frame"),
        (np.zeros((2, 5), dtype=bool), {}, "dimensions at least 3"),
        (np.zeros((5, 5, 1), dtype=bool), {}, "2-D"),
        (np.full((5, 5), 0.5), {}, "exact numeric 0/1"),
        (np.array([[0, 0, 0], [0, np.nan, 0], [0, 1, 0]]), {}, "finite"),
        (_rect(), {"pixel_spacing": (0.0, 1.0)}, "strictly positive"),
        (_rect(), {"pixel_spacing": (1.0, np.inf)}, "finite"),
        (_rect(), {"mask_id": ""}, "mask_id"),
        (_rect(), {"segmentation_provenance": " "}, "segmentation_provenance"),
        (_rect(), {"frame_provenance": ""}, "frame_provenance"),
    ],
)
def test_bad_inputs_fail_closed(mask, override, message):
    kwargs = {**KW, **override}
    with pytest.raises(ValueError, match=message):
        measure_presegmented_mask_morphology(mask, **kwargs)


def test_one_pixel_wide_axis_fails_closed():
    mask = np.zeros((7, 7), dtype=bool)
    mask[2:5, 3] = True
    with pytest.raises(ValueError, match="span at least two distinct"):
        measure_presegmented_mask_morphology(mask, **KW)
