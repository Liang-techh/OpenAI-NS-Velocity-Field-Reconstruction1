import numpy as np
import pytest

from openai_ns_reconstruction.constrained_candidate_projection_mask import (
    rasterize_candidate_projection_points,
)


KW = dict(
    projection="xz",
    shape=(11, 11),
    transverse_bounds=(-1.0, 1.0),
    axial_bounds=(-1.0, 1.0),
    point_radius_pixels=0,
    candidate_id="candidate-a",
    point_selection_provenance="fixed vorticity selection protocol v1",
    frame_provenance="fixed orthographic public-comparison frame v1",
)


def test_fixed_xz_mapping_and_axis_direction():
    points = np.array([
        [-0.8, 99.0, 0.8],
        [0.0, -50.0, 0.0],
        [0.8, 12.0, -0.8],
    ])
    result = rasterize_candidate_projection_points(points, **KW)
    occupied = np.argwhere(result.raw_occupancy)
    assert occupied.tolist() == [[1, 1], [5, 5], [9, 9]]
    assert result.point_count == 3
    assert result.occupied_pixels_before_dilation == 3
    assert not result.touches_image_frame
    assert not result.mask.flags.writeable
    assert not result.raw_occupancy.flags.writeable


def test_yz_projection_uses_y_not_x():
    points = np.array([
        [100.0, -0.8, 0.0],
        [-100.0, 0.8, 0.0],
    ])
    result = rasterize_candidate_projection_points(
        points,
        **{**KW, "projection": "yz"},
    )
    assert np.argwhere(result.raw_occupancy).tolist() == [[5, 1], [5, 9]]


def test_declared_radius_one_uses_local_disk_dilation():
    points = np.array([[0.0, 0.0, 0.0]])
    result = rasterize_candidate_projection_points(
        points,
        **{**KW, "point_radius_pixels": 1},
    )
    assert result.occupied_pixels_before_dilation == 1
    assert result.occupied_pixels_after_dilation == 5
    assert result.mask[5, 5]
    assert result.mask[4, 5]
    assert result.mask[6, 5]
    assert result.mask[5, 4]
    assert result.mask[5, 6]
    assert not result.mask[4, 4]


def test_no_hole_filling_is_performed():
    # Eight points around the central pixel. Radius zero must preserve the hole.
    xy = [-0.18, 0.0, 0.18]
    points = []
    for z in xy:
        for x in xy:
            if x == 0.0 and z == 0.0:
                continue
            points.append([x, 0.0, z])
    result = rasterize_candidate_projection_points(np.array(points), **KW)
    assert not result.mask[5, 5]
    assert result.occupied_pixels_after_dilation == 8


def test_input_permutation_is_deterministic():
    points = np.array([
        [-0.3, 0.0, 0.4],
        [0.1, 0.0, 0.2],
        [0.5, 0.0, -0.7],
        [-0.9, 0.0, -0.1],
    ])
    a = rasterize_candidate_projection_points(points, **KW)
    b = rasterize_candidate_projection_points(points[[2, 0, 3, 1]], **KW)
    np.testing.assert_array_equal(a.mask, b.mask)


def test_frame_contact_is_reported_not_interpreted_as_support():
    points = np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
    result = rasterize_candidate_projection_points(points, **KW)
    assert result.touches_image_frame
    assert result.truth_boundary["visualization_ready"] is False
    assert result.truth_boundary["pde_validated"] is False
    assert result.truth_boundary["openai_field_identified"] is False


def test_out_of_frame_points_fail_closed_instead_of_cropping():
    with pytest.raises(ValueError, match="silent cropping is forbidden"):
        rasterize_candidate_projection_points(
            np.array([[1.01, 0.0, 0.0], [0.0, 0.0, 0.0]]),
            **KW,
        )


@pytest.mark.parametrize(
    "points, override, message",
    [
        (np.empty((0, 3)), {}, "N >= 1"),
        (np.zeros((2, 2)), {}, "shape \\(N, 3\\)"),
        (np.array([[np.nan, 0.0, 0.0]]), {}, "finite"),
        (np.array([[0.0, 0.0, 0.0]]), {"projection": "xy"}, "projection"),
        (np.array([[0.0, 0.0, 0.0]]), {"shape": (2, 11)}, "at least 3"),
        (np.array([[0.0, 0.0, 0.0]]), {"transverse_bounds": (1.0, -1.0)}, "strictly increasing"),
        (np.array([[0.0, 0.0, 0.0]]), {"point_radius_pixels": -1}, "non-negative integer"),
        (np.array([[0.0, 0.0, 0.0]]), {"candidate_id": ""}, "candidate_id"),
        (np.array([[0.0, 0.0, 0.0]]), {"point_selection_provenance": ""}, "point_selection_provenance"),
        (np.array([[0.0, 0.0, 0.0]]), {"frame_provenance": ""}, "frame_provenance"),
    ],
)
def test_bad_inputs_fail_closed(points, override, message):
    kwargs = {**KW, **override}
    with pytest.raises(ValueError, match=message):
        rasterize_candidate_projection_points(points, **kwargs)


def test_frame_filling_dilation_fails_closed():
    with pytest.raises(ValueError, match="filled the entire frame"):
        rasterize_candidate_projection_points(
            np.array([[0.0, 0.0, 0.0]]),
            **{**KW, "shape": (5, 5), "point_radius_pixels": 4},
        )
