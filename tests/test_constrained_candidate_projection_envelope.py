import numpy as np
import pytest

from openai_ns_reconstruction.constrained_candidate_projection_envelope import (
    project_candidate_meridional_envelope,
)


def _rectangle_prism_points():
    xz = np.array([
        [-2.0, -1.0],
        [2.0, -1.0],
        [2.0, 3.0],
        [-2.0, 3.0],
        [0.0, 0.0],
        [0.5, 1.0],
    ])
    points = []
    for y in (-7.0, 0.0, 11.0):
        for x, z in xz:
            points.append([x, y, z])
    return np.asarray(points)


def test_xz_projection_ignores_depth_and_recovers_outer_rectangle():
    result = project_candidate_meridional_envelope(
        _rectangle_prism_points(),
        candidate_id="manufactured-prism",
        observable="outer_cloud",
        projection="xz",
    )

    expected = np.array([
        [-2.0, -1.0],
        [2.0, -1.0],
        [2.0, 3.0],
        [-2.0, 3.0],
    ])
    np.testing.assert_allclose(result.boundary_points, expected)
    assert result.projected_area == pytest.approx(16.0)
    assert result.projected_perimeter == pytest.approx(16.0)
    assert result.source_point_count == 18
    assert result.unique_projected_point_count == 6
    assert result.boundary_points.flags.writeable is False
    assert result.camera_fitted is False
    assert result.rotation_fitted is False
    assert result.velocity_changed is False
    assert result.pde_validated is False
    assert result.openai_field_identified is False


def test_projection_choice_is_explicit_and_not_camera_fitted():
    points = np.array([
        [x, y, z]
        for x in (-4.0, 4.0)
        for y in (-0.5, 0.5)
        for z in (-1.0, 1.0)
    ] + [[0.0, 0.0, 0.0]])

    xz = project_candidate_meridional_envelope(
        points, candidate_id="anisotropic", observable="cloud", projection="xz"
    )
    yz = project_candidate_meridional_envelope(
        points, candidate_id="anisotropic", observable="cloud", projection="yz"
    )

    assert xz.projected_area == pytest.approx(16.0)
    assert yz.projected_area == pytest.approx(2.0)
    assert xz.projection == "xz"
    assert yz.projection == "yz"


def test_boundary_order_is_deterministic_under_input_permutation():
    points = _rectangle_prism_points()
    a = project_candidate_meridional_envelope(
        points, candidate_id="same", observable="same", projection="xz"
    )
    b = project_candidate_meridional_envelope(
        points[::-1], candidate_id="same", observable="same", projection="xz"
    )
    np.testing.assert_array_equal(a.boundary_points, b.boundary_points)


@pytest.mark.parametrize(
    "points, kwargs, message",
    [
        (np.zeros((3, 2)), {}, "shape"),
        (np.array([[0.0, 0.0, 0.0], [1.0, 0.0, np.nan], [0.0, 0.0, 1.0]]), {}, "finite"),
        (np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 1.0], [2.0, 0.0, 2.0]]), {}, "degenerate"),
        (_rectangle_prism_points(), {"projection": "xy"}, "projection"),
    ],
)
def test_invalid_geometry_fails_closed(points, kwargs, message):
    with pytest.raises(ValueError, match=message):
        project_candidate_meridional_envelope(
            points,
            candidate_id="candidate",
            observable="observable",
            **kwargs,
        )


def test_empty_identity_and_observable_fail_closed():
    points = _rectangle_prism_points()
    with pytest.raises(ValueError, match="candidate_id"):
        project_candidate_meridional_envelope(points, candidate_id=" ", observable="cloud")
    with pytest.raises(ValueError, match="observable"):
        project_candidate_meridional_envelope(points, candidate_id="candidate", observable="")
