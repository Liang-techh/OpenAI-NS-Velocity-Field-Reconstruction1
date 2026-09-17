import numpy as np
import pytest

from openai_ns_reconstruction.constrained_streamline_collection import (
    make_matplotlib_line3d_collection,
    pack_streamline_collection,
)


def _helices(n_lines=200, n_points=81):
    theta = np.linspace(0.0, 6.0 * np.pi, n_points)
    lines = []
    scalars = []
    for index in range(n_lines):
        phase = 2.0 * np.pi * index / n_lines
        radius = 0.25 + 0.55 * (index + 0.5) / n_lines
        z = np.linspace(-1.0, 1.0, n_points)
        line = np.column_stack(
            (
                radius * np.cos(theta + phase),
                radius * np.sin(theta + phase),
                z,
            )
        )
        lines.append(line)
        scalars.append(np.sqrt(radius * radius + (0.4 + 0.1 * z) ** 2))
    return lines, scalars


def test_packs_200_colored_lines_into_one_collection_artist(monkeypatch):
    lines, scalars = _helices()
    pack = pack_streamline_collection(lines, scalars)

    assert pack.line_count == 200
    assert pack.segment_count == 200 * 80
    assert pack.segments.shape == (16000, 2, 3)
    assert pack.segment_scalars.shape == (16000,)
    assert pack.line_segment_offsets.shape == (201,)
    assert np.array_equal(pack.line_segment_offsets, np.arange(201) * 80)
    assert not pack.segments.flags.writeable
    assert not pack.matlab_x.flags.writeable

    import sys
    import types

    class FakeLine3DCollection:
        def __init__(self, segments, *, linewidths, cmap):
            self.segments = np.asarray(segments)
            self.linewidths = linewidths
            self.cmap = cmap
            self.array = None
            self.clim = None

        def set_array(self, array):
            self.array = np.asarray(array)

        def set_clim(self, low, high):
            self.clim = (low, high)

    fake_art3d = types.ModuleType("mpl_toolkits.mplot3d.art3d")
    fake_art3d.Line3DCollection = FakeLine3DCollection
    monkeypatch.setitem(sys.modules, "mpl_toolkits.mplot3d.art3d", fake_art3d)

    collection = make_matplotlib_line3d_collection(pack, linewidth=0.3)
    assert isinstance(collection, FakeLine3DCollection)
    assert collection.segments.shape == pack.segments.shape
    assert np.array_equal(collection.array, pack.segment_scalars)
    assert collection.clim == (pack.scalar_min, pack.scalar_max)

    assert pack.claim_scope == "visualization_render_packing_only"
    assert pack.visualization_ready is False
    assert pack.visual_correspondence_verified is False
    assert pack.pde_validated is False
    assert pack.paper_exact is False
    assert pack.openai_field_identified is False
    assert pack.blowup_proved is False


def test_matlab_surface_bundle_is_nan_separated_and_geometry_preserving():
    lines = [
        np.array([[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [2.0, 0.0, 1.0]]),
        np.array([[0.0, 1.0, 2.0], [1.0, 1.0, 3.0]]),
    ]
    scalars = [np.array([10.0, 20.0, 30.0]), np.array([40.0, 50.0])]
    pack = pack_streamline_collection(lines, scalars)

    assert pack.matlab_x.shape == (2, 6)
    assert np.allclose(pack.matlab_x[0, :3], [0.0, 1.0, 2.0])
    assert np.isnan(pack.matlab_x[0, 3])
    assert np.allclose(pack.matlab_x[0, 4:], [0.0, 1.0])
    assert np.array_equal(pack.matlab_x[0], pack.matlab_x[1], equal_nan=True)
    assert np.array_equal(pack.matlab_c[0], pack.matlab_c[1], equal_nan=True)
    assert np.allclose(pack.matlab_c[0, :3], scalars[0])
    assert np.isnan(pack.matlab_c[0, 3])
    assert np.allclose(pack.matlab_c[0, 4:], scalars[1])

    assert np.allclose(pack.segments[0], [lines[0][0], lines[0][1]])
    assert np.allclose(pack.segments[1], [lines[0][1], lines[0][2]])
    assert np.allclose(pack.segments[2], [lines[1][0], lines[1][1]])
    assert np.allclose(pack.segment_scalars, [15.0, 25.0, 45.0])


def test_height_default_and_fail_closed_inputs():
    line = np.array([[0.0, 0.0, -0.5], [0.2, 0.1, 0.0], [0.4, 0.3, 0.5]])
    pack = pack_streamline_collection([line])
    assert np.allclose(pack.segment_scalars, [-0.25, 0.25])

    with pytest.raises(ValueError):
        pack_streamline_collection([])
    with pytest.raises(ValueError):
        pack_streamline_collection([np.zeros((1, 3))])
    with pytest.raises(ValueError):
        pack_streamline_collection([np.zeros((3, 3))])
    with pytest.raises(ValueError):
        pack_streamline_collection([np.array([[0.0, 0.0, 0.0], [np.nan, 0.0, 1.0]])])
    with pytest.raises(ValueError):
        pack_streamline_collection([line], [np.array([1.0, 2.0])])
    with pytest.raises(ValueError):
        pack_streamline_collection([line], [np.array([1.0, np.inf, 3.0])])

    pack = pack_streamline_collection([line])
    with pytest.raises(ValueError):
        make_matplotlib_line3d_collection(pack, linewidth=0.0)
    with pytest.raises(ValueError):
        make_matplotlib_line3d_collection(pack, cmap="")
