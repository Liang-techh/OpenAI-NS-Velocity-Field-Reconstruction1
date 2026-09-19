import numpy as np

import agent7_st052m_tip_odd_poloidal_screen as screen


def _radial_component(points, velocity):
    points = np.asarray(points, dtype=float)
    velocity = np.asarray(velocity, dtype=float)
    r = np.hypot(points[:, 0], points[:, 1])
    return (points[:, 0] * velocity[:, 0] + points[:, 1] * velocity[:, 1]) / r


def test_frozen_activation_reuses_linear_time_law():
    assert screen.activation(0.25) == 0.0
    assert screen.activation(0.50) == 0.5
    assert screen.activation(0.75) == 1.0


def test_tip_orientation_is_inward_on_both_axial_signs():
    pts = np.array(
        [
            [0.6, 0.0, 1.15],
            [0.9, 0.0, 1.15],
            [1.2, 0.0, 1.15],
            [0.6, 0.0, -1.15],
            [0.9, 0.0, -1.15],
            [1.2, 0.0, -1.15],
        ],
        dtype=float,
    )
    correction = screen.tip_odd_poloidal_correction(pts)
    ur = _radial_component(pts, correction)
    assert np.all(ur < 0.0)


def test_shoulder_and_declared_exterior_are_exactly_zero():
    shoulder = np.array(
        [[r, 0.0, z] for z in (-0.85, 0.85) for r in (0.6, 0.9, 1.2)],
        dtype=float,
    )
    outside = np.array(
        [
            [1.61, 0.0, 1.2],
            [0.9, 0.0, 0.99],
            [0.9, 0.0, -0.99],
            [0.9, 0.0, 1.81],
            [0.9, 0.0, -1.81],
        ],
        dtype=float,
    )
    assert np.max(np.abs(screen.tip_odd_poloidal_correction(shoulder))) == 0.0
    assert np.max(np.abs(screen.tip_odd_poloidal_correction(outside))) == 0.0


def test_analytic_curl_is_numerically_divergence_free():
    rng = np.random.default_rng(9177003)
    pts = np.column_stack(
        (
            rng.uniform(-1.45, 1.45, 40),
            rng.uniform(-1.45, 1.45, 40),
            rng.uniform(-1.75, 1.75, 40),
        )
    )
    h = 1.0e-5
    div = np.zeros(len(pts), dtype=float)
    for axis in range(3):
        shift = np.zeros(3, dtype=float)
        shift[axis] = h
        up = screen.tip_odd_poloidal_correction(pts + shift)
        um = screen.tip_odd_poloidal_correction(pts - shift)
        div += (up[:, axis] - um[:, axis]) / (2.0 * h)
    assert np.max(np.abs(div)) < 2.0e-7


def test_truth_boundary_and_single_new_channel_are_frozen():
    assert screen.PREREG_ISSUE == 700
    assert screen.SOURCE_PARENT_PR == 692
    assert screen.SOURCE_PARENT_HEAD == "e3262d4c03dc2d0806eb46639a6e14b2edb5fcb7"
    assert screen.TRUTH["screened_extra_tip_odd_poloidal_basis"] is True
    assert screen.TRUTH["new_spatial_basis_promoted"] is False
    assert screen.TRUTH["new_temporal_basis_added"] is False
    assert screen.TRUTH["parameter_grid_scan_performed"] is False
    assert screen.TRUTH["optimization_performed"] is False
    assert screen.TRUTH["held_out_pde_residual_evaluated"] is False
    assert screen.TRUTH["visual_correspondence_verified"] is False
    assert screen.TRUTH["pde_validated"] is False
    assert screen.TRUTH["openai_field_identified"] is False
