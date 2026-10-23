import numpy as np
import pytest

from openai_ns_reconstruction.constrained_pathlines import integrate_pathlines


def test_constant_time_dependent_translation_matches_closed_form():
    def velocity(points, t):
        out = np.empty_like(points)
        out[:, 0] = 1.0
        out[:, 1] = 2.0 * t
        out[:, 2] = 0.0
        return out

    seeds = np.array([[0.0, 0.0, 0.0], [1.0, -2.0, 3.0]])
    times = np.linspace(0.25, 0.75, 21)

    bundle = integrate_pathlines(
        velocity,
        seeds,
        times,
        allowed_time_interval=(0.25, 0.75),
        rtol=1e-10,
        atol=1e-12,
        max_step=0.01,
    )

    dt = times - times[0]
    expected = np.empty_like(bundle.positions)
    expected[:, :, 0] = seeds[None, :, 0] + dt[:, None]
    expected[:, :, 1] = seeds[None, :, 1] + (times**2 - times[0] ** 2)[:, None]
    expected[:, :, 2] = seeds[None, :, 2]

    np.testing.assert_allclose(bundle.positions, expected, rtol=0, atol=2e-11)
    assert bundle.claim_scope == "visualization_kinematics_only"
    assert bundle.pde_validated is False
    assert bundle.paper_exact is False
    assert bundle.openai_field_identified is False
    assert bundle.blowup_proved is False
    assert bundle.positions.flags.writeable is False


def test_solid_rotation_preserves_radius_and_tracks_angle():
    omega = 3.0

    def velocity(points, t):
        del t
        out = np.zeros_like(points)
        out[:, 0] = -omega * points[:, 1]
        out[:, 1] = omega * points[:, 0]
        return out

    seeds = np.array([[1.0, 0.0, 0.2], [0.0, 0.5, -0.4]])
    times = np.linspace(0.0, 1.0, 51)
    bundle = integrate_pathlines(
        velocity, seeds, times, rtol=1e-10, atol=1e-12, max_step=0.01
    )

    radii = np.linalg.norm(bundle.positions[:, :, :2], axis=2)
    expected_radii = np.linalg.norm(seeds[:, :2], axis=1)
    np.testing.assert_allclose(
        radii, np.broadcast_to(expected_radii, radii.shape), atol=2e-10, rtol=0
    )

    expected0 = np.column_stack([np.cos(omega * times), np.sin(omega * times)])
    np.testing.assert_allclose(bundle.positions[:, 0, :2], expected0, atol=2e-10, rtol=0)


def test_fail_closed_inputs_and_velocity_output():
    good_times = np.array([0.25, 0.5, 0.75])
    seeds = np.array([[0.0, 0.0, 0.0]])

    with pytest.raises(ValueError):
        integrate_pathlines(lambda p, t: p, seeds, [0.25, 0.25])
    with pytest.raises(ValueError):
        integrate_pathlines(
            lambda p, t: p,
            seeds,
            good_times,
            allowed_time_interval=(0.3, 0.7),
        )
    with pytest.raises(ValueError):
        integrate_pathlines(lambda p, t: np.zeros((2, 3)), seeds, good_times)
    with pytest.raises(ValueError):
        integrate_pathlines(lambda p, t: np.full_like(p, np.nan), seeds, good_times)
