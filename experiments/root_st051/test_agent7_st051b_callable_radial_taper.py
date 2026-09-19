import numpy as np

import agent7_st051b_callable_radial_taper as a7
import replay_st051


def test_frozen_contract_and_truth_boundary():
    assert a7.PREREG_ISSUE == 519
    assert a7.TAPER_TAU == 0.05
    assert not any(a7.TRUTH.values())
    assert a7.CRITERIA["tip_radial_rms_reduction_min"] == 0.01
    assert a7.CRITERIA["response_condition_max"] == 5.0
    assert a7.CRITERIA["response_abs_cosine_max"] == 0.90


def test_taper_profile_is_identity_at_midplane_and_support_and_unit_normalized():
    z_star = 2.0 / np.sqrt(3.0)
    z = np.array([-2.0, 0.0, 2.0, z_star], dtype=float)
    q, dq = a7.taper_profile_and_derivative(z)
    assert q[0] == 0.0
    assert q[1] == 0.0
    assert q[2] == 0.0
    assert abs(q[3] - 1.0) < 2.0e-15
    assert np.isfinite(dq).all()

    grid = np.linspace(-2.0, 2.0, 8001)
    qg, _ = a7.taper_profile_and_derivative(grid)
    radial_scale = 1.0 + a7.TAPER_TAU * qg
    assert np.min(radial_scale) >= 1.0
    assert np.max(radial_scale) <= 1.0 + a7.TAPER_TAU + 1.0e-12
    assert radial_scale[0] == 1.0
    assert radial_scale[len(radial_scale) // 2] == 1.0
    assert radial_scale[-1] == 1.0
    assert np.min(radial_scale * radial_scale) > 0.0


def test_tau_zero_callable_replay_is_exact_on_reconstructed_st051b():
    f, raw = replay_st051.reconstruct("ST051-B")
    points = np.array(
        [[0.31, 0.17, 0.22], [0.72, -0.18, -0.61], [1.18, 0.29, 0.83]],
        dtype=float,
    )
    control = a7.base.redistributed_velocity(
        f,
        raw,
        points,
        0.5,
        gain=a7.base.SOURCE_REDISTRIBUTION_GAIN,
        scale=1.0,
    )
    replay = a7.taper_velocity(
        f,
        raw,
        points,
        0.5,
        redistribution_scale=1.0,
        tau=0.0,
        taper_scale=1.0,
    )
    assert np.array_equal(control, replay)
    assert np.isfinite(replay).all()
    assert np.linalg.norm(replay) > 0.0


def test_taper_map_is_identity_on_midplane_and_support_endpoints():
    points = np.array(
        [[0.6, -0.2, 0.0], [0.6, -0.2, 2.0], [0.6, -0.2, -2.0]],
        dtype=float,
    )
    mapped, a, ap = a7.taper_map(points, a7.TAPER_TAU)
    assert np.array_equal(mapped, points)
    assert np.array_equal(a, np.ones(3))
    assert np.array_equal(ap, np.zeros(3))
