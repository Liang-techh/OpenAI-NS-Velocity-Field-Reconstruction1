import numpy as np

import agent7_st051b_callable_piola_transfer as a7
import replay_st051


def test_frozen_contract_and_truth_boundary():
    assert a7.PREREG_ISSUE == 509
    assert a7.PARENT_ID == "ST051-B"
    assert a7.SOURCE_REDISTRIBUTION_GAIN == 0.025
    assert a7.PIOLA_BETA == 0.075
    assert not any(a7.TRUTH.values())
    assert a7.CRITERIA["axial_vorticity_rms_gain_min"] == 0.01
    assert a7.CRITERIA["response_condition_max"] == 5.0


def test_piola_coordinate_map_is_orientation_preserving_and_endpoint_fixed():
    z = np.linspace(-2.0, 2.0, 4001)
    mapped, jac = a7.warp_z_and_jacobian(z, a7.PIOLA_BETA)
    assert np.all(jac > 0.0)
    assert np.all(np.diff(mapped) > 0.0)
    assert mapped[0] == -2.0
    assert mapped[-1] == 2.0
    assert np.min(jac) >= 0.924999999999
    assert np.max(jac) < 1.05


def test_frozen_redistribution_has_inner_mid_gain_and_outer_compensation():
    vals = a7.h_profile(np.array([0.6, 0.9, 1.2]))
    assert vals[0] > 0.0
    assert vals[1] > 0.0
    assert vals[2] < 0.0


def test_beta_zero_callable_replay_is_exact_on_reconstructed_st051b():
    f, raw = replay_st051.reconstruct("ST051-B")
    points = np.array(
        [[0.31, 0.17, 0.22], [0.72, -0.18, -0.61], [1.18, 0.29, 0.83]],
        dtype=float,
    )
    base = a7.redistributed_velocity(
        f, raw, points, 0.5,
        gain=a7.SOURCE_REDISTRIBUTION_GAIN,
        scale=1.0,
    )
    replay = a7.piola_velocity(
        f, raw, points, 0.5,
        redistribution_gain=a7.SOURCE_REDISTRIBUTION_GAIN,
        redistribution_scale=1.0,
        beta=0.0,
        piola_scale=1.0,
    )
    assert np.array_equal(base, replay)
    assert np.isfinite(replay).all()
    assert np.linalg.norm(replay) > 0.0
