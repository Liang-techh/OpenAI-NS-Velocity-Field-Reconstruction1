import numpy as np

import agent7_st048s_piola_warp_screen as screen


class _LinearDivergenceFreeField:
    def at_points(self, points, time):
        points = np.asarray(points, dtype=float)
        x, y, _z = points.T
        return np.column_stack((x, -y, np.zeros_like(x)))


def _fd_divergence(field, points, time=0.5, h=1e-6):
    points = np.asarray(points, dtype=float)
    out = np.zeros(len(points))
    for axis in range(3):
        delta = np.zeros_like(points)
        delta[:, axis] = h
        plus = field.at_points(points + delta, time)
        minus = field.at_points(points - delta, time)
        out += (plus[:, axis] - minus[:, axis]) / (2 * h)
    return out


def test_piola_warp_beta_zero_replay_and_divergence_identity():
    parent = _LinearDivergenceFreeField()
    points = np.array([[.2, .3, .1], [.7, -.4, .8], [-.5, .6, -1.1]])
    zero = screen._PiolaWarpField(parent, 0.0)
    assert np.array_equal(zero.at_points(points, .5), parent.at_points(points, .5))

    warped = screen._PiolaWarpField(parent, .15)
    assert np.max(np.abs(_fd_divergence(warped, points))) < 2e-9
    mapped, jac = screen._warp_z_and_jacobian(np.linspace(-2, 2, 2001), .2)
    assert np.all(np.diff(mapped) > 0)
    assert np.min(jac) > 0
    assert mapped[0] == -2 and mapped[-1] == 2


def test_frozen_st048s_warp_screen_smoke_contract():
    report = screen.audit_st048s_piola_warp_screen(
        betas=(0.0, .05),
        times=(.5,),
        energy_times=(.25, .75),
        energy_orders=(16, 20),
        grid_size=17,
    )
    assert report["task_id"] == screen.TASK_ID
    assert report["parent_id"] == "ST048-S"
    assert report["upstream_st048_residual_receipt"]["paired_results"]
    assert len(report["rows"]) == 2
    assert report["rows"][0]["beta"] == 0.0
    assert report["rows"][1]["parameter_count"] == 1
    assert report["truth_boundary"]["new_spatial_basis_added"] is False
    assert report["truth_boundary"]["held_out_pde_residual_evaluated"] is False
    assert "not evaluated" in report["held_out_full_momentum"]
    assert np.isfinite(report["energy_quadrature_max_relative_refinement_change"])
