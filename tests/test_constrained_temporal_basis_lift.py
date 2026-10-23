import numpy as np
import pytest

from openai_ns_reconstruction.constrained_temporal_basis_lift import (
    diagnose_affine_temporal_lift,
)


def test_temporal_lift_resolves_opposite_time_requirements_for_same_channel():
    times = np.array([0.25, 0.5, 0.75])
    # Two spatial directions; only the swirl-like second direction needs to
    # reverse across time. Static coefficients cannot do that.
    jr = np.zeros((3, 2, 2))
    jv = np.zeros((3, 2, 2))
    jr[:, 0, 0] = 1.0
    jv[:, 0, 0] = 1.0
    jr[:, 1, 1] = 1.0
    jv[:, 1, 1] = 1.0

    tau = np.array([-1.0, 0.0, 1.0])
    rt = np.column_stack((np.ones(3), 0.8 * tau))
    vt = np.column_stack((np.ones(3), 0.6 * tau))

    report = diagnose_affine_temporal_lift(
        times,
        jr,
        rt,
        jv,
        vt,
        residual_scales=[1, 1],
        visual_scales=[1, 1],
        lower_bounds=[-2, -2],
        upper_bounds=[2, 2],
        lifted_parameter_index=1,
        slope_bound=1.0,
        parameter_labels=["poloidal", "swirl"],
        material_gain=0.2,
    )

    assert report.static.residual_remaining_ratio > 0.45
    assert report.static.visual_remaining_ratio > 0.35
    assert report.affine_lift.residual_remaining_ratio < 0.1
    assert report.affine_lift.visual_remaining_ratio < 0.1
    assert report.lift_materially_improves_both
    assert report.lifted_parameter_label == "swirl:time_slope"
    assert not report.pde_validated
    assert not report.visual_correspondence_verified


def test_lift_does_not_help_when_target_is_static():
    times = [0.25, 0.5, 0.75]
    jr = np.ones((3, 1, 1))
    jv = np.ones((3, 1, 1))
    rt = np.ones((3, 1)) * 0.4
    vt = np.ones((3, 1)) * 0.7

    report = diagnose_affine_temporal_lift(
        times,
        jr,
        rt,
        jv,
        vt,
        residual_scales=[1],
        visual_scales=[1],
        lower_bounds=[-1],
        upper_bounds=[1],
        lifted_parameter_index=0,
        slope_bound=1.0,
        parameter_labels=["F0"],
        material_gain=0.05,
    )
    assert report.affine_lift.residual_remaining_ratio == pytest.approx(
        report.static.residual_remaining_ratio, abs=1e-10
    )
    assert report.affine_lift.visual_remaining_ratio == pytest.approx(
        report.static.visual_remaining_ratio, abs=1e-10
    )
    assert not report.lift_materially_improves_both


def test_bound_saturation_is_reported_and_prevents_fake_capacity():
    times = [0.25, 0.5, 0.75]
    jr = np.ones((3, 1, 1))
    jv = np.ones((3, 1, 1))
    tau = np.array([-1.0, 0.0, 1.0])
    rt = (3.0 * tau)[:, None]
    vt = (3.0 * tau)[:, None]

    report = diagnose_affine_temporal_lift(
        times,
        jr,
        rt,
        jv,
        vt,
        residual_scales=[1],
        visual_scales=[1],
        lower_bounds=[-1],
        upper_bounds=[1],
        lifted_parameter_index=0,
        slope_bound=0.25,
        parameter_labels=["swirl"],
        material_gain=0.5,
    )
    assert "swirl:time_slope" in report.affine_lift.active_upper_bounds
    assert report.affine_lift.residual_remaining_ratio > 0.8
    assert not report.lift_materially_improves_both


def test_fail_closed_inputs():
    jr = np.ones((3, 1, 1))
    jv = np.ones((3, 1, 1))
    rt = np.ones((3, 1))
    vt = np.ones((3, 1))
    kwargs = dict(
        residual_scales=[1],
        visual_scales=[1],
        lower_bounds=[-1],
        upper_bounds=[1],
        lifted_parameter_index=0,
        slope_bound=1.0,
    )
    with pytest.raises(ValueError):
        diagnose_affine_temporal_lift(
            [0.5, 0.25, 0.75], jr, rt, jv, vt, **kwargs
        )
    with pytest.raises(ValueError):
        diagnose_affine_temporal_lift(
            [0.25, 0.5, 0.75],
            jr,
            rt,
            jv,
            vt,
            **{**kwargs, "slope_bound": 0.0},
        )
    bad = jr.copy()
    bad[0, 0, 0] = np.nan
    with pytest.raises(ValueError):
        diagnose_affine_temporal_lift(
            [0.25, 0.5, 0.75], bad, rt, jv, vt, **kwargs
        )
