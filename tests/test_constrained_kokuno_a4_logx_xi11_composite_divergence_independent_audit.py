from __future__ import annotations

import inspect
import math

import numpy as np
import pytest

from openai_ns_reconstruction import (
    kokuno_a4_logx_xi11_composite_divergence_independent_audit as m,
)


class _Leading:
    D = 1.0
    lambda_value = 0.05
    log_X_p = 2.0


class _Field:
    leading_backend = _Leading()


class _LinearSolenoidal:
    def velocity(self, x, y, z, t):
        x, y, z, t = np.broadcast_arrays(
            np.asarray(x, float),
            np.asarray(y, float),
            np.asarray(z, float),
            np.asarray(t, float),
        )
        return np.stack((y, -x, np.zeros_like(z)), axis=-1)


def test_frozen_protocol_and_truth_boundary():
    assert m.TASK == "K4-VAL-117"
    assert m.UPSTREAM_PR == 1117
    assert m.UPSTREAM_HEAD == "27741d9c0a27262f7fabf61eebaa2fbd507e9f03"
    assert m.AGENT1_PR == 1107
    assert m.SPATIAL_STEPS == (0.02, 0.01, 0.005)
    assert m.DIVERGENCE_GATE == 1e-5
    assert m.SEED == 9173891
    assert m._TRUTH_BOUNDARY["pde_validated"] is False
    assert m._TRUTH_BOUNDARY["leading_plus_oscillatory_ns_residual_assessed"] is False
    assert m._TRUTH_BOUNDARY["newer_agent1_1133_pulse_end_leading_consumed"] is False


def test_logx_cartesian_map_avoids_explicit_X_overflow():
    field = _Field()
    log_x = np.asarray((10.0, 400.0))
    eta = np.asarray((0.0, 0.2))
    theta = np.asarray((0.0, 1.0))
    time = np.asarray((0.5, 0.5))
    xyz = m._source_to_cartesian_logx(field, log_x, eta, theta, time)
    assert xyz.shape == (2, 3)
    assert np.all(np.isfinite(xyz))
    assert xyz[1, 0] != 0.0 or xyz[1, 1] != 0.0


def test_log_weighted_metrics_remain_stable_when_volume_overflows():
    values = np.asarray((1e-6, 2e-6, 3e-6))
    base = np.asarray((0.0, 1.0, 2.0))
    a = m._stable_weighted_metrics(values, base)
    b = m._stable_weighted_metrics(values, base + 1000.0)
    assert a["weighted_rms"] == pytest.approx(b["weighted_rms"], rel=1e-14)
    assert b["estimated_volume"] is None
    assert b["volume_l2_estimate"] is None
    assert math.isfinite(b["log_estimated_volume"])
    assert math.isfinite(b["log_volume_l2_estimate"])


def test_fd_representability_accepts_normal_coordinates_and_rejects_huge_tail():
    ordinary = np.asarray(((0.2, -0.3, 0.1), (1.0, 2.0, -1.0)))
    mask, report = m._fd_representability(ordinary, 0.005)
    assert np.all(mask)
    assert report["all_probes_representable"] is True

    huge = np.asarray(((1e200, 0.2, 0.1),))
    mask, report = m._fd_representability(huge, 0.005)
    assert not np.all(mask)
    assert report["all_probes_representable"] is False
    assert report["per_axis"][0]["nonrepresentable_count"] == 1


def test_independent_fd2_recovers_zero_divergence_on_linear_field():
    points = np.asarray(((0.2, -0.3, 0.1), (-0.4, 0.1, -0.2)))
    times = np.asarray((0.31, 0.63))
    jac = m.independent_fd2_jacobian(
        _LinearSolenoidal(), points, times, m.SPATIAL_STEPS[-1]
    )
    div = np.trace(jac, axis1=1, axis2=2)
    assert np.max(np.abs(div)) <= 1e-12
    assert m._manufactured_calibration() is True


def _synthetic_metrics(*, representable: bool = True):
    rep = {"all_probes_representable": representable}
    return {
        "sampled_max": 1e-7,
        "pooled_weighted_rms": 2e-7,
        "pulse_entry_seam_sampled_max": 1e-7,
        "axis_axis_near_sampled_max": 1e-7,
        "late_tail_sampled_max": 1e-7,
        "speed_rms": 1.0,
        "main_fd_representability": rep,
        "pulse_entry_seam_fd_representability": rep,
        "axis_fd_representability": rep,
        "late_tail_fd_representability": rep,
    }


def _synthetic_receipt(*, representable: bool = True):
    metric = _synthetic_metrics(representable=representable)
    return {
        "resolutions": [
            {"leading_plus_oscillatory": metric},
            {"leading_plus_oscillatory": metric},
            {"leading_plus_oscillatory": metric},
        ],
        "inner_oscillatory_signal_rms": 1e-5,
        "checks": {
            "all_resolutions_fd_representable": representable,
            "save_reload_semantic_identity": True,
            "save_reload_velocity_replay_max_abs": 0.0,
            "velocity_mutation_detected": True,
            "offgrid_order_invariant_all_components": True,
            "rehash_consistent_leading_mutation_rejected": True,
            "rehash_consistent_oscillatory_mutation_rejected": True,
            "post_xi11_fails_closed": True,
            "manufactured_solenoidal_calibration_pass": True,
            "source_parameter_perturbation": {"all_finite": True},
        },
    }


def test_scientific_gate_fails_closed_on_fd_roundoff_collapse():
    assert m._derive_audit_pass(_synthetic_receipt()) is True
    assert m._derive_audit_pass(
        _synthetic_receipt(representable=False)
    ) is False


def test_public_scientific_entrypoints_have_no_tuning_knobs():
    forbidden = {
        "threshold",
        "residual",
        "forcing",
        "pressure",
        "viscosity",
        "gain",
        "step",
        "steps",
        "rtol",
        "atol",
        "seed",
    }
    for fn in (
        m.materialize_receipt,
        m.enforce_preregistered_gates,
    ):
        assert not (set(inspect.signature(fn).parameters) & forbidden)


def test_logx_tail_definition_is_strictly_inside_xi11():
    assert all(0.0 < lo < hi < m.SOURCE_XI_END for lo, hi in m.XI_BANDS)
    assert max(m.TAIL_XI) < m.SOURCE_XI_END
    assert m.PRE_PULSE_XI_BAND[1] < 0.0
