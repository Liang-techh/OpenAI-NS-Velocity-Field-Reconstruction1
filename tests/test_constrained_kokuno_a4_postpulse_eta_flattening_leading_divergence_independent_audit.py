from __future__ import annotations

import copy
import inspect

import numpy as np
import pytest

from openai_ns_reconstruction import (
    kokuno_a4_postpulse_eta_flattening_leading_divergence_independent_audit as m,
)
from openai_ns_reconstruction.kokuno_pa16_current_cartesian_postpulse_eta_flattening import (
    KokunoPA16CurrentCartesianPostPulseEtaFlattening,
)


def test_preregistered_identity_and_thresholds_are_frozen():
    assert m.TASK == "K4-VAL-118"
    assert m.UPSTREAM_HEAD == "c5442c11165d7f0976889b826bf58dce38a3d3dd"
    assert m.SEED == 9173901
    assert m.SPATIAL_STEPS == (0.02, 0.01, 0.005)
    assert m.DIVERGENCE_GATE == 1e-5
    assert m.FINAL_PROJECT_MOMENTUM_GATE == 1e-3
    assert m.FINAL_PROJECT_DIVERGENCE_GATE == 1e-5


def test_heldout_postpulse_probe_geometry_is_deterministic_and_finite():
    field = m.default_field()
    a = m.make_heldout_probes(field)
    b = m.make_heldout_probes(field)
    assert a.points.shape == (60, 3)
    assert a.seam_points.shape == (16, 3)
    assert a.axis_points.shape == (6, 3)
    assert a.late_points.shape == (3, 3)
    assert np.array_equal(a.points, b.points)
    assert np.array_equal(a.times, b.times)
    assert np.array_equal(a.region, b.region)
    assert np.all(np.isfinite(a.points))
    assert np.all(np.isfinite(a.log_weights))
    assert np.allclose(a.late_s, (0.92, 0.96, 0.985))


def test_save_reload_is_identity_bound_and_provenance_mutation_fails_closed(tmp_path):
    field = m.default_field()
    path = tmp_path / "candidate.json"
    payload = field.save_configuration(path)
    loaded = KokunoPA16CurrentCartesianPostPulseEtaFlattening.load_configuration(path)
    assert loaded.semantic_sha256 == field.semantic_sha256

    mutated = copy.deepcopy(payload)
    mutated["bound_scope"]["T_f_role"] = "source_exact"
    with pytest.raises(ValueError):
        KokunoPA16CurrentCartesianPostPulseEtaFlattening.from_configuration(mutated)


def test_post_stage_fails_closed():
    field = m.default_field()
    with pytest.raises(ValueError):
        field.similarity_profile_values_logX(field.log_X_flatten_end + 0.02, 0.0)


class _Solenoidal:
    def velocity(self, x, y, z, t):
        x, y, z, t = np.broadcast_arrays(
            np.asarray(x, float),
            np.asarray(y, float),
            np.asarray(z, float),
            np.asarray(t, float),
        )
        return np.stack(((1.0 + t) * y, -(1.0 + t) * x, np.zeros_like(z)), axis=-1)


class _Corrupt:
    def __init__(self, base):
        self.base = base

    def velocity(self, x, y, z, t):
        out = np.asarray(self.base.velocity(x, y, z, t), float).copy()
        out[..., 0] += 1e-3 * np.asarray(x, float)
        return out


def test_independent_fd2_calibration_and_mutation_detection():
    points = np.asarray(((0.2, -0.3, 0.1), (-0.4, 0.1, -0.2), (0.6, 0.2, 0.3)))
    times = np.asarray((0.31, 0.47, 0.71))
    for h in m.SPATIAL_STEPS:
        jac = m.independent_fd2_jacobian(_Solenoidal(), points, times, h)
        div = np.trace(jac, axis1=1, axis2=2)
        assert np.max(np.abs(div)) <= 1e-12
    jac = m.independent_fd2_jacobian(_Corrupt(_Solenoidal()), points, times, 0.005)
    div = np.trace(jac, axis1=1, axis2=2)
    assert np.min(np.abs(div)) >= 5e-4


def test_binary64_representability_firewall_rejects_collapsed_absolute_step():
    ordinary = np.asarray(((0.2, -0.3, 0.1), (1.0, 2.0, -3.0)))
    ok, receipt = m.fd_representability(ordinary, 0.005)
    assert np.all(ok)
    assert receipt["all_probes_representable"] is True

    huge = np.asarray(((1.0e200, 0.0, 0.0),))
    ok, receipt = m.fd_representability(huge, 0.005)
    assert not np.all(ok)
    assert receipt["all_probes_representable"] is False
    assert receipt["per_axis"][0]["nonrepresentable_count"] == 1


def _passing_report():
    keys = (
        "fine_sampled_max_pass",
        "fine_weighted_rms_pass",
        "pulse_end_seam_pass",
        "axis_near_pass",
        "late_stage_pass",
        "nontrivial_speed_pass",
        "medium_to_fine_sampled_max_stable",
        "medium_to_fine_weighted_rms_stable",
        "all_canonical_fd_perturbations_representable",
        "save_reload_replay_pass",
        "manufactured_solenoidal_calibration_pass",
        "mutation_detection_pass",
        "ordering_invariance_pass",
        "configuration_mutation_rejected",
        "post_stage_fail_closed",
        "public_api_tuning_knobs_absent",
    )
    gates = {key: True for key in keys}
    gates.update(
        divergence_gate=1e-5,
        final_project_momentum_gate_unchanged=1e-3,
        final_project_divergence_gate_unchanged=1e-5,
    )
    return {
        "schema": m.SCHEMA,
        "task": m.TASK,
        "upstream_head": m.UPSTREAM_HEAD,
        "gates": gates,
        "truth_boundary": {
            "leading_only_ns_residual_assessed": False,
            "leading_plus_oscillatory_ns_residual_assessed": False,
            "after_correction_ns_residual_assessed": False,
            "pde_validated": False,
        },
    }


def test_gate_enforcer_accepts_only_frozen_contract():
    report = _passing_report()
    m.enforce_preregistered_gates(report)
    bad = copy.deepcopy(report)
    bad["gates"]["divergence_gate"] = 1e-4
    with pytest.raises(AssertionError):
        m.enforce_preregistered_gates(bad)
    bad = copy.deepcopy(report)
    bad["truth_boundary"]["pde_validated"] = True
    with pytest.raises(AssertionError):
        m.enforce_preregistered_gates(bad)


def test_public_scientific_entrypoints_expose_no_tuning_knobs():
    forbidden = ("threshold", "forcing", "pressure", "viscosity", "step", "residual", "target")
    for fn in (m.materialize_receipt, m.enforce_preregistered_gates):
        params = tuple(inspect.signature(fn).parameters)
        assert not any(any(token in name.lower() for token in forbidden) for name in params)
