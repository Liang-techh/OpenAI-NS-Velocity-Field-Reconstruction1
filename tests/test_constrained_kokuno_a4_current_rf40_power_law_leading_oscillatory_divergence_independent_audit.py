from __future__ import annotations

import copy
import inspect
import math

import numpy as np
import pytest

from openai_ns_reconstruction import (
    kokuno_a4_current_rf40_power_law_leading_oscillatory_divergence_independent_audit as mod,
)


class _LinearField:
    def __init__(self, matrix: np.ndarray):
        self.matrix = np.asarray(matrix, float)

    def velocity(self, x, y, z, t):
        x, y, z, t = np.broadcast_arrays(
            np.asarray(x, float),
            np.asarray(y, float),
            np.asarray(z, float),
            np.asarray(t, float),
        )
        p = np.stack((x, y, z), axis=-1)
        return np.einsum("ij,...j->...i", self.matrix, p) * (1.0 + 0.1 * t[..., None])


class _ProbeGeometry:
    D = 0.5
    X_3 = 2.0
    T_w = 4.0
    X_4 = X_3 * math.exp(T_w)


def test_fixed_protocol_and_truth_boundary():
    assert mod.TASK == "K4-VAL-102"
    assert mod.UPSTREAM_PR == 1010
    assert mod.UPSTREAM_HEAD == "e36d9da4b7f037e998f5b1658f8c0ea291a76b80"
    assert mod.AGENT1_PR == 1005
    assert mod.AGENT1_HEAD == "2c76ebdc41d6c566f43a1305034ba2ff9dce410b"
    assert mod.SEED == 9173741
    assert mod.STEPS == (0.02, 0.01, 0.005)
    assert mod.DIVERGENCE_GATE == 1.0e-5
    assert mod.TRUTH_BOUNDARY["RF40_power_law_scoped_divergence_assessed"] is True
    assert mod.TRUTH_BOUNDARY["leading_only_ns_residual_assessed"] is False
    assert mod.TRUTH_BOUNDARY["leading_plus_oscillatory_ns_residual_assessed"] is False
    assert mod.TRUTH_BOUNDARY["after_correction_ns_residual_assessed"] is False
    assert mod.TRUTH_BOUNDARY["pde_validated"] is False


def test_fd2_jacobian_is_exact_for_linear_public_velocity():
    A = np.asarray(
        [
            [1.2, -0.4, 0.7],
            [0.3, -0.8, 1.1],
            [-0.5, 0.2, -0.4],
        ]
    )
    field = _LinearField(A)
    points = np.asarray(
        [
            [0.3, -0.2, 0.1],
            [-0.4, 0.5, 0.2],
            [0.1, 0.2, -0.3],
        ]
    )
    times = np.asarray([0.31, 0.47, 0.71])
    J = mod.fd2_jacobian(field, points, times, 0.005)
    expected = np.stack([(1.0 + 0.1 * t) * A for t in times], axis=0)
    np.testing.assert_allclose(J, expected, rtol=0.0, atol=2.0e-13)


def test_manufactured_solenoidal_calibration_and_mutation_detection():
    points = np.asarray(
        [
            [0.31, -0.17, 0.08],
            [-0.42, 0.21, -0.11],
            [0.19, 0.33, 0.14],
        ]
    )
    times = np.asarray([0.31, 0.47, 0.71])
    field = mod._ManufacturedSolenoidal()
    J = mod.fd2_jacobian(field, points, times, mod.STEPS[-1])
    div = np.trace(J, axis1=1, axis2=2)
    assert float(np.max(np.abs(div))) <= 1.0e-11

    mutated = mod._DivergenceMutation(field)
    Jm = mod.fd2_jacobian(mutated, points, times, mod.STEPS[-1])
    delta = np.trace(Jm, axis1=1, axis2=2) - div
    np.testing.assert_allclose(
        delta,
        np.full(len(points), mod.MUTATION_EPSILON),
        rtol=0.0,
        atol=2.0e-13,
    )
    assert float(np.max(np.abs(delta))) >= mod.MUTATION_DETECTION_FLOOR


def test_power_law_probe_protocol_is_fresh_weighted_and_strict():
    probes = mod.make_probes(_ProbeGeometry())
    assert probes.inner.points.shape == (32, 3)
    assert probes.power_law.points.shape == (72, 3)
    assert probes.seam_points.shape == (16, 3)
    assert probes.axis_points.shape == (5, 3)
    assert np.all(np.isfinite(probes.power_law.points))
    assert np.all(probes.power_law.weights > 0.0)
    assert np.all(np.isfinite(probes.power_law.weights))
    assert len(np.unique(probes.power_law.group_index)) == 12

    probes2 = mod.make_probes(_ProbeGeometry())
    np.testing.assert_array_equal(probes.power_law.points, probes2.power_law.points)
    np.testing.assert_array_equal(probes.power_law.weights, probes2.power_law.weights)


def test_power_law_weights_include_Tw_change_of_variables():
    a = _ProbeGeometry()
    p = mod.make_probes(a)
    assert float(np.sum(p.power_law.weights)) > 0.0

    class _Longer:
        D = 0.5
        X_3 = 2.0
        T_w = 5.0
        X_4 = X_3 * math.exp(T_w)

    q = mod.make_probes(_Longer())
    assert not math.isclose(
        float(np.sum(p.power_law.weights)),
        float(np.sum(q.power_law.weights)),
        rel_tol=1.0e-8,
        abs_tol=0.0,
    )


def test_staged_attribution_uses_two_independent_black_box_jacobians():
    leading_A = np.diag([1.0, -1.0, 0.0])
    osc_A = np.asarray(
        [
            [0.2, 0.1, 0.0],
            [0.0, -0.1, 0.3],
            [0.4, 0.0, -0.1],
        ]
    )
    leading = _LinearField(leading_A)
    total = _LinearField(leading_A + osc_A)
    points = np.asarray([[0.31, -0.14, 0.09], [-0.27, 0.18, -0.07]])
    times = np.asarray([0.31, 0.63])
    staged, mats = mod._staged(total, leading, points, times, 0.005)
    jt, jl, jo = mats
    np.testing.assert_allclose(jo, jt - jl, rtol=0.0, atol=0.0)
    assert staged["leading"]["sampled_max"] <= 2.0e-13
    assert staged["oscillatory_increment"]["sampled_max"] <= 2.0e-13


def test_recursive_parameter_perturbation_is_fail_closed_and_local():
    payload = {
        "outer": {
            "nested": [{"x": 1.0}, {"eta_fd_step": 2.0e-5}],
            "unchanged": 7,
        }
    }
    before = copy.deepcopy(payload)
    assert mod._recursive_eta_fd_perturb(payload) is True
    assert payload["outer"]["nested"][1]["eta_fd_step"] == pytest.approx(2.002e-5)
    assert payload["outer"]["unchanged"] == before["outer"]["unchanged"]


def test_public_fd_operator_does_not_call_candidate_derivative_helpers():
    text = inspect.getsource(mod.fd2_jacobian)
    assert ".velocity" not in text
    for forbidden in (
        "similarity_radial_derivatives",
        "primitive",
        ".jacobian(",
        ".divergence(",
        ".vorticity(",
        "pressure",
        "forcing",
        "residual",
    ):
        assert forbidden not in text.lower()


def test_gate_rejects_posthoc_threshold_laundering():
    good = {
        "inner": {
            k: {"sampled_max": 0.5e-5, "weighted_rms": 0.5e-5}
            for k in ("total", "leading", "oscillatory_increment")
        }
    }
    mod._stage_gate(good, "inner")
    bad = copy.deepcopy(good)
    bad["inner"]["total"]["sampled_max"] = 1.0001e-5
    with pytest.raises(AssertionError):
        mod._stage_gate(bad, "inner")


def test_source_xyz_rejects_chart_escape():
    field = _ProbeGeometry()
    with pytest.raises(ValueError):
        mod._source_xyz(field, [1.0], [1.0], [0.0], [0.5])
    with pytest.raises(ValueError):
        mod._source_xyz(field, [-1.0], [0.0], [0.0], [0.5])
    with pytest.raises(ValueError):
        mod._source_xyz(field, [1.0], [0.0], [0.0], [1.0])
