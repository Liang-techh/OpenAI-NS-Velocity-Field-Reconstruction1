import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_joint_gain_margin import (
    connected_sublevel_interval,
    intersect_gain_constraints,
    quartic_relative_rms_sublevel_coefficients,
    relative_rms,
)


def test_quartic_matches_direct_relative_rms_sublevel_value():
    base = np.array([[1.0, -0.5], [0.75, 0.25], [-0.3, 0.8]])
    linear = np.array([[-0.4, 0.1], [-0.2, -0.3], [0.5, -0.1]])
    quadratic = np.array([[0.08, -0.04], [0.02, 0.06], [-0.03, 0.01]])
    bound = 0.92
    coefficients = quartic_relative_rms_sublevel_coefficients(
        base,
        linear,
        quadratic,
        relative_bound=bound,
    )
    polynomial = np.polynomial.Polynomial(coefficients)
    base_norm_sq = float(np.sum(base * base))
    for damping in (0.0, 0.2, 0.51, 0.8, 1.0):
        direct = base + damping * linear + damping * damping * quadratic
        expected = float(np.sum(direct * direct) - bound * bound * base_norm_sq)
        assert float(polynomial(damping)) == pytest.approx(expected, rel=2e-14, abs=2e-14)


def test_connected_sublevel_interval_recovers_known_linear_residual_window():
    # residual = 1-lambda and inherited relative bound = 3/4, hence on [0,1]
    # the sublevel is exactly lambda >= 1/4.
    coefficients = quartic_relative_rms_sublevel_coefficients(
        np.array([1.0]),
        np.array([-1.0]),
        np.array([0.0]),
        relative_bound=0.75,
    )
    result = connected_sublevel_interval(coefficients, point=0.5)
    assert result["sublevel_interval"][0] == pytest.approx(0.25, abs=1e-10)
    assert result["sublevel_interval"][1] == pytest.approx(1.0, abs=1e-12)
    assert result["left_margin"] == pytest.approx(0.25, abs=1e-10)
    assert result["right_margin"] == pytest.approx(0.5, abs=1e-12)


def test_connected_sublevel_interval_handles_quartic_two_sided_window():
    # residual = 1 - 2 lambda + 2 lambda^2.  Relative bound 0.8 gives a
    # nontrivial connected interval around lambda=0.5, bounded on both sides.
    base = np.array([1.0])
    linear = np.array([-2.0])
    quadratic = np.array([2.0])
    coefficients = quartic_relative_rms_sublevel_coefficients(
        base, linear, quadratic, relative_bound=0.8
    )
    result = connected_sublevel_interval(coefficients, point=0.5)
    left, right = result["sublevel_interval"]
    assert 0.0 < left < 0.5 < right < 1.0
    assert relative_rms(base, linear, quadratic, left) == pytest.approx(0.8, abs=2e-9)
    assert relative_rms(base, linear, quadratic, right) == pytest.approx(0.8, abs=2e-9)


def test_intersection_reports_binding_constraints_without_reoptimizing_point():
    # First constraint admits [0.2,1]; second admits a symmetric interior window.
    first = quartic_relative_rms_sublevel_coefficients(
        np.array([1.0]), np.array([-1.0]), np.array([0.0]), relative_bound=0.8
    )
    second = quartic_relative_rms_sublevel_coefficients(
        np.array([1.0]), np.array([-2.0]), np.array([2.0]), relative_bound=0.75
    )
    result = intersect_gain_constraints(
        [
            {"name": "first", "quartic_coefficients": first.tolist()},
            {"name": "second", "quartic_coefficients": second.tolist()},
        ],
        frozen_damping=0.5,
    )
    left, right = result["joint_admissible_interval"]
    assert left < 0.5 < right
    assert result["frozen_damping"] == 0.5
    assert result["validation_used_for_damping_selection"] is False
    assert result["frozen_damping_changed"] is False
    assert result["frozen_damping_strictly_inside_joint_interval"] is True
    assert result["minimum_two_sided_margin"] > 0.0
    assert "second" in result["binding_right_constraints"]


def test_sublevel_fails_closed_when_frozen_damping_is_not_admissible():
    coefficients = quartic_relative_rms_sublevel_coefficients(
        np.array([1.0]), np.array([-1.0]), np.array([0.0]), relative_bound=0.5
    )
    with pytest.raises(ValueError, match="outside this gain sublevel"):
        connected_sublevel_interval(coefficients, point=0.1)


def test_quartic_rejects_zero_base_and_shape_mismatch():
    with pytest.raises(ValueError, match="base norm"):
        quartic_relative_rms_sublevel_coefficients(
            np.zeros(3), np.ones(3), np.ones(3), relative_bound=1.0
        )
    with pytest.raises(ValueError, match="common nonempty shape"):
        quartic_relative_rms_sublevel_coefficients(
            np.ones(3), np.ones(4), np.ones(3), relative_bound=1.0
        )
