from __future__ import annotations

import copy
import inspect

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_a4_current_rf40_power_law_nonlinear_mean_independent_audit import (
    AGENT1_HEAD,
    AGENT1_PR,
    AGENT2_HEAD,
    AGENT2_PR,
    FROZEN_ANGULAR_ORDERS,
    FROZEN_SEED,
    FROZEN_SPATIAL_STEPS,
    PARENT_BACKEND_KIND,
    PIECE_CLOSURE_RELATIVE_GATE,
    RingSpec,
    _validate_a3_truth_boundary,
    _validate_receipt,
    enforce_scoped_gates,
    frozen_axis_near_specs,
    frozen_inner_ring_specs,
    independent_ring_means,
    scoped_gates_pass,
    truth_boundary,
)


class _Leading:
    def velocity(self, x, y, z, t):
        x, y, z, t = np.broadcast_arrays(x, y, z, t)
        return np.stack((x, -y, np.zeros_like(z)), axis=-1)


class _Total:
    def velocity(self, x, y, z, t):
        x, y, z, t = np.broadcast_arrays(x, y, z, t)
        # leading [x,-y,0] + solid rotation [-y,x,0]
        return np.stack((x - y, -y + x, np.zeros_like(z)), axis=-1)


def _production_receipt(spec: RingSpec) -> dict:
    zeros = [[0.0, 0.0, 0.0]]
    return {
        "backend_kind": PARENT_BACKEND_KIND,
        "radius": [spec.radius],
        "z": [spec.z],
        "t": [spec.t],
        "mean_inner_advects_oscillation_cylindrical": zeros,
        "mean_oscillation_advects_inner_cylindrical": zeros,
        "mean_mixed_cross_cylindrical": zeros,
        "mean_oscillatory_self_advection_cylindrical": zeros,
        "mean_aggregate_nonlinear_cylindrical": zeros,
        "pointwise_three_piece_closure_absolute_max": 0.0,
        "projected_three_piece_closure_absolute_max": 0.0,
        "agent2_backend": {
            "agent2_composite_pr": AGENT2_PR,
            "agent2_composite_head": AGENT2_HEAD,
            "agent1_leading_pr": AGENT1_PR,
            "agent1_leading_head": AGENT1_HEAD,
        },
    }


def _a3_boundary() -> dict:
    return {
        "current_leading_plus_oscillation_through_RF40_power_law_consumed": True,
        "RF40_power_law_velocity_materialized": True,
        "current_RF40_power_law_mixed_nonlinear_mean_materialized": True,
        "current_RF40_power_law_quadratic_nonlinear_mean_materialized": True,
        "current_RF40_power_law_aggregate_nonlinear_mean_materialized": True,
        "velocity_after_RF40_power_law_materialized": False,
        "full_post_XR_RF40_current_lineage_materialized": False,
        "outer_global_leading_velocity_materialized": False,
        "radial_inverse_performed_in_this_increment": False,
        "pressure_gradient_included": False,
        "restricted_forcing_included": False,
        "complete_ns_defect": False,
        "scoped_nonlinear_mean_authorized_as_correction_target": False,
        "mean_correction_velocity_materialized": False,
        "real_candidate_finite_correction_cycle_run": False,
        "heldout_normalized_ns_residual_assessed": False,
        "residual_reduction_claimed": False,
        "same_protocol_comparable_to_st006": False,
        "pde_validated": False,
        "final_normalized_momentum_gate": 1e-3,
        "final_normalized_divergence_gate": 1e-5,
    }


def _passing_report() -> dict:
    fine = {
        p: {"relative_rms": 0.01, "relative_weighted_l2": 0.01, "relative_max": 0.02}
        for p in ("mixed", "quadratic", "aggregate")
    }
    resolution = {
        "pieces": {
            p: {
                "coarse_to_medium_relative_rms": 0.02,
                "medium_to_fine_relative_rms": 0.01,
            }
            for p in ("mixed", "quadratic", "aggregate")
        }
    }
    axis = {
        p: {"normalized_error": 0.01}
        for p in ("mixed", "quadratic", "aggregate")
    }
    power = {
        p: {
            "production_abs_max": 0.0,
            "independent_abs_max": 0.0,
            "difference_abs_max": 0.0,
        }
        for p in ("mixed", "quadratic", "aggregate")
    }
    return {
        "fine_production_vs_independent": fine,
        "spatial_resolution": copy.deepcopy(resolution),
        "angular_resolution": copy.deepcopy(resolution),
        "independent_piece_closure_relative_rms": 0.0,
        "independent_A_plus_B_closure_relative_rms": 0.0,
        "axis_near": axis,
        "power_law_zero_support": power,
        "inner_production_quadratic_rms": 1e-3,
        "inner_independent_quadratic_rms": 1e-3,
        "inner_production_aggregate_rms": 1e-3,
        "inner_independent_aggregate_rms": 1e-3,
        "scientific_boundary": {
            "complete_ns_defect": False,
            "matched_pressure_materialized": False,
            "restricted_forcing_materialized": False,
            "correction_velocity_materialized": False,
            "after_correction_residual_assessed": False,
            "heldout_normalized_ns_residual_assessed": False,
            "same_protocol_comparable_to_st006": False,
            "pde_validated": False,
            "final_normalized_momentum_gate": 1e-3,
            "final_normalized_divergence_gate": 1e-5,
        },
    }


def test_frozen_protocol_is_deterministic_and_not_agent3_angular_ladder():
    first = frozen_inner_ring_specs()
    second = frozen_inner_ring_specs()
    assert first == second
    assert len(first) == 6
    assert FROZEN_SEED == 9173771
    assert FROZEN_SPATIAL_STEPS == (0.012, 0.006, 0.003)
    assert FROZEN_ANGULAR_ORDERS == (20, 40, 80)
    assert FROZEN_ANGULAR_ORDERS != (32, 64, 128)
    assert tuple(s.radius for s in frozen_axis_near_specs()) == (1e-4, 0.005, 0.02)


def test_fd4_plus_gauss_legendre_recovers_manufactured_rotating_mean():
    r = 0.47
    result = independent_ring_means(
        _Total(),
        _Leading(),
        RingSpec("manufactured", "inner", r, -0.08, 0.39),
        spatial_step=0.01,
        angular_order=80,
    )
    np.testing.assert_allclose(result["mixed"], np.zeros(3), atol=3e-13, rtol=0.0)
    np.testing.assert_allclose(result["quadratic"], np.asarray([-r, 0.0, 0.0]), atol=3e-13, rtol=0.0)
    np.testing.assert_allclose(result["aggregate"], np.asarray([-r, 0.0, 0.0]), atol=3e-13, rtol=0.0)
    np.testing.assert_allclose(result["A"] + result["B"], result["mixed"], atol=3e-13, rtol=0.0)


def test_public_receipt_identity_and_coordinates_are_fail_closed():
    spec = frozen_inner_ring_specs()[0]
    parsed = _validate_receipt(_production_receipt(spec), spec)
    assert set(parsed) == {"A", "B", "mixed", "quadratic", "aggregate"}
    bad = _production_receipt(spec)
    bad["radius"] = [spec.radius + 1e-3]
    with pytest.raises(ValueError, match="radius drifted"):
        _validate_receipt(bad, spec)


def test_public_receipt_backend_mutation_is_rejected():
    spec = frozen_inner_ring_specs()[0]
    bad = _production_receipt(spec)
    bad["agent2_backend"]["agent2_composite_head"] = "0" * 40
    with pytest.raises(ValueError, match="A2 identity drifted"):
        _validate_receipt(bad, spec)


def test_a3_truth_boundary_cannot_promote_pde_state():
    _validate_a3_truth_boundary(_a3_boundary())
    bad = _a3_boundary()
    bad["pde_validated"] = True
    with pytest.raises(ValueError, match="pde_validated"):
        _validate_a3_truth_boundary(bad)


def test_preregistered_gate_accepts_consistent_scoped_report():
    report = _passing_report()
    assert scoped_gates_pass(report)
    enforce_scoped_gates(report)


def test_mean_mutation_is_detected_without_threshold_retiming():
    report = _passing_report()
    report["fine_production_vs_independent"]["quadratic"]["relative_rms"] = 0.051
    assert not scoped_gates_pass(report)
    with pytest.raises(AssertionError):
        enforce_scoped_gates(report)


def test_resolution_degradation_and_piece_closure_mutations_are_detected():
    report = _passing_report()
    report["spatial_resolution"]["pieces"]["aggregate"] = {
        "coarse_to_medium_relative_rms": 0.021,
        "medium_to_fine_relative_rms": 0.049,
    }
    assert not scoped_gates_pass(report)

    report = _passing_report()
    report["independent_piece_closure_relative_rms"] = PIECE_CLOSURE_RELATIVE_GATE * 1.01
    assert not scoped_gates_pass(report)


def test_power_law_zero_support_mutation_is_detected():
    report = _passing_report()
    report["power_law_zero_support"]["aggregate"]["independent_abs_max"] = 1.01e-10
    assert not scoped_gates_pass(report)


def test_truth_boundary_exposes_no_caller_threshold_knob():
    boundary = truth_boundary()
    assert boundary["pde_validated"] is False
    assert boundary["complete_ns_defect"] is False
    assert boundary["final_normalized_momentum_gate"] == 1e-3
    assert boundary["final_normalized_divergence_gate"] == 1e-5
    params = inspect.signature(
        __import__(
            "openai_ns_reconstruction.kokuno_a4_current_rf40_power_law_nonlinear_mean_independent_audit",
            fromlist=["audit_current_rf40_power_law_nonlinear_mean"],
        ).audit_current_rf40_power_law_nonlinear_mean
    ).parameters
    assert "threshold" not in params
    assert "forcing" not in params
    assert "residual" not in params
