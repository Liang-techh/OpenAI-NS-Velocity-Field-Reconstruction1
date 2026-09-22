from __future__ import annotations

import copy
import inspect

import numpy as np
import pytest

import openai_ns_reconstruction.kokuno_a4_current_i4_nonlinear_mean_independent_audit as audit


class ZeroField:
    def velocity(self, x, y, z, t):
        x, y, z, t = np.broadcast_arrays(
            np.asarray(x, float), np.asarray(y, float), np.asarray(z, float), np.asarray(t, float)
        )
        return np.zeros(x.shape + (3,), dtype=float)


class SolidRotation:
    def velocity(self, x, y, z, t):
        x, y, z, t = np.broadcast_arrays(
            np.asarray(x, float), np.asarray(y, float), np.asarray(z, float), np.asarray(t, float)
        )
        return np.stack((-y, x, np.zeros_like(z)), axis=-1)


class I4Geometry:
    X_I4_start = 10.0
    X_I4_end = 100.0
    D = 0.3


def _truth() -> dict[str, object]:
    out: dict[str, object] = {
        "forbidden_scientific_controls_exposed": False,
        "final_normalized_momentum_gate": 1.0e-3,
        "final_normalized_divergence_gate": 1.0e-5,
    }
    for key in (
        "current_I4_leading_plus_oscillatory_identity_consumed",
        "current_I4_mixed_nonlinear_mean_materialized",
        "current_I4_quadratic_oscillatory_mean_materialized",
        "current_I4_aggregate_nonlinear_mean_materialized",
        "source_I4_reserved_mean_correction_interval_recorded",
    ):
        out[key] = True
    for key in (
        "source_I3_positive_order_correction_materialized",
        "source_I4_five_row_mean_correction_materialized",
        "current_I4_correction_velocity_materialized",
        "post_I4_or_pulse_leading_identity_consumed",
        "outer_global_velocity_consumed",
        "radial_inverse_performed_in_this_increment",
        "compact_radial_stress_materialized_in_this_increment",
        "radial_force_materialized_in_this_increment",
        "scoped_nonlinear_mean_authorized_as_complete_ns_defect",
        "scoped_nonlinear_mean_authorized_as_correction_target",
        "matched_global_pressure_materialized",
        "restricted_forcing_materialized",
        "complete_ns_defect",
        "current_real_ns_correction_velocity_materialized",
        "real_candidate_finite_correction_cycle_run",
        "heldout_normalized_ns_residual_assessed",
        "residual_reduction_claimed",
        "same_protocol_comparable_to_st006",
        "pde_validated",
    ):
        out[key] = False
    return out


def _receipt(spec: audit.RingSpec) -> dict[str, object]:
    z = [0.0, 0.0, 0.0]
    return {
        "backend_kind": audit.PARENT_BACKEND_KIND,
        "radius": [spec.radius],
        "z": [spec.z],
        "t": [spec.t],
        "agent2_backend": {
            "agent2_composite_pr": audit.AGENT2_PR,
            "agent2_composite_head": audit.AGENT2_HEAD,
            "agent1_leading_pr": audit.AGENT1_PR,
            "agent1_leading_head": audit.AGENT1_HEAD,
        },
        "pointwise_three_piece_closure_absolute_max": 0.0,
        "projected_three_piece_closure_absolute_max": 0.0,
        "mean_inner_advects_oscillation_cylindrical": z,
        "mean_oscillation_advects_inner_cylindrical": z,
        "mean_mixed_cross_cylindrical": z,
        "mean_oscillatory_self_advection_cylindrical": z,
        "mean_aggregate_nonlinear_cylindrical": z,
    }


def test_manufactured_solid_rotation_recovers_quadratic_radial_mean():
    spec = audit.RingSpec("manufactured", "inner", 0.43, -0.07, 0.39)
    result = audit.independent_ring_means(
        SolidRotation(), ZeroField(), spec, spatial_step=0.003, angular_order=80
    )
    assert result["quadratic"][0] == pytest.approx(-spec.radius, abs=3.0e-12)
    assert result["quadratic"][1] == pytest.approx(0.0, abs=3.0e-12)
    assert result["quadratic"][2] == pytest.approx(0.0, abs=3.0e-12)
    assert np.max(np.abs(result["A"])) < 3.0e-12
    assert np.max(np.abs(result["B"])) < 3.0e-12


def test_frozen_probe_sets_are_deterministic_and_i4_strict():
    a = audit.frozen_inner_ring_specs()
    b = audit.frozen_inner_ring_specs()
    assert a == b
    assert len(a) == audit.FROZEN_INNER_RING_COUNT
    assert tuple(s.radius for s in audit.frozen_axis_near_specs()) == audit.FROZEN_AXIS_RADII
    i4 = audit.frozen_i4_specs(I4Geometry())
    assert len(i4) == 3
    for spec in i4:
        assert spec.radius > 0.0
        assert np.isfinite(spec.radius)
        assert np.isfinite(spec.z)


def test_receipt_provenance_mutation_fails_closed():
    spec = audit.RingSpec("r", "inner", 0.41, -0.08, 0.39)
    good = _receipt(spec)
    parsed = audit._validate_receipt(good, spec)
    assert set(parsed) == set(audit.PIECES)
    bad = copy.deepcopy(good)
    bad["agent2_backend"]["agent2_composite_head"] = "0" * 40
    with pytest.raises(ValueError, match="identity drifted"):
        audit._validate_receipt(bad, spec)


def test_truth_promotion_and_gate_drift_fail_closed():
    audit._validate_a3_truth_boundary(_truth())
    promoted = _truth()
    promoted["pde_validated"] = True
    with pytest.raises(ValueError, match="pde_validated"):
        audit._validate_a3_truth_boundary(promoted)
    drifted = _truth()
    drifted["final_normalized_momentum_gate"] = 2.0e-3
    with pytest.raises(ValueError, match="momentum gate"):
        audit._validate_a3_truth_boundary(drifted)


def test_low_signal_rule_uses_absolute_gate_not_vacuous_relative_ratio():
    good = {
        "signal_scale": 1.0e-14,
        "absolute_max_error": 0.5 * audit.LOW_SIGNAL_ABSOLUTE_GATE,
        "relative_rms": 1.0,
        "relative_weighted_l2": 1.0,
        "relative_max": 1.0,
    }
    assert audit._comparison_pass(good)
    bad = dict(good)
    bad["absolute_max_error"] = 1.01 * audit.LOW_SIGNAL_ABSOLUTE_GATE
    assert not audit._comparison_pass(bad)


def test_resolution_regression_is_rejected_above_floor():
    good = {
        "fine_signal_scale": 1.0,
        "coarse_to_medium_relative_rms": 0.03,
        "medium_to_fine_relative_rms": 0.02,
        "coarse_to_medium_absolute_max": 0.03,
        "medium_to_fine_absolute_max": 0.02,
    }
    assert audit._resolution_piece_pass(good)
    bad = dict(good)
    bad["coarse_to_medium_relative_rms"] = 0.02
    bad["medium_to_fine_relative_rms"] = 0.04
    assert not audit._resolution_piece_pass(bad)


def test_public_contract_exposes_no_scientific_tuning_knob():
    contract = audit.public_contract()
    assert contract["forbidden_scientific_controls_exposed"] is False
    assert contract["pde_validated"] is False
    assert tuple(inspect.signature(audit.audit_current_i4_nonlinear_mean).parameters) == (
        "field", "production_receipts", "a3_truth_boundary"
    )
