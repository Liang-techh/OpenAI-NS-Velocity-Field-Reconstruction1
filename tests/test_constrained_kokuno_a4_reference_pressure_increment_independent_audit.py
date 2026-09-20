from __future__ import annotations

import inspect
import json

import numpy as np
import pytest

from openai_ns_reconstruction import kokuno_a4_reference_pressure_increment_independent_audit as audit
from openai_ns_reconstruction.kokuno_pa10_reference_pressure_increment_x100 import (
    KokunoPA10ReferencePressureIncrementToX100,
)


class _ManufacturedPressureIncrement:
    semantic_sha256 = "manufactured-pressure-increment-v1"

    @staticmethod
    def pressure_increment(X, eta):
        X, eta = np.broadcast_arrays(np.asarray(X, dtype=float), np.asarray(eta, dtype=float))
        return X * (1.0 + 0.20 * eta + 0.05 * eta * eta) + 1.0e-3 * X * X

    @staticmethod
    def radial_derivative(X, eta):
        X, eta = np.broadcast_arrays(np.asarray(X, dtype=float), np.asarray(eta, dtype=float))
        return 1.0 + 0.20 * eta + 0.05 * eta * eta + 2.0e-3 * X

    @staticmethod
    def eta_derivative(X, eta):
        X, eta = np.broadcast_arrays(np.asarray(X, dtype=float), np.asarray(eta, dtype=float))
        return X * (0.20 + 0.10 * eta)


class _AxisShiftedPressureIncrement(_ManufacturedPressureIncrement):
    semantic_sha256 = "manufactured-axis-shifted-pressure-increment-v1"

    @staticmethod
    def pressure_increment(X, eta):
        return _ManufacturedPressureIncrement.pressure_increment(X, eta) + 1.0e-6


def test_public_scientific_entry_has_no_caller_protocol_or_threshold_knobs():
    signature = inspect.signature(audit.materialize_reference_pressure_increment_audit)
    assert list(signature.parameters) == ["pressure"]
    contract = audit.public_contract()
    assert contract["protocol"]["seed"] == 9173571
    assert contract["protocol"]["x_steps"] == (8.0e-3, 4.0e-3, 2.0e-3)
    assert contract["protocol"]["eta_steps"] == (8.0e-4, 4.0e-4, 2.0e-4)
    assert contract["frozen_final_gates"]["normalized_momentum_max"] == 1.0e-3
    assert contract["frozen_final_gates"]["normalized_momentum_volume_l2"] == 1.0e-3
    assert contract["frozen_final_gates"]["divergence_max"] == 1.0e-5
    assert contract["frozen_final_gates"]["divergence_volume_l2"] == 1.0e-5


def test_manufactured_blackbox_pressure_derivatives_pass_frozen_audit():
    receipt = audit.materialize_reference_pressure_increment_audit(_ManufacturedPressureIncrement())
    assert receipt["passed"] is True
    assert receipt["sample_counts"] == {
        "fresh_random_offgrid": 128,
        "axis_near": 3,
        "exact_axis": 5,
    }
    assert receipt["checks"]["axis_increment"] is True
    assert receipt["checks"]["radial_relative_rms"] is True
    assert receipt["checks"]["eta_relative_rms"] is True
    assert receipt["checks"]["radial_resolution_stability"] is True
    assert receipt["checks"]["eta_resolution_stability"] is True
    assert receipt["truth_boundary"]["reference_pressure_increment_independently_audited"] is True
    assert receipt["truth_boundary"]["absolute_reference_pressure_materialized"] is False
    assert receipt["truth_boundary"]["cartesian_pressure_gradient_assessed"] is False
    assert receipt["truth_boundary"]["heldout_complete_ns_residual_assessed"] is False
    assert receipt["truth_boundary"]["pde_validated"] is False


def test_mutations_are_detected_without_changing_original_gates():
    receipt = audit.materialize_reference_pressure_increment_audit(_ManufacturedPressureIncrement())
    assert receipt["mutation_checks"] == {
        "production_radial_derivative_times_0p99_detected": True,
        "production_eta_derivative_times_0p90_detected": True,
        "pressure_increment_plus_1e_minus3_X_detected": True,
    }
    assert receipt["checks"]["production_radial_derivative_times_0p99_detected"] is True
    assert receipt["checks"]["production_eta_derivative_times_0p90_detected"] is True
    assert receipt["checks"]["pressure_increment_plus_1e_minus3_X_detected"] is True


def test_exact_axis_increment_is_fail_closed():
    receipt = audit.materialize_reference_pressure_increment_audit(_AxisShiftedPressureIncrement())
    assert receipt["axis_increment_max_abs"] >= 1.0e-6
    assert receipt["checks"]["axis_increment"] is False
    assert receipt["passed"] is False
    assert receipt["truth_boundary"]["reference_pressure_increment_independently_audited"] is False
    assert receipt["truth_boundary"]["pde_validated"] is False


def test_agent1_configuration_mutation_is_rejected_after_serialization(tmp_path):
    pressure = KokunoPA10ReferencePressureIncrementToX100()
    path = tmp_path / "pressure.json"
    pressure.save_configuration(path)
    payload = json.loads(path.read_text())
    payload["eta_fd_step"] *= 1.001
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    with pytest.raises(ValueError, match="configuration is frozen"):
        KokunoPA10ReferencePressureIncrementToX100.load_configuration(path)


def test_real_agent1_serialized_pressure_increment_runs_frozen_scientific_audit(tmp_path):
    original = KokunoPA10ReferencePressureIncrementToX100()
    config_path = tmp_path / "reference-pressure-increment.json"
    original.save_configuration(config_path)
    rebound = KokunoPA10ReferencePressureIncrementToX100.load_configuration(config_path)
    assert rebound.semantic_sha256 == original.semantic_sha256

    receipt = audit.materialize_reference_pressure_increment_audit(rebound)
    audit.save_audit_receipt(tmp_path / "receipt.json", receipt)
    assert receipt["agent1_exact_head"] == "3f43931cec1fc762678a5a8134adfceaf7edfa81"
    assert receipt["pressure_semantic_sha256"] == original.semantic_sha256
    assert receipt["truth_boundary"]["absolute_reference_pressure_materialized"] is False
    assert receipt["truth_boundary"]["cartesian_pressure_gradient_assessed"] is False
    assert receipt["truth_boundary"]["pde_validated"] is False
    assert receipt["passed"] is True, json.dumps(
        {
            "radial_three_resolution_match": receipt["radial_three_resolution_match"],
            "eta_three_resolution_match": receipt["eta_three_resolution_match"],
            "radial_resolution_stability": receipt["radial_resolution_stability"],
            "eta_resolution_stability": receipt["eta_resolution_stability"],
            "worst_radial_match": receipt["worst_radial_match"],
            "worst_eta_match": receipt["worst_eta_match"],
            "checks": receipt["checks"],
        },
        sort_keys=True,
    )


def test_parent_identity_and_truth_boundary_are_exactly_scoped():
    contract = audit.public_contract()
    assert contract["agent1_pr"] == 907
    assert contract["agent1_exact_head"] == "3f43931cec1fc762678a5a8134adfceaf7edfa81"
    assert contract["agent1_source_blob"] == "4704fed4271da2c6070110a6935d2f33f7932f56"
    assert contract["scientific_input_surface"] == "pressure_increment(X,eta)"
    assert "Richardson" in contract["independent_operator"]
    truth = contract["truth_boundary"]
    assert truth["absolute_axis_pressure_Pi0_materialized"] is False
    assert truth["matched_global_pressure_materialized"] is False
    assert truth["leading_only_ns_residual_assessed"] is False
    assert truth["leading_plus_oscillatory_ns_residual_assessed"] is False
    assert truth["after_correction_ns_residual_assessed"] is False
    assert truth["pde_validated"] is False
