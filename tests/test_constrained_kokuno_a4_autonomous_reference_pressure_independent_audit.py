from __future__ import annotations

import inspect
import json

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_a4_autonomous_reference_pressure_independent_audit import (
    A1_AUTONOMOUS_PRESSURE_HEAD,
    A1_AUTONOMOUS_PRESSURE_SOURCE_BLOB,
    FROZEN_FINAL_GATES,
    PROTOCOL,
    materialize_autonomous_reference_pressure_audit,
    public_contract,
    save_real_artifact_audit,
)
from openai_ns_reconstruction.kokuno_pa10_autonomous_axis_pressure_seed_x100 import (
    KokunoPA10AutonomousAxisPressureSeedToX100,
)


class _PolynomialPressure:
    """Analytic mechanics fixture; not Kokuno candidate evidence."""

    semantic_sha256 = "a" * 64

    @staticmethod
    def axis_pressure(eta):
        eta = np.asarray(eta, dtype=float)
        return 1.0 + 0.2 * eta + 0.05 * eta * eta

    @staticmethod
    def axis_pressure_eta(eta):
        eta = np.asarray(eta, dtype=float)
        return 0.2 + 0.1 * eta

    @classmethod
    def pressure(cls, X, eta):
        X = np.asarray(X, dtype=float)
        eta = np.asarray(eta, dtype=float)
        return cls.axis_pressure(eta) + 0.4 * X + 0.03 * X * X

    @staticmethod
    def radial_derivative(X, eta):
        del eta
        X = np.asarray(X, dtype=float)
        return 0.4 + 0.06 * X

    @classmethod
    def eta_derivative(cls, X, eta):
        del X
        return cls.axis_pressure_eta(eta)


def test_public_contract_freezes_upstream_identity_and_final_gates() -> None:
    contract = public_contract()
    assert contract["agent1_pr"] == 913
    assert contract["agent1_exact_head"] == A1_AUTONOMOUS_PRESSURE_HEAD
    assert A1_AUTONOMOUS_PRESSURE_HEAD == "a5903ce594bf1c1f775d6fbba8a5eeca226622fa"
    assert contract["agent1_source_blob"] == A1_AUTONOMOUS_PRESSURE_SOURCE_BLOB
    assert A1_AUTONOMOUS_PRESSURE_SOURCE_BLOB == "68dddc4f6ad1fb5bb6955498099b380240be6750"
    assert FROZEN_FINAL_GATES["normalized_momentum_max"] == 1.0e-3
    assert FROZEN_FINAL_GATES["normalized_momentum_volume_l2"] == 1.0e-3
    assert FROZEN_FINAL_GATES["divergence_max"] == 1.0e-5
    assert FROZEN_FINAL_GATES["divergence_volume_l2"] == 1.0e-5

    truth = contract["truth_boundary"]
    assert truth["source_prepared_appendixA_Pi0_materialized"] is False
    assert truth["matched_global_pressure_materialized"] is False
    assert truth["cartesian_pressure_gradient_assessed"] is False
    assert truth["heldout_complete_ns_residual_assessed"] is False
    assert truth["pde_validated"] is False


def test_scientific_entry_point_exposes_no_protocol_or_threshold_knobs() -> None:
    params = inspect.signature(materialize_autonomous_reference_pressure_audit).parameters
    assert list(params) == ["pressure"]
    assert PROTOCOL.seed == 9173581
    assert PROTOCOL.x_steps == (4.0e-3, 2.0e-3, 1.0e-3)
    assert PROTOCOL.eta_steps == (4.0e-4, 2.0e-4, 1.0e-4)
    assert PROTOCOL.refinement_ratio_gate == 20.0


def test_polynomial_mechanics_fixture_passes_and_mutations_are_detected() -> None:
    receipt = materialize_autonomous_reference_pressure_audit(_PolynomialPressure())
    if not receipt["passed"]:
        pytest.fail(json.dumps(receipt, indent=2, sort_keys=True))

    assert receipt["checks"]["radial_relative_rms"]
    assert receipt["checks"]["eta_relative_rms"]
    assert receipt["checks"]["axis_eta_relative_rms"]
    assert receipt["checks"]["exact_axis_composition"]
    assert receipt["checks"]["radial_resolution_stability"]
    assert receipt["checks"]["eta_resolution_stability"]
    assert receipt["checks"]["axis_eta_resolution_stability"]
    assert all(receipt["mutation_checks"].values())
    assert receipt["truth_boundary"]["autonomous_reference_pressure_independently_audited"]
    assert receipt["truth_boundary"]["pde_validated"] is False
    assert receipt["truth_boundary"]["cartesian_pressure_gradient_assessed"] is False


def test_real_agent1_artifact_save_reload_and_independent_audit(tmp_path) -> None:
    original = KokunoPA10AutonomousAxisPressureSeedToX100()
    config_path = tmp_path / "a1_pressure.json"
    original.save_configuration(config_path)
    rebound = KokunoPA10AutonomousAxisPressureSeedToX100.load_configuration(config_path)
    assert rebound.semantic_sha256 == original.semantic_sha256

    receipt = materialize_autonomous_reference_pressure_audit(rebound)
    if not receipt["passed"]:
        pytest.fail(json.dumps(receipt, indent=2, sort_keys=True))

    assert receipt["pressure_semantic_sha256"] == original.semantic_sha256
    assert receipt["sample_counts"] == {
        "fresh_random_offgrid": 128,
        "axis_near": 3,
        "exact_axis": 5,
    }
    assert receipt["exact_axis_composition_max_abs"] <= PROTOCOL.axis_composition_abs_gate
    assert all(receipt["mutation_checks"].values())
    assert receipt["truth_boundary"]["autonomous_reference_pressure_independently_audited"]
    assert receipt["truth_boundary"]["source_prepared_appendixA_Pi0_materialized"] is False
    assert receipt["truth_boundary"]["matched_global_pressure_materialized"] is False
    assert receipt["truth_boundary"]["pde_validated"] is False


def test_real_receipt_writer_binds_save_reload_identity(tmp_path) -> None:
    path = tmp_path / "receipt.json"
    payload = save_real_artifact_audit(path)
    loaded = json.loads(path.read_text())
    assert loaded == payload
    assert payload["save_reload_semantic_identity_exact"] is True
    assert len(payload["receipt_sha256"]) == 64
    if not payload["passed"]:
        pytest.fail(json.dumps(payload, indent=2, sort_keys=True))


def test_missing_public_scalar_or_derivative_surface_fails_closed() -> None:
    class MissingPressure:
        semantic_sha256 = "b" * 64

        @staticmethod
        def radial_derivative(X, eta):
            return np.asarray(X) * 0.0 + np.asarray(eta) * 0.0

        @staticmethod
        def eta_derivative(X, eta):
            return np.asarray(X) * 0.0 + np.asarray(eta) * 0.0

        @staticmethod
        def axis_pressure(eta):
            return np.asarray(eta) * 0.0

        @staticmethod
        def axis_pressure_eta(eta):
            return np.asarray(eta) * 0.0

    with pytest.raises(TypeError, match="pressure"):
        materialize_autonomous_reference_pressure_audit(MissingPressure())


def test_scope_never_promotes_autonomous_seed_to_matched_pressure_or_ns() -> None:
    receipt = materialize_autonomous_reference_pressure_audit(_PolynomialPressure())
    truth = receipt["truth_boundary"]
    assert truth["repository_autonomous_axis_pressure_seed_materialized"] is True
    assert truth["source_prepared_appendixA_Pi0_materialized"] is False
    assert truth["source_exact_absolute_reference_pressure_materialized"] is False
    assert truth["matched_global_pressure_materialized"] is False
    assert truth["cartesian_pressure_gradient_assessed"] is False
    assert truth["leading_only_ns_residual_assessed"] is False
    assert truth["leading_plus_oscillatory_ns_residual_assessed"] is False
    assert truth["after_correction_ns_residual_assessed"] is False
    assert truth["heldout_complete_ns_residual_assessed"] is False
    assert truth["same_protocol_comparable_to_st006"] is False
    assert truth["pde_validated"] is False
