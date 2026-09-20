from __future__ import annotations

import json

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_a4_autonomous_reference_ns_independent_audit import (
    A1_AUTONOMOUS_REFERENCE_NS_HEAD,
    A1_AUTONOMOUS_REFERENCE_NS_PR,
    A1_AUTONOMOUS_REFERENCE_NS_SOURCE_BLOB,
    FROZEN_FINAL_GATES,
    PROTOCOL,
    materialize_autonomous_reference_stress_audit,
    public_contract,
    save_real_artifact_audit,
)
from openai_ns_reconstruction.kokuno_pa10_autonomous_reference_ns_x100 import (
    KokunoPA10AutonomousReferenceNsToX100,
)


class _ManufacturedStress:
    X_0 = 0.2
    X_b_ref = 100.0
    h = 0.1
    semantic_sha256 = "a" * 64

    @staticmethod
    def _L(eta):
        eta = np.asarray(eta, dtype=float)
        return 1.0 - 0.2 * eta * eta

    @staticmethod
    def _N(X, eta):
        X = np.asarray(X, dtype=float)
        eta = np.asarray(eta, dtype=float)
        return (1.0 + 0.15 * eta * eta) * (1.0 + 0.2 * X + 0.03 * X * X)

    @staticmethod
    def _N_X(X, eta):
        X = np.asarray(X, dtype=float)
        eta = np.asarray(eta, dtype=float)
        return (1.0 + 0.15 * eta * eta) * (0.2 + 0.06 * X)

    def values(self, X, eta):
        N = self._N(X, eta)
        n = N / self._L(eta)
        return {
            "N_s_reference_autonomous": N,
            "n_s_reference_autonomous": n,
            "p2_reference_autonomous": 0.5 + np.asarray(X, dtype=float) * n,
        }

    def radial_derivatives(self, X, eta):
        N_X = self._N_X(X, eta)
        return {
            "N_s_reference_autonomous_X": N_X,
            "n_s_reference_autonomous_X": N_X / self._L(eta),
        }

    def initial_values(self, eta):
        eta = np.asarray(eta, dtype=float)
        X = np.full_like(eta, self.X_0)
        N = self._N(X, eta)
        return {
            "N_s_initial": N,
            "n_s_initial": N / self._L(eta),
        }


class _MutatedProductionStress(_ManufacturedStress):
    def radial_derivatives(self, X, eta):
        out = super().radial_derivatives(X, eta)
        return {
            "N_s_reference_autonomous_X": 0.97 * out["N_s_reference_autonomous_X"],
            "n_s_reference_autonomous_X": 1.04 * out["n_s_reference_autonomous_X"],
        }


def test_public_contract_freezes_upstream_identity_and_final_gates() -> None:
    contract = public_contract()
    assert contract["agent1_pr"] == A1_AUTONOMOUS_REFERENCE_NS_PR == 918
    assert contract["agent1_exact_head"] == A1_AUTONOMOUS_REFERENCE_NS_HEAD
    assert contract["agent1_source_blob"] == A1_AUTONOMOUS_REFERENCE_NS_SOURCE_BLOB
    assert PROTOCOL.seed == 9173591
    assert PROTOCOL.relative_x_steps == (2.0e-3, 1.0e-3, 5.0e-4)
    assert FROZEN_FINAL_GATES["normalized_momentum_max"] == 1.0e-3
    assert FROZEN_FINAL_GATES["normalized_momentum_volume_l2"] == 1.0e-3
    assert FROZEN_FINAL_GATES["divergence_max"] == 1.0e-5
    assert FROZEN_FINAL_GATES["divergence_volume_l2"] == 1.0e-5
    assert contract["truth_boundary"]["pde_validated"] is False


def test_manufactured_public_values_pass_implementation_distinct_fd6_audit() -> None:
    receipt = materialize_autonomous_reference_stress_audit(_ManufacturedStress())
    assert receipt["passed"] is True
    assert all(receipt["checks"].values())
    assert receipt["truth_boundary"]["autonomous_reference_stress_independently_audited"] is True
    assert receipt["truth_boundary"]["pde_validated"] is False
    assert receipt["N_three_resolution_match"][-1]["relative_rms"] < 1.0e-8
    assert receipt["n_three_resolution_match"][-1]["relative_rms"] < 1.0e-8


def test_production_derivative_mutation_is_not_laundered() -> None:
    receipt = materialize_autonomous_reference_stress_audit(_MutatedProductionStress())
    assert receipt["passed"] is False
    assert (
        receipt["checks"]["N_derivative_relative_rms"] is False
        or receipt["checks"]["N_derivative_relative_sampled_max"] is False
    )
    assert (
        receipt["checks"]["n_derivative_relative_rms"] is False
        or receipt["checks"]["n_derivative_relative_sampled_max"] is False
    )
    assert receipt["truth_boundary"]["pde_validated"] is False


def test_real_upstream_configuration_drift_is_rejected(tmp_path) -> None:
    obj = KokunoPA10AutonomousReferenceNsToX100()
    path = tmp_path / "stress.json"
    cfg = obj.save_configuration(path)
    bad = json.loads(json.dumps(cfg))
    bad["Ns_quadrature_order"] = int(bad["Ns_quadrature_order"]) + 1
    path.write_text(json.dumps(bad))
    with pytest.raises(ValueError):
        KokunoPA10AutonomousReferenceNsToX100.load_configuration(path)


def test_real_save_reload_artifact_runs_frozen_audit(tmp_path) -> None:
    path = tmp_path / "a4_reference_stress_receipt.json"
    receipt = save_real_artifact_audit(path)
    assert path.exists()
    assert receipt["save_reload_semantic_identity_exact"] is True
    assert receipt["stress_semantic_sha256"] == KokunoPA10AutonomousReferenceNsToX100().semantic_sha256
    if not receipt["passed"]:
        pytest.fail(
            "real autonomous-reference-stress audit failed: "
            + json.dumps(
                {
                    "checks": receipt["checks"],
                    "N_three_resolution_match": receipt["N_three_resolution_match"],
                    "n_three_resolution_match": receipt["n_three_resolution_match"],
                    "N_resolution_stability": receipt["N_resolution_stability"],
                    "n_resolution_stability": receipt["n_resolution_stability"],
                    "worst_N_derivative_match": receipt["worst_N_derivative_match"],
                    "worst_n_derivative_match": receipt["worst_n_derivative_match"],
                },
                sort_keys=True,
            )
        )
    assert receipt["truth_boundary"]["pde_validated"] is False
    assert receipt["truth_boundary"]["heldout_complete_ns_residual_assessed"] is False
    assert receipt["truth_boundary"]["same_protocol_comparable_to_st006"] is False
