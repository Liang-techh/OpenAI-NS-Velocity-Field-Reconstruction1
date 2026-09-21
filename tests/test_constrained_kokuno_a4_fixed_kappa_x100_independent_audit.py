from __future__ import annotations

import copy
from pathlib import Path

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_a4_fixed_kappa_x100_independent_audit import (
    A1_FIXED_KAPPA_HEAD,
    A1_FIXED_KAPPA_PR,
    FROZEN_FINAL_GATES,
    PROTOCOL,
    TRUTH_BOUNDARY,
    centered_fd8_first,
    materialize_fixed_kappa_x100_audit,
    materialize_real_fixed_kappa_x100_audit,
    public_contract,
)
from openai_ns_reconstruction.kokuno_pa10_fixed_kappa_x100 import (
    KokunoPA10FixedKappaContinuationToX100,
)


class _AnalyticField:
    X_1 = 2.0
    X_b_ref = 10.0
    eta_interval = (-0.8, 0.8)
    kappa_0 = 0.1
    semantic_sha256 = "analytic-fixed-kappa-fixture"

    def values(self, X, eta):
        X, eta = np.broadcast_arrays(np.asarray(X, dtype=float), np.asarray(eta, dtype=float))
        F = np.exp(-0.03 * X) * (1.0 + 0.02 * eta * eta)
        U = 0.4 + 0.01 * X + 0.03 * eta
        E = np.sqrt(2.0 * X) * F
        return {
            "F_fixed_kappa": F,
            "U_fixed_kappa": U,
            "E_fixed_kappa": E,
            "delta_log_F_from_activation_endpoint": np.log(F / self.values_at_x1(eta)[0]),
            "delta_U_from_activation_endpoint": U - self.values_at_x1(eta)[1],
        }

    def values_at_x1(self, eta):
        eta = np.asarray(eta, dtype=float)
        F = np.exp(-0.03 * self.X_1) * (1.0 + 0.02 * eta * eta)
        U = 0.4 + 0.01 * self.X_1 + 0.03 * eta
        return F, U

    def radial_derivatives(self, X, eta):
        X, eta = np.broadcast_arrays(np.asarray(X, dtype=float), np.asarray(eta, dtype=float))
        values = self.values(X, eta)
        F = values["F_fixed_kappa"]
        F_X = -0.03 * F
        U_X = np.full_like(X, 0.01)
        E_X = F / np.sqrt(2.0 * X) + np.sqrt(2.0 * X) * F_X
        return {
            "F_fixed_kappa_X": F_X,
            "U_fixed_kappa_X": U_X,
            "E_fixed_kappa_X": E_X,
        }

    def endpoint_values(self, eta):
        F, U = self.values_at_x1(eta)
        return {"F_activation": F, "U_activation": U}

    def save_configuration(self, path):
        Path(path).write_text("{}\n")
        return {}


def test_frozen_protocol_and_final_scientific_gates_are_not_relaxed() -> None:
    assert A1_FIXED_KAPPA_PR == 925
    assert A1_FIXED_KAPPA_HEAD == "1f3dec711f8ce2052276cd951b2f44d2849579bf"
    assert PROTOCOL.seed == 9173601
    assert PROTOCOL.random_offgrid_count == 128
    assert PROTOCOL.relative_x_steps == (1.6e-3, 8.0e-4, 4.0e-4)
    assert FROZEN_FINAL_GATES == {
        "normalized_momentum_max": 1e-3,
        "normalized_momentum_volume_l2": 1e-3,
        "divergence_max": 1e-5,
        "divergence_volume_l2": 1e-5,
    }
    assert TRUTH_BOUNDARY["pde_validated"] is False
    assert TRUTH_BOUNDARY["leading_only_ns_residual_assessed"] is False
    assert TRUTH_BOUNDARY["after_correction_ns_residual_assessed"] is False


def test_public_contract_is_black_box_profile_only() -> None:
    contract = public_contract()
    assert contract["scientific_reference_surface"] == "save/reloaded values(X,eta) -> F,U,E"
    assert contract["production_surfaces_under_audit"] == [
        "radial_derivatives(X,eta) -> F_X,U_X,E_X"
    ]
    text = str(contract).lower()
    for forbidden in ("training loss", "free forcing", "residual-defined forcing"):
        assert forbidden not in text
    assert contract["truth_boundary"]["outer_global_leading_velocity_materialized"] is False
    assert contract["truth_boundary"]["matched_global_pressure_materialized"] is False


def test_fd8_sign_and_order_on_degree_seven_polynomial() -> None:
    x = np.asarray([1.3, 2.1, 3.7])
    eta = np.zeros_like(x)

    def fn(xx, ee):
        del ee
        return xx**7 - 2.0 * xx**4 + 0.3 * xx**2 - 5.0

    observed = centered_fd8_first(fn, x, eta, 2.0e-4)
    expected = 7.0 * x**6 - 8.0 * x**3 + 0.6 * x
    np.testing.assert_allclose(observed, expected, rtol=3.0e-10, atol=2.0e-9)


def test_analytic_fixture_passes_three_resolution_audit_and_mutations() -> None:
    report = materialize_fixed_kappa_x100_audit(_AnalyticField())
    assert report["passed"] is True, report
    assert report["gate_map"]["fine_derivative_match"] is True
    assert report["gate_map"]["three_resolution_stability"] is True
    assert report["gate_map"]["public_profile_identities"] is True
    assert report["gate_map"]["mutation_firewall"] is True
    assert all(report["mutation_checks"].values())
    assert report["truth_boundary"][
        "candidate_side_fixed_kappa_continuation_independently_audited"
    ] is True
    assert report["truth_boundary"]["pde_validated"] is False


def test_upstream_configuration_drift_is_rejected(tmp_path) -> None:
    field = KokunoPA10FixedKappaContinuationToX100()
    saved = field.configuration()

    kappa_mutation = copy.deepcopy(saved)
    kappa_mutation["selected_kappa0"] = float(saved["selected_kappa0"]) * 1.001
    with pytest.raises(ValueError, match="configuration is frozen"):
        KokunoPA10FixedKappaContinuationToX100.from_configuration(kappa_mutation)

    quadrature_mutation = copy.deepcopy(saved)
    quadrature_mutation["continuation_quadrature_order"] = int(
        saved["continuation_quadrature_order"]
    ) + 1
    with pytest.raises(ValueError, match="configuration is frozen"):
        KokunoPA10FixedKappaContinuationToX100.from_configuration(quadrature_mutation)

    path = tmp_path / "fixed-kappa.json"
    field.save_configuration(path)
    rebound = KokunoPA10FixedKappaContinuationToX100.load_configuration(path)
    assert rebound.semantic_sha256 == field.semantic_sha256


def test_real_saved_reloaded_fixed_kappa_artifact_passes_or_localizes_failure() -> None:
    report = materialize_real_fixed_kappa_x100_audit()
    assert report["save_reload_exact_replay"] is True
    assert report["passed"] is True, {
        "gate_map": report["gate_map"],
        "resolution_metrics": report["resolution_metrics"],
        "stability": report["stability"],
        "fine_worst_witness": report["fine_worst_witness"],
        "profile_identities": report["profile_identities"],
        "nontriviality": report["nontriviality"],
        "mutation_checks": report["mutation_checks"],
    }
    assert report["truth_boundary"]["pde_validated"] is False
    assert report["truth_boundary"]["same_protocol_comparable_to_st006"] is False
