from __future__ import annotations

import inspect
import json
import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_a4_reference_stress_xi110_independent_audit import (
    A1_REFERENCE_STRESS_HEAD,
    A1_REFERENCE_STRESS_PR,
    FROZEN_FINAL_GATES,
    LOCAL_POLY_NODES,
    PROTOCOL,
    TRUTH_BOUNDARY,
    local_polynomial_first,
    local_polynomial_first_weights,
    public_contract,
    run_real_artifact_audit,
)
from openai_ns_reconstruction.kokuno_pa10_autonomous_reference_stress_xi110 import (
    KokunoPA10AutonomousReferenceStressToXi110,
)


def test_contract_freezes_upstream_identity_and_final_gates() -> None:
    contract = public_contract()
    assert contract["agent1_pr"] == A1_REFERENCE_STRESS_PR == 933
    assert (
        contract["agent1_exact_head"]
        == A1_REFERENCE_STRESS_HEAD
        == "15791891fbbcfdaa614930791a843983afde2b4d"
    )
    assert FROZEN_FINAL_GATES == {
        "normalized_momentum_max": 1.0e-3,
        "normalized_momentum_volume_l2": 1.0e-3,
        "divergence_max": 1.0e-5,
        "divergence_volume_l2": 1.0e-5,
    }
    assert contract["protocol"]["seed"] == 9173611
    assert tuple(contract["protocol"]["x_halfwidths"]) == (0.24, 0.12, 0.06)
    assert contract["truth_boundary"]["pde_validated"] is False


def test_scientific_entry_point_exposes_no_threshold_resolution_or_seed_knob() -> None:
    from openai_ns_reconstruction import (
        kokuno_a4_reference_stress_xi110_independent_audit as module,
    )

    params = inspect.signature(module.materialize_reference_stress_xi_audit).parameters
    assert tuple(params) == ("field",)
    forbidden = {
        "seed",
        "resolution",
        "halfwidth",
        "threshold",
        "momentum_threshold",
        "divergence_threshold",
        "forcing",
        "pressure",
        "residual",
    }
    assert forbidden.isdisjoint(params)


def test_nonuniform_local_polynomial_derivative_is_degree_six_exact() -> None:
    weights = local_polynomial_first_weights()
    for degree in range(7):
        observed = float(np.sum(weights * LOCAL_POLY_NODES**degree))
        expected = 1.0 if degree == 1 else 0.0
        assert observed == pytest.approx(expected, abs=2.0e-12)

    x = np.asarray([0.27, 0.43, 0.71])
    eta = np.asarray([-0.2, 0.0, 0.31])

    def poly(xx: np.ndarray, ee: np.ndarray) -> np.ndarray:
        return (
            0.3
            + 0.7 * xx
            - 0.4 * xx**2
            + 0.2 * xx**3
            - 0.03 * xx**4
            + 0.01 * xx**5
            - 0.002 * xx**6
            + 0.05 * ee
        )

    exact = (
        0.7
        - 0.8 * x
        + 0.6 * x**2
        - 0.12 * x**3
        + 0.05 * x**4
        - 0.012 * x**5
    )
    observed = local_polynomial_first(poly, x, eta, 0.09)
    assert np.max(np.abs(observed - exact)) <= 3.0e-12


def test_upstream_loader_rejects_checksum_valid_scientific_configuration_drift(
    tmp_path,
) -> None:
    obj = KokunoPA10AutonomousReferenceStressToXi110()
    path = tmp_path / "candidate.json"
    payload = obj.save_configuration(path)
    mutated = dict(payload)
    mutated["stress_quadrature_order"] = int(mutated["stress_quadrature_order"]) + 1
    path.write_text(json.dumps(mutated, indent=2, sort_keys=True) + "\n")
    with pytest.raises(ValueError, match="configuration is frozen"):
        KokunoPA10AutonomousReferenceStressToXi110.load_configuration(path)


def test_real_save_reload_reference_stress_xi_audit_passes_scoped_gates() -> None:
    receipt = run_real_artifact_audit()

    assert receipt["save_reload_exact_replay"] is True
    assert receipt["passed"] is True, json.dumps(receipt, indent=2, sort_keys=True)
    assert receipt["sample_count"] == PROTOCOL.random_offgrid_count + 10
    assert all(receipt["gates"].values())
    assert all(receipt["mutation_checks"].values())

    for key in (
        "E_reference_X",
        "Pi_reference_autonomous_X",
        "p1_reference_X",
        "N_s_reference_autonomous_X",
        "n_s_reference_autonomous_X",
    ):
        fine = receipt["resolution_metrics"][key][-1]
        assert fine["relative_rms"] <= PROTOCOL.derivative_relative_rms_gate
        assert fine["relative_sampled_max"] <= PROTOCOL.derivative_relative_max_gate
        assert receipt["resolution_stability"][key]["passed"] is True
        assert math.isfinite(receipt["worst_witness"][key]["X"])
        assert math.isfinite(receipt["worst_witness"][key]["eta"])

    assert (
        receipt["algebraic_closure_relative"]
        <= PROTOCOL.algebraic_closure_relative_gate
    )
    assert (
        receipt["independent_derivative_closure_relative"]
        <= PROTOCOL.derivative_closure_relative_gate
    )
    assert receipt["plateau_F_abs_max"] <= PROTOCOL.plateau_abs_gate
    assert receipt["plateau_U_abs_max"] <= PROTOCOL.plateau_abs_gate
    assert receipt["Xi_public_handoff_abs_max"] <= PROTOCOL.xi_handoff_abs_gate

    truth = receipt["truth_boundary"]
    assert truth["autonomous_reference_stress_to_Xi_independently_audited"] is True
    assert truth["actual_final_interpolation_100_to_Xi_materialized"] is False
    assert truth["outer_global_leading_velocity_materialized"] is False
    assert truth["matched_global_pressure_materialized"] is False
    assert truth["restricted_forcing_materialized"] is False
    assert truth["correction_velocity_materialized"] is False
    assert truth["leading_only_ns_residual_assessed"] is False
    assert truth["leading_plus_oscillatory_ns_residual_assessed"] is False
    assert truth["after_correction_ns_residual_assessed"] is False
    assert truth["pde_validated"] is False


def test_module_truth_boundary_never_preclaims_full_ns() -> None:
    assert TRUTH_BOUNDARY["pde_validated"] is False
    assert TRUTH_BOUNDARY["heldout_complete_ns_residual_assessed"] is False
    assert TRUTH_BOUNDARY["same_protocol_comparable_to_st006"] is False
