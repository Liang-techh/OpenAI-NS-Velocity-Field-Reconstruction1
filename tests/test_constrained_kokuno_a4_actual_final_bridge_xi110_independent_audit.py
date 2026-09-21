from __future__ import annotations

import inspect
import json
import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_a4_actual_final_bridge_xi110_independent_audit import (
    A1_FINAL_BRIDGE_HEAD,
    A1_FINAL_BRIDGE_PR,
    BACKWARD_ENDPOINT_NODES,
    FORWARD_ENDPOINT_NODES,
    FROZEN_FINAL_GATES,
    INTERIOR_NODES,
    PROTOCOL,
    TRUTH_BOUNDARY,
    first_derivative_weights,
    polynomial_first,
    public_contract,
    run_real_artifact_audit,
)
from openai_ns_reconstruction.kokuno_pa10_actual_final_bridge_xi110 import (
    KokunoPA10ActualFinalBridgeToXi110,
)


def test_contract_freezes_upstream_identity_and_final_gates() -> None:
    contract = public_contract()
    assert contract["agent1_pr"] == A1_FINAL_BRIDGE_PR == 940
    assert (
        contract["agent1_exact_head"]
        == A1_FINAL_BRIDGE_HEAD
        == "f430f8f97df3bfafce56bf094047c769ecec922d"
    )
    assert FROZEN_FINAL_GATES == {
        "normalized_momentum_max": 1.0e-3,
        "normalized_momentum_volume_l2": 1.0e-3,
        "divergence_max": 1.0e-5,
        "divergence_volume_l2": 1.0e-5,
    }
    assert contract["protocol"]["seed"] == 9173621
    assert tuple(contract["protocol"]["interior_halfwidths"]) == (0.12, 0.06, 0.03)
    assert tuple(contract["protocol"]["endpoint_steps"]) == (0.08, 0.04, 0.02)
    assert contract["truth_boundary"]["pde_validated"] is False


def test_scientific_entry_point_has_no_threshold_resolution_or_seed_knob() -> None:
    from openai_ns_reconstruction import (
        kokuno_a4_actual_final_bridge_xi110_independent_audit as module,
    )

    params = inspect.signature(module.materialize_actual_final_bridge_audit).parameters
    assert tuple(params) == ("field",)
    forbidden = {
        "seed",
        "resolution",
        "step",
        "halfwidth",
        "threshold",
        "momentum_threshold",
        "divergence_threshold",
        "forcing",
        "pressure",
        "residual",
    }
    assert forbidden.isdisjoint(params)


@pytest.mark.parametrize(
    "nodes",
    [INTERIOR_NODES, FORWARD_ENDPOINT_NODES, BACKWARD_ENDPOINT_NODES],
)
def test_polynomial_derivative_weights_are_degree_six_exact(nodes: np.ndarray) -> None:
    weights = first_derivative_weights(nodes)
    for degree in range(7):
        observed = float(np.sum(weights * nodes**degree))
        expected = 1.0 if degree == 1 else 0.0
        assert observed == pytest.approx(expected, abs=5.0e-12)

    x = np.asarray([0.31, 0.57, 0.83])
    eta = np.asarray([-0.2, 0.0, 0.27])

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
    observed = polynomial_first(poly, x, eta, 0.017, nodes, weights)
    np.testing.assert_allclose(observed, exact, rtol=0.0, atol=2.0e-9)


def test_upstream_loader_rejects_scientific_configuration_drift(tmp_path) -> None:
    obj = KokunoPA10ActualFinalBridgeToXi110()
    path = tmp_path / "candidate.json"
    payload = obj.save_configuration(path)
    mutated = dict(payload)
    mutated["angular_settle_log_width"] = float(
        mutated["angular_settle_log_width"]
    ) + 1.0e-4
    path.write_text(json.dumps(mutated, indent=2, sort_keys=True) + "\n")
    with pytest.raises(ValueError, match="configuration is frozen"):
        KokunoPA10ActualFinalBridgeToXi110.load_configuration(path)


def test_real_save_reload_actual_bridge_audit_materializes_scoped_receipt() -> None:
    receipt = run_real_artifact_audit()

    assert receipt["save_reload_exact_replay"] is True
    assert receipt["sample_count"] == 3 * PROTOCOL.random_per_region + 10
    assert receipt["endpoint_sample_count"] == 7
    assert all(receipt["mutation_checks"].values())
    assert isinstance(receipt["passed"], bool)

    for key in (
        "F_final_bridge_X",
        "U_final_bridge_X",
        "E_final_bridge_X",
    ):
        assert len(receipt["interior_resolution_metrics"][key]) == 3
        assert math.isfinite(
            receipt["interior_worst_witness"][key]["absolute_error"]
        )
        for endpoint in ("X100", "Xi"):
            assert len(receipt["endpoint_results"][endpoint]["metrics"][key]) == 3
            assert math.isfinite(
                receipt["endpoint_results"][endpoint]["worst_witness"][key][
                    "absolute_error"
                ]
            )

    assert receipt["truth_boundary"]["outer_global_leading_velocity_materialized"] is False
    assert receipt["truth_boundary"]["matched_global_pressure_materialized"] is False
    assert receipt["truth_boundary"]["restricted_forcing_materialized"] is False
    assert receipt["truth_boundary"]["correction_velocity_materialized"] is False
    assert receipt["truth_boundary"]["heldout_complete_ns_residual_assessed"] is False
    assert receipt["truth_boundary"]["pde_validated"] is False


def test_module_truth_boundary_never_preclaims_full_ns() -> None:
    assert TRUTH_BOUNDARY["actual_final_interpolation_100_to_Xi_materialized"] is True
    assert TRUTH_BOUNDARY["actual_final_bridge_independently_audited"] is False
    assert TRUTH_BOUNDARY["Xi_endpoint_radial_derivatives_independently_audited"] is False
    assert TRUTH_BOUNDARY["leading_only_ns_residual_assessed"] is False
    assert TRUTH_BOUNDARY["leading_plus_oscillatory_ns_residual_assessed"] is False
    assert TRUTH_BOUNDARY["after_correction_ns_residual_assessed"] is False
    assert TRUTH_BOUNDARY["pde_validated"] is False
