from __future__ import annotations

from dataclasses import replace
import inspect

import numpy as np
import pytest

from openai_ns_reconstruction import kokuno_a4_blackbox_full_ns_validator as parent
from openai_ns_reconstruction import kokuno_a4_canonical_quadrature_evidence as gate


def _contract():
    return parent.build_fixed_physical_contract(
        contract_id="a4-quadrature-test-contract",
        matched_pressure_identity_sha256="pressure-v1",
        restricted_forcing_parameters_sha256="force-params-v1",
        restricted_forcing_preregistration_sha256="force-preregistered-v1",
    )


def _levels(
    *,
    coarse_momentum=6.0e-4,
    medium_momentum=5.0e-4,
    fine_momentum=4.9e-4,
    coarse_divergence=6.0e-6,
    medium_divergence=5.0e-6,
    fine_divergence=4.9e-6,
):
    values = (
        (24, coarse_momentum, coarse_divergence),
        (48, medium_momentum, medium_divergence),
        (96, fine_momentum, fine_divergence),
    )
    out = []
    for order, momentum, divergence in values:
        rows = tuple(
            gate.QuadratureTimeMetrics(
                time=float(time_value),
                momentum_volume_l2=float(momentum),
                divergence_volume_l2=float(divergence),
                kinetic_energy=1.0,
            )
            for time_value in parent.VALIDATION_TIMES
        )
        out.append(gate.QuadratureOrderMetrics(order_per_axis=order, by_time=rows))
    return tuple(out)


def _evidence(*, candidate_sha="candidate-a", contract=None, levels=None):
    if contract is None:
        contract = _contract()
    if levels is None:
        levels = _levels()
    return gate.build_canonical_quadrature_evidence(
        candidate_sha256=candidate_sha,
        stage="leading_only",
        physical_contract_sha256=contract.physical_contract_sha256,
        levels=levels,
    )


class SolidRotationCandidate:
    def __init__(self, contract_sha: str):
        self.contract_sha = contract_sha

    def velocity(self, x, y, z, t):
        xb, yb, zb, _ = np.broadcast_arrays(x, y, z, t)
        return np.stack((-yb, xb, np.zeros_like(zb)), axis=-1).astype(float)

    def pressure(self, x, y, z, t):
        xb, yb, _, _ = np.broadcast_arrays(x, y, z, t)
        return (0.5 * (xb * xb + yb * yb)).astype(float)

    def forcing(self, x, y, z, t):
        xb, _, _, _ = np.broadcast_arrays(x, y, z, t)
        return np.zeros(xb.shape + (3,), dtype=float)

    def validation_metadata(self):
        return {
            "candidate_id": "solid-rotation-mechanics",
            "candidate_sha256": "solid-rotation-sha",
            "physical_contract_sha256": self.contract_sha,
            "stage": "leading_only",
            "global_leading_velocity_materialized": True,
            "outer_join_materialized": True,
            "matched_pressure_included": True,
            "restricted_forcing_included": True,
            "restricted_forcing_preregistered": True,
            "residual_defined_forcing": False,
            "whole_domain_support_materialized": True,
            "leading_included": True,
            "oscillatory_included": False,
            "correction_included": False,
        }


def test_typed_receipt_binds_exact_candidate_contract_operator_protocol_and_levels():
    contract = _contract()
    evidence = _evidence(contract=contract)
    assert evidence.receipt_sha256 == evidence.recomputed_sha256()
    assert evidence.validator_operator_sha256 == gate.VALIDATOR_OPERATOR_SHA256
    assert evidence.heldout_protocol_sha256 == gate.HELDOUT_PROTOCOL_SHA256
    assert evidence.derivative_steps == parent.DERIVATIVE_STEPS
    assert evidence.validation_times == parent.VALIDATION_TIMES
    assert evidence.quadrature_orders_per_axis == (24, 48, 96)
    assert tuple(level.order_per_axis for level in evidence.levels) == (24, 48, 96)
    assert evidence.max_scaled_change_48_to_96 <= gate.QUADRATURE_STABILITY_GATE
    assert gate.canonical_quadrature_failures(
        evidence,
        candidate_sha256="candidate-a",
        stage="leading_only",
        physical_contract_sha256=contract.physical_contract_sha256,
    ) == []


def test_checksum_valid_candidate_and_physical_contract_reuse_are_rejected():
    contract = _contract()
    evidence = _evidence(contract=contract)

    changed_candidate = replace(evidence, candidate_sha256="candidate-b", receipt_sha256="")
    changed_candidate = replace(
        changed_candidate, receipt_sha256=changed_candidate.recomputed_sha256()
    )
    failures = gate.canonical_quadrature_failures(
        changed_candidate,
        candidate_sha256="candidate-a",
        stage="leading_only",
        physical_contract_sha256=contract.physical_contract_sha256,
    )
    assert "quadrature_candidate_identity" in failures
    assert "quadrature_receipt_digest" not in failures

    changed_contract = replace(evidence, physical_contract_sha256="other-contract", receipt_sha256="")
    changed_contract = replace(
        changed_contract, receipt_sha256=changed_contract.recomputed_sha256()
    )
    failures = gate.canonical_quadrature_failures(
        changed_contract,
        candidate_sha256="candidate-a",
        stage="leading_only",
        physical_contract_sha256=contract.physical_contract_sha256,
    )
    assert "quadrature_physical_contract_identity" in failures
    assert "quadrature_receipt_digest" not in failures


def test_operator_threshold_and_receipt_identity_mutations_fail_closed():
    contract = _contract()
    evidence = _evidence(contract=contract)

    operator_mutation = replace(evidence, validator_operator_sha256="0" * 64, receipt_sha256="")
    operator_mutation = replace(
        operator_mutation, receipt_sha256=operator_mutation.recomputed_sha256()
    )
    failures = gate.canonical_quadrature_failures(
        operator_mutation,
        candidate_sha256="candidate-a",
        stage="leading_only",
        physical_contract_sha256=contract.physical_contract_sha256,
    )
    assert "quadrature_operator_identity" in failures

    loosened = replace(
        evidence,
        thresholds=replace(evidence.thresholds, momentum_volume_l2=1.001e-3),
        receipt_sha256="",
    )
    loosened = replace(loosened, receipt_sha256=loosened.recomputed_sha256())
    failures = gate.canonical_quadrature_failures(
        loosened,
        candidate_sha256="candidate-a",
        stage="leading_only",
        physical_contract_sha256=contract.physical_contract_sha256,
    )
    assert "quadrature_threshold_identity" in failures

    malformed_digest = replace(evidence, receipt_sha256="f" * 64)
    failures = gate.canonical_quadrature_failures(
        malformed_digest,
        candidate_sha256="candidate-a",
        stage="leading_only",
        physical_contract_sha256=contract.physical_contract_sha256,
    )
    assert "quadrature_receipt_digest" in failures


def test_finest_canonical_volume_l2_gates_remain_exactly_1e3_and_1e5():
    contract = _contract()
    bad_momentum = _levels(
        coarse_momentum=1.001e-3,
        medium_momentum=1.001e-3,
        fine_momentum=1.001e-3,
    )
    evidence = _evidence(contract=contract, levels=bad_momentum)
    failures = gate.canonical_quadrature_failures(
        evidence,
        candidate_sha256="candidate-a",
        stage="leading_only",
        physical_contract_sha256=contract.physical_contract_sha256,
    )
    assert any(item.startswith("quadrature_momentum_volume_l2@") for item in failures)

    bad_divergence = _levels(
        coarse_divergence=1.001e-5,
        medium_divergence=1.001e-5,
        fine_divergence=1.001e-5,
    )
    evidence = _evidence(contract=contract, levels=bad_divergence)
    failures = gate.canonical_quadrature_failures(
        evidence,
        candidate_sha256="candidate-a",
        stage="leading_only",
        physical_contract_sha256=contract.physical_contract_sha256,
    )
    assert any(item.startswith("quadrature_divergence_volume_l2@") for item in failures)
    assert gate.FrozenQuadratureThresholds().momentum_volume_l2 == 1e-3
    assert gate.FrozenQuadratureThresholds().divergence_volume_l2 == 1e-5


def test_48_to_96_quadrature_stability_is_an_independent_fail_closed_gate():
    contract = _contract()
    unstable = _levels(
        coarse_momentum=7.0e-4,
        medium_momentum=5.0e-4,
        fine_momentum=1.0e-4,
    )
    evidence = _evidence(contract=contract, levels=unstable)
    assert evidence.max_scaled_change_48_to_96 > gate.QUADRATURE_STABILITY_GATE
    failures = gate.canonical_quadrature_failures(
        evidence,
        candidate_sha256="candidate-a",
        stage="leading_only",
        physical_contract_sha256=contract.physical_contract_sha256,
    )
    assert "quadrature_resolution_stability" in failures


def test_successor_validator_has_no_bare_quadrature_boolean_and_missing_evidence_cannot_promote():
    signature = inspect.signature(gate.validate_complete_candidate_with_canonical_quadrature)
    assert "quadrature_ladder_assessed" not in signature.parameters
    contract = _contract()
    candidate = SolidRotationCandidate(contract.physical_contract_sha256)
    x = np.asarray((0.1, -0.2, 0.3), dtype=float)
    y = np.asarray((0.2, 0.1, -0.1), dtype=float)
    z = np.asarray((0.0, 0.1, -0.2), dtype=float)
    t = np.asarray((0.25, 0.50, 0.75), dtype=float)
    weights = np.ones_like(x)

    receipt = gate.validate_complete_candidate_with_canonical_quadrature(
        candidate,
        contract,
        x,
        y,
        z,
        t,
        expected_stage="leading_only",
        volume_weights=weights,
    )
    assert receipt.canonical_quadrature_assessed is False
    assert "canonical_quadrature_evidence_missing" in receipt.quadrature_gate_failures
    assert receipt.pde_validated is False

    with pytest.raises(TypeError):
        gate.validate_complete_candidate_with_canonical_quadrature(
            candidate,
            contract,
            x,
            y,
            z,
            t,
            expected_stage="leading_only",
            volume_weights=weights,
            quadrature_ladder_assessed=True,
        )


def test_public_contract_records_cr002_response_and_current_candidate_remains_ineligible():
    contract = gate.public_contract()
    assert contract["parent_a4_pr"] == 896
    assert contract["cr002_audit_pr"] == 898
    assert contract["quadrature_orders_per_axis"] == [24, 48, 96]
    assert contract["derivative_steps"] == [0.02, 0.01, 0.005]
    assert contract["frozen_gates"]["momentum_volume_l2"] == 1e-3
    assert contract["frozen_gates"]["divergence_volume_l2"] == 1e-5
    assert contract["bare_quadrature_ladder_assessed_parameter_present"] is False
    assert contract["caller_sample_weights_are_canonical_quadrature"] is False
    assert contract["typed_canonical_quadrature_evidence_required_for_pde_promotion"] is True
    assert contract["current_ineligibility"]["canonical_quadrature_evidence_available"] is False
    assert contract["current_ineligibility"]["complete_blackbox_candidate_available"] is False
    assert contract["truth_boundary"]["real_canonical_quadrature_run_completed"] is False
    assert contract["truth_boundary"]["pde_validated"] is False
