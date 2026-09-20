from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.audit_kokuno_quadrature_evidence_identity_scope import (
    EXPECTED_OBSERVED_A4,
    audit,
    observed_agent4_semantics,
    quadrature_receipt_identity_failures,
)


ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "configs/kokuno_quadrature_evidence_identity_scope.json"
AGENT4_SOURCE = (
    ROOT / "src/openai_ns_reconstruction/kokuno_a4_blackbox_full_ns_validator.py"
)


def _identity_complete_receipt() -> dict[str, object]:
    return {
        "receipt_schema": "canonical-cr001-quadrature-receipt-v1",
        "receipt_sha256": "a" * 64,
        "candidate_sha256": "b" * 64,
        "physical_contract_sha256": "c" * 64,
        "validator_operator_identity": "independent-blackbox-fd4-v1",
        "derivative_step_ladder": [0.02, 0.01, 0.005],
        "validation_times": [0.25, 0.3125, 0.4375, 0.5625, 0.6875, 0.75],
        "quadrature_orders_per_axis": [24, 48, 96],
        "momentum_volume_l2_by_time_and_order": {"mechanics_only": []},
        "divergence_volume_l2_by_time_and_order": {"mechanics_only": []},
        "energy_by_time_and_order": {"mechanics_only": []},
        "convergence_assessment": {"mechanics_only": True},
        "thresholds_bound": {
            "momentum_volume_l2": 0.001,
            "divergence_volume_l2": 1e-05,
            "reference_energy_abs_tolerance": 0.001,
        },
        "independent_of_training_and_held_in_selection": True,
        "residual_defined_forcing_forbidden": True,
    }


def test_audit_passes_on_exact_agent4_parent() -> None:
    result = audit(ROOT)
    assert result["agent4_source_blob_sha"] == (
        "a9e46daf60d909a63261dd10821e1a08101dd872"
    )
    assert result["constraints_blob_sha"] == (
        "6c559e42895a606e2ef025ade4cb448966d75814"
    )
    assert result["bare_boolean_receipt_failures"] == ["typed_receipt_required"]
    assert result["canonical_quadrature_receipt_exists"] is False
    assert result["pde_validated"] is False
    assert result["velocity_or_scientific_state_changed"] is False


def test_agent4_source_exposes_the_scoped_boolean_and_weight_seam() -> None:
    source = AGENT4_SOURCE.read_text(encoding="utf-8")
    assert observed_agent4_semantics(source) == EXPECTED_OBSERVED_A4


def test_bare_boolean_is_not_a_quadrature_receipt() -> None:
    failures = quadrature_receipt_identity_failures(
        True,
        expected_candidate_sha256="b" * 64,
        expected_physical_contract_sha256="c" * 64,
    )
    assert failures == ("typed_receipt_required",)


@pytest.mark.parametrize("value", [False, None, 0, 1, "assessed", [24, 48, 96]])
def test_untyped_quadrature_assertions_fail_closed(value: object) -> None:
    failures = quadrature_receipt_identity_failures(
        value,
        expected_candidate_sha256="b" * 64,
        expected_physical_contract_sha256="c" * 64,
    )
    assert failures == ("typed_receipt_required",)


def test_identity_complete_mechanics_receipt_only_passes_identity_checks() -> None:
    receipt = _identity_complete_receipt()
    failures = quadrature_receipt_identity_failures(
        receipt,
        expected_candidate_sha256="b" * 64,
        expected_physical_contract_sha256="c" * 64,
    )
    assert failures == ()
    # This helper intentionally checks identity/provenance shape only.  The
    # synthetic values above are mechanics placeholders and are not scientific
    # evidence or a PDE-validation receipt.
    assert receipt["convergence_assessment"] == {"mechanics_only": True}


def test_receipt_cannot_be_reused_across_candidate_identity() -> None:
    receipt = _identity_complete_receipt()
    failures = quadrature_receipt_identity_failures(
        receipt,
        expected_candidate_sha256="d" * 64,
        expected_physical_contract_sha256="c" * 64,
    )
    assert "candidate_sha256_mismatch" in failures


def test_receipt_cannot_be_reused_across_physical_contract_identity() -> None:
    receipt = _identity_complete_receipt()
    failures = quadrature_receipt_identity_failures(
        receipt,
        expected_candidate_sha256="b" * 64,
        expected_physical_contract_sha256="d" * 64,
    )
    assert "physical_contract_sha256_mismatch" in failures


@pytest.mark.parametrize(
    ("field", "value", "failure"),
    [
        ("derivative_step_ladder", [0.02, 0.01], "derivative_step_ladder_drift"),
        ("validation_times", [0.25, 0.75], "validation_times_drift"),
        ("quadrature_orders_per_axis", [24, 48], "quadrature_ladder_drift"),
        (
            "thresholds_bound",
            {
                "momentum_volume_l2": 0.002,
                "divergence_volume_l2": 1e-05,
                "reference_energy_abs_tolerance": 0.001,
            },
            "threshold_binding_drift",
        ),
        (
            "independent_of_training_and_held_in_selection",
            False,
            "independence_not_bound",
        ),
        (
            "residual_defined_forcing_forbidden",
            False,
            "residual_defined_forcing_not_forbidden",
        ),
    ],
)
def test_receipt_rejects_protocol_or_governance_drift(
    field: str, value: object, failure: str
) -> None:
    receipt = _identity_complete_receipt()
    receipt[field] = value
    failures = quadrature_receipt_identity_failures(
        receipt,
        expected_candidate_sha256="b" * 64,
        expected_physical_contract_sha256="c" * 64,
    )
    assert failure in failures


def test_receipt_requires_all_identity_and_metric_surfaces() -> None:
    receipt = _identity_complete_receipt()
    for key in (
        "receipt_sha256",
        "candidate_sha256",
        "physical_contract_sha256",
        "validator_operator_identity",
        "momentum_volume_l2_by_time_and_order",
        "divergence_volume_l2_by_time_and_order",
        "energy_by_time_and_order",
        "convergence_assessment",
    ):
        mutated = copy.deepcopy(receipt)
        mutated.pop(key)
        failures = quadrature_receipt_identity_failures(
            mutated,
            expected_candidate_sha256="b" * 64,
            expected_physical_contract_sha256="c" * 64,
        )
        assert f"missing:{key}" in failures


def test_machine_contract_keeps_four_provenance_classes_and_false_promotions() -> None:
    config = json.loads(CONFIG.read_text(encoding="utf-8"))
    assert set(config["provenance"]) == {
        "user_requirement",
        "public_source_fact",
        "autonomous_design",
        "pending_unknown",
    }
    assert config["provenance"]["public_source_fact"] == []
    assert config["mechanics_witness"]["candidate_scientific_evidence"] is False
    assert config["mechanics_witness"]["canonical_quadrature_was_run"] is False
    assert config["promotion_guards"] == {
        "bare_boolean_is_canonical_quadrature_evidence": False,
        "caller_supplied_sample_weights_are_canonical_24_48_96_quadrature": False,
        "sampled_or_monte_carlo_l2_may_replace_canonical_volume_l2": False,
        "quadrature_ladder_assessed_true_alone_may_authorize_pde_validated": False,
        "quadrature_evidence_may_be_reused_across_candidate_sha": False,
        "quadrature_evidence_may_be_reused_across_physical_contract_sha": False,
        "quadrature_evidence_may_omit_operator_and_derivative_identity": False,
        "quadrature_evidence_may_omit_per_level_metrics": False,
        "pde_pending_blocks_callable_velocity_delivery": False,
        "validator_contract_registration_implies_pde_validated": False,
    }
    truth = config["truth_boundary"]
    assert truth["candidate_bytes_changed"] is False
    assert truth["velocity_changed"] is False
    assert truth["pressure_changed"] is False
    assert truth["forcing_changed"] is False
    assert truth["threshold_changed"] is False
    assert truth["canonical_quadrature_receipt_exists"] is False
    assert truth["real_kokuno_complete_ns_residual_assessed"] is False
    assert truth["velocity_export_ready_promoted"] is False
    assert truth["visual_correspondence_verified"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False
