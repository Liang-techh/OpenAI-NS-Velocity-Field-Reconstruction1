from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.audit_kokuno_pressure_increment_gradient_scope import (
    SCOPE_PATH,
    audit_registered_scope,
    mechanics_witness,
    validate_scope,
)


ROOT = Path(__file__).resolve().parents[1]


def _scope() -> dict:
    return json.loads((ROOT / SCOPE_PATH).read_text(encoding="utf-8"))


def test_registered_pressure_increment_gradient_scope_is_fail_closed() -> None:
    receipt = audit_registered_scope(repo_root=ROOT)
    assert receipt["ok"] is True
    assert receipt["agent5_source_blob"] == "c8c23490c62c981440e7e725f8345335df3e59f9"
    assert receipt["canonical_constraints_blob"] == "6c559e42895a606e2ef025ade4cb448966d75814"


def test_autonomous_witness_same_increment_but_different_cartesian_gradient() -> None:
    witness = mechanics_witness()
    assert witness["max_pressure_increment_gap"] == 0.0
    assert witness["max_pressure_increment_X_derivative_gap"] == 0.0
    assert witness["max_pressure_increment_eta_derivative_gap"] == 0.0
    assert witness["max_cartesian_gradient_gap"] == pytest.approx(0.3, abs=1.0e-15)
    assert witness["max_spatial_constant_gradient_gap"] == 0.0
    assert witness["candidate_evidence"] is False
    assert witness["kokuno_or_openai_parameter_evidence"] is False


@pytest.mark.parametrize(
    ("key", "value"),
    [
        ("pressure_increment_identity_is_absolute_pressure_identity", True),
        ("pressure_increment_derivative_audit_authorizes_cartesian_grad_p", True),
        ("eta_dependent_axis_datum_is_gradient_gauge", True),
        ("Pi0_eta_required_for_full_eta_pressure_derivative", False),
        ("coordinate_pullback_required_for_cartesian_grad_p", False),
        ("matched_global_pressure_required_before_complete_ns_residual", False),
        ("reference_increment_failure_blocks_existing_callable_eq45_delivery", True),
    ],
)
def test_representation_claim_promotions_are_rejected(key: str, value: bool) -> None:
    mutated = _scope()
    mutated["representation_rules"][key] = value
    with pytest.raises(ValueError, match="representation rule"):
        validate_scope(mutated, repo_root=ROOT)


def test_autonomous_mechanics_cannot_be_laundered_into_public_source_facts() -> None:
    mutated = _scope()
    item = mutated["provenance"]["autonomous_design"].pop(1)
    mutated["provenance"]["public_source_fact"].append(item)
    with pytest.raises(ValueError, match="public-source fact classification drift"):
        validate_scope(mutated, repo_root=ROOT)


def test_missing_pi0_cannot_be_removed_from_pending_unknown() -> None:
    mutated = _scope()
    mutated["provenance"]["pending_unknown"].pop(0)
    with pytest.raises(ValueError, match="pending/unknown classification drift"):
        validate_scope(mutated, repo_root=ROOT)


def test_mechanics_witness_cannot_be_promoted_to_candidate_or_source_parameter_evidence() -> None:
    mutated = _scope()
    mutated["mechanics_witness"]["candidate_evidence"] = True
    with pytest.raises(ValueError, match="candidate evidence"):
        validate_scope(mutated, repo_root=ROOT)

    mutated = _scope()
    mutated["mechanics_witness"]["kokuno_or_openai_parameter_evidence"] = True
    with pytest.raises(ValueError, match="source parameter evidence"):
        validate_scope(mutated, repo_root=ROOT)


def test_cr001_threshold_and_free_force_mutations_are_rejected() -> None:
    mutated = _scope()
    mutated["cr001_invariants"]["momentum_volume_l2_gate"] = 2.0e-3
    with pytest.raises(ValueError, match="momentum L2 gate drift"):
        validate_scope(mutated, repo_root=ROOT)

    mutated = _scope()
    mutated["cr001_invariants"]["divergence_volume_l2_gate"] = 2.0e-5
    with pytest.raises(ValueError, match="divergence L2 gate drift"):
        validate_scope(mutated, repo_root=ROOT)

    mutated = _scope()
    mutated["cr001_invariants"]["residual_defined_pointwise_free_force_forbidden"] = False
    with pytest.raises(ValueError, match="free-force firewall disabled"):
        validate_scope(mutated, repo_root=ROOT)


def test_scoped_pressure_increment_audit_cannot_promote_pde_or_exactness() -> None:
    for key in (
        "kokuno_pde_validated",
        "visual_correspondence_verified",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        mutated = _scope()
        mutated["truth_state"][key] = True
        with pytest.raises(ValueError, match="truth-state promotion"):
            validate_scope(mutated, repo_root=ROOT)


def test_upstream_identity_mutation_is_rejected() -> None:
    mutated = _scope()
    mutated["audited_upstream"]["agent5_source_blob"] = "0" * 40
    with pytest.raises(ValueError, match="A5 source-blob registration drift"):
        validate_scope(mutated, repo_root=ROOT)
