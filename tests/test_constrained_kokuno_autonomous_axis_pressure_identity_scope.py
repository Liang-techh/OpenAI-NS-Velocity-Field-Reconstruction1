from __future__ import annotations

import copy
import json
from pathlib import Path

from openai_ns_reconstruction.audit_kokuno_autonomous_axis_pressure_identity_scope import (
    A1_REL,
    A4_REL,
    CONFIG_REL,
    audit_scope,
    mechanics_witness,
)


ROOT = Path(__file__).resolve().parents[1]


def _config() -> dict:
    return json.loads((ROOT / CONFIG_REL).read_text(encoding="utf-8"))


def test_exact_registered_stack_passes_scope_audit() -> None:
    assert audit_scope(ROOT) == []


def test_nonuniqueness_witness_has_two_distinct_admissible_profiles() -> None:
    witness = mechanics_witness(epsilon=0.1)
    assert witness["alternative_obeys_seed_envelope"] is True
    assert witness["seed_derivative_sign_passes"] is True
    assert witness["alternative_derivative_sign_passes"] is True
    assert witness["seed_even_on_probe"] is True
    assert witness["alternative_even_on_probe"] is True
    assert witness["profiles_distinct"] is True
    assert witness["derivatives_distinct"] is True
    assert witness["max_profile_gap"] > 0.0
    assert witness["max_derivative_gap"] > 0.0


def test_autonomous_seed_cannot_be_laundered_into_public_source_fact() -> None:
    cfg = copy.deepcopy(_config())
    cfg["provenance"]["public_source_fact"].append(
        "Pi0_seed is the source-prepared datum because it saturates the envelope."
    )
    errors = audit_scope(ROOT, config_override=cfg)
    assert any("autonomous choice leaked into public_source_fact" in e for e in errors)


def test_admissibility_cannot_be_promoted_to_source_identity() -> None:
    cfg = copy.deepcopy(_config())
    cfg["truth_boundary"]["source_prepared_appendixA_Pi0_materialized"] = True
    cfg["truth_boundary"]["source_prepared_appendixA_Pi0_identified_with_autonomous_seed"] = True
    errors = audit_scope(ROOT, config_override=cfg)
    assert any("source_prepared_appendixA_Pi0_materialized" in e for e in errors)
    assert any("source_prepared_appendixA_Pi0_identified_with_autonomous_seed" in e for e in errors)


def test_representation_contract_cannot_weaken_identity_boundary() -> None:
    cfg = copy.deepcopy(_config())
    cfg["representation_contract"][
        "admissibility_conditions_satisfied_does_not_imply_source_identity"
    ] = False
    errors = audit_scope(ROOT, config_override=cfg)
    assert any("representation contract weakened" in e for e in errors)


def test_independent_audit_cannot_promote_pde_or_visual_claims() -> None:
    cfg = copy.deepcopy(_config())
    cfg["truth_boundary"]["pde_validated"] = True
    cfg["truth_boundary"]["visual_correspondence_verified"] = True
    cfg["truth_boundary"]["paper_exact"] = True
    cfg["truth_boundary"]["openai_field_identified"] = True
    errors = audit_scope(ROOT, config_override=cfg)
    for key in (
        "pde_validated",
        "visual_correspondence_verified",
        "paper_exact",
        "openai_field_identified",
    ):
        assert any(key in e for e in errors)


def test_cr001_gates_cannot_be_relaxed_on_this_seam() -> None:
    cfg = copy.deepcopy(_config())
    cfg["cr001_invariants"]["momentum_max"] = 0.01
    cfg["cr001_invariants"]["momentum_volume_l2"] = 0.01
    cfg["cr001_invariants"]["divergence_max"] = 0.001
    cfg["cr001_invariants"]["divergence_volume_l2"] = 0.001
    errors = audit_scope(ROOT, config_override=cfg)
    assert "CR001 frozen invariants changed" in errors


def test_free_force_and_amplitude_collapse_remain_forbidden() -> None:
    cfg = copy.deepcopy(_config())
    cfg["cr001_invariants"]["residual_defined_pointwise_free_force_allowed"] = True
    cfg["cr001_invariants"]["amplitude_collapse_allowed"] = True
    errors = audit_scope(ROOT, config_override=cfg)
    assert "CR001 frozen invariants changed" in errors


def test_pending_kokuno_pressure_cannot_block_eq45_velocity_delivery() -> None:
    cfg = copy.deepcopy(_config())
    cfg["truth_boundary"]["canonical_eq45_velocity_export_ready_remains_independent"] = False
    errors = audit_scope(ROOT, config_override=cfg)
    assert any("Eq45 callable velocity delivery" in e for e in errors)


def test_agent1_source_disclaimer_and_false_states_are_fail_closed() -> None:
    text = (ROOT / A1_REL).read_text(encoding="utf-8")
    mutated = text.replace(
        "does **not** pretend to recover it",
        "claims to recover it exactly",
        1,
    ).replace(
        '"source_prepared_appendixA_Pi0_materialized": False',
        '"source_prepared_appendixA_Pi0_materialized": True',
        1,
    )
    errors = audit_scope(ROOT, a1_text_override=mutated)
    assert any("source_prepared_appendixA_Pi0_materialized" in e for e in errors)
    assert any("identity disclaimer" in e for e in errors)


def test_agent4_audit_cannot_upgrade_autonomous_seed_provenance() -> None:
    text = (ROOT / A4_REL).read_text(encoding="utf-8")
    mutated = text.replace(
        "repository-autonomous rather than the source-prepared",
        "verified identical to the source-prepared",
        1,
    ).replace(
        '"pde_validated": False',
        '"pde_validated": True',
        1,
    )
    errors = audit_scope(ROOT, a4_text_override=mutated)
    assert any("pde_validated" in e for e in errors)
    assert any("identity disclaimer" in e for e in errors)
