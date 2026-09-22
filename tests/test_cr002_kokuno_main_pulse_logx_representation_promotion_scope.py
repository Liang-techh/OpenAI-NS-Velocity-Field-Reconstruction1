from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.audit_kokuno_main_pulse_logx_representation_promotion_scope import (
    audit_contract,
    representation_witness,
)

ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "configs/kokuno_main_pulse_logx_representation_promotion_scope.json"


def _contract():
    return json.loads(CONTRACT.read_text(encoding="utf-8"))


def test_repository_live_contract_is_consistent():
    assert audit_contract(ROOT) == []


def test_four_way_provenance_remains_exactly_partitioned():
    provenance = _contract()["four_way_provenance"]
    assert set(provenance) == {
        "user_requirements",
        "public_source_facts",
        "autonomous_design",
        "pending_unknown",
    }
    assert all(provenance[key] for key in provenance)


def test_representation_mechanics_witness_is_nonvacuous_and_not_scientific_evidence():
    witness = representation_witness()
    assert witness["classification"] == "autonomous_mechanics_only"
    assert witness["explicit_X_float64_representable"] is False
    assert witness["physical_radius_float64_finite"] is True
    assert witness["log_X"] > witness["log_float64_max"]


@pytest.mark.parametrize(
    ("section", "key", "bad_value"),
    [
        ("representation_transition", "parent_and_child_same_semantic_candidate_identity", True),
        ("representation_transition", "new_representation_identity_required", False),
        ("representation_transition", "source_exact_main_pulse_recovered", True),
        ("evidence_transfer", "parent_finite_x_receipts_prove_new_logx_only_interval", True),
        ("evidence_transfer", "a2_pr1108_composite_consumes_child_pr1107", True),
        ("evidence_transfer", "a4_pr1110_divergence_audits_child_pr1107", True),
        ("evidence_transfer", "a3_pr1109_correction_evidence_belongs_to_child_pr1107", True),
        ("evidence_transfer", "fresh_matching_a2_composition_required_for_child", False),
        ("evidence_transfer", "fresh_implementation_distinct_a4_audit_required_for_child", False),
        ("truth_states", "source_exact_amplitude_root_materialized", True),
        ("truth_states", "source_pulse_end_MJ_corrections_materialized", True),
        ("truth_states", "outer_global_leading_velocity_materialized", True),
        ("truth_states", "kokuno_velocity_export_ready", True),
        ("truth_states", "visual_correspondence_verified", True),
        ("truth_states", "pde_validated", True),
        ("truth_states", "paper_exact", True),
        ("truth_states", "openai_field_identified", True),
        ("truth_states", "blowup_proved", True),
        ("canonical_cr001", "residual_defined_free_forcing_forbidden", False),
        ("canonical_cr001", "candidate_collapse_forbidden", False),
        ("canonical_cr001", "post_hoc_threshold_relaxation_forbidden", False),
        ("canonical_delivery", "velocity_export_ready", False),
        ("canonical_delivery", "pde_validated", True),
    ],
)
def test_forbidden_scope_or_truth_mutations_fail_closed(section, key, bad_value):
    mutated = copy.deepcopy(_contract())
    mutated[section][key] = bad_value
    assert audit_contract(ROOT, mutated), (section, key, bad_value)


def test_threshold_mutation_fails_closed():
    mutated = copy.deepcopy(_contract())
    mutated["canonical_cr001"]["momentum_max"] = 0.01
    assert any("momentum_max" in item for item in audit_contract(ROOT, mutated))


def test_base_identity_mutation_fails_closed():
    mutated = copy.deepcopy(_contract())
    mutated["exact_stack"]["base_head"] = "0" * 40
    assert any("#1107 base head" in item for item in audit_contract(ROOT, mutated))


def test_representation_promotion_does_not_promote_scientific_states():
    c = _contract()
    assert c["representation_transition"]["parent_finite_x_gap_closed_by_child_logx_through_public_principal_xi11"]
    assert c["truth_states"]["logx_leading_only_public_principal_xi11_materialized"]
    assert not c["truth_states"]["source_exact_amplitude_root_materialized"]
    assert not c["truth_states"]["outer_global_leading_velocity_materialized"]
    assert not c["truth_states"]["pde_validated"]
    assert not c["truth_states"]["paper_exact"]
