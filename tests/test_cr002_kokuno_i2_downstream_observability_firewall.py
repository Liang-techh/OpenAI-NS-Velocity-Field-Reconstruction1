from __future__ import annotations

import copy
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "src/openai_ns_reconstruction/audit_kokuno_i2_downstream_observability_firewall.py"
SPEC = importlib.util.spec_from_file_location("audit_i2_downstream", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)


def _contract() -> dict:
    return json.loads((ROOT / "configs/kokuno_i2_downstream_observability_firewall.json").read_text(encoding="utf-8"))


def test_contract_baseline_is_fail_closed_clean() -> None:
    assert AUDIT.audit_contract(_contract()) == []


def test_mechanics_witness_is_nonvacuous_and_not_candidate_data() -> None:
    witness = AUDIT.descendant_nonisolation_witness()
    assert witness["classification"] == "autonomous_mechanics_only_not_kokuno_or_openai_candidate_data"
    assert witness["i2_high_precision_delta_nonzero"] is True
    assert witness["i2_public_float64_unchanged"] is True
    assert witness["descendant_public_float64_changed"] is True
    assert witness["descendant_difference_does_not_isolate_i2"] is True


@pytest.mark.parametrize(
    "mutator",
    [
        lambda p: p.update(contract_id="wrong"),
        lambda p: p["target"].update(base_a5_head="0" * 40),
        lambda p: p["target"].update(origin_i2_head="1" * 40),
        lambda p: p["target"].update(origin_i1_parent_head="2" * 40),
        lambda p: p["target"].update(prior_cr002_head="3" * 40),
        lambda p: p["target"].update(prior_cr002_contract_blob="4" * 40),
        lambda p: p["downstream_lineage"]["agent1_i4"].update(head="5" * 40),
        lambda p: p["downstream_lineage"]["agent2_i4_composite"].update(head="6" * 40),
        lambda p: p["downstream_lineage"]["agent4_i4_divergence_audit"].update(head="7" * 40),
        lambda p: p["downstream_lineage"]["agent3_i2_stress"].update(head="8" * 40),
        lambda p: p["representation_truth"].update(i2_public_float64_delta_from_i1_parent_verified=True),
        lambda p: p["representation_truth"].update(i2_public_float64_delta_probe_protocol_executed=True),
        lambda p: p["representation_truth"].update(i2_public_float64_visual_effect_verified=True),
        lambda p: p["representation_truth"].update(downstream_semantic_progress_implies_i2_public_delta=True),
        lambda p: p["representation_truth"].update(downstream_save_load_implies_i2_public_delta=True),
        lambda p: p["representation_truth"].update(downstream_scoped_divergence_implies_i2_public_delta=True),
        lambda p: p["representation_truth"].update(downstream_i2_stress_implies_i2_public_delta=True),
        lambda p: p["representation_truth"].update(descendant_public_difference_from_i1_parent_isolates_i2_repair=True),
        lambda p: p["representation_truth"].update(i2_public_delta_implies_visual_correspondence=True),
        lambda p: p["representation_truth"].update(i2_public_delta_implies_pde_validation=True),
        lambda p: p["forbidden_implications"].remove("I4 scoped divergence audit -> I2 public float64 delta verified"),
        lambda p: p["provenance_classes"].pop("pending_or_unknown"),
        lambda p: p["cr001_freeze"].update(momentum_max=0.01),
        lambda p: p["cr001_freeze"].update(divergence_L2=0.001),
        lambda p: p["cr001_freeze"].update(residual_defined_free_force_forbidden=False),
        lambda p: p["cr001_freeze"].update(candidate_collapse_forbidden=False),
        lambda p: p["cr001_freeze"].update(posthoc_threshold_relaxation_forbidden=False),
        lambda p: p["independent_project_states"].update(canonical_eq45_velocity_export_ready=False),
        lambda p: p["independent_project_states"].update(canonical_eq45_pde_validated=True),
        lambda p: p["independent_project_states"].update(kokuno_current_route_velocity_export_ready=True),
        lambda p: p["independent_project_states"].update(kokuno_current_route_visual_correspondence_verified=True),
        lambda p: p["independent_project_states"].update(kokuno_current_route_pde_validated=True),
        lambda p: p["mechanics_witness"].update(classification="public_source_fact"),
        lambda p: p.update(truth_boundary="This changes velocity."),
    ],
)
def test_forbidden_promotions_and_scope_drift_fail_closed(mutator) -> None:
    payload = copy.deepcopy(_contract())
    mutator(payload)
    assert AUDIT.audit_contract(payload), payload


def test_promotion_contract_cannot_substitute_descendants_for_direct_i2_parent_audit() -> None:
    payload = copy.deepcopy(_contract())
    requirements = payload["promotion_requirements"]["i2_public_float64_delta_from_i1_parent_verified"]
    payload["promotion_requirements"]["i2_public_float64_delta_from_i1_parent_verified"] = [
        item.replace("#1061", "#1079").replace("#1051", "#1080") for item in requirements
    ]
    assert AUDIT.audit_contract(payload)


def test_promotion_contract_requires_ulp_and_negative_evidence() -> None:
    payload = copy.deepcopy(_contract())
    requirements = payload["promotion_requirements"]["i2_public_float64_delta_from_i1_parent_verified"]
    payload["promotion_requirements"]["i2_public_float64_delta_from_i1_parent_verified"] = [
        item.replace("ULP/spacing", "precision context").replace("zero/sub-ULP output as negative evidence", "small output")
        for item in requirements
    ]
    assert AUDIT.audit_contract(payload)


def test_repository_live_replay() -> None:
    violations = AUDIT.repository_live_violations(ROOT)
    assert violations == []
