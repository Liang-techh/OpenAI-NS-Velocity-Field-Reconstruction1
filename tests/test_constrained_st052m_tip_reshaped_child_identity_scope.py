from __future__ import annotations

from copy import deepcopy

import pytest

from openai_ns_reconstruction.constrained_st052m_tip_reshaped_child_identity_scope import (
    GovernanceError,
    audit,
    load_contract,
    repository_root,
)


def test_live_scope_audits() -> None:
    result = audit()
    assert result["audited_pr"] == 740
    assert result["audited_head"] == "21019df839557e64b8d8689f68cdf8fc685a3a36"
    assert result["procedural_velocity_evaluable"] is True
    assert result["standalone_saved_loadable_child_ready"] is False
    assert result["fresh_evidence_requires_frozen_child_identity"] is True
    assert result["pde_validated"] is False


@pytest.mark.parametrize(
    ("section", "key", "value"),
    [
        ("procedural_child_scope", "standalone_candidate_payload_emitted", True),
        ("procedural_child_scope", "standalone_candidate_load_api_provided", True),
        ("procedural_child_scope", "standalone_candidate_identity_digest_bound", True),
        ("claim_states", "standalone_saved_loadable_child_ready", True),
        ("claim_states", "velocity_export_ready_promoted_by_740", True),
        ("claim_states", "visual_correspondence_verified", True),
        ("claim_states", "pde_validated", True),
        ("claim_states", "paper_exact", True),
        ("claim_states", "openai_field_identified", True),
        (
            "identity_rule",
            "740_targeted_and_historical_development_evidence_is_independent_visual_validation",
            True,
        ),
        (
            "identity_rule",
            "fresh_path_or_render_evidence_may_promote_visual_correspondence_without_public_target_protocol",
            True,
        ),
    ],
)
def test_claim_and_identity_promotions_fail_closed(section: str, key: str, value: bool) -> None:
    contract = deepcopy(load_contract())
    contract[section][key] = value
    with pytest.raises(GovernanceError):
        audit(contract=contract)


def test_fresh_path_must_not_be_reclassified_as_alpha_selection_data() -> None:
    contract = deepcopy(load_contract())
    contract["procedural_child_scope"]["fresh_714_data_used_for_alpha"] = True
    with pytest.raises(GovernanceError):
        audit(contract=contract)


def test_source_mutation_fresh_data_leak_fails_closed() -> None:
    root = repository_root()
    source = (root / "experiments/root_st052/agent7_st052m_tip_reshaped_nonlinear_replay.py").read_text(
        encoding="utf-8"
    )
    mutated = source.replace('"fresh_714_data_used": False', '"fresh_714_data_used": True')
    assert mutated != source
    with pytest.raises(GovernanceError):
        audit(replay_source=mutated)


def test_source_mutation_candidate_formula_binding_fails_closed() -> None:
    root = repository_root()
    source = (root / "experiments/root_st052/agent7_st052m_tip_reshaped_nonlinear_replay.py").read_text(
        encoding="utf-8"
    )
    mutated = source.replace(
        "candidate_fn = lambda p, t: reshaped_velocity(child_fn, p, t, alpha_94)",
        "candidate_fn = lambda p, t: child_fn(p, t)",
    )
    assert mutated != source
    with pytest.raises(GovernanceError):
        audit(replay_source=mutated)


def test_cr001_free_force_and_threshold_drift_fail_closed() -> None:
    contract = deepcopy(load_contract())
    contract["cr001_lock"]["residual_defined_or_pointwise_free_force_allowed"] = True
    with pytest.raises(GovernanceError):
        audit(contract=contract)

    contract = deepcopy(load_contract())
    contract["cr001_lock"]["pde_residual_max"] = 0.01
    with pytest.raises(GovernanceError):
        audit(contract=contract)


def test_frozen_identity_minimum_fields_cannot_be_weakened() -> None:
    contract = deepcopy(load_contract())
    contract["identity_rule"]["minimum_frozen_identity_fields"].remove("numeric alpha_94")
    with pytest.raises(GovernanceError):
        audit(contract=contract)
