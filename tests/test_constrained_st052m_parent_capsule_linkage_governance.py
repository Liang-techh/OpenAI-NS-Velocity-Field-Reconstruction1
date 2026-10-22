from __future__ import annotations

import copy
import hashlib
import json

import pytest

from openai_ns_reconstruction.constrained_st052m_parent_capsule_linkage_governance import (
    audit_contract,
    load_contract,
    verify_candidate_validation_binding,
)


def test_governance_contract_passes() -> None:
    receipt = audit_contract()
    assert receipt["status"] == "pass"
    assert receipt["upstream_pr"] == 619
    assert receipt["candidate_validation_sha_bound_in_exact_artifact"] is True
    assert receipt["builder_cross_file_binding_gap_recorded"] is True
    assert receipt["complete_candidate_save_load_ready"] is False
    assert receipt["pde_validated"] is False


def _write_pair(tmp_path):
    candidate = tmp_path / "candidate.json"
    candidate.write_text('{"candidate":"synthetic-binding-test"}\n', encoding="utf-8")
    candidate_sha = hashlib.sha256(candidate.read_bytes()).hexdigest()
    validation = tmp_path / "validation.json"
    validation.write_text(
        json.dumps(
            {
                "candidate_sha256": candidate_sha,
                "pde_validated": False,
                "all_numeric_gates_pass": False,
                "gates": {"momentum_max": False, "momentum_L2": False},
            }
        ),
        encoding="utf-8",
    )
    return candidate, validation, candidate_sha


def test_candidate_validation_pair_binding_passes_then_rejects_substitution(tmp_path) -> None:
    candidate, validation, candidate_sha = _write_pair(tmp_path)
    receipt = verify_candidate_validation_binding(candidate, validation)
    assert receipt["candidate_sha256"] == candidate_sha
    assert receipt["validation_candidate_sha256"] == candidate_sha
    assert receipt["candidate_validation_sha_pair_matches"] is True

    candidate.write_text('{"candidate":"substituted-after-validation"}\n', encoding="utf-8")
    with pytest.raises(ValueError, match="candidate/validation SHA mismatch"):
        verify_candidate_validation_binding(candidate, validation)


@pytest.mark.parametrize(
    ("mutator", "message"),
    [
        (
            lambda c: c["governance_finding"].__setitem__(
                "candidate_validation_sha_binding_required_before_manifest_claims_saved_file_validation",
                False,
            ),
            "binding prerequisite removed",
        ),
        (
            lambda c: c["delivery_state"].__setitem__(
                "experimental_st052_temporal_velocity_export_ready", True
            ),
            "premature promotion",
        ),
        (
            lambda c: c["exact_artifact_pair"].__setitem__("validation_seed", 914027),
            "diagnostic seed drift",
        ),
        (
            lambda c: c["cr001_snapshot"].__setitem__("momentum_max", 0.01),
            "momentum threshold drift",
        ),
        (
            lambda c: c["exact_artifact_pair"].__setitem__(
                "replay_identity_scope", "complete_candidate_file_identity"
            ),
            "replay identity scope laundering",
        ),
    ],
)
def test_governance_mutations_fail_closed(mutator, message) -> None:
    contract = copy.deepcopy(load_contract())
    mutator(contract)
    with pytest.raises(ValueError, match=message):
        audit_contract(contract)
