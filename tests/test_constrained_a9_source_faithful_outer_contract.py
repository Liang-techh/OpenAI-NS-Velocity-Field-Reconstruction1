from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
CONTRACT_PATH = ROOT / "configs" / "external_openai_source_faithful_outer_v1.json"


def _load() -> dict:
    return json.loads(CONTRACT_PATH.read_text(encoding="utf-8"))


def _validate(contract: dict) -> None:
    assert contract["schema"] == "a9_external_source_outer_crosswalk_v1"
    assert contract["task_id"] == "CR-A9-099"

    source = contract["source_identity"]
    assert source["workbench_repository"] == "KokunoYumeto/yang-mills-interacting-workbench"
    assert source["workbench_commit"] == "9fd8c271372b5f53bb7f5613118794d6ecd0250f"
    assert source["source_record_path"] == "navier-stokes/SOURCE_READING_20260920.json"
    assert source["source_record_schema"] == "pinned-external-mathematical-source-v1"
    assert source["primary_paper_pages"] == 166
    assert source["primary_paper_sha256"] == "0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f"
    assert source["source_faithful_reconstruction_commit"] == "5e162f34cd2d3581f890660e81fbf063509085d0"
    assert source["source_faithful_pdf_sha256"] == "67bca97d638868c2fcb34b0ef1acdc4918cfcbb1a210afe98f0019a3c93fe8d9"
    assert source["workbench_mirror_zip_sha256"] == "997e825e34fa825d189bd301edaa6251e5505710fe294bec99c83a9712558e4c"
    assert source["independent_reconstruction"] is True
    assert source["official_source_release"] is False
    assert source["independent_mathematical_proof_verification"] is False
    assert source["distinct_from_208_page_analytical_reader"] is True

    migration = contract["license_and_migration"]
    assert migration["license_status"] == "no_identified_license_for_paper_or_source_faithful_reconstruction"
    assert migration["source_record_assigns_new_license"] is False
    assert migration["direct_source_code_or_text_copy_allowed_by_this_screen"] is False
    assert migration["classification"] == "reimplement_math_only"
    assert migration["migration_scope"] == "equations_and_structural_schedule_only_no_direct_source_code_or_text_copy"

    crosswalk = contract["directly_useful_crosswalk"]
    assert [entry["source_locator"] for entry in crosswalk] == [
        "Section 4, Eq. 4.5",
        "Section 4, Eq. 4.7",
        "Section 4, Lemma 4.8",
        "Appendix A, Eq. A.5 and Eq. A.9-A.13",
    ]
    assert crosswalk[-1]["migration_class"] == "priority_reimplementation"
    assert crosswalk[-1]["role"] == "post_power_law_outer_completion_schedule"

    governance = contract["constraint_governance"]
    assert governance == {
        "canonical_constraints_path": "configs/constraints.json",
        "canonical_constraints_unchanged": True,
        "forcing_family_changed": False,
        "normalization_changed": False,
        "validation_split_changed": False,
        "residual_thresholds_changed": False,
        "visual_observables_used_as_pde_evidence": False,
    }

    truth = contract["truth_boundary"]
    assert truth
    assert all(value is False for value in truth.values())


def test_source_faithful_outer_contract_is_pinned_and_fail_closed() -> None:
    _validate(_load())


@pytest.mark.parametrize(
    ("section", "key", "replacement"),
    [
        ("source_identity", "workbench_commit", "0" * 40),
        ("source_identity", "official_source_release", True),
        ("license_and_migration", "direct_source_code_or_text_copy_allowed_by_this_screen", True),
        ("license_and_migration", "classification", "direct_migration"),
        ("constraint_governance", "residual_thresholds_changed", True),
        ("truth_boundary", "pde_validated", True),
        ("truth_boundary", "source_correspondence_verified", True),
    ],
)
def test_contract_rejects_identity_license_or_truth_promotion(
    section: str, key: str, replacement: object
) -> None:
    mutated = copy.deepcopy(_load())
    mutated[section][key] = replacement
    with pytest.raises(AssertionError):
        _validate(mutated)
