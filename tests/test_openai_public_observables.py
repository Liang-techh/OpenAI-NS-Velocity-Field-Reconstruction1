from __future__ import annotations

from copy import deepcopy

import pytest

from openai_ns_reconstruction.openai_public_observables import (
    contract_sha256,
    load_openai_public_observable_contract,
    observable_contract_report,
    validate_openai_public_observable_contract,
)


EXPECTED_IDS = [
    "vortex_swirl_presence",
    "inward_spiraling_trajectories",
    "axial_stretching_trajectories",
    "shrinking_accelerating_core",
    "angular_rotation_variation",
    "radius_dependent_circulation",
]
EXPECTED_EXTERNAL_REPO = "scikit-image/scikit-image"
EXPECTED_EXTERNAL_COMMIT = "cc0a4b16ebf655873cd6f770e897460abdade288"
EXPECTED_EXTERNAL_LICENSE = "BSD-3-Clause (project default; individual vendored files may differ)"


def test_public_observable_contract_is_qualitative_and_fail_closed() -> None:
    contract = load_openai_public_observable_contract()
    assert [item["id"] for item in contract["observables"]] == EXPECTED_IDS
    assert all(item["numerical_target"] is None for item in contract["observables"])
    screening = contract["external_method_screening"][0]
    assert screening["source_repo"] == EXPECTED_EXTERNAL_REPO
    assert screening["source_commit"] == EXPECTED_EXTERNAL_COMMIT
    assert screening["license"] == EXPECTED_EXTERNAL_LICENSE
    assert screening["candidate_scope"].strip()
    assert screening["classification"] == "screened_not_adopted"
    assert screening["migration_scope"] == "none"
    assert screening["reason"].strip()
    assert all(value is False for value in contract["truth_boundary"].values())

    report = observable_contract_report()
    assert report["observable_ids"] == EXPECTED_IDS
    assert report["external_method_source_repo"] == EXPECTED_EXTERNAL_REPO
    assert report["external_method_source_commit"] == EXPECTED_EXTERNAL_COMMIT
    assert report["external_method_license"] == EXPECTED_EXTERNAL_LICENSE
    assert report["external_method_candidate_scope"] == screening["candidate_scope"]
    assert report["external_method_classification"] == "screened_not_adopted"
    assert report["external_method_migration_scope"] == "none"
    assert report["external_method_reason"] == screening["reason"]
    assert report["numeric_visual_target_defined"] is False
    assert report["pde_validated"] is False
    assert report["visual_correspondence_verified"] is False
    assert report["openai_field_identified"] is False
    assert report["contract_sha256"] == contract_sha256(contract)
    assert len(report["contract_sha256"]) == 64


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        ("numeric_target", "numerical targets"),
        ("pde_promotion", "truth flags"),
        ("visual_promotion", "truth flags"),
        ("pixel_tool_promotion", "silently promoted"),
        ("source_drift", "source identity drift"),
        ("external_license_drift", "license drift"),
        ("external_scope_missing", "candidate scope missing"),
        ("external_reason_missing", "rationale missing"),
    ],
)
def test_public_observable_contract_rejects_claim_laundering(mutation: str, message: str) -> None:
    contract = deepcopy(load_openai_public_observable_contract())
    if mutation == "numeric_target":
        contract["observables"][1]["numerical_target"] = 0.75
    elif mutation == "pde_promotion":
        contract["truth_boundary"]["pde_validated"] = True
    elif mutation == "visual_promotion":
        contract["truth_boundary"]["visual_correspondence_verified"] = True
    elif mutation == "pixel_tool_promotion":
        contract["external_method_screening"][0]["classification"] = "direct_migration"
    elif mutation == "source_drift":
        contract["source"]["published_date"] = "2026-09-09"
    elif mutation == "external_license_drift":
        contract["external_method_screening"][0]["license"] = "unknown"
    elif mutation == "external_scope_missing":
        contract["external_method_screening"][0]["candidate_scope"] = ""
    elif mutation == "external_reason_missing":
        contract["external_method_screening"][0]["reason"] = ""
    else:  # pragma: no cover - parametrization is closed above
        raise AssertionError(mutation)

    with pytest.raises(ValueError, match=message):
        validate_openai_public_observable_contract(contract)
