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


def test_public_observable_contract_is_qualitative_and_fail_closed() -> None:
    contract = load_openai_public_observable_contract()
    assert [item["id"] for item in contract["observables"]] == EXPECTED_IDS
    assert all(item["numerical_target"] is None for item in contract["observables"])
    assert contract["external_method_screening"][0]["classification"] == "screened_not_adopted"
    assert contract["external_method_screening"][0]["migration_scope"] == "none"
    assert all(value is False for value in contract["truth_boundary"].values())

    report = observable_contract_report()
    assert report["observable_ids"] == EXPECTED_IDS
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
    else:  # pragma: no cover - parametrization is closed above
        raise AssertionError(mutation)

    with pytest.raises(ValueError, match=message):
        validate_openai_public_observable_contract(contract)
