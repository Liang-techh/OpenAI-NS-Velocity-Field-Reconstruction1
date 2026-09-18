import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def _load(path: str) -> dict:
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def test_axial_cap_routing_scope_matches_live_project_route():
    status = _load("project_status.json")
    scope = _load("configs/bipolar_axial_cap_poloidal_representation_scope.json")
    routing = scope["integration_routing_scope"]

    assert routing["scope_kind"] == (
        "active_child_materialization_route_with_prematerialization_evidence"
    )
    assert routing["live_integrated_mode"] == "AXIAL_CAP_BANDED_C4_ODD_Z_POLOIDAL"
    assert routing["live_route"] == status["active_scientific_route"]
    assert routing["live_mode_nonzero_child_materialized"] is False
    assert routing["axial_cap_capacity_is_integrated_prematerialization_evidence"] is True
    assert routing["additional_capacity_growth_allowed_before_live_materialization_screen"] is False
    assert routing["route_guard"] == (
        "axial_cap_prematerialization_evidence_does_not_itself_constitute_"
        "production_materialization_promotion_or_pde_validation"
    )
    assert routing["promotion_requires"] == [
        "explicit_autonomous_coefficient_bound",
        "explicit_nonzero_coefficient_value",
        "new_representation_family_identity",
        "new_candidate_sha256",
        "fresh_full_candidate_validation",
    ]

    assert status["next_integration_task"] == (
        "materialize_bounded_nonzero_AXIAL_CAP_BANDED_C4_ODD_Z_POLOIDAL_child_then_validate"
    )
