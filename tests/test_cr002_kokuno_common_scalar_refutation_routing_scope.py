from __future__ import annotations

import inspect
import json
import math
from pathlib import Path

import pytest

from openai_ns_reconstruction.audit_kokuno_common_scalar_refutation_routing_scope import (
    audit,
    route_refutation_report,
)
from openai_ns_reconstruction.kokuno_current_bridge_common_scalar_refutation import (
    KokunoCurrentBridgeCommonScalarRefutation,
)


def _truth_boundary() -> dict[str, bool]:
    return {
        "current_cartesian_terminal_multiplier_composed": False,
        "outer_global_leading_velocity_materialized": False,
        "unified_global_cartesian_velocity_export_ready": False,
        "heldout_ns_residual_assessed": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
    }


def _common_false() -> dict[str, object]:
    return {
        "finite_witness_is_sufficient_only_for_refutation": True,
        "absence_of_refutation_does_not_establish_common_scalar_bridge": True,
        "full_eta_interval_entry_above_qp_established": False,
        "full_eta_interval_target_totality_established": False,
        "full_eta_interval_root_uniqueness_established": False,
        "full_eta_interval_transversality_established": False,
        "smooth_eta_target_time_map_established": False,
        "full_eta_interval_common_scalar_bridge_length_established": False,
        "current_l_minus_h_matching_bridge_materialized": False,
        "eta_dependent_cartesian_matching_boundary_materialized": False,
        "truth_boundary": _truth_boundary(),
    }


def _positive_report() -> dict[str, object]:
    report = _common_false()
    report.update(
        {
            "common_scalar_bridge_refuted_by_disjoint_unique_root_brackets": True,
            "disjoint_root_pair": {
                "earlier_root_bracket": [0.3, 0.31],
                "later_root_bracket": [0.6, 0.61],
                "disjoint_gap_lower_bound": 0.29,
            },
            "disjoint_gap_lower_bound": 0.29,
        }
    )
    return report


def _negative_report() -> dict[str, object]:
    report = _common_false()
    report.update(
        {
            "common_scalar_bridge_refuted_by_disjoint_unique_root_brackets": False,
            "disjoint_root_pair": None,
            "disjoint_gap_lower_bound": 0.0,
        }
    )
    return report


def test_contract_audit_passes_and_preserves_canonical_delivery_independence() -> None:
    report = audit()
    assert report["dependency_pr"] == 1242
    assert report["dependency_exact_head"] == "7e987ff0a3bea743fe3f0fe9d02d5b3c8c593ea3"
    assert report["dependency_module_blob_sha1"] == "409031a36425befefa7038ebda8ca848cdde239a"
    truth = report["truth_state"]
    assert truth["current_common_scalar_bridge_witness_check_materialized"] is True
    assert truth["common_scalar_bridge_established"] is False
    assert truth["current_l_minus_h_matching_bridge_materialized"] is False
    assert truth["eta_dependent_cartesian_matching_boundary_materialized"] is False
    assert truth["kokuno_unified_global_velocity_export_ready"] is False
    assert truth["pde_validated"] is False
    canonical = report["canonical_states"]
    assert canonical["velocity_export_ready"] is True
    assert canonical["visual_correspondence_verified"] is False
    assert canonical["pde_validated"] is False
    assert canonical["paper_exact"] is False
    assert canonical["openai_field_identified"] is False


def test_positive_refutation_routes_to_new_identity_without_scientific_promotion() -> None:
    routed = route_refutation_report(_positive_report())
    assert routed["dependency_positive_refutation"] is True
    assert routed["common_scalar_bridge_status"] == "refuted_for_exact_current_candidate_identity"
    assert routed["current_common_scalar_bridge_established"] is False
    assert routed["current_l_minus_h_matching_bridge_materialized"] is False
    assert routed["eta_dependent_cartesian_matching_boundary_materialized"] is False
    assert routed["current_cartesian_terminal_multiplier_composed"] is False
    assert routed["kokuno_unified_global_velocity_export_ready"] is False
    assert routed["pde_validated"] is False
    assert routed["paper_exact"] is False
    assert routed["openai_field_identified"] is False
    assert "new_upstream_candidate_identity_then_retest_scalar_compatibility" in routed[
        "allowed_next_representation_routes"
    ]
    assert "new_eta_dependent_autonomous_bridge_identity_after_continuum_gates" in routed[
        "allowed_next_representation_routes"
    ]


def test_negative_refutation_outcome_is_inconclusive_not_scalar_bridge_evidence() -> None:
    routed = route_refutation_report(_negative_report())
    assert routed["dependency_positive_refutation"] is False
    assert routed["common_scalar_bridge_status"] == "unresolved_not_established"
    assert routed["current_common_scalar_bridge_established"] is False
    assert routed["current_l_minus_h_matching_bridge_materialized"] is False
    assert "continuum_common_scalar_proof_before_promotion" in routed[
        "allowed_next_representation_routes"
    ]


def test_malformed_positive_refutation_fails_closed() -> None:
    report = _positive_report()
    report["disjoint_root_pair"] = None
    with pytest.raises(AssertionError, match="requires a witness pair"):
        route_refutation_report(report)

    report = _positive_report()
    report["disjoint_root_pair"]["later_root_bracket"] = [0.30, 0.305]
    with pytest.raises(AssertionError, match="ordered disjoint brackets"):
        route_refutation_report(report)


def test_dependency_truth_promotion_fails_closed() -> None:
    report = _positive_report()
    report["truth_boundary"]["pde_validated"] = True
    with pytest.raises(AssertionError, match="pde_validated"):
        route_refutation_report(report)

    report = _negative_report()
    report["full_eta_interval_target_totality_established"] = True
    with pytest.raises(AssertionError, match="target_totality"):
        route_refutation_report(report)


def test_exact_current_dependency_runtime_maps_through_either_legal_outcome() -> None:
    dependency = KokunoCurrentBridgeCommonScalarRefutation()
    source_report = dependency.refutation_report()
    routed = route_refutation_report(source_report)
    assert routed["dependency_positive_refutation"] is source_report[
        "common_scalar_bridge_refuted_by_disjoint_unique_root_brackets"
    ]
    assert routed["current_common_scalar_bridge_established"] is False
    assert routed["current_l_minus_h_matching_bridge_materialized"] is False
    assert routed["eta_dependent_cartesian_matching_boundary_materialized"] is False
    assert routed["pde_validated"] is False


def test_autonomous_mechanics_witness_has_distinct_unique_target_times() -> None:
    report = audit()
    witness = report["mechanics_witness"]
    assert witness["role"] == "autonomous_representation_mechanics_only"
    assert witness["target"] == 1.0
    assert witness["root_curve_a"] == pytest.approx(math.log(2.0))
    assert witness["root_curve_b"] == pytest.approx(0.5 * math.log(2.0))
    assert witness["absolute_root_gap"] > 0.0
    assert witness["common_scalar_time_exists_for_both_curves"] is False


def test_contract_separates_four_provenance_classes_and_forbids_transfer() -> None:
    root = Path(__file__).resolve().parents[1]
    contract = json.loads(
        (root / "configs/kokuno_common_scalar_refutation_routing_scope.json").read_text(
            encoding="utf-8"
        )
    )
    assert set(contract["source_classification"]) == {
        "user_requirement",
        "public_source_fact",
        "autonomous_design",
        "pending_unknown",
    }
    assert all(contract["source_classification"][key] for key in contract["source_classification"])
    assert all(value is False for value in contract["evidence_transfer"].values())
    identity = contract["identity_policy"]
    assert identity["upstream_current_state_revision_requires_new_candidate_identity"] is True
    assert identity["eta_dependent_bridge_requires_new_representation_identity"] is True
    assert identity["prior_refutation_receipt_transfers_to_revised_candidate"] is False
    assert identity["source_scalar_stage_label_transfers_to_eta_dependent_autonomous_design"] is False


def test_route_api_exposes_no_scientific_or_bridge_selection_tuning_knobs() -> None:
    sig = inspect.signature(route_refutation_report)
    assert tuple(sig.parameters) == ("report",)
    forbidden = {
        "residual",
        "forcing",
        "pressure",
        "optimizer",
        "threshold",
        "tolerance",
        "gain",
        "h",
        "bridge_length",
        "representative_eta",
        "witness_count",
    }
    assert forbidden.isdisjoint(sig.parameters)
