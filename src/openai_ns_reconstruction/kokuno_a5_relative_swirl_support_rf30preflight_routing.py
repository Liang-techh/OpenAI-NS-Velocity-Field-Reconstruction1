"""Kokuno Agent-5 fail-closed routing for the newest A1/A2 seams.

This module is integration/provenance glue only. It does not alter Agent-1--4
mathematics or transfer scoped evidence across candidate identities.

Fresh facts bound here:

* A1 #1159 materializes only the current-minus-imported contribution to the
  late relative-swirl angular target. The imported/base absolute target is not
  materialized, so the total target and Cartesian relative-swirl composition
  remain unavailable. The newest bounded Cartesian leading therefore remains
  exact A1 #1148.
* A2 #1158 adds identity-bound conservative outer-support execution
  certificates to the bounded complete-curl multi-harmonic family. The source
  input provider remains external; this is not a self-contained project field.
* A3's newest delivered frontier remains #1150. Central coordination has a
  newer claimed RF30 fixed-Q preflight, but a claim is not a delivery and is
  recorded only as pending/not evidence.
* A4 #1152 remains an implementation-distinct, leading-only divergence audit
  of exact A1 #1148. It does not audit #1159 and is not complete-NS evidence.

Project gates remain frozen: normalized momentum sampled-max and volume-L2
<=1e-3; normalized divergence sampled-max and volume-L2 <=1e-5; canonical
quadrature [24,48,96]; residual-defined/free forcing forbidden.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from . import kokuno_a5_holdprefix_multiharmonic_rf49_audit_routing as parent_a5

SCHEMA = "kokuno-a5-relative-swirl-support-rf30preflight-routing-v1"
TASK = "KOKUNO-A5-RELATIVE-SWIRL-SUPPORT-RF30PREFLIGHT-ROUTING-118"

PARENT_A5 = {
    "pr": 1153,
    "head": "eb15089a312b9ea1a3c4632145d501d244ec1324",
    "branch": "kokuno-agent5/holdprefix-multiharmonic-rf49-audit-routing-117",
    "source_path": "src/openai_ns_reconstruction/kokuno_a5_holdprefix_multiharmonic_rf49_audit_routing.py",
    "task": "KOKUNO-A5-HOLDPREFIX-MULTIHARMONIC-RF49-AUDIT-ROUTING-117",
    "schema": "kokuno-a5-holdprefix-multiharmonic-rf49-audit-routing-v1",
    "dedicated_run": 35702318127,
    "repository_tests_run": 35702318118,
    "actions_status_at_registration": "queued",
}

AGENT1_CURRENT_ANGULAR_CORRECTION = {
    "task": "KOKUNO-A1-CURRENT-ANGULAR-DEFECT-120",
    "pr": 1159,
    "head": "6ddf35b9dc726c720033ee44161e29efc77d755e",
    "branch": "codex/kokuno-a1-current-angular-defect",
    "source_path": "src/openai_ns_reconstruction/kokuno_current_relative_swirl_angular_correction.py",
    "source_blob": "37e6b0814fbad14084a8c8c4815122111252ac9b",
    "parent_pr": 1154,
    "parent_head": "d0e855e37c8f5d6d58f1a85b1dcb874651d20838",
    "stage": "current-minus-imported-relative-swirl-angular-target-correction-only",
    "current_lineage_angular_target_correction_materialized": True,
    "imported_base_absolute_angular_target_materialized": False,
    "current_lineage_angular_entry_I_materialized": False,
    "current_relative_swirl_total_target_materialized": False,
    "current_cartesian_relative_swirl_composed": False,
    "terminal_global_leading_velocity_materialized": False,
    "complete_ns_residual_assessed": False,
    "pde_validated": False,
    "dedicated_run": 35706739448,
    "repository_tests_run": 35706739342,
    "actions_status_at_registration": "queued",
}

AGENT2_SUPPORT_CERTIFIED_MULTIHARMONIC = {
    "task": "K2-OSC-106",
    "pr": 1158,
    "head": "f2b7ddb909b8fc2bcf575942ed63ec4a1fbc98bf",
    "branch": "kokuno-agent2-support-certified-multiharmonic-106",
    "source_path": "src/openai_ns_reconstruction/kokuno_source_support_certified_multiharmonic_velocity.py",
    "source_blob": "5b2bdab4e7fff24ecec8f8cb0fbd8d118a862381",
    "parent_pr": 1149,
    "parent_head": "db7a9a4efb4faed2fa46d321fdf5a66b8686de1f",
    "stage": "provider-driven-support-certified-bounded-complete-curl-multiharmonic",
    "provider_certified_outer_spatial_support_materialized": True,
    "outer_support_exact_zero_skips_provider_evaluation": True,
    "support_applied_as_execution_skip_not_velocity_hard_mask": True,
    "source_exact_outer_support_geometry_recovered": False,
    "pulse_time_support_certificate_materialized": False,
    "source_provider_self_contained": False,
    "self_contained_velocity_xyzt_provider": False,
    "complete_ns_residual_assessed": False,
    "residual_reduction_claimed": False,
    "pde_validated": False,
    "dedicated_run": 35707176374,
    "repository_tests_run": 35707176333,
    "actions_status_at_registration": "queued",
}

AGENT3_DELIVERED_RF49 = copy.deepcopy(parent_a5.AGENT3_RF49_GAIN_FIREWALL)
AGENT3_PENDING_RF30_PREFLIGHT_CLAIM = {
    "task": "KOKUNO-A3-CURRENT-I4-RF30-FIXEDQ-PREFLIGHT-126",
    "central_issue": 15,
    "claim_comment_id": 5773639838,
    "exact_base_head": "e0b48ccfa964db0caa9d657961d451a21e1d53bf",
    "state": "claimed_not_delivered_at_registration",
    "pr": None,
    "head": None,
    "materialized": False,
    "candidate_evidence": False,
    "may_promote_truth_state": False,
}

LATEST_CARTESIAN_LEADING = copy.deepcopy(parent_a5.AGENT1_HOLD_PREFIX)
AGENT4_HOLD_PREFIX_AUDIT = copy.deepcopy(parent_a5.AGENT4_HOLD_PREFIX_LEADING_AUDIT)
LATEST_SELF_CONTAINED_PROJECT_COMPOSITE = copy.deepcopy(parent_a5.LATEST_SELF_CONTAINED_PROJECT_COMPOSITE)
AGENT4_XI11_COMPOSITE_AUDIT = copy.deepcopy(parent_a5.AGENT4_XI11_COMPOSITE_DIVERGENCE_AUDIT)
CURRENT_I4_REAL_RADIAL_FORCE = copy.deepcopy(parent_a5.CURRENT_I4_REAL_RADIAL_FORCE)
AGENT4_CURRENT_I4_RADIAL_FORCE_AUDIT = copy.deepcopy(parent_a5.AGENT4_CURRENT_I4_RADIAL_FORCE_AUDIT)
CORRECTED_SOURCE = copy.deepcopy(parent_a5.CORRECTED_SOURCE)
FROZEN_GATES = copy.deepcopy(parent_a5.FROZEN_GATES)
CORE_STATE = copy.deepcopy(parent_a5.CORE_STATE)


def _truth_boundary() -> dict[str, Any]:
    return {
        "latest_cartesian_leading_is_agent1_1148": True,
        "current_relative_swirl_angular_correction_materialized": True,
        "imported_base_absolute_angular_target_materialized": False,
        "current_lineage_angular_entry_I_materialized": False,
        "current_relative_swirl_total_target_materialized": False,
        "current_cartesian_relative_swirl_composed": False,
        "terminal_global_leading_velocity_materialized": False,
        "matching_latest_cartesian_leading_agent4_audit_present": True,
        "matching_latest_cartesian_leading_agent4_audit_scoped_gate_passed": None,
        "agent4_1152_audits_agent1_1159": False,
        "support_certified_bounded_multiharmonic_materialized": True,
        "support_certified_multiharmonic_source_provider_self_contained": False,
        "support_certified_multiharmonic_complete_ns_residual_assessed": False,
        "source_exact_outer_support_geometry_recovered": False,
        "pulse_time_support_certificate_materialized": False,
        "latest_self_contained_project_composite_is_agent2_1117": True,
        "agent4_matching_xi11_composite_audit_present": True,
        "agent4_matching_xi11_composite_audit_scoped_gate_passed": None,
        "rf49_source_exponent_gain_firewall_materialized": True,
        "rf49_gains_are_raw_residual_contraction_factors": False,
        "rf30_fixedq_preflight_claim_present": True,
        "rf30_fixedq_preflight_delivered": False,
        "current_i4_source_chart_backend_materialized": False,
        "current_i4_rf30_defect_materialized": False,
        "current_i4_source_specific_gain_evidence": False,
        "current_i4_radial_force_materialized": True,
        "agent4_matching_current_i4_radial_force_audit_present": True,
        "cartesian_correction_velocity_materialized": False,
        "matched_cartesian_pressure_gradient_materialized": False,
        "preregistered_restricted_forcing_materialized": False,
        "restricted_forcing_proved_not_residual_defined": False,
        "complete_candidate_api_ready": False,
        "complete_identity_bound_ns_defect_materialized": False,
        "finite_correction_cycle_run": False,
        "heldout_complete_ns_residual_assessed": False,
        "canonical_whole_domain_admission_run": False,
        "same_protocol_comparable_to_st006": False,
        "scientific_admission": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
        "pde_validated": False,
    }


def build_identity_firewall() -> dict[str, bool]:
    return {
        "agent1_1159_correction_not_absolute_angular_target": True,
        "agent1_1159_not_cartesian_relative_swirl_composition": True,
        "agent4_1152_evidence_bound_to_agent1_1148_only": True,
        "agent4_1152_not_evidence_for_agent1_1159": True,
        "agent2_1158_support_metadata_not_source_exact_support": True,
        "agent2_1158_provider_family_not_self_contained_project_candidate": True,
        "agent3_pending_claim_not_delivery_or_candidate_evidence": True,
        "agent3_rf49_exponent_gain_not_raw_residual_contraction": True,
        "latest_self_contained_project_composite_stays_agent2_1117": True,
        "cross_identity_evidence_transfer_allowed": False,
    }


def _canonical_bytes(value: Mapping[str, Any]) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")


def _without_digest(registration: Mapping[str, Any]) -> dict[str, Any]:
    payload = copy.deepcopy(dict(registration))
    payload.pop("registration_sha256", None)
    return payload


def registration_sha256(registration: Mapping[str, Any]) -> str:
    return hashlib.sha256(_canonical_bytes(_without_digest(registration))).hexdigest()


def build_registration() -> dict[str, Any]:
    registration: dict[str, Any] = {
        "schema": SCHEMA,
        "task": TASK,
        "role": "Agent-5 integration/provenance/artifact/CI routing only",
        "parent_a5": copy.deepcopy(PARENT_A5),
        "corrected_source": copy.deepcopy(CORRECTED_SOURCE),
        "frontiers": {
            "latest_cartesian_leading": copy.deepcopy(LATEST_CARTESIAN_LEADING),
            "current_angular_target_correction": copy.deepcopy(AGENT1_CURRENT_ANGULAR_CORRECTION),
            "support_certified_multiharmonic": copy.deepcopy(AGENT2_SUPPORT_CERTIFIED_MULTIHARMONIC),
            "latest_self_contained_project_composite": copy.deepcopy(LATEST_SELF_CONTAINED_PROJECT_COMPOSITE),
            "xi11_composite_validator": copy.deepcopy(AGENT4_XI11_COMPOSITE_AUDIT),
            "rf49_delivered": copy.deepcopy(AGENT3_DELIVERED_RF49),
            "rf30_fixedq_preflight_pending_claim": copy.deepcopy(AGENT3_PENDING_RF30_PREFLIGHT_CLAIM),
            "holdprefix_leading_validator": copy.deepcopy(AGENT4_HOLD_PREFIX_AUDIT),
            "current_i4_real_radial_force": copy.deepcopy(CURRENT_I4_REAL_RADIAL_FORCE),
            "current_i4_radial_force_validator": copy.deepcopy(AGENT4_CURRENT_I4_RADIAL_FORCE_AUDIT),
        },
        "identity_firewall": build_identity_firewall(),
        "truth_boundary": _truth_boundary(),
        "core_state": copy.deepcopy(CORE_STATE),
        "frozen_gates": copy.deepcopy(FROZEN_GATES),
        "shortest_closure": [
            "materialize the missing imported/base absolute angular target before adding A1 #1159 current correction into a total relative-swirl target",
            "compose the resulting two relative-swirl bumps into the exact A1 current Cartesian identity, then finish terminal/exterior/global leading",
            "make A2 consume that exact global identity into one self-contained leading-plus-oscillatory candidate; provider support certificates alone do not close this seam",
            "deliver the claimed candidate-bound fixed-Q RF30 preflight and then the real source auxiliary-T2/Haar covariance backend before applying RF34-RF49 mechanics",
            "bind matched pressure and preregistered non-residual-defined restricted forcing before constructing a complete NS defect",
            "only then run the real finite correction cycle and Agent-4 held-out whole-domain 1e-3 momentum gate",
        ],
    }
    registration["registration_sha256"] = registration_sha256(registration)
    validate_registration(registration)
    return registration


def validate_registration(registration: Mapping[str, Any]) -> None:
    if registration.get("schema") != SCHEMA or registration.get("task") != TASK:
        raise ValueError("unexpected Agent-5 routing schema/task")
    if registration.get("registration_sha256") != registration_sha256(registration):
        raise ValueError("routing registration digest mismatch")

    frontiers = registration["frontiers"]
    a1lead = frontiers["latest_cartesian_leading"]
    a1corr = frontiers["current_angular_target_correction"]
    a2 = frontiers["support_certified_multiharmonic"]
    a3pending = frontiers["rf30_fixedq_preflight_pending_claim"]
    a4 = frontiers["holdprefix_leading_validator"]
    truth = registration["truth_boundary"]
    firewall = registration["identity_firewall"]

    if a1lead["pr"] != 1148 or a1corr["pr"] != 1159:
        raise ValueError("A1 frontier identity drift")
    if a1corr["imported_base_absolute_angular_target_materialized"]:
        raise ValueError("A1 current correction laundered into absolute angular target")
    if a1corr["current_cartesian_relative_swirl_composed"]:
        raise ValueError("A1 algebra/correction laundered into Cartesian composition")
    if a2["pr"] != 1158 or a2["source_provider_self_contained"]:
        raise ValueError("A2 provider-driven support family misclassified")
    if a2["source_exact_outer_support_geometry_recovered"]:
        raise ValueError("A2 provider support metadata laundered into source truth")
    if a3pending["materialized"] or a3pending["candidate_evidence"] or a3pending["pr"] is not None:
        raise ValueError("pending A3 claim laundered into delivery/evidence")
    if a4["pr"] != 1152 or a4["audited_agent1_pr"] != 1148:
        raise ValueError("A4 hold-prefix audit identity drift")
    if a4["scoped_gate_passed"] is not None or a4["scientifically_admitted"]:
        raise ValueError("unresolved A4 audit promoted to PASS/admission")

    required_false = (
        "imported_base_absolute_angular_target_materialized",
        "current_relative_swirl_total_target_materialized",
        "current_cartesian_relative_swirl_composed",
        "terminal_global_leading_velocity_materialized",
        "rf30_fixedq_preflight_delivered",
        "current_i4_source_chart_backend_materialized",
        "current_i4_rf30_defect_materialized",
        "cartesian_correction_velocity_materialized",
        "matched_cartesian_pressure_gradient_materialized",
        "preregistered_restricted_forcing_materialized",
        "complete_candidate_api_ready",
        "complete_identity_bound_ns_defect_materialized",
        "finite_correction_cycle_run",
        "heldout_complete_ns_residual_assessed",
        "same_protocol_comparable_to_st006",
        "scientific_admission",
        "pde_validated",
    )
    if any(truth[key] for key in required_false):
        raise ValueError("fail-closed project truth boundary promoted")
    if firewall["cross_identity_evidence_transfer_allowed"]:
        raise ValueError("cross-identity evidence transfer enabled")

    expected_state = {
        "leading_ready": False,
        "oscillatory_ready": True,
        "correction_ready": False,
        "velocity_export_ready": False,
        "pde_validated": False,
    }
    if registration["core_state"] != expected_state:
        raise ValueError("core state drift")
    if registration["frozen_gates"] != FROZEN_GATES:
        raise ValueError("frozen project gates changed")
    if FROZEN_GATES["normalized_momentum_sampled_max"] != 1.0e-3:
        raise ValueError("momentum max gate drift")
    if FROZEN_GATES["normalized_momentum_volume_l2"] != 1.0e-3:
        raise ValueError("momentum L2 gate drift")
    if FROZEN_GATES["normalized_divergence_sampled_max"] != 1.0e-5:
        raise ValueError("divergence max gate drift")
    if FROZEN_GATES["normalized_divergence_volume_l2"] != 1.0e-5:
        raise ValueError("divergence L2 gate drift")
    if FROZEN_GATES["canonical_quadrature"] != [24, 48, 96]:
        raise ValueError("canonical quadrature drift")
    if FROZEN_GATES["residual_defined_free_forcing_allowed"]:
        raise ValueError("residual-defined/free forcing unexpectedly allowed")


def write_registration(path: str | Path) -> dict[str, Any]:
    registration = build_registration()
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(registration, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return registration


def load_registration(path: str | Path) -> dict[str, Any]:
    registration = json.loads(Path(path).read_text(encoding="utf-8"))
    validate_registration(registration)
    return registration


def _main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    registration = write_registration(args.output)
    print(json.dumps(registration, indent=2, sort_keys=True))


if __name__ == "__main__":
    _main()
