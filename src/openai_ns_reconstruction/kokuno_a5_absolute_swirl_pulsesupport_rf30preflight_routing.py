"""Kokuno Agent-5 fail-closed routing for the absolute-swirl/pulse-support seam.

This module is integration/provenance glue only.  It does not copy or alter
Agent-1--4 mathematics and it never transfers scoped evidence across candidate
identities.

Fresh bounded deliveries registered here:
* Agent 1 #1168 materializes an axis-anchored absolute current r_I at the
  eta-flattening endpoint together with the imported/base and total late
  relative-swirl angular targets.  It still does not compose the two bumps into
  Cartesian velocity, and GitHub reported the open PR non-mergeable at freeze.
* Agent 2 #1164 adds provider-certified source pulse-coordinate support to the
  bounded complete-curl family.  The v-map/support interval remain repository
  provider metadata; v is not identified with chart T or physical time t and
  the family is still not a self-contained project candidate.
* Agent 3 #1161 remains the current-I4 fixed-Q observable preflight and refuses
  to substitute physical azimuthal theta averaging for the missing source
  auxiliary-T^2 Haar mean required by RF30.
* Agent 4 #1160 remains an implementation-distinct audit of Agent-1 #1154's
  low-dimensional compensator algebra only.  It does not audit #1168's new
  absolute/total target and is not complete-NS evidence.

Frozen project gates are inherited unchanged from the previous Agent-5
registration: normalized momentum max/L2 <= 1e-3; divergence max/L2 <= 1e-5;
canonical quadrature [24,48,96]; residual-defined/free forcing forbidden.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from . import kokuno_a5_relative_swirl_support_rf30preflight_routing as parent_a5

SCHEMA = "kokuno-a5-absolute-swirl-pulsesupport-rf30preflight-routing-v1"
TASK = "KOKUNO-A5-ABSOLUTE-SWIRL-PULSESUPPORT-RF30PREFLIGHT-ROUTING-119"

PARENT_A5 = {
    "pr": 1162,
    "head": "ab8501a3f2e89e2f2e6f47c2ff67b53ead8eaa8f",
    "branch": "kokuno-agent5/relative-swirl-support-rf30preflight-routing-118",
    "source_path": "src/openai_ns_reconstruction/kokuno_a5_relative_swirl_support_rf30preflight_routing.py",
    "source_blob": "c9d0f0c9f8ac3007f11a9add94eb04afbbdfb1dd",
    "task": "KOKUNO-A5-RELATIVE-SWIRL-SUPPORT-RF30PREFLIGHT-ROUTING-118",
    "schema": "kokuno-a5-relative-swirl-support-rf30preflight-routing-v1",
    "dedicated_run": 35708528672,
    "repository_tests_run": 35708528631,
    "actions_status_at_registration": "queued",
}

AGENT1_ABSOLUTE_RELATIVE_SWIRL_TARGET = {
    "task": "KOKUNO-A1-ABSOLUTE-RELATIVE-SWIRL-TARGET",
    "pr": 1168,
    "head": "bb246140c6c2451b22b59aad9ba90604bc2acb21",
    "branch": "codex/kokuno-agent1-absolute-relative-swirl-target",
    "base_pr": 1159,
    "base_head": "6ddf35b9dc726c720033ee44161e29efc77d755e",
    "source_path": "src/openai_ns_reconstruction/kokuno_current_relative_swirl_absolute_target.py",
    "source_blob": "86b4804a77b2d17f1a881579dfbed7a17264fb49",
    "stage": "absolute-current-base-total-relative-swirl-angular-target",
    "exact_compare_ahead": 3,
    "exact_compare_behind": 0,
    "open_pr_mergeable_at_registration": False,
    "absolute_current_r_I_at_flattening_endpoint_materialized": True,
    "current_relative_swirl_total_target_materialized": True,
    "imported_base_absolute_angular_target_materialized": True,
    "eta_jet_materialized": True,
    "eta_jet_is_numerical_not_source_analytic": True,
    "current_cartesian_relative_swirl_composed": False,
    "terminal_global_leading_velocity_materialized": False,
    "matched_global_pressure_materialized": False,
    "restricted_forcing_materialized": False,
    "heldout_ns_residual_assessed": False,
    "scientifically_admitted": False,
    "pde_validated": False,
    "dedicated_run": 35709299481,
    "repository_tests_run": 35709299432,
    "actions_status_at_registration": "queued",
}

AGENT2_PULSE_COORDINATE_SUPPORT = {
    "task": "K2-OSC-107",
    "pr": 1164,
    "head": "0cbf74972af6589f7eb9a2edfe21c957d7512a5c",
    "branch": "kokuno-agent2-pulse-coordinate-support-107",
    "base_pr": 1158,
    "base_head": "f2b7ddb909b8fc2bcf575942ed63ec4a1fbc98bf",
    "source_path": "src/openai_ns_reconstruction/kokuno_source_pulse_coordinate_support_multiharmonic_velocity.py",
    "source_blob": "d95265098a9a9eb95ee673fecb3e7b115f2866e6",
    "stage": "provider-certified-pulse-coordinate-support-bounded-complete-curl-multiharmonic",
    "exact_compare_ahead": 3,
    "exact_compare_behind": 0,
    "open_pr_mergeable_at_registration": True,
    "provider_certified_pulse_coordinate_support_materialized": True,
    "pulse_coordinate_mapping_identity_bound": True,
    "pulse_coordinate_exact_zero_skips_harmonic_provider_evaluation": True,
    "support_applied_as_execution_skip_not_velocity_hard_mask": True,
    "source_exact_pulse_coordinate_mapping_recovered": False,
    "source_exact_pulse_support_interval_recovered": False,
    "source_pulse_coordinate_identified_with_chart_T": False,
    "source_pulse_coordinate_identified_with_physical_t": False,
    "physical_time_support_certificate_materialized": False,
    "source_provider_self_contained": False,
    "self_contained_velocity_xyzt_provider": False,
    "complete_ns_residual_assessed": False,
    "residual_reduction_claimed": False,
    "scientifically_admitted": False,
    "pde_validated": False,
    "dedicated_run": 35708875850,
    "repository_tests_run": 35708875701,
    "actions_status_at_registration": "queued",
}

_PARENT_REGISTRATION = parent_a5.build_registration()
PARENT_FRONTIERS = copy.deepcopy(_PARENT_REGISTRATION["frontiers"])
CORRECTED_SOURCE = copy.deepcopy(_PARENT_REGISTRATION["corrected_source"])
FROZEN_GATES = copy.deepcopy(_PARENT_REGISTRATION["frozen_gates"])
CORE_STATE = copy.deepcopy(_PARENT_REGISTRATION["core_state"])


def _truth_boundary() -> dict[str, Any]:
    truth = copy.deepcopy(_PARENT_REGISTRATION["truth_boundary"])
    # New A1 truths are target/cumulative truths only; Cartesian composition is
    # deliberately not inferred from them.
    truth.update(
        {
            "absolute_current_r_I_at_flattening_endpoint_materialized": True,
            "imported_base_absolute_angular_target_materialized": True,
            "current_relative_swirl_total_target_materialized": True,
            "current_cartesian_relative_swirl_composed": False,
            "terminal_global_leading_velocity_materialized": False,
            "agent1_1168_open_pr_mergeable_at_registration": False,
            "agent1_1168_scientifically_admitted": False,
            "agent4_1160_audits_agent1_1168_absolute_target": False,
            "provider_certified_pulse_coordinate_support_materialized": True,
            "source_exact_pulse_coordinate_mapping_recovered": False,
            "source_exact_pulse_support_interval_recovered": False,
            "source_pulse_coordinate_identified_with_chart_T": False,
            "source_pulse_coordinate_identified_with_physical_t": False,
            "physical_time_support_certificate_materialized": False,
            "pulse_support_source_provider_self_contained": False,
            "latest_self_contained_project_composite_is_agent2_1117": True,
            "current_i4_source_auxiliary_t2_provider_materialized": False,
            "current_i4_source_auxiliary_t2_haar_mean_materialized": False,
            "current_i4_rf30_repository_candidate_state_authorized": False,
            "current_i4_rf30_defect_materialized": False,
            "matched_cartesian_pressure_gradient_materialized": False,
            "preregistered_restricted_forcing_materialized": False,
            "restricted_forcing_proved_not_residual_defined": False,
            "complete_candidate_api_ready": False,
            "complete_identity_bound_ns_defect_materialized": False,
            "cartesian_correction_velocity_materialized": False,
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
    )
    return truth


def build_identity_firewall() -> dict[str, bool]:
    firewall = copy.deepcopy(_PARENT_REGISTRATION["identity_firewall"])
    firewall.update(
        {
            "agent1_1168_absolute_target_not_cartesian_relative_swirl_composition": True,
            "agent1_1168_nonmergeable_open_pr_not_scientific_admission": True,
            "agent4_1160_evidence_bound_to_agent1_1154_algebra_only": True,
            "agent4_1160_not_evidence_for_agent1_1168_absolute_target": True,
            "agent2_1164_pulse_coordinate_provider_metadata_not_source_exact_map": True,
            "agent2_1164_pulse_coordinate_not_chart_T_or_physical_t": True,
            "agent2_1164_provider_family_not_self_contained_project_candidate": True,
            "agent3_1161_physical_theta_mean_not_source_auxiliary_t2_haar_mean": True,
            "latest_self_contained_project_composite_stays_agent2_1117": True,
            "cross_identity_evidence_transfer_allowed": False,
        }
    )
    return firewall


SHORTEST_CLOSURE = [
    "compose exact Agent-1 #1168 total target through the existing two-bump solver into the same current Cartesian leading identity",
    "obtain implementation-distinct Agent-4 evidence for that exact composed relative-swirl Cartesian identity before transferring any scoped claim",
    "finish terminal/exterior/global leading and matched pressure, then build an exact-identity self-contained Agent-2 leading-plus-oscillatory project candidate",
    "materialize the source auxiliary-T2 RF30 Haar backend for the exact current correction identity; never substitute physical theta averaging",
    "bind preregistered non-residual-defined restricted forcing, complete NS defect, real Cartesian correction and finite cycle, then run Agent-4 held-out and canonical [24,48,96] admission",
]


def _canonical_bytes(value: Mapping[str, Any]) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")


def _without_digest(registration: Mapping[str, Any]) -> dict[str, Any]:
    payload = copy.deepcopy(dict(registration))
    payload.pop("registration_sha256", None)
    return payload


def registration_sha256(registration: Mapping[str, Any]) -> str:
    return hashlib.sha256(_canonical_bytes(_without_digest(registration))).hexdigest()


def build_registration() -> dict[str, Any]:
    frontiers = copy.deepcopy(PARENT_FRONTIERS)
    frontiers["absolute_relative_swirl_target"] = copy.deepcopy(
        AGENT1_ABSOLUTE_RELATIVE_SWIRL_TARGET
    )
    frontiers["pulse_coordinate_support_multiharmonic"] = copy.deepcopy(
        AGENT2_PULSE_COORDINATE_SUPPORT
    )
    registration: dict[str, Any] = {
        "schema": SCHEMA,
        "task": TASK,
        "role": "Agent-5 integration/provenance/artifact/CI routing only",
        "parent_a5": copy.deepcopy(PARENT_A5),
        "corrected_source": copy.deepcopy(CORRECTED_SOURCE),
        "frontiers": frontiers,
        "identity_firewall": build_identity_firewall(),
        "truth_boundary": _truth_boundary(),
        "core_state": copy.deepcopy(CORE_STATE),
        "frozen_gates": copy.deepcopy(FROZEN_GATES),
        "shortest_closure": list(SHORTEST_CLOSURE),
    }
    registration["registration_sha256"] = registration_sha256(registration)
    return registration


def validate_registration(registration: Mapping[str, Any]) -> None:
    reg = copy.deepcopy(dict(registration))
    if reg.get("schema") != SCHEMA or reg.get("task") != TASK:
        raise ValueError("Agent-5 routing schema/task drifted")
    if reg.get("registration_sha256") != registration_sha256(reg):
        raise ValueError("Agent-5 registration digest mismatch")
    if reg.get("parent_a5") != PARENT_A5:
        raise ValueError("Agent-5 parent identity drifted")
    if reg.get("corrected_source") != CORRECTED_SOURCE:
        raise ValueError("corrected-source provenance drifted")
    if reg.get("core_state") != CORE_STATE:
        raise ValueError("core readiness state was promoted or changed")
    if reg.get("frozen_gates") != FROZEN_GATES:
        raise ValueError("frozen scientific gates changed")
    if reg.get("identity_firewall") != build_identity_firewall():
        raise ValueError("identity firewall drifted")
    if reg.get("truth_boundary") != _truth_boundary():
        raise ValueError("truth boundary drifted")
    if reg.get("shortest_closure") != SHORTEST_CLOSURE:
        raise ValueError("shortest-closure routing drifted")

    frontiers = reg.get("frontiers")
    if not isinstance(frontiers, dict):
        raise ValueError("frontiers must be a mapping")
    expected_parent_keys = set(PARENT_FRONTIERS)
    if not expected_parent_keys.issubset(frontiers):
        raise ValueError("parent Agent-5 frontier was dropped")
    for key in expected_parent_keys:
        if frontiers[key] != PARENT_FRONTIERS[key]:
            raise ValueError(f"inherited frontier drifted: {key}")
    if frontiers.get("absolute_relative_swirl_target") != AGENT1_ABSOLUTE_RELATIVE_SWIRL_TARGET:
        raise ValueError("Agent-1 #1168 absolute-target identity drifted")
    if frontiers.get("pulse_coordinate_support_multiharmonic") != AGENT2_PULSE_COORDINATE_SUPPORT:
        raise ValueError("Agent-2 #1164 pulse-support identity drifted")

    # Explicit fail-closed semantic checks beyond exact dictionary equality.
    a1 = frontiers["absolute_relative_swirl_target"]
    if a1["open_pr_mergeable_at_registration"] is not False:
        raise ValueError("Agent-1 #1168 mergeability truth was laundered")
    if a1["current_cartesian_relative_swirl_composed"] is not False:
        raise ValueError("absolute angular target was laundered into Cartesian composition")
    if a1["scientifically_admitted"] is not False:
        raise ValueError("open sibling delivery was laundered into scientific admission")

    a2 = frontiers["pulse_coordinate_support_multiharmonic"]
    forbidden_true = (
        "source_exact_pulse_coordinate_mapping_recovered",
        "source_exact_pulse_support_interval_recovered",
        "source_pulse_coordinate_identified_with_chart_T",
        "source_pulse_coordinate_identified_with_physical_t",
        "physical_time_support_certificate_materialized",
        "source_provider_self_contained",
        "self_contained_velocity_xyzt_provider",
        "complete_ns_residual_assessed",
        "residual_reduction_claimed",
        "scientifically_admitted",
        "pde_validated",
    )
    if any(a2[name] is not False for name in forbidden_true):
        raise ValueError("pulse-coordinate provider boundary was promoted")

    if frontiers["relative_swirl_compensator_validator"]["audited_agent1_pr"] != 1154:
        raise ValueError("Agent-4 #1160 evidence transferred across Agent-1 identities")
    if frontiers["current_i4_rf30_fixedq_preflight"]["source_auxiliary_t2_provider_available"] is not False:
        raise ValueError("missing source auxiliary-T2 provider was fabricated")
    if frontiers["current_i4_rf30_fixedq_preflight"]["rf30_repository_candidate_state_authorized"] is not False:
        raise ValueError("RF30 candidate state was authorized without source Haar backend")
    if frontiers["latest_self_contained_project_composite"]["pr"] != 1117:
        raise ValueError("provider-driven pulse family was laundered into self-contained candidate")

    required_false_truth = (
        "current_cartesian_relative_swirl_composed",
        "terminal_global_leading_velocity_materialized",
        "agent4_1160_audits_agent1_1168_absolute_target",
        "physical_time_support_certificate_materialized",
        "pulse_support_source_provider_self_contained",
        "current_i4_source_auxiliary_t2_provider_materialized",
        "current_i4_source_auxiliary_t2_haar_mean_materialized",
        "current_i4_rf30_repository_candidate_state_authorized",
        "current_i4_rf30_defect_materialized",
        "matched_cartesian_pressure_gradient_materialized",
        "preregistered_restricted_forcing_materialized",
        "complete_identity_bound_ns_defect_materialized",
        "cartesian_correction_velocity_materialized",
        "finite_correction_cycle_run",
        "heldout_complete_ns_residual_assessed",
        "canonical_whole_domain_admission_run",
        "same_protocol_comparable_to_st006",
        "scientific_admission",
        "pde_validated",
    )
    truth = reg["truth_boundary"]
    if any(truth[name] is not False for name in required_false_truth):
        raise ValueError("final or scoped truth boundary was promoted")


def write_registration(path: str | Path) -> dict[str, Any]:
    registration = build_registration()
    validate_registration(registration)
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(registration, indent=2, sort_keys=True) + "\n")
    return registration


def load_registration(path: str | Path) -> dict[str, Any]:
    payload = json.loads(Path(path).read_text())
    if not isinstance(payload, dict):
        raise ValueError("registration must decode to a JSON object")
    validate_registration(payload)
    return payload


def _main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, help="Path for deterministic routing JSON")
    args = parser.parse_args(argv)
    write_registration(args.output)


if __name__ == "__main__":
    _main()
