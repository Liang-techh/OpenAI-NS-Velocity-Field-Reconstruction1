"""Kokuno Agent-5 fail-closed routing for split-delta/orientation integration.

Integration/provenance glue only. No Agent-1--4 mathematics is copied or
changed and scoped evidence is never transferred across semantic identities.

Registered deliveries:
* A1 #1171 preserves the mathematically nonzero late relative-swirl correction
  as representable base+delta channels. The correction is below ordinary
  float64 relative resolution, so a unified ordinary-float64 Cartesian field
  is deliberately not claimed composed.
* A2 #1170 adds one bounded repository-autonomous rigid-z azimuth to the
  provider-driven complete-curl family. Source-exact frame orientation is not
  recovered and the provider remains non-self-contained.
* A3 #1172 makes the corrected-source normalized auxiliary-T^2 Haar covariance
  operator executable from raw wave samples, but deliberately leaves the
  repository provider blob unpinned. Therefore RF30 repository candidate state
  and defect remain unauthorized/unmaterialized.
* A4 #1174 independently audits exact A1 #1171's save/reloaded public split
  Cartesian evaluator. It is scoped split/divergence evidence only; its exact
  Actions are unresolved at registration and it is not complete-NS evidence.

Frozen project gates remain unchanged: normalized momentum max/L2 <=1e-3,
divergence max/L2 <=1e-5, canonical quadrature [24,48,96], and no
residual-defined/free forcing.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from . import kokuno_a5_absolute_swirl_pulsesupport_rf30preflight_routing as parent_a5

SCHEMA = "kokuno-a5-splitdelta-orientation-rf30preflight-routing-v1"
TASK = "KOKUNO-A5-SPLITDELTA-ORIENTATION-RF30PREFLIGHT-ROUTING-120"

PARENT_A5 = {
    "pr": 1169,
    "head": "a035784fcba6d373f5dffc062872b632c6a26309",
    "branch": "kokuno-agent5/absolute-swirl-pulsesupport-rf30preflight-routing-119",
    "source_path": "src/openai_ns_reconstruction/kokuno_a5_absolute_swirl_pulsesupport_rf30preflight_routing.py",
    "source_blob": "83ff6dfb3092e20382e05b6fe689491512e9eff6",
    "task": "KOKUNO-A5-ABSOLUTE-SWIRL-PULSESUPPORT-RF30PREFLIGHT-ROUTING-119",
    "schema": "kokuno-a5-absolute-swirl-pulsesupport-rf30preflight-routing-v1",
    "dedicated_run": 35710456661,
    "repository_tests_run": 35710456445,
    "actions_status_at_registration": "queued",
}

AGENT1_RELATIVE_SWIRL_SPLIT = {
    "task": "KOKUNO-A1-CURRENT-RELATIVE-SWIRL-SPLIT",
    "pr": 1171,
    "head": "657818dd83e119e6091bd2490d3804c04c3ef723",
    "branch": "codex/kokuno-agent1-relative-swirl-split",
    "base_pr": 1168,
    "base_head": "bb246140c6c2451b22b59aad9ba90604bc2acb21",
    "source_path": "src/openai_ns_reconstruction/kokuno_pa16_current_cartesian_relative_swirl_split.py",
    "source_blob": "c6d7d93ee33e5551513bdc817b2f118e4c544298",
    "stage": "relative-swirl-base-plus-sub-epsilon-delta-representation",
    "exact_compare_ahead": 4,
    "exact_compare_behind": 0,
    "open_pr_mergeable_at_registration": True,
    "relative_swirl_profile_delta_channel_materialized": True,
    "relative_swirl_cartesian_delta_channel_materialized": True,
    "binary64_total_relative_swirl_sum_is_resolved": False,
    "current_cartesian_relative_swirl_composed": False,
    "terminal_global_leading_velocity_materialized": False,
    "matched_global_pressure_materialized": False,
    "restricted_forcing_materialized": False,
    "heldout_ns_residual_assessed": False,
    "scientifically_admitted": False,
    "pde_validated": False,
    "dedicated_run": 35712062236,
    "repository_tests_run": 35712062144,
    "actions_status_at_registration": "queued",
}

AGENT2_AZIMUTHAL_FRAME_ORIENTATION = {
    "task": "K2-OSC-108",
    "pr": 1170,
    "head": "59d3e6d64fcdfeaf145cbee67e68b8a25e89489c",
    "branch": "kokuno-agent2-azimuthal-frame-orientation-108",
    "base_pr": 1164,
    "base_head": "0cbf74972af6589f7eb9a2edfe21c957d7512a5c",
    "source_path": "src/openai_ns_reconstruction/kokuno_source_azimuthal_frame_oriented_multiharmonic_velocity.py",
    "source_blob": "ca8acae2c3a53b79528cb55f45502443f6ffb2b6",
    "stage": "provider-driven-bounded-rigid-z-frame-oriented-complete-curl-family",
    "exact_compare_ahead": 3,
    "exact_compare_behind": 0,
    "open_pr_mergeable_at_registration": True,
    "bounded_azimuthal_frame_orientation_materialized": True,
    "same_rigid_rotation_applied_to_vector_potential_and_complete_curl": True,
    "rigid_rotation_preserves_divergence_free_contract_by_covariance": True,
    "orientation_classification": "repository_autonomous_candidate_design",
    "source_exact_frame_vectors_recovered": False,
    "source_exact_orientation_recovered": False,
    "source_provider_self_contained": False,
    "self_contained_velocity_xyzt_provider": False,
    "complete_ns_residual_assessed": False,
    "residual_reduction_claimed": False,
    "scientifically_admitted": False,
    "pde_validated": False,
    "dedicated_run": 35711831391,
    "repository_tests_run": 35711831373,
    "actions_status_at_registration": "queued",
}

AGENT3_AUXILIARY_T2_HAAR_BRIDGE = {
    "task": "KOKUNO-A3-AUXT2-HAAR-RF30-BRIDGE-127",
    "pr": 1172,
    "head": "9b0cecdf7eb9195bd468c3ddcd026c9d28d8e8d7",
    "branch": "agent3/rf30-auxiliary-t2-haar-bridge-20260922",
    "base_pr": 1161,
    "base_head": "cafab2098f00c3dc2ed5668e1b12bdb38a1f31e8",
    "source_path": "src/openai_ns_reconstruction/kokuno_rf30_auxiliary_t2_haar_bridge.py",
    "source_blob": "5b20d4945f3aeb177e327e8d6b5709a6726c885a",
    "stage": "normalized-auxiliary-t2-haar-covariance-operator-bridge",
    "exact_compare_ahead": 3,
    "exact_compare_behind": 0,
    "open_pr_mergeable_at_registration": True,
    "normalized_source_auxiliary_t2_haar_operator_executable": True,
    "haar_measure_total_mass_one": True,
    "covariance_formed_in_bridge_from_raw_wave_samples": True,
    "preaveraged_covariance_input_exposed": False,
    "repository_provider_blob_pinned": False,
    "rf30_repository_candidate_state_authorized": False,
    "rf30_defect_materialized": False,
    "rf31_or_rf34_rf39_run_on_real_candidate": False,
    "heldout_ns_residual_assessed": False,
    "residual_reduction_claimed": False,
    "scientifically_admitted": False,
    "pde_validated": False,
    "dedicated_run": 35712706101,
    "repository_tests_run": 35712706063,
    "actions_status_at_registration": "queued",
}

AGENT4_RELATIVE_SWIRL_SPLIT_AUDIT = {
    "task": "K4-VAL-121",
    "pr": 1174,
    "head": "8c9149e9ce9d198b23e8fcaf214d5d21f4b592ec",
    "branch": "kokuno-agent4/relative-swirl-split-audit-121",
    "base_pr": 1171,
    "base_head": "657818dd83e119e6091bd2490d3804c04c3ef723",
    "source_path": "src/openai_ns_reconstruction/kokuno_a4_relative_swirl_split_independent_audit.py",
    "source_blob": "307eb9f1f4f002c690bba55c18fe804ef0b6c24e",
    "stage": "implementation-distinct-save-reload-split-cartesian-divergence-audit",
    "exact_compare_ahead": 3,
    "exact_compare_behind": 0,
    "open_pr_mergeable_at_registration": True,
    "audited_agent1_pr": 1171,
    "audited_agent1_head": "657818dd83e119e6091bd2490d3804c04c3ef723",
    "scientific_path_consumes_public_velocity_split_only_after_save_load": True,
    "split_relative_swirl_cartesian_channel_independently_audited": True,
    "canonical_absolute_fd_transfer_ready": False,
    "scoped_gate_passed": None,
    "complete_ns_residual_assessed": False,
    "final_project_admission_ready": False,
    "scientifically_admitted": False,
    "pde_validated": False,
    "dedicated_run": 35713472612,
    "repository_tests_run": 35713472711,
    "actions_status_at_registration": "queued",
}

_PARENT_REGISTRATION = parent_a5.build_registration()
PARENT_FRONTIERS = copy.deepcopy(_PARENT_REGISTRATION["frontiers"])
CORRECTED_SOURCE = copy.deepcopy(_PARENT_REGISTRATION["corrected_source"])
FROZEN_GATES = copy.deepcopy(_PARENT_REGISTRATION["frozen_gates"])
CORE_STATE = copy.deepcopy(_PARENT_REGISTRATION["core_state"])


def _truth_boundary() -> dict[str, Any]:
    truth = copy.deepcopy(_PARENT_REGISTRATION["truth_boundary"])
    truth.update(
        {
            "relative_swirl_profile_delta_channel_materialized": True,
            "relative_swirl_cartesian_delta_channel_materialized": True,
            "binary64_total_relative_swirl_sum_is_resolved": False,
            "current_cartesian_relative_swirl_composed": False,
            "terminal_global_leading_velocity_materialized": False,
            "agent1_1171_scientifically_admitted": False,
            "agent4_matching_agent1_1171_split_audit_present": True,
            "agent4_matching_agent1_1171_split_audit_registered": True,
            "agent4_matching_agent1_1171_split_scoped_gate_passed": None,
            "agent4_matching_agent1_1171_split_canonical_fd_transfer_ready": False,
            "bounded_azimuthal_frame_orientation_materialized": True,
            "source_exact_frame_vectors_recovered": False,
            "source_exact_orientation_recovered": False,
            "oriented_source_provider_self_contained": False,
            "oriented_self_contained_velocity_xyzt_provider": False,
            "agent4_matching_agent2_1170_orientation_audit_present": False,
            "latest_self_contained_project_composite_is_agent2_1117": True,
            "current_i4_source_auxiliary_t2_haar_operator_executable": True,
            "current_i4_source_auxiliary_t2_provider_blob_pinned": False,
            "current_i4_source_auxiliary_t2_haar_mean_materialized_for_repository_candidate": False,
            "current_i4_rf30_repository_candidate_state_authorized": False,
            "current_i4_rf30_defect_materialized": False,
            "agent4_matching_agent3_1172_bridge_audit_present": False,
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
            "agent1_1171_split_delta_not_unified_float64_cartesian_composition": True,
            "agent1_1171_split_delta_not_terminal_or_global_leading": True,
            "agent4_1174_evidence_bound_to_exact_agent1_1171_split_identity": True,
            "agent4_1174_relative_fd_ladder_not_replacement_for_canonical_project_fd": True,
            "agent4_1174_split_divergence_evidence_not_complete_ns_validation": True,
            "agent2_1170_orientation_is_repository_autonomous_not_source_exact": True,
            "agent2_1170_provider_family_not_self_contained_project_candidate": True,
            "agent2_1170_covariance_contract_not_independent_pde_validation": True,
            "agent3_1172_haar_operator_not_repository_candidate_without_pinned_provider": True,
            "agent3_1172_raw_wave_bridge_not_rf30_defect_or_correction": True,
            "agent4_1174_not_evidence_for_agent2_1170_or_agent3_1172": True,
            "latest_self_contained_project_composite_stays_agent2_1117": True,
            "cross_identity_evidence_transfer_allowed": False,
        }
    )
    return firewall


SHORTEST_CLOSURE = [
    "let exact Agent-4 #1174 finish its preregistered split-channel audit, but do not treat it as unified-field or full-NS admission; retain the A1 #1171 split with a representation that does not round away the delta",
    "finish terminal/exterior/global leading and matched pressure on the same semantic identity",
    "turn Agent-2 provider-driven orientation/support assets into an exact-identity self-contained leading-plus-oscillatory project candidate; autonomous orientation must remain provenance-labeled",
    "supply and checksum-pin an implementation-distinct same-identity Agent-2 raw auxiliary-T2 provider to Agent-3 #1172; only then authorize RF30 candidate covariance/defect and run RF31/RF34-RF49",
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
    frontiers["relative_swirl_split_delta"] = copy.deepcopy(AGENT1_RELATIVE_SWIRL_SPLIT)
    frontiers["azimuthal_frame_oriented_multiharmonic"] = copy.deepcopy(
        AGENT2_AZIMUTHAL_FRAME_ORIENTATION
    )
    frontiers["auxiliary_t2_haar_rf30_bridge"] = copy.deepcopy(AGENT3_AUXILIARY_T2_HAAR_BRIDGE)
    frontiers["relative_swirl_split_validator"] = copy.deepcopy(AGENT4_RELATIVE_SWIRL_SPLIT_AUDIT)
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
    for key in PARENT_FRONTIERS:
        if frontiers.get(key) != PARENT_FRONTIERS[key]:
            raise ValueError(f"inherited frontier drifted: {key}")
    if frontiers.get("relative_swirl_split_delta") != AGENT1_RELATIVE_SWIRL_SPLIT:
        raise ValueError("Agent-1 #1171 identity drifted")
    if frontiers.get("azimuthal_frame_oriented_multiharmonic") != AGENT2_AZIMUTHAL_FRAME_ORIENTATION:
        raise ValueError("Agent-2 #1170 identity drifted")
    if frontiers.get("auxiliary_t2_haar_rf30_bridge") != AGENT3_AUXILIARY_T2_HAAR_BRIDGE:
        raise ValueError("Agent-3 #1172 identity drifted")
    if frontiers.get("relative_swirl_split_validator") != AGENT4_RELATIVE_SWIRL_SPLIT_AUDIT:
        raise ValueError("Agent-4 #1174 identity drifted")

    a1 = frontiers["relative_swirl_split_delta"]
    if a1["relative_swirl_cartesian_delta_channel_materialized"] is not True:
        raise ValueError("Agent-1 split delta delivery was dropped")
    if a1["binary64_total_relative_swirl_sum_is_resolved"] is not False:
        raise ValueError("sub-epsilon representation blocker was laundered")
    if a1["current_cartesian_relative_swirl_composed"] is not False:
        raise ValueError("split delta was laundered into unified Cartesian composition")
    if a1["scientifically_admitted"] is not False:
        raise ValueError("A1 split representation was laundered into admission")

    a2 = frontiers["azimuthal_frame_oriented_multiharmonic"]
    if a2["bounded_azimuthal_frame_orientation_materialized"] is not True:
        raise ValueError("A2 orientation delivery was dropped")
    if a2["source_exact_orientation_recovered"] is not False:
        raise ValueError("autonomous orientation was laundered into source-exact data")
    if a2["self_contained_velocity_xyzt_provider"] is not False:
        raise ValueError("provider-driven A2 family was laundered into project candidate")
    if a2["complete_ns_residual_assessed"] is not False:
        raise ValueError("A2 orientation covariance was laundered into NS validation")

    a3 = frontiers["auxiliary_t2_haar_rf30_bridge"]
    if a3["normalized_source_auxiliary_t2_haar_operator_executable"] is not True:
        raise ValueError("A3 normalized-Haar operator delivery was dropped")
    if a3["preaveraged_covariance_input_exposed"] is not False:
        raise ValueError("pre-averaged covariance bypass was introduced")
    if a3["repository_provider_blob_pinned"] is not False:
        raise ValueError("unavailable A2 auxiliary-T2 provider was fabricated")
    if a3["rf30_repository_candidate_state_authorized"] is not False:
        raise ValueError("RF30 candidate state was authorized without pinned provider")
    if a3["rf30_defect_materialized"] is not False:
        raise ValueError("A3 Haar mechanics was laundered into RF30 defect")

    a4 = frontiers["relative_swirl_split_validator"]
    if a4["audited_agent1_pr"] != 1171 or a4["audited_agent1_head"] != AGENT1_RELATIVE_SWIRL_SPLIT["head"]:
        raise ValueError("Agent-4 split audit is not bound to exact A1 #1171")
    if a4["scientific_path_consumes_public_velocity_split_only_after_save_load"] is not True:
        raise ValueError("Agent-4 implementation-distinct public path was weakened")
    if a4["scoped_gate_passed"] is not None:
        raise ValueError("queued Agent-4 audit was laundered into a gate conclusion")
    if a4["complete_ns_residual_assessed"] is not False or a4["pde_validated"] is not False:
        raise ValueError("scoped split audit was laundered into PDE validation")

    if frontiers["latest_self_contained_project_composite"]["pr"] != 1117:
        raise ValueError("provider-driven orientation was laundered into self-contained candidate")

    expected_core = {
        "leading_ready": False,
        "oscillatory_ready": True,
        "correction_ready": False,
        "velocity_export_ready": False,
        "pde_validated": False,
    }
    if reg["core_state"] != expected_core:
        raise ValueError("core state no longer matches frozen Agent-5 readiness")

    gates = reg["frozen_gates"]
    if gates["st006_momentum_sampled_max"] != 0.1082289305112118:
        raise ValueError("ST006 sampled-max baseline changed")
    if gates["st006_momentum_volume_l2"] != 0.10758432876230622:
        raise ValueError("ST006 volume-L2 baseline changed")
    if gates["normalized_momentum_sampled_max"] != 1e-3:
        raise ValueError("momentum sampled-max gate changed")
    if gates["normalized_momentum_volume_l2"] != 1e-3:
        raise ValueError("momentum volume-L2 gate changed")
    if gates["normalized_divergence_sampled_max"] != 1e-5:
        raise ValueError("divergence sampled-max gate changed")
    if gates["normalized_divergence_volume_l2"] != 1e-5:
        raise ValueError("divergence volume-L2 gate changed")
    if gates["canonical_quadrature"] != [24, 48, 96]:
        raise ValueError("canonical quadrature changed")
    if gates["residual_defined_free_forcing_allowed"] is not False:
        raise ValueError("free-forcing prohibition was weakened")

    truth = reg["truth_boundary"]
    if truth["agent4_matching_agent1_1171_split_audit_present"] is not True:
        raise ValueError("matching A4 #1174 audit delivery was dropped")
    if truth["agent4_matching_agent1_1171_split_scoped_gate_passed"] is not None:
        raise ValueError("unresolved A4 #1174 gate was pre-claimed")
    if truth["current_i4_source_auxiliary_t2_haar_operator_executable"] is not True:
        raise ValueError("A3 executable Haar operator truth was dropped")

    required_false_truth = (
        "binary64_total_relative_swirl_sum_is_resolved",
        "current_cartesian_relative_swirl_composed",
        "terminal_global_leading_velocity_materialized",
        "agent1_1171_scientifically_admitted",
        "agent4_matching_agent1_1171_split_canonical_fd_transfer_ready",
        "source_exact_frame_vectors_recovered",
        "source_exact_orientation_recovered",
        "oriented_source_provider_self_contained",
        "oriented_self_contained_velocity_xyzt_provider",
        "agent4_matching_agent2_1170_orientation_audit_present",
        "current_i4_source_auxiliary_t2_provider_blob_pinned",
        "current_i4_source_auxiliary_t2_haar_mean_materialized_for_repository_candidate",
        "current_i4_rf30_repository_candidate_state_authorized",
        "current_i4_rf30_defect_materialized",
        "agent4_matching_agent3_1172_bridge_audit_present",
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
    if any(truth[name] is not False for name in required_false_truth):
        raise ValueError("truth boundary was promoted")


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
