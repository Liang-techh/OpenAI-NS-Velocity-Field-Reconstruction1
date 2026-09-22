"""Fail-closed Kokuno Agent-5 routing after post-pulse/A2-multiharmonic/A4-xi11-audit advances.

Integration/provenance glue only. This module does not rewrite Agent-1--4
mathematics and never transfers evidence across candidate identities.

Fresh frontiers frozen here:

* Agent 1 #1140 carries the current Cartesian leading field from the public
  xi=13 pulse endpoint through one bounded post-pulse eta-flattening stage.
  It is callable/save-loadable on that bounded stage, but later hold,
  relative-swirl, exterior-heat and global completion remain absent.
* Agent 2 #1141 composes up to four identity-bound provider-driven
  complete-curl harmonics with bounded repository-autonomous complex
  coefficients. The source/background/support providers remain external, so
  this is not a new self-contained project candidate.
* Agent 3 #1134 remains source-specific RF34--RF39 mechanics only until a
  real same-identity current-I4 fixed-Q RF30 backend exists.
* Agent 4 #1136 now independently audits exact self-contained A2 #1117 /
  A1 #1107 xi=11 composite divergence after save/load. Its exact-head CI is
  still queued at registration, so this records presence/registration only,
  never PASS or scientific admission.

The final project gates remain fixed: normalized momentum sampled-max and
volume-L2 <= 1e-3; normalized divergence sampled-max and volume-L2 <= 1e-5;
canonical quadrature [24,48,96]; residual-defined/free forcing forbidden.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from . import kokuno_a5_xi13_axissafe_rf39_routing as parent_a5

SCHEMA = "kokuno-a5-postpulse-multiharmonic-xi11-audit-routing-v1"
TASK = "KOKUNO-A5-POSTPULSE-MULTIHARMONIC-XI11-AUDIT-ROUTING-116"

PARENT_A5 = {
    "pr": 1135,
    "head": "a6976c654926c3c96364c66bb18eac461688bde7",
    "branch": "kokuno-agent5/xi13-axissafe-rf39-routing-115",
    "source_path": "src/openai_ns_reconstruction/kokuno_a5_xi13_axissafe_rf39_routing.py",
    "source_blob": "a04996365234f33e81172a2d6a4f6bee8a7da993",
    "task": "KOKUNO-A5-XI13-AXISSAFE-RF39-ROUTING-115",
    "schema": "kokuno-a5-xi13-axissafe-rf39-routing-v1",
}

AGENT1_POSTPULSE_ETA_FLATTENING = {
    "pr": 1140,
    "head": "c5442c11165d7f0976889b826bf58dce38a3d3dd",
    "branch": "codex/kokuno-a1-postpulse-eta-flattening",
    "source_path": "src/openai_ns_reconstruction/kokuno_pa16_current_cartesian_postpulse_eta_flattening.py",
    "source_blob": "32cbfc00406aa8d47904e2492ad8916730541d28",
    "parent_pr": 1133,
    "stage": "current-cartesian-leading-bounded-postpulse-eta-flattening",
    "public_parent_pulse_endpoint_xi": 13.0,
    "autonomous_T_f": 100.0,
    "source_exact_T_f_recovered": False,
    "save_load_materialized": True,
    "bounded_velocity_xyzt_materialized": True,
    "later_terminal_hold_relative_swirl_exterior_global_materialized": False,
    "terminal_global_leading_velocity_materialized": False,
    "repository_tests_run": 35695446759,
    "dedicated_run": 35695446964,
    "actions_status_at_registration": "queued",
}

AGENT2_BOUNDED_MULTIHARMONIC = {
    "pr": 1141,
    "head": "4241b0b4fcdc36e750798133101baaf803ba8c44",
    "branch": "codex/k2-osc-104-bounded-multiharmonic",
    "source_path": "src/openai_ns_reconstruction/kokuno_source_bounded_multiharmonic_velocity.py",
    "source_blob": "06bc24722b4219419e91cd2d5b08cb1294665f24",
    "parent_pr": 1132,
    "stage": "provider-driven-bounded-multiharmonic-complete-curl-family",
    "max_terms": 4,
    "amplitude_l1_cap": 1.0,
    "nontriviality_floor": 1.0e-12,
    "complete_curl_family_materialized": True,
    "provider_driven_velocity_xyzt_materialized": True,
    "source_exact_amplitude_phase_recovered": False,
    "project_domain_source_input_provider_materialized": False,
    "self_contained_velocity_xyzt_provider": False,
    "repository_tests_run": 35695804431,
    "dedicated_run": 35695804435,
    "actions_status_at_registration": "queued",
}

LATEST_SELF_CONTAINED_PROJECT_COMPOSITE = {
    "pr": 1117,
    "head": "27741d9c0a27262f7fabf61eebaa2fbd507e9f03",
    "consumed_agent1_pr": 1107,
    "consumed_agent1_head": "45da043dd2b4cd067f005a72c7e21fd0f2bcf309",
    "stage": "project-domain-logX-leading-plus-frozen-complete-curl-through-xi11",
    "public_pulse_endpoint_xi": 11.0,
    "save_load_materialized": True,
    "matching_agent4_audit_present": True,
    "matching_agent4_audit_pr": 1136,
    "matching_agent4_audit_scoped_gate_passed": None,
    "matching_postpulse_agent2_composite": False,
}

AGENT4_XI11_COMPOSITE_DIVERGENCE_AUDIT = {
    "task": "K4-VAL-117",
    "pr": 1136,
    "head": "03269e7ba2053e0bc35813089577d274d9ddf8f9",
    "branch": "kokuno-agent4/logx-xi11-composite-divergence-117",
    "source_path": "src/openai_ns_reconstruction/kokuno_a4_logx_xi11_composite_divergence_independent_audit.py",
    "source_blob": "bd605bfb650b7eadbebb8d860e9720b434540ba3",
    "audited_agent2_pr": 1117,
    "audited_agent2_head": "27741d9c0a27262f7fabf61eebaa2fbd507e9f03",
    "audited_agent1_pr": 1107,
    "audited_agent1_head": "45da043dd2b4cd067f005a72c7e21fd0f2bcf309",
    "stage": "xi11-self-contained-composite-independent-divergence-audit",
    "fd_steps": [0.02, 0.01, 0.005],
    "seed": 9173891,
    "save_load_roundtrip": True,
    "implementation_distinct": True,
    "present": True,
    "registered": True,
    "scoped_gate_passed": None,
    "scientifically_admitted": False,
    "complete_ns_momentum_residual_assessed": False,
    "pde_validated": False,
    "repository_tests_run": 35692557189,
    "dedicated_run": 35692557216,
    "actions_status_at_registration": "queued",
}

AGENT3_RF34_RF39_MECHANICS = copy.deepcopy(parent_a5.AGENT3_RF34_RF39_MECHANICS)
CURRENT_I4_REAL_RADIAL_FORCE = copy.deepcopy(parent_a5.CURRENT_I4_REAL_RADIAL_FORCE)
AGENT4_CURRENT_I4_RADIAL_FORCE_AUDIT = copy.deepcopy(parent_a5.AGENT4_CURRENT_I4_RADIAL_FORCE_AUDIT)

CORRECTED_SOURCE = copy.deepcopy(parent_a5.CORRECTED_SOURCE)
FROZEN_GATES = copy.deepcopy(parent_a5.FROZEN_GATES)
CORE_STATE = copy.deepcopy(parent_a5.CORE_STATE)


def _truth_boundary() -> dict[str, Any]:
    return {
        "postpulse_eta_flattening_leading_materialized": True,
        "postpulse_eta_flattening_save_load_materialized": True,
        "source_exact_postpulse_T_f_recovered": False,
        "terminal_global_leading_velocity_materialized": False,
        "matching_postpulse_agent2_project_composite_materialized": False,
        "matching_postpulse_agent4_audit_present": False,
        "bounded_multiharmonic_complete_curl_family_materialized": True,
        "bounded_multiharmonic_provider_driven_velocity_xyzt_materialized": True,
        "bounded_multiharmonic_self_contained_velocity_xyzt_provider": False,
        "source_exact_multiharmonic_amplitude_phase_recovered": False,
        "project_domain_source_input_provider_materialized": False,
        "latest_self_contained_project_composite_is_agent2_1117": True,
        "agent4_matching_xi11_composite_audit_present": True,
        "agent4_matching_xi11_composite_audit_registered": True,
        "agent4_matching_xi11_composite_audit_scoped_gate_passed": None,
        "agent4_matching_xi11_composite_audit_admitted": False,
        "rf34_rf39_compact_mechanics_materialized": True,
        "current_i4_source_chart_backend_materialized": False,
        "current_i4_rf30_defect_materialized": False,
        "current_i4_rf34_rf39_correction_materialized": False,
        "cartesian_correction_velocity_materialized": False,
        "current_i4_radial_force_materialized": True,
        "agent4_matching_current_i4_radial_force_audit_present": True,
        "agent4_matching_current_i4_radial_force_audit_admitted": False,
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
        "agent1_1140_postpulse_not_consumed_by_agent2_1117_xi11": True,
        "agent2_1141_provider_family_not_self_contained_project_composite": True,
        "agent4_1136_xi11_audit_not_postpulse_or_global_or_full_ns_evidence": True,
        "agent3_1134_mechanics_not_current_i4_candidate_residual_evidence": True,
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
            "leading": copy.deepcopy(AGENT1_POSTPULSE_ETA_FLATTENING),
            "bounded_multiharmonic_sibling": copy.deepcopy(AGENT2_BOUNDED_MULTIHARMONIC),
            "latest_self_contained_project_composite": copy.deepcopy(LATEST_SELF_CONTAINED_PROJECT_COMPOSITE),
            "xi11_composite_validator": copy.deepcopy(AGENT4_XI11_COMPOSITE_DIVERGENCE_AUDIT),
            "rf34_rf39_mechanics_sibling": copy.deepcopy(AGENT3_RF34_RF39_MECHANICS),
            "current_i4_real_radial_force": copy.deepcopy(CURRENT_I4_REAL_RADIAL_FORCE),
            "current_i4_radial_force_validator": copy.deepcopy(AGENT4_CURRENT_I4_RADIAL_FORCE_AUDIT),
        },
        "identity_firewall": build_identity_firewall(),
        "truth_boundary": _truth_boundary(),
        "core_state": copy.deepcopy(CORE_STATE),
        "frozen_gates": copy.deepcopy(FROZEN_GATES),
        "shortest_closure": [
            "finish the remaining A1 terminal/relative-swirl/exterior/global leading continuation without relabeling the bounded post-pulse stage as global",
            "bind the eventual stable global A1 identity into one self-contained A2 leading+oscillatory candidate; A2 #1141 remains provider-driven until its source inputs are materialized",
            "materialize the real same-identity fixed-Q RF30 covariance/derivative backend before A3 #1134 can become candidate-specific correction evidence",
            "bind matched pressure and preregistered non-residual-defined restricted forcing before constructing a complete NS defect",
            "only then run the real finite correction cycle and Agent-4 held-out full-NS 1e-3 gate",
        ],
    }
    registration["registration_sha256"] = registration_sha256(registration)
    validate_registration(registration)
    return registration


def validate_registration(registration: Mapping[str, Any]) -> None:
    reg = dict(registration)
    if reg.get("registration_sha256") != registration_sha256(reg):
        raise ValueError("registration digest mismatch")
    if reg.get("schema") != SCHEMA or reg.get("task") != TASK:
        raise ValueError("schema/task identity drift")
    if reg.get("parent_a5") != PARENT_A5:
        raise ValueError("parent A5 identity drift")
    if reg.get("corrected_source") != CORRECTED_SOURCE:
        raise ValueError("corrected-source provenance drift")

    frontiers = reg.get("frontiers")
    if not isinstance(frontiers, Mapping):
        raise ValueError("frontiers missing")
    expected_frontiers = {
        "leading": AGENT1_POSTPULSE_ETA_FLATTENING,
        "bounded_multiharmonic_sibling": AGENT2_BOUNDED_MULTIHARMONIC,
        "latest_self_contained_project_composite": LATEST_SELF_CONTAINED_PROJECT_COMPOSITE,
        "xi11_composite_validator": AGENT4_XI11_COMPOSITE_DIVERGENCE_AUDIT,
        "rf34_rf39_mechanics_sibling": AGENT3_RF34_RF39_MECHANICS,
        "current_i4_real_radial_force": CURRENT_I4_REAL_RADIAL_FORCE,
        "current_i4_radial_force_validator": AGENT4_CURRENT_I4_RADIAL_FORCE_AUDIT,
    }
    if dict(frontiers) != expected_frontiers:
        raise ValueError("frontier identity drift")
    if reg.get("identity_firewall") != build_identity_firewall():
        raise ValueError("identity firewall drift")
    if reg.get("core_state") != CORE_STATE:
        raise ValueError("core readiness promotion/drift")
    if reg.get("frozen_gates") != FROZEN_GATES:
        raise ValueError("frozen project gates changed")
    if reg.get("truth_boundary") != _truth_boundary():
        raise ValueError("truth boundary promotion/drift")


def write_registration(path: str | Path) -> dict[str, Any]:
    registration = build_registration()
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(registration, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return registration


def load_registration(path: str | Path) -> dict[str, Any]:
    registration = json.loads(Path(path).read_text(encoding="utf-8"))
    validate_registration(registration)
    return registration


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    registration = write_registration(args.output)
    print(json.dumps({"registration_sha256": registration["registration_sha256"], "output": str(args.output)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
