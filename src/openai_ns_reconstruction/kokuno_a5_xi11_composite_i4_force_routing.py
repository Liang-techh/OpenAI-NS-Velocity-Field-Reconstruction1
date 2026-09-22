"""Fail-closed Agent-5 routing checkpoint for the xi=11 composite / I4 force seam.

Integration/provenance glue only.  The current frontiers are deliberately kept
separate:

* A1 #1116 materializes public low-dimensional c1/c2 pulse-end compensator
  algebra, but does not compose it into the current Cartesian field.
* A2 #1117 composes exact A1 #1107's overflow-safe log-X leading field through
  xi=11 with the frozen complete-curl oscillation and supplies deterministic
  save/load for that project-domain candidate identity.
* A3 #1118 advances the older exact current-I4 identity from compact axial
  stress to the source-required radial force (div T)_r = partial_z sigma_1.
* A4 #1119 independently audits exact A3 #1109 current-I4 compact stress, but
  not the newer #1118 radial force.  A4 #1110 remains finite-X A2 #1108
  divergence evidence and does not audit A2 #1117.

No scoped evidence is transferred across candidate identities.  No pressure,
forcing, residual, correction coefficient, optimization sample, or scientific
threshold is accepted by this module.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import re
from pathlib import Path
from typing import Any, Mapping

from . import kokuno_a5_xi11_pulse_i4_stress_routing as parent_a5

SCHEMA = "kokuno-a5-xi11-composite-i4-force-routing-v1"
TASK = "KOKUNO-A5-XI11-COMPOSITE-I4-FORCE-ROUTING-113"
_HEX40 = re.compile(r"^[0-9a-f]{40}$")
_HEX64 = re.compile(r"^[0-9a-f]{64}$")

PARENT_A5 = {
    "pr": 1111,
    "head": "f35bbc55455bb7bd6773c4e273adad896e744132",
    "source_path": "src/openai_ns_reconstruction/kokuno_a5_xi11_pulse_i4_stress_routing.py",
    "source_blob": "5572b3ffc7fa73d51c5a15d9dc76eb9eb6dc72af",
    "task": "KOKUNO-A5-XI11-PULSE-I4-STRESS-ROUTING-112",
    "schema": "kokuno-a5-xi11-pulse-i4-stress-routing-v1",
}

AGENT1_PULSE_END_ALGEBRA = {
    "pr": 1116,
    "head": "c144d00fd81a938e497ca5d841fed6d1fe448a8b",
    "source_path": "src/openai_ns_reconstruction/kokuno_public_pulse_end_compensator.py",
    "source_blob": "3c4d7a2d51cbccad7322a9bf6825a69e2b45e486",
    "parent_pr": 1107,
    "parent_head": "45da043dd2b4cd067f005a72c7e21fd0f2bcf309",
    "stage": "public-low-dimensional-pulse-end-compensator-algebra",
    "corrected_two_row_c1_c2_system_materialized": True,
    "public_bump_centers_width_bound": True,
    "source_exact_bump_shape_recovered": False,
    "source_exact_amplitude_root_materialized": False,
    "current_lineage_J_materialized": False,
    "cartesian_c1_c2_composition_materialized": False,
    "terminal_global_leading_velocity_materialized": False,
    "consumed_by_latest_agent2_composite": False,
    "exact_head_tests_run": 35683731026,
    "exact_head_dedicated_run": 35683731031,
    "actions_status_at_registration": "queued",
}

AGENT2_XI11_LOGX_COMPOSITE = {
    "pr": 1117,
    "head": "27741d9c0a27262f7fabf61eebaa2fbd507e9f03",
    "source_path": "src/openai_ns_reconstruction/kokuno_current_main_pulse_logx_leading_oscillatory_identity.py",
    "source_blob": "5504c2cde2eb69bfcdb4e4a87cc48f541bb9f8d3",
    "parent_pr": 1108,
    "parent_head": "9f221dd57ef4b5e2d2e80e24c8cf531e991a5387",
    "consumed_agent1_pr": 1107,
    "consumed_agent1_head": "45da043dd2b4cd067f005a72c7e21fd0f2bcf309",
    "candidate_stage": "current-main-pulse-logX-through-xi11",
    "self_contained_project_cartesian_velocity": True,
    "leading_plus_frozen_complete_curl_oscillation_materialized": True,
    "identity_preserving_save_load_available": True,
    "public_principal_xi11_interval_materialized": True,
    "repository_autonomous_principal_amplitude_consumed": True,
    "source_exact_amplitude_root_materialized": False,
    "pulse_end_compensator_algebra_1116_consumed": False,
    "terminal_global_leading_velocity_materialized": False,
    "complete_ns_residual_assessed": False,
    "exact_head_tests_run": 35684187222,
    "exact_head_dedicated_run": 35684187205,
    "actions_status_at_registration": "queued",
}

AGENT3_CURRENT_I4_RADIAL_FORCE = {
    "pr": 1118,
    "head": "922f7aa10460ded44af212d313eceff33a2ac647",
    "source_path": "src/openai_ns_reconstruction/kokuno_current_i4_nonlinear_radial_force.py",
    "source_blob": "e28eac6f7fb7f3191a024b052bcf677cff8d74d2",
    "parent_pr": 1109,
    "parent_head": "49590ff311fef1dec4bd850b3013fd989485c87a",
    "agent2_composite_pr": 1080,
    "agent2_composite_head": "c40d8ddecd2971544a6e07dab093436b423cf326",
    "agent1_leading_pr": 1079,
    "agent1_leading_head": "b06742ca6e189499192ede3cce40f62cdc1e35ca",
    "candidate_stage": "current-I4",
    "nonlinear_mean_attribution_consumed": True,
    "compact_theta_e2_stress_materialized": True,
    "compact_axial_e1_stress_materialized": True,
    "radial_force_partial_z_sigma_1_materialized": True,
    "physical_z_derivative_steps": [0.02, 0.01, 0.005],
    "production_derivative_operator": "centered-FD2",
    "source_i4_five_row_mean_correction_materialized": False,
    "matching_agent4_radial_force_audit_present": False,
    "authorized_as_complete_ns_correction_target": False,
    "exact_head_tests_run": 35684255309,
    "exact_head_dedicated_run": 35684255510,
    "actions_status_at_registration": "queued",
}

AGENT4_CURRENT_I4_STRESS_AUDIT = {
    "pr": 1119,
    "head": "3f3fd169ce79a2e44b90b50353415d057119e63f",
    "source_path": "src/openai_ns_reconstruction/kokuno_a4_current_i4_radial_stress_independent_audit.py",
    "source_blob": "97387bcab05098ecc622372ec3aa08150d79b89d",
    "audited_agent3_pr": 1109,
    "audited_agent3_head": "49590ff311fef1dec4bd850b3013fd989485c87a",
    "audited_agent2_pr": 1080,
    "audited_agent2_head": "c40d8ddecd2971544a6e07dab093436b423cf326",
    "audited_agent1_pr": 1079,
    "audited_agent1_head": "b06742ca6e189499192ede3cce40f62cdc1e35ca",
    "audited_agent3_1118_radial_force": False,
    "seed": 9173871,
    "operator": "local-cubic interpolation plus order-8 Gauss-Legendre cell integration",
    "implementation_distinct": True,
    "registered": True,
    "scientific_admission": False,
    "complete_ns_residual_evidence": False,
    "exact_head_tests_run": 35684778994,
    "exact_head_dedicated_run": 35684779061,
    "actions_status_at_registration": "queued",
}

AGENT4_FINITE_PREFIX_DIVERGENCE_AUDIT = {
    "pr": 1110,
    "head": "71432cec314c2ab11b1bf38c83609ad7d64a6654",
    "source_path": "src/openai_ns_reconstruction/kokuno_a4_current_main_pulse_composite_divergence_independent_audit.py",
    "source_blob": "50fc8557f3e0d06621c5b801221712a2c77d091a",
    "audited_agent2_pr": 1108,
    "audited_agent2_head": "9f221dd57ef4b5e2d2e80e24c8cf531e991a5387",
    "audited_agent1_pr": 1100,
    "audited_agent1_head": "2a6c58ff6bf6d1ba6b875bca0e6abf9d1118e156",
    "audited_agent2_1117": False,
    "audited_agent3_1118_radial_force": False,
    "seed": 9173861,
    "spatial_steps": [0.02, 0.01, 0.005],
    "operator": "public-Cartesian centered-FD2 after real candidate save/reload",
    "implementation_distinct": True,
    "registered": True,
    "scientific_admission": False,
    "complete_ns_residual_evidence": False,
    "exact_head_tests_run": 35681237729,
    "exact_head_dedicated_run": 35681237794,
    "actions_status_at_registration": "queued",
}

FINAL_GATES = copy.deepcopy(parent_a5.FINAL_GATES)
ST006_BASELINE = copy.deepcopy(parent_a5.ST006_BASELINE)
READINESS = copy.deepcopy(parent_a5.READINESS)

TRUTH_BOUNDARY = {
    "parent_a5_checkpoint_preserved": True,
    "leading_frontier_stage": "public-pulse-end-compensator-algebra-sibling",
    "latest_self_contained_composite_frontier_stage": "current-main-pulse-logX-through-xi11",
    "correction_frontier_stage": "current-I4-radial-force-partial-z-sigma-1",
    "agent1_public_pulse_end_compensator_algebra_materialized": True,
    "agent1_public_pulse_end_compensator_algebra_consumed_in_cartesian_velocity": False,
    "source_exact_main_pulse_amplitude_materialized": False,
    "current_lineage_J_materialized": False,
    "cartesian_pulse_end_compensation_materialized": False,
    "terminal_global_leading_completion_materialized": False,
    "matching_xi11_agent2_project_composite_materialized": True,
    "matching_xi11_agent2_project_composite_save_load_available": True,
    "matching_xi11_agent2_project_composite_consumes_exact_agent1_1107": True,
    "matching_xi11_agent2_project_composite_consumes_agent1_1116": False,
    "agent4_matching_xi11_composite_audit_present": False,
    "agent4_matching_xi11_composite_audit_registered": False,
    "agent4_1110_finite_prefix_audit_transferred_to_xi11_identity": False,
    "current_i4_nonlinear_mean_attribution_preserved": True,
    "current_i4_compact_radial_stress_materialized": True,
    "agent4_matching_current_i4_radial_stress_audit_present": True,
    "agent4_matching_current_i4_radial_stress_audit_registered": True,
    "agent4_matching_current_i4_radial_stress_audit_admitted": False,
    "current_i4_radial_force_materialized": True,
    "agent4_matching_current_i4_radial_force_audit_present": False,
    "agent4_matching_current_i4_radial_force_audit_registered": False,
    "source_i4_five_row_mean_correction_materialized": False,
    "current_i4_radial_force_authorized_as_complete_ns_correction_target": False,
    "current_i4_stress_audit_authorized_as_complete_ns_correction_target": False,
    "current_i4_force_evidence_transferred_to_xi11_composite": False,
    "xi11_composite_evidence_transferred_to_current_i4_force": False,
    "source_positive_order_i3_correction_materialized": False,
    "global_leading_plus_oscillatory_velocity_materialized": False,
    "matched_cartesian_pressure_materialized": False,
    "cartesian_matched_pressure_gradient_materialized": False,
    "preregistered_restricted_forcing_materialized": False,
    "forcing_proven_not_residual_defined": False,
    "complete_identity_bound_ns_defect_materialized": False,
    "current_cartesian_correction_velocity_materialized": False,
    "finite_correction_cycle_admitted": False,
    "real_candidate_finite_correction_cycle_run": False,
    "heldout_normalized_full_ns_residual_assessed": False,
    "canonical_24_48_96_whole_domain_admission_assessed": False,
    "same_protocol_comparable_to_st006": False,
    "scientific_admission": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}

PIPELINE_POSITION = {
    "leading_frontier": (
        "A1 #1116 exposes the public c1/c2 end-compensator algebra as a low-dimensional sibling, "
        "but source-exact Amp, current-lineage J and Cartesian composition remain absent"
    ),
    "latest_self_contained_composite_frontier": (
        "A2 #1117 consumes exact A1 #1107 and provides a save/loadable project-domain leading+"
        "frozen-complete-curl candidate through the public principal xi=11 endpoint"
    ),
    "correction_frontier": (
        "A3 #1118 advances the distinct older current-I4 A2 #1080 identity through compact stress "
        "to radial force partial_z sigma_1; source five-row I4 correction remains absent"
    ),
    "validator_frontier": (
        "A4 #1119 independently audits exact A3 #1109 current-I4 compact stress but not #1118 "
        "radial force; A4 #1110 independently audits only A2 #1108 finite-X divergence, not #1117"
    ),
    "identity_firewall": (
        "A1 #1116, A2 #1117, A3 #1118, A4 #1119 and A4 #1110 retain distinct semantic scopes; "
        "no scoped divergence/stress/force evidence is transferred across those boundaries"
    ),
    "next_shortest_blocker": (
        "obtain an implementation-distinct A4 audit for exact A2 #1117 and a separate A4 radial-"
        "force audit for exact A3 #1118; in parallel A1 must supply current-lineage J and Cartesian "
        "pulse-end composition before terminal/global leading completion. Global velocity, matched "
        "pressure and preregistered restricted forcing remain mandatory before a complete NS defect, "
        "real finite correction cycle, or final held-out 1e-3 PDE gate can be admitted"
    ),
}


def _canonical_json(payload: Mapping[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _sha256(payload: Mapping[str, Any]) -> str:
    return hashlib.sha256(_canonical_json(payload).encode()).hexdigest()


def _require_hex40(value: Any, label: str) -> None:
    if not isinstance(value, str) or _HEX40.fullmatch(value) is None:
        raise ValueError(f"{label} is not an exact git/blob identity")


def _validate_parent() -> None:
    if (parent_a5.SCHEMA, parent_a5.TASK) != (PARENT_A5["schema"], PARENT_A5["task"]):
        raise ValueError("parent A5 schema/task drifted")
    if parent_a5.FINAL_GATES != FINAL_GATES:
        raise ValueError("frozen final gates drifted")
    if parent_a5.ST006_BASELINE != ST006_BASELINE:
        raise ValueError("ST006 baseline drifted")
    if parent_a5.READINESS != READINESS:
        raise ValueError("core readiness drifted")


def build_registration() -> dict[str, Any]:
    _validate_parent()
    payload = {
        "parent_a5": copy.deepcopy(PARENT_A5),
        "agent1_pulse_end_algebra": copy.deepcopy(AGENT1_PULSE_END_ALGEBRA),
        "agent2_xi11_logx_composite": copy.deepcopy(AGENT2_XI11_LOGX_COMPOSITE),
        "agent3_current_i4_radial_force": copy.deepcopy(AGENT3_CURRENT_I4_RADIAL_FORCE),
        "agent4_current_i4_stress_audit": copy.deepcopy(AGENT4_CURRENT_I4_STRESS_AUDIT),
        "agent4_finite_prefix_divergence_audit": copy.deepcopy(AGENT4_FINITE_PREFIX_DIVERGENCE_AUDIT),
        "final_gates": copy.deepcopy(FINAL_GATES),
        "st006_baseline": copy.deepcopy(ST006_BASELINE),
        "readiness": copy.deepcopy(READINESS),
        "truth_boundary": copy.deepcopy(TRUTH_BOUNDARY),
        "pipeline_position": copy.deepcopy(PIPELINE_POSITION),
    }
    registration = {"schema": SCHEMA, "task": TASK, **payload, "digest": _sha256(payload)}
    validate_registration(registration)
    return registration


def validate_registration(registration: Mapping[str, Any]) -> None:
    _validate_parent()
    if registration.get("schema") != SCHEMA or registration.get("task") != TASK:
        raise ValueError("A5 registration identity drifted")

    exact_blocks = {
        "parent_a5": PARENT_A5,
        "agent1_pulse_end_algebra": AGENT1_PULSE_END_ALGEBRA,
        "agent2_xi11_logx_composite": AGENT2_XI11_LOGX_COMPOSITE,
        "agent3_current_i4_radial_force": AGENT3_CURRENT_I4_RADIAL_FORCE,
        "agent4_current_i4_stress_audit": AGENT4_CURRENT_I4_STRESS_AUDIT,
        "agent4_finite_prefix_divergence_audit": AGENT4_FINITE_PREFIX_DIVERGENCE_AUDIT,
    }
    for name, frozen in exact_blocks.items():
        if registration.get(name) != frozen:
            raise ValueError(f"{name} provenance drifted")
        _require_hex40(frozen["head"], f"{name}.head")
        _require_hex40(frozen["source_blob"], f"{name}.source_blob")

    if registration.get("final_gates") != FINAL_GATES:
        raise ValueError("final scientific gates drifted")
    if registration.get("st006_baseline") != ST006_BASELINE:
        raise ValueError("ST006 baseline drifted")
    if registration.get("readiness") != READINESS:
        raise ValueError("core readiness drifted")
    if registration.get("truth_boundary") != TRUTH_BOUNDARY:
        raise ValueError("truth boundary drifted")
    if registration.get("pipeline_position") != PIPELINE_POSITION:
        raise ValueError("pipeline routing drifted")

    if AGENT2_XI11_LOGX_COMPOSITE["consumed_agent1_head"] != AGENT1_PULSE_END_ALGEBRA["parent_head"]:
        raise ValueError("A2 xi11 composite is not bound to expected A1 #1107")
    if AGENT2_XI11_LOGX_COMPOSITE["pulse_end_compensator_algebra_1116_consumed"]:
        raise ValueError("A1 #1116 algebra was incorrectly promoted into A2 #1117 field")
    if AGENT4_FINITE_PREFIX_DIVERGENCE_AUDIT["audited_agent2_head"] == AGENT2_XI11_LOGX_COMPOSITE["head"]:
        raise ValueError("finite-prefix A4 evidence was laundered onto xi11 composite")
    if AGENT4_CURRENT_I4_STRESS_AUDIT["audited_agent3_head"] != AGENT3_CURRENT_I4_RADIAL_FORCE["parent_head"]:
        raise ValueError("A4 #1119 does not bind exact A3 #1109 stress parent")
    if AGENT4_CURRENT_I4_STRESS_AUDIT["audited_agent3_1118_radial_force"]:
        raise ValueError("A4 #1119 was incorrectly relabelled as #1118 radial-force audit")
    if AGENT3_CURRENT_I4_RADIAL_FORCE["authorized_as_complete_ns_correction_target"]:
        raise ValueError("scoped current-I4 force was incorrectly authorized as complete-NS target")

    truth = registration["truth_boundary"]
    forbidden_promotions = (
        "agent4_matching_xi11_composite_audit_present",
        "agent4_matching_current_i4_radial_force_audit_present",
        "current_i4_radial_force_authorized_as_complete_ns_correction_target",
        "current_i4_stress_audit_authorized_as_complete_ns_correction_target",
        "global_leading_plus_oscillatory_velocity_materialized",
        "matched_cartesian_pressure_materialized",
        "preregistered_restricted_forcing_materialized",
        "complete_identity_bound_ns_defect_materialized",
        "finite_correction_cycle_admitted",
        "heldout_normalized_full_ns_residual_assessed",
        "scientific_admission",
        "pde_validated",
    )
    if any(bool(truth[key]) for key in forbidden_promotions):
        raise ValueError("unsupported scientific-state promotion detected")
    if truth["agent4_matching_current_i4_radial_stress_audit_admitted"]:
        raise ValueError("scoped A4 stress audit cannot be scientific admission")
    if truth["agent4_1110_finite_prefix_audit_transferred_to_xi11_identity"]:
        raise ValueError("A4 finite-prefix evidence transfer is forbidden")
    if truth["current_i4_force_evidence_transferred_to_xi11_composite"]:
        raise ValueError("current-I4 correction evidence transfer is forbidden")
    if truth["xi11_composite_evidence_transferred_to_current_i4_force"]:
        raise ValueError("xi11 composite evidence transfer is forbidden")

    payload_keys = (
        "parent_a5",
        "agent1_pulse_end_algebra",
        "agent2_xi11_logx_composite",
        "agent3_current_i4_radial_force",
        "agent4_current_i4_stress_audit",
        "agent4_finite_prefix_divergence_audit",
        "final_gates",
        "st006_baseline",
        "readiness",
        "truth_boundary",
        "pipeline_position",
    )
    payload = {key: copy.deepcopy(registration[key]) for key in payload_keys}
    digest = registration.get("digest")
    if not isinstance(digest, str) or _HEX64.fullmatch(digest) is None or digest != _sha256(payload):
        raise ValueError("registration digest mismatch")


def save_registration(path: str | Path) -> dict[str, Any]:
    registration = build_registration()
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(registration, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return registration


def load_registration(path: str | Path) -> dict[str, Any]:
    registration = json.loads(Path(path).read_text(encoding="utf-8"))
    validate_registration(registration)
    return registration


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)
    registration = save_registration(args.output)
    print(json.dumps({"schema": SCHEMA, "task": TASK, "digest": registration["digest"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
