"""Fail-closed Agent-5 routing checkpoint for xi=11 pulse and current-I4 stress.

This module is integration/provenance glue only. It harvests the latest bounded
Agent-1/2/3/4 deliveries while preserving their exact candidate identities:
A1 #1107 reaches the public main-pulse endpoint xi=11 in log-X; A2 #1108 still
consumes A1 #1100 and therefore remains only a finite-X main-pulse-prefix
leading+oscillatory composite; A3 #1109 advances the older current-I4 identity
from nonlinear mean to compact radial stress; A4 #1110 independently audits
only the exact #1108 composite divergence. No scoped evidence is transferred
across those identities and no scoped audit is promoted to Navier--Stokes
validation.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import re
from pathlib import Path
from typing import Any, Mapping

from . import kokuno_a5_current_i4_mean_audit_pulse_routing as parent_a5

SCHEMA = "kokuno-a5-xi11-pulse-i4-stress-routing-v1"
TASK = "KOKUNO-A5-XI11-PULSE-I4-STRESS-ROUTING-112"
_HEX40 = re.compile(r"^[0-9a-f]{40}$")
_HEX64 = re.compile(r"^[0-9a-f]{64}$")

PARENT_A5 = {
    "pr": 1103,
    "head": "68535271cdbb79cead9431109b5311d6047cfd0e",
    "source_path": "src/openai_ns_reconstruction/kokuno_a5_current_i4_mean_audit_pulse_routing.py",
    "source_blob": "07f7818984d42d7249eee100bbd485c9791d40a1",
    "task": "KOKUNO-A5-CURRENT-I4-MEAN-AUDIT-PULSE-ROUTING-111",
    "schema": "kokuno-a5-current-i4-mean-audit-pulse-routing-v1",
}

AGENT1_XI11_MAIN_PULSE = {
    "pr": 1107,
    "head": "45da043dd2b4cd067f005a72c7e21fd0f2bcf309",
    "source_path": "src/openai_ns_reconstruction/kokuno_pa16_current_cartesian_main_pulse_logx.py",
    "source_blob": "477bda899694e4c75786d59850a2eddada669f07",
    "parent_pr": 1100,
    "parent_head": "2a6c58ff6bf6d1ba6b875bca0e6abf9d1118e156",
    "candidate_stage": "current-main-pulse-leading-only-logX-through-xi11",
    "overflow_safe_logx_similarity_materialized": True,
    "full_source_xi_11_current_cartesian_materialized": True,
    "identity_preserving_save_load_available": True,
    "repository_autonomous_principal_amplitude_consumed": True,
    "source_exact_amplitude_root_materialized": False,
    "source_pulse_end_mj_corrections_materialized": False,
    "terminal_global_leading_velocity_materialized": False,
    "matching_agent2_composite_consuming_this_exact_identity_materialized": False,
    "exact_head_tests_run": 35680025231,
    "exact_head_dedicated_run": 35680025259,
    "actions_status_at_registration": "queued",
}

AGENT2_MAIN_PULSE_PREFIX_COMPOSITE = {
    "pr": 1108,
    "head": "9f221dd57ef4b5e2d2e80e24c8cf531e991a5387",
    "source_path": "src/openai_ns_reconstruction/kokuno_current_main_pulse_leading_oscillatory_identity.py",
    "source_blob": "d8bde2336dc4a8f8dd7cf0ac6fecd184ea693390",
    "parent_pr": 1099,
    "parent_head": "3f0988ff016f41126ec546681cba696c5fa8001d",
    "consumed_agent1_pr": 1100,
    "consumed_agent1_head": "2a6c58ff6bf6d1ba6b875bca0e6abf9d1118e156",
    "newer_agent1_pr": 1107,
    "newer_agent1_head": "45da043dd2b4cd067f005a72c7e21fd0f2bcf309",
    "newer_agent1_1107_consumed": False,
    "candidate_stage": "current-main-pulse-finite-X-prefix",
    "self_contained_project_cartesian_velocity": True,
    "leading_plus_frozen_complete_curl_oscillation_materialized": True,
    "identity_preserving_save_load_available": True,
    "full_source_xi_11_current_cartesian_materialized": False,
    "terminal_global_leading_velocity_materialized": False,
    "complete_ns_residual_assessed": False,
    "exact_head_tests_run": 35680128360,
    "exact_head_dedicated_run": 35680128408,
    "actions_status_at_registration": "queued",
}

AGENT3_CURRENT_I4_RADIAL_STRESS = {
    "pr": 1109,
    "head": "49590ff311fef1dec4bd850b3013fd989485c87a",
    "source_path": "src/openai_ns_reconstruction/kokuno_current_i4_nonlinear_radial_stress.py",
    "source_blob": "b8d9928934d4df74b52e8cd0d2c27c468257cdc9",
    "parent_pr": 1101,
    "parent_head": "6c9a4e7b0dec4b3669fb1dc60297b18f10163e0b",
    "agent2_composite_pr": 1080,
    "agent2_composite_head": "c40d8ddecd2971544a6e07dab093436b423cf326",
    "agent1_leading_pr": 1079,
    "agent1_leading_head": "b06742ca6e189499192ede3cce40f62cdc1e35ca",
    "candidate_stage": "current-I4",
    "nonlinear_mean_attribution_consumed": True,
    "compact_theta_e2_stress_materialized": True,
    "compact_axial_e1_stress_materialized": True,
    "current_i4_compact_radial_stress_materialized": True,
    "source_i4_five_row_mean_correction_materialized": False,
    "radial_force_partial_z_sigma_1_materialized": False,
    "matching_agent4_stress_audit_present": False,
    "authorized_as_complete_ns_correction_target": False,
    "exact_head_tests_run": 35680519473,
    "exact_head_dedicated_run": 35680519595,
    "actions_status_at_registration": "queued",
}

AGENT4_MAIN_PULSE_PREFIX_DIVERGENCE_AUDIT = {
    "pr": 1110,
    "head": "8a9fb9ed1a67ec5ad7d9d856be71ee23467654bf",
    "source_path": "src/openai_ns_reconstruction/kokuno_a4_current_main_pulse_composite_divergence_independent_audit.py",
    "source_blob": "6dce05c1306a855b0ccf52b3415d7b7774279d09",
    "audited_agent2_pr": 1108,
    "audited_agent2_head": "9f221dd57ef4b5e2d2e80e24c8cf531e991a5387",
    "audited_agent1_pr": 1100,
    "audited_agent1_head": "2a6c58ff6bf6d1ba6b875bca0e6abf9d1118e156",
    "newer_agent1_1107_consumed": False,
    "seed": 9173861,
    "spatial_steps": [0.02, 0.01, 0.005],
    "operator": "public-Cartesian centered-FD2 after real candidate save/reload",
    "implementation_distinct": True,
    "registered": True,
    "scientific_admission": False,
    "complete_ns_residual_evidence": False,
    "exact_head_tests_run": 35680904633,
    "exact_head_dedicated_run": 35680904602,
    "actions_status_at_registration": "queued",
}

FINAL_GATES = copy.deepcopy(parent_a5.FINAL_GATES)
ST006_BASELINE = copy.deepcopy(parent_a5.ST006_BASELINE)
READINESS = copy.deepcopy(parent_a5.READINESS)

TRUTH_BOUNDARY = {
    "parent_a5_checkpoint_preserved": True,
    "leading_frontier_stage": "current-main-pulse-leading-only-logX-through-xi11",
    "latest_self_contained_composite_frontier_stage": "current-main-pulse-finite-X-prefix",
    "correction_frontier_stage": "current-I4-compact-radial-stress",
    "agent1_xi11_logx_main_pulse_leading_materialized": True,
    "agent1_xi11_logx_save_load_available": True,
    "agent1_xi11_uses_repository_autonomous_amplitude_proxy": True,
    "source_exact_main_pulse_amplitude_materialized": False,
    "source_pulse_end_mj_corrections_materialized": False,
    "terminal_global_leading_completion_materialized": False,
    "main_pulse_prefix_agent2_project_composite_materialized": True,
    "main_pulse_prefix_agent2_composite_consumes_agent1_1100": True,
    "main_pulse_prefix_agent2_composite_consumes_newer_agent1_1107": False,
    "matching_xi11_agent2_project_composite_materialized": False,
    "agent4_matching_main_pulse_prefix_divergence_audit_present": True,
    "agent4_matching_main_pulse_prefix_divergence_audit_registered": True,
    "agent4_matching_main_pulse_prefix_divergence_audit_admitted": False,
    "agent4_matching_main_pulse_prefix_divergence_audit_is_complete_ns_evidence": False,
    "agent4_main_pulse_prefix_audit_transferred_to_xi11_identity": False,
    "current_i4_nonlinear_mean_attribution_preserved": True,
    "current_i4_compact_radial_stress_materialized": True,
    "current_i4_radial_force_materialized": False,
    "agent4_matching_current_i4_radial_stress_audit_present": False,
    "agent4_matching_current_i4_radial_stress_audit_registered": False,
    "source_i4_five_row_mean_correction_materialized": False,
    "current_i4_stress_authorized_as_complete_ns_correction_target": False,
    "current_i4_stress_evidence_transferred_to_main_pulse_prefix": False,
    "current_i4_stress_evidence_transferred_to_xi11_pulse": False,
    "main_pulse_prefix_divergence_evidence_transferred_to_current_i4_stress": False,
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
        "A1 #1107 now evaluates the current main-pulse leading-only field through "
        "the public principal endpoint xi=11 in an overflow-safe log-X/log-F chart; "
        "it still lacks source-exact amplitude, M/J end closure and terminal/global completion"
    ),
    "latest_self_contained_composite_frontier": (
        "A2 #1108 is a save/loadable project-domain leading+frozen-complete-curl "
        "composite, but it consumes exact A1 #1100 and therefore remains a finite-X "
        "main-pulse-prefix identity rather than the newer #1107 xi=11 identity"
    ),
    "composite_validator_frontier": (
        "A4 #1110 independently audits only the exact A2 #1108 finite-X composite "
        "with public Cartesian velocity and centered FD2; it is scoped divergence evidence only"
    ),
    "correction_frontier": (
        "A3 #1109 advances the distinct older current-I4 A2 #1080 identity from "
        "nonlinear mean to compact theta/axial radial stress; partial_z sigma_1 and "
        "a matching A4 stress audit remain absent"
    ),
    "identity_firewall": (
        "A1 #1107, A2 #1108/A4 #1110, and A3 #1109 are distinct semantic identities; "
        "none of their scoped evidence is transferred across those identity boundaries"
    ),
    "next_shortest_blocker": (
        "recompose A2 on exact A1 #1107 to obtain a self-contained xi=11 composite "
        "and matching A4 scoped audit; on the correction lane carry exact A3 #1109 "
        "axial stress through partial_z sigma_1 and obtain matching A4 stress/force audits. "
        "Global velocity, matched pressure and preregistered restricted forcing are still "
        "mandatory before a complete NS defect or finite correction cycle can be admitted"
    ),
}


def _canonical_json(payload: Mapping[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _sha256(payload: Mapping[str, Any]) -> str:
    return hashlib.sha256(_canonical_json(payload).encode()).hexdigest()


def _hex40(value: Any, label: str) -> None:
    if not isinstance(value, str) or _HEX40.fullmatch(value) is None:
        raise ValueError(f"{label} is not an exact git identity")


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
        "agent1_xi11_main_pulse": copy.deepcopy(AGENT1_XI11_MAIN_PULSE),
        "agent2_main_pulse_prefix_composite": copy.deepcopy(AGENT2_MAIN_PULSE_PREFIX_COMPOSITE),
        "agent3_current_i4_radial_stress": copy.deepcopy(AGENT3_CURRENT_I4_RADIAL_STRESS),
        "agent4_main_pulse_prefix_divergence_audit": copy.deepcopy(AGENT4_MAIN_PULSE_PREFIX_DIVERGENCE_AUDIT),
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

    expected = {
        "parent_a5": PARENT_A5,
        "agent1_xi11_main_pulse": AGENT1_XI11_MAIN_PULSE,
        "agent2_main_pulse_prefix_composite": AGENT2_MAIN_PULSE_PREFIX_COMPOSITE,
        "agent3_current_i4_radial_stress": AGENT3_CURRENT_I4_RADIAL_STRESS,
        "agent4_main_pulse_prefix_divergence_audit": AGENT4_MAIN_PULSE_PREFIX_DIVERGENCE_AUDIT,
    }
    for name, frozen in expected.items():
        block = registration.get(name)
        if not isinstance(block, Mapping):
            raise ValueError(f"missing {name}")
        _hex40(block.get("head"), f"{name}.head")
        _hex40(block.get("source_blob"), f"{name}.source_blob")
        if dict(block) != frozen:
            raise ValueError(f"{name} drifted")

    if registration.get("final_gates") != FINAL_GATES or registration.get("st006_baseline") != ST006_BASELINE:
        raise ValueError("fixed science gates or ST006 baseline drifted")
    if registration.get("readiness") != READINESS or registration.get("truth_boundary") != TRUTH_BOUNDARY:
        raise ValueError("readiness/truth boundary drifted")
    if registration.get("pipeline_position") != PIPELINE_POSITION:
        raise ValueError("pipeline position drifted")

    a1 = registration["agent1_xi11_main_pulse"]
    a2 = registration["agent2_main_pulse_prefix_composite"]
    a3 = registration["agent3_current_i4_radial_stress"]
    a4 = registration["agent4_main_pulse_prefix_divergence_audit"]

    if a1["full_source_xi_11_current_cartesian_materialized"] is not True:
        raise ValueError("A1 xi=11 materialization disappeared")
    if a1["source_exact_amplitude_root_materialized"] is not False:
        raise ValueError("autonomous pulse amplitude laundered as source exact")
    if a1["source_pulse_end_mj_corrections_materialized"] is not False:
        raise ValueError("missing M/J end closure was promoted")
    if a1["terminal_global_leading_velocity_materialized"] is not False:
        raise ValueError("xi=11 leading promoted to terminal/global")

    if a2["consumed_agent1_head"] != a1["parent_head"]:
        raise ValueError("A2 #1108 is not bound to exact A1 #1100")
    if a2["newer_agent1_head"] != a1["head"]:
        raise ValueError("A2 newer-A1 firewall identity drifted")
    if a2["newer_agent1_1107_consumed"] is not False:
        raise ValueError("A2 #1108 laundered onto newer A1 #1107")
    if a2["full_source_xi_11_current_cartesian_materialized"] is not False:
        raise ValueError("finite-X composite promoted to xi=11")

    if a4["audited_agent2_head"] != a2["head"]:
        raise ValueError("A4 #1110 detached from exact A2 #1108")
    if a4["audited_agent1_head"] != a2["consumed_agent1_head"]:
        raise ValueError("A4 #1110 audited a different A1 identity")
    if a4["newer_agent1_1107_consumed"] is not False:
        raise ValueError("A4 #1110 evidence transplanted to A1 #1107")
    if a4["implementation_distinct"] is not True or a4["registered"] is not True:
        raise ValueError("A4 main-pulse divergence audit not registered")
    if a4["scientific_admission"] is not False or a4["complete_ns_residual_evidence"] is not False:
        raise ValueError("A4 scoped divergence audit promoted to PDE evidence")

    if a3["parent_head"] != parent_a5.AGENT3_CURRENT_I4_NONLINEAR_MEAN["head"]:
        raise ValueError("A3 #1109 detached from exact A3 #1101")
    if a3["current_i4_compact_radial_stress_materialized"] is not True:
        raise ValueError("current-I4 radial stress disappeared")
    if a3["radial_force_partial_z_sigma_1_materialized"] is not False:
        raise ValueError("A3 stress promoted to radial force")
    if a3["matching_agent4_stress_audit_present"] is not False:
        raise ValueError("missing A4 stress audit promoted")
    if a3["authorized_as_complete_ns_correction_target"] is not False:
        raise ValueError("scoped current-I4 stress prematurely authorized")

    if registration["readiness"] != {
        "leading_ready": False,
        "oscillatory_ready": True,
        "correction_ready": False,
        "velocity_export_ready": False,
        "pde_validated": False,
    }:
        raise ValueError("core readiness promoted")

    if FINAL_GATES["normalized_momentum_sampled_max"] != 1e-3:
        raise ValueError("momentum sampled-max gate drifted")
    if FINAL_GATES["normalized_momentum_volume_l2"] != 1e-3:
        raise ValueError("momentum volume-L2 gate drifted")
    if FINAL_GATES["normalized_divergence_sampled_max"] != 1e-5:
        raise ValueError("divergence sampled-max gate drifted")
    if FINAL_GATES["normalized_divergence_volume_l2"] != 1e-5:
        raise ValueError("divergence volume-L2 gate drifted")
    if FINAL_GATES["canonical_quadrature"] != [24, 48, 96]:
        raise ValueError("canonical quadrature drifted")
    if FINAL_GATES["residual_defined_free_forcing_allowed"] is not False:
        raise ValueError("residual-defined free forcing was enabled")

    hard_false = (
        "source_exact_main_pulse_amplitude_materialized",
        "source_pulse_end_mj_corrections_materialized",
        "terminal_global_leading_completion_materialized",
        "main_pulse_prefix_agent2_composite_consumes_newer_agent1_1107",
        "matching_xi11_agent2_project_composite_materialized",
        "agent4_matching_main_pulse_prefix_divergence_audit_admitted",
        "agent4_matching_main_pulse_prefix_divergence_audit_is_complete_ns_evidence",
        "agent4_main_pulse_prefix_audit_transferred_to_xi11_identity",
        "current_i4_radial_force_materialized",
        "agent4_matching_current_i4_radial_stress_audit_present",
        "agent4_matching_current_i4_radial_stress_audit_registered",
        "source_i4_five_row_mean_correction_materialized",
        "current_i4_stress_authorized_as_complete_ns_correction_target",
        "current_i4_stress_evidence_transferred_to_main_pulse_prefix",
        "current_i4_stress_evidence_transferred_to_xi11_pulse",
        "main_pulse_prefix_divergence_evidence_transferred_to_current_i4_stress",
        "source_positive_order_i3_correction_materialized",
        "global_leading_plus_oscillatory_velocity_materialized",
        "matched_cartesian_pressure_materialized",
        "cartesian_matched_pressure_gradient_materialized",
        "preregistered_restricted_forcing_materialized",
        "forcing_proven_not_residual_defined",
        "complete_identity_bound_ns_defect_materialized",
        "current_cartesian_correction_velocity_materialized",
        "finite_correction_cycle_admitted",
        "real_candidate_finite_correction_cycle_run",
        "heldout_normalized_full_ns_residual_assessed",
        "canonical_24_48_96_whole_domain_admission_assessed",
        "same_protocol_comparable_to_st006",
        "scientific_admission",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    )
    for key in hard_false:
        if registration["truth_boundary"][key] is not False:
            raise ValueError(f"truth state promoted at {key}")

    digest = registration.get("digest")
    if not isinstance(digest, str) or _HEX64.fullmatch(digest) is None:
        raise ValueError("invalid registration digest")
    payload = {
        key: copy.deepcopy(value)
        for key, value in registration.items()
        if key not in {"schema", "task", "digest"}
    }
    if digest != _sha256(payload):
        raise ValueError("registration digest mismatch")


def save_registration(path: str | Path) -> Path:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(build_registration(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return destination


def load_registration(path: str | Path) -> dict[str, Any]:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("registration must be a JSON object")
    validate_registration(payload)
    return payload


def _main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    output = save_registration(args.output)
    receipt = load_registration(output)
    print(json.dumps({
        "path": str(output),
        "schema": receipt["schema"],
        "task": receipt["task"],
        "digest": receipt["digest"],
    }, sort_keys=True))


if __name__ == "__main__":
    _main()
