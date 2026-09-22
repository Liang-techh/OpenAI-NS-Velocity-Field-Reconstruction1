"""Fail-closed Agent-5 routing checkpoint for current-I4 mean evidence.

This module is integration/provenance glue only.  It advances the correction
frontier to the exact current-I4 nonlinear mean plus a matching
implementation-distinct Agent-4 audit, while recording newer pulse-leading and
source-shell Agent-1/2 work as non-consumed siblings.  It does not transfer
older current-I2 stress/force evidence across candidate identities and does not
promote any scoped audit to complete Navier--Stokes validation.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import re
from pathlib import Path
from typing import Any, Mapping

from . import kokuno_a5_post_i4_routing_i2_force_audit_ingest as parent_a5

SCHEMA = "kokuno-a5-current-i4-mean-audit-pulse-routing-v1"
TASK = "KOKUNO-A5-CURRENT-I4-MEAN-AUDIT-PULSE-ROUTING-111"
_HEX40 = re.compile(r"^[0-9a-f]{40}$")
_HEX64 = re.compile(r"^[0-9a-f]{64}$")

PARENT_A5 = {
    "pr": 1093,
    "head": "666d5dad83dee8a9313c3a59bf102b10bb2917cc",
    "source_path": "src/openai_ns_reconstruction/kokuno_a5_post_i4_routing_i2_force_audit_ingest.py",
    "source_blob": "59be383c2f266ed83cc5518ce4b43746f52161e1",
    "task": "KOKUNO-A5-POST-I4-ROUTING-I2-FORCE-AUDIT-INGEST-110",
    "schema": "kokuno-a5-post-i4-routing-i2-force-audit-ingest-v1",
}

AGENT1_CURRENT_MAIN_PULSE = {
    "pr": 1100,
    "head": "2a6c58ff6bf6d1ba6b875bca0e6abf9d1118e156",
    "source_path": "src/openai_ns_reconstruction/kokuno_pa16_current_cartesian_main_pulse.py",
    "source_blob": "efcce36612ffc1cb3a74f5eb9c9ae6c86748f60d",
    "schema": "kokuno-pa16-current-cartesian-main-pulse-v1",
    "parent_agent1_pr": 1094,
    "parent_agent1_head": "dc33914d3dd60c1dc4c8dfaab0760ef071385d54",
    "candidate_stage": "active-main-pulse-leading-only-finite-X-prefix",
    "current_cartesian_main_pulse_prefix_materialized": True,
    "identity_preserving_save_load_available": True,
    "repository_autonomous_principal_amplitude_consumed": True,
    "source_exact_amplitude_root_materialized": False,
    "source_pulse_end_mj_corrections_materialized": False,
    "full_source_xi_11_current_cartesian_materialized": False,
    "terminal_global_leading_velocity_materialized": False,
    "matching_agent2_project_cartesian_composite_materialized": False,
    "exact_head_tests_run": 35676517652,
    "exact_head_dedicated_run": 35676517667,
    "actions_status_at_registration": "queued",
}

AGENT2_SOURCE_SHELL_EDGE = {
    "pr": 1099,
    "head": "3f0988ff016f41126ec546681cba696c5fa8001d",
    "source_path": "src/openai_ns_reconstruction/kokuno_source_shell_edge_localized_harmonic.py",
    "source_blob": "a16fba5d2c77fe952d39ae0dae9be2c04b2031fd",
    "task": "K2-OSC-099",
    "parent_agent2_pr": 1089,
    "parent_agent2_head": "4bc9408d10b2f3c636abd4f6ce33cefab391093d",
    "source_shell_edge_sqrt_zeta_materialized": True,
    "analytic_shell_edge_x_derivative_materialized": True,
    "complete_curl_localization_retains_support_gradient_terms": True,
    "caller_supplies_source_partition_and_pulse_cutoff": True,
    "caller_supplies_corrected_background_and_forcing": True,
    "project_domain_coordinate_map_applied": False,
    "self_contained_velocity_xyzt_provider": False,
    "consumed_by_agent1_1100_identity": False,
    "exact_head_tests_run": 35676502166,
    "exact_head_dedicated_run": 35676502159,
    "actions_status_at_registration": "queued",
}

AGENT3_CURRENT_I4_NONLINEAR_MEAN = {
    "pr": 1101,
    "head": "6c9a4e7b0dec4b3669fb1dc60297b18f10163e0b",
    "source_path": "src/openai_ns_reconstruction/kokuno_current_i4_nonlinear_mean_attribution.py",
    "source_blob": "10d067570ddba834db3f19c9859c2efdad16b3e8",
    "task": "KOKUNO-A3-CURRENT-I4-NONLINEAR-MEAN-119",
    "schema": "kokuno-a3-current-i4-nonlinear-mean-v1",
    "parent_agent3_pr": 1090,
    "parent_agent3_head": "5e8512392df683375fabd132f35577e1b94773ff",
    "agent2_composite_pr": 1080,
    "agent2_composite_head": "c40d8ddecd2971544a6e07dab093436b423cf326",
    "agent1_leading_pr": 1079,
    "agent1_leading_head": "b06742ca6e189499192ede3cce40f62cdc1e35ca",
    "candidate_stage": "current-I4",
    "nonlinear_mean_attribution_materialized": True,
    "source_i4_five_row_mean_correction_materialized": False,
    "current_i4_compact_radial_stress_materialized": False,
    "current_i4_radial_force_materialized": False,
    "authorized_as_complete_ns_correction_target": False,
    "exact_head_tests_run": 35676609539,
    "exact_head_dedicated_run": 35676609594,
    "actions_status_at_registration": "queued",
}

AGENT4_CURRENT_I4_NONLINEAR_MEAN_AUDIT = {
    "pr": 1102,
    "head": "d3999c0e4e0e3483d4199f009e936060c42848c8",
    "source_path": "src/openai_ns_reconstruction/kokuno_a4_current_i4_nonlinear_mean_independent_audit.py",
    "source_blob": "5f50c05f5fc106c0a9e9a59be1b0f2ff687c9b9c",
    "task": "K4-VAL-113",
    "schema": "kokuno-a4-current-i4-nonlinear-mean-independent-audit-v1",
    "audited_agent3_pr": 1101,
    "audited_agent3_head": "6c9a4e7b0dec4b3669fb1dc60297b18f10163e0b",
    "agent2_composite_pr": 1080,
    "agent2_composite_head": "c40d8ddecd2971544a6e07dab093436b423cf326",
    "operator": "public-Cartesian centered five-point FD4 plus nonuniform Gauss-Legendre cylindrical m=0 projection",
    "seed": 9173851,
    "spatial_steps": [0.012, 0.006, 0.003],
    "angular_orders": [20, 40, 80],
    "implementation_distinct": True,
    "registered": True,
    "scientific_admission": False,
    "complete_ns_residual_evidence": False,
    "exact_head_tests_run": 35677321703,
    "exact_head_dedicated_run": 35677321707,
    "actions_status_at_registration": "queued",
}

FINAL_GATES = copy.deepcopy(parent_a5.FINAL_GATES)
ST006_BASELINE = copy.deepcopy(parent_a5.ST006_BASELINE)
READINESS = copy.deepcopy(parent_a5.READINESS)

TRUTH_BOUNDARY = {
    "parent_a5_checkpoint_preserved": True,
    "leading_frontier_stage": "active-main-pulse-leading-only-finite-X-prefix",
    "latest_self_contained_composite_frontier_stage": "current-I4",
    "correction_frontier_stage": "current-I4-nonlinear-mean",
    "current_main_pulse_leading_prefix_materialized": True,
    "current_main_pulse_leading_save_load_available": True,
    "current_main_pulse_uses_repository_autonomous_amplitude_proxy": True,
    "source_exact_main_pulse_amplitude_materialized": False,
    "full_xi_11_current_cartesian_pulse_materialized": False,
    "terminal_global_leading_completion_materialized": False,
    "matching_main_pulse_agent2_project_cartesian_composite_materialized": False,
    "agent2_source_shell_edge_adapter_present": True,
    "agent2_source_shell_edge_adapter_self_contained_project_velocity": False,
    "agent2_source_shell_edge_consumed_by_main_pulse_candidate": False,
    "current_i4_self_contained_leading_plus_oscillatory_composite_preserved": True,
    "current_i4_nonlinear_mean_attribution_materialized": True,
    "agent4_matching_current_i4_nonlinear_mean_audit_present": True,
    "agent4_matching_current_i4_nonlinear_mean_audit_registered": True,
    "agent4_matching_current_i4_nonlinear_mean_audit_admitted": False,
    "agent4_matching_current_i4_nonlinear_mean_audit_is_complete_ns_evidence": False,
    "source_positive_order_i3_correction_materialized": False,
    "source_i4_five_row_mean_correction_materialized": False,
    "current_i4_compact_radial_stress_materialized": False,
    "current_i4_radial_force_materialized": False,
    "current_i4_mean_authorized_as_complete_ns_correction_target": False,
    "older_current_i2_stress_force_evidence_preserved": True,
    "older_current_i2_stress_force_evidence_transferred_to_current_i4": False,
    "older_current_i2_stress_force_evidence_transferred_to_main_pulse": False,
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
        "A1 #1100 extends the exact current leading-only Cartesian field into a "
        "finite-X main-pulse prefix; it is not terminal/global and has no matching "
        "self-contained A2 project-domain composite"
    ),
    "latest_self_contained_composite_frontier": (
        "A1 #1079 -> A2 #1080 remains the latest self-contained current-I4 "
        "leading+frozen-complete-curl Cartesian candidate"
    ),
    "source_oscillatory_frontier": (
        "A2 #1099 materializes the corrected source sqrt(zeta) shell edge inside "
        "complete-curl localization but remains caller-dependent/source-chart only"
    ),
    "correction_frontier": (
        "A3 #1101 now binds the exact current-I4 A2 #1080 identity to nonlinear "
        "m=0 mean attribution; A4 #1102 independently audits that mean and is "
        "registered but unresolved, while current-I4 stress/force remain absent"
    ),
    "stage_firewall": (
        "current-I2 stress/force/audit receipts remain historical scoped evidence "
        "and are not transferred onto current-I4 or main-pulse identities"
    ),
    "next_shortest_blocker": (
        "carry exact A3 #1101 current-I4 mean through compact radial stress and "
        "source-required partial_z sigma_1 radial force, then obtain matching A4 "
        "audits; in parallel finish the pulse/terminal-global leading route and "
        "build a matching self-contained A2 composite before pressure/forcing/full-defect work"
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
    if parent_a5.AGENT3_CURRENT_I2_RADIAL_FORCE["head"] != AGENT3_CURRENT_I4_NONLINEAR_MEAN["parent_agent3_head"]:
        raise ValueError("A3 #1101 detached from the exact A3 #1090 parent")


def build_registration() -> dict[str, Any]:
    _validate_parent()
    payload = {
        "parent_a5": copy.deepcopy(PARENT_A5),
        "agent1_current_main_pulse": copy.deepcopy(AGENT1_CURRENT_MAIN_PULSE),
        "agent2_source_shell_edge": copy.deepcopy(AGENT2_SOURCE_SHELL_EDGE),
        "agent3_current_i4_nonlinear_mean": copy.deepcopy(AGENT3_CURRENT_I4_NONLINEAR_MEAN),
        "agent4_current_i4_nonlinear_mean_audit": copy.deepcopy(AGENT4_CURRENT_I4_NONLINEAR_MEAN_AUDIT),
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
        "agent1_current_main_pulse": AGENT1_CURRENT_MAIN_PULSE,
        "agent2_source_shell_edge": AGENT2_SOURCE_SHELL_EDGE,
        "agent3_current_i4_nonlinear_mean": AGENT3_CURRENT_I4_NONLINEAR_MEAN,
        "agent4_current_i4_nonlinear_mean_audit": AGENT4_CURRENT_I4_NONLINEAR_MEAN_AUDIT,
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

    a1 = registration["agent1_current_main_pulse"]
    a2 = registration["agent2_source_shell_edge"]
    a3 = registration["agent3_current_i4_nonlinear_mean"]
    a4 = registration["agent4_current_i4_nonlinear_mean_audit"]

    if a1["source_exact_amplitude_root_materialized"] is not False:
        raise ValueError("autonomous pulse amplitude laundered as source exact")
    if a1["terminal_global_leading_velocity_materialized"] is not False or a1["matching_agent2_project_cartesian_composite_materialized"] is not False:
        raise ValueError("main-pulse leading promoted to global/composite")
    if a2["self_contained_velocity_xyzt_provider"] is not False or a2["project_domain_coordinate_map_applied"] is not False:
        raise ValueError("source shell adapter laundered into project velocity")
    if a2["consumed_by_agent1_1100_identity"] is not False:
        raise ValueError("unconsumed A2 sibling was attached to A1 #1100")
    if a3["agent2_composite_head"] != a4["agent2_composite_head"]:
        raise ValueError("A3/A4 current-I4 candidate identity mismatch")
    if a4["audited_agent3_head"] != a3["head"]:
        raise ValueError("A4 audit detached from exact A3 #1101")
    if a4["implementation_distinct"] is not True or a4["registered"] is not True:
        raise ValueError("A4 independent audit not registered")
    if a4["scientific_admission"] is not False or a4["complete_ns_residual_evidence"] is not False:
        raise ValueError("A4 scoped mean audit promoted to PDE evidence")
    if a3["current_i4_compact_radial_stress_materialized"] is not False or a3["current_i4_radial_force_materialized"] is not False:
        raise ValueError("current-I4 correction frontier overstated")
    if a3["authorized_as_complete_ns_correction_target"] is not False:
        raise ValueError("current-I4 nonlinear mean prematurely authorized as NS correction target")

    if registration["readiness"] != {
        "leading_ready": False,
        "oscillatory_ready": True,
        "correction_ready": False,
        "velocity_export_ready": False,
        "pde_validated": False,
    }:
        raise ValueError("core readiness promoted")
    if FINAL_GATES["normalized_momentum_sampled_max"] != 1e-3 or FINAL_GATES["normalized_momentum_volume_l2"] != 1e-3:
        raise ValueError("momentum gate drifted")
    if FINAL_GATES["normalized_divergence_sampled_max"] != 1e-5 or FINAL_GATES["normalized_divergence_volume_l2"] != 1e-5:
        raise ValueError("divergence gate drifted")
    if FINAL_GATES["canonical_quadrature"] != [24, 48, 96]:
        raise ValueError("canonical quadrature drifted")
    if FINAL_GATES["residual_defined_free_forcing_allowed"] is not False:
        raise ValueError("residual-defined free forcing was enabled")

    hard_false = (
        "source_exact_main_pulse_amplitude_materialized",
        "full_xi_11_current_cartesian_pulse_materialized",
        "terminal_global_leading_completion_materialized",
        "matching_main_pulse_agent2_project_cartesian_composite_materialized",
        "agent2_source_shell_edge_adapter_self_contained_project_velocity",
        "agent2_source_shell_edge_consumed_by_main_pulse_candidate",
        "agent4_matching_current_i4_nonlinear_mean_audit_admitted",
        "agent4_matching_current_i4_nonlinear_mean_audit_is_complete_ns_evidence",
        "source_positive_order_i3_correction_materialized",
        "source_i4_five_row_mean_correction_materialized",
        "current_i4_compact_radial_stress_materialized",
        "current_i4_radial_force_materialized",
        "current_i4_mean_authorized_as_complete_ns_correction_target",
        "older_current_i2_stress_force_evidence_transferred_to_current_i4",
        "older_current_i2_stress_force_evidence_transferred_to_main_pulse",
        "global_leading_plus_oscillatory_velocity_materialized",
        "matched_cartesian_pressure_materialized",
        "preregistered_restricted_forcing_materialized",
        "complete_identity_bound_ns_defect_materialized",
        "current_cartesian_correction_velocity_materialized",
        "finite_correction_cycle_admitted",
        "heldout_normalized_full_ns_residual_assessed",
        "scientific_admission",
        "pde_validated",
    )
    for key in hard_false:
        if registration["truth_boundary"][key] is not False:
            raise ValueError(f"truth state promoted at {key}")

    digest = registration.get("digest")
    if not isinstance(digest, str) or _HEX64.fullmatch(digest) is None:
        raise ValueError("invalid registration digest")
    payload = {
        k: copy.deepcopy(v)
        for k, v in registration.items()
        if k not in {"schema", "task", "digest"}
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
    print(
        json.dumps(
            {
                "path": str(output),
                "schema": receipt["schema"],
                "task": receipt["task"],
                "digest": receipt["digest"],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    _main()
