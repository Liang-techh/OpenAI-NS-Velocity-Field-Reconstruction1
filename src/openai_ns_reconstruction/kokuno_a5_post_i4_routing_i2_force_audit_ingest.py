"""Fail-closed Agent-5 routing checkpoint after current I4.

Integration only: preserve the latest self-contained current-I4 composite,
register newer post-I4 leading-only/source-harmonic siblings, and advance the
older current-I2 correction route from stress to radial force plus matching
implementation-distinct A4 audit.  No evidence is transferred across candidate
identities and no full-NS/PDE state is promoted.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import re
from pathlib import Path
from typing import Any, Mapping

from . import kokuno_a5_current_i4_routing_i2_stress_ingest as parent_a5

SCHEMA = "kokuno-a5-post-i4-routing-i2-force-audit-ingest-v1"
TASK = "KOKUNO-A5-POST-I4-ROUTING-I2-FORCE-AUDIT-INGEST-110"
_HEX40 = re.compile(r"^[0-9a-f]{40}$")
_HEX64 = re.compile(r"^[0-9a-f]{64}$")

PARENT_A5 = {
    "pr": 1083,
    "head": "2886d888cafa80a30e9b21ffee5011e079d944a0",
    "source_path": "src/openai_ns_reconstruction/kokuno_a5_current_i4_routing_i2_stress_ingest.py",
    "source_blob": "505fc333463db23030a4eee72da30755317982fb",
    "task": "KOKUNO-A5-CURRENT-I4-ROUTING-I2-STRESS-INGEST-109",
    "schema": "kokuno-a5-current-i4-routing-i2-stress-ingest-v1",
}
AGENT1_POST_I4_LEADING = {
    "pr": 1088,
    "head": "bd85ed48323feb3b33b9c404a078630245dfe783",
    "source_path": "src/openai_ns_reconstruction/kokuno_pa16_current_cartesian_post_i4_pulse_entry_bridge.py",
    "source_blob": "6826085ee248a47f445d91031219bc2e5b655fd3",
    "schema": "kokuno-pa16-current-cartesian-post-i4-pulse-entry-bridge-v1",
    "parent_agent1_pr": 1079,
    "parent_agent1_head": "b06742ca6e189499192ede3cce40f62cdc1e35ca",
    "candidate_stage": "post-I4/pre-pulse-leading-only",
    "current_leading_to_pulse_entry_materialized": True,
    "identity_preserving_save_load_available": True,
    "matching_agent2_project_cartesian_composite_materialized": False,
    "source_axial_pulse_materialized": False,
    "source_positive_order_i3_correction_materialized": False,
    "source_i4_mean_correction_materialized": False,
    "outer_global_leading_velocity_materialized": False,
}
AGENT2_SOURCE_LOCALIZED_HARMONIC = {
    "pr": 1089,
    "head": "4bc9408d10b2f3c636abd4f6ce33cefab391093d",
    "source_path": "src/openai_ns_reconstruction/kokuno_source_zero_data_localized_harmonic.py",
    "source_blob": "8087e619bb1812df6dcb616e1a51f1ab4dd728d4",
    "task": "K2-OSC-098",
    "schema": "kokuno-source-zero-data-localized-harmonic-v1",
    "parent_agent2_pr": 1080,
    "parent_agent2_head": "c40d8ddecd2971544a6e07dab093436b423cf326",
    "potential_level_product_rule_localization_executable": True,
    "complete_curl_retains_support_gradient_terms": True,
    "real_pair_velocity_materialized_in_source_chart": True,
    "caller_supplies_corrected_background": True,
    "caller_supplies_mode_forcing_directional_jet": True,
    "caller_supplies_source_partition_jet": True,
    "caller_supplies_pulse_cutoff_value": True,
    "project_domain_coordinate_map_applied": False,
    "self_contained_velocity_xyzt_provider": False,
    "consumed_by_agent1_1088_identity": False,
}
AGENT3_CURRENT_I2_RADIAL_FORCE = {
    "pr": 1090,
    "head": "5e8512392df683375fabd132f35577e1b94773ff",
    "source_path": "src/openai_ns_reconstruction/kokuno_current_i2_nonlinear_radial_force.py",
    "source_blob": "9259ecdf4aa3912fb521acdf7e62060272395e38",
    "task": "KOKUNO-A3-CURRENT-I2-NONLINEAR-RADIAL-FORCE-118",
    "schema": "kokuno-a3-current-i2-nonlinear-radial-force-v1",
    "parent_agent3_pr": 1081,
    "parent_agent3_head": "cb9ee25d0c466d42ec7ca30ee89cc10517d4b5c4",
    "candidate_stage": "current-I2",
    "source_formula": "(div T)_r = partial_z sigma_1",
    "production_operator": "centered physical-z FD2",
    "z_steps": [0.02, 0.01, 0.005],
    "radial_force_materialized": True,
    "authorized_as_complete_ns_correction_target": False,
    "consumed_by_current_i4_candidate": False,
    "consumed_by_post_i4_candidate": False,
}
AGENT4_CURRENT_I2_RADIAL_FORCE_AUDIT = {
    "pr": 1092,
    "head": "a01bc1f088b6b0d6556d81cf417c2b73bd7c99e8",
    "source_path": "src/openai_ns_reconstruction/kokuno_a4_current_i2_radial_force_independent_audit.py",
    "source_blob": "6c217186fb57a6fbe5c0e20ebc3734897d7f58ac",
    "task": "K4-VAL-112",
    "schema": "kokuno-a4-current-i2-radial-force-independent-audit-v1",
    "audited_agent3_pr": 1090,
    "audited_agent3_head": "5e8512392df683375fabd132f35577e1b94773ff",
    "stress_agent3_pr": 1081,
    "stress_agent3_head": "cb9ee25d0c466d42ec7ca30ee89cc10517d4b5c4",
    "operator": "centered five-point FD4 from JSON-round-tripped public stress receipts",
    "seed": 9173841,
    "z_centers": [-0.149, -0.061, 0.053],
    "z_steps": [0.02, 0.01, 0.005],
    "implementation_distinct": True,
    "registered": True,
    "scientific_admission": False,
    "complete_ns_residual_evidence": False,
    "applies_to_current_i4_or_post_i4_identity": False,
}

FINAL_GATES = copy.deepcopy(parent_a5.FINAL_GATES)
ST006_BASELINE = copy.deepcopy(parent_a5.ST006_BASELINE)
READINESS = copy.deepcopy(parent_a5.READINESS)

TRUTH_BOUNDARY = {
    "parent_a5_checkpoint_preserved": True,
    "leading_frontier_stage": "post-I4/pre-pulse-leading-only",
    "latest_self_contained_composite_frontier_stage": "current-I4",
    "correction_frontier_stage": "current-I2-radial-force",
    "post_i4_pre_pulse_leading_materialized": True,
    "post_i4_pre_pulse_leading_save_load_available": True,
    "matching_post_i4_agent2_project_cartesian_composite_materialized": False,
    "current_i4_self_contained_leading_plus_oscillatory_composite_preserved": True,
    "current_i4_matching_a4_scoped_divergence_audit_preserved": True,
    "agent2_source_localized_harmonic_adapter_present": True,
    "agent2_source_localized_harmonic_adapter_self_contained_project_velocity": False,
    "agent2_source_localized_harmonic_consumed_by_post_i4_candidate": False,
    "current_i2_compact_radial_stress_materialized": True,
    "current_i2_radial_force_materialized": True,
    "agent4_matching_current_i2_radial_force_audit_present": True,
    "agent4_matching_current_i2_radial_force_audit_registered": True,
    "agent4_matching_current_i2_radial_force_audit_admitted": False,
    "current_i2_radial_force_consumed_by_current_i4_candidate": False,
    "current_i2_radial_force_consumed_by_post_i4_candidate": False,
    "current_i2_radial_force_authorized_as_complete_ns_correction_target": False,
    "agent4_current_i2_force_audit_transferred_to_current_i4_or_post_i4": False,
    "source_positive_order_i3_correction_materialized": False,
    "source_i4_mean_correction_materialized": False,
    "source_axial_pulse_materialized": False,
    "terminal_global_leading_completion_materialized": False,
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
    "leading_frontier": "A1 #1088 extends exact current leading-only velocity from I4 exit to pulse entry; no matching A2 project-domain composite exists yet",
    "latest_self_contained_composite_frontier": "parent A5 #1083 preserves A1 #1079 -> A2 #1080 current-I4 Cartesian composite plus A4 #1082 scoped divergence audit",
    "source_oscillatory_frontier": "A2 #1089 localizes zero-data source harmonics at potential level but still requires caller source data and is not a self-contained project velocity",
    "correction_frontier": "A1 #1061 -> A2 #1071 -> A3 #1073 -> A3 #1081 -> A3 #1090 current-I2 mean -> stress -> radial force, with matching A4 #1092 independent audit registered",
    "stage_firewall": "current-I2 force/audit evidence is not consumed by current-I4 or post-I4 identities; A2 #1089 source-chart output is not promoted to the missing post-I4 project-domain composite",
    "next_shortest_blocker": "materialize source-role-correct axial pulse and terminal/global leading velocity, then a matching self-contained A2 composite; matched pressure plus preregistered restricted forcing must exist before a complete NS defect can authorize a real finite correction cycle",
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
    if parent_a5.AGENT1_CURRENT_I4_LEADING["head"] != AGENT1_POST_I4_LEADING["parent_agent1_head"]:
        raise ValueError("A1 #1088 detached from exact current-I4 parent")
    if parent_a5.AGENT2_CURRENT_I4_COMPOSITE["head"] != AGENT2_SOURCE_LOCALIZED_HARMONIC["parent_agent2_head"]:
        raise ValueError("A2 #1089 detached from exact current-I4 A2 parent")
    if parent_a5.AGENT3_CURRENT_I2_STRESS["head"] != AGENT3_CURRENT_I2_RADIAL_FORCE["parent_agent3_head"]:
        raise ValueError("A3 #1090 detached from exact current-I2 stress")
    if parent_a5.FINAL_GATES != FINAL_GATES or parent_a5.ST006_BASELINE != ST006_BASELINE or parent_a5.READINESS != READINESS:
        raise ValueError("frozen parent science/readiness drifted")


def build_registration() -> dict[str, Any]:
    _validate_parent()
    payload = {
        "parent_a5": copy.deepcopy(PARENT_A5),
        "agent1_post_i4_leading": copy.deepcopy(AGENT1_POST_I4_LEADING),
        "agent2_source_localized_harmonic": copy.deepcopy(AGENT2_SOURCE_LOCALIZED_HARMONIC),
        "agent3_current_i2_radial_force": copy.deepcopy(AGENT3_CURRENT_I2_RADIAL_FORCE),
        "agent4_current_i2_radial_force_audit": copy.deepcopy(AGENT4_CURRENT_I2_RADIAL_FORCE_AUDIT),
        "final_gates": copy.deepcopy(FINAL_GATES),
        "st006_baseline": copy.deepcopy(ST006_BASELINE),
        "readiness": copy.deepcopy(READINESS),
        "truth_boundary": copy.deepcopy(TRUTH_BOUNDARY),
        "pipeline_position": copy.deepcopy(PIPELINE_POSITION),
    }
    out = {"schema": SCHEMA, "task": TASK, **payload, "digest": _sha256(payload)}
    validate_registration(out)
    return out


def validate_registration(registration: Mapping[str, Any]) -> None:
    _validate_parent()
    if registration.get("schema") != SCHEMA or registration.get("task") != TASK:
        raise ValueError("A5 registration identity drifted")
    expected = {
        "parent_a5": PARENT_A5,
        "agent1_post_i4_leading": AGENT1_POST_I4_LEADING,
        "agent2_source_localized_harmonic": AGENT2_SOURCE_LOCALIZED_HARMONIC,
        "agent3_current_i2_radial_force": AGENT3_CURRENT_I2_RADIAL_FORCE,
        "agent4_current_i2_radial_force_audit": AGENT4_CURRENT_I2_RADIAL_FORCE_AUDIT,
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
        raise ValueError("fixed gates/ST006 drifted")
    if registration.get("readiness") != READINESS or registration.get("truth_boundary") != TRUTH_BOUNDARY:
        raise ValueError("readiness/truth boundary drifted")
    if registration.get("pipeline_position") != PIPELINE_POSITION:
        raise ValueError("pipeline position drifted")

    a1, a2 = registration["agent1_post_i4_leading"], registration["agent2_source_localized_harmonic"]
    a3, a4 = registration["agent3_current_i2_radial_force"], registration["agent4_current_i2_radial_force_audit"]
    if a1["parent_agent1_head"] != parent_a5.AGENT1_CURRENT_I4_LEADING["head"]:
        raise ValueError("post-I4 leading lineage mismatch")
    if a2["parent_agent2_head"] != parent_a5.AGENT2_CURRENT_I4_COMPOSITE["head"]:
        raise ValueError("source-harmonic lineage mismatch")
    if a3["parent_agent3_head"] != parent_a5.AGENT3_CURRENT_I2_STRESS["head"]:
        raise ValueError("current-I2 force lineage mismatch")
    if a4["audited_agent3_head"] != a3["head"] or a4["stress_agent3_head"] != a3["parent_agent3_head"]:
        raise ValueError("A4 independent audit lineage mismatch")
    if a1["matching_agent2_project_cartesian_composite_materialized"] is not False:
        raise ValueError("post-I4 matching composite invented")
    if a2["self_contained_velocity_xyzt_provider"] is not False or a2["consumed_by_agent1_1088_identity"] is not False:
        raise ValueError("source adapter laundered into project candidate")
    if a3["authorized_as_complete_ns_correction_target"] is not False or a3["consumed_by_current_i4_candidate"] is not False or a3["consumed_by_post_i4_candidate"] is not False:
        raise ValueError("current-I2 force evidence transplanted/promoted")
    if a4["implementation_distinct"] is not True or a4["registered"] is not True or a4["scientific_admission"] is not False or a4["complete_ns_residual_evidence"] is not False or a4["applies_to_current_i4_or_post_i4_identity"] is not False:
        raise ValueError("A4 scoped evidence was promoted or transferred")

    if registration["readiness"] != {"leading_ready": False, "oscillatory_ready": True, "correction_ready": False, "velocity_export_ready": False, "pde_validated": False}:
        raise ValueError("core readiness promoted")
    if FINAL_GATES["normalized_momentum_sampled_max"] != 1e-3 or FINAL_GATES["normalized_momentum_volume_l2"] != 1e-3:
        raise ValueError("momentum gate drifted")
    if FINAL_GATES["normalized_divergence_sampled_max"] != 1e-5 or FINAL_GATES["normalized_divergence_volume_l2"] != 1e-5:
        raise ValueError("divergence gate drifted")
    if FINAL_GATES["canonical_quadrature"] != [24, 48, 96] or FINAL_GATES["residual_defined_free_forcing_allowed"] is not False:
        raise ValueError("canonical quadrature/free-forcing rule drifted")

    for key in (
        "matching_post_i4_agent2_project_cartesian_composite_materialized",
        "agent2_source_localized_harmonic_adapter_self_contained_project_velocity",
        "agent4_matching_current_i2_radial_force_audit_admitted",
        "current_i2_radial_force_consumed_by_current_i4_candidate",
        "current_i2_radial_force_consumed_by_post_i4_candidate",
        "current_i2_radial_force_authorized_as_complete_ns_correction_target",
        "source_positive_order_i3_correction_materialized", "source_i4_mean_correction_materialized",
        "source_axial_pulse_materialized", "terminal_global_leading_completion_materialized",
        "global_leading_plus_oscillatory_velocity_materialized", "matched_cartesian_pressure_materialized",
        "preregistered_restricted_forcing_materialized", "complete_identity_bound_ns_defect_materialized",
        "current_cartesian_correction_velocity_materialized", "finite_correction_cycle_admitted",
        "heldout_normalized_full_ns_residual_assessed", "pde_validated",
    ):
        if registration["truth_boundary"][key] is not False:
            raise ValueError(f"truth state promoted at {key}")

    digest = registration.get("digest")
    if not isinstance(digest, str) or _HEX64.fullmatch(digest) is None:
        raise ValueError("invalid registration digest")
    payload = {k: copy.deepcopy(v) for k, v in registration.items() if k not in {"schema", "task", "digest"}}
    if digest != _sha256(payload):
        raise ValueError("registration digest mismatch")


def save_registration(path: str | Path) -> Path:
    dst = Path(path)
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text(json.dumps(build_registration(), indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return dst


def load_registration(path: str | Path) -> dict[str, Any]:
    obj = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(obj, dict):
        raise ValueError("registration must be a JSON object")
    validate_registration(obj)
    return obj


def _main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    path = save_registration(args.output)
    receipt = load_registration(path)
    print(json.dumps({"path": str(path), "schema": receipt["schema"], "task": receipt["task"], "digest": receipt["digest"]}, sort_keys=True))


if __name__ == "__main__":
    _main()
