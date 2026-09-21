"""Agent-5 integration contract for the current Cartesian leading-velocity seam.

This module is provenance/integration glue only. It registers the first
current-lineage candidate-side Cartesian leading ``velocity(x,y,z,t)`` surface
from Agent 1 PR #965 together with Agent 4 PR #966's implementation-distinct,
public-velocity-only divergence audit protocol. Agent 2 #960 and Agent 3 #961
remain exact sibling execution state.

Truth boundary: the Agent-1 velocity exists only through the currently
materialized ``X_h`` region. It is not a global compactly supported leading
field, does not include matched pressure or preregistered restricted forcing,
and has not been composed with the oscillatory/correction lanes into a complete
Navier--Stokes candidate. Agent-4 #966 assesses scoped leading-only divergence;
it does not assess momentum residual or canonical whole-domain quadrature.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
import re
from typing import Any, Mapping

SCHEMA_NAME = "kokuno-agent5-current-cartesian-leading-ingest-v1"
TASK_ID = "KOKUNO-A5-CURRENT-CARTESIAN-LEADING-INGEST-092"
_HEX40 = re.compile(r"^[0-9a-f]{40}$")

PARENT_A5 = {
    "pr": 963,
    "head": "0a17e4b8d2ca90997bc9397d548cef94d792ccbe",
    "branch": "codex/kokuno-a5-current-pa16-profile-application-ingest-091",
    "source_path": "src/openai_ns_reconstruction/kokuno_a5_current_pa16_profile_application_ingest_contract.py",
    "observed_ci": {
        "dedicated": {"run_id": 35566807300, "status": "queued", "conclusion": None},
        "repository_tests": {"run_id": 35566807506, "status": "queued", "conclusion": None},
    },
}

AGENT1_CARTESIAN_LEADING = {
    "pr": 965,
    "head": "7dc09c59584a6b1bfaabcbb8df00283d9f9a40c1",
    "role": "current-lineage candidate-side Cartesian leading velocity through current X_h",
    "source_path": "src/openai_ns_reconstruction/kokuno_pa16_current_cartesian_leading_velocity.py",
    "source_blob": "ca8b80b0451be1a8f31deaf620f542d8c2e92c0c",
    "test_path": "tests/test_kokuno_pa16_current_cartesian_leading_velocity.py",
    "test_blob": "c8c4322789df1b24b1499ff574571243ffc073e4",
    "workflow_path": ".github/workflows/kokuno-agent1-current-cartesian-leading.yml",
    "workflow_blob": "794d85ec77d1a2cfae121bf4d11cbe89017df44a",
    "parent_joined_profile_pr": 959,
    "parent_joined_profile_head": "61b8730fba5623c5acfbd3ec42d54f7cc4f45748",
    "public_api": ["velocity", "save_configuration", "load_configuration"],
    "vectorized_cartesian_velocity": True,
    "candidate_save_reload": True,
    "reuses_governed_q_inverse": True,
    "velocity_materialized_through_xh": True,
    "velocity_beyond_xh_materialized": False,
    "global_compact_supported_leading_velocity": False,
    "matched_pressure": False,
    "restricted_forcing": False,
    "complete_ns_residual": False,
    "observed_ci": {
        "dedicated": {"run_id": 35568968563, "status": "queued", "conclusion": None},
        "repository_tests": {"run_id": 35568968517, "status": "queued", "conclusion": None},
    },
}

AGENT4_CARTESIAN_DIVERGENCE_AUDIT = {
    "pr": 966,
    "head": "50485166808ac9d942ab017708c23438dfdb01a6",
    "role": "implementation-distinct public-velocity-only divergence audit of A1 #965 leading field",
    "source_path": "src/openai_ns_reconstruction/kokuno_a4_current_cartesian_leading_divergence_independent_audit.py",
    "source_blob": "c0773b7c30a99216147dcc605a9d61b88ebbccb5",
    "test_path": "tests/test_constrained_kokuno_a4_current_cartesian_leading_divergence_independent_audit.py",
    "test_blob": "2780fdb995a344052c8801ee5e101940ed245fa0",
    "workflow_path": ".github/workflows/kokuno-agent4-current-cartesian-leading-divergence-audit.yml",
    "workflow_blob": "a68106945d939a49a4edf2c7ceb9d7ba61aed765",
    "audited_agent1_pr": 965,
    "audited_agent1_head": "7dc09c59584a6b1bfaabcbb8df00283d9f9a40c1",
    "public_velocity_only_scientific_path": True,
    "candidate_save_reload_required": True,
    "derivative_operator": "centered Cartesian FD2",
    "protocol": {
        "seed": 9173661,
        "times": [0.31, 0.47, 0.63, 0.71],
        "source_x_interval": [0.02, 1.20],
        "source_eta_interval": [-0.65, 0.65],
        "points_per_time": 8,
        "integration_probe_count": 32,
        "axis_probe_count": 5,
        "spatial_step_ladder": [0.02, 0.01, 0.005],
        "finest_sampled_max_gate": 1.0e-5,
        "finest_weighted_rms_gate": 1.0e-5,
        "finest_axis_max_gate": 1.0e-5,
        "nontrivial_speed_rms_floor": 1.0e-10,
        "medium_to_fine_degradation_factor": 1.25,
        "stability_floor": 2.0e-8,
    },
    "audit_registered": True,
    "audit_admitted": False,
    "scoped_leading_only_divergence": True,
    "canonical_whole_domain_divergence_l2": False,
    "momentum_residual_assessed": False,
    "complete_ns_residual_audit": False,
    "same_protocol_comparable_to_st006": False,
    "observed_ci": {
        "dedicated": {"run_id": 35570193367, "status": "queued", "conclusion": None},
        "repository_tests": {"run_id": 35570193372, "status": "queued", "conclusion": None},
    },
}

AGENT2_SIBLING = {
    "pr": 960,
    "head": "6d2fb1f701a34f783dca15a267ae2ce0734ba741",
    "role": "axis-safe batch oscillatory velocity differentials",
    "source_path": "src/openai_ns_reconstruction/kokuno_oscillatory_batch_differentials.py",
    "source_blob": "12df3bf6c949baeebf609b366f2973e7f515a157",
    "coupled_into_current_cartesian_leading": False,
    "complete_ns_residual": False,
    "observed_ci": {
        "dedicated": {"run_id": 35565320969, "status": "queued", "conclusion": None},
        "repository_tests": {"run_id": 35565320955, "status": "queued", "conclusion": None},
    },
}

AGENT3_SIBLING = {
    "pr": 961,
    "head": "670f7b71d72735d441dcf48f2391c271ee5c163c",
    "role": "current PA.16 source-profile repair application",
    "source_path": "src/openai_ns_reconstruction/kokuno_current_pa16_profile_application.py",
    "source_blob": "cb3dcd25e0c55a7935daee1b77ffb1741320f83b",
    "source_profile_application_materialized": True,
    "real_ns_correction_velocity_materialized": False,
    "discrepancy_from_complete_ns_defect": False,
    "authorized_for_gain_gated_ns_stage": False,
    "observed_ci": {
        "dedicated": {"run_id": 35565700824, "status": "queued", "conclusion": None},
        "repository_tests": {"run_id": 35565700813, "status": "queued", "conclusion": None},
    },
}

FROZEN_SCIENCE = {
    "viscosity": 0.01,
    "evaluation_box": [[-2.0, 2.0], [-2.0, 2.0], [-2.0, 2.0]],
    "time_interval": [0.25, 0.75],
    "forcing_family": "restricted two-parameter curl forcing",
    "residual_defined_free_forcing_forbidden": True,
}

FINAL_GATE = {
    "normalized_momentum_sampled_max": 1.0e-3,
    "normalized_momentum_volume_l2": 1.0e-3,
    "divergence_sampled_max": 1.0e-5,
    "divergence_volume_l2": 1.0e-5,
    "canonical_volume_quadrature_ladder": [24, 48, 96],
}

ST006_BASELINE = {
    "same_protocol_full_pde_baseline": True,
    "momentum_sampled_max": 0.1082289305112118,
    "momentum_volume_l2": 0.10758432876230622,
}

READINESS = {
    "leading_ready": False,
    "oscillatory_ready": True,
    "correction_ready": False,
    "velocity_export_ready": False,
    "pde_validated": False,
}

TRUTH_BOUNDARY = {
    "current_joined_source_profile_materialized": True,
    "current_pa16_source_profile_application_materialized": True,
    "current_cartesian_leading_velocity_materialized_through_xh": True,
    "current_cartesian_leading_save_load_materialized": True,
    "current_cartesian_leading_a4_divergence_audit_registered": True,
    "current_cartesian_leading_a4_divergence_audit_admitted": False,
    "current_cartesian_leading_divergence_scoped_assessed": False,
    "global_cartesian_leading_velocity_materialized": False,
    "velocity_beyond_xh_materialized": False,
    "leading_plus_oscillatory_cartesian_velocity_materialized": False,
    "complete_velocity_candidate_materialized": False,
    "matched_cartesian_pressure_materialized": False,
    "cartesian_matched_pressure_gradient_materialized": False,
    "restricted_forcing_materialized": False,
    "discrepancy_from_complete_ns_defect": False,
    "authorized_for_gain_gated_ns_stage": False,
    "real_agent3_ns_correction_velocity_materialized": False,
    "real_candidate_finite_correction_cycle_run": False,
    "complete_candidate_api_ready": False,
    "canonical_whole_domain_divergence_l2_assessed": False,
    "same_protocol_full_ns_residual_available": False,
    "same_protocol_st006_comparison_available_now": False,
    "current_candidate_eligible_for_full_ns_validation": False,
    "upstream_ci_admitted_as_pass": False,
    "scientific_admission": False,
    "paper_exact": False,
}

PIPELINE_POSITION = {
    "stage": "current Cartesian leading-velocity registration through X_h",
    "input": "A5 #963 joined-profile/application seam plus A1 #965 Cartesian leading assembly",
    "new_output": "checksum-bound registration of callable/saveable Cartesian leading velocity through X_h and A4 #966 scoped divergence audit protocol",
    "independent_audit_scope": "A4 #966 save/reloads A1 #965 and differentiates only public velocity with centered Cartesian FD2 on a frozen held-out source-cell protocol",
    "not_output": "global velocity beyond X_h, leading+oscillatory composition, matched pressure/forcing, momentum residual, correction velocity, complete candidate artifact, canonical whole-domain admission, or PDE validation",
    "next_shortest_blocker": "materialize the global/outer leading completion beyond X_h or the exact composition seam needed to form one complete velocity candidate, then add matched pressure/restricted forcing before complete residual work",
}


def _canonical_json(payload: Mapping[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _sha256(payload: Mapping[str, Any]) -> str:
    return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()


def _require_hex40(value: Any, label: str) -> None:
    if not isinstance(value, str) or _HEX40.fullmatch(value) is None:
        raise ValueError(f"{label} must be a lowercase 40-hex SHA")


def build_contract(exact_head: str) -> dict[str, Any]:
    _require_hex40(exact_head, "exact_head")
    payload: dict[str, Any] = {
        "schema": SCHEMA_NAME,
        "task_id": TASK_ID,
        "agent5_exact_head": exact_head,
        "parent_a5": copy.deepcopy(PARENT_A5),
        "agent1_cartesian_leading": copy.deepcopy(AGENT1_CARTESIAN_LEADING),
        "agent4_cartesian_divergence_audit": copy.deepcopy(AGENT4_CARTESIAN_DIVERGENCE_AUDIT),
        "agent2_sibling": copy.deepcopy(AGENT2_SIBLING),
        "agent3_sibling": copy.deepcopy(AGENT3_SIBLING),
        "frozen_science": copy.deepcopy(FROZEN_SCIENCE),
        "final_gate": copy.deepcopy(FINAL_GATE),
        "st006_baseline": copy.deepcopy(ST006_BASELINE),
        "readiness": copy.deepcopy(READINESS),
        "truth_boundary": copy.deepcopy(TRUTH_BOUNDARY),
        "pipeline_position": copy.deepcopy(PIPELINE_POSITION),
        "registration_only": True,
    }
    payload["contract_sha256"] = _sha256(payload)
    return payload


def validate_contract(payload: Mapping[str, Any]) -> list[str]:
    errors: list[str] = []
    observed = copy.deepcopy(dict(payload))
    supplied_sha = observed.pop("contract_sha256", None)
    if supplied_sha != _sha256(observed):
        errors.append("contract_sha256_mismatch")

    expected = {
        "schema": SCHEMA_NAME,
        "task_id": TASK_ID,
        "parent_a5": PARENT_A5,
        "agent1_cartesian_leading": AGENT1_CARTESIAN_LEADING,
        "agent4_cartesian_divergence_audit": AGENT4_CARTESIAN_DIVERGENCE_AUDIT,
        "agent2_sibling": AGENT2_SIBLING,
        "agent3_sibling": AGENT3_SIBLING,
        "frozen_science": FROZEN_SCIENCE,
        "final_gate": FINAL_GATE,
        "st006_baseline": ST006_BASELINE,
        "readiness": READINESS,
        "truth_boundary": TRUTH_BOUNDARY,
        "pipeline_position": PIPELINE_POSITION,
        "registration_only": True,
    }
    for key, value in expected.items():
        if observed.get(key) != value:
            errors.append(f"{key}_drift")

    try:
        _require_hex40(observed.get("agent5_exact_head"), "agent5_exact_head")
    except ValueError:
        errors.append("agent5_exact_head_invalid")

    a1 = observed.get("agent1_cartesian_leading", {})
    a4 = observed.get("agent4_cartesian_divergence_audit", {})
    a2 = observed.get("agent2_sibling", {})
    a3 = observed.get("agent3_sibling", {})
    truth = observed.get("truth_boundary", {})
    readiness = observed.get("readiness", {})
    science = observed.get("frozen_science", {})
    gate = observed.get("final_gate", {})

    if a4.get("audited_agent1_head") != a1.get("head") or a4.get("audited_agent1_pr") != a1.get("pr"):
        errors.append("a1_a4_lineage_mismatch")
    if a1.get("velocity_materialized_through_xh") is not True or truth.get("current_cartesian_leading_velocity_materialized_through_xh") is not True:
        errors.append("cartesian_leading_materialization_lost")
    if a1.get("velocity_beyond_xh_materialized") is not False or truth.get("velocity_beyond_xh_materialized") is not False:
        errors.append("global_leading_laundering")
    if a1.get("global_compact_supported_leading_velocity") is not False or truth.get("global_cartesian_leading_velocity_materialized") is not False:
        errors.append("global_leading_laundering")
    if a4.get("momentum_residual_assessed") is not False or a4.get("complete_ns_residual_audit") is not False:
        errors.append("a4_scope_laundering")
    if a4.get("canonical_whole_domain_divergence_l2") is not False or truth.get("canonical_whole_domain_divergence_l2_assessed") is not False:
        errors.append("canonical_norm_laundering")
    if truth.get("current_cartesian_leading_a4_divergence_audit_admitted") is not False or truth.get("current_cartesian_leading_divergence_scoped_assessed") is not False:
        errors.append("queued_evidence_laundering")
    if a2.get("coupled_into_current_cartesian_leading") is not False:
        errors.append("oscillatory_composition_laundering")
    if truth.get("leading_plus_oscillatory_cartesian_velocity_materialized") is not False:
        errors.append("oscillatory_composition_laundering")
    if a3.get("real_ns_correction_velocity_materialized") is not False or truth.get("real_agent3_ns_correction_velocity_materialized") is not False:
        errors.append("correction_velocity_laundering")
    if a3.get("discrepancy_from_complete_ns_defect") is not False or truth.get("discrepancy_from_complete_ns_defect") is not False:
        errors.append("complete_defect_laundering")
    if a3.get("authorized_for_gain_gated_ns_stage") is not False or truth.get("authorized_for_gain_gated_ns_stage") is not False:
        errors.append("gain_stage_laundering")
    if truth.get("matched_cartesian_pressure_materialized") is not False or truth.get("restricted_forcing_materialized") is not False:
        errors.append("pressure_forcing_laundering")
    if truth.get("same_protocol_full_ns_residual_available") is not False or truth.get("same_protocol_st006_comparison_available_now") is not False:
        errors.append("residual_comparison_laundering")
    if science.get("residual_defined_free_forcing_forbidden") is not True:
        errors.append("free_forcing_firewall_weakened")
    if gate.get("normalized_momentum_sampled_max") != 1.0e-3 or gate.get("normalized_momentum_volume_l2") != 1.0e-3:
        errors.append("momentum_gate_drift")
    if gate.get("divergence_sampled_max") != 1.0e-5 or gate.get("divergence_volume_l2") != 1.0e-5:
        errors.append("divergence_gate_drift")
    if gate.get("canonical_volume_quadrature_ladder") != [24, 48, 96]:
        errors.append("canonical_quadrature_drift")
    if any(readiness.get(k) is not False for k in ("leading_ready", "correction_ready", "velocity_export_ready", "pde_validated")):
        errors.append("readiness_premature_promotion")
    if readiness.get("oscillatory_ready") is not True:
        errors.append("oscillatory_readiness_drift")
    if truth.get("scientific_admission") is not False or truth.get("upstream_ci_admitted_as_pass") is not False:
        errors.append("scientific_admission_laundering")
    if truth.get("paper_exact") is not False:
        errors.append("paper_exact_laundering")
    return sorted(set(errors))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exact-head", required=True)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    payload = build_contract(args.exact_head)
    errors = validate_contract(payload)
    payload["validation_errors"] = errors
    rendered = json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if args.out is not None:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
