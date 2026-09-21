"""Agent-5 integration contract for the current Cartesian leading field through X_R.

This module is integration/provenance glue only.  It registers Agent-1 PR #980,
which extends the current candidate-side Cartesian leading velocity from the
previous X_h boundary through the public exterior-preservation interval ending
at X_R.  It does not reimplement Agent-1 mathematics and it does not promote
that scoped extension into a global compactly-supported velocity field.

Current registered seam:

    A5 #978 current correction-side/radial-stress registration
      -> A1 #980 current Cartesian leading velocity through X_R

Agent-2 #981 (identity-preserving save/load for the older partial
leading+oscillatory composite through X_h) and Agent-3 #982 (current partial
radial force from d_z sigma_1) are recorded as siblings only.  Neither is
silently retargeted to the new X_R leading lineage here.  No Agent-4 audit of
#980 exists at this freeze, so independent exterior validation stays false.

The fixed full-PDE gates remain unchanged.  In particular, no scoped profile,
divergence, stress, or delivery receipt is equivalent to held-out normalized
Navier--Stokes validation.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import re
from pathlib import Path
from typing import Any, Mapping

SCHEMA_NAME = "kokuno-agent5-current-cartesian-exterior-xr-ingest-v1"
TASK_ID = "KOKUNO-A5-CURRENT-CARTESIAN-EXTERIOR-XR-INGEST-095"
_HEX40 = re.compile(r"^[0-9a-f]{40}$")

PARENT_A5 = {
    "pr": 978,
    "head": "cb67fed42c459d9d753ad17cbc339744e51a9b84",
    "branch": "codex/kokuno-a5-current-partial-radial-stress-ingest-094",
    "role": "current partial nonlinear radial-stress registration",
}

AGENT1_EXTERIOR_TO_XR = {
    "pr": 980,
    "head": "d3c971f2c62e272333e124e532212d23cca4908d",
    "branch": "agent-kokuno-1/current-cartesian-exterior-to-xr-096",
    "base_pr": 965,
    "base_head": "7dc09c59584a6b1bfaabcbb8df00283d9f9a40c1",
    "source_path": "src/openai_ns_reconstruction/kokuno_pa16_current_cartesian_exterior_to_xr.py",
    "source_blob": "2ba28636a56b252fa485719e7b3e8da754d7e588",
    "test_path": "tests/test_kokuno_pa16_current_cartesian_exterior_to_xr.py",
    "test_blob": "ecbab9822c54f440f5555f25cd5e7e1648d84bd2",
    "workflow_path": ".github/workflows/kokuno-agent1-current-cartesian-exterior-to-xr.yml",
    "workflow_blob": "5c2d68211e8fea526e06087fc4475f7bab9eb563",
    "public_class": "KokunoPA16CurrentCartesianExteriorToXR",
    "velocity_api": "velocity(x,y,z,t)->[...,3]",
    "profile_api": "similarity_profile_values(X,eta)",
    "radial_derivative_api": "similarity_radial_derivatives(X,eta)",
    "save_load": True,
    "semantic_identity": True,
    "materialized_domain": "0<=X<=X_R",
    "fails_closed_after_X_R": True,
    "post_XR_RF40_materialized": False,
    "global_compact_support_completed": False,
    "matched_pressure_materialized": False,
    "restricted_forcing_materialized": False,
    "heldout_complete_ns_residual_assessed": False,
}

AGENT2_IDENTITY_SAVE_LOAD_SIBLING = {
    "pr": 981,
    "head": "9997fc55455d126f935643da36bf17eaa0491aa4",
    "branch": "codex/kokuno-a2-partial-composite-identity-save-load-082",
    "role": "identity-preserving save/load for the existing partial u_lead(#965 through X_h)+u_osc composite",
    "full_concrete_oscillatory_runtime_digest_bound": True,
    "identity_preserving_partial_composite_save_load_available": True,
    "composite_extended_to_agent1_980_X_R": False,
    "consumed_by_this_increment": False,
}

AGENT3_RADIAL_FORCE_SIBLING = {
    "pr": 982,
    "head": "dd11348e8eecf298656d7da58aee4526c6102306",
    "branch": "codex/kokuno-a3-current-partial-nonlinear-radial-force-104",
    "role": "current partial radial force from physical-z derivative of axial/e=1 compact stress",
    "source_relation": "(div T)_r = partial_z sigma_1",
    "complete_ns_defect": False,
    "authorized_as_correction_target": False,
    "consumes_agent1_980_exterior": False,
    "consumed_by_this_increment": False,
}

LATEST_MATCHING_A4 = {
    "pr": 977,
    "head": "2b015feccb5d55735ced8763fcb05a563a4cf12e",
    "role": "implementation-distinct audit of A3 #976 theta/axial scoped radial stresses",
    "audits_agent1_980_exterior": False,
    "agent1_980_independent_a4_audit_available": False,
}

OBSERVED_CI_AT_FREEZE = {
    "parent_a5_978": {
        "repository_tests_run": 35580677745,
        "dedicated_run": 35580677804,
        "status": "queued",
        "conclusion": None,
    },
    "agent1_980": {
        "repository_tests_run": 35583982998,
        "dedicated_run": 35583983106,
        "status": "queued",
        "conclusion": None,
    },
    "agent2_981": {
        "repository_tests_run": 35584060033,
        "dedicated_run": 35584060095,
        "status": "queued",
        "conclusion": None,
    },
    "agent3_982": {
        "repository_tests_run": 35584768010,
        "dedicated_run": 35584768136,
        "status": "queued",
        "conclusion": None,
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
    "current_cartesian_leading_velocity_materialized_through_xh": True,
    "current_cartesian_leading_velocity_materialized_through_xr": True,
    "velocity_beyond_xh_materialized": True,
    "velocity_beyond_xr_materialized": False,
    "agent1_980_save_load_materialized": True,
    "agent1_980_semantic_identity_materialized": True,
    "agent1_980_independent_a4_audit_registered": False,
    "agent1_980_independent_a4_audit_admitted": False,
    "post_xr_rf40_current_lineage_materialized": False,
    "cone_i1_i2_i3_i4_outer_overlays_completed": False,
    "outer_global_leading_velocity_materialized": False,
    "global_compact_support_completed": False,
    "current_partial_composite_save_load_materialized_through_xh": True,
    "current_partial_composite_extended_through_xr": False,
    "matched_cartesian_pressure_materialized": False,
    "cartesian_matched_pressure_gradient_materialized": False,
    "restricted_forcing_materialized": False,
    "complete_ns_defect_materialized": False,
    "real_agent3_ns_correction_velocity_materialized": False,
    "real_candidate_finite_correction_cycle_run": False,
    "complete_candidate_api_ready": False,
    "canonical_whole_domain_divergence_l2_assessed": False,
    "heldout_normalized_ns_residual_assessed": False,
    "same_protocol_st006_comparison_available_now": False,
    "residual_reduction_claimed": False,
    "scientific_admission": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}

PIPELINE_POSITION = {
    "stage": "current Cartesian leading extension from X_h through X_R",
    "input": "A5 #978 plus exact A1 #980; A2 #981 and A3 #982 are non-consumed siblings",
    "new_output": "checksum-bound registration of a callable/save-load Cartesian leading velocity through X_R",
    "not_output": "post-X_R/global velocity, X_R-matched leading+oscillatory composite, pressure, forcing, correction velocity, full residual, or PDE validation",
    "next_shortest_blocker": "post-X_R RF40 plus cone/I1-I4 outer completion, then recompose the exact oscillatory runtime on that extended/global leading identity before matched pressure/restricted forcing and complete-defect correction",
}


def _canonical_json(payload: Mapping[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _sha256(payload: Mapping[str, Any]) -> str:
    return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()


def _require_hex40(value: Any, label: str) -> None:
    if not isinstance(value, str) or _HEX40.fullmatch(value) is None:
        raise ValueError(f"{label} must be a lowercase 40-hex SHA")


def build_contract(exact_head: str) -> dict[str, Any]:
    """Build the immutable Agent-5 registration receipt for this increment."""
    _require_hex40(exact_head, "exact_head")
    payload: dict[str, Any] = {
        "schema": SCHEMA_NAME,
        "task_id": TASK_ID,
        "agent5_exact_head": exact_head,
        "parent_a5": copy.deepcopy(PARENT_A5),
        "agent1_exterior_to_xr": copy.deepcopy(AGENT1_EXTERIOR_TO_XR),
        "agent2_identity_save_load_sibling": copy.deepcopy(AGENT2_IDENTITY_SAVE_LOAD_SIBLING),
        "agent3_radial_force_sibling": copy.deepcopy(AGENT3_RADIAL_FORCE_SIBLING),
        "latest_matching_a4": copy.deepcopy(LATEST_MATCHING_A4),
        "observed_ci_at_freeze": copy.deepcopy(OBSERVED_CI_AT_FREEZE),
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
    """Return fail-closed validation errors for a serialized registration receipt."""
    errors: list[str] = []
    observed = copy.deepcopy(dict(payload))
    supplied_sha = observed.pop("contract_sha256", None)
    if supplied_sha != _sha256(observed):
        errors.append("contract_sha256_mismatch")

    expected = {
        "schema": SCHEMA_NAME,
        "task_id": TASK_ID,
        "parent_a5": PARENT_A5,
        "agent1_exterior_to_xr": AGENT1_EXTERIOR_TO_XR,
        "agent2_identity_save_load_sibling": AGENT2_IDENTITY_SAVE_LOAD_SIBLING,
        "agent3_radial_force_sibling": AGENT3_RADIAL_FORCE_SIBLING,
        "latest_matching_a4": LATEST_MATCHING_A4,
        "observed_ci_at_freeze": OBSERVED_CI_AT_FREEZE,
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

    truth = observed.get("truth_boundary", {})
    readiness = observed.get("readiness", {})
    a1 = observed.get("agent1_exterior_to_xr", {})
    a2 = observed.get("agent2_identity_save_load_sibling", {})
    a3 = observed.get("agent3_radial_force_sibling", {})
    a4 = observed.get("latest_matching_a4", {})

    if truth.get("current_cartesian_leading_velocity_materialized_through_xr") is not True:
        errors.append("xr_materialization_lost")
    if truth.get("velocity_beyond_xr_materialized") is not False:
        errors.append("post_xr_scope_promoted")
    if truth.get("outer_global_leading_velocity_materialized") is not False:
        errors.append("global_leading_promoted")
    if truth.get("current_partial_composite_extended_through_xr") is not False:
        errors.append("composite_lineage_laundered_to_xr")
    if a1.get("fails_closed_after_X_R") is not True:
        errors.append("agent1_xr_fail_closed_lost")
    if a2.get("consumed_by_this_increment") is not False:
        errors.append("agent2_sibling_silently_consumed")
    if a2.get("composite_extended_to_agent1_980_X_R") is not False:
        errors.append("agent2_partial_composite_silently_extended")
    if a3.get("consumed_by_this_increment") is not False:
        errors.append("agent3_sibling_silently_consumed")
    if a3.get("consumes_agent1_980_exterior") is not False:
        errors.append("agent3_lineage_silently_retargeted")
    if a4.get("agent1_980_independent_a4_audit_available") is not False:
        errors.append("missing_a4_audit_promoted")

    for key in ("leading_ready", "correction_ready", "velocity_export_ready", "pde_validated"):
        if readiness.get(key) is not False:
            errors.append(f"readiness_{key}_promoted")
    if readiness.get("oscillatory_ready") is not True:
        errors.append("oscillatory_ready_regressed")

    forbidden_true = (
        "global_compact_support_completed",
        "matched_cartesian_pressure_materialized",
        "restricted_forcing_materialized",
        "complete_ns_defect_materialized",
        "real_agent3_ns_correction_velocity_materialized",
        "heldout_normalized_ns_residual_assessed",
        "same_protocol_st006_comparison_available_now",
        "residual_reduction_claimed",
        "scientific_admission",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    )
    for key in forbidden_true:
        if truth.get(key) is not False:
            errors.append(f"truth_{key}_promoted")

    return errors


def write_contract(path: str | Path, exact_head: str) -> dict[str, Any]:
    payload = build_contract(exact_head)
    errors = validate_contract(payload)
    if errors:
        raise ValueError("invalid registration contract: " + ", ".join(errors))
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return payload


def _main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--head", required=True, help="exact Agent-5 commit SHA")
    parser.add_argument("--out", required=True, help="JSON receipt output path")
    args = parser.parse_args()
    payload = write_contract(args.out, args.head)
    print(json.dumps({
        "schema": payload["schema"],
        "task_id": payload["task_id"],
        "contract_sha256": payload["contract_sha256"],
        "readiness": payload["readiness"],
        "truth_boundary": payload["truth_boundary"],
    }, sort_keys=True))


if __name__ == "__main__":
    _main()
