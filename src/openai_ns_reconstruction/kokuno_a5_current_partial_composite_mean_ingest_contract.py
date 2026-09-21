"""Agent-5 integration contract for the current partial composite/mean seam.

This module is integration/provenance glue only.  It registers two fresh sibling
artifacts without reimplementing either mathematical lane:

* Agent 2 PR #970: current Cartesian leading velocity through ``X_h`` plus the
  frozen complete-curl oscillatory field;
* Agent 3 PR #971: mixed, quadratic and aggregate cylindrical ``m=0`` nonlinear
  mean attribution for that exact partial candidate.

Agent 4 PR #966 remains a *leading-only* scoped divergence audit.  No independent
A4 audit currently covers the A2 composite or A3 nonlinear mean, and this module
must not launder that absence into scientific admission.

The registered candidate is still partial-domain through the current ``X_h``.
There is no outer/global completion, matched pressure, preregistered restricted
forcing, complete NS defect, correction velocity, same-protocol full residual or
PDE validation in this increment.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
import re
from typing import Any, Mapping

SCHEMA_NAME = "kokuno-agent5-current-partial-composite-mean-ingest-v1"
TASK_ID = "KOKUNO-A5-CURRENT-PARTIAL-COMPOSITE-MEAN-INGEST-093"
_HEX40 = re.compile(r"^[0-9a-f]{40}$")

PARENT_A5 = {
    "pr": 968,
    "head": "5558d0859d142914ae8921a4286f56077fe2232c",
    "branch": "codex/kokuno-a5-current-cartesian-leading-ingest-092",
    "role": "current Cartesian leading registration through current X_h",
    "observed_ci": {
        "dedicated": {"run_id": 35570786470, "status": "queued", "conclusion": None},
        "repository_tests": {"run_id": 35570786541, "status": "queued", "conclusion": None},
    },
}

AGENT1_LEADING = {
    "pr": 965,
    "head": "7dc09c59584a6b1bfaabcbb8df00283d9f9a40c1",
    "source_path": "src/openai_ns_reconstruction/kokuno_pa16_current_cartesian_leading_velocity.py",
    "source_blob": "ca8b80b0451be1a8f31deaf620f542d8c2e92c0c",
    "role": "current-lineage Cartesian leading velocity through current X_h",
    "velocity_beyond_xh_materialized": False,
}

AGENT2_PARTIAL_COMPOSITE = {
    "pr": 970,
    "head": "3a6405bbd3d10b8c3b38078f0c989c45e5d407b4",
    "branch": "codex/kokuno-a2-current-partial-leading-oscillatory-080",
    "source_path": "src/openai_ns_reconstruction/kokuno_current_partial_leading_oscillatory_velocity.py",
    "source_blob": "532435705bd5db341e82371180da8c04c3dd6c14",
    "parent_oscillatory_differentials_pr": 960,
    "parent_oscillatory_differentials_head": "6d2fb1f701a34f783dca15a267ae2ce0734ba741",
    "parent_oscillatory_differentials_source_blob": "12df3bf6c949baeebf609b366f2973e7f515a157",
    "consumed_agent1_pr": 965,
    "consumed_agent1_head": "7dc09c59584a6b1bfaabcbb8df00283d9f9a40c1",
    "public_api": "velocity(x,y,z,t)",
    "composition": "u_partial = u_lead_current_through_Xh + u_osc_frozen_complete_curl",
    "partial_domain_through_xh": True,
    "velocity_beyond_xh_materialized": False,
    "outer_global_leading_velocity_materialized": False,
    "correction_velocity_composed": False,
    "matched_pressure": False,
    "restricted_forcing": False,
    "complete_ns_residual": False,
    "observed_ci": {
        "dedicated": {"run_id": 35573770388, "status": "queued", "conclusion": None},
        "repository_tests": {"run_id": 35573770361, "status": "queued", "conclusion": None},
    },
}

AGENT3_PARTIAL_NONLINEAR_MEAN = {
    "pr": 971,
    "head": "bdb02ad487a175ab62748687c3172127d09aeb9c",
    "branch": "codex/kokuno-a3-current-partial-nonlinear-mean-102",
    "source_path": "src/openai_ns_reconstruction/kokuno_current_partial_nonlinear_mean_attribution.py",
    "source_blob": "9c62ae3c3b4b3ca2d9082f822d57feef700ea6cf",
    "test_path": "tests/test_constrained_kokuno_current_partial_nonlinear_mean_attribution.py",
    "test_blob": "4e70f5b7b51105e00ff44b32a2fe16473854c810",
    "workflow_path": ".github/workflows/kokuno-agent3-current-partial-nonlinear-mean.yml",
    "workflow_blob": "2b8514c25dc74e7b2aa4df43cf3cf23d00993572",
    "consumed_agent2_composite_pr": 970,
    "consumed_agent2_composite_head": "3a6405bbd3d10b8c3b38078f0c989c45e5d407b4",
    "consumed_agent2_differentials_head": "6d2fb1f701a34f783dca15a267ae2ce0734ba741",
    "consumed_agent1_leading_head": "7dc09c59584a6b1bfaabcbb8df00283d9f9a40c1",
    "fixed_leading_cartesian_fd_step": 1.0e-3,
    "nonlinear_pieces": [
        "A=(u_lead.grad)u_osc",
        "B=(u_osc.grad)u_lead",
        "Q=(u_osc.grad)u_osc",
        "N=A+B+Q",
    ],
    "outputs": ["mixed_m0_mean", "quadratic_m0_mean", "aggregate_m0_mean"],
    "frozen_scoped_receipt_gates": {
        "pointwise_three_piece_closure_absolute_max": 5.0e-13,
        "projected_three_piece_closure_absolute_max": 1.0e-12,
        "quadratic_mean_rms_min": 1.0e-12,
    },
    "radial_inverse_performed": False,
    "pressure_gradient_included": False,
    "restricted_forcing_included": False,
    "complete_ns_defect": False,
    "authorized_as_correction_target": False,
    "correction_velocity_materialized": False,
    "observed_ci": {
        "dedicated": {"run_id": 35574362920, "status": "queued", "conclusion": None},
        "repository_tests": {"run_id": 35574362841, "status": "queued", "conclusion": None},
    },
}

AGENT4_LEADING_ONLY_AUDIT = {
    "pr": 966,
    "head": "50485166808ac9d942ab017708c23438dfdb01a6",
    "source_path": "src/openai_ns_reconstruction/kokuno_a4_current_cartesian_leading_divergence_independent_audit.py",
    "source_blob": "c0773b7c30a99216147dcc605a9d61b88ebbccb5",
    "role": "independent public-velocity-only scoped divergence audit of A1 #965 leading field",
    "audits_agent1_leading": True,
    "audits_agent2_partial_composite": False,
    "audits_agent3_partial_nonlinear_mean": False,
    "independent_composite_audit_available": False,
    "independent_nonlinear_mean_audit_available": False,
    "canonical_whole_domain_admission": False,
    "momentum_residual_assessed": False,
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
    "current_partial_leading_plus_oscillatory_velocity_materialized_through_xh": True,
    "current_partial_mixed_nonlinear_mean_materialized": True,
    "current_partial_quadratic_nonlinear_mean_materialized": True,
    "current_partial_aggregate_nonlinear_mean_materialized": True,
    "velocity_beyond_xh_materialized": False,
    "outer_global_leading_velocity_materialized": False,
    "global_compact_support_completed": False,
    "independent_a4_composite_audit_available": False,
    "independent_a4_nonlinear_mean_audit_available": False,
    "radial_inverse_performed_on_current_partial_mean": False,
    "scoped_nonlinear_mean_authorized_as_correction_target": False,
    "real_agent3_ns_correction_velocity_materialized": False,
    "matched_cartesian_pressure_materialized": False,
    "cartesian_matched_pressure_gradient_materialized": False,
    "restricted_forcing_materialized": False,
    "complete_ns_defect_materialized": False,
    "complete_candidate_api_ready": False,
    "canonical_whole_domain_divergence_l2_assessed": False,
    "same_protocol_full_ns_residual_available": False,
    "same_protocol_st006_comparison_available_now": False,
    "current_candidate_eligible_for_full_ns_validation": False,
    "residual_reduction_claimed": False,
    "scientific_admission": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}

PIPELINE_POSITION = {
    "stage": "current partial Cartesian leading+oscillatory composition and nonlinear-mean registration through X_h",
    "input": "A5 #968 leading seam + A2 #970 partial composite + A3 #971 nonlinear mean attribution",
    "new_output": "checksum-bound registration of the first current-lineage leading+oscillatory callable and its mixed/quadratic/aggregate m=0 nonlinear means",
    "a4_scope": "A4 #966 remains leading-only; no A4 audit currently covers the composite or nonlinear means",
    "not_output": "velocity beyond X_h, global compact support, radial inverse, complete NS defect, correction velocity, matched pressure/forcing, same-protocol residual, canonical whole-domain admission, or PDE validation",
    "next_shortest_blocker": "harvest outer/global completion beyond X_h and/or a scoped radial-inverse adapter for the current mean; do not authorize a correction stage until matched pressure plus preregistered restricted forcing define a complete NS defect",
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
        "agent1_leading": copy.deepcopy(AGENT1_LEADING),
        "agent2_partial_composite": copy.deepcopy(AGENT2_PARTIAL_COMPOSITE),
        "agent3_partial_nonlinear_mean": copy.deepcopy(AGENT3_PARTIAL_NONLINEAR_MEAN),
        "agent4_leading_only_audit": copy.deepcopy(AGENT4_LEADING_ONLY_AUDIT),
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
        "agent1_leading": AGENT1_LEADING,
        "agent2_partial_composite": AGENT2_PARTIAL_COMPOSITE,
        "agent3_partial_nonlinear_mean": AGENT3_PARTIAL_NONLINEAR_MEAN,
        "agent4_leading_only_audit": AGENT4_LEADING_ONLY_AUDIT,
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

    a1 = observed.get("agent1_leading", {})
    a2 = observed.get("agent2_partial_composite", {})
    a3 = observed.get("agent3_partial_nonlinear_mean", {})
    a4 = observed.get("agent4_leading_only_audit", {})
    truth = observed.get("truth_boundary", {})
    readiness = observed.get("readiness", {})
    gate = observed.get("final_gate", {})
    science = observed.get("frozen_science", {})

    if a2.get("consumed_agent1_head") != a1.get("head"):
        errors.append("a1_a2_lineage_mismatch")
    if a3.get("consumed_agent2_composite_head") != a2.get("head"):
        errors.append("a2_a3_lineage_mismatch")
    if a3.get("consumed_agent1_leading_head") != a1.get("head"):
        errors.append("a1_a3_lineage_mismatch")
    if a4.get("audits_agent2_partial_composite") is not False:
        errors.append("a4_composite_scope_laundered")
    if a4.get("audits_agent3_partial_nonlinear_mean") is not False:
        errors.append("a4_mean_scope_laundered")

    for key in (
        "current_cartesian_leading_velocity_materialized_through_xh",
        "current_partial_leading_plus_oscillatory_velocity_materialized_through_xh",
        "current_partial_mixed_nonlinear_mean_materialized",
        "current_partial_quadratic_nonlinear_mean_materialized",
        "current_partial_aggregate_nonlinear_mean_materialized",
    ):
        if truth.get(key) is not True:
            errors.append(f"missing_narrow_truth_{key}")

    for key in (
        "velocity_beyond_xh_materialized",
        "outer_global_leading_velocity_materialized",
        "global_compact_support_completed",
        "independent_a4_composite_audit_available",
        "independent_a4_nonlinear_mean_audit_available",
        "radial_inverse_performed_on_current_partial_mean",
        "scoped_nonlinear_mean_authorized_as_correction_target",
        "real_agent3_ns_correction_velocity_materialized",
        "matched_cartesian_pressure_materialized",
        "cartesian_matched_pressure_gradient_materialized",
        "restricted_forcing_materialized",
        "complete_ns_defect_materialized",
        "complete_candidate_api_ready",
        "canonical_whole_domain_divergence_l2_assessed",
        "same_protocol_full_ns_residual_available",
        "same_protocol_st006_comparison_available_now",
        "current_candidate_eligible_for_full_ns_validation",
        "residual_reduction_claimed",
        "scientific_admission",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        if truth.get(key) is not False:
            errors.append(f"premature_promotion_{key}")

    if readiness != READINESS:
        errors.append("readiness_promoted")
    if gate != FINAL_GATE:
        errors.append("final_gate_drift")
    if science.get("residual_defined_free_forcing_forbidden") is not True:
        errors.append("free_forcing_firewall_removed")
    if a3.get("radial_inverse_performed") is not False or a3.get("complete_ns_defect") is not False:
        errors.append("a3_scope_laundered")
    if a3.get("authorized_as_correction_target") is not False:
        errors.append("a3_correction_authorization_laundered")
    return errors


def write_contract(path: str | Path, exact_head: str) -> dict[str, Any]:
    payload = build_contract(exact_head)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return payload


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exact-head")
    parser.add_argument("--out", type=Path)
    parser.add_argument("--check", type=Path)
    args = parser.parse_args(argv)
    if args.check is not None:
        payload = json.loads(args.check.read_text(encoding="utf-8"))
        errors = validate_contract(payload)
        if errors:
            raise SystemExit("invalid A5 contract: " + ", ".join(errors))
        return 0
    if args.exact_head is None or args.out is None:
        parser.error("--exact-head and --out are required unless --check is used")
    write_contract(args.out, args.exact_head)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
