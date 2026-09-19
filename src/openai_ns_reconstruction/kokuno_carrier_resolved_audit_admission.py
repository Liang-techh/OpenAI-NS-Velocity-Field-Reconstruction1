"""Agent-3 admission for the independently audited carrier-resolved oscillatory field.

Agent 4 PR #553 audited Agent 2 PR #551 through the frozen public black-box
protocol.  That execution independently clears support and the phase-mean
covariance-rank seam, but it still rejects the field on divergence RMS and
refinement stability.  This module pins the exact execution/receipt identity,
recomputes the frozen scientific inequalities, and exposes only the truthful
partial promotion needed by downstream mean-correction work.

In particular, covariance availability alone never authorizes requested-stress
materialization, ``Delta C / epsilon`` inversion, or a finite correction cycle.
A future overall oscillatory PASS must be independently audited and receive a
new Agent-3 execution-identity pin.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

TASK = "KOKUNO-A3-CARRIER-RESOLVED-AUDIT-IDENTITY-037"
SCHEMA = "kokuno-agent3-carrier-resolved-audit-admission-v1"
EXPECTED_AUDIT_SCHEMA = "kokuno-agent4-carrier-resolved-independent-audit-v1"

ACTUAL_AUDITED_AGENT2_PR = 551
ACTUAL_AUDITED_AGENT2_HEAD = "92852046e6e1ec0a5889b53f989df3ea5d79cb3b"
ACTUAL_AUDITOR_AGENT4_PR = 553
ACTUAL_AUDITOR_AGENT4_HEAD = "2a18a87db0ef8c6682f10fc5a6355c90b4168c7c"
PROTOCOL_ORIGIN_AGENT4_HEAD = "9498209376fd9dc5abf228248967958759e9f64a"
AUDIT_WORKFLOW_RUN_ID = 35419518424
AUDIT_ARTIFACT_ID = 10576976589
AUDIT_ARTIFACT_ZIP_DIGEST = (
    "sha256:3ee35376dd1fa8cb1417639ebe40513b74cc6d8fe686d3fe380fa669ac6239cf"
)
AUDIT_RECEIPT_FILENAME = "kokuno_agent4_carrier_resolved_independent_audit.json"
AUDIT_RECEIPT_RAW_SHA256 = "c316f8045cdc4388396e5c76a16e004d5e8055903c528398a4a7406c7be7d127"
AUDIT_RECEIPT_CANONICAL_SHA256 = "3ebe93ca8f49240d08d8f8f97142e6bdd3e9fc6ed73171519a8b4d59854ce114"

SEED = 9173241
FD4_STEPS = (0.02, 0.01, 0.005)
PHASE_RESOLUTIONS = (32, 64, 128)
FROZEN_GUARDS = {
    "finest_relative_divergence_rms": 2.0e-5,
    "finest_relative_divergence_max": 1.0e-4,
    "minimum_divergence_refinement_ratio": 3.0,
    "minimum_covariance_rank_ratio": 2.0e-2,
    "maximum_covariance_resolution_drift": 1.0e-6,
    "axis_near_absolute_max": 1.0e-14,
    "radial_exterior_absolute_max": 1.0e-14,
    "project_axial_exterior_absolute_max": 1.0e-12,
    "minimum_nontrivial_velocity_rms": 1.0e-8,
    "minimum_parameter_perturbation_relative_change": 1.0e-4,
    "minimum_divergence_mutation_relative_rms": 5.0e-2,
    "maximum_duplicated_covariance_rank_ratio": 1.0e-8,
}


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _canonical_payload_bytes(payload: Mapping[str, Any]) -> bytes:
    return json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    ).encode("utf-8")


def _parse_and_bind_receipt(raw_receipt: bytes) -> Mapping[str, Any]:
    if not isinstance(raw_receipt, bytes):
        raise TypeError("raw_receipt must be bytes")
    if _sha256_bytes(raw_receipt) != AUDIT_RECEIPT_RAW_SHA256:
        raise ValueError("audit receipt bytes do not match pinned Agent-4 #553 artifact member")
    try:
        payload = json.loads(raw_receipt.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError("audit receipt is not valid UTF-8 JSON") from exc
    if not isinstance(payload, Mapping):
        raise ValueError("audit receipt JSON must be an object")
    if _sha256_bytes(_canonical_payload_bytes(payload)) != AUDIT_RECEIPT_CANONICAL_SHA256:
        raise ValueError("audit receipt canonical payload does not match pinned execution")
    if payload.get("schema") != EXPECTED_AUDIT_SCHEMA:
        raise ValueError("pinned audit schema changed")
    if payload.get("audited_agent2_head") != ACTUAL_AUDITED_AGENT2_HEAD:
        raise ValueError("pinned audit Agent-2 candidate identity changed")
    if payload.get("frozen_protocol_origin_agent4_head") != PROTOCOL_ORIGIN_AGENT4_HEAD:
        raise ValueError("pinned audit protocol-origin Agent-4 head changed")
    return payload


def _require_frozen_protocol(audit: Mapping[str, Any]) -> None:
    protocol = audit["protocol"]
    if audit.get("seed") != SEED:
        raise ValueError("independent audit seed changed")
    if tuple(protocol.get("fd4_steps", ())) != FD4_STEPS:
        raise ValueError("independent FD4 ladder changed")
    if tuple(protocol.get("phase_resolutions", ())) != PHASE_RESOLUTIONS:
        raise ValueError("phase-mean covariance resolutions changed")
    if protocol.get("heldout_offgrid_points") != 24:
        raise ValueError("held-out point count changed")
    if protocol.get("candidate_internal_derivative_oracle_used") is not False:
        raise ValueError("candidate-internal derivative oracle is forbidden")
    if protocol.get("training_tensor_or_loss_used") is not False:
        raise ValueError("training tensor/loss reuse is forbidden")
    if protocol.get("pressure_or_forcing_fit") is not False:
        raise ValueError("pressure/forcing fit is forbidden in this audit")
    if protocol.get("scientific_guards_changed_from_agent4_543") is not False:
        raise ValueError("scientific guards changed from the frozen protocol")
    if audit.get("guards_frozen_before_actions") != FROZEN_GUARDS:
        raise ValueError("frozen scientific guard values changed")


def _recompute_scientific_verdict(audit: Mapping[str, Any]) -> dict[str, Any]:
    guards = FROZEN_GUARDS
    div = audit["divergence"]
    div_rows = div["rows"]
    if tuple(float(row["step"]) for row in div_rows) != FD4_STEPS:
        raise ValueError("divergence row ladder changed")
    refinement = tuple(float(v) for v in div["relative_rms_refinement_ratios"])
    if len(refinement) != 2:
        raise ValueError("divergence refinement ladder must have two ratios")
    finest = div_rows[-1]
    divergence_rms_passed = bool(
        float(finest["relative_rms"]) <= guards["finest_relative_divergence_rms"]
    )
    divergence_max_passed = bool(
        float(finest["relative_max"]) <= guards["finest_relative_divergence_max"]
    )
    divergence_refinement_passed = bool(
        min(refinement) >= guards["minimum_divergence_refinement_ratio"]
    )
    local_divergence_passed = bool(
        divergence_rms_passed and divergence_max_passed and divergence_refinement_passed
    )

    cov = audit["phase_mean_covariance"]
    cov_rows = cov["rows"]
    if tuple(int(row["phase_resolution"]) for row in cov_rows) != PHASE_RESOLUTIONS:
        raise ValueError("covariance phase-resolution ladder changed")
    min_covariance_rank_ratio = min(float(row["rank_ratio"]) for row in cov_rows)
    max_covariance_drift = max(float(v) for v in cov["relative_drift_to_finest"])
    duplicated_rank_ratio = float(cov["duplicated_first_column_rank_ratio"])
    covariance_passed = bool(
        min_covariance_rank_ratio >= guards["minimum_covariance_rank_ratio"]
        and max_covariance_drift <= guards["maximum_covariance_resolution_drift"]
        and duplicated_rank_ratio <= guards["maximum_duplicated_covariance_rank_ratio"]
    )

    support = audit["support"]
    support_passed = bool(
        float(support["axis_near_absolute_max"]) <= guards["axis_near_absolute_max"]
        and float(support["radial_exterior_absolute_max"])
        <= guards["radial_exterior_absolute_max"]
        and float(support["project_axial_exterior_absolute_max"])
        <= guards["project_axial_exterior_absolute_max"]
    )
    nontriviality_passed = bool(
        float(audit["nontriviality"]["heldout_velocity_rms"])
        >= guards["minimum_nontrivial_velocity_rms"]
    )
    parameter_perturbation_passed = bool(
        min(float(v) for k, v in audit["parameter_perturbation"].items() if k.endswith("relative_change"))
        >= guards["minimum_parameter_perturbation_relative_change"]
    )
    divergence_mutation_passed = bool(
        float(audit["divergence_mutation"]["relative_rms"])
        >= guards["minimum_divergence_mutation_relative_rms"]
    )

    failed: list[str] = []
    if not divergence_rms_passed:
        failed.append("finest_relative_divergence_rms")
    if not divergence_max_passed:
        failed.append("finest_relative_divergence_max")
    if not divergence_refinement_passed:
        failed.append("minimum_divergence_refinement_ratio")
    if not covariance_passed:
        failed.append("phase_mean_covariance_rank_or_resolution")
    if not support_passed:
        failed.append("project_support")
    if not nontriviality_passed:
        failed.append("nontriviality")
    if not parameter_perturbation_passed:
        failed.append("parameter_perturbation")
    if not divergence_mutation_passed:
        failed.append("divergence_mutation_detection")

    overall = bool(
        local_divergence_passed
        and covariance_passed
        and support_passed
        and nontriviality_passed
        and parameter_perturbation_passed
        and divergence_mutation_passed
    )
    return {
        "finest_relative_divergence_rms": float(finest["relative_rms"]),
        "finest_relative_divergence_max": float(finest["relative_max"]),
        "divergence_refinement_ratios": refinement,
        "divergence_rms_passed": divergence_rms_passed,
        "divergence_max_passed": divergence_max_passed,
        "divergence_refinement_passed": divergence_refinement_passed,
        "local_divergence_passed": local_divergence_passed,
        "minimum_covariance_rank_ratio": min_covariance_rank_ratio,
        "maximum_covariance_resolution_drift": max_covariance_drift,
        "duplicated_covariance_rank_ratio": duplicated_rank_ratio,
        "covariance_rank_preflight_passed": covariance_passed,
        "project_support_preflight_passed": support_passed,
        "nontriviality_passed": nontriviality_passed,
        "parameter_perturbation_passed": parameter_perturbation_passed,
        "divergence_mutation_passed": divergence_mutation_passed,
        "failed_guards": tuple(failed),
        "public_oscillatory_preflight_passed": overall,
    }


def evaluate_carrier_resolved_audit_admission(raw_receipt: bytes) -> dict[str, Any]:
    """Bind #553 and expose covariance PASS without opening correction ingestion."""
    audit = _parse_and_bind_receipt(raw_receipt)
    _require_frozen_protocol(audit)
    scientific = _recompute_scientific_verdict(audit)

    declared_pairs = (
        ("local_divergence_passed", audit["divergence"]["local_divergence_passed"]),
        ("covariance_rank_preflight_passed", audit["phase_mean_covariance"]["covariance_rank_preflight_passed"]),
        ("project_support_preflight_passed", audit["project_support_preflight_passed"]),
        ("public_oscillatory_preflight_passed", audit["public_oscillatory_preflight_passed"]),
    )
    for key, declared in declared_pairs:
        if bool(declared) is not bool(scientific[key]):
            raise ValueError(f"Agent-4 declared {key} disagrees with recomputed frozen guards")

    truth = audit["truth_boundary"]
    if float(truth["formal_momentum_normalized_max_l2_gate"]) != 1.0e-3:
        raise ValueError("formal momentum gate changed")
    if float(truth["formal_divergence_max_l2_gate"]) != 1.0e-5:
        raise ValueError("formal divergence gate changed")
    for key in (
        "formal_full_domain_pde_gate_assessed",
        "global_leading_velocity_pressure_available",
        "heldout_ns_residual_assessed",
        "materialized_correction_available",
        "pde_validated",
    ):
        if truth.get(key) is not False:
            raise ValueError(f"pinned component audit cannot promote {key}")

    if scientific["covariance_rank_preflight_passed"] is not True:
        raise ValueError("pinned #553 receipt unexpectedly lost the independently cleared covariance seam")
    if scientific["project_support_preflight_passed"] is not True:
        raise ValueError("pinned #553 receipt unexpectedly lost the independently cleared support seam")
    if scientific["public_oscillatory_preflight_passed"] is not False:
        raise ValueError("pinned #553 receipt unexpectedly became an overall oscillatory PASS")
    if scientific["failed_guards"] != (
        "finest_relative_divergence_rms",
        "minimum_divergence_refinement_ratio",
    ):
        raise ValueError("pinned #553 scientific failure boundary changed")

    return {
        "task": TASK,
        "schema": SCHEMA,
        "execution_identity_bound": True,
        "actual_audited_agent2_pr": ACTUAL_AUDITED_AGENT2_PR,
        "actual_audited_agent2_head": ACTUAL_AUDITED_AGENT2_HEAD,
        "actual_auditor_agent4_pr": ACTUAL_AUDITOR_AGENT4_PR,
        "actual_auditor_agent4_head": ACTUAL_AUDITOR_AGENT4_HEAD,
        "protocol_origin_agent4_head": PROTOCOL_ORIGIN_AGENT4_HEAD,
        "audit_workflow_run_id": AUDIT_WORKFLOW_RUN_ID,
        "audit_artifact_id": AUDIT_ARTIFACT_ID,
        "audit_artifact_zip_digest": AUDIT_ARTIFACT_ZIP_DIGEST,
        "audit_receipt_filename": AUDIT_RECEIPT_FILENAME,
        "audit_receipt_raw_sha256": AUDIT_RECEIPT_RAW_SHA256,
        "audit_receipt_canonical_sha256": AUDIT_RECEIPT_CANONICAL_SHA256,
        **scientific,
        "independent_candidate_covariance_rank_preflight_passed": True,
        "candidate_second_covariance_direction_independently_observed": True,
        "covariance_input_structurally_available": True,
        "correction_ingest_allowed": False,
        "same_cycle_requested_stress_materialization_allowed": False,
        "candidate_finite_head_mean_debt_materialization_allowed": False,
        "requested_stress_actual_state_values_materialized": False,
        "finite_head_mean_debt_materialized": False,
        "real_candidate_defect_consumed": False,
        "signed_mean_inverse_input_ready": False,
        "public_velocity_correction_materialized": False,
        "finite_correction_cycle_rerun_allowed": False,
        "finite_correction_cycle_run": False,
        "heldout_ns_residual_assessed": False,
        "residual_reduction_claimed": False,
        "pde_validated": False,
        "surrogate_defect_used": False,
        "future_overall_pass_requires_new_agent3_identity_pin": True,
        "scope": (
            "exact Agent-4 #553 carrier-resolved execution: independently clears support and "
            "phase-mean covariance rank, but divergence RMS/refinement still reject the "
            "oscillatory handoff; no mean defect or correction is materialized"
        ),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    report = evaluate_carrier_resolved_audit_admission(args.audit.read_bytes())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
