"""Agent-3 admission for Agent-4's independent compact radial-stress audit.

Agent 4 PR #595 independently rebuilt the admitted oscillatory self-defect with
Cartesian FD6, reconstructed the compact moment-complement radial stress with a
separate Simpson/cumulative-Simpson implementation, and checked the radial
operator with a new FD6 derivative. This module pins that exact execution and
recomputes every frozen local guard from the exact artifact receipt.

A PASS here establishes that the compact/radial inverse used by Agent 3 has an
independent black-box cross-audit on the real oscillatory component. It does
not create the missing Agent-1 leading/cross defect, does not materialize a
full same-cycle requestedStress or finite-head target, and does not authorize a
signed inverse or finite correction-cycle run yet.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping

TASK = "KOKUNO-A3-INDEPENDENT-RADIAL-STRESS-ADMISSION-042"
SCHEMA = "kokuno-agent3-independent-radial-stress-admission-v1"
EXPECTED_AUDIT_SCHEMA = "kokuno-a4-independent-radial-stress-audit-v1"

PARENT_AGENT3_PR = 589
PARENT_AGENT3_HEAD = "2271c52412922c8255cede8aee9a01e68ba55a55"
ACTUAL_AUDITOR_AGENT4_PR = 595
ACTUAL_AUDITOR_AGENT4_HEAD = "19296acad4f84ba05d3b095bb8c01bf5a2c95891"
ADMITTED_AGENT2_PR = 561
ADMITTED_AGENT2_HEAD = "732800ce4990464b49c8aa32d0dff4580f6684d4"
INDEPENDENT_AGENT4_VELOCITY_PR = 563
INDEPENDENT_AGENT4_VELOCITY_HEAD = "9b0f86012c53fa8e32a19f766dbc150931870425"
INDEPENDENT_AGENT4_DT_PR = 580
INDEPENDENT_AGENT4_DT_HEAD = "6f7edb4d66dcb162a8510a40e4227f48f3c57e28"
AUDIT_WORKFLOW_RUN_ID = 35431599990
AUDIT_ARTIFACT_ID = 10579939940
AUDIT_ARTIFACT_ZIP_DIGEST = (
    "sha256:295980b23ff31877f341bbf3d7866718678dcd716fc4bd5d4efb3d0ff8c3845a"
)
AUDIT_RECEIPT_FILENAME = "agent4_595_independent_radial_stress_audit.json"
AUDIT_RECEIPT_RAW_SHA256 = "f1c8448d0eb640ed6c3c1eb9f105412aa25a0d4ff5639930d6224ad3eb198259"
AUDIT_RECEIPT_CANONICAL_SHA256 = "794482cfd48fb7369d2a3150c028ec64ac482fbdea93f02ea866a39ed11dfb27"

RADIAL_COUNTS = (25, 49, 97)
ANGULAR_COUNT = 16
TIME = 0.5
Z = 0.08
NU = 0.01
CARTESIAN_FD6_STEP = 0.003
FROZEN_GUARDS = {
    "edge_relative_max": 1.0e-8,
    "finest_relative_max_max": 5.0e-2,
    "finest_relative_rms_max": 2.0e-2,
    "minimum_refinement_ratio": 2.0,
    "moment_complement_relative_max": 1.0e-10,
    "project_divergence_max_l2_gate": 1.0e-5,
    "project_momentum_max_l2_gate": 1.0e-3,
    "sign_flip_mutation_relative_rms_min": 0.5,
}


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _canonical_payload_bytes(payload: Mapping[str, Any]) -> bytes:
    return json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")


def _approx_equal(
    left: float,
    right: float,
    *,
    rel: float = 1.0e-12,
    abs_: float = 1.0e-15,
) -> bool:
    return math.isclose(float(left), float(right), rel_tol=rel, abs_tol=abs_)


def _parse_and_bind_receipt(raw_receipt: bytes) -> Mapping[str, Any]:
    if not isinstance(raw_receipt, bytes):
        raise TypeError("raw_receipt must be bytes")
    if _sha256_bytes(raw_receipt) != AUDIT_RECEIPT_RAW_SHA256:
        raise ValueError("audit receipt bytes do not match pinned Agent-4 #595 artifact member")
    try:
        payload = json.loads(raw_receipt.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError("audit receipt is not valid UTF-8 JSON") from exc
    if not isinstance(payload, Mapping):
        raise ValueError("audit receipt JSON must be an object")
    if _sha256_bytes(_canonical_payload_bytes(payload)) != AUDIT_RECEIPT_CANONICAL_SHA256:
        raise ValueError("audit receipt canonical payload does not match pinned execution")
    if payload.get("schema") != EXPECTED_AUDIT_SCHEMA:
        raise ValueError("pinned independent radial-stress audit schema changed")
    provenance = payload.get("provenance")
    if not isinstance(provenance, Mapping):
        raise ValueError("pinned audit provenance is missing")
    expected = {
        "parent_agent3_pr": PARENT_AGENT3_PR,
        "parent_agent3_head": PARENT_AGENT3_HEAD,
        "admitted_agent2_pr": ADMITTED_AGENT2_PR,
        "admitted_agent2_head": ADMITTED_AGENT2_HEAD,
        "independent_agent4_velocity_pr": INDEPENDENT_AGENT4_VELOCITY_PR,
        "independent_agent4_velocity_head": INDEPENDENT_AGENT4_VELOCITY_HEAD,
        "independent_agent4_dt_pr": INDEPENDENT_AGENT4_DT_PR,
        "independent_agent4_dt_head": INDEPENDENT_AGENT4_DT_HEAD,
    }
    for key, value in expected.items():
        if provenance.get(key) != value:
            raise ValueError(f"pinned audit provenance changed: {key}")
    if provenance.get("velocity_provider") != (
        "openai_ns_reconstruction.kokuno_public_z_pullback_velocity:velocity_osc"
    ):
        raise ValueError("pinned public velocity provider changed")
    if provenance.get("time_derivative_provider") != (
        "openai_ns_reconstruction.kokuno_public_oscillatory_time_derivative:velocity_osc_dt"
    ):
        raise ValueError("pinned public time-derivative provider changed")
    return payload


def _require_frozen_protocol(audit: Mapping[str, Any]) -> None:
    state = audit.get("frozen_state")
    if not isinstance(state, Mapping):
        raise ValueError("pinned audit frozen_state is missing")
    if tuple(int(v) for v in state.get("radial_counts", ())) != RADIAL_COUNTS:
        raise ValueError("pinned radial-resolution ladder changed")
    if int(state.get("angular_count", -1)) != ANGULAR_COUNT:
        raise ValueError("pinned angular count changed")
    if not _approx_equal(float(state.get("time", math.nan)), TIME):
        raise ValueError("pinned time changed")
    if not _approx_equal(float(state.get("z", math.nan)), Z):
        raise ValueError("pinned z changed")
    if not _approx_equal(float(state.get("nu", math.nan)), NU):
        raise ValueError("pinned viscosity changed")
    if not _approx_equal(
        float(state.get("cartesian_fd6_step", math.nan)), CARTESIAN_FD6_STEP
    ):
        raise ValueError("pinned Cartesian FD6 step changed")

    guards = audit.get("frozen_guards")
    if not isinstance(guards, Mapping):
        raise ValueError("pinned audit frozen_guards is missing")
    if set(guards) != set(FROZEN_GUARDS):
        raise ValueError("pinned frozen guard set changed")
    for key, expected in FROZEN_GUARDS.items():
        if not _approx_equal(float(guards[key]), expected):
            raise ValueError(f"pinned frozen guard changed: {key}")


def _recompute_channel(audit: Mapping[str, Any], channel: str) -> dict[str, Any]:
    levels = audit.get("levels")
    if not isinstance(levels, list) or len(levels) != len(RADIAL_COUNTS):
        raise ValueError("pinned audit levels changed")
    counts = tuple(int(level.get("radial_count", -1)) for level in levels)
    if counts != RADIAL_COUNTS:
        raise ValueError("pinned audit level radial counts changed")

    rows: list[Mapping[str, Any]] = []
    for level in levels:
        row = level.get(channel)
        if not isinstance(row, Mapping):
            raise ValueError(f"pinned audit channel missing: {channel}")
        rows.append(row)

    errors = tuple(float(row["operator_relative_rms"]) for row in rows)
    ratios = tuple(
        errors[i] / max(errors[i + 1], 1.0e-300)
        for i in range(len(errors) - 1)
    )
    finest = rows[-1]
    result = {
        "relative_rms_by_radial_count": errors,
        "refinement_ratios": ratios,
        "finest_relative_rms": float(finest["operator_relative_rms"]),
        "finest_relative_max": float(finest["operator_relative_max"]),
        "finest_moment_complement_relative": float(finest["moment_complement_relative"]),
        "finest_edge_relative": float(finest["edge_relative"]),
        "finest_sign_flip_mutation_relative_rms": float(
            finest["sign_flip_mutation_relative_rms"]
        ),
    }

    convergence = audit.get("convergence")
    declared = convergence.get(channel) if isinstance(convergence, Mapping) else None
    if not isinstance(declared, Mapping):
        raise ValueError(f"pinned audit convergence summary missing: {channel}")
    for key in (
        "finest_relative_rms",
        "finest_relative_max",
        "finest_moment_complement_relative",
        "finest_edge_relative",
        "finest_sign_flip_mutation_relative_rms",
    ):
        if not _approx_equal(float(declared[key]), float(result[key])):
            raise ValueError(f"pinned audit declared {channel}.{key} disagrees with levels")
    declared_errors = tuple(float(v) for v in declared["relative_rms_by_radial_count"])
    declared_ratios = tuple(float(v) for v in declared["refinement_ratios"])
    if len(declared_errors) != len(errors) or any(
        not _approx_equal(a, b) for a, b in zip(declared_errors, errors)
    ):
        raise ValueError(f"pinned audit declared {channel} RMS ladder disagrees with levels")
    if len(declared_ratios) != len(ratios) or any(
        not _approx_equal(a, b) for a, b in zip(declared_ratios, ratios)
    ):
        raise ValueError(f"pinned audit declared {channel} refinement ratios disagree with levels")

    failed: list[str] = []
    if result["finest_relative_rms"] > FROZEN_GUARDS["finest_relative_rms_max"]:
        failed.append("finest_relative_rms")
    if result["finest_relative_max"] > FROZEN_GUARDS["finest_relative_max_max"]:
        failed.append("finest_relative_max")
    if min(result["refinement_ratios"]) < FROZEN_GUARDS["minimum_refinement_ratio"]:
        failed.append("minimum_refinement_ratio")
    if result["finest_moment_complement_relative"] > FROZEN_GUARDS[
        "moment_complement_relative_max"
    ]:
        failed.append("moment_complement_relative")
    if result["finest_edge_relative"] > FROZEN_GUARDS["edge_relative_max"]:
        failed.append("edge_relative")
    if result["finest_sign_flip_mutation_relative_rms"] < FROZEN_GUARDS[
        "sign_flip_mutation_relative_rms_min"
    ]:
        failed.append("sign_flip_mutation_relative_rms")
    result["failed_guards"] = tuple(failed)
    result["passed"] = not failed
    return result


def evaluate_independent_radial_stress_admission(raw_receipt: bytes) -> dict[str, Any]:
    """Pin Agent-4 #595 and admit only the radial-inverse operator seam."""
    audit = _parse_and_bind_receipt(raw_receipt)
    _require_frozen_protocol(audit)
    theta = _recompute_channel(audit, "theta_e2")
    axial = _recompute_channel(audit, "axial_e1")
    failed = tuple(
        [f"theta_e2:{name}" for name in theta["failed_guards"]]
        + [f"axial_e1:{name}" for name in axial["failed_guards"]]
    )

    declared_failed = tuple(str(v) for v in audit.get("failed_guards", ()))
    if declared_failed != failed:
        raise ValueError("Agent-4 declared failed_guards disagrees with recomputed frozen guards")
    declared_pass = audit.get("independent_radial_stress_preflight_passed")
    if bool(declared_pass) is not (not failed):
        raise ValueError("Agent-4 declared radial-stress PASS disagrees with recomputed guards")
    if failed:
        raise ValueError("pinned Agent-4 #595 radial-stress audit no longer passes frozen guards")

    truth = audit.get("truth_boundary")
    if not isinstance(truth, Mapping):
        raise ValueError("pinned audit truth_boundary is missing")
    required_true = (
        "independent_compact_radial_stress_reconstruction_audited",
        "independent_real_oscillatory_self_defect_consumed",
        "public_black_box_time_derivative_consumed",
        "public_black_box_velocity_consumed",
    )
    for key in required_true:
        if truth.get(key) is not True:
            raise ValueError(f"pinned independent audit lost required evidence: {key}")
    required_false = (
        "agent1_leading_cross_terms_included",
        "agent3_compact_stress_constructor_used",
        "agent3_fd4_self_defect_helper_used",
        "agent3_radial_derivative_used",
        "candidate_finite_head_mean_debt_materialized",
        "finite_correction_cycle_run",
        "formal_full_domain_pde_gate_assessed",
        "full_same_cycle_composite_requested_stress_materialized",
        "heldout_ns_momentum_residual_assessed",
        "matched_pressure_included",
        "public_velocity_correction_materialized",
        "restricted_forcing_included",
        "signed_mean_inverse_input_ready",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    )
    for key in required_false:
        if truth.get(key) is not False:
            raise ValueError(f"component radial-stress audit cannot promote {key}")

    return {
        "task": TASK,
        "schema": SCHEMA,
        "execution_identity_bound": True,
        "parent_agent3_pr": PARENT_AGENT3_PR,
        "parent_agent3_head": PARENT_AGENT3_HEAD,
        "actual_auditor_agent4_pr": ACTUAL_AUDITOR_AGENT4_PR,
        "actual_auditor_agent4_head": ACTUAL_AUDITOR_AGENT4_HEAD,
        "admitted_agent2_pr": ADMITTED_AGENT2_PR,
        "admitted_agent2_head": ADMITTED_AGENT2_HEAD,
        "audit_workflow_run_id": AUDIT_WORKFLOW_RUN_ID,
        "audit_artifact_id": AUDIT_ARTIFACT_ID,
        "audit_artifact_zip_digest": AUDIT_ARTIFACT_ZIP_DIGEST,
        "audit_receipt_filename": AUDIT_RECEIPT_FILENAME,
        "audit_receipt_raw_sha256": AUDIT_RECEIPT_RAW_SHA256,
        "audit_receipt_canonical_sha256": AUDIT_RECEIPT_CANONICAL_SHA256,
        "theta_e2": theta,
        "axial_e1": axial,
        "failed_guards": failed,
        "independent_compact_radial_stress_cross_audit_passed": True,
        "compact_radial_stress_operator_independently_admitted": True,
        "full_composite_radial_stress_execution_allowed_when_actual_defect_available": True,
        "real_oscillatory_component_defect_independently_consumed": True,
        "full_same_cycle_composite_requested_stress_materialized": False,
        "agent1_leading_cross_terms_included": False,
        "matched_pressure_included": False,
        "restricted_forcing_included": False,
        "candidate_finite_head_mean_debt_materialized": False,
        "real_full_candidate_defect_consumed": False,
        "signed_mean_inverse_input_ready": False,
        "public_velocity_correction_materialized": False,
        "finite_correction_cycle_rerun_allowed": False,
        "finite_correction_cycle_run": False,
        "heldout_ns_residual_assessed": False,
        "residual_reduction_claimed": False,
        "pde_validated": False,
        "paper_exact": False,
        "blowup_proved": False,
        "surrogate_defect_used": False,
        "scope": (
            "Agent-4 #595 independently validates the compact radial inverse on the real admitted "
            "oscillatory self-defect; the operator seam is admitted for reuse once the actual full "
            "same-cycle composite defect exists, while signed inversion, correction materialization, "
            "finite-cycle execution and PDE claims remain closed"
        ),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    result = evaluate_independent_radial_stress_admission(args.audit.read_bytes())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
