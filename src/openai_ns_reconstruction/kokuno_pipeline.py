"""Kokuno integration preflight and reproducible receipt.

This module is the integration seam for the current Agent 1--4 deliveries.  It
intentionally stops before a full Kokuno 3-D leading field exists.  The existing
capped bipolar candidate is used only as a temporary engineering bridge for the
Agent-3 real phase-mean defect measurement; it is not relabeled as the Kokuno
leading field and cannot make the final velocity or PDE gates pass.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

from .kokuno_complete_curl import KokunoCompleteCurlCorrection
from .kokuno_independent_preflight import REGISTERED_PDE_THRESHOLD, run_independent_preflight
from .kokuno_leading_axis_profiles import KokunoLeadingAxisProfile
from .kokuno_mean_defect import generate_current_candidate_report


SCHEMA = "kokuno-agent5-integration-preflight-v1"
UPSTREAM_HEADS = {
    "agent1_leading_axis_profile": "357f3f05af2bd20d0c092b62531cc890d370c82f",
    "agent2_complete_curl": "7e3a21e8f3b1950e28490bc35124571f285abe8a",
    "agent3_real_mean_defect": "98acbcde149f71faea8f3367760b7f7dd71fcbe2",
    "agent4_independent_preflight": "9c8e2a33baa0b8476f29e99e992cfbebe9bc04fa",
}
BRIDGE_ROLE = "temporary_existing_constrained_candidate_not_kokuno_leading"


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _receipt_sha(payload: dict[str, Any]) -> str:
    bare = dict(payload)
    bare.pop("sha256", None)
    return hashlib.sha256(_canonical_json(bare).encode("utf-8")).hexdigest()


def _preflight_summary(report: dict[str, Any]) -> dict[str, Any]:
    finest_step = min(float(step) for step in report["steps"])
    finest = [row for row in report["rows"] if float(row["step"]) == finest_step]
    return {
        "task_id": report["task_id"],
        "held_out_seed": report["held_out_seed"],
        "point_count": report["point_count"],
        "steps": report["steps"],
        "times": report["times"],
        "structural_preflight_passed": report["structural_preflight_passed"],
        "outside_support_max": report["outside_support_max"],
        "finest_step": finest_step,
        "finest_curl_error_rms_max_over_times": max(float(row["curl_error_rms"]) for row in finest),
        "finest_divergence_rms_max_over_times": max(float(row["divergence_rms"]) for row in finest),
        "registered_full_ns_gate": report["registered_full_ns_gate"],
        "pde_validated": report["pde_validated"],
    }


def build_current_pipeline_receipt(
    *,
    candidate_path: str | Path = "artifacts/bipolar_joint_capped/candidate.json",
    output_dir: str | Path = "artifacts/kokuno_agent5_pipeline",
    mean_defect_points: int = 32,
) -> dict[str, Any]:
    """Run the currently available Kokuno integration stages once.

    The resulting receipt records executable source/profile, oscillatory,
    mean-defect, and independent-preflight evidence while fail-closing the
    still-missing full leading/correction/PDE stages.
    """
    if mean_defect_points <= 0:
        raise ValueError("mean_defect_points must be positive")
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    leading = KokunoLeadingAxisProfile()
    leading_path = leading.save_json(output_dir / "leading_axis_profile.json")
    correction = KokunoCompleteCurlCorrection()

    mean_report_path = output_dir / "mean_defect_report.json"
    mean_report = generate_current_candidate_report(
        candidate_path=candidate_path,
        output=mean_report_path,
        point_count=mean_defect_points,
    )
    preflight = run_independent_preflight(correction)

    receipt: dict[str, Any] = {
        "schema": SCHEMA,
        "upstream_heads": dict(UPSTREAM_HEADS),
        "source": {
            "reader": "KokunoYumeto corrected 208-page reconstruction",
            "edition": "2026.09.09-consolidated",
            "doi": "10.5281/zenodo.22678406",
            "truth_boundary": "source-guided reconstruction; not paper-exact",
        },
        "leading": {
            "axis_profile_path": str(leading_path),
            "axis_profile_sha256": leading.sha256,
            "axis_profile_payload": leading.to_payload(),
            "axis_profile_ready": True,
            "full_3d_leading_ready": False,
            "coordinate_mapping_status": leading.to_payload()["coordinates"]["mapping_status"],
        },
        "oscillatory": {
            "metadata": correction.metadata(),
            "structural_preflight": _preflight_summary(preflight),
            "ready": True,
        },
        "mean_radial_correction": {
            "mean_defect_report_path": str(mean_report_path),
            "bridge_candidate_sha256": mean_report["candidate_sha256"],
            "mean_defect_rows": mean_report["rows"],
            "projector_refinement": mean_report["projector_refinement"],
            "quadratic_scaling_check": mean_report["quadratic_scaling_check"],
            "real_candidate_defect_consumed": True,
            "mean_defect_ready": True,
            "radial_stress_inverse_applied": False,
            "finite_correction_cycle_run": False,
            "correction_ready": False,
        },
        "bridge_base": {
            "role": BRIDGE_ROLE,
            "candidate_path": str(candidate_path),
            "candidate_sha256": mean_report["candidate_sha256"],
            "may_be_promoted_to_kokuno_leading": False,
        },
        "target_api": {
            "signature": "velocity(x,y,z,t)->[...,3] Cartesian [u,v,w]",
            "python_export_required": True,
            "matlab_export_required": True,
            "current_full_kokuno_velocity_available": False,
        },
        "registered_pde_gate": {
            "held_out_normalized_ns_residual_threshold": REGISTERED_PDE_THRESHOLD,
            "threshold_changed": False,
            "assessed_on_full_kokuno_candidate": False,
            "passed": False,
        },
        "status": {
            "leading_axis_profile_ready": True,
            "leading_ready": False,
            "oscillatory_ready": True,
            "mean_defect_ready": True,
            "correction_ready": False,
            "independent_structural_preflight_ready": True,
            "velocity_export_ready": False,
            "pde_validated": False,
        },
        "blockers": [
            "reconcile the source native q/eta convention with the repository Eq45 convention",
            "materialize the full 3-D Kokuno leading velocity through the public velocity API",
            "apply the typed mean/radial inverse and run at least one finite correction cycle",
            "attach compatible pressure and fixed/restricted forcing contracts to the composite candidate",
            "run Agent-4 held-out full normalized NS residual validation at the unchanged 1e-3 gate",
            "only after a frozen passing candidate exists, export the same identity to MATLAB/Python visualization",
        ],
        "truth_boundary": {
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
            "visual_success_is_pde_validation": False,
        },
    }
    receipt["sha256"] = _receipt_sha(receipt)
    (output_dir / "pipeline_receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    validate_pipeline_receipt(receipt)
    return receipt


def validate_pipeline_receipt(receipt: dict[str, Any]) -> bool:
    if not isinstance(receipt, dict) or receipt.get("schema") != SCHEMA:
        raise ValueError("unsupported Kokuno integration receipt schema")
    if receipt.get("upstream_heads") != UPSTREAM_HEADS:
        raise ValueError("upstream Kokuno agent identity drift")
    if receipt.get("sha256") != _receipt_sha(receipt):
        raise ValueError("Kokuno integration receipt SHA mismatch")

    leading = receipt["leading"]
    profile = KokunoLeadingAxisProfile.from_payload(leading["axis_profile_payload"])
    if profile.sha256 != leading["axis_profile_sha256"]:
        raise ValueError("leading-axis profile identity drift")
    if leading["full_3d_leading_ready"] is not False:
        raise ValueError("full 3-D Kokuno leading field is not yet available")
    if leading["coordinate_mapping_status"] != "pending_explicit_convention_reconciliation":
        raise ValueError("coordinate-convention blocker was laundered away")

    oscillatory = receipt["oscillatory"]
    parameters = oscillatory["metadata"]["parameters"]
    correction = KokunoCompleteCurlCorrection(**parameters)
    if correction.metadata()["family"] != oscillatory["metadata"]["family"]:
        raise ValueError("oscillatory correction family drift")
    structural = oscillatory["structural_preflight"]
    if structural["structural_preflight_passed"] is not True:
        raise ValueError("independent oscillatory structural preflight is not passing")
    if structural["registered_full_ns_gate"]["assessed"] is not False:
        raise ValueError("oscillatory preflight cannot assess the full NS gate")
    if structural["pde_validated"] is not False:
        raise ValueError("oscillatory structural preflight cannot validate PDE")

    mean = receipt["mean_radial_correction"]
    if mean["real_candidate_defect_consumed"] is not True or mean["mean_defect_ready"] is not True:
        raise ValueError("real Agent-3 mean defect must remain bound")
    if mean["radial_stress_inverse_applied"] or mean["finite_correction_cycle_run"] or mean["correction_ready"]:
        raise ValueError("mean/radial correction stage is not complete")
    if receipt["bridge_base"]["role"] != BRIDGE_ROLE or receipt["bridge_base"]["may_be_promoted_to_kokuno_leading"]:
        raise ValueError("temporary bridge base was relabeled as Kokuno leading data")

    gate = receipt["registered_pde_gate"]
    if float(gate["held_out_normalized_ns_residual_threshold"]) != 1.0e-3:
        raise ValueError("registered 1e-3 PDE gate changed")
    if gate["threshold_changed"] or gate["assessed_on_full_kokuno_candidate"] or gate["passed"]:
        raise ValueError("full Kokuno PDE gate cannot be promoted yet")

    expected_status = {
        "leading_axis_profile_ready": True,
        "leading_ready": False,
        "oscillatory_ready": True,
        "mean_defect_ready": True,
        "correction_ready": False,
        "independent_structural_preflight_ready": True,
        "velocity_export_ready": False,
        "pde_validated": False,
    }
    if receipt["status"] != expected_status:
        raise ValueError("Kokuno pipeline status drift")
    if receipt["truth_boundary"] != {
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
        "visual_success_is_pde_validation": False,
    }:
        raise ValueError("Kokuno truth boundary drift")
    return True


def load_pipeline_receipt(path: str | Path) -> dict[str, Any]:
    receipt = json.loads(Path(path).read_text(encoding="utf-8"))
    validate_pipeline_receipt(receipt)
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the current Kokuno Agent 1-4 integration preflight")
    parser.add_argument("--candidate", default="artifacts/bipolar_joint_capped/candidate.json")
    parser.add_argument("--output-dir", default="artifacts/kokuno_agent5_pipeline")
    parser.add_argument("--mean-defect-points", type=int, default=32)
    args = parser.parse_args()
    receipt = build_current_pipeline_receipt(
        candidate_path=args.candidate,
        output_dir=args.output_dir,
        mean_defect_points=args.mean_defect_points,
    )
    print(json.dumps({
        "sha256": receipt["sha256"],
        "status": receipt["status"],
        "registered_pde_gate": receipt["registered_pde_gate"],
        "blockers": receipt["blockers"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
