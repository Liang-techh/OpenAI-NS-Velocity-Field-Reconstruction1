"""Kokuno Agent-5 fail-closed routing for the support-localized curl seam.

This integration increment adds no Kokuno mathematics.  It binds Agent 2's
support-localized complete-curl implementation to Agent 4's independent
Cartesian audit and preserves the *failed* preregistered small-radius
divergence guard as a routing blocker.

The local audit is not the formal Navier--Stokes gate.  The project gates stay
held-out normalized momentum max/L2 <= 1e-3 and divergence max/L2 <= 1e-5.
No free residual-defined forcing is introduced and no failed guard is weakened.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

from .kokuno_support_localized_cartesian_independent import build_report as build_agent4_report

SCHEMA = "kokuno-agent5-support-localization-routing-checkpoint-v17"
TASK_ID = "KOKUNO-A5-SUPPORT-LOCALIZATION-REJECTION-ROUTING-017"
BASE_AGENT4_PR = 353
BASE_AGENT4_HEAD = "a9f24388245996eb7b0e2cc66e4847eae2833f3d"
AGENT2_PR = 351
AGENT2_HEAD = "26112140d8a57873d20c6104ea19f32d8b9afa1d"

FIXED_GATES = {
    "held_out_normalized_full_momentum_max": 1.0e-3,
    "held_out_normalized_full_momentum_l2": 1.0e-3,
    "divergence_max": 1.0e-5,
    "divergence_l2": 1.0e-5,
    "changed": False,
}

ST006_REFERENCE = {
    "candidate_sha256": "6b4d84b48ab9dbcd2ee1a1858d3e56ef81523f5864369d7e96c6431fccf107a3",
    "validation_seed": 9172801,
    "held_out_cartesian_points": 4096,
    "validation_times": 6,
    "finest_spatial_step": 0.005,
    "momentum_sampled_max": 0.1082289305112118,
    "volume_l2": 0.10758432876230622,
    "pde_validated": False,
}

UPSTREAM = {
    "agent1": {
        "pr": 345,
        "head_sha": "57df8375f134277ebea509847f6b319aa84acfd7",
        "segmented_leading_velocity_router_ready": True,
        "global_radial_coverage_complete": False,
        "global_matched_pressure_ready": False,
        "leading_ready": False,
        "consumed_in_executable_ancestry": False,
    },
    "agent2": {
        "pr": AGENT2_PR,
        "head_sha": AGENT2_HEAD,
        "dedicated_run": 35322804017,
        "standard_run": 35322804023,
        "support_localized_vector_potential_curl_ready": True,
        "actual_source_partition_values_ready": False,
        "actual_positive_order_source_background_ready": False,
        "public_source_oscillatory_xyz_t_velocity_ready": False,
        "genuinely_independent_second_covariance_column_ready": False,
        "consumed_in_executable_ancestry": True,
    },
    "agent3": {
        "pr": 346,
        "head_sha": "61c0f91dc4b0175cd398a7b27f941ba8b2dd1cd4",
        "dedicated_run": 35320473725,
        "standard_run": 35320473638,
        "artifact_id": 10537575118,
        "artifact_digest": "sha256:290268203086d941a99798a418af2461b373caf351b834f0477dd07ec22d74fd",
        "radial_counts": [33, 65, 129],
        "fine_pair_max_guarded_relative_change": 0.00029353357568215327,
        "finest_missing_relative_vector_rms": 0.2101949905750188,
        "finest_nodes_requiring_second_direction": 99,
        "duplicate_control_rank2_required_nodes": 0,
        "finite_correction_cycle_rerun_allowed": False,
        "consumed_in_executable_ancestry": False,
    },
    "agent4": {
        "pr": BASE_AGENT4_PR,
        "head_sha": BASE_AGENT4_HEAD,
        "dedicated_run": 35323537621,
        "artifact_id": 10537559838,
        "artifact_digest": "sha256:aada7589a6067c91057b447785d0914203f9016c85708bdcf07c936d2ce36f8b",
        "independent_cartesian_support_localization_audit_ready": True,
        "consumed_in_executable_ancestry": True,
    },
    "previous_agent5": {
        "pr": 350,
        "head_sha": "0eb7ca53a93020e419b14d2475c032129bd9d0d5",
        "segmented_leading_frontier_recorded": True,
        "consumed_in_executable_ancestry": False,
    },
}

STATES = {
    "segmented_leading_velocity_router_ready": True,
    "support_localized_vector_potential_curl_ready": True,
    "support_localized_cartesian_preflight_passed": False,
    "support_localized_small_radius_divergence_guard_passed": False,
    "leading_ready": False,
    "oscillatory_ready": True,
    "public_source_oscillatory_xyz_t_velocity_ready": False,
    "genuinely_independent_second_covariance_column_ready": False,
    "correction_ready": False,
    "complete_kokuno_composite_velocity_ready": False,
    "velocity_export_ready": False,
    "formal_full_domain_pde_gate_assessed": False,
    "pde_validated": False,
}


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _sha(payload: dict[str, Any]) -> str:
    return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()


def _validate_agent4_report(report: dict[str, Any]) -> None:
    if report.get("task_id") != "KOKUNO-A4-SUPPORT-LOCALIZED-CARTESIAN-AUDIT-017":
        raise ValueError("unexpected Agent-4 support-localization audit task")
    dependency = report.get("dependency", {})
    if dependency.get("agent2_pr") != AGENT2_PR or dependency.get("exact_head") != AGENT2_HEAD:
        raise ValueError("Agent-4 audit is not bound to the expected Agent-2 exact head")
    formal = report.get("formal_project_gates", {})
    if formal.get("normalized_momentum_max_l2") != 1.0e-3:
        raise ValueError("formal momentum gate changed")
    if formal.get("divergence_max_l2") != 1.0e-5:
        raise ValueError("formal divergence gate changed")
    if formal.get("formal_full_domain_pde_gate_assessed") is not False:
        raise ValueError("local support audit cannot be relabeled as the formal PDE gate")
    truth = report.get("truth_boundary", {})
    if truth.get("pde_validated") is not False:
        raise ValueError("local support audit cannot promote pde_validated")

    guards = report.get("local_guards", {})
    thresholds = guards.get("thresholds", {})
    if thresholds.get("finest_curl_relative_rms_max") != 2.0e-4:
        raise ValueError("Agent-4 local curl guard changed")
    if thresholds.get("finest_divergence_max_abs_max") != 2.0e-4:
        raise ValueError("Agent-4 local divergence guard changed")
    checks = guards.get("checks", {})
    # The exact audited head is deliberately a scientific local rejection.
    if guards.get("structural_preflight_passed") is not False:
        raise ValueError("expected Agent-4 local structural rejection disappeared")
    if checks.get("finest_cartesian_divergence_max_abs") is not False:
        raise ValueError("expected small-radius divergence guard failure disappeared")
    for key in (
        "finest_cartesian_curl_relative_rms",
        "curl_refinement_ratio",
        "divergence_refinement_ratio",
        "missing_support_gradient_curl_detected",
        "missing_support_gradient_divergence_detected",
    ):
        if checks.get(key) is not True:
            raise ValueError(f"unexpected additional Agent-4 local guard failure: {key}")


def _finest(report: dict[str, Any]) -> dict[str, Any]:
    ladder = report.get("resolution_ladder")
    if not isinstance(ladder, list) or len(ladder) != 3:
        raise ValueError("unexpected Agent-4 resolution ladder")
    row = ladder[-1]
    return {
        "step": float(row["step"]),
        "sample_count": int(row["sample_count"]),
        "curl_relative_rms": float(row["curl_relative_rms"]),
        "divergence_max_abs": float(row["divergence_max_abs"]),
        "divergence_rms_abs": float(row["divergence_rms_abs"]),
        "regions": row["regions"],
        "parameter_cases": row["parameter_cases"],
    }


def build_checkpoint(report: dict[str, Any] | None = None) -> dict[str, Any]:
    if report is None:
        report = build_agent4_report()
    _validate_agent4_report(report)
    finest = _finest(report)
    refinement = report["refinement"]
    mutation = report["mutation"]

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "base": {"agent4_pr": BASE_AGENT4_PR, "agent4_head": BASE_AGENT4_HEAD},
        "fixed_gates": FIXED_GATES,
        "st006_reference": ST006_REFERENCE,
        "upstream": UPSTREAM,
        "agent4_independent_support_localization_audit": {
            "seed": int(report["seed"]),
            "operator": report["operator"],
            "sampling": report["sampling"],
            "finest": finest,
            "refinement": refinement,
            "mutation": mutation,
            "local_guards": report["local_guards"],
            "structural_preflight_passed": False,
            "formal_full_domain_verdict": None,
            "directly_comparable_to_ST006": False,
        },
        "routing": {
            "accepted_fact": (
                "Agent-2 vector-potential localization retains the support-gradient curl terms; "
                "Agent-4 Cartesian curl agreement and mutation sensitivity pass their frozen guards."
            ),
            "blocking_fact": (
                "Agent-4 finest Cartesian divergence max fails the frozen 2e-4 local guard, "
                "with the failure concentrated in the small-radius probe set."
            ),
            "do_not_do": [
                "do not weaken the 2e-4 local divergence guard",
                "do not promote this manufactured local seam to actual-source xyz,t velocity",
                "do not rerun Agent-3 finite correction before a genuine second public column exists",
                "do not treat this local audit as ST006-comparable or as the formal PDE gate",
            ],
            "shortest_next_closure": [
                "Agent 2/4: isolate the small-radius divergence discrepancy with an axis-regular or smaller-step independent check while keeping the same frozen guard",
                "Agent 2: after that seam passes, bind concrete source partition values/derivatives and the actual positive-order background into the existing localized complete-curl chain",
                "Agent 3: screen the resulting genuinely independent public covariance column through the frozen 3x3/radial guards before any two-column finite cycle",
                "Agent 1: independently continue closing segmented-leading radial gaps and matched pressure/forcing; no gap may be invented by Agent 5",
                "Agent 4: run the fixed held-out full-domain 1e-3/1e-5 gate only after a global leading+oscillatory+corrected candidate exists",
            ],
        },
        "baseline_vs_kokuno": [
            {
                "name": "ST006",
                "scope": "held-out full-domain normalized NS momentum",
                "momentum_sampled_max": ST006_REFERENCE["momentum_sampled_max"],
                "volume_l2": ST006_REFERENCE["volume_l2"],
                "pde_validated": False,
            },
            {
                "name": "Kokuno support-localized complete-curl seam",
                "scope": "manufactured local Cartesian curl/divergence implementation audit",
                "momentum_sampled_max": None,
                "volume_l2": None,
                "directly_comparable_to_ST006": False,
                "local_structural_preflight_passed": False,
            },
        ],
        "states": STATES,
        "truth_boundary": {
            "threshold_relaxed": False,
            "surrogate_promoted_to_candidate": False,
            "free_residual_defined_forcing_used": False,
            "kokuno_replay_used_as_independent_pde_validation": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
            "formal_full_domain_pde_gate_assessed": False,
            "pde_validated": False,
        },
    }
    payload["checkpoint_sha256"] = _sha(payload)
    return payload


def validate_checkpoint(payload: dict[str, Any]) -> dict[str, Any]:
    if payload.get("schema") != SCHEMA or payload.get("task_id") != TASK_ID:
        raise ValueError("unexpected Agent-5 checkpoint schema/task")
    signature = payload.get("checkpoint_sha256")
    unsigned = dict(payload)
    unsigned.pop("checkpoint_sha256", None)
    if signature != _sha(unsigned):
        raise ValueError("checkpoint SHA mismatch")
    if payload.get("fixed_gates") != FIXED_GATES:
        raise ValueError("fixed gates changed")
    if payload.get("st006_reference") != ST006_REFERENCE:
        raise ValueError("ST006 reference changed")
    if payload.get("states") != STATES:
        raise ValueError("routing states changed")
    audit = payload.get("agent4_independent_support_localization_audit", {})
    if audit.get("structural_preflight_passed") is not False:
        raise ValueError("failed Agent-4 structural preflight was promoted")
    if audit.get("directly_comparable_to_ST006") is not False:
        raise ValueError("local audit cannot be relabeled as ST006-comparable")
    truth = payload.get("truth_boundary", {})
    for key in (
        "threshold_relaxed",
        "surrogate_promoted_to_candidate",
        "free_residual_defined_forcing_used",
        "kokuno_replay_used_as_independent_pde_validation",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
        "formal_full_domain_pde_gate_assessed",
        "pde_validated",
    ):
        if truth.get(key) is not False:
            raise ValueError(f"truth boundary promoted: {key}")
    return payload


def write_bundle(output_dir: str | Path) -> dict[str, Path]:
    root = Path(output_dir)
    if root.exists():
        raise FileExistsError(f"refusing to overwrite existing output directory: {root}")
    root.mkdir(parents=True)
    report = build_agent4_report()
    checkpoint = build_checkpoint(report)
    report_path = root / "agent4_support_localized_cartesian_audit.json"
    checkpoint_path = root / "support_localization_routing_checkpoint.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    checkpoint_path.write_text(json.dumps(checkpoint, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    loaded = json.loads(checkpoint_path.read_text(encoding="utf-8"))
    validate_checkpoint(loaded)
    return {"agent4_report": report_path, "checkpoint": checkpoint_path}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        default="artifacts/kokuno_agent5/support_localization_routing_checkpoint_v17",
    )
    args = parser.parse_args(argv)
    paths = write_bundle(args.output_dir)
    print(paths["checkpoint"].read_text(encoding="utf-8"), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
