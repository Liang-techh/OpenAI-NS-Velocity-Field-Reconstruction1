"""Agent-5 integration checkpoint for Agent-4's independent signed-curl audit.

This module adds no new Kokuno mathematics.  It keeps the v6 closure-frontier
routing unchanged, replays Agent 4's independent black-box FD2 audit of the
materialized Agent-3 signed complete-curl correction, and binds that evidence
into one fail-closed receipt.  The signed correction remains rejected, the
formal full-domain PDE gate remains unassessed, and the registered thresholds
remain exactly 1e-3 for normalized momentum and 1e-5 for divergence.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_closure_frontier_checkpoint import (
    REGISTERED_DIVERGENCE_THRESHOLD,
    REGISTERED_PDE_THRESHOLD,
    STATES as PARENT_STATES,
    build_checkpoint as build_parent_checkpoint,
    validate_checkpoint as validate_parent_checkpoint,
    write_bundle as write_parent_bundle,
)
from .kokuno_independent_signed_curl_cycle import (
    ANNULUS_COUNT,
    AXIS_COUNT,
    RESIDUAL_REFERENCE,
    DIVERGENCE_REFERENCE,
    STEPS,
    TASK as AGENT4_TASK,
    generate_report as generate_agent4_report,
)


SCHEMA = "kokuno-agent5-independent-routing-checkpoint-v7"
TASK_ID = "KOKUNO-A5-INDEPENDENT-ROUTING-CHECKPOINT-007"
PARENT_AGENT5_HEAD = "852f2d38d5e106997bc88588f58df3854852a0ef"

AGENT4_UPSTREAM = {
    "pr": 264,
    "head_sha": "ea18ec24d3f9251eb60dc4722e98ecef44d53453",
    "dedicated_run": 35284021382,
    "standard_run": 35284021330,
    "artifact_id": 10523007283,
    "artifact_digest": (
        "sha256:6af52a86d8907780ce610bbfde562bbd27285e2a77260f7302925ee8a1c8ba9c"
    ),
    "focused_tests_passed": 7,
    "constrained_tests_passed": 204,
    "reported_finest_step": 0.001,
    "reported_oscillatory_over_leading_rms": 0.98601051785,
    "reported_corrected_over_oscillatory_rms": 1.00008706040,
    "reported_corrected_over_leading_rms": 0.98609636032,
    "reported_corrected_sample_rms": 10.3075981486,
    "reported_corrected_sample_max": 22.7787616638,
    "reported_correction_annulus_nontrivial_rms": 0.0042096510,
    "formal_full_domain_gate_assessed": False,
    "pde_validated": False,
}

STATES = dict(PARENT_STATES)
STATES.update(
    {
        "independent_signed_curl_validation_ready": True,
        "signed_curl_correction_independently_promoted": False,
        "leading_ready": False,
        "correction_ready": False,
        "velocity_export_ready": False,
        "formal_full_domain_pde_gate_assessed": False,
        "pde_validated": False,
    }
)


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _sha(payload: dict[str, Any]) -> str:
    return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()


def summarize_agent4_report(report: dict[str, Any]) -> dict[str, Any]:
    """Extract only routing-relevant evidence from Agent 4's independent report."""
    if not isinstance(report, dict) or report.get("task") != AGENT4_TASK:
        raise ValueError("unexpected Agent-4 report/task")
    validator = report.get("validator", {})
    if validator.get("agent3_cr006_residual_operator_used") is not False:
        raise ValueError("Agent-4 report is not independent of Agent-3 residual operator")
    if validator.get("training_loss_or_tensor_used") is not False:
        raise ValueError("Agent-4 report reused training data/loss")
    if float(validator.get("fixed_residual_reference", np.nan)) != REGISTERED_PDE_THRESHOLD:
        raise ValueError("Agent-4 residual reference differs from registered gate")
    if float(validator.get("fixed_divergence_reference", np.nan)) != REGISTERED_DIVERGENCE_THRESHOLD:
        raise ValueError("Agent-4 divergence reference differs from registered gate")
    if RESIDUAL_REFERENCE != REGISTERED_PDE_THRESHOLD or DIVERGENCE_REFERENCE != REGISTERED_DIVERGENCE_THRESHOLD:
        raise RuntimeError("imported Agent-4 constants drifted from registered gates")

    finest_key = str(float(STEPS[-1]))
    finest = report.get("aggregates", {}).get(finest_key, {})
    required_stages = {
        "leading_only",
        "leading_plus_oscillatory",
        "after_signed_correction",
    }
    if set(finest) != required_stages:
        raise ValueError("Agent-4 finest aggregate stages are incomplete")

    ratios = report.get("finest_stage_ratios", {})
    support = report.get("signed_correction_support", {})
    truth = report.get("truth_boundary", {})
    local = report.get("local_reference_checks", {})
    sampling = report.get("sampling", {})

    values = {
        "finest_step": float(STEPS[-1]),
        "annulus_count": int(sampling.get("annulus_count", -1)),
        "axis_count": int(sampling.get("axis_count", -1)),
        "leading_sample_rms": float(finest["leading_only"]["sample_rms"]),
        "oscillatory_sample_rms": float(finest["leading_plus_oscillatory"]["sample_rms"]),
        "corrected_sample_rms": float(finest["after_signed_correction"]["sample_rms"]),
        "corrected_sample_max": float(finest["after_signed_correction"]["sample_max"]),
        "corrected_divergence_rms": float(finest["after_signed_correction"]["divergence_rms"]),
        "corrected_divergence_max": float(finest["after_signed_correction"]["divergence_max"]),
        "oscillatory_over_leading_rms": float(ratios["oscillatory_over_leading_rms"]),
        "corrected_over_oscillatory_rms": float(ratios["corrected_over_oscillatory_rms"]),
        "corrected_over_leading_rms": float(ratios["corrected_over_leading_rms"]),
        "correction_outside_probe_exact_zero": bool(support["outside_probe_exact_zero"]),
        "correction_annulus_nontrivial_rms": float(support["annulus_nontrivial_rms"]),
        "local_residual_rms_reference_pass": bool(
            local["finest_corrected_normalized_sample_rms_le_1e_minus_3"]
        ),
        "local_residual_max_reference_pass": bool(
            local["finest_corrected_normalized_sample_max_le_1e_minus_3"]
        ),
        "local_divergence_rms_reference_pass": bool(
            local["finest_corrected_divergence_rms_le_1e_minus_5"]
        ),
        "local_divergence_max_reference_pass": bool(
            local["finest_corrected_divergence_max_le_1e_minus_5"]
        ),
        "formal_full_domain_gate_assessed": bool(truth["formal_full_domain_gate_assessed"]),
        "pde_validated": bool(truth["pde_validated"]),
    }
    if not all(np.isfinite(v) for k, v in values.items() if isinstance(v, float)):
        raise ValueError("Agent-4 routing evidence contains nonfinite values")
    if values["annulus_count"] != ANNULUS_COUNT or values["axis_count"] != AXIS_COUNT:
        raise ValueError("Agent-4 sampling contract changed")
    if values["oscillatory_over_leading_rms"] >= 1.0:
        raise ValueError("expected independently observed oscillatory improvement is absent")
    if values["corrected_over_oscillatory_rms"] <= 1.0:
        raise ValueError("signed correction no longer has the independently observed adverse routing signal")
    if values["corrected_sample_rms"] <= REGISTERED_PDE_THRESHOLD:
        raise ValueError("unexpected local residual crossing requires a fresh scientific review")
    if values["correction_annulus_nontrivial_rms"] <= 0.0:
        raise ValueError("signed correction became trivial")
    if not values["correction_outside_probe_exact_zero"]:
        raise ValueError("signed correction support contract changed")
    if values["formal_full_domain_gate_assessed"] or values["pde_validated"]:
        raise ValueError("Agent-4 local replay cannot be promoted to formal PDE validation")
    return values


def _check_replay_matches_agent4_reference(summary: dict[str, Any]) -> None:
    checks = {
        "finest_step": AGENT4_UPSTREAM["reported_finest_step"],
        "oscillatory_over_leading_rms": AGENT4_UPSTREAM["reported_oscillatory_over_leading_rms"],
        "corrected_over_oscillatory_rms": AGENT4_UPSTREAM["reported_corrected_over_oscillatory_rms"],
        "corrected_over_leading_rms": AGENT4_UPSTREAM["reported_corrected_over_leading_rms"],
        "corrected_sample_rms": AGENT4_UPSTREAM["reported_corrected_sample_rms"],
        "corrected_sample_max": AGENT4_UPSTREAM["reported_corrected_sample_max"],
        "correction_annulus_nontrivial_rms": AGENT4_UPSTREAM[
            "reported_correction_annulus_nontrivial_rms"
        ],
    }
    for key, expected in checks.items():
        actual = float(summary[key])
        if not np.isclose(actual, float(expected), rtol=2.0e-8, atol=2.0e-10):
            raise RuntimeError(
                f"independent Agent-4 replay drift for {key}: {actual} != {expected}"
            )


def build_checkpoint(agent4_report: dict[str, Any]) -> dict[str, Any]:
    parent = validate_parent_checkpoint(build_parent_checkpoint())
    summary = summarize_agent4_report(agent4_report)
    _check_replay_matches_agent4_reference(summary)

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "parent_agent5": {
            "head_sha": PARENT_AGENT5_HEAD,
            "checkpoint_sha256": parent["checkpoint_sha256"],
            "routing_preserved": True,
        },
        "agent4_independent_signed_curl_audit": {
            "provenance": AGENT4_UPSTREAM,
            "replayed_under_agent5": True,
            "replay_matches_exact_head_reference": True,
            "evidence": summary,
            "used_for_parameter_selection": False,
            "used_as_formal_full_domain_pde_gate": False,
        },
        "fixed_gates": {
            "held_out_normalized_full_momentum_residual": REGISTERED_PDE_THRESHOLD,
            "divergence_max": REGISTERED_DIVERGENCE_THRESHOLD,
            "changed": False,
        },
        "routing": {
            "apply_agent3_signed_curl_correction": False,
            "agent3_rejection_independently_supported": True,
            "reason": (
                "Agent 4 independently finds the signed correction slightly worsens the "
                "oscillatory-stage raw RMS while Agent 3 finds the intended mean-defect "
                "target worsens materially"
            ),
            "priority_blocker": "axis-near axial leading momentum and global leading matching",
            "retune_single_signed_column_now": False,
            "run_formal_full_domain_gate": False,
            "formal_gate_blocker": (
                "global matched core-to-heat leading field and compatible final forcing "
                "contract are still absent"
            ),
        },
        "states": STATES,
        "truth_boundary": {
            "source_version_and_provenance_preserved": True,
            "free_residual_defined_forcing_used": False,
            "pressure_or_forcing_refit_in_agent5": False,
            "agent4_replay_called_formal_full_domain_validation": False,
            "surrogate_promoted_to_final_candidate": False,
            "threshold_relaxed": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
        "next_shortest_closure": [
            "Agent 1: reduce the dominant axis-near axial/global leading mismatch and complete matched core-to-heat assembly",
            "Agent 2: compose the source phase/covector only after the Agent-1/base V/G mapping exists",
            "Agent 3: do not widen the rejected one-column signed correction; revisit a source-motivated independent covariance/stress column only after the leading blocker is reduced",
            "Agent 4: run the unchanged full-domain 1e-3 gate only after the global composite and compatible forcing contract exist",
        ],
    }
    payload["checkpoint_sha256"] = _sha(payload)
    return payload


def validate_checkpoint(payload: dict[str, Any]) -> dict[str, Any]:
    if not isinstance(payload, dict):
        raise ValueError("checkpoint must be a JSON object")
    unsigned = dict(payload)
    claimed = unsigned.pop("checkpoint_sha256", None)
    if claimed != _sha(unsigned):
        raise ValueError("checkpoint SHA mismatch")
    if unsigned.get("schema") != SCHEMA or unsigned.get("task_id") != TASK_ID:
        raise ValueError("checkpoint schema/task mismatch")
    if unsigned.get("fixed_gates") != {
        "held_out_normalized_full_momentum_residual": REGISTERED_PDE_THRESHOLD,
        "divergence_max": REGISTERED_DIVERGENCE_THRESHOLD,
        "changed": False,
    }:
        raise ValueError("registered gates changed")
    if unsigned.get("states") != STATES:
        raise ValueError("scientific states changed")
    parent = unsigned.get("parent_agent5", {})
    if parent.get("head_sha") != PARENT_AGENT5_HEAD or parent.get("routing_preserved") is not True:
        raise ValueError("parent Agent-5 identity/routing changed")
    audit = unsigned.get("agent4_independent_signed_curl_audit", {})
    if audit.get("provenance") != AGENT4_UPSTREAM:
        raise ValueError("Agent-4 provenance changed")
    if audit.get("used_for_parameter_selection") is not False:
        raise ValueError("independent validation cannot be reused for selection")
    if audit.get("used_as_formal_full_domain_pde_gate") is not False:
        raise ValueError("local independent audit cannot be relabeled as full-domain gate")
    routing = unsigned.get("routing", {})
    if routing.get("apply_agent3_signed_curl_correction") is not False:
        raise ValueError("rejected Agent-3 signed correction cannot be routed")
    if routing.get("agent3_rejection_independently_supported") is not True:
        raise ValueError("independent routing conclusion changed")
    if routing.get("retune_single_signed_column_now") is not False:
        raise ValueError("validation result cannot silently authorize retuning")
    if routing.get("run_formal_full_domain_gate") is not False:
        raise ValueError("formal full-domain gate remains blocked")
    truth = unsigned.get("truth_boundary", {})
    required_truth = {
        "source_version_and_provenance_preserved": True,
        "free_residual_defined_forcing_used": False,
        "pressure_or_forcing_refit_in_agent5": False,
        "agent4_replay_called_formal_full_domain_validation": False,
        "surrogate_promoted_to_final_candidate": False,
        "threshold_relaxed": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
    }
    if truth != required_truth:
        raise ValueError("truth boundary changed")
    if STATES["pde_validated"] is not False:
        raise ValueError("PDE state cannot be promoted")
    return payload


def load_checkpoint(path: str | Path) -> dict[str, Any]:
    return validate_checkpoint(json.loads(Path(path).read_text(encoding="utf-8")))


def write_bundle(output_dir: str | Path) -> dict[str, Any]:
    target = Path(output_dir)
    if target.exists() and any(target.iterdir()):
        raise FileExistsError("output directory must be absent or empty")
    target.mkdir(parents=True, exist_ok=True)

    parent_dir = target / "closure_frontier_v6"
    parent = write_parent_bundle(parent_dir)
    if parent["checkpoint_sha256"] != build_parent_checkpoint()["checkpoint_sha256"]:
        raise RuntimeError("parent Agent-5 checkpoint drifted during bundle generation")

    report_path = target / "independent_signed_curl_report.json"
    agent4_report = generate_agent4_report(
        output=report_path,
        annulus_count=ANNULUS_COUNT,
        axis_count=AXIS_COUNT,
    )
    checkpoint = build_checkpoint(agent4_report)
    checkpoint_path = target / "independent_routing_checkpoint.json"
    checkpoint_path.write_text(
        json.dumps(checkpoint, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    load_checkpoint(checkpoint_path)
    return checkpoint


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        default="artifacts/kokuno_agent5/independent_routing_checkpoint_v7",
    )
    args = parser.parse_args(argv)
    checkpoint = write_bundle(args.output_dir)
    print(json.dumps(checkpoint, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
