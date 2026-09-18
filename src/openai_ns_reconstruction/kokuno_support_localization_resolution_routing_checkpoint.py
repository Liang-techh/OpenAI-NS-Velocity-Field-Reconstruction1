"""Kokuno Agent-5 routing after the support-localization axis refinement.

This increment adds no new Kokuno mathematics. It binds Agent 4 PR #360's
independent FD4 resolution of the earlier manufactured small-radius divergence
rejection, while keeping Agent 2 PR #358's source-wave-shell wrapper and Agent 3
PR #357's radial-force completeness result as explicit sibling dependencies.

The manufactured local pass is not a source-wave-shell pass and is not the
formal Navier--Stokes gate. Formal held-out normalized momentum max/L2 remain
<=1e-3 and divergence max/L2 <=1e-5. No free residual-defined forcing is used.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

from .kokuno_support_localized_axis_refinement import build_report as build_agent4_report

SCHEMA = "kokuno-agent5-support-localization-resolution-routing-checkpoint-v18"
TASK_ID = "KOKUNO-A5-SUPPORT-LOCALIZATION-RESOLUTION-ROUTING-018"
BASE_AGENT4_PR = 360
BASE_AGENT4_HEAD = "948200cfac9c95e5bed2cf08371c5d113f00188f"

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
        "pr": 358,
        "head_sha": "e24d918bb0cc16b1f1006710183058296fbf760f",
        "dedicated_run": 35325772034,
        "standard_run": 35325771949,
        "source_wave_shell_wrapper_ready": True,
        "source_angular_grid_enforced": True,
        "source_inner_annulus_guard_ready": True,
        "source_actual_lambda_recovered": False,
        "concrete_source_partition_bumps_reconstructed": False,
        "source_actual_background_path_instantiated": False,
        "public_source_oscillatory_xyz_t_velocity_ready": False,
        "genuinely_independent_second_covariance_column_ready": False,
        "consumed_in_executable_ancestry": False,
    },
    "agent3": {
        "pr": 357,
        "head_sha": "169558a6994a255119d6af97430e717a221c034c",
        "dedicated_run": 35325747834,
        "standard_run": 35325747777,
        "artifact_id": 10539133848,
        "artifact_digest": "sha256:a3022ccca9a66d07515ea67e8f14fc18064804cf8eb6da3aba7f783c61e874b6",
        "radial_force_rms": 8.718341895613187e-05,
        "radial_force_max": 1.7292926497973056e-04,
        "radial_to_tangential_force_rms_ratio": 0.4638871624620116,
        "fine_pair_radial_force_relative_change": 0.0010797618861827264,
        "radial_force_derivative_stability_preflight_passed": True,
        "finite_correction_cycle_rerun_allowed": False,
        "consumed_in_executable_ancestry": False,
    },
    "agent4": {
        "pr": BASE_AGENT4_PR,
        "head_sha": BASE_AGENT4_HEAD,
        "dedicated_run": 35326140220,
        "standard_run": 35326140228,
        "artifact_id": 10538894527,
        "artifact_digest": "sha256:6cb6317bcc58ba4c65ec48a9e7eb614e5552f1ff34929a57681e31bd9821189d",
        "manufactured_support_localization_axis_refinement_ready": True,
        "consumed_in_executable_ancestry": True,
    },
    "previous_agent5": {
        "pr": 356,
        "head_sha": "f8617b6bc14ab4d3053ac16eff89f4a279bcb163",
        "recorded_parent_fd2_rejection": True,
        "consumed_in_executable_ancestry": False,
    },
}

STATES = {
    "segmented_leading_velocity_router_ready": True,
    "support_localized_vector_potential_curl_ready": True,
    "manufactured_support_localization_axis_refinement_passed": True,
    "prior_fd2_small_radius_rejection_resolved_as_persistent_leak": False,
    "source_wave_shell_wrapper_ready": True,
    "source_wave_shell_independent_cartesian_preflight_passed": False,
    "leading_ready": False,
    "oscillatory_ready": True,
    "public_source_oscillatory_xyz_t_velocity_ready": False,
    "genuinely_independent_second_covariance_column_ready": False,
    "radial_force_side_effect_retained": True,
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
    if report.get("task_id") != "KOKUNO-A4-SUPPORT-LOCALIZED-AXIS-REFINEMENT-018":
        raise ValueError("unexpected Agent-4 axis-refinement task")
    dependency = report.get("dependency", {})
    if dependency.get("agent4_parent_pr") != 353:
        raise ValueError("Agent-4 refinement is not bound to expected parent PR")
    formal = report.get("formal_project_gates", {})
    if formal.get("normalized_momentum_max_l2") != 1.0e-3:
        raise ValueError("formal momentum gate changed")
    if formal.get("divergence_max_l2") != 1.0e-5:
        raise ValueError("formal divergence gate changed")
    if formal.get("formal_full_domain_pde_gate_assessed") is not False:
        raise ValueError("local refinement cannot be relabeled as formal PDE gate")
    truth = report.get("truth_boundary", {})
    if truth.get("pde_validated") is not False:
        raise ValueError("local refinement cannot promote pde_validated")

    guards = report.get("local_guards", {})
    thresholds = guards.get("thresholds", {})
    if thresholds.get("finest_curl_relative_rms_max") != 2.0e-4:
        raise ValueError("Agent-4 local curl guard changed")
    if thresholds.get("finest_divergence_max_abs_max") != 2.0e-4:
        raise ValueError("Agent-4 local divergence guard changed")
    if guards.get("same_thresholds_as_parent_rejection") is not True:
        raise ValueError("parent rejection guard was not preserved")
    if guards.get("structural_preflight_passed") is not True:
        raise ValueError("expected Agent-4 manufactured local pass is absent")
    checks = guards.get("checks", {})
    if not checks or not all(value is True for value in checks.values()):
        raise ValueError("one or more frozen Agent-4 local checks failed")


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
        "agent4_axis_refinement": {
            "seed": int(report["seed"]),
            "operator": report["operator"],
            "sampling": report["sampling"],
            "finest": finest,
            "refinement": refinement,
            "mutation": mutation,
            "local_guards": report["local_guards"],
            "manufactured_structural_preflight_passed": True,
            "source_wave_shell_preflight_assessed": False,
            "formal_full_domain_verdict": None,
            "directly_comparable_to_ST006": False,
        },
        "routing": {
            "resolved_fact": (
                "The prior FD2 small-radius divergence rejection does not persist under Agent-4's fresh Cartesian FD4/smaller-step audit with the same 2e-4 guard; the manufactured complete-curl/product-rule seam passes."
            ),
            "remaining_oscillatory_blocker": (
                "Agent-2 PR #358 adds the source angular grid and away-from-axis source wave shell, but that exact source-domain wrapper has not yet been independently audited by Agent 4."
            ),
            "correction_completeness_fact": (
                "Agent-3 PR #357 shows the source-required radial d_z sigma_1 force is stable and about 46.39% of reconstructed tangential-force RMS, so a future two-column theta/z fit must retain the radial side effect."
            ),
            "do_not_do": [
                "do not treat the manufactured FD4 pass as validation of the source-wave-shell wrapper",
                "do not weaken the 2e-4 local support-localization guards",
                "do not drop Agent-3 radial d_z sigma_1 when materializing a correction",
                "do not rerun the finite correction cycle before a genuinely independent public second column passes the frozen screens",
                "do not treat any local support/covariance diagnostic as ST006-comparable or as the formal PDE gate",
            ],
            "shortest_next_closure": [
                "Agent 4: consume Agent-2 PR #358 exact head and rerun the unchanged Cartesian curl/divergence guard on a source-grid phase inside the source wave shell",
                "Agent 2: only after that independent source-domain pass, bind concrete source partition values/derivatives plus the actual positive-order background and expose direct Q-scaled velocity(x,y,z,t)",
                "Agent 3: screen that exact public column through the frozen spacetime/radial bounded-inverse guards and retain radial d_z sigma_1 before any finite cycle",
                "Agent 1: close the segmented-leading radial gaps and globally matched pressure/admissible forcing independently of this oscillatory seam",
                "Agent 4: run the fixed held-out full-domain 1e-3/1e-5 gate only after the complete global leading+oscillatory+corrected candidate exists",
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
                "name": "Kokuno support-localized complete-curl manufactured seam",
                "scope": "local Cartesian FD4 curl/divergence implementation audit",
                "momentum_sampled_max": None,
                "volume_l2": None,
                "directly_comparable_to_ST006": False,
                "local_structural_preflight_passed": True,
                "source_wave_shell_preflight_assessed": False,
            },
        ],
        "states": STATES,
        "truth_boundary": {
            "threshold_relaxed": False,
            "surrogate_promoted_to_candidate": False,
            "free_residual_defined_forcing_used": False,
            "kokuno_replay_used_as_independent_pde_validation": False,
            "manufactured_local_pass_promoted_to_source_domain": False,
            "radial_force_side_effect_dropped": False,
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
    audit = payload.get("agent4_axis_refinement", {})
    if audit.get("manufactured_structural_preflight_passed") is not True:
        raise ValueError("Agent-4 manufactured local pass disappeared")
    if audit.get("source_wave_shell_preflight_assessed") is not False:
        raise ValueError("source-wave-shell audit was falsely promoted")
    if audit.get("directly_comparable_to_ST006") is not False:
        raise ValueError("local audit cannot be relabeled as ST006-comparable")
    truth = payload.get("truth_boundary", {})
    for key in (
        "threshold_relaxed",
        "surrogate_promoted_to_candidate",
        "free_residual_defined_forcing_used",
        "kokuno_replay_used_as_independent_pde_validation",
        "manufactured_local_pass_promoted_to_source_domain",
        "radial_force_side_effect_dropped",
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
    report_path = root / "agent4_support_localized_axis_refinement.json"
    checkpoint_path = root / "support_localization_resolution_routing_checkpoint.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    checkpoint_path.write_text(json.dumps(checkpoint, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    validate_checkpoint(json.loads(checkpoint_path.read_text(encoding="utf-8")))
    return {"agent4_report": report_path, "checkpoint": checkpoint_path}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        default="artifacts/kokuno_agent5/support_localization_resolution_routing_checkpoint_v18",
    )
    args = parser.parse_args(argv)
    paths = write_bundle(args.output_dir)
    print(json.dumps({key: str(value) for key, value in paths.items()}, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
