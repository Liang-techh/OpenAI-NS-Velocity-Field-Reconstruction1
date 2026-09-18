"""Kokuno Agent-5 segmented-frontier integration/routing checkpoint.

This module adds no new Kokuno mathematics.  It harvests Agent 1's first
candidate-facing segmented leading velocity router and Agent 4's independent
black-box audit, while binding the newest Agent 2/3 lane state as sibling
provenance only.

The checkpoint deliberately distinguishes "a unified public router exists for
the reconstructed segments" from "a global leading candidate exists".  Large
radial gaps and the missing matched I2/global pressure remain machine-visible,
so ``leading_ready`` and every formal PDE state stay false.

The fixed project gates are unchanged: held-out normalized momentum <= 1e-3
and divergence <= 1e-5.  No free residual-defined forcing is introduced.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

from .kokuno_segmented_leading_assembly import KokunoSegmentedLeadingAssembly
from .kokuno_segmented_leading_independent import run_audit as run_agent4_audit


SCHEMA = "kokuno-agent5-segmented-frontier-routing-checkpoint-v16"
TASK_ID = "KOKUNO-A5-SEGMENTED-FRONTIER-ROUTING-CHECKPOINT-016"
BASE_AGENT4_PR = 348
BASE_AGENT4_HEAD = "c66ab597278bb6a35274c99a2f9db7035a0a7926"

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
        "standard_run": 35319669834,
        "segmented_candidate_velocity_router_ready": True,
        "reference_velocity_pressure_stage_ready": True,
        "audited_repaired_i2_velocity_stage_ready": True,
        "unsupported_radial_gaps_fail_closed": True,
        "global_matched_pressure_ready": False,
        "global_leading_profile_reconstructed": False,
        "consumed_in_executable_ancestry": True,
    },
    "agent2": {
        "pr": 347,
        "head_sha": "4675a0a3b04fe8891357619a7eb5398500539f64",
        "source_phase_label_run": 35320962052,
        "native_phase_contract_run": 35320962056,
        "standard_run": 35320962079,
        "source_compatible_phase_label_selector_ready": True,
        "homogeneous_pulse_spatial_sensitivity_ready": True,
        "actual_positive_order_source_background_path_ready": False,
        "public_source_oscillatory_xyz_t_velocity_ready": False,
        "genuinely_independent_second_covariance_column_ready": False,
        "consumed_in_executable_ancestry": False,
    },
    "agent3": {
        "pr": 346,
        "head_sha": "61c0f91dc4b0175cd398a7b27f941ba8b2dd1cd4",
        "dedicated_run": 35320473725,
        "standard_run": 35320473638,
        "agent5_reference_run": 35320473641,
        "artifact_id": 10537575118,
        "artifact_digest": "sha256:290268203086d941a99798a418af2461b373caf351b834f0477dd07ec22d74fd",
        "radial_counts": [33, 65, 129],
        "fine_pair_stability_tolerance": 0.02,
        "fine_pair_max_guarded_relative_change": 0.00029353357568215327,
        "finest_missing_relative_vector_rms": 0.2101949905750188,
        "finest_nodes_requiring_second_direction": 99,
        "duplicate_control_rank2_required_nodes": 0,
        "duplicate_control_required_nodes": 99,
        "radial_target_stable_for_future_column_screen": True,
        "second_public_covariance_column_available": False,
        "finite_correction_cycle_rerun_allowed": False,
        "consumed_in_executable_ancestry": False,
    },
    "agent4": {
        "pr": BASE_AGENT4_PR,
        "head_sha": BASE_AGENT4_HEAD,
        "standard_run": 35322435268,
        "independent_segmented_public_contract_audit_executable": True,
        "conservative_i2_attribution_boundary_pinned": True,
        "consumed_in_executable_ancestry": True,
    },
    "prior_agent4_heat_audit": {
        "pr": 339,
        "head_sha": "abf4152b981fba12a493f1fc9db0495661995505",
        "dedicated_run": 35316852392,
        "standard_run": 35316876623,
        "continuous_high_precision_three_moment_i2_repair_independently_preflighted": True,
        "consumed_in_executable_ancestry": True,
    },
    "previous_agent5": {
        "pr": 340,
        "head_sha": "49ab9879cb86142e50760d96e8874103ec8e8d17",
        "high_precision_heat_handoff_frontier_recorded": True,
        "consumed_in_executable_ancestry": True,
    },
}

STATES = {
    "segmented_leading_velocity_router_ready": True,
    "reference_velocity_pressure_stage_ready": True,
    "repaired_i2_velocity_stage_ready": True,
    "global_radial_coverage_complete": False,
    "global_matched_pressure_ready": False,
    "leading_ready": False,
    # "oscillatory_ready" means the source/complete-curl machinery lane is usable;
    # the stronger candidate-facing public source velocity flag remains false.
    "oscillatory_ready": True,
    "public_source_oscillatory_xyz_t_velocity_ready": False,
    "genuinely_independent_second_covariance_column_ready": False,
    "correction_ready": False,
    "complete_kokuno_composite_velocity_ready": False,
    "velocity_export_ready": False,
    "reference_local_zero_force_residual_gate_met": False,
    "reference_local_zero_force_divergence_gate_met": False,
    "formal_full_domain_pde_gate_assessed": False,
    "pde_validated": False,
}


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _sha(payload: dict[str, Any]) -> str:
    return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()


def _finest_region_summary(audit: dict[str, Any], region: str) -> dict[str, Any]:
    stage = audit["reference_stage"][region]
    ladder = stage["resolution_ladder"]
    if not isinstance(ladder, list) or len(ladder) != 3:
        raise ValueError(f"unexpected Agent-4 resolution ladder for {region}")
    finest = ladder[-1]
    return {
        "step": float(finest["step"]),
        "sample_count": int(finest["sample_count"]),
        "max_vector_residual": float(finest["max_vector_residual"]),
        "sample_l2_vector_residual": float(finest["sample_l2_vector_residual"]),
        "divergence_max_abs": float(finest["divergence_max_abs"]),
        "divergence_sample_l2": float(finest["divergence_sample_l2"]),
        "local_residual_gate_met": bool(stage["finest_local_residual_gate_met"]),
        "local_divergence_gate_met": bool(stage["finest_local_divergence_gate_met"]),
        "fine_to_previous_relative_change": stage["fine_to_previous_relative_change"],
    }


def _validate_agent4_audit(audit: dict[str, Any]) -> None:
    if audit.get("task_id") != "KOKUNO-A4-SEGMENTED-PUBLIC-CONTRACT-AUDIT-016":
        raise ValueError("unexpected Agent-4 audit task")
    if audit.get("base", {}).get("head") != UPSTREAM["agent1"]["head_sha"]:
        raise ValueError("Agent-4 audit is not bound to the expected Agent-1 segmented head")
    contract = audit.get("fixed_project_contract", {})
    expected = {
        "nu": 0.01,
        "normalized_momentum_max_gate": 1.0e-3,
        "normalized_momentum_l2_gate": 1.0e-3,
        "divergence_max_gate": 1.0e-5,
        "divergence_l2_gate": 1.0e-5,
        "thresholds_changed": False,
    }
    if contract != expected:
        raise ValueError("Agent-4 fixed project contract changed")
    truth = audit.get("truth_boundary", {})
    if truth.get("formal_full_domain_pde_gate_assessed") is not False:
        raise ValueError("segmented audit cannot be relabeled as a full-domain PDE gate")
    if truth.get("pde_validated") is not False:
        raise ValueError("segmented audit cannot promote pde_validated")
    if audit.get("preregistered_local_guards", {}).get(
        "all_local_implementation_guards_passed"
    ) is not True:
        raise ValueError("Agent-4 local implementation guards did not pass")


def build_checkpoint(audit: dict[str, Any] | None = None) -> dict[str, Any]:
    """Build one fail-closed routing receipt from the current Agent-4 public audit."""

    if audit is None:
        audit = run_agent4_audit()
    _validate_agent4_audit(audit)

    assembly = KokunoSegmentedLeadingAssembly()
    coverage = assembly.coverage_report()
    if coverage.get("global_coverage_complete") is not False:
        raise RuntimeError("segmented assembly unexpectedly claims global radial coverage")
    gaps = coverage.get("unreconstructed_log_X_gaps")
    if not isinstance(gaps, list) or len(gaps) < 1:
        raise RuntimeError("unreconstructed radial gaps disappeared without a new source stage")

    i2 = audit["i2_public_contract"]
    i2_visible = bool(i2["public_float64_heat_repair_visible"])

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "base": {"agent4_pr": BASE_AGENT4_PR, "agent4_head": BASE_AGENT4_HEAD},
        "fixed_gates": FIXED_GATES,
        "st006_reference": ST006_REFERENCE,
        "upstream": UPSTREAM,
        "segmented_leading_contract": {
            "assembly_sha256": assembly.sha256,
            "coverage": coverage,
            "candidate_facing_velocity_router": "KokunoSegmentedLeadingAssembly.velocity",
            "candidate_facing_pressure_router": "KokunoSegmentedLeadingAssembly.pressure",
            "unsupported_regions_fail_closed": True,
            "global_candidate": False,
        },
        "agent4_independent_audit": {
            "seed": int(audit["seed"]),
            "operator": audit["independent_operator"],
            "off_grid_finest": _finest_region_summary(audit, "off_grid"),
            "axis_near_finest": _finest_region_summary(audit, "axis_near"),
            "mutation_calibration": audit["mutation_calibration"],
            "local_implementation_guards": audit["preregistered_local_guards"],
            "i2_public_contract": i2,
            "formal_full_domain_verdict": None,
            "directly_comparable_to_ST006": False,
        },
        "baseline_vs_kokuno": [
            {
                "name": "ST006",
                "scope": "held-out full-domain normalized NS momentum",
                "protocol": "seed 9172801 / 4096 Cartesian points / six times / finest h=.005",
                "momentum_sampled_max": ST006_REFERENCE["momentum_sampled_max"],
                "volume_l2": ST006_REFERENCE["volume_l2"],
                "pde_validated": False,
            },
            {
                "name": "Kokuno segmented leading public contract",
                "scope": "reference-stage local FD2 audit plus repaired-I2 representation seam",
                "momentum_sampled_max": None,
                "volume_l2": None,
                "directly_comparable_to_ST006": False,
                "reason": (
                    "the Kokuno object is still segmented, has unreconstructed radial gaps, "
                    "lacks matched I2/global pressure, and is not a complete global candidate"
                ),
                "pde_validated": False,
            },
        ],
        "states": {
            **STATES,
            "i2_high_precision_heat_repair_prior_independently_preflighted": True,
            "i2_float64_public_repair_observable_at_agent4_probes": i2_visible,
            "i2_float64_public_repair_robustly_attributable": bool(
                audit["truth_boundary"].get(
                    "i2_float64_public_repair_robustly_attributable", False
                )
            ),
            "agent3_radial_target_stability_preflight_passed": True,
            "reference_local_zero_force_residual_gate_met": bool(
                audit["reference_stage"]["off_grid"]["finest_local_residual_gate_met"]
                and audit["reference_stage"]["axis_near"]["finest_local_residual_gate_met"]
            ),
            "reference_local_zero_force_divergence_gate_met": bool(
                audit["reference_stage"]["off_grid"]["finest_local_divergence_gate_met"]
                and audit["reference_stage"]["axis_near"]["finest_local_divergence_gate_met"]
            ),
        },
        "routing": {
            "run_formal_full_domain_gate_now": False,
            "export_segmented_router_as_complete_candidate": False,
            "hide_or_bridge_unreconstructed_radial_gaps_autonomously": False,
            "rerun_one_or_duplicate_column_finite_cycle": False,
            "replace_ST006_as_repository_pde_baseline": False,
            "agent1_next": (
                "first bind the admissible source/reference forcing contract (or prove zero forcing "
                "is the intended local contract) to explain the independent reference-stage raw "
                "zero-force FD2 residual, which is far above the project gate; if it remains large "
                "under the correct allowed forcing, repair the reference stage. In parallel, "
                "reconstruct the missing radial source stages and matched non-reference/global "
                "pressure. The I2 repair changes a subset of float64 probes only at rounding scale, "
                "so keep module-level after-repair attribution fail-closed until a "
                "precision-aware/rescaled public contract exists"
            ),
            "agent2_next": (
                "bind the selected source phase label to an actual positive-order/background "
                "path and materialize a direct public Q-scaled velocity(x,y,z,t); provide a "
                "genuinely independent second covariance column"
            ),
            "agent3_next": (
                "screen that exact Agent-2 second column through the frozen spacetime guard and "
                "the now resolution-stable radial target; only if all guards pass, materialize "
                "a bounded two-column correction and rerun the frozen finite correction cycle"
            ),
            "agent4_next": (
                "after a global leading+oscillatory+correction velocity/pressure/allowed-forcing "
                "candidate exists, run the unchanged held-out <=1e-3 momentum and <=1e-5 "
                "divergence gate plus a same-protocol ST006 comparison"
            ),
            "agent5_next": (
                "keep the segmented schema/provenance as an intermediate artifact; promote to a "
                "single save/load/export candidate schema only when the global pressure/velocity "
                "and the Agent-2/3 composite are actually available"
            ),
        },
        "truth_boundary": {
            "source_version_and_provenance_preserved": True,
            "agent5_new_kokuno_mathematics_added": False,
            "agent2_or_agent3_sibling_code_laundered_into_executable_ancestry": False,
            "segmented_velocity_router_called_global": False,
            "local_reference_fd2_called_formal_full_domain_validation": False,
            "local_zero_force_reference_failure_hidden": False,
            "unreconstructed_radial_gaps_hidden": False,
            "i2_float64_visibility_reported_exactly_as_measured": True,
            "i2_float64_module_level_repair_attribution_promoted": False,
            "free_residual_defined_forcing_used": False,
            "pressure_or_forcing_refit_in_agent5": False,
            "threshold_relaxed": False,
            "st006_used_as_kokuno_source_truth": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
        "next_shortest_closure": [
            "Agent 1/4: bind the admissible source/reference forcing contract and rerun the same independent local operator; if O(10) residual persists, repair the reference stage",
            "Agent 1: missing radial source stages + matched non-reference/global pressure; preserve explicit I2 attribution",
            "Agent 2: actual positive-order/background source path + public xyz,t oscillatory velocity + independent second covariance column",
            "Agent 3: resolution-stable radial/spacetime column guard -> bounded two-column finite cycle only after a real second column passes",
            "Agent 4: formal full-domain held-out 1e-3 / 1e-5 gate only after the global composite exists",
            "Agent 5: final unified candidate artifact, ST006 same-protocol comparison, Python/MATLAB smoke only after that gateable global API exists",
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
    if unsigned.get("base") != {
        "agent4_pr": BASE_AGENT4_PR,
        "agent4_head": BASE_AGENT4_HEAD,
    }:
        raise ValueError("base Agent-4 provenance changed")
    if unsigned.get("fixed_gates") != FIXED_GATES:
        raise ValueError("fixed gates changed")
    if unsigned.get("st006_reference") != ST006_REFERENCE:
        raise ValueError("ST006 reference changed")
    if unsigned.get("upstream") != UPSTREAM:
        raise ValueError("upstream provenance changed")
    states = unsigned.get("states", {})
    for key in (
        "leading_ready",
        "public_source_oscillatory_xyz_t_velocity_ready",
        "genuinely_independent_second_covariance_column_ready",
        "correction_ready",
        "complete_kokuno_composite_velocity_ready",
        "velocity_export_ready",
        "formal_full_domain_pde_gate_assessed",
        "pde_validated",
    ):
        if states.get(key) is not False:
            raise ValueError(f"{key} must remain fail-closed")
    if states.get("segmented_leading_velocity_router_ready") is not True:
        raise ValueError("segmented velocity router unexpectedly unavailable")
    if states.get("reference_local_zero_force_residual_gate_met") is not False:
        raise ValueError("reference local zero-force residual failure cannot be promoted")
    if states.get("reference_local_zero_force_divergence_gate_met") is not False:
        raise ValueError("reference local zero-force divergence failure cannot be promoted")
    if states.get("i2_float64_public_repair_robustly_attributable") is not False:
        raise ValueError("I2 float64 module-level repair attribution must remain fail-closed")
    coverage = unsigned.get("segmented_leading_contract", {}).get("coverage", {})
    if coverage.get("global_coverage_complete") is not False:
        raise ValueError("global coverage cannot be promoted")
    if not coverage.get("unreconstructed_log_X_gaps"):
        raise ValueError("unreconstructed radial gaps must remain explicit")
    truth = unsigned.get("truth_boundary", {})
    if truth.get("free_residual_defined_forcing_used") is not False:
        raise ValueError("free residual-defined forcing is forbidden")
    if truth.get("threshold_relaxed") is not False:
        raise ValueError("threshold relaxation is forbidden")
    if truth.get("paper_exact") is not False or truth.get("openai_field_identified") is not False:
        raise ValueError("source-truth boundary was promoted")
    return payload


def write_bundle(output_dir: str | Path) -> dict[str, Any]:
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    audit = run_agent4_audit()
    _validate_agent4_audit(audit)
    checkpoint = build_checkpoint(audit)
    validate_checkpoint(checkpoint)
    audit_path = out / "agent4_segmented_public_contract_audit.json"
    checkpoint_path = out / "segmented_frontier_routing_checkpoint.json"
    audit_path.write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    checkpoint_path.write_text(
        json.dumps(checkpoint, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    reloaded = json.loads(checkpoint_path.read_text(encoding="utf-8"))
    validate_checkpoint(reloaded)
    return checkpoint


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        default="artifacts/kokuno_agent5/segmented_frontier_routing_checkpoint_v16",
    )
    args = parser.parse_args(argv)
    checkpoint = write_bundle(args.output_dir)
    print(json.dumps(checkpoint, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
