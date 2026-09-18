"""Kokuno Agent-5 integration checkpoint after the reserved-patch audit.

This module adds no new Kokuno mathematics. It binds two independently checked
local prerequisites into one fail-closed routing artifact:

* Agent 1 PR #300 now supplies source-derived ``X_star/e_star`` scales and a
  serializable local I2 reserved-patch velocity component.
* Agent 4 PR #303 independently reconstructs those scales and the local source
  swirl through a different numerical path.
* The corrected heat-discrepancy V2 / three-bump repair machinery inherited
  through Agent-5 PR #295 remains locally validated.

The missing seam is now narrower: the terminal-tail schedule still has to be
bound to the same outer schedule and the corrected V2 discrepancy still has to
be applied to the I2 patch before a globally matched leading velocity/pressure
candidate exists. Agent-2 PR #301 and Agent-3 PR #302 are sibling evidence and
are recorded as provenance only; they are not laundered into this executable
ancestry.

Nothing in this checkpoint assesses the full-domain Navier--Stokes gate.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_heat_v2_routing_checkpoint import run_bound_heat_v2_audit
from .kokuno_outer_patch_independent_audit import run_audit as run_outer_patch_audit
from .kokuno_outer_reserved_patch_schedule import KokunoOuterReservedPatchSchedule

SCHEMA = "kokuno-agent5-outer-patch-routing-checkpoint-v11"
TASK_ID = "KOKUNO-A5-OUTER-PATCH-ROUTING-CHECKPOINT-011"
BASE_AGENT4_HEAD = "dc1ae2a111d0ea49c75d42fdd191d87b0abb989c"

FIXED_GATES = {
    "held_out_normalized_full_momentum_residual": 1.0e-3,
    "divergence_max": 1.0e-5,
    "changed": False,
}

ST006_REFERENCE = {
    "candidate_sha256": "6b4d84b48ab9dbcd2ee1a1858d3e56ef81523f5864369d7e96c6431fccf107a3",
    "validation_seed": 9172801,
    "held_out_cartesian_points": 4096,
    "momentum_sampled_max": 0.1082289305112118,
    "volume_l2": 0.10758432876230622,
    "pde_validated": False,
}

UPSTREAM = {
    "agent1": {
        "pr": 300,
        "head_sha": "8aab23c07e0d239624d7c72bb9c3e6577674d984",
        "standard_run": 35299727733,
        "reserved_patch_schedule_ready": True,
        "X_star_e_star_binding_ready": True,
        "local_I2_patch_velocity_ready": True,
        "terminal_tail_schedule_bound": False,
        "actual_heat_discrepancy_applied_to_patch": False,
        "global_leading_profile_reconstructed": False,
        "consumed_in_executable_ancestry": True,
    },
    "agent2": {
        "pr": 301,
        "head_sha": "1f79a308bc06f0456e1f756c210df2296b8cab22",
        "standard_run": 35300128056,
        "projected_pulse_sensitivity_contract_ready": True,
        "source_actual_D_r_path_instantiated": False,
        "source_actual_D_z_path_instantiated": False,
        "second_public_covariance_column_realized": False,
        "consumed_in_executable_ancestry": False,
    },
    "agent3": {
        "pr": 302,
        "head_sha": "7da46de54d8643b982e9f98e1e1a0eda082fbf14",
        "dedicated_run": 35300304572,
        "standard_run": 35300304569,
        "second_column_bounded_inverse_preflight_ready": True,
        "missing_relative_vector_rms": 0.2101945021227428,
        "second_public_covariance_column_available": False,
        "finite_correction_cycle_run": False,
        "consumed_in_executable_ancestry": False,
    },
    "agent4": {
        "pr": 303,
        "head_sha": BASE_AGENT4_HEAD,
        "dedicated_run": 35301062354,
        "standard_run": 35301062359,
        "artifact_id": 10529956826,
        "artifact_digest": "sha256:dbfe22c6df90fbecfb608d35dd09c8ed2b05b172ddaca251539562532b0bfe3a",
        "outer_reserved_patch_independently_validated": True,
        "consumed_in_executable_ancestry": True,
    },
    "previous_agent5": {
        "pr": 295,
        "head_sha": "4b206974e9912dbd31ddcb2b4dfe279ff49389f9",
        "corrected_heat_v2_routing_ready": True,
        "consumed_in_executable_ancestry": True,
    },
}

STATES = {
    "corrected_heat_discrepancy_independently_validated": True,
    "reserved_patch_schedule_ready": True,
    "X_star_e_star_binding_ready": True,
    "local_I2_patch_velocity_ready": True,
    "local_I2_patch_velocity_independently_validated": True,
    "terminal_tail_schedule_bound_to_outer_schedule": False,
    "physical_three_bump_target_ready": False,
    "heat_compensation_completed": False,
    "global_leading_profile_reconstructed": False,
    "source_pulse_sensitivity_contract_upstream_ready": True,
    "source_actual_pulse_derivative_paths_ready": False,
    "second_column_bounded_inverse_preflight_upstream_ready": True,
    "second_public_covariance_column_realized": False,
    "leading_ready": False,
    "oscillatory_ready": True,
    "correction_ready": False,
    "velocity_export_ready": False,
    "formal_full_domain_pde_gate_assessed": False,
    "pde_validated": False,
}


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _sha(payload: dict[str, Any]) -> str:
    return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()


def bind_outer_patch_audit() -> dict[str, Any]:
    """Rerun Agent 4's independent audit and retain only routing-level evidence."""

    report = run_outer_patch_audit()
    if report.get("structural_preflight_passed") is not True:
        raise RuntimeError("outer reserved-patch independent preflight regressed")
    formal = report.get("formal_project_gate", {})
    if formal.get("assessed") is not False:
        raise ValueError("local outer-patch audit cannot become the formal PDE gate")
    if formal.get("normalized_full_momentum_threshold") != 1.0e-3:
        raise ValueError("formal momentum threshold changed")
    if formal.get("divergence_max_threshold") != 1.0e-5:
        raise ValueError("formal divergence threshold changed")
    st006 = report.get("st006_cross_route_boundary", {})
    if st006.get("directly_comparable") is not False:
        raise ValueError("local source-map audit cannot be relabeled ST006-comparable")

    cases = report.get("cases", [])
    if len(cases) != 3 or not all(case.get("passed") is True for case in cases):
        raise RuntimeError("unexpected outer-patch independent case matrix")
    max_velocity_relative = max(
        float(case["max_public_vs_independent_velocity_relative_error"]) for case in cases
    )
    max_log_scale_abs = max(
        float(case["max_public_vs_independent_log_scale_abs_error"]) for case in cases
    )
    finest_key = "2.0e-04"
    max_finest_divergence = max(
        float(case["normalized_divergence_by_relative_step"][finest_key]["max"])
        for case in cases
    )
    min_mutation_divergence = min(
        float(case["mutation_finest_normalized_divergence_min"]) for case in cases
    )
    if max_velocity_relative >= 2.0e-10:
        raise RuntimeError("independent public/source velocity agreement regressed")
    if max_log_scale_abs >= 2.0e-11:
        raise RuntimeError("independent source-scale agreement regressed")
    if max_finest_divergence >= 2.0e-5:
        raise RuntimeError("local solenoidality guard regressed")
    if min_mutation_divergence < 5.0e-4:
        raise RuntimeError("independent divergence mutation is no longer detected")

    return {
        "task_id": report["task_id"],
        "seed": report["seed"],
        "parameter_case_count": len(cases),
        "max_public_vs_independent_velocity_relative_error": max_velocity_relative,
        "max_public_vs_independent_log_scale_abs_error": max_log_scale_abs,
        "max_finest_local_normalized_divergence": max_finest_divergence,
        "min_mutation_normalized_divergence": min_mutation_divergence,
        "structural_preflight_passed": True,
        "formal_full_domain_pde_gate_assessed": False,
        "st006_directly_comparable": False,
    }


def build_checkpoint(
    schedule: KokunoOuterReservedPatchSchedule,
    outer_audit: dict[str, Any],
    heat_audit: dict[str, Any],
) -> dict[str, Any]:
    scales = schedule.log_scale_report()
    repair = schedule.repair_scales()
    if not (np.isfinite(repair["X_star"]) and repair["X_star"] > 0.0):
        raise ValueError("X_star is not a positive finite physical scale")
    if not (np.isfinite(repair["e_star"]) and repair["e_star"] > 0.0):
        raise ValueError("e_star is not a positive finite physical scale")

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "base_agent4_head": BASE_AGENT4_HEAD,
        "upstream": UPSTREAM,
        "fixed_gates": FIXED_GATES,
        "st006_reference": ST006_REFERENCE,
        "typed_component": {
            "schema": "kokuno-leading-component-capsule-v1",
            "component_kind": "source_reserved_patch_I2",
            "schedule_schema": schedule.to_payload()["schema"],
            "schedule_sha256": schedule.sha256,
            "velocity_api": "KokunoOuterReservedPatchSchedule.patch_velocity(x,y,z,t)->[...,3]",
            "serialization_api": "KokunoOuterReservedPatchSchedule.save/load",
            "pressure_api_ready": False,
            "forcing_api_ready": False,
            "domain_scope": "open source reserved interval I2 only",
            "full_domain_candidate": False,
            "matlab_export_ready": False,
        },
        "derived_default_scales": {
            "log_X_star": float(scales["log_X_star"]),
            "log_e_star": float(scales["log_e_star"]),
            "X_star": float(repair["X_star"]),
            "e_star": float(repair["e_star"]),
            "I2_log_interval": list(schedule.reserved_log_intervals()["I2"]),
            "autonomous_choice": True,
            "hidden_source_numbers_recovered": False,
        },
        "independent_outer_patch_audit": outer_audit,
        "independent_heat_v2_audit": heat_audit,
        "baseline_vs_kokuno": [
            {
                "name": "ST006",
                "scope": "held-out full-domain normalized momentum",
                "momentum_sampled_max": ST006_REFERENCE["momentum_sampled_max"],
                "volume_l2": ST006_REFERENCE["volume_l2"],
                "pde_validated": False,
            },
            {
                "name": "Kokuno local I2 reserved patch",
                "scope": "local source-map/solenoidality implementation preflight",
                "momentum_sampled_max": None,
                "volume_l2": None,
                "pde_validated": False,
                "directly_comparable_to_ST006": False,
                "reason": "no globally matched velocity/pressure/forcing candidate exists",
            },
        ],
        "states": STATES,
        "routing": {
            "repeat_X_star_e_star_search": False,
            "repeat_outer_patch_source_map_validation": False,
            "apply_heat_v2_to_patch_before_terminal_tail_binding": False,
            "promote_local_I2_patch_to_global_leading": False,
            "run_formal_full_domain_gate_now": False,
            "rerun_one_column_correction_cycle": False,
            "replace_ST006_as_repository_pde_baseline": False,
            "priority_blocker": (
                "Agent 1 must bind the source terminal-tail schedule to this same outer schedule, "
                "then feed corrected V2 through the existing three-bump repair using the now-bound "
                "X_star/e_star scales and materialize a globally matched leading velocity/pressure"
            ),
            "agent2_blocker": (
                "instantiate actual source D_r/D_z pulse/background paths and materialize a genuine "
                "second public covariance velocity column through the tested complete-curl machinery"
            ),
            "agent3_blocker": (
                "wait for that physical second column, run the bounded-inverse preflight, and only "
                "then materialize/run the two-column finite correction cycle"
            ),
        },
        "truth_boundary": {
            "source_version_and_provenance_preserved": True,
            "agent5_new_kokuno_mathematics_added": False,
            "sibling_evidence_laundered_into_executable_ancestry": False,
            "local_component_promoted_to_global_candidate": False,
            "free_residual_defined_forcing_used": False,
            "pressure_or_forcing_refit_in_agent5": False,
            "threshold_relaxed": False,
            "st006_used_as_kokuno_source_truth": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
        "next_shortest_closure": [
            "Agent 1: bind terminal-tail scales to the same outer schedule; apply corrected V2 through the preflighted three-bump repair; emit globally matched leading velocity/pressure",
            "Agent 2: instantiate actual source D_r/D_z pulse/background paths; materialize a genuine second public covariance column",
            "Agent 3: preflight that physical second column, then run a frozen two-column finite correction cycle",
            "Agent 4: once the global composite exists, run independent held-out normalized momentum <=1e-3 and divergence <=1e-5 validation plus same-protocol ST006 comparison",
            "Agent 5: serialize the complete velocity/pressure/forcing candidate and run Python/MATLAB export smoke only after the global interface exists",
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
    if unsigned.get("base_agent4_head") != BASE_AGENT4_HEAD:
        raise ValueError("base Agent-4 head changed")
    if unsigned.get("fixed_gates") != FIXED_GATES:
        raise ValueError("fixed PDE/divergence gates changed")
    if unsigned.get("st006_reference") != ST006_REFERENCE:
        raise ValueError("ST006 reference changed")
    if unsigned.get("upstream") != UPSTREAM:
        raise ValueError("upstream provenance changed")
    if unsigned.get("states") != STATES:
        raise ValueError("scientific state vector changed")

    typed = unsigned.get("typed_component", {})
    if typed.get("full_domain_candidate") is not False:
        raise ValueError("local I2 component was promoted to a global candidate")
    if typed.get("pressure_api_ready") is not False or typed.get("forcing_api_ready") is not False:
        raise ValueError("missing pressure/forcing APIs were relabeled ready")

    outer = unsigned.get("independent_outer_patch_audit", {})
    if outer.get("structural_preflight_passed") is not True:
        raise ValueError("outer-patch independent preflight is not bound")
    if outer.get("formal_full_domain_pde_gate_assessed") is not False:
        raise ValueError("local outer audit was relabeled as the full PDE gate")
    if outer.get("st006_directly_comparable") is not False:
        raise ValueError("local outer audit was relabeled ST006-comparable")

    heat = unsigned.get("independent_heat_v2_audit", {})
    if heat.get("structural_preflight_passed") is not True:
        raise ValueError("heat-v2 independent preflight is not bound")
    if heat.get("used_as_formal_full_domain_pde_gate") is not False:
        raise ValueError("local heat audit was relabeled as the full PDE gate")

    routing = unsigned.get("routing", {})
    for key in (
        "repeat_X_star_e_star_search",
        "repeat_outer_patch_source_map_validation",
        "apply_heat_v2_to_patch_before_terminal_tail_binding",
        "promote_local_I2_patch_to_global_leading",
        "run_formal_full_domain_gate_now",
        "rerun_one_column_correction_cycle",
        "replace_ST006_as_repository_pde_baseline",
    ):
        if routing.get(key) is not False:
            raise ValueError(f"blocked routing step was promoted: {key}")

    truth = unsigned.get("truth_boundary", {})
    if truth.get("agent5_new_kokuno_mathematics_added") is not False:
        raise ValueError("Agent 5 cannot claim new Kokuno mathematics in this checkpoint")
    if truth.get("free_residual_defined_forcing_used") is not False:
        raise ValueError("free residual-defined forcing is prohibited")
    if truth.get("threshold_relaxed") is not False:
        raise ValueError("fixed project threshold was relaxed")
    if truth.get("paper_exact") is not False or truth.get("openai_field_identified") is not False:
        raise ValueError("truth boundary was promoted")
    return payload


def load_checkpoint(path: str | Path) -> dict[str, Any]:
    return validate_checkpoint(json.loads(Path(path).read_text(encoding="utf-8")))


def write_bundle(output_dir: str | Path) -> dict[str, Any]:
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    schedule = KokunoOuterReservedPatchSchedule()
    schedule_path = out / "outer_reserved_patch_schedule.json"
    schedule.save(schedule_path)
    reloaded = KokunoOuterReservedPatchSchedule.load(schedule_path)
    if reloaded.sha256 != schedule.sha256:
        raise RuntimeError("outer reserved-patch schedule round-trip changed its SHA")

    outer_audit = bind_outer_patch_audit()
    heat_audit = run_bound_heat_v2_audit()
    checkpoint = build_checkpoint(reloaded, outer_audit, heat_audit)
    validate_checkpoint(checkpoint)

    (out / "outer_patch_audit_summary.json").write_text(
        json.dumps(outer_audit, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (out / "heat_v2_audit_summary.json").write_text(
        json.dumps(heat_audit, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    checkpoint_path = out / "outer_patch_routing_checkpoint.json"
    checkpoint_path.write_text(
        json.dumps(checkpoint, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    load_checkpoint(checkpoint_path)
    return checkpoint


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser


def main() -> int:
    args = _parser().parse_args()
    checkpoint = write_bundle(args.output_dir)
    print(json.dumps(checkpoint, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
