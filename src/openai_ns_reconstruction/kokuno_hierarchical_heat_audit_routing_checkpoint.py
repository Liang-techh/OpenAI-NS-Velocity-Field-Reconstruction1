"""Kokuno Agent-5 fail-closed routing checkpoint after the hierarchical heat audit.

This module adds no new Kokuno mathematics.  It binds Agent 1's high-precision
hierarchical I2 repair to Agent 4's independently integrated rejection receipt,
and records the current Agent 2/3 handoff state.  The purpose is to prevent a
frozen-discrete moment closure from being promoted to continuous heat
compensation or to a full-domain Navier--Stokes candidate.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

from .kokuno_hierarchical_heat_repair import KokunoHierarchicalHeatRepair
from .kokuno_hierarchical_heat_repair_independent import build_report as build_independent_report

SCHEMA = "kokuno-agent5-hierarchical-heat-audit-routing-checkpoint-v14"
TASK_ID = "KOKUNO-A5-HIERARCHICAL-HEAT-AUDIT-ROUTING-CHECKPOINT-014"
BASE_AGENT4_HEAD = "749b2fdae1d76787e88448ed15fe44098cd9f723"

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
        "pr": 328,
        "head_sha": "3afbf5a019e09f9a351a49e09803e6d6ac25fc61",
        "standard_run": 35311841481,
        "hierarchical_discrete_three_bump_repair_ready": True,
        "continuous_source_moment_compensation_certified": False,
        "global_leading_profile_reconstructed": False,
        "consumed_in_executable_ancestry": True,
    },
    "agent2": {
        "pr": 330,
        "head_sha": "9f077f206eb7abf73ec24c35e8214863dc250417",
        "dedicated_run": 35311856780,
        "standard_run": 35311856777,
        "source_homogeneous_primary_pulse_ready": True,
        "caller_supplied_source_path_geometry_still_required": True,
        "public_xyz_t_velocity_correction_materialized": False,
        "genuinely_independent_second_covariance_column_ready": False,
        "consumed_in_executable_ancestry": False,
    },
    "agent3": {
        "pr": 329,
        "head_sha": "42cd3ae4f70ad1260930bb24b94e756a70c487b9",
        "dedicated_run": 35311643130,
        "standard_run": 35311643070,
        "held_out_times": [0.375, 0.5, 0.625],
        "node_time_requirements_total": 75,
        "rank_two_requirements_passed": 0,
        "missing_relative_vector_rms": 0.21019450214,
        "finite_cycle_rerun_allowed": False,
        "consumed_in_executable_ancestry": False,
    },
    "agent4": {
        "pr": 331,
        "head_sha": BASE_AGENT4_HEAD,
        "dedicated_run": 35312449473,
        "artifact_id": 10533489918,
        "artifact_digest": "sha256:5555c444cce5d66295d46bae02bca6642a1ea67206aa6a273b7f18d423d3ce51",
        "independent_hierarchical_heat_audit_executed": True,
        "continuous_source_moment_compensation_certified": False,
        "consumed_in_executable_ancestry": True,
    },
    "previous_agent5": {
        "pr": 324,
        "head_sha": "381c0cfb9b2f4076b107d03b9eb45c2df5cf8f31",
        "signed_log_target_routed": True,
        "consumed_in_executable_ancestry": False,
    },
}
STATES = {
    "signed_log_heat_target_ready": True,
    "hierarchical_discrete_three_bump_repair_ready": True,
    "independent_hierarchical_heat_audit_executed": True,
    "continuous_source_moment_compensation_certified": False,
    "heat_compensation_completed": False,
    "global_leading_profile_reconstructed": False,
    "source_homogeneous_primary_pulse_upstream_ready": True,
    "public_source_oscillatory_xyz_t_velocity_ready": False,
    "second_public_covariance_column_realized": False,
    "finite_cycle_rerun_allowed": False,
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


def bind_hierarchical_heat_audit() -> dict[str, Any]:
    """Rerun Agent-4's independent audit and preserve its scientific rejection."""

    report = build_independent_report()
    summary = report["summary"]
    guards = report["local_guards"]
    if report["schema"] != "kokuno-agent4-hierarchical-heat-repair-independent-v1":
        raise RuntimeError("unexpected Agent-4 audit schema")
    if summary["independent_continuous_moment_closure"] is not False:
        raise RuntimeError("Agent-4 rejection unexpectedly disappeared")
    if summary["frozen_discrete_map_only"] is not True:
        raise RuntimeError("frozen-discrete truth boundary changed")
    if report["local_audit_completed"] is not False:
        raise RuntimeError("failed independent closure guard was silently promoted")
    if guards["dominant_I_sub_channel_relative_error_le_1e-10"] is not False:
        raise RuntimeError("I_sub independent precision barrier unexpectedly vanished")
    if float(summary["max_I_sub_relative_error"]) <= 1.0e-10:
        raise RuntimeError("independent I_sub mismatch no longer demonstrates the blocker")
    if float(summary["minimum_Cp_log10_relative_error"]) < 50.0:
        raise RuntimeError("C_p negative-control separation regressed")
    if float(summary["minimum_S_log10_relative_error"]) < 50.0:
        raise RuntimeError("S negative-control separation regressed")
    formal = report["formal_project_gates"]
    if formal["normalized_momentum_max"] != 1.0e-3 or formal["divergence_max"] != 1.0e-5:
        raise ValueError("formal project gates changed")
    if formal["assessed_here"] is not False or formal["pde_validated"] is not False:
        raise ValueError("local heat audit cannot promote the full-domain PDE gate")
    return {
        "schema": report["schema"],
        "summary": summary,
        "local_guards": guards,
        "local_audit_completed": False,
        "scientific_result": "REJECT_CONTINUOUS_MOMENT_CLOSURE",
        "formal_full_domain_pde_gate_assessed": False,
        "pde_validated": False,
    }


def build_checkpoint(audit: dict[str, Any]) -> dict[str, Any]:
    repair = KokunoHierarchicalHeatRepair()
    solution = repair.solve(0.2)
    if tuple(solution.target_sign) != (1, -1, 1):
        raise RuntimeError("Agent-1 signed-log target signs changed")
    if solution.precision_digits < 480:
        raise RuntimeError("Agent-1 adaptive Decimal precision unexpectedly regressed")
    if not math.isfinite(solution.max_relative_residual) or solution.max_relative_residual >= 1.0e-40:
        raise RuntimeError("Agent-1 frozen-discrete nonlinear closure regressed")

    independent_i_sub_relative = float(audit["summary"]["max_I_sub_relative_error"])
    if independent_i_sub_relative <= 0.0:
        raise RuntimeError("independent mismatch must be positive")
    i_sub_guard_gap_orders = math.log10(independent_i_sub_relative / 1.0e-10)

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "base_agent4_head": BASE_AGENT4_HEAD,
        "upstream": UPSTREAM,
        "fixed_gates": FIXED_GATES,
        "st006_reference": ST006_REFERENCE,
        "agent1_discrete_hierarchical_receipt": {
            "eta": solution.eta,
            "solution_sha256": solution.sha256,
            "precision_digits": solution.precision_digits,
            "target_sign": list(solution.target_sign),
            "max_relative_residual": solution.max_relative_residual,
            "newton_steps": solution.newton_steps,
            "scope": "repository frozen discrete three-bump moment tensor only",
            "continuous_source_moment_compensation_certified": False,
        },
        "agent4_independent_audit": audit,
        "construction_vs_independent_precision_gap": {
            "frozen_discrete_relative_residual_lt": 1.0e-40,
            "independent_I_sub_relative_error": independent_i_sub_relative,
            "independent_I_sub_guard": 1.0e-10,
            "guard_excess_orders": i_sub_guard_gap_orders,
            "minimum_Cp_log10_relative_error": float(audit["summary"]["minimum_Cp_log10_relative_error"]),
            "minimum_S_log10_relative_error": float(audit["summary"]["minimum_S_log10_relative_error"]),
            "interpretation": (
                "high-precision coefficient arithmetic is no longer the blocker; the underlying "
                "moment tensor/integration accuracy is. The current coefficients close the frozen "
                "float64 tensor but fail an independently integrated continuous-moment approximation."
            ),
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
                "name": "Kokuno hierarchical I2 heat repair",
                "scope": "local source heat-moment precision audit",
                "momentum_sampled_max": None,
                "volume_l2": None,
                "directly_comparable_to_ST006": False,
                "pde_validated": False,
            },
        ],
        "states": STATES,
        "routing": {
            "repeat_signed_log_target_derivation": False,
            "repeat_float64_tensor_decimal_newton_only": False,
            "promote_current_hierarchical_repair_to_heat_compensation": False,
            "run_formal_full_domain_gate_now": False,
            "rerun_one_or_duplicate_column_finite_cycle": False,
            "replace_ST006_as_repository_pde_baseline": False,
            "priority_blocker": (
                "Agent 1: rebuild the three-bump moment tensor/integrals at precision commensurate "
                "with the signed-log channels, or reformulate the hierarchy to avoid cancellation "
                "through float64 tensor coefficients; then rerun Agent 4's independent alternate-"
                "quadrature audit before any global leading promotion"
            ),
            "agent2_blocker": (
                "instantiate actual source path/background geometry for the homogeneous primary "
                "pulse, propagate the needed source sensitivities, expose direct public Q-scaled "
                "velocity(x,y,z,t), and provide a genuinely independent second covariance column"
            ),
            "agent3_blocker": (
                "screen the next genuinely independent Agent-2 physical column on all three frozen "
                "held-out slices; current duplicate control is 0/75 rank-two, so no finite cycle"
            ),
            "agent4_next": (
                "repeat the independent alternate-quadrature heat audit only after Agent 1 changes "
                "the underlying moment representation; run the fixed full-domain PDE gate only when "
                "a global composite velocity/pressure/forcing candidate exists"
            ),
            "agent5_next": (
                "bind the independently accepted heat repair into the unified typed candidate schema "
                "only after continuous moment closure passes; do not export a partial local I2 field "
                "as the global Kokuno candidate"
            ),
        },
        "truth_boundary": {
            "source_version_and_provenance_preserved": True,
            "agent5_new_kokuno_mathematics_added": False,
            "scientific_rejection_preserved_in_green_CI": True,
            "sibling_evidence_laundered_into_executable_ancestry": False,
            "local_heat_repair_promoted_to_global_candidate": False,
            "free_residual_defined_forcing_used": False,
            "pressure_or_forcing_refit_in_agent5": False,
            "threshold_relaxed": False,
            "st006_used_as_kokuno_source_truth": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
        "next_shortest_closure": [
            "Agent 1: high-precision/precision-stable moment tensor -> hierarchical repair -> Agent 4 independent closure pass",
            "Agent 1: only after that pass, assemble repaired I2 into global leading velocity/pressure",
            "Agent 2: actual source path/background -> public real-pair xyz,t velocity -> independent second covariance column",
            "Agent 3: multi-slice second-column preflight -> bounded two-column finite cycle only if all guards pass",
            "Agent 4: independent global held-out <=1e-3 momentum and <=1e-5 divergence gate",
            "Agent 5: final candidate serialization plus minimal Python/MATLAB export smoke after global API exists",
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
    if unsigned.get("fixed_gates") != FIXED_GATES or unsigned.get("st006_reference") != ST006_REFERENCE:
        raise ValueError("fixed gates or ST006 reference changed")
    if unsigned.get("upstream") != UPSTREAM or unsigned.get("states") != STATES:
        raise ValueError("upstream provenance or scientific state changed")
    audit = unsigned.get("agent4_independent_audit", {})
    if audit.get("scientific_result") != "REJECT_CONTINUOUS_MOMENT_CLOSURE":
        raise ValueError("independent scientific rejection missing")
    if audit.get("local_audit_completed") is not False:
        raise ValueError("failed independent audit was promoted")
    if unsigned.get("agent1_discrete_hierarchical_receipt", {}).get("continuous_source_moment_compensation_certified") is not False:
        raise ValueError("frozen-discrete closure was promoted")
    rows = unsigned.get("baseline_vs_kokuno", [{}, {}])
    if len(rows) != 2 or rows[1].get("directly_comparable_to_ST006") is not False:
        raise ValueError("local heat audit was relabeled ST006-comparable")
    routing = unsigned.get("routing", {})
    if routing.get("run_formal_full_domain_gate_now") is not False:
        raise ValueError("formal PDE gate routed before a global candidate exists")
    if routing.get("promote_current_hierarchical_repair_to_heat_compensation") is not False:
        raise ValueError("rejected heat repair was routed for promotion")
    if routing.get("rerun_one_or_duplicate_column_finite_cycle") is not False:
        raise ValueError("rank-deficient correction route was reopened")
    truth = unsigned.get("truth_boundary", {})
    if truth.get("threshold_relaxed") is not False or truth.get("free_residual_defined_forcing_used") is not False:
        raise ValueError("scientific contract weakened")
    if any(unsigned["states"][key] for key in (
        "leading_ready",
        "correction_ready",
        "velocity_export_ready",
        "formal_full_domain_pde_gate_assessed",
        "pde_validated",
    )):
        raise ValueError("premature scientific state promotion")
    return payload


def write_bundle(output_dir: str | Path) -> dict[str, Any]:
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    audit = bind_hierarchical_heat_audit()
    checkpoint = build_checkpoint(audit)
    validate_checkpoint(checkpoint)
    path = out / "hierarchical_heat_audit_routing_checkpoint.json"
    path.write_text(json.dumps(checkpoint, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    return checkpoint


def load_checkpoint(path: str | Path) -> dict[str, Any]:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    return validate_checkpoint(payload)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Build the Kokuno Agent-5 hierarchical heat audit routing checkpoint")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    checkpoint = write_bundle(args.output_dir)
    summary = {
        "checkpoint_sha256": checkpoint["checkpoint_sha256"],
        "scientific_result": checkpoint["agent4_independent_audit"]["scientific_result"],
        "independent_I_sub_relative_error": checkpoint["construction_vs_independent_precision_gap"]["independent_I_sub_relative_error"],
        "guard_excess_orders": checkpoint["construction_vs_independent_precision_gap"]["guard_excess_orders"],
        "leading_ready": checkpoint["states"]["leading_ready"],
        "pde_validated": checkpoint["states"]["pde_validated"],
    }
    print(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
