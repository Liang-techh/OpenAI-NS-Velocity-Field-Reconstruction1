"""Kokuno Agent-5 routing checkpoint after Agent-1's high-precision moment rebuild.

This module adds no Kokuno mathematics.  It harvests the new Agent-1 I2
three-bump moment backend, retains Agent-4's previous independent rejection as
historical evidence, performs only the ordinary-precision dominant-channel
cross-check that the new construction is already capable of supporting, and
fails closed on the still-missing independent *high-precision* alternate-
quadrature certification of the tiny C_p/S channels.

The checkpoint is integration/routing evidence, not a Navier--Stokes gate.
"""
from __future__ import annotations

import argparse
from decimal import Decimal
import hashlib
import json
import math
from pathlib import Path
from typing import Any

from .kokuno_high_precision_heat_repair import KokunoHighPrecisionHierarchicalHeatRepair
from .kokuno_hierarchical_heat_repair import KokunoHierarchicalHeatRepair
from .kokuno_hierarchical_heat_repair_independent import (
    _map_decimal,
    build_report as build_legacy_independent_report,
    independent_moment_tensors,
)

SCHEMA = "kokuno-agent5-high-precision-heat-handoff-routing-checkpoint-v15"
TASK_ID = "KOKUNO-A5-HIGHPREC-HEAT-HANDOFF-ROUTING-CHECKPOINT-015"
BASE_AGENT1_HEAD = "887fe6939a3d6fd09b40a090b671136a6d2e8b13"

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
        "pr": 337,
        "head_sha": BASE_AGENT1_HEAD,
        "dedicated_run": 35315916416,
        "standard_run": 35315916378,
        "high_precision_continuous_moment_tensor_used": True,
        "float64_moment_tensor_used": False,
        "independent_high_precision_moment_audit_completed": False,
        "continuous_source_moment_compensation_certified": False,
        "global_leading_profile_reconstructed": False,
        "consumed_in_executable_ancestry": True,
    },
    "agent2": {
        "pr": 336,
        "head_sha": "acfa0153cf4d2a8deb080302688d72895fe06180",
        "dedicated_run": 35315632374,
        "standard_run": 35315632353,
        "homogeneous_pulse_spatial_sensitivity_ready": True,
        "actual_positive_order_source_path_instantiated": False,
        "public_xyz_t_velocity_correction_materialized": False,
        "genuinely_independent_second_covariance_column_ready": False,
        "consumed_in_executable_ancestry": False,
    },
    "agent3": {
        "pr": 338,
        "head_sha": "12f34deb99260b0457484ee0d38cd41e6428a549",
        "spacetime_screen_times": [0.375, 0.5, 0.625],
        "spacetime_screen_z": [0.06, 0.08, 0.10],
        "actual_second_public_covariance_column_available": False,
        "finite_cycle_rerun_allowed": False,
        "consumed_in_executable_ancestry": False,
    },
    "agent4": {
        "pr": 331,
        "head_sha": "749b2fdae1d76787e88448ed15fe44098cd9f723",
        "dedicated_run": 35312449473,
        "standard_run": 35312452405,
        "legacy_float64_tensor_rejection_preserved": True,
        "new_agent1_high_precision_backend_independently_audited": False,
        "consumed_in_executable_ancestry": True,
    },
    "previous_agent5": {
        "pr": 332,
        "head_sha": "9fc32ca30bcf9813b672ba977a7ced6ccd4fbbed",
        "routed_high_precision_tensor_rebuild": True,
        "consumed_in_executable_ancestry": True,
    },
}
STATES = {
    "signed_log_heat_target_ready": True,
    "high_precision_three_bump_moment_tensor_ready": True,
    "high_precision_hierarchical_repair_executable": True,
    "dominant_I_sub_ordinary_precision_crosscheck_passed": True,
    "independent_high_precision_moment_audit_completed": False,
    "continuous_source_moment_compensation_certified": False,
    "heat_compensation_completed": False,
    "global_leading_profile_reconstructed": False,
    "homogeneous_pulse_spatial_sensitivity_ready": True,
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


def _relative_decimal(actual: Decimal, target: Decimal) -> Decimal:
    if target == 0:
        return abs(actual)
    return abs(actual - target) / abs(target)


def build_dominant_channel_crosscheck(eta: float = 0.2) -> dict[str, Any]:
    """Check the newly rebuilt tensor only where float64 alternate quadrature can speak.

    Agent-4's existing Simpson tensor is ordinary precision.  It is adequate to
    detect the old dominant I_sub error but cannot certify C_p/S cancellation
    hundreds of decades below float64.  Those channels remain deliberately
    pending an Agent-4 high-precision alternate operator.
    """

    repair = KokunoHighPrecisionHierarchicalHeatRepair()
    solution = repair.solve(float(eta))
    if tuple(solution.target_sign) != (1, -1, 1):
        raise RuntimeError("signed-log target signs changed")
    if solution.precision_digits < 480:
        raise RuntimeError("high-precision solve lost the required hierarchy")
    if not math.isfinite(solution.max_relative_residual) or solution.max_relative_residual >= 1.0e-40:
        raise RuntimeError("construction-side high-precision closure regressed")

    independent = independent_moment_tensors(
        lambda_outer=repair.outer_schedule.lambda_outer,
        f_eta=float(repair.outer_schedule.source_f(float(eta))),
        resolution=4096,
    )
    coefficients = solution.coefficient_decimals()
    achieved = _map_decimal(coefficients, independent, solution.precision_digits)
    target = tuple(Decimal(value) for value in solution.target)
    i_sub_relative = _relative_decimal(achieved[2], target[2])

    legacy = KokunoHierarchicalHeatRepair()
    legacy_solution = legacy.solve(float(eta))
    legacy_achieved = _map_decimal(
        legacy_solution.coefficient_decimals(), independent, legacy_solution.precision_digits
    )
    legacy_target = tuple(Decimal(value) for value in legacy_solution.target)
    legacy_i_sub_relative = _relative_decimal(legacy_achieved[2], legacy_target[2])

    if i_sub_relative >= Decimal("1e-10"):
        raise RuntimeError("new high-precision tensor still fails the dominant I_sub cross-check")
    if i_sub_relative >= legacy_i_sub_relative * Decimal("1e-3"):
        raise RuntimeError("new tensor did not improve the legacy I_sub mismatch by >=1000x")

    precision_report = repair.precision_report(float(eta))
    if precision_report.get("independent_high_precision_moment_audit_completed") is not False:
        raise RuntimeError("construction incorrectly claims an independent high-precision audit")
    if precision_report.get("continuous_moment_compensation_certified") is not False:
        raise RuntimeError("construction incorrectly claims continuous compensation certification")

    return {
        "eta": float(eta),
        "solution_sha256": solution.sha256,
        "precision_digits": int(solution.precision_digits),
        "target_sign": list(solution.target_sign),
        "construction_max_relative_residual": float(solution.max_relative_residual),
        "alternate_operator": "Agent-4 4096-interval composite-Simpson moment tensor (ordinary precision)",
        "new_I_sub_relative_error": float(i_sub_relative),
        "legacy_I_sub_relative_error": float(legacy_i_sub_relative),
        "improvement_factor": float(legacy_i_sub_relative / i_sub_relative),
        "dominant_I_sub_guard": 1.0e-10,
        "dominant_I_sub_guard_passed": True,
        "tiny_Cp_S_independently_high_precision_certified": False,
        "independent_high_precision_moment_audit_completed": False,
        "continuous_source_moment_compensation_certified": False,
        "scope": (
            "ordinary-precision alternate-operator cross-check of the dominant I_sub channel only; "
            "not certification of the hundreds-of-decades C_p/S cancellations"
        ),
    }


def bind_previous_agent4_rejection() -> dict[str, Any]:
    """Preserve the independent rejection that motivated Agent-1 #337."""

    report = build_legacy_independent_report()
    summary = report["summary"]
    if summary["independent_continuous_moment_closure"] is not False:
        raise RuntimeError("previous Agent-4 rejection unexpectedly disappeared")
    old_error = float(summary["max_I_sub_relative_error"])
    if old_error <= 1.0e-10:
        raise RuntimeError("previous dominant-channel rejection no longer reproduces")
    return {
        "schema": report["schema"],
        "legacy_max_I_sub_relative_error": old_error,
        "legacy_dominant_guard": 1.0e-10,
        "legacy_continuous_source_moment_compensation_certified": False,
        "scientific_result": "REJECT_LEGACY_FLOAT64_MOMENT_TENSOR",
    }


def build_checkpoint() -> dict[str, Any]:
    dominant = build_dominant_channel_crosscheck()
    previous = bind_previous_agent4_rejection()
    if dominant["new_I_sub_relative_error"] >= previous["legacy_max_I_sub_relative_error"]:
        raise RuntimeError("new moment backend did not improve the previous independent mismatch")

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "base_agent1_head": BASE_AGENT1_HEAD,
        "upstream": UPSTREAM,
        "fixed_gates": FIXED_GATES,
        "st006_reference": ST006_REFERENCE,
        "previous_agent4_rejection": previous,
        "agent1_high_precision_handoff": dominant,
        "baseline_vs_kokuno": [
            {
                "name": "ST006",
                "scope": "held-out full-domain normalized NS momentum",
                "momentum_sampled_max": ST006_REFERENCE["momentum_sampled_max"],
                "volume_l2": ST006_REFERENCE["volume_l2"],
                "pde_validated": False,
            },
            {
                "name": "Kokuno high-precision I2 heat repair handoff",
                "scope": "local source heat-moment construction/cross-check",
                "momentum_sampled_max": None,
                "volume_l2": None,
                "directly_comparable_to_ST006": False,
                "pde_validated": False,
            },
        ],
        "states": STATES,
        "routing": {
            "repeat_float64_tensor_decimal_newton": False,
            "repeat_high_precision_tensor_rebuild_without_new_evidence": False,
            "promote_repaired_I2_to_completed_heat_compensation_now": False,
            "run_formal_full_domain_gate_now": False,
            "rerun_one_or_duplicate_column_finite_cycle": False,
            "replace_ST006_as_repository_pde_baseline": False,
            "priority_blocker": (
                "Agent 4 must independently audit Agent-1 #337 with a genuinely high-precision "
                "alternate quadrature/operator over all three signed-log channels, especially the "
                "C_p/S cancellations that ordinary float64 Simpson cannot certify"
            ),
            "agent1_after_agent4_pass": (
                "only after the independent high-precision all-channel audit passes, assemble the "
                "repaired I2 into the global matched leading velocity/pressure and emit a typed, "
                "serializable public API"
            ),
            "agent2_blocker": (
                "bind the actual positive-order/source-background path and derivatives to the "
                "homogeneous pulse machinery, expose direct Q-scaled velocity(x,y,z,t), and supply "
                "a genuinely independent second covariance velocity column"
            ),
            "agent3_blocker": (
                "screen that real second column through the existing t x z spacetime guard before "
                "any bounded two-column finite correction cycle; no actual second column exists yet"
            ),
            "agent4_global_gate_after_composite": (
                "when the complete global leading+oscillatory+correction candidate exists, run the "
                "fixed held-out <=1e-3 normalized momentum and <=1e-5 divergence gate and a "
                "same-protocol ST006 comparison"
            ),
            "agent5_next": (
                "after Agent-4 high-precision heat certification, bind the repaired leading API into "
                "the unified candidate schema; final MATLAB/Python export remains gated on a global field"
            ),
        },
        "truth_boundary": {
            "source_version_and_provenance_preserved": True,
            "agent5_new_kokuno_mathematics_added": False,
            "ordinary_precision_I_sub_crosscheck_is_not_full_high_precision_certification": True,
            "tiny_Cp_S_channels_still_pending_independent_high_precision_audit": True,
            "sibling_evidence_laundered_into_executable_ancestry": False,
            "partial_I2_field_promoted_to_global_candidate": False,
            "free_residual_defined_forcing_used": False,
            "pressure_or_forcing_refit_in_agent5": False,
            "threshold_relaxed": False,
            "st006_used_as_kokuno_source_truth": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
        "next_shortest_closure": [
            "Agent 4: independent high-precision alternate-quadrature audit of Agent-1 #337 over C_p/S/I_sub",
            "Agent 1: after that pass only, global repaired leading velocity/pressure assembly and serialization",
            "Agent 2: actual source background/path -> direct public xyz,t oscillatory velocity -> independent second covariance column",
            "Agent 3: t x z spacetime second-column guard -> bounded two-column finite cycle only if all cells pass",
            "Agent 4: independent full-domain held-out <=1e-3 momentum / <=1e-5 divergence validation",
            "Agent 5: unified final candidate artifact plus minimal Python/MATLAB export smoke after the global API exists",
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
    if unsigned.get("base_agent1_head") != BASE_AGENT1_HEAD:
        raise ValueError("base Agent-1 head changed")
    if unsigned.get("upstream") != UPSTREAM or unsigned.get("states") != STATES:
        raise ValueError("upstream provenance or states changed")
    if unsigned.get("fixed_gates") != FIXED_GATES or unsigned.get("st006_reference") != ST006_REFERENCE:
        raise ValueError("fixed gates or ST006 reference changed")
    handoff = unsigned.get("agent1_high_precision_handoff", {})
    if handoff.get("dominant_I_sub_guard_passed") is not True:
        raise ValueError("dominant I_sub handoff guard missing")
    if handoff.get("tiny_Cp_S_independently_high_precision_certified") is not False:
        raise ValueError("tiny channels were prematurely certified")
    if handoff.get("independent_high_precision_moment_audit_completed") is not False:
        raise ValueError("Agent-4 high-precision audit was invented")
    if handoff.get("continuous_source_moment_compensation_certified") is not False:
        raise ValueError("continuous heat compensation was prematurely certified")
    rows = unsigned.get("baseline_vs_kokuno", [{}, {}])
    if len(rows) != 2 or rows[1].get("directly_comparable_to_ST006") is not False:
        raise ValueError("local heat handoff was relabeled ST006-comparable")
    routing = unsigned.get("routing", {})
    if routing.get("promote_repaired_I2_to_completed_heat_compensation_now") is not False:
        raise ValueError("I2 repair promoted before independent high-precision audit")
    if routing.get("run_formal_full_domain_gate_now") is not False:
        raise ValueError("formal PDE gate routed before a global candidate exists")
    if routing.get("rerun_one_or_duplicate_column_finite_cycle") is not False:
        raise ValueError("rank-deficient correction route was reopened")
    truth = unsigned.get("truth_boundary", {})
    if truth.get("threshold_relaxed") is not False or truth.get("free_residual_defined_forcing_used") is not False:
        raise ValueError("scientific contract weakened")
    for key in (
        "leading_ready",
        "correction_ready",
        "velocity_export_ready",
        "formal_full_domain_pde_gate_assessed",
        "pde_validated",
    ):
        if unsigned["states"][key] is not False:
            raise ValueError("premature scientific state promotion")
    return payload


def write_bundle(output_dir: str | Path) -> dict[str, Any]:
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    checkpoint = build_checkpoint()
    validate_checkpoint(checkpoint)
    path = out / "high_precision_heat_handoff_routing_checkpoint.json"
    path.write_text(
        json.dumps(checkpoint, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return checkpoint


def load_checkpoint(path: str | Path) -> dict[str, Any]:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    return validate_checkpoint(payload)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Build the Kokuno Agent-5 high-precision heat handoff routing checkpoint"
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    checkpoint = write_bundle(args.output_dir)
    handoff = checkpoint["agent1_high_precision_handoff"]
    print(
        json.dumps(
            {
                "checkpoint_sha256": checkpoint["checkpoint_sha256"],
                "new_I_sub_relative_error": handoff["new_I_sub_relative_error"],
                "legacy_I_sub_relative_error": handoff["legacy_I_sub_relative_error"],
                "improvement_factor": handoff["improvement_factor"],
                "independent_high_precision_moment_audit_completed": handoff[
                    "independent_high_precision_moment_audit_completed"
                ],
                "leading_ready": checkpoint["states"]["leading_ready"],
                "pde_validated": checkpoint["states"]["pde_validated"],
            },
            indent=2,
            sort_keys=True,
            allow_nan=False,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
