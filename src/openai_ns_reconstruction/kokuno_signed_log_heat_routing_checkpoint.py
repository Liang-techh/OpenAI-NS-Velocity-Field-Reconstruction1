"""Kokuno Agent-5 fail-closed routing checkpoint for the signed-log heat frontier.

This module introduces no new Kokuno mathematics.  It binds Agent 1 PR #320,
Agent 4 PR #323's independent audit, and the current Agent 2/3 sibling handoff
state into one deterministic, serializable routing receipt.  The full-domain
Navier--Stokes gates remain unassessed until a global composite candidate exists.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_independent_signed_log_heat_target_audit import run_audit
from .kokuno_log_rescaled_heat_repair_target import KokunoLogRescaledHeatRepairTarget

SCHEMA = "kokuno-agent5-signed-log-heat-routing-checkpoint-v13"
TASK_ID = "KOKUNO-A5-SIGNED-LOG-HEAT-ROUTING-CHECKPOINT-013"
BASE_AGENT4_HEAD = "ff81bf49bc45d3575bd43a52113562354f381b3e"

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
        "pr": 320,
        "head_sha": "2a7d51a05bcf46a2071ed7e767e8f9f467a79f86",
        "standard_run": 35307255084,
        "signed_log_heat_target_ready": True,
        "float64_three_bump_inverse_certified": False,
        "heat_compensation_completed": False,
        "global_leading_profile_reconstructed": False,
        "consumed_in_executable_ancestry": True,
    },
    "agent2": {
        "pr": 321,
        "head_sha": "1c5cf93b7468441a0f735336837a4ecdcb8e5afd",
        "dedicated_run": 35307264120,
        "standard_run": 35307264110,
        "source_real_conjugate_pair_kernel_ready": True,
        "actual_source_pulse_background_instantiated": False,
        "public_xyz_t_velocity_correction_materialized": False,
        "genuinely_independent_second_covariance_column_ready": False,
        "consumed_in_executable_ancestry": False,
    },
    "agent3": {
        "pr": 322,
        "head_sha": "2168b4fe8a7b2c6a177e6f9c4fd9a10428f61f95",
        "dedicated_run": 35307772992,
        "standard_run": 35307772950,
        "velocity_column_covariance_adapter_executable": True,
        "active_target_nodes": 27,
        "nodes_requiring_second_direction": 25,
        "rank_two_required_nodes": 0,
        "duplicate_two_column_relative_stress_residual": 0.21019450212274277,
        "finite_cycle_rerun_allowed": False,
        "consumed_in_executable_ancestry": False,
    },
    "agent4": {
        "pr": 323,
        "head_sha": BASE_AGENT4_HEAD,
        "dedicated_run": 35308426745,
        "standard_run": 35308426643,
        "artifact_id": 10531804129,
        "artifact_digest": "sha256:567e5dc407c6f9ee521530814b8cf88b5ded1919a94fc13533f3c01300fb17f8",
        "signed_log_heat_target_independently_preflighted": True,
        "consumed_in_executable_ancestry": True,
    },
    "previous_agent5": {
        "pr": 315,
        "head_sha": "3d20508865b91de21f11a069f7c6cbc5b26bba12",
        "terminal_tail_schedule_independently_preflighted": True,
        "consumed_in_executable_ancestry": False,
    },
}
STATES = {
    "terminal_tail_schedule_independently_preflighted": True,
    "signed_log_heat_target_ready": True,
    "signed_log_heat_target_independently_preflighted": True,
    "precision_preserving_three_bump_representation_ready": False,
    "heat_compensation_completed": False,
    "global_leading_profile_reconstructed": False,
    "source_real_conjugate_pair_kernel_upstream_ready": True,
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


def bind_independent_heat_audit() -> dict[str, Any]:
    report = run_audit()
    if report["structural_preflight_passed"] is not True or not all(report["guards"].values()):
        raise RuntimeError("independent signed-log heat-target preflight regressed")
    if report["formal_momentum_gate"] != 1.0e-3 or report["formal_divergence_gate"] != 1.0e-5:
        raise ValueError("formal gates changed")
    if report["formal_full_domain_pde_gate_assessed"] is not False or report["pde_validated"] is not False:
        raise ValueError("local heat audit cannot promote the full-domain PDE gate")
    if report["st006_context"]["directly_comparable_to_this_local_audit"] is not False:
        raise ValueError("local heat audit cannot be relabeled ST006-comparable")
    cases = {case["name"]: case for case in report["cases"]}
    if set(cases) != {"default_extreme", "moderate_materializable", "perturbed_outer"}:
        raise RuntimeError("unexpected signed-log audit matrix")
    default = cases["default_extreme"]["finest_eta_0p2"]
    if default["sign"] != [1, -1, 1]:
        raise RuntimeError("default signed-log target signs changed")
    if default["underflow_channels"] != ["C_p", "S"]:
        raise RuntimeError("default underflow hierarchy changed")
    if float(default["span_decades"]) < 400.0:
        raise RuntimeError("default signed-log channel span unexpectedly collapsed")
    return {
        "schema": report["schema"],
        "metrics": report["metrics"],
        "guards": report["guards"],
        "structural_preflight_passed": True,
        "default_eta_0p2": default,
        "formal_full_domain_pde_gate_assessed": False,
        "st006_directly_comparable": False,
    }


def build_checkpoint(audit: dict[str, Any]) -> dict[str, Any]:
    target = KokunoLogRescaledHeatRepairTarget()
    encoded = target.signed_log_target(0.2)
    signs = np.asarray(encoded["sign"], dtype=int).reshape(3)
    logs = np.asarray(encoded["log_abs"], dtype=float).reshape(3)
    log10_abs = logs / math.log(10.0)
    precision = target.precision_report(0.2)
    if signs.tolist() != [1, -1, 1]:
        raise RuntimeError("Agent-1 default signed-log target signs changed")
    if list(precision["underflow_channels"]) != ["C_p", "S"]:
        raise RuntimeError("Agent-1 default precision barrier changed")
    try:
        target.materialize_target(0.2)
    except OverflowError:
        pass
    else:
        raise RuntimeError("default target unexpectedly became float64 materializable")

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "base_agent4_head": BASE_AGENT4_HEAD,
        "upstream": UPSTREAM,
        "fixed_gates": FIXED_GATES,
        "st006_reference": ST006_REFERENCE,
        "typed_heat_target_component": {
            "schema": target.to_payload()["schema"],
            "sha256": target.sha256,
            "serialization_api": "KokunoLogRescaledHeatRepairTarget.save_json/load_json",
            "signed_log_api": "KokunoLogRescaledHeatRepairTarget.signed_log_target",
            "rescaled_api": "KokunoLogRescaledHeatRepairTarget.rescaled_target",
            "precision_api": "KokunoLogRescaledHeatRepairTarget.precision_report",
            "float64_materialization_api": "KokunoLogRescaledHeatRepairTarget.materialize_target",
            "velocity_api_ready": False,
            "pressure_api_ready": False,
            "forcing_api_ready": False,
            "matlab_export_ready": False,
        },
        "default_eta_0p2_precision_frontier": {
            "sign": signs.tolist(),
            "log_abs": logs.tolist(),
            "log10_abs": log10_abs.tolist(),
            "underflow_channels": list(precision["underflow_channels"]),
            "channel_span_decades": float(np.max(log10_abs) - np.min(log10_abs)),
            "fully_float64_materializable": False,
            "existing_float64_three_bump_inverse_faithful_for_default": False,
        },
        "independent_signed_log_heat_audit": audit,
        "baseline_vs_kokuno": [
            {
                "name": "ST006",
                "scope": "held-out full-domain normalized momentum",
                "momentum_sampled_max": ST006_REFERENCE["momentum_sampled_max"],
                "volume_l2": ST006_REFERENCE["volume_l2"],
                "pde_validated": False,
            },
            {
                "name": "Kokuno signed-log heat target",
                "scope": "local source heat-target algebra / precision preflight",
                "momentum_sampled_max": None,
                "volume_l2": None,
                "pde_validated": False,
                "directly_comparable_to_ST006": False,
                "reason": "no globally matched Kokuno velocity/pressure/forcing candidate exists",
            },
        ],
        "historical_stage_context": {
            "scope": "older core-local zero-force diagnostics; not ST006-comparable",
            "leading_only_raw_sample_rms": 10.4529319480,
            "leading_plus_oscillatory_raw_sample_rms": 10.3067008431,
            "after_signed_correction_raw_sample_rms": 10.3075981486,
            "oscillatory_over_leading_rms_ratio": 0.98601051785,
            "signed_corrected_over_oscillatory_rms_ratio": 1.00008706040,
            "source_pr": 264,
            "directly_comparable_to_ST006": False,
        },
        "states": STATES,
        "routing": {
            "repeat_terminal_scale_derivation": False,
            "repeat_signed_log_target_derivation": False,
            "feed_underflowed_zero_channels_to_float64_three_bump_inverse": False,
            "rerun_duplicate_or_one_column_correction_cycle": False,
            "run_formal_full_domain_gate_now": False,
            "replace_ST006_as_repository_pde_baseline": False,
            "priority_blocker": (
                "Agent 1: implement a precision-preserving/hierarchical three-bump coefficient "
                "representation, apply the independently preflighted signed-log target to I2, "
                "then emit the first global leading velocity/pressure contract"
            ),
            "agent2_blocker": (
                "bind actual source pulse/background data into the real conjugate-pair kernel and "
                "materialize a direct public velocity(x,y,z,t) correction; a genuinely independent "
                "second physical pair/column is still required for Agent 3"
            ),
            "agent3_blocker": (
                "the velocity-column adapter is ready, but the duplicate column is rank deficient; "
                "screen the next genuinely independent Agent-2 physical column before any finite cycle"
            ),
            "agent4_next": (
                "independently audit the hierarchical heat repair when materialized; return to the "
                "fixed held-out full-domain PDE gate as soon as a global composite exists"
            ),
        },
        "truth_boundary": {
            "source_version_and_provenance_preserved": True,
            "agent5_new_kokuno_mathematics_added": False,
            "sibling_evidence_laundered_into_executable_ancestry": False,
            "local_heat_target_promoted_to_global_candidate": False,
            "free_residual_defined_forcing_used": False,
            "pressure_or_forcing_refit_in_agent5": False,
            "threshold_relaxed": False,
            "st006_used_as_kokuno_source_truth": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
        "next_shortest_closure": [
            "Agent 1: hierarchical/high-precision three-bump repair -> apply to I2 -> global leading velocity/pressure",
            "Agent 2: actual source pulse/background -> public Q-scaled real-pair velocity -> genuinely independent second column",
            "Agent 3: screen second physical column -> bounded two-column finite correction cycle only if rank/budget guard passes",
            "Agent 4: independent repair audit, then held-out <=1e-3 momentum and <=1e-5 divergence on the global composite",
            "Agent 5: serialize the complete candidate and add minimal Python/MATLAB export smoke only after the global interface exists",
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
    typed = unsigned.get("typed_heat_target_component", {})
    if typed.get("velocity_api_ready") is not False or typed.get("pressure_api_ready") is not False:
        raise ValueError("local heat target was promoted to a field API")
    rows = unsigned.get("baseline_vs_kokuno", [{}, {}])
    if len(rows) != 2 or rows[1].get("directly_comparable_to_ST006") is not False:
        raise ValueError("local heat target was relabeled ST006-comparable")
    routing = unsigned.get("routing", {})
    if routing.get("run_formal_full_domain_gate_now") is not False:
        raise ValueError("formal PDE gate routed before a global candidate exists")
    if routing.get("rerun_duplicate_or_one_column_correction_cycle") is not False:
        raise ValueError("rank-deficient correction path was rerouted")
    truth = unsigned.get("truth_boundary", {})
    if truth.get("threshold_relaxed") is not False or truth.get("free_residual_defined_forcing_used") is not False:
        raise ValueError("scientific contract weakened")
    if any(unsigned["states"][key] for key in (
        "leading_ready", "correction_ready", "velocity_export_ready",
        "formal_full_domain_pde_gate_assessed", "pde_validated",
    )):
        raise ValueError("premature scientific promotion")
    return payload


def write_bundle(output_dir: str | Path) -> dict[str, Any]:
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    target = KokunoLogRescaledHeatRepairTarget()
    target_path = output / "signed_log_heat_target.json"
    target.save_json(target_path)
    replayed = KokunoLogRescaledHeatRepairTarget.load_json(target_path)
    if replayed.to_payload() != target.to_payload():
        raise RuntimeError("signed-log heat target round-trip changed")
    audit = bind_independent_heat_audit()
    checkpoint = build_checkpoint(audit)
    validate_checkpoint(checkpoint)
    (output / "independent_signed_log_heat_audit_summary.json").write_text(
        json.dumps(audit, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    path = output / "signed_log_heat_routing_checkpoint.json"
    path.write_text(
        json.dumps(checkpoint, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    if load_checkpoint(path) != checkpoint:
        raise RuntimeError("checkpoint round-trip changed")
    return checkpoint


def load_checkpoint(path: str | Path) -> dict[str, Any]:
    return validate_checkpoint(json.loads(Path(path).read_text(encoding="utf-8")))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("artifacts/kokuno_agent5/signed_log_heat_routing_checkpoint_v13"),
    )
    args = parser.parse_args()
    checkpoint = write_bundle(args.output_dir)
    print(json.dumps({
        "checkpoint_sha256": checkpoint["checkpoint_sha256"],
        "states": checkpoint["states"],
        "precision_frontier": checkpoint["default_eta_0p2_precision_frontier"],
        "independent_metrics": checkpoint["independent_signed_log_heat_audit"]["metrics"],
    }, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
