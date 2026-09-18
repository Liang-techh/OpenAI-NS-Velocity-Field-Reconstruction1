"""Kokuno Agent-5 integration checkpoint after the terminal-tail audit.

This module adds no new Kokuno mathematics. It binds Agent 1's executable
terminal release/tail schedule to Agent 4's independent schedule audit and the
already independently checked corrected heat-V2 machinery. Agent 2/3 newest
work is recorded as sibling provenance only.

The routing consequence is narrow but important: source-scale discovery through
``X_star/e_star`` and ``X_tail/c_inf/rho_o`` is no longer the leading blocker.
The default terminal schedule is intentionally kept in log space because its
absolute ``X_tail`` exceeds IEEE-754 range. The next Agent-1 step is therefore
to apply the corrected V2 discrepancy / three-bump repair in log-normalized
variables and materialize a globally matched leading velocity/pressure.

This is a local source-schedule integration receipt, not a Navier--Stokes gate.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

from .kokuno_heat_v2_routing_checkpoint import run_bound_heat_v2_audit
from .kokuno_independent_terminal_tail_audit import run_audit as run_terminal_tail_audit
from .kokuno_terminal_tail_schedule import KokunoTerminalTailSchedule

SCHEMA = "kokuno-agent5-terminal-tail-routing-checkpoint-v12"
TASK_ID = "KOKUNO-A5-TERMINAL-TAIL-ROUTING-CHECKPOINT-012"
BASE_AGENT4_HEAD = "f5a3697c5a0fa03000785856feab90b5a491890b"

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
        "pr": 310,
        "head_sha": "3ffff01961840b4fcecf1fe278af5b1cc7301789",
        "standard_run": 35303791325,
        "terminal_tail_schedule_ready": True,
        "terminal_tail_schedule_bound_to_outer_schedule": True,
        "terminal_scales_ready_in_log_space": True,
        "default_heat_inputs_materializable_as_float": False,
        "corrected_heat_discrepancy_applied_to_patch": False,
        "heat_compensation_completed": False,
        "global_leading_profile_reconstructed": False,
        "consumed_in_executable_ancestry": True,
    },
    "agent2": {
        "pr": 312,
        "head_sha": "09eb62601bce831fd1490b1dbd18b97c845d8358",
        "source_normalized_Dz_pulse_sensitivity_ready": True,
        "actual_source_pulse_background_instantiated": False,
        "second_public_covariance_velocity_column_realized": False,
        "consumed_in_executable_ancestry": False,
    },
    "agent3": {
        "pr": 313,
        "head_sha": "ca274e9f41f7b1cac669e78c6007fc0f33cfde78",
        "dedicated_run": 35304061342,
        "standard_run": 35304061303,
        "public_covariance_tangent_bridge_ready": True,
        "amplitude_tangent_is_independent_second_direction": False,
        "active_nodes_requiring_second_direction": 25,
        "rank_two_nodes": 0,
        "missing_relative_vector_rms": 0.21019450212274277,
        "finite_cycle_rerun_allowed": False,
        "consumed_in_executable_ancestry": False,
    },
    "agent4": {
        "pr": 314,
        "head_sha": BASE_AGENT4_HEAD,
        "dedicated_run": 35304724349,
        "standard_run": 35304724348,
        "artifact_id": 10530773583,
        "artifact_digest": "sha256:e0febd271cd06b240746422e1624dd561c227efdddfa5c0c2898caaa2c77c776",
        "terminal_tail_schedule_independently_preflighted": True,
        "consumed_in_executable_ancestry": True,
    },
    "previous_agent5_outer_patch": {
        "pr": 305,
        "head_sha": "5bab222ff722f9e86fb254ae7fc2576dfac50687",
        "outer_reserved_patch_independently_validated": True,
        "corrected_heat_v2_routing_ready": True,
        "consumed_in_executable_ancestry": False,
        "reason": "sibling Agent-5 receipt; its mathematical prerequisites are present but the receipt module is not in this ancestry",
    },
}

STATES = {
    "corrected_heat_discrepancy_independently_validated": True,
    "reserved_patch_schedule_ready": True,
    "X_star_e_star_binding_ready": True,
    "terminal_tail_schedule_bound_to_outer_schedule": True,
    "terminal_tail_schedule_independently_preflighted": True,
    "terminal_scales_ready_in_log_space": True,
    "default_terminal_heat_inputs_float_materializable": False,
    "physical_three_bump_target_ready": False,
    "heat_compensation_completed": False,
    "global_leading_profile_reconstructed": False,
    "source_normalized_Dz_sensitivity_upstream_ready": True,
    "actual_source_pulse_background_ready": False,
    "public_covariance_tangent_bridge_upstream_ready": True,
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


def bind_terminal_audit() -> dict[str, Any]:
    report = run_terminal_tail_audit()
    summary = report.get("summary", {})
    if summary.get("structural_preflight_passed") is not True:
        raise RuntimeError("terminal-tail independent preflight regressed")
    if not all(summary.get("guard_results", {}).values()):
        raise RuntimeError("terminal-tail independent guard matrix regressed")
    truth = report.get("truth_boundary", {})
    if truth.get("terminal_tail_schedule_independently_preflighted") is not True:
        raise RuntimeError("terminal-tail source schedule was not independently preflighted")
    if truth.get("formal_full_domain_pde_gate_assessed") is not False:
        raise ValueError("local terminal audit cannot become the formal PDE gate")
    if truth.get("pde_validated") is not False:
        raise ValueError("local terminal audit cannot promote PDE validation")
    formal = report.get("formal_gates_unchanged", {})
    if formal.get("normalized_momentum") != 1.0e-3 or formal.get("divergence_max") != 1.0e-5:
        raise ValueError("formal project gates changed")
    baseline = report.get("cross_route_baseline", {})
    if baseline.get("directly_comparable") is not False:
        raise ValueError("local terminal audit cannot be relabeled ST006-comparable")

    cases = report.get("cases", [])
    if {case.get("name") for case in cases} != {
        "default_extreme",
        "moderate_materializable",
        "perturbed_outer",
    }:
        raise RuntimeError("unexpected terminal-tail audit case matrix")
    case_map = {case["name"]: case for case in cases}
    if case_map["default_extreme"].get("public_materialized") is not False:
        raise RuntimeError("default extreme schedule unexpectedly materialized as floats")
    if case_map["moderate_materializable"].get("public_materialized") is not True:
        raise RuntimeError("moderate terminal schedule stopped materializing")

    return {
        "summary": summary,
        "default_extreme_public_materialized": False,
        "moderate_public_materialized": True,
        "parameter_case_count": len(cases),
        "formal_full_domain_pde_gate_assessed": False,
        "st006_directly_comparable": False,
    }


def build_checkpoint(
    schedule: KokunoTerminalTailSchedule,
    terminal_audit: dict[str, Any],
    heat_audit: dict[str, Any],
) -> dict[str, Any]:
    scales = schedule.log_scale_report()
    ratio = schedule.log_ratio_report()
    if scales["X_tail"] is not None or scales["c_inf"] is not None:
        raise RuntimeError("default extreme schedule must remain log-space only")

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "base_agent4_head": BASE_AGENT4_HEAD,
        "upstream": UPSTREAM,
        "fixed_gates": FIXED_GATES,
        "st006_reference": ST006_REFERENCE,
        "typed_schedule_component": {
            "schema": schedule.to_payload()["schema"],
            "sha256": schedule.sha256,
            "serialization_api": "KokunoTerminalTailSchedule.save_json/load_json",
            "log_scale_api": "KokunoTerminalTailSchedule.log_scale_report",
            "heat_materialization_api": "KokunoTerminalTailSchedule.materialize_heat_inputs",
            "default_heat_materialization_ready": False,
            "full_domain_velocity_api_ready": False,
            "pressure_api_ready": False,
            "forcing_api_ready": False,
            "matlab_export_ready": False,
        },
        "default_terminal_scales": {
            "Q_release": float(scales["Q_release"]),
            "Q_in": float(scales["Q_in"]),
            "Q_p": float(scales["Q_p"]),
            "release_hold_length": float(scales["release_hold_length"]),
            "rho_o": float(scales["rho_o"]),
            "log_X_tail": float(scales["log_X_tail"]),
            "log_c_inf": float(scales["log_c_inf"]),
            "log_X_K": float(scales["log_X_K"]),
            "log_X_b": float(scales["log_X_b"]),
            "log_X_star_over_X_K": float(ratio["log_X_star_over_X_K"]),
            "log_e_star_over_e_K": float(ratio["log_e_star_over_e_K"]),
            "X_tail": None,
            "c_inf": None,
            "source_hidden_numeric_choices_recovered": False,
            "autonomous_T_f_and_c_o": True,
        },
        "independent_terminal_tail_audit": terminal_audit,
        "independent_corrected_heat_v2_audit": heat_audit,
        "baseline_vs_kokuno": [
            {
                "name": "ST006",
                "scope": "held-out full-domain normalized momentum",
                "momentum_sampled_max": ST006_REFERENCE["momentum_sampled_max"],
                "volume_l2": ST006_REFERENCE["volume_l2"],
                "pde_validated": False,
            },
            {
                "name": "Kokuno terminal-tail schedule",
                "scope": "local source schedule / scale-chain implementation preflight",
                "momentum_sampled_max": None,
                "volume_l2": None,
                "pde_validated": False,
                "directly_comparable_to_ST006": False,
                "reason": "no globally matched Kokuno velocity/pressure/forcing candidate exists",
            },
        ],
        "states": STATES,
        "routing": {
            "repeat_X_star_e_star_search": False,
            "repeat_terminal_tail_scale_derivation": False,
            "repeat_terminal_tail_independent_audit": False,
            "materialize_default_absolute_tail_scales_in_float": False,
            "rerun_one_column_correction_cycle": False,
            "run_formal_full_domain_gate_now": False,
            "replace_ST006_as_repository_pde_baseline": False,
            "priority_blocker": (
                "Agent 1 should carry the corrected V2 heat discrepancy and existing three-bump repair "
                "in log-normalized/rescaled variables across the extreme terminal scales, apply it to "
                "the I2/heat matching construction, and emit the first globally matched leading velocity/pressure"
            ),
            "agent2_blocker": (
                "turn the source-normalized Dz sensitivity contract into an actual source pulse/background "
                "perturbation and materialize a genuinely independent second public [u,v,w] covariance column"
            ),
            "agent3_blocker": (
                "the public covariance-tangent bridge is ready, but the tested amplitude tangent is rank-degenerate; "
                "wait for Agent 2's physical second direction before rerunning the finite correction cycle"
            ),
        },
        "truth_boundary": {
            "source_version_and_provenance_preserved": True,
            "agent5_new_kokuno_mathematics_added": False,
            "sibling_evidence_laundered_into_executable_ancestry": False,
            "local_schedule_promoted_to_global_candidate": False,
            "free_residual_defined_forcing_used": False,
            "pressure_or_forcing_refit_in_agent5": False,
            "threshold_relaxed": False,
            "st006_used_as_kokuno_source_truth": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
        "next_shortest_closure": [
            "Agent 1: apply corrected V2 + three-bump heat repair in log-normalized variables and emit globally matched leading velocity/pressure",
            "Agent 2: materialize an actual source pulse/background perturbation and genuinely independent second public covariance [u,v,w] column",
            "Agent 3: evaluate that physical second column through the public covariance-tangent bridge, then run a frozen two-column finite correction cycle only if rank/budget guards pass",
            "Agent 4: after the global composite exists, run independent held-out normalized momentum <=1e-3 and divergence <=1e-5 plus same-protocol ST006 comparison",
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

    typed = unsigned.get("typed_schedule_component", {})
    if typed.get("full_domain_velocity_api_ready") is not False:
        raise ValueError("terminal schedule was promoted to a global velocity API")
    if typed.get("pressure_api_ready") is not False or typed.get("forcing_api_ready") is not False:
        raise ValueError("missing pressure/forcing APIs were relabeled ready")
    if unsigned.get("baseline_vs_kokuno", [{}, {}])[1].get("directly_comparable_to_ST006") is not False:
        raise ValueError("local terminal schedule was relabeled ST006-comparable")
    if unsigned.get("routing", {}).get("run_formal_full_domain_gate_now") is not False:
        raise ValueError("formal PDE gate routed before a global candidate exists")
    if unsigned.get("routing", {}).get("rerun_one_column_correction_cycle") is not False:
        raise ValueError("rank-degenerate one-column correction was rerouted")
    truth = unsigned.get("truth_boundary", {})
    if truth.get("threshold_relaxed") is not False:
        raise ValueError("validation threshold was relaxed")
    if truth.get("free_residual_defined_forcing_used") is not False:
        raise ValueError("free residual-defined forcing was introduced")
    if any(unsigned["states"][key] for key in (
        "leading_ready",
        "correction_ready",
        "velocity_export_ready",
        "formal_full_domain_pde_gate_assessed",
        "pde_validated",
    )):
        raise ValueError("premature scientific promotion")
    return payload


def write_bundle(output_dir: str | Path) -> dict[str, Any]:
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)

    schedule = KokunoTerminalTailSchedule()
    schedule_path = output / "terminal_tail_schedule.json"
    schedule.save_json(schedule_path)
    replayed = KokunoTerminalTailSchedule.load_json(schedule_path)
    if replayed.to_payload() != schedule.to_payload():
        raise RuntimeError("terminal-tail schedule round-trip changed")

    terminal_audit = bind_terminal_audit()
    heat_audit = run_bound_heat_v2_audit()
    checkpoint = build_checkpoint(replayed, terminal_audit, heat_audit)
    validate_checkpoint(checkpoint)

    (output / "terminal_tail_independent_summary.json").write_text(
        json.dumps(terminal_audit, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    (output / "corrected_heat_v2_independent_summary.json").write_text(
        json.dumps(heat_audit, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    checkpoint_path = output / "terminal_tail_routing_checkpoint.json"
    checkpoint_path.write_text(
        json.dumps(checkpoint, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    if load_checkpoint(checkpoint_path) != checkpoint:
        raise RuntimeError("terminal-tail routing checkpoint round-trip changed")
    return checkpoint


def load_checkpoint(path: str | Path) -> dict[str, Any]:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    return validate_checkpoint(payload)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("artifacts/kokuno_agent5/terminal_tail_routing_checkpoint_v12"),
    )
    args = parser.parse_args()
    checkpoint = write_bundle(args.output_dir)
    print(json.dumps({
        "checkpoint_sha256": checkpoint["checkpoint_sha256"],
        "states": checkpoint["states"],
        "default_terminal_scales": checkpoint["default_terminal_scales"],
        "terminal_audit_summary": checkpoint["independent_terminal_tail_audit"]["summary"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
