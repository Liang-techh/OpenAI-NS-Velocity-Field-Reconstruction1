"""Kokuno Agent-5 fail-closed routing checkpoint for the terminal-tail stage.

No new Kokuno mathematics is introduced here. The checkpoint binds Agent 1's
terminal release/tail schedule, Agent 4's independent audit, and the inherited
corrected heat-V2 audit. It records the newest Agent 2/3 work as sibling
provenance and keeps the formal full-domain PDE gate unassessed.
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
        "terminal_scales_ready_in_log_space": True,
        "global_leading_profile_reconstructed": False,
        "consumed_in_executable_ancestry": True,
    },
    "agent2": {
        "pr": 312,
        "head_sha": "09eb62601bce831fd1490b1dbd18b97c845d8358",
        "dedicated_run": 35304744001,
        "standard_run": 35304744009,
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
    "previous_agent5": {
        "pr": 305,
        "head_sha": "5bab222ff722f9e86fb254ae7fc2576dfac50687",
        "outer_reserved_patch_independently_validated": True,
        "corrected_heat_v2_routing_ready": True,
        "consumed_in_executable_ancestry": False,
    },
}
STATES = {
    "corrected_heat_discrepancy_independently_validated": True,
    "X_star_e_star_binding_ready": True,
    "terminal_tail_schedule_bound_to_outer_schedule": True,
    "terminal_tail_schedule_independently_preflighted": True,
    "terminal_scales_ready_in_log_space": True,
    "default_terminal_heat_inputs_float_materializable": False,
    "physical_three_bump_target_ready": False,
    "heat_compensation_completed": False,
    "global_leading_profile_reconstructed": False,
    "source_normalized_Dz_sensitivity_upstream_ready": True,
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
    summary = report["summary"]
    if summary["structural_preflight_passed"] is not True:
        raise RuntimeError("terminal-tail independent preflight regressed")
    if not all(summary["guard_results"].values()):
        raise RuntimeError("terminal-tail guard matrix regressed")
    truth = report["truth_boundary"]
    if truth["terminal_tail_schedule_independently_preflighted"] is not True:
        raise RuntimeError("terminal-tail schedule not independently preflighted")
    if truth["formal_full_domain_pde_gate_assessed"] is not False or truth["pde_validated"] is not False:
        raise ValueError("local audit cannot promote the full-domain PDE gate")
    if report["formal_gates_unchanged"] != {"normalized_momentum": 1.0e-3, "divergence_max": 1.0e-5}:
        raise ValueError("formal gates changed")
    if report["cross_route_baseline"]["directly_comparable"] is not False:
        raise ValueError("local audit cannot be relabeled ST006-comparable")
    cases = {case["name"]: case for case in report["cases"]}
    if set(cases) != {"default_extreme", "moderate_materializable", "perturbed_outer"}:
        raise RuntimeError("unexpected terminal-tail audit matrix")
    if cases["default_extreme"]["public_materialized"] is not False:
        raise RuntimeError("extreme schedule unexpectedly materialized as complete heat inputs")
    if cases["moderate_materializable"]["public_materialized"] is not True:
        raise RuntimeError("moderate schedule stopped materializing")
    return {
        "summary": summary,
        "parameter_case_count": len(cases),
        "default_extreme_public_materialized": False,
        "moderate_public_materialized": True,
        "formal_full_domain_pde_gate_assessed": False,
        "st006_directly_comparable": False,
    }


def build_checkpoint(
    schedule: KokunoTerminalTailSchedule,
    terminal_audit: dict[str, Any],
    heat_audit: dict[str, Any],
) -> dict[str, Any]:
    scales = schedule.log_scale_report()
    ratios = schedule.log_ratio_report()
    # The default X_tail is intentionally beyond float range, while c_inf itself
    # remains finite. materialize_heat_inputs() still fails closed because the
    # complete absolute input tuple cannot be represented safely.
    if scales["X_tail"] is not None:
        raise RuntimeError("default X_tail unexpectedly fits in float")
    if scales["c_inf"] is None or not scales["c_inf"] > 0.0:
        raise RuntimeError("default c_inf should remain a finite positive scalar")
    try:
        schedule.materialize_heat_inputs()
    except OverflowError:
        pass
    else:
        raise RuntimeError("default complete heat input tuple unexpectedly materialized")

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
            "default_complete_heat_materialization_ready": False,
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
            "log_X_star_over_X_K": float(ratios["log_X_star_over_X_K"]),
            "log_e_star_over_e_K": float(ratios["log_e_star_over_e_K"]),
            "X_tail": None,
            "c_inf": float(scales["c_inf"]),
            "complete_heat_inputs_materializable": False,
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
                "scope": "local source schedule / scale-chain preflight",
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
            "materialize_default_absolute_X_tail_in_float": False,
            "rerun_one_column_correction_cycle": False,
            "run_formal_full_domain_gate_now": False,
            "replace_ST006_as_repository_pde_baseline": False,
            "priority_blocker": (
                "Agent 1 should carry corrected V2 and the existing three-bump repair in "
                "log-normalized/rescaled variables across the extreme tail scales, apply the "
                "repair to the I2/heat match, and emit global leading velocity/pressure"
            ),
            "agent2_blocker": (
                "materialize the actual source pulse/background perturbation and a genuinely "
                "independent second public [u,v,w] covariance column"
            ),
            "agent3_blocker": (
                "the tangent bridge is ready but the amplitude tangent is rank-degenerate; "
                "wait for Agent 2's physical second direction before the finite cycle"
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
            "Agent 1: apply corrected V2 + three-bump repair in log-normalized variables; emit global leading velocity/pressure",
            "Agent 2: materialize actual source pulse/background perturbation plus an independent second covariance [u,v,w] column",
            "Agent 3: pass that physical second column through the tangent/rank guards, then run a frozen two-column finite correction cycle",
            "Agent 4: once global composite exists, run independent held-out <=1e-3 momentum and <=1e-5 divergence plus same-protocol ST006 comparison",
            "Agent 5: serialize complete velocity/pressure/forcing candidate and run Python/MATLAB export smoke only after the global interface exists",
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
        raise ValueError("fixed gate or ST006 reference changed")
    if unsigned.get("upstream") != UPSTREAM:
        raise ValueError("upstream provenance changed")
    if unsigned.get("states") != STATES:
        raise ValueError("scientific state vector changed")
    typed = unsigned.get("typed_schedule_component", {})
    if typed.get("full_domain_velocity_api_ready") is not False:
        raise ValueError("terminal schedule was promoted to a global velocity API")
    if typed.get("pressure_api_ready") is not False or typed.get("forcing_api_ready") is not False:
        raise ValueError("missing pressure/forcing APIs were relabeled ready")
    rows = unsigned.get("baseline_vs_kokuno", [{}, {}])
    if len(rows) != 2 or rows[1].get("directly_comparable_to_ST006") is not False:
        raise ValueError("local terminal schedule was relabeled ST006-comparable")
    routing = unsigned.get("routing", {})
    if routing.get("run_formal_full_domain_gate_now") is not False:
        raise ValueError("formal PDE gate routed before a global candidate exists")
    if routing.get("rerun_one_column_correction_cycle") is not False:
        raise ValueError("rank-degenerate one-column correction was rerouted")
    truth = unsigned.get("truth_boundary", {})
    if truth.get("threshold_relaxed") is not False:
        raise ValueError("validation threshold was relaxed")
    if truth.get("free_residual_defined_forcing_used") is not False:
        raise ValueError("free residual-defined forcing was introduced")
    if any(unsigned["states"][key] for key in (
        "leading_ready", "correction_ready", "velocity_export_ready",
        "formal_full_domain_pde_gate_assessed", "pde_validated",
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
        raise RuntimeError("terminal schedule round-trip changed")
    terminal_audit = bind_terminal_audit()
    heat_audit = run_bound_heat_v2_audit()
    checkpoint = build_checkpoint(replayed, terminal_audit, heat_audit)
    validate_checkpoint(checkpoint)
    (output / "terminal_tail_independent_summary.json").write_text(
        json.dumps(terminal_audit, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (output / "corrected_heat_v2_independent_summary.json").write_text(
        json.dumps(heat_audit, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    path = output / "terminal_tail_routing_checkpoint.json"
    path.write_text(json.dumps(checkpoint, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if load_checkpoint(path) != checkpoint:
        raise RuntimeError("checkpoint round-trip changed")
    return checkpoint


def load_checkpoint(path: str | Path) -> dict[str, Any]:
    return validate_checkpoint(json.loads(Path(path).read_text(encoding="utf-8")))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir", type=Path,
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
