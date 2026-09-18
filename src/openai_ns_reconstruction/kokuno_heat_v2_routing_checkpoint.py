"""Kokuno Agent-5 routing checkpoint after the corrected heat-tail audit.

This module adds no new Kokuno mathematics.  It consumes Agent 1's corrected
heat-replacement discrepancy (PR #289) through Agent 4's independent robustness
audit (PR #294), binds the newest Agent 2 / Agent 3 sibling state, and emits a
deterministic fail-closed routing receipt.

The important routing change relative to Agent-5 PR #284 is narrow: the local
first-order ``Delta C_p`` terminal-tail factor blocker is now independently
resolved.  That does *not* make a global leading field ready.  The source
reserved-patch scales (notably ``X_star/e_star``) and remaining outer schedule
still have to be instantiated before the corrected discrepancy can become the
physical target for the three-bump repair/global core-to-heat assembly.

This local source-moment audit is not directly comparable with the retained
ST006 full-domain momentum benchmark and cannot set ``pde_validated=true``.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_heat_discrepancy_v2_robustness import run_audit

SCHEMA = "kokuno-agent5-heat-v2-routing-checkpoint-v10"
TASK_ID = "KOKUNO-A5-HEAT-V2-ROUTING-CHECKPOINT-010"
BASE_AGENT4_HEAD = "365a70b911d45b7b3e105abfce758b4a95d9efdb"

AGENT1_UPSTREAM = {
    "pr": 289,
    "head_sha": "eb736740373792f51b261549e958420d4b37b843",
    "standard_run": 35295107876,
    "corrected_heat_discrepancy_v2_ready": True,
    "delta_cp_terminal_tail_factor_fixed": True,
    "earlier_reserved_patch_scales_recovered": False,
    "global_leading_profile_reconstructed": False,
    "consumed_in_executable_ancestry": True,
}
AGENT2_UPSTREAM = {
    "pr": 291,
    "head_sha": "2655d90d9fa561eb49e982245ed550cf635271df",
    "standard_run": 35295932605,
    "projected_pulse_inverse_contract_ready": True,
    "source_actual_pulse_forcing_instantiated": False,
    "source_actual_background_path_instantiated": False,
    "public_velocity_correction_materialized": False,
    "consumed_in_executable_ancestry": False,
}
AGENT3_UPSTREAM = {
    "pr": 293,
    "head_sha": "78c89b42f4ea5ad45884a292040339088111f9cf",
    "dedicated_run": 35296279431,
    "standard_run": 35296279396,
    "missing_second_column_target_measured": True,
    "one_column_cycle_accepted": False,
    "held_out_total_defect_ratio": 1.1506888629778744,
    "held_out_theta_defect_ratio": 1.0876983964807052,
    "second_public_covariance_column_realized": False,
    "consumed_in_executable_ancestry": False,
}
AGENT4_UPSTREAM = {
    "pr": 294,
    "head_sha": BASE_AGENT4_HEAD,
    "dedicated_run": 35296597060,
    "standard_run": 35296619078,
    "corrected_heat_discrepancy_independently_robust": True,
    "consumed_in_executable_ancestry": True,
}
PREVIOUS_AGENT5 = {
    "pr": 284,
    "head_sha": "e6eba01735ca3d378d36f56bff110a1c12fc4578",
    "previous_tail_factor_blocker_active": True,
    "superseded_by_agent1_289_and_agent4_294": True,
    "consumed_in_executable_ancestry": False,
}
COORDINATION = {
    "pr": 290,
    "head_sha": "4cedee7a1aec267d994e00747aa3e9853c7e3ab0",
    "st006_cross_route_rule_ready": True,
    "consumed_in_executable_ancestry": False,
}

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
    "directly_comparable_to_this_heat_moment_audit": False,
    "comparability_reason": (
        "this checkpoint validates a local source heat-moment map; ST006 is a "
        "frozen full-domain normalized momentum protocol"
    ),
}

STATES = {
    "delta_cp_tail_factor_blocker_resolved": True,
    "corrected_heat_discrepancy_formula_ready": True,
    "corrected_heat_discrepancy_independently_validated": True,
    "local_three_bump_repair_operator_ready": True,
    "physical_three_bump_target_ready": False,
    "heat_compensation_completed": False,
    "global_leading_profile_reconstructed": False,
    "source_complete_curl_operator_upstream_ready": True,
    "source_actual_pulse_inverse_ready": False,
    "source_oscillatory_public_velocity_ready": False,
    "missing_second_covariance_column_target_upstream_ready": True,
    "second_covariance_column_realized": False,
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


def run_bound_heat_v2_audit() -> dict[str, Any]:
    """Rerun Agent 4's independent robustness audit and bind only routed facts."""

    report = run_audit()
    if report.get("structural_preflight_passed") is not True:
        raise RuntimeError("corrected heat discrepancy no longer passes Agent-4 robustness preflight")
    if report.get("st006_directly_comparable") is not False:
        raise ValueError("local heat audit cannot be relabeled as directly comparable to ST006")
    if report.get("formal_full_domain_pde_gate_assessed") is not False:
        raise ValueError("local heat audit cannot become the formal PDE gate")
    if report.get("pde_validated") is not False:
        raise ValueError("local heat audit cannot promote PDE validity")

    summary = report.get("summary", {})
    component = np.asarray(summary.get("finest_max_relative_by_component", ()), dtype=float)
    if component.shape != (3,) or not np.all(np.isfinite(component)):
        raise ValueError("unexpected Agent-4 heat-v2 robustness summary")
    if not (component[0] < 2.0e-5 and component[1] < 2.0e-6 and component[2] < 2.0e-5):
        raise RuntimeError("corrected heat-v2 independent component agreement regressed")
    if not float(summary.get("medium_to_fine_max_relative", np.inf)) < 5.0e-5:
        raise RuntimeError("corrected heat-v2 production quadrature did not stabilize")
    if not float(summary.get("tail_fix_vs_independent_half_first_order_max_relative", np.inf)) < 5.0e-5:
        raise RuntimeError("Delta C_p restored-copy attribution no longer closes")
    if float(summary.get("unchanged_delta_s_i_max_abs", np.inf)) != 0.0:
        raise RuntimeError("heat-v2 fix leaked into Delta S or Delta I_sub")
    if not float(summary.get("missing_copy_mutation_min_cp_relative_error", -np.inf)) >= 1.0e-4:
        raise RuntimeError("independent audit no longer detects the removed-copy mutation")

    return {
        "task_id": report["task_id"],
        "seed": report["seed"],
        "summary": summary,
        "checks": report["checks"],
        "structural_preflight_passed": True,
        "parameter_case_count": len(report.get("parameter_cases", ())),
        "training_loss_used": False,
        "used_for_parameter_selection": False,
        "used_as_formal_full_domain_pde_gate": False,
        "st006_directly_comparable": False,
    }


def build_checkpoint(audit_summary: dict[str, Any]) -> dict[str, Any]:
    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "base_agent4_head": BASE_AGENT4_HEAD,
        "upstream": {
            "agent1": AGENT1_UPSTREAM,
            "agent2": AGENT2_UPSTREAM,
            "agent3": AGENT3_UPSTREAM,
            "agent4": AGENT4_UPSTREAM,
            "previous_agent5": PREVIOUS_AGENT5,
            "coordination": COORDINATION,
        },
        "independent_heat_v2_robustness": audit_summary,
        "fixed_gates": FIXED_GATES,
        "st006_reference": ST006_REFERENCE,
        "routing": {
            "repeat_delta_cp_tail_factor_work": False,
            "accept_corrected_heat_v2_as_local_source_moment_formula": True,
            "use_corrected_heat_v2_as_physical_three_bump_target_now": False,
            "run_three_bump_repair_without_explicit_X_star_e_star": False,
            "promote_global_leading_field_from_current_heat_data": False,
            "materialize_source_complete_curl_before_actual_background_and_pulse_inverse": False,
            "rerun_agent3_one_column_finite_cycle": False,
            "run_formal_full_domain_gate": False,
            "replace_st006_as_repository_pde_baseline": False,
            "priority_blocker": (
                "Agent 1 must bind the explicit reserved-patch/outer schedule scales, including "
                "X_star/e_star, then feed the independently validated corrected heat discrepancy "
                "through the existing three-bump repair and materialize a globally matched leading "
                "velocity/pressure candidate"
            ),
            "agent2_blocker": (
                "the projected pulse inverse and source complete-curl algebra are executable, but "
                "the actual source pulse forcing/background path and a second public covariance "
                "column are not materialized"
            ),
            "agent3_blocker": (
                "the real one-column finite cycle is rejected; wait for a genuinely independent "
                "second public covariance column that passes the measured rank/coverage contract"
            ),
        },
        "states": STATES,
        "truth_boundary": {
            "source_version_and_provenance_preserved": True,
            "agent5_new_kokuno_mathematics_added": False,
            "corrected_heat_result_promoted_beyond_local_scope": False,
            "sibling_evidence_laundered_into_executable_ancestry": False,
            "free_residual_defined_forcing_used": False,
            "pressure_or_forcing_refit_in_agent5": False,
            "threshold_relaxed": False,
            "st006_used_as_kokuno_source_truth": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
        "next_shortest_closure": [
            "Agent 1: bind X_star/e_star and remaining outer schedule scales; use corrected V2 heat discrepancy with the already preflighted three-bump repair; materialize globally matched leading velocity/pressure",
            "Agent 2: instantiate the actual source pulse/background path and materialize a genuinely independent second public covariance column through the tested complete-curl operator",
            "Agent 3: evaluate that second column against the measured missing-column contract; only then run a real two-column finite correction cycle",
            "Agent 4: once a global candidate exists, run the unchanged held-out full-domain normalized momentum <=1e-3 and divergence <=1e-5 gate with an independent operator",
            "Agent 5: serialize the complete candidate, bind the independent verdict, compare against ST006 on a directly comparable protocol, and run Python/MATLAB export smoke",
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
    if unsigned.get("states") != STATES:
        raise ValueError("scientific state vector changed")
    if unsigned.get("st006_reference") != ST006_REFERENCE:
        raise ValueError("ST006 reference/comparability boundary changed")
    expected_upstream = {
        "agent1": AGENT1_UPSTREAM,
        "agent2": AGENT2_UPSTREAM,
        "agent3": AGENT3_UPSTREAM,
        "agent4": AGENT4_UPSTREAM,
        "previous_agent5": PREVIOUS_AGENT5,
        "coordination": COORDINATION,
    }
    if unsigned.get("upstream") != expected_upstream:
        raise ValueError("upstream provenance changed")

    routing = unsigned.get("routing", {})
    blocked = (
        "repeat_delta_cp_tail_factor_work",
        "use_corrected_heat_v2_as_physical_three_bump_target_now",
        "run_three_bump_repair_without_explicit_X_star_e_star",
        "promote_global_leading_field_from_current_heat_data",
        "materialize_source_complete_curl_before_actual_background_and_pulse_inverse",
        "rerun_agent3_one_column_finite_cycle",
        "run_formal_full_domain_gate",
        "replace_st006_as_repository_pde_baseline",
    )
    if any(routing.get(key) is not False for key in blocked):
        raise ValueError("a blocked routing step was promoted")
    if routing.get("accept_corrected_heat_v2_as_local_source_moment_formula") is not True:
        raise ValueError("resolved local heat-v2 routing was reverted")

    audit = unsigned.get("independent_heat_v2_robustness", {})
    if audit.get("structural_preflight_passed") is not True:
        raise ValueError("Agent-4 heat-v2 robustness preflight is not bound")
    if audit.get("st006_directly_comparable") is not False:
        raise ValueError("local heat audit was relabeled as ST006-comparable")
    if audit.get("used_as_formal_full_domain_pde_gate") is not False:
        raise ValueError("local heat audit was relabeled as a formal PDE gate")

    truth = unsigned.get("truth_boundary", {})
    expected_truth = {
        "source_version_and_provenance_preserved": True,
        "agent5_new_kokuno_mathematics_added": False,
        "corrected_heat_result_promoted_beyond_local_scope": False,
        "sibling_evidence_laundered_into_executable_ancestry": False,
        "free_residual_defined_forcing_used": False,
        "pressure_or_forcing_refit_in_agent5": False,
        "threshold_relaxed": False,
        "st006_used_as_kokuno_source_truth": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
    }
    if truth != expected_truth:
        raise ValueError("truth boundary changed")
    return payload


def load_checkpoint(path: str | Path) -> dict[str, Any]:
    return validate_checkpoint(json.loads(Path(path).read_text(encoding="utf-8")))


def write_bundle(output_dir: str | Path) -> dict[str, Any]:
    target = Path(output_dir)
    if target.exists() and any(target.iterdir()):
        raise FileExistsError("output directory must be absent or empty")
    target.mkdir(parents=True, exist_ok=True)

    full_report = run_audit()
    if full_report.get("structural_preflight_passed") is not True:
        raise RuntimeError("Agent-4 heat-v2 robustness audit failed")
    report_path = target / "heat_v2_robustness_report.json"
    report_path.write_text(
        json.dumps(full_report, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )

    audit_summary = run_bound_heat_v2_audit()
    checkpoint = build_checkpoint(audit_summary)
    checkpoint_path = target / "heat_v2_routing_checkpoint.json"
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
        default="artifacts/kokuno_agent5/heat_v2_routing_checkpoint_v10",
    )
    args = parser.parse_args(argv)
    checkpoint = write_bundle(args.output_dir)
    print(json.dumps(checkpoint, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
