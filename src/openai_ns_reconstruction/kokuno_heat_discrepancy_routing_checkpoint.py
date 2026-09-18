"""Fail-closed Kokuno Agent-5 checkpoint for the current heat-discrepancy frontier.

This module adds no new Kokuno mathematics.  It reruns Agent 4's independent
heat-replacement discrepancy audit against Agent 1's public evaluator and binds
the latest Agent 2 / Agent 3 handoff state into one deterministic receipt.

The current Agent-4 audit detects a concrete first-order ``Delta C_p`` tail
factor defect in Agent 1 #280.  Consequently the public heat discrepancy is
*not* eligible to drive the three-bump repair or a global leading-field
assembly until Agent 1 fixes the factor and reruns its own deterministic
checkpoint.  This local routing decision does not assess the full Navier--Stokes
gate and does not alter the fixed ``1e-3`` residual or ``1e-5`` divergence
references.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_heat_replacement_discrepancy import KokunoHeatReplacementDiscrepancy
from .kokuno_heat_replacement_discrepancy_independent import (
    IndependentHeatDiscrepancyAudit,
)

SCHEMA = "kokuno-agent5-heat-discrepancy-routing-checkpoint-v9"
TASK_ID = "KOKUNO-A5-HEAT-DISCREPANCY-ROUTING-CHECKPOINT-009"
BASE_AGENT4_HEAD = "ecf1be298e880afe62b1f568a8f42bd27bc54d4a"

AGENT1_UPSTREAM = {
    "pr": 280,
    "head_sha": "51076990abdd36bcb7dc598a7abd85d677dee53d",
    "standard_run": 35291190773,
    "actual_heat_replacement_discrepancy_formula_executable": True,
    "source_outer_schedule_numeric_choices_recovered": False,
    "earlier_reserved_patch_scales_recovered": False,
    "consumed_in_executable_ancestry": True,
}
AGENT2_UPSTREAM = {
    "pr": 281,
    "head_sha": "4f9be3a159f02d0c3d4d30b07f847d02c3304f43",
    "standard_run": 35291927653,
    "source_complete_curl_operator_ready": True,
    "source_pulse_inverse_reconstructed": False,
    "actual_background_V_G_mapping_completed": False,
    "public_velocity_correction_materialized": False,
    "consumed_in_executable_ancestry": False,
}
AGENT3_UPSTREAM = {
    "pr": 282,
    "head_sha": "56501232583f7e08f02c223e9b9e5b9d6ca94c27",
    "dedicated_run": 35292279365,
    "standard_run": 35292279284,
    "artifact_id": 10527151558,
    "artifact_digest": "sha256:56107f975bfcad3157eda8b14fc59fdeba3d2227a26a1930a8637814f6197874",
    "missing_second_column_target_measured": True,
    "active_nodes": 27,
    "active_nodes_requiring_second_direction": 25,
    "pointwise_rank_one_relative_rms_remainder": 0.2101945021,
    "new_oscillatory_column_constructed": False,
    "consumed_in_executable_ancestry": False,
}
AGENT4_UPSTREAM = {
    "pr": 283,
    "head_sha": BASE_AGENT4_HEAD,
    "standard_run": 35292637019,
    "independent_heat_discrepancy_audit_ready": True,
    "consumed_in_executable_ancestry": True,
}
PREVIOUS_AGENT5 = {
    "pr": 276,
    "head_sha": "97f7f4e9067b3cc5f3e0a479b6bcbb1a1a0e5bac",
    "heat_discrepancy_repair_ready": True,
    "independent_heat_repair_preflight_ready": True,
    "consumed_in_executable_ancestry": False,
}

FIXED_GATES = {
    "held_out_normalized_full_momentum_residual": 1.0e-3,
    "divergence_max": 1.0e-5,
    "changed": False,
}

STATES = {
    "actual_heat_discrepancy_formula_executable": True,
    "independent_heat_discrepancy_audit_ready": True,
    "heat_discrepancy_ready": False,
    "heat_discrepancy_promotion_blocked": True,
    "source_complete_curl_operator_upstream_ready": True,
    "missing_second_covariance_column_target_upstream_ready": True,
    "source_oscillatory_public_velocity_ready": False,
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


def run_independent_heat_audit() -> dict[str, Any]:
    """Rerun Agent 4's public-vs-independent heat discrepancy comparison."""

    model = KokunoHeatReplacementDiscrepancy()
    audit = IndependentHeatDiscrepancyAudit()
    report = audit.compare_public(model)

    rel = np.asarray(report["max_relative_by_component"], dtype=float)
    ratios = np.asarray(report["cp_missing_half_factor_ratio"], dtype=float)
    public = np.asarray(report["public"], dtype=float)
    reference = np.asarray(report["independent_reference"], dtype=float)
    if rel.shape != (3,) or public.shape != reference.shape or public.shape[-1] != 3:
        raise ValueError("unexpected Agent-4 independent heat-audit shape")
    if not np.all(np.isfinite(rel)) or not np.all(np.isfinite(ratios)):
        raise ValueError("independent heat audit contains nonfinite values")
    if not (2.3e-3 < rel[0] < 2.6e-3):
        raise RuntimeError("Delta C_p discrepancy no longer matches the independently diagnosed factor bug")
    if not (rel[1] < 2.0e-8 and rel[2] < 2.0e-6):
        raise RuntimeError("independent heat audit drifted in Delta S or Delta I_sub")
    if not np.allclose(ratios, 1.0, rtol=0.0, atol=3.0e-5):
        raise RuntimeError("missing-half-tail attribution no longer closes")
    if report.get("formal_full_domain_pde_gate_assessed") is not False:
        raise ValueError("local independent heat audit cannot become the formal PDE gate")
    if report.get("pde_validated") is not False:
        raise ValueError("local independent heat audit cannot promote PDE validity")

    return {
        "etas": report["etas"],
        "public": report["public"],
        "independent_reference": report["independent_reference"],
        "max_relative_by_component": report["max_relative_by_component"],
        "cp_missing_half_factor_ratio": report["cp_missing_half_factor_ratio"],
        "diagnosis": report["diagnosis"],
        "agent1_heat_discrepancy_eligible_for_repair_target": False,
        "used_for_parameter_selection": False,
        "used_as_formal_full_domain_pde_gate": False,
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
        },
        "independent_heat_discrepancy_audit": audit_summary,
        "fixed_gates": FIXED_GATES,
        "routing": {
            "use_agent1_280_heat_discrepancy_as_repair_target": False,
            "run_three_bump_repair_from_agent1_280": False,
            "promote_global_leading_field_from_current_heat_data": False,
            "materialize_source_complete_curl_before_actual_background_and_pulse_inverse": False,
            "rerun_agent3_one_column_finite_cycle": False,
            "run_formal_full_domain_gate": False,
            "priority_blocker": (
                "Agent 1 must correct the missing factor in the first-order y>=3 Delta C_p tail, "
                "rerun its deterministic heat-discrepancy checkpoint, then supply the still-missing "
                "outer schedule/reserved-patch scales before the three-bump repair can be instantiated"
            ),
            "agent2_blocker": (
                "source complete-curl algebra is executable, but actual-background V/G positive-order "
                "corrections and the source pulse inverse are still absent"
            ),
            "agent3_blocker": (
                "the measured target still needs a genuinely independent second public covariance column; "
                "25 of 27 active nodes require a second direction"
            ),
        },
        "states": STATES,
        "truth_boundary": {
            "source_version_and_provenance_preserved": True,
            "agent5_new_kokuno_mathematics_added": False,
            "agent1_detected_heat_bug_hidden": False,
            "sibling_evidence_laundered_into_executable_ancestry": False,
            "free_residual_defined_forcing_used": False,
            "pressure_or_forcing_refit_in_agent5": False,
            "threshold_relaxed": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
        "next_shortest_closure": [
            "Agent 1: fix the Delta C_p first-order tail factor, rerun the independent comparison, and only then expose a repair target tied to explicit X_star/e_star schedule scales",
            "Agent 2: keep the source complete-curl operator ready, but do not materialize the source oscillatory public velocity until actual-background V/G corrections and the pulse inverse exist",
            "Agent 3: evaluate a future genuinely independent second public column against the measured missing covariance target before any new finite correction cycle",
            "Agent 4: independently re-audit the corrected heat discrepancy and run the unchanged full-domain PDE gate only after a globally matched candidate exists",
            "Agent 5: integrate the corrected heat target, source oscillatory column and accepted finite correction only after those upstream gates pass",
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
    if unsigned.get("upstream") != {
        "agent1": AGENT1_UPSTREAM,
        "agent2": AGENT2_UPSTREAM,
        "agent3": AGENT3_UPSTREAM,
        "agent4": AGENT4_UPSTREAM,
        "previous_agent5": PREVIOUS_AGENT5,
    }:
        raise ValueError("upstream provenance changed")

    routing = unsigned.get("routing", {})
    blocked = (
        "use_agent1_280_heat_discrepancy_as_repair_target",
        "run_three_bump_repair_from_agent1_280",
        "promote_global_leading_field_from_current_heat_data",
        "materialize_source_complete_curl_before_actual_background_and_pulse_inverse",
        "rerun_agent3_one_column_finite_cycle",
        "run_formal_full_domain_gate",
    )
    if any(routing.get(key) is not False for key in blocked):
        raise ValueError("a blocked routing step was promoted")

    audit = unsigned.get("independent_heat_discrepancy_audit", {})
    rel = np.asarray(audit.get("max_relative_by_component", ()), dtype=float)
    ratios = np.asarray(audit.get("cp_missing_half_factor_ratio", ()), dtype=float)
    if rel.shape != (3,) or not (2.3e-3 < rel[0] < 2.6e-3):
        raise ValueError("independent Delta C_p diagnosis changed")
    if not (rel[1] < 2.0e-8 and rel[2] < 2.0e-6):
        raise ValueError("independent non-Cp channels changed")
    if not np.allclose(ratios, 1.0, rtol=0.0, atol=3.0e-5):
        raise ValueError("independent missing-half-tail attribution changed")
    if audit.get("agent1_heat_discrepancy_eligible_for_repair_target") is not False:
        raise ValueError("buggy Agent-1 heat discrepancy was promoted")
    if audit.get("used_as_formal_full_domain_pde_gate") is not False:
        raise ValueError("local heat audit was relabeled as a formal PDE gate")

    truth = unsigned.get("truth_boundary", {})
    expected_truth = {
        "source_version_and_provenance_preserved": True,
        "agent5_new_kokuno_mathematics_added": False,
        "agent1_detected_heat_bug_hidden": False,
        "sibling_evidence_laundered_into_executable_ancestry": False,
        "free_residual_defined_forcing_used": False,
        "pressure_or_forcing_refit_in_agent5": False,
        "threshold_relaxed": False,
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

    audit_summary = run_independent_heat_audit()
    audit_path = target / "independent_heat_discrepancy_audit.json"
    audit_path.write_text(
        json.dumps(audit_summary, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )

    checkpoint = build_checkpoint(audit_summary)
    checkpoint_path = target / "heat_discrepancy_routing_checkpoint.json"
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
        default="artifacts/kokuno_agent5/heat_discrepancy_routing_checkpoint_v9",
    )
    args = parser.parse_args(argv)
    checkpoint = write_bundle(args.output_dir)
    print(json.dumps(checkpoint, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
