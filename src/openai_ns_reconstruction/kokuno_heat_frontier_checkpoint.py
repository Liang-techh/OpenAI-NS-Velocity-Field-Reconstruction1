"""Agent-5 frontier checkpoint for the latest Kokuno heat-repair evidence.

This module adds no new Kokuno mathematics. It consumes Agent 1's executable
three-row heat-discrepancy repair through Agent 4's independent Simpson/finite-
difference audit, binds the current Agent-2/3/5 routing state, and writes one
fail-closed integration receipt. The actual source heat discrepancy/scales and
a globally matched core-to-heat velocity remain absent, so the formal held-out
Navier--Stokes gate remains unassessed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_independent_heat_repair import (
    DIVERGENCE_REFERENCE,
    LEVELS,
    RESIDUAL_REFERENCE,
    SEED,
    TASK as AGENT4_TASK,
    UPSTREAM_HEAD as AGENT1_HEAD,
    run_audit,
)

SCHEMA = "kokuno-agent5-heat-frontier-checkpoint-v8"
TASK_ID = "KOKUNO-A5-HEAT-FRONTIER-CHECKPOINT-008"
BASE_AGENT4_HEAD = "f7359ce64a0e1719fb63d2d61ea9ec1b7415f689"

AGENT1_UPSTREAM = {
    "pr": 271,
    "head_sha": "d2b20eeb52e1914e2aaa99be9ef6437d819a3e53",
    "heat_discrepancy_repair_executable": True,
    "actual_heat_discrepancy_supplied": False,
    "source_outer_scale_instantiated": False,
    "core_to_heat_matching_completed": False,
}
AGENT2_UPSTREAM = {
    "pr": 272,
    "head_sha": "8dcae5276ebda4279021fddf17cb838ab04e1db2",
    "standard_run": 35287586426,
    "leading_chart_phase_bridge_ready": True,
    "actual_background_V_G_mapping_completed": False,
    "positive_order_background_corrections_included": False,
    "complete_curl_velocity_changed": False,
    "consumed_in_executable_ancestry": False,
}
AGENT3_UPSTREAM = {
    "pr": 263,
    "head_sha": "04f1335dd73eecffb95d18c523f1afb81ee89195",
    "dedicated_run": 35283445927,
    "standard_run": 35283445968,
    "held_out_mean_defect_ratio": 1.15068886,
    "held_out_theta_ratio": 1.08769840,
    "accepted_for_next_cycle": False,
    "consumed_in_executable_ancestry": False,
}
AGENT4_UPSTREAM = {
    "pr": 274,
    "head_sha": BASE_AGENT4_HEAD,
    "dedicated_run": 35288512615,
    "standard_run": 35288512643,
    "artifact_id": 10525560769,
    "artifact_digest": "sha256:8859d6b109e269049801e71141178cad43ce5918c39648e4f1da9f828f9e1b35",
    "consumed_in_executable_ancestry": True,
}
PREVIOUS_AGENT5 = {
    "pr": 268,
    "head_sha": "6e78e34e0eb1b30cc19ad09ef47ea8e6318ac7e1",
    "checkpoint_sha256": "717aa4d9a6ef42d3fa9b1f7cc439301c6d72091fc02f3d7a8278b25aa9fb53be",
    "dedicated_run": 35285465199,
    "standard_run": 35285465095,
    "consumed_in_executable_ancestry": False,
}

REFERENCE_METRICS = {
    "max_finest_public_vs_independent_relative_error": 3.6519317704e-6,
    "max_medium_to_fine_relative_change": 5.6453779296e-14,
    "max_analytic_vs_independent_jacobian_relative_error": 2.3655584468e-6,
    "finest_profile_derivative_relative_error": 1.7732530982e-7,
    "min_parameter_mutation_abs_detection": 4.4131554746e-4,
    "min_independent_increment_norm": 8.8318753261e-4,
}

STATES = {
    "heat_discrepancy_repair_ready": True,
    "independent_heat_repair_preflight_ready": True,
    "leading_chart_phase_bridge_upstream_ready": True,
    "source_phase_composed_into_complete_curl": False,
    "agent3_signed_correction_accepted": False,
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


def summarize_heat_audit(report: dict[str, Any]) -> dict[str, Any]:
    if not isinstance(report, dict) or report.get("task") != AGENT4_TASK:
        raise ValueError("unexpected Agent-4 heat-audit task")
    if report.get("upstream_agent1_head") != AGENT1_HEAD:
        raise ValueError("Agent-1 heat-repair identity changed")
    if int(report.get("seed", -1)) != SEED:
        raise ValueError("Agent-4 heat-audit seed changed")
    if tuple(report.get("levels", ())) != tuple(LEVELS):
        raise ValueError("Agent-4 heat-audit resolution ladder changed")
    if report.get("fixed_references") != {
        "normalized_ns_residual": RESIDUAL_REFERENCE,
        "divergence": DIVERGENCE_REFERENCE,
        "changed": False,
    }:
        raise ValueError("registered references changed")

    summary = report.get("summary", {})
    values = {key: float(summary[key]) for key in REFERENCE_METRICS}
    if not all(np.isfinite(value) for value in values.values()):
        raise ValueError("Agent-4 heat-audit summary contains nonfinite values")
    if summary.get("structural_preflight_passed") is not True:
        raise ValueError("Agent-4 heat-repair structural preflight did not pass")
    for key, expected in REFERENCE_METRICS.items():
        if not np.isclose(values[key], expected, rtol=2.0e-8, atol=2.0e-12):
            raise RuntimeError(
                f"Agent-4 exact-head heat-audit replay drift for {key}: "
                f"{values[key]} != {expected}"
            )

    truth = report.get("truth_boundary", {})
    required_false = (
        "actual_heat_discrepancy_supplied",
        "global_leading_profile_reconstructed",
        "complete_kokuno_composite_velocity",
        "formal_full_domain_pde_gate_assessed",
        "normalized_ns_residual_le_1e-3_claimed",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    )
    if any(truth.get(key) is not False for key in required_false):
        raise ValueError("Agent-4 heat audit crossed its declared truth boundary")

    return {
        **values,
        "structural_preflight_passed": True,
        "sample_count": int(report.get("sample_count", -1)),
        "seed": SEED,
        "levels": list(LEVELS),
    }


def build_checkpoint(agent4_report: dict[str, Any]) -> dict[str, Any]:
    evidence = summarize_heat_audit(agent4_report)
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
        "agent4_independent_heat_repair_audit": {
            "replayed_under_agent5": True,
            "replay_matches_exact_head_reference": True,
            "evidence": evidence,
            "used_for_parameter_selection": False,
            "used_as_formal_full_domain_pde_gate": False,
        },
        "fixed_gates": {
            "held_out_normalized_full_momentum_residual": RESIDUAL_REFERENCE,
            "divergence_max": DIVERGENCE_REFERENCE,
            "changed": False,
        },
        "routing": {
            "instantiate_heat_repair_into_global_velocity": False,
            "compose_agent2_leading_bridge_into_complete_curl": False,
            "apply_agent3_signed_curl_correction": False,
            "run_formal_full_domain_gate": False,
            "priority_blocker": (
                "Agent 1 actual source outer/heat discrepancy and scales, followed by "
                "a matched core-to-heat public velocity/pressure assembly"
            ),
            "heat_repair_blocker": (
                "the local repair operator is independently preflighted, but the actual "
                "source discrepancy, outer scales and eta-dependent coefficient family are absent"
            ),
            "agent2_blocker": (
                "the leading-only V/G chart bridge is executable, but positive-order actual "
                "background corrections are absent and the source phase has not changed the complete curl"
            ),
            "agent3_blocker": (
                "the materialized one-column signed correction worsens held-out mean-defect "
                "and theta targets, so it remains rejected"
            ),
        },
        "states": STATES,
        "truth_boundary": {
            "source_version_and_provenance_preserved": True,
            "agent5_new_kokuno_mathematics_added": False,
            "actual_heat_discrepancy_synthesized": False,
            "source_outer_scales_synthesized": False,
            "free_residual_defined_forcing_used": False,
            "pressure_or_forcing_refit_in_agent5": False,
            "independent_preflight_called_formal_pde_validation": False,
            "threshold_relaxed": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
        "next_shortest_closure": [
            "Agent 1: supply the actual source outer/heat discrepancy and scales, solve the now-preflighted local repair, and materialize a matched core-to-heat public velocity/pressure candidate",
            "Agent 2: keep the leading chart bridge as provenance-ready glue, but compose a source-oriented complete curl only after the actual-background V/G corrections exist",
            "Agent 3: keep the rejected one-column signed correction out of routing until the leading/global blocker is reduced and a source-motivated independent correction column exists",
            "Agent 4: once a globally matched candidate exists, run the unchanged held-out full-domain 1e-3 residual and 1e-5 divergence gate",
            "Agent 5: then serialize the complete candidate and bind Agent-4's verdict into the final reproducible artifact/report",
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
    if unsigned.get("fixed_gates") != {
        "held_out_normalized_full_momentum_residual": RESIDUAL_REFERENCE,
        "divergence_max": DIVERGENCE_REFERENCE,
        "changed": False,
    }:
        raise ValueError("registered gates changed")
    if unsigned.get("states") != STATES:
        raise ValueError("scientific states changed")
    if unsigned.get("upstream") != {
        "agent1": AGENT1_UPSTREAM,
        "agent2": AGENT2_UPSTREAM,
        "agent3": AGENT3_UPSTREAM,
        "agent4": AGENT4_UPSTREAM,
        "previous_agent5": PREVIOUS_AGENT5,
    }:
        raise ValueError("upstream provenance changed")

    audit = unsigned.get("agent4_independent_heat_repair_audit", {})
    if audit.get("replayed_under_agent5") is not True:
        raise ValueError("Agent-4 audit was not replayed")
    if audit.get("replay_matches_exact_head_reference") is not True:
        raise ValueError("Agent-4 replay no longer matches exact-head evidence")
    if audit.get("used_for_parameter_selection") is not False:
        raise ValueError("independent audit cannot be used for parameter selection")
    if audit.get("used_as_formal_full_domain_pde_gate") is not False:
        raise ValueError("local structural preflight cannot be relabeled as a PDE gate")
    evidence = audit.get("evidence", {})
    if evidence.get("structural_preflight_passed") is not True:
        raise ValueError("heat-repair structural preflight state changed")
    for key, expected in REFERENCE_METRICS.items():
        if not np.isclose(float(evidence.get(key, np.nan)), expected, rtol=2.0e-8, atol=2.0e-12):
            raise ValueError(f"heat-repair evidence changed for {key}")

    routing = unsigned.get("routing", {})
    for key in (
        "instantiate_heat_repair_into_global_velocity",
        "compose_agent2_leading_bridge_into_complete_curl",
        "apply_agent3_signed_curl_correction",
        "run_formal_full_domain_gate",
    ):
        if routing.get(key) is not False:
            raise ValueError(f"blocked routing step was promoted: {key}")

    truth = unsigned.get("truth_boundary", {})
    required_truth = {
        "source_version_and_provenance_preserved": True,
        "agent5_new_kokuno_mathematics_added": False,
        "actual_heat_discrepancy_synthesized": False,
        "source_outer_scales_synthesized": False,
        "free_residual_defined_forcing_used": False,
        "pressure_or_forcing_refit_in_agent5": False,
        "independent_preflight_called_formal_pde_validation": False,
        "threshold_relaxed": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
    }
    if truth != required_truth:
        raise ValueError("truth boundary changed")
    return payload


def load_checkpoint(path: str | Path) -> dict[str, Any]:
    return validate_checkpoint(json.loads(Path(path).read_text(encoding="utf-8")))


def write_bundle(output_dir: str | Path) -> dict[str, Any]:
    target = Path(output_dir)
    if target.exists() and any(target.iterdir()):
        raise FileExistsError("output directory must be absent or empty")
    target.mkdir(parents=True, exist_ok=True)

    report = run_audit()
    report_path = target / "independent_heat_repair_report.json"
    report_path.write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )

    checkpoint = build_checkpoint(report)
    checkpoint_path = target / "heat_frontier_checkpoint.json"
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
        default="artifacts/kokuno_agent5/heat_frontier_checkpoint_v8",
    )
    args = parser.parse_args(argv)
    checkpoint = write_bundle(args.output_dir)
    print(json.dumps(checkpoint, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
