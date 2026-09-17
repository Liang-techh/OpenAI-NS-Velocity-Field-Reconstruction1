from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Mapping

_ALLOWED_STATES = {"pass", "fail", "pending"}


def _finest_rows(validation: Mapping[str, Any], steps: list[float]) -> tuple[float, list[Mapping[str, Any]]]:
    if not steps:
        raise ValueError("validation derivative_steps must be non-empty")
    finest = min(float(step) for step in steps)
    rows = validation.get("rows")
    if not isinstance(rows, list) or not rows:
        raise ValueError("validation rows must be a non-empty list")
    selected = [row for row in rows if abs(float(row["step"]) - finest) <= 1e-15]
    if not selected:
        raise ValueError(f"no validation rows at finest step {finest}")
    return finest, selected


def _entry(
    evidence_id: str,
    evidence_class: str,
    status: str,
    source: str,
    detail: str,
    observed: Any = None,
    threshold: Any = None,
) -> dict[str, Any]:
    if status not in _ALLOWED_STATES:
        raise ValueError(f"invalid evidence status for {evidence_id}: {status}")
    return {
        "id": evidence_id,
        "evidence_class": evidence_class,
        "status": status,
        "source": source,
        "detail": detail,
        "observed": observed,
        "threshold": threshold,
    }


def _metric_entry(evidence_id: str, observed: float, threshold: float, finest: float) -> dict[str, Any]:
    return _entry(
        evidence_id,
        "independent_numeric_evidence",
        "pass" if observed <= threshold else "fail",
        "independent held-out validation rows",
        f"maximum over finest-step rows h={finest}",
        observed,
        threshold,
    )


def build_acceptance_ledger(
    config: Mapping[str, Any],
    validation: Mapping[str, Any],
    structure: Mapping[str, Any],
    *,
    candidate_path: str,
    task_states: Mapping[str, str],
    boundary_support_status: str = "pending",
    generalization_status: str = "pending",
    clean_replay_status: str = "pending",
) -> dict[str, Any]:
    thresholds = config["validation"]["thresholds"]
    finest, rows = _finest_rows(validation, list(config["validation"]["derivative_steps"]))
    metrics = {
        "pde_residual_max": max(float(row["residual_sampled_max"]) for row in rows),
        "pde_residual_L2": max(float(row["residual_L2_estimate"]) for row in rows),
        "divergence_max": max(float(row["divergence_sampled_max"]) for row in rows),
        "divergence_L2": max(float(row["divergence_L2_estimate"]) for row in rows),
    }
    structure_block = validation["structure"]
    core_drift = float(structure_block["core_drift"])
    energy_range = [float(value) for value in structure_block["energy_range"]]
    if len(energy_range) != 2:
        raise ValueError("energy_range must contain [min,max]")
    identities = structure.get("identities", [])
    identity_pass = (
        structure.get("status") == "five_symbolic_identities_verified"
        and isinstance(identities, list)
        and len(identities) >= 3
        and all(item.get("result") is not None for item in identities)
    )

    evidence = [
        _metric_entry("pde_residual_max", metrics["pde_residual_max"], float(thresholds["pde_residual_max"]), finest),
        _metric_entry("pde_residual_L2", metrics["pde_residual_L2"], float(thresholds["pde_residual_L2"]), finest),
        _metric_entry("divergence_max", metrics["divergence_max"], float(thresholds["divergence_max"]), finest),
        _metric_entry("divergence_L2", metrics["divergence_L2"], float(thresholds["divergence_L2"]), finest),
        _entry(
            "sampled_structure",
            "independent_numeric_evidence",
            "pass" if bool(structure_block.get("sampled_constraints_pass")) else "fail",
            "candidate validation structure block",
            "sampled structural checks reported by the validation artifact",
            bool(structure_block.get("sampled_constraints_pass")),
            True,
        ),
        _entry(
            "core_drift",
            "independent_numeric_evidence",
            "pass" if core_drift <= float(thresholds["scaled_core_profile_relative_drift"]) else "fail",
            "candidate validation structure block",
            "scaled-core drift compared with preregistered tolerance",
            core_drift,
            float(thresholds["scaled_core_profile_relative_drift"]),
        ),
        _entry(
            "energy_range",
            "independent_numeric_evidence",
            "pass"
            if energy_range[0] >= float(config["nontriviality"]["minimum_energy_each_validation_time"])
            and energy_range[1] <= float(config["nontriviality"]["maximum_energy_each_validation_time"])
            else "fail",
            "candidate validation structure block",
            "validation-time energy range compared with preregistered nontriviality bounds",
            f"[{energy_range[0]}, {energy_range[1]}]",
            f"[{config['nontriviality']['minimum_energy_each_validation_time']}, {config['nontriviality']['maximum_energy_each_validation_time']}]",
        ),
        _entry(
            "structure_identities",
            "strict_local_proof",
            "pass" if identity_pass else "fail",
            "artifacts/constrained/structure_identities.json",
            "at least three checked symbolic ansatz identities; limited to their stated assumptions",
            len(identities),
            ">=3",
        ),
        _entry(
            "boundary_support_current_candidate",
            "unresolved" if boundary_support_status == "pending" else "independent_numeric_evidence",
            boundary_support_status,
            "current-candidate boundary/support evidence",
            "candidate-specific independent shell/support result, not merely an operator calibration",
        ),
        _entry(
            "cr009_generalization",
            "unresolved" if generalization_status == "pending" else "independent_numeric_evidence",
            generalization_status,
            "CR009 perturbation/generalization evidence",
            "seed/initialization/parameter perturbation generalization remains separately gated",
        ),
        _entry(
            "clean_environment_replay",
            "unresolved" if clean_replay_status == "pending" else "independent_numeric_evidence",
            clean_replay_status,
            "CR011 clean-environment replay",
            "generate/load -> independent validator -> diagnostics replay",
        ),
    ]

    metric_ids = {"pde_residual_max", "pde_residual_L2", "divergence_max", "divergence_L2"}
    independent_fail = any(item["status"] == "fail" for item in evidence if item["id"] in metric_ids)
    declared = validation.get("status")
    if declared not in {"failed_validation", "passed_validation"}:
        raise ValueError(f"unrecognized validation status: {declared}")
    if declared == "passed_validation" and independent_fail:
        raise ValueError("validation artifact claims pass while preregistered finest-step metrics fail")

    tasks = {key: task_states.get(key, "MISSING") for key in ("CR005", "CR009", "CR011", "CR012")}
    dependency_gate_open = all(tasks[key] == "DONE" for key in ("CR005", "CR009", "CR011"))
    acceptance_ready = dependency_gate_open and all(item["status"] == "pass" for item in evidence)
    return {
        "schema_version": 1,
        "candidate": candidate_path,
        "finest_derivative_step": finest,
        "validation_seed": validation.get("seed"),
        "validation_points": validation.get("points"),
        "task_states": tasks,
        "dependency_gate_open": dependency_gate_open,
        "evidence": evidence,
        "acceptance_ready": acceptance_ready,
        "truth_boundary": {
            "paper_exact": False,
            "full_blowup_proof": False,
            "small_residual_implies_singularity": False,
        },
    }


def _task_states(items: list[str]) -> dict[str, str]:
    out: dict[str, str] = {}
    for item in items:
        if "=" not in item:
            raise ValueError(f"task state must be KEY=STATE, got {item!r}")
        key, value = item.split("=", 1)
        out[key] = value
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description="Build a truth-preserving constrained acceptance ledger")
    parser.add_argument("--config", required=True)
    parser.add_argument("--validation", required=True)
    parser.add_argument("--structure", required=True)
    parser.add_argument("--candidate", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--task-state", action="append", default=[])
    parser.add_argument("--boundary-support-status", choices=sorted(_ALLOWED_STATES), default="pending")
    parser.add_argument("--generalization-status", choices=sorted(_ALLOWED_STATES), default="pending")
    parser.add_argument("--clean-replay-status", choices=sorted(_ALLOWED_STATES), default="pending")
    args = parser.parse_args()
    load = lambda path: json.loads(Path(path).read_text(encoding="utf-8"))
    ledger = build_acceptance_ledger(
        load(args.config),
        load(args.validation),
        load(args.structure),
        candidate_path=args.candidate,
        task_states=_task_states(args.task_state),
        boundary_support_status=args.boundary_support_status,
        generalization_status=args.generalization_status,
        clean_replay_status=args.clean_replay_status,
    )
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(ledger, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
