"""Fail-closed accounting audit for the integrated ST052-M whole-child runtime.

This module deliberately separates executable facts from repository status prose.
At the CR-A9-070 integration head the ST052-M linear-temporal whole child is
materialized and reloadable when its authenticated exact historical source
runtime is supplied, while ``project_status.json`` still contains older wording
that says the child has not been materialized.  Neither side may be used to
launder a stronger delivery, visual, PDE, or exact-field claim.
"""
from __future__ import annotations

import inspect
import json
from pathlib import Path
from typing import Any

from . import st052_linear_temporal_capsule as capsule

ROOT = Path(__file__).resolve().parents[2]
CONTRACT_PATH = ROOT / "configs" / "st052m_live_delivery_state_accounting_contract.json"

_FALSE_ST052_STATES = (
    "st052_velocity_export_ready",
    "st052_production_candidate_selected",
    "st052_visualization_ready",
    "st052_visual_correspondence_verified",
    "st052_held_out_temporal_child_pde_residual_evaluated",
    "st052_pde_validated",
    "st052_source_correspondence_verified",
    "st052_paper_exact",
    "st052_openai_field_identified",
    "st052_blowup_proved",
)
_EXPECTED_CLASSES = {
    "user_requirement",
    "public_source_fact",
    "autonomous_design",
    "pending_unknown",
}


def _read_json(path: Path) -> dict[str, Any]:
    obj = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(obj, dict):
        raise ValueError(f"expected JSON object: {path}")
    return obj


def load_contract(path: str | Path = CONTRACT_PATH) -> dict[str, Any]:
    return _read_json(Path(path))


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def validate_contract(contract: dict[str, Any]) -> None:
    """Validate the governance object independently of repository files."""
    _require(contract.get("schema_version") == 1, "unsupported accounting schema")
    _require(
        contract.get("task_id") == "CR002-ST052M-LIVE-DELIVERY-STATE-ACCOUNTING-072",
        "accounting task identity drifted",
    )

    integrated = contract.get("integrated_st052_runtime")
    _require(isinstance(integrated, dict), "integrated ST052 runtime section missing")
    _require(
        integrated.get("candidate_id") == "ST052-M-linear-temporal-child-v1",
        "ST052 whole-child candidate identity drifted",
    )
    for key in (
        "whole_child_bundle_materialized",
        "whole_child_save_load_ready_with_authenticated_exact_source_runtime",
        "exact_source_runtime_identity_closed",
        "historical_module_cache_isolated",
    ):
        _require(integrated.get(key) is True, f"integrated positive fact missing: {key}")
    _require(
        integrated.get("standalone_package_parent_runtime_ready") is False,
        "standalone package parent runtime must remain false in this accounting version",
    )

    observation = contract.get("status_reporting_observation")
    _require(isinstance(observation, dict), "status-reporting observation missing")
    _require(
        observation.get("project_status_st052_materialization_fields_stale") is True,
        "v1 must explicitly record the known project-status materialization staleness",
    )

    states = contract.get("claim_states")
    _require(isinstance(states, dict), "claim state section missing")
    _require(
        states.get("canonical_eq45_velocity_export_ready") is True,
        "canonical Eq45 delivery readiness must remain independent and true",
    )
    _require(
        states.get("st052_source_runtime_backed_callable_save_load_ready") is True,
        "integrated source-runtime-backed ST052 callable/save-load fact missing",
    )
    for key in _FALSE_ST052_STATES:
        _require(states.get(key) is False, f"premature ST052 claim promotion: {key}")

    source_classes = contract.get("source_classification")
    _require(isinstance(source_classes, list) and len(source_classes) == 4, "source classification must contain four governed classes")
    classes = {entry.get("classification") for entry in source_classes if isinstance(entry, dict)}
    _require(classes == _EXPECTED_CLASSES, "source classification classes drifted")

    cr001 = contract.get("cr001_nonmutation")
    _require(isinstance(cr001, dict), "CR001 nonmutation section missing")
    _require(cr001.get("thresholds_changed") is False, "CR001 threshold relaxation is forbidden")
    _require(cr001.get("forcing_contract_changed") is False, "CR001 forcing-contract mutation is forbidden")
    _require(cr001.get("validation_sample_changed") is False, "CR001 validation-sample mutation is forbidden")


def _audit_cr001(contract: dict[str, Any], constraints: dict[str, Any]) -> None:
    expected = contract["cr001_nonmutation"]
    domain = constraints["domain"]
    forcing = constraints["forcing"]
    nontriviality = constraints["nontriviality"]
    validation = constraints["validation"]
    thresholds = validation["thresholds"]

    _require(constraints.get("nu") == expected["nu"], "nu drifted")
    _require(domain.get("physical") == expected["physical_domain"], "physical domain drifted")
    _require(domain.get("evaluation_box") == expected["evaluation_box"], "evaluation box drifted")
    _require(domain.get("support") == expected["support"], "support contract drifted")
    _require(domain.get("time_interval") == expected["time_interval"], "time interval drifted")
    _require(forcing.get("mode") == expected["forcing_mode"], "forcing mode drifted")
    _require(forcing["parameters"].get("a") == expected["forcing_parameter_bounds"]["a"], "forcing a bounds drifted")
    _require(forcing["parameters"].get("c") == expected["forcing_parameter_bounds"]["c"], "forcing c bounds drifted")
    _require(nontriviality.get("reference_energy") == expected["reference_energy"], "reference energy drifted")
    _require(
        nontriviality.get("reference_energy_abs_tolerance") == expected["reference_energy_abs_tolerance"],
        "reference energy tolerance drifted",
    )
    _require(validation.get("seed") == expected["validation_seed"], "validation seed drifted")
    _require(validation.get("held_out_points") == expected["held_out_points"], "held-out sample count drifted")
    _require(validation.get("times") == expected["validation_times"], "validation times drifted")
    _require(validation.get("derivative_steps") == expected["derivative_steps"], "derivative ladder drifted")
    _require(
        validation.get("quadrature_orders_per_axis") == expected["quadrature_orders_per_axis"],
        "quadrature ladder drifted",
    )
    for key in ("divergence_max", "divergence_L2", "pde_residual_max", "pde_residual_L2"):
        _require(thresholds.get(key) == expected[key], f"CR001 threshold drifted: {key}")
    _require("No residual-dependent basis or pointwise free force" in forcing.get("restriction", ""), "free residual-defined forcing prohibition missing")


def audit_repository(
    *,
    root: str | Path = ROOT,
    contract: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Audit live executable facts, stale status wording, and immutable CR001 gates."""
    root = Path(root)
    contract = load_contract(root / "configs" / CONTRACT_PATH.name) if contract is None else contract
    validate_contract(contract)

    status = _read_json(root / "project_status.json")
    delivery = _read_json(root / "configs" / "velocity_delivery_contract.json")
    constraints = _read_json(root / "configs" / "constraints.json")

    integrated = contract["integrated_st052_runtime"]
    _require(capsule.SCHEMA == integrated["whole_candidate_schema"], "live whole-candidate schema drifted")
    _require(capsule.TASK_ID == integrated["task_id"], "live whole-candidate task identity drifted")

    build_source = inspect.getsource(capsule.build_bundle)
    load_source = inspect.getsource(capsule.load_bundle_runtime)
    parent_source = inspect.getsource(capsule._load_exact_source_parent)
    for token in (
        '"candidate_id": "ST052-M-linear-temporal-child-v1"',
        '"whole_child_bundle_materialized": True',
        '"whole_child_save_load_ready_with_exact_source_runtime": True',
        '"exact_source_runtime_identity_closed": True',
        '"historical_module_cache_isolated": True',
        '"velocity_export_ready": False',
        '"pde_validated": False',
    ):
        _require(token in build_source, f"live ST052 capsule truth boundary changed: {token}")
    _require("_load_exact_source_parent" in load_source, "whole-child loader no longer binds the parent runtime")
    _require("authenticate_source_runtime" in parent_source, "exact-source runtime authentication disappeared")
    _require("verify_source_module_cache" in parent_source, "historical module-cache verification disappeared")

    observation = contract["status_reporting_observation"]
    latest = status.get("latest_st052m_visual_candidate_evidence")
    _require(isinstance(latest, dict), "project status ST052 evidence block missing")
    _require(
        latest.get("status") == observation["observed_latest_st052_status"],
        "project_status ST052 materialization wording changed; revise this transitional accounting contract",
    )
    _require(
        str(status.get("next_delivery_task", "")).startswith(observation["observed_next_delivery_task_prefix"]),
        "project_status next-delivery wording changed; revise this transitional accounting contract",
    )
    for key in ("velocity_export_ready", "visualization_ready", "visual_correspondence_verified", "pde_validated"):
        _require(latest.get(key) is False, f"stale ST052 status block unexpectedly promoted {key}")

    primary = delivery.get("primary_deliverable")
    claims = delivery.get("claim_status")
    _require(isinstance(primary, dict) and isinstance(claims, dict), "canonical delivery contract malformed")
    _require(primary.get("candidate_family") == "eq45_supported_velocity_candidate_v1", "canonical delivery family drifted")
    _require(primary.get("api") == "openai_ns_reconstruction.eq45_supported_delivery:velocity", "canonical velocity API drifted")
    _require(claims.get("velocity_export_ready") is True, "canonical Eq45 delivery readiness drifted")
    for key in ("visualization_ready", "visual_correspondence_verified", "pde_validated", "paper_exact", "openai_field_identified", "blowup_proved"):
        _require(claims.get(key) is False, f"canonical claim unexpectedly promoted: {key}")

    _audit_cr001(contract, constraints)
    return {
        "task_id": contract["task_id"],
        "integration_base": contract["integration_base"],
        "st052_candidate_id": integrated["candidate_id"],
        "st052_source_runtime_backed_callable_save_load_ready": True,
        "project_status_st052_materialization_fields_stale": True,
        "canonical_eq45_velocity_export_ready": True,
        "st052_velocity_export_ready": False,
        "visual_correspondence_verified": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "cr001_unchanged": True,
    }


def main() -> None:
    print(json.dumps(audit_repository(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
