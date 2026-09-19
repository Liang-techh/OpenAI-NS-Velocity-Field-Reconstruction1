"""Fail-closed governance for the canonical function-first velocity deliverable.

This module does not evaluate Navier--Stokes residuals. It keeps the current
callable/save-load delivery identity separate from legacy compatibility,
visual correspondence, PDE acceptance, paper-exact reconstruction, and a
blow-up theorem.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping

_ALLOWED_CLASSES = {
    "user_requirement",
    "public_source_fact",
    "autonomous_design",
    "pending_unknown",
}
_REQUIRED_CLASSES = _ALLOWED_CLASSES
_REQUIRED_FALSE_CLAIMS = (
    "visualization_ready",
    "visual_correspondence_verified",
    "pde_validated",
    "paper_exact",
    "openai_field_identified",
    "blowup_proved",
)
_CANONICAL_API = "openai_ns_reconstruction.eq45_supported_delivery:velocity"
_CANONICAL_FAMILY = "eq45_supported_velocity_candidate_v1"
_LEGACY_API = "openai_ns_reconstruction.velocity_components:velocity"


def _read_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def audit_velocity_delivery_contract(
    contract: Mapping[str, Any], *, repo_root: str | Path
) -> dict[str, Any]:
    """Validate delivery identity and claim boundaries against live artifacts.

    A pass means only that the canonical delivery metadata is internally
    consistent. It is not a PDE or visual-correspondence result.
    """

    root = Path(repo_root)
    _require(contract.get("schema_version") == 2, "unsupported schema_version")

    deliverable = contract.get("primary_deliverable")
    _require(isinstance(deliverable, Mapping), "missing primary_deliverable")
    _require(deliverable.get("canonical") is True, "primary deliverable is not canonical")
    _require(deliverable.get("api") == _CANONICAL_API, "unexpected primary velocity API")
    _require(
        deliverable.get("candidate_family") == _CANONICAL_FAMILY,
        "unexpected primary candidate family",
    )
    _require(deliverable.get("component_order") == ["u", "v", "w"], "component order drift")

    documentation_path = root / str(deliverable.get("documentation", ""))
    capsule_path = root / str(deliverable.get("delivery_capsule", ""))
    _require(documentation_path.is_file(), "velocity API documentation is missing")
    _require(capsule_path.is_file(), "delivery capsule is missing")
    _require((root / "project_status.json").is_file(), "project status is missing")

    capsule = _read_json(capsule_path)
    project_status = _read_json(root / "project_status.json")
    candidate = capsule.get("candidate", {})
    public_interface = capsule.get("public_interface", {})

    _require(
        deliverable.get("candidate_family")
        == project_status.get("candidate_family")
        == candidate.get("family"),
        "canonical candidate family identity drift",
    )
    _require(
        deliverable.get("candidate_sha256")
        == project_status.get("candidate_sha256")
        == candidate.get("child_sha256"),
        "canonical candidate SHA identity drift",
    )
    _require(
        deliverable.get("api")
        == project_status.get("velocity_api")
        == public_interface.get("velocity"),
        "canonical velocity API identity drift",
    )
    _require(
        deliverable.get("grid_api")
        == project_status.get("grid_api")
        == public_interface.get("grid"),
        "canonical grid API identity drift",
    )
    _require(
        deliverable.get("load_api")
        == project_status.get("candidate_save_load_api")
        == public_interface.get("load_candidate"),
        "canonical load API identity drift",
    )

    legacy = contract.get("legacy_compatibility_delivery")
    _require(isinstance(legacy, Mapping), "missing legacy compatibility delivery")
    _require(legacy.get("canonical") is False, "legacy delivery cannot be canonical")
    _require(legacy.get("api") == _LEGACY_API, "legacy velocity API identity drift")
    _require(legacy.get("candidate_family") == "coupled_velocity_v1", "legacy candidate family drift")
    _require(legacy.get("cli") == "ns-velocity", "legacy CLI identity drift")

    source_rows = contract.get("source_classification")
    _require(isinstance(source_rows, list) and source_rows, "source_classification must be nonempty")
    observed_classes: set[str] = set()
    for row in source_rows:
        _require(isinstance(row, Mapping), "source classification row must be an object")
        classification = row.get("classification")
        _require(classification in _ALLOWED_CLASSES, f"invalid source classification: {classification!r}")
        observed_classes.add(str(classification))
        source = row.get("source")
        _require(isinstance(source, str) and source, "source classification row lacks source")
        _require((root / source).is_file(), f"source path does not exist: {source}")
    _require(observed_classes == _REQUIRED_CLASSES, "source classification coverage is incomplete")

    gates = contract.get("claim_gates")
    status = contract.get("claim_status")
    _require(isinstance(gates, Mapping), "missing claim_gates")
    _require(isinstance(status, Mapping), "missing claim_status")
    _require(status.get("velocity_export_ready") is True, "velocity export is not marked ready")

    project_states = project_status.get("states", {})
    capsule_states = capsule.get("states", {})
    _require(project_states.get("velocity_export_ready") is True, "project velocity export is not ready")
    _require(capsule_states.get("velocity_export_ready") is True, "capsule velocity export is not ready")

    visual_gate = gates.get("visual_correspondence")
    _require(isinstance(visual_gate, Mapping), "missing visual_correspondence gate")
    _require(
        visual_gate.get("blocking_for_velocity_delivery") is False,
        "visual correspondence must not block callable velocity delivery",
    )
    if status.get("visual_correspondence_verified") is True:
        _require(
            visual_gate.get("time_mapping_status") == "verified",
            "verified correspondence lacks time mapping",
        )
        evidence = visual_gate.get("comparison_evidence")
        _require(
            isinstance(evidence, list) and evidence,
            "verified correspondence lacks comparison evidence",
        )

    pde_gate = gates.get("pde_validation")
    _require(isinstance(pde_gate, Mapping), "missing pde_validation gate")
    _require(
        pde_gate.get("blocking_for_velocity_delivery") is False,
        "PDE validation must not block callable velocity delivery",
    )
    _require(pde_gate.get("blocking_for_pde_claim") is True, "PDE evidence must gate PDE claims")
    pde_config_path = root / str(pde_gate.get("config", ""))
    _require(pde_config_path.is_file(), "active PDE config is missing")
    pde_config = _read_json(pde_config_path)

    _require(
        pde_config.get("domain", {}).get("time_interval") == deliverable.get("time_interval"),
        "velocity delivery time interval drifted from active config",
    )
    _require(
        pde_config.get("domain", {}).get("support") == deliverable.get("support"),
        "velocity delivery support drifted from active config",
    )
    if status.get("pde_validated") is True:
        evidence = pde_gate.get("evidence")
        _require(isinstance(evidence, list) and evidence, "PDE-valid status lacks independent evidence")

    for key in _REQUIRED_FALSE_CLAIMS:
        _require(status.get(key) is False, f"unsupported claim promotion: {key}")
        _require(project_states.get(key) is False, f"project status promoted unsupported claim: {key}")
        _require(capsule_states.get(key) is False, f"delivery capsule promoted unsupported claim: {key}")

    guard = contract.get("cr001_nonmutation")
    _require(isinstance(guard, Mapping), "missing CR001 non-mutation guard")
    for key in (
        "scientific_parameters_changed",
        "thresholds_changed",
        "forcing_contract_changed",
        "validation_sample_changed",
    ):
        _require(guard.get(key) is False, f"CR001 mutation flag must remain false: {key}")

    return {
        "contract_pass": True,
        "velocity_export_ready": True,
        "visualization_ready": False,
        "visual_correspondence_verified": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
        "candidate_family": deliverable.get("candidate_family"),
        "candidate_sha256": deliverable.get("candidate_sha256"),
        "velocity_api": deliverable.get("api"),
        "time_interval": deliverable.get("time_interval"),
        "support": deliverable.get("support"),
        "source_classes": sorted(observed_classes),
    }


def audit_velocity_delivery_contract_file(path: str | Path, *, repo_root: str | Path) -> dict[str, Any]:
    return audit_velocity_delivery_contract(_read_json(Path(path)), repo_root=repo_root)
