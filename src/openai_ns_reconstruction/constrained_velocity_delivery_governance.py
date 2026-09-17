"""Fail-closed governance for the function-first velocity deliverable.

This module does not evaluate Navier--Stokes residuals. It only prevents the
callable u/v/w deliverable from being conflated with visual correspondence,
PDE acceptance, paper-exact reconstruction, or a blow-up theorem.
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
    "openai_numerical_field_identified",
    "openai_correspondence_verified",
    "pde_validated",
    "paper_exact",
    "blowup_proof",
)


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
    """Validate the delivery/claim boundary against repository artifacts.

    A pass means only that scope metadata is internally consistent with the
    packaged velocity API and active constrained config. It is not a PDE or
    visual-correspondence result.
    """

    root = Path(repo_root)
    _require(contract.get("schema_version") == 1, "unsupported schema_version")

    deliverable = contract.get("primary_deliverable")
    _require(isinstance(deliverable, Mapping), "missing primary_deliverable")
    _require(
        deliverable.get("api") == "openai_ns_reconstruction.velocity_components:velocity",
        "unexpected primary velocity API",
    )
    _require(deliverable.get("component_order") == ["u", "v", "w"], "component order drift")
    _require(deliverable.get("component_apis") == ["u", "v", "w"], "component API drift")

    candidate_path = root / str(deliverable.get("bundled_candidate", ""))
    documentation_path = root / str(deliverable.get("documentation", ""))
    _require(candidate_path.is_file(), "bundled candidate is missing")
    _require(documentation_path.is_file(), "velocity API documentation is missing")
    candidate = _read_json(candidate_path)
    _require(
        candidate.get("family") == deliverable.get("bundled_candidate_family"),
        "bundled candidate family does not match contract",
    )

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
    _require(status.get("public_target_asset_identified") is True, "public target asset is not identified")

    visual_gate = gates.get("visual_correspondence")
    _require(isinstance(visual_gate, Mapping), "missing visual_correspondence gate")
    _require(
        visual_gate.get("blocking_for_velocity_delivery") is False,
        "visual correspondence must not block callable velocity delivery",
    )
    if status.get("openai_correspondence_verified") is True:
        _require(visual_gate.get("time_mapping_status") == "verified", "verified correspondence lacks time mapping")
        evidence = visual_gate.get("comparison_evidence")
        _require(isinstance(evidence, list) and evidence, "verified correspondence lacks comparison evidence")

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

    return {
        "contract_pass": True,
        "velocity_export_ready": True,
        "visual_correspondence_verified": False,
        "pde_validated": False,
        "paper_exact": False,
        "blowup_proof": False,
        "candidate_family": candidate.get("family"),
        "time_interval": deliverable.get("time_interval"),
        "support": deliverable.get("support"),
        "source_classes": sorted(observed_classes),
    }


def audit_velocity_delivery_contract_file(path: str | Path, *, repo_root: str | Path) -> dict[str, Any]:
    return audit_velocity_delivery_contract(_read_json(Path(path)), repo_root=repo_root)
