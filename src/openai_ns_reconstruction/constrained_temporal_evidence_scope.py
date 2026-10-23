"""Fail-closed governance for temporal evidence scopes.

Equal numeric time windows do not make different evidence roles interchangeable.
This module separates the CR001 PDE-validation window, callable velocity domain,
autonomous visualization snapshot times, and the still-pending mapping to public
OpenAI visualization frames. Passing this audit is governance only: it does
not validate the PDE or visual correspondence.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping

_ALLOWED_CLASSES = {"autonomous", "pending"}


def _read_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _ordered_interval(value: Any, *, label: str) -> list[float]:
    _require(
        isinstance(value, list)
        and len(value) == 2
        and all(isinstance(item, (int, float)) for item in value),
        f"{label} must be a two-number interval",
    )
    interval = [float(value[0]), float(value[1])]
    _require(interval[0] < interval[1], f"{label} must be strictly ordered")
    return interval


def _times(value: Any, *, label: str) -> list[float]:
    _require(
        isinstance(value, list)
        and value
        and all(isinstance(item, (int, float)) for item in value),
        f"{label} must be a nonempty numeric list",
    )
    times = [float(item) for item in value]
    _require(all(a < b for a, b in zip(times, times[1:])), f"{label} must be strictly increasing")
    return times


def audit_temporal_evidence_scope(
    contract: Mapping[str, Any], *, repo_root: str | Path
) -> dict[str, Any]:
    """Audit temporal source classification and fail-closed evidence transfer."""

    root = Path(repo_root)
    _require(contract.get("schema_version") == 1, "unsupported schema_version")
    scopes = contract.get("scopes")
    _require(isinstance(scopes, Mapping), "missing temporal scopes")

    pde = scopes.get("pde_validation")
    delivery = scopes.get("velocity_delivery")
    visual = scopes.get("visualization_reference")
    public = scopes.get("public_visual_time_mapping")
    for label, scope in (
        ("pde_validation", pde),
        ("velocity_delivery", delivery),
        ("visualization_reference", visual),
        ("public_visual_time_mapping", public),
    ):
        _require(isinstance(scope, Mapping), f"missing {label} scope")
        classification = scope.get("classification")
        _require(classification in _ALLOWED_CLASSES, f"invalid {label} classification")

    _require(pde.get("classification") == "autonomous", "PDE finite window is an autonomous choice")
    _require(delivery.get("classification") == "autonomous", "delivery time window is autonomous")
    _require(visual.get("classification") == "autonomous", "visualization reference times are autonomous")
    _require(public.get("classification") == "pending", "public visualization time mapping is pending")

    pde_path = root / str(pde.get("source", ""))
    delivery_path = root / str(delivery.get("source", ""))
    visual_path = root / str(visual.get("source", ""))
    public_path = root / str(public.get("source", ""))
    for label, path in (
        ("PDE source", pde_path),
        ("delivery source", delivery_path),
        ("visualization source", visual_path),
        ("public mapping source", public_path),
    ):
        _require(path.is_file(), f"{label} is missing: {path}")

    pde_config = _read_json(pde_path)
    delivery_contract = _read_json(delivery_path)

    pde_interval = _ordered_interval(pde.get("time_interval"), label="PDE time interval")
    configured_pde_interval = _ordered_interval(
        pde_config.get("domain", {}).get("time_interval"),
        label="configured PDE time interval",
    )
    _require(pde_interval == configured_pde_interval, "PDE time interval drifted from preregistration")

    validation_times = _times(pde.get("validation_times"), label="PDE validation times")
    configured_validation_times = _times(
        pde_config.get("validation", {}).get("times"),
        label="configured PDE validation times",
    )
    _require(validation_times == configured_validation_times, "PDE validation times drifted from preregistration")
    _require(
        all(pde_interval[0] <= time <= pde_interval[1] for time in validation_times),
        "PDE validation time lies outside its preregistered interval",
    )
    _require(pde.get("evidence_role") == "pde_validation_only", "PDE evidence role drift")

    delivery_interval = _ordered_interval(delivery.get("time_interval"), label="delivery time interval")
    configured_delivery_interval = _ordered_interval(
        delivery_contract.get("primary_deliverable", {}).get("time_interval"),
        label="configured delivery time interval",
    )
    _require(delivery_interval == configured_delivery_interval, "delivery time interval drifted from velocity contract")
    _require(
        delivery.get("blocking_on_pde_validation") is False,
        "PDE validation must not block callable velocity delivery",
    )
    _require(delivery.get("evidence_role") == "callable_export_domain", "delivery evidence role drift")

    visual_times = _times(visual.get("times"), label="visualization reference times")
    _require(
        all(delivery_interval[0] <= time <= delivery_interval[1] for time in visual_times),
        "visualization reference time lies outside callable velocity interval",
    )
    _require(
        visual.get("evidence_role") == "visualization_diagnostic_only",
        "visualization snapshots must remain diagnostic-only evidence",
    )

    _require(public.get("status") == "pending", "public visualization time mapping is not independently verified")
    _require(public.get("mapped_times") == [], "pending public visualization mapping must not contain mapped times")

    transfer = contract.get("evidence_transfer")
    _require(isinstance(transfer, Mapping), "missing evidence_transfer policy")
    _require(
        transfer.get("visualization_snapshot_to_pde_validation") is False,
        "visualization snapshots cannot be transferred as PDE-validation evidence",
    )
    _require(
        transfer.get("pde_failure_blocks_velocity_delivery") is False,
        "PDE failure must not be converted into a callable-velocity delivery blocker",
    )
    _require(
        transfer.get("snapshot_stability_implies_public_correspondence") is False,
        "snapshot stability cannot establish public visual correspondence",
    )

    status = contract.get("claim_status")
    _require(isinstance(status, Mapping), "missing claim_status")
    for key in ("visual_correspondence_verified", "pde_validated", "paper_exact", "openai_field_identified"):
        _require(status.get(key) is False, f"unsupported claim promotion: {key}")

    return {
        "contract_pass": True,
        "pde_time_interval": pde_interval,
        "pde_validation_times": validation_times,
        "delivery_time_interval": delivery_interval,
        "visualization_reference_times": visual_times,
        "public_visual_time_mapping": "pending",
        "visual_to_pde_transfer_allowed": False,
        "pde_failure_blocks_velocity_delivery": False,
        "visual_correspondence_verified": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
    }


def audit_temporal_evidence_scope_file(
    path: str | Path, *, repo_root: str | Path
) -> dict[str, Any]:
    return audit_temporal_evidence_scope(_read_json(Path(path)), repo_root=repo_root)
