"""CR012 acceptance-state snapshot for the support-connected Eq. (4.5) child.

This module is deliberately narrower than the historical coupled-field acceptance
ledger.  It binds the *current supported child identity* to the post-support
energy evidence already present in this branch and keeps delivery, visual, and
PDE states orthogonal.

Sibling PR evidence is not copied into this branch: until such evidence is in
this ancestry (or otherwise consumed by a checked receipt), those gates remain
pending.  In particular, an exportable field may coexist with a failed CR001
energy-normalization gate and with ``pde_validated == False``.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping

from .constrained_eq45_candidate import Eq45VelocityCandidate
from .constrained_eq45_supported_candidate import Eq45SupportedVelocityCandidate

SCHEMA = "eq45_supported_acceptance_state_v1"
TASK_ID = "CR012-EQ45-SUPPORTED-ACCEPTANCE-STATE-011"
EXPECTED_PARENT_SHA256 = "48f1845fcfd71ec95495c98bb7aac3fca4653e748a8856a5421234e9601525a7"
EXPECTED_CHILD_SHA256 = "2fdff812c22131d56eac1b7d6e455187ff3207500c6385aa9abf282a8e3d7b1d"


def _require_close(name: str, left: float, right: float, tol: float = 1e-14) -> None:
    if abs(float(left) - float(right)) > tol:
        raise ValueError(f"{name} drifted: {left!r} != {right!r}")


def _evidence(
    evidence_id: str,
    status: str,
    evidence_class: str,
    source: str,
    detail: str,
    *,
    observed: Any = None,
    threshold: Any = None,
) -> dict[str, Any]:
    if status not in {"pass", "fail", "pending"}:
        raise ValueError(f"invalid status for {evidence_id}: {status!r}")
    return {
        "id": evidence_id,
        "status": status,
        "evidence_class": evidence_class,
        "source": source,
        "detail": detail,
        "observed": observed,
        "threshold": threshold,
    }


def build_supported_acceptance_state(
    child: Eq45SupportedVelocityCandidate,
    constraints: Mapping[str, Any],
    energy_audit: Mapping[str, Any],
    delivery_contract: Mapping[str, Any],
) -> dict[str, Any]:
    """Build the checked CR012 state for the current supported Eq45 child."""
    if not isinstance(child, Eq45SupportedVelocityCandidate):
        raise TypeError("child must be Eq45SupportedVelocityCandidate")
    if constraints.get("experiment_id") != "compact_axisymmetric_window_v1":
        raise ValueError("unexpected CR001 experiment identity")

    if child.parent_sha256 != EXPECTED_PARENT_SHA256:
        raise ValueError("supported child parent identity drifted")
    if child.sha256 != EXPECTED_CHILD_SHA256:
        raise ValueError("supported child identity drifted")
    if energy_audit.get("candidate_sha256") != child.sha256:
        raise ValueError("energy audit is not bound to this supported child")
    if energy_audit.get("parent_sha256") != child.parent_sha256:
        raise ValueError("energy audit parent identity drifted")

    nontriviality = constraints["nontriviality"]
    reference_time = float(nontriviality["reference_time"])
    target = float(nontriviality["reference_energy"])
    tolerance = float(nontriviality["reference_energy_abs_tolerance"])
    observed = float(energy_audit["reference_observed_finest_energy"])
    error = abs(observed - target)
    energy_pass = bool(error <= tolerance)

    _require_close("reference time", energy_audit["reference_time"], reference_time)
    _require_close("reference target", energy_audit["reference_energy"], target)
    _require_close(
        "reference tolerance",
        energy_audit["reference_energy_abs_tolerance"],
        tolerance,
    )
    _require_close(
        "reported reference-energy error",
        energy_audit["reference_energy_abs_error"],
        error,
    )
    if bool(energy_audit["reference_energy_gate_pass"]) != energy_pass:
        raise ValueError("energy audit reference gate is inconsistent with CR001")

    thresholds = constraints["validation"]["thresholds"]
    _require_close(
        "energy quadrature threshold",
        energy_audit["quadrature_relative_change_threshold"],
        float(thresholds["energy_quadrature_relative_change"]),
    )

    truth = child.to_dict()["truth_boundary"]
    required_truth = {
        "callable_serializable": True,
        "velocity_export_ready": True,
        "physical_support_connection_implemented": True,
        "physical_support_validated": False,
        "visualization_ready": False,
        "visual_correspondence_verified": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
    }
    for key, expected in required_truth.items():
        if truth.get(key) is not expected:
            raise ValueError(f"supported-child truth boundary drifted at {key}")

    state_definitions = delivery_contract.get("states")
    if not isinstance(state_definitions, Mapping):
        raise ValueError("delivery state contract is missing states")
    for key in (
        "velocity_export_ready",
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
    ):
        if key not in state_definitions:
            raise ValueError(f"delivery state contract is missing {key}")

    forbidden = delivery_contract.get("forbidden_inferences", [])
    if not any(
        item.get("from") == ["pde_validation_failed"]
        and item.get("to") == "velocity_export_not_allowed"
        for item in forbidden
        if isinstance(item, Mapping)
    ):
        raise ValueError("delivery contract lost the PDE-failure/export separation")

    quadrature_pass = bool(energy_audit["all_quadrature_gates_pass"])
    range_pass = bool(energy_audit["all_energy_range_gates_pass"])
    evidence = [
        _evidence(
            "candidate_identity",
            "pass",
            "strict_artifact_identity",
            "Eq45SupportedVelocityCandidate.sha256",
            "exact supported-child and parent identities are pinned",
            observed=child.sha256,
            threshold=EXPECTED_CHILD_SHA256,
        ),
        _evidence(
            "velocity_export_ready",
            "pass",
            "representation_delivery_evidence",
            "Eq45SupportedVelocityCandidate truth boundary",
            "candidate is callable, serializable, griddable, and explicitly export-ready",
            observed=True,
            threshold=True,
        ),
        _evidence(
            "physical_support_connection_implemented",
            "pass",
            "strict_representation_property",
            "Eq45SupportedVelocityCandidate truth boundary",
            "the governed physical-space support transform is part of this child identity",
            observed=True,
            threshold=True,
        ),
        _evidence(
            "cr001_reference_energy",
            "pass" if energy_pass else "fail",
            "independent_numeric_evidence",
            "artifacts/constrained/eq45_supported_energy_audit.json",
            "post-support kinetic energy at the preregistered reference time",
            observed=observed,
            threshold={"target": target, "abs_tolerance": tolerance},
        ),
        _evidence(
            "energy_quadrature_convergence",
            "pass" if quadrature_pass else "fail",
            "independent_numeric_evidence",
            "artifacts/constrained/eq45_supported_energy_audit.json",
            "all registered energy quadrature convergence gates",
            observed=quadrature_pass,
            threshold=True,
        ),
        _evidence(
            "validation_time_energy_range",
            "pass" if range_pass else "fail",
            "independent_numeric_evidence",
            "artifacts/constrained/eq45_supported_energy_audit.json",
            "all checked times remain within the preregistered nontriviality range",
            observed=range_pass,
            threshold=True,
        ),
        _evidence(
            "standalone_committed_child_artifact",
            "pending",
            "unconsumed_sibling_evidence",
            "PR #123",
            "persistent child JSON exists on a sibling PR but is not copied into this ancestry",
        ),
        _evidence(
            "supported_child_divergence_validation",
            "pending",
            "unconsumed_sibling_evidence",
            "PR #121",
            "three-level child divergence evidence remains outside this branch ancestry",
        ),
        _evidence(
            "physical_support_validated",
            "pending",
            "unresolved",
            "independent child support/boundary acceptance",
            "representation-level exact support does not by itself close the validation state",
        ),
        _evidence(
            "supported_child_pde_validation",
            "pending",
            "unconsumed_sibling_evidence",
            "PR #129",
            "independent supported-child PDE/vorticity evidence is not inherited across sibling PRs",
        ),
        _evidence(
            "visualization_ready",
            "pending",
            "unresolved",
            "supported visualization entry point + resolution sanity evidence",
            "callable velocity alone does not establish a supported reproducible visualization path",
        ),
        _evidence(
            "visual_correspondence_verified",
            "pending",
            "unresolved",
            "public-reference visual comparison",
            "no public-observable morphology/time-evolution acceptance is consumed here",
        ),
    ]

    acceptance_ready = all(item["status"] == "pass" for item in evidence)
    return {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "candidate_sha256": child.sha256,
        "parent_sha256": child.parent_sha256,
        "candidate_scope": "support_connected_eq45_visualization_candidate",
        "cr001_experiment_id": constraints["experiment_id"],
        "reference_time": reference_time,
        "evidence": evidence,
        "states": {
            "velocity_export_ready": True,
            "visualization_ready": False,
            "visual_correspondence_verified": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
        "acceptance_ready": acceptance_ready,
        "interpretation": {
            "energy_failure_blocks_velocity_export": False,
            "energy_failure_blocks_current_cr001_acceptance": not energy_pass,
            "sibling_pr_evidence_promoted_without_consumption": False,
            "visual_similarity_would_not_imply_pde_validity": True,
            "pde_validity_would_not_imply_visual_correspondence": True,
        },
    }


def audit_repository(root: str | Path) -> dict[str, Any]:
    root = Path(root)
    parent = Eq45VelocityCandidate.load_json(
        root / "artifacts" / "constrained" / "eq45_velocity_candidate_seed.json"
    )
    child = Eq45SupportedVelocityCandidate(parent=parent)
    load = lambda path: json.loads((root / path).read_text(encoding="utf-8"))
    return build_supported_acceptance_state(
        child,
        load("configs/constraints.json"),
        load("artifacts/constrained/eq45_supported_energy_audit.json"),
        load("configs/delivery_state_contract.json"),
    )
