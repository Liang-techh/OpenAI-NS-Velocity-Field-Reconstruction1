"""Fail-closed CR002 audit for ST052-M temporal serialization semantics.

The live temporal module can save/load its transform *specification*.  That is a
useful reproducibility primitive, but it is not a serialized velocity candidate
because the ST052-M parent evaluator remains an externally injected callable.
This auditor keeps that distinction machine-checkable while preserving CR001
and the canonical Eq45 delivery truth states.
"""
from __future__ import annotations

import inspect
import json
from pathlib import Path
from typing import Any

from .st052_linear_temporal_transform import (
    DEFAULT_TEMPORAL_SPEC,
    St052LinearTemporalAdapter,
    TRUTH_BOUNDARY as TEMPORAL_TRUTH_BOUNDARY,
)

EXPECTED_TASK_ID = "CR002-ST052M-TEMPORAL-SERIALIZATION-SCOPE-091"
EXPECTED_TEMPORAL_SPEC_SHA256 = (
    "003ca8434ae9ec59ced410c4bd5204abe24a40ba8d5edde5d31ec132a5d4f791"
)
EXPECTED_STATIC_SPEC_SHA256 = (
    "470103b68fd5ccece5d43d054711c894e40e4cb2f7ce890d52d04c4623bf7c25"
)
EXPECTED_CANONICAL_FAMILY = "eq45_supported_velocity_candidate_v1"
EXPECTED_CANONICAL_SHA256 = (
    "2fdff812c22131d56eac1b7d6e455187ff3207500c6385aa9abf282a8e3d7b1d"
)
EXPECTED_CANONICAL_API = "openai_ns_reconstruction.eq45_supported_delivery:velocity"


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _classifications(contract: dict[str, Any]) -> dict[str, str]:
    return {
        str(entry["item"]): str(entry["classification"])
        for entry in contract["source_classification"]
    }


def audit_contract(
    contract: dict[str, Any], *, root: Path | None = None
) -> dict[str, Any]:
    """Audit one in-memory governance contract against the live repository."""
    root = _repo_root() if root is None else Path(root)
    constraints = _load_json(root / "configs" / "constraints.json")
    project_status = _load_json(root / "project_status.json")
    delivery = _load_json(root / "configs" / "velocity_delivery_contract.json")

    _require(contract.get("schema_version") == 1, "governance schema drift")
    _require(contract.get("task_id") == EXPECTED_TASK_ID, "governance task drift")

    temporal = contract["live_temporal_transform"]
    _require(
        DEFAULT_TEMPORAL_SPEC.sha256() == EXPECTED_TEMPORAL_SPEC_SHA256,
        "live temporal transform spec SHA drift",
    )
    _require(
        temporal["transform_spec_sha256"] == EXPECTED_TEMPORAL_SPEC_SHA256,
        "contract temporal transform SHA drift",
    )
    _require(
        temporal["static_transform_spec_sha256"] == EXPECTED_STATIC_SPEC_SHA256,
        "contract static transform SHA drift",
    )
    _require(
        temporal["transform_id"] == DEFAULT_TEMPORAL_SPEC.transform_id,
        "temporal transform identity drift",
    )
    _require(
        temporal["parent_candidate_id"] == DEFAULT_TEMPORAL_SPEC.parent_candidate_id,
        "temporal parent candidate metadata drift",
    )
    _require(
        temporal["parent_source_head"] == DEFAULT_TEMPORAL_SPEC.parent_source_head,
        "temporal parent source metadata drift",
    )
    _require(temporal["spec_save_load_ready"] is True, "spec save/load regressed")
    _require(
        temporal["spec_roundtrip_scope"] == "temporal_transform_metadata_only",
        "spec round-trip scope laundering",
    )
    for key in (
        "base_velocity_serialized",
        "parent_artifact_sha_bound",
        "composite_candidate_sha_exists",
        "complete_candidate_save_load_ready",
        "experimental_velocity_export_ready",
    ):
        _require(temporal[key] is False, f"premature temporal delivery promotion: {key}")

    semantics = contract["identity_semantics"]
    for key in (
        "transform_spec_sha_is_complete_candidate_sha",
        "spec_save_load_is_complete_candidate_save_load",
        "parent_metadata_strings_authenticate_parent_callable",
        "static_transform_binding_authenticates_st052m_parent",
        "endpoint_equivalence_defines_interior_time_field",
        "complete_temporal_child_materialized",
    ):
        _require(semantics[key] is False, f"identity laundering: {key}")
    _require(
        semantics["temporal_transform_kernel_materialized"] is True,
        "integrated temporal kernel must remain recorded as materialized",
    )

    # Runtime evidence for the serialization boundary: the saved payload is only
    # the transform spec, while constructing the evaluator still requires a
    # separately supplied parent callable.
    payload_keys = set(DEFAULT_TEMPORAL_SPEC.payload())
    _require("parent_candidate_sha256" not in payload_keys, "unexpected parent SHA binding")
    _require("complete_candidate_sha256" not in payload_keys, "unexpected child SHA binding")
    signature = inspect.signature(St052LinearTemporalAdapter.__init__)
    base_velocity = signature.parameters.get("base_velocity")
    _require(base_velocity is not None, "temporal adapter no longer exposes parent callable")
    _require(
        base_velocity.default is inspect.Parameter.empty,
        "temporal adapter parent callable unexpectedly became optional",
    )
    _require(
        TEMPORAL_TRUTH_BOUNDARY["full_st052_parent_materialized_here"] is False,
        "temporal module truth boundary promoted parent materialization",
    )
    _require(
        TEMPORAL_TRUTH_BOUNDARY["complete_candidate_save_load_ready"] is False,
        "temporal module truth boundary promoted candidate save/load",
    )

    prerequisites = set(contract["requirements_before_temporal_child_export_ready"])
    required_prerequisites = {
        "materialize deterministic ST052-M parent evaluator",
        "assign immutable parent artifact or candidate SHA and bind it at runtime",
        "compose parent identity, static transform identity, and temporal transform identity into a versioned child manifest",
        "assign a complete child candidate ID and SHA",
        "save and reload the complete child so the same parent plus transforms are reconstructed without an externally supplied unauthenticated callable",
        "replay the unified velocity(x,y,z,t) API after reload and verify finite pointwise equality",
        "keep visual correspondence and PDE validation as independent claim states",
    }
    _require(
        required_prerequisites.issubset(prerequisites),
        "complete-child delivery prerequisite removed",
    )

    canonical = contract["canonical_delivery"]
    _require(
        canonical["candidate_family"] == EXPECTED_CANONICAL_FAMILY,
        "canonical family relabelled",
    )
    _require(
        canonical["candidate_sha256"] == EXPECTED_CANONICAL_SHA256,
        "canonical candidate SHA relabelled",
    )
    _require(canonical["velocity_api"] == EXPECTED_CANONICAL_API, "canonical API relabelled")
    _require(canonical["velocity_export_ready"] is True, "canonical export readiness regressed")
    for key in (
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        _require(canonical[key] is False, f"canonical truth state promoted: {key}")

    states = project_status["states"]
    _require(project_status["candidate_family"] == EXPECTED_CANONICAL_FAMILY, "project family drift")
    _require(project_status["candidate_sha256"] == EXPECTED_CANONICAL_SHA256, "project SHA drift")
    _require(project_status["velocity_api"] == EXPECTED_CANONICAL_API, "project API drift")
    _require(states["velocity_export_ready"] is True, "project canonical export readiness drift")
    _require(
        project_status["latest_st052m_visual_candidate_evidence"]["velocity_export_ready"]
        is False,
        "experimental ST052-M export readiness was prematurely promoted",
    )
    for key in (
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        _require(states[key] is False, f"project truth state promoted: {key}")

    primary = delivery["primary_deliverable"]
    _require(primary["candidate_family"] == EXPECTED_CANONICAL_FAMILY, "delivery family drift")
    _require(primary["candidate_sha256"] == EXPECTED_CANONICAL_SHA256, "delivery SHA drift")
    _require(primary["api"] == EXPECTED_CANONICAL_API, "delivery API drift")
    _require(delivery["claim_status"]["velocity_export_ready"] is True, "delivery readiness drift")
    for key in (
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        _require(delivery["claim_status"][key] is False, f"delivery truth state promoted: {key}")

    frozen = contract["cr001_nonmutation"]
    _require(frozen["nu"] == constraints["nu"] == 0.01, "CR001 nu drift")
    _require(frozen["evaluation_box"] == constraints["domain"]["evaluation_box"], "CR001 box drift")
    _require(frozen["support"] == constraints["domain"]["support"], "CR001 support drift")
    _require(frozen["time_interval"] == constraints["domain"]["time_interval"], "CR001 time drift")
    _require(frozen["force_mode"] == constraints["forcing"]["mode"], "CR001 force mode drift")
    _require(
        frozen["force_parameter_bounds"] == constraints["forcing"]["parameters"],
        "CR001 force bounds drift",
    )
    _require(frozen["reference_energy"] == constraints["nontriviality"]["reference_energy"], "CR001 energy drift")
    _require(
        frozen["reference_energy_abs_tolerance"]
        == constraints["nontriviality"]["reference_energy_abs_tolerance"]
        == 0.001,
        "CR001 energy tolerance drift",
    )
    validation = constraints["validation"]
    _require(frozen["validation_seed"] == validation["seed"] == 914027, "CR001 validation seed drift")
    _require(frozen["held_out_points"] == validation["held_out_points"] == 4096, "CR001 sample drift")
    _require(frozen["validation_times"] == validation["times"], "CR001 validation-time drift")
    _require(frozen["derivative_steps"] == validation["derivative_steps"], "CR001 FD ladder drift")
    _require(
        frozen["quadrature_orders_per_axis"] == validation["quadrature_orders_per_axis"],
        "CR001 quadrature ladder drift",
    )
    thresholds = validation["thresholds"]
    for key in ("divergence_max", "divergence_L2", "pde_residual_max", "pde_residual_L2"):
        _require(frozen[key] == thresholds[key], f"CR001 threshold drift: {key}")

    classes = _classifications(contract)
    expected_classes = {
        "callable and saveable/loadable velocity(x,y,z,t)->[u,v,w] delivery": "user_requirement",
        "qualitative observable structure of the published OpenAI velocity-field visualization": "public_source_fact",
        "ST052-M local-swirl transform and linear temporal activation, including serialization layout": "autonomous_design",
        "exact OpenAI numerical field, hidden coefficients, exact frame/time map, camera and trajectory seeds": "pending_unknown",
    }
    _require(classes == expected_classes, "source classification drift")

    return {
        "status": "pass",
        "temporal_transform_spec_sha256": DEFAULT_TEMPORAL_SPEC.sha256(),
        "transform_spec_save_load_ready": True,
        "complete_candidate_save_load_ready": False,
        "experimental_velocity_export_ready": False,
        "canonical_velocity_export_ready": True,
        "pde_validated": False,
    }


def run_audit(
    contract_path: str | Path | None = None, *, root: Path | None = None
) -> dict[str, Any]:
    root = _repo_root() if root is None else Path(root)
    path = (
        root / "configs" / "st052m_temporal_serialization_governance_contract.json"
        if contract_path is None
        else Path(contract_path)
    )
    return audit_contract(_load_json(path), root=root)


if __name__ == "__main__":
    print(json.dumps(run_audit(), indent=2, sort_keys=True))
