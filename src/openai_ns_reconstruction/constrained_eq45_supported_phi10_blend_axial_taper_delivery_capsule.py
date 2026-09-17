"""Replayable delivery capsule for the supported Eq45 blend/axial-taper family.

The caller declares both the compact--quartic blend weight and the bounded
axial taper plateau.  This module packages that exact public ``[u,v,w]``
identity for downstream validators and visualization code.  It deliberately
selects neither parameter and does not promote visualization or PDE status.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .constrained_eq45_supported_collar_vorticity import governed_supported_seed
from .constrained_eq45_supported_phi10_blend_axial_taper_mode import (
    AXIAL_PLATEAU_Q_BOUNDS,
    BASE_AXIAL_PLATEAU_Q,
    MAX_AXIAL_PLATEAU_Q,
    Eq45SupportedPhi10BlendAxialTaperCandidate,
)
from .constrained_eq45_supported_phi10_compact_quartic_blend_temporal_mode import (
    Eq45SupportedPhi10CompactQuarticBlendTemporalCandidate,
)

SCHEMA = "eq45_supported_phi10_blend_axial_taper_delivery_recipe_v1"
CAPSULE_SCHEMA = "eq45_supported_phi10_blend_axial_taper_delivery_capsule_v1"
TASK_ID = "CR011-EQ45-BLEND-AXIAL-TAPER-DELIVERY-CAPSULE-019"
DEPENDENCY_HEAD = "e50a536aab3f8c65215ee28290571ebc1ed45f77"
EXPECTED_SUPPORTED_SHA256 = (
    "2fdff812c22131d56eac1b7d6e455187ff3207500c6385aa9abf282a8e3d7b1d"
)
EXPECTED_REFERENCE_TIMES = (0.25, 0.3125, 0.5, 0.6875, 0.75)
EXPECTED_EARLY_DELTA = -1.4

_RECIPE_PATH = (
    Path(__file__).resolve().parents[2]
    / "artifacts/constrained/eq45_supported_phi10_blend_axial_taper_delivery_recipe.json"
)

_EXPECTED_STATUS = {
    "axial_taper_control_materialized": True,
    "axial_taper_value_selected": False,
    "blend_weight_selected": False,
    "blowup_proved": False,
    "callable_serializable": True,
    "caller_declared_axial_taper": True,
    "caller_declared_blend_weight": True,
    "candidate_selection_resolved": False,
    "canonical_velocity_changed": False,
    "force_or_pressure_fitted": False,
    "openai_field_identified": False,
    "paper_exact": False,
    "pde_validated": False,
    "physical_support_connection_implemented": True,
    "physical_support_validated": False,
    "public_image_fitted": False,
    "velocity_export_ready": True,
    "visual_correspondence_verified": False,
    "visualization_candidate_only": True,
    "visualization_ready": False,
}


def _canonical_sha256(payload: Mapping[str, Any]) -> str:
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _validate_recipe(payload: Mapping[str, Any]) -> dict[str, Any]:
    if not isinstance(payload, Mapping):
        raise ValueError("axial-taper delivery recipe must be an object")
    data = dict(payload)
    if data.get("schema") != SCHEMA:
        raise ValueError(f"recipe schema must be {SCHEMA!r}")
    if data.get("task_id") != TASK_ID:
        raise ValueError("delivery task identity drifted")
    if data.get("base_dependency_head") != DEPENDENCY_HEAD:
        raise ValueError("delivery dependency head drifted")
    if data.get("base_supported_sha256") != EXPECTED_SUPPORTED_SHA256:
        raise ValueError("supported parent identity drifted")

    representation = data.get("representation")
    if not isinstance(representation, Mapping):
        raise ValueError("representation must be an object")
    if representation.get("family") != (
        "Phi(1,0) compact-quartic blend with autonomous axial support taper"
    ):
        raise ValueError("representation family drifted")
    if representation.get("blend_weight_bounds") != [0.0, 1.0]:
        raise ValueError("blend-weight bounds drifted")
    if representation.get("blend_weight_selected") is not False:
        raise ValueError("recipe must not select a blend weight")
    if representation.get("axial_taper_parameter") != (
        "AxisymmetricPhysicalTaper.axial_plateau_q"
    ):
        raise ValueError("axial-taper parameter identity drifted")
    if representation.get("axial_plateau_q_bounds") != list(AXIAL_PLATEAU_Q_BOUNDS):
        raise ValueError("axial-taper bounds drifted")
    if representation.get("axial_plateau_q_selected") is not False:
        raise ValueError("recipe must not select an axial taper")
    if representation.get("axial_identity_half_height_formula") != (
        "2*sqrt(axial_plateau_q)"
    ):
        raise ValueError("axial identity-height formula drifted")
    if float(representation.get("axial_support_half_height")) != 2.0:
        raise ValueError("axial support half-height drifted")
    if float(representation.get("radial_support")) != 2.0:
        raise ValueError("radial support drifted")
    if float(representation.get("radial_plateau_q")) != 0.64:
        raise ValueError("radial taper must remain unchanged")

    delivery = data.get("delivery")
    if not isinstance(delivery, Mapping):
        raise ValueError("delivery must be an object")
    if delivery.get("spatial_box") != [[-2.0, 2.0], [-2.0, 2.0], [-2.0, 2.0]]:
        raise ValueError("delivery spatial box drifted")
    if delivery.get("time_interval") != [0.25, 0.75]:
        raise ValueError("delivery time interval drifted")
    if tuple(float(x) for x in delivery.get("reference_times", [])) != EXPECTED_REFERENCE_TIMES:
        raise ValueError("delivery reference times drifted")
    if set(delivery.get("velocity_interfaces", [])) != {
        "velocity",
        "at_points",
        "velocity_xyz",
        "grid",
    }:
        raise ValueError("public velocity interface inventory drifted")
    if set(delivery.get("serialization_interfaces", [])) != {"save_json", "load_json"}:
        raise ValueError("serialization interface inventory drifted")
    if delivery.get("component_order") != ["u", "v", "w"]:
        raise ValueError("Cartesian component order drifted")
    if delivery.get("grid_layout") != "time,x,y,z,component":
        raise ValueError("grid layout drifted")

    siblings = data.get("sibling_evidence")
    if not isinstance(siblings, Mapping):
        raise ValueError("sibling_evidence must be an object")
    expected_siblings = {
        "restricted_force_pr": 195,
        "restricted_force_status": "open_unconsumed",
        "divergence_pr": 196,
        "divergence_status": "open_unconsumed",
        "morphology_pr": 197,
        "morphology_status": "open_unconsumed_calibration_pending",
        "numeric_summary_bound_here": False,
        "may_select_blend_weight": False,
        "may_select_axial_taper": False,
        "may_promote_visualization_ready": False,
        "may_promote_pde_validated": False,
    }
    if dict(siblings) != expected_siblings:
        raise ValueError("sibling evidence provenance or consumption state drifted")

    if data.get("status") != _EXPECTED_STATUS:
        raise ValueError("delivery truth-boundary status drifted")

    pending = data.get("pending")
    if not isinstance(pending, Mapping):
        raise ValueError("pending must be an object")
    for key in (
        "axial_taper_selection",
        "blend_weight_selection",
        "candidate_specific_energy_acceptance",
        "formal_pde_gate",
        "public_observable_comparison",
        "three_dimensional_render_review",
    ):
        if pending.get(key) is not True:
            raise ValueError(f"{key} must remain pending")
    return data


def load_recipe(path: str | Path | None = None) -> dict[str, Any]:
    """Load and fail-closed validate the committed family recipe."""
    recipe_path = Path(path) if path is not None else _RECIPE_PATH
    return _validate_recipe(json.loads(recipe_path.read_text(encoding="utf-8")))


def build_candidate(
    blend_weight: float,
    axial_plateau_q: float,
    recipe: Mapping[str, Any] | None = None,
) -> Eq45SupportedPhi10BlendAxialTaperCandidate:
    """Build one caller-declared family member without selecting either control."""
    payload = load_recipe() if recipe is None else _validate_recipe(recipe)
    weight = float(blend_weight)
    q_value = float(axial_plateau_q)
    if not np.isfinite(weight) or not (0.0 <= weight <= 1.0):
        raise ValueError("blend_weight must be finite and lie in [0,1]")
    if not np.isfinite(q_value) or not (
        BASE_AXIAL_PLATEAU_Q <= q_value <= MAX_AXIAL_PLATEAU_Q
    ):
        raise ValueError(
            "axial_plateau_q must be finite and lie in "
            f"[{BASE_AXIAL_PLATEAU_Q},{MAX_AXIAL_PLATEAU_Q}]"
        )

    supported = governed_supported_seed()
    if supported.sha256 != EXPECTED_SUPPORTED_SHA256:
        raise RuntimeError("governed supported parent identity drifted")
    if supported.sha256 != payload["base_supported_sha256"]:
        raise RuntimeError("recipe no longer binds the governed supported parent")

    blend = Eq45SupportedPhi10CompactQuarticBlendTemporalCandidate(
        base=supported,
        blend_weight=weight,
        early_delta=EXPECTED_EARLY_DELTA,
    )
    candidate = Eq45SupportedPhi10BlendAxialTaperCandidate(
        base=blend,
        axial_plateau_q=q_value,
    )
    truth = candidate.to_dict()["truth_boundary"]
    if truth.get("velocity_export_ready") is not True:
        raise RuntimeError("reconstructed candidate lost velocity export readiness")
    for key in (
        "physical_support_validated",
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        if truth.get(key) is not False:
            raise RuntimeError(f"reconstructed candidate over-promoted {key}")
    return candidate


def delivery_capsule(
    blend_weight: float,
    axial_plateau_q: float,
    recipe_path: str | Path | None = None,
) -> dict[str, Any]:
    """Return checked family provenance plus one explicit member identity."""
    recipe = load_recipe(recipe_path)
    candidate = build_candidate(blend_weight, axial_plateau_q, recipe)
    return {
        "schema": CAPSULE_SCHEMA,
        "task_id": TASK_ID,
        "recipe_sha256": _canonical_sha256(recipe),
        "candidate_sha256": candidate.sha256,
        "base_blend_sha256": candidate.base_sha256,
        "base_supported_sha256": candidate.supported_base_sha256,
        "declared_blend_weight": float(candidate.base.blend_weight),
        "declared_axial_plateau_q": float(candidate.axial_plateau_q),
        "axial_identity_half_height": float(candidate.axial_identity_half_height),
        "blend_weight_selected_by_capsule": False,
        "axial_taper_selected_by_capsule": False,
        "delivery": dict(recipe["delivery"]),
        "representation": dict(recipe["representation"]),
        "sibling_evidence": dict(recipe["sibling_evidence"]),
        "status": dict(recipe["status"]),
        "pending": dict(recipe["pending"]),
    }


def write_delivery_bundle(
    output_dir: str | Path,
    blend_weight: float,
    axial_plateau_q: float,
    recipe_path: str | Path | None = None,
) -> dict[str, Path]:
    """Write ``candidate.json`` + ``capsule.json`` and prove exact SHA replay."""
    recipe = load_recipe(recipe_path)
    candidate = build_candidate(blend_weight, axial_plateau_q, recipe)
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)

    candidate_path = output / "candidate.json"
    capsule_path = output / "capsule.json"
    candidate.save_json(candidate_path)
    reloaded = Eq45SupportedPhi10BlendAxialTaperCandidate.load_json(candidate_path)
    if reloaded.sha256 != candidate.sha256:
        raise RuntimeError("saved axial-taper candidate failed exact SHA replay")

    capsule = delivery_capsule(blend_weight, axial_plateau_q, recipe_path)
    capsule["candidate_file"] = candidate_path.name
    capsule["candidate_file_sha256"] = hashlib.sha256(candidate_path.read_bytes()).hexdigest()
    capsule_path.write_text(
        json.dumps(capsule, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return {"candidate": candidate_path, "capsule": capsule_path}
