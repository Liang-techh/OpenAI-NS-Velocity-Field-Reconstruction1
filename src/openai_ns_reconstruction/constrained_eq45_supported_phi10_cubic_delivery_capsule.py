"""Replayable delivery capsule for the cubic-localized supported Phi(1,0) field.

This packages an already-materialized velocity candidate and already-completed
public-output temporal-derivative audit.  It does not alter the velocity field,
fit pressure or forcing, select a public-image target, or promote scientific
readiness.  Export, visualization acceptance, and PDE acceptance remain
independent states.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .constrained_eq45_candidate import Eq45VelocityCandidate
from .constrained_eq45_supported_candidate import Eq45SupportedVelocityCandidate
from .constrained_eq45_supported_phi10_cubic_temporal_mode import (
    Eq45SupportedPhi10CubicLocalizedTemporalCandidate,
)

SCHEMA = "eq45_supported_phi10_cubic_delivery_recipe_v1"
CAPSULE_SCHEMA = "eq45_supported_phi10_cubic_delivery_capsule_v1"
TASK_ID = "CR011-EQ45-PHI10-CUBIC-DELIVERY-CAPSULE-015"
DEPENDENCY_HEAD = "b5c9d0cb3145e014911bffd71420148c2e16c4f4"
EXPECTED_BASE_SHA256 = "2fdff812c22131d56eac1b7d6e455187ff3207500c6385aa9abf282a8e3d7b1d"
EXPECTED_CANDIDATE_SHA256 = "310fc2ad1bf3d721860369b9a58c7599202182964158708867691f4445f21576"
EXPECTED_EARLY_DELTA = -1.4
EXPECTED_TIMES = (0.25, 0.375, 0.5, 0.625, 0.6875, 0.75)
EXPECTED_COEFFICIENTS = (-1.7, -0.65, -0.3, -0.3, -0.321875, -0.3)

_RECIPE_PATH = (
    Path(__file__).resolve().parents[2]
    / "artifacts/constrained/eq45_supported_phi10_cubic_recipe.json"
)

_EXPECTED_STATUS = {
    "callable_serializable": True,
    "velocity_export_ready": True,
    "visualization_candidate_only": True,
    "canonical_velocity_changed": False,
    "screened_temporal_shape_materialized": True,
    "production_temporal_shape_promoted": False,
    "physical_support_connection_implemented": True,
    "physical_support_validated": False,
    "spatial_basis_grown": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _canonical_sha256(payload: Mapping[str, Any]) -> str:
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _validate_recipe(payload: Mapping[str, Any]) -> dict[str, Any]:
    if not isinstance(payload, Mapping):
        raise ValueError("cubic delivery recipe must be an object")
    data = dict(payload)
    if data.get("schema") != SCHEMA:
        raise ValueError(f"recipe schema must be {SCHEMA!r}")
    if data.get("task_id") != TASK_ID:
        raise ValueError("delivery task identity drifted")
    if data.get("base_dependency_head") != DEPENDENCY_HEAD:
        raise ValueError("delivery dependency head drifted")
    if data.get("base_supported_sha256") != EXPECTED_BASE_SHA256:
        raise ValueError("supported parent identity drifted")
    if data.get("candidate_sha256") != EXPECTED_CANDIDATE_SHA256:
        raise ValueError("cubic candidate identity drifted")

    mode = data.get("temporal_mode")
    if not isinstance(mode, Mapping):
        raise ValueError("temporal_mode must be an object")
    if mode.get("family") != "phi" or mode.get("index") != [1, 0]:
        raise ValueError("temporal mode must remain exactly Phi(1,0)")
    if float(mode.get("early_delta")) != EXPECTED_EARLY_DELTA:
        raise ValueError("early_delta drifted from the screened -1.4 value")
    if mode.get("schedule") != "early_delta*tau*(tau-0.5)*(tau-1)/(-3)":
        raise ValueError("cubic temporal schedule drifted")
    if mode.get("coefficient_bound") != [-4.0, 4.0]:
        raise ValueError("profile coefficient bound drifted")
    if tuple(float(x) for x in mode.get("reference_times", [])) != EXPECTED_TIMES:
        raise ValueError("temporal reference times drifted")
    if tuple(float(x) for x in mode.get("reference_coefficients", [])) != EXPECTED_COEFFICIENTS:
        raise ValueError("registered cubic coefficients drifted")

    delivery = data.get("delivery")
    if not isinstance(delivery, Mapping):
        raise ValueError("delivery must be an object")
    if delivery.get("time_interval") != [0.25, 0.75]:
        raise ValueError("delivery time interval drifted")
    if tuple(delivery.get("reference_times", [])) != EXPECTED_TIMES:
        raise ValueError("delivery reference times drifted")
    if delivery.get("spatial_box") != [[-2.0, 2.0], [-2.0, 2.0], [-2.0, 2.0]]:
        raise ValueError("delivery spatial box drifted")
    if set(delivery.get("velocity_interfaces", [])) != {
        "velocity", "at_points", "velocity_xyz", "grid"
    }:
        raise ValueError("public velocity interface inventory drifted")
    if set(delivery.get("serialization_interfaces", [])) != {"save_json", "load_json"}:
        raise ValueError("serialization interface inventory drifted")
    if delivery.get("grid_layout") != "time,x,y,z,component":
        raise ValueError("grid layout drifted")
    if delivery.get("component_order") != ["u", "v", "w"]:
        raise ValueError("Cartesian component order drifted")

    representation = data.get("representation_evidence")
    if not isinstance(representation, Mapping):
        raise ValueError("representation_evidence must be an object")
    if representation.get("capacity_pr") != 158 or representation.get("materialization_pr") != 159:
        raise ValueError("cubic representation provenance drifted")
    if representation.get("static_return_times") != [0.5, 0.625, 0.75]:
        raise ValueError("static-return keyframes drifted")
    if representation.get("off_keyframe_nontrivial_times") != [0.375, 0.6875]:
        raise ValueError("off-keyframe inventory drifted")
    if representation.get("public_openai_reference_used_for_selection") is not False:
        raise ValueError("target-free cubic selection cannot be relabeled as public-image fitting")
    if representation.get("visual_correspondence_established") is not False:
        raise ValueError("capacity evidence cannot establish public correspondence")
    cubic_excursion = float(representation.get("cubic_late_half_max_coefficient_excursion"))
    quadratic_excursion = float(representation.get("quadratic_late_half_max_coefficient_excursion"))
    if not (0.0 < cubic_excursion < quadratic_excursion):
        raise ValueError("registered cubic late-half localization advantage drifted")

    temporal = data.get("temporal_derivative_evidence")
    if not isinstance(temporal, Mapping):
        raise ValueError("temporal_derivative_evidence must be an object")
    if temporal.get("audit_pr") != 161 or temporal.get("audit_exact_head") != DEPENDENCY_HEAD:
        raise ValueError("temporal-derivative audit provenance drifted")
    if temporal.get("probe_seed") != 914509 or temporal.get("probe_count") != 24:
        raise ValueError("temporal-derivative held-out probe contract drifted")
    if tuple(temporal.get("audit_times", [])) != EXPECTED_TIMES:
        raise ValueError("temporal-derivative audit times drifted")
    if temporal.get("second_order_steps") != [0.01, 0.005, 0.0025]:
        raise ValueError("time-derivative refinement ladder drifted")
    if float(temporal.get("fourth_order_reference_step")) != 0.00125:
        raise ValueError("time-derivative reference step drifted")
    orders = tuple(float(x) for x in temporal.get("observed_rms_orders", []))
    if len(orders) != 2 or min(orders) <= 1.9:
        raise ValueError("registered public-u_t second-order convergence evidence drifted")
    normalized = tuple(float(x) for x in temporal.get("normalized_rms_errors", []))
    if len(normalized) != 3 or not (normalized[2] < normalized[1] < normalized[0]):
        raise ValueError("public-u_t normalized errors no longer refine monotonically")
    if temporal.get("snapshot_identity_not_derivative_identity") is not True:
        raise ValueError("snapshot/derivative distinction must remain explicit")
    if temporal.get("formal_pde_gate_assessed") is not False or temporal.get("pde_validated") is not False:
        raise ValueError("temporal convergence cannot be promoted to formal PDE validation")
    if float(temporal.get("sign_mutation_bad_rms_error")) <= 1.0:
        raise ValueError("time-derivative sign mutation calibration drifted")

    sibling = data.get("sibling_evidence")
    if not isinstance(sibling, Mapping):
        raise ValueError("sibling_evidence must be an object")
    if sibling.get("restricted_force_crosscheck_pr") != 160:
        raise ValueError("restricted-force sibling provenance drifted")
    if sibling.get("ancestry_status") != "open_unconsumed_sibling_evidence":
        raise ValueError("concurrent force/PDE sibling must remain unconsumed")
    if sibling.get("numeric_summary_bound_here") is not False:
        raise ValueError("unconsumed sibling numerics cannot be laundered into this capsule")
    if sibling.get("may_promote_pde_validated") is not False:
        raise ValueError("sibling evidence cannot promote PDE validity here")

    if data.get("status") != _EXPECTED_STATUS:
        raise ValueError("delivery truth-boundary status drifted")
    pending = data.get("pending")
    if not isinstance(pending, Mapping):
        raise ValueError("pending must be an object")
    for key in (
        "public_observable_comparison",
        "three_dimensional_render_review",
        "formal_pde_gate",
        "candidate_specific_energy_acceptance",
        "restricted_force_sibling_integration",
    ):
        if pending.get(key) is not True:
            raise ValueError(f"{key} must remain pending")
    if pending.get("canonical_promotion") is not False:
        raise ValueError("candidate must not be marked canonically promoted")
    return data


def load_recipe(path: str | Path | None = None) -> dict[str, Any]:
    """Load and fail-closed validate the committed cubic delivery recipe."""
    recipe_path = Path(path) if path is not None else _RECIPE_PATH
    return _validate_recipe(json.loads(recipe_path.read_text(encoding="utf-8")))


def build_candidate(
    recipe: Mapping[str, Any] | None = None,
) -> Eq45SupportedPhi10CubicLocalizedTemporalCandidate:
    """Reconstruct exactly the already-materialized cubic temporal candidate."""
    payload = load_recipe() if recipe is None else _validate_recipe(recipe)
    base = Eq45SupportedVelocityCandidate(parent=Eq45VelocityCandidate.seed())
    if base.sha256 != EXPECTED_BASE_SHA256 or base.sha256 != payload["base_supported_sha256"]:
        raise RuntimeError("governed supported parent no longer matches the delivery recipe")
    candidate = Eq45SupportedPhi10CubicLocalizedTemporalCandidate(
        base=base,
        early_delta=EXPECTED_EARLY_DELTA,
    )
    if candidate.sha256 != EXPECTED_CANDIDATE_SHA256:
        raise RuntimeError("reconstructed cubic candidate identity drifted")
    if candidate.sha256 != payload["candidate_sha256"]:
        raise RuntimeError("recipe candidate identity does not match reconstructed field")
    actual = tuple(float(candidate.coefficient_at(t)) for t in EXPECTED_TIMES)
    if not np.allclose(actual, EXPECTED_COEFFICIENTS, rtol=0.0, atol=1.0e-14):
        raise RuntimeError("reconstructed cubic coefficient schedule drifted")
    truth = candidate.to_dict()["truth_boundary"]
    if truth.get("velocity_export_ready") is not True:
        raise RuntimeError("reconstructed candidate lost export readiness")
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


def delivery_capsule(recipe_path: str | Path | None = None) -> dict[str, Any]:
    """Return checked identity, provenance, validation and truth-boundary metadata."""
    recipe = load_recipe(recipe_path)
    candidate = build_candidate(recipe)
    return {
        "schema": CAPSULE_SCHEMA,
        "task_id": TASK_ID,
        "recipe_sha256": _canonical_sha256(recipe),
        "candidate_sha256": candidate.sha256,
        "base_supported_sha256": candidate.base_sha256,
        "temporal_mode": dict(recipe["temporal_mode"]),
        "delivery": dict(recipe["delivery"]),
        "representation_evidence": dict(recipe["representation_evidence"]),
        "temporal_derivative_evidence": dict(recipe["temporal_derivative_evidence"]),
        "sibling_evidence": dict(recipe["sibling_evidence"]),
        "status": dict(recipe["status"]),
        "pending": dict(recipe["pending"]),
    }


def write_delivery_bundle(
    directory: str | Path,
    *,
    recipe_path: str | Path | None = None,
) -> dict[str, Path]:
    """Write replayable ``candidate.json`` and ``capsule.json`` into a directory."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    recipe = load_recipe(recipe_path)
    candidate = build_candidate(recipe)
    capsule = delivery_capsule(recipe_path)
    candidate_path = directory / "candidate.json"
    capsule_path = directory / "capsule.json"
    candidate.save_json(candidate_path)
    capsule_path.write_text(
        json.dumps(capsule, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    loaded = Eq45SupportedPhi10CubicLocalizedTemporalCandidate.load_json(candidate_path)
    if loaded.sha256 != capsule["candidate_sha256"]:
        raise RuntimeError("written candidate no longer matches the delivery capsule identity")
    return {"candidate": candidate_path, "capsule": capsule_path}
