"""Replayable delivery capsule for the derivative-balanced quartic Phi(1,0) field.

This packages the already-materialized quartic support-connected Eq45 velocity
candidate together with the already-completed fixed train/disjoint-holdout
restricted-force diagnostic. It does not change the field, fit a new force or
pressure, consume concurrent morphology/fresh-seed siblings, or promote
visualization/PDE/canonical readiness.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .constrained_eq45_candidate import Eq45VelocityCandidate
from .constrained_eq45_supported_candidate import Eq45SupportedVelocityCandidate
from .constrained_eq45_supported_phi10_quartic_temporal_mode import (
    Eq45SupportedPhi10QuarticDerivativeBalancedTemporalCandidate,
)

SCHEMA = "eq45_supported_phi10_quartic_delivery_recipe_v1"
CAPSULE_SCHEMA = "eq45_supported_phi10_quartic_delivery_capsule_v1"
TASK_ID = "CR011-EQ45-PHI10-QUARTIC-DELIVERY-CAPSULE-016"
DEPENDENCY_HEAD = "d3a4c1d6895c6ef08c73bc27d933c0dc1e74ea90"
EXPECTED_BASE_SHA256 = "2fdff812c22131d56eac1b7d6e455187ff3207500c6385aa9abf282a8e3d7b1d"
EXPECTED_CANDIDATE_SHA256 = "03fdae73170469ae1160297b489b1b83a27adddfc941119e116a21b98bd15049"
EXPECTED_EARLY_DELTA = -1.4
EXPECTED_NULLSPACE_COEFFICIENT = -0.2831460674157303
EXPECTED_TIMES = (0.25, 0.375, 0.5, 0.5625, 0.625, 0.6875, 0.75)
EXPECTED_COEFFICIENTS = (
    -1.7,
    -0.5438202247191011,
    -0.3,
    -0.29471558988764046,
    -0.3,
    -0.2986481741573033,
    -0.3,
)

_RECIPE_PATH = (
    Path(__file__).resolve().parents[2]
    / "artifacts/constrained/eq45_supported_phi10_quartic_recipe.json"
)

_EXPECTED_STATUS = {
    "callable_serializable": True,
    "velocity_export_ready": True,
    "visualization_candidate_only": True,
    "canonical_velocity_changed": False,
    "screened_temporal_shape_materialized": True,
    "derivative_balanced_temporal_shape": True,
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
        raise ValueError("quartic delivery recipe must be an object")
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
        raise ValueError("quartic candidate identity drifted")

    mode = data.get("temporal_mode")
    if not isinstance(mode, Mapping):
        raise ValueError("temporal_mode must be an object")
    if mode.get("family") != "phi" or mode.get("index") != [1, 0]:
        raise ValueError("temporal mode must remain exactly Phi(1,0)")
    if float(mode.get("early_delta")) != EXPECTED_EARLY_DELTA:
        raise ValueError("early_delta drifted from the screened -1.4 value")
    if float(mode.get("nullspace_coefficient")) != EXPECTED_NULLSPACE_COEFFICIENT:
        raise ValueError("balanced quartic nullspace coefficient drifted")
    if (
        mode.get("nullspace_rule")
        != "least_squares_minimize_return_node_coefficient_slope"
    ):
        raise ValueError("balanced-nullspace rule drifted")
    if (
        mode.get("schedule")
        != "cubic_phi10_delta(tau,early_delta)+a*(tau+1)*tau*(tau-0.5)*(tau-1)"
    ):
        raise ValueError("quartic temporal schedule drifted")
    if mode.get("coefficient_bound") != [-4.0, 4.0]:
        raise ValueError("profile coefficient bound drifted")
    if tuple(float(x) for x in mode.get("reference_times", [])) != EXPECTED_TIMES:
        raise ValueError("temporal reference times drifted")
    coefficients = tuple(float(x) for x in mode.get("reference_coefficients", []))
    if len(coefficients) != len(EXPECTED_COEFFICIENTS) or not np.allclose(
        coefficients, EXPECTED_COEFFICIENTS, rtol=0.0, atol=1.0e-14
    ):
        raise ValueError("registered quartic coefficients drifted")

    delivery = data.get("delivery")
    if not isinstance(delivery, Mapping):
        raise ValueError("delivery must be an object")
    if delivery.get("time_interval") != [0.25, 0.75]:
        raise ValueError("delivery time interval drifted")
    if tuple(float(x) for x in delivery.get("reference_times", [])) != EXPECTED_TIMES:
        raise ValueError("delivery reference times drifted")
    if delivery.get("spatial_box") != [[-2.0, 2.0], [-2.0, 2.0], [-2.0, 2.0]]:
        raise ValueError("delivery spatial box drifted")
    if set(delivery.get("velocity_interfaces", [])) != {
        "velocity",
        "at_points",
        "velocity_xyz",
        "grid",
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
    if representation.get("capacity_pr") != 166:
        raise ValueError("quartic capacity provenance drifted")
    if representation.get("materialization_pr") != 168:
        raise ValueError("quartic materialization provenance drifted")
    if representation.get("early_identity") != "quartic_equals_cubic_at_t0p25":
        raise ValueError("early snapshot identity drifted")
    if representation.get("static_return_times") != [0.5, 0.625, 0.75]:
        raise ValueError("static-return keyframes drifted")
    if representation.get("off_keyframe_nontrivial_times") != [0.375, 0.5625, 0.6875]:
        raise ValueError("off-keyframe inventory drifted")
    if float(representation.get("return_node_slope_l2_ratio_vs_cubic")) >= 0.31:
        raise ValueError("registered derivative-balancing advantage drifted")
    if float(representation.get("late_half_excursion_ratio_vs_cubic")) >= 0.29:
        raise ValueError("registered quartic late-half localization advantage drifted")
    if float(representation.get("start_endpoint_slope_ratio_vs_cubic")) <= 1.27:
        raise ValueError("quartic early-time sharpening warning drifted")
    if representation.get("public_openai_reference_used_for_selection") is not False:
        raise ValueError("target-free quartic selection cannot be relabeled as image fitting")
    if representation.get("visual_correspondence_established") is not False:
        raise ValueError("representation evidence cannot establish public correspondence")

    force = data.get("restricted_force_evidence")
    if not isinstance(force, Mapping):
        raise ValueError("restricted_force_evidence must be an object")
    if force.get("audit_pr") != 169 or force.get("audit_exact_head") != DEPENDENCY_HEAD:
        raise ValueError("quartic force-crosscheck provenance drifted")
    if force.get("force_family") != "RestrictedForce(a,c)":
        raise ValueError("force family drifted")
    if force.get("bounds") != {"a": [0.0, 10.0], "c": [0.0, 10.0]}:
        raise ValueError("restricted-force bounds drifted")
    if force.get("training_probe_count") != 16 or force.get("holdout_probe_count") != 16:
        raise ValueError("train/holdout probe contract drifted")
    if force.get("spatial_derivative_steps") != [0.02, 0.01, 0.005]:
        raise ValueError("spatial derivative ladder drifted")
    if float(force.get("time_step")) != 0.0025 or float(force.get("nu")) != 0.01:
        raise ValueError("quartic force diagnostic discretization drifted")
    coefficients_force = force.get("frozen_force_coefficients")
    if not isinstance(coefficients_force, Mapping):
        raise ValueError("frozen_force_coefficients must be an object")
    if not np.isclose(float(coefficients_force.get("a")), 2.1615475834844833e-15):
        raise ValueError("fitted restricted-force a coefficient drifted")
    if not np.isclose(float(coefficients_force.get("c")), 0.07770084654400304):
        raise ValueError("fitted restricted-force c coefficient drifted")
    before = tuple(float(x) for x in force.get("heldout_rms_before", []))
    after = tuple(float(x) for x in force.get("heldout_rms_after", []))
    if before != (
        3.678322090166512,
        3.6969949228504446,
        3.7016892575139453,
    ):
        raise ValueError("registered quartic held-out zero-force RMS drifted")
    if after != (
        3.677015083694525,
        3.6957137765888333,
        3.700414843717181,
    ):
        raise ValueError("registered quartic held-out projected RMS drifted")
    if force.get("pressure_fitted") is not False:
        raise ValueError("pressure-free diagnostic cannot be relabeled as pressure fit")
    if force.get("formal_pde_gate_assessed") is not False:
        raise ValueError("restricted-force diagnostic is not the formal PDE gate")
    if force.get("pde_validated") is not False:
        raise ValueError("restricted-force diagnostic cannot promote PDE validity")

    siblings = data.get("sibling_evidence")
    if not isinstance(siblings, Mapping):
        raise ValueError("sibling_evidence must be an object")
    if siblings.get("quartic_seed_generalization_task") != (
        "CR009-EQ45-SUPPORTED-PHI10-QUARTIC-SEED-GENERALIZATION-028"
    ):
        raise ValueError("quartic seed-generalization sibling identity drifted")
    if siblings.get("seed_generalization_status") != "claimed_unconsumed":
        raise ValueError("fresh-seed sibling must remain unconsumed in this ancestry")
    if siblings.get("late_morphology_pr") != 170:
        raise ValueError("quartic morphology sibling provenance drifted")
    if siblings.get("late_morphology_status") != "open_unconsumed_calibration_pending":
        raise ValueError("quartic morphology sibling state drifted")
    if siblings.get("numeric_summary_bound_here") is not False:
        raise ValueError("unconsumed sibling numerics cannot be laundered into this capsule")
    if siblings.get("may_promote_visualization_ready") is not False:
        raise ValueError("open sibling evidence cannot promote visualization readiness")
    if siblings.get("may_promote_pde_validated") is not False:
        raise ValueError("open sibling evidence cannot promote PDE validity")

    if data.get("status") != _EXPECTED_STATUS:
        raise ValueError("delivery truth-boundary status drifted")

    pending = data.get("pending")
    if not isinstance(pending, Mapping):
        raise ValueError("pending must be an object")
    for key in (
        "quartic_specific_temporal_derivative_audit",
        "fresh_seed_generalization_integration",
        "late_morphology_integration",
        "public_observable_comparison",
        "three_dimensional_render_review",
        "formal_pde_gate",
        "candidate_specific_energy_acceptance",
    ):
        if pending.get(key) is not True:
            raise ValueError(f"{key} must remain pending")
    if pending.get("canonical_promotion") is not False:
        raise ValueError("candidate must not be marked canonically promoted")
    return data


def load_recipe(path: str | Path | None = None) -> dict[str, Any]:
    """Load and fail-closed validate the committed quartic delivery recipe."""
    recipe_path = Path(path) if path is not None else _RECIPE_PATH
    return _validate_recipe(json.loads(recipe_path.read_text(encoding="utf-8")))


def build_candidate(
    recipe: Mapping[str, Any] | None = None,
) -> Eq45SupportedPhi10QuarticDerivativeBalancedTemporalCandidate:
    """Reconstruct exactly the already-materialized quartic temporal candidate."""
    payload = load_recipe() if recipe is None else _validate_recipe(recipe)
    base = Eq45SupportedVelocityCandidate(parent=Eq45VelocityCandidate.seed())
    if base.sha256 != EXPECTED_BASE_SHA256 or base.sha256 != payload["base_supported_sha256"]:
        raise RuntimeError("governed supported parent no longer matches the delivery recipe")
    candidate = Eq45SupportedPhi10QuarticDerivativeBalancedTemporalCandidate(
        base=base,
        early_delta=EXPECTED_EARLY_DELTA,
    )
    if candidate.sha256 != EXPECTED_CANDIDATE_SHA256:
        raise RuntimeError("reconstructed quartic candidate identity drifted")
    if candidate.sha256 != payload["candidate_sha256"]:
        raise RuntimeError("recipe candidate identity does not match reconstructed field")
    if candidate.nullspace_coefficient != EXPECTED_NULLSPACE_COEFFICIENT:
        raise RuntimeError("reconstructed quartic balancing coefficient drifted")
    actual = tuple(float(candidate.coefficient_at(t)) for t in EXPECTED_TIMES)
    if not np.allclose(actual, EXPECTED_COEFFICIENTS, rtol=0.0, atol=1.0e-14):
        raise RuntimeError("reconstructed quartic coefficient schedule drifted")
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
    """Return checked identity, provenance, diagnostics and truth-boundary metadata."""
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
        "restricted_force_evidence": dict(recipe["restricted_force_evidence"]),
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
    loaded = Eq45SupportedPhi10QuarticDerivativeBalancedTemporalCandidate.load_json(
        candidate_path
    )
    if loaded.sha256 != capsule["candidate_sha256"]:
        raise RuntimeError("written candidate no longer matches the delivery capsule identity")
    return {"candidate": candidate_path, "capsule": capsule_path}
