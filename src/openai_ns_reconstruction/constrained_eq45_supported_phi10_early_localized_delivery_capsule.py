"""Reproducible delivery capsule for the early-localized supported Phi(1,0) field.

This module packages the already-screened quadratic temporal candidate from the
active constrained stack.  It does not select another representation, rerun an
optimizer, refit pressure/forcing, or recompute PDE/morphology diagnostics.

The capsule intentionally keeps three states independent: callable/export-ready
velocity, visualization acceptance, and PDE acceptance.  A replayable [u,v,w]
artifact is useful even while the latter two remain false.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .constrained_eq45_candidate import Eq45VelocityCandidate
from .constrained_eq45_supported_candidate import Eq45SupportedVelocityCandidate
from .constrained_eq45_supported_phi10_early_localized_temporal_mode import (
    Eq45SupportedPhi10EarlyLocalizedTemporalCandidate,
)

SCHEMA = "eq45_supported_phi10_early_localized_delivery_recipe_v1"
CAPSULE_SCHEMA = "eq45_supported_phi10_early_localized_delivery_capsule_v1"
TASK_ID = "CR011-EQ45-PHI10-EARLY-LOCALIZED-DELIVERY-CAPSULE-014"
DEPENDENCY_HEAD = "73429a605c750707cbcc5e6e26ddcc5511894ac4"
EXPECTED_BASE_SHA256 = "2fdff812c22131d56eac1b7d6e455187ff3207500c6385aa9abf282a8e3d7b1d"
EXPECTED_CANDIDATE_SHA256 = "5b08abbfdce0e96d39a5fc394ceedea0d0ff549cb921ff50e84004355c0e2bfe"
EXPECTED_EARLY_DELTA = -1.4
EXPECTED_TIMES = (0.25, 0.5, 0.625, 0.75)
EXPECTED_COEFFICIENTS = (-1.7, -0.3, -0.125, -0.3)

_RECIPE_PATH = (
    Path(__file__).resolve().parents[2]
    / "artifacts/constrained/eq45_supported_phi10_early_localized_recipe.json"
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
        raise ValueError("early-localized delivery recipe must be an object")
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
        raise ValueError("early-localized candidate identity drifted")

    mode = data.get("temporal_mode")
    if not isinstance(mode, Mapping):
        raise ValueError("temporal_mode must be an object")
    if mode.get("family") != "phi" or mode.get("index") != [1, 0]:
        raise ValueError("temporal mode must remain exactly Phi(1,0)")
    if float(mode.get("early_delta")) != EXPECTED_EARLY_DELTA:
        raise ValueError("early_delta drifted from the screened -1.4 value")
    if mode.get("schedule") != "0.5*early_delta*tau*(tau-1)":
        raise ValueError("quadratic temporal schedule drifted")
    if mode.get("coefficient_bound") != [-4.0, 4.0]:
        raise ValueError("profile coefficient bound drifted")
    coefficients = tuple(
        float(mode.get(key))
        for key in (
            "coefficient_start",
            "coefficient_midpoint",
            "coefficient_intermediate_t0625",
            "coefficient_end",
        )
    )
    if coefficients != EXPECTED_COEFFICIENTS:
        raise ValueError("registered temporal coefficients drifted")

    delivery = data.get("delivery")
    if not isinstance(delivery, Mapping):
        raise ValueError("delivery must be an object")
    if delivery.get("time_interval") != [0.25, 0.75]:
        raise ValueError("delivery time interval drifted")
    if delivery.get("reference_times") != list(EXPECTED_TIMES):
        raise ValueError("delivery reference times drifted")
    if delivery.get("spatial_box") != [[-2.0, 2.0], [-2.0, 2.0], [-2.0, 2.0]]:
        raise ValueError("delivery spatial box drifted")
    required_velocity = {"velocity", "at_points", "velocity_xyz", "grid"}
    if set(delivery.get("velocity_interfaces", [])) != required_velocity:
        raise ValueError("public velocity interface inventory drifted")
    if set(delivery.get("serialization_interfaces", [])) != {"save_json", "load_json"}:
        raise ValueError("serialization interface inventory drifted")
    if delivery.get("component_order") != ["u", "v", "w"]:
        raise ValueError("Cartesian component order drifted")

    selection = data.get("selection_evidence")
    if not isinstance(selection, Mapping):
        raise ValueError("selection_evidence must be an object")
    if selection.get("capacity_pr") != 150 or selection.get("whole_domain_morphology_pr") != 146:
        raise ValueError("morphology evidence provenance drifted")
    if selection.get("public_openai_reference_used_for_selection") is not False:
        raise ValueError("target-free morphology selection cannot be relabeled as public-image fitting")
    if selection.get("visual_correspondence_established") is not False:
        raise ValueError("internal morphology evidence cannot establish public correspondence")
    if float(selection.get("early_radial_q99_trial")) >= float(selection.get("early_radial_q99_baseline")):
        raise ValueError("registered early radial-envelope improvement drifted")
    if float(selection.get("early_radial_outer_vorticity2_fraction_trial")) >= float(
        selection.get("early_radial_outer_vorticity2_fraction_baseline")
    ):
        raise ValueError("registered early collar-energy improvement drifted")

    pde = data.get("pde_evidence")
    if not isinstance(pde, Mapping):
        raise ValueError("pde_evidence must be an object")
    if pde.get("fixed_holdout_pr") != 152 or pde.get("fresh_seed_generalization_pr") != 153:
        raise ValueError("PDE diagnostic provenance drifted")
    if pde.get("fresh_seed_generalization_exact_head") != DEPENDENCY_HEAD:
        raise ValueError("fresh-seed dependency identity drifted")
    if pde.get("fresh_seed_validation_seeds") != [914457, 914471, 914483]:
        raise ValueError("fresh-seed validation set drifted")
    if pde.get("fresh_seed_ci_success") is not True:
        raise ValueError("capsule requires successful exact-head fresh-seed CI")
    if pde.get("fresh_seed_numeric_summary_bound_here") is not False:
        raise ValueError("unharvested fresh-seed numerics must not be invented in this capsule")
    if pde.get("formal_pde_gate_assessed") is not False or pde.get("pde_validated") is not False:
        raise ValueError("fixed/fresh probe diagnostics cannot be promoted to formal PDE validation")
    if float(pde.get("early_localized_zero_force_rms")) <= float(pde.get("static_zero_force_rms")):
        raise ValueError("registered fixed-holdout PDE warning drifted")

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
    ):
        if pending.get(key) is not True:
            raise ValueError(f"{key} must remain pending")
    if pending.get("canonical_promotion") is not False:
        raise ValueError("candidate must not be marked canonically promoted")
    return data


def load_recipe(path: str | Path | None = None) -> dict[str, Any]:
    """Load and fail-closed validate the committed delivery recipe."""
    recipe_path = Path(path) if path is not None else _RECIPE_PATH
    return _validate_recipe(json.loads(recipe_path.read_text(encoding="utf-8")))


def build_candidate(
    recipe: Mapping[str, Any] | None = None,
) -> Eq45SupportedPhi10EarlyLocalizedTemporalCandidate:
    """Reconstruct exactly the already-screened callable temporal candidate."""
    payload = load_recipe() if recipe is None else _validate_recipe(recipe)
    base = Eq45SupportedVelocityCandidate(parent=Eq45VelocityCandidate.seed())
    if base.sha256 != EXPECTED_BASE_SHA256 or base.sha256 != payload["base_supported_sha256"]:
        raise RuntimeError("governed supported parent no longer matches the delivery recipe")
    candidate = Eq45SupportedPhi10EarlyLocalizedTemporalCandidate(
        base=base,
        early_delta=EXPECTED_EARLY_DELTA,
    )
    if candidate.sha256 != EXPECTED_CANDIDATE_SHA256:
        raise RuntimeError("reconstructed early-localized candidate identity drifted")
    if candidate.sha256 != payload["candidate_sha256"]:
        raise RuntimeError("recipe candidate identity does not match reconstructed field")

    actual_coefficients = tuple(float(candidate.coefficient_at(t)) for t in EXPECTED_TIMES)
    if not np.allclose(actual_coefficients, EXPECTED_COEFFICIENTS, rtol=0.0, atol=1.0e-14):
        raise RuntimeError("reconstructed temporal coefficient schedule drifted")
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
    """Return the checked identity/evidence/status capsule for downstream consumers."""
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
        "selection_evidence": dict(recipe["selection_evidence"]),
        "pde_evidence": dict(recipe["pde_evidence"]),
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
        json.dumps(capsule, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    loaded = Eq45SupportedPhi10EarlyLocalizedTemporalCandidate.load_json(candidate_path)
    if loaded.sha256 != capsule["candidate_sha256"]:
        raise RuntimeError("written candidate no longer matches the delivery capsule identity")
    return {"candidate": candidate_path, "capsule": capsule_path}
