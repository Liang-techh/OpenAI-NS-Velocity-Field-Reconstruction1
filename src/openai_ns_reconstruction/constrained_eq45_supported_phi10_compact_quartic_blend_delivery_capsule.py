"""Replayable delivery capsule for the compact--quartic Phi(1,0) blend family.

This packages the already-materialized support-connected Eq45 temporal blend
together with the already-completed public-velocity time-derivative audit.
The caller must declare the blend weight explicitly.  This module does not
select that weight, change the field, fit pressure/forcing, consume concurrent
force/spectral siblings, or promote visualization/PDE/OpenAI-field readiness.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .constrained_eq45_candidate import Eq45VelocityCandidate
from .constrained_eq45_supported_candidate import Eq45SupportedVelocityCandidate
from .constrained_eq45_supported_phi10_compact_quartic_blend_temporal_mode import (
    Eq45SupportedPhi10CompactQuarticBlendTemporalCandidate,
)

SCHEMA = "eq45_supported_phi10_compact_quartic_blend_delivery_recipe_v1"
CAPSULE_SCHEMA = "eq45_supported_phi10_compact_quartic_blend_delivery_capsule_v1"
TASK_ID = "CR011-EQ45-PHI10-COMPACT-QUARTIC-BLEND-DELIVERY-CAPSULE-018"
DEPENDENCY_HEAD = "37be0d6ee963cd9b98b8853df00be731d8e28bf4"
EXPECTED_BASE_SHA256 = "2fdff812c22131d56eac1b7d6e455187ff3207500c6385aa9abf282a8e3d7b1d"
EXPECTED_EARLY_DELTA = -1.4
REPRESENTATIVE_SHA256 = {
    "0.000000": "6eb71af7ebafc84b246e1855248582c7f84dca37abc263bf3295fae788b1cd68",
    "0.500000": "3fc3fa3a1770931709a4b76fc13942b776386aa7bf40edc992dc8b2514af19d4",
    "1.000000": "99373abdf5c4cb453d36bf1a8086b620e7dda7eea5bfbfbed8947aacea9cedad",
}
EXPECTED_REFERENCE_TIMES = (0.25, 0.3125, 0.375, 0.4375, 0.5, 0.75)

_RECIPE_PATH = (
    Path(__file__).resolve().parents[2]
    / "artifacts/constrained/eq45_supported_phi10_compact_quartic_blend_recipe.json"
)

_EXPECTED_STATUS = {
    "callable_serializable": True,
    "velocity_export_ready": True,
    "visualization_candidate_only": True,
    "caller_declared_blend_weight": True,
    "blend_weight_selected": False,
    "candidate_selection_resolved": False,
    "canonical_velocity_changed": False,
    "physical_support_connection_implemented": True,
    "physical_support_validated": False,
    "force_or_pressure_fitted": False,
    "public_image_fitted": False,
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
        raise ValueError("blend delivery recipe must be an object")
    data = dict(payload)
    if data.get("schema") != SCHEMA:
        raise ValueError(f"recipe schema must be {SCHEMA!r}")
    if data.get("task_id") != TASK_ID:
        raise ValueError("delivery task identity drifted")
    if data.get("base_dependency_head") != DEPENDENCY_HEAD:
        raise ValueError("delivery dependency head drifted")
    if data.get("base_supported_sha256") != EXPECTED_BASE_SHA256:
        raise ValueError("supported parent identity drifted")

    representation = data.get("representation")
    if not isinstance(representation, Mapping):
        raise ValueError("representation must be an object")
    if representation.get("family") != "phi" or representation.get("index") != [1, 0]:
        raise ValueError("temporal mode must remain exactly Phi(1,0)")
    if float(representation.get("early_delta")) != EXPECTED_EARLY_DELTA:
        raise ValueError("early_delta drifted from the screened -1.4 value")
    if (
        representation.get("blend_formula")
        != "(1-lambda)*delta_quartic(tau)+lambda*delta_compact(tau)"
    ):
        raise ValueError("blend formula drifted")
    if representation.get("blend_weight_interval") != [0.0, 1.0]:
        raise ValueError("blend-weight interval drifted")
    if representation.get("blend_weight_selected") is not False:
        raise ValueError("recipe must not select a blend weight")
    if representation.get("endpoint_families") != [
        "derivative_balanced_quartic",
        "slope_capped_c2_compact",
    ]:
        raise ValueError("endpoint family inventory drifted")
    if representation.get("tau_definition") != (
        "delivery_interval_endpoints_map_to_minus1_plus1"
    ):
        raise ValueError("tau definition drifted")
    if representation.get("compact_window_rule") != (
        "earliest_c2_return_under_balanced_quartic_peak_slope_cap"
    ):
        raise ValueError("compact-window rule drifted")
    if representation.get("representative_candidate_sha256") != REPRESENTATIVE_SHA256:
        raise ValueError("representative candidate identities drifted")

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
    if delivery.get("grid_layout") != "time,x,y,z,component":
        raise ValueError("grid layout drifted")
    if delivery.get("component_order") != ["u", "v", "w"]:
        raise ValueError("Cartesian component order drifted")

    evidence = data.get("temporal_derivative_evidence")
    if not isinstance(evidence, Mapping):
        raise ValueError("temporal_derivative_evidence must be an object")
    if evidence.get("audit_pr") != 187 or evidence.get("audit_exact_head") != DEPENDENCY_HEAD:
        raise ValueError("temporal derivative provenance drifted")
    if evidence.get("actions_run") != 35234995835:
        raise ValueError("temporal derivative CI provenance drifted")
    if evidence.get("probe_seed") != 914617 or evidence.get("probe_count") != 24:
        raise ValueError("temporal derivative probe contract drifted")
    if evidence.get("probe_region_counts") != {
        "axial_collar": 6,
        "corner_collar": 4,
        "plateau": 8,
        "radial_collar": 6,
    }:
        raise ValueError("temporal derivative probe regions drifted")
    if evidence.get("blend_weights") != [0.0, 0.5, 1.0]:
        raise ValueError("audited blend-weight inventory drifted")
    if evidence.get("time_levels") != [0.01, 0.005, 0.0025]:
        raise ValueError("time-difference ladder drifted")
    if float(evidence.get("reference_step")) != 0.00125:
        raise ValueError("time-derivative reference step drifted")
    if evidence.get("finest_normalized_rms") != {
        "0.000000": 0.00019000133460634884,
        "0.500000": 0.0006986742153697458,
        "1.000000": 0.0012568302260969515,
    }:
        raise ValueError("registered finest normalized RMS drifted")
    if evidence.get("finest_sampled_max_vector") != {
        "0.000000": 0.0006105780615444523,
        "0.500000": 0.0016384499832866476,
        "1.000000": 0.0032176812930085137,
    }:
        raise ValueError("registered finest max-vector errors drifted")
    if evidence.get("observed_orders") != {
        "0.000000": [1.9814367127354424, 1.993893149355883],
        "0.500000": [1.8015385320652726, 1.90443075105992],
        "1.000000": [1.787112457217247, 1.8980297142790958],
    }:
        raise ValueError("registered temporal convergence orders drifted")
    if float(evidence.get("midpoint_instantaneous_affinity_max_error")) > 1.2e-16:
        raise ValueError("midpoint instantaneous affine identity drifted")
    if float(evidence.get("midpoint_derivative_affinity_normalized_rms_max")) > 5.94e-13:
        raise ValueError("midpoint derivative affine identity drifted")
    if evidence.get("endpoint_derivative_span_rms_at_static_anchors") != {
        "0.500000": 0.0064438239000152,
        "0.750000": 0.0027456278999935003,
    }:
        raise ValueError("registered endpoint derivative spans drifted")
    if not np.isclose(
        float(evidence.get("mutation_good_sampled_rms_error")),
        1.8665839559811222e-15,
        rtol=0.0,
        atol=1.0e-28,
    ):
        raise ValueError("good derivative mutation calibration drifted")
    if not np.isclose(
        float(evidence.get("mutation_sign_flipped_sampled_rms_error")),
        3.464101615137756,
        rtol=0.0,
        atol=1.0e-15,
    ):
        raise ValueError("sign-mutated derivative calibration drifted")
    if float(evidence.get("mutation_sign_flipped_max_vector_error")) != 6.000000000000003:
        raise ValueError("sign-mutated max-vector calibration drifted")
    if evidence.get("public_velocity_only") is not True:
        raise ValueError("temporal derivative evidence must remain public-velocity only")
    if evidence.get("serialized_candidates_reloaded") is not True:
        raise ValueError("temporal derivative evidence lost reload provenance")
    if evidence.get("blend_weight_selected") is not False:
        raise ValueError("temporal derivative audit cannot select the blend weight")
    if evidence.get("formal_pde_gate_assessed") is not False:
        raise ValueError("time-derivative audit is not the formal PDE gate")
    if evidence.get("pde_validated") is not False:
        raise ValueError("time-derivative audit cannot promote PDE validity")

    siblings = data.get("sibling_evidence")
    if not isinstance(siblings, Mapping):
        raise ValueError("sibling_evidence must be an object")
    if siblings.get("restricted_force_pr") != 188:
        raise ValueError("restricted-force sibling provenance drifted")
    if siblings.get("restricted_force_status") != "open_unconsumed":
        raise ValueError("restricted-force sibling must remain unconsumed here")
    if siblings.get("spectral_fingerprint_pr") != 189:
        raise ValueError("spectral sibling provenance drifted")
    if siblings.get("spectral_fingerprint_status") != "open_unconsumed":
        raise ValueError("spectral sibling must remain unconsumed here")
    for key in (
        "numeric_summary_bound_here",
        "may_select_blend_weight",
        "may_promote_visualization_ready",
        "may_promote_pde_validated",
    ):
        if siblings.get(key) is not False:
            raise ValueError(f"{key} must remain false for unconsumed sibling evidence")

    if data.get("status") != _EXPECTED_STATUS:
        raise ValueError("delivery truth-boundary status drifted")

    pending = data.get("pending")
    if not isinstance(pending, Mapping):
        raise ValueError("pending must be an object")
    for key in (
        "restricted_force_tradeoff_integration",
        "spectral_fingerprint_integration",
        "public_observable_comparison",
        "three_dimensional_render_review",
        "formal_pde_gate",
        "candidate_specific_energy_acceptance",
        "blend_weight_selection",
    ):
        if pending.get(key) is not True:
            raise ValueError(f"{key} must remain pending")
    if pending.get("canonical_promotion") is not False:
        raise ValueError("blend family must not be marked canonically promoted")
    return data


def load_recipe(path: str | Path | None = None) -> dict[str, Any]:
    """Load and fail-closed validate the committed blend-family recipe."""
    recipe_path = Path(path) if path is not None else _RECIPE_PATH
    return _validate_recipe(json.loads(recipe_path.read_text(encoding="utf-8")))


def build_candidate(
    blend_weight: float,
    recipe: Mapping[str, Any] | None = None,
) -> Eq45SupportedPhi10CompactQuarticBlendTemporalCandidate:
    """Build one caller-declared member of the governed blend family."""
    payload = load_recipe() if recipe is None else _validate_recipe(recipe)
    weight = float(blend_weight)
    if not np.isfinite(weight) or not (0.0 <= weight <= 1.0):
        raise ValueError("blend_weight must be finite and lie in [0,1]")

    base = Eq45SupportedVelocityCandidate(parent=Eq45VelocityCandidate.seed())
    if base.sha256 != EXPECTED_BASE_SHA256 or base.sha256 != payload["base_supported_sha256"]:
        raise RuntimeError("governed supported parent no longer matches the delivery recipe")
    candidate = Eq45SupportedPhi10CompactQuarticBlendTemporalCandidate(
        base=base,
        blend_weight=weight,
        early_delta=EXPECTED_EARLY_DELTA,
    )
    representative_key = f"{weight:.6f}"
    if representative_key in REPRESENTATIVE_SHA256:
        if candidate.sha256 != REPRESENTATIVE_SHA256[representative_key]:
            raise RuntimeError("representative blend candidate identity drifted")
    truth = candidate.to_dict()["truth_boundary"]
    if truth.get("velocity_export_ready") is not True:
        raise RuntimeError("reconstructed blend candidate lost export readiness")
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
            raise RuntimeError(f"reconstructed blend candidate over-promoted {key}")
    return candidate


def delivery_capsule(
    blend_weight: float,
    recipe_path: str | Path | None = None,
) -> dict[str, Any]:
    """Return checked family provenance plus one explicit member identity."""
    recipe = load_recipe(recipe_path)
    candidate = build_candidate(blend_weight, recipe)
    return {
        "schema": CAPSULE_SCHEMA,
        "task_id": TASK_ID,
        "recipe_sha256": _canonical_sha256(recipe),
        "candidate_sha256": candidate.sha256,
        "base_supported_sha256": candidate.base_sha256,
        "declared_blend_weight": float(candidate.blend_weight),
        "blend_weight_selected_by_capsule": False,
        "representation": dict(recipe["representation"]),
        "delivery": dict(recipe["delivery"]),
        "temporal_derivative_evidence": dict(recipe["temporal_derivative_evidence"]),
        "sibling_evidence": dict(recipe["sibling_evidence"]),
        "status": dict(recipe["status"]),
        "pending": dict(recipe["pending"]),
    }


def write_delivery_bundle(
    output_dir: str | Path,
    blend_weight: float,
    recipe_path: str | Path | None = None,
) -> dict[str, Path]:
    """Write candidate.json + capsule.json and verify exact public artifact replay."""
    directory = Path(output_dir)
    directory.mkdir(parents=True, exist_ok=True)
    candidate = build_candidate(blend_weight, load_recipe(recipe_path))
    capsule = delivery_capsule(blend_weight, recipe_path)

    candidate_path = directory / "candidate.json"
    capsule_path = directory / "capsule.json"
    candidate.save_json(candidate_path)
    reloaded = Eq45SupportedPhi10CompactQuarticBlendTemporalCandidate.load_json(
        candidate_path
    )
    if reloaded.sha256 != candidate.sha256:
        raise RuntimeError("serialized blend candidate changed identity after reload")
    if capsule["candidate_sha256"] != reloaded.sha256:
        raise RuntimeError("capsule identity does not match reloaded blend candidate")
    if capsule["declared_blend_weight"] != reloaded.blend_weight:
        raise RuntimeError("capsule blend weight does not match reloaded candidate")

    capsule_path.write_text(
        json.dumps(capsule, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return {"candidate": candidate_path, "capsule": capsule_path}
