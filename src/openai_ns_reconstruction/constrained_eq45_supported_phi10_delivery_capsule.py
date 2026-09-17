"""Reproducible delivery capsule for the selected supported Phi(1,0) temporal trial.

This module does not select another representation and does not rerun any PDE or
morphology optimizer.  It binds the already-selected target-free slope-1.4 trial
to a checked recipe and provides a small bundle writer so downstream MATLAB /
Python visualization and validation lanes can consume the *same* public
``[u,v,w]`` object.

The recipe deliberately keeps export readiness, visual acceptance, and PDE
acceptance independent.  In particular, a serializable velocity candidate is not
a claim of OpenAI-image correspondence or Navier--Stokes validity.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .constrained_eq45_supported_phi10_temporal_replay import build_phi10_temporal_trial
from .constrained_eq45_supported_temporal_mode import (
    Eq45SupportedAffineTemporalModeCandidate,
)

SCHEMA = "eq45_supported_phi10_temporal_trial_recipe_v1"
CAPSULE_SCHEMA = "eq45_supported_phi10_temporal_delivery_capsule_v1"
TASK_ID = "CR011-EQ45-PHI10-TEMPORAL-TRIAL-CAPSULE-013"
DEPENDENCY_HEAD = "83eaa359e63023de7d43abb60988808c1c5d8378"
EXPECTED_BASE_SHA256 = "2fdff812c22131d56eac1b7d6e455187ff3207500c6385aa9abf282a8e3d7b1d"
EXPECTED_MODE = ("phi", 1, 0)
EXPECTED_SLOPE = 1.4
EXPECTED_TIMES = (0.25, 0.5, 0.75)

_RECIPE_PATH = (
    Path(__file__).resolve().parents[2]
    / "artifacts/constrained/eq45_supported_phi10_temporal_trial_recipe.json"
)

_EXPECTED_STATUS = {
    "velocity_export_ready": True,
    "visualization_candidate_only": True,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "physical_support_validated": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "canonical_velocity_changed": False,
    "production_slope_promoted": False,
}


def _canonical_sha256(payload: Mapping[str, Any]) -> str:
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def load_trial_recipe(path: str | Path | None = None) -> dict[str, Any]:
    """Load and fail-closed validate the committed trial recipe."""
    recipe_path = Path(path) if path is not None else _RECIPE_PATH
    payload = json.loads(recipe_path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("temporal trial recipe must be a JSON object")
    if payload.get("schema") != SCHEMA:
        raise ValueError(f"recipe schema must be {SCHEMA!r}")
    if payload.get("task_id") != TASK_ID:
        raise ValueError("temporal trial task identity drifted")
    if payload.get("base_dependency_head") != DEPENDENCY_HEAD:
        raise ValueError("temporal trial dependency head drifted")
    if payload.get("base_supported_sha256") != EXPECTED_BASE_SHA256:
        raise ValueError("temporal trial base supported identity drifted")

    mode = payload.get("temporal_mode")
    if not isinstance(mode, dict):
        raise ValueError("temporal_mode must be an object")
    index = mode.get("index")
    if mode.get("family") != EXPECTED_MODE[0] or index != [EXPECTED_MODE[1], EXPECTED_MODE[2]]:
        raise ValueError("temporal trial mode drifted from Phi(1,0)")
    if float(mode.get("slope")) != EXPECTED_SLOPE:
        raise ValueError("temporal trial slope drifted from the selected 1.4 value")
    if [float(mode.get(key)) for key in ("coefficient_start", "coefficient_midpoint", "coefficient_end")] != [-1.7, -0.3, 1.1]:
        raise ValueError("temporal trial endpoint coefficients drifted")
    if mode.get("coefficient_bound") != [-4.0, 4.0]:
        raise ValueError("temporal trial coefficient bound drifted")

    delivery = payload.get("delivery")
    if not isinstance(delivery, dict):
        raise ValueError("delivery must be an object")
    if delivery.get("time_interval") != [0.25, 0.75]:
        raise ValueError("delivery time interval drifted")
    if delivery.get("reference_times") != list(EXPECTED_TIMES):
        raise ValueError("delivery reference times drifted")
    if delivery.get("spatial_box") != [[-2.0, 2.0], [-2.0, 2.0], [-2.0, 2.0]]:
        raise ValueError("delivery spatial box drifted")

    if payload.get("status") != _EXPECTED_STATUS:
        raise ValueError("temporal trial truth-boundary status drifted")
    pde = payload.get("pde_evidence")
    if not isinstance(pde, dict):
        raise ValueError("pde_evidence must be an object")
    if pde.get("uniform_generalization") is not False or pde.get("formal_pde_gate_assessed") is not False:
        raise ValueError("PDE diagnostic must remain mixed and non-formal")
    pending = payload.get("pending")
    if not isinstance(pending, dict) or pending.get("whole_domain_morphology_acceptance") != "PR146_open_unconsumed":
        raise ValueError("whole-domain morphology acceptance must remain pending")
    if pending.get("canonical_promotion") is not False:
        raise ValueError("temporal trial must not be marked canonically promoted")
    return payload


def build_trial_candidate(
    recipe: Mapping[str, Any] | None = None,
) -> Eq45SupportedAffineTemporalModeCandidate:
    """Reconstruct exactly the selected callable support-connected trial."""
    payload = load_trial_recipe() if recipe is None else dict(recipe)
    if recipe is not None:
        # Reuse the same fail-closed checks for caller-supplied recipe data.
        candidate_path = None
        if payload.get("schema") != SCHEMA or payload.get("status") != _EXPECTED_STATUS:
            raise ValueError("caller-supplied recipe failed the frozen schema/status contract")
        if payload.get("base_dependency_head") != DEPENDENCY_HEAD:
            raise ValueError("caller-supplied recipe dependency head drifted")
        if payload.get("base_supported_sha256") != EXPECTED_BASE_SHA256:
            raise ValueError("caller-supplied recipe base identity drifted")
        mode = payload.get("temporal_mode", {})
        if mode.get("family") != "phi" or mode.get("index") != [1, 0] or float(mode.get("slope")) != EXPECTED_SLOPE:
            raise ValueError("caller-supplied recipe mode/slope drifted")
        del candidate_path

    trial = build_phi10_temporal_trial(slope=EXPECTED_SLOPE)
    if trial.base_sha256 != EXPECTED_BASE_SHA256:
        raise RuntimeError("governed supported parent no longer matches the frozen recipe")
    if (trial.family, trial.mode_i, trial.mode_j) != EXPECTED_MODE:
        raise RuntimeError("reconstructed temporal mode no longer matches Phi(1,0)")
    if float(trial.slope) != EXPECTED_SLOPE:
        raise RuntimeError("reconstructed temporal slope no longer matches 1.4")

    expected_coefficients = (-1.7, -0.3, 1.1)
    actual_coefficients = tuple(float(trial.coefficient_at(time)) for time in EXPECTED_TIMES)
    if not np.allclose(actual_coefficients, expected_coefficients, rtol=0.0, atol=1.0e-14):
        raise RuntimeError("reconstructed temporal endpoint coefficients drifted")

    truth = trial.to_dict()["truth_boundary"]
    if truth.get("velocity_export_ready") is not True:
        raise RuntimeError("reconstructed temporal trial lost export readiness")
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
            raise RuntimeError(f"reconstructed temporal trial over-promoted {key}")
    return trial


def delivery_capsule(recipe_path: str | Path | None = None) -> dict[str, Any]:
    """Return the checked identity/status capsule used by downstream consumers."""
    recipe = load_trial_recipe(recipe_path)
    trial = build_trial_candidate(recipe)
    return {
        "schema": CAPSULE_SCHEMA,
        "task_id": TASK_ID,
        "recipe_sha256": _canonical_sha256(recipe),
        "candidate_sha256": trial.sha256,
        "base_supported_sha256": trial.base_sha256,
        "temporal_mode": {
            "family": trial.family,
            "index": [trial.mode_i, trial.mode_j],
            "slope": float(trial.slope),
        },
        "delivery": dict(recipe["delivery"]),
        "selection_evidence": dict(recipe["selection_evidence"]),
        "pde_evidence": dict(recipe["pde_evidence"]),
        "status": dict(recipe["status"]),
        "pending": dict(recipe["pending"]),
    }


def write_trial_bundle(
    directory: str | Path,
    *,
    recipe_path: str | Path | None = None,
) -> dict[str, Path]:
    """Write replayable candidate JSON plus its checked capsule to ``directory``."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    trial = build_trial_candidate(load_trial_recipe(recipe_path))
    capsule = delivery_capsule(recipe_path)
    candidate_path = directory / "candidate.json"
    capsule_path = directory / "capsule.json"
    trial.save_json(candidate_path)
    capsule_path.write_text(
        json.dumps(capsule, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    loaded = Eq45SupportedAffineTemporalModeCandidate.load_json(candidate_path)
    if loaded.sha256 != capsule["candidate_sha256"]:
        raise RuntimeError("written candidate no longer matches the delivery capsule identity")
    return {"candidate": candidate_path, "capsule": capsule_path}
