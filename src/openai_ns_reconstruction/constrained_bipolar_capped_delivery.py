"""Replayable delivery capsule for the energy-defect-capped bipolar Eq45 candidate.

This module freezes provenance and public-delivery semantics for the already
optimized candidate. It does not rerun optimization or independent PDE
validation and cannot promote scientific readiness states.
"""
from __future__ import annotations

import hashlib
import inspect
import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from . import constrained_bipolar_capped_fresh as fresh_validation
from . import constrained_bipolar_joint_energy as joint_energy
from .eq45_supported_delivery import Eq45SupportedDeliveryField


SCHEMA = "eq45_bipolar_joint_capped_delivery_recipe_v1"
CAPSULE_SCHEMA = "eq45_bipolar_joint_capped_delivery_capsule_v1"
TASK_ID = "CR011-BIPOLAR-JOINT-CAPPED-DELIVERY-CAPSULE-020"
DEPENDENCY_HEAD = "c358a66d3ba74067c16010be8399fc01eea26453"
EXPECTED_CANDIDATE_SHA256 = (
    "7875214ad65f5e1ba4f7c5e362217f05ef47bc3d893e8641e8e40a85fc4e7609"
)
EXPECTED_SUPPORTED_PARENT_SHA256 = (
    "c0e27269adfb6f305c2c0b5a2483f9eddd54f79691f69cacd7bb735ba592a702"
)
EXPECTED_EQ45_PARENT_SHA256 = (
    "39ae3141214e5b91a2df39a1a277610b19afe424d681ccb9d84ff031e151a81d"
)
EXPECTED_OPTIMIZER_SEED = 9172609
EXPECTED_FRESH_SEED = 9172631
EXPECTED_FRESH_POINTS = 2048
EXPECTED_ENERGY_CAP = 0.14
EXPECTED_REFERENCE_TIMES = (0.25, 0.5, 0.75)

_ROOT = Path(__file__).resolve().parents[2]
_RECIPE_PATH = _ROOT / "artifacts/bipolar_joint_capped/delivery_recipe.json"
_CANDIDATE_PATH = _ROOT / "artifacts/bipolar_joint_capped/candidate.json"
_TRAINING_REPORT_PATH = _ROOT / "artifacts/bipolar_joint_capped/report.json"
_CONFIG_PATH = _ROOT / "configs/constraints.json"

_EXPECTED_STATUS = {
    "callable_serializable": True,
    "velocity_export_ready": True,
    "visualization_candidate_only": True,
    "candidate_selection_resolved": False,
    "physical_support_connection_implemented": True,
    "physical_support_validated": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _read_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text())
    if not isinstance(payload, dict):
        raise ValueError(f"{path} must contain a JSON object")
    return payload


def _canonical_sha256(payload: Mapping[str, Any]) -> str:
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _require_equal(actual: Any, expected: Any, message: str) -> None:
    if actual != expected:
        raise ValueError(message)


def validate_recipe(payload: Mapping[str, Any]) -> dict[str, Any]:
    """Validate provenance, registered thresholds, and independent status flags."""
    if not isinstance(payload, Mapping):
        raise ValueError("delivery recipe must be an object")
    data = dict(payload)
    _require_equal(data.get("schema"), SCHEMA, "delivery recipe schema drifted")
    _require_equal(data.get("task_id"), TASK_ID, "delivery task identity drifted")

    dependency = data.get("base_dependency")
    if not isinstance(dependency, Mapping):
        raise ValueError("base_dependency must be an object")
    _require_equal(dependency.get("pr"), 206, "fresh-validation PR identity drifted")
    _require_equal(
        dependency.get("head"), DEPENDENCY_HEAD, "dependency head identity drifted"
    )

    candidate = data.get("candidate")
    if not isinstance(candidate, Mapping):
        raise ValueError("candidate must be an object")
    _require_equal(
        candidate.get("path"),
        "artifacts/bipolar_joint_capped/candidate.json",
        "candidate path drifted",
    )
    _require_equal(
        candidate.get("sha256"), EXPECTED_CANDIDATE_SHA256, "candidate SHA drifted"
    )
    _require_equal(
        candidate.get("supported_parent_field_sha256"),
        EXPECTED_SUPPORTED_PARENT_SHA256,
        "supported parent field identity drifted",
    )
    _require_equal(
        candidate.get("underlying_eq45_parent_sha256"),
        EXPECTED_EQ45_PARENT_SHA256,
        "underlying Eq45 parent identity drifted",
    )

    optimizer = data.get("optimizer")
    if not isinstance(optimizer, Mapping):
        raise ValueError("optimizer must be an object")
    _require_equal(
        optimizer.get("entrypoint"),
        "openai_ns_reconstruction.constrained_bipolar_joint_energy:run",
        "optimizer entrypoint drifted",
    )
    _require_equal(
        optimizer.get("multistart_seed"),
        EXPECTED_OPTIMIZER_SEED,
        "optimizer multistart seed drifted",
    )
    _require_equal(optimizer.get("start_count"), 12, "optimizer start count drifted")
    _require_equal(
        optimizer.get("selected_start_index"), 6, "selected optimizer start drifted"
    )
    _require_equal(
        optimizer.get("odd_phi03_extension"), True, "odd Phi(0,3) scope drifted"
    )
    _require_equal(
        optimizer.get("training_quadrature_order"),
        96,
        "training quadrature order drifted",
    )
    _require_equal(
        optimizer.get("training_times"),
        [0.3125, 0.5, 0.6875],
        "training times drifted",
    )
    _require_equal(
        optimizer.get("training_spatial_derivative_step"),
        0.005,
        "training derivative step drifted",
    )
    _require_equal(
        optimizer.get("theta_objective_enabled"),
        True,
        "theta objective classification drifted",
    )
    _require_equal(
        float(optimizer.get("energy_defect_cap")),
        EXPECTED_ENERGY_CAP,
        "optimization energy-defect cap drifted",
    )
    _require_equal(
        optimizer.get("energy_defect_cap_role"),
        "optimization_only_not_acceptance_threshold",
        "optimization cap was relabeled as an acceptance threshold",
    )
    _require_equal(
        optimizer.get("acceptance_threshold_changed"),
        False,
        "acceptance threshold must remain unchanged",
    )

    source = inspect.getsource(joint_energy.run)
    if "np.random.default_rng(9172609)" not in source:
        raise ValueError("optimizer multistart seed no longer matches implementation")

    training = _read_json(_TRAINING_REPORT_PATH)
    _require_equal(
        training.get("candidate_sha256"),
        EXPECTED_CANDIDATE_SHA256,
        "training report candidate SHA drifted",
    )
    _require_equal(
        training.get("parent_sha256"),
        EXPECTED_SUPPORTED_PARENT_SHA256,
        "training report parent SHA drifted",
    )
    _require_equal(
        training.get("optimization_energy_cap"),
        EXPECTED_ENERGY_CAP,
        "training report energy cap drifted",
    )
    _require_equal(
        training.get("cap_is_acceptance_threshold"),
        False,
        "training report falsely promotes optimization cap",
    )
    _require_equal(
        training.get("theta_objective_enabled"),
        True,
        "training report theta objective drifted",
    )
    _require_equal(
        training.get("odd_phi03_extension"),
        True,
        "training report odd extension drifted",
    )
    _require_equal(
        training.get("training_quadrature_order"),
        96,
        "training report quadrature order drifted",
    )
    trials = training.get("multistart_trials")
    if not isinstance(trials, list) or len(trials) != 12:
        raise ValueError("training report must retain all 12 multistart trials")
    selected = trials[6]
    if selected.get("optimizer_success") is not True:
        raise ValueError("selected optimizer start is no longer successful")
    if selected.get("physical_constraints_satisfied") is not True:
        raise ValueError("selected optimizer start lost its training constraints")
    if not np.allclose(
        np.asarray(selected.get("parameters"), dtype=float),
        np.asarray(training.get("parameters"), dtype=float),
        rtol=0.0,
        atol=5e-10,
    ):
        raise ValueError("selected start no longer reconstructs reported parameters")
    if training.get("pde_validated") is not False:
        raise ValueError("training report must not claim PDE validation")

    config = _read_json(_CONFIG_PATH)
    contract = data.get("validation_contract")
    if not isinstance(contract, Mapping):
        raise ValueError("validation_contract must be an object")
    _require_equal(
        contract.get("registered_seed"),
        int(config["validation"]["seed"]),
        "registered validation seed drifted",
    )
    _require_equal(
        contract.get("registered_held_out_points"),
        int(config["validation"]["held_out_points"]),
        "registered validation sample count drifted",
    )
    _require_equal(
        contract.get("times"),
        [float(value) for value in config["validation"]["times"]],
        "registered validation times drifted",
    )
    _require_equal(
        contract.get("derivative_steps"),
        [float(value) for value in config["validation"]["derivative_steps"]],
        "registered derivative ladder drifted",
    )
    expected_thresholds = {
        key: float(config["validation"]["thresholds"][key])
        for key in (
            "divergence_max",
            "divergence_L2",
            "pde_residual_max",
            "pde_residual_L2",
        )
    }
    _require_equal(
        contract.get("thresholds"),
        expected_thresholds,
        "registered PDE/divergence thresholds drifted",
    )

    receipt = data.get("fresh_validation_receipt")
    if not isinstance(receipt, Mapping):
        raise ValueError("fresh_validation_receipt must be an object")
    _require_equal(receipt.get("pr"), 206, "fresh-validation PR drifted")
    _require_equal(receipt.get("head"), DEPENDENCY_HEAD, "fresh-validation head drifted")
    _require_equal(
        receipt.get("seed"), EXPECTED_FRESH_SEED, "fresh validation seed drifted"
    )
    _require_equal(
        receipt.get("held_out_points"),
        EXPECTED_FRESH_POINTS,
        "fresh validation point count drifted",
    )
    _require_equal(
        tuple(receipt.get("derivative_steps", ())),
        tuple(float(value) for value in fresh_validation.DEFAULT_STEPS),
        "fresh derivative ladder drifted",
    )
    _require_equal(
        tuple(receipt.get("times", ())),
        tuple(float(value) for value in fresh_validation.DEFAULT_TIMES),
        "fresh validation times drifted",
    )
    _require_equal(
        fresh_validation.DEFAULT_SEED,
        EXPECTED_FRESH_SEED,
        "fresh validation implementation seed drifted",
    )
    if receipt.get("fresh_sample_divergence_thresholds_passed") is not True:
        raise ValueError("fresh divergence evidence was relabeled")
    if receipt.get("formal_registered_divergence_gate_assessed") is not False:
        raise ValueError("2048-point fresh audit cannot masquerade as registered gate")
    if receipt.get("necessary_axisymmetric_pressure_PDE_condition_passed") is not False:
        raise ValueError("theta obstruction failure was relabeled")
    if receipt.get("formal_full_momentum_gate_assessed") is not False:
        raise ValueError("fresh audit cannot masquerade as formal full momentum gate")
    if receipt.get("pde_validated") is not False:
        raise ValueError("fresh diagnostic cannot promote PDE validation")

    status = data.get("status")
    if status != _EXPECTED_STATUS:
        raise ValueError("delivery/scientific status flags drifted")

    sibling = data.get("sibling_evidence")
    if not isinstance(sibling, Mapping):
        raise ValueError("sibling_evidence must be an object")
    _require_equal(
        sibling.get("spatial_redistribution_pr"), 207, "spatial sibling PR drifted"
    )
    _require_equal(
        sibling.get("spatial_redistribution_status"),
        "open_unconsumed",
        "unmerged spatial evidence was falsely consumed",
    )
    if sibling.get("may_promote_visualization_ready") is not False:
        raise ValueError("target-free sibling evidence cannot promote visualization")
    if sibling.get("may_promote_pde_validated") is not False:
        raise ValueError("spatial redistribution cannot promote PDE validation")

    field = Eq45SupportedDeliveryField.load_candidate(_CANDIDATE_PATH)
    _require_equal(field.sha256, EXPECTED_CANDIDATE_SHA256, "runtime candidate SHA drifted")
    _require_equal(
        field.candidate.parent_sha256,
        EXPECTED_EQ45_PARENT_SHA256,
        "runtime underlying Eq45 parent SHA drifted",
    )
    truth = field.metadata()["truth_boundary"]
    if truth.get("velocity_export_ready") is not True:
        raise ValueError("candidate unexpectedly lost export readiness")
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
            raise ValueError(f"candidate truth boundary unexpectedly promoted {key}")

    return data


def load_recipe() -> dict[str, Any]:
    return validate_recipe(_read_json(_RECIPE_PATH))


def load_field() -> Eq45SupportedDeliveryField:
    load_recipe()
    field = Eq45SupportedDeliveryField.load_candidate(_CANDIDATE_PATH)
    if field.sha256 != EXPECTED_CANDIDATE_SHA256:
        raise ValueError("runtime candidate identity changed after validation")
    return field


def write_delivery_bundle(directory: str | Path) -> dict[str, Path]:
    """Write candidate.json + capsule.json and verify exact public-velocity replay."""
    recipe = load_recipe()
    source = load_field()
    out = Path(directory)
    out.mkdir(parents=True, exist_ok=True)
    candidate_path = out / "candidate.json"
    capsule_path = out / "capsule.json"
    source.save_candidate(candidate_path)
    replay = Eq45SupportedDeliveryField.load_candidate(candidate_path)
    if replay.sha256 != EXPECTED_CANDIDATE_SHA256:
        raise ValueError("serialized candidate SHA does not replay")

    probes = np.array(
        [
            [0.12, 0.00, 0.15],
            [0.45, -0.25, 0.30],
            [-0.60, 0.35, -0.40],
            [1.10, 0.20, 0.75],
            [0.30, 0.40, -1.10],
        ],
        dtype=float,
    )
    original = np.stack(
        [source.at_points(probes, time) for time in EXPECTED_REFERENCE_TIMES], axis=0
    )
    restored = np.stack(
        [replay.at_points(probes, time) for time in EXPECTED_REFERENCE_TIMES], axis=0
    )
    if not np.array_equal(original, restored):
        raise ValueError("serialized candidate changed public off-grid [u,v,w] output")

    grid = replay.grid(
        np.array([-0.5, 0.5]),
        np.array([-0.5, 0.5]),
        np.array([-0.5, 0.5]),
        EXPECTED_REFERENCE_TIMES,
    )
    if grid.shape != (3, 2, 2, 2, 3) or not np.all(np.isfinite(grid)):
        raise ValueError("public grid replay failed")

    fingerprint = hashlib.sha256(
        np.ascontiguousarray(restored, dtype="<f8").tobytes()
    ).hexdigest()
    capsule = {
        "schema": CAPSULE_SCHEMA,
        "task_id": TASK_ID,
        "recipe_sha256": _canonical_sha256(recipe),
        "candidate_sha256": replay.sha256,
        "candidate_file": candidate_path.name,
        "public_probe_fingerprint_sha256": fingerprint,
        "reference_times": list(EXPECTED_REFERENCE_TIMES),
        "grid_layout": "time,x,y,z,component",
        "status": dict(_EXPECTED_STATUS),
        "fresh_validation_receipt": dict(recipe["fresh_validation_receipt"]),
    }
    capsule_path.write_text(json.dumps(capsule, indent=2) + "\n")
    return {"candidate": candidate_path, "capsule": capsule_path}
