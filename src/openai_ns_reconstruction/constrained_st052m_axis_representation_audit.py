"""Conditional axis-regularity audit for the integrated ST052-M transforms.

This audit distinguishes numerical axis safety from an analytic smoothness
claim.  The migrated ST052 transforms avoid division by zero at the axis and
preserve transverse O(r) scaling for an axis-regular parent, but they use a
finite ``r > 1e-14`` implementation guard.  That guard is autonomous numerical
machinery and is not, by itself, a C1/C-infinity regularity proof across the
finite-radius switch surface.

No candidate, pressure, forcing, coefficient, or scientific threshold is
changed here.
"""
from __future__ import annotations

from copy import deepcopy
import inspect
import json
from pathlib import Path
from typing import Any, Callable

import numpy as np

from . import st052_linear_temporal_capsule as capsule
from . import st052_linear_temporal_transform as temporal
from . import st052_local_swirl_transform as static

ROOT = Path(__file__).resolve().parents[2]
CONTRACT_PATH = ROOT / "configs" / "st052m_axis_representation_contract.json"

_FALSE_CLAIMS = (
    "velocity_formula_changed",
    "candidate_bytes_changed",
    "parent_evaluator_changed",
    "pressure_or_force_changed",
    "scientific_threshold_changed",
    "production_candidate_selected",
    "visual_correspondence_verified",
    "held_out_temporal_child_pde_residual_evaluated",
    "pde_validated",
    "paper_exact",
    "openai_field_identified",
    "blowup_proved",
)
_EXPECTED_CLASSES = {
    "user_requirement",
    "public_source_fact",
    "autonomous_design",
    "pending_unknown",
}

BaseVelocity = Callable[[np.ndarray, float], np.ndarray]


def _read_json(path: Path) -> dict[str, Any]:
    obj = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(obj, dict):
        raise ValueError(f"expected JSON object: {path}")
    return obj


def load_contract(path: str | Path = CONTRACT_PATH) -> dict[str, Any]:
    return _read_json(Path(path))


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def validate_contract(contract: dict[str, Any]) -> None:
    _require(contract.get("schema_version") == 1, "unsupported axis-audit schema")
    _require(
        contract.get("task_id") == "CR002-ST052M-AXIS-REPRESENTATION-073",
        "axis-audit task identity drifted",
    )
    axis = contract.get("axis_representation")
    _require(isinstance(axis, dict), "axis representation section missing")
    _require(axis.get("exact_axis_division_safe") is True, "axis division-safety fact missing")
    _require(
        axis.get("conditional_transverse_O_r_preserved") is True,
        "conditional O(r) preservation fact missing",
    )
    _require(axis.get("static_radius_guard") == 1.0e-14, "static radius guard drifted")
    _require(axis.get("hard_radius_guard_present") is True, "hard radius guard must be recorded")
    _require(
        axis.get("hard_radius_guard_is_autonomous_numerical_design") is True,
        "radius guard must remain autonomous design",
    )
    for key in (
        "analytic_C1_or_smoother_axis_regularity_proved",
        "whole_parent_axis_regularity_independently_proved_here",
        "whole_candidate_axis_regularity_theorem_proved",
    ):
        _require(axis.get(key) is False, f"premature analytic axis claim: {key}")

    states = contract.get("claim_states")
    _require(isinstance(states, dict), "claim state section missing")
    for key in _FALSE_CLAIMS:
        _require(states.get(key) is False, f"premature claim promotion: {key}")

    source_classes = contract.get("source_classification")
    _require(
        isinstance(source_classes, list) and len(source_classes) == 4,
        "source classification must contain four governed entries",
    )
    classes = {
        entry.get("classification")
        for entry in source_classes
        if isinstance(entry, dict)
    }
    _require(classes == _EXPECTED_CLASSES, "source classification classes drifted")


def _audit_cr001(contract: dict[str, Any], constraints: dict[str, Any]) -> None:
    expected = contract["cr001_nonmutation"]
    domain = constraints["domain"]
    forcing = constraints["forcing"]
    nontriviality = constraints["nontriviality"]
    validation = constraints["validation"]
    thresholds = validation["thresholds"]

    _require(constraints.get("nu") == expected["nu"], "nu drifted")
    _require(domain.get("physical") == expected["physical_domain"], "physical domain drifted")
    _require(domain.get("evaluation_box") == expected["evaluation_box"], "evaluation box drifted")
    _require(domain.get("support") == expected["support"], "support contract drifted")
    _require(domain.get("time_interval") == expected["time_interval"], "time interval drifted")
    _require(forcing.get("mode") == expected["forcing_mode"], "forcing mode drifted")
    _require(
        forcing["parameters"].get("a") == expected["forcing_parameter_bounds"]["a"],
        "forcing a bounds drifted",
    )
    _require(
        forcing["parameters"].get("c") == expected["forcing_parameter_bounds"]["c"],
        "forcing c bounds drifted",
    )
    _require(
        "No residual-dependent basis or pointwise free force" in forcing.get("restriction", ""),
        "free residual-defined forcing prohibition missing",
    )
    _require(nontriviality.get("reference_energy") == expected["reference_energy"], "reference energy drifted")
    _require(
        nontriviality.get("reference_energy_abs_tolerance")
        == expected["reference_energy_abs_tolerance"],
        "reference energy tolerance drifted",
    )
    _require(validation.get("seed") == expected["validation_seed"], "validation seed drifted")
    _require(validation.get("held_out_points") == expected["held_out_points"], "held-out count drifted")
    _require(validation.get("times") == expected["validation_times"], "validation times drifted")
    _require(validation.get("derivative_steps") == expected["derivative_steps"], "derivative ladder drifted")
    _require(
        validation.get("quadrature_orders_per_axis") == expected["quadrature_orders_per_axis"],
        "quadrature ladder drifted",
    )
    for key in ("divergence_max", "divergence_L2", "pde_residual_max", "pde_residual_L2"):
        _require(thresholds.get(key) == expected[key], f"CR001 threshold drifted: {key}")


def _regular_parent(points: np.ndarray, time: float) -> np.ndarray:
    """Smooth manufactured parent with transverse velocity exactly O(r)."""
    points = np.asarray(points, dtype=float)
    x, y, z = points[:, 0], points[:, 1], points[:, 2]
    time = float(time)
    contraction = 0.20 + 0.05 * time
    swirl = 0.70 + 0.10 * time
    ux = -contraction * x - swirl * y
    uy = -contraction * y + swirl * x
    uz = 0.30 + 0.05 * z + 0.02 * (x * x + y * y) + 0.01 * time
    return np.column_stack((ux, uy, uz))


def _irregular_parent(points: np.ndarray, time: float) -> np.ndarray:
    """Manufactured counterexample that violates the parent axis prerequisite."""
    points = np.asarray(points, dtype=float)
    return np.column_stack(
        (
            np.ones(len(points), dtype=float),
            np.zeros(len(points), dtype=float),
            np.full(len(points), 0.25 + 0.0 * float(time), dtype=float),
        )
    )


def check_parent_axis_precondition(
    parent: BaseVelocity,
    *,
    ratio_upper: float,
) -> dict[str, float | bool]:
    """Check the explicit prerequisite used by this conditional transform audit."""
    axis_points = np.array([[0.0, 0.0, z] for z in (0.0, 0.8, 1.2)], dtype=float)
    near_points = np.array(
        [[r, 0.0, z] for z in (0.0, 0.8, 1.2) for r in (1.0e-12, 1.0e-9, 1.0e-6)],
        dtype=float,
    )
    max_axis_transverse = 0.0
    max_ratio = 0.0
    for time in (0.25, 0.5, 0.75):
        axis_value = np.asarray(parent(axis_points, time), dtype=float)
        near_value = np.asarray(parent(near_points, time), dtype=float)
        _require(axis_value.shape == axis_points.shape, "parent axis shape mismatch")
        _require(near_value.shape == near_points.shape, "parent near-axis shape mismatch")
        _require(np.isfinite(axis_value).all() and np.isfinite(near_value).all(), "parent returned non-finite values")
        axis_transverse = np.hypot(axis_value[:, 0], axis_value[:, 1])
        max_axis_transverse = max(max_axis_transverse, float(np.max(axis_transverse)))
        radii = np.hypot(near_points[:, 0], near_points[:, 1])
        ratios = np.hypot(near_value[:, 0], near_value[:, 1]) / radii
        max_ratio = max(max_ratio, float(np.max(ratios)))
    _require(max_axis_transverse <= 1.0e-14, "parent transverse velocity is nonzero on axis")
    _require(max_ratio <= float(ratio_upper), "parent transverse O(r) engineering bound failed")
    return {
        "passed": True,
        "max_axis_transverse": max_axis_transverse,
        "max_transverse_over_r": max_ratio,
    }


def _evaluate_static(points: np.ndarray, time: float) -> np.ndarray:
    return static.transformed_velocity_points(_regular_parent, points, time, static.DEFAULT_SPEC)


def _evaluate_temporal(points: np.ndarray, time: float) -> np.ndarray:
    return temporal.temporal_transformed_velocity_points(
        _regular_parent, points, time, temporal.DEFAULT_TEMPORAL_SPEC
    )


def _synthetic_axis_receipt(contract: dict[str, Any]) -> dict[str, Any]:
    spec = contract["synthetic_regression"]
    ratio_upper = float(spec["engineering_transverse_over_r_upper"])
    parent_receipt = check_parent_axis_precondition(_regular_parent, ratio_upper=ratio_upper)

    max_ratio: dict[str, float] = {"static": 0.0, "temporal": 0.0}
    max_axis_transverse: dict[str, float] = {"static": 0.0, "temporal": 0.0}
    evaluators = {"static": _evaluate_static, "temporal": _evaluate_temporal}
    for label, evaluator in evaluators.items():
        for time in spec["times"]:
            for z in spec["z_values"]:
                for angle in spec["angles"]:
                    ca, sa = float(np.cos(angle)), float(np.sin(angle))
                    for radius in spec["radii"]:
                        radius = float(radius)
                        point = np.array([[radius * ca, radius * sa, float(z)]], dtype=float)
                        value = np.asarray(evaluator(point, float(time)), dtype=float)
                        _require(value.shape == (1, 3), f"{label} transformed shape mismatch")
                        _require(np.isfinite(value).all(), f"{label} transform returned non-finite values")
                        transverse = float(np.hypot(value[0, 0], value[0, 1]))
                        if radius == 0.0:
                            max_axis_transverse[label] = max(max_axis_transverse[label], transverse)
                        else:
                            max_ratio[label] = max(max_ratio[label], transverse / radius)
        _require(max_axis_transverse[label] <= 1.0e-14, f"{label} transform is not finite/zero-transverse on axis")
        _require(max_ratio[label] <= ratio_upper, f"{label} transform lost conditional O(r) scaling")

    z = float(spec["guard_probe_z"])
    time = float(spec["guard_probe_time"])
    r_below = float(spec["guard_probe_r_below"])
    r_above = float(spec["guard_probe_r_above"])
    probe = np.array([[r_below, 0.0, z], [r_above, 0.0, z]], dtype=float)
    static_values = _evaluate_static(probe, time)
    temporal_values = _evaluate_temporal(probe, time)
    static_tangential_over_r = static_values[:, 1] / probe[:, 0]
    temporal_tangential_over_r = temporal_values[:, 1] / probe[:, 0]
    static_guard_effect = float(abs(static_tangential_over_r[1] - static_tangential_over_r[0]))
    temporal_guard_effect = float(abs(temporal_tangential_over_r[1] - temporal_tangential_over_r[0]))
    minimum = float(spec["guard_effect_minimum"])
    _require(static_guard_effect > minimum, "static finite-radius branch effect was not observed; re-audit implementation")
    _require(temporal_guard_effect > minimum, "temporal finite-radius branch effect was not observed; re-audit implementation")

    return {
        "parent_precondition": parent_receipt,
        "max_axis_transverse": max_axis_transverse,
        "max_transverse_over_r": max_ratio,
        "static_guard_effect_in_tangential_over_r": static_guard_effect,
        "temporal_guard_effect_in_tangential_over_r": temporal_guard_effect,
        "conditional_transverse_O_r_preserved": True,
        "analytic_C1_or_smoother_axis_regularity_proved": False,
    }


def _audit_live_implementation() -> None:
    redistribution_source = inspect.getsource(static.redistributed_control)
    shoulder_source = inspect.getsource(static._apply_shoulder_swirl)
    taper_source = inspect.getsource(static.transformed_velocity_points)
    temporal_source = inspect.getsource(temporal.effective_static_spec)

    _require("mask = radius > 1.0e-14" in redistribution_source, "redistribution radius guard changed")
    _require("mask = radius > 1.0e-14" in shoulder_source, "shoulder radius guard changed")
    _require("x * uz" in taper_source and "y * uz" in taper_source, "taper axis-limit structure changed")
    _require("replace(" in temporal_source, "temporal transform no longer reuses static kernel")
    _require("taper_tau=" in temporal_source and "shoulder_beta=" in temporal_source, "temporal activation parameters changed")

    _require(capsule.TRUTH_BOUNDARY.get("whole_child_bundle_materialized") is True, "whole-child materialization fact disappeared")
    _require(capsule.TRUTH_BOUNDARY.get("pde_validated") is False, "ST052 PDE state unexpectedly promoted")
    _require(capsule.TRUTH_BOUNDARY.get("visual_correspondence_verified") is False, "ST052 visual correspondence unexpectedly promoted")


def audit_repository(
    *,
    root: str | Path = ROOT,
    contract: dict[str, Any] | None = None,
) -> dict[str, Any]:
    root = Path(root)
    contract = load_contract(root / "configs" / CONTRACT_PATH.name) if contract is None else contract
    validate_contract(contract)
    _audit_live_implementation()
    constraints = _read_json(root / "configs" / "constraints.json")
    _audit_cr001(contract, constraints)
    synthetic = _synthetic_axis_receipt(contract)
    return {
        "task_id": contract["task_id"],
        "integration_base": contract["integration_base"],
        "exact_axis_division_safe": True,
        "conditional_transverse_O_r_preserved": True,
        "hard_radius_guard_present": True,
        "hard_radius_guard": contract["axis_representation"]["static_radius_guard"],
        "analytic_C1_or_smoother_axis_regularity_proved": False,
        "whole_parent_axis_regularity_independently_proved_here": False,
        "whole_candidate_axis_regularity_theorem_proved": False,
        "synthetic_regression": synthetic,
        "pde_validated": False,
        "visual_correspondence_verified": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "cr001_unchanged": True,
    }


def mutation_rejects_irregular_parent() -> bool:
    """Expose the conditional nature of the audit as an executable negative control."""
    contract = load_contract()
    ratio_upper = float(contract["synthetic_regression"]["engineering_transverse_over_r_upper"])
    try:
        check_parent_axis_precondition(_irregular_parent, ratio_upper=ratio_upper)
    except ValueError:
        return True
    return False


def main() -> None:
    print(json.dumps(audit_repository(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
