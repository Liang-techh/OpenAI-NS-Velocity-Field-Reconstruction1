"""Independent numerical audit of exact A1 #1154 relative-swirl compensator algebra.

This validator is deliberately implementation-distinct from the construction path:

* the candidate configuration is save/reloaded before scientific use;
* the declared autonomous C-infinity bump is reimplemented locally;
* the angular / pressure moments are recomputed with composite Simpson on a
  frozen three-resolution ladder rather than the production GL96 path;
* the small branch is reconstructed by continuation + Newton rather than the
  production quadratic-formula root solver.

The result is scoped algebra/operator-consistency evidence only.  #1154 does
not materialize the current-lineage entering I/r_I, does not compose the
relative-swirl edit into Cartesian velocity, and has no matched pressure,
restricted forcing, or complete Navier--Stokes defect.  Therefore this module
cannot promote ``pde_validated``.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
import tempfile
from typing import Any, Mapping

import numpy as np

from .kokuno_public_relative_swirl_compensator import (
    KokunoPublicRelativeSwirlCompensator,
)

TASK = "K4-VAL-120"
SCHEMA = "kokuno-a4-relative-swirl-compensator-independent-audit-v1"
UPSTREAM_PR = 1154
UPSTREAM_HEAD = "d0e855e37c8f5d6d58f1a85b1dcb874651d20838"
UPSTREAM_SOURCE_BLOB = "7715dbf9fadb5d9ed37213677703b45f21d29724"
EXPECTED_PARENT_HEAD = "5bafa4199ef127bf50bce30a58cc69ded99e1aa9"
EXPECTED_SCHEMA = "kokuno-public-relative-swirl-compensator-v1"
EXPECTED_LAMBDA = 0.05
EXPECTED_HOLD_MULTIPLIER = 30.0
EXPECTED_WIDTH = 0.3
EXPECTED_FIRST_CENTER_FROM_END = 3.0
EXPECTED_SECOND_CENTER_FROM_END = 1.0
EXPECTED_PRODUCTION_QUADRATURE = 96

SIMPSON_NODES = (2049, 4097, 8193)
ANGULAR_TARGETS = (-0.0187, -0.0113, -0.0049, 0.0, 0.0061, 0.0137, 0.0193)
CONTINUATION_STEPS = 32
NEWTON_MAX_STEPS = 24
NEWTON_TOL = 2.0e-14

ROW_RELATIVE_MAX_GATE = 2.0e-6
COEFFICIENT_RELATIVE_RMS_GATE = 5.0e-5
COEFFICIENT_RELATIVE_MAX_GATE = 2.0e-4
INDEPENDENT_CLOSURE_GATE = 2.0e-10
SAVE_RELOAD_COEFFICIENT_GATE = 2.0e-12
NONZERO_COEFFICIENT_L2_FLOOR = 1.0e-5
STABILITY_FACTOR = 1.25
STABILITY_FLOOR = 2.0e-10

FINAL_PROJECT_MOMENTUM_GATE = 1.0e-3
FINAL_PROJECT_DIVERGENCE_GATE = 1.0e-5

_TRUTH_BOUNDARY = {
    "exact_a1_1154_relative_swirl_algebra_consumed": True,
    "candidate_configuration_save_reload_required": True,
    "independent_bump_reimplementation_used": True,
    "independent_composite_simpson_used": True,
    "independent_continuation_newton_branch_solver_used": True,
    "relative_swirl_algebra_scoped_assessed": True,
    "current_lineage_angular_entry_I_materialized": False,
    "current_cartesian_relative_swirl_composed": False,
    "complete_terminal_global_leading_materialized": False,
    "leading_only_ns_residual_assessed": False,
    "leading_plus_oscillatory_ns_residual_assessed": False,
    "after_correction_ns_residual_assessed": False,
    "matched_global_pressure_materialized": False,
    "restricted_forcing_materialized": False,
    "complete_ns_residual_assessed": False,
    "same_protocol_comparable_to_st006": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def default_compensator() -> KokunoPublicRelativeSwirlCompensator:
    return KokunoPublicRelativeSwirlCompensator()


def _canonical_json(payload: Mapping[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _assert_exact_contract(candidate: Any) -> None:
    config = candidate.configuration()
    checks = {
        "schema": EXPECTED_SCHEMA,
        "parent_exact_head": EXPECTED_PARENT_HEAD,
        "lambda_value": EXPECTED_LAMBDA,
        "hold_multiplier": EXPECTED_HOLD_MULTIPLIER,
        "bump_width": EXPECTED_WIDTH,
        "first_center_from_end": EXPECTED_FIRST_CENTER_FROM_END,
        "second_center_from_end": EXPECTED_SECOND_CENTER_FROM_END,
        "quadrature_order": EXPECTED_PRODUCTION_QUADRATURE,
    }
    for key, expected in checks.items():
        if config.get(key) != expected:
            raise ValueError(f"exact #1154 contract drift in {key}")
    truth = config.get("truth_boundary")
    if not isinstance(truth, dict):
        raise ValueError("missing #1154 truth boundary")
    for key in (
        "current_lineage_angular_entry_I_materialized",
        "current_cartesian_relative_swirl_composed",
        "matched_global_pressure_materialized",
        "restricted_forcing_materialized",
        "heldout_ns_residual_assessed",
        "pde_validated",
    ):
        if truth.get(key) is not False:
            raise ValueError(f"#1154 truth boundary unexpectedly promotes {key}")
    if not math.isclose(
        float(candidate.hold_length),
        EXPECTED_HOLD_MULTIPLIER * math.log(1.0 / EXPECTED_LAMBDA),
        rel_tol=0.0,
        abs_tol=2.0e-14,
    ):
        raise ValueError("hold length drift")
    if not math.isclose(
        float(candidate.y1),
        float(candidate.hold_length) - EXPECTED_FIRST_CENTER_FROM_END,
        rel_tol=0.0,
        abs_tol=2.0e-14,
    ):
        raise ValueError("first bump center drift")
    if not math.isclose(
        float(candidate.y2),
        float(candidate.hold_length) - EXPECTED_SECOND_CENTER_FROM_END,
        rel_tol=0.0,
        abs_tol=2.0e-14,
    ):
        raise ValueError("second bump center drift")


def _a4_bump(y: Any, center: float, width: float) -> np.ndarray:
    """Independent implementation of the declared autonomous bump realization."""
    yy = np.asarray(y, dtype=float)
    if np.any(~np.isfinite(yy)):
        raise ValueError("nonfinite independent bump coordinates")
    z = 2.0 * (yy - float(center)) / float(width)
    out = np.zeros_like(z, dtype=float)
    inside = np.abs(z) < 1.0
    if np.any(inside):
        q = z[inside]
        out[inside] = np.exp(1.0 - 1.0 / (1.0 - q * q))
    return out


def _simpson_integral(func: Any, a: float, b: float, nodes: int) -> float:
    if nodes < 3 or nodes % 2 != 1:
        raise ValueError("composite Simpson requires an odd node count >=3")
    if not (math.isfinite(a) and math.isfinite(b) and b >= a):
        raise ValueError("invalid Simpson bounds")
    if a == b:
        return 0.0
    x = np.linspace(a, b, nodes, dtype=float)
    values = np.asarray(func(x), dtype=float)
    if values.shape != x.shape or np.any(~np.isfinite(values)):
        raise RuntimeError("invalid independent quadrature integrand")
    h = (b - a) / float(nodes - 1)
    return float(
        (h / 3.0)
        * (
            values[0]
            + values[-1]
            + 4.0 * np.sum(values[1:-1:2])
            + 2.0 * np.sum(values[2:-2:2])
        )
    )


def _reference_rows(candidate: Any, nodes: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    _assert_exact_contract(candidate)
    y1 = float(candidate.y1)
    y2 = float(candidate.y2)
    width = float(candidate.width)
    angular_slope = 1.0 - EXPECTED_LAMBDA
    pressure_slope = -1.0 - 2.0 * EXPECTED_LAMBDA
    centers = (y1, y2)
    half = 0.5 * width

    angular = np.empty(2, dtype=float)
    pressure = np.empty(2, dtype=float)
    for j, center in enumerate(centers):
        lo = center - half
        hi = center + half
        angular[j] = _simpson_integral(
            lambda y, c=center: (
                np.exp(angular_slope * (y - y1)) * _a4_bump(y, c, width)
            ),
            lo,
            hi,
            nodes,
        )
        pressure[j] = _simpson_integral(
            lambda y, c=center: (
                np.exp(pressure_slope * (y - y1)) * _a4_bump(y, c, width)
            ),
            lo,
            hi,
            nodes,
        )

    quadratic = np.zeros((2, 2), dtype=float)
    for i, ci in enumerate(centers):
        for j, cj in enumerate(centers):
            lo = max(ci - half, cj - half)
            hi = min(ci + half, cj + half)
            if hi <= lo:
                continue
            quadratic[i, j] = _simpson_integral(
                lambda y, a=ci, b=cj: (
                    np.exp(pressure_slope * (y - y1))
                    * _a4_bump(y, a, width)
                    * _a4_bump(y, b, width)
                ),
                lo,
                hi,
                nodes,
            )
    return angular, pressure, quadratic


def _production_rows(candidate: Any) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    angular = np.asarray(candidate.angular_row_scaled, dtype=float)
    pressure = np.asarray(candidate.pressure_linear_row_scaled, dtype=float)
    quadratic = np.asarray(candidate.pressure_quadratic_matrix_scaled, dtype=float)
    if angular.shape != (2,) or pressure.shape != (2,) or quadratic.shape != (2, 2):
        raise RuntimeError("invalid public production row shape")
    if (
        np.any(~np.isfinite(angular))
        or np.any(~np.isfinite(pressure))
        or np.any(~np.isfinite(quadratic))
    ):
        raise RuntimeError("nonfinite public production row")
    return angular, pressure, quadratic


def _flatten_rows(rows: tuple[np.ndarray, np.ndarray, np.ndarray]) -> np.ndarray:
    a, b, q = rows
    return np.concatenate(
        (
            np.asarray(a).reshape(-1),
            np.asarray(b).reshape(-1),
            np.asarray(q).reshape(-1),
        )
    )


def _row_metrics(
    production: tuple[np.ndarray, np.ndarray, np.ndarray],
    reference: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> dict[str, float]:
    p = _flatten_rows(production)
    r = _flatten_rows(reference)
    delta = r - p
    scale = np.maximum(np.abs(p), 1.0e-12)
    return {
        "relative_max": float(np.max(np.abs(delta) / scale)),
        "relative_rms": float(np.sqrt(np.mean((delta / scale) ** 2))),
        "absolute_max": float(np.max(np.abs(delta))),
        "reference_norm_l2": float(np.linalg.norm(r)),
    }


def _system_value(
    coeff: np.ndarray,
    target: float,
    rows: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> np.ndarray:
    a, b, q = rows
    c = np.asarray(coeff, dtype=float)
    return np.asarray(
        (
            float(a @ c - target),
            float(2.0 * b @ c + c @ q @ c),
        ),
        dtype=float,
    )


def _system_jacobian(
    coeff: np.ndarray,
    rows: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> np.ndarray:
    a, b, q = rows
    c = np.asarray(coeff, dtype=float)
    return np.vstack((a, 2.0 * b + 2.0 * (q @ c)))


def _reference_solution(
    target: float,
    rows: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> np.ndarray:
    """Follow the zero-connected branch with continuation + Newton."""
    if not math.isfinite(target):
        raise ValueError("nonfinite angular target")
    coeff = np.zeros(2, dtype=float)
    if target == 0.0:
        return coeff
    for tau in np.linspace(0.0, target, CONTINUATION_STEPS + 1, dtype=float)[1:]:
        converged = False
        for _ in range(NEWTON_MAX_STEPS):
            value = _system_value(coeff, float(tau), rows)
            if float(np.max(np.abs(value))) <= NEWTON_TOL:
                converged = True
                break
            jac = _system_jacobian(coeff, rows)
            if not np.all(np.isfinite(jac)) or abs(float(np.linalg.det(jac))) <= 1.0e-14:
                raise RuntimeError("independent branch Newton Jacobian is singular")
            delta = np.linalg.solve(jac, -value)
            coeff = coeff + delta
            if np.any(~np.isfinite(coeff)):
                raise RuntimeError("nonfinite independent branch iterate")
            if float(np.max(np.abs(delta))) <= NEWTON_TOL:
                value = _system_value(coeff, float(tau), rows)
                if float(np.max(np.abs(value))) <= 5.0 * NEWTON_TOL:
                    converged = True
                    break
        if not converged:
            raise RuntimeError("independent continuation/Newton did not converge")
    return coeff


def _production_coefficients(candidate: Any) -> np.ndarray:
    target = np.asarray(ANGULAR_TARGETS, dtype=float)
    solution = candidate.solve(target)
    coeff = np.stack(
        (
            np.asarray(solution.c1, dtype=float),
            np.asarray(solution.c2, dtype=float),
        ),
        axis=-1,
    )
    if coeff.shape != (len(ANGULAR_TARGETS), 2) or np.any(~np.isfinite(coeff)):
        raise RuntimeError("invalid public production coefficients")
    return coeff


def _reference_coefficients(
    rows: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> np.ndarray:
    return np.asarray(
        [_reference_solution(float(target), rows) for target in ANGULAR_TARGETS],
        dtype=float,
    )


def _coefficient_metrics(
    production: np.ndarray,
    reference: np.ndarray,
    rows: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> dict[str, Any]:
    per_target_relative: list[float] = []
    closure_angular: list[float] = []
    closure_pressure: list[float] = []
    nonzero_norms: list[float] = []
    target_rows: list[dict[str, Any]] = []
    for idx, target in enumerate(ANGULAR_TARGETS):
        pc = np.asarray(production[idx], dtype=float)
        rc = np.asarray(reference[idx], dtype=float)
        abs_delta = float(np.linalg.norm(rc - pc))
        denom = max(float(np.linalg.norm(pc)), 1.0e-15)
        rel = abs_delta / denom if target != 0.0 else abs_delta
        if target != 0.0:
            per_target_relative.append(rel)
            nonzero_norms.append(float(np.linalg.norm(rc)))
        value = _system_value(rc, float(target), rows)
        closure_angular.append(abs(float(value[0])))
        closure_pressure.append(abs(float(value[1])))
        target_rows.append(
            {
                "target": float(target),
                "production": pc.tolist(),
                "independent": rc.tolist(),
                "coefficient_l2_absolute_difference": abs_delta,
                "coefficient_relative_difference": float(rel),
                "independent_angular_abs_closure": abs(float(value[0])),
                "independent_pressure_abs_closure": abs(float(value[1])),
                "independent_coefficient_l2": float(np.linalg.norm(rc)),
            }
        )
    rel_arr = np.asarray(per_target_relative, dtype=float)
    return {
        "relative_rms": float(np.sqrt(np.mean(rel_arr * rel_arr))),
        "relative_max": float(np.max(rel_arr)),
        "zero_target_absolute_max": float(np.max(np.abs(reference[3] - production[3]))),
        "independent_angular_abs_closure_max": float(np.max(closure_angular)),
        "independent_pressure_abs_closure_max": float(np.max(closure_pressure)),
        "independent_nonzero_coefficient_l2_min": float(np.min(nonzero_norms)),
        "targets": target_rows,
    }


def _resolution_report(candidate: Any, nodes: int) -> dict[str, Any]:
    production_rows = _production_rows(candidate)
    reference_rows = _reference_rows(candidate, nodes)
    production_coeff = _production_coefficients(candidate)
    reference_coeff = _reference_coefficients(reference_rows)
    return {
        "nodes": int(nodes),
        "row_metrics": _row_metrics(production_rows, reference_rows),
        "coefficient_metrics": _coefficient_metrics(
            production_coeff, reference_coeff, reference_rows
        ),
        "reference_rows": {
            "angular": reference_rows[0].tolist(),
            "pressure_linear": reference_rows[1].tolist(),
            "pressure_quadratic": reference_rows[2].tolist(),
        },
    }


def _replay_coefficient_error(candidate: Any, pre_serialization_reference: Any) -> float:
    a = _production_coefficients(candidate)
    b = _production_coefficients(pre_serialization_reference)
    return float(np.max(np.abs(a - b)))


def _stability_ok(fine: float, medium: float) -> bool:
    if medium <= STABILITY_FLOOR:
        return fine <= STABILITY_FLOOR
    return fine <= STABILITY_FACTOR * medium


def audit_loaded_compensator(
    candidate: Any,
    pre_serialization_reference: Any,
) -> dict[str, Any]:
    """Run the frozen K4-VAL-120 audit with no caller-tunable scientific knobs."""
    _assert_exact_contract(candidate)
    _assert_exact_contract(pre_serialization_reference)
    if candidate.semantic_sha256() != pre_serialization_reference.semantic_sha256():
        raise ValueError("save/reload semantic identity mismatch")
    if _canonical_json(candidate.configuration()) != _canonical_json(
        pre_serialization_reference.configuration()
    ):
        raise ValueError("save/reload configuration mismatch")

    resolution = [_resolution_report(candidate, n) for n in SIMPSON_NODES]
    medium = resolution[1]
    fine = resolution[2]
    replay_error = _replay_coefficient_error(candidate, pre_serialization_reference)

    fine_row = float(fine["row_metrics"]["relative_max"])
    fine_coeff_rms = float(fine["coefficient_metrics"]["relative_rms"])
    fine_coeff_max = float(fine["coefficient_metrics"]["relative_max"])
    fine_ang_closure = float(
        fine["coefficient_metrics"]["independent_angular_abs_closure_max"]
    )
    fine_pressure_closure = float(
        fine["coefficient_metrics"]["independent_pressure_abs_closure_max"]
    )
    nontrivial = float(
        fine["coefficient_metrics"]["independent_nonzero_coefficient_l2_min"]
    )
    medium_row = float(medium["row_metrics"]["relative_max"])
    medium_coeff = float(medium["coefficient_metrics"]["relative_max"])

    gates = {
        "fine_row_relative_max_le_2e-6": fine_row <= ROW_RELATIVE_MAX_GATE,
        "fine_coefficient_relative_rms_le_5e-5": (
            fine_coeff_rms <= COEFFICIENT_RELATIVE_RMS_GATE
        ),
        "fine_coefficient_relative_max_le_2e-4": (
            fine_coeff_max <= COEFFICIENT_RELATIVE_MAX_GATE
        ),
        "independent_angular_closure_le_2e-10": (
            fine_ang_closure <= INDEPENDENT_CLOSURE_GATE
        ),
        "independent_pressure_closure_le_2e-10": (
            fine_pressure_closure <= INDEPENDENT_CLOSURE_GATE
        ),
        "save_reload_coefficient_replay_le_2e-12": (
            replay_error <= SAVE_RELOAD_COEFFICIENT_GATE
        ),
        "nonzero_target_coefficient_l2_ge_1e-5": (
            nontrivial >= NONZERO_COEFFICIENT_L2_FLOOR
        ),
        "medium_to_fine_row_not_worse": _stability_ok(fine_row, medium_row),
        "medium_to_fine_coefficient_not_worse": _stability_ok(
            fine_coeff_max, medium_coeff
        ),
        "final_project_momentum_gate_unchanged": FINAL_PROJECT_MOMENTUM_GATE,
        "final_project_divergence_gate_unchanged": FINAL_PROJECT_DIVERGENCE_GATE,
    }
    boolean_gate_values = [value for value in gates.values() if isinstance(value, bool)]
    failures = [
        name for name, value in gates.items() if isinstance(value, bool) and not value
    ]

    return {
        "schema": SCHEMA,
        "task": TASK,
        "upstream_pr": UPSTREAM_PR,
        "upstream_head": UPSTREAM_HEAD,
        "upstream_source_blob": UPSTREAM_SOURCE_BLOB,
        "candidate_semantic_sha256": candidate.semantic_sha256(),
        "protocol": {
            "simpson_nodes": list(SIMPSON_NODES),
            "angular_targets": list(ANGULAR_TARGETS),
            "continuation_steps": CONTINUATION_STEPS,
            "newton_max_steps": NEWTON_MAX_STEPS,
            "newton_tolerance": NEWTON_TOL,
            "production_quadrature_excluded_from_independent_reference": True,
            "production_quadratic_formula_solver_excluded_from_independent_reference": True,
        },
        "resolution_reports": resolution,
        "save_reload_coefficient_absolute_max": replay_error,
        "gates": gates,
        "audit_pass": bool(all(boolean_gate_values)),
        "failures": failures,
        "truth_boundary": dict(_TRUTH_BOUNDARY),
    }


def materialize_receipt() -> dict[str, Any]:
    original = default_compensator()
    with tempfile.TemporaryDirectory(prefix="kokuno_a4_rel_swirl_") as td:
        path = Path(td) / "candidate.json"
        original.save_configuration(path)
        loaded = KokunoPublicRelativeSwirlCompensator.load_configuration(path)
        return audit_loaded_compensator(loaded, original)


def enforce_preregistered_gates(report: Mapping[str, Any]) -> None:
    if report.get("schema") != SCHEMA or report.get("task") != TASK:
        raise AssertionError("wrong K4-VAL-120 report identity")
    truth = report.get("truth_boundary")
    if not isinstance(truth, dict) or truth != _TRUTH_BOUNDARY:
        raise AssertionError("truth boundary drift")
    gates = report.get("gates")
    if not isinstance(gates, dict):
        raise AssertionError("missing preregistered gates")
    expected_boolean = (
        "fine_row_relative_max_le_2e-6",
        "fine_coefficient_relative_rms_le_5e-5",
        "fine_coefficient_relative_max_le_2e-4",
        "independent_angular_closure_le_2e-10",
        "independent_pressure_closure_le_2e-10",
        "save_reload_coefficient_replay_le_2e-12",
        "nonzero_target_coefficient_l2_ge_1e-5",
        "medium_to_fine_row_not_worse",
        "medium_to_fine_coefficient_not_worse",
    )
    for key in expected_boolean:
        if gates.get(key) is not True:
            raise AssertionError(f"preregistered gate failed: {key}")
    if gates.get("final_project_momentum_gate_unchanged") != FINAL_PROJECT_MOMENTUM_GATE:
        raise AssertionError("final project momentum gate drift")
    if gates.get("final_project_divergence_gate_unchanged") != FINAL_PROJECT_DIVERGENCE_GATE:
        raise AssertionError("final project divergence gate drift")
    if report.get("audit_pass") is not True or report.get("failures") != []:
        raise AssertionError("K4-VAL-120 scoped audit did not pass")


__all__ = [
    "TASK",
    "SCHEMA",
    "UPSTREAM_PR",
    "UPSTREAM_HEAD",
    "SIMPSON_NODES",
    "ANGULAR_TARGETS",
    "default_compensator",
    "audit_loaded_compensator",
    "materialize_receipt",
    "enforce_preregistered_gates",
]
