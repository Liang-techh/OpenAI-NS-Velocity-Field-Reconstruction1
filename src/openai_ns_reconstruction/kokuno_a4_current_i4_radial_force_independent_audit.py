"""Independent A4 audit of the exact current-I4 correction-side radial force.

The production Agent-3 path in PR #1118 evaluates

    (div T)_r = partial_z sigma_1

with the inherited centered-FD2 physical-z derivative ladder.  This audit sees
only JSON-round-tripped public #1118 force receipts and public #1109 shifted
axial-stress receipts.  It reconstructs the same derivative with a centered
five-point FD4 stencil and never calls the production differentiation helper.

This is scoped correction-side operator-consistency evidence only.  It is not
an RF30--RF49 five-row correction, complete NS defect, correction velocity,
pressure/forcing validation, residual-reduction claim, or PDE admission.
"""
from __future__ import annotations

import inspect
from typing import Any, Mapping, Sequence

import numpy as np

TASK = "K4-VAL-116"
SCHEMA = "kokuno-a4-current-i4-radial-force-independent-audit-v1"

PARENT_AGENT3_PR = 1118
PARENT_AGENT3_HEAD = "922f7aa10460ded44af212d313eceff33a2ac647"
PARENT_AGENT3_SOURCE_BLOB = "e28eac6f7fb7f3191a024b052bcf677cff8d74d2"
PARENT_FORCE_SCHEMA = "kokuno-a3-current-i4-nonlinear-radial-force-v1"

STRESS_AGENT3_PR = 1109
STRESS_AGENT3_HEAD = "49590ff311fef1dec4bd850b3013fd989485c87a"
STRESS_AGENT3_SOURCE_BLOB = "b8d9928934d4df74b52e8cd0d2c27c468257cdc9"
STRESS_SCHEMA = "kokuno-a3-current-i4-nonlinear-radial-stress-v1"
STRESS_PARENT_AGENT3_1101_HEAD = "6c9a4e7b0dec4b3669fb1dc60297b18f10163e0b"
STRESS_PARENT_AGENT3_1101_SOURCE_BLOB = "10d067570ddba834db3f19c9859c2efdad16b3e8"

AGENT2_1080_HEAD = "c40d8ddecd2971544a6e07dab093436b423cf326"
AGENT1_1079_HEAD = "b06742ca6e189499192ede3cce40f62cdc1e35ca"

FROZEN_SEED = 9173881
FROZEN_TIME = 0.39
FROZEN_Z_CENTERS = (-0.137, -0.071, 0.061)
FROZEN_Z_STEPS = (0.02, 0.01, 0.005)
FROZEN_R_MIN = 0.005
FROZEN_R_MAX = 0.44
FROZEN_R_COUNT = 169
FROZEN_BUMP_CENTER = 0.30
FROZEN_BUMP_HALFWIDTH = 0.10
FROZEN_AXIS_NEAR_MAX_R = 0.02

FINE_FORCE_REL_RMS_GATE = 5.0e-2
FINE_FORCE_REL_MAX_GATE = 1.5e-1
FINE_FORCE_REL_INTEGRAL_L2_GATE = 5.0e-2
AXIS_NEAR_NORMALIZED_ERROR_GATE = 1.5e-1
INDEPENDENT_FINE_PAIR_REL_RMS_GATE = 5.0e-2
RESOLUTION_DEGRADATION_FACTOR = 1.25
RESOLUTION_NUMERICAL_FLOOR = 2.0e-5
NONTRIVIAL_FORCE_RMS_FLOOR = 1.0e-12
PIECE_CLOSURE_RELATIVE_GATE = 5.0e-10

FINAL_NORMALIZED_MOMENTUM_GATE = 1.0e-3
FINAL_NORMALIZED_DIVERGENCE_GATE = 1.0e-5


def frozen_radial_grid() -> np.ndarray:
    """Return the preregistered deterministic positive-radius off-grid ladder."""
    r = np.linspace(FROZEN_R_MIN, FROZEN_R_MAX, FROZEN_R_COUNT, dtype=float)
    h = (FROZEN_R_MAX - FROZEN_R_MIN) / (FROZEN_R_COUNT - 1)
    rng = np.random.default_rng(FROZEN_SEED)
    r[1:-1] += rng.uniform(-0.22, 0.22, size=FROZEN_R_COUNT - 2) * h
    if np.any(np.diff(r) <= 0.0):
        raise RuntimeError("frozen radial jitter lost monotonicity")
    return r


def required_z_offsets() -> tuple[float, ...]:
    return tuple(
        sorted(
            {
                sign * factor * h
                for h in FROZEN_Z_STEPS
                for factor in (1.0, 2.0)
                for sign in (-1.0, 1.0)
            }
        )
    )


def _rms(values: np.ndarray) -> float:
    a = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(a * a)))


def _relative_rms(a: np.ndarray, b: np.ndarray) -> float:
    x = np.asarray(a, dtype=float)
    y = np.asarray(b, dtype=float)
    if x.shape != y.shape or not (np.all(np.isfinite(x)) and np.all(np.isfinite(y))):
        raise ValueError("arrays must have matching finite shapes")
    return _rms(x - y) / max(_rms(x), _rms(y), np.finfo(float).tiny)


def _radial_integral_l2(values: np.ndarray, radii: np.ndarray) -> float:
    v = np.asarray(values, dtype=float)
    r = np.asarray(radii, dtype=float)
    if v.shape != r.shape or v.ndim != 1 or np.any(np.diff(r) <= 0.0):
        raise ValueError("radial L2 requires matching increasing 1D arrays")
    density = r * v * v
    integral = float(np.sum(0.5 * (density[:-1] + density[1:]) * np.diff(r)))
    return float(np.sqrt(max(integral, 0.0)))


def _comparison_metrics(
    production: np.ndarray, independent: np.ndarray, radii: np.ndarray
) -> dict[str, float]:
    p = np.asarray(production, dtype=float)
    q = np.asarray(independent, dtype=float)
    if p.shape != q.shape or p.shape != radii.shape:
        raise ValueError("force arrays must match preregistered radial grid")
    if not (np.all(np.isfinite(p)) and np.all(np.isfinite(q))):
        raise ValueError("force arrays must be finite")
    delta = p - q
    rms_scale = max(_rms(p), _rms(q), np.finfo(float).tiny)
    max_scale = max(
        float(np.max(np.abs(p))), float(np.max(np.abs(q))), np.finfo(float).tiny
    )
    lp = _radial_integral_l2(p, radii)
    lq = _radial_integral_l2(q, radii)
    ld = _radial_integral_l2(delta, radii)
    axis = radii <= FROZEN_AXIS_NEAR_MAX_R
    if not np.any(axis):
        raise ValueError("preregistered grid lost positive-radius axis-near samples")
    worst = int(np.argmax(np.abs(delta)))
    return {
        "relative_rms": _rms(delta) / rms_scale,
        "relative_max": float(np.max(np.abs(delta))) / max_scale,
        "relative_integral_l2": ld / max(lp, lq, np.finfo(float).tiny),
        "axis_near_normalized_error": float(np.max(np.abs(delta[axis]))) / max_scale,
        "production_rms": _rms(p),
        "independent_rms": _rms(q),
        "production_integral_l2": lp,
        "independent_integral_l2": lq,
        "worst_radius": float(radii[worst]),
        "worst_signed_error": float(delta[worst]),
    }


def _validate_truth_boundary(boundary: Mapping[str, Any]) -> None:
    for key in (
        "complete_ns_defect",
        "current_real_ns_correction_velocity_materialized",
        "real_candidate_finite_correction_cycle_run",
        "heldout_normalized_ns_residual_assessed",
        "residual_reduction_claimed",
        "same_protocol_comparable_to_st006",
        "pde_validated",
    ):
        if boundary.get(key) is not False:
            raise ValueError(f"scientific truth boundary drifted at {key}")
    if float(boundary.get("final_normalized_momentum_gate")) != FINAL_NORMALIZED_MOMENTUM_GATE:
        raise ValueError("final momentum gate drifted")
    if float(boundary.get("final_normalized_divergence_gate")) != FINAL_NORMALIZED_DIVERGENCE_GATE:
        raise ValueError("final divergence gate drifted")


def _validate_common_geometry(
    receipt: Mapping[str, Any], *, expected_z: float, radii: np.ndarray
) -> None:
    geometry = receipt.get("geometry")
    if not isinstance(geometry, Mapping):
        raise ValueError("receipt geometry is missing")
    rr = np.asarray(geometry.get("radii"), dtype=float)
    if rr.shape != radii.shape or not np.allclose(rr, radii, rtol=0.0, atol=2e-15):
        raise ValueError("receipt radial grid drifted")
    if abs(float(geometry.get("time")) - FROZEN_TIME) > 1e-15:
        raise ValueError("receipt time drifted")
    if abs(float(geometry.get("axial_z")) - expected_z) > 2e-14:
        raise ValueError("receipt axial z drifted")
    if abs(float(geometry.get("bump_center")) - FROZEN_BUMP_CENTER) > 1e-15:
        raise ValueError("receipt bump center drifted")
    if abs(float(geometry.get("bump_halfwidth")) - FROZEN_BUMP_HALFWIDTH) > 1e-15:
        raise ValueError("receipt bump halfwidth drifted")


def _validate_force_receipt(
    receipt: Mapping[str, Any], *, center_z: float, radii: np.ndarray
) -> None:
    if receipt.get("schema") != PARENT_FORCE_SCHEMA:
        raise ValueError("unexpected A3 #1118 force receipt schema")
    if receipt.get("parent_agent3_1109_source_blob") != STRESS_AGENT3_SOURCE_BLOB:
        raise ValueError("A3 #1109 source identity drifted")
    if receipt.get("source_radial_force_formula") != "(div T)_r = partial_z sigma_1":
        raise ValueError("source radial-force formula drifted")
    if tuple(float(v) for v in receipt.get("z_derivative_step_ladder", ())) != FROZEN_Z_STEPS:
        raise ValueError("production z-derivative ladder drifted")
    provenance = receipt.get("provenance")
    if not isinstance(provenance, Mapping):
        raise ValueError("force receipt provenance is missing")
    expected = {
        "parent_agent3_head": STRESS_AGENT3_HEAD,
        "agent2_composite_head": AGENT2_1080_HEAD,
        "agent1_leading_head": AGENT1_1079_HEAD,
    }
    for key, value in expected.items():
        if provenance.get(key) != value:
            raise ValueError(f"force provenance drifted at {key}")
    boundary = receipt.get("truth_boundary")
    if not isinstance(boundary, Mapping):
        raise ValueError("force receipt truth boundary is missing")
    _validate_truth_boundary(boundary)
    _validate_common_geometry(receipt, expected_z=center_z, radii=radii)
    closure = tuple(float(v) for v in receipt.get("piece_closure_relative_max_levels", ()))
    if len(closure) != len(FROZEN_Z_STEPS) or max(closure) > PIECE_CLOSURE_RELATIVE_GATE:
        raise ValueError("production force piece attribution no longer closes")


def _extract_stress_piece(
    receipt: Mapping[str, Any], piece: str, radii: np.ndarray
) -> np.ndarray:
    axial = receipt.get("axial_e1")
    if not isinstance(axial, Mapping) or not isinstance(axial.get(piece), Mapping):
        raise ValueError(f"stress receipt missing axial_e1/{piece}")
    arr = np.asarray(axial[piece].get("stress"), dtype=float)
    if arr.shape != radii.shape or not np.all(np.isfinite(arr)):
        raise ValueError(f"invalid axial stress array for {piece}")
    return arr


def _validate_stress_receipt(
    receipt: Mapping[str, Any], *, expected_z: float, radii: np.ndarray
) -> dict[str, np.ndarray]:
    if receipt.get("schema") != STRESS_SCHEMA:
        raise ValueError("unexpected A3 #1109 stress receipt schema")
    if receipt.get("parent_agent3_1101_source_blob") != STRESS_PARENT_AGENT3_1101_SOURCE_BLOB:
        raise ValueError("A3 #1101 source identity drifted")
    provenance = receipt.get("provenance")
    if not isinstance(provenance, Mapping):
        raise ValueError("stress receipt provenance is missing")
    expected = {
        "parent_agent3_head": STRESS_PARENT_AGENT3_1101_HEAD,
        "agent2_composite_head": AGENT2_1080_HEAD,
        "agent1_leading_head": AGENT1_1079_HEAD,
    }
    for key, value in expected.items():
        if provenance.get(key) != value:
            raise ValueError(f"stress provenance drifted at {key}")
    boundary = receipt.get("truth_boundary")
    if not isinstance(boundary, Mapping):
        raise ValueError("stress receipt truth boundary is missing")
    _validate_truth_boundary(boundary)
    _validate_common_geometry(receipt, expected_z=expected_z, radii=radii)
    pieces = {
        piece: _extract_stress_piece(receipt, piece, radii)
        for piece in ("quadratic", "mixed", "aggregate")
    }
    closure = pieces["aggregate"] - pieces["mixed"] - pieces["quadratic"]
    scale = max(
        float(np.max(np.abs(pieces["aggregate"]))),
        float(np.max(np.abs(pieces["mixed"]))),
        float(np.max(np.abs(pieces["quadratic"]))),
        np.finfo(float).tiny,
    )
    if float(np.max(np.abs(closure))) / scale > 2.0e-11:
        raise ValueError("public axial stress piece attribution no longer closes")
    return pieces


def _index_stress_receipts(
    receipts: Sequence[Mapping[str, Any]], *, center_z: float, radii: np.ndarray
) -> dict[float, dict[str, np.ndarray]]:
    expected_offsets = required_z_offsets()
    if len(receipts) != len(expected_offsets):
        raise ValueError("A4 requires exactly eight shifted stress receipts per held-out center")
    indexed: dict[float, dict[str, np.ndarray]] = {}
    used: set[float] = set()
    for receipt in receipts:
        geometry = receipt.get("geometry")
        if not isinstance(geometry, Mapping):
            raise ValueError("stress receipt geometry is missing")
        z = float(geometry.get("axial_z"))
        matches = [off for off in expected_offsets if abs(z - (center_z + off)) <= 2e-14]
        if len(matches) != 1 or matches[0] in used:
            raise ValueError("shifted stress receipt z does not match unique preregistered offset")
        off = matches[0]
        used.add(off)
        indexed[off] = _validate_stress_receipt(
            receipt, expected_z=center_z + off, radii=radii
        )
    if used != set(expected_offsets):
        raise ValueError("shifted stress receipt set is incomplete")
    return indexed


def _fd4(
    indexed: Mapping[float, Mapping[str, np.ndarray]], piece: str, h: float
) -> np.ndarray:
    return (
        np.asarray(indexed[-2.0 * h][piece], dtype=float)
        - 8.0 * np.asarray(indexed[-h][piece], dtype=float)
        + 8.0 * np.asarray(indexed[h][piece], dtype=float)
        - np.asarray(indexed[2.0 * h][piece], dtype=float)
    ) / (12.0 * h)


def _production_force_levels(
    receipt: Mapping[str, Any], piece: str, radii: np.ndarray
) -> dict[float, np.ndarray]:
    block = receipt.get(piece)
    if not isinstance(block, Mapping):
        raise ValueError(f"force receipt missing {piece}")
    levels = block.get("levels")
    if not isinstance(levels, Sequence) or len(levels) != len(FROZEN_Z_STEPS):
        raise ValueError(f"force receipt has invalid {piece} levels")
    result: dict[float, np.ndarray] = {}
    for level in levels:
        if not isinstance(level, Mapping):
            raise ValueError("invalid force level")
        h = float(level.get("z_step"))
        if h not in FROZEN_Z_STEPS or h in result:
            raise ValueError("force receipt step drifted")
        arr = np.asarray(level.get("radial_force"), dtype=float)
        if arr.shape != radii.shape or not np.all(np.isfinite(arr)):
            raise ValueError("invalid production radial-force array")
        result[h] = arr
    return result


def _audit_case(
    case: Mapping[str, Any], *, expected_center: float, radii: np.ndarray
) -> dict[str, Any]:
    force = case.get("force_receipt")
    stresses = case.get("stress_receipts")
    if not isinstance(force, Mapping) or not isinstance(stresses, Sequence):
        raise ValueError("case must contain public force_receipt and stress_receipts")
    _validate_force_receipt(force, center_z=expected_center, radii=radii)
    indexed = _index_stress_receipts(stresses, center_z=expected_center, radii=radii)

    independent_levels: dict[str, dict[float, np.ndarray]] = {}
    pieces_out: dict[str, Any] = {}
    for piece in ("quadratic", "mixed", "aggregate"):
        independent_levels[piece] = {h: _fd4(indexed, piece, h) for h in FROZEN_Z_STEPS}
        production_levels = _production_force_levels(force, piece, radii)
        levels_out = [
            {
                "z_step": h,
                **_comparison_metrics(production_levels[h], independent_levels[piece][h], radii),
            }
            for h in FROZEN_Z_STEPS
        ]
        pieces_out[piece] = {
            "levels": levels_out,
            "coarse_to_medium_independent_relative_rms": _relative_rms(
                independent_levels[piece][0.02], independent_levels[piece][0.01]
            ),
            "medium_to_fine_independent_relative_rms": _relative_rms(
                independent_levels[piece][0.01], independent_levels[piece][0.005]
            ),
        }

    closure_levels = []
    for h in FROZEN_Z_STEPS:
        a = independent_levels["aggregate"][h]
        m = independent_levels["mixed"][h]
        q = independent_levels["quadratic"][h]
        scale = max(_rms(a), _rms(m), _rms(q), np.finfo(float).tiny)
        closure_levels.append(_rms(a - m - q) / scale)

    return {
        "center_z": float(expected_center),
        "pieces": pieces_out,
        "independent_piece_closure_relative_rms_levels": closure_levels,
        "production_aggregate_derivative_stability_preflight_passed": bool(
            force.get("aggregate_derivative_stability_preflight_passed")
        ),
    }


def audit_serialized_cases(cases: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    """Audit preregistered held-out z centers using public JSON receipts only."""
    if len(cases) != len(FROZEN_Z_CENTERS):
        raise ValueError("A4 requires one case for each preregistered held-out z center")
    radii = frozen_radial_grid()
    ordered: list[Mapping[str, Any]] = []
    unused = list(cases)
    for center in FROZEN_Z_CENTERS:
        matches = []
        for case in unused:
            force = case.get("force_receipt") if isinstance(case, Mapping) else None
            geometry = force.get("geometry") if isinstance(force, Mapping) else None
            if isinstance(geometry, Mapping) and abs(float(geometry.get("axial_z")) - center) <= 2e-14:
                matches.append(case)
        if len(matches) != 1:
            raise ValueError("held-out force receipts do not match unique preregistered centers")
        ordered.append(matches[0])
        unused.remove(matches[0])

    audited = [
        _audit_case(case, expected_center=center, radii=radii)
        for case, center in zip(ordered, FROZEN_Z_CENTERS, strict=True)
    ]
    failures: list[str] = []
    finest_metrics: list[Mapping[str, Any]] = []
    max_nontrivial = 0.0
    for case in audited:
        center = case["center_z"]
        if not case["production_aggregate_derivative_stability_preflight_passed"]:
            failures.append(f"z={center}:production_preflight")
        if max(case["independent_piece_closure_relative_rms_levels"]) > PIECE_CLOSURE_RELATIVE_GATE:
            failures.append(f"z={center}:independent_piece_closure")
        for piece in ("quadratic", "mixed", "aggregate"):
            block = case["pieces"][piece]
            fine = block["levels"][-1]
            finest_metrics.append(fine)
            if piece == "aggregate":
                max_nontrivial = max(max_nontrivial, float(fine["independent_rms"]))
            prefix = f"z={center}:{piece}"
            if float(fine["relative_rms"]) > FINE_FORCE_REL_RMS_GATE:
                failures.append(f"{prefix}:relative_rms")
            if float(fine["relative_max"]) > FINE_FORCE_REL_MAX_GATE:
                failures.append(f"{prefix}:relative_max")
            if float(fine["relative_integral_l2"]) > FINE_FORCE_REL_INTEGRAL_L2_GATE:
                failures.append(f"{prefix}:relative_integral_l2")
            if float(fine["axis_near_normalized_error"]) > AXIS_NEAR_NORMALIZED_ERROR_GATE:
                failures.append(f"{prefix}:axis_near")
            cm = float(block["coarse_to_medium_independent_relative_rms"])
            mf = float(block["medium_to_fine_independent_relative_rms"])
            if mf > INDEPENDENT_FINE_PAIR_REL_RMS_GATE:
                failures.append(f"{prefix}:fine_pair_stability")
            if mf > max(RESOLUTION_DEGRADATION_FACTOR * cm, RESOLUTION_NUMERICAL_FLOOR):
                failures.append(f"{prefix}:resolution_degradation")

    if max_nontrivial < NONTRIVIAL_FORCE_RMS_FLOOR:
        failures.append("aggregate_force_nontriviality")

    return {
        "schema": SCHEMA,
        "task": TASK,
        "parent_agent3_pr": PARENT_AGENT3_PR,
        "parent_agent3_head": PARENT_AGENT3_HEAD,
        "stress_agent3_pr": STRESS_AGENT3_PR,
        "stress_agent3_head": STRESS_AGENT3_HEAD,
        "seed": FROZEN_SEED,
        "time": FROZEN_TIME,
        "heldout_z_centers": list(FROZEN_Z_CENTERS),
        "z_step_ladder": list(FROZEN_Z_STEPS),
        "radial_count": int(radii.size),
        "radial_min": float(radii[0]),
        "radial_max": float(radii[-1]),
        "cases": audited,
        "fine_max_force_relative_rms": max(float(m["relative_rms"]) for m in finest_metrics),
        "fine_max_force_relative_max": max(float(m["relative_max"]) for m in finest_metrics),
        "fine_max_force_relative_integral_l2": max(
            float(m["relative_integral_l2"]) for m in finest_metrics
        ),
        "fine_max_axis_near_normalized_error": max(
            float(m["axis_near_normalized_error"]) for m in finest_metrics
        ),
        "max_aggregate_independent_force_rms": max_nontrivial,
        "failures": failures,
        "audit_pass": not failures,
        "truth_boundary": {
            "complete_ns_defect": False,
            "authorized_correction_target": False,
            "cartesian_correction_velocity_materialized": False,
            "after_correction_ns_residual_assessed": False,
            "same_protocol_comparable_to_st006": False,
            "pde_validated": False,
            "final_normalized_momentum_gate": FINAL_NORMALIZED_MOMENTUM_GATE,
            "final_normalized_divergence_gate": FINAL_NORMALIZED_DIVERGENCE_GATE,
        },
    }


def enforce_scientific_gates(report: Mapping[str, Any]) -> None:
    if report.get("schema") != SCHEMA or report.get("task") != TASK:
        raise RuntimeError("unexpected A4 report identity")
    if report.get("parent_agent3_head") != PARENT_AGENT3_HEAD:
        raise RuntimeError("A3 #1118 identity drifted")
    if report.get("audit_pass") is not True:
        raise RuntimeError(
            "independent current-I4 radial-force audit failed: "
            f"{report.get('failures')}"
        )
    boundary = report.get("truth_boundary")
    if not isinstance(boundary, Mapping) or boundary.get("pde_validated") is not False:
        raise RuntimeError("scientific truth boundary drifted")


def truth_boundary() -> dict[str, Any]:
    signature = inspect.signature(audit_serialized_cases)
    forbidden = {
        "threshold", "residual", "defect", "pressure", "forcing", "viscosity",
        "nu", "correction", "gain", "damping", "seed", "resolution", "z_step",
    }
    return {
        "schema": SCHEMA,
        "task": TASK,
        "parent_agent3_pr": PARENT_AGENT3_PR,
        "parent_agent3_head": PARENT_AGENT3_HEAD,
        "parent_agent3_source_blob": PARENT_AGENT3_SOURCE_BLOB,
        "stress_agent3_pr": STRESS_AGENT3_PR,
        "stress_agent3_head": STRESS_AGENT3_HEAD,
        "stress_agent3_source_blob": STRESS_AGENT3_SOURCE_BLOB,
        "agent2_composite_head": AGENT2_1080_HEAD,
        "agent1_leading_head": AGENT1_1079_HEAD,
        "implementation_distinct_from_parent_fd2": True,
        "independent_derivative_rule": "centered_fd4_five_point",
        "public_json_receipts_only": True,
        "caller_scientific_tuning_absent": forbidden.isdisjoint(signature.parameters),
        "heldout_z_centers": list(FROZEN_Z_CENTERS),
        "z_step_ladder": list(FROZEN_Z_STEPS),
        "current_i4_identity_consumed": True,
        "source_I4_five_row_mean_correction_materialized": False,
        "complete_ns_defect": False,
        "authorized_correction_target": False,
        "cartesian_correction_velocity_materialized": False,
        "after_correction_ns_residual_assessed": False,
        "same_protocol_comparable_to_st006": False,
        "pde_validated": False,
        "final_normalized_momentum_gate": FINAL_NORMALIZED_MOMENTUM_GATE,
        "final_normalized_divergence_gate": FINAL_NORMALIZED_DIVERGENCE_GATE,
    }
