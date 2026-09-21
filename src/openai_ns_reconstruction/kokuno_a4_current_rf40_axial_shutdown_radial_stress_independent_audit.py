"""Independent A4 audit of the current RF40 axial-shutdown radial stress.

Consumes serialized/public receipts from A3 PR #1015 only.  The numerical
reference does not call the inherited A3 #920/#875 radial inverse.  Instead it
reconstructs the compact moment-complement stresses from public cylindrical
mean samples with local cubic interpolation and order-8 Gauss--Legendre cell
integration.

This is scoped correction-side operator-consistency evidence.  It is not a
complete Navier--Stokes defect, not an authorized correction target, not a
Cartesian correction velocity, and not held-out momentum-residual evidence.
"""
from __future__ import annotations

import math
from typing import Any, Mapping, Sequence

import numpy as np

TASK = "K4-VAL-103"
SCHEMA = "kokuno-a4-current-rf40-axial-shutdown-radial-stress-independent-audit-v1"

PARENT_AGENT3_PR = 1015
PARENT_AGENT3_HEAD = "2e7683f06bea810881afbe39d8f0cf64815ddd56"
PARENT_AGENT3_SOURCE_BLOB = "9c6082e6dc5edacf7b6a64e87549cbec142b3d4a"
PARENT_SCHEMA = "kokuno-a3-current-rf40-axial-shutdown-nonlinear-radial-stress-v1"
PARENT_AGENT3_1011_SOURCE_BLOB = "1ccef44607487d7c3cbdfa6dbce0067fdd302ea5"
AGENT2_999_HEAD = "93b99292fcf141f24b7c6d7e4fbf95775e6a07e5"
AGENT2_999_SOURCE_BLOB = "255e65ce0c37652858509c33e7c9fad40b73ca97"
AGENT1_993_HEAD = "2ac6460b483efb1c07f2fa65e7fed781a32f2718"
AGENT1_993_SOURCE_BLOB = "057a0514c8a941c3e59922b8158f480434b4441e"
RADIAL_OPERATOR_SOURCE_BLOB = "19f8cf9fa909fef7be2cf53c84df269e696b647c"

FROZEN_SEED = 9173751
FROZEN_TIME = 0.39
FROZEN_Z = -0.08
FROZEN_R_MIN = 0.005
FROZEN_R_MAX = 0.44
FROZEN_BUMP_CENTER = 0.30
FROZEN_BUMP_HALFWIDTH = 0.10
FROZEN_GRID_COUNTS = (43, 85, 169)
FROZEN_GL_ORDER = 8
FROZEN_AXIS_NEAR_MAX_R = 0.02

# Preserved unchanged from A4 PR #977.  These are scoped stress gates, not PDE gates.
FINE_STRESS_REL_RMS_GATE = 5.0e-2
FINE_STRESS_REL_MAX_GATE = 1.5e-1
FINE_MOMENT_REL_GATE = 3.0e-2
AXIS_NEAR_NORMALIZED_ERROR_GATE = 1.5e-1
OUTER_EDGE_NORMALIZED_GATE = 1.0e-8
RESOLUTION_DEGRADATION_FACTOR = 1.25
RESOLUTION_NUMERICAL_FLOOR = 2.0e-5
NONTRIVIAL_STRESS_RMS_FLOOR = 1.0e-12

FINAL_NORMALIZED_MOMENTUM_GATE = 1.0e-3
FINAL_NORMALIZED_DIVERGENCE_GATE = 1.0e-5

_GL_X, _GL_W = np.polynomial.legendre.leggauss(FROZEN_GL_ORDER)


def frozen_radial_grids() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return frozen coarse/medium/fine grids and one shifted fine off-grid grid."""
    regular = tuple(
        np.linspace(FROZEN_R_MIN, FROZEN_R_MAX, n, dtype=float)
        for n in FROZEN_GRID_COUNTS
    )
    fine = np.array(regular[-1], copy=True)
    h = (FROZEN_R_MAX - FROZEN_R_MIN) / (fine.size - 1)
    rng = np.random.default_rng(FROZEN_SEED)
    shifted = np.array(fine, copy=True)
    shifted[1:-1] += rng.uniform(-0.24, 0.24, size=fine.size - 2) * h
    if np.any(np.diff(shifted) <= 0.0):
        raise RuntimeError("frozen off-grid jitter lost monotonicity")
    return (*regular, shifted)


def _rms(values: np.ndarray) -> float:
    a = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(a * a)))


def _compact_cos8_bump(x: np.ndarray, center: float, halfwidth: float) -> np.ndarray:
    s = (np.asarray(x, dtype=float) - float(center)) / float(halfwidth)
    out = np.zeros_like(s)
    mask = np.abs(s) < 1.0
    out[mask] = np.cos(0.5 * math.pi * s[mask]) ** 8
    return out


def _piecewise_cubic_weighted_cumulative(
    radii: np.ndarray,
    values: np.ndarray,
    *,
    exponent: int,
) -> np.ndarray:
    """Integrate r**exponent * values by local cubic interpolation + GL8 cells."""
    r = np.asarray(radii, dtype=float)
    y = np.asarray(values, dtype=float)
    if r.ndim != 1 or y.shape != r.shape or r.size < 9:
        raise ValueError("independent radial integral needs matching 1D arrays with >=9 nodes")
    if exponent not in (1, 2):
        raise ValueError("only e=1/e=2 are supported")
    if np.any(r <= 0.0) or np.any(np.diff(r) <= 0.0):
        raise ValueError("radii must be positive and strictly increasing")
    if not (np.all(np.isfinite(r)) and np.all(np.isfinite(y))):
        raise ValueError("radii and values must be finite")

    edges = np.concatenate(([0.0], r))
    cumulative = np.empty_like(r)
    total = 0.0
    for cell in range(r.size):
        a = float(edges[cell])
        b = float(edges[cell + 1])
        start = 0 if cell == 0 else max(0, min(cell - 2, r.size - 4))
        nodes = r[start : start + 4]
        vals = y[start : start + 4]
        mid = 0.5 * (a + b)
        half = 0.5 * (b - a)
        coeff = np.polynomial.polynomial.polyfit(nodes - mid, vals, deg=3)
        sample = mid + half * _GL_X
        interp = np.polynomial.polynomial.polyval(sample - mid, coeff)
        total += half * float(np.sum(_GL_W * (sample**exponent) * interp))
        cumulative[cell] = total
    return cumulative


def _direct_bump_weighted_cumulative(
    radii: np.ndarray,
    *,
    exponent: int,
    center: float,
    halfwidth: float,
) -> np.ndarray:
    r = np.asarray(radii, dtype=float)
    edges = np.concatenate(([0.0], r))
    cumulative = np.empty_like(r)
    total = 0.0
    for cell in range(r.size):
        a = float(edges[cell])
        b = float(edges[cell + 1])
        mid = 0.5 * (a + b)
        half = 0.5 * (b - a)
        sample = mid + half * _GL_X
        bump = _compact_cos8_bump(sample, center, halfwidth)
        total += half * float(np.sum(_GL_W * (sample**exponent) * bump))
        cumulative[cell] = total
    return cumulative


def _independent_compact_stress(
    radii: np.ndarray,
    source: np.ndarray,
    *,
    exponent: int,
    bump_center: float,
    bump_halfwidth: float,
) -> dict[str, Any]:
    source_cumulative = _piecewise_cubic_weighted_cumulative(
        radii, source, exponent=exponent
    )
    bump_cumulative = _direct_bump_weighted_cumulative(
        radii,
        exponent=exponent,
        center=bump_center,
        halfwidth=bump_halfwidth,
    )
    moment = float(source_cumulative[-1])
    bump_integral = float(bump_cumulative[-1])
    if not math.isfinite(bump_integral) or bump_integral <= 0.0:
        raise ValueError("independent compact bump lost positive normalization")
    complement_cumulative = source_cumulative - (moment / bump_integral) * bump_cumulative
    r = np.asarray(radii, dtype=float)
    stress = -complement_cumulative / (r**exponent)
    return {
        "weighted_moment": moment,
        "bump_weighted_integral": bump_integral,
        "stress": stress,
        "stress_rms": _rms(stress),
        "stress_outer_edge": float(stress[-1]),
    }


def _extract_sources(receipt: Mapping[str, Any]) -> dict[str, np.ndarray]:
    mean = receipt.get("current_rf40_axial_shutdown_mean_witness")
    if not isinstance(mean, Mapping):
        raise ValueError("receipt is missing current_rf40_axial_shutdown_mean_witness")
    names = {
        "quadratic": "mean_oscillatory_self_advection_cylindrical",
        "mixed": "mean_mixed_cross_cylindrical",
        "aggregate": "mean_aggregate_nonlinear_cylindrical",
    }
    out: dict[str, np.ndarray] = {}
    for piece, key in names.items():
        arr = np.asarray(mean.get(key), dtype=float)
        if arr.ndim != 2 or arr.shape[1] != 3 or not np.all(np.isfinite(arr)):
            raise ValueError(f"invalid public {piece} cylindrical mean array")
        out[piece] = arr
    return out


def _relative_metrics(production: np.ndarray, independent: np.ndarray) -> tuple[float, float]:
    p = np.asarray(production, dtype=float)
    q = np.asarray(independent, dtype=float)
    if p.shape != q.shape or not (np.all(np.isfinite(p)) and np.all(np.isfinite(q))):
        raise ValueError("stress arrays must have matching finite shapes")
    delta = p - q
    rms_scale = max(_rms(p), _rms(q), np.finfo(float).tiny)
    max_scale = max(
        float(np.max(np.abs(p))),
        float(np.max(np.abs(q))),
        np.finfo(float).tiny,
    )
    return _rms(delta) / rms_scale, float(np.max(np.abs(delta))) / max_scale


def _audit_one_receipt(
    receipt: Mapping[str, Any],
    *,
    expected_radii: np.ndarray,
) -> dict[str, Any]:
    if receipt.get("schema") != PARENT_SCHEMA:
        raise ValueError("unexpected A3 #1015 receipt schema")
    if receipt.get("parent_agent3_1011_source_blob") != PARENT_AGENT3_1011_SOURCE_BLOB:
        raise ValueError("A3 #1011 provenance drifted")

    provenance = receipt.get("provenance")
    if not isinstance(provenance, Mapping):
        raise ValueError("receipt provenance is missing")
    expected_provenance = {
        "parent_agent3_head": "6aabf7b6dd28ae683a3776d518e9a977507911e9",
        "parent_agent3_source_blob": PARENT_AGENT3_1011_SOURCE_BLOB,
        "radial_operator_source_blob": RADIAL_OPERATOR_SOURCE_BLOB,
        "agent2_composite_head": AGENT2_999_HEAD,
        "agent2_composite_source_blob": AGENT2_999_SOURCE_BLOB,
        "agent1_leading_head": AGENT1_993_HEAD,
        "agent1_leading_source_blob": AGENT1_993_SOURCE_BLOB,
    }
    for key, expected in expected_provenance.items():
        if provenance.get(key) != expected:
            raise ValueError(f"A3 #1015 provenance drifted at {key}")

    geometry = receipt.get("geometry")
    if not isinstance(geometry, Mapping):
        raise ValueError("receipt geometry is missing")
    radii = np.asarray(geometry.get("radii"), dtype=float)
    if (
        radii.shape != expected_radii.shape
        or not np.allclose(radii, expected_radii, rtol=0.0, atol=2e-15)
    ):
        raise ValueError("receipt radial grid does not match preregistered A4 grid")
    if abs(float(geometry.get("time")) - FROZEN_TIME) > 1e-15:
        raise ValueError("receipt time drifted")
    if abs(float(geometry.get("axial_z")) - FROZEN_Z) > 1e-15:
        raise ValueError("receipt axial z drifted")
    center = float(geometry.get("bump_center"))
    halfwidth = float(geometry.get("bump_halfwidth"))
    if (
        abs(center - FROZEN_BUMP_CENTER) > 1e-15
        or abs(halfwidth - FROZEN_BUMP_HALFWIDTH) > 1e-15
    ):
        raise ValueError("receipt compact-bump geometry drifted")

    boundary = receipt.get("truth_boundary")
    if not isinstance(boundary, Mapping):
        raise ValueError("receipt truth boundary is missing")
    false_keys = (
        "complete_ns_defect",
        "scoped_current_RF40_axial_shutdown_radial_stress_authorized_as_correction_target",
        "current_real_ns_correction_velocity_materialized",
        "real_candidate_finite_correction_cycle_run",
        "heldout_normalized_ns_residual_assessed",
        "residual_reduction_claimed",
        "same_protocol_comparable_to_st006",
        "pde_validated",
    )
    for key in false_keys:
        if boundary.get(key) is not False:
            raise ValueError(f"A3 #1015 scientific boundary drifted at {key}")
    if boundary.get("final_normalized_momentum_gate") != FINAL_NORMALIZED_MOMENTUM_GATE:
        raise ValueError("final momentum gate drifted")
    if boundary.get("final_normalized_divergence_gate") != FINAL_NORMALIZED_DIVERGENCE_GATE:
        raise ValueError("final divergence gate drifted")

    sources = _extract_sources(receipt)
    n = radii.size
    if any(arr.shape[0] != n for arr in sources.values()):
        raise ValueError("mean sample count does not match radial grid")
    mean_closure = float(
        np.max(np.abs(sources["aggregate"] - sources["mixed"] - sources["quadratic"]))
    )
    if mean_closure > 2.0e-11:
        raise ValueError("public mean attribution no longer closes")

    output: dict[str, Any] = {
        "count": int(n),
        "mean_piece_closure_absolute_max": mean_closure,
        "channels": {},
    }
    for channel, exponent, component in (("theta_e2", 2, 1), ("axial_e1", 1, 2)):
        parent_channel = receipt.get(channel)
        if not isinstance(parent_channel, Mapping):
            raise ValueError(f"receipt is missing {channel}")
        channel_out: dict[str, Any] = {}
        for piece in ("quadratic", "mixed", "aggregate"):
            parent_piece = parent_channel.get(piece)
            if not isinstance(parent_piece, Mapping):
                raise ValueError(f"receipt is missing {channel}/{piece}")
            production = np.asarray(parent_piece.get("stress"), dtype=float)
            if production.shape != radii.shape or not np.all(np.isfinite(production)):
                raise ValueError(f"invalid production stress for {channel}/{piece}")
            independent = _independent_compact_stress(
                radii,
                sources[piece][:, component],
                exponent=exponent,
                bump_center=center,
                bump_halfwidth=halfwidth,
            )
            rel_rms, rel_max = _relative_metrics(production, independent["stress"])
            prod_moment = float(parent_piece.get("weighted_moment"))
            ind_moment = float(independent["weighted_moment"])
            moment_rel = abs(prod_moment - ind_moment) / max(
                abs(prod_moment), abs(ind_moment), np.finfo(float).tiny
            )
            stress_scale = max(
                float(np.max(np.abs(production))),
                float(np.max(np.abs(independent["stress"]))),
                np.finfo(float).tiny,
            )
            axis_mask = radii <= FROZEN_AXIS_NEAR_MAX_R
            if not np.any(axis_mask):
                raise ValueError("preregistered grid lost positive-radius axis-near samples")
            axis_error = float(
                np.max(np.abs(production[axis_mask] - independent["stress"][axis_mask]))
                / stress_scale
            )
            outer_norm = abs(float(production[-1])) / stress_scale
            channel_out[piece] = {
                "stress_relative_rms": float(rel_rms),
                "stress_relative_max": float(rel_max),
                "weighted_moment_relative_error": float(moment_rel),
                "axis_near_normalized_error": float(axis_error),
                "outer_edge_normalized_stress": float(outer_norm),
                "production_stress_rms": _rms(production),
                "independent_stress_rms": float(independent["stress_rms"]),
                "independent_weighted_moment": ind_moment,
            }
        output["channels"][channel] = channel_out
    return output


def _gate_level(level: Mapping[str, Any]) -> list[str]:
    failures: list[str] = []
    channels = level["channels"]
    for channel in ("theta_e2", "axial_e1"):
        for piece in ("quadratic", "mixed", "aggregate"):
            m = channels[channel][piece]
            prefix = f"{channel}/{piece}"
            if float(m["stress_relative_rms"]) > FINE_STRESS_REL_RMS_GATE:
                failures.append(f"{prefix}:stress_relative_rms")
            if float(m["stress_relative_max"]) > FINE_STRESS_REL_MAX_GATE:
                failures.append(f"{prefix}:stress_relative_max")
            if float(m["weighted_moment_relative_error"]) > FINE_MOMENT_REL_GATE:
                failures.append(f"{prefix}:weighted_moment_relative_error")
            if float(m["axis_near_normalized_error"]) > AXIS_NEAR_NORMALIZED_ERROR_GATE:
                failures.append(f"{prefix}:axis_near_normalized_error")
            if float(m["outer_edge_normalized_stress"]) > OUTER_EDGE_NORMALIZED_GATE:
                failures.append(f"{prefix}:outer_edge_normalized_stress")
    aggregate_rms = max(
        float(channels[channel]["aggregate"]["production_stress_rms"])
        for channel in ("theta_e2", "axial_e1")
    )
    if aggregate_rms < NONTRIVIAL_STRESS_RMS_FLOOR:
        failures.append("aggregate_stress_nontriviality")
    return failures


def audit_serialized_receipts(receipts: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    """Audit three nested-resolution receipts plus one shifted fine-grid receipt."""
    grids = frozen_radial_grids()
    if len(receipts) != len(grids):
        raise ValueError("A4 audit requires coarse/medium/fine/off-grid receipts")
    levels = [
        _audit_one_receipt(receipt, expected_radii=grid)
        for receipt, grid in zip(receipts, grids, strict=True)
    ]

    failures: list[str] = []
    # Scientific acceptance is enforced on fine and independently shifted off-grid data.
    for label, index in (("fine", 2), ("offgrid", 3)):
        failures.extend(f"{label}:{failure}" for failure in _gate_level(levels[index]))

    # Resolution trend must not degrade outside the frozen numerical floor.
    for channel in ("theta_e2", "axial_e1"):
        for piece in ("quadratic", "mixed", "aggregate"):
            vals = [
                float(levels[i]["channels"][channel][piece]["stress_relative_rms"])
                for i in range(3)
            ]
            if vals[1] > max(
                RESOLUTION_DEGRADATION_FACTOR * vals[0], RESOLUTION_NUMERICAL_FLOOR
            ):
                failures.append(f"resolution:{channel}/{piece}:coarse_to_medium")
            if vals[2] > max(
                RESOLUTION_DEGRADATION_FACTOR * vals[1], RESOLUTION_NUMERICAL_FLOOR
            ):
                failures.append(f"resolution:{channel}/{piece}:medium_to_fine")

    fine_metrics = [
        levels[2]["channels"][channel][piece]
        for channel in ("theta_e2", "axial_e1")
        for piece in ("quadratic", "mixed", "aggregate")
    ]
    offgrid_metrics = [
        levels[3]["channels"][channel][piece]
        for channel in ("theta_e2", "axial_e1")
        for piece in ("quadratic", "mixed", "aggregate")
    ]
    return {
        "schema": SCHEMA,
        "task": TASK,
        "parent_agent3_pr": PARENT_AGENT3_PR,
        "parent_agent3_head": PARENT_AGENT3_HEAD,
        "parent_agent3_source_blob": PARENT_AGENT3_SOURCE_BLOB,
        "seed": FROZEN_SEED,
        "grid_counts": list(FROZEN_GRID_COUNTS),
        "offgrid_count": int(grids[-1].size),
        "gauss_legendre_order": FROZEN_GL_ORDER,
        "levels": levels,
        "fine_max_stress_relative_rms": max(float(m["stress_relative_rms"]) for m in fine_metrics),
        "fine_max_stress_relative_max": max(float(m["stress_relative_max"]) for m in fine_metrics),
        "fine_max_weighted_moment_relative_error": max(
            float(m["weighted_moment_relative_error"]) for m in fine_metrics
        ),
        "fine_max_axis_near_normalized_error": max(
            float(m["axis_near_normalized_error"]) for m in fine_metrics
        ),
        "fine_max_outer_edge_normalized_stress": max(
            float(m["outer_edge_normalized_stress"]) for m in fine_metrics
        ),
        "offgrid_max_stress_relative_rms": max(
            float(m["stress_relative_rms"]) for m in offgrid_metrics
        ),
        "offgrid_max_stress_relative_max": max(
            float(m["stress_relative_max"]) for m in offgrid_metrics
        ),
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
    """Fail closed on any preregistered scoped stress-consistency gate."""
    if report.get("schema") != SCHEMA or report.get("task") != TASK:
        raise RuntimeError("unexpected A4 report identity")
    if report.get("parent_agent3_head") != PARENT_AGENT3_HEAD:
        raise RuntimeError("A3 #1015 identity drifted")
    if report.get("audit_pass") is not True:
        failures = report.get("failures")
        raise RuntimeError(f"independent RF40 axial-shutdown stress audit failed: {failures}")
    boundary = report.get("truth_boundary")
    if not isinstance(boundary, Mapping) or boundary.get("pde_validated") is not False:
        raise RuntimeError("scientific truth boundary drifted")
