"""Independent A4 audit of exact current-I4 nonlinear radial stress.

Consumes only JSON-round-tripped public receipts from Agent-3 PR #1109.  The
reference never calls the inherited Agent-3 #920/#875 radial inverse.  Instead
it reconstructs the compact moment-complement stress from public cylindrical
mean samples with local cubic interpolation and order-8 Gauss-Legendre cell
integration, including the formal 0->r first-cell primitive.

This is scoped correction-side operator-consistency evidence only.  It is not
a complete Navier-Stokes defect, not an authorized correction target, and not
held-out momentum-residual evidence.
"""
from __future__ import annotations

import inspect
import math
from typing import Any, Mapping, Sequence

import numpy as np

TASK = "K4-VAL-115"
SCHEMA = "kokuno-a4-current-i4-radial-stress-independent-audit-v1"

PARENT_AGENT3_PR = 1109
PARENT_AGENT3_HEAD = "49590ff311fef1dec4bd850b3013fd989485c87a"
PARENT_AGENT3_SOURCE_BLOB = "b8d9928934d4df74b52e8cd0d2c27c468257cdc9"
PARENT_SCHEMA = "kokuno-a3-current-i4-nonlinear-radial-stress-v1"
PARENT_MEAN_PR = 1101
PARENT_MEAN_HEAD = "6c9a4e7b0dec4b3669fb1dc60297b18f10163e0b"
PARENT_MEAN_SOURCE_BLOB = "10d067570ddba834db3f19c9859c2efdad16b3e8"
AGENT2_COMPOSITE_PR = 1080
AGENT2_COMPOSITE_HEAD = "c40d8ddecd2971544a6e07dab093436b423cf326"
AGENT2_COMPOSITE_SOURCE_BLOB = "2a0a5aa5966b02da856bdcf51940f3c186042802"
AGENT2_DIFFERENTIAL_PR = 960
AGENT2_DIFFERENTIAL_HEAD = "6d2fb1f701a34f783dca15a267ae2ce0734ba741"
AGENT2_DIFFERENTIAL_SOURCE_BLOB = "12df3bf6c949baeebf609b366f2973e7f515a157"
AGENT1_LEADING_PR = 1079
AGENT1_LEADING_HEAD = "b06742ca6e189499192ede3cce40f62cdc1e35ca"
AGENT1_LEADING_SOURCE_BLOB = "6f04ce0a856b44430402576dad88438da90d1ebb"
RADIAL_OPERATOR_SOURCE_BLOB = "19f8cf9fa909fef7be2cf53c84df269e696b647c"

FROZEN_SEED = 9173871
FROZEN_TIME = 0.39
FROZEN_Z = -0.08
FROZEN_R_MIN = 0.005
FROZEN_R_MAX = 0.44
FROZEN_BUMP_CENTER = 0.30
FROZEN_BUMP_HALFWIDTH = 0.10
FROZEN_GRID_COUNTS = (43, 85, 169)
FROZEN_GL_ORDER = 8
FROZEN_AXIS_NEAR_MAX_R = 0.02

FINE_STRESS_REL_RMS_GATE = 5.0e-2
FINE_STRESS_REL_MAX_GATE = 1.5e-1
FINE_MOMENT_REL_GATE = 3.0e-2
AXIS_NEAR_NORMALIZED_ERROR_GATE = 1.5e-1
OUTER_EDGE_NORMALIZED_GATE = 1.0e-8
RESOLUTION_DEGRADATION_FACTOR = 1.25
RESOLUTION_NUMERICAL_FLOOR = 2.0e-5
NONTRIVIAL_STRESS_RMS_FLOOR = 1.0e-12
MEAN_CLOSURE_ABSOLUTE_GATE = 2.0e-11
STRESS_CLOSURE_RELATIVE_GATE = 5.0e-10

FINAL_NORMALIZED_MOMENTUM_GATE = 1.0e-3
FINAL_NORMALIZED_DIVERGENCE_GATE = 1.0e-5

_GL_X, _GL_W = np.polynomial.legendre.leggauss(FROZEN_GL_ORDER)


def frozen_radial_grids() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return frozen coarse/medium/fine grids and one independent off-grid ladder."""
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
    radii: np.ndarray, values: np.ndarray, *, exponent: int
) -> np.ndarray:
    """Integrate r**e * values from the formal axis using local cubic + GL8."""
    r = np.asarray(radii, dtype=float)
    y = np.asarray(values, dtype=float)
    if r.ndim != 1 or y.shape != r.shape or r.size < 9:
        raise ValueError("independent radial integral requires matching 1D arrays with >=9 nodes")
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
    radii: np.ndarray, *, exponent: int, center: float, halfwidth: float
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
        radii, exponent=exponent, center=bump_center, halfwidth=bump_halfwidth
    )
    moment = float(source_cumulative[-1])
    bump_integral = float(bump_cumulative[-1])
    if not math.isfinite(bump_integral) or bump_integral <= 0.0:
        raise ValueError("independent compact bump lost positive normalization")
    complement = source_cumulative - (moment / bump_integral) * bump_cumulative
    r = np.asarray(radii, dtype=float)
    stress = -complement / (r**exponent)
    return {
        "weighted_moment": moment,
        "bump_weighted_integral": bump_integral,
        "stress": stress,
        "stress_rms": _rms(stress),
        "stress_outer_edge": float(stress[-1]),
    }


def _extract_sources(receipt: Mapping[str, Any]) -> dict[str, np.ndarray]:
    mean = receipt.get("current_i4_mean_witness")
    if not isinstance(mean, Mapping):
        raise ValueError("receipt is missing current_i4_mean_witness")
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
    max_scale = max(float(np.max(np.abs(p))), float(np.max(np.abs(q))), np.finfo(float).tiny)
    return _rms(delta) / rms_scale, float(np.max(np.abs(delta))) / max_scale


def _validate_receipt(
    receipt: Mapping[str, Any], *, expected_radii: np.ndarray
) -> tuple[np.ndarray, dict[str, np.ndarray], float, float]:
    if receipt.get("schema") != PARENT_SCHEMA:
        raise ValueError("unexpected A3 #1109 receipt schema")
    if receipt.get("parent_agent3_1101_source_blob") != PARENT_MEAN_SOURCE_BLOB:
        raise ValueError("A3 #1101 source provenance drifted")

    provenance = receipt.get("provenance")
    if not isinstance(provenance, Mapping):
        raise ValueError("receipt provenance is missing")
    expected = {
        "parent_agent3_head": PARENT_MEAN_HEAD,
        "radial_operator_source_blob": RADIAL_OPERATOR_SOURCE_BLOB,
        "agent2_composite_head": AGENT2_COMPOSITE_HEAD,
        "agent2_composite_source_blob": AGENT2_COMPOSITE_SOURCE_BLOB,
        "agent2_differential_head": AGENT2_DIFFERENTIAL_HEAD,
        "agent2_differential_source_blob": AGENT2_DIFFERENTIAL_SOURCE_BLOB,
        "agent1_leading_head": AGENT1_LEADING_HEAD,
        "agent1_leading_source_blob": AGENT1_LEADING_SOURCE_BLOB,
    }
    for key, value in expected.items():
        if provenance.get(key) != value:
            raise ValueError(f"A3 #1109 provenance drifted at {key}")

    geometry = receipt.get("geometry")
    if not isinstance(geometry, Mapping):
        raise ValueError("receipt geometry is missing")
    radii = np.asarray(geometry.get("radii"), dtype=float)
    if radii.shape != expected_radii.shape or not np.allclose(
        radii, expected_radii, rtol=0.0, atol=2e-15
    ):
        raise ValueError("receipt radial grid does not match preregistered A4 grid")
    if abs(float(geometry.get("time")) - FROZEN_TIME) > 1e-15:
        raise ValueError("receipt time drifted")
    if abs(float(geometry.get("axial_z")) - FROZEN_Z) > 1e-15:
        raise ValueError("receipt axial z drifted")
    center = float(geometry.get("bump_center"))
    halfwidth = float(geometry.get("bump_halfwidth"))
    if abs(center - FROZEN_BUMP_CENTER) > 1e-15 or abs(
        halfwidth - FROZEN_BUMP_HALFWIDTH
    ) > 1e-15:
        raise ValueError("receipt compact-bump geometry drifted")

    boundary = receipt.get("truth_boundary")
    if not isinstance(boundary, Mapping):
        raise ValueError("receipt truth boundary is missing")
    for key in (
        "source_I4_five_row_mean_correction_materialized",
        "current_I4_correction_velocity_materialized",
        "radial_force_materialized_in_this_increment",
        "complete_ns_defect",
        "current_real_ns_correction_velocity_materialized",
        "real_candidate_finite_correction_cycle_run",
        "heldout_normalized_ns_residual_assessed",
        "residual_reduction_claimed",
        "same_protocol_comparable_to_st006",
        "pde_validated",
    ):
        if boundary.get(key) is not False:
            raise ValueError(f"A3 #1109 scientific boundary drifted at {key}")
    if boundary.get("final_normalized_momentum_gate") != FINAL_NORMALIZED_MOMENTUM_GATE:
        raise ValueError("final momentum gate drifted")
    if boundary.get("final_normalized_divergence_gate") != FINAL_NORMALIZED_DIVERGENCE_GATE:
        raise ValueError("final divergence gate drifted")

    sources = _extract_sources(receipt)
    if any(arr.shape[0] != radii.size for arr in sources.values()):
        raise ValueError("mean sample count does not match radial grid")
    closure = float(np.max(np.abs(sources["aggregate"] - sources["mixed"] - sources["quadratic"])))
    if closure > MEAN_CLOSURE_ABSOLUTE_GATE:
        raise ValueError("public mean attribution no longer closes")
    return radii, sources, center, halfwidth


def _audit_one_receipt(receipt: Mapping[str, Any], *, expected_radii: np.ndarray) -> dict[str, Any]:
    radii, sources, center, halfwidth = _validate_receipt(
        receipt, expected_radii=expected_radii
    )
    output: dict[str, Any] = {"count": int(radii.size), "channels": {}}
    for channel, exponent, component in (("theta_e2", 2, 1), ("axial_e1", 1, 2)):
        parent_channel = receipt.get(channel)
        if not isinstance(parent_channel, Mapping):
            raise ValueError(f"receipt is missing {channel}")
        channel_out: dict[str, Any] = {}
        prod_arrays: dict[str, np.ndarray] = {}
        ref_arrays: dict[str, np.ndarray] = {}
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
            reference = np.asarray(independent["stress"], dtype=float)
            rel_rms, rel_max = _relative_metrics(production, reference)
            prod_moment = float(parent_piece.get("weighted_moment"))
            ref_moment = float(independent["weighted_moment"])
            moment_scale = max(abs(prod_moment), abs(ref_moment), np.finfo(float).tiny)
            moment_rel = abs(prod_moment - ref_moment) / moment_scale
            axis = radii <= FROZEN_AXIS_NEAR_MAX_R
            if not np.any(axis):
                raise ValueError("preregistered grid lost axis-near samples")
            global_scale = max(
                float(np.max(np.abs(production))),
                float(np.max(np.abs(reference))),
                np.finfo(float).tiny,
            )
            axis_error = float(np.max(np.abs(production[axis] - reference[axis]))) / global_scale
            outer = max(abs(float(production[-1])), abs(float(reference[-1]))) / max(
                _rms(production), _rms(reference), np.finfo(float).tiny
            )
            delta = np.abs(production - reference)
            worst = int(np.argmax(delta))
            channel_out[piece] = {
                "production_stress_rms": _rms(production),
                "independent_stress_rms": _rms(reference),
                "relative_rms": float(rel_rms),
                "relative_max": float(rel_max),
                "weighted_moment_relative_error": float(moment_rel),
                "axis_near_normalized_error": float(axis_error),
                "outer_edge_normalized_stress": float(outer),
                "worst_r": float(radii[worst]),
                "worst_absolute_delta": float(delta[worst]),
                "worst_production": float(production[worst]),
                "worst_independent": float(reference[worst]),
            }
            prod_arrays[piece] = production
            ref_arrays[piece] = reference

        prod_scale = max(
            _rms(prod_arrays["aggregate"]),
            _rms(prod_arrays["mixed"]),
            _rms(prod_arrays["quadratic"]),
            np.finfo(float).tiny,
        )
        ref_scale = max(
            _rms(ref_arrays["aggregate"]),
            _rms(ref_arrays["mixed"]),
            _rms(ref_arrays["quadratic"]),
            np.finfo(float).tiny,
        )
        channel_out["production_piece_closure_relative_max"] = float(
            np.max(np.abs(prod_arrays["aggregate"] - prod_arrays["mixed"] - prod_arrays["quadratic"]))
            / prod_scale
        )
        channel_out["independent_piece_closure_relative_max"] = float(
            np.max(np.abs(ref_arrays["aggregate"] - ref_arrays["mixed"] - ref_arrays["quadratic"]))
            / ref_scale
        )
        output["channels"][channel] = channel_out
    return output


def _gate_level(level: Mapping[str, Any]) -> list[str]:
    failures: list[str] = []
    for channel in ("theta_e2", "axial_e1"):
        c = level["channels"][channel]
        if c["production_piece_closure_relative_max"] > STRESS_CLOSURE_RELATIVE_GATE:
            failures.append(f"{channel}:production_piece_closure")
        if c["independent_piece_closure_relative_max"] > STRESS_CLOSURE_RELATIVE_GATE:
            failures.append(f"{channel}:independent_piece_closure")
        for piece in ("quadratic", "mixed", "aggregate"):
            p = c[piece]
            if p["relative_rms"] > FINE_STRESS_REL_RMS_GATE:
                failures.append(f"{channel}/{piece}:relative_rms")
            if p["relative_max"] > FINE_STRESS_REL_MAX_GATE:
                failures.append(f"{channel}/{piece}:relative_max")
            if p["weighted_moment_relative_error"] > FINE_MOMENT_REL_GATE:
                failures.append(f"{channel}/{piece}:weighted_moment")
            if p["axis_near_normalized_error"] > AXIS_NEAR_NORMALIZED_ERROR_GATE:
                failures.append(f"{channel}/{piece}:axis_near")
            if p["outer_edge_normalized_stress"] > OUTER_EDGE_NORMALIZED_GATE:
                failures.append(f"{channel}/{piece}:outer_edge")
        if c["aggregate"]["independent_stress_rms"] < NONTRIVIAL_STRESS_RMS_FLOOR:
            failures.append(f"{channel}:aggregate_nontriviality")
    return failures


def audit_serialized_receipts(receipts: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    """Audit four JSON-round-tripped A3 receipts under one frozen protocol."""
    grids = frozen_radial_grids()
    if len(receipts) != len(grids):
        raise ValueError("A4 audit requires coarse/medium/fine/off-grid receipts")
    levels = [
        _audit_one_receipt(receipt, expected_radii=grid)
        for receipt, grid in zip(receipts, grids)
    ]

    failures: list[str] = []
    for label, index in (("fine", 2), ("offgrid", 3)):
        failures.extend(f"{label}:{item}" for item in _gate_level(levels[index]))

    trends: dict[str, Any] = {}
    for channel in ("theta_e2", "axial_e1"):
        trends[channel] = {}
        for piece in ("quadratic", "mixed", "aggregate"):
            seq = [float(levels[i]["channels"][channel][piece]["relative_rms"]) for i in range(3)]
            checks = []
            for previous, current in zip(seq[:-1], seq[1:]):
                allowed = max(
                    RESOLUTION_NUMERICAL_FLOOR,
                    RESOLUTION_DEGRADATION_FACTOR * previous,
                )
                checks.append(bool(current <= allowed))
            trends[channel][piece] = {
                "coarse_medium_fine_relative_rms": seq,
                "adjacent_nondegradation": checks,
            }
            if not all(checks):
                failures.append(f"resolution:{channel}/{piece}")

    return {
        "schema": SCHEMA,
        "task": TASK,
        "frozen_protocol": {
            "seed": FROZEN_SEED,
            "time": FROZEN_TIME,
            "z": FROZEN_Z,
            "r_min": FROZEN_R_MIN,
            "r_max": FROZEN_R_MAX,
            "grid_counts": list(FROZEN_GRID_COUNTS),
            "offgrid_count": FROZEN_GRID_COUNTS[-1],
            "gl_order": FROZEN_GL_ORDER,
            "axis_near_max_r": FROZEN_AXIS_NEAR_MAX_R,
            "bump_center": FROZEN_BUMP_CENTER,
            "bump_halfwidth": FROZEN_BUMP_HALFWIDTH,
        },
        "levels": levels,
        "resolution_trends": trends,
        "failures": failures,
        "audit_pass": not failures,
        "truth_boundary": truth_boundary(),
    }


def enforce_scientific_gates(report: Mapping[str, Any]) -> None:
    if report.get("schema") != SCHEMA or report.get("task") != TASK:
        raise RuntimeError("unexpected A4 audit report identity")
    if report.get("audit_pass") is not True:
        failures = report.get("failures", ())
        raise RuntimeError(f"preregistered current-I4 radial-stress audit failed: {failures}")


def truth_boundary() -> dict[str, Any]:
    signature = inspect.signature(audit_serialized_receipts)
    forbidden = {
        "threshold", "scientific_threshold", "residual", "defect", "forcing",
        "pressure", "viscosity", "nu", "derivative_step", "gain", "damping",
    }
    return {
        "parent_agent3_pr": PARENT_AGENT3_PR,
        "parent_agent3_head": PARENT_AGENT3_HEAD,
        "parent_agent3_source_blob": PARENT_AGENT3_SOURCE_BLOB,
        "parent_mean_pr": PARENT_MEAN_PR,
        "parent_mean_head": PARENT_MEAN_HEAD,
        "agent2_composite_pr": AGENT2_COMPOSITE_PR,
        "agent2_composite_head": AGENT2_COMPOSITE_HEAD,
        "agent1_leading_pr": AGENT1_LEADING_PR,
        "agent1_leading_head": AGENT1_LEADING_HEAD,
        "independent_radial_operator": "local_cubic_plus_GL8_formal_axis_primitive",
        "production_radial_operator_reused_by_reference": False,
        "serialized_public_receipts_only": True,
        "caller_supplied_thresholds_absent": forbidden.isdisjoint(signature.parameters),
        "source_I4_five_row_mean_correction_materialized": False,
        "radial_force_partial_z_sigma1_materialized": False,
        "cartesian_correction_velocity_materialized": False,
        "matched_pressure_materialized": False,
        "restricted_forcing_materialized": False,
        "complete_ns_defect": False,
        "authorized_correction_target": False,
        "after_correction_ns_residual_assessed": False,
        "same_protocol_comparable_to_st006": False,
        "pde_validated": False,
        "final_normalized_momentum_gate": FINAL_NORMALIZED_MOMENTUM_GATE,
        "final_normalized_divergence_gate": FINAL_NORMALIZED_DIVERGENCE_GATE,
    }
