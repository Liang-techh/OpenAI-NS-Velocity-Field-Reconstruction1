"""Independent Agent-4 audit of the current partial nonlinear radial stress.

This module consumes only the serialized/public receipt emitted by Agent-3 PR
#976.  It does not call the inherited #920/#875 radial inverse when forming its
numerical reference.  Instead, each public cylindrical mean profile is
re-integrated with a distinct piecewise-cubic / Gauss--Legendre path and the
compact moment-complement stress is reconstructed from those public samples.

The audit is deliberately scoped.  It checks consistency of the current
partial-domain nonlinear mean -> radial-stress step only.  It is not a complete
Navier--Stokes defect, a correction target, a correction velocity, or residual
evidence.
"""
from __future__ import annotations

import copy
import math
from typing import Any, Mapping, Sequence

import numpy as np

TASK = "K4-VAL-096"
SCHEMA = "kokuno-a4-current-partial-radial-stress-independent-audit-v1"

PARENT_AGENT3_PR = 976
PARENT_AGENT3_HEAD = "0a777668050ff9e73d6c8db1f72e5fabfa9b04a6"
PARENT_AGENT3_SOURCE_BLOB = "9afa50c573186b86cd86fa7ffe44a5986a808b6d"
PARENT_SCHEMA = "kokuno-a3-current-partial-nonlinear-radial-stress-v1"

FROZEN_SEED = 9173681
FROZEN_TIME = 0.39
FROZEN_Z = -0.08
FROZEN_R_MIN = 0.005
FROZEN_R_MAX = 0.44
FROZEN_BUMP_CENTER = 0.30
FROZEN_BUMP_HALFWIDTH = 0.10
FROZEN_GRID_COUNTS = (43, 85, 169)
FROZEN_GL_ORDER = 8
FROZEN_AXIS_NEAR_MAX_R = 0.02

# Scoped stress-consistency gates.  These are not CR001 PDE gates.
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
    """Return the preregistered coarse/medium/fine plus shifted off-grid grids."""
    regular = tuple(
        np.linspace(FROZEN_R_MIN, FROZEN_R_MAX, n, dtype=float)
        for n in FROZEN_GRID_COUNTS
    )
    fine = np.array(regular[-1], copy=True)
    h = (FROZEN_R_MAX - FROZEN_R_MIN) / (fine.size - 1)
    rng = np.random.default_rng(FROZEN_SEED)
    jitter = rng.uniform(-0.24, 0.24, size=fine.size - 2) * h
    shifted = np.array(fine, copy=True)
    shifted[1:-1] += jitter
    if np.any(np.diff(shifted) <= 0.0):
        raise RuntimeError("frozen off-grid radial jitter lost monotonicity")
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
    """Integrate r**e * values with local cubic interpolation and GL8 cells.

    The first cell [0,r_min] uses a cubic extrapolant through the first four
    public positive-radius samples.  Every later cell uses four nearby public
    samples.  This is intentionally distinct from the production quadratic
    first-cell + positive-node trapezoid rule.
    """
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
        if cell == 0:
            start = 0
        else:
            start = max(0, min(cell - 2, r.size - 4))
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
    mean = receipt.get("current_partial_mean_witness")
    if not isinstance(mean, Mapping):
        raise ValueError("receipt is missing current_partial_mean_witness")
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


def _audit_one_receipt(receipt: Mapping[str, Any], *, expected_radii: np.ndarray) -> dict[str, Any]:
    if receipt.get("schema") != PARENT_SCHEMA:
        raise ValueError("unexpected A3 #976 receipt schema")
    if receipt.get("parent_agent3_971_source_blob") != "9c62ae3c3b4b3ca2d9082f822d57feef700ea6cf":
        raise ValueError("A3 #971 provenance drifted")
    geometry = receipt.get("geometry")
    if not isinstance(geometry, Mapping):
        raise ValueError("receipt geometry is missing")
    radii = np.asarray(geometry.get("radii"), dtype=float)
    if radii.shape != expected_radii.shape or not np.allclose(radii, expected_radii, rtol=0.0, atol=2e-15):
        raise ValueError("receipt radial grid does not match preregistered A4 grid")
    if abs(float(geometry.get("time")) - FROZEN_TIME) > 1e-15:
        raise ValueError("receipt time drifted")
    if abs(float(geometry.get("axial_z")) - FROZEN_Z) > 1e-15:
        raise ValueError("receipt axial z drifted")
    center = float(geometry.get("bump_center"))
    halfwidth = float(geometry.get("bump_halfwidth"))
    if abs(center - FROZEN_BUMP_CENTER) > 1e-15 or abs(halfwidth - FROZEN_BUMP_HALFWIDTH) > 1e-15:
        raise ValueError("receipt compact-bump geometry drifted")

    boundary = receipt.get("truth_boundary")
    if not isinstance(boundary, Mapping):
        raise ValueError("receipt truth boundary is missing")
    for key in (
        "complete_ns_defect",
        "scoped_current_partial_radial_stress_authorized_as_correction_target",
        "current_real_ns_correction_velocity_materialized",
        "heldout_normalized_ns_residual_assessed",
        "residual_reduction_claimed",
        "same_protocol_comparable_to_st006",
        "pde_validated",
    ):
        if boundary.get(key) is not False:
            raise ValueError(f"A3 #976 scientific boundary drifted at {key}")

    sources = _extract_sources(receipt)
    n = radii.size
    if any(arr.shape[0] != n for arr in sources.values()):
        raise ValueError("mean sample count does not match radial grid")
    mean_closure = float(np.max(np.abs(sources["aggregate"] - sources["mixed"] - sources["quadratic"])))
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
            if production.shape != radii.shape:
                raise ValueError(f"invalid production stress shape for {channel}/{piece}")
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
            axis_error = float(
                np.max(np.abs(production[axis_mask] - independent["stress"][axis_mask]), initial=0.0)
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


def audit_serialized_receipts(receipts: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    """Audit exactly three nested grids plus one shifted fine-grid receipt."""
    grids = frozen_radial_grids()
    if len(receipts) != len(grids):
        raise ValueError("A4 audit requires coarse/medium/fine/off-grid receipts")
    levels = [
        _audit_one_receipt(receipt, expected_radii=grid)
        for receipt, grid in zip(receipts, grids, strict=True)
    ]
    report = {
        "schema": SCHEMA,
        "task": TASK,
        "parent_agent3_pr": PARENT_AGENT3_PR,
        "parent_agent3_head": PARENT_AGENT3_HEAD,
        "parent_agent3_source_blob": PARENT_AGENT3_SOURCE_BLOB,
        "seed": FROZEN_SEED,
        "grid_counts": list(FROZEN_GRID_COUNTS),
        "offgrid_count": int(grids[-1].size),
        "gauss_legendre_order": FROZEN_GL_ORDER,
        "levels": levels[:3],
        "offgrid": levels[3],
        "axis_evidence_scope": "positive-radius axis-near only; not exact Cartesian-axis evidence",
        "candidate_residual_evidence": False,
        "truth_boundary": truth_boundary(),
    }
    return report


def enforce_audit(report: Mapping[str, Any]) -> None:
    levels = report.get("levels")
    offgrid = report.get("offgrid")
    if not isinstance(levels, Sequence) or len(levels) != 3 or not isinstance(offgrid, Mapping):
        raise ValueError("malformed A4 radial-stress audit report")

    for channel in ("theta_e2", "axial_e1"):
        for piece in ("quadratic", "mixed", "aggregate"):
            series = [float(level["channels"][channel][piece]["stress_relative_rms"]) for level in levels]
            fine = levels[-1]["channels"][channel][piece]
            shifted = offgrid["channels"][channel][piece]
            if float(fine["stress_relative_rms"]) > FINE_STRESS_REL_RMS_GATE:
                raise AssertionError(f"fine stress RMS gate failed for {channel}/{piece}")
            if float(fine["stress_relative_max"]) > FINE_STRESS_REL_MAX_GATE:
                raise AssertionError(f"fine stress max gate failed for {channel}/{piece}")
            if float(fine["weighted_moment_relative_error"]) > FINE_MOMENT_REL_GATE:
                raise AssertionError(f"fine moment gate failed for {channel}/{piece}")
            if float(fine["axis_near_normalized_error"]) > AXIS_NEAR_NORMALIZED_ERROR_GATE:
                raise AssertionError(f"axis-near gate failed for {channel}/{piece}")
            if float(fine["outer_edge_normalized_stress"]) > OUTER_EDGE_NORMALIZED_GATE:
                raise AssertionError(f"outer-edge gate failed for {channel}/{piece}")
            if float(shifted["stress_relative_rms"]) > FINE_STRESS_REL_RMS_GATE:
                raise AssertionError(f"off-grid stress RMS gate failed for {channel}/{piece}")
            if float(shifted["stress_relative_max"]) > FINE_STRESS_REL_MAX_GATE:
                raise AssertionError(f"off-grid stress max gate failed for {channel}/{piece}")
            if float(shifted["weighted_moment_relative_error"]) > FINE_MOMENT_REL_GATE:
                raise AssertionError(f"off-grid moment gate failed for {channel}/{piece}")
            if series[1] > max(RESOLUTION_DEGRADATION_FACTOR * series[0], RESOLUTION_NUMERICAL_FLOOR):
                raise AssertionError(f"coarse->medium stability failed for {channel}/{piece}")
            if series[2] > max(RESOLUTION_DEGRADATION_FACTOR * series[1], RESOLUTION_NUMERICAL_FLOOR):
                raise AssertionError(f"medium->fine stability failed for {channel}/{piece}")

    q_rms = max(
        float(levels[-1]["channels"][channel]["quadratic"]["production_stress_rms"])
        for channel in ("theta_e2", "axial_e1")
    )
    if q_rms <= NONTRIVIAL_STRESS_RMS_FLOOR:
        raise AssertionError("current quadratic radial stress is numerically trivial")
    boundary = report.get("truth_boundary")
    if not isinstance(boundary, Mapping) or boundary.get("pde_validated") is not False:
        raise AssertionError("A4 truth boundary drifted")
    if report.get("candidate_residual_evidence") is not False:
        raise AssertionError("scoped stress audit was laundered into residual evidence")


def mutate_fine_stress_tangent(receipt: Mapping[str, Any]) -> dict[str, Any]:
    """Return an outer-edge-value-preserving stress mutation for firewall tests."""
    out = copy.deepcopy(dict(receipt))
    r = np.asarray(out["geometry"]["radii"], dtype=float)
    shape = r * (r[-1] - r)
    shape /= max(float(np.max(np.abs(shape))), np.finfo(float).tiny)
    for channel in ("theta_e2", "axial_e1"):
        piece = out[channel]["aggregate"]
        stress = np.asarray(piece["stress"], dtype=float)
        scale = max(_rms(stress), 1.0e-6)
        piece["stress"] = (stress + 0.50 * scale * shape).tolist()
    return out


def mutate_bump_geometry(receipt: Mapping[str, Any]) -> dict[str, Any]:
    """Perturb the public bump center without recomputing the serialized stress."""
    out = copy.deepcopy(dict(receipt))
    out["geometry"]["bump_center"] = FROZEN_BUMP_CENTER + 0.02
    return out


def truth_boundary() -> dict[str, Any]:
    return {
        "parent_agent3_pr": PARENT_AGENT3_PR,
        "parent_agent3_head": PARENT_AGENT3_HEAD,
        "independent_reference_uses_parent_radial_inverse": False,
        "independent_reference_uses_public_serialized_mean_samples": True,
        "independent_reference_integrator": "piecewise_cubic_plus_gauss_legendre_8",
        "three_resolution_radial_stress_audit": True,
        "shifted_offgrid_radial_stress_audit": True,
        "positive_radius_axis_near_audit": True,
        "exact_cartesian_axis_audit": False,
        "partial_domain_through_current_Xh_only": True,
        "complete_ns_defect": False,
        "scoped_radial_stress_authorized_as_correction_target": False,
        "matched_global_pressure_materialized": False,
        "restricted_forcing_materialized": False,
        "current_real_ns_correction_velocity_materialized": False,
        "after_correction_residual_assessed": False,
        "heldout_normalized_ns_residual_assessed": False,
        "same_protocol_comparable_to_st006": False,
        "pde_validated": False,
        "final_normalized_momentum_gate": FINAL_NORMALIZED_MOMENTUM_GATE,
        "final_normalized_divergence_gate": FINAL_NORMALIZED_DIVERGENCE_GATE,
    }
