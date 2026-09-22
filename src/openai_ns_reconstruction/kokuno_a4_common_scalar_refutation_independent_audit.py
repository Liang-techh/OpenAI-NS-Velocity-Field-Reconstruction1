"""Independent A4 audit of A1 #1242's one-way common-scalar bridge witness.

K4-VAL-130 is deliberately narrower than Navier--Stokes admission.  A1 #1242
adds a sufficient analytic certificate which may refute one *common scalar*
``l=-h`` matching time if two eta slices have individually certified unique
``Q_s=Q_p`` roots with disjoint brackets.

A4 does not use A1's ``unique_root_certificate``, ``unique_root_bracket`` or
96-bisection path to build its scientific reference.  After a real save/load,
this module first builds an implementation-distinct finite-sample reference
from the public parent ``Q_s(y,eta)`` only:

* deterministic held-out/off-grid eta probes;
* safeguarded regula-falsi root finding, rather than production bisection;
* centered five-point FD4 root-slope diagnostics at three frozen resolutions;
* fixed refinement and target-closure checks.

Only after that independent reference exists is A1's production refutation
report read.  If A1 asserts a disjoint unique-root pair, A4 independently
resolves both eta roots and checks target closure, negative local slope,
production-vs-independent root agreement and independent interval disjointness.
If A1 does not assert a pair, the outcome remains inconclusive: this audit does
not turn absence of a finite witness into evidence for a common bridge.

This is finite scalar-geometry evidence only.  It does not establish full-eta
target totality, global uniqueness/transversality, a smooth eta->T map, an
eta-dependent Cartesian terminal surface, terminal/global velocity, pressure,
forcing, a correction cycle, or a complete NS residual.
"""

from __future__ import annotations

import inspect
import math
from typing import Any, Callable, Mapping

import numpy as np

from .kokuno_current_bridge_common_scalar_refutation import (
    KokunoCurrentBridgeCommonScalarRefutation,
)
from .kokuno_current_exterior_lminus_h_bridge import (
    KokunoCurrentExteriorLMinusHBridge,
)

TASK = "K4-VAL-130"
UPSTREAM_HEAD = "7e987ff0a3bea743fe3f0fe9d02d5b3c8c593ea3"
SEED = 9174011

FD4_Y_STEPS = (2.0**-6, 2.0**-7, 2.0**-8)
ROOT_SNAPSHOT_ITERATIONS = (32, 48, 64)
ROOT_EXPANSIONS = 24
ROOT_EXPANSION_FACTOR = 1.7
ROOT_TARGET_NORM_GATE = 2.0e-10
PRODUCTION_ROOT_DISAGREEMENT_GATE = 2.0e-8
FINE_NORMALIZED_SLOPE_MIN = 1.0e-8
REFINEMENT_WORSEN_FACTOR = 1.25
REFINEMENT_FLOOR = 5.0e-10

FINAL_MOMENTUM_GATE = 1.0e-3
FINAL_DIVERGENCE_GATE = 1.0e-5
FINAL_QUADRATURE = (24, 48, 96)

_FORBIDDEN_PUBLIC_KNOBS = {
    "step",
    "steps",
    "threshold",
    "tol",
    "tolerance",
    "residual",
    "forcing",
    "pressure",
    "viscosity",
    "nu",
    "gain",
    "optimizer",
    "seed",
    "resolution",
    "quadrature",
    "target",
    "bridge_length",
    "root_iterations",
}


def _semantic_value(obj: Any) -> str:
    value = getattr(obj, "semantic_sha256")
    return str(value() if callable(value) else value)


def default_candidate() -> KokunoCurrentBridgeCommonScalarRefutation:
    return KokunoCurrentBridgeCommonScalarRefutation()


def heldout_eta(candidate: KokunoCurrentBridgeCommonScalarRefutation) -> np.ndarray:
    """Return deterministic strict-interior eta probes distinct from A1's grid."""
    lo, hi = (float(v) for v in candidate.eta_interval)
    if not (math.isfinite(lo) and math.isfinite(hi) and lo < hi):
        raise ValueError("candidate eta interval must be finite and nonempty")
    center = 0.5 * (lo + hi)
    half = 0.5 * (hi - lo)
    anchors = np.asarray(
        [-0.947, -0.811, -0.633, -0.421, -0.173, 0.071, 0.287, 0.503, 0.739, 0.913],
        dtype=float,
    )
    rng = np.random.default_rng(SEED)
    random_fraction = rng.uniform(-0.965, 0.965, size=31)
    eta = center + half * np.concatenate([anchors, random_fraction])
    if eta.shape != (41,):
        raise RuntimeError("unexpected held-out eta count")
    if np.any((eta <= lo) | (eta >= hi)):
        raise RuntimeError("held-out eta generation left the strict interior")
    if len(np.unique(eta)) != eta.size:
        raise RuntimeError("held-out eta generation produced duplicates")
    return np.asarray(eta, dtype=float)


def _public_residual(
    parent: KokunoCurrentExteriorLMinusHBridge,
    eta: float,
    y: float,
    q_p: float,
) -> float:
    value = float(np.asarray(parent.state(float(y), float(eta))["Q_s"]))
    out = value - float(q_p)
    if not math.isfinite(out):
        raise RuntimeError("public Q_s target residual became non-finite")
    return out


def _bracket_public_root(
    parent: KokunoCurrentExteriorLMinusHBridge,
    eta: float,
    q_p: float,
) -> tuple[float, float, float, float] | None:
    """Bracket one positive target root using non-production expansion mechanics."""
    lo = 0.0
    f_lo = _public_residual(parent, eta, lo, q_p)
    if f_lo <= 0.0:
        return None
    hi = 0.75
    f_hi = _public_residual(parent, eta, hi, q_p)
    for _ in range(ROOT_EXPANSIONS):
        if f_hi <= 0.0:
            return lo, hi, f_lo, f_hi
        hi *= ROOT_EXPANSION_FACTOR
        f_hi = _public_residual(parent, eta, hi, q_p)
    return None


def _regula_falsi_snapshots(
    residual: Callable[[float], float],
    lo: float,
    hi: float,
    f_lo: float,
    f_hi: float,
) -> dict[int, dict[str, float]]:
    """Safeguarded Illinois regula-falsi snapshots at fixed iteration counts."""
    if not (math.isfinite(lo) and math.isfinite(hi) and lo < hi):
        raise ValueError("invalid root bracket")
    if not (math.isfinite(f_lo) and math.isfinite(f_hi) and f_lo > 0.0 and f_hi <= 0.0):
        raise ValueError("root bracket must be sign-oriented positive-to-nonpositive")

    left, right = float(lo), float(hi)
    fl, fr = float(f_lo), float(f_hi)
    weighted_fl, weighted_fr = fl, fr
    last_side = 0
    snapshots: dict[int, dict[str, float]] = {}
    max_iter = max(ROOT_SNAPSHOT_ITERATIONS)

    for iteration in range(1, max_iter + 1):
        denom = weighted_fr - weighted_fl
        if denom == 0.0 or not math.isfinite(denom):
            x = 0.5 * (left + right)
        else:
            x = (left * weighted_fr - right * weighted_fl) / denom
            if not (left < x < right) or not math.isfinite(x):
                x = 0.5 * (left + right)
            guard = 1.0e-12 * max(1.0, abs(left), abs(right))
            if min(x - left, right - x) < guard:
                x = 0.5 * (left + right)
        fx = float(residual(x))
        if not math.isfinite(fx):
            raise RuntimeError("independent root residual became non-finite")

        if fx > 0.0:
            left, fl = x, fx
            weighted_fl = fl
            if last_side == 1:
                weighted_fr *= 0.5
            else:
                weighted_fr = fr
            last_side = 1
        else:
            right, fr = x, fx
            weighted_fr = fr
            if last_side == -1:
                weighted_fl *= 0.5
            else:
                weighted_fl = fl
            last_side = -1

        if iteration in ROOT_SNAPSHOT_ITERATIONS:
            mid = 0.5 * (left + right)
            f_mid = float(residual(mid))
            snapshots[iteration] = {
                "lower": float(left),
                "upper": float(right),
                "width": float(right - left),
                "midpoint": float(mid),
                "midpoint_residual": float(f_mid),
                "lower_residual": float(fl),
                "upper_residual": float(fr),
            }

    if set(snapshots) != set(ROOT_SNAPSHOT_ITERATIONS):
        raise RuntimeError("independent root solver did not materialize all snapshots")
    return snapshots


def _fd4_slope(
    parent: KokunoCurrentExteriorLMinusHBridge,
    eta: float,
    root: float,
    step: float,
) -> float:
    h = float(step)
    if not (math.isfinite(h) and h > 0.0 and root > 2.0 * h):
        raise ValueError("root is unsafe for frozen FD4 slope stencil")

    def q(offset: float) -> float:
        return float(np.asarray(parent.state(root + offset, eta)["Q_s"]))

    out = (q(-2.0 * h) - 8.0 * q(-h) + 8.0 * q(h) - q(2.0 * h)) / (12.0 * h)
    if not math.isfinite(out):
        raise RuntimeError("independent FD4 slope became non-finite")
    return float(out)


def _normalized_target_closure(value: float, q_p: float) -> float:
    scale = max(1.0, abs(float(value)), abs(float(q_p)))
    return abs(float(value) - float(q_p)) / scale


def _root_record(
    parent: KokunoCurrentExteriorLMinusHBridge,
    eta: float,
    q_p: float,
) -> dict[str, Any] | None:
    bracket = _bracket_public_root(parent, float(eta), float(q_p))
    if bracket is None:
        return None
    lo, hi, f_lo, f_hi = bracket
    residual = lambda y: _public_residual(parent, float(eta), float(y), float(q_p))
    snapshots = _regula_falsi_snapshots(residual, lo, hi, f_lo, f_hi)
    root_by_resolution = [
        snapshots[it]["midpoint"] for it in ROOT_SNAPSHOT_ITERATIONS
    ]
    root = float(root_by_resolution[-1])
    q_at_root = float(np.asarray(parent.state(root, float(eta))["Q_s"]))
    slopes = [
        _fd4_slope(parent, float(eta), root, step) for step in FD4_Y_STEPS
    ]
    normalized_slopes = [
        float(s / max(1.0, abs(q_at_root), abs(q_p))) for s in slopes
    ]
    root_changes = [
        abs(root_by_resolution[1] - root_by_resolution[0]),
        abs(root_by_resolution[2] - root_by_resolution[1]),
    ]
    slope_changes = [
        abs(normalized_slopes[1] - normalized_slopes[0]),
        abs(normalized_slopes[2] - normalized_slopes[1]),
    ]
    return {
        "eta": float(eta),
        "initial_bracket": [float(lo), float(hi)],
        "root_snapshots": {str(k): dict(v) for k, v in snapshots.items()},
        "root_by_resolution": [float(v) for v in root_by_resolution],
        "root_changes": [float(v) for v in root_changes],
        "root_fine": root,
        "q_at_root": q_at_root,
        "target_normalized_closure": _normalized_target_closure(q_at_root, q_p),
        "fd4_steps": list(FD4_Y_STEPS),
        "fd4_slopes": [float(v) for v in slopes],
        "normalized_slopes": [float(v) for v in normalized_slopes],
        "slope_changes": [float(v) for v in slope_changes],
        "fine_normalized_slope": float(normalized_slopes[-1]),
        "final_interval": [
            float(snapshots[ROOT_SNAPSHOT_ITERATIONS[-1]]["lower"]),
            float(snapshots[ROOT_SNAPSHOT_ITERATIONS[-1]]["upper"]),
        ],
        "final_interval_width": float(
            snapshots[ROOT_SNAPSHOT_ITERATIONS[-1]]["width"]
        ),
    }


def _independent_reference(
    loaded: KokunoCurrentBridgeCommonScalarRefutation,
) -> dict[str, Any]:
    """Form the A4 reference without reading any A1 witness/certificate method."""
    eta = heldout_eta(loaded)
    parent = loaded.parent
    q_p = float(loaded.q_p)
    records: list[dict[str, Any]] = []
    missing: list[float] = []
    for value in eta:
        record = _root_record(parent, float(value), q_p)
        if record is None:
            missing.append(float(value))
        else:
            records.append(record)

    roots = np.asarray([r["root_fine"] for r in records], dtype=float)
    spread = float(np.max(roots) - np.min(roots)) if roots.size >= 2 else 0.0
    return {
        "heldout_seed": SEED,
        "eta_count": int(eta.size),
        "root_found_count": int(len(records)),
        "eta_without_positive_bracket": missing,
        "records": records,
        "finite_sample_root_time_min": float(np.min(roots)) if roots.size else None,
        "finite_sample_root_time_max": float(np.max(roots)) if roots.size else None,
        "finite_sample_root_time_spread": spread,
        "finite_sample_spread_is_diagnostic_only": True,
        "does_not_establish_global_uniqueness": True,
    }


def _pair_audit_after_reference(
    loaded: KokunoCurrentBridgeCommonScalarRefutation,
    production_report: Mapping[str, Any],
) -> dict[str, Any]:
    pair = production_report.get("disjoint_root_pair")
    claimed = bool(
        production_report.get(
            "common_scalar_bridge_refuted_by_disjoint_unique_root_brackets", False
        )
    )
    if not claimed:
        if pair is not None:
            raise RuntimeError("production report exposes a pair while refutation flag is false")
        return {
            "production_refutation_claim_present": False,
            "independent_pair_audit_required": False,
            "absence_of_production_pair_is_inconclusive": True,
        }
    if not isinstance(pair, Mapping):
        raise RuntimeError("production refutation flag lacks a disjoint-root pair")

    eta_a = float(pair["eta_earlier"])
    eta_b = float(pair["eta_later"])
    rec_a = _root_record(loaded.parent, eta_a, float(loaded.q_p))
    rec_b = _root_record(loaded.parent, eta_b, float(loaded.q_p))
    if rec_a is None or rec_b is None:
        raise RuntimeError("independent solver failed to bracket a production witness root")

    prod_a = [float(v) for v in pair["earlier_root_bracket"]]
    prod_b = [float(v) for v in pair["later_root_bracket"]]
    prod_mid_a = 0.5 * (prod_a[0] + prod_a[1])
    prod_mid_b = 0.5 * (prod_b[0] + prod_b[1])
    disagree_a = abs(float(rec_a["root_fine"]) - prod_mid_a)
    disagree_b = abs(float(rec_b["root_fine"]) - prod_mid_b)

    ind_a = [float(v) for v in rec_a["final_interval"]]
    ind_b = [float(v) for v in rec_b["final_interval"]]
    if rec_a["root_fine"] <= rec_b["root_fine"]:
        early, late = rec_a, rec_b
        early_interval, late_interval = ind_a, ind_b
    else:
        early, late = rec_b, rec_a
        early_interval, late_interval = ind_b, ind_a
    independent_gap = float(late_interval[0] - early_interval[1])

    return {
        "production_refutation_claim_present": True,
        "independent_pair_audit_required": True,
        "eta_pair": [eta_a, eta_b],
        "first_root": rec_a,
        "second_root": rec_b,
        "production_root_midpoints": [prod_mid_a, prod_mid_b],
        "production_vs_independent_root_abs_disagreement": [
            float(disagree_a),
            float(disagree_b),
        ],
        "independent_final_intervals_disjoint": bool(independent_gap > 0.0),
        "independent_disjoint_gap": independent_gap,
        "production_disjoint_gap_lower_bound": float(
            pair["disjoint_gap_lower_bound"]
        ),
        "finite_pair_supports_only_one_way_refutation": True,
    }


def _refinement_ok(changes: list[float]) -> bool:
    if len(changes) != 2:
        return False
    coarse_medium, medium_fine = (float(v) for v in changes)
    return medium_fine <= max(
        REFINEMENT_FLOOR,
        REFINEMENT_WORSEN_FACTOR * coarse_medium,
    )


def audit_loaded_common_scalar_refutation(
    loaded: KokunoCurrentBridgeCommonScalarRefutation,
    pre_serialization_reference: KokunoCurrentBridgeCommonScalarRefutation,
) -> dict[str, Any]:
    """Audit exact A1 #1242 after save/load with an independent numerical path."""
    if not isinstance(loaded, KokunoCurrentBridgeCommonScalarRefutation):
        raise TypeError("loaded must be KokunoCurrentBridgeCommonScalarRefutation")
    if not isinstance(pre_serialization_reference, KokunoCurrentBridgeCommonScalarRefutation):
        raise TypeError("pre_serialization_reference must match the candidate type")

    # Scientific reference first.  This function deliberately does not call
    # A1 witness/certificate methods until _independent_reference is complete.
    independent = _independent_reference(loaded)

    production = loaded.refutation_report()
    pair_audit = _pair_audit_after_reference(loaded, production)

    probe_eta = heldout_eta(loaded)[:7]
    loaded_q = np.asarray(loaded.parent.entry_state(probe_eta)["Q_s_entry"], dtype=float)
    pre_q = np.asarray(
        pre_serialization_reference.parent.entry_state(probe_eta)["Q_s_entry"],
        dtype=float,
    )
    save_load = {
        "semantic_exact": _semantic_value(loaded) == _semantic_value(pre_serialization_reference),
        "configuration_exact": loaded.configuration() == pre_serialization_reference.configuration(),
        "entry_q_exact": bool(np.array_equal(loaded_q, pre_q)),
    }

    truth = dict(loaded.truth_boundary)
    truth.update(
        {
            "current_common_scalar_refutation_independent_audit_executed": True,
            "full_eta_interval_entry_above_qp_established": False,
            "full_eta_interval_target_totality_established": False,
            "full_eta_interval_root_uniqueness_established": False,
            "full_eta_interval_transversality_established": False,
            "smooth_eta_target_time_map_established": False,
            "full_eta_interval_common_scalar_bridge_length_established": False,
            "current_l_minus_h_matching_bridge_materialized": False,
            "eta_dependent_cartesian_matching_boundary_materialized": False,
            "cartesian_terminal_multiplier_composed": False,
            "leading_only_ns_residual_assessed": False,
            "leading_plus_oscillatory_ns_residual_assessed": False,
            "after_correction_ns_residual_assessed": False,
            "complete_ns_admission_ready": False,
            "pde_validated": False,
        }
    )

    return {
        "schema": "kokuno-a4-common-scalar-refutation-independent-audit-v1",
        "task": TASK,
        "upstream_head": UPSTREAM_HEAD,
        "candidate_semantic_sha256": _semantic_value(loaded),
        "parent_semantic_sha256": _semantic_value(loaded.parent),
        "q_p": float(loaded.q_p),
        "heldout_seed": SEED,
        "fd4_y_steps": list(FD4_Y_STEPS),
        "root_snapshot_iterations": list(ROOT_SNAPSHOT_ITERATIONS),
        "independent_reference": independent,
        "production_report_read_only_after_reference": True,
        "production_report": production,
        "production_pair_independent_audit": pair_audit,
        "save_load": save_load,
        "gates": {
            "root_normalized_target_closure": ROOT_TARGET_NORM_GATE,
            "production_root_disagreement": PRODUCTION_ROOT_DISAGREEMENT_GATE,
            "fine_normalized_slope_minimum_magnitude": FINE_NORMALIZED_SLOPE_MIN,
            "refinement_worsen_factor": REFINEMENT_WORSEN_FACTOR,
            "refinement_floor": REFINEMENT_FLOOR,
            "final_project_momentum_gate_unchanged": FINAL_MOMENTUM_GATE,
            "final_project_divergence_gate_unchanged": FINAL_DIVERGENCE_GATE,
            "final_project_quadrature_unchanged": list(FINAL_QUADRATURE),
        },
        "truth_boundary": truth,
        "final_project_admission_ready": False,
    }


def _enforce_root_record(record: Mapping[str, Any]) -> None:
    assert float(record["target_normalized_closure"]) <= ROOT_TARGET_NORM_GATE
    fine_slope = float(record["fine_normalized_slope"])
    assert fine_slope < 0.0
    assert abs(fine_slope) >= FINE_NORMALIZED_SLOPE_MIN
    assert _refinement_ok(list(record["root_changes"]))
    assert _refinement_ok(list(record["slope_changes"]))
    interval = list(record["final_interval"])
    assert len(interval) == 2 and float(interval[0]) < float(interval[1])
    assert float(record["final_interval_width"]) >= 0.0


def enforce_preregistered_gates(report: Mapping[str, Any]) -> None:
    """Enforce K4-VAL-130 mechanics without promoting continuum/PDE truth."""
    assert report["task"] == TASK
    assert report["upstream_head"] == UPSTREAM_HEAD
    assert report["heldout_seed"] == SEED
    assert report["fd4_y_steps"] == list(FD4_Y_STEPS)
    assert report["root_snapshot_iterations"] == list(ROOT_SNAPSHOT_ITERATIONS)

    save_load = report["save_load"]
    assert save_load["semantic_exact"] is True
    assert save_load["configuration_exact"] is True
    assert save_load["entry_q_exact"] is True

    independent = report["independent_reference"]
    assert int(independent["eta_count"]) == 41
    assert int(independent["root_found_count"]) >= 1
    for record in independent["records"]:
        _enforce_root_record(record)

    pair = report["production_pair_independent_audit"]
    if bool(pair["production_refutation_claim_present"]):
        _enforce_root_record(pair["first_root"])
        _enforce_root_record(pair["second_root"])
        for value in pair["production_vs_independent_root_abs_disagreement"]:
            assert float(value) <= PRODUCTION_ROOT_DISAGREEMENT_GATE
        assert pair["independent_final_intervals_disjoint"] is True
        assert float(pair["independent_disjoint_gap"]) > 0.0
    else:
        assert pair["independent_pair_audit_required"] is False
        assert pair["absence_of_production_pair_is_inconclusive"] is True

    truth = report["truth_boundary"]
    assert truth["current_common_scalar_refutation_independent_audit_executed"] is True
    assert truth["full_eta_interval_entry_above_qp_established"] is False
    assert truth["full_eta_interval_target_totality_established"] is False
    assert truth["full_eta_interval_root_uniqueness_established"] is False
    assert truth["full_eta_interval_transversality_established"] is False
    assert truth["smooth_eta_target_time_map_established"] is False
    assert truth["full_eta_interval_common_scalar_bridge_length_established"] is False
    assert truth["current_l_minus_h_matching_bridge_materialized"] is False
    assert truth["eta_dependent_cartesian_matching_boundary_materialized"] is False
    assert truth["cartesian_terminal_multiplier_composed"] is False
    assert truth["complete_ns_admission_ready"] is False
    assert truth["pde_validated"] is False, "PDE truth promotion is forbidden"
    assert report["final_project_admission_ready"] is False
    assert report["gates"]["final_project_momentum_gate_unchanged"] == 1.0e-3
    assert report["gates"]["final_project_divergence_gate_unchanged"] == 1.0e-5
    assert report["gates"]["final_project_quadrature_unchanged"] == [24, 48, 96]


def public_api_has_no_scientific_tuning_knobs() -> bool:
    sig = inspect.signature(audit_loaded_common_scalar_refutation)
    names = {name.lower() for name in sig.parameters}
    return not any(any(token in name for token in _FORBIDDEN_PUBLIC_KNOBS) for name in names)
