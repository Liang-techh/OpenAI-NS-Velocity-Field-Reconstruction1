"""Independent A4 audit of the public/source-ideal exterior Q_s schedule.

K4-VAL-127 is intentionally a scoped numerical-reference check, not a
Navier--Stokes admission.  The upstream A1 #1217 surface explicitly labels its
Q_s transport as SOURCE-IDEAL because it assumes I=XH/(1-lambda), U=M=J=0,
while the current repository candidate carries a nonzero M/M_eta history.

The A4 reference below therefore checks only that source-ideal scalar schedule.
It reimplements the public flat step and solves

    Q' + (1+l) Q = -l-h

with classical RK4 on three frozen resolutions.  It does not call A1's
_sigma/_sigma_integral/_integral_0_upper or its GL192 construction to build the
reference.  Production values are read only after all independent values have
been formed, for comparison.

Passing this audit must never be used to choose a current-candidate matching
bridge, extend Cartesian velocity, or promote pde_validated.
"""

from __future__ import annotations

import inspect
import math
from typing import Any, Callable, Mapping

import numpy as np

from .kokuno_public_ideal_exterior_qs_schedule import (
    KokunoPublicIdealExteriorQsSchedule,
)

TASK = "K4-VAL-127"
UPSTREAM_HEAD = "3978560078104bff7c9d5556bd9e658e4874c963"
SEED = 9173981
RK4_STEPS_PER_UNIT = (2048, 4096, 8192)
ENDPOINT_PUBLIC_REL_GATE = 1.0e-8
OFFGRID_PUBLIC_REL_GATE = 1.0e-8
MEDIUM_FINE_ENDPOINT_REL_GATE = 2.0e-9
HOLD_LINEAR_REL_GATE = 2.0e-12
BRIDGE_ENDPOINT_REL_GATE = 2.0e-12
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
}


def _semantic_value(obj: Any) -> str:
    value = getattr(obj, "semantic_sha256")
    return str(value() if callable(value) else value)


def _sigma_ref(s: Any) -> np.ndarray:
    """A4-owned stable realization of the public C-infinity flat step."""
    arr = np.asarray(s, dtype=float)
    if np.any(~np.isfinite(arr)):
        raise ValueError("sigma argument must be finite")
    out = np.zeros_like(arr, dtype=float)
    out[arr >= 1.0] = 1.0
    inside = (arr > 0.0) & (arr < 1.0)
    if np.any(inside):
        x = arr[inside]
        logit = -1.0 / (x * x) + 1.0 / ((1.0 - x) ** 2)
        vals = np.empty_like(logit)
        positive = logit >= 0.0
        if np.any(positive):
            e = np.exp(-logit[positive])
            vals[positive] = 1.0 / (1.0 + e)
        if np.any(~positive):
            e = np.exp(logit[~positive])
            vals[~positive] = e / (1.0 + e)
        out[inside] = vals
    return out


def _release1_rhs_ref(s: float, q: float, lam: float, h: float) -> float:
    sig = float(_sigma_ref(np.asarray(s)))
    ell = -lam - (1.0 - lam) * sig
    return -ell - h - (1.0 + ell) * q


def _release2_rhs_ref(s: float, q: float, lam: float, h: float) -> float:
    del lam
    sig = float(_sigma_ref(np.asarray(s)))
    ell = -1.0 + (1.0 - h) * sig
    return -ell - h - (1.0 + ell) * q


def _rk4_to(
    rhs: Callable[[float, float, float, float], float],
    q0: float,
    upper: float,
    n_per_unit: int,
    lam: float,
    h: float,
) -> float:
    """Integrate to one arbitrary point with a frozen max step 1/n_per_unit."""
    if n_per_unit <= 0:
        raise ValueError("n_per_unit must be positive")
    if not (math.isfinite(upper) and 0.0 <= upper <= 1.0):
        raise ValueError("release coordinate must lie in [0,1]")
    q = float(q0)
    if not math.isfinite(q):
        raise ValueError("initial Q must be finite")
    dt = 1.0 / float(n_per_unit)
    full_steps = int(math.floor(upper / dt + 1.0e-13))
    # Protect exact upper=1 from a one-step overshoot caused by roundoff.
    full_steps = min(full_steps, n_per_unit)
    s = 0.0
    for _ in range(full_steps):
        k1 = rhs(s, q, lam, h)
        k2 = rhs(s + 0.5 * dt, q + 0.5 * dt * k1, lam, h)
        k3 = rhs(s + 0.5 * dt, q + 0.5 * dt * k2, lam, h)
        k4 = rhs(s + dt, q + dt * k3, lam, h)
        q += (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        s += dt
    remainder = float(upper) - s
    if remainder > 16.0 * np.finfo(float).eps:
        d = remainder
        k1 = rhs(s, q, lam, h)
        k2 = rhs(s + 0.5 * d, q + 0.5 * d * k1, lam, h)
        k3 = rhs(s + 0.5 * d, q + 0.5 * d * k2, lam, h)
        k4 = rhs(s + d, q + d * k3, lam, h)
        q += (d / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
    if not (math.isfinite(q) and q > 0.0):
        raise RuntimeError("independent RK4 source-ideal Q became invalid")
    return q


def _relative_difference(a: float, b: float) -> float:
    scale = max(abs(float(a)), abs(float(b)), np.finfo(float).tiny)
    return abs(float(a) - float(b)) / scale


def _max_relative(a: np.ndarray, b: np.ndarray) -> float:
    aa = np.asarray(a, dtype=float)
    bb = np.asarray(b, dtype=float)
    if aa.shape != bb.shape:
        raise ValueError("relative-comparison shapes differ")
    scale = np.maximum(np.maximum(np.abs(aa), np.abs(bb)), np.finfo(float).tiny)
    return float(np.max(np.abs(aa - bb) / scale)) if aa.size else 0.0


def _heldout_coordinates() -> dict[str, np.ndarray]:
    """Deterministic off-grid points, fixed before any production comparison."""
    rng = np.random.default_rng(SEED)
    release1 = np.sort(0.035 + 0.93 * rng.random(10))
    release2 = np.sort(0.035 + 0.93 * rng.random(10))
    hold_fraction = np.sort(0.025 + 0.95 * rng.random(8))
    return {
        "release1": release1,
        "release2": release2,
        "hold_fraction": hold_fraction,
    }


def _independent_realization(
    lam: float,
    h: float,
    q_p: float,
    n_per_unit: int,
    heldout: Mapping[str, np.ndarray],
) -> dict[str, Any]:
    q0 = (lam - h) / (1.0 - lam)
    if not (math.isfinite(q0) and q0 > 0.0):
        raise RuntimeError("independent source-ideal exterior entry Q is invalid")

    q1 = _rk4_to(_release1_rhs_ref, q0, 1.0, n_per_unit, lam, h)
    hold_length = 4.0 * math.log(1.0 / h)
    qh = q1 + (1.0 - h) * hold_length
    qin = _rk4_to(_release2_rhs_ref, qh, 1.0, n_per_unit, lam, h)
    tmatch = math.log(qin / q_p) / (1.0 - h)
    if not (math.isfinite(tmatch) and tmatch > 0.0):
        raise RuntimeError("independent source-ideal matching length is invalid")

    r1_points = np.asarray(heldout["release1"], dtype=float)
    r2_points = np.asarray(heldout["release2"], dtype=float)
    hold_fraction = np.asarray(heldout["hold_fraction"], dtype=float)
    r1_values = np.asarray(
        [_rk4_to(_release1_rhs_ref, q0, float(s), n_per_unit, lam, h) for s in r1_points]
    )
    r2_values = np.asarray(
        [_rk4_to(_release2_rhs_ref, qh, float(s), n_per_unit, lam, h) for s in r2_points]
    )
    hold_s = hold_fraction * hold_length
    hold_values = q1 + (1.0 - h) * hold_s
    bridge_endpoint = qin * math.exp(-(1.0 - h) * tmatch)

    return {
        "n_per_unit": int(n_per_unit),
        "q0": float(q0),
        "q_release1_end": float(q1),
        "hold_length": float(hold_length),
        "q_hold_end": float(qh),
        "q_in": float(qin),
        "matching_length": float(tmatch),
        "bridge_endpoint": float(bridge_endpoint),
        "release1_points": r1_points,
        "release1_values": r1_values,
        "release2_points": r2_points,
        "release2_values": r2_values,
        "hold_s": hold_s,
        "hold_values": hold_values,
    }


def audit_loaded_public_ideal_schedule(
    loaded: KokunoPublicIdealExteriorQsSchedule,
    pre_serialization_reference: KokunoPublicIdealExteriorQsSchedule,
) -> dict[str, Any]:
    """Run K4-VAL-127 after a real A1 save/load round trip."""
    if not isinstance(loaded, KokunoPublicIdealExteriorQsSchedule):
        raise TypeError("loaded must be KokunoPublicIdealExteriorQsSchedule")
    if not isinstance(pre_serialization_reference, KokunoPublicIdealExteriorQsSchedule):
        raise TypeError("pre_serialization_reference must be KokunoPublicIdealExteriorQsSchedule")

    loaded_semantic = _semantic_value(loaded)
    reference_semantic = _semantic_value(pre_serialization_reference)
    if loaded_semantic != reference_semantic:
        raise ValueError("source-ideal schedule semantic identity changed across save/load")
    config = loaded.configuration()
    if config != pre_serialization_reference.configuration():
        raise ValueError("source-ideal schedule configuration changed across save/load")
    truth = config.get("truth_boundary", {})
    required_false = (
        "current_q_s_release2_endpoint_materialized",
        "current_l_minus_h_matching_bridge_materialized",
        "current_cartesian_terminal_multiplier_composed",
        "outer_global_leading_velocity_materialized",
        "heldout_ns_residual_assessed",
        "pde_validated",
    )
    for key in required_false:
        if truth.get(key) is not False:
            raise ValueError(f"upstream source-ideal schedule illegally promotes {key}")

    lam = float(config["lambda_current"])
    h = float(config["h_current"])
    if not (0.0 < h < 0.01 and 0.0 < lam < 1.0 and h < lam):
        raise ValueError("serialized lambda/h leave the frozen public construction range")

    # Q_p is an already-public boundary value audited separately by K4-VAL-126.
    # It is not used to build release1/release2 Q; only the final matching length.
    q_p = float(loaded.q_p)
    if not (0.0 < q_p < h):
        raise ValueError("public terminal target must satisfy 0<Q_p<h")

    heldout = _heldout_coordinates()

    # Form ALL independent values before reading any A1 production Q_s values.
    independent = {
        str(n): _independent_realization(lam, h, q_p, n, heldout)
        for n in RK4_STEPS_PER_UNIT
    }

    # Production paths are comparison targets only after the A4 reference exists.
    public = loaded.schedule_report()
    public_r1 = np.asarray(loaded.q_release1(heldout["release1"]), dtype=float)
    public_r2 = np.asarray(loaded.q_release2(heldout["release2"]), dtype=float)
    public_hold_s = heldout["hold_fraction"] * float(public["hold_length"])
    public_hold = np.asarray(loaded.q_hold(public_hold_s), dtype=float)

    by_resolution: dict[str, Any] = {}
    endpoint_names = (
        ("q_release1_end", "q_release1_end_source_ideal"),
        ("q_hold_end", "q_hold_end_source_ideal"),
        ("q_in", "q_in_source_ideal"),
        ("matching_length", "ideal_matching_length"),
    )
    for n in RK4_STEPS_PER_UNIT:
        ref = independent[str(n)]
        endpoint_rel = {
            local: _relative_difference(float(ref[local]), float(public[upstream]))
            for local, upstream in endpoint_names
        }
        by_resolution[str(n)] = {
            "endpoints": {
                "q0": float(ref["q0"]),
                "q_release1_end": float(ref["q_release1_end"]),
                "q_hold_end": float(ref["q_hold_end"]),
                "q_in": float(ref["q_in"]),
                "matching_length": float(ref["matching_length"]),
                "bridge_endpoint": float(ref["bridge_endpoint"]),
            },
            "endpoint_relative_differences_to_public": endpoint_rel,
            "endpoint_relative_max_to_public": max(endpoint_rel.values()),
            "release1_offgrid_relative_max_to_public": _max_relative(
                ref["release1_values"], public_r1
            ),
            "release2_offgrid_relative_max_to_public": _max_relative(
                ref["release2_values"], public_r2
            ),
            "hold_linear_relative_max_to_public": _max_relative(
                ref["hold_values"], public_hold
            ),
            "bridge_endpoint_relative_error_to_q_p": _relative_difference(
                ref["bridge_endpoint"], q_p
            ),
        }

    medium = independent[str(RK4_STEPS_PER_UNIT[1])]
    fine = independent[str(RK4_STEPS_PER_UNIT[2])]
    medium_fine = {
        name: _relative_difference(float(medium[name]), float(fine[name]))
        for name in ("q_release1_end", "q_hold_end", "q_in", "matching_length")
    }

    fine_metrics = by_resolution[str(RK4_STEPS_PER_UNIT[-1])]
    fine_offgrid_max = max(
        float(fine_metrics["release1_offgrid_relative_max_to_public"]),
        float(fine_metrics["release2_offgrid_relative_max_to_public"]),
    )

    return {
        "schema": "kokuno-a4-public-ideal-qs-independent-audit-v1",
        "task": TASK,
        "upstream_head": UPSTREAM_HEAD,
        "candidate_semantic_sha256": loaded_semantic,
        "save_load_semantic_exact": True,
        "independent_operator": {
            "flat_step": "A4-owned stable public sigma reimplementation",
            "ode_solver": "classical explicit RK4",
            "steps_per_unit": list(RK4_STEPS_PER_UNIT),
            "seed": SEED,
            "production_paths_excluded_from_reference": [
                "A1 _sigma",
                "A1 _sigma_integral",
                "A1 _integral_0_upper",
                "A1 q_release1",
                "A1 q_hold",
                "A1 q_release2",
                "A1 q_matching_bridge_source_ideal",
                "A1 GL192 nodes/weights",
            ],
        },
        "parameters": {"lambda": lam, "h": h, "q_p": q_p},
        "heldout": {
            "release1_points": [float(v) for v in heldout["release1"]],
            "release2_points": [float(v) for v in heldout["release2"]],
            "hold_fractions": [float(v) for v in heldout["hold_fraction"]],
        },
        "public_endpoints": {
            "q0": float(public["q_exterior_entry_source_ideal"]),
            "q_release1_end": float(public["q_release1_end_source_ideal"]),
            "q_hold_end": float(public["q_hold_end_source_ideal"]),
            "q_in": float(public["q_in_source_ideal"]),
            "matching_length": float(public["ideal_matching_length"]),
        },
        "by_resolution": by_resolution,
        "medium_to_fine_endpoint_relative_changes": medium_fine,
        "medium_to_fine_endpoint_relative_max": max(medium_fine.values()),
        "fine_endpoint_relative_max_to_public": float(
            fine_metrics["endpoint_relative_max_to_public"]
        ),
        "fine_offgrid_release_relative_max_to_public": fine_offgrid_max,
        "fine_hold_linear_relative_max_to_public": float(
            fine_metrics["hold_linear_relative_max_to_public"]
        ),
        "fine_bridge_endpoint_relative_error_to_q_p": float(
            fine_metrics["bridge_endpoint_relative_error_to_q_p"]
        ),
        "truth_boundary": {
            "source_ideal_exterior_qs_numerically_audited": True,
            "current_candidate_q_s_audited": False,
            "current_q_s_release2_endpoint_materialized": False,
            "current_l_minus_h_matching_bridge_materialized": False,
            "current_cartesian_bridge_or_terminal_velocity_materialized": False,
            "leading_only_ns_residual_assessed": False,
            "leading_plus_oscillatory_ns_residual_assessed": False,
            "after_correction_ns_residual_assessed": False,
            "complete_ns_residual_assessed": False,
            "pde_validated": False,
        },
        "gates": {
            "fine_endpoint_public_relative": ENDPOINT_PUBLIC_REL_GATE,
            "fine_offgrid_public_relative": OFFGRID_PUBLIC_REL_GATE,
            "medium_fine_endpoint_relative": MEDIUM_FINE_ENDPOINT_REL_GATE,
            "hold_linear_relative": HOLD_LINEAR_REL_GATE,
            "bridge_endpoint_relative": BRIDGE_ENDPOINT_REL_GATE,
            "final_project_momentum_gate_unchanged": FINAL_MOMENTUM_GATE,
            "final_project_divergence_gate_unchanged": FINAL_DIVERGENCE_GATE,
            "final_project_quadrature_unchanged": list(FINAL_QUADRATURE),
        },
        "final_project_admission_ready": False,
    }


def enforce_preregistered_gates(report: Mapping[str, Any]) -> None:
    assert float(report["fine_endpoint_relative_max_to_public"]) <= ENDPOINT_PUBLIC_REL_GATE
    assert float(report["fine_offgrid_release_relative_max_to_public"]) <= OFFGRID_PUBLIC_REL_GATE
    assert float(report["medium_to_fine_endpoint_relative_max"]) <= MEDIUM_FINE_ENDPOINT_REL_GATE
    assert float(report["fine_hold_linear_relative_max_to_public"]) <= HOLD_LINEAR_REL_GATE
    assert float(report["fine_bridge_endpoint_relative_error_to_q_p"]) <= BRIDGE_ENDPOINT_REL_GATE
    assert report["save_load_semantic_exact"] is True
    truth = report["truth_boundary"]
    assert truth["source_ideal_exterior_qs_numerically_audited"] is True
    assert truth["current_candidate_q_s_audited"] is False
    assert truth["current_q_s_release2_endpoint_materialized"] is False
    assert truth["current_l_minus_h_matching_bridge_materialized"] is False
    assert truth["complete_ns_residual_assessed"] is False
    assert truth["pde_validated"] is False
    assert report["final_project_admission_ready"] is False
    gates = report["gates"]
    assert float(gates["final_project_momentum_gate_unchanged"]) == FINAL_MOMENTUM_GATE
    assert float(gates["final_project_divergence_gate_unchanged"]) == FINAL_DIVERGENCE_GATE
    assert tuple(gates["final_project_quadrature_unchanged"]) == FINAL_QUADRATURE


def public_api_has_no_scientific_tuning_knobs() -> bool:
    params = inspect.signature(audit_loaded_public_ideal_schedule).parameters
    lowered = {name.lower() for name in params}
    return not any(name in _FORBIDDEN_PUBLIC_KNOBS for name in lowered)


def default_schedule() -> KokunoPublicIdealExteriorQsSchedule:
    return KokunoPublicIdealExteriorQsSchedule()
