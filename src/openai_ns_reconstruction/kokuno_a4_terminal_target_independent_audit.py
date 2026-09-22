"""Kokuno Agent 4 independent audit of the A1 terminal target algebra.

K4-VAL-126 is deliberately narrower than PDE validation.  It consumes the
save/reloaded public A1 #1213 target configuration, reconstructs the public
flat-step formulas on an A4-owned code path, and recomputes Q_p with a
composite-Simpson integration-by-parts reference distinct from A1's GL192
direct/reduced implementations.

No Cartesian velocity, pressure, forcing, residual, or optimization quantity
is introduced here.  Passing this audit cannot promote pde_validated.
"""

from __future__ import annotations

import inspect
import math
from typing import Any, Mapping

import numpy as np

from .kokuno_public_terminal_multiplier_target import (
    KokunoPublicTerminalMultiplierTarget,
)

TASK = "K4-VAL-126"
UPSTREAM_HEAD = "fb81b8e90272caa9e8bbd075f05b94dbb1195310"
SEED = 9173971
SIMPSON_SUBINTERVALS = (4096, 8192, 16384)
TERMINAL_LENGTH = 3.0
Q_FINE_PUBLIC_REL_GATE = 2.0e-10
Q_MEDIUM_FINE_REL_GATE = 5.0e-11
PUBLIC_DIRECT_REDUCED_REL_GATE = 2.0e-12
SLOPE_NEGATIVE_FLOOR = 1.0e-15
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


def _sigma_ref(z: Any) -> np.ndarray:
    """Independent stable realization of the public C-infinity flat step."""
    arr = np.asarray(z, dtype=float)
    if np.any(~np.isfinite(arr)):
        raise ValueError("sigma argument must be finite")
    out = np.zeros_like(arr, dtype=float)
    out[arr >= 1.0] = 1.0
    inside = (arr > 0.0) & (arr < 1.0)
    if np.any(inside):
        x = arr[inside]
        logit = -1.0 / (x * x) + 1.0 / ((1.0 - x) ** 2)
        vals = np.empty_like(logit)
        nonnegative = logit >= 0.0
        if np.any(nonnegative):
            e = np.exp(-logit[nonnegative])
            vals[nonnegative] = 1.0 / (1.0 + e)
        if np.any(~nonnegative):
            e = np.exp(logit[~nonnegative])
            vals[~nonnegative] = e / (1.0 + e)
        out[inside] = vals
    return out


def _sigma_prime_ref(z: Any) -> np.ndarray:
    """A4-owned analytic derivative of the public flat step."""
    arr = np.asarray(z, dtype=float)
    if np.any(~np.isfinite(arr)):
        raise ValueError("sigma derivative argument must be finite")
    out = np.zeros_like(arr, dtype=float)
    inside = (arr > 0.0) & (arr < 1.0)
    if np.any(inside):
        x = arr[inside]
        sig = _sigma_ref(x)
        out[inside] = sig * (1.0 - sig) * (
            2.0 / (x ** 3) + 2.0 / ((1.0 - x) ** 3)
        )
    return out


def _f_ref(y: Any, rho: float) -> np.ndarray:
    yy = np.asarray(y, dtype=float)
    z = 0.5 * (yy - 1.0)
    return 1.0 - float(rho) * (1.0 - _sigma_ref(z))


def _f_prime_ref(y: Any, rho: float) -> np.ndarray:
    yy = np.asarray(y, dtype=float)
    z = 0.5 * (yy - 1.0)
    return 0.5 * float(rho) * _sigma_prime_ref(z)


def _simpson_integral(values: np.ndarray, length: float, n: int) -> float:
    if n <= 0 or n % 2:
        raise ValueError("composite Simpson requires a positive even subinterval count")
    vals = np.asarray(values, dtype=float)
    if vals.shape != (n + 1,):
        raise ValueError("Simpson value array has the wrong shape")
    dx = float(length) / float(n)
    return (dx / 3.0) * float(
        vals[0]
        + vals[-1]
        + 4.0 * np.sum(vals[1:-1:2])
        + 2.0 * np.sum(vals[2:-2:2])
    )


def _q_p_ibp_reference(h: float, rho: float, n: int) -> float:
    """Independent Q_p via integration by parts, integrating f rather than f'."""
    x = np.linspace(0.0, TERMINAL_LENGTH, n + 1, dtype=float)
    f = _f_ref(x, rho)
    a = 1.0 - float(h)
    weighted = np.exp(a * x) * f
    integral = _simpson_integral(weighted, TERMINAL_LENGTH, n)
    f0 = float(_f_ref(np.asarray(0.0), rho))
    f3 = float(_f_ref(np.asarray(TERMINAL_LENGTH), rho))
    q = (math.exp(a * TERMINAL_LENGTH) * f3 - f0 - a * integral) / f0
    if not (math.isfinite(q) and q > 0.0):
        raise RuntimeError("independent Q_p reference must be finite and positive")
    return q


def _relative_difference(a: float, b: float) -> float:
    scale = max(abs(float(a)), abs(float(b)), np.finfo(float).tiny)
    return abs(float(a) - float(b)) / scale


def _slope_audit(h: float, rho: float) -> dict[str, Any]:
    """Deterministic off-grid audit of 0 <= f'/f < h/4."""
    rng = np.random.default_rng(SEED)
    random_points = 1.0 + 2.0 * rng.random(8192)
    # Include the symmetry center, but keep the bulk of the set independently off-grid.
    y = np.concatenate((random_points, np.asarray([2.0], dtype=float)))
    f = _f_ref(y, rho)
    fp = _f_prime_ref(y, rho)
    slope = fp / f
    if np.any(~np.isfinite(slope)):
        raise RuntimeError("independent slope audit became non-finite")
    index = int(np.argmax(slope))
    return {
        "seed": SEED,
        "point_count": int(y.size),
        "minimum": float(np.min(slope)),
        "maximum": float(slope[index]),
        "maximum_y": float(y[index]),
        "h_over_4": float(h / 4.0),
        "strict_upper_margin": float(h / 4.0 - slope[index]),
    }


def audit_loaded_terminal_target(
    loaded: KokunoPublicTerminalMultiplierTarget,
    pre_serialization_reference: KokunoPublicTerminalMultiplierTarget,
) -> dict[str, Any]:
    """Run K4-VAL-126 after a real A1 save/load round trip."""
    if not isinstance(loaded, KokunoPublicTerminalMultiplierTarget):
        raise TypeError("loaded must be KokunoPublicTerminalMultiplierTarget")
    if not isinstance(pre_serialization_reference, KokunoPublicTerminalMultiplierTarget):
        raise TypeError("pre_serialization_reference must be KokunoPublicTerminalMultiplierTarget")

    loaded_semantic = _semantic_value(loaded)
    reference_semantic = _semantic_value(pre_serialization_reference)
    save_load_exact = loaded_semantic == reference_semantic
    if not save_load_exact:
        raise ValueError("terminal target semantic identity changed across save/load")

    cfg = loaded.configuration()
    if cfg != pre_serialization_reference.configuration():
        raise ValueError("terminal target configuration changed across save/load")
    if cfg.get("truth_boundary", {}).get("pde_validated") is not False:
        raise ValueError("upstream target illegally promotes pde_validated")
    if cfg.get("truth_boundary", {}).get("heldout_ns_residual_assessed") is not False:
        raise ValueError("upstream target illegally promotes held-out NS residual")

    h = float(cfg["h_current"])
    c_o = float(cfg["c_o_autonomous"])
    rho = float(cfg["rho_o"])
    if not math.isclose(rho, c_o * h, rel_tol=0.0, abs_tol=0.0):
        raise ValueError("serialized rho_o != c_o*h")
    if not (0.0 < h < 0.01 and 0.0 < c_o < 1.0 and 0.0 < rho < 1.0):
        raise ValueError("serialized terminal target parameters are outside frozen bounds")

    # Form all independent scientific values BEFORE asking A1 for its production receipt.
    q_by_resolution = {
        str(n): _q_p_ibp_reference(h, rho, n) for n in SIMPSON_SUBINTERVALS
    }
    q_coarse = q_by_resolution[str(SIMPSON_SUBINTERVALS[0])]
    q_medium = q_by_resolution[str(SIMPSON_SUBINTERVALS[1])]
    q_fine = q_by_resolution[str(SIMPSON_SUBINTERVALS[2])]
    slope = _slope_audit(h, rho)

    # Production values are comparison targets only; they are not used to build the reference.
    public = loaded.target_report()
    q_public = float(public["q_p_reduced"])
    q_public_direct = float(public["q_p_direct"])
    q_public_internal_rel = float(public["q_p_relative_replay_error"])

    report = {
        "schema": "kokuno-a4-terminal-target-independent-audit-v1",
        "task": TASK,
        "upstream_head": UPSTREAM_HEAD,
        "candidate_semantic_sha256": loaded_semantic,
        "save_load_semantic_exact": save_load_exact,
        "independent_operator": {
            "flat_step": "A4-owned stable public sigma reimplementation",
            "quadrature": "composite Simpson",
            "identity": (
                "Q_p=[exp((1-h)v)f(v)]_0^3/f(0)"
                "-(1-h)/f(0)*int_0^3 exp((1-h)v)f(v)dv"
            ),
            "subintervals": list(SIMPSON_SUBINTERVALS),
            "production_paths_excluded": [
                "A1 _sigma",
                "A1 _sigma_prime",
                "A1 q_p_direct",
                "A1 q_p_reduced",
                "A1 GL192 nodes/weights",
            ],
        },
        "parameters": {"h": h, "c_o": c_o, "rho_o": rho},
        "q_p": {
            "independent_by_resolution": q_by_resolution,
            "coarse_to_medium_relative_change": _relative_difference(q_coarse, q_medium),
            "medium_to_fine_relative_change": _relative_difference(q_medium, q_fine),
            "public_q_p_reduced": q_public,
            "public_q_p_direct": q_public_direct,
            "fine_independent_vs_public_relative_difference": _relative_difference(
                q_fine, q_public
            ),
            "public_direct_vs_reduced_relative_difference": _relative_difference(
                q_public_direct, q_public
            ),
            "public_reported_direct_reduced_relative_error": q_public_internal_rel,
        },
        "slope": slope,
        "truth_boundary": {
            "terminal_target_independent_algebra_audited": True,
            "cartesian_terminal_velocity_materialized": False,
            "leading_only_ns_residual_assessed": False,
            "leading_plus_oscillatory_ns_residual_assessed": False,
            "after_correction_ns_residual_assessed": False,
            "matched_pressure_assessed": False,
            "restricted_forcing_assessed": False,
            "complete_ns_residual_assessed": False,
            "pde_validated": False,
        },
        "gates": {
            "fine_public_q_rel": Q_FINE_PUBLIC_REL_GATE,
            "medium_fine_q_rel": Q_MEDIUM_FINE_REL_GATE,
            "public_direct_reduced_rel": PUBLIC_DIRECT_REDUCED_REL_GATE,
            "slope_negative_floor": SLOPE_NEGATIVE_FLOOR,
            "final_project_momentum_gate_unchanged": FINAL_MOMENTUM_GATE,
            "final_project_divergence_gate_unchanged": FINAL_DIVERGENCE_GATE,
            "final_project_quadrature_unchanged": list(FINAL_QUADRATURE),
        },
        "final_project_admission_ready": False,
    }
    return report


def enforce_preregistered_gates(report: Mapping[str, Any]) -> None:
    q = report["q_p"]
    slope = report["slope"]
    params = report["parameters"]
    truth = report["truth_boundary"]

    assert float(q["fine_independent_vs_public_relative_difference"]) <= Q_FINE_PUBLIC_REL_GATE
    assert float(q["medium_to_fine_relative_change"]) <= Q_MEDIUM_FINE_REL_GATE
    assert float(q["public_direct_vs_reduced_relative_difference"]) <= PUBLIC_DIRECT_REDUCED_REL_GATE
    assert float(q["public_reported_direct_reduced_relative_error"]) <= PUBLIC_DIRECT_REDUCED_REL_GATE
    assert float(slope["minimum"]) >= -SLOPE_NEGATIVE_FLOOR
    assert float(slope["maximum"]) < float(slope["h_over_4"])
    q_fine = float(q["independent_by_resolution"][str(SIMPSON_SUBINTERVALS[-1])])
    assert 0.0 < q_fine < float(params["h"])
    assert report["save_load_semantic_exact"] is True
    assert truth["complete_ns_residual_assessed"] is False
    assert truth["pde_validated"] is False
    assert report["final_project_admission_ready"] is False
    gates = report["gates"]
    assert float(gates["final_project_momentum_gate_unchanged"]) == FINAL_MOMENTUM_GATE
    assert float(gates["final_project_divergence_gate_unchanged"]) == FINAL_DIVERGENCE_GATE
    assert tuple(gates["final_project_quadrature_unchanged"]) == FINAL_QUADRATURE


def public_api_has_no_scientific_tuning_knobs() -> bool:
    params = inspect.signature(audit_loaded_terminal_target).parameters
    lowered = {name.lower() for name in params}
    return not any(name in _FORBIDDEN_PUBLIC_KNOBS for name in lowered)


def default_target() -> KokunoPublicTerminalMultiplierTarget:
    return KokunoPublicTerminalMultiplierTarget()
