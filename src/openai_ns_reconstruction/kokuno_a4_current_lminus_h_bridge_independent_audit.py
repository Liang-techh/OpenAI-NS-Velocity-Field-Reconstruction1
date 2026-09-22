"""Independent A4 audit of the current inhomogeneous ``l=-h`` Q_s bridge.

K4-VAL-129 is a bounded scalar-transport audit, not Navier--Stokes
admission.  Upstream A1 #1233 materializes the current eta-resolved bridge

    Q_s' + (1-h) Q_s = h (W-1),

with the carried current M/M_eta state.  On this stage

    P_0(eta) = 1 - W(0,eta),
    W(y,eta)-1 = -P_0(eta) exp(-y),

so an equivalent public ODE is

    Q_s' = -(1-h) Q_s - h P_0 exp(-y).

A4 deliberately does not use A1's analytic ``D_y_Q_s`` or its pointwise
bisection path to construct its scientific reference.  After a real
save/load, this module:

* integrates the public bridge ODE independently with classical RK4 on
  three frozen resolutions;
* differentiates public ``Q_s(y,eta)`` independently with a centered
  five-point FD4 stencil on three frozen y steps;
* reports transport disagreement, the independently reconstructed ODE
  defect, refinement behavior, carried-state nontriviality and pointwise
  target-time spread.

Production analytic derivatives and production pointwise target times are
read only after the independent references have been formed, as diagnostics.
Passing this scoped audit cannot establish one full-eta scalar bridge length,
an eta-dependent Cartesian terminal boundary, a Cartesian velocity extension,
or any complete Navier--Stokes residual.
"""

from __future__ import annotations

import inspect
import math
from typing import Any, Mapping

import numpy as np

from .kokuno_current_exterior_lminus_h_bridge import (
    KokunoCurrentExteriorLMinusHBridge,
)

TASK = "K4-VAL-129"
UPSTREAM_HEAD = "6ae44fdac17483b1fd042dca6f5bdf06701ae068"
SEED = 9174001

RK4_STEPS_PER_UNIT = (64, 128, 256)
FD4_Y_STEPS = (2.0**-7, 2.0**-8, 2.0**-9)

TRANSPORT_NORM_MAX_GATE = 2.0e-8
TRANSPORT_NORM_RMS_GATE = 2.0e-8
ODE_DEFECT_NORM_MAX_GATE = 2.0e-8
ODE_DEFECT_NORM_RMS_GATE = 2.0e-8
REFINEMENT_WORSEN_FACTOR = 1.25
REFINEMENT_FLOOR = 5.0e-10
INHOMOGENEOUS_P_MIN = 1.0e-12

FINAL_MOMENTUM_GATE = 1.0e-3
FINAL_DIVERGENCE_GATE = 1.0e-5
FINAL_QUADRATURE = (24, 48, 96)

_TARGET_SEARCH_MAX_Y = 32.0

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
    "h",
    "lambda_value",
    "target",
    "bridge_length",
}


def _semantic_value(obj: Any) -> str:
    value = getattr(obj, "semantic_sha256")
    return str(value() if callable(value) else value)


def default_candidate() -> KokunoCurrentExteriorLMinusHBridge:
    return KokunoCurrentExteriorLMinusHBridge()


def heldout_eta(candidate: KokunoCurrentExteriorLMinusHBridge) -> np.ndarray:
    """Return 33 deterministic off-grid eta probes on the public bridge chart."""
    lo, hi = (float(v) for v in candidate.eta_interval)
    if not (math.isfinite(lo) and math.isfinite(hi) and lo < hi):
        raise ValueError("candidate eta interval must be finite and nonempty")
    center = 0.5 * (lo + hi)
    half = 0.5 * (hi - lo)
    anchors = np.asarray(
        [-0.91, -0.74, -0.55, -0.31, 0.0, 0.23, 0.46, 0.71, 0.91],
        dtype=float,
    )
    rng = np.random.default_rng(SEED)
    random_fraction = rng.uniform(-0.93, 0.93, size=24)
    eta = center + half * np.concatenate([anchors, random_fraction])
    if np.any((eta <= lo) | (eta >= hi)):
        raise RuntimeError("held-out eta generation left the strict interior")
    return np.asarray(eta, dtype=float)


def heldout_y() -> np.ndarray:
    """Return frozen bridge-y probes, separate from A1's root-search path."""
    anchors = np.asarray(
        [0.125, 0.375, 0.9, 1.8, 3.4, 5.6, 8.2, 10.4, 12.0],
        dtype=float,
    )
    rng = np.random.default_rng(SEED + 1)
    random_values = rng.uniform(0.14, 11.86, size=15)
    y = np.sort(np.concatenate([anchors, random_values]))
    if np.any(y <= 2.0 * max(FD4_Y_STEPS)):
        raise RuntimeError("held-out y values are unsafe for frozen FD4 stencil")
    return np.asarray(y, dtype=float)


def _entry_channels(
    candidate: KokunoCurrentExteriorLMinusHBridge,
    eta: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read public entry primitives and independently rebuild P_0."""
    ee = np.asarray(eta, dtype=float)
    entry = candidate.entry_state(ee)
    q0 = np.asarray(entry["Q_s_entry"], dtype=float)
    m0 = np.asarray(entry["M_over_X_entry"], dtype=float)
    me0 = np.asarray(entry["M_eta_over_X_entry"], dtype=float)
    w0 = np.asarray(entry["W_entry"], dtype=float)
    d = 1.0 - ee * ee
    p0 = 2.0 * float(candidate.D) * ee * m0 + d * me0
    arrays = (q0, p0, w0)
    if any(np.any(~np.isfinite(a)) for a in arrays):
        raise RuntimeError("bridge entry channels became non-finite")
    return q0, p0, w0


def _rhs(
    q: np.ndarray,
    p0: np.ndarray,
    y: float,
    h_value: float,
) -> np.ndarray:
    return -(1.0 - h_value) * q - h_value * p0 * math.exp(-float(y))


def _rk4_step(
    q: np.ndarray,
    p0: np.ndarray,
    y: float,
    dt: float,
    h_value: float,
) -> np.ndarray:
    k1 = _rhs(q, p0, y, h_value)
    k2 = _rhs(q + 0.5 * dt * k1, p0, y + 0.5 * dt, h_value)
    k3 = _rhs(q + 0.5 * dt * k2, p0, y + 0.5 * dt, h_value)
    k4 = _rhs(q + dt * k3, p0, y + dt, h_value)
    out = q + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
    if np.any(~np.isfinite(out)):
        raise RuntimeError("independent RK4 bridge integration became non-finite")
    return out


def _rk4_at_targets(
    q0: np.ndarray,
    p0: np.ndarray,
    targets: np.ndarray,
    steps_per_unit: int,
    h_value: float,
) -> np.ndarray:
    """Integrate all eta channels simultaneously to sorted y targets."""
    yy = np.asarray(targets, dtype=float)
    if yy.ndim != 1 or np.any(~np.isfinite(yy)) or np.any(yy < 0.0):
        raise ValueError("targets must be a finite nonnegative 1-D array")
    order = np.argsort(yy)
    sorted_y = yy[order]
    max_dt = 1.0 / int(steps_per_unit)
    if not (math.isfinite(max_dt) and max_dt > 0.0):
        raise ValueError("steps_per_unit must be positive")

    q = np.asarray(q0, dtype=float).copy()
    p = np.asarray(p0, dtype=float)
    if q.shape != p.shape:
        raise ValueError("q0 and p0 must have identical shape")
    out_sorted = np.empty((q.size, sorted_y.size), dtype=float)
    y = 0.0
    for j, target in enumerate(sorted_y):
        target_f = float(target)
        while y + 0.5 * max_dt < target_f:
            dt = min(max_dt, target_f - y)
            q = _rk4_step(q, p, y, dt, h_value)
            y += dt
        if y < target_f:
            dt = target_f - y
            q = _rk4_step(q, p, y, dt, h_value)
            y = target_f
        out_sorted[:, j] = q

    inverse = np.empty_like(order)
    inverse[order] = np.arange(order.size)
    return out_sorted[:, inverse]


def _fd4_public_q(
    candidate: KokunoCurrentExteriorLMinusHBridge,
    y: np.ndarray,
    eta: np.ndarray,
    step: float,
) -> np.ndarray:
    """Differentiate only public Q_s values with an A4-owned five-point FD4."""
    yy = np.asarray(y, dtype=float)
    ee = np.asarray(eta, dtype=float)
    h = float(step)
    if not (math.isfinite(h) and h > 0.0):
        raise ValueError("FD4 y step must be positive and finite")
    if np.any(yy - 2.0 * h < 0.0):
        raise ValueError("FD4 stencil left nonnegative bridge-y domain")

    def q_at(offset: float) -> np.ndarray:
        state = candidate.state(yy + offset, ee)
        value = np.asarray(state["Q_s"], dtype=float)
        if np.any(~np.isfinite(value)):
            raise RuntimeError("public Q_s became non-finite")
        return value

    fm2 = q_at(-2.0 * h)
    fm1 = q_at(-h)
    fp1 = q_at(h)
    fp2 = q_at(2.0 * h)
    derivative = (fm2 - 8.0 * fm1 + 8.0 * fp1 - fp2) / (12.0 * h)
    if np.any(~np.isfinite(derivative)):
        raise RuntimeError("independent FD4 derivative became non-finite")
    return derivative


def _normalized_error(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    aa = np.asarray(a, dtype=float)
    bb = np.asarray(b, dtype=float)
    scale = np.maximum.reduce([np.ones_like(aa), np.abs(aa), np.abs(bb)])
    return np.abs(aa - bb) / scale


def _normalized_defect(
    dq: np.ndarray,
    q: np.ndarray,
    p: np.ndarray,
    h_value: float,
) -> np.ndarray:
    term_q = (1.0 - h_value) * np.asarray(q, dtype=float)
    term_p = h_value * np.asarray(p, dtype=float)
    defect = np.asarray(dq, dtype=float) + term_q + term_p
    scale = np.maximum.reduce(
        [
            np.ones_like(defect),
            np.abs(dq),
            np.abs(term_q),
            np.abs(term_p),
        ]
    )
    return np.abs(defect) / scale


def _resolution_report(
    candidate: KokunoCurrentExteriorLMinusHBridge,
    eta: np.ndarray,
    y: np.ndarray,
    steps_per_unit: int,
    fd_step: float,
) -> dict[str, Any]:
    """Build independent RK4 and FD4 references before reading A1 derivative jets."""
    q0, p0, _ = _entry_channels(candidate, eta)
    q_independent = _rk4_at_targets(
        q0,
        p0,
        y,
        steps_per_unit,
        float(candidate.h_value),
    )

    eta_grid, y_grid = np.meshgrid(eta, y, indexing="ij")
    public_state = candidate.state(y_grid, eta_grid)
    q_public = np.asarray(public_state["Q_s"], dtype=float)

    transport_norm = _normalized_error(q_independent, q_public)

    dq_fd = _fd4_public_q(candidate, y_grid, eta_grid, fd_step)
    p_public_formula = p0[:, None] * np.exp(-y_grid)
    ode_norm = _normalized_defect(
        dq_fd,
        q_public,
        p_public_formula,
        float(candidate.h_value),
    )

    # Production analytic derivative is intentionally consumed only after the
    # two independent references above have already been formed.
    dq_production = np.asarray(public_state["D_y_Q_s"], dtype=float)
    derivative_diag = _normalized_error(dq_fd, dq_production)

    worst_transport = np.unravel_index(
        int(np.argmax(transport_norm)), transport_norm.shape
    )
    worst_ode = np.unravel_index(int(np.argmax(ode_norm)), ode_norm.shape)

    return {
        "rk4_steps_per_unit": int(steps_per_unit),
        "fd4_y_step": float(fd_step),
        "transport_normalized_max": float(np.max(transport_norm)),
        "transport_normalized_rms": float(
            np.sqrt(np.mean(transport_norm * transport_norm))
        ),
        "transport_absolute_max": float(np.max(np.abs(q_independent - q_public))),
        "ode_defect_normalized_max": float(np.max(ode_norm)),
        "ode_defect_normalized_rms": float(np.sqrt(np.mean(ode_norm * ode_norm))),
        "analytic_derivative_diagnostic_normalized_max": float(
            np.max(derivative_diag)
        ),
        "analytic_derivative_diagnostic_normalized_rms": float(
            np.sqrt(np.mean(derivative_diag * derivative_diag))
        ),
        "worst_transport_witness": {
            "eta": float(eta[worst_transport[0]]),
            "y": float(y[worst_transport[1]]),
            "q_public": float(q_public[worst_transport]),
            "q_independent_rk4": float(q_independent[worst_transport]),
            "normalized_error": float(transport_norm[worst_transport]),
        },
        "worst_ode_witness": {
            "eta": float(eta[worst_ode[0]]),
            "y": float(y[worst_ode[1]]),
            "q_public": float(q_public[worst_ode]),
            "dq_fd4": float(dq_fd[worst_ode]),
            "normalized_defect": float(ode_norm[worst_ode]),
        },
    }


def _independent_target_times(
    candidate: KokunoCurrentExteriorLMinusHBridge,
    eta: np.ndarray,
) -> np.ndarray:
    """Find pointwise Q_p first hits by independent fixed-step RK4 only."""
    q, p0, _ = _entry_channels(candidate, eta)
    q_p = float(candidate.q_p)
    if np.any(q <= q_p):
        raise RuntimeError("bridge entry must lie strictly above Q_p")
    dt = 1.0 / RK4_STEPS_PER_UNIT[-1]
    y = 0.0
    times = np.full(q.shape, np.nan, dtype=float)
    active = np.ones(q.shape, dtype=bool)

    while y < _TARGET_SEARCH_MAX_Y and np.any(active):
        old = q.copy()
        q = _rk4_step(q, p0, y, dt, float(candidate.h_value))
        new_y = y + dt
        crossed = active & (q <= q_p)
        if np.any(crossed):
            denom = old[crossed] - q[crossed]
            frac = np.where(
                denom != 0.0,
                (old[crossed] - q_p) / denom,
                1.0,
            )
            frac = np.clip(frac, 0.0, 1.0)
            times[crossed] = y + dt * frac
            active[crossed] = False
        y = new_y

    if np.any(active) or np.any(~np.isfinite(times)):
        raise RuntimeError("independent RK4 failed to reach Q_p on all eta probes")
    return times


def audit_loaded_current_bridge(
    loaded: KokunoCurrentExteriorLMinusHBridge,
    pre_serialization_reference: KokunoCurrentExteriorLMinusHBridge,
) -> dict[str, Any]:
    """Run K4-VAL-129 using only preregistered internal choices."""
    if not isinstance(loaded, KokunoCurrentExteriorLMinusHBridge):
        raise TypeError("loaded must be KokunoCurrentExteriorLMinusHBridge")
    if not isinstance(pre_serialization_reference, KokunoCurrentExteriorLMinusHBridge):
        raise TypeError("pre_serialization_reference must match the candidate type")

    eta = heldout_eta(loaded)
    y = heldout_y()

    before = pre_serialization_reference.state(
        np.asarray([0.125, 3.4, 10.4])[:, None],
        eta[None, :],
    )
    after = loaded.state(
        np.asarray([0.125, 3.4, 10.4])[:, None],
        eta[None, :],
    )
    semantic_exact = _semantic_value(loaded) == _semantic_value(pre_serialization_reference)
    configuration_exact = loaded.configuration() == pre_serialization_reference.configuration()
    q_replay_exact = np.array_equal(
        np.asarray(before["Q_s"], dtype=float),
        np.asarray(after["Q_s"], dtype=float),
    )

    resolutions = [
        _resolution_report(loaded, eta, y, rk4, fd)
        for rk4, fd in zip(RK4_STEPS_PER_UNIT, FD4_Y_STEPS)
    ]

    _, p0, w0 = _entry_channels(loaded, eta)
    p_from_w = 1.0 - w0
    p_identity_norm = _normalized_error(p0, p_from_w)

    target_independent = _independent_target_times(loaded, eta)
    # Production pointwise bisection is read only after independent target
    # times are fully formed.
    target_production = np.asarray(loaded.pointwise_target_time(eta), dtype=float)
    target_delta = target_independent - target_production
    target_spread = float(np.max(target_independent) - np.min(target_independent))

    return {
        "schema": "kokuno-a4-current-lminus-h-bridge-independent-audit-v1",
        "task": TASK,
        "upstream_head": UPSTREAM_HEAD,
        "seed": SEED,
        "candidate_semantic_sha256": _semantic_value(loaded),
        "heldout_eta": eta.tolist(),
        "heldout_y": y.tolist(),
        "rk4_steps_per_unit": list(RK4_STEPS_PER_UNIT),
        "fd4_y_steps": list(FD4_Y_STEPS),
        "resolutions": resolutions,
        "save_load": {
            "semantic_exact": bool(semantic_exact),
            "configuration_exact": bool(configuration_exact),
            "q_state_exact": bool(q_replay_exact),
        },
        "inhomogeneous_entry": {
            "max_abs_P_entry": float(np.max(np.abs(p0))),
            "rms_P_entry": float(np.sqrt(np.mean(p0 * p0))),
            "P_equals_one_minus_W_normalized_max": float(np.max(p_identity_norm)),
            "P_equals_one_minus_W_normalized_rms": float(
                np.sqrt(np.mean(p_identity_norm * p_identity_norm))
            ),
        },
        "target_time_diagnostic": {
            "independent_rk4_target_time_min": float(np.min(target_independent)),
            "independent_rk4_target_time_max": float(np.max(target_independent)),
            "independent_rk4_target_time_spread": target_spread,
            "independent_vs_production_target_time_max_abs": float(
                np.max(np.abs(target_delta))
            ),
            "independent_vs_production_target_time_rms": float(
                np.sqrt(np.mean(target_delta * target_delta))
            ),
            "diagnostic_only_not_a_scalar_bridge_construction": True,
        },
        "gates": {
            "fine_transport_normalized_max": TRANSPORT_NORM_MAX_GATE,
            "fine_transport_normalized_rms": TRANSPORT_NORM_RMS_GATE,
            "fine_ode_defect_normalized_max": ODE_DEFECT_NORM_MAX_GATE,
            "fine_ode_defect_normalized_rms": ODE_DEFECT_NORM_RMS_GATE,
            "medium_to_fine_worsen_factor": REFINEMENT_WORSEN_FACTOR,
            "refinement_floor": REFINEMENT_FLOOR,
            "inhomogeneous_P_min": INHOMOGENEOUS_P_MIN,
            "final_project_momentum_gate_unchanged": FINAL_MOMENTUM_GATE,
            "final_project_divergence_gate_unchanged": FINAL_DIVERGENCE_GATE,
            "final_project_quadrature_unchanged": list(FINAL_QUADRATURE),
        },
        "truth_boundary": {
            "current_l_minus_h_bridge_independent_audit_executed": True,
            "full_eta_interval_common_scalar_bridge_length_established": False,
            "current_l_minus_h_matching_bridge_materialized": False,
            "eta_dependent_cartesian_matching_boundary_materialized": False,
            "cartesian_terminal_multiplier_composed": False,
            "leading_only_ns_residual_assessed": False,
            "leading_plus_oscillatory_ns_residual_assessed": False,
            "after_correction_ns_residual_assessed": False,
            "complete_ns_admission_ready": False,
            "pde_validated": False,
        },
        "final_project_admission_ready": False,
    }


def _check_refinement(
    resolutions: list[Mapping[str, Any]],
    key: str,
    label: str,
) -> None:
    medium = float(resolutions[-2][key])
    fine = float(resolutions[-1][key])
    if fine > max(REFINEMENT_FLOOR, REFINEMENT_WORSEN_FACTOR * medium):
        raise AssertionError(
            f"medium-to-fine {label} worsened beyond preregistered allowance"
        )


def enforce_preregistered_gates(report: Mapping[str, Any]) -> None:
    """Fail closed on the exact K4-VAL-129 gates preregistered in Issue #15."""
    resolutions = list(report["resolutions"])
    if len(resolutions) != 3:
        raise AssertionError("K4-VAL-129 requires exactly three resolutions")
    fine = resolutions[-1]

    if fine["transport_normalized_max"] > TRANSPORT_NORM_MAX_GATE:
        raise AssertionError("fine independent bridge transport normalized max exceeds gate")
    if fine["transport_normalized_rms"] > TRANSPORT_NORM_RMS_GATE:
        raise AssertionError("fine independent bridge transport normalized RMS exceeds gate")
    if fine["ode_defect_normalized_max"] > ODE_DEFECT_NORM_MAX_GATE:
        raise AssertionError("fine independent bridge ODE-defect normalized max exceeds gate")
    if fine["ode_defect_normalized_rms"] > ODE_DEFECT_NORM_RMS_GATE:
        raise AssertionError("fine independent bridge ODE-defect normalized RMS exceeds gate")

    _check_refinement(resolutions, "transport_normalized_max", "transport metric")
    _check_refinement(resolutions, "ode_defect_normalized_max", "ODE-defect metric")

    replay = report["save_load"]
    if not (
        replay["semantic_exact"]
        and replay["configuration_exact"]
        and replay["q_state_exact"]
    ):
        raise AssertionError("save/load replay is not exact")

    if report["inhomogeneous_entry"]["max_abs_P_entry"] < INHOMOGENEOUS_P_MIN:
        raise AssertionError("carried inhomogeneous bridge term fell below preregistered floor")

    truth = report["truth_boundary"]
    if truth["pde_validated"] or truth["complete_ns_admission_ready"]:
        raise AssertionError("scoped scalar bridge audit attempted a PDE truth promotion")
    if truth["full_eta_interval_common_scalar_bridge_length_established"]:
        raise AssertionError("finite eta probes cannot establish a full-interval scalar bridge")
    if truth["current_l_minus_h_matching_bridge_materialized"]:
        raise AssertionError("A4 audit cannot materialize the A1 matching bridge")
    if truth["eta_dependent_cartesian_matching_boundary_materialized"]:
        raise AssertionError("A4 audit cannot materialize an eta-dependent Cartesian boundary")
    if report["final_project_admission_ready"]:
        raise AssertionError("scoped scalar bridge audit cannot be final project admission")


def public_api_has_no_scientific_tuning_knobs() -> bool:
    sig = inspect.signature(audit_loaded_current_bridge)
    return _FORBIDDEN_PUBLIC_KNOBS.isdisjoint(sig.parameters)
