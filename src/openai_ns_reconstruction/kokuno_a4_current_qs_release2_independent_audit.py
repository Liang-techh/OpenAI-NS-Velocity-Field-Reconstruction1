"""Independent A4 audit of the current-candidate release2-end ``Q_s`` state.

K4-VAL-128 is a bounded scalar-state audit, not Navier--Stokes admission.
Upstream A1 #1225 materializes the current release2-end radial state using the
full general identity

    W = 1 - 2 D eta M/X - d M_eta/X,
    Q_s = -W + (1-h) I/(XH) - D eta I_eta/(XH)
          - d J_eta/(XH) + 2(h-D) eta J/(XH).

The production implementation carries analytic/current eta jets.  A4 does not
reuse those jets as its numerical reference.  After a real save/load, this
module differentiates only the public *base* channels ``M/X``, ``I/(XH)`` and
``J/(XH)`` with an A4-owned centered FD4 stencil on three frozen eta steps,
then rebuilds ``Q_s`` independently.  The production derivative channels are
read only afterwards as comparison targets.

The same receipt also reports, without gating or constructing it, the spread
of the implied ``l=-h`` bridge length to the already-fixed public ``Q_p``.
A nonzero spread would mean one scalar bridge length cannot exactly match all
eta probes and must be resolved by A1 before Cartesian continuation.

Passing this scoped audit cannot promote a matching bridge, velocity, pressure,
forcing, complete residual, or ``pde_validated``.
"""

from __future__ import annotations

import inspect
import math
from typing import Any, Callable, Mapping

import numpy as np

from .kokuno_current_exterior_qs_release2_state import (
    KokunoCurrentExteriorQsRelease2State,
)

TASK = "K4-VAL-128"
UPSTREAM_HEAD = "8503f2ede2dd604392198b31ffdeb79321822a30"
SEED = 9173991
ETA_FD_STEPS = (2.0**-8, 2.0**-9, 2.0**-10)
Q_RECON_NORM_MAX_GATE = 1.0e-7
Q_RECON_NORM_RMS_GATE = 1.0e-7
DERIVATIVE_CHANNEL_NORM_MAX_GATE = 2.0e-6
REFINEMENT_WORSEN_FACTOR = 1.25
REFINEMENT_FLOOR = 2.0e-9
CURRENT_MINUS_IDEAL_MIN = 1.0e-12
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
    "h",
    "lambda_value",
}


def _semantic_value(obj: Any) -> str:
    value = getattr(obj, "semantic_sha256")
    return str(value() if callable(value) else value)


def default_candidate() -> KokunoCurrentExteriorQsRelease2State:
    return KokunoCurrentExteriorQsRelease2State()


def heldout_eta(candidate: KokunoCurrentExteriorQsRelease2State) -> np.ndarray:
    """Deterministic off-grid eta probes with margin for the coarsest FD4 stencil."""
    lo, hi = (float(v) for v in candidate.eta_interval)
    if not (math.isfinite(lo) and math.isfinite(hi) and lo < hi):
        raise ValueError("candidate eta interval must be finite and nonempty")
    center = 0.5 * (lo + hi)
    half = 0.5 * (hi - lo)
    # Keep ten percent of the chart on both sides for +/-2h at the frozen ladder.
    if 0.10 * half <= 2.0 * max(ETA_FD_STEPS):
        raise ValueError("eta chart is too narrow for preregistered FD4 ladder")
    anchors = np.asarray([-0.90, -0.73, -0.51, -0.29, 0.0, 0.24, 0.47, 0.69, 0.90])
    rng = np.random.default_rng(SEED)
    random_fraction = rng.uniform(-0.86, 0.86, size=24)
    fraction = np.concatenate([anchors, random_fraction])
    eta = center + half * fraction
    return np.asarray(eta, dtype=float)


def _fd4(
    evaluator: Callable[[np.ndarray], np.ndarray],
    x: np.ndarray,
    step: float,
) -> np.ndarray:
    """A4-owned centered five-point first derivative."""
    h = float(step)
    if not (math.isfinite(h) and h > 0.0):
        raise ValueError("FD4 step must be positive and finite")
    values = np.asarray(x, dtype=float)
    fm2 = np.asarray(evaluator(values - 2.0 * h), dtype=float)
    fm1 = np.asarray(evaluator(values - h), dtype=float)
    fp1 = np.asarray(evaluator(values + h), dtype=float)
    fp2 = np.asarray(evaluator(values + 2.0 * h), dtype=float)
    out = (fm2 - 8.0 * fm1 + 8.0 * fp1 - fp2) / (12.0 * h)
    if np.any(~np.isfinite(out)):
        raise RuntimeError("independent FD4 derivative became non-finite")
    return out


def _channel(candidate: KokunoCurrentExteriorQsRelease2State, eta: np.ndarray, key: str) -> np.ndarray:
    state = candidate.state(np.asarray(eta, dtype=float))
    value = np.asarray(state[key], dtype=float)
    if np.any(~np.isfinite(value)):
        raise RuntimeError(f"public base channel {key} became non-finite")
    return value


def _normalized_error(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    aa = np.asarray(a, dtype=float)
    bb = np.asarray(b, dtype=float)
    scale = np.maximum.reduce([np.ones_like(aa), np.abs(aa), np.abs(bb)])
    return np.abs(aa - bb) / scale


def _resolution_report(
    candidate: KokunoCurrentExteriorQsRelease2State,
    eta: np.ndarray,
    step: float,
) -> dict[str, Any]:
    center = candidate.state(eta)
    M = np.asarray(center["M_over_X"], dtype=float)
    r_I = np.asarray(center["I_over_XH"], dtype=float)
    J = np.asarray(center["J_over_XH"], dtype=float)
    q_public = np.asarray(center["Q_s_current_release2_end"], dtype=float)

    M_eta_fd = _fd4(lambda e: _channel(candidate, e, "M_over_X"), eta, step)
    I_eta_fd = _fd4(lambda e: _channel(candidate, e, "I_over_XH"), eta, step)
    J_eta_fd = _fd4(lambda e: _channel(candidate, e, "J_over_XH"), eta, step)

    d = 1.0 - eta * eta
    W_fd = 1.0 - 2.0 * candidate.D * eta * M - d * M_eta_fd
    q_fd = (
        -W_fd
        + (1.0 - candidate.h_value) * r_I
        - candidate.D * eta * I_eta_fd
        - d * J_eta_fd
        + 2.0 * (candidate.h_value - candidate.D) * eta * J
    )
    if np.any(~np.isfinite(q_fd)):
        raise RuntimeError("independent Q_s reconstruction became non-finite")

    q_norm = _normalized_error(q_fd, q_public)
    deriv_pairs = {
        "M_eta_over_X": (M_eta_fd, np.asarray(center["M_eta_over_X"], dtype=float)),
        "I_eta_over_XH": (I_eta_fd, np.asarray(center["I_eta_over_XH"], dtype=float)),
        "J_eta_over_XH": (J_eta_fd, np.asarray(center["J_eta_over_XH"], dtype=float)),
    }
    deriv_report: dict[str, dict[str, float]] = {}
    all_deriv_norm = []
    for name, (independent, production) in deriv_pairs.items():
        norm = _normalized_error(independent, production)
        all_deriv_norm.append(norm)
        deriv_report[name] = {
            "normalized_max": float(np.max(norm)),
            "normalized_rms": float(np.sqrt(np.mean(norm * norm))),
            "absolute_max": float(np.max(np.abs(independent - production))),
        }
    stacked = np.concatenate(all_deriv_norm)
    worst = int(np.argmax(q_norm))
    return {
        "eta_fd_step": float(step),
        "q_reconstruction_normalized_max": float(np.max(q_norm)),
        "q_reconstruction_normalized_rms": float(np.sqrt(np.mean(q_norm * q_norm))),
        "q_reconstruction_absolute_max": float(np.max(np.abs(q_fd - q_public))),
        "derivative_channel_normalized_max": float(np.max(stacked)),
        "derivative_channel_normalized_rms": float(np.sqrt(np.mean(stacked * stacked))),
        "derivative_channels": deriv_report,
        "worst_q_reconstruction_witness": {
            "eta": float(eta[worst]),
            "q_public": float(q_public[worst]),
            "q_independent": float(q_fd[worst]),
            "normalized_error": float(q_norm[worst]),
        },
        "q_independent": q_fd.tolist(),
    }


def _bridge_diagnostic(
    candidate: KokunoCurrentExteriorQsRelease2State,
    q_independent: np.ndarray,
) -> dict[str, Any]:
    q = np.asarray(q_independent, dtype=float)
    q_p = float(candidate.parent.q_p)
    rate = 1.0 - float(candidate.h_value)
    valid = np.isfinite(q) & (q > q_p) & (q_p > 0.0) & (rate > 0.0)
    lengths = np.full_like(q, np.nan, dtype=float)
    lengths[valid] = np.log(q[valid] / q_p) / rate
    finite = lengths[np.isfinite(lengths)]
    if finite.size:
        minimum = float(np.min(finite))
        maximum = float(np.max(finite))
        spread = maximum - minimum
        rms_about_mean = float(np.sqrt(np.mean((finite - np.mean(finite)) ** 2)))
    else:
        minimum = maximum = spread = rms_about_mean = math.nan
    return {
        "q_p": q_p,
        "eligible_probe_count": int(np.sum(valid)),
        "ineligible_probe_count": int(q.size - np.sum(valid)),
        "implied_bridge_lengths": [None if not math.isfinite(v) else float(v) for v in lengths],
        "implied_bridge_length_min": minimum,
        "implied_bridge_length_max": maximum,
        "implied_bridge_length_spread": spread,
        "implied_bridge_length_rms_about_mean": rms_about_mean,
        "diagnostic_only_not_a_bridge_construction": True,
    }


def audit_loaded_current_qs_state(
    loaded: KokunoCurrentExteriorQsRelease2State,
    pre_serialization_reference: KokunoCurrentExteriorQsRelease2State,
) -> dict[str, Any]:
    """Run K4-VAL-128 using only frozen internal protocol choices."""
    if not isinstance(loaded, KokunoCurrentExteriorQsRelease2State):
        raise TypeError("loaded must be KokunoCurrentExteriorQsRelease2State")
    if not isinstance(pre_serialization_reference, KokunoCurrentExteriorQsRelease2State):
        raise TypeError("pre_serialization_reference must match the candidate type")

    eta = heldout_eta(loaded)
    before = pre_serialization_reference.state(eta)
    after = loaded.state(eta)
    semantic_exact = _semantic_value(loaded) == _semantic_value(pre_serialization_reference)
    configuration_exact = loaded.configuration() == pre_serialization_reference.configuration()
    q_replay_exact = np.array_equal(
        np.asarray(before["Q_s_current_release2_end"], dtype=float),
        np.asarray(after["Q_s_current_release2_end"], dtype=float),
    )

    resolutions = [
        _resolution_report(loaded, eta, step)
        for step in ETA_FD_STEPS
    ]
    fine_q = np.asarray(resolutions[-1]["q_independent"], dtype=float)
    delta = np.asarray(after["Q_s_current_minus_source_ideal"], dtype=float)
    q_public = np.asarray(after["Q_s_current_release2_end"], dtype=float)
    source_ideal = np.asarray(after["Q_s_source_ideal_release2_end"], dtype=float)
    bridge = _bridge_diagnostic(loaded, fine_q)

    return {
        "schema": "kokuno-a4-current-qs-release2-independent-audit-v1",
        "task": TASK,
        "upstream_head": UPSTREAM_HEAD,
        "seed": SEED,
        "candidate_semantic_sha256": _semantic_value(loaded),
        "heldout_eta": eta.tolist(),
        "eta_fd_steps": list(ETA_FD_STEPS),
        "resolutions": resolutions,
        "save_load": {
            "semantic_exact": bool(semantic_exact),
            "configuration_exact": bool(configuration_exact),
            "q_state_exact": bool(q_replay_exact),
        },
        "current_vs_source_ideal": {
            "max_abs_delta": float(np.max(np.abs(delta))),
            "rms_delta": float(np.sqrt(np.mean(delta * delta))),
            "q_current_min": float(np.min(q_public)),
            "q_current_max": float(np.max(q_public)),
            "q_source_ideal_min": float(np.min(source_ideal)),
            "q_source_ideal_max": float(np.max(source_ideal)),
        },
        "bridge_diagnostic": bridge,
        "gates": {
            "fine_q_reconstruction_normalized_max": Q_RECON_NORM_MAX_GATE,
            "fine_q_reconstruction_normalized_rms": Q_RECON_NORM_RMS_GATE,
            "fine_derivative_channel_normalized_max": DERIVATIVE_CHANNEL_NORM_MAX_GATE,
            "medium_to_fine_worsen_factor": REFINEMENT_WORSEN_FACTOR,
            "refinement_floor": REFINEMENT_FLOOR,
            "current_minus_source_ideal_min": CURRENT_MINUS_IDEAL_MIN,
            "final_project_momentum_gate_unchanged": FINAL_MOMENTUM_GATE,
            "final_project_divergence_gate_unchanged": FINAL_DIVERGENCE_GATE,
            "final_project_quadrature_unchanged": list(FINAL_QUADRATURE),
        },
        "truth_boundary": {
            "current_q_s_release2_independent_audit_executed": True,
            "current_l_minus_h_matching_bridge_materialized": False,
            "cartesian_terminal_multiplier_composed": False,
            "leading_only_ns_residual_assessed": False,
            "leading_plus_oscillatory_ns_residual_assessed": False,
            "after_correction_ns_residual_assessed": False,
            "complete_ns_admission_ready": False,
            "pde_validated": False,
        },
        "final_project_admission_ready": False,
    }


def enforce_preregistered_gates(report: Mapping[str, Any]) -> None:
    """Fail closed on the exact gates preregistered in Issue #15."""
    resolutions = list(report["resolutions"])
    if len(resolutions) != 3:
        raise AssertionError("K4-VAL-128 requires exactly three eta resolutions")
    fine = resolutions[-1]
    medium = resolutions[-2]
    if fine["q_reconstruction_normalized_max"] > Q_RECON_NORM_MAX_GATE:
        raise AssertionError("fine independent Q_s normalized max exceeds preregistered gate")
    if fine["q_reconstruction_normalized_rms"] > Q_RECON_NORM_RMS_GATE:
        raise AssertionError("fine independent Q_s normalized RMS exceeds preregistered gate")
    if fine["derivative_channel_normalized_max"] > DERIVATIVE_CHANNEL_NORM_MAX_GATE:
        raise AssertionError("fine derivative-channel disagreement exceeds preregistered gate")
    m = float(medium["q_reconstruction_normalized_max"])
    f = float(fine["q_reconstruction_normalized_max"])
    if f > max(REFINEMENT_FLOOR, REFINEMENT_WORSEN_FACTOR * m):
        raise AssertionError("medium-to-fine Q_s reconstruction worsened beyond preregistered allowance")
    replay = report["save_load"]
    if not (replay["semantic_exact"] and replay["configuration_exact"] and replay["q_state_exact"]):
        raise AssertionError("save/load replay is not exact")
    if report["current_vs_source_ideal"]["max_abs_delta"] < CURRENT_MINUS_IDEAL_MIN:
        raise AssertionError("current-vs-source-ideal state difference fell below preregistered nontriviality floor")
    truth = report["truth_boundary"]
    if truth["pde_validated"] or truth["complete_ns_admission_ready"]:
        raise AssertionError("scoped scalar audit attempted a PDE truth promotion")
    if report["final_project_admission_ready"]:
        raise AssertionError("scoped scalar audit cannot be final project admission")


def public_api_has_no_scientific_tuning_knobs() -> bool:
    sig = inspect.signature(audit_loaded_current_qs_state)
    return _FORBIDDEN_PUBLIC_KNOBS.isdisjoint(sig.parameters)
