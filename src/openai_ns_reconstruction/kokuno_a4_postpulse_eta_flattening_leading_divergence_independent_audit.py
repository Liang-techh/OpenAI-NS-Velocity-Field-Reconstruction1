"""Independent divergence audit for exact A1 #1140 post-pulse leading field.

The scientific path deliberately consumes only the save/reloaded public Cartesian
``velocity(x,y,z,t)`` surface.  Cartesian derivatives are reconstructed here
with an A4-owned centered FD2 operator at the preregistered absolute steps.
Analytic/source-coordinate derivatives exposed by the candidate are not used to
compute scientific divergence.

This is scoped leading-only divergence evidence.  It is not a complete
Navier--Stokes momentum residual and cannot promote ``pde_validated``.
"""
from __future__ import annotations

import copy
from dataclasses import dataclass
import inspect
import json
import math
from pathlib import Path
import tempfile
from typing import Any, Mapping

import numpy as np

from .kokuno_pa16_current_cartesian_postpulse_eta_flattening import (
    KokunoPA16CurrentCartesianPostPulseEtaFlattening,
)

TASK = "K4-VAL-118"
SCHEMA = "kokuno-a4-postpulse-leading-divergence-audit-v1"
UPSTREAM_PR = 1140
UPSTREAM_HEAD = "c5442c11165d7f0976889b826bf58dce38a3d3dd"
PARENT_PR = 1133
PARENT_HEAD = "3b4af71c4ab547bc11f9f3e320afa4b62c25b930"
SEED = 9173901
SPATIAL_STEPS = (0.02, 0.01, 0.005)
TIMES = (0.31, 0.47, 0.63, 0.71)
POSTPULSE_S_BANDS = ((0.08, 0.24), (0.38, 0.58), (0.72, 0.88))
LATE_S = (0.92, 0.96, 0.985)
ETA_INTERVAL = (-0.40, 0.40)
DIVERGENCE_GATE = 1.0e-5
STABILITY_FACTOR = 1.25
STABILITY_FLOOR = 2.0e-8
NONTRIVIAL_SPEED_RMS = 1.0e-10
SAVE_RELOAD_VELOCITY_TOLERANCE = 2.0e-12
ORDER_INVARIANCE_ATOL = 1.0e-13
MUTATION_EPSILON = 1.0e-3
MUTATION_DETECTION_FLOOR = 5.0e-4
FINAL_PROJECT_MOMENTUM_GATE = 1.0e-3
FINAL_PROJECT_DIVERGENCE_GATE = 1.0e-5
_LOG2 = math.log(2.0)
_LOG_MAX = math.log(np.finfo(float).max)

_TRUTH_BOUNDARY = {
    "exact_a1_1140_postpulse_leading_consumed": True,
    "candidate_save_reload_required": True,
    "public_cartesian_velocity_only_for_scientific_derivatives": True,
    "independent_cartesian_fd2_operator_used": True,
    "canonical_absolute_fd_ladder_preserved": True,
    "fd_representability_fail_closed": True,
    "leading_only_divergence_scoped_assessed": True,
    "leading_only_ns_residual_assessed": False,
    "leading_plus_oscillatory_ns_residual_assessed": False,
    "after_correction_ns_residual_assessed": False,
    "matching_postpulse_oscillatory_composite_materialized": False,
    "current_lineage_correction_velocity_materialized": False,
    "terminal_global_leading_completion_materialized": False,
    "matched_global_pressure_materialized": False,
    "restricted_forcing_materialized": False,
    "canonical_whole_domain_admission_assessed": False,
    "same_protocol_comparable_to_st006": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def default_field() -> KokunoPA16CurrentCartesianPostPulseEtaFlattening:
    return KokunoPA16CurrentCartesianPostPulseEtaFlattening()


def _D(field: Any) -> float:
    value = float(getattr(field, "D"))
    if not math.isfinite(value):
        raise RuntimeError("invalid similarity D")
    return value


def _source_to_cartesian_logx(
    field: Any, log_X: Any, eta: Any, theta: Any, t: Any
) -> np.ndarray:
    log_X, eta, theta, t = np.broadcast_arrays(
        np.asarray(log_X, float),
        np.asarray(eta, float),
        np.asarray(theta, float),
        np.asarray(t, float),
    )
    if (
        np.any(~np.isfinite(log_X))
        or np.any(~np.isfinite(eta))
        or np.any(~np.isfinite(theta))
        or np.any(~np.isfinite(t))
        or np.any(np.abs(eta) >= 1.0)
        or np.any(t >= 1.0)
    ):
        raise ValueError("invalid log-X similarity probe")
    q = (1.0 - t) / (1.0 - eta * eta)
    if np.any(q <= 0.0) or np.any(~np.isfinite(q)):
        raise ValueError("invalid q")
    log_r = 0.5 * (_LOG2 + np.log(q) + log_X)
    if np.any(log_r >= _LOG_MAX):
        raise RuntimeError("physical radius is not binary64 representable")
    r = np.exp(log_r)
    z = np.power(q, _D(field)) * eta
    return np.stack((r * np.cos(theta), r * np.sin(theta), z), axis=-1)


def _log_physical_volume_jacobian(field: Any, eta: Any, t: Any) -> np.ndarray:
    eta, t = np.broadcast_arrays(np.asarray(eta, float), np.asarray(t, float))
    q = (1.0 - t) / (1.0 - eta * eta)
    D = _D(field)
    correction = 1.0 + 2.0 * D * eta * eta / (1.0 - eta * eta)
    if (
        np.any(q <= 0.0)
        or np.any(correction <= 0.0)
        or np.any(~np.isfinite(q))
        or np.any(~np.isfinite(correction))
    ):
        raise RuntimeError("invalid physical-volume Jacobian")
    return (D + 1.0) * np.log(q) + np.log(correction)


def _lhs(rng: np.random.Generator, n: int) -> np.ndarray:
    out = (np.arange(n, dtype=float) + rng.random(n)) / n
    rng.shuffle(out)
    return out


@dataclass(frozen=True)
class HeldoutProbeSet:
    points: np.ndarray
    times: np.ndarray
    log_weights: np.ndarray
    region: np.ndarray
    time_index: np.ndarray
    seam_points: np.ndarray
    seam_times: np.ndarray
    axis_points: np.ndarray
    axis_times: np.ndarray
    late_points: np.ndarray
    late_times: np.ndarray
    late_s: np.ndarray


def _region(
    field: Any,
    rng: np.random.Generator,
    time: float,
    tid: int,
    label: str,
    s_lo: float,
    s_hi: float,
    n: int = 5,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    us, ue, ut = _lhs(rng, n), _lhs(rng, n), _lhs(rng, n)
    s = s_lo + (s_hi - s_lo) * us
    local_y = field.T_f * s
    log_X = field.log_X_pulse_end + local_y
    eta = ETA_INTERVAL[0] + (ETA_INTERVAL[1] - ETA_INTERVAL[0]) * ue
    theta = 2.0 * math.pi * ut
    tt = np.full(n, time)
    xyz = _source_to_cartesian_logx(field, log_X, eta, theta, tt)
    dlogX = field.T_f * (s_hi - s_lo)
    cell = dlogX * (ETA_INTERVAL[1] - ETA_INTERVAL[0]) * 2.0 * math.pi / n
    log_weight = math.log(cell) + log_X + _log_physical_volume_jacobian(field, eta, tt)
    return xyz, tt, log_weight, np.full(n, label, object), np.full(n, tid, int)


def make_heldout_probes(field: Any) -> HeldoutProbeSet:
    rng = np.random.default_rng(SEED)
    blocks = []
    for tid, time in enumerate(TIMES):
        for zid, (lo, hi) in enumerate(POSTPULSE_S_BANDS):
            blocks.append(
                _region(
                    field,
                    rng,
                    time,
                    tid,
                    f"postpulse_flattening_zone_{zid}",
                    lo,
                    hi,
                )
            )
    points = np.concatenate([b[0] for b in blocks])
    times = np.concatenate([b[1] for b in blocks])
    log_weights = np.concatenate([b[2] for b in blocks])
    regions = np.concatenate([b[3] for b in blocks])
    tids = np.concatenate([b[4] for b in blocks])

    seam_points, seam_times = [], []
    for time in (0.47, 0.63):
        for eta0 in (-0.22, 0.22):
            local_y = np.asarray((-0.08, -0.02, 0.02, 0.08), dtype=float)
            log_X = field.log_X_pulse_end + local_y
            eta = np.full(4, eta0)
            theta = np.asarray((0.37, 1.11, 2.03, 2.81))
            tt = np.full(4, time)
            seam_points.append(_source_to_cartesian_logx(field, log_X, eta, theta, tt))
            seam_times.append(tt)

    axis_points = np.asarray(
        (
            (0.0, 0.0, 0.0),
            (0.0, 0.0, 0.10),
            (1.0e-8, 0.0, -0.10),
            (0.0, -1.0e-6, 0.05),
            (1.0e-4, 1.0e-4, -0.03),
            (1.0e-2, 0.0, 0.02),
        ),
        dtype=float,
    )
    axis_times = np.asarray((0.31, 0.47, 0.63, 0.71, 0.39, 0.55), dtype=float)

    late_s = np.asarray(LATE_S, dtype=float)
    late_log_X = field.log_X_pulse_end + field.T_f * late_s
    late_eta = np.asarray((-0.17, 0.08, 0.19))
    late_theta = np.asarray((0.41, 1.37, 2.41))
    late_times = np.asarray((0.43, 0.55, 0.67))
    late_points = _source_to_cartesian_logx(
        field, late_log_X, late_eta, late_theta, late_times
    )
    return HeldoutProbeSet(
        points=points,
        times=times,
        log_weights=log_weights,
        region=regions,
        time_index=tids,
        seam_points=np.concatenate(seam_points),
        seam_times=np.concatenate(seam_times),
        axis_points=axis_points,
        axis_times=axis_times,
        late_points=late_points,
        late_times=late_times,
        late_s=late_s,
    )


def _velocity(field: Any, points: np.ndarray, times: np.ndarray) -> np.ndarray:
    p, tt = np.asarray(points, float), np.asarray(times, float)
    out = np.asarray(field.velocity(p[:, 0], p[:, 1], p[:, 2], tt), float)
    if out.shape != (len(p), 3) or np.any(~np.isfinite(out)):
        raise RuntimeError("invalid public velocity output")
    return out


def fd_representability(points: np.ndarray, step: float) -> tuple[np.ndarray, dict[str, Any]]:
    p = np.asarray(points, dtype=float)
    h = float(step)
    if p.ndim != 2 or p.shape[1] != 3 or np.any(~np.isfinite(p)) or h <= 0.0:
        raise ValueError("invalid representability input")
    ok = np.ones(len(p), dtype=bool)
    per_axis = []
    worst_ratio = 0.0
    worst = None
    for axis in range(3):
        center = p[:, axis]
        plus, minus = center + h, center - h
        changed = (plus != center) & (minus != center) & (plus != minus)
        ok &= changed
        spacing = np.spacing(np.abs(center))
        ratio = np.divide(spacing, h, out=np.zeros_like(spacing), where=np.isfinite(spacing))
        if len(ratio):
            iw = int(np.argmax(ratio))
            if float(ratio[iw]) >= worst_ratio:
                worst_ratio = float(ratio[iw])
                worst = {
                    "point_index": iw,
                    "axis": axis,
                    "coordinate": float(center[iw]),
                    "spacing_over_step": worst_ratio,
                }
        per_axis.append({"axis": axis, "nonrepresentable_count": int(np.count_nonzero(~changed))})
    return ok, {
        "step": h,
        "probe_count": len(p),
        "representable_count": int(np.count_nonzero(ok)),
        "all_probes_representable": bool(np.all(ok)),
        "per_axis": per_axis,
        "worst_spacing_over_step": worst,
    }


def independent_fd2_jacobian(field: Any, points: np.ndarray, times: np.ndarray, step: float) -> np.ndarray:
    p, tt = np.asarray(points, float), np.asarray(times, float)
    if p.ndim != 2 or p.shape[1] != 3 or tt.shape != (len(p),):
        raise ValueError("points/times shape mismatch")
    mask, _ = fd_representability(p, step)
    if not np.all(mask):
        raise RuntimeError("canonical Cartesian FD perturbation is not representable")
    h = float(step)
    jac = np.empty((len(p), 3, 3), dtype=float)
    for axis in range(3):
        plus, minus = p.copy(), p.copy()
        plus[:, axis] += h
        minus[:, axis] -= h
        jac[:, :, axis] = (_velocity(field, plus, tt) - _velocity(field, minus, tt)) / (2.0 * h)
    return jac


def _divergence_on_mask(
    field: Any, points: np.ndarray, times: np.ndarray, step: float, mask: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    use = np.asarray(mask, dtype=bool)
    if use.shape != (len(points),):
        raise ValueError("mask shape mismatch")
    if not np.any(use):
        return np.empty(0), np.empty((0, 3, 3))
    jac = independent_fd2_jacobian(field, np.asarray(points)[use], np.asarray(times)[use], step)
    return np.trace(jac, axis1=1, axis2=2), jac


def _stable_weighted_metrics(values: np.ndarray, log_weights: np.ndarray) -> dict[str, Any]:
    values, lw = np.asarray(values, float), np.asarray(log_weights, float)
    if values.shape != lw.shape or values.size == 0:
        raise ValueError("weighted metric shape mismatch")
    shift = float(np.max(lw))
    w = np.exp(lw - shift)
    sw = float(np.sum(w))
    sumsq = float(np.sum(w * values * values))
    if not (sw > 0.0 and math.isfinite(sw) and math.isfinite(sumsq)):
        raise RuntimeError("invalid scaled physical weights")
    rms = math.sqrt(sumsq / sw)
    log_volume = shift + math.log(sw)
    log_l2 = -math.inf if sumsq == 0.0 else 0.5 * (shift + math.log(sumsq))
    return {
        "estimated_volume": math.exp(log_volume) if log_volume < _LOG_MAX else None,
        "log_estimated_volume": log_volume,
        "volume_l2_estimate": 0.0 if sumsq == 0.0 else (math.exp(log_l2) if log_l2 < _LOG_MAX else None),
        "log_volume_l2_estimate": log_l2,
        "weighted_rms": rms,
    }


def _metrics(field: Any, probes: HeldoutProbeSet, step: float) -> dict[str, Any]:
    representable, rep = fd_representability(probes.points, step)
    div, jac = _divergence_on_mask(field, probes.points, probes.times, step, representable)
    if div.size:
        valid_points = probes.points[representable]
        valid_times = probes.times[representable]
        valid_regions = probes.region[representable]
        frob = np.linalg.norm(jac, axis=(1, 2))
        normed = np.abs(div) / np.maximum(1.0, frob)
        weighted = _stable_weighted_metrics(div, probes.log_weights[representable])
        norm_weighted = _stable_weighted_metrics(normed, probes.log_weights[representable])
        speed = np.linalg.norm(_velocity(field, valid_points, valid_times), axis=1)
        iw = int(np.argmax(np.abs(div)))
        worst = {
            "point": valid_points[iw].tolist(),
            "time": float(valid_times[iw]),
            "region": str(valid_regions[iw]),
            "divergence": float(div[iw]),
            "jacobian_frobenius": float(frob[iw]),
            "normalized_divergence": float(normed[iw]),
        }
        sampled_max = float(np.max(np.abs(div)))
        normalized_max = float(np.max(normed))
        speed_rms = float(np.sqrt(np.mean(speed**2)))
    else:
        weighted = {"estimated_volume": None, "log_estimated_volume": None, "volume_l2_estimate": None, "log_volume_l2_estimate": None, "weighted_rms": None}
        norm_weighted = copy.deepcopy(weighted)
        worst = None
        sampled_max = normalized_max = speed_rms = None

    by_region: dict[str, Any] = {}
    for label in sorted(set(map(str, probes.region.tolist()))):
        all_mask = probes.region == label
        mask = all_mask & representable
        item = {
            "probe_count": int(np.count_nonzero(all_mask)),
            "representable_count": int(np.count_nonzero(mask)),
        }
        if np.any(mask):
            rd, _ = _divergence_on_mask(field, probes.points, probes.times, step, mask)
            item["sampled_max"] = float(np.max(np.abs(rd)))
            item.update(_stable_weighted_metrics(rd, probes.log_weights[mask]))
        else:
            item.update({"sampled_max": None, "weighted_rms": None, "volume_l2_estimate": None, "log_volume_l2_estimate": None})
        by_region[label] = item

    by_time = []
    for tid, time in enumerate(TIMES):
        all_mask = probes.time_index == tid
        mask = all_mask & representable
        item = {"time": time, "probe_count": int(np.count_nonzero(all_mask)), "representable_count": int(np.count_nonzero(mask))}
        if np.any(mask):
            td, _ = _divergence_on_mask(field, probes.points, probes.times, step, mask)
            item["sampled_max"] = float(np.max(np.abs(td)))
            item.update(_stable_weighted_metrics(td, probes.log_weights[mask]))
        else:
            item.update({"sampled_max": None, "weighted_rms": None, "log_volume_l2_estimate": None})
        by_time.append(item)

    seam_mask, seam_rep = fd_representability(probes.seam_points, step)
    seam_div, _ = _divergence_on_mask(field, probes.seam_points, probes.seam_times, step, seam_mask)
    axis_mask, axis_rep = fd_representability(probes.axis_points, step)
    axis_div, _ = _divergence_on_mask(field, probes.axis_points, probes.axis_times, step, axis_mask)
    late_mask, late_rep = fd_representability(probes.late_points, step)
    late_div, _ = _divergence_on_mask(field, probes.late_points, probes.late_times, step, late_mask)
    return {
        "sampled_max": sampled_max,
        "pooled_weighted_rms": weighted["weighted_rms"],
        "pooled_volume_l2_estimate": weighted["volume_l2_estimate"],
        "pooled_log_volume_l2_estimate": weighted["log_volume_l2_estimate"],
        "normalized_sampled_max": normalized_max,
        "normalized_weighted_rms": norm_weighted["weighted_rms"],
        "speed_rms": speed_rms,
        "main_fd_representability": rep,
        "pulse_end_seam_fd_representability": seam_rep,
        "axis_fd_representability": axis_rep,
        "late_stage_fd_representability": late_rep,
        "pulse_end_seam_sampled_max": float(np.max(np.abs(seam_div))) if seam_div.size else None,
        "axis_axis_near_sampled_max": float(np.max(np.abs(axis_div))) if axis_div.size else None,
        "late_stage_sampled_max": float(np.max(np.abs(late_div))) if late_div.size else None,
        "per_region": by_region,
        "per_time": by_time,
        "worst_witness": worst,
    }


class _XLinearMutation:
    def __init__(self, base: Any, epsilon: float = MUTATION_EPSILON):
        self.base, self.epsilon = base, float(epsilon)

    def velocity(self, x, y, z, t):
        out = np.asarray(self.base.velocity(x, y, z, t), float).copy()
        out[..., 0] += self.epsilon * np.asarray(x, float)
        return out


class _ManufacturedSolenoidalField:
    def velocity(self, x, y, z, t):
        x, y, z, t = np.broadcast_arrays(np.asarray(x, float), np.asarray(y, float), np.asarray(z, float), np.asarray(t, float))
        return np.stack(((1.0 + t) * y, -(1.0 + t) * x, np.zeros_like(z)), axis=-1)


def _manufactured_calibration() -> dict[str, Any]:
    points = np.asarray(((0.2, -0.3, 0.1), (-0.4, 0.1, -0.2), (0.6, 0.2, 0.3)), dtype=float)
    times = np.asarray((0.31, 0.47, 0.71), dtype=float)
    maxima = []
    for h in SPATIAL_STEPS:
        jac = independent_fd2_jacobian(_ManufacturedSolenoidalField(), points, times, h)
        maxima.append(float(np.max(np.abs(np.trace(jac, axis1=1, axis2=2)))))
    return {"max_abs_divergence_by_step": maxima, "passed": bool(max(maxima) <= 1.0e-12)}


def _mutation_detection() -> dict[str, Any]:
    base = _ManufacturedSolenoidalField()
    mutated = _XLinearMutation(base)
    points = np.asarray(((0.2, -0.3, 0.1), (-0.4, 0.1, -0.2), (0.6, 0.2, 0.3)), dtype=float)
    times = np.asarray((0.31, 0.47, 0.71), dtype=float)
    jac = independent_fd2_jacobian(mutated, points, times, SPATIAL_STEPS[-1])
    div = np.trace(jac, axis1=1, axis2=2)
    signal = float(np.min(np.abs(div)))
    return {"minimum_abs_detected_divergence": signal, "passed": bool(signal >= MUTATION_DETECTION_FLOOR)}


def _configuration_mutation_rejected(saved: Mapping[str, Any]) -> bool:
    raw = copy.deepcopy(dict(saved))
    raw["bound_scope"]["T_f_role"] = "source_exact"
    try:
        KokunoPA16CurrentCartesianPostPulseEtaFlattening.from_configuration(raw)
    except (ValueError, RuntimeError):
        return True
    return False


def _post_stage_fails_closed(field: Any) -> bool:
    try:
        field.similarity_profile_values_logX(field.log_X_flatten_end + 0.02, 0.0)
    except (ValueError, RuntimeError):
        return True
    return False


def _stability_ok(medium: Any, fine: Any) -> bool:
    if medium is None or fine is None:
        return False
    return bool(float(fine) <= STABILITY_FACTOR * float(medium) + STABILITY_FLOOR)


def _numeric_leq(value: Any, threshold: float) -> bool:
    return value is not None and math.isfinite(float(value)) and float(value) <= threshold


def _all_fd_representable(metrics: Mapping[str, Any]) -> bool:
    return all(
        bool(metrics[key]["all_probes_representable"])
        for key in (
            "main_fd_representability",
            "pulse_end_seam_fd_representability",
            "axis_fd_representability",
            "late_stage_fd_representability",
        )
    )


def _order_invariance(field: Any, probes: HeldoutProbeSet) -> bool:
    h = SPATIAL_STEPS[-1]
    mask, _ = fd_representability(probes.points, h)
    if not np.any(mask):
        return False
    p, tt = probes.points[mask], probes.times[mask]
    d0 = np.trace(independent_fd2_jacobian(field, p, tt, h), axis1=1, axis2=2)
    d1 = np.trace(independent_fd2_jacobian(field, p[::-1], tt[::-1], h), axis1=1, axis2=2)[::-1]
    return bool(np.max(np.abs(d0 - d1)) <= ORDER_INVARIANCE_ATOL)


def _public_api_has_no_tuning_knobs() -> bool:
    forbidden = {"threshold", "residual", "forcing", "pressure", "viscosity", "step", "tolerance", "gain", "target"}
    for fn in (materialize_receipt, enforce_preregistered_gates):
        names = {name.lower() for name in inspect.signature(fn).parameters}
        if any(any(token in name for token in forbidden) for name in names):
            return False
    return True


def materialize_receipt() -> dict[str, Any]:
    original = default_field()
    with tempfile.TemporaryDirectory() as td:
        path = Path(td) / "candidate.json"
        saved = original.save_configuration(path)
        field = KokunoPA16CurrentCartesianPostPulseEtaFlattening.load_configuration(path)
        probes = make_heldout_probes(field)
        replay_points = np.concatenate((probes.axis_points, probes.seam_points[:4]))
        replay_times = np.concatenate((probes.axis_times, probes.seam_times[:4]))
        replay_error = float(np.max(np.abs(_velocity(original, replay_points, replay_times) - _velocity(field, replay_points, replay_times))))
        mutation_rejected = _configuration_mutation_rejected(saved)

    levels = {f"h={h:g}": _metrics(field, probes, h) for h in SPATIAL_STEPS}
    coarse = levels[f"h={SPATIAL_STEPS[0]:g}"]
    medium = levels[f"h={SPATIAL_STEPS[1]:g}"]
    fine = levels[f"h={SPATIAL_STEPS[2]:g}"]
    calibration = _manufactured_calibration()
    mutation = _mutation_detection()
    gates = {
        "divergence_gate": DIVERGENCE_GATE,
        "fine_sampled_max_pass": _numeric_leq(fine["sampled_max"], DIVERGENCE_GATE),
        "fine_weighted_rms_pass": _numeric_leq(fine["pooled_weighted_rms"], DIVERGENCE_GATE),
        "pulse_end_seam_pass": _numeric_leq(fine["pulse_end_seam_sampled_max"], DIVERGENCE_GATE),
        "axis_near_pass": _numeric_leq(fine["axis_axis_near_sampled_max"], DIVERGENCE_GATE),
        "late_stage_pass": _numeric_leq(fine["late_stage_sampled_max"], DIVERGENCE_GATE),
        "nontrivial_speed_pass": fine["speed_rms"] is not None and float(fine["speed_rms"]) >= NONTRIVIAL_SPEED_RMS,
        "medium_to_fine_sampled_max_stable": _stability_ok(medium["sampled_max"], fine["sampled_max"]),
        "medium_to_fine_weighted_rms_stable": _stability_ok(medium["pooled_weighted_rms"], fine["pooled_weighted_rms"]),
        "all_canonical_fd_perturbations_representable": all(_all_fd_representable(levels[f"h={h:g}"]) for h in SPATIAL_STEPS),
        "save_reload_replay_pass": replay_error <= SAVE_RELOAD_VELOCITY_TOLERANCE,
        "manufactured_solenoidal_calibration_pass": calibration["passed"],
        "mutation_detection_pass": mutation["passed"],
        "ordering_invariance_pass": _order_invariance(field, probes),
        "configuration_mutation_rejected": mutation_rejected,
        "post_stage_fail_closed": _post_stage_fails_closed(field),
        "public_api_tuning_knobs_absent": _public_api_has_no_tuning_knobs(),
        "final_project_momentum_gate_unchanged": FINAL_PROJECT_MOMENTUM_GATE,
        "final_project_divergence_gate_unchanged": FINAL_PROJECT_DIVERGENCE_GATE,
    }
    boolean_gate_names = [key for key, value in gates.items() if isinstance(value, bool)]
    audit_pass = all(bool(gates[key]) for key in boolean_gate_names)
    return {
        "schema": SCHEMA,
        "task": TASK,
        "upstream_pr": UPSTREAM_PR,
        "upstream_head": UPSTREAM_HEAD,
        "parent_pr": PARENT_PR,
        "parent_head": PARENT_HEAD,
        "seed": SEED,
        "spatial_steps": list(SPATIAL_STEPS),
        "times": list(TIMES),
        "postpulse_s_bands": [list(v) for v in POSTPULSE_S_BANDS],
        "late_s": list(LATE_S),
        "candidate_semantic_sha256": field.semantic_sha256,
        "save_reload_velocity_max_abs_error": replay_error,
        "probe_counts": {
            "main": len(probes.points),
            "seam": len(probes.seam_points),
            "axis_near": len(probes.axis_points),
            "late_stage": len(probes.late_points),
        },
        "resolution_metrics": levels,
        "manufactured_calibration": calibration,
        "mutation_detection": mutation,
        "gates": gates,
        "audit_pass": audit_pass,
        "truth_boundary": dict(_TRUTH_BOUNDARY),
        "limitations": [
            "scoped leading-only divergence audit; no complete NS momentum residual is computed",
            "A2 has no matching self-contained post-pulse composite on exact A1 #1140",
            "A3 has no same-identity Cartesian correction velocity on this post-pulse candidate",
            "matched pressure, preregistered restricted forcing, terminal/global leading completion, and whole-domain admission remain absent",
            "binary64 failure of any canonical Cartesian perturbation is a scientific failure, not permission to retune the FD ladder",
        ],
    }


def enforce_preregistered_gates(report: Mapping[str, Any]) -> None:
    if report.get("schema") != SCHEMA or report.get("task") != TASK:
        raise AssertionError("unexpected A4 report identity")
    if report.get("upstream_head") != UPSTREAM_HEAD:
        raise AssertionError("wrong upstream identity")
    gates = report.get("gates")
    if not isinstance(gates, Mapping):
        raise AssertionError("missing preregistered gates")
    required = (
        "fine_sampled_max_pass",
        "fine_weighted_rms_pass",
        "pulse_end_seam_pass",
        "axis_near_pass",
        "late_stage_pass",
        "nontrivial_speed_pass",
        "medium_to_fine_sampled_max_stable",
        "medium_to_fine_weighted_rms_stable",
        "all_canonical_fd_perturbations_representable",
        "save_reload_replay_pass",
        "manufactured_solenoidal_calibration_pass",
        "mutation_detection_pass",
        "ordering_invariance_pass",
        "configuration_mutation_rejected",
        "post_stage_fail_closed",
        "public_api_tuning_knobs_absent",
    )
    failed = [name for name in required if gates.get(name) is not True]
    if float(gates.get("divergence_gate", math.nan)) != DIVERGENCE_GATE:
        failed.append("divergence_gate_mutated")
    if float(gates.get("final_project_momentum_gate_unchanged", math.nan)) != FINAL_PROJECT_MOMENTUM_GATE:
        failed.append("momentum_gate_mutated")
    if float(gates.get("final_project_divergence_gate_unchanged", math.nan)) != FINAL_PROJECT_DIVERGENCE_GATE:
        failed.append("project_divergence_gate_mutated")
    truth = report.get("truth_boundary", {})
    for key in (
        "leading_only_ns_residual_assessed",
        "leading_plus_oscillatory_ns_residual_assessed",
        "after_correction_ns_residual_assessed",
        "pde_validated",
    ):
        if truth.get(key) is not False:
            failed.append(f"truth_boundary:{key}")
    if failed:
        raise AssertionError("preregistered A4 gate failure: " + ", ".join(failed))


__all__ = [
    "TASK",
    "SCHEMA",
    "SEED",
    "SPATIAL_STEPS",
    "DIVERGENCE_GATE",
    "default_field",
    "make_heldout_probes",
    "fd_representability",
    "independent_fd2_jacobian",
    "materialize_receipt",
    "enforce_preregistered_gates",
]
