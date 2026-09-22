"""Independent divergence audit for exact A2 #1117 log-X main-pulse composite.

Scientific path:
1. exact candidate save/load;
2. consume only public Cartesian velocity(x,y,z,t) from total and authenticated
   leading backend;
3. reconstruct Cartesian derivatives with an A4-owned centered FD2 ladder.

The late log-X tail is intentionally audited with the canonical absolute
Cartesian steps.  If binary64 cannot represent x +/- h, the scientific gate
fails closed rather than interpreting a rounded zero derivative as evidence
of incompressibility.
"""
from __future__ import annotations

import argparse
import copy
from dataclasses import dataclass
import json
import math
from pathlib import Path
import tempfile
from typing import Any, Mapping

import numpy as np

from . import kokuno_current_main_pulse_logx_leading_oscillatory_identity as _candidate

CurrentMainPulseLogXLeadingOscillatoryField = (
    _candidate.CurrentMainPulseLogXLeadingOscillatoryField
)
default_field = _candidate.default_field

TASK = "K4-VAL-117"
SCHEMA = "kokuno-a4-logx-xi11-composite-divergence-audit-v1"
UPSTREAM_PR = 1117
UPSTREAM_HEAD = "27741d9c0a27262f7fabf61eebaa2fbd507e9f03"
AGENT1_PR = 1107
AGENT1_HEAD = "45da043dd2b4cd067f005a72c7e21fd0f2bcf309"
NEWER_AGENT1_PR = 1133
NEWER_AGENT1_HEAD = "3b4af71c4ab547bc11f9f3e320afa4b62c25b930"

SEED = 9173891
SPATIAL_STEPS = (0.02, 0.01, 0.005)
TIMES = (0.31, 0.47, 0.63, 0.71)
PRE_PULSE_XI_BAND = (-0.55, -0.12)
XI_BANDS = ((0.7, 2.1), (4.2, 5.8), (8.3, 10.2))
ETA_INTERVAL = (-0.40, 0.40)
TAIL_XI = (9.7, 10.3, 10.8)
DIVERGENCE_GATE = 1.0e-5
STABILITY_FACTOR = 1.25
STABILITY_FLOOR = 2.0e-8
NONTRIVIAL_SPEED_RMS = 1.0e-10
OSCILLATORY_SIGNAL_RMS_MIN = 1.0e-12
MUTATION_EPSILON = 1.0e-3
MUTATION_DETECTION_FLOOR = 5.0e-4
SAVE_RELOAD_VELOCITY_TOLERANCE = 2.0e-12
ORDER_INVARIANCE_ATOL = 1.0e-13
SOURCE_XI_END = 11.0
_LOG2 = math.log(2.0)
_LOG_MAX = math.log(np.finfo(float).max)

_INNER_POINTS = np.asarray(
    [
        (0.32, 0.11, -0.20),
        (0.41, -0.17, -0.08),
        (-0.36, 0.24, 0.05),
        (-0.52, -0.16, 0.16),
        (0.58, 0.21, 0.22),
        (0.47, -0.31, -0.14),
    ],
    dtype=float,
)
_INNER_TIMES = np.asarray((0.31, 0.39, 0.47, 0.55, 0.63, 0.71), dtype=float)

_TRUTH_BOUNDARY = {
    "current_logx_xi11_leading_plus_oscillatory_velocity_consumed": True,
    "candidate_save_reload_required": True,
    "independent_cartesian_fd_operator_used": True,
    "canonical_absolute_fd_ladder_preserved": True,
    "fd_representability_fail_closed": True,
    "staged_leading_total_increment_divergence_scoped_assessed": True,
    "source_exact_main_pulse_amplitude_materialized": False,
    "source_pulse_end_MJ_corrections_materialized": False,
    "newer_agent1_1133_pulse_end_leading_consumed": False,
    "terminal_global_leading_completion_materialized": False,
    "leading_only_ns_residual_assessed": False,
    "leading_plus_oscillatory_ns_residual_assessed": False,
    "after_correction_ns_residual_assessed": False,
    "current_lineage_correction_velocity_materialized": False,
    "matched_global_pressure_materialized": False,
    "restricted_forcing_materialized": False,
    "canonical_whole_domain_admission_assessed": False,
    "same_protocol_comparable_to_st006": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _leading(field: Any) -> Any:
    backend = getattr(field, "leading_backend", None)
    if backend is None:
        raise RuntimeError("candidate has no authenticated leading backend")
    return backend


def _D(field: Any) -> float:
    value = float(getattr(_leading(field), "D"))
    if not math.isfinite(value):
        raise RuntimeError("invalid similarity D")
    return value


def _lambda(field: Any) -> float:
    value = float(getattr(_leading(field), "lambda_value"))
    if not (math.isfinite(value) and 0.0 < value < 1.0):
        raise RuntimeError("invalid current lambda")
    return value


def _log_X_p(field: Any) -> float:
    value = float(getattr(_leading(field), "log_X_p"))
    if not math.isfinite(value):
        raise RuntimeError("invalid log X_p")
    return value


def _log_X_from_xi(field: Any, xi: Any) -> np.ndarray:
    values = np.asarray(xi, dtype=float)
    if np.any(~np.isfinite(values)):
        raise ValueError("xi must be finite")
    return _log_X_p(field) + values / _lambda(field)


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
    tail_points: np.ndarray
    tail_times: np.ndarray
    tail_xi: np.ndarray


def _region(
    field: Any,
    rng: np.random.Generator,
    time: float,
    tid: int,
    label: str,
    xi_lo: float,
    xi_hi: float,
    n: int = 5,
):
    ux, ue, ut = _lhs(rng, n), _lhs(rng, n), _lhs(rng, n)
    xi = xi_lo + (xi_hi - xi_lo) * ux
    log_X = _log_X_from_xi(field, xi)
    eta = ETA_INTERVAL[0] + (ETA_INTERVAL[1] - ETA_INTERVAL[0]) * ue
    theta = 2.0 * math.pi * ut
    tt = np.full(n, time)
    xyz = _source_to_cartesian_logx(field, log_X, eta, theta, tt)

    # d log X = d xi / lambda. Store logarithmic weights so the physical
    # volume integral remains representable even when X itself would overflow.
    cell = (
        (xi_hi - xi_lo)
        / _lambda(field)
        * (ETA_INTERVAL[1] - ETA_INTERVAL[0])
        * 2.0
        * math.pi
        / n
    )
    log_weight = (
        math.log(cell)
        + log_X
        + _log_physical_volume_jacobian(field, eta, tt)
    )
    return (
        xyz,
        tt,
        log_weight,
        np.full(n, label, object),
        np.full(n, tid, int),
    )


def make_heldout_probes(field: Any) -> HeldoutProbeSet:
    rng = np.random.default_rng(SEED)
    blocks = []
    for tid, time in enumerate(TIMES):
        blocks.append(
            _region(
                field,
                rng,
                time,
                tid,
                "pre_pulse_control",
                PRE_PULSE_XI_BAND[0],
                PRE_PULSE_XI_BAND[1],
            )
        )
        for zid, (lo, hi) in enumerate(XI_BANDS):
            blocks.append(
                _region(
                    field,
                    rng,
                    time,
                    tid,
                    f"logx_main_pulse_zone_{zid}",
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
            xi = np.asarray((-0.08, -0.02, 0.02, 0.08))
            log_X = _log_X_from_xi(field, xi)
            eta = np.full(4, eta0)
            theta = np.asarray((0.37, 1.11, 2.03, 2.81))
            tt = np.full(4, time)
            seam_points.append(
                _source_to_cartesian_logx(field, log_X, eta, theta, tt)
            )
            seam_times.append(tt)

    # Axis and near-axis checks remain in ordinary physical coordinates. They
    # are deliberately independent of the enormous log-X tail.
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

    tail_xi = np.asarray(TAIL_XI, dtype=float)
    tail_log_X = _log_X_from_xi(field, tail_xi)
    tail_eta = np.asarray((-0.17, 0.08, 0.19))
    tail_theta = np.asarray((0.41, 1.37, 2.41))
    tail_times = np.asarray((0.43, 0.55, 0.67))
    tail_points = _source_to_cartesian_logx(
        field, tail_log_X, tail_eta, tail_theta, tail_times
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
        tail_points=tail_points,
        tail_times=tail_times,
        tail_xi=tail_xi,
    )


def _velocity(field: Any, points: np.ndarray, times: np.ndarray) -> np.ndarray:
    p, tt = np.asarray(points, float), np.asarray(times, float)
    out = np.asarray(field.velocity(p[:, 0], p[:, 1], p[:, 2], tt), float)
    if out.shape != (len(p), 3) or np.any(~np.isfinite(out)):
        raise RuntimeError("invalid public velocity output")
    return out


def _fd_representability(
    points: np.ndarray, step: float
) -> tuple[np.ndarray, dict[str, Any]]:
    p = np.asarray(points, dtype=float)
    h = float(step)
    if (
        p.ndim != 2
        or p.shape[1] != 3
        or np.any(~np.isfinite(p))
        or not math.isfinite(h)
        or h <= 0.0
    ):
        raise ValueError("invalid representability input")
    per_axis = []
    ok = np.ones(len(p), dtype=bool)
    worst_ratio = 0.0
    worst = None
    for axis in range(3):
        center = p[:, axis]
        plus = center + h
        minus = center - h
        changed = (plus != center) & (minus != center) & (plus != minus)
        ok &= changed
        spacing = np.spacing(np.abs(center))
        ratio = np.divide(
            spacing,
            h,
            out=np.zeros_like(spacing),
            where=np.isfinite(spacing),
        )
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
        per_axis.append(
            {
                "axis": axis,
                "nonrepresentable_count": int(np.count_nonzero(~changed)),
            }
        )
    return ok, {
        "step": h,
        "probe_count": len(p),
        "representable_count": int(np.count_nonzero(ok)),
        "all_probes_representable": bool(np.all(ok)),
        "per_axis": per_axis,
        "worst_spacing_over_step": worst,
    }


def independent_fd2_jacobian(
    field: Any, points: np.ndarray, times: np.ndarray, step: float
) -> np.ndarray:
    p, tt = np.asarray(points, float), np.asarray(times, float)
    if p.ndim != 2 or p.shape[1] != 3 or tt.shape != (len(p),):
        raise ValueError("points/times shape mismatch")
    mask, _ = _fd_representability(p, step)
    if not np.all(mask):
        raise RuntimeError("canonical Cartesian FD perturbation is not representable")
    h = float(step)
    jac = np.empty((len(p), 3, 3))
    for axis in range(3):
        plus, minus = p.copy(), p.copy()
        plus[:, axis] += h
        minus[:, axis] -= h
        jac[:, :, axis] = (
            _velocity(field, plus, tt) - _velocity(field, minus, tt)
        ) / (2.0 * h)
    return jac


def _divergence_on_mask(
    field: Any,
    points: np.ndarray,
    times: np.ndarray,
    step: float,
    mask: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    p, tt = np.asarray(points, float), np.asarray(times, float)
    use = np.asarray(mask, dtype=bool)
    if use.shape != (len(p),):
        raise ValueError("mask shape mismatch")
    if not np.any(use):
        return np.empty(0), np.empty((0, 3, 3))
    jac = independent_fd2_jacobian(field, p[use], tt[use], step)
    return np.trace(jac, axis1=1, axis2=2), jac


def _stable_weighted_metrics(
    values: np.ndarray, log_weights: np.ndarray
) -> dict[str, Any]:
    values = np.asarray(values, float)
    lw = np.asarray(log_weights, float)
    if values.shape != lw.shape or values.size == 0:
        raise ValueError("weighted metric shape mismatch")
    shift = float(np.max(lw))
    w = np.exp(lw - shift)
    sw = float(np.sum(w))
    if not (sw > 0.0 and math.isfinite(sw)):
        raise RuntimeError("invalid scaled physical weights")
    sumsq = float(np.sum(w * values * values))
    rms = math.sqrt(sumsq / sw)
    log_volume = shift + math.log(sw)
    if sumsq == 0.0:
        log_l2 = -math.inf
        l2 = 0.0
    else:
        log_l2 = 0.5 * (shift + math.log(sumsq))
        l2 = math.exp(log_l2) if log_l2 < _LOG_MAX else None
    volume = math.exp(log_volume) if log_volume < _LOG_MAX else None
    return {
        "estimated_volume": volume,
        "log_estimated_volume": log_volume,
        "volume_l2_estimate": l2,
        "log_volume_l2_estimate": log_l2,
        "weighted_rms": rms,
    }


def _metrics(
    field: Any,
    probes: HeldoutProbeSet,
    step: float,
) -> dict[str, Any]:
    representable, rep = _fd_representability(probes.points, step)
    div, jac = _divergence_on_mask(
        field, probes.points, probes.times, step, representable
    )
    valid_logw = probes.log_weights[representable]
    if div.size:
        frob = np.linalg.norm(jac, axis=(1, 2))
        normed = np.abs(div) / np.maximum(1.0, frob)
        weighted = _stable_weighted_metrics(div, valid_logw)
        norm_weighted = _stable_weighted_metrics(normed, valid_logw)
        valid_points = probes.points[representable]
        valid_times = probes.times[representable]
        valid_regions = probes.region[representable]
        speed = np.linalg.norm(
            _velocity(field, valid_points, valid_times), axis=1
        )
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
        weighted = {
            "estimated_volume": None,
            "log_estimated_volume": None,
            "volume_l2_estimate": None,
            "log_volume_l2_estimate": None,
            "weighted_rms": None,
        }
        norm_weighted = copy.deepcopy(weighted)
        worst = None
        sampled_max = normalized_max = speed_rms = None

    by_region: dict[str, Any] = {}
    for label in sorted(set(map(str, probes.region.tolist()))):
        region_all = probes.region == label
        region_valid = region_all & representable
        item = {
            "probe_count": int(np.count_nonzero(region_all)),
            "representable_count": int(np.count_nonzero(region_valid)),
            "representable_fraction": float(
                np.count_nonzero(region_valid) / np.count_nonzero(region_all)
            ),
        }
        if np.any(region_valid):
            rd, _ = _divergence_on_mask(
                field, probes.points, probes.times, step, region_valid
            )
            item["sampled_max"] = float(np.max(np.abs(rd)))
            item.update(
                _stable_weighted_metrics(rd, probes.log_weights[region_valid])
            )
        else:
            item.update(
                {
                    "sampled_max": None,
                    "estimated_volume": None,
                    "log_estimated_volume": None,
                    "volume_l2_estimate": None,
                    "log_volume_l2_estimate": None,
                    "weighted_rms": None,
                }
            )
        by_region[label] = item

    by_time = []
    for tid, time in enumerate(TIMES):
        time_all = probes.time_index == tid
        time_valid = time_all & representable
        item = {
            "time": time,
            "probe_count": int(np.count_nonzero(time_all)),
            "representable_count": int(np.count_nonzero(time_valid)),
        }
        if np.any(time_valid):
            td, _ = _divergence_on_mask(
                field, probes.points, probes.times, step, time_valid
            )
            item["sampled_max"] = float(np.max(np.abs(td)))
            item.update(
                _stable_weighted_metrics(td, probes.log_weights[time_valid])
            )
        else:
            item["sampled_max"] = None
            item["weighted_rms"] = None
            item["log_volume_l2_estimate"] = None
        by_time.append(item)

    seam_mask, seam_rep = _fd_representability(probes.seam_points, step)
    seam_div, _ = _divergence_on_mask(
        field, probes.seam_points, probes.seam_times, step, seam_mask
    )
    axis_mask, axis_rep = _fd_representability(probes.axis_points, step)
    axis_div, _ = _divergence_on_mask(
        field, probes.axis_points, probes.axis_times, step, axis_mask
    )
    tail_mask, tail_rep = _fd_representability(probes.tail_points, step)
    tail_div, _ = _divergence_on_mask(
        field, probes.tail_points, probes.tail_times, step, tail_mask
    )

    return {
        "sampled_max": sampled_max,
        "pooled_weighted_rms": weighted["weighted_rms"],
        "pooled_volume_l2_estimate": weighted["volume_l2_estimate"],
        "pooled_log_volume_l2_estimate": weighted["log_volume_l2_estimate"],
        "estimated_volume": weighted["estimated_volume"],
        "log_estimated_volume": weighted["log_estimated_volume"],
        "normalized_sampled_max": normalized_max,
        "normalized_weighted_rms": norm_weighted["weighted_rms"],
        "speed_rms": speed_rms,
        "main_fd_representability": rep,
        "pulse_entry_seam_fd_representability": seam_rep,
        "axis_fd_representability": axis_rep,
        "late_tail_fd_representability": tail_rep,
        "pulse_entry_seam_sampled_max": (
            float(np.max(np.abs(seam_div))) if seam_div.size else None
        ),
        "axis_axis_near_sampled_max": (
            float(np.max(np.abs(axis_div))) if axis_div.size else None
        ),
        "late_tail_sampled_max": (
            float(np.max(np.abs(tail_div))) if tail_div.size else None
        ),
        "per_region": by_region,
        "per_time": by_time,
        "worst_witness": worst,
    }


class _DifferenceField:
    def __init__(self, total: Any, leading: Any):
        self.total = total
        self.leading = leading

    def velocity(self, x, y, z, t):
        return np.asarray(self.total.velocity(x, y, z, t), float) - np.asarray(
            self.leading.velocity(x, y, z, t), float
        )


class _XLinearMutation:
    def __init__(self, base: Any, epsilon: float = MUTATION_EPSILON):
        self.base = base
        self.epsilon = float(epsilon)

    def velocity(self, x, y, z, t):
        out = np.asarray(self.base.velocity(x, y, z, t), float).copy()
        out[..., 0] += self.epsilon * np.asarray(x, float)
        return out


def _stability_ok(medium: Any, fine: Any) -> bool:
    if medium is None or fine is None:
        return False
    return bool(
        float(fine)
        <= STABILITY_FACTOR * float(medium) + STABILITY_FLOOR
    )


def _numeric_leq(value: Any, threshold: float) -> bool:
    return (
        value is not None
        and math.isfinite(float(value))
        and float(value) <= threshold
    )


def _all_fd_representable(metrics: Mapping[str, Any]) -> bool:
    return all(
        bool(metrics[key]["all_probes_representable"])
        for key in (
            "main_fd_representability",
            "pulse_entry_seam_fd_representability",
            "axis_fd_representability",
            "late_tail_fd_representability",
        )
    )


def _mutate_numeric(obj: Any) -> bool:
    if isinstance(obj, dict):
        for key in sorted(obj):
            value = obj[key]
            if (
                not isinstance(value, bool)
                and isinstance(value, (int, float))
                and math.isfinite(float(value))
            ):
                obj[key] = float(value) + 1.0e-7
                return True
            if _mutate_numeric(value):
                return True
    elif isinstance(obj, list):
        for index, value in enumerate(obj):
            if (
                not isinstance(value, bool)
                and isinstance(value, (int, float))
                and math.isfinite(float(value))
            ):
                obj[index] = float(value) + 1.0e-7
                return True
            if _mutate_numeric(value):
                return True
    return False


def _save_json(path: Path, obj: Mapping[str, Any]) -> None:
    path.write_text(
        json.dumps(obj, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def _leading_mutation_rejected(
    saved: Mapping[str, Any], tmp: Path
) -> bool:
    raw = copy.deepcopy(saved)
    cfg = raw["configuration"]
    lead = cfg["agent1"]
    nested = lead["configuration"]
    nested["schema"] = str(nested.get("schema", "missing")) + "-a4-mutation"
    lead["configuration_sha256"] = _candidate._id._sha256(nested)
    lead["semantic_sha256"] = _candidate._id._sha256(nested)
    raw["semantic_sha256"] = _candidate._id._sha256(cfg)
    path = tmp / "mutated-leading.json"
    _save_json(path, raw)
    try:
        CurrentMainPulseLogXLeadingOscillatoryField.load_candidate(path)
    except (ValueError, RuntimeError):
        return True
    return False


def _osc_mutation_rejected(
    saved: Mapping[str, Any], tmp: Path
) -> bool:
    raw = copy.deepcopy(saved)
    cfg = raw["configuration"]
    osc = cfg["oscillatory_runtime"]
    if not _mutate_numeric(osc["payload"]):
        raise RuntimeError("no numeric oscillatory value to mutate")
    osc["payload_sha256"] = _candidate._id._sha256(osc["payload"])
    raw["semantic_sha256"] = _candidate._id._sha256(cfg)
    path = tmp / "mutated-oscillatory.json"
    _save_json(path, raw)
    try:
        CurrentMainPulseLogXLeadingOscillatoryField.load_candidate(path)
    except (ValueError, RuntimeError):
        return True
    return False


def _post_xi11_fails_closed(field: Any) -> bool:
    xi = SOURCE_XI_END + 0.02
    log_X = _log_X_from_xi(field, [xi])
    point = _source_to_cartesian_logx(
        field, log_X, [0.0], [0.71], [0.47]
    )
    try:
        _velocity(field, point, np.asarray([0.47]))
    except (ValueError, RuntimeError):
        return True
    return False


def _component_order_invariant(
    a: Mapping[str, Any], b: Mapping[str, Any]
) -> bool:
    for key in ("sampled_max", "pooled_weighted_rms"):
        av, bv = a.get(key), b.get(key)
        if av is None or bv is None:
            if av is not bv:
                return False
        elif abs(float(av) - float(bv)) > ORDER_INVARIANCE_ATOL:
            return False
    return True


class _ManufacturedSolenoidalField:
    def velocity(self, x, y, z, t):
        x, y, z, t = np.broadcast_arrays(
            np.asarray(x, float),
            np.asarray(y, float),
            np.asarray(z, float),
            np.asarray(t, float),
        )
        return np.stack(
            ((1.0 + t) * y, -(1.0 + t) * x, np.zeros_like(z)),
            axis=-1,
        )


def _manufactured_calibration() -> bool:
    points = np.asarray(
        ((0.2, -0.3, 0.1), (-0.4, 0.1, -0.2), (0.6, 0.2, 0.3))
    )
    times = np.asarray((0.31, 0.47, 0.63))
    jac = independent_fd2_jacobian(
        _ManufacturedSolenoidalField(),
        points,
        times,
        SPATIAL_STEPS[-1],
    )
    div = np.trace(jac, axis1=1, axis2=2)
    return bool(np.max(np.abs(div)) <= 1.0e-12)


def _mutation_detection(field: Any) -> tuple[float, bool]:
    points = _INNER_POINTS
    times = _INNER_TIMES
    jac0 = independent_fd2_jacobian(
        field, points, times, SPATIAL_STEPS[-1]
    )
    jac1 = independent_fd2_jacobian(
        _XLinearMutation(field), points, times, SPATIAL_STEPS[-1]
    )
    signal = float(
        np.max(
            np.abs(
                np.trace(jac1, axis1=1, axis2=2)
                - np.trace(jac0, axis1=1, axis2=2)
            )
        )
    )
    return signal, bool(signal >= MUTATION_DETECTION_FLOOR)


def _source_parameter_perturbation_report(field: Any) -> dict[str, Any]:
    # Small source-coordinate perturbations are evaluation probes only: they do
    # not modify candidate parameters. The report detects discontinuous or
    # nonfinite public evaluation around a fresh off-grid point.
    xi0, eta0, theta0, t0 = 1.37, 0.173, 0.83, 0.527
    eps = 1.0e-6
    values = []
    for deta in (-eps, 0.0, eps):
        point = _source_to_cartesian_logx(
            field,
            _log_X_from_xi(field, [xi0]),
            [eta0 + deta],
            [theta0],
            [t0],
        )
        values.append(_velocity(field, point, np.asarray([t0]))[0])
    arr = np.asarray(values)
    second = arr[2] - 2.0 * arr[1] + arr[0]
    return {
        "eta_perturbation": eps,
        "all_finite": bool(np.all(np.isfinite(arr))),
        "centered_second_difference_norm": float(np.linalg.norm(second)),
    }


def _derive_audit_pass(receipt: Mapping[str, Any]) -> bool:
    medium = receipt["resolutions"][1]["leading_plus_oscillatory"]
    fine = receipt["resolutions"][2]["leading_plus_oscillatory"]
    checks = receipt["checks"]
    keys = (
        "sampled_max",
        "pooled_weighted_rms",
        "pulse_entry_seam_sampled_max",
        "axis_axis_near_sampled_max",
        "late_tail_sampled_max",
    )
    return bool(
        all(_numeric_leq(fine.get(key), DIVERGENCE_GATE) for key in keys)
        and fine.get("speed_rms") is not None
        and float(fine["speed_rms"]) >= NONTRIVIAL_SPEED_RMS
        and float(receipt["inner_oscillatory_signal_rms"])
        >= OSCILLATORY_SIGNAL_RMS_MIN
        and all(
            _stability_ok(medium.get(key), fine.get(key))
            for key in keys
        )
        and _all_fd_representable(fine)
        and checks["all_resolutions_fd_representable"]
        and checks["save_reload_semantic_identity"]
        and checks["save_reload_velocity_replay_max_abs"]
        <= SAVE_RELOAD_VELOCITY_TOLERANCE
        and checks["velocity_mutation_detected"]
        and checks["offgrid_order_invariant_all_components"]
        and checks["rehash_consistent_leading_mutation_rejected"]
        and checks["rehash_consistent_oscillatory_mutation_rejected"]
        and checks["post_xi11_fails_closed"]
        and checks["manufactured_solenoidal_calibration_pass"]
        and checks["source_parameter_perturbation"]["all_finite"]
    )


def materialize_receipt() -> dict[str, Any]:
    original = default_field()
    with tempfile.TemporaryDirectory(prefix="k4-logx-xi11-") as directory:
        tmp = Path(directory)
        path = tmp / "candidate.json"
        saved = original.save_candidate(path)
        field = CurrentMainPulseLogXLeadingOscillatoryField.load_candidate(path)
        lead_mut = _leading_mutation_rejected(saved, tmp)
        osc_mut = _osc_mutation_rejected(saved, tmp)

    leading = _leading(field)
    oscillatory = _DifferenceField(field, leading)
    probes = make_heldout_probes(field)

    resolutions = []
    for h in SPATIAL_STEPS:
        resolutions.append(
            {
                "step": h,
                "leading": _metrics(leading, probes, h),
                "leading_plus_oscillatory": _metrics(field, probes, h),
                "oscillatory_increment": _metrics(oscillatory, probes, h),
            }
        )

    n = min(8, len(probes.points))
    replay = float(
        np.max(
            np.abs(
                _velocity(original, probes.points[:n], probes.times[:n])
                - _velocity(field, probes.points[:n], probes.times[:n])
            )
        )
    )
    inner_osc = _velocity(field, _INNER_POINTS, _INNER_TIMES) - _velocity(
        leading, _INNER_POINTS, _INNER_TIMES
    )
    inner_rms = float(np.sqrt(np.mean(np.sum(inner_osc**2, axis=1))))

    reverse = np.arange(len(probes.points) - 1, -1, -1)
    reversed_probes = HeldoutProbeSet(
        points=probes.points[reverse],
        times=probes.times[reverse],
        log_weights=probes.log_weights[reverse],
        region=probes.region[reverse],
        time_index=probes.time_index[reverse],
        seam_points=probes.seam_points,
        seam_times=probes.seam_times,
        axis_points=probes.axis_points,
        axis_times=probes.axis_times,
        tail_points=probes.tail_points,
        tail_times=probes.tail_times,
        tail_xi=probes.tail_xi,
    )
    reverse_fine = {
        "leading": _metrics(leading, reversed_probes, SPATIAL_STEPS[-1]),
        "leading_plus_oscillatory": _metrics(
            field, reversed_probes, SPATIAL_STEPS[-1]
        ),
        "oscillatory_increment": _metrics(
            oscillatory, reversed_probes, SPATIAL_STEPS[-1]
        ),
    }
    order_ok = all(
        _component_order_invariant(
            resolutions[-1][name], reverse_fine[name]
        )
        for name in (
            "leading",
            "leading_plus_oscillatory",
            "oscillatory_increment",
        )
    )

    mutation_signal, mutation_ok = _mutation_detection(field)
    all_rep = all(
        _all_fd_representable(res[stage])
        for res in resolutions
        for stage in (
            "leading",
            "leading_plus_oscillatory",
            "oscillatory_increment",
        )
    )

    receipt = {
        "schema": SCHEMA,
        "task": TASK,
        "upstream_pr": UPSTREAM_PR,
        "upstream_head": UPSTREAM_HEAD,
        "agent1_pr": AGENT1_PR,
        "agent1_head": AGENT1_HEAD,
        "newer_unconsumed_agent1_pr": NEWER_AGENT1_PR,
        "newer_unconsumed_agent1_head": NEWER_AGENT1_HEAD,
        "candidate_semantic_sha256": field.semantic_sha256,
        "oscillatory_runtime_sha256": field.oscillatory_runtime_sha256,
        "domain": {
            "finite_parent_xi_max": float(field.finite_xi_max),
            "source_main_xi_end": SOURCE_XI_END,
            "log_X_p": _log_X_p(field),
            "log_X_source_end": float(field.log_X_source_end),
            "full_source_xi_11_current_cartesian_leading_materialized": True,
            "terminal_global_leading_velocity_materialized": False,
        },
        "protocol": {
            "seed": SEED,
            "spatial_steps": list(SPATIAL_STEPS),
            "times": list(TIMES),
            "pre_pulse_xi_band": list(PRE_PULSE_XI_BAND),
            "main_pulse_xi_bands": [list(v) for v in XI_BANDS],
            "eta_interval": list(ETA_INTERVAL),
            "late_tail_xi": list(TAIL_XI),
            "integration_probe_count": len(probes.points),
            "pulse_entry_seam_probe_count": len(probes.seam_points),
            "axis_axis_near_probe_count": len(probes.axis_points),
            "late_tail_probe_count": len(probes.tail_points),
            "derivative_operator": (
                "centered Cartesian FD2 from public velocity only"
            ),
            "fd_perturbation_representability_required": True,
            "staged_components": [
                "leading",
                "leading_plus_oscillatory",
                "oscillatory_increment",
            ],
            "physical_volume_weighting": (
                "log-weighted independent fixed-time chart Jacobian "
                "with dX=X*dlogX"
            ),
        },
        "gates": {
            "total_finest_sampled_max": DIVERGENCE_GATE,
            "total_finest_weighted_rms": DIVERGENCE_GATE,
            "total_pulse_entry_seam_sampled_max": DIVERGENCE_GATE,
            "total_axis_axis_near_sampled_max": DIVERGENCE_GATE,
            "total_late_tail_sampled_max": DIVERGENCE_GATE,
            "total_nontrivial_speed_rms_min": NONTRIVIAL_SPEED_RMS,
            "inner_oscillatory_signal_rms_min": OSCILLATORY_SIGNAL_RMS_MIN,
            "stability_factor": STABILITY_FACTOR,
            "stability_floor": STABILITY_FLOOR,
            "mutation_detection_floor": MUTATION_DETECTION_FLOOR,
            "save_reload_velocity_tolerance": SAVE_RELOAD_VELOCITY_TOLERANCE,
            "final_project_momentum_gate_unchanged": 1.0e-3,
            "final_project_divergence_gate_unchanged": 1.0e-5,
        },
        "inner_oscillatory_signal_rms": inner_rms,
        "resolutions": resolutions,
        "checks": {
            "save_reload_semantic_identity": (
                field.semantic_sha256 == original.semantic_sha256
            ),
            "save_reload_velocity_replay_max_abs": replay,
            "offgrid_order_invariant_all_components": order_ok,
            "velocity_mutation_signal": mutation_signal,
            "velocity_mutation_detected": mutation_ok,
            "all_resolutions_fd_representable": all_rep,
            "rehash_consistent_leading_mutation_rejected": lead_mut,
            "rehash_consistent_oscillatory_mutation_rejected": osc_mut,
            "post_xi11_fails_closed": _post_xi11_fails_closed(field),
            "manufactured_solenoidal_calibration_pass": (
                _manufactured_calibration()
            ),
            "source_parameter_perturbation": (
                _source_parameter_perturbation_report(field)
            ),
        },
        "truth_boundary": copy.deepcopy(_TRUTH_BOUNDARY),
        "limitations": [
            (
                "scoped log-X xi<=11 divergence evidence, not a complete "
                "Navier-Stokes momentum residual"
            ),
            (
                "canonical absolute Cartesian FD steps are frozen; any "
                "binary64 perturbation collapse fails closed rather than "
                "accepting zero-by-roundoff derivatives"
            ),
            (
                "ordinary physical volume-L2 may exceed binary64 range; "
                "the receipt therefore also records logarithmic volume-L2"
            ),
            (
                "A1 #1133 pulse-end xi<=13 leading is newer but has no "
                "matching A2 composite and is not transferred"
            ),
            (
                "no Cartesian correction velocity, matched pressure, "
                "preregistered restricted forcing, or terminal/global field"
            ),
            "same-protocol ST006 comparison remains unavailable",
        ],
    }
    receipt["audit_pass"] = _derive_audit_pass(receipt)
    return receipt


def enforce_preregistered_gates(receipt: Mapping[str, Any]) -> None:
    if (
        receipt.get("schema") != SCHEMA
        or receipt.get("task") != TASK
        or receipt.get("upstream_head") != UPSTREAM_HEAD
        or receipt.get("agent1_head") != AGENT1_HEAD
    ):
        raise AssertionError("audit identity drifted")
    if receipt.get("truth_boundary") != _TRUTH_BOUNDARY:
        raise AssertionError("truth boundary drifted")
    if not _derive_audit_pass(receipt):
        fine = receipt["resolutions"][-1]["leading_plus_oscillatory"]
        raise AssertionError(
            "preregistered log-X xi11 divergence audit failed: "
            f"max={fine.get('sampled_max')!r}, "
            f"rms={fine.get('pooled_weighted_rms')!r}, "
            f"seam={fine.get('pulse_entry_seam_sampled_max')!r}, "
            f"axis={fine.get('axis_axis_near_sampled_max')!r}, "
            f"tail={fine.get('late_tail_sampled_max')!r}, "
            f"fd_all={_all_fd_representable(fine)}"
        )


def _main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("--enforce-existing", action="store_true")
    args = parser.parse_args()
    path = Path(args.output)
    if args.enforce_existing:
        enforce_preregistered_gates(
            json.loads(path.read_text(encoding="utf-8"))
        )
        return
    receipt = materialize_receipt()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    _main()
