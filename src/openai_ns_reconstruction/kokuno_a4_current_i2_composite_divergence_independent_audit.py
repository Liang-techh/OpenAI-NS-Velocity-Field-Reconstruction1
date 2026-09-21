"""Independent divergence audit for exact Kokuno current-I2 composite #1071.

Scientific path: real save/reload, then public Cartesian velocity only.  All
Cartesian derivatives in this module are centered FD2 and independent of the
candidate construction/Jacobian path.  This is scoped divergence evidence, not
complete Navier--Stokes momentum validation.
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

from . import kokuno_current_i2_leading_oscillatory_identity as _candidate

CurrentI2LeadingOscillatoryField = _candidate.CurrentI2LeadingOscillatoryField
default_field = _candidate.default_field

SCHEMA = "kokuno-a4-current-i2-composite-divergence-audit-v1"
UPSTREAM_PR = 1071
UPSTREAM_HEAD = "48d69e37e78e7f7f0e4e9936f28ff7974719288d"
AGENT1_PR = 1061
AGENT1_HEAD = "9c2abbc69ba1ef8c8a9e8d1699a466a75d0293c3"
SEED = 9173821
SPATIAL_STEPS = (0.02, 0.01, 0.005)
TIMES = (0.31, 0.47, 0.63, 0.71)
I2_ZONE_FRACTIONS = ((0.16, 0.30), (0.39, 0.55), (0.68, 0.84))
PRE_I2_FRACTIONS = (-0.24, -0.06)
ETA_INTERVAL = (-0.45, 0.45)
DIVERGENCE_GATE = 1.0e-5
STABILITY_FACTOR = 1.25
STABILITY_FLOOR = 2.0e-8
NONTRIVIAL_SPEED_RMS = 1.0e-10
OSCILLATORY_SIGNAL_RMS_MIN = 1.0e-12
MUTATION_EPSILON = 1.0e-3
MUTATION_DETECTION_FLOOR = 5.0e-4
SAVE_RELOAD_VELOCITY_TOLERANCE = 2.0e-12
ORDER_INVARIANCE_ATOL = 1.0e-13

_INNER_POINTS = np.asarray([
    (0.32, 0.11, -0.20), (0.41, -0.17, -0.08), (-0.36, 0.24, 0.05),
    (-0.52, -0.16, 0.16), (0.58, 0.21, 0.22), (0.47, -0.31, -0.14),
], dtype=float)
_INNER_TIMES = np.asarray((0.31, 0.39, 0.47, 0.55, 0.63, 0.71), dtype=float)

_TRUTH_BOUNDARY = {
    "current_i2_leading_plus_oscillatory_velocity_consumed": True,
    "candidate_save_reload_required": True,
    "independent_cartesian_fd_operator_used": True,
    "staged_leading_total_increment_divergence_scoped_assessed": True,
    "public_float64_i2_delta_verified_by_this_audit": False,
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


def _source_to_cartesian(field: Any, X: Any, eta: Any, theta: Any, t: Any) -> np.ndarray:
    X, eta, theta, t = np.broadcast_arrays(
        np.asarray(X, float), np.asarray(eta, float), np.asarray(theta, float), np.asarray(t, float)
    )
    if np.any(X < 0.0) or np.any(np.abs(eta) >= 1.0) or np.any(t >= 1.0):
        raise ValueError("probe outside similarity chart")
    q = (1.0 - t) / (1.0 - eta * eta)
    r = np.sqrt(2.0 * q * X)
    z = np.power(q, float(field.D)) * eta
    return np.stack((r * np.cos(theta), r * np.sin(theta), z), axis=-1)


def _physical_volume_jacobian(field: Any, eta: Any, t: Any) -> np.ndarray:
    eta, t = np.broadcast_arrays(np.asarray(eta, float), np.asarray(t, float))
    q = (1.0 - t) / (1.0 - eta * eta)
    out = np.power(q, float(field.D) + 1.0) * (
        1.0 + 2.0 * float(field.D) * eta * eta / (1.0 - eta * eta)
    )
    if np.any(~np.isfinite(out)) or np.any(out <= 0.0):
        raise RuntimeError("invalid physical-volume Jacobian")
    return out


def _lhs(rng: np.random.Generator, n: int) -> np.ndarray:
    v = (np.arange(n, dtype=float) + rng.random(n)) / n
    rng.shuffle(v)
    return v


@dataclass(frozen=True)
class HeldoutProbeSet:
    points: np.ndarray
    times: np.ndarray
    weights: np.ndarray
    region: np.ndarray
    time_index: np.ndarray
    seam_points: np.ndarray
    seam_times: np.ndarray
    axis_points: np.ndarray
    axis_times: np.ndarray


def _region(field: Any, rng: np.random.Generator, time: float, tid: int,
            label: str, lo: float, hi: float, n: int = 4):
    uX, ue, ut = (_lhs(rng, n), _lhs(rng, n), _lhs(rng, n))
    logX = lo + (hi - lo) * uX
    X = np.exp(logX)
    eta = ETA_INTERVAL[0] + (ETA_INTERVAL[1] - ETA_INTERVAL[0]) * ue
    theta = 2.0 * math.pi * ut
    tt = np.full(n, time)
    xyz = _source_to_cartesian(field, X, eta, theta, tt)
    jac = _physical_volume_jacobian(field, eta, tt)
    cell = (hi - lo) * (ETA_INTERVAL[1] - ETA_INTERVAL[0]) * 2.0 * math.pi
    weight = cell * X * jac / n
    return xyz, tt, weight, np.full(n, label, object), np.full(n, tid, int)


def make_heldout_probes(field: Any) -> HeldoutProbeSet:
    rng = np.random.default_rng(SEED)
    s0, s1 = math.log(float(field.X_I2_start)), math.log(float(field.X_I2_end))
    span = s1 - s0
    if span <= 0.0:
        raise RuntimeError("invalid I2 interval")
    blocks = []
    for tid, time in enumerate(TIMES):
        blocks.append(_region(field, rng, time, tid, "pre_i2_control",
                              s0 + PRE_I2_FRACTIONS[0] * span,
                              s0 + PRE_I2_FRACTIONS[1] * span))
        for zid, (a, b) in enumerate(I2_ZONE_FRACTIONS):
            blocks.append(_region(field, rng, time, tid, f"i2_zone_{zid}",
                                  s0 + a * span, s0 + b * span))
    points = np.concatenate([b[0] for b in blocks])
    times = np.concatenate([b[1] for b in blocks])
    weights = np.concatenate([b[2] for b in blocks])
    regions = np.concatenate([b[3] for b in blocks])
    tids = np.concatenate([b[4] for b in blocks])

    seam_points, seam_times = [], []
    for time in (0.47, 0.63):
        off = np.asarray((-0.025, -0.008, 0.008, 0.025)) * span
        for eta0 in (-0.22, 0.22):
            X = np.exp(s0 + off)
            eta = np.full(4, eta0)
            th = np.asarray((0.37, 1.11, 2.03, 2.81))
            tt = np.full(4, time)
            seam_points.append(_source_to_cartesian(field, X, eta, th, tt))
            seam_times.append(tt)
    axis_X = np.asarray((0.0, 0.0, 1e-10, 1e-8, 1e-6))
    axis_eta = np.asarray((0.0, 0.40, -0.40, 0.22, 0.0))
    axis_th = np.asarray((0.0, 0.9, 1.4, 2.0, 2.7))
    axis_t = np.asarray((0.50, 0.47, 0.63, 0.31, 0.71))
    return HeldoutProbeSet(
        points, times, weights, regions, tids,
        np.concatenate(seam_points), np.concatenate(seam_times),
        _source_to_cartesian(field, axis_X, axis_eta, axis_th, axis_t), axis_t,
    )


def _velocity(field: Any, points: np.ndarray, times: np.ndarray) -> np.ndarray:
    p, tt = np.asarray(points, float), np.asarray(times, float)
    out = np.asarray(field.velocity(p[:, 0], p[:, 1], p[:, 2], tt), float)
    if out.shape != (len(p), 3) or np.any(~np.isfinite(out)):
        raise RuntimeError("invalid public velocity output")
    return out


def independent_fd2_jacobian(field: Any, points: np.ndarray, times: np.ndarray,
                             step: float) -> np.ndarray:
    p, tt = np.asarray(points, float), np.asarray(times, float)
    if p.ndim != 2 or p.shape[1] != 3 or tt.shape != (len(p),):
        raise ValueError("points/times shape mismatch")
    h = float(step)
    if not math.isfinite(h) or h <= 0.0:
        raise ValueError("invalid step")
    J = np.empty((len(p), 3, 3))
    for axis in range(3):
        plus, minus = p.copy(), p.copy()
        plus[:, axis] += h
        minus[:, axis] -= h
        J[:, :, axis] = (_velocity(field, plus, tt) - _velocity(field, minus, tt)) / (2.0 * h)
    return J


def _divergence(field: Any, points: np.ndarray, times: np.ndarray, step: float):
    J = independent_fd2_jacobian(field, points, times, step)
    return np.trace(J, axis1=1, axis2=2), J


def _weighted_metrics(div: np.ndarray, weights: np.ndarray):
    vol = float(np.sum(weights))
    if vol <= 0.0:
        raise RuntimeError("nonpositive probe volume")
    l2 = float(np.sqrt(np.sum(weights * div * div)))
    return vol, l2, l2 / math.sqrt(vol)


def _metrics(field: Any, probes: HeldoutProbeSet, h: float) -> dict[str, Any]:
    div, J = _divergence(field, probes.points, probes.times, h)
    frob = np.linalg.norm(J, axis=(1, 2))
    normed = np.abs(div) / np.maximum(1.0, frob)
    speed = np.linalg.norm(_velocity(field, probes.points, probes.times), axis=1)
    vol, l2, rms = _weighted_metrics(div, probes.weights)
    by_region = {}
    for label in sorted(set(map(str, probes.region.tolist()))):
        m = probes.region == label
        rv, rl2, rrms = _weighted_metrics(div[m], probes.weights[m])
        by_region[label] = {"estimated_volume": rv, "sampled_max": float(np.max(np.abs(div[m]))),
                            "volume_l2_estimate": rl2, "weighted_rms": rrms}
    by_time = []
    for tid, time in enumerate(TIMES):
        m = probes.time_index == tid
        tv, tl2, trms = _weighted_metrics(div[m], probes.weights[m])
        by_time.append({"time": time, "estimated_volume": tv,
                        "sampled_max": float(np.max(np.abs(div[m]))),
                        "volume_l2_estimate": tl2, "weighted_rms": trms})
    seam, _ = _divergence(field, probes.seam_points, probes.seam_times, h)
    axis, _ = _divergence(field, probes.axis_points, probes.axis_times, h)
    iw = int(np.argmax(np.abs(div)))
    return {
        "sampled_max": float(np.max(np.abs(div))), "estimated_volume": vol,
        "pooled_volume_l2_estimate": l2, "pooled_weighted_rms": rms,
        "normalized_sampled_max": float(np.max(normed)),
        "normalized_weighted_rms": float(np.sqrt(np.sum(probes.weights * normed**2) / vol)),
        "i2_entry_seam_sampled_max": float(np.max(np.abs(seam))),
        "axis_axis_near_sampled_max": float(np.max(np.abs(axis))),
        "speed_rms": float(np.sqrt(np.mean(speed**2))),
        "per_region": by_region, "per_time": by_time,
        "worst_witness": {"point": probes.points[iw].tolist(), "time": float(probes.times[iw]),
                          "region": str(probes.region[iw]), "divergence": float(div[iw]),
                          "jacobian_frobenius": float(frob[iw]),
                          "normalized_divergence": float(normed[iw])},
    }


class _DifferenceField:
    def __init__(self, total: Any, leading: Any):
        self.total, self.leading, self.D = total, leading, total.D
    def velocity(self, x, y, z, t):
        return np.asarray(self.total.velocity(x, y, z, t), float) - np.asarray(
            self.leading.velocity(x, y, z, t), float)


class _XLinearMutation:
    def __init__(self, base: Any, epsilon: float = MUTATION_EPSILON):
        self.base, self.epsilon, self.D = base, float(epsilon), base.D
    def velocity(self, x, y, z, t):
        out = np.asarray(self.base.velocity(x, y, z, t), float).copy()
        out[..., 0] += self.epsilon * np.asarray(x, float)
        return out


def _stability_ok(medium: float, fine: float) -> bool:
    return bool(fine <= STABILITY_FACTOR * medium + STABILITY_FLOOR)


def _mutate_numeric(obj: Any) -> bool:
    if isinstance(obj, dict):
        for k in sorted(obj):
            v = obj[k]
            if not isinstance(v, bool) and isinstance(v, (int, float)) and math.isfinite(float(v)):
                obj[k] = float(v) + 1e-7
                return True
            if _mutate_numeric(v):
                return True
    if isinstance(obj, list):
        for i, v in enumerate(obj):
            if not isinstance(v, bool) and isinstance(v, (int, float)) and math.isfinite(float(v)):
                obj[i] = float(v) + 1e-7
                return True
            if _mutate_numeric(v):
                return True
    return False


def _save_json(path: Path, obj: Mapping[str, Any]):
    path.write_text(json.dumps(obj, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _leading_mutation_rejected(saved: Mapping[str, Any], tmp: Path) -> bool:
    raw = copy.deepcopy(saved)
    cfg = raw["configuration"]
    lead = cfg["agent1_leading"]
    nested = lead["configuration"]
    nested["schema"] = str(nested.get("schema", "missing")) + "-mutation"
    lead["configuration_sha256"] = _candidate._identity._sha256(nested)
    lead["semantic_sha256"] = _candidate._identity._sha256(nested)
    raw["semantic_sha256"] = _candidate._identity._sha256(cfg)
    path = tmp / "mutated-leading.json"
    _save_json(path, raw)
    try:
        CurrentI2LeadingOscillatoryField.load_candidate(path)
    except (ValueError, RuntimeError):
        return True
    return False


def _osc_mutation_rejected(saved: Mapping[str, Any], tmp: Path) -> bool:
    raw = copy.deepcopy(saved)
    cfg = raw["configuration"]
    osc = cfg["oscillatory_runtime"]
    if not _mutate_numeric(osc["payload"]):
        raise RuntimeError("no numeric oscillatory value to mutate")
    osc["payload_sha256"] = _candidate._identity._sha256(osc["payload"])
    raw["semantic_sha256"] = _candidate._identity._sha256(cfg)
    path = tmp / "mutated-osc.json"
    _save_json(path, raw)
    try:
        CurrentI2LeadingOscillatoryField.load_candidate(path)
    except (ValueError, RuntimeError):
        return True
    return False


def _post_i2_fails_closed(field: Any) -> bool:
    p = _source_to_cartesian(field, [1.001 * float(field.X_I2_end)], [0.0], [0.71], [0.47])
    try:
        _velocity(field, p, np.asarray([0.47]))
    except (ValueError, RuntimeError):
        return True
    return False


def _component_order_invariant(a: Mapping[str, Any], b: Mapping[str, Any]) -> bool:
    return all(abs(float(a[k]) - float(b[k])) <= ORDER_INVARIANCE_ATOL
               for k in ("sampled_max", "pooled_weighted_rms", "pooled_volume_l2_estimate"))


def _derive_audit_pass(receipt: Mapping[str, Any]) -> bool:
    med = receipt["resolutions"][1]["leading_plus_oscillatory"]
    fine = receipt["resolutions"][2]["leading_plus_oscillatory"]
    c = receipt["checks"]
    keys = ("sampled_max", "pooled_weighted_rms", "i2_entry_seam_sampled_max",
            "axis_axis_near_sampled_max")
    return bool(
        all(fine[k] <= DIVERGENCE_GATE for k in keys)
        and fine["speed_rms"] >= NONTRIVIAL_SPEED_RMS
        and receipt["inner_oscillatory_signal_rms"] >= OSCILLATORY_SIGNAL_RMS_MIN
        and all(_stability_ok(med[k], fine[k]) for k in keys)
        and c["save_reload_semantic_identity"]
        and c["save_reload_velocity_replay_max_abs"] <= SAVE_RELOAD_VELOCITY_TOLERANCE
        and c["velocity_mutation_detected"] and c["offgrid_order_invariant_all_components"]
        and c["rehash_consistent_leading_mutation_rejected"]
        and c["rehash_consistent_oscillatory_mutation_rejected"]
        and c["post_i2_fails_closed"] and c["manufactured_solenoidal_calibration_pass"]
    )


class _ManufacturedSolenoidalField:
    D = 1.0
    def velocity(self, x, y, z, t):
        x, y, z, t = np.broadcast_arrays(np.asarray(x, float), np.asarray(y, float),
                                          np.asarray(z, float), np.asarray(t, float))
        return np.stack(((1 + t) * y, -(1 + t) * x, np.zeros_like(z)), axis=-1)


def _manufactured_calibration() -> bool:
    p = np.asarray(((0.2, -0.3, 0.1), (-0.4, 0.1, -0.2), (0.6, 0.2, 0.3)))
    tt = np.asarray((0.31, 0.47, 0.63))
    div, _ = _divergence(_ManufacturedSolenoidalField(), p, tt, SPATIAL_STEPS[-1])
    return bool(np.max(np.abs(div)) <= 1e-12)


def materialize_receipt() -> dict[str, Any]:
    original = default_field()
    with tempfile.TemporaryDirectory(prefix="k4-i2-") as d:
        tmp = Path(d)
        path = tmp / "candidate.json"
        saved = original.save_candidate(path)
        field = CurrentI2LeadingOscillatoryField.load_candidate(path)
        lead_mut = _leading_mutation_rejected(saved, tmp)
        osc_mut = _osc_mutation_rejected(saved, tmp)

    leading = field.leading_backend
    osc = _DifferenceField(field, leading)
    probes = make_heldout_probes(field)
    resolutions = []
    for h in SPATIAL_STEPS:
        resolutions.append({"step": h, "leading": _metrics(leading, probes, h),
                            "leading_plus_oscillatory": _metrics(field, probes, h),
                            "oscillatory_increment": _metrics(osc, probes, h)})

    n = min(8, len(probes.points))
    replay = float(np.max(np.abs(
        _velocity(original, probes.points[:n], probes.times[:n]) -
        _velocity(field, probes.points[:n], probes.times[:n]))))
    inner_osc = _velocity(field, _INNER_POINTS, _INNER_TIMES) - _velocity(
        leading, _INNER_POINTS, _INNER_TIMES)
    inner_rms = float(np.sqrt(np.mean(np.sum(inner_osc**2, axis=1))))

    rev = np.arange(len(probes.points) - 1, -1, -1)
    rp = HeldoutProbeSet(probes.points[rev], probes.times[rev], probes.weights[rev],
                         probes.region[rev], probes.time_index[rev], probes.seam_points,
                         probes.seam_times, probes.axis_points, probes.axis_times)
    reverse_fine = {"leading": _metrics(leading, rp, SPATIAL_STEPS[-1]),
                    "leading_plus_oscillatory": _metrics(field, rp, SPATIAL_STEPS[-1]),
                    "oscillatory_increment": _metrics(osc, rp, SPATIAL_STEPS[-1])}
    order_ok = all(_component_order_invariant(resolutions[-1][name], reverse_fine[name])
                   for name in ("leading", "leading_plus_oscillatory", "oscillatory_increment"))
    base_div, _ = _divergence(field, probes.points, probes.times, SPATIAL_STEPS[-1])
    mut_div, _ = _divergence(_XLinearMutation(field), probes.points, probes.times,
                             SPATIAL_STEPS[-1])
    mut_signal = float(np.max(np.abs(mut_div - base_div)))

    receipt = {
        "schema": SCHEMA, "upstream_pr": UPSTREAM_PR, "upstream_head": UPSTREAM_HEAD,
        "agent1_pr": AGENT1_PR, "agent1_head": AGENT1_HEAD,
        "candidate_semantic_sha256": field.semantic_sha256,
        "oscillatory_runtime_sha256": field.oscillatory_runtime_sha256,
        "protocol": {"seed": SEED, "spatial_steps": list(SPATIAL_STEPS), "times": list(TIMES),
                     "pre_i2_fractions_of_log_span": list(PRE_I2_FRACTIONS),
                     "i2_zone_fractions": [list(v) for v in I2_ZONE_FRACTIONS],
                     "eta_interval": list(ETA_INTERVAL), "integration_probe_count": len(probes.points),
                     "i2_entry_seam_probe_count": len(probes.seam_points),
                     "axis_axis_near_probe_count": len(probes.axis_points),
                     "derivative_operator": "centered Cartesian FD2 from public velocity only",
                     "staged_components": ["leading", "leading_plus_oscillatory", "oscillatory_increment"],
                     "physical_volume_weighting": "independent fixed-time chart Jacobian with dX=X*dlogX"},
        "gates": {"total_finest_sampled_max": DIVERGENCE_GATE,
                  "total_finest_weighted_rms": DIVERGENCE_GATE,
                  "total_i2_entry_seam_sampled_max": DIVERGENCE_GATE,
                  "total_axis_axis_near_sampled_max": DIVERGENCE_GATE,
                  "total_nontrivial_speed_rms_min": NONTRIVIAL_SPEED_RMS,
                  "inner_oscillatory_signal_rms_min": OSCILLATORY_SIGNAL_RMS_MIN,
                  "stability_factor": STABILITY_FACTOR, "stability_floor": STABILITY_FLOOR,
                  "mutation_detection_floor": MUTATION_DETECTION_FLOOR,
                  "save_reload_velocity_tolerance": SAVE_RELOAD_VELOCITY_TOLERANCE,
                  "final_project_momentum_gate_unchanged": 1e-3,
                  "final_project_divergence_gate_unchanged": 1e-5},
        "inner_oscillatory_signal_rms": inner_rms, "resolutions": resolutions,
        "checks": {"save_reload_semantic_identity": field.semantic_sha256 == original.semantic_sha256,
                   "save_reload_velocity_replay_max_abs": replay,
                   "offgrid_order_invariant_all_components": order_ok,
                   "velocity_mutation_signal": mut_signal,
                   "velocity_mutation_detected": mut_signal >= MUTATION_DETECTION_FLOOR,
                   "rehash_consistent_leading_mutation_rejected": lead_mut,
                   "rehash_consistent_oscillatory_mutation_rejected": osc_mut,
                   "post_i2_fails_closed": _post_i2_fails_closed(field),
                   "manufactured_solenoidal_calibration_pass": _manufactured_calibration()},
        "truth_boundary": copy.deepcopy(_TRUTH_BOUNDARY),
        "limitations": [
            "scoped current-I2 divergence evidence, not complete NS momentum residual",
            "high-precision I2 repair may be sub-ULP at public float64 boundary; no float64-delta promotion",
            "no current-lineage Cartesian correction velocity, matched pressure, or restricted forcing",
            "I3/I4/global completion and same-protocol ST006 comparison remain unavailable",
        ],
    }
    receipt["audit_pass"] = _derive_audit_pass(receipt)
    return receipt


def enforce_preregistered_gates(receipt: Mapping[str, Any]) -> None:
    if receipt.get("schema") != SCHEMA or receipt.get("upstream_head") != UPSTREAM_HEAD:
        raise AssertionError("audit identity drifted")
    if receipt.get("truth_boundary") != _TRUTH_BOUNDARY:
        raise AssertionError("truth boundary drifted")
    if not _derive_audit_pass(receipt):
        f = receipt["resolutions"][-1]["leading_plus_oscillatory"]
        raise AssertionError("preregistered current-I2 divergence audit failed: "
                             f"max={f['sampled_max']:.6e}, rms={f['pooled_weighted_rms']:.6e}, "
                             f"seam={f['i2_entry_seam_sampled_max']:.6e}, "
                             f"axis={f['axis_axis_near_sampled_max']:.6e}")


def _main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--output", required=True)
    p.add_argument("--enforce-existing", action="store_true")
    a = p.parse_args()
    path = Path(a.output)
    if a.enforce_existing:
        enforce_preregistered_gates(json.loads(path.read_text(encoding="utf-8")))
        return
    receipt = materialize_receipt()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    _main()
