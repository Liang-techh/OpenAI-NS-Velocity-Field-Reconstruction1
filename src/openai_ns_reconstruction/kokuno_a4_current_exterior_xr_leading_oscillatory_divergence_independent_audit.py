"""Independent staged Cartesian divergence audit for A2 #987 through X_R.

The scientific path save/reloads the identity-bound A2 composite and then calls only
public ``velocity(x,y,z,t)``.  For staged attribution it reconstructs the exact A1 #980
leading field from the leading configuration stored in that artifact and likewise calls
only public velocity.  No production Jacobian/divergence, source derivative, pressure,
forcing, training tensor, or residual implementation is reused.
"""
from __future__ import annotations

import copy
from dataclasses import dataclass
import hashlib
import importlib
import inspect
import json
import math
from pathlib import Path
import tempfile
from typing import Any, Mapping

import numpy as np

from . import kokuno_current_exterior_xr_leading_oscillatory_identity as _upstream

SCHEMA = "kokuno-a4-current-exterior-xr-leading-oscillatory-divergence-audit-v1"
UPSTREAM_PR = 987
UPSTREAM_HEAD = "4be2c9ee898c44dd1ad2217e90601b161fe81964"
UPSTREAM_SOURCE_BLOB = "93a3da5c3aed4be8310cc510edfab6af940b64dc"
AGENT1_PR = 980
AGENT1_HEAD = "d3c971f2c62e272333e124e532212d23cca4908d"
AGENT1_SOURCE_BLOB = "2ba28636a56b252fa485719e7b3e8da754d7e588"
PARENT_A4_PR = 983
PARENT_A4_HEAD = "5fdb59afbbe89525e43164f92f0856aad2949600"

SEED = 9173701
STEPS = (0.02, 0.01, 0.005)
TIMES = (0.31, 0.47, 0.63, 0.71)
INNER_RADIAL_INTERVAL = (0.30, 0.62)
INNER_AXIAL_INTERVAL = (-0.24, 0.24)
INNER_POINTS_PER_TIME = 8
ETA_INTERVAL = (-0.60, 0.60)
EXTERIOR_ZONES = (
    ("near_seam_exterior", -4.80, -4.10),
    ("middle_exterior", -3.45, -1.85),
    ("near_XR_exterior", -1.35, -0.35),
)
EXTERIOR_POINTS_PER_ZONE_TIME = 4
SEAM_OFFSETS = (-0.06, -0.02, 0.02, 0.06)
DIVERGENCE_GATE = 1.0e-5
STABILITY_FACTOR = 1.25
STABILITY_FLOOR = 2.0e-8
TOTAL_SPEED_FLOOR = 1.0e-10
OSCILLATORY_SPEED_FLOOR = 1.0e-12
MUTATION_EPSILON = 1.0e-3
MUTATION_DETECTION_FLOOR = 5.0e-4
REPLAY_ATOL = 2.0e-12
ORDER_ATOL = 2.0e-13

TRUTH_BOUNDARY = {
    "current_leading_plus_oscillatory_velocity_through_XR_consumed": True,
    "candidate_save_reload_required": True,
    "independent_cartesian_fd_operator_used": True,
    "leading_divergence_scoped_assessed": True,
    "leading_plus_oscillatory_divergence_scoped_assessed": True,
    "oscillatory_divergence_increment_scoped_assessed": True,
    "leading_only_ns_residual_assessed": False,
    "leading_plus_oscillatory_ns_residual_assessed": False,
    "after_correction_ns_residual_assessed": False,
    "post_XR_velocity_materialized": False,
    "outer_global_leading_velocity_materialized": False,
    "matched_global_pressure_materialized": False,
    "restricted_forcing_materialized": False,
    "agent3_correction_velocity_materialized": False,
    "canonical_24_48_96_volume_admission_assessed": False,
    "same_protocol_comparable_to_st006": False,
    "residual_reduction_claimed": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
}


def _canonical_json(payload: Mapping[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _git_blob_sha1(path: str | Path) -> str:
    raw = Path(path).read_bytes()
    return hashlib.sha1(f"blob {len(raw)}\0".encode("ascii") + raw).hexdigest()


def _latin_hypercube(rng: np.random.Generator, count: int) -> np.ndarray:
    out = (np.arange(count, dtype=float) + rng.random(count)) / float(count)
    rng.shuffle(out)
    return out


def _source_xyz(field: Any, X: Any, eta: Any, theta: Any, t: Any) -> np.ndarray:
    X, eta, theta, t = np.broadcast_arrays(
        np.asarray(X, float), np.asarray(eta, float), np.asarray(theta, float), np.asarray(t, float)
    )
    q = (1.0 - t) / (1.0 - eta * eta)
    radius = np.sqrt(2.0 * q * X)
    z = np.power(q, float(field.D)) * eta
    return np.stack((radius * np.cos(theta), radius * np.sin(theta), z), axis=-1)


def _source_volume_jacobian(field: Any, eta: Any, t: Any) -> np.ndarray:
    eta, t = np.broadcast_arrays(np.asarray(eta, float), np.asarray(t, float))
    D = float(field.D)
    q = (1.0 - t) / (1.0 - eta * eta)
    return np.power(q, D + 1.0) * (1.0 + 2.0 * D * eta * eta / (1.0 - eta * eta))


@dataclass(frozen=True)
class ProbeStratum:
    points: np.ndarray
    times: np.ndarray
    weights: np.ndarray
    group_index: np.ndarray


@dataclass(frozen=True)
class ProbeSet:
    inner: ProbeStratum
    exterior: ProbeStratum
    seam_points: np.ndarray
    seam_times: np.ndarray
    axis_points: np.ndarray
    axis_times: np.ndarray


def make_probes(field: Any) -> ProbeSet:
    rng = np.random.default_rng(SEED)
    ip, it, iw, ig = [], [], [], []
    r0, r1 = INNER_RADIAL_INTERVAL
    z0, z1 = INNER_AXIAL_INTERVAL
    measure = (r1 - r0) * (z1 - z0) * 2.0 * math.pi
    for time_id, time in enumerate(TIMES):
        r = r0 + (r1 - r0) * _latin_hypercube(rng, INNER_POINTS_PER_TIME)
        z = z0 + (z1 - z0) * _latin_hypercube(rng, INNER_POINTS_PER_TIME)
        th = 2.0 * math.pi * _latin_hypercube(rng, INNER_POINTS_PER_TIME)
        ip.append(np.stack((r * np.cos(th), r * np.sin(th), z), axis=-1))
        it.append(np.full(INNER_POINTS_PER_TIME, time))
        iw.append(measure * r / INNER_POINTS_PER_TIME)
        ig.append(np.full(INNER_POINTS_PER_TIME, time_id, dtype=int))

    ep, et, ew, eg = [], [], [], []
    eta0, eta1 = ETA_INTERVAL
    for time_id, time in enumerate(TIMES):
        for zone_id, (_, s0, s1) in enumerate(EXTERIOR_ZONES):
            n = EXTERIOR_POINTS_PER_ZONE_TIME
            s = s0 + (s1 - s0) * _latin_hypercube(rng, n)
            X = float(field.X_R) * np.exp(s)
            eta = eta0 + (eta1 - eta0) * _latin_hypercube(rng, n)
            th = 2.0 * math.pi * _latin_hypercube(rng, n)
            tv = np.full(n, time)
            ep.append(_source_xyz(field, X, eta, th, tv))
            ew.append((s1 - s0) * (eta1 - eta0) * 2.0 * math.pi * X * _source_volume_jacobian(field, eta, tv) / n)
            et.append(tv)
            eg.append(np.full(n, time_id * len(EXTERIOR_ZONES) + zone_id, dtype=int))

    seam, seam_t = [], []
    for time_id, time in enumerate(TIMES):
        for j, offset in enumerate(SEAM_OFFSETS):
            X = float(field.X_h) * math.exp(offset)
            eta = -0.45 + 0.30 * ((time_id + j) % 4)
            theta = 0.37 + 0.61 * time_id + 0.29 * j
            seam.append(_source_xyz(field, [X], [eta], [theta], [time])[0])
            seam_t.append(time)

    ar = np.asarray([0.0, 0.0, 1.0e-10, 1.0e-8, 1.0e-6])
    ath = np.asarray([0.0, 1.1, 0.7, 1.7, 2.3])
    az = np.asarray([-0.20, 0.0, 0.18, -0.12, 0.15])
    at = np.asarray([0.31, 0.47, 0.63, 0.71, 0.55])
    axis = np.stack((ar * np.cos(ath), ar * np.sin(ath), az), axis=-1)
    return ProbeSet(
        ProbeStratum(np.concatenate(ip), np.concatenate(it), np.concatenate(iw), np.concatenate(ig)),
        ProbeStratum(np.concatenate(ep), np.concatenate(et), np.concatenate(ew), np.concatenate(eg)),
        np.asarray(seam), np.asarray(seam_t), axis, at,
    )


def _velocity(field: Any, points: np.ndarray, times: np.ndarray) -> np.ndarray:
    points, times = np.asarray(points, float), np.asarray(times, float)
    out = np.asarray(field.velocity(points[:, 0], points[:, 1], points[:, 2], times), float)
    if out.shape != (len(points), 3) or not np.all(np.isfinite(out)):
        raise RuntimeError("bad public velocity")
    return out


def fd2_jacobian(field: Any, points: np.ndarray, times: np.ndarray, step: float) -> np.ndarray:
    points, times = np.asarray(points, float), np.asarray(times, float)
    h = float(step)
    if points.ndim != 2 or points.shape[1] != 3 or times.shape != (len(points),):
        raise ValueError("invalid point/time shape")
    if not math.isfinite(h) or h <= 0.0:
        raise ValueError("step must be positive")
    J = np.empty((len(points), 3, 3), float)
    for axis in range(3):
        plus, minus = points.copy(), points.copy()
        plus[:, axis] += h
        minus[:, axis] -= h
        J[:, :, axis] = (_velocity(field, plus, times) - _velocity(field, minus, times)) / (2.0 * h)
    if not np.all(np.isfinite(J)):
        raise RuntimeError("nonfinite independent Jacobian")
    return J


def _stats(J: np.ndarray, weights: np.ndarray | None = None) -> dict[str, Any]:
    div = np.trace(J, axis1=1, axis2=2)
    frob = np.linalg.norm(J, axis=(1, 2))
    norm = np.abs(div) / np.maximum(1.0, frob)
    worst = int(np.argmax(np.abs(div)))
    out = {
        "sampled_max": float(np.max(np.abs(div))),
        "rms": float(np.sqrt(np.mean(div * div))),
        "normalized_sampled_max": float(np.max(norm)),
        "normalized_rms": float(np.sqrt(np.mean(norm * norm))),
        "worst_index": worst,
        "worst_divergence": float(div[worst]),
        "worst_normalized_divergence": float(norm[worst]),
    }
    if weights is not None:
        w = np.asarray(weights, float)
        volume = float(np.sum(w))
        out.update({
            "estimated_volume": volume,
            "weighted_rms": float(np.sqrt(np.sum(w * div * div) / volume)),
            "volume_l2_estimate": float(np.sqrt(np.sum(w * div * div))),
            "normalized_weighted_rms": float(np.sqrt(np.sum(w * norm * norm) / volume)),
        })
    return out


def _staged(total_field: Any, leading_field: Any, points: np.ndarray, times: np.ndarray, step: float, weights=None) -> dict[str, Any]:
    jt = fd2_jacobian(total_field, points, times, step)
    jl = fd2_jacobian(leading_field, points, times, step)
    jo = jt - jl
    uv = _velocity(total_field, points, times)
    lv = _velocity(leading_field, points, times)
    result = {
        "total": _stats(jt, weights),
        "leading": _stats(jl, weights),
        "oscillatory_increment": _stats(jo, weights),
        "total_speed_rms": float(np.sqrt(np.mean(np.sum(uv * uv, axis=1)))),
        "oscillatory_speed_rms": float(np.sqrt(np.mean(np.sum((uv - lv) ** 2, axis=1)))),
        "oscillatory_speed_abs_max": float(np.max(np.linalg.norm(uv - lv, axis=1))),
    }
    return result


def metrics_for_step(total: Any, leading: Any, probes: ProbeSet, step: float) -> dict[str, Any]:
    inner = _staged(total, leading, probes.inner.points, probes.inner.times, step, probes.inner.weights)
    exterior = _staged(total, leading, probes.exterior.points, probes.exterior.times, step, probes.exterior.weights)
    seam = _staged(total, leading, probes.seam_points, probes.seam_times, step)
    axis = _staged(total, leading, probes.axis_points, probes.axis_times, step)
    wi = int(inner["total"].pop("worst_index"))
    wo = int(inner["oscillatory_increment"].pop("worst_index"))
    inner["leading"].pop("worst_index")
    for group in (exterior, seam, axis):
        for stage in ("total", "leading", "oscillatory_increment"):
            group[stage].pop("worst_index")
    return {
        "step": float(step),
        "inner": inner,
        "exterior": exterior,
        "seam": seam,
        "axis_axis_near": axis,
        "worst_inner_total_witness": {
            "point": [float(v) for v in probes.inner.points[wi]],
            "time": float(probes.inner.times[wi]),
            "divergence": float(inner["total"]["worst_divergence"]),
            "normalized_divergence": float(inner["total"]["worst_normalized_divergence"]),
        },
        "worst_inner_oscillatory_witness": {
            "point": [float(v) for v in probes.inner.points[wo]],
            "time": float(probes.inner.times[wo]),
            "divergence_increment": float(inner["oscillatory_increment"]["worst_divergence"]),
            "normalized_divergence": float(inner["oscillatory_increment"]["worst_normalized_divergence"]),
        },
    }


class XLinearMutation:
    def __init__(self, base: Any):
        self.base = base
    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        out = np.asarray(self.base.velocity(x, y, z, t), float).copy()
        out[..., 0] += MUTATION_EPSILON * np.asarray(x, float)
        return out


def _load_leading(field: Any) -> Any:
    lead = field.configuration().get("agent1_leading")
    if not isinstance(lead, Mapping) or lead.get("pr") != AGENT1_PR or lead.get("head") != AGENT1_HEAD:
        raise RuntimeError("saved composite Agent-1 identity drifted")
    if lead.get("source_blob_sha1") != AGENT1_SOURCE_BLOB:
        raise RuntimeError("saved composite Agent-1 source identity drifted")
    cfg = lead.get("configuration")
    if not isinstance(cfg, Mapping):
        raise RuntimeError("saved Agent-1 configuration missing")
    module = importlib.import_module(_upstream.AGENT1_MODULE)
    cls = getattr(module, _upstream.AGENT1_CLASS)
    obj = cls.from_configuration(cfg)
    if str(obj.semantic_sha256) != str(lead.get("semantic_sha256")):
        raise RuntimeError("reconstructed Agent-1 semantic identity drifted")
    return obj


def _upstream_truth_ok(field: Any) -> bool:
    truth = getattr(field, "truth_boundary", None)
    req = {
        "current_leading_plus_oscillatory_velocity_through_XR_materialized": True,
        "identity_preserving_through_XR_composite_save_load_available": True,
        "post_XR_velocity_materialized": False,
        "outer_global_leading_velocity_materialized": False,
        "agent3_correction_velocity_composed": False,
        "matched_global_pressure_materialized": False,
        "restricted_forcing_materialized": False,
        "heldout_ns_residual_assessed": False,
        "pde_validated": False,
    }
    return isinstance(truth, Mapping) and all(truth.get(k) is v for k, v in req.items())


def _stable(medium: float, fine: float) -> bool:
    return fine <= STABILITY_FACTOR * medium + STABILITY_FLOOR


def _metric_keys():
    for region in ("inner", "exterior"):
        for stage in ("total", "leading", "oscillatory_increment"):
            for name in ("sampled_max", "weighted_rms"):
                yield region, stage, name
    for region in ("seam", "axis_axis_near"):
        for stage in ("total", "leading", "oscillatory_increment"):
            yield region, stage, "sampled_max"


def _derive_audit_pass(receipt: Mapping[str, Any]) -> bool:
    medium, fine = receipt["resolutions"][1], receipt["resolutions"][2]
    gated = all(float(fine[r][s][n]) <= DIVERGENCE_GATE for r, s, n in _metric_keys())
    stable = all(_stable(float(medium[r][s][n]), float(fine[r][s][n])) for r, s, n in _metric_keys())
    return bool(
        gated and stable
        and float(fine["inner"]["total_speed_rms"]) >= TOTAL_SPEED_FLOOR
        and float(fine["inner"]["oscillatory_speed_rms"]) >= OSCILLATORY_SPEED_FLOOR
        and all(bool(v) for v in receipt["checks"].values())
    )


def _post_xr_fail_closed(field: Any) -> bool:
    p = _source_xyz(field, [float(field.X_R) * 1.001], [0.2], [0.7], [0.5])[0]
    try:
        field.velocity(p[0], p[1], p[2], 0.5)
    except ValueError:
        return True
    return False


def _leading_mutation_changes_identity(leading: Any) -> bool:
    cfg = copy.deepcopy(leading.configuration())
    current = cfg.get("current")
    if not isinstance(current, dict) or "eta_fd_step" not in current:
        return False
    current["eta_fd_step"] = float(current["eta_fd_step"]) * 1.001
    mutated = type(leading).from_configuration(cfg)
    return str(mutated.semantic_sha256) != str(leading.semantic_sha256)


def _osc_mutation_rejected(saved: Mapping[str, Any]) -> bool:
    mutated = copy.deepcopy(dict(saved))
    config = mutated["configuration"]
    osc = config["oscillatory_runtime"]
    osc["payload"]["parameters"]["mode_imaginary_ratio"] = float(osc["payload"]["parameters"]["mode_imaginary_ratio"]) + 0.125
    osc["payload_sha256"] = _upstream._identity._sha256(osc["payload"])
    mutated["semantic_sha256"] = _upstream._identity._sha256(config)
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "mutated.json"
        path.write_text(json.dumps(mutated, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        try:
            _upstream.CurrentExteriorXRLeadingOscillatoryField.load_candidate(path)
        except ValueError:
            return True
    return False


def materialize_receipt() -> dict[str, Any]:
    original = _upstream.default_field()
    source_path = inspect.getsourcefile(type(original))
    if source_path is None:
        raise RuntimeError("cannot locate A2 #987 source")
    source_blob = _git_blob_sha1(source_path)
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "candidate.json"
        saved = original.save_candidate(path)
        field = _upstream.CurrentExteriorXRLeadingOscillatoryField.load_candidate(path)
        leading = _load_leading(field)
        probes = make_probes(field)
        replay_p = np.concatenate((probes.inner.points[:8], probes.exterior.points[:8]))
        replay_t = np.concatenate((probes.inner.times[:8], probes.exterior.times[:8]))
        replay = float(np.max(np.abs(_velocity(original, replay_p, replay_t) - _velocity(field, replay_p, replay_t))))
        resolutions = [metrics_for_step(field, leading, probes, h) for h in STEPS]

        ri = np.arange(len(probes.inner.points) - 1, -1, -1)
        re = np.arange(len(probes.exterior.points) - 1, -1, -1)
        reverse = ProbeSet(
            ProbeStratum(probes.inner.points[ri], probes.inner.times[ri], probes.inner.weights[ri], probes.inner.group_index[ri]),
            ProbeStratum(probes.exterior.points[re], probes.exterior.times[re], probes.exterior.weights[re], probes.exterior.group_index[re]),
            probes.seam_points[::-1], probes.seam_times[::-1], probes.axis_points[::-1], probes.axis_times[::-1]
        )
        rev = metrics_for_step(field, leading, reverse, STEPS[-1])
        fine = resolutions[-1]
        order_ok = (
            abs(rev["inner"]["total"]["sampled_max"] - fine["inner"]["total"]["sampled_max"]) <= ORDER_ATOL
            and abs(rev["inner"]["total"]["weighted_rms"] - fine["inner"]["total"]["weighted_rms"]) <= ORDER_ATOL
            and abs(rev["exterior"]["total"]["weighted_rms"] - fine["exterior"]["total"]["weighted_rms"]) <= ORDER_ATOL
            and abs(rev["inner"]["oscillatory_increment"]["weighted_rms"] - fine["inner"]["oscillatory_increment"]["weighted_rms"]) <= ORDER_ATOL
        )
        mut = metrics_for_step(XLinearMutation(field), leading, probes, STEPS[-1])
        mutation_detected = max(mut["inner"]["total"]["sampled_max"], mut["exterior"]["total"]["sampled_max"]) >= MUTATION_DETECTION_FLOOR
        checks = {
            "upstream_source_blob_exact": source_blob == UPSTREAM_SOURCE_BLOB,
            "save_reload_semantic_identity": str(field.semantic_sha256) == str(original.semantic_sha256),
            "save_reload_velocity_replay": replay <= REPLAY_ATOL,
            "independent_leading_identity_matches_saved_composite": str(leading.semantic_sha256) == str(field.configuration()["agent1_leading"]["semantic_sha256"]),
            "upstream_truth_boundary_fail_closed": _upstream_truth_ok(field),
            "velocity_mutation_detected": bool(mutation_detected),
            "offgrid_order_invariant": bool(order_ok),
            "leading_parameter_perturbation_changes_semantic_identity": _leading_mutation_changes_identity(leading),
            "oscillatory_runtime_parameter_mutation_rejected": _osc_mutation_rejected(saved),
            "post_XR_fail_closed": _post_xr_fail_closed(field),
        }
        receipt = {
            "schema": SCHEMA,
            "upstream": {
                "pr": UPSTREAM_PR, "head": UPSTREAM_HEAD, "source_blob": source_blob,
                "expected_source_blob": UPSTREAM_SOURCE_BLOB, "semantic_sha256": str(field.semantic_sha256),
                "agent1_pr": AGENT1_PR, "agent1_head": AGENT1_HEAD, "agent1_source_blob": AGENT1_SOURCE_BLOB,
                "agent1_semantic_sha256": str(leading.semantic_sha256), "parent_a4_pr": PARENT_A4_PR, "parent_a4_head": PARENT_A4_HEAD,
            },
            "protocol": {
                "seed": SEED, "steps": list(STEPS), "times": list(TIMES),
                "inner_radial_interval": list(INNER_RADIAL_INTERVAL), "inner_axial_interval": list(INNER_AXIAL_INTERVAL),
                "inner_points_per_time": INNER_POINTS_PER_TIME, "eta_interval": list(ETA_INTERVAL),
                "exterior_zones": [list(v) for v in EXTERIOR_ZONES], "exterior_points_per_zone_time": EXTERIOR_POINTS_PER_ZONE_TIME,
                "seam_offsets": list(SEAM_OFFSETS), "inner_probe_count": len(probes.inner.points), "exterior_probe_count": len(probes.exterior.points),
                "seam_probe_count": len(probes.seam_points), "axis_probe_count": len(probes.axis_points),
                "operator": "centered Cartesian FD2 on save/reloaded public total and independently reconstructed public leading velocity only",
                "divergence_gate": DIVERGENCE_GATE, "stability_factor": STABILITY_FACTOR, "stability_floor": STABILITY_FLOOR,
                "total_speed_floor": TOTAL_SPEED_FLOOR, "oscillatory_speed_floor": OSCILLATORY_SPEED_FLOOR,
                "mutation_epsilon": MUTATION_EPSILON, "mutation_detection_floor": MUTATION_DETECTION_FLOOR, "replay_atol": REPLAY_ATOL,
                "scope_note": "scoped physical inner/exterior volume estimates, not canonical 24/48/96 whole-domain admission",
                "momentum_note": "matched Cartesian pressure/restricted forcing absent; complete NS momentum residual not assessed",
            },
            "resolutions": resolutions,
            "mutation": {"kind": "total u_x += 1e-3*x", "finest_inner_total_sampled_max": mut["inner"]["total"]["sampled_max"], "finest_exterior_total_sampled_max": mut["exterior"]["total"]["sampled_max"]},
            "checks": checks,
            "truth_boundary": copy.deepcopy(TRUTH_BOUNDARY),
        }
        receipt["audit_pass"] = _derive_audit_pass(receipt)
        receipt["receipt_sha256"] = hashlib.sha256(_canonical_json(receipt).encode("utf-8")).hexdigest()
        return receipt


def enforce_receipt(receipt: Mapping[str, Any]) -> None:
    if receipt.get("schema") != SCHEMA:
        raise ValueError("unexpected A4 receipt schema")
    protocol = receipt.get("protocol", {})
    expected = {
        "seed": SEED, "steps": list(STEPS), "times": list(TIMES),
        "inner_radial_interval": list(INNER_RADIAL_INTERVAL), "inner_axial_interval": list(INNER_AXIAL_INTERVAL),
        "inner_points_per_time": INNER_POINTS_PER_TIME, "eta_interval": list(ETA_INTERVAL),
        "exterior_zones": [list(v) for v in EXTERIOR_ZONES], "exterior_points_per_zone_time": EXTERIOR_POINTS_PER_ZONE_TIME,
        "seam_offsets": list(SEAM_OFFSETS), "divergence_gate": DIVERGENCE_GATE, "stability_factor": STABILITY_FACTOR,
        "stability_floor": STABILITY_FLOOR, "total_speed_floor": TOTAL_SPEED_FLOOR, "oscillatory_speed_floor": OSCILLATORY_SPEED_FLOOR,
        "mutation_epsilon": MUTATION_EPSILON, "mutation_detection_floor": MUTATION_DETECTION_FLOOR, "replay_atol": REPLAY_ATOL,
    }
    for key, value in expected.items():
        if protocol.get(key) != value:
            raise ValueError(f"frozen protocol drift at {key}")
    upstream = receipt.get("upstream", {})
    if upstream.get("head") != UPSTREAM_HEAD or upstream.get("expected_source_blob") != UPSTREAM_SOURCE_BLOB or upstream.get("source_blob") != UPSTREAM_SOURCE_BLOB:
        raise ValueError("upstream exact identity drift")
    if upstream.get("agent1_head") != AGENT1_HEAD or upstream.get("agent1_source_blob") != AGENT1_SOURCE_BLOB:
        raise ValueError("Agent-1 exact identity drift")
    for key, value in TRUTH_BOUNDARY.items():
        if receipt.get("truth_boundary", {}).get(key) is not value:
            raise ValueError(f"truth-boundary drift at {key}")
    derived = _derive_audit_pass(receipt)
    if bool(receipt.get("audit_pass")) != derived:
        raise ValueError("audit_pass laundering")
    payload = copy.deepcopy(dict(receipt))
    observed = payload.pop("receipt_sha256", None)
    if observed != hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest():
        raise ValueError("receipt checksum mismatch")
    if not derived:
        fine = receipt["resolutions"][-1]
        raise RuntimeError(
            "through-XR composite divergence audit failed: "
            f"inner_total_max={fine['inner']['total']['sampled_max']:.6e}, "
            f"inner_total_rms={fine['inner']['total']['weighted_rms']:.6e}, "
            f"inner_osc_max={fine['inner']['oscillatory_increment']['sampled_max']:.6e}, "
            f"exterior_total_max={fine['exterior']['total']['sampled_max']:.6e}, "
            f"seam_total={fine['seam']['total']['sampled_max']:.6e}, axis_total={fine['axis_axis_near']['total']['sampled_max']:.6e}"
        )


def _main() -> None:
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--verify-receipt", type=Path)
    args = parser.parse_args()
    if args.verify_receipt is not None:
        enforce_receipt(json.loads(args.verify_receipt.read_text(encoding="utf-8")))
        print("A4 through-XR leading+oscillatory divergence receipt passes frozen gates")
        return
    receipt = materialize_receipt()
    text = json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False)
    if args.output is None:
        print(text)
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + "\n", encoding="utf-8")
        print(text)


if __name__ == "__main__":
    _main()
