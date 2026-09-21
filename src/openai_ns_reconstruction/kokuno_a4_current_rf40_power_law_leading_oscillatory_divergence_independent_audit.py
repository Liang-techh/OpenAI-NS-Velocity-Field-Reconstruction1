"""Independent Cartesian divergence audit for the current Kokuno RF40 power-law composite.

K4-VAL-102 consumes exact Agent-2 PR #1010.  The scientific path performs a
real candidate save/reload and then differentiates only the reloaded public
``velocity(x,y,z,t)``.  For staged attribution it reconstructs exact Agent-1
PR #1005 from the configuration carried by that saved artifact and likewise
uses only its public Cartesian velocity.

No Agent-1 source-coordinate derivative/primitive helper, Agent-2 production
Jacobian/divergence/vorticity, training tensor, pressure/forcing helper, or
candidate-side residual routine is used by the independent numerical operator.
This is preregistered scoped incompressibility evidence on ``X_3 < X < X_4``;
it is not complete Navier--Stokes momentum validation.
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

from . import kokuno_current_rf40_power_law_leading_oscillatory_identity as _upstream

SCHEMA = "kokuno-a4-current-rf40-power-law-leading-oscillatory-divergence-audit-v1"
TASK = "K4-VAL-102"

UPSTREAM_PR = 1010
UPSTREAM_HEAD = "e36d9da4b7f037e998f5b1658f8c0ea291a76b80"
UPSTREAM_TASK = "K2-OSC-087"
UPSTREAM_SCHEMA = "kokuno-a2-current-rf40-power-law-leading-oscillatory-identity-v1"

AGENT1_PR = 1005
AGENT1_HEAD = "2c76ebdc41d6c566f43a1305034ba2ff9dce410b"
AGENT1_SOURCE_BLOB = "36a0183352a2005bd0b5f477594e05a2e4fe3868"
PARENT_A4_PR = 1007
PARENT_A4_HEAD = "8b5b337940edd568bf30ddba7afe68815bce1fd9"

SEED = 9173741
STEPS = (0.02, 0.01, 0.005)
TIMES = (0.31, 0.47, 0.63, 0.71)
INNER_RADIAL_INTERVAL = (0.30, 0.62)
INNER_AXIAL_INTERVAL = (-0.24, 0.24)
INNER_POINTS_PER_TIME = 8
ETA_INTERVAL = (-0.60, 0.60)
POWER_LAW_ZONES = (
    ("power_law_early", 0.08, 0.28),
    ("power_law_middle", 0.38, 0.62),
    ("power_law_late", 0.72, 0.92),
)
POWER_LAW_POINTS_PER_ZONE_TIME = 6
X3_SEAM_LOG_OFFSETS = (-0.06, -0.02, 0.02, 0.06)

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
    "current_leading_plus_oscillatory_velocity_through_RF40_power_law_consumed": True,
    "candidate_save_reload_required": True,
    "independent_cartesian_fd_operator_used": True,
    "RF40_power_law_scoped_divergence_assessed": True,
    "leading_divergence_scoped_assessed": True,
    "leading_plus_oscillatory_divergence_scoped_assessed": True,
    "oscillatory_divergence_increment_scoped_assessed": True,
    "leading_only_ns_residual_assessed": False,
    "leading_plus_oscillatory_ns_residual_assessed": False,
    "after_correction_ns_residual_assessed": False,
    "velocity_after_RF40_power_law_materialized": False,
    "current_lineage_cone_I1_I2_I3_I4_outer_overlays_completed": False,
    "outer_global_leading_velocity_materialized": False,
    "global_compact_support_completed": False,
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


def _sha256(payload: Mapping[str, Any]) -> str:
    return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()


def _git_blob_sha1(path: str | Path) -> str:
    raw = Path(path).read_bytes()
    return hashlib.sha1(f"blob {len(raw)}\0".encode("ascii") + raw).hexdigest()


def _latin_hypercube(rng: np.random.Generator, count: int) -> np.ndarray:
    out = (np.arange(count, dtype=float) + rng.random(count)) / float(count)
    rng.shuffle(out)
    return out


def _source_xyz(field: Any, X: Any, eta: Any, theta: Any, t: Any) -> np.ndarray:
    X, eta, theta, t = np.broadcast_arrays(
        np.asarray(X, float),
        np.asarray(eta, float),
        np.asarray(theta, float),
        np.asarray(t, float),
    )
    if not all(np.all(np.isfinite(a)) for a in (X, eta, theta, t)):
        raise ValueError("source probes must be finite")
    if np.any(X < 0.0) or np.any(np.abs(eta) >= 1.0) or np.any((t <= 0.0) | (t >= 1.0)):
        raise ValueError("source probes escaped the registered chart")
    q = (1.0 - t) / (1.0 - eta * eta)
    radius = np.sqrt(2.0 * q * X)
    z = np.power(q, float(field.D)) * eta
    return np.stack((radius * np.cos(theta), radius * np.sin(theta), z), axis=-1)


def _source_volume_jacobian(field: Any, eta: Any, t: Any) -> np.ndarray:
    eta, t = np.broadcast_arrays(np.asarray(eta, float), np.asarray(t, float))
    D = float(field.D)
    q = (1.0 - t) / (1.0 - eta * eta)
    return np.power(q, D + 1.0) * (
        1.0 + 2.0 * D * eta * eta / (1.0 - eta * eta)
    )


@dataclass(frozen=True)
class ProbeStratum:
    points: np.ndarray
    times: np.ndarray
    weights: np.ndarray
    group_index: np.ndarray


@dataclass(frozen=True)
class ProbeSet:
    inner: ProbeStratum
    power_law: ProbeStratum
    seam_points: np.ndarray
    seam_times: np.ndarray
    axis_points: np.ndarray
    axis_times: np.ndarray


def make_probes(field: Any) -> ProbeSet:
    """Build preregistered fresh held-out probes without construction tensors."""
    rng = np.random.default_rng(SEED)

    ip, it, iw, ig = [], [], [], []
    r0, r1 = INNER_RADIAL_INTERVAL
    z0, z1 = INNER_AXIAL_INTERVAL
    inner_measure = (r1 - r0) * (z1 - z0) * 2.0 * math.pi
    for time_id, time in enumerate(TIMES):
        r = r0 + (r1 - r0) * _latin_hypercube(rng, INNER_POINTS_PER_TIME)
        z = z0 + (z1 - z0) * _latin_hypercube(rng, INNER_POINTS_PER_TIME)
        theta = 2.0 * math.pi * _latin_hypercube(rng, INNER_POINTS_PER_TIME)
        ip.append(np.stack((r * np.cos(theta), r * np.sin(theta), z), axis=-1))
        it.append(np.full(INNER_POINTS_PER_TIME, time))
        iw.append(inner_measure * r / INNER_POINTS_PER_TIME)
        ig.append(np.full(INNER_POINTS_PER_TIME, time_id, dtype=int))

    X3, X4, Tw = float(field.X_3), float(field.X_4), float(field.T_w)
    if not (X3 > 0.0 and X4 > X3 and Tw > 0.0):
        raise RuntimeError("RF40 power-law geometry is invalid")
    if not math.isclose(math.log(X4 / X3), Tw, rel_tol=0.0, abs_tol=2.0e-12):
        raise RuntimeError("RF40 power-law logarithmic length drifted")

    pp, pt, pw, pg = [], [], [], []
    eta0, eta1 = ETA_INTERVAL
    for time_id, time in enumerate(TIMES):
        for zone_id, (_, rho0, rho1) in enumerate(POWER_LAW_ZONES):
            n = POWER_LAW_POINTS_PER_ZONE_TIME
            rho = rho0 + (rho1 - rho0) * _latin_hypercube(rng, n)
            X = X3 * np.exp(Tw * rho)
            if not np.all((X > X3) & (X < X4)):
                raise RuntimeError("power-law sampler escaped X_3 < X < X_4")
            eta = eta0 + (eta1 - eta0) * _latin_hypercube(rng, n)
            theta = 2.0 * math.pi * _latin_hypercube(rng, n)
            tv = np.full(n, time)
            pp.append(_source_xyz(field, X, eta, theta, tv))
            # Uniform rho = log(X/X3)/T_w, hence dX = X*T_w*d rho.
            pw.append(
                (rho1 - rho0)
                * (eta1 - eta0)
                * 2.0
                * math.pi
                * X
                * Tw
                * _source_volume_jacobian(field, eta, tv)
                / n
            )
            pt.append(tv)
            pg.append(
                np.full(
                    n,
                    time_id * len(POWER_LAW_ZONES) + zone_id,
                    dtype=int,
                )
            )

    seam, seam_t = [], []
    for time_id, time in enumerate(TIMES):
        for j, offset in enumerate(X3_SEAM_LOG_OFFSETS):
            X = X3 * math.exp(offset)
            if offset > 0.0 and not X < X4:
                raise RuntimeError("positive X3 seam probe escaped power-law stage")
            eta = -0.45 + 0.30 * ((time_id + j) % 4)
            theta = 0.37 + 0.61 * time_id + 0.29 * j
            seam.append(_source_xyz(field, [X], [eta], [theta], [time])[0])
            seam_t.append(time)

    ar = np.asarray([0.0, 0.0, 1.0e-10, 1.0e-8, 1.0e-6])
    atheta = np.asarray([0.0, 1.1, 0.7, 1.7, 2.3])
    az = np.asarray([-0.20, 0.0, 0.18, -0.12, 0.15])
    axis_t = np.asarray([0.31, 0.47, 0.63, 0.71, 0.55])
    axis = np.stack((ar * np.cos(atheta), ar * np.sin(atheta), az), axis=-1)

    return ProbeSet(
        ProbeStratum(
            np.concatenate(ip),
            np.concatenate(it),
            np.concatenate(iw),
            np.concatenate(ig),
        ),
        ProbeStratum(
            np.concatenate(pp),
            np.concatenate(pt),
            np.concatenate(pw),
            np.concatenate(pg),
        ),
        np.asarray(seam),
        np.asarray(seam_t),
        axis,
        axis_t,
    )


def _velocity(field: Any, points: np.ndarray, times: np.ndarray) -> np.ndarray:
    points, times = np.asarray(points, float), np.asarray(times, float)
    out = np.asarray(
        field.velocity(points[:, 0], points[:, 1], points[:, 2], times),
        float,
    )
    if out.shape != (len(points), 3) or not np.all(np.isfinite(out)):
        raise RuntimeError("bad public velocity")
    return out


def fd2_jacobian(
    field: Any,
    points: np.ndarray,
    times: np.ndarray,
    step: float,
) -> np.ndarray:
    """Implementation-distinct centered Cartesian FD2 Jacobian."""
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
        J[:, :, axis] = (
            _velocity(field, plus, times) - _velocity(field, minus, times)
        ) / (2.0 * h)
    if not np.all(np.isfinite(J)):
        raise RuntimeError("nonfinite independent Jacobian")
    return J


def _stats(J: np.ndarray, weights: np.ndarray | None = None) -> dict[str, Any]:
    div = np.trace(J, axis1=1, axis2=2)
    frob = np.linalg.norm(J, axis=(1, 2))
    normalized = np.abs(div) / np.maximum(1.0, frob)
    worst = int(np.argmax(np.abs(div)))
    out: dict[str, Any] = {
        "sampled_max": float(np.max(np.abs(div))),
        "rms": float(np.sqrt(np.mean(div * div))),
        "normalized_sampled_max": float(np.max(normalized)),
        "normalized_rms": float(np.sqrt(np.mean(normalized * normalized))),
        "worst_index": worst,
        "worst_divergence": float(div[worst]),
        "worst_normalized_divergence": float(normalized[worst]),
    }
    if weights is not None:
        w = np.asarray(weights, float)
        volume = float(np.sum(w))
        if not (volume > 0.0 and np.all(w > 0.0) and np.all(np.isfinite(w))):
            raise RuntimeError("invalid physical-volume weights")
        out.update(
            {
                "estimated_volume": volume,
                "weighted_rms": float(np.sqrt(np.sum(w * div * div) / volume)),
                "volume_l2_estimate": float(np.sqrt(np.sum(w * div * div))),
                "normalized_weighted_rms": float(
                    np.sqrt(np.sum(w * normalized * normalized) / volume)
                ),
            }
        )
    return out


def _staged(
    total: Any,
    leading: Any,
    points: np.ndarray,
    times: np.ndarray,
    step: float,
    weights: np.ndarray | None = None,
) -> tuple[dict[str, Any], tuple[np.ndarray, np.ndarray, np.ndarray]]:
    jt = fd2_jacobian(total, points, times, step)
    jl = fd2_jacobian(leading, points, times, step)
    jo = jt - jl
    uv, lv = _velocity(total, points, times), _velocity(leading, points, times)
    return (
        {
            "total": _stats(jt, weights),
            "leading": _stats(jl, weights),
            "oscillatory_increment": _stats(jo, weights),
            "total_speed_rms": float(
                np.sqrt(np.mean(np.sum(uv * uv, axis=1)))
            ),
            "oscillatory_speed_rms": float(
                np.sqrt(np.mean(np.sum((uv - lv) ** 2, axis=1)))
            ),
            "oscillatory_speed_abs_max": float(
                np.max(np.linalg.norm(uv - lv, axis=1))
            ),
        },
        (jt, jl, jo),
    )


def _grouped_stats(
    J: np.ndarray,
    weights: np.ndarray,
    groups: np.ndarray,
) -> dict[str, Any]:
    out: dict[str, Any] = {}
    zone_count = len(POWER_LAW_ZONES)
    for tid, time in enumerate(TIMES):
        for zid, (name, _, _) in enumerate(POWER_LAW_ZONES):
            gid = tid * zone_count + zid
            mask = groups == gid
            row = _stats(J[mask], weights[mask])
            row.pop("worst_index", None)
            out[f"t={time:.2f}/{name}"] = row
    return out


def metrics_for_step(
    total: Any,
    leading: Any,
    probes: ProbeSet,
    step: float,
) -> dict[str, Any]:
    inner, ji = _staged(
        total,
        leading,
        probes.inner.points,
        probes.inner.times,
        step,
        probes.inner.weights,
    )
    power, jp = _staged(
        total,
        leading,
        probes.power_law.points,
        probes.power_law.times,
        step,
        probes.power_law.weights,
    )
    seam, _ = _staged(
        total,
        leading,
        probes.seam_points,
        probes.seam_times,
        step,
    )
    axis, _ = _staged(
        total,
        leading,
        probes.axis_points,
        probes.axis_times,
        step,
    )
    power["per_time_zone"] = {
        name: _grouped_stats(
            J,
            probes.power_law.weights,
            probes.power_law.group_index,
        )
        for name, J in zip(
            ("total", "leading", "oscillatory_increment"),
            jp,
        )
    }
    for region, points, times, mats in (
        ("inner", probes.inner.points, probes.inner.times, ji),
        ("power_law", probes.power_law.points, probes.power_law.times, jp),
    ):
        target = inner if region == "inner" else power
        for name, J in zip(("total", "leading", "oscillatory_increment"), mats):
            idx = int(_stats(J)["worst_index"])
            target[name]["worst_point"] = points[idx].tolist()
            target[name]["worst_time"] = float(times[idx])
    return {
        "step": float(step),
        "inner": inner,
        "power_law": power,
        "X3_seam": seam,
        "axis": axis,
    }


def _check_upstream_identity() -> None:
    if getattr(_upstream, "TASK", None) != UPSTREAM_TASK:
        raise RuntimeError("A2 #1010 task identity drifted")
    if getattr(_upstream, "SCHEMA", None) != UPSTREAM_SCHEMA:
        raise RuntimeError("A2 #1010 schema identity drifted")
    if getattr(_upstream, "AGENT1_PR", None) != AGENT1_PR:
        raise RuntimeError("A2 #1010 Agent-1 binding drifted")
    if getattr(_upstream, "AGENT1_HEAD", None) != AGENT1_HEAD:
        raise RuntimeError("A2 #1010 Agent-1 head binding drifted")
    if getattr(_upstream, "AGENT1_SOURCE_BLOB_SHA1", None) != AGENT1_SOURCE_BLOB:
        raise RuntimeError("A2 #1010 Agent-1 source binding drifted")
    cls = getattr(_upstream, "CurrentRF40PowerLawLeadingOscillatoryField", None)
    if cls is None or cls.__module__ != _upstream.__name__:
        raise RuntimeError("A2 #1010 public candidate class drifted")


def _reconstruct_leading_from_saved(total: Any) -> Any:
    cfg = total.configuration()
    block = cfg.get("agent1_leading")
    if (
        not isinstance(block, Mapping)
        or block.get("pr") != AGENT1_PR
        or block.get("head") != AGENT1_HEAD
        or block.get("source_blob_sha1") != AGENT1_SOURCE_BLOB
    ):
        raise RuntimeError("saved artifact does not bind exact A1 #1005")
    configuration = block.get("configuration")
    if not isinstance(configuration, Mapping):
        raise RuntimeError("saved artifact lacks A1 #1005 configuration")
    module = importlib.import_module(_upstream.AGENT1_MODULE)
    cls = getattr(module, _upstream.AGENT1_CLASS)
    source_path = inspect.getsourcefile(cls)
    if source_path is None or _git_blob_sha1(source_path) != AGENT1_SOURCE_BLOB:
        raise RuntimeError("A1 #1005 source identity drifted")
    leading = cls.from_configuration(configuration)
    if str(leading.semantic_sha256) != str(block.get("semantic_sha256")):
        raise RuntimeError("A1 #1005 semantic identity drifted during reload")
    return leading


class _DivergenceMutation:
    def __init__(self, field: Any, epsilon: float = MUTATION_EPSILON):
        self.field, self.epsilon = field, float(epsilon)

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        out = np.asarray(self.field.velocity(x, y, z, t), float).copy()
        out[..., 0] += self.epsilon * np.asarray(x, float)
        return out


class _ManufacturedSolenoidal:
    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        x, y, z, t = np.broadcast_arrays(
            np.asarray(x, float),
            np.asarray(y, float),
            np.asarray(z, float),
            np.asarray(t, float),
        )
        s = 1.0 + 0.2 * t
        return np.stack(
            (s * (x + 2.0 * y), s * (-y + 3.0 * z), s * (4.0 * x)),
            axis=-1,
        )


def _recursive_eta_fd_perturb(payload: Any) -> bool:
    if isinstance(payload, dict):
        for key, value in payload.items():
            if key == "eta_fd_step" and isinstance(value, (int, float)):
                payload[key] = float(value) * 1.001
                return True
            if _recursive_eta_fd_perturb(value):
                return True
    elif isinstance(payload, list):
        for value in payload:
            if _recursive_eta_fd_perturb(value):
                return True
    return False


def _stage_gate(metrics: Mapping[str, Any], region: str) -> None:
    block = metrics[region]
    for name in ("total", "leading", "oscillatory_increment"):
        if float(block[name]["sampled_max"]) > DIVERGENCE_GATE:
            raise AssertionError(f"{region}/{name} sampled divergence gate failed")
        if float(block[name]["weighted_rms"]) > DIVERGENCE_GATE:
            raise AssertionError(f"{region}/{name} weighted RMS divergence gate failed")


def _point_gate(metrics: Mapping[str, Any], region: str) -> None:
    for name in ("total", "leading", "oscillatory_increment"):
        if float(metrics[region][name]["sampled_max"]) > DIVERGENCE_GATE:
            raise AssertionError(f"{region}/{name} sampled divergence gate failed")


def _stability_gate(
    medium: Mapping[str, Any],
    fine: Mapping[str, Any],
) -> None:
    for region in ("inner", "power_law"):
        for name in ("total", "leading", "oscillatory_increment"):
            for metric in ("sampled_max", "weighted_rms"):
                m = float(medium[region][name][metric])
                f = float(fine[region][name][metric])
                if (
                    m > STABILITY_FLOOR
                    and f > max(STABILITY_FLOOR, STABILITY_FACTOR * m)
                ):
                    raise AssertionError(
                        f"resolution instability at {region}/{name}/{metric}"
                    )


def _aggregate_signature(metrics: Mapping[str, Any]) -> tuple[float, ...]:
    vals: list[float] = []
    for region in ("inner", "power_law"):
        for name in ("total", "leading", "oscillatory_increment"):
            vals.extend(
                (
                    float(metrics[region][name]["sampled_max"]),
                    float(metrics[region][name]["weighted_rms"]),
                )
            )
    return tuple(vals)


def materialize_receipt() -> dict[str, Any]:
    """Run the frozen independent audit and return a JSON-serializable receipt."""
    _check_upstream_identity()
    field = _upstream.default_field()

    with tempfile.TemporaryDirectory() as td:
        path = Path(td) / "candidate.json"
        saved = field.save_candidate(path)
        loaded = _upstream.CurrentRF40PowerLawLeadingOscillatoryField.load_candidate(path)
        raw_saved = json.loads(path.read_text(encoding="utf-8"))

    if (
        loaded.semantic_sha256 != field.semantic_sha256
        or loaded.semantic_sha256 != saved["semantic_sha256"]
    ):
        raise RuntimeError("A2 #1010 semantic identity did not survive save/load")

    leading = _reconstruct_leading_from_saved(loaded)
    probes = make_probes(loaded)
    ladder = [metrics_for_step(loaded, leading, probes, h) for h in STEPS]

    replay_points = np.concatenate(
        (probes.inner.points[:8], probes.power_law.points[:8])
    )
    replay_times = np.concatenate(
        (probes.inner.times[:8], probes.power_law.times[:8])
    )
    replay = float(
        np.max(
            np.abs(
                _velocity(field, replay_points, replay_times)
                - _velocity(loaded, replay_points, replay_times)
            )
        )
    )

    fine_total = fd2_jacobian(
        loaded,
        probes.inner.points,
        probes.inner.times,
        STEPS[-1],
    )
    fine_mut = fd2_jacobian(
        _DivergenceMutation(loaded),
        probes.inner.points,
        probes.inner.times,
        STEPS[-1],
    )
    mutation_detected = float(
        np.max(
            np.abs(
                np.trace(fine_mut, axis1=1, axis2=2)
                - np.trace(fine_total, axis1=1, axis2=2)
            )
        )
    )

    reverse = ProbeSet(
        ProbeStratum(
            probes.inner.points[::-1],
            probes.inner.times[::-1],
            probes.inner.weights[::-1],
            probes.inner.group_index[::-1],
        ),
        ProbeStratum(
            probes.power_law.points[::-1],
            probes.power_law.times[::-1],
            probes.power_law.weights[::-1],
            probes.power_law.group_index[::-1],
        ),
        probes.seam_points[::-1],
        probes.seam_times[::-1],
        probes.axis_points[::-1],
        probes.axis_times[::-1],
    )
    reordered = metrics_for_step(loaded, leading, reverse, STEPS[-1])
    ordering_error = float(
        max(
            abs(a - b)
            for a, b in zip(
                _aggregate_signature(ladder[-1]),
                _aggregate_signature(reordered),
            )
        )
    )

    manufactured = _ManufacturedSolenoidal()
    manufactured_J = fd2_jacobian(
        manufactured,
        probes.inner.points,
        probes.inner.times,
        STEPS[-1],
    )
    manufactured_max = float(
        np.max(np.abs(np.trace(manufactured_J, axis1=1, axis2=2)))
    )

    cfg_mut = copy.deepcopy(raw_saved["configuration"])
    a1cfg = cfg_mut["agent1_leading"]["configuration"]
    if not _recursive_eta_fd_perturb(a1cfg):
        raise RuntimeError("A1 configuration lacks preregistered eta_fd_step target")
    mutated_leading_module = importlib.import_module(_upstream.AGENT1_MODULE)
    mutated_leading_cls = getattr(mutated_leading_module, _upstream.AGENT1_CLASS)
    mutated_leading = mutated_leading_cls.from_configuration(a1cfg)
    a1_semantic_mutation_detected = str(mutated_leading.semantic_sha256) != str(
        raw_saved["configuration"]["agent1_leading"]["semantic_sha256"]
    )

    osc_mut = copy.deepcopy(raw_saved)
    osc = osc_mut["configuration"]["oscillatory_runtime"]
    osc["payload"]["parameters"]["mode_imaginary_ratio"] += 0.125
    osc["payload_sha256"] = _upstream._identity._sha256(osc["payload"])
    osc_mut["semantic_sha256"] = _upstream._identity._sha256(
        osc_mut["configuration"]
    )
    oscillatory_mutation_rejected = False
    with tempfile.TemporaryDirectory() as td:
        p = Path(td) / "mutated.json"
        p.write_text(
            json.dumps(osc_mut, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        try:
            _upstream.CurrentRF40PowerLawLeadingOscillatoryField.load_candidate(p)
        except ValueError:
            oscillatory_mutation_rejected = True

    beyond = _source_xyz(
        loaded,
        [1.001 * loaded.X_4],
        [0.2],
        [0.37],
        [0.47],
    )[0]
    beyond_fail_closed = False
    try:
        loaded.velocity(beyond[0], beyond[1], beyond[2], 0.47)
    except (ValueError, RuntimeError):
        beyond_fail_closed = True

    return {
        "schema": SCHEMA,
        "task": TASK,
        "seed": SEED,
        "steps": list(STEPS),
        "times": list(TIMES),
        "upstream": {
            "pr": UPSTREAM_PR,
            "head": UPSTREAM_HEAD,
            "task": UPSTREAM_TASK,
            "schema": UPSTREAM_SCHEMA,
            "semantic_sha256": loaded.semantic_sha256,
        },
        "agent1": {
            "pr": AGENT1_PR,
            "head": AGENT1_HEAD,
            "source_blob_sha1": AGENT1_SOURCE_BLOB,
            "semantic_sha256": str(leading.semantic_sha256),
        },
        "probe_counts": {
            "inner": int(len(probes.inner.points)),
            "power_law": int(len(probes.power_law.points)),
            "X3_seam": int(len(probes.seam_points)),
            "axis": int(len(probes.axis_points)),
        },
        "resolution_ladder": ladder,
        "firewalls": {
            "save_load_velocity_replay_max_abs": replay,
            "divergence_mutation_detected": mutation_detected,
            "ordering_aggregate_max_abs": ordering_error,
            "manufactured_solenoidal_divergence_max_abs": manufactured_max,
            "agent1_semantic_mutation_detected": bool(
                a1_semantic_mutation_detected
            ),
            "oscillatory_runtime_mutation_rejected": bool(
                oscillatory_mutation_rejected
            ),
            "beyond_X4_fail_closed": bool(beyond_fail_closed),
        },
        "frozen_gates": {
            "divergence_sampled_max": DIVERGENCE_GATE,
            "divergence_weighted_rms": DIVERGENCE_GATE,
            "stability_factor": STABILITY_FACTOR,
            "stability_floor": STABILITY_FLOOR,
            "total_speed_rms_min": TOTAL_SPEED_FLOOR,
            "oscillatory_speed_rms_min": OSCILLATORY_SPEED_FLOOR,
            "mutation_detection_min": MUTATION_DETECTION_FLOOR,
            "save_load_replay_max_abs": REPLAY_ATOL,
            "ordering_max_abs": ORDER_ATOL,
            "formal_momentum_max_L2_gate": 1.0e-3,
            "formal_divergence_max_L2_gate": 1.0e-5,
        },
        "scientific_scope_note": (
            "scoped divergence only; physical-volume L2 is an estimate on the "
            "registered inner/power-law strata and is not canonical [24,48,96] "
            "whole-domain admission or momentum/full-NS evidence"
        ),
        "truth_boundary": dict(TRUTH_BOUNDARY),
        "receipt_sha256": "",
    }


def enforce_receipt(receipt: Mapping[str, Any]) -> None:
    if receipt.get("schema") != SCHEMA or receipt.get("task") != TASK:
        raise AssertionError("receipt schema/task drifted")
    if (
        receipt.get("seed") != SEED
        or tuple(receipt.get("steps", ())) != STEPS
        or tuple(receipt.get("times", ())) != TIMES
    ):
        raise AssertionError("frozen protocol drifted")
    if receipt.get("truth_boundary") != TRUTH_BOUNDARY:
        raise AssertionError("truth boundary drifted")

    upstream = receipt.get("upstream", {})
    if (
        upstream.get("pr") != UPSTREAM_PR
        or upstream.get("head") != UPSTREAM_HEAD
        or upstream.get("task") != UPSTREAM_TASK
        or upstream.get("schema") != UPSTREAM_SCHEMA
    ):
        raise AssertionError("A2 #1010 identity drifted")

    agent1 = receipt.get("agent1", {})
    if (
        agent1.get("pr") != AGENT1_PR
        or agent1.get("head") != AGENT1_HEAD
        or agent1.get("source_blob_sha1") != AGENT1_SOURCE_BLOB
    ):
        raise AssertionError("A1 #1005 identity drifted")

    counts = receipt.get("probe_counts", {})
    if counts != {
        "inner": 32,
        "power_law": 72,
        "X3_seam": 16,
        "axis": 5,
    }:
        raise AssertionError("held-out probe counts drifted")

    ladder = receipt.get("resolution_ladder")
    if not isinstance(ladder, list) or len(ladder) != 3:
        raise AssertionError("three-resolution ladder missing")
    for got, expected in zip(ladder, STEPS):
        if not math.isclose(
            float(got["step"]),
            expected,
            rel_tol=0.0,
            abs_tol=0.0,
        ):
            raise AssertionError("resolution step drifted")

    fine, medium = ladder[-1], ladder[-2]
    _stage_gate(fine, "inner")
    _stage_gate(fine, "power_law")
    _point_gate(fine, "X3_seam")
    _point_gate(fine, "axis")
    if float(fine["inner"]["total_speed_rms"]) < TOTAL_SPEED_FLOOR:
        raise AssertionError("candidate became trivial")
    if (
        float(fine["inner"]["oscillatory_speed_rms"])
        < OSCILLATORY_SPEED_FLOOR
    ):
        raise AssertionError("oscillatory signal became trivial")
    _stability_gate(medium, fine)

    fw = receipt.get("firewalls", {})
    if float(fw.get("save_load_velocity_replay_max_abs", math.inf)) > REPLAY_ATOL:
        raise AssertionError("save/load public velocity replay failed")
    if (
        float(fw.get("divergence_mutation_detected", -math.inf))
        < MUTATION_DETECTION_FLOOR
    ):
        raise AssertionError("divergence mutation was not detected")
    if float(fw.get("ordering_aggregate_max_abs", math.inf)) > ORDER_ATOL:
        raise AssertionError("held-out ordering changed aggregate metrics")
    if (
        float(fw.get("manufactured_solenoidal_divergence_max_abs", math.inf))
        > 1.0e-11
    ):
        raise AssertionError("manufactured solenoidal calibration failed")
    for key in (
        "agent1_semantic_mutation_detected",
        "oscillatory_runtime_mutation_rejected",
        "beyond_X4_fail_closed",
    ):
        if fw.get(key) is not True:
            raise AssertionError(f"firewall failed: {key}")


def write_receipt(path: str | Path) -> dict[str, Any]:
    receipt = materialize_receipt()
    core = dict(receipt)
    core.pop("receipt_sha256", None)
    receipt["receipt_sha256"] = _sha256(core)
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return receipt


__all__ = [
    "SCHEMA",
    "TASK",
    "SEED",
    "STEPS",
    "TIMES",
    "TRUTH_BOUNDARY",
    "make_probes",
    "fd2_jacobian",
    "metrics_for_step",
    "materialize_receipt",
    "enforce_receipt",
    "write_receipt",
]
