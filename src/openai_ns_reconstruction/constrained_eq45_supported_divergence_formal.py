"""Registered held-out divergence metrics for the support-connected Eq. (4.5) child.

This module is a narrow extension of the independent public-output divergence
operator in :mod:`constrained_eq45_supported_divergence`.  It closes one metric
scope left deliberately open there: the preregistered CR001 validation contract
uses 4096 held-out Cartesian points, six fixed validation times, a three-level
spatial derivative ladder, a sampled maximum, and a volume-weighted spatial L2
norm.

The field is loaded from a serialized ``Eq45SupportedVelocityCandidate`` and all
velocity values enter only through ``at_points(points, time)``.  The held-out
spatial points are generated deterministically from the preregistered validation
seed and reused at each fixed time so time-to-time changes are not confounded by
a different point cloud.  The volume-weighted L2 is the Monte-Carlo estimate
``sqrt(V * mean(div(u)^2))`` on the declared evaluation box.

This is still a finite numerical acceptance check, not a continuum proof.  In
particular, the 4096-point maximum is not claimed to be the supremum over R^3,
and passing divergence alone would not establish the Navier--Stokes momentum
equation, visualization correspondence, paper exactness, or identification of
an OpenAI hidden field.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .constrained_eq45_supported_candidate import Eq45SupportedVelocityCandidate
from .constrained_eq45_supported_divergence import _gradient, _validated_steps


CLAIM_SCOPE = "supported_child_registered_held_out_divergence_metrics"


def _load_contract(path: str | Path) -> dict[str, Any]:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(payload, Mapping):
        raise ValueError("constraints payload must be an object")

    domain = payload.get("domain")
    validation = payload.get("validation")
    if not isinstance(domain, Mapping) or not isinstance(validation, Mapping):
        raise ValueError("constraints payload is missing domain/validation contract")

    raw_box = np.asarray(domain.get("evaluation_box"), dtype=float)
    if raw_box.shape != (3, 2) or not np.all(np.isfinite(raw_box)):
        raise ValueError("evaluation_box must contain three finite [lo, hi] axes")
    lower = raw_box[:, 0]
    upper = raw_box[:, 1]
    if np.any(upper <= lower):
        raise ValueError("evaluation_box upper bounds must exceed lower bounds")
    volume = float(np.prod(upper - lower))

    seed = int(validation.get("seed"))
    point_count = int(validation.get("held_out_points"))
    if point_count <= 0:
        raise ValueError("held_out_points must be positive")

    times = np.asarray(validation.get("times"), dtype=float)
    if times.ndim != 1 or times.size == 0 or not np.all(np.isfinite(times)):
        raise ValueError("validation times must be a nonempty finite 1D sequence")
    if np.any(np.diff(times) <= 0.0):
        raise ValueError("validation times must be strictly increasing")

    time_interval = np.asarray(domain.get("time_interval"), dtype=float)
    if time_interval.shape != (2,) or not np.all(np.isfinite(time_interval)):
        raise ValueError("time_interval must contain two finite values")
    if np.any(times < time_interval[0]) or np.any(times > time_interval[1]):
        raise ValueError("validation times must lie inside the declared time interval")

    steps = _validated_steps(validation.get("derivative_steps", ()))
    thresholds = validation.get("thresholds")
    if not isinstance(thresholds, Mapping):
        raise ValueError("validation thresholds must be an object")
    divergence_max = float(thresholds.get("divergence_max"))
    divergence_l2 = float(thresholds.get("divergence_L2"))
    if not np.isfinite(divergence_max) or not np.isfinite(divergence_l2):
        raise ValueError("divergence thresholds must be finite")
    if divergence_max <= 0.0 or divergence_l2 <= 0.0:
        raise ValueError("divergence thresholds must be positive")

    norms = validation.get("norms")
    if not isinstance(norms, list) or "volume-weighted L2 spatial norm at each time" not in norms:
        raise ValueError("registered volume-weighted spatial L2 norm is missing")

    return {
        "seed": seed,
        "point_count": point_count,
        "times": tuple(float(value) for value in times),
        "steps": steps,
        "lower": lower,
        "upper": upper,
        "volume": volume,
        "divergence_max": divergence_max,
        "divergence_L2": divergence_l2,
    }


def registered_held_out_points(contract: Mapping[str, Any]) -> np.ndarray:
    """Generate the deterministic preregistered uniform Cartesian point cloud."""

    rng = np.random.default_rng(int(contract["seed"]))
    points = rng.uniform(
        np.asarray(contract["lower"], dtype=float),
        np.asarray(contract["upper"], dtype=float),
        size=(int(contract["point_count"]), 3),
    )
    if points.shape != (int(contract["point_count"]), 3):
        raise RuntimeError("held-out point generation returned the wrong shape")
    if not np.all(np.isfinite(points)):
        raise RuntimeError("held-out point generation returned nonfinite values")
    return points


def _spatial_metrics(field, points: np.ndarray, time: float, step: float, volume: float) -> dict[str, float]:
    times = np.full(points.shape[0], float(time), dtype=float)
    gradient = _gradient(field, points, times, float(step))
    divergence = gradient[:, 0, 0] + gradient[:, 1, 1] + gradient[:, 2, 2]
    squared_mean = float(np.mean(divergence * divergence))
    rms = float(np.sqrt(squared_mean))
    return {
        "sampled_max_abs": float(np.max(np.abs(divergence))),
        "sampled_rms": rms,
        "volume_weighted_L2": float(np.sqrt(float(volume) * squared_mean)),
    }


def audit_supported_candidate_registered_divergence(
    candidate_path: str | Path,
    constraints_path: str | Path,
) -> dict[str, Any]:
    """Assess the registered finite held-out divergence contract for one artifact."""

    field = Eq45SupportedVelocityCandidate.load_json(candidate_path)
    contract = _load_contract(constraints_path)
    points = registered_held_out_points(contract)

    levels: list[dict[str, Any]] = []
    for step in contract["steps"]:
        per_time = [
            {
                "time": float(time),
                **_spatial_metrics(field, points, time, step, contract["volume"]),
            }
            for time in contract["times"]
        ]
        levels.append(
            {
                "step": float(step),
                "per_time": per_time,
                "worst_sampled_max_abs": float(
                    max(item["sampled_max_abs"] for item in per_time)
                ),
                "worst_volume_weighted_L2": float(
                    max(item["volume_weighted_L2"] for item in per_time)
                ),
            }
        )

    orders_by_time: list[dict[str, Any]] = []
    for time_index, time in enumerate(contract["times"]):
        l2_values = [level["per_time"][time_index]["volume_weighted_L2"] for level in levels]
        orders: list[float | None] = []
        for coarse_index in range(len(levels) - 1):
            coarse = l2_values[coarse_index]
            fine = l2_values[coarse_index + 1]
            if coarse <= 0.0 or fine <= 0.0:
                orders.append(None)
            else:
                orders.append(
                    float(
                        np.log(coarse / fine)
                        / np.log(levels[coarse_index]["step"] / levels[coarse_index + 1]["step"])
                    )
                )
        orders_by_time.append({"time": float(time), "volume_L2_orders": orders})

    finest = levels[-1]
    max_pass = finest["worst_sampled_max_abs"] <= contract["divergence_max"]
    l2_pass = finest["worst_volume_weighted_L2"] <= contract["divergence_L2"]

    return {
        "claim_scope": CLAIM_SCOPE,
        "candidate_sha256": field.sha256,
        "parent_sha256": field.parent_sha256,
        "artifact_reloaded": True,
        "velocity_access": "at_points_only",
        "sampling": {
            "distribution": "uniform_cartesian_box",
            "seed": int(contract["seed"]),
            "point_count": int(contract["point_count"]),
            "same_spatial_cloud_reused_at_each_time": True,
            "evaluation_box": [
                [float(lo), float(hi)]
                for lo, hi in zip(contract["lower"], contract["upper"])
            ],
            "box_volume": float(contract["volume"]),
            "times": list(contract["times"]),
        },
        "operator": "independent_centered_second_order_cartesian_space",
        "norm_definition": "sqrt(box_volume * mean(divergence^2)) at each fixed time",
        "levels": levels,
        "observed_volume_L2_orders_by_time": orders_by_time,
        "registered_thresholds": {
            "divergence_max": float(contract["divergence_max"]),
            "divergence_L2": float(contract["divergence_L2"]),
        },
        "finest_level": {
            "step": finest["step"],
            "worst_sampled_max_abs": finest["worst_sampled_max_abs"],
            "worst_volume_weighted_L2": finest["worst_volume_weighted_L2"],
            "registered_sample_max_pass": bool(max_pass),
            "registered_volume_L2_pass": bool(l2_pass),
        },
        "registered_finite_divergence_gate_assessed": True,
        "registered_finite_divergence_gate_passed": bool(max_pass and l2_pass),
        "continuum_supremum_proved": False,
        "full_momentum_residual_assessed": False,
        "pde_validated": False,
        "physical_support_validated": False,
        "visualization_ready": False,
        "visual_correspondence_verified": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
    }
