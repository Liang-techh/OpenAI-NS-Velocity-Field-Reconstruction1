"""Independent divergence revalidation for a support-transformed Eq. (4.5) child.

The support transform changes ``[u,v,w]`` in the exterior collar, so divergence
evidence from the untapered parent cannot be inherited. This module treats a
serialized :class:`Eq45SupportedVelocityCandidate` as a black-box public field:
every sample is obtained through ``at_points(points, time)`` and Cartesian
centered differences are reconstructed outside the candidate implementation.

The audit reports the preregistered three-level derivative ladder and preserves
the unchanged CR001 divergence thresholds as reference metadata. The production
probe statistic here is a held-out sampled RMS, not the preregistered
volume-weighted spatial L2 norm, and its sampled maximum is not a global-domain
maximum. Therefore this module deliberately does *not* assess the formal CR001
divergence acceptance gate. Passing sampled divergence is only an
incompressibility diagnostic; it does not establish momentum balance,
visualization correspondence, paper exactness, or identification of an OpenAI
hidden field.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Iterable

import numpy as np

from .constrained_eq45_supported_candidate import Eq45SupportedVelocityCandidate


CLAIM_SCOPE = "supported_child_independent_public_velocity_divergence_only"
MIN_ACTIVITY_RMS = 1e-12


def _validated_points_times_labels(points, times, labels):
    pts = np.asarray(points, dtype=float)
    if pts.ndim != 2 or pts.shape[1] != 3 or pts.shape[0] == 0:
        raise ValueError("points must have shape (n, 3) with n > 0")
    if not np.all(np.isfinite(pts)):
        raise ValueError("points must be finite")

    t = np.asarray(times, dtype=float)
    if not np.all(np.isfinite(t)):
        raise ValueError("times must be finite")
    try:
        t = np.broadcast_to(t, (pts.shape[0],)).copy()
    except ValueError as exc:
        raise ValueError("times must be scalar or broadcast to one value per point") from exc

    region = np.asarray(labels, dtype=object)
    if region.shape != (pts.shape[0],):
        raise ValueError("labels must provide one region label per point")
    if any(not isinstance(value, str) or not value for value in region.tolist()):
        raise ValueError("region labels must be nonempty strings")
    return pts, t, region


def _validated_steps(steps: Iterable[float]) -> tuple[float, ...]:
    values = tuple(float(step) for step in steps)
    if len(values) < 3:
        raise ValueError("at least three spatial derivative steps are required")
    if not all(np.isfinite(step) and step > 0.0 for step in values):
        raise ValueError("spatial derivative steps must be finite and positive")
    if not all(values[index + 1] < values[index] for index in range(len(values) - 1)):
        raise ValueError("spatial derivative steps must be strictly decreasing")
    return values


def _public_velocity(field, points: np.ndarray, times: np.ndarray) -> np.ndarray:
    evaluator = getattr(field, "at_points", None)
    if not callable(evaluator):
        raise TypeError("field must expose callable at_points(points, time)")
    values = np.asarray(evaluator(points, times), dtype=float)
    if values.shape != points.shape:
        raise ValueError("public velocity output must have shape (n, 3)")
    if not np.all(np.isfinite(values)):
        raise ValueError("public velocity output must be finite")
    return values


def _gradient(field, points: np.ndarray, times: np.ndarray, step: float) -> np.ndarray:
    gradient = np.empty((points.shape[0], 3, 3), dtype=float)
    for axis in range(3):
        offset = np.zeros(3, dtype=float)
        offset[axis] = step
        plus = _public_velocity(field, points + offset, times)
        minus = _public_velocity(field, points - offset, times)
        gradient[:, :, axis] = (plus - minus) / (2.0 * step)
    return gradient


def _metrics(divergence: np.ndarray, gradient: np.ndarray) -> dict[str, float]:
    rms = float(np.sqrt(np.mean(divergence * divergence)))
    gradient_rms = float(np.sqrt(np.mean(np.sum(gradient * gradient, axis=(1, 2)))))
    if not np.isfinite(gradient_rms) or gradient_rms <= 0.0:
        raise ValueError("sampled velocity-gradient norm is not usable")
    return {
        "max_abs": float(np.max(np.abs(divergence))),
        "rms": rms,
        "gradient_frobenius_rms": gradient_rms,
        "normalized_rms": rms / gradient_rms,
    }


def audit_public_velocity_divergence(
    field,
    points,
    times,
    labels,
    *,
    steps=(0.02, 0.01, 0.005),
) -> dict[str, Any]:
    """Audit one public field at fixed probes while varying only derivative step."""

    pts, t, region = _validated_points_times_labels(points, times, labels)
    derivative_steps = _validated_steps(steps)

    base = _public_velocity(field, pts, t)
    velocity_rms = float(np.sqrt(np.mean(np.sum(base * base, axis=1))))
    if not np.isfinite(velocity_rms) or velocity_rms <= MIN_ACTIVITY_RMS:
        raise ValueError("reference public velocity is numerically inactive")

    unique_regions = tuple(dict.fromkeys(region.tolist()))
    levels: list[dict[str, Any]] = []
    for step in derivative_steps:
        gradient = _gradient(field, pts, t, step)
        divergence = gradient[:, 0, 0] + gradient[:, 1, 1] + gradient[:, 2, 2]
        by_region = {}
        for name in unique_regions:
            mask = region == name
            by_region[name] = _metrics(divergence[mask], gradient[mask])
        levels.append(
            {
                "step": step,
                **_metrics(divergence, gradient),
                "by_region": by_region,
            }
        )

    observed_orders: list[float | None] = []
    for coarse, fine in zip(levels[:-1], levels[1:]):
        if coarse["rms"] <= 0.0 or fine["rms"] <= 0.0:
            observed_orders.append(None)
        else:
            observed_orders.append(
                float(
                    np.log(coarse["rms"] / fine["rms"])
                    / np.log(coarse["step"] / fine["step"])
                )
            )

    return {
        "claim_scope": CLAIM_SCOPE,
        "operator": "centered_second_order_cartesian_space",
        "velocity_access": "at_points_only",
        "point_count": int(pts.shape[0]),
        "regions": list(unique_regions),
        "velocity_rms": velocity_rms,
        "levels": levels,
        "observed_rms_orders": observed_orders,
        "sampled_rms_only": True,
        "global_domain_max_assessed": False,
        "volume_weighted_L2_assessed": False,
        "pde_validated": False,
        "visualization_ready": False,
        "visual_correspondence_verified": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
    }


def _load_validation_contract(constraints_path: str | Path) -> dict[str, Any]:
    payload = json.loads(Path(constraints_path).read_text(encoding="utf-8"))
    validation = payload.get("validation")
    if not isinstance(validation, dict):
        raise ValueError("constraints payload is missing validation contract")
    steps = _validated_steps(validation.get("derivative_steps", ()))
    thresholds = validation.get("thresholds")
    if not isinstance(thresholds, dict):
        raise ValueError("constraints payload is missing validation thresholds")
    divergence_max = float(thresholds["divergence_max"])
    divergence_l2 = float(thresholds["divergence_L2"])
    if divergence_max <= 0.0 or divergence_l2 <= 0.0:
        raise ValueError("divergence thresholds must be positive")
    return {
        "steps": steps,
        "divergence_max": divergence_max,
        "divergence_L2": divergence_l2,
    }


def audit_supported_candidate_artifact_divergence(
    candidate_path: str | Path,
    constraints_path: str | Path,
    points,
    times,
    labels,
) -> dict[str, Any]:
    """Reload a supported child and compare samples to the unchanged contract."""

    field = Eq45SupportedVelocityCandidate.load_json(candidate_path)
    contract = _load_validation_contract(constraints_path)
    report = audit_public_velocity_divergence(
        field,
        points,
        times,
        labels,
        steps=contract["steps"],
    )
    finest = report["levels"][-1]
    report.update(
        {
            "candidate_sha256": field.sha256,
            "parent_sha256": field.parent_sha256,
            "artifact_reloaded": True,
            "registered_thresholds": {
                "divergence_max": contract["divergence_max"],
                "divergence_L2": contract["divergence_L2"],
            },
            "sampled_threshold_indicators": {
                "sampled_max_below_registered_max": bool(
                    finest["max_abs"] <= contract["divergence_max"]
                ),
                "sampled_rms_below_registered_L2_number": bool(
                    finest["rms"] <= contract["divergence_L2"]
                ),
            },
            "cr001_divergence_gate_assessed": False,
            "full_momentum_residual_assessed": False,
            "physical_support_validated": False,
        }
    )
    return report


class PublicDivergenceMutation:
    """Public-output mutation adding ``epsilon*x`` to the x velocity component."""

    def __init__(self, field, epsilon: float = 0.03):
        self.field = field
        self.epsilon = float(epsilon)
        if not np.isfinite(self.epsilon) or self.epsilon == 0.0:
            raise ValueError("epsilon must be finite and nonzero")

    def at_points(self, points, time):
        pts = np.asarray(points, dtype=float)
        values = _public_velocity(
            self.field,
            pts,
            np.broadcast_to(np.asarray(time, dtype=float), (pts.shape[0],)),
        )
        out = values.copy()
        out[:, 0] += self.epsilon * pts[:, 0]
        return out
