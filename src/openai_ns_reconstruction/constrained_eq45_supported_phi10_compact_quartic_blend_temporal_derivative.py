"""Independent public-output time-derivative audit for the compact--quartic blend.

Agent 1 PR #186 materializes a bounded one-scalar interpolation between the
already-audited derivative-balanced quartic and slope-capped compact-C2
``Phi(1,0)`` temporal schedules.  This module does not choose that scalar.  It
asks only whether representative serialized/reloaded blend members have a
resolved public ``u_t`` and whether the 50/50 member preserves the endpoint
interpolation at the derivative level.

The derivative operator is the independent public-velocity finite-difference
methodology already introduced for the cubic candidate.  No training loss,
pressure, force, vorticity residual, hidden tensor, or analytic internal
velocity derivative is used here.  This remains a CR009 convergence diagnostic
for a visualization/research candidate, not the formal momentum PDE gate.
"""
from __future__ import annotations

import json
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Iterable

import numpy as np

from .constrained_eq45_candidate import Eq45VelocityCandidate
from .constrained_eq45_supported_candidate import Eq45SupportedVelocityCandidate
from .constrained_eq45_supported_phi10_compact_quartic_blend_temporal_mode import (
    Eq45SupportedPhi10CompactQuarticBlendTemporalCandidate,
)
from .constrained_eq45_supported_phi10_cubic_temporal_derivative import (
    REFERENCE_STEP,
    TIME_LEVELS,
    _fourth_order_time_derivative,
    _fresh_probe_points,
    _manufactured_sign_mutation_calibration,
    _metrics,
    _second_order_time_derivative,
)


SCHEMA = "eq45_supported_phi10_compact_quartic_blend_temporal_derivative_audit_v1"
PROBE_SEED = 914617
BLEND_WEIGHTS = (0.0, 0.5, 1.0)
AUDIT_TIMES = (0.25, 0.3125, 0.375, 0.4375, 0.50, 0.75)

_TRUTH_BOUNDARY = {
    "serialized_candidates_reloaded": True,
    "public_velocity_only": True,
    "blend_weight_selected": False,
    "training_loss_reused": False,
    "pressure_fitted": False,
    "forcing_fitted": False,
    "vorticity_or_momentum_residual_evaluated": False,
    "formal_pde_gate_assessed": False,
    "visualization_candidate_only": True,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _validated_weights(weights: Iterable[float]) -> tuple[float, ...]:
    values = tuple(float(value) for value in weights)
    if len(values) < 3:
        raise ValueError("at least three blend weights are required")
    if not all(np.isfinite(value) and 0.0 <= value <= 1.0 for value in values):
        raise ValueError("blend weights must be finite and lie in [0,1]")
    if len(set(values)) != len(values):
        raise ValueError("blend weights must be unique")
    if 0.0 not in values or 0.5 not in values or 1.0 not in values:
        raise ValueError("audit requires frozen endpoint and midpoint weights 0, 0.5, 1")
    return values


def _serialize_reload(weight: float):
    base = Eq45SupportedVelocityCandidate(parent=Eq45VelocityCandidate.seed())
    candidate = Eq45SupportedPhi10CompactQuarticBlendTemporalCandidate(
        base=base,
        blend_weight=float(weight),
    )
    with TemporaryDirectory() as tmp:
        path = Path(tmp) / f"blend_{weight:.6f}.json"
        candidate.save_json(path)
        loaded = Eq45SupportedPhi10CompactQuarticBlendTemporalCandidate.load_json(path)
    if loaded.sha256 != candidate.sha256:
        raise RuntimeError("blend candidate serialization changed identity")
    return loaded


def _rms(values: np.ndarray) -> float:
    array = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(array * array)))


def _convergence(candidate, points: np.ndarray) -> dict[str, object]:
    reference = np.stack(
        [
            _fourth_order_time_derivative(candidate, points, float(time), REFERENCE_STEP)
            for time in AUDIT_TIMES
        ],
        axis=0,
    )
    levels: list[dict[str, float]] = []
    errors: list[float] = []
    for step in TIME_LEVELS:
        estimate = np.stack(
            [
                _second_order_time_derivative(candidate, points, float(time), float(step))
                for time in AUDIT_TIMES
            ],
            axis=0,
        )
        metrics = _metrics(estimate - reference, reference)
        levels.append({"step": float(step), **metrics.to_dict()})
        errors.append(float(metrics.sampled_rms))

    observed_orders = [
        float(np.log(left / right) / np.log(2.0))
        for left, right in zip(errors[:-1], errors[1:])
    ]
    return {
        "levels": levels,
        "observed_orders": observed_orders,
        "reference_step": float(REFERENCE_STEP),
        "reference_operator": "fourth_order_public_velocity_time_difference",
    }


def _affinity_rows(candidates: dict[float, object], points: np.ndarray) -> list[dict[str, float]]:
    quartic = candidates[0.0]
    midpoint = candidates[0.5]
    compact = candidates[1.0]
    rows: list[dict[str, float]] = []
    for time in AUDIT_TIMES:
        q_value = np.asarray(quartic.at_points(points, float(time)), dtype=float)
        m_value = np.asarray(midpoint.at_points(points, float(time)), dtype=float)
        c_value = np.asarray(compact.at_points(points, float(time)), dtype=float)
        expected_value = 0.5 * (q_value + c_value)
        value_error = m_value - expected_value

        q_dt = _fourth_order_time_derivative(quartic, points, float(time), REFERENCE_STEP)
        m_dt = _fourth_order_time_derivative(midpoint, points, float(time), REFERENCE_STEP)
        c_dt = _fourth_order_time_derivative(compact, points, float(time), REFERENCE_STEP)
        expected_dt = 0.5 * (q_dt + c_dt)
        dt_error = m_dt - expected_dt
        expected_dt_rms = max(_rms(expected_dt), np.finfo(float).tiny)
        rows.append(
            {
                "time": float(time),
                "instantaneous_velocity_affine_rms_error": _rms(value_error),
                "instantaneous_velocity_affine_max_vector_error": float(
                    np.max(np.linalg.norm(value_error, axis=1))
                ),
                "time_derivative_affine_rms_error": _rms(dt_error),
                "time_derivative_affine_max_vector_error": float(
                    np.max(np.linalg.norm(dt_error, axis=1))
                ),
                "time_derivative_affine_normalized_rms_error": float(
                    _rms(dt_error) / expected_dt_rms
                ),
                "endpoint_time_derivative_span_rms": _rms(c_dt - q_dt),
            }
        )
    return rows


def audit_blend_temporal_derivative(
    *,
    weights: Iterable[float] = BLEND_WEIGHTS,
    probe_seed: int = PROBE_SEED,
) -> dict[str, object]:
    """Audit public ``u_t`` convergence for frozen blend endpoints and midpoint."""
    checked_weights = _validated_weights(weights)
    if not isinstance(probe_seed, (int, np.integer)) or int(probe_seed) < 0:
        raise ValueError("probe_seed must be a nonnegative integer")

    points, regions = _fresh_probe_points(int(probe_seed))
    candidates = {weight: _serialize_reload(weight) for weight in checked_weights}
    convergence = {
        f"{weight:.6f}": _convergence(candidates[weight], points)
        for weight in checked_weights
    }
    affinity = _affinity_rows(candidates, points)

    return {
        "schema": SCHEMA,
        "probe_seed": int(probe_seed),
        "probe_count": int(len(points)),
        "probe_region_counts": {
            name: int(sum(region == name for region in regions))
            for name in sorted(set(regions))
        },
        "blend_weights": [float(weight) for weight in checked_weights],
        "audit_times": [float(time) for time in AUDIT_TIMES],
        "time_levels": [float(step) for step in TIME_LEVELS],
        "reference_step": float(REFERENCE_STEP),
        "candidate_sha256": {
            f"{weight:.6f}": candidates[weight].sha256 for weight in checked_weights
        },
        "convergence": convergence,
        "midpoint_endpoint_affinity": affinity,
        "mutation_calibration": _manufactured_sign_mutation_calibration(),
        "truth_boundary": dict(_TRUTH_BOUNDARY),
    }


def write_report(path: str | Path, **kwargs) -> dict[str, object]:
    report = audit_blend_temporal_derivative(**kwargs)
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report


if __name__ == "__main__":
    print(json.dumps(audit_blend_temporal_derivative(), indent=2, sort_keys=True))
