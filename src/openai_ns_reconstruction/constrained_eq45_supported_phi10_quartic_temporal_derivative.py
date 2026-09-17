"""Independent public-output time-derivative audit for the quartic Phi(1,0) trial.

The derivative-balanced quartic supported Eq45 candidate was selected target-free
to reduce coefficient-slope collateral at static-return keyframes.  This module
checks whether that representation-level improvement survives at the actual
public Cartesian ``[u,v,w]`` interface and whether its time derivative is
numerically resolved.

The candidate is serialized/reloaded before validation.  Every derivative is
reconstructed only from public ``at_points(...)->[u,v,w]`` evaluations.  The
three-level second-order time ladder is compared against a separate fourth-order
finer public-velocity reference.  This is a CR009 convergence/generalization
diagnostic, not the preregistered full momentum PDE gate.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import tempfile

import numpy as np

from .constrained_eq45_candidate import Eq45VelocityCandidate
from .constrained_eq45_supported_candidate import Eq45SupportedVelocityCandidate
from .constrained_eq45_supported_phi10_cubic_temporal_derivative import (
    AUDIT_TIMES,
    REFERENCE_STEP,
    TIME_LEVELS,
    _convergence_report,
    _evaluate_public,
    _fourth_order_time_derivative,
    _fresh_probe_points,
    _manufactured_sign_mutation_calibration,
    _rms,
)
from .constrained_eq45_supported_phi10_cubic_temporal_mode import (
    Eq45SupportedPhi10CubicLocalizedTemporalCandidate,
)
from .constrained_eq45_supported_phi10_quartic_temporal_mode import (
    Eq45SupportedPhi10QuarticDerivativeBalancedTemporalCandidate,
)


SCHEMA = "eq45_supported_phi10_quartic_temporal_derivative_audit_v1"
PROBE_SEED = 914593
PROBE_COUNT = 24
STATIC_RETURN_TIMES = (0.50, 0.625, 0.75)

_TRUTH_BOUNDARY = {
    "serialized_candidate_reloaded": True,
    "public_velocity_only": True,
    "training_loss_reused": False,
    "pressure_fitted": False,
    "forcing_fitted": False,
    "formal_pde_gate_assessed": False,
    "visualization_candidate_only": True,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _reference_derivative(candidate, points: np.ndarray, time: float) -> np.ndarray:
    return _fourth_order_time_derivative(candidate, points, float(time), REFERENCE_STEP)


def _anchor_dynamic_comparisons(quartic, cubic, base, points: np.ndarray) -> dict[str, object]:
    """Measure dynamic differences at identical public-velocity snapshots."""
    early_velocity_difference = (
        _evaluate_public(quartic, points, 0.25) - _evaluate_public(cubic, points, 0.25)
    )
    quartic_early_dt = _reference_derivative(quartic, points, 0.25)
    cubic_early_dt = _reference_derivative(cubic, points, 0.25)
    early_dt_difference = quartic_early_dt - cubic_early_dt

    returns: list[dict[str, float]] = []
    quartic_return_differences: list[np.ndarray] = []
    cubic_return_differences: list[np.ndarray] = []
    for time in STATIC_RETURN_TIMES:
        quartic_velocity_difference = (
            _evaluate_public(quartic, points, time) - _evaluate_public(base, points, time)
        )
        cubic_velocity_difference = (
            _evaluate_public(cubic, points, time) - _evaluate_public(base, points, time)
        )
        quartic_dt_difference = (
            _reference_derivative(quartic, points, time)
            - _reference_derivative(base, points, time)
        )
        cubic_dt_difference = (
            _reference_derivative(cubic, points, time)
            - _reference_derivative(base, points, time)
        )
        quartic_rms = _rms(quartic_dt_difference)
        cubic_rms = _rms(cubic_dt_difference)
        returns.append(
            {
                "time": float(time),
                "quartic_static_velocity_rms_difference": _rms(quartic_velocity_difference),
                "quartic_static_velocity_max_vector_difference": float(
                    np.max(np.linalg.norm(quartic_velocity_difference, axis=1))
                ),
                "cubic_static_velocity_rms_difference": _rms(cubic_velocity_difference),
                "cubic_static_velocity_max_vector_difference": float(
                    np.max(np.linalg.norm(cubic_velocity_difference, axis=1))
                ),
                "quartic_static_time_derivative_rms_difference": quartic_rms,
                "cubic_static_time_derivative_rms_difference": cubic_rms,
                "quartic_to_cubic_dynamic_collateral_ratio": float(
                    quartic_rms / max(cubic_rms, np.finfo(float).tiny)
                ),
            }
        )
        quartic_return_differences.append(quartic_dt_difference)
        cubic_return_differences.append(cubic_dt_difference)

    quartic_stack = np.stack(quartic_return_differences, axis=0)
    cubic_stack = np.stack(cubic_return_differences, axis=0)
    quartic_aggregate = _rms(quartic_stack)
    cubic_aggregate = _rms(cubic_stack)
    return {
        "early_same_snapshot": {
            "time": 0.25,
            "quartic_cubic_velocity_rms_difference": _rms(early_velocity_difference),
            "quartic_cubic_velocity_max_vector_difference": float(
                np.max(np.linalg.norm(early_velocity_difference, axis=1))
            ),
            "quartic_cubic_time_derivative_rms_difference": _rms(early_dt_difference),
            "quartic_cubic_time_derivative_difference_over_cubic_rms": float(
                _rms(early_dt_difference)
                / max(_rms(cubic_early_dt), np.finfo(float).tiny)
            ),
        },
        "static_returns": returns,
        "aggregate_static_return_dynamic_collateral": {
            "quartic_rms": quartic_aggregate,
            "cubic_rms": cubic_aggregate,
            "quartic_to_cubic_ratio": float(
                quartic_aggregate / max(cubic_aggregate, np.finfo(float).tiny)
            ),
        },
    }


def _off_keyframe_derivative_comparisons(quartic, cubic, base, points: np.ndarray) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    for time in (0.375, 0.5625, 0.6875):
        quartic_dt = _reference_derivative(quartic, points, time)
        cubic_dt = _reference_derivative(cubic, points, time)
        base_dt = _reference_derivative(base, points, time)
        rows.append(
            {
                "time": float(time),
                "quartic_cubic_velocity_rms_difference": _rms(
                    _evaluate_public(quartic, points, time)
                    - _evaluate_public(cubic, points, time)
                ),
                "quartic_cubic_time_derivative_rms_difference": _rms(quartic_dt - cubic_dt),
                "quartic_static_time_derivative_rms_difference": _rms(quartic_dt - base_dt),
                "cubic_static_time_derivative_rms_difference": _rms(cubic_dt - base_dt),
            }
        )
    return rows


def audit_reloaded_candidate(
    candidate: Eq45SupportedPhi10QuarticDerivativeBalancedTemporalCandidate,
) -> dict[str, object]:
    if not isinstance(candidate, Eq45SupportedPhi10QuarticDerivativeBalancedTemporalCandidate):
        raise TypeError("candidate must be a reloaded derivative-balanced quartic candidate")
    points, regions = _fresh_probe_points(PROBE_SEED)
    base = candidate.base
    cubic = Eq45SupportedPhi10CubicLocalizedTemporalCandidate(
        base=base, early_delta=float(candidate.early_delta)
    )

    return {
        "schema": SCHEMA,
        "task_id": "CR009-EQ45-SUPPORTED-PHI10-QUARTIC-TEMPORAL-DERIVATIVE-CONVERGENCE-031",
        "candidate_sha256": candidate.sha256,
        "base_sha256": base.sha256,
        "cubic_sha256": cubic.sha256,
        "probe_seed": PROBE_SEED,
        "probe_count": int(points.shape[0]),
        "probe_region_counts": {
            name: int(sum(region == name for region in regions))
            for name in sorted(set(regions))
        },
        "audit_times": [float(value) for value in AUDIT_TIMES],
        "time_levels": [float(value) for value in TIME_LEVELS],
        "reference_step": float(REFERENCE_STEP),
        "quartic_convergence": _convergence_report(candidate, points),
        "cubic_convergence": _convergence_report(cubic, points),
        "base_convergence": _convergence_report(base, points),
        "anchor_dynamic_comparisons": _anchor_dynamic_comparisons(
            candidate, cubic, base, points
        ),
        "off_keyframe_derivative_comparisons": _off_keyframe_derivative_comparisons(
            candidate, cubic, base, points
        ),
        "manufactured_sign_mutation": _manufactured_sign_mutation_calibration(),
        "truth_boundary": dict(_TRUTH_BOUNDARY),
    }


def build_default_report() -> dict[str, object]:
    """Serialize/reload the exact quartic candidate before auditing public output."""
    base = Eq45SupportedVelocityCandidate(parent=Eq45VelocityCandidate.seed())
    candidate = Eq45SupportedPhi10QuarticDerivativeBalancedTemporalCandidate(base=base)
    with tempfile.TemporaryDirectory(prefix="eq45-quartic-ut-") as directory:
        path = Path(directory) / "candidate.json"
        candidate.save_json(path)
        loaded = Eq45SupportedPhi10QuarticDerivativeBalancedTemporalCandidate.load_json(path)
        if loaded.sha256 != candidate.sha256 or loaded.to_dict() != candidate.to_dict():
            raise RuntimeError("serialized/reloaded quartic candidate identity drifted")
        return audit_reloaded_candidate(loaded)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    text = json.dumps(build_default_report(), indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(text, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
        print(args.output)


if __name__ == "__main__":
    main()
