"""Small held-out post-core diagnostic for the native-X complete-curl adapter.

K2-OSC-005 showed that the retained Cartesian-box curl has essentially no
activity on Agent-1's reference-continuation samples.  This diagnostic changes
only the localization coordinate: the new correction uses the predeclared
native-X annulus ``X0<X<2*X0``.  Amplitude, affine phase, wave vector,
polarization and frequency remain frozen at the retained Agent-2 values.

The report is deliberately small and is not used for parameter selection.  It
compares one post-core held-out slice with Agent-1 pressure unchanged, nu=.01,
and zero forcing.  Both leading and composite residuals are evaluated by the
same independent FD4 black-box operator so this is a preliminary engineering
check, not the formal full-domain normalized PDE gate.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_complete_curl import KokunoCompleteCurlCorrection
from .kokuno_independent_leading_core import evaluate_fd4
from .kokuno_native_x_curl import KokunoNativeXCompleteCurlCorrection
from .kokuno_reference_continuation import KokunoReferenceContinuationCandidate
from .kokuno_similarity_coordinates import KokunoNativeSimilarityCoordinates


SCHEMA = "kokuno-native-x-complete-curl-screen-v1"
TASK_ID = "K2-OSC-006"
VALIDATION_SEED = 9_173_001
VALIDATION_TIME = 0.5
VALIDATION_STEP = 0.01
NU = 0.01
AMPLITUDE = 0.125
PHASE = 0.0


def _sample_post_core_native_x(*, count: int, seed: int, time: float) -> np.ndarray:
    if isinstance(count, bool) or not isinstance(count, (int, np.integer)):
        raise TypeError("count must be an integer")
    count = int(count)
    if count < 4:
        raise ValueError("count must be at least 4")
    rng = np.random.default_rng(int(seed))
    coordinates = KokunoNativeSimilarityCoordinates()
    z = rng.uniform(-0.02, 0.02, count)
    target_X = rng.uniform(0.52, 0.68, count)
    angle = rng.uniform(-np.pi, np.pi, count)
    q = coordinates.solve_q(z, float(time))
    radius = np.sqrt(2.0 * q * target_X)
    return np.column_stack((radius * np.cos(angle), radius * np.sin(angle), z))


def _metrics(result: dict[str, np.ndarray]) -> dict[str, float]:
    residual_norm = np.linalg.norm(np.asarray(result["residual"], dtype=float), axis=1)
    divergence = np.asarray(result["divergence"], dtype=float)
    speed = np.linalg.norm(np.asarray(result["velocity"], dtype=float), axis=1)
    return {
        "residual_rms": float(np.sqrt(np.mean(residual_norm * residual_norm))),
        "residual_max": float(np.max(residual_norm)),
        "divergence_rms": float(np.sqrt(np.mean(divergence * divergence))),
        "divergence_max": float(np.max(np.abs(divergence))),
        "velocity_rms": float(np.sqrt(np.mean(speed * speed))),
        "velocity_max": float(np.max(speed)),
    }


def _rms_speed(values: np.ndarray) -> float:
    norm = np.linalg.norm(np.asarray(values, dtype=float), axis=1)
    return float(np.sqrt(np.mean(norm * norm)))


def run_screen(
    *,
    count: int = 4,
    seed: int = VALIDATION_SEED,
    time: float = VALIDATION_TIME,
    step: float = VALIDATION_STEP,
) -> dict[str, Any]:
    if not np.isfinite(time) or not (0.25 <= float(time) <= 0.75):
        raise ValueError("time must be finite and lie in [0.25,0.75]")
    if not np.isfinite(step) or not (0.0 < float(step) <= 0.02):
        raise ValueError("step must be finite and lie in (0,.02]")

    leading = KokunoReferenceContinuationCandidate()
    native = KokunoNativeXCompleteCurlCorrection(
        amplitude=AMPLITUDE,
        h=leading.h,
        Lambda=leading.Lambda,
        phase=PHASE,
    )
    old_box = KokunoCompleteCurlCorrection(amplitude=AMPLITUDE, phase=PHASE)
    points = _sample_post_core_native_x(count=count, seed=seed, time=float(time))
    source_coordinates = leading.coordinates(
        points[:, 0], points[:, 1], points[:, 2], float(time)
    )
    point_X = np.asarray(source_coordinates["X"], dtype=float)

    native_velocity = native.at_points(points, float(time))
    old_velocity = old_box.at_points(points, float(time))
    composite_velocity = native.compose_velocity(leading.velocity)

    leading_result = evaluate_fd4(
        leading.velocity,
        leading.pressure,
        points,
        float(time),
        float(step),
        nu=NU,
    )
    composite_result = evaluate_fd4(
        composite_velocity,
        leading.pressure,
        points,
        float(time),
        float(step),
        nu=NU,
    )
    leading_metrics = _metrics(leading_result)
    composite_metrics = _metrics(composite_result)
    native_rms = _rms_speed(native_velocity)
    old_rms = _rms_speed(old_velocity)

    try:
        leading.core.at_points(points, float(time))
    except (ValueError, RuntimeError):
        old_core_rejected = True
    else:
        old_core_rejected = False

    residual_ratio = float(
        composite_metrics["residual_rms"] / leading_metrics["residual_rms"]
    )
    activity_ratio = float(native_rms / old_rms) if old_rms > 0.0 else float("inf")

    return {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "purpose": (
            "one frozen post-core engineering screen of source-native-X localization "
            "at the vector-potential level; no parameter selection"
        ),
        "source_stage": {
            "candidate_schema": leading.to_payload()["schema"],
            "candidate_sha256": leading.sha256,
            "source_commit": leading.to_payload()["source"]["commit"],
            "corrected_release": leading.to_payload()["source"]["corrected_release"],
            "reference_continuation_is_final_leading": False,
        },
        "localization": {
            "source_native_coordinate": "X=(x^2+y^2)/(2q)",
            "source_reference_start": "X0=4/Lambda",
            "autonomous_support_rule": "X0<X<2*X0",
            "X_support": list(native.X_support),
            "source_exact_support_claimed": False,
            "potential_level_localization": True,
            "complete_curl_remainder_retained": True,
        },
        "oscillation": {
            "amplitude": AMPLITUDE,
            "phase": PHASE,
            "wave_vector": list(native.wave_vector),
            "polarization": list(native.polarization),
            "omega": native.omega,
            "affine_cartesian_phase_frame_is_autonomous": True,
            "parameter_selection_performed": False,
        },
        "validation": {
            "seed": int(seed),
            "points": int(count),
            "time": float(time),
            "step": float(step),
            "nu": NU,
            "X_min": float(np.min(point_X)),
            "X_max": float(np.max(point_X)),
            "old_core_rejected": old_core_rejected,
            "operator": "same independent fourth-order centered Cartesian FD4 for leading and composite",
            "pressure": "Agent-1 reference-continuation pressure unchanged",
            "forcing": "zero raw diagnostic only; no fit",
            "training_samples_used": False,
            "validation_used_for_parameter_selection": False,
        },
        "activity": {
            "old_cartesian_box_correction_velocity_rms": old_rms,
            "native_x_correction_velocity_rms": native_rms,
            "native_to_old_activity_ratio": activity_ratio,
        },
        "leading": leading_metrics,
        "leading_plus_native_x_curl": composite_metrics,
        "residual_rms_ratio": residual_ratio,
        "routing": {
            "candidate_promoted": False,
            "routing_changed": False,
            "interpretation": (
                "diagnostic only: native-X localization tests support transfer; "
                "the autonomous support interval and affine phase are not fitted here"
            ),
        },
        "truth_boundary": {
            "source_native_localization_coordinate": True,
            "complete_curl_correction_executable": True,
            "source_exact_support": False,
            "source_exact_phase_frame": False,
            "global_leading_profile_reconstructed": False,
            "formal_full_domain_pde_gate_assessed": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
        },
    }


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--count", type=int, default=4)
    parser.add_argument("--seed", type=int, default=VALIDATION_SEED)
    parser.add_argument("--time", type=float, default=VALIDATION_TIME)
    parser.add_argument("--step", type=float, default=VALIDATION_STEP)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    report = run_screen(count=args.count, seed=args.seed, time=args.time, step=args.step)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
