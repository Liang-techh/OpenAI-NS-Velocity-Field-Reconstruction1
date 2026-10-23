"""Truth-bounded oscillatory screen on Kokuno's executable reference continuation.

This Agent-2 increment deliberately does not add an oscillatory degree of freedom.
It freezes the complete-curl surrogate already retained by the routed core
checkpoint (amplitude=0.125, phase=0) and asks one narrower question: does its
small core-only momentum improvement survive Agent-1's first source-derived
outer/reference continuation stage?

The leading stage is ``KokunoReferenceContinuationCandidate``.  It implements
Kokuno's public flat logarithmic-X reference continuation but intentionally stops
before cone modulation, moment restoration, heat compensation and the final
matched core-to-heat splice.  Consequently every residual below is a local
engineering diagnostic, not the project's formal full-domain PDE gate.

The oscillatory correction remains the autonomous Cartesian complete-curl
surrogate from ``kokuno_complete_curl.py``.  Source-derived structure is limited
to exact-curl/potential-level localization; its Cartesian phase, frame, compact
box support and numerical wave parameters are not Kokuno/OpenAI exact data.

For cost and numerical separation, the leading residual/Jacobian are measured by
the independent FD4 black-box operator, while the already independently checked
Agent-2 analytic ``w_t``, ``grad w`` and ``Delta w`` are used in the exact
momentum increment

    N(b+w)-N(b) = w_t + (b.grad)w + (w.grad)b + (w.grad)w - nu Delta w.

This avoids numerically differentiating the high-frequency correction and avoids
evaluating the expensive reference-continuation leading field twice per stencil.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_complete_curl import KokunoCompleteCurlCorrection
from .kokuno_independent_leading_core import evaluate_fd4
from .kokuno_reference_continuation import KokunoReferenceContinuationCandidate


SCHEMA = "kokuno-reference-continuation-oscillatory-screen-v1"
TASK_ID = "K2-OSC-005"
VALIDATION_SEED = 9_172_981
DEFAULT_TIMES = (0.375, 0.5, 0.625)
DEFAULT_STEPS = (0.02, 0.01, 0.005)
NU = 0.01
AMPLITUDE = 0.125
PHASE = 0.0
REGISTERED_PDE_THRESHOLD = 1.0e-3
REGISTERED_DIVERGENCE_THRESHOLD = 1.0e-5


def _sample_points(*, count: int, seed: int) -> dict[str, np.ndarray]:
    """Return disjoint inner and post-core-X cylindrical strata.

    The outer radial interval is chosen so its unshifted points are beyond the
    old finite core's ``Lambda*X<=4.1`` domain at all registered report times.
    FD stencils are allowed to cross that historical seam because the new
    reference continuation is precisely the public stage being tested.
    """

    if isinstance(count, bool) or not isinstance(count, (int, np.integer)):
        raise TypeError("count must be an integer")
    count = int(count)
    if count < 4:
        raise ValueError("count must be at least 4 per stratum")
    rng = np.random.default_rng(int(seed))

    def ring(radius_low: float, radius_high: float, z_half: float) -> np.ndarray:
        radius = rng.uniform(radius_low, radius_high, count)
        angle = rng.uniform(-np.pi, np.pi, count)
        z = rng.uniform(-z_half, z_half, count)
        return np.column_stack((radius * np.cos(angle), radius * np.sin(angle), z))

    inner = ring(0.10, 0.24, 0.08)
    post_core = ring(0.74, 0.79, 0.04)
    return {"inner": inner, "post_core": post_core}


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


def _ratio(numerator: float, denominator: float) -> float:
    if denominator <= 0.0 or not np.isfinite(denominator):
        raise ArithmeticError("diagnostic denominator must be positive and finite")
    return float(numerator / denominator)


def _old_core_rejects(
    leading: KokunoReferenceContinuationCandidate,
    points: np.ndarray,
    times: tuple[float, ...],
) -> bool:
    for time in times:
        try:
            leading.core.at_points(points, float(time))
        except (ValueError, RuntimeError):
            continue
        return False
    return True


def _analytic_composite_from_base(
    base: dict[str, np.ndarray],
    *,
    correction_velocity: np.ndarray,
    correction_time_derivative: np.ndarray,
    correction_jacobian: np.ndarray,
    correction_laplacian: np.ndarray,
    correction_divergence: np.ndarray,
    nu: float,
) -> tuple[dict[str, np.ndarray], np.ndarray]:
    """Add the exact complete-curl momentum increment to one FD4 base result."""

    b = np.asarray(base["velocity"], dtype=float)
    jac_b = np.asarray(base["jacobian"], dtype=float)
    w = np.asarray(correction_velocity, dtype=float)
    jac_w = np.asarray(correction_jacobian, dtype=float)

    b_grad_w = np.einsum("nij,nj->ni", jac_w, b)
    w_grad_b = np.einsum("nij,nj->ni", jac_b, w)
    w_grad_w = np.einsum("nij,nj->ni", jac_w, w)
    increment = (
        np.asarray(correction_time_derivative, dtype=float)
        + b_grad_w
        + w_grad_b
        + w_grad_w
        - float(nu) * np.asarray(correction_laplacian, dtype=float)
    )
    composite = {
        "velocity": b + w,
        "divergence": np.asarray(base["divergence"], dtype=float)
        + np.asarray(correction_divergence, dtype=float),
        "residual": np.asarray(base["residual"], dtype=float) + increment,
    }
    return composite, increment


def run_screen(
    *,
    count: int = 8,
    seed: int = VALIDATION_SEED,
    times: tuple[float, ...] = DEFAULT_TIMES,
    steps: tuple[float, ...] = DEFAULT_STEPS,
) -> dict[str, Any]:
    """Compare reference-leading versus reference-leading + frozen complete curl.

    No quantity returned by this function is used to optimize or select an
    oscillatory parameter.  The amplitude, phase, wave vector, polarization,
    frequency and compact support are frozen before the held-out evaluation.
    """

    times = tuple(float(value) for value in times)
    steps = tuple(float(value) for value in steps)
    if not times or not steps:
        raise ValueError("times and steps must be nonempty")
    if any(not np.isfinite(value) for value in times + steps):
        raise ValueError("times and steps must be finite")
    if any(value <= 0.0 for value in steps):
        raise ValueError("steps must be positive")

    leading = KokunoReferenceContinuationCandidate()
    correction = KokunoCompleteCurlCorrection(amplitude=AMPLITUDE, phase=PHASE)
    samples = _sample_points(count=count, seed=seed)

    rows: list[dict[str, Any]] = []
    for stratum, points in samples.items():
        for time in times:
            coordinates = leading.coordinates(
                points[:, 0], points[:, 1], points[:, 2], float(time)
            )
            point_X = np.asarray(coordinates["X"], dtype=float)
            w = np.asarray(correction.at_points(points, float(time)), dtype=float)
            w_t = np.asarray(
                correction.time_derivative(
                    points[:, 0], points[:, 1], points[:, 2], float(time)
                ),
                dtype=float,
            )
            jac_w = np.asarray(
                correction.spatial_jacobian(
                    points[:, 0], points[:, 1], points[:, 2], float(time)
                ),
                dtype=float,
            )
            lap_w = np.asarray(
                correction.laplacian(
                    points[:, 0], points[:, 1], points[:, 2], float(time)
                ),
                dtype=float,
            )
            div_w = np.asarray(
                correction.divergence(
                    points[:, 0], points[:, 1], points[:, 2], float(time)
                ),
                dtype=float,
            )
            correction_speed = np.linalg.norm(w, axis=1)

            for step in steps:
                base = evaluate_fd4(
                    leading.velocity,
                    leading.pressure,
                    points,
                    float(time),
                    float(step),
                    nu=NU,
                )
                composite, increment = _analytic_composite_from_base(
                    base,
                    correction_velocity=w,
                    correction_time_derivative=w_t,
                    correction_jacobian=jac_w,
                    correction_laplacian=lap_w,
                    correction_divergence=div_w,
                    nu=NU,
                )
                base_metrics = _metrics(base)
                composite_metrics = _metrics(composite)
                increment_norm = np.linalg.norm(increment, axis=1)
                rows.append(
                    {
                        "stratum": stratum,
                        "time": float(time),
                        "step": float(step),
                        "X_min": float(np.min(point_X)),
                        "X_max": float(np.max(point_X)),
                        "correction_velocity_rms": float(
                            np.sqrt(np.mean(correction_speed * correction_speed))
                        ),
                        "correction_velocity_max": float(np.max(correction_speed)),
                        "analytic_correction_divergence_max": float(np.max(np.abs(div_w))),
                        "analytic_momentum_increment_rms": float(
                            np.sqrt(np.mean(increment_norm * increment_norm))
                        ),
                        "leading": base_metrics,
                        "leading_plus_oscillation": composite_metrics,
                        "residual_rms_ratio": _ratio(
                            composite_metrics["residual_rms"],
                            base_metrics["residual_rms"],
                        ),
                        "residual_max_ratio": _ratio(
                            composite_metrics["residual_max"],
                            base_metrics["residual_max"],
                        ),
                    }
                )

    finest = min(steps)
    finest_rows = [row for row in rows if row["step"] == finest]
    summaries: dict[str, dict[str, float]] = {}
    for stratum in samples:
        selected = [row for row in finest_rows if row["stratum"] == stratum]
        base_mean = float(np.mean([row["leading"]["residual_rms"] for row in selected]))
        composite_mean = float(
            np.mean([row["leading_plus_oscillation"]["residual_rms"] for row in selected])
        )
        summaries[stratum] = {
            "mean_leading_residual_rms": base_mean,
            "mean_composite_residual_rms": composite_mean,
            "mean_residual_rms_ratio": _ratio(composite_mean, base_mean),
            "worst_leading_residual_max": float(
                max(row["leading"]["residual_max"] for row in selected)
            ),
            "worst_composite_residual_max": float(
                max(row["leading_plus_oscillation"]["residual_max"] for row in selected)
            ),
            "worst_composite_divergence_max": float(
                max(row["leading_plus_oscillation"]["divergence_max"] for row in selected)
            ),
            "max_analytic_correction_divergence": float(
                max(row["analytic_correction_divergence_max"] for row in selected)
            ),
        }

    all_base = float(np.mean([row["leading"]["residual_rms"] for row in finest_rows]))
    all_composite = float(
        np.mean([row["leading_plus_oscillation"]["residual_rms"] for row in finest_rows])
    )

    post_core_old_core_rejected = _old_core_rejects(
        leading, samples["post_core"], times
    )
    post_core_X_min = float(
        min(
            row["X_min"]
            for row in finest_rows
            if row["stratum"] == "post_core"
        )
    )

    report: dict[str, Any] = {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "purpose": (
            "held-out transfer check of the retained Agent-2 exact-curl surrogate "
            "from the finite leading core to Agent-1's first source reference-continuation stage"
        ),
        "source_stage": {
            "candidate_schema": leading.to_payload()["schema"],
            "candidate_sha256": leading.sha256,
            "source_commit": leading.to_payload()["source"]["commit"],
            "corrected_release": leading.to_payload()["source"]["corrected_release"],
            "reference_continuation_is_final_leading": False,
            "cone_moment_heat_matching_complete": False,
        },
        "oscillation": {
            "candidate_class": "KokunoCompleteCurlCorrection",
            "amplitude": AMPLITUDE,
            "phase": PHASE,
            "wave_vector": list(correction.wave_vector),
            "polarization": list(correction.polarization),
            "omega": correction.omega,
            "center": list(correction.center),
            "half_widths": list(correction.half_widths),
            "source_derived_contract": "complete curl with potential-level localization",
            "phase_frame_support_parameters": "autonomous_cartesian_surrogate",
            "parameter_selection_performed": False,
        },
        "validation": {
            "seed": int(seed),
            "points_per_stratum": int(count),
            "strata": list(samples),
            "times": list(times),
            "steps": list(steps),
            "nu": NU,
            "leading_operator": "independent_fourth_order_centered_cartesian_fd4",
            "oscillatory_derivatives": "Agent-2 analytic complete-curl derivatives checked independently in K2-OSC-002/003",
            "composite_residual": "FD4 leading residual plus exact analytic complete-curl momentum increment",
            "pressure": "Agent-1 reference-continuation pressure unchanged",
            "forcing": "zero_raw_local_diagnostic_only_no_fit",
            "training_samples_used": False,
            "validation_used_for_parameter_selection": False,
            "post_core_stratum_old_core_rejected": post_core_old_core_rejected,
            "post_core_X_min": post_core_X_min,
            "old_core_nominal_X_limit": 4.1 / leading.Lambda,
        },
        "rows": rows,
        "finest_step": finest,
        "finest_summaries": summaries,
        "finest_all_strata": {
            "mean_leading_residual_rms": all_base,
            "mean_composite_residual_rms": all_composite,
            "mean_residual_rms_ratio": _ratio(all_composite, all_base),
        },
        "routing": {
            "previous_core_routing": {"amplitude": AMPLITUDE, "phase": PHASE},
            "routing_changed_by_this_screen": False,
            "candidate_promoted": False,
            "interpretation": (
                "transfer diagnostic only; result may support or weaken the retained probe, "
                "but no parameter is selected on these held-out points"
            ),
        },
        "fixed_gates": {
            "normalized_full_momentum_residual": REGISTERED_PDE_THRESHOLD,
            "divergence_max": REGISTERED_DIVERGENCE_THRESHOLD,
            "changed": False,
        },
        "truth_boundary": {
            "complete_curl_correction_executable": True,
            "reference_continuation_executable": True,
            "post_core_reference_stage_audited": True,
            "global_leading_profile_reconstructed": False,
            "core_to_heat_matching_completed": False,
            "mean_correction_applied": False,
            "pressure_or_forcing_fitted": False,
            "formal_full_domain_pde_gate_assessed": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
    }
    if not post_core_old_core_rejected:
        raise AssertionError("post-core stratum failed to exercise the new continuation domain")
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--count", type=int, default=8)
    parser.add_argument("--seed", type=int, default=VALIDATION_SEED)
    args = parser.parse_args(argv)
    report = run_screen(count=args.count, seed=args.seed)
    text = json.dumps(report, indent=2, sort_keys=True, allow_nan=False)
    if args.output is None:
        print(text)
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
