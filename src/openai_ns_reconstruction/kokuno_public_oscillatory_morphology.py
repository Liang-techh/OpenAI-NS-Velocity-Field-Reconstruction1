"""Target-free three-resolution morphology audit for the admitted Kokuno A2 field.

This module does not change the oscillatory velocity candidate.  It consumes the
public provider introduced by Agent 2 #561 and measures velocity/vorticity
statistics on a deterministic off-grid sample cloud.  Vorticity is reconstructed
with an independent Cartesian fourth-order finite-difference curl at three
spatial steps, so the report adds a resolution/stability diagnostic without
reusing the source complete-curl derivatives as an oracle.

The metrics are descriptive only.  They are not a visual correspondence score,
not an NS momentum residual, and not a paper-exact/source-recovery certificate.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Callable

import numpy as np

from .kokuno_public_z_pullback_velocity import velocity_osc

SCHEMA = "kokuno-a2-public-oscillatory-morphology-v1"
PARENT_AGENT2_HEAD = "732800ce4990464b49c8aa32d0dff4580f6684d4"
INDEPENDENT_AGENT4_HEAD = "9b0f86012c53fa8e32a19f766dbc150931870425"
INDEPENDENT_AGENT4_WORKFLOW = 35422203621
INDEPENDENT_AGENT4_ARTIFACT = 10577644750
INDEPENDENT_AGENT4_ARTIFACT_DIGEST = (
    "sha256:7205e041c31d2a68b345235bbd7b0f3316a50ec6c41fa3012c107789c6a07029"
)
FROZEN_FD4_STEPS = (0.02, 0.01, 0.005)

VelocityEvaluator = Callable[[Any, Any, Any, Any], np.ndarray]


def frozen_offgrid_points() -> dict[str, np.ndarray]:
    """Return a deterministic support-interior cloud, independent of A4 points.

    The radii and z values stay at least 0.19 away from the public radial support
    edge and at least 0.71 away from the public axial support edge even after the
    largest +/-2h FD4 displacement.  Angles are irrationally-spaced rather than
    aligned to Cartesian axes.  No point is an optimization/training sample.
    """

    radii = np.asarray((0.42, 0.59, 0.77, 0.95, 1.12), dtype=float)
    z_values = np.asarray((-1.25, -0.52, 0.37, 1.11), dtype=float)
    times = np.asarray((0.31, 0.50, 0.69), dtype=float)
    rows: list[tuple[float, float, float, float]] = []
    golden = math.pi * (3.0 - math.sqrt(5.0))
    index = 0
    for time in times:
        for z in z_values:
            for radius in radii:
                theta = 0.173 + (index + 1) * golden
                rows.append(
                    (
                        radius * math.cos(theta),
                        radius * math.sin(theta),
                        float(z),
                        float(time),
                    )
                )
                index += 1
    array = np.asarray(rows, dtype=float)
    return {
        "x": array[:, 0],
        "y": array[:, 1],
        "z": array[:, 2],
        "t": array[:, 3],
    }


def _as_velocity(values: Any, expected_size: int) -> np.ndarray:
    velocity = np.asarray(values, dtype=float)
    if velocity.shape != (expected_size, 3):
        raise ValueError(
            f"velocity evaluator must return shape {(expected_size, 3)}, got {velocity.shape}"
        )
    if not np.all(np.isfinite(velocity)):
        raise ValueError("velocity evaluator returned non-finite values")
    return velocity


def _fd4_velocity_derivative(
    evaluator: VelocityEvaluator,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
    *,
    axis: int,
    step: float,
) -> np.ndarray:
    """Independent centered FD4 derivative of the public Cartesian velocity."""

    if axis not in (0, 1, 2):
        raise ValueError("axis must be 0, 1, or 2")
    if not math.isfinite(step) or step <= 0.0:
        raise ValueError("step must be positive and finite")
    coordinates = [np.asarray(x, dtype=float), np.asarray(y, dtype=float), np.asarray(z, dtype=float)]
    t = np.asarray(t, dtype=float)
    n = t.size
    if any(value.shape != (n,) for value in coordinates) or t.shape != (n,):
        raise ValueError("x, y, z, t must be aligned one-dimensional sample arrays")

    def shifted(multiplier: float) -> np.ndarray:
        shifted_coordinates = [value.copy() for value in coordinates]
        shifted_coordinates[axis] = shifted_coordinates[axis] + multiplier * step
        return _as_velocity(
            evaluator(
                shifted_coordinates[0],
                shifted_coordinates[1],
                shifted_coordinates[2],
                t,
            ),
            n,
        )

    return (
        -shifted(2.0)
        + 8.0 * shifted(1.0)
        - 8.0 * shifted(-1.0)
        + shifted(-2.0)
    ) / (12.0 * step)


def fd4_vorticity(
    evaluator: VelocityEvaluator,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    *,
    step: float,
) -> np.ndarray:
    """Compute Cartesian curl with an evaluator-only FD4 oracle."""

    arrays = [np.asarray(value, dtype=float) for value in (x, y, z, t)]
    if not arrays or arrays[0].ndim != 1:
        raise ValueError("x, y, z, t must be one-dimensional")
    n = arrays[0].size
    if any(value.shape != (n,) or not np.all(np.isfinite(value)) for value in arrays):
        raise ValueError("x, y, z, t must be finite aligned one-dimensional arrays")
    d_dx = _fd4_velocity_derivative(evaluator, *arrays, axis=0, step=step)
    d_dy = _fd4_velocity_derivative(evaluator, *arrays, axis=1, step=step)
    d_dz = _fd4_velocity_derivative(evaluator, *arrays, axis=2, step=step)
    return np.column_stack(
        (
            d_dy[:, 2] - d_dz[:, 1],
            d_dz[:, 0] - d_dx[:, 2],
            d_dx[:, 1] - d_dy[:, 0],
        )
    )


def _rms(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    return float(math.sqrt(float(np.mean(values * values))))


def _vector_rms(vectors: np.ndarray) -> float:
    vectors = np.asarray(vectors, dtype=float)
    return float(math.sqrt(float(np.mean(np.sum(vectors * vectors, axis=-1)))))


def _morphology_metrics(velocity: np.ndarray, vorticity: np.ndarray) -> dict[str, Any]:
    speed = np.linalg.norm(velocity, axis=1)
    omega_norm = np.linalg.norm(vorticity, axis=1)
    total_enstrophy = float(np.mean(omega_norm * omega_norm))
    if not math.isfinite(total_enstrophy) or total_enstrophy <= 0.0:
        raise RuntimeError("public field has zero/non-finite sampled vorticity")
    component_rms = [_rms(vorticity[:, index]) for index in range(3)]
    axial_enstrophy = float(np.mean(vorticity[:, 2] ** 2))
    radial_enstrophy = float(np.mean(vorticity[:, 0] ** 2 + vorticity[:, 1] ** 2))
    nonzero = omega_norm > np.finfo(float).tiny
    axis_alignment = np.zeros_like(omega_norm)
    axis_alignment[nonzero] = np.abs(vorticity[nonzero, 2]) / omega_norm[nonzero]
    return {
        "velocity_rms": _vector_rms(velocity),
        "velocity_speed_max": float(np.max(speed)),
        "vorticity_rms": _vector_rms(vorticity),
        "vorticity_abs_max": float(np.max(omega_norm)),
        "vorticity_component_rms_xyz": component_rms,
        "axial_enstrophy_fraction": axial_enstrophy / total_enstrophy,
        "radial_enstrophy_fraction": radial_enstrophy / total_enstrophy,
        "mean_abs_vorticity_axis_alignment": float(np.mean(axis_alignment)),
    }


def build_morphology_receipt(
    evaluator: VelocityEvaluator = velocity_osc,
    *,
    steps: tuple[float, float, float] = FROZEN_FD4_STEPS,
) -> dict[str, Any]:
    """Evaluate one target-free three-resolution morphology/stability receipt."""

    steps = tuple(float(step) for step in steps)
    if len(steps) != 3 or any(not math.isfinite(step) or step <= 0.0 for step in steps):
        raise ValueError("exactly three positive finite differentiation steps are required")
    if not (steps[0] > steps[1] > steps[2]):
        raise ValueError("steps must be strictly decreasing coarse -> fine")

    points = frozen_offgrid_points()
    x, y, z, t = (points[key] for key in ("x", "y", "z", "t"))
    velocity = _as_velocity(evaluator(x, y, z, t), x.size)
    vorticities: list[np.ndarray] = []
    by_step: dict[str, Any] = {}
    for step in steps:
        omega = fd4_vorticity(evaluator, x, y, z, t, step=step)
        if not np.all(np.isfinite(omega)):
            raise RuntimeError("FD4 vorticity contains non-finite values")
        vorticities.append(omega)
        by_step[f"{step:.6g}"] = _morphology_metrics(velocity, omega)

    coarse_medium = _vector_rms(vorticities[0] - vorticities[1])
    medium_fine = _vector_rms(vorticities[1] - vorticities[2])
    fine_scale = max(_vector_rms(vorticities[2]), np.finfo(float).tiny)
    pairwise = {
        "coarse_to_medium_relative_vorticity_rms_change": coarse_medium
        / max(_vector_rms(vorticities[1]), np.finfo(float).tiny),
        "medium_to_fine_relative_vorticity_rms_change": medium_fine / fine_scale,
        "coarse_medium_to_medium_fine_change_ratio": coarse_medium
        / max(medium_fine, np.finfo(float).tiny),
    }
    finest_metrics = by_step[f"{steps[-1]:.6g}"]

    return {
        "schema": SCHEMA,
        "candidate": {
            "provider": "openai_ns_reconstruction.kokuno_public_z_pullback_velocity:velocity_osc",
            "agent2_head": PARENT_AGENT2_HEAD,
            "candidate_provenance": "repository-autonomous/source-compatible",
            "paper_exact": False,
        },
        "independent_agent4_preflight_reference": {
            "agent4_head": INDEPENDENT_AGENT4_HEAD,
            "workflow": INDEPENDENT_AGENT4_WORKFLOW,
            "artifact_id": INDEPENDENT_AGENT4_ARTIFACT,
            "artifact_digest": INDEPENDENT_AGENT4_ARTIFACT_DIGEST,
            "public_oscillatory_preflight_passed": True,
            "note": "external independent receipt; not recomputed or promoted by this morphology module",
        },
        "protocol": {
            "sample_count": int(x.size),
            "sample_design": "deterministic off-grid support-interior cloud independent of Agent-4 held-out points",
            "times": sorted({float(value) for value in t}),
            "fd4_steps": list(steps),
            "derivative_oracle": "independent Cartesian centered FD4 using only public velocity evaluations",
            "target_image_or_visual_score_used": False,
        },
        "results_by_step": by_step,
        "offgrid_resolution_stability": pairwise,
        "finest_summary": {
            "vorticity_rms": finest_metrics["vorticity_rms"],
            "vorticity_abs_max": finest_metrics["vorticity_abs_max"],
            "axial_enstrophy_fraction": finest_metrics["axial_enstrophy_fraction"],
            "radial_enstrophy_fraction": finest_metrics["radial_enstrophy_fraction"],
            "mean_abs_vorticity_axis_alignment": finest_metrics[
                "mean_abs_vorticity_axis_alignment"
            ],
        },
        "truth_boundary": {
            "velocity_candidate_changed": False,
            "source_formula_changed": False,
            "agent4_scientific_guards_changed": False,
            "three_resolution_vorticity_morphology_assessed": True,
            "offgrid_vorticity_stability_assessed": True,
            "morphology_acceptance_threshold_preregistered": False,
            "visual_correspondence_assessed": False,
            "heldout_ns_momentum_residual_assessed": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
        },
    }


def _main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    receipt = build_morphology_receipt()
    text = json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(text, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
        print(text, end="")


if __name__ == "__main__":
    _main()
