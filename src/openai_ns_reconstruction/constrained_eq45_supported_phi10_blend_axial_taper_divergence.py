"""Independent public-output divergence audit for the bounded axial-taper family.

This CR006/CR009 increment validates the new
``Eq45SupportedPhi10BlendAxialTaperCandidate`` without using any internal
streamfunction derivatives.  Each frozen taper member is serialized, reloaded,
and sampled only through the public ``at_points(...)->[u,v,w]`` interface.
Cartesian first derivatives are reconstructed independently with centered finite
differences on three fixed spatial steps.

The report is intentionally a finite held-out convergence diagnostic rather than
the preregistered full-domain PDE acceptance gate.  It asks a narrower question:
does moving the already-governed axial identity plateau introduce a persistent
incompressibility defect, or does the observed divergence continue to decay like
finite-difference truncation error?  No taper value is selected here.
"""
from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
import tempfile
from typing import Any

import numpy as np

from .constrained_eq45_supported_collar_vorticity import governed_supported_seed
from .constrained_eq45_supported_phi10_blend_axial_taper_mode import (
    Eq45SupportedPhi10BlendAxialTaperCandidate,
)
from .constrained_eq45_supported_phi10_compact_quartic_blend_temporal_mode import (
    Eq45SupportedPhi10CompactQuarticBlendTemporalCandidate,
)


TASK_ID = "CR006-EQ45-SUPPORTED-BLEND-AXIAL-TAPER-DIVERGENCE-CONVERGENCE-035"
SEED = 914641
BLEND_WEIGHT = 0.5
AXIAL_PLATEAU_VALUES = (0.64, 0.725, 0.81)
TIMES = (0.3125, 0.6875)
SPATIAL_STEPS = (0.02, 0.01, 0.005)
MUTATION_DIVERGENCE = 0.03
REGION_SIZE = 8
REGIONS = ("plateau", "radial_collar", "axial_collar", "corner_collar")


@dataclass(frozen=True)
class _PublicVelocityMutation:
    """Calibration wrapper that adds a known +epsilon divergence to public u."""

    base: Eq45SupportedPhi10BlendAxialTaperCandidate
    epsilon: float = MUTATION_DIVERGENCE

    def at_points(self, points, time):
        points_arr = np.asarray(points, dtype=float)
        values = np.asarray(self.base.at_points(points_arr, time), dtype=float).copy()
        values[..., 0] += float(self.epsilon) * points_arr[..., 0]
        return values


def governed_candidate(axial_plateau_q: float) -> Eq45SupportedPhi10BlendAxialTaperCandidate:
    """Construct one frozen member of the bounded taper family."""

    blend = Eq45SupportedPhi10CompactQuarticBlendTemporalCandidate(
        base=governed_supported_seed(),
        blend_weight=BLEND_WEIGHT,
    )
    return Eq45SupportedPhi10BlendAxialTaperCandidate(
        base=blend,
        axial_plateau_q=float(axial_plateau_q),
    )


def _polar_points(rng: np.random.Generator, radii: np.ndarray, z: np.ndarray) -> np.ndarray:
    theta = rng.uniform(-np.pi, np.pi, size=radii.shape[0])
    return np.column_stack((radii * np.cos(theta), radii * np.sin(theta), z))


def fresh_stratified_points(seed: int = SEED) -> tuple[np.ndarray, np.ndarray]:
    """Return fresh deterministic probes buffered away from support faces/knots."""

    rng = np.random.default_rng(int(seed))
    chunks: list[np.ndarray] = []
    labels: list[str] = []

    # Interior plateau: comfortably inside every axial/radial identity region.
    r = rng.uniform(0.20, 1.35, size=REGION_SIZE)
    z = rng.uniform(-1.35, 1.35, size=REGION_SIZE)
    chunks.append(_polar_points(rng, r, z))
    labels.extend(["plateau"] * REGION_SIZE)

    # Radial-only collar: z stays in the common axial plateau.
    r = rng.uniform(1.70, 1.88, size=REGION_SIZE)
    z = rng.uniform(-1.20, 1.20, size=REGION_SIZE)
    chunks.append(_polar_points(rng, r, z))
    labels.extend(["radial_collar"] * REGION_SIZE)

    # Axial collar: |z| > 1.84 is inside the taper region for every q in the
    # audited family, while remaining at least 0.06 from the outer support face.
    r = rng.uniform(0.25, 1.30, size=REGION_SIZE)
    abs_z = rng.uniform(1.84, 1.92, size=REGION_SIZE)
    signs = rng.choice(np.asarray([-1.0, 1.0]), size=REGION_SIZE)
    chunks.append(_polar_points(rng, r, signs * abs_z))
    labels.extend(["axial_collar"] * REGION_SIZE)

    # Combined radial+axial collar.
    r = rng.uniform(1.70, 1.88, size=REGION_SIZE)
    abs_z = rng.uniform(1.84, 1.92, size=REGION_SIZE)
    signs = rng.choice(np.asarray([-1.0, 1.0]), size=REGION_SIZE)
    chunks.append(_polar_points(rng, r, signs * abs_z))
    labels.extend(["corner_collar"] * REGION_SIZE)

    points = np.concatenate(chunks, axis=0)
    label_arr = np.asarray(labels, dtype="U32")
    if points.shape != (REGION_SIZE * len(REGIONS), 3):
        raise RuntimeError("unexpected stratified probe shape")
    if not np.all(np.isfinite(points)):
        raise RuntimeError("stratified probes became nonfinite")
    return points, label_arr


def _gradient(field, points: np.ndarray, time: float, step: float) -> np.ndarray:
    """Independent centered Cartesian gradient of the public velocity output."""

    h = float(step)
    if not np.isfinite(h) or h <= 0.0:
        raise ValueError("step must be positive and finite")
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError("points must have shape (n,3)")

    gradient = np.empty((points.shape[0], 3, 3), dtype=float)
    for axis in range(3):
        delta = np.zeros(3, dtype=float)
        delta[axis] = h
        plus = np.asarray(field.at_points(points + delta, float(time)), dtype=float)
        minus = np.asarray(field.at_points(points - delta, float(time)), dtype=float)
        gradient[:, :, axis] = (plus - minus) / (2.0 * h)
    if not np.all(np.isfinite(gradient)):
        raise RuntimeError("finite-difference gradient became nonfinite")
    return gradient


def _metrics(gradient: np.ndarray) -> dict[str, float]:
    divergence = gradient[:, 0, 0] + gradient[:, 1, 1] + gradient[:, 2, 2]
    div_rms = float(np.sqrt(np.mean(divergence * divergence)))
    grad_rms = float(np.sqrt(np.mean(np.sum(gradient * gradient, axis=(1, 2)))))
    return {
        "sampled_max_abs": float(np.max(np.abs(divergence))),
        "sampled_rms": div_rms,
        "gradient_rms": grad_rms,
        "gradient_normalized_rms": float(div_rms / max(grad_rms, np.finfo(float).tiny)),
    }


def _region_metrics(gradient: np.ndarray, labels: np.ndarray) -> dict[str, dict[str, float]]:
    result: dict[str, dict[str, float]] = {}
    for region in REGIONS:
        mask = labels == region
        result[region] = _metrics(gradient[mask])
    return result


def _observed_orders(levels: list[dict[str, Any]]) -> list[float | None]:
    orders: list[float | None] = []
    for coarse, fine in zip(levels[:-1], levels[1:]):
        coarse_error = float(coarse["sampled_rms"])
        fine_error = float(fine["sampled_rms"])
        if coarse_error <= 0.0 or fine_error <= 0.0:
            orders.append(None)
            continue
        orders.append(
            float(
                np.log(coarse_error / fine_error)
                / np.log(float(coarse["step"]) / float(fine["step"]))
            )
        )
    return orders


def _velocity_response(
    baseline: Eq45SupportedPhi10BlendAxialTaperCandidate,
    candidate: Eq45SupportedPhi10BlendAxialTaperCandidate,
    points: np.ndarray,
    time: float,
) -> dict[str, float]:
    reference = np.asarray(baseline.at_points(points, time), dtype=float)
    values = np.asarray(candidate.at_points(points, time), dtype=float)
    delta = values - reference
    reference_rms = float(np.sqrt(np.mean(np.sum(reference * reference, axis=1))))
    delta_rms = float(np.sqrt(np.mean(np.sum(delta * delta, axis=1))))
    return {
        "rms_vector_change": delta_rms,
        "relative_rms_vector_change": float(delta_rms / max(reference_rms, np.finfo(float).tiny)),
        "sampled_max_vector_change": float(np.max(np.linalg.norm(delta, axis=1))),
    }


def audit_axial_taper_divergence() -> dict[str, Any]:
    """Run the deterministic three-level public-output divergence audit."""

    points, labels = fresh_stratified_points()
    candidates: dict[float, Eq45SupportedPhi10BlendAxialTaperCandidate] = {}

    with tempfile.TemporaryDirectory(prefix="eq45-axial-taper-divergence-") as directory:
        root = Path(directory)
        for q in AXIAL_PLATEAU_VALUES:
            candidate = governed_candidate(q)
            path = root / f"candidate-{q:.6f}.json"
            candidate.save_json(path)
            loaded = Eq45SupportedPhi10BlendAxialTaperCandidate.load_json(path)
            if loaded.sha256 != candidate.sha256:
                raise RuntimeError("serialized axial-taper candidate identity changed on reload")
            candidates[float(q)] = loaded

        baseline = candidates[AXIAL_PLATEAU_VALUES[0]]
        members: list[dict[str, Any]] = []
        for q in AXIAL_PLATEAU_VALUES:
            field = candidates[q]
            time_reports: list[dict[str, Any]] = []
            for time in TIMES:
                levels: list[dict[str, Any]] = []
                for step in SPATIAL_STEPS:
                    gradient = _gradient(field, points, time, step)
                    levels.append(
                        {
                            "step": float(step),
                            **_metrics(gradient),
                            "regions": _region_metrics(gradient, labels),
                        }
                    )
                time_reports.append(
                    {
                        "time": float(time),
                        "levels": levels,
                        "observed_rms_orders": _observed_orders(levels),
                        "velocity_response_vs_q064": _velocity_response(
                            baseline, field, points, time
                        ),
                    }
                )
            members.append(
                {
                    "axial_plateau_q": float(q),
                    "axial_identity_half_height": float(field.axial_identity_half_height),
                    "candidate_sha256": field.sha256,
                    "artifact_reloaded": True,
                    "time_reports": time_reports,
                }
            )

        mutant = _PublicVelocityMutation(candidates[AXIAL_PLATEAU_VALUES[-1]])
        mutation_gradient = _gradient(mutant, points, TIMES[0], SPATIAL_STEPS[-1])
        mutation_metrics = _metrics(mutation_gradient)
        base_gradient = _gradient(
            candidates[AXIAL_PLATEAU_VALUES[-1]], points, TIMES[0], SPATIAL_STEPS[-1]
        )
        base_metrics = _metrics(base_gradient)

    return {
        "task_id": TASK_ID,
        "candidate_family": "Eq45SupportedPhi10BlendAxialTaperCandidate",
        "blend_weight_frozen": BLEND_WEIGHT,
        "axial_plateau_values_frozen": list(AXIAL_PLATEAU_VALUES),
        "times_frozen": list(TIMES),
        "spatial_steps": list(SPATIAL_STEPS),
        "probe_contract": {
            "seed": SEED,
            "point_count": int(points.shape[0]),
            "regions": {region: REGION_SIZE for region in REGIONS},
            "off_grid": True,
            "buffered_from_outer_support_faces": True,
        },
        "velocity_access": "serialized_reloaded_public_at_points_only",
        "operator": "independent_centered_second_order_cartesian_space",
        "members": members,
        "mutation_calibration": {
            "mutation": "u <- u + 0.03*x",
            "expected_added_divergence": MUTATION_DIVERGENCE,
            "baseline_finest_rms": base_metrics["sampled_rms"],
            "mutant_finest_rms": mutation_metrics["sampled_rms"],
            "mutant_finest_max_abs": mutation_metrics["sampled_max_abs"],
        },
        "axial_taper_value_selected": False,
        "registered_full_domain_divergence_gate_assessed": False,
        "full_momentum_or_vorticity_residual_assessed": False,
        "visualization_candidate_only": True,
        "visualization_ready": False,
        "visual_correspondence_verified": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
    }


def main() -> None:
    print(json.dumps(audit_axial_taper_divergence(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
