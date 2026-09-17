"""Local axial/end-shape capacity audit for the materialized Phi10 blend family.

The compact--quartic ``Phi(1,0)`` blend already supplies one bounded temporal
control that primarily changes early radial/collar morphology.  This Agent-7
increment asks whether the already-existing ``Phi(1,2)`` spatial coefficient
adds a genuinely independent axial/end-shape direction before anybody grows
the spatial basis further.

Only centered local responses of the public support-connected ``[u,v,w]``
field are measured.  No blend weight or profile coefficient is selected, no
public image/PDE residual/force/pressure is fitted, and no saved candidate is
modified or promoted.
"""
from __future__ import annotations

from dataclasses import replace
from typing import Iterable

import numpy as np

from .constrained_eq45_supported_candidate import Eq45SupportedVelocityCandidate
from .constrained_eq45_supported_collar_vorticity import governed_supported_seed
from .constrained_eq45_supported_phi10_compact_quartic_blend_temporal_mode import (
    Eq45SupportedPhi10CompactQuarticBlendTemporalCandidate,
)


TASK_ID = "CR003-EQ45-SUPPORTED-BLEND-PHI12-AXIAL-CAPACITY-033"
EXPECTED_BASE_SHA256 = "2fdff812c22131d56eac1b7d6e455187ff3207500c6385aa9abf282a8e3d7b1d"
PHI12_MODE = (1, 2)
PARAMETER_LABELS = ("compact_quartic_blend_weight", "Phi(1,2)")
DEFAULT_BLEND_WEIGHT = 0.5
DEFAULT_TIME = 0.3125
DEFAULT_BLEND_STEPS = (0.10, 0.05)
DEFAULT_PHI12_STEPS = (0.04, 0.02)


def _with_phi12_delta(
    base: Eq45SupportedVelocityCandidate,
    *,
    delta: float,
) -> Eq45SupportedVelocityCandidate:
    if not isinstance(base, Eq45SupportedVelocityCandidate):
        raise TypeError("base must be Eq45SupportedVelocityCandidate")
    if not np.isscalar(delta) or not np.isfinite(delta):
        raise ValueError("Phi(1,2) delta must be finite")

    parent = base.parent
    basis = parent.profile_basis
    try:
        index = basis.mode_indices.index(PHI12_MODE)
    except ValueError as exc:
        raise ValueError("base profile basis must contain Phi(1,2)") from exc

    coefficients = list(basis.phi_coefficients)
    shifted = float(coefficients[index]) + float(delta)
    if abs(shifted) > float(basis.coefficient_limit):
        raise ValueError("Phi(1,2) perturbation exceeds the registered coefficient bound")
    coefficients[index] = shifted
    shifted_basis = replace(basis, phi_coefficients=tuple(coefficients))
    shifted_parent = replace(parent, profile_basis=shifted_basis)
    return replace(base, parent=shifted_parent)


def _cylindrical_points(
    radii: Iterable[float],
    z_values: Iterable[float],
    *,
    n_angles: int,
) -> np.ndarray:
    if not isinstance(n_angles, int) or n_angles < 4:
        raise ValueError("n_angles must be an integer >= 4")
    angles = np.linspace(0.0, 2.0 * np.pi, n_angles, endpoint=False)
    points: list[tuple[float, float, float]] = []
    for radius in radii:
        for z_value in z_values:
            for angle in angles:
                points.append(
                    (
                        float(radius) * float(np.cos(angle)),
                        float(radius) * float(np.sin(angle)),
                        float(z_value),
                    )
                )
    return np.asarray(points, dtype=float)


def _probe_sets(n_angles: int) -> dict[str, np.ndarray]:
    # Radial probes emphasize the shoulder / radial support collar while keeping
    # |z| moderate.  Axial probes keep radius moderate and approach the upper
    # and lower support caps.  Interfaces and the exact r=2, |z|=2 faces are
    # deliberately avoided.
    radial = _cylindrical_points(
        (1.55, 1.72, 1.84),
        (-0.65, 0.0, 0.65),
        n_angles=n_angles,
    )
    axial = _cylindrical_points(
        (0.35, 0.80, 1.25),
        (-1.82, -1.68, 1.68, 1.82),
        n_angles=n_angles,
    )
    combined = np.concatenate((radial, axial), axis=0)
    return {"radial": radial, "axial": axial, "combined": combined}


def _blend_response(
    base: Eq45SupportedVelocityCandidate,
    *,
    blend_weight: float,
    step: float,
    points: np.ndarray,
    time: float,
) -> np.ndarray:
    plus_weight = float(blend_weight) + float(step)
    minus_weight = float(blend_weight) - float(step)
    if not (0.0 <= minus_weight < plus_weight <= 1.0):
        raise ValueError("blend centered difference must remain inside [0,1]")
    plus = Eq45SupportedPhi10CompactQuarticBlendTemporalCandidate(
        base=base,
        blend_weight=plus_weight,
    )
    minus = Eq45SupportedPhi10CompactQuarticBlendTemporalCandidate(
        base=base,
        blend_weight=minus_weight,
    )
    derivative = (
        plus.at_points(points, time) - minus.at_points(points, time)
    ) / (2.0 * float(step))
    out = np.asarray(derivative, dtype=float).reshape(-1)
    if not np.all(np.isfinite(out)):
        raise RuntimeError("blend-weight response became nonfinite")
    return out


def _phi12_response(
    base: Eq45SupportedVelocityCandidate,
    *,
    blend_weight: float,
    step: float,
    points: np.ndarray,
    time: float,
) -> np.ndarray:
    plus_base = _with_phi12_delta(base, delta=step)
    minus_base = _with_phi12_delta(base, delta=-step)
    plus = Eq45SupportedPhi10CompactQuarticBlendTemporalCandidate(
        base=plus_base,
        blend_weight=blend_weight,
    )
    minus = Eq45SupportedPhi10CompactQuarticBlendTemporalCandidate(
        base=minus_base,
        blend_weight=blend_weight,
    )
    derivative = (
        plus.at_points(points, time) - minus.at_points(points, time)
    ) / (2.0 * float(step))
    out = np.asarray(derivative, dtype=float).reshape(-1)
    if not np.all(np.isfinite(out)):
        raise RuntimeError("Phi(1,2) response became nonfinite")
    return out


def _response_matrix(
    base: Eq45SupportedVelocityCandidate,
    *,
    blend_weight: float,
    blend_step: float,
    phi12_step: float,
    points: np.ndarray,
    time: float,
) -> np.ndarray:
    return np.column_stack(
        (
            _blend_response(
                base,
                blend_weight=blend_weight,
                step=blend_step,
                points=points,
                time=time,
            ),
            _phi12_response(
                base,
                blend_weight=blend_weight,
                step=phi12_step,
                points=points,
                time=time,
            ),
        )
    )


def _relative_change(coarse: np.ndarray, fine: np.ndarray) -> float:
    denominator = max(float(np.linalg.norm(fine)), np.finfo(float).tiny)
    return float(np.linalg.norm(coarse - fine) / denominator)


def _rms_columns(matrix: np.ndarray) -> np.ndarray:
    return np.sqrt(np.mean(np.square(matrix), axis=0))


def audit_blend_phi12_axial_capacity(
    base: Eq45SupportedVelocityCandidate | None = None,
    *,
    blend_weight: float = DEFAULT_BLEND_WEIGHT,
    time: float = DEFAULT_TIME,
    blend_steps: tuple[float, ...] = DEFAULT_BLEND_STEPS,
    phi12_steps: tuple[float, ...] = DEFAULT_PHI12_STEPS,
    n_angles: int = 8,
) -> dict[str, object]:
    """Measure whether existing Phi(1,2) complements the Phi10 blend direction."""
    if base is None:
        base = governed_supported_seed()
    if not isinstance(base, Eq45SupportedVelocityCandidate):
        raise TypeError("base must be Eq45SupportedVelocityCandidate")
    if base.sha256 != EXPECTED_BASE_SHA256:
        raise ValueError("audit requires the frozen governed supported Eq45 base")

    weight = float(blend_weight)
    checked_time = float(time)
    if not np.isfinite(weight) or not (0.0 < weight < 1.0):
        raise ValueError("blend_weight must lie strictly inside (0,1)")
    if not np.isfinite(checked_time) or not (base.time_start <= checked_time <= base.time_end):
        raise ValueError("time must lie inside the supported candidate delivery interval")

    blend_steps = tuple(float(value) for value in blend_steps)
    phi12_steps = tuple(float(value) for value in phi12_steps)
    if len(blend_steps) < 2 or len(phi12_steps) < 2 or len(blend_steps) != len(phi12_steps):
        raise ValueError("blend_steps and phi12_steps must have equal length >= 2")
    for name, values in (("blend_steps", blend_steps), ("phi12_steps", phi12_steps)):
        if not all(np.isfinite(value) and value > 0.0 for value in values):
            raise ValueError(f"{name} must contain positive finite values")
        if not all(values[index] > values[index + 1] for index in range(len(values) - 1)):
            raise ValueError(f"{name} must be strictly decreasing")
    if weight - max(blend_steps) < 0.0 or weight + max(blend_steps) > 1.0:
        raise ValueError("largest blend step leaves [0,1]")

    probes = _probe_sets(n_angles)
    combined_by_level = [
        _response_matrix(
            base,
            blend_weight=weight,
            blend_step=blend_step,
            phi12_step=phi12_step,
            points=probes["combined"],
            time=checked_time,
        )
        for blend_step, phi12_step in zip(blend_steps, phi12_steps)
    ]
    fine = combined_by_level[-1]
    radial = _response_matrix(
        base,
        blend_weight=weight,
        blend_step=blend_steps[-1],
        phi12_step=phi12_steps[-1],
        points=probes["radial"],
        time=checked_time,
    )
    axial = _response_matrix(
        base,
        blend_weight=weight,
        blend_step=blend_steps[-1],
        phi12_step=phi12_steps[-1],
        points=probes["axial"],
        time=checked_time,
    )

    singular_values = np.linalg.svd(fine, compute_uv=False)
    largest = float(singular_values[0])
    smallest = float(singular_values[-1])
    tolerance = max(fine.shape) * np.finfo(float).eps * max(largest, np.finfo(float).tiny)
    numerical_rank = int(np.count_nonzero(singular_values > tolerance))
    condition_number = float(largest / smallest) if smallest > 0.0 else float("inf")

    blend_direction = fine[:, 0]
    phi12_direction = fine[:, 1]
    blend_norm_sq = float(np.dot(blend_direction, blend_direction))
    projection = blend_direction * (
        float(np.dot(blend_direction, phi12_direction))
        / max(blend_norm_sq, np.finfo(float).tiny)
    )
    phi12_novelty = float(
        np.linalg.norm(phi12_direction - projection)
        / max(float(np.linalg.norm(phi12_direction)), np.finfo(float).tiny)
    )
    column_cosine = float(
        np.dot(blend_direction, phi12_direction)
        / max(
            float(np.linalg.norm(blend_direction) * np.linalg.norm(phi12_direction)),
            np.finfo(float).tiny,
        )
    )

    radial_rms = _rms_columns(radial)
    axial_rms = _rms_columns(axial)
    combined_rms = _rms_columns(fine)
    axial_to_radial = axial_rms / np.maximum(radial_rms, np.finfo(float).tiny)

    return {
        "task_id": TASK_ID,
        "base_sha256": base.sha256,
        "parameter_labels": list(PARAMETER_LABELS),
        "blend_weight": weight,
        "time": checked_time,
        "blend_steps": list(blend_steps),
        "phi12_steps": list(phi12_steps),
        "n_angles": int(n_angles),
        "probe_counts": {name: int(points.shape[0]) for name, points in probes.items()},
        "singular_values": [float(value) for value in singular_values],
        "numerical_rank": numerical_rank,
        "condition_number": condition_number,
        "phi12_novelty_outside_blend_span": phi12_novelty,
        "response_column_cosine": column_cosine,
        "step_refinement_relative_changes": [
            _relative_change(combined_by_level[index], combined_by_level[index + 1])
            for index in range(len(combined_by_level) - 1)
        ],
        "combined_response_rms": [float(value) for value in combined_rms],
        "radial_response_rms": [float(value) for value in radial_rms],
        "axial_response_rms": [float(value) for value in axial_rms],
        "axial_to_radial_response_rms_ratio": [float(value) for value in axial_to_radial],
        "truth_boundary": {
            "velocity_changed": False,
            "canonical_velocity_changed": False,
            "production_coefficients_changed": False,
            "new_basis_added": False,
            "blend_weight_selected": False,
            "phi12_coefficient_selected": False,
            "force_or_pressure_fitted": False,
            "pde_objective_used": False,
            "public_image_fitted": False,
            "visualization_ready": False,
            "visual_correspondence_verified": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
    }


if __name__ == "__main__":
    import json

    print(json.dumps(audit_blend_phi12_axial_capacity(), indent=2, sort_keys=True))
