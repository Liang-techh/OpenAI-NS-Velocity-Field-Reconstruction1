"""Candidate-specific Eq. (4.5) basis-growth screen for visualization morphology.

This module deliberately consumes the frozen public-style
``Eq45VelocityCandidate.at_points(points, time)`` path. It does not fit a new
candidate or infer hidden OpenAI profile data. The purpose is narrower: test
whether four minimal symmetry-preserving spatial correction directions provide
independent control of observable velocity morphology before adding a larger
basis.

The proposed directions extend the current (radial_degree=1, eta_degree=2)
Phi/F basis to (2, 4) while keeping every existing coefficient unchanged:

* Phi (2, 0): one extra radial/poloidal direction.
* Phi (0, 4): one extra even-in-eta axial/poloidal direction.
* F   (2, 0): one extra radial/swirl direction.
* F   (0, 4): one extra even-in-eta axial/swirl direction.

All reported features are dimensionless finite-sample diagnostics derived from
[u,v,w]. They are not acceptance thresholds and are not PDE evidence.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np

from .constrained_eq45_candidate import Eq45VelocityCandidate
from .constrained_eq45_profile_basis import Eq45CompactProfileBasis


SCHEMA = "eq45_targeted_mode_screen_v1"
REFERENCE_TIMES = (0.25, 0.5, 0.75)
FEATURE_NAMES = (
    "radial_thickness_ratio",
    "axial_reach_ratio",
    "shoulder_swirl_energy_fraction",
    "shoulder_axial_poloidal_fraction",
)
MODE_SPECS = (
    ("phi", (2, 0), "phi_radial_i2_j0"),
    ("phi", (0, 4), "phi_axial_i0_j4"),
    ("swirl", (2, 0), "swirl_radial_i2_j0"),
    ("swirl", (0, 4), "swirl_axial_i0_j4"),
)
TARGET_RADIAL_DEGREE = 2
TARGET_ETA_DEGREE = 4

TRUTH_BOUNDARY = {
    "claim_scope": "local_visualization_expression_capacity_only",
    "velocity_changed": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


@dataclass(frozen=True)
class TargetedModeScreen:
    times: tuple[float, ...]
    feature_names: tuple[str, ...]
    mode_labels: tuple[str, ...]
    baseline_features: np.ndarray
    response_matrix: np.ndarray
    response_norms: np.ndarray
    singular_values: np.ndarray
    numerical_rank: int
    condition_number: float
    max_abs_column_cosine: float
    coefficient_delta: float
    refinement_delta: float
    derivative_refinement_relative_change: float

    def to_dict(self) -> dict:
        return {
            "schema": SCHEMA,
            "times": list(self.times),
            "feature_names": list(self.feature_names),
            "mode_labels": list(self.mode_labels),
            "baseline_features": self.baseline_features.tolist(),
            "response_matrix": self.response_matrix.tolist(),
            "response_norms": self.response_norms.tolist(),
            "singular_values": self.singular_values.tolist(),
            "numerical_rank": int(self.numerical_rank),
            "condition_number": float(self.condition_number),
            "max_abs_column_cosine": float(self.max_abs_column_cosine),
            "coefficient_delta": float(self.coefficient_delta),
            "refinement_delta": float(self.refinement_delta),
            "derivative_refinement_relative_change": float(
                self.derivative_refinement_relative_change
            ),
            "truth_boundary": dict(TRUTH_BOUNDARY),
        }


def extend_eq45_basis(
    basis: Eq45CompactProfileBasis,
    *,
    radial_degree: int = TARGET_RADIAL_DEGREE,
    eta_degree: int = TARGET_ETA_DEGREE,
    perturbation: tuple[str, tuple[int, int], float] | None = None,
) -> Eq45CompactProfileBasis:
    """Embed a smaller Phi/F basis exactly into a larger tensor-product basis."""
    if not isinstance(basis, Eq45CompactProfileBasis):
        raise TypeError("basis must be Eq45CompactProfileBasis")
    if radial_degree < basis.radial_degree or eta_degree < basis.eta_degree:
        raise ValueError("target degrees may not discard existing basis modes")

    target_modes = tuple(
        (i, j)
        for i in range(radial_degree + 1)
        for j in range(eta_degree + 1)
    )
    target_index = {mode: index for index, mode in enumerate(target_modes)}
    phi = np.zeros(len(target_modes), dtype=float)
    swirl = np.zeros(len(target_modes), dtype=float)

    for source_index, mode in enumerate(basis.mode_indices):
        destination_index = target_index[mode]
        phi[destination_index] = basis.phi_coefficients[source_index]
        swirl[destination_index] = basis.swirl_coefficients[source_index]

    if perturbation is not None:
        channel, mode, delta = perturbation
        if channel not in {"phi", "swirl"}:
            raise ValueError("perturbation channel must be 'phi' or 'swirl'")
        if mode not in target_index:
            raise ValueError("perturbation mode is outside the target basis")
        if mode in basis.mode_indices:
            raise ValueError("targeted growth perturbation must be a newly added mode")
        if not np.isfinite(delta):
            raise ValueError("perturbation delta must be finite")
        values = phi if channel == "phi" else swirl
        values[target_index[mode]] += float(delta)

    return Eq45CompactProfileBasis(
        radial_degree=radial_degree,
        eta_degree=eta_degree,
        phi_coefficients=tuple(float(value) for value in phi),
        swirl_coefficients=tuple(float(value) for value in swirl),
        x_cut=basis.x_cut,
        eta_cut=basis.eta_cut,
        cutoff_power=basis.cutoff_power,
        coefficient_limit=basis.coefficient_limit,
    )


def _candidate_with_perturbation(
    candidate: Eq45VelocityCandidate,
    perturbation: tuple[str, tuple[int, int], float] | None,
) -> Eq45VelocityCandidate:
    return Eq45VelocityCandidate(
        profile_basis=extend_eq45_basis(
            candidate.profile_basis,
            perturbation=perturbation,
        ),
        h=candidate.h,
        time_start=candidate.time_start,
        time_end=candidate.time_end,
    )


def _ring_points(radius: float, z: float, n_angles: int) -> np.ndarray:
    theta = np.linspace(0.0, 2.0 * np.pi, n_angles, endpoint=False)
    return np.column_stack(
        (
            radius * np.cos(theta),
            radius * np.sin(theta),
            np.full(n_angles, z, dtype=float),
        )
    )


def _velocity(field: Eq45VelocityCandidate, points: np.ndarray, time: float) -> np.ndarray:
    values = np.asarray(field.at_points(points, time), dtype=float)
    if values.shape != points.shape:
        raise ValueError(f"velocity must return shape {points.shape}, got {values.shape}")
    if not np.all(np.isfinite(values)):
        raise ValueError("velocity returned nonfinite values")
    return values


def _cylindrical_components(points: np.ndarray, velocity: np.ndarray):
    radius = np.hypot(points[:, 0], points[:, 1])
    if np.any(radius <= 0.0):
        raise ValueError("ring probes must stay away from the symmetry axis")
    er_x = points[:, 0] / radius
    er_y = points[:, 1] / radius
    radial = velocity[:, 0] * er_x + velocity[:, 1] * er_y
    azimuthal = -velocity[:, 0] * er_y + velocity[:, 1] * er_x
    return radial, azimuthal, velocity[:, 2]


def _rms(values: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.asarray(values, dtype=float) ** 2)))


def dimensionless_morphology_features(
    field: Eq45VelocityCandidate,
    time: float,
    *,
    n_angles: int = 24,
    inner_radius: float = 0.45,
    mid_radius: float = 0.65,
    outer_radius: float = 1.05,
    shoulder_z: float = 0.65,
) -> np.ndarray:
    """Return four dimensionless public-velocity morphology features.

    The first two ratios probe radial thickness and axial reach. The latter two
    split shoulder-ring kinetic energy into swirl and axial-poloidal content.
    Symmetric +/-z shoulder rings are used so the diagnostic does not assume a
    sign convention for the north/south poloidal circulation.
    """
    if not isinstance(field, Eq45VelocityCandidate):
        raise TypeError("field must be Eq45VelocityCandidate")
    if not isinstance(n_angles, (int, np.integer)) or n_angles < 8:
        raise ValueError("n_angles must be an integer >= 8")
    geometry = np.asarray(
        (inner_radius, mid_radius, outer_radius, shoulder_z), dtype=float
    )
    if not np.all(np.isfinite(geometry)) or np.any(geometry <= 0.0):
        raise ValueError("probe radii and shoulder_z must be finite and positive")
    if not (inner_radius < mid_radius < outer_radius):
        raise ValueError("probe radii must satisfy inner < mid < outer")
    if not np.isfinite(time) or not (field.time_start <= time <= field.time_end):
        raise ValueError("time must lie inside the candidate interval")

    inner = _ring_points(inner_radius, 0.0, n_angles)
    mid = _ring_points(mid_radius, 0.0, n_angles)
    outer = _ring_points(outer_radius, 0.0, n_angles)
    shoulder_plus = _ring_points(mid_radius, shoulder_z, n_angles)
    shoulder_minus = _ring_points(mid_radius, -shoulder_z, n_angles)

    velocity_inner = _velocity(field, inner, time)
    velocity_mid = _velocity(field, mid, time)
    velocity_outer = _velocity(field, outer, time)
    velocity_shoulder = np.concatenate(
        (
            _velocity(field, shoulder_plus, time),
            _velocity(field, shoulder_minus, time),
        ),
        axis=0,
    )
    shoulder_points = np.concatenate((shoulder_plus, shoulder_minus), axis=0)

    speed_inner = np.linalg.norm(velocity_inner, axis=1)
    speed_mid = np.linalg.norm(velocity_mid, axis=1)
    speed_outer = np.linalg.norm(velocity_outer, axis=1)
    speed_shoulder = np.linalg.norm(velocity_shoulder, axis=1)

    inner_rms = _rms(speed_inner)
    mid_rms = _rms(speed_mid)
    if inner_rms <= 1e-14 or mid_rms <= 1e-14:
        raise ValueError("candidate is numerically inactive on morphology reference rings")

    radial, azimuthal, axial = _cylindrical_components(
        shoulder_points, velocity_shoulder
    )
    total_energy = float(np.mean(np.sum(velocity_shoulder**2, axis=1)))
    poloidal_energy = float(np.mean(radial**2 + axial**2))
    if total_energy <= 1e-28 or poloidal_energy <= 1e-28:
        raise ValueError("candidate is numerically inactive on shoulder probes")

    return np.asarray(
        (
            _rms(speed_outer) / inner_rms,
            _rms(speed_shoulder) / mid_rms,
            float(np.mean(azimuthal**2)) / total_energy,
            float(np.mean(axial**2)) / poloidal_energy,
        ),
        dtype=float,
    )


def _feature_vector(
    candidate: Eq45VelocityCandidate,
    times: Iterable[float],
    *,
    n_angles: int,
) -> np.ndarray:
    return np.concatenate(
        [
            dimensionless_morphology_features(candidate, float(time), n_angles=n_angles)
            for time in times
        ]
    )


def _response_matrix(
    candidate: Eq45VelocityCandidate,
    times: tuple[float, ...],
    *,
    coefficient_delta: float,
    n_angles: int,
) -> np.ndarray:
    columns = []
    for channel, mode, _label in MODE_SPECS:
        plus = _candidate_with_perturbation(
            candidate, (channel, mode, coefficient_delta)
        )
        minus = _candidate_with_perturbation(
            candidate, (channel, mode, -coefficient_delta)
        )
        response = (
            _feature_vector(plus, times, n_angles=n_angles)
            - _feature_vector(minus, times, n_angles=n_angles)
        ) / (2.0 * coefficient_delta)
        columns.append(response)
    return np.column_stack(columns)


def screen_eq45_targeted_modes(
    candidate: Eq45VelocityCandidate,
    *,
    times: tuple[float, ...] = REFERENCE_TIMES,
    coefficient_delta: float = 1e-3,
    refinement_delta: float = 5e-4,
    n_angles: int = 24,
    rank_rtol: float = 1e-8,
) -> TargetedModeScreen:
    """Screen four minimal symmetry-preserving growth directions on [u,v,w]."""
    if not isinstance(candidate, Eq45VelocityCandidate):
        raise TypeError("candidate must be Eq45VelocityCandidate")
    times = tuple(float(time) for time in times)
    if len(times) < 3 or not np.all(np.isfinite(times)):
        raise ValueError("times must contain at least three finite values")
    if any(left >= right for left, right in zip(times, times[1:])):
        raise ValueError("times must be strictly increasing")
    if times[0] < candidate.time_start or times[-1] > candidate.time_end:
        raise ValueError("times must lie within the candidate interval")
    if not np.isfinite(coefficient_delta) or coefficient_delta <= 0.0:
        raise ValueError("coefficient_delta must be finite and positive")
    if not np.isfinite(refinement_delta) or not (
        0.0 < refinement_delta < coefficient_delta
    ):
        raise ValueError(
            "refinement_delta must satisfy 0 < refinement_delta < coefficient_delta"
        )
    if not np.isfinite(rank_rtol) or not (0.0 < rank_rtol < 1.0):
        raise ValueError("rank_rtol must lie in (0,1)")

    baseline = np.vstack(
        [
            dimensionless_morphology_features(candidate, time, n_angles=n_angles)
            for time in times
        ]
    )
    response = _response_matrix(
        candidate,
        times,
        coefficient_delta=coefficient_delta,
        n_angles=n_angles,
    )
    refined = _response_matrix(
        candidate,
        times,
        coefficient_delta=refinement_delta,
        n_angles=n_angles,
    )

    singular_values = np.linalg.svd(response, compute_uv=False)
    if singular_values.size == 0 or singular_values[0] <= 0.0:
        raise ValueError("targeted mode response is numerically inactive")
    threshold = rank_rtol * singular_values[0]
    numerical_rank = int(np.count_nonzero(singular_values > threshold))
    if numerical_rank == 0:
        raise ValueError("targeted mode response has zero numerical rank")
    if numerical_rank < response.shape[1]:
        condition_number = float("inf")
    else:
        condition_number = float(singular_values[0] / singular_values[-1])

    response_norms = np.linalg.norm(response, axis=0)
    if np.any(response_norms <= 0.0):
        raise ValueError("one or more targeted directions are locally dead")
    cosine = (response.T @ response) / (
        response_norms[:, None] * response_norms[None, :]
    )
    off_diagonal = np.abs(cosine - np.eye(cosine.shape[0]))
    max_abs_column_cosine = float(np.max(off_diagonal))

    refined_norm = float(np.linalg.norm(refined))
    if refined_norm <= 0.0:
        raise ValueError("refined targeted response is numerically inactive")
    derivative_refinement_relative_change = float(
        np.linalg.norm(response - refined) / refined_norm
    )

    arrays = (baseline, response, response_norms, singular_values)
    if not all(np.all(np.isfinite(array)) for array in arrays):
        raise ValueError("targeted mode screen produced nonfinite diagnostics")

    baseline.setflags(write=False)
    response.setflags(write=False)
    response_norms.setflags(write=False)
    singular_values.setflags(write=False)
    return TargetedModeScreen(
        times=times,
        feature_names=FEATURE_NAMES,
        mode_labels=tuple(spec[2] for spec in MODE_SPECS),
        baseline_features=baseline,
        response_matrix=response,
        response_norms=response_norms,
        singular_values=singular_values,
        numerical_rank=numerical_rank,
        condition_number=condition_number,
        max_abs_column_cosine=max_abs_column_cosine,
        coefficient_delta=float(coefficient_delta),
        refinement_delta=float(refinement_delta),
        derivative_refinement_relative_change=derivative_refinement_relative_change,
    )
