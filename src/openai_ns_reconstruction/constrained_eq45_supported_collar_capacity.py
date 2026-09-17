"""Local basis-capacity audit after the Eq45 physical support transform.

This Agent-7 diagnostic asks a narrow representation question: after the
streamfunction-level physical-space support connection is applied, do two
already-existing even Phi channels still provide distinct controls over the
callable ``[u,v,w]`` field?

The compared directions are ``Phi(1,0)`` (primarily radial) and ``Phi(1,2)``
(mixed radial x axial/end-shape).  No coefficient is fitted or promoted.  The
audit uses centered coefficient perturbations only to form a local velocity
Jacobian on deterministic plateau and support-collar probes.
"""
from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from typing import Iterable

import numpy as np

from .constrained_eq45_candidate import Eq45VelocityCandidate
from .constrained_eq45_supported_candidate import Eq45SupportedVelocityCandidate


TASK_ID = "CR003-EQ45-SUPPORTED-COLLAR-CAPACITY-017"
PARAMETER_LABELS = ("Phi(1,0)", "Phi(1,2)")
MODES = ((1, 0), (1, 2))
EXPECTED_PARENT_SHA256 = "48f1845fcfd71ec95495c98bb7aac3fca4653e748a8856a5421234e9601525a7"
DEFAULT_COEFFICIENT_STEPS = (0.04, 0.02)
DEFAULT_TIMES = (0.25, 0.50, 0.75)

ROOT = Path(__file__).resolve().parents[2]
SEED_PATH = ROOT / "artifacts" / "constrained" / "eq45_velocity_candidate_seed.json"


def _with_phi_delta(
    candidate: Eq45VelocityCandidate,
    *,
    mode: tuple[int, int],
    delta: float,
) -> Eq45VelocityCandidate:
    """Return a bounded child parent-candidate with one Phi coefficient shifted."""
    if mode not in MODES:
        raise ValueError(f"mode must be one of {MODES}")
    if not np.isscalar(delta) or not np.isfinite(delta):
        raise ValueError("delta must be finite")

    basis = candidate.profile_basis
    try:
        index = basis.mode_indices.index(mode)
    except ValueError as exc:
        raise ValueError(f"candidate basis does not contain mode {mode}") from exc

    values = list(basis.phi_coefficients)
    shifted = float(values[index]) + float(delta)
    if abs(shifted) > float(basis.coefficient_limit):
        raise ValueError("perturbation exceeds the registered coefficient bound")
    values[index] = shifted
    return replace(candidate, profile_basis=replace(basis, phi_coefficients=tuple(values)))


def _cylindrical_points(
    radii: Iterable[float],
    z_values: Iterable[float],
    *,
    n_angles: int,
) -> np.ndarray:
    if not isinstance(n_angles, int) or n_angles < 4:
        raise ValueError("n_angles must be an integer >= 4")
    angles = np.linspace(0.0, 2.0 * np.pi, n_angles, endpoint=False)
    points = []
    for radius in radii:
        for z in z_values:
            for angle in angles:
                points.append(
                    (
                        float(radius) * float(np.cos(angle)),
                        float(radius) * float(np.sin(angle)),
                        float(z),
                    )
                )
    return np.asarray(points, dtype=float)


def _probe_sets(n_angles: int) -> dict[str, np.ndarray]:
    # All plateau probes are comfortably inside the default r<=1.6, |z|<=1.6
    # identity region.  Collar probes sit strictly between plateau and support
    # faces, avoiding the nonsmooth bookkeeping issue of sampling an interface.
    plateau = _cylindrical_points((0.35, 0.85, 1.35), (-1.10, 0.0, 1.10), n_angles=n_angles)
    radial_collar = _cylindrical_points(
        (1.66, 1.74, 1.84), (-0.85, 0.0, 0.85), n_angles=n_angles
    )
    axial_collar = _cylindrical_points(
        (0.35, 0.85, 1.35), (-1.84, -1.72, 1.72, 1.84), n_angles=n_angles
    )
    collar = np.concatenate((radial_collar, axial_collar), axis=0)
    return {
        "plateau": plateau,
        "radial_collar": radial_collar,
        "axial_collar": axial_collar,
        "collar": collar,
    }


def _response_vector(
    parent: Eq45VelocityCandidate,
    *,
    mode: tuple[int, int],
    step: float,
    points: np.ndarray,
    times: tuple[float, ...],
    supported: bool,
) -> np.ndarray:
    plus_parent = _with_phi_delta(parent, mode=mode, delta=step)
    minus_parent = _with_phi_delta(parent, mode=mode, delta=-step)
    if supported:
        plus = Eq45SupportedVelocityCandidate(parent=plus_parent)
        minus = Eq45SupportedVelocityCandidate(parent=minus_parent)
    else:
        plus = plus_parent
        minus = minus_parent

    blocks = []
    for time in times:
        plus_values = plus.at_points(points, time)
        minus_values = minus.at_points(points, time)
        derivative = (plus_values - minus_values) / (2.0 * step)
        blocks.append(np.asarray(derivative, dtype=float).reshape(-1))
    out = np.concatenate(blocks)
    if not np.all(np.isfinite(out)):
        raise RuntimeError("coefficient response became nonfinite")
    return out


def _response_matrix(
    parent: Eq45VelocityCandidate,
    *,
    step: float,
    points: np.ndarray,
    times: tuple[float, ...],
    supported: bool,
) -> np.ndarray:
    return np.column_stack(
        [
            _response_vector(
                parent,
                mode=mode,
                step=step,
                points=points,
                times=times,
                supported=supported,
            )
            for mode in MODES
        ]
    )


def _relative_change(a: np.ndarray, b: np.ndarray) -> float:
    denominator = max(float(np.linalg.norm(b)), np.finfo(float).tiny)
    return float(np.linalg.norm(a - b) / denominator)


def audit_eq45_supported_collar_capacity(
    parent: Eq45VelocityCandidate | None = None,
    *,
    coefficient_steps: tuple[float, ...] = DEFAULT_COEFFICIENT_STEPS,
    times: tuple[float, ...] = DEFAULT_TIMES,
    n_angles: int = 8,
) -> dict[str, object]:
    """Audit two existing Phi controls after the physical support connection.

    The returned Jacobian is a local representation diagnostic only.  It does
    not fit a target, modify the saved candidate, or establish visual/PDE
    acceptance.
    """
    if parent is None:
        parent = Eq45VelocityCandidate.load_json(SEED_PATH)
    if not isinstance(parent, Eq45VelocityCandidate):
        raise TypeError("parent must be Eq45VelocityCandidate")
    if parent.sha256 != EXPECTED_PARENT_SHA256:
        raise ValueError("audit requires the frozen governed Eq45 parent identity")

    if len(coefficient_steps) < 2:
        raise ValueError("coefficient_steps must contain at least two levels")
    steps = tuple(float(value) for value in coefficient_steps)
    if not all(np.isfinite(value) and value > 0.0 for value in steps):
        raise ValueError("coefficient_steps must be positive and finite")
    if not all(steps[index] > steps[index + 1] for index in range(len(steps) - 1)):
        raise ValueError("coefficient_steps must be strictly decreasing")

    checked_times = tuple(float(value) for value in times)
    if not checked_times or not all(np.isfinite(value) for value in checked_times):
        raise ValueError("times must be finite and nonempty")
    if min(checked_times) < parent.time_start or max(checked_times) > parent.time_end:
        raise ValueError("times must stay inside the parent delivery interval")

    probes = _probe_sets(n_angles)
    supported_collar_by_step = []
    for step in steps:
        supported_collar_by_step.append(
            _response_matrix(
                parent,
                step=step,
                points=probes["collar"],
                times=checked_times,
                supported=True,
            )
        )

    finest_step = steps[-1]
    supported_collar = supported_collar_by_step[-1]
    parent_collar = _response_matrix(
        parent,
        step=finest_step,
        points=probes["collar"],
        times=checked_times,
        supported=False,
    )
    supported_plateau = _response_matrix(
        parent,
        step=finest_step,
        points=probes["plateau"],
        times=checked_times,
        supported=True,
    )
    parent_plateau = _response_matrix(
        parent,
        step=finest_step,
        points=probes["plateau"],
        times=checked_times,
        supported=False,
    )
    supported_radial = _response_matrix(
        parent,
        step=finest_step,
        points=probes["radial_collar"],
        times=checked_times,
        supported=True,
    )
    supported_axial = _response_matrix(
        parent,
        step=finest_step,
        points=probes["axial_collar"],
        times=checked_times,
        supported=True,
    )

    singular_values = np.linalg.svd(supported_collar, compute_uv=False)
    largest = float(singular_values[0])
    smallest = float(singular_values[-1])
    tolerance = max(supported_collar.shape) * np.finfo(float).eps * largest
    numerical_rank = int(np.count_nonzero(singular_values > tolerance))
    condition_number = float(largest / smallest) if smallest > 0.0 else float("inf")

    radial_direction = supported_collar[:, 0]
    mixed_direction = supported_collar[:, 1]
    projection = radial_direction * (
        float(np.dot(radial_direction, mixed_direction))
        / max(float(np.dot(radial_direction, radial_direction)), np.finfo(float).tiny)
    )
    mixed_novelty = float(
        np.linalg.norm(mixed_direction - projection)
        / max(float(np.linalg.norm(mixed_direction)), np.finfo(float).tiny)
    )

    refinement_changes = [
        _relative_change(supported_collar_by_step[index], supported_collar_by_step[index + 1])
        for index in range(len(supported_collar_by_step) - 1)
    ]

    response_norms = np.linalg.norm(supported_collar, axis=0)
    radial_region_norms = np.linalg.norm(supported_radial, axis=0)
    axial_region_norms = np.linalg.norm(supported_axial, axis=0)

    return {
        "task_id": TASK_ID,
        "parent_sha256": parent.sha256,
        "supported_child_sha256": Eq45SupportedVelocityCandidate(parent=parent).sha256,
        "parameter_labels": list(PARAMETER_LABELS),
        "coefficient_steps": list(steps),
        "times": list(checked_times),
        "n_angles": int(n_angles),
        "probe_counts": {name: int(value.shape[0]) for name, value in probes.items()},
        "finest_collar_response_norms": [float(value) for value in response_norms],
        "radial_collar_response_norms": [float(value) for value in radial_region_norms],
        "axial_collar_response_norms": [float(value) for value in axial_region_norms],
        "axial_response_share": [
            float(axial_region_norms[index] / max(response_norms[index], np.finfo(float).tiny))
            for index in range(len(PARAMETER_LABELS))
        ],
        "singular_values": [float(value) for value in singular_values],
        "numerical_rank": numerical_rank,
        "condition_number": condition_number,
        "phi_1_2_novelty_outside_phi_1_0_span": mixed_novelty,
        "step_refinement_relative_changes": refinement_changes,
        "plateau_supported_vs_parent_response_relative_change": [
            _relative_change(supported_plateau[:, index], parent_plateau[:, index])
            for index in range(len(PARAMETER_LABELS))
        ],
        "collar_supported_vs_parent_response_relative_change": [
            _relative_change(supported_collar[:, index], parent_collar[:, index])
            for index in range(len(PARAMETER_LABELS))
        ],
        "truth_boundary": {
            "velocity_changed": False,
            "production_coefficients_changed": False,
            "new_basis_added": False,
            "forcing_or_pressure_refit": False,
            "physical_support_transform_changed": False,
            "visualization_ready": False,
            "visual_correspondence_verified": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
    }
