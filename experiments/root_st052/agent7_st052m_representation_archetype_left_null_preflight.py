"""Screen minimal representation archetypes against the existing ST052-M defect-space left null.

Preregistered in issue #1146 and stacked exactly on Agent-7 PR #1139 head.
The probes are analysis-only divergence-free tangents.  They do not change the
selected/canonical/saved velocity, add a basis to the candidate, select a
coefficient, or establish public-source visual correspondence.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Callable

import numpy as np

import agent7_st052m_defect_space_projection as parent

TASK_ID = "CR003-ST052M-REPRESENTATION-ARCHETYPE-LEFT-NULL-PREFLIGHT-136"
PREREG_ISSUE = 1146
STACK_BASE_PR = 1139
STACK_BASE_HEAD = "dc801576b7037667642db49665f97b3dead29033"

CHANNELS = tuple(parent.CHANNELS)
DEFECT_COORDINATES = tuple(parent.DEFECT_COORDINATES)
PROBE_NAMES = (
    "poloidal_scalar_space_extension",
    "compact_toroidal_envelope_extension",
    "time_skewed_poloidal_extension",
)

SUPPORT_RADIUS = 1.55
SUPPORT_ABS_Z = 1.75
RANK_TOL = 1.0e-10
ROW_RESPONSE_FLOOR = 1.0e-12
DERIVATIVE_DRIFT_MAX = 0.05
DERIVATIVE_COSINE_MIN = 0.999
CONDITION_MAX = 25.0
COSINE_MAX_ABS = 0.995
DIVERGENCE_FD_STEP = 1.0e-5
DIVERGENCE_FD_MAX = 5.0e-7
OUTSIDE_SUPPORT_MAX = 1.0e-12
AXIS_TRANSVERSE_MAX = 1.0e-12

TRUTH = {
    "candidate_velocity_changed": False,
    "canonical_velocity_changed": False,
    "saved_velocity_changed": False,
    "basis_dimension_changed": False,
    "new_spatial_basis_added_to_candidate": False,
    "new_temporal_basis_added_to_candidate": False,
    "analysis_only_probe_tangents_added": 3,
    "coefficient_selected": False,
    "parameter_scan_performed": False,
    "optimization_performed": False,
    "post_result_retuning_performed": False,
    "pressure_or_force_changed": False,
    "held_out_pde_residual_evaluated": False,
    "st006_comparison_performed": False,
    "public_image_numeric_target_used": False,
    "pixel_similarity_objective_used": False,
    "renderer_or_camera_fit_used": False,
    "direct_visualization_fingerprint_improvement": 0.0,
    "closer_visualization_delivery_established": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "source_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _scientific_parent():
    return parent._scientific_parent()


def _assert_source_lock() -> None:
    parent._assert_source_lock()
    if parent.TASK_ID != "CR003-ST052M-DEFECT-SPACE-PROJECTION-135":
        raise RuntimeError("#1139 parent task identity drifted")
    if parent.PREREG_ISSUE != 1138:
        raise RuntimeError("#1139 preregistration identity drifted")
    if parent.STACK_BASE_PR != 1131:
        raise RuntimeError("#1139 stack lineage drifted")
    if tuple(parent.CHANNELS) != CHANNELS:
        raise RuntimeError("existing-control order drifted")
    if tuple(parent.DEFECT_COORDINATES) != DEFECT_COORDINATES:
        raise RuntimeError("defect-coordinate order drifted")
    scientific = _scientific_parent()
    if scientific.TASK_ID != "CR003-ST052M-DEFECT-COORDINATE-CONTROLLABILITY-128":
        raise RuntimeError("#1078 scientific-parent identity drifted")
    if tuple(scientific.EPSILONS) != (1.0e-3, 5.0e-4):
        raise RuntimeError("#1078 finite-difference steps drifted")
    if scientific.PRIMARY_EPSILON != 5.0e-4:
        raise RuntimeError("#1078 primary finite-difference step drifted")


def _activation_g1(time: float) -> float:
    return 2.0 * (float(time) - 0.25)


def _activation_skew(time: float) -> float:
    tau = (float(time) - 0.25) / 0.50
    return 4.0 * tau * (1.0 - tau) * (2.0 * tau - 1.0)


def _inside(points: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    pts = np.asarray(points, dtype=float)
    if pts.ndim != 2 or pts.shape[1] != 3:
        raise ValueError("points must have shape (N,3)")
    x, y, z = pts[:, 0], pts[:, 1], pts[:, 2]
    q = x * x + y * y
    mask = (q < SUPPORT_RADIUS * SUPPORT_RADIUS) & (np.abs(z) < SUPPORT_ABS_Z)
    return x, y, z, mask


def _poloidal_unit(points: np.ndarray) -> np.ndarray:
    """curl(-y f, x f, 0) for the frozen compact odd scalar f."""
    x, y, z, mask = _inside(points)
    out = np.zeros((x.size, 3), dtype=float)
    if not np.any(mask):
        return out

    xm, ym, zm = x[mask], y[mask], z[mask]
    q = xm * xm + ym * ym
    r2 = SUPPORT_RADIUS * SUPPORT_RADIUS
    z2 = SUPPORT_ABS_Z * SUPPORT_ABS_Z
    a = 1.0 - q / r2
    b = 1.0 - (zm * zm) / z2
    radial = a**4
    axial = b**4
    radial_q = -(4.0 / r2) * a**3
    axial_z = -(8.0 * zm / z2) * b**3

    f = zm * radial * axial
    f_z = radial * (axial + zm * axial_z)
    f_q = zm * radial_q * axial

    out[mask, 0] = -xm * f_z
    out[mask, 1] = -ym * f_z
    out[mask, 2] = 2.0 * (f + q * f_q)
    return out


def _toroidal_unit(points: np.ndarray) -> np.ndarray:
    """curl(0,0,phi) for the frozen compact even scalar phi."""
    x, y, z, mask = _inside(points)
    out = np.zeros((x.size, 3), dtype=float)
    if not np.any(mask):
        return out

    xm, ym, zm = x[mask], y[mask], z[mask]
    q = xm * xm + ym * ym
    r2 = SUPPORT_RADIUS * SUPPORT_RADIUS
    z2 = SUPPORT_ABS_Z * SUPPORT_ABS_Z
    a = 1.0 - q / r2
    b = 1.0 - (zm * zm) / z2
    phi_q = -(5.0 / r2) * a**4 * b**5

    out[mask, 0] = 2.0 * ym * phi_q
    out[mask, 1] = -2.0 * xm * phi_q
    return out


def _probe_tangents() -> dict[str, Callable[[np.ndarray, float], np.ndarray]]:
    def poloidal(points: np.ndarray, time: float) -> np.ndarray:
        return _activation_g1(time) * _poloidal_unit(points)

    def toroidal(points: np.ndarray, time: float) -> np.ndarray:
        return _activation_g1(time) * _toroidal_unit(points)

    def time_skewed(points: np.ndarray, time: float) -> np.ndarray:
        return _activation_skew(time) * _poloidal_unit(points)

    return {
        "poloidal_scalar_space_extension": poloidal,
        "compact_toroidal_envelope_extension": toroidal,
        "time_skewed_poloidal_extension": time_skewed,
    }


def _divergence_fd(
    tangent: Callable[[np.ndarray, float], np.ndarray],
    points: np.ndarray,
    time: float,
) -> np.ndarray:
    pts = np.asarray(points, dtype=float)
    h = DIVERGENCE_FD_STEP
    div = np.zeros(pts.shape[0], dtype=float)
    for axis in range(3):
        plus = pts.copy()
        minus = pts.copy()
        plus[:, axis] += h
        minus[:, axis] -= h
        vp = np.asarray(tangent(plus, float(time)), dtype=float)
        vm = np.asarray(tangent(minus, float(time)), dtype=float)
        div += (vp[:, axis] - vm[:, axis]) / (2.0 * h)
    return div


def _structural_guard(
    tangent: Callable[[np.ndarray, float], np.ndarray]
) -> dict[str, Any]:
    interior = np.asarray(
        [
            [0.20, 0.10, 0.40],
            [-0.35, 0.25, -0.70],
            [0.55, -0.40, 0.90],
            [0.0, 0.0, 0.50],
            [0.75, 0.20, -1.10],
        ],
        dtype=float,
    )
    divergence_values = []
    for time in (0.375, 0.625):
        divergence_values.extend(
            np.abs(_divergence_fd(tangent, interior, time)).tolist()
        )
    divergence_max = float(max(divergence_values))

    outside = np.asarray(
        [
            [SUPPORT_RADIUS + 0.05, 0.0, 0.0],
            [0.0, SUPPORT_RADIUS + 0.05, 0.0],
            [0.2, 0.1, SUPPORT_ABS_Z + 0.05],
            [-0.2, 0.1, -SUPPORT_ABS_Z - 0.05],
        ],
        dtype=float,
    )
    outside_max = 0.0
    for time in (0.25, 0.50, 0.75):
        outside_max = max(
            outside_max,
            float(np.max(np.abs(np.asarray(tangent(outside, time), dtype=float)))),
        )

    axis = np.asarray(
        [[0.0, 0.0, -1.0], [0.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
        dtype=float,
    )
    axis_transverse_max = 0.0
    axis_finite = True
    for time in (0.375, 0.625):
        values = np.asarray(tangent(axis, time), dtype=float)
        axis_finite = axis_finite and bool(np.all(np.isfinite(values)))
        axis_transverse_max = max(
            axis_transverse_max, float(np.max(np.abs(values[:, :2])))
        )

    passes = bool(
        divergence_max <= DIVERGENCE_FD_MAX
        and outside_max <= OUTSIDE_SUPPORT_MAX
        and axis_finite
        and axis_transverse_max <= AXIS_TRANSVERSE_MAX
    )
    return {
        "centered_divergence_fd_step": DIVERGENCE_FD_STEP,
        "centered_divergence_max_abs": divergence_max,
        "centered_divergence_gate": DIVERGENCE_FD_MAX,
        "outside_support_max_abs_velocity": outside_max,
        "outside_support_gate": OUTSIDE_SUPPORT_MAX,
        "axis_transverse_max_abs_velocity": axis_transverse_max,
        "axis_transverse_gate": AXIS_TRANSVERSE_MAX,
        "axis_values_finite": axis_finite,
        "passes": passes,
    }


def _column_normalize(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    m = np.asarray(matrix, dtype=float)
    norms = np.linalg.norm(m, axis=0)
    out = np.zeros_like(m)
    good = norms > np.finfo(float).tiny
    out[:, good] = m[:, good] / norms[good]
    return out, norms


def _score_probe(
    existing_balanced: np.ndarray,
    projector: np.ndarray,
    null_projector: np.ndarray,
    row_norms: np.ndarray,
    coarse: np.ndarray,
    fine: np.ndarray,
    structural: dict[str, Any],
) -> dict[str, Any]:
    existing = np.asarray(existing_balanced, dtype=float)
    p = np.asarray(projector, dtype=float)
    q = np.asarray(null_projector, dtype=float)
    rn = np.asarray(row_norms, dtype=float)
    vc = np.asarray(coarse, dtype=float)
    vf = np.asarray(fine, dtype=float)
    if existing.shape != (len(DEFECT_COORDINATES), len(CHANNELS)):
        raise RuntimeError("unexpected existing balanced-Jacobian shape")
    if p.shape != (len(DEFECT_COORDINATES),) * 2 or q.shape != p.shape:
        raise RuntimeError("unexpected projector shape")
    if rn.shape != (len(DEFECT_COORDINATES),) or np.any(rn <= ROW_RESPONSE_FLOOR):
        raise RuntimeError("invalid frozen row scales")
    if vf.shape != rn.shape or vc.shape != rn.shape:
        raise RuntimeError("unexpected probe derivative shape")

    fine_raw_norm = float(np.linalg.norm(vf))
    coarse_raw_norm = float(np.linalg.norm(vc))
    if fine_raw_norm <= ROW_RESPONSE_FLOOR:
        drift = None
        derivative_cosine = None
        stable = False
    else:
        drift = float(np.linalg.norm(vf - vc) / fine_raw_norm)
        derivative_cosine = (
            None
            if coarse_raw_norm <= np.finfo(float).tiny
            else float(np.dot(vf, vc) / (fine_raw_norm * coarse_raw_norm))
        )
        stable = bool(
            drift <= DERIVATIVE_DRIFT_MAX
            and derivative_cosine is not None
            and derivative_cosine >= DERIVATIVE_COSINE_MIN
        )

    balanced = vf / rn
    balanced_norm = float(np.linalg.norm(balanced))
    if balanced_norm <= ROW_RESPONSE_FLOOR:
        unit = np.zeros_like(balanced)
        span_fraction = 0.0
        null_fraction = 0.0
    else:
        unit = balanced / balanced_norm
        span_fraction = float(np.linalg.norm(p @ unit))
        null_fraction = float(np.linalg.norm(q @ unit))

    existing_normed, existing_norms = _column_normalize(existing)
    augmented = np.column_stack((existing, balanced))
    augmented_normed, _augmented_norms = _column_normalize(augmented)
    existing_rank = int(np.linalg.matrix_rank(existing_normed, tol=RANK_TOL))
    augmented_rank = int(np.linalg.matrix_rank(augmented_normed, tol=RANK_TOL))
    singular = np.linalg.svd(augmented_normed, compute_uv=False)
    condition = (
        None
        if singular.size == 0 or singular[-1] <= np.finfo(float).tiny
        else float(singular[0] / singular[-1])
    )
    if balanced_norm <= np.finfo(float).tiny:
        cosines = {channel: None for channel in CHANNELS}
    else:
        cosines = {
            channel: (
                None
                if existing_norms[i] <= np.finfo(float).tiny
                else float(np.dot(existing_normed[:, i], unit))
            )
            for i, channel in enumerate(CHANNELS)
        }
    finite_cos = [abs(v) for v in cosines.values() if v is not None]
    max_abs_cosine = None if not finite_cos else float(max(finite_cos))

    order = np.argsort(-np.abs(unit))
    dominant = [
        {
            "coordinate": DEFECT_COORDINATES[int(i)],
            "row_balanced_unit_derivative": float(unit[int(i)]),
            "fine_raw_derivative": float(vf[int(i)]),
        }
        for i in order[:4]
    ]

    eligible = bool(
        stable
        and structural["passes"]
        and augmented_rank == existing_rank + 1
        and condition is not None
        and condition <= CONDITION_MAX
        and max_abs_cosine is not None
        and max_abs_cosine < COSINE_MAX_ABS
    )
    return {
        "fine_raw_derivative_l2": fine_raw_norm,
        "coarse_raw_derivative_l2": coarse_raw_norm,
        "relative_derivative_drift": drift,
        "coarse_fine_cosine": derivative_cosine,
        "derivative_stable": stable,
        "fine_raw_derivative_by_coordinate": {
            name: float(vf[i]) for i, name in enumerate(DEFECT_COORDINATES)
        },
        "row_balanced_unit_derivative_by_coordinate": {
            name: float(unit[i]) for i, name in enumerate(DEFECT_COORDINATES)
        },
        "existing_span_fraction_l2": span_fraction,
        "left_null_novelty_fraction_l2": null_fraction,
        "existing_rank": existing_rank,
        "augmented_rank": augmented_rank,
        "rank_gain": int(augmented_rank - existing_rank),
        "augmented_normalized_column_singular_values": [float(v) for v in singular],
        "augmented_normalized_column_condition": condition,
        "probe_cosine_with_existing_normalized_controls": cosines,
        "max_abs_probe_existing_cosine": max_abs_cosine,
        "dominant_balanced_defect_coordinates": dominant,
        "structural_guard": structural,
        "independent_preflight_eligible": eligible,
    }


def build_report() -> dict[str, Any]:
    _assert_source_lock()
    scientific = _scientific_parent()
    observable_parent = scientific.parent
    base_fn, _existing_tangents, base_metadata = observable_parent._build_field_and_tangents()
    baseline_speed = observable_parent._baseline_speed_scale(base_fn)

    scientific_audit = scientific.controllability_audit()
    projection = parent.projection_audit_from_parent(scientific_audit)
    raw_existing = parent._raw_jacobian(scientific_audit)
    row_norms = np.asarray(
        [projection["raw_row_l2"][name] for name in DEFECT_COORDINATES], dtype=float
    )
    existing_balanced = raw_existing / row_norms[:, None]
    projector = np.asarray(projection["existing_control_projector"], dtype=float)
    null_projector = np.asarray(projection["left_null_projector"], dtype=float)

    tangents = _probe_tangents()
    probes: dict[str, Any] = {}
    for name in PROBE_NAMES:
        tangent = tangents[name]
        derivatives = {
            float(epsilon): scientific._coordinate_derivative(
                base_fn, tangent, float(epsilon), baseline_speed
            )
            for epsilon in scientific.EPSILONS
        }
        structural = _structural_guard(tangent)
        probes[name] = _score_probe(
            existing_balanced,
            projector,
            null_projector,
            row_norms,
            derivatives[float(scientific.EPSILONS[0])],
            derivatives[float(scientific.PRIMARY_EPSILON)],
            structural,
        )

    eligible = [name for name in PROBE_NAMES if probes[name]["independent_preflight_eligible"]]
    ranked = sorted(
        eligible,
        key=lambda name: (
            probes[name]["left_null_novelty_fraction_l2"],
            -float(probes[name]["augmented_normalized_column_condition"]),
        ),
        reverse=True,
    )

    return {
        "task_id": TASK_ID,
        "prereg_issue": PREREG_ISSUE,
        "stack_base": {"pr": STACK_BASE_PR, "head": STACK_BASE_HEAD},
        "scientific_parent": {
            "pr": parent.SCIENTIFIC_PARENT_PR,
            "task": parent.SCIENTIFIC_PARENT_TASK,
        },
        "frozen_protocol": {
            "probe_names": list(PROBE_NAMES),
            "streamfunction_vector_potential_same_poloidal_class": True,
            "streamfunction_vector_potential_prior_reference_pr": 828,
            "support_radius": SUPPORT_RADIUS,
            "support_abs_z": SUPPORT_ABS_Z,
            "poloidal_scalar": (
                "f=z*(1-r^2/R^2)^4*(1-z^2/Z^2)^4; "
                "A=(-y*f,x*f,0); psi=r^2*f"
            ),
            "toroidal_scalar": (
                "phi=(1-r^2/R^2)^5*(1-z^2/Z^2)^5; A=(0,0,phi)"
            ),
            "linear_activation": "g1(t)=2*(t-.25)",
            "time_skew_activation": (
                "h(t)=4*tau*(1-tau)*(2*tau-1), tau=(t-.25)/.5"
            ),
            "finite_difference_epsilons": [float(v) for v in scientific.EPSILONS],
            "primary_epsilon": float(scientific.PRIMARY_EPSILON),
            "derivative_drift_max": DERIVATIVE_DRIFT_MAX,
            "derivative_cosine_min": DERIVATIVE_COSINE_MIN,
            "condition_max": CONDITION_MAX,
            "cosine_max_abs": COSINE_MAX_ABS,
            "rank_tolerance": RANK_TOL,
            "row_response_floor": ROW_RESPONSE_FLOOR,
            "divergence_fd_max": DIVERGENCE_FD_MAX,
            "outside_support_max": OUTSIDE_SUPPORT_MAX,
            "axis_transverse_max": AXIS_TRANSVERSE_MAX,
            "public_openai_numeric_target": None,
            "coefficient_selected": None,
            "candidate_basis_dimension_added": 0,
        },
        "base_metadata": base_metadata,
        "existing_defect_tangent": {
            "rank": projection["rank"],
            "left_nullity": projection["left_nullity"],
            "row_balancing": "reuse #1139 raw-row L2 balancing",
        },
        "probe_results": probes,
        "eligible_probe_archetypes": eligible,
        "eligible_probe_archetypes_ranked_by_left_null_novelty": ranked,
        "decision": {
            "basis_candidate_promoted": None,
            "candidate_mutation_authorized": False,
            "coefficient_selection_authorized": False,
            "interpretation": (
                "eligible probes only identify structurally admissible local directions "
                "that add morphology-tangent capacity. A real candidate-bound public-"
                "observable discrepancy must first project materially onto the same "
                "left-null direction before a one-basis child is justified."
            ),
            "if_none_eligible": (
                "redesign or reparameterize one probe against the exposed left-null "
                "defect combination; do not perform unrestricted basis growth"
            ),
            "actual_velocity_changed": False,
            "direct_visualization_fingerprint_improvement": 0.0,
            "closer_visualization_delivery_established": False,
        },
        "held_out_pde_residual": {
            "evaluated": False,
            "reason": (
                "analysis-only tangent preflight; candidate velocity and compatible "
                "pressure/restricted forcing are unchanged"
            ),
        },
        **TRUTH,
    }


def run(out: Path) -> dict[str, Any]:
    report = build_report()
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    report = run(args.out)
    summary = {
        "existing_defect_tangent": report["existing_defect_tangent"],
        "eligible_probe_archetypes": report["eligible_probe_archetypes"],
        "eligible_probe_archetypes_ranked_by_left_null_novelty": report[
            "eligible_probe_archetypes_ranked_by_left_null_novelty"
        ],
        "probe_results": {
            name: {
                "left_null_novelty_fraction_l2": value[
                    "left_null_novelty_fraction_l2"
                ],
                "rank_gain": value["rank_gain"],
                "augmented_normalized_column_condition": value[
                    "augmented_normalized_column_condition"
                ],
                "max_abs_probe_existing_cosine": value[
                    "max_abs_probe_existing_cosine"
                ],
                "derivative_stable": value["derivative_stable"],
                "structural_guard_passes": value["structural_guard"]["passes"],
                "independent_preflight_eligible": value[
                    "independent_preflight_eligible"
                ],
            }
            for name, value in report["probe_results"].items()
        },
        "decision": report["decision"],
    }
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
