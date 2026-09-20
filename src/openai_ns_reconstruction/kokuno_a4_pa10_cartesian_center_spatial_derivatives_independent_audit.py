"""Independent Agent-4 audit of PA.10 Cartesian-center spatial derivatives.

The object under audit is Agent-1 #819's public analytic
``velocity_jacobian``, ``divergence`` and ``vorticity`` for the executable
PA.10 inner contraction-center velocity.  The independent numerical reference
never calls those analytic spatial-derivative paths.  It differentiates only
the public ``velocity(x,y,z,t)`` callable with a centered sixth-order Cartesian
finite-difference operator on three preregistered spatial steps.

This is deliberately an inner-field audit, not a complete Navier--Stokes
validation.  No corrected fixed point, outer/global join, matched pressure,
restricted forcing, oscillatory/correction composite, or full momentum
residual is inferred here.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_pa10_cartesian_center_spatial_derivatives import (
    KokunoPA10CartesianCenterSpatialDerivatives,
)


SCHEMA = "kokuno-a4-pa10-cartesian-center-spatial-derivatives-independent-audit-v1"
SEED = 9173481
SAMPLE_COUNT = 384
SPACE_STEPS = (1.2e-3, 6.0e-4, 3.0e-4)
JACOBIAN_FINE_RELATIVE_RMS_GATE = 2.0e-6
JACOBIAN_FINE_RELATIVE_MAX_GATE = 2.0e-5
VORTICITY_FINE_RELATIVE_RMS_GATE = 3.0e-6
VORTICITY_FINE_RELATIVE_MAX_GATE = 3.0e-5
REFINEMENT_RATIO_GATE = 8.0
REFINEMENT_FLOOR = 2.0e-9
DIVERGENCE_MAX_GATE = 1.0e-5
DIVERGENCE_L2_GATE = 1.0e-5
NONTRIVIAL_JACOBIAN_RMS_FLOOR = 1.0e-8
NONTRIVIAL_VORTICITY_RMS_FLOOR = 1.0e-8
PUBLIC_CLOSURE_GATE = 1.0e-13
FINAL_MOMENTUM_GATE = 1.0e-3
FINAL_DIVERGENCE_GATE = 1.0e-5

_BASE_TRUTH_BOUNDARY = {
    "agent4_public_velocity_only_fd6_spatial_derivative_used": True,
    "agent4_fresh_offgrid_fixed_cartesian_samples_used": True,
    "agent4_three_spatial_resolutions_used": True,
    "inner_cartesian_center_spatial_derivative_protocol_executed": True,
    "inner_cartesian_center_spatial_derivatives_independently_admitted": False,
    "source_complex_C_normalization_independently_admitted": False,
    "source_center_is_final_corrected_fixed_point": False,
    "fixed_point_correction_materialized": False,
    "global_cartesian_spacetime_leading_velocity_materialized": False,
    "outer_join_localization_materialized": False,
    "matched_global_pressure_materialized": False,
    "complete_restricted_forcing_materialized": False,
    "complete_kokuno_composite_velocity": False,
    "leading_only_ns_residual_assessed": False,
    "leading_plus_oscillatory_ns_residual_assessed": False,
    "after_correction_ns_residual_assessed": False,
    "heldout_ns_momentum_residual_assessed": False,
    "formal_full_domain_pde_gate_assessed": False,
    "pde_validated": False,
}


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _relative_rms(actual: np.ndarray, reference: np.ndarray) -> float:
    actual_arr = np.asarray(actual, dtype=float)
    reference_arr = np.asarray(reference, dtype=float)
    numerator = float(np.sqrt(np.mean(np.square(actual_arr - reference_arr))))
    denominator = max(float(np.sqrt(np.mean(np.square(reference_arr)))), 1.0e-14)
    return numerator / denominator


def _relative_max(actual: np.ndarray, reference: np.ndarray) -> float:
    actual_arr = np.asarray(actual, dtype=float)
    reference_arr = np.asarray(reference, dtype=float)
    numerator = float(np.max(np.abs(actual_arr - reference_arr)))
    denominator = max(float(np.max(np.abs(reference_arr))), 1.0e-14)
    return numerator / denominator


def _curl_from_jacobian(jacobian: np.ndarray) -> np.ndarray:
    jac = np.asarray(jacobian, dtype=float)
    return np.stack(
        (
            jac[..., 2, 1] - jac[..., 1, 2],
            jac[..., 0, 2] - jac[..., 2, 0],
            jac[..., 1, 0] - jac[..., 0, 1],
        ),
        axis=-1,
    )


def _sample_fixed_cartesian_points(
    provider: KokunoPA10CartesianCenterSpatialDerivatives,
    rng: np.random.Generator,
) -> dict[str, np.ndarray]:
    """Generate fresh source-interior points without using the inverse map."""
    _, xmax = provider.source_X_interval
    X = rng.uniform(0.10 * xmax, 0.28 * xmax, size=SAMPLE_COUNT)
    eta = rng.uniform(-0.34, 0.34, size=SAMPLE_COUNT)
    t = rng.uniform(0.38, 0.62, size=SAMPLE_COUNT)
    theta = rng.uniform(-math.pi, math.pi, size=SAMPLE_COUNT)

    q = (1.0 - t) / (1.0 - eta * eta)
    z = np.power(q, provider.D) * eta
    radius = np.sqrt(2.0 * q * X)
    x = radius * np.cos(theta)
    y = radius * np.sin(theta)

    # Three exact-axis and two axis-near probes.  They are diagnostics for
    # regularity; no nonzero-radius value is incorrectly required to be zero.
    axis_eta = np.asarray([-0.30, -1.0e-12, 0.0, 1.0e-12, 0.30], dtype=float)
    axis_t = np.asarray([0.40, 0.45, 0.50, 0.55, 0.60], dtype=float)
    axis_q = (1.0 - axis_t) / (1.0 - axis_eta * axis_eta)
    axis_z = np.power(axis_q, provider.D) * axis_eta
    axis_radius = np.asarray([0.0, 1.0e-10, 0.0, 1.0e-10, 0.0], dtype=float)
    axis_theta = np.asarray([0.0, 0.37, 0.0, -0.81, 0.0], dtype=float)
    axis_x = axis_radius * np.cos(axis_theta)
    axis_y = axis_radius * np.sin(axis_theta)

    return {
        "x": np.concatenate((x, axis_x)),
        "y": np.concatenate((y, axis_y)),
        "z": np.concatenate((z, axis_z)),
        "t": np.concatenate((t, axis_t)),
        "axis_probe_radius": axis_radius,
    }


def _fd6_axis(
    provider: Any,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
    axis: int,
    step: float,
) -> np.ndarray:
    """Centered sixth-order derivative of the public velocity only."""
    coords = [
        np.asarray(x, dtype=float),
        np.asarray(y, dtype=float),
        np.asarray(z, dtype=float),
    ]

    def shifted(multiplier: float) -> np.ndarray:
        moved = [value.copy() for value in coords]
        moved[axis] = moved[axis] + multiplier * step
        return np.asarray(
            provider.velocity(moved[0], moved[1], moved[2], t), dtype=float
        )

    return (
        -shifted(-3.0)
        + 9.0 * shifted(-2.0)
        - 45.0 * shifted(-1.0)
        + 45.0 * shifted(1.0)
        - 9.0 * shifted(2.0)
        + shifted(3.0)
    ) / (60.0 * step)


def _fd6_jacobian(
    provider: Any,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
    step: float,
) -> np.ndarray:
    columns = [
        _fd6_axis(provider, x, y, z, t, axis=axis, step=step)
        for axis in range(3)
    ]
    return np.stack(columns, axis=-1)


def _component_errors(actual: np.ndarray, reference: np.ndarray) -> dict[str, Any]:
    component_names = ("u", "v", "w")
    axis_names = ("x", "y", "z")
    result: dict[str, Any] = {}
    for component in range(3):
        for axis in range(3):
            key = f"d{component_names[component]}_d{axis_names[axis]}"
            result[key] = {
                "relative_rms": _relative_rms(
                    actual[:, component, axis], reference[:, component, axis]
                ),
                "relative_max": _relative_max(
                    actual[:, component, axis], reference[:, component, axis]
                ),
                "absolute_max": float(
                    np.max(
                        np.abs(
                            actual[:, component, axis]
                            - reference[:, component, axis]
                        )
                    )
                ),
            }
    return result


@dataclass(frozen=True)
class AuditResult:
    payload: dict[str, Any]

    @property
    def passed(self) -> bool:
        return bool(self.payload["passed"])


def run_independent_audit() -> AuditResult:
    provider = KokunoPA10CartesianCenterSpatialDerivatives()
    rng = np.random.default_rng(SEED)
    samples = _sample_fixed_cartesian_points(provider, rng)
    x = np.asarray(samples["x"], dtype=float)
    y = np.asarray(samples["y"], dtype=float)
    z = np.asarray(samples["z"], dtype=float)
    t = np.asarray(samples["t"], dtype=float)

    public_jacobian = np.asarray(provider.velocity_jacobian(x, y, z, t), dtype=float)
    public_divergence = np.asarray(provider.divergence(x, y, z, t), dtype=float)
    public_vorticity = np.asarray(provider.vorticity(x, y, z, t), dtype=float)
    if public_jacobian.shape != (x.size, 3, 3):
        raise RuntimeError("public velocity_jacobian returned malformed data")
    if public_divergence.shape != (x.size,):
        raise RuntimeError("public divergence returned malformed data")
    if public_vorticity.shape != (x.size, 3):
        raise RuntimeError("public vorticity returned malformed data")
    if not (
        np.all(np.isfinite(public_jacobian))
        and np.all(np.isfinite(public_divergence))
        and np.all(np.isfinite(public_vorticity))
    ):
        raise RuntimeError("public spatial derivative surface became nonfinite")

    public_trace = np.trace(public_jacobian, axis1=-2, axis2=-1)
    public_curl = _curl_from_jacobian(public_jacobian)
    public_divergence_closure = float(np.max(np.abs(public_divergence - public_trace)))
    public_vorticity_closure = float(np.max(np.abs(public_vorticity - public_curl)))

    estimates: list[np.ndarray] = []
    by_step: dict[str, Any] = {}
    for step in SPACE_STEPS:
        estimate = _fd6_jacobian(provider, x, y, z, t, step)
        estimates.append(estimate)
        fd_divergence = np.trace(estimate, axis1=-2, axis2=-1)
        fd_vorticity = _curl_from_jacobian(estimate)
        by_step[f"{step:.7f}"] = {
            "jacobian_relative_rms": _relative_rms(public_jacobian, estimate),
            "jacobian_relative_max": _relative_max(public_jacobian, estimate),
            "vorticity_relative_rms": _relative_rms(public_vorticity, fd_vorticity),
            "vorticity_relative_max": _relative_max(public_vorticity, fd_vorticity),
            "divergence_sampled_max": float(np.max(np.abs(fd_divergence))),
            "divergence_sampled_l2_rms": float(
                np.sqrt(np.mean(np.square(fd_divergence)))
            ),
            "component_errors": _component_errors(public_jacobian, estimate),
        }

    coarse, medium, fine = estimates
    coarse_medium_difference_rms = float(
        np.sqrt(np.mean(np.square(coarse - medium)))
    )
    medium_fine_difference_rms = float(
        np.sqrt(np.mean(np.square(medium - fine)))
    )
    refinement_ratio = coarse_medium_difference_rms / max(
        medium_fine_difference_rms, np.finfo(float).tiny
    )
    refinement_stable = bool(
        medium_fine_difference_rms <= REFINEMENT_FLOOR
        or refinement_ratio >= REFINEMENT_RATIO_GATE
    )

    fine_vorticity = _curl_from_jacobian(fine)
    fine_divergence = np.trace(fine, axis1=-2, axis2=-1)
    jacobian_relative_rms = _relative_rms(public_jacobian, fine)
    jacobian_relative_max = _relative_max(public_jacobian, fine)
    vorticity_relative_rms = _relative_rms(public_vorticity, fine_vorticity)
    vorticity_relative_max = _relative_max(public_vorticity, fine_vorticity)
    divergence_max = float(np.max(np.abs(fine_divergence)))
    divergence_l2 = float(np.sqrt(np.mean(np.square(fine_divergence))))
    jacobian_rms = float(np.sqrt(np.mean(np.square(public_jacobian))))
    vorticity_rms = float(np.sqrt(np.mean(np.square(public_vorticity))))

    jacobian_error = np.abs(public_jacobian - fine)
    flat_index = int(np.argmax(jacobian_error))
    point_index, component_index, axis_index = np.unravel_index(
        flat_index, jacobian_error.shape
    )
    worst = {
        "point_index": int(point_index),
        "component": ("u", "v", "w")[component_index],
        "axis": ("x", "y", "z")[axis_index],
        "x": float(x[point_index]),
        "y": float(y[point_index]),
        "z": float(z[point_index]),
        "t": float(t[point_index]),
        "public": float(public_jacobian[point_index, component_index, axis_index]),
        "independent_fd6": float(fine[point_index, component_index, axis_index]),
        "absolute_error": float(jacobian_error[point_index, component_index, axis_index]),
    }

    axis_start = SAMPLE_COUNT
    axis_jacobian = public_jacobian[axis_start:]
    axis_vorticity = public_vorticity[axis_start:]
    axis_all_finite = bool(
        np.all(np.isfinite(axis_jacobian)) and np.all(np.isfinite(axis_vorticity))
    )

    # Frozen negative controls.  The same gates must reject these mutations.
    scaled_jacobian_relative_max = _relative_max(0.99 * public_jacobian, fine)
    flipped_vorticity_relative_max = _relative_max(-public_vorticity, fine_vorticity)
    divergence_mutation = np.array(public_jacobian, copy=True)
    divergence_mutation[:, 0, 0] += 1.0e-3
    divergence_mutation_max = float(
        np.max(np.abs(np.trace(divergence_mutation, axis1=-2, axis2=-1)))
    )
    mutations_detected = bool(
        scaled_jacobian_relative_max > JACOBIAN_FINE_RELATIVE_MAX_GATE
        and flipped_vorticity_relative_max > VORTICITY_FINE_RELATIVE_MAX_GATE
        and divergence_mutation_max > DIVERGENCE_MAX_GATE
    )

    gates = {
        "jacobian_fine_relative_rms": (
            jacobian_relative_rms <= JACOBIAN_FINE_RELATIVE_RMS_GATE
        ),
        "jacobian_fine_relative_max": (
            jacobian_relative_max <= JACOBIAN_FINE_RELATIVE_MAX_GATE
        ),
        "vorticity_fine_relative_rms": (
            vorticity_relative_rms <= VORTICITY_FINE_RELATIVE_RMS_GATE
        ),
        "vorticity_fine_relative_max": (
            vorticity_relative_max <= VORTICITY_FINE_RELATIVE_MAX_GATE
        ),
        "three_level_refinement_stable": refinement_stable,
        "independent_divergence_sampled_max": divergence_max <= DIVERGENCE_MAX_GATE,
        "independent_divergence_sampled_l2_rms": divergence_l2 <= DIVERGENCE_L2_GATE,
        "jacobian_nontrivial": jacobian_rms >= NONTRIVIAL_JACOBIAN_RMS_FLOOR,
        "vorticity_nontrivial": vorticity_rms >= NONTRIVIAL_VORTICITY_RMS_FLOOR,
        "public_divergence_trace_closure": public_divergence_closure <= PUBLIC_CLOSURE_GATE,
        "public_vorticity_curl_closure": public_vorticity_closure <= PUBLIC_CLOSURE_GATE,
        "axis_spatial_derivatives_finite": axis_all_finite,
        "mutations_detected": mutations_detected,
    }
    passed = bool(all(gates.values()))

    truth_boundary = dict(_BASE_TRUTH_BOUNDARY)
    truth_boundary["inner_cartesian_center_spatial_derivatives_independently_admitted"] = passed

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "seed": SEED,
        "sample_count_random": SAMPLE_COUNT,
        "sample_count_total": int(x.size),
        "space_steps": list(SPACE_STEPS),
        "frozen_gates": {
            "jacobian_fine_relative_rms": JACOBIAN_FINE_RELATIVE_RMS_GATE,
            "jacobian_fine_relative_max": JACOBIAN_FINE_RELATIVE_MAX_GATE,
            "vorticity_fine_relative_rms": VORTICITY_FINE_RELATIVE_RMS_GATE,
            "vorticity_fine_relative_max": VORTICITY_FINE_RELATIVE_MAX_GATE,
            "refinement_ratio_or_floor": {
                "ratio_min": REFINEMENT_RATIO_GATE,
                "fine_difference_floor": REFINEMENT_FLOOR,
            },
            "independent_divergence_sampled_max": DIVERGENCE_MAX_GATE,
            "independent_divergence_sampled_l2_rms": DIVERGENCE_L2_GATE,
            "final_normalized_momentum_max_l2": FINAL_MOMENTUM_GATE,
            "final_divergence_max_l2": FINAL_DIVERGENCE_GATE,
        },
        "measurements": {
            "jacobian_rms": jacobian_rms,
            "vorticity_rms": vorticity_rms,
            "by_step": by_step,
            "jacobian_fine_relative_rms": jacobian_relative_rms,
            "jacobian_fine_relative_max": jacobian_relative_max,
            "vorticity_fine_relative_rms": vorticity_relative_rms,
            "vorticity_fine_relative_max": vorticity_relative_max,
            "independent_divergence_sampled_max": divergence_max,
            "independent_divergence_sampled_l2_rms": divergence_l2,
            "coarse_medium_difference_rms": coarse_medium_difference_rms,
            "medium_fine_difference_rms": medium_fine_difference_rms,
            "refinement_ratio": refinement_ratio,
            "public_divergence_trace_closure_abs_max": public_divergence_closure,
            "public_vorticity_curl_closure_abs_max": public_vorticity_closure,
            "axis_probe_radius": np.asarray(samples["axis_probe_radius"]).tolist(),
            "worst_jacobian_error": worst,
            "scaled_0p99_jacobian_mutation_relative_max": scaled_jacobian_relative_max,
            "flipped_vorticity_mutation_relative_max": flipped_vorticity_relative_max,
            "plus_1e_minus_3_du_dx_divergence_mutation_max": divergence_mutation_max,
        },
        "gates": gates,
        "truth_boundary": truth_boundary,
        "limitations": [
            "This is an inner PA.10 contraction-center spatial-derivative audit, not a global leading-field admission.",
            "The independent reference differentiates only public velocity values; it does not independently replay Agent-1 coordinate/profile derivative algebra.",
            "The reported divergence L2 quantity is a held-out sampled RMS, not a whole-domain volume integral norm.",
            "No pressure, restricted forcing, outer join, oscillatory/correction composite, or full NS momentum residual is evaluated.",
            "The source-complex C normalization is not independently admitted by this audit.",
            "No source parameter, finite-difference step, sample set, or gate may be retuned from these measurements.",
        ],
        "passed": passed,
    }
    payload["receipt_sha256"] = hashlib.sha256(
        _canonical_json(payload).encode("utf-8")
    ).hexdigest()
    return AuditResult(payload)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)

    result = run_independent_audit()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result.payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(result.payload, sort_keys=True, allow_nan=False))
    return 0 if result.passed else 1


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
