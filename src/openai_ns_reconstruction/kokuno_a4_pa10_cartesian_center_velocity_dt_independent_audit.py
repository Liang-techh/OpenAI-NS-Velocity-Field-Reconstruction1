"""Independent Agent-4 audit of the PA.10 Cartesian-center time derivative.

The object under audit is Agent-1 #811's public
``KokunoPA10CartesianCenterVelocityTimeDerivative.velocity_dt``.  This module
never calls Agent-1's coordinate-time-derivative helper, physical-profile
analytic derivatives, ``v0_eta`` helper, or report as a numerical oracle.

Instead, fresh physical points are generated from the displayed source forward
coordinates and the derivative is reconstructed only from the public
``velocity(x,y,z,t)`` callable with a fixed-Cartesian centered FD4 operator on
three preregistered time steps.  This is therefore an implementation-distinct
check of the new ``u_t`` seam needed by a future independent NS residual.

Scope remains deliberately narrow: the audited velocity is only the PA.10
inner contraction center.  No corrected fixed point, outer/global join,
matched pressure, restricted forcing, complete leading+oscillatory+correction
candidate, or full Navier--Stokes momentum residual is inferred here.
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

from .kokuno_pa10_cartesian_center_velocity_dt import (
    KokunoPA10CartesianCenterVelocityTimeDerivative,
)


SCHEMA = "kokuno-a4-pa10-cartesian-center-velocity-dt-independent-audit-v1"
SEED = 9173471
SAMPLE_COUNT = 1024
TIME_STEPS = (4.0e-4, 2.0e-4, 1.0e-4)
FINE_RELATIVE_RMS_GATE = 5.0e-6
FINE_RELATIVE_MAX_GATE = 3.0e-5
REFINEMENT_RATIO_GATE = 6.0
REFINEMENT_FLOOR = 2.0e-10
NONTRIVIAL_DT_RMS_FLOOR = 1.0e-8
EXACT_AXIS_TRANSVERSE_DT_GATE = 1.0e-13
FINAL_MOMENTUM_GATE = 1.0e-3
FINAL_DIVERGENCE_GATE = 1.0e-5

_TRUTH_BOUNDARY = {
    "agent4_public_velocity_only_fd4_time_derivative_used": True,
    "agent4_fresh_offgrid_fixed_cartesian_samples_used": True,
    "agent4_three_time_resolutions_used": True,
    "inner_cartesian_center_velocity_dt_independently_audited": True,
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


def _regularized_relative_rms(actual: np.ndarray, reference: np.ndarray) -> float:
    actual_arr = np.asarray(actual, dtype=float)
    reference_arr = np.asarray(reference, dtype=float)
    numerator = float(np.sqrt(np.mean(np.square(actual_arr - reference_arr))))
    denominator = 1.0 + float(np.sqrt(np.mean(np.square(reference_arr))))
    return numerator / denominator


def _regularized_relative_max(actual: np.ndarray, reference: np.ndarray) -> float:
    actual_arr = np.asarray(actual, dtype=float)
    reference_arr = np.asarray(reference, dtype=float)
    scale = 1.0 + np.abs(reference_arr)
    return float(np.max(np.abs(actual_arr - reference_arr) / scale))


def _sample_fixed_cartesian_points(
    provider: KokunoPA10CartesianCenterVelocityTimeDerivative,
    rng: np.random.Generator,
) -> dict[str, np.ndarray]:
    """Generate fresh points from source forward coordinates, not the inverse helper."""
    _, xmax = provider.source_X_interval
    X = rng.uniform(0.06 * xmax, 0.30 * xmax, size=SAMPLE_COUNT)
    eta = rng.uniform(-0.42, 0.42, size=SAMPLE_COUNT)
    t = rng.uniform(0.34, 0.66, size=SAMPLE_COUNT)
    theta = rng.uniform(-math.pi, math.pi, size=SAMPLE_COUNT)

    q = (1.0 - t) / (1.0 - eta * eta)
    z = np.power(q, provider.D) * eta
    radius = np.sqrt(2.0 * q * X)
    x = radius * np.cos(theta)
    y = radius * np.sin(theta)

    # Three exact-axis probes and two radius=1e-12 axis-near probes.  Exact
    # transverse-zero is gated only on the exact axis; the near-axis values are
    # retained separately as finite diagnostics rather than falsely required to
    # be exactly zero.
    axis_eta = np.asarray([-0.38, -1.0e-12, 0.0, 1.0e-12, 0.38], dtype=float)
    axis_t = np.asarray([0.37, 0.44, 0.50, 0.56, 0.63], dtype=float)
    axis_q = (1.0 - axis_t) / (1.0 - axis_eta * axis_eta)
    axis_z = np.power(axis_q, provider.D) * axis_eta

    near_radius = np.asarray([0.0, 1.0e-12, 0.0, 1.0e-12, 0.0], dtype=float)
    axis_x = near_radius
    axis_y = np.zeros_like(near_radius)

    return {
        "X": X,
        "eta": eta,
        "t": np.concatenate((t, axis_t)),
        "x": np.concatenate((x, axis_x)),
        "y": np.concatenate((y, axis_y)),
        "z": np.concatenate((z, axis_z)),
        "random_count": np.asarray([SAMPLE_COUNT], dtype=int),
        "axis_probe_radius": near_radius,
    }


def _fd4_velocity_dt(
    provider: KokunoPA10CartesianCenterVelocityTimeDerivative,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
    step: float,
) -> np.ndarray:
    """Differentiate only the public velocity callable at fixed Cartesian xyz."""
    fm2 = np.asarray(provider.velocity(x, y, z, t - 2.0 * step), dtype=float)
    fm1 = np.asarray(provider.velocity(x, y, z, t - step), dtype=float)
    fp1 = np.asarray(provider.velocity(x, y, z, t + step), dtype=float)
    fp2 = np.asarray(provider.velocity(x, y, z, t + 2.0 * step), dtype=float)
    return (fm2 - 8.0 * fm1 + 8.0 * fp1 - fp2) / (12.0 * step)


def _component_errors(actual: np.ndarray, reference: np.ndarray) -> dict[str, Any]:
    names = ("u_t", "v_t", "w_t")
    result: dict[str, Any] = {}
    for index, name in enumerate(names):
        result[name] = {
            "relative_rms": _regularized_relative_rms(actual[:, index], reference[:, index]),
            "relative_max": _regularized_relative_max(actual[:, index], reference[:, index]),
            "absolute_max": float(np.max(np.abs(actual[:, index] - reference[:, index]))),
        }
    return result


@dataclass(frozen=True)
class AuditResult:
    payload: dict[str, Any]

    @property
    def passed(self) -> bool:
        return bool(self.payload["passed"])


def run_independent_audit() -> AuditResult:
    provider = KokunoPA10CartesianCenterVelocityTimeDerivative()
    rng = np.random.default_rng(SEED)
    samples = _sample_fixed_cartesian_points(provider, rng)
    x = np.asarray(samples["x"], dtype=float)
    y = np.asarray(samples["y"], dtype=float)
    z = np.asarray(samples["z"], dtype=float)
    t = np.asarray(samples["t"], dtype=float)

    public_dt = np.asarray(provider.velocity_dt(x, y, z, t), dtype=float)
    if public_dt.shape != (x.size, 3) or np.any(~np.isfinite(public_dt)):
        raise RuntimeError("public velocity_dt returned malformed/nonfinite data")

    estimates: list[np.ndarray] = []
    by_step: dict[str, Any] = {}
    for step in TIME_STEPS:
        estimate = _fd4_velocity_dt(provider, x, y, z, t, step)
        estimates.append(estimate)
        by_step[f"{step:.7f}"] = {
            "relative_rms": _regularized_relative_rms(public_dt, estimate),
            "relative_max": _regularized_relative_max(public_dt, estimate),
            "component_errors": _component_errors(public_dt, estimate),
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

    fine_relative_rms = _regularized_relative_rms(public_dt, fine)
    fine_relative_max = _regularized_relative_max(public_dt, fine)
    dt_rms = float(np.sqrt(np.mean(np.square(public_dt))))
    dt_max = float(np.max(np.abs(public_dt)))

    error = np.abs(public_dt - fine)
    flat_index = int(np.argmax(error))
    point_index, component_index = np.unravel_index(flat_index, error.shape)
    worst = {
        "point_index": int(point_index),
        "component": ("u_t", "v_t", "w_t")[component_index],
        "x": float(x[point_index]),
        "y": float(y[point_index]),
        "z": float(z[point_index]),
        "t": float(t[point_index]),
        "public": float(public_dt[point_index, component_index]),
        "independent_fd4": float(fine[point_index, component_index]),
        "absolute_error": float(error[point_index, component_index]),
    }

    axis_start = SAMPLE_COUNT
    probe_radii = np.asarray(samples["axis_probe_radius"], dtype=float)
    exact_axis_mask = probe_radii == 0.0
    near_axis_mask = ~exact_axis_mask
    axis_probe_dt = public_dt[axis_start:]
    exact_axis_transverse_max = float(
        np.max(np.abs(axis_probe_dt[exact_axis_mask, :2]))
    )
    near_axis_transverse_max = float(
        np.max(np.abs(axis_probe_dt[near_axis_mask, :2]))
    )
    near_axis_finite = bool(np.all(np.isfinite(axis_probe_dt[near_axis_mask])))

    # Pre-registered mutations.  They must not be able to pass the frozen gate.
    scaled_mutation_rel_max = _regularized_relative_max(0.99 * public_dt, fine)
    sign_mutation = np.array(public_dt, copy=True)
    sign_mutation[:, 2] *= -1.0
    sign_mutation_rel_max = _regularized_relative_max(sign_mutation, fine)
    mutations_detected = bool(
        scaled_mutation_rel_max > FINE_RELATIVE_MAX_GATE
        and sign_mutation_rel_max > FINE_RELATIVE_MAX_GATE
    )

    gates = {
        "fine_relative_rms": fine_relative_rms <= FINE_RELATIVE_RMS_GATE,
        "fine_relative_max": fine_relative_max <= FINE_RELATIVE_MAX_GATE,
        "three_level_refinement_stable": refinement_stable,
        "velocity_dt_nontrivial": dt_rms >= NONTRIVIAL_DT_RMS_FLOOR,
        "exact_axis_transverse_velocity_dt": (
            exact_axis_transverse_max <= EXACT_AXIS_TRANSVERSE_DT_GATE
        ),
        "axis_near_velocity_dt_finite": near_axis_finite,
        "mutations_detected": mutations_detected,
    }
    passed = bool(all(gates.values()))

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "seed": SEED,
        "sample_count_random": SAMPLE_COUNT,
        "sample_count_total": int(x.size),
        "time_steps": list(TIME_STEPS),
        "frozen_gates": {
            "fine_relative_rms": FINE_RELATIVE_RMS_GATE,
            "fine_relative_max": FINE_RELATIVE_MAX_GATE,
            "refinement_ratio_or_floor": {
                "ratio_min": REFINEMENT_RATIO_GATE,
                "fine_difference_floor": REFINEMENT_FLOOR,
            },
            "nontrivial_velocity_dt_rms": NONTRIVIAL_DT_RMS_FLOOR,
            "exact_axis_transverse_velocity_dt_max": EXACT_AXIS_TRANSVERSE_DT_GATE,
            "final_normalized_momentum_max_l2": FINAL_MOMENTUM_GATE,
            "final_divergence_max_l2": FINAL_DIVERGENCE_GATE,
        },
        "measurements": {
            "velocity_dt_rms": dt_rms,
            "velocity_dt_max": dt_max,
            "by_step": by_step,
            "fine_relative_rms": fine_relative_rms,
            "fine_relative_max": fine_relative_max,
            "coarse_medium_difference_rms": coarse_medium_difference_rms,
            "medium_fine_difference_rms": medium_fine_difference_rms,
            "refinement_ratio": refinement_ratio,
            "exact_axis_transverse_velocity_dt_max": exact_axis_transverse_max,
            "axis_near_radius": probe_radii[near_axis_mask].tolist(),
            "axis_near_transverse_velocity_dt_max": near_axis_transverse_max,
            "worst_error": worst,
            "scaled_0p99_mutation_relative_max": scaled_mutation_rel_max,
            "axial_sign_mutation_relative_max": sign_mutation_rel_max,
        },
        "gates": gates,
        "truth_boundary": dict(_TRUTH_BOUNDARY),
        "limitations": [
            "This is an inner contraction-center velocity_dt audit, not a global leading-field admission.",
            "No pressure, restricted forcing, outer join, oscillatory/correction composite, or full NS momentum residual is evaluated.",
            "The source-complex C certificate is not independently admitted by this audit.",
            "Exact transverse zero is asserted only on the exact axis; radius=1e-12 probes are retained as finite axis-near diagnostics.",
            "No source parameter is retuned from these held-out measurements.",
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
