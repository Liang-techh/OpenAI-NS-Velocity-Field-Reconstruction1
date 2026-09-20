"""Independent Agent-4 audit for the PA.10 Cartesian contraction-center velocity.

This validator treats :class:`KokunoPA10CartesianCenterVelocity` as a public
candidate-facing object.  It does not call its coordinate inverse, forward
coordinate helper, or report as a numerical oracle.

The independent path uses a safeguarded Newton solve for the public implicit
q-coordinate and then reconstructs the displayed Cartesian velocity formula
from the upstream public physical profile values.  Divergence is checked only
through the public ``velocity(x,y,z,t)`` callable with a separate centered FD4
operator on three frozen spatial step sizes.

Scope is intentionally narrow: the audited object is still only the PA.10
inner contraction center.  No outer/global join, fixed-point correction,
matched pressure, restricted forcing, complete Kokuno composite, or full NS
residual exists here, so none is inferred by this audit.
"""

from __future__ import annotations

import argparse
import copy
from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_pa10_cartesian_center_velocity import (
    KokunoPA10CartesianCenterVelocity,
)


SCHEMA = "kokuno-a4-pa10-cartesian-center-independent-audit-v1"
SEED = 9173461
MAPPING_SAMPLE_COUNT = 2048
DIVERGENCE_SAMPLE_COUNT = 48
FD_STEPS = (0.004, 0.002, 0.001)
VELOCITY_RELATIVE_MAX_GATE = 5.0e-10
Q_RELATIVE_MAX_GATE = 5.0e-12
DIVERGENCE_MAX_GATE = 1.0e-5
DIVERGENCE_L2_GATE = 1.0e-5
NONTRIVIAL_RMS_FLOOR = 1.0e-8

_TRUTH_BOUNDARY = {
    "agent4_independent_coordinate_inverse_used": True,
    "agent4_independent_cartesian_formula_reconstruction_used": True,
    "agent4_public_velocity_only_fd4_divergence_used": True,
    "inner_cartesian_center_velocity_independently_audited": True,
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


def _relative_max(actual: np.ndarray, expected: np.ndarray) -> float:
    actual_arr = np.asarray(actual, dtype=float)
    expected_arr = np.asarray(expected, dtype=float)
    scale = 1.0 + np.abs(expected_arr)
    return float(np.max(np.abs(actual_arr - expected_arr) / scale))


def _relative_rms(actual: np.ndarray, expected: np.ndarray) -> float:
    actual_arr = np.asarray(actual, dtype=float)
    expected_arr = np.asarray(expected, dtype=float)
    num = float(np.sqrt(np.mean(np.square(actual_arr - expected_arr))))
    den = 1.0 + float(np.sqrt(np.mean(np.square(expected_arr))))
    return num / den


def _independent_q_newton(
    field: KokunoPA10CartesianCenterVelocity,
    z: np.ndarray,
    t: np.ndarray,
    *,
    h_override: float | None = None,
) -> np.ndarray:
    """Safeguarded Newton solve independent of Agent-1's fixed bisection path."""
    z_ld = np.asarray(z, dtype=np.longdouble)
    t_ld = np.asarray(t, dtype=np.longdouble)
    h = np.longdouble(field.h if h_override is None else h_override)
    D = np.longdouble(0.5) - h
    tau = np.longdouble(1.0) - t_ld
    derivative_floor = np.longdouble(1.0) - np.longdouble(2.0) * h
    if derivative_floor <= 0:
        raise RuntimeError("independent q derivative floor is nonpositive")

    q_star = np.power(np.abs(z_ld), np.longdouble(1.0) / D)
    low = q_star.copy()
    high = q_star + tau / derivative_floor
    q = low + np.longdouble(0.5) * (high - low)

    def g(value: np.ndarray) -> np.ndarray:
        return value - z_ld * z_ld * np.power(value, np.longdouble(2.0) * h) - tau

    for _ in range(28):
        gv = g(q)
        derivative_term = np.zeros_like(q)
        nonzero = np.abs(z_ld) > 0
        derivative_term[nonzero] = (
            z_ld[nonzero]
            * z_ld[nonzero]
            * np.power(q[nonzero], np.longdouble(2.0) * h - np.longdouble(1.0))
        )
        dg = np.longdouble(1.0) - np.longdouble(2.0) * h * derivative_term
        newton = q - gv / dg
        midpoint = low + np.longdouble(0.5) * (high - low)
        valid = np.isfinite(newton) & (newton > low) & (newton < high)
        trial = np.where(valid, newton, midpoint)
        g_trial = g(trial)
        move_low = g_trial <= 0
        low = np.where(move_low, trial, low)
        high = np.where(move_low, high, trial)
        q = trial

    # Finish with a few bracket contractions so the answer does not depend on
    # the last Newton step when binary64 is already the dominant error.
    for _ in range(12):
        midpoint = low + np.longdouble(0.5) * (high - low)
        move_low = g(midpoint) <= 0
        low = np.where(move_low, midpoint, low)
        high = np.where(move_low, high, midpoint)
    return np.asarray(low + np.longdouble(0.5) * (high - low), dtype=float)


def _sample_similarity_points(
    field: KokunoPA10CartesianCenterVelocity,
    rng: np.random.Generator,
    count: int,
) -> dict[str, np.ndarray]:
    _, xmax = field.source_X_interval
    X = rng.uniform(0.02 * xmax, 0.62 * xmax, size=count)
    eta = rng.uniform(-0.68, 0.68, size=count)
    t = rng.uniform(0.29, 0.71, size=count)
    theta = rng.uniform(-math.pi, math.pi, size=count)
    tau = 1.0 - t
    q = tau / (1.0 - eta * eta)
    z = np.power(q, field.D) * eta
    r = np.sqrt(2.0 * q * X)
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return {
        "X": X,
        "eta": eta,
        "t": t,
        "theta": theta,
        "tau": tau,
        "q": q,
        "x": x,
        "y": y,
        "z": z,
    }


def _independent_cartesian_velocity(
    field: KokunoPA10CartesianCenterVelocity,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
    *,
    h_override: float | None = None,
) -> dict[str, np.ndarray]:
    q = _independent_q_newton(field, z, t, h_override=h_override)
    h = field.h if h_override is None else float(h_override)
    D = 0.5 - h
    A = 0.5 + h
    eta = z / np.power(q, D)
    X = (x * x + y * y) / (2.0 * q)
    profiles = field.source_center.values(X, eta)
    radial = profiles["v_0"] / (2.0 * q)
    swirl = np.power(q, -A - 0.5) * profiles["F_0"]
    u = radial * x - swirl * y
    v = radial * y + swirl * x
    w = np.power(q, -A) * profiles["U_0"]
    return {
        "q": q,
        "X": X,
        "eta": eta,
        "velocity": np.stack((u, v, w), axis=-1),
    }


def _fd4_derivative_component(
    field: KokunoPA10CartesianCenterVelocity,
    xyz: np.ndarray,
    t: np.ndarray,
    *,
    axis: int,
    component: int,
    step: float,
) -> np.ndarray:
    values: list[np.ndarray] = []
    for offset in (-2.0, -1.0, 1.0, 2.0):
        moved = np.array(xyz, copy=True)
        moved[:, axis] += offset * step
        velocity = field.velocity(moved[:, 0], moved[:, 1], moved[:, 2], t)
        values.append(np.asarray(velocity[:, component], dtype=float))
    fm2, fm1, fp1, fp2 = values
    return (fm2 - 8.0 * fm1 + 8.0 * fp1 - fp2) / (12.0 * step)


def _fd4_divergence(
    field: KokunoPA10CartesianCenterVelocity,
    xyz: np.ndarray,
    t: np.ndarray,
    step: float,
) -> np.ndarray:
    return sum(
        _fd4_derivative_component(
            field,
            xyz,
            t,
            axis=axis,
            component=axis,
            step=step,
        )
        for axis in range(3)
    )


def _divergence_cloud(
    field: KokunoPA10CartesianCenterVelocity,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    _, xmax = field.source_X_interval
    count = DIVERGENCE_SAMPLE_COUNT
    X = rng.uniform(0.08 * xmax, 0.34 * xmax, size=count)
    eta = rng.uniform(-0.48, 0.48, size=count)
    t = rng.uniform(0.32, 0.68, size=count)
    theta = rng.uniform(-math.pi, math.pi, size=count)
    tau = 1.0 - t
    q = tau / (1.0 - eta * eta)
    z = np.power(q, field.D) * eta
    r = np.sqrt(2.0 * q * X)
    xyz = np.column_stack((r * np.cos(theta), r * np.sin(theta), z))

    # Axis and axis-near samples are appended explicitly.  These stay safely
    # inside the source X interval under the frozen FD4 stencil.
    axis_eta = np.asarray([-0.42, -1e-12, 0.0, 1e-12, 0.42], dtype=float)
    axis_t = np.asarray([0.34, 0.43, 0.50, 0.57, 0.66], dtype=float)
    axis_q = (1.0 - axis_t) / (1.0 - axis_eta * axis_eta)
    axis_z = np.power(axis_q, field.D) * axis_eta
    axis_xyz = np.column_stack(
        (np.zeros(axis_eta.size), np.zeros(axis_eta.size), axis_z)
    )
    return np.vstack((xyz, axis_xyz)), np.concatenate((t, axis_t))


@dataclass(frozen=True)
class AuditResult:
    payload: dict[str, Any]

    @property
    def passed(self) -> bool:
        return bool(self.payload["passed"])


def run_independent_audit() -> AuditResult:
    field = KokunoPA10CartesianCenterVelocity()
    rng = np.random.default_rng(SEED)
    samples = _sample_similarity_points(field, rng, MAPPING_SAMPLE_COUNT)

    public_velocity = np.asarray(
        field.velocity(samples["x"], samples["y"], samples["z"], samples["t"]),
        dtype=float,
    )
    independent = _independent_cartesian_velocity(
        field,
        samples["x"],
        samples["y"],
        samples["z"],
        samples["t"],
    )

    q_rel = _relative_max(independent["q"], samples["q"])
    X_rel = _relative_max(independent["X"], samples["X"])
    eta_rel = _relative_max(independent["eta"], samples["eta"])
    velocity_rel_max = _relative_max(public_velocity, independent["velocity"])
    velocity_rel_rms = _relative_rms(public_velocity, independent["velocity"])

    velocity_rms = float(np.sqrt(np.mean(np.square(public_velocity))))
    velocity_max = float(np.max(np.abs(public_velocity)))

    # Explicit axis regularity/nontriviality probes through the public callable.
    axis_eta = np.asarray([-0.4, 0.0, 0.4], dtype=float)
    axis_t = np.asarray([0.35, 0.50, 0.65], dtype=float)
    axis_q = (1.0 - axis_t) / (1.0 - axis_eta * axis_eta)
    axis_z = np.power(axis_q, field.D) * axis_eta
    axis_velocity = np.asarray(
        field.velocity(np.zeros(3), np.zeros(3), axis_z, axis_t), dtype=float
    )
    axis_transverse_max = float(np.max(np.abs(axis_velocity[:, :2])))

    xyz_div, t_div = _divergence_cloud(field, rng)
    divergence_by_step: dict[str, dict[str, float]] = {}
    divergence_arrays: list[np.ndarray] = []
    for step in FD_STEPS:
        div = _fd4_divergence(field, xyz_div, t_div, step)
        divergence_arrays.append(div)
        divergence_by_step[f"{step:.6f}"] = {
            "max_abs": float(np.max(np.abs(div))),
            "l2": float(np.sqrt(np.mean(np.square(div)))),
        }

    fine = divergence_arrays[-1]
    fine_max = float(np.max(np.abs(fine)))
    fine_l2 = float(np.sqrt(np.mean(np.square(fine))))
    coarse_max = divergence_by_step[f"{FD_STEPS[0]:.6f}"]["max_abs"]
    medium_max = divergence_by_step[f"{FD_STEPS[1]:.6f}"]["max_abs"]
    refinement_stable = bool(
        medium_max <= max(1.25 * coarse_max, 1.0e-8)
        and fine_max <= max(1.25 * medium_max, 1.0e-8)
    )

    # Pre-registered mutation controls.  They are evaluated after the real
    # measurements but use no adaptive threshold or parameter selection.
    scaled_public_rel = _relative_max(0.999 * public_velocity, independent["velocity"])
    q_mut = 0.999 * independent["q"]
    eta_mut = samples["z"] / np.power(q_mut, field.D)
    X_mut = (samples["x"] ** 2 + samples["y"] ** 2) / (2.0 * q_mut)
    q_coordinate_mutation_error = max(
        _relative_max(q_mut, samples["q"]),
        _relative_max(eta_mut, samples["eta"]),
        _relative_max(X_mut, samples["X"]),
    )

    # Diagnostic-only sensitivity to the source h parameter; no retuning uses it.
    plus = _independent_cartesian_velocity(
        field,
        samples["x"],
        samples["y"],
        samples["z"],
        samples["t"],
        h_override=field.h * 1.001,
    )["velocity"]
    minus = _independent_cartesian_velocity(
        field,
        samples["x"],
        samples["y"],
        samples["z"],
        samples["t"],
        h_override=field.h * 0.999,
    )["velocity"]
    h_perturbation_relative_rms = _relative_rms(plus, minus)

    guards = {
        "independent_q_relative_max": q_rel <= Q_RELATIVE_MAX_GATE,
        "independent_X_relative_max": X_rel <= Q_RELATIVE_MAX_GATE,
        "independent_eta_relative_max": eta_rel <= Q_RELATIVE_MAX_GATE,
        "public_velocity_relative_max": velocity_rel_max <= VELOCITY_RELATIVE_MAX_GATE,
        "public_velocity_nontrivial": velocity_rms >= NONTRIVIAL_RMS_FLOOR,
        "axis_transverse_exact_zero": axis_transverse_max == 0.0,
        "fd4_divergence_finest_max": fine_max <= DIVERGENCE_MAX_GATE,
        "fd4_divergence_finest_l2": fine_l2 <= DIVERGENCE_L2_GATE,
        "fd4_three_level_resolution_stable": refinement_stable,
        "scaled_velocity_mutation_detected": scaled_public_rel > VELOCITY_RELATIVE_MAX_GATE,
        "q_coordinate_mutation_detected": q_coordinate_mutation_error > Q_RELATIVE_MAX_GATE,
    }
    failed = sorted(name for name, passed in guards.items() if not passed)

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "candidate": {
            "module": "openai_ns_reconstruction.kokuno_pa10_cartesian_center_velocity",
            "class": "KokunoPA10CartesianCenterVelocity",
            "field_sha256": field.field_sha256,
            "scope": "PA.10 inner contraction-center Cartesian velocity only",
        },
        "protocol": {
            "seed": SEED,
            "mapping_sample_count": MAPPING_SAMPLE_COUNT,
            "divergence_random_sample_count": DIVERGENCE_SAMPLE_COUNT,
            "divergence_axis_probe_count": 5,
            "fd4_spatial_steps": list(FD_STEPS),
            "independent_q_solver": "longdouble safeguarded Newton plus bracket contraction",
            "candidate_q_solver": "fixed source-monotone bisection",
            "velocity_relative_max_gate": VELOCITY_RELATIVE_MAX_GATE,
            "q_relative_max_gate": Q_RELATIVE_MAX_GATE,
            "divergence_max_gate": DIVERGENCE_MAX_GATE,
            "divergence_l2_gate": DIVERGENCE_L2_GATE,
            "final_project_momentum_gate": 1.0e-3,
            "final_project_divergence_gate": 1.0e-5,
            "free_residual_defined_forcing_forbidden": True,
        },
        "measurements": {
            "q_relative_max": q_rel,
            "X_relative_max": X_rel,
            "eta_relative_max": eta_rel,
            "velocity_relative_max": velocity_rel_max,
            "velocity_relative_rms": velocity_rel_rms,
            "velocity_rms": velocity_rms,
            "velocity_max_abs": velocity_max,
            "axis_transverse_max_abs": axis_transverse_max,
            "divergence_by_step": divergence_by_step,
            "h_plus_minus_0p1pct_velocity_relative_rms": h_perturbation_relative_rms,
        },
        "mutation_controls": {
            "public_velocity_times_0p999_relative_max": scaled_public_rel,
            "independent_q_times_0p999_coordinate_error": q_coordinate_mutation_error,
        },
        "guards": guards,
        "failed_guards": failed,
        "passed": not failed,
        "truth_boundary": copy.deepcopy(_TRUTH_BOUNDARY),
    }
    payload["receipt_sha256"] = hashlib.sha256(
        _canonical_json(payload).encode("utf-8")
    ).hexdigest()
    return AuditResult(payload)


def save_independent_audit(path: str | Path) -> dict[str, Any]:
    result = run_independent_audit()
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(result.payload, indent=2, sort_keys=True) + "\n")
    return result.payload


def _main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = save_independent_audit(args.output)
    print("passed=", payload["passed"])
    print("failed_guards=", payload["failed_guards"])
    print("measurements=", payload["measurements"])
    print("receipt_sha256=", payload["receipt_sha256"])
    raise SystemExit(0 if payload["passed"] else 1)


if __name__ == "__main__":
    _main()
