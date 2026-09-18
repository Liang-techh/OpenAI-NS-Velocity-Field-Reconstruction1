"""Independent black-box audit of the Kokuno reserved-patch leading velocity.

This module validates :class:`KokunoOuterReservedPatchSchedule` through a
numerical path deliberately different from its construction path.  In
particular, it

* reconstructs all outer logarithmic scales with high-precision ``Decimal``
  arithmetic from the raw public parameters;
* inverts the source similarity relation with monotone bisection rather than
  the production fixed-point solver;
* evaluates the source power-law swirl formula independently and compares it
  only with the public ``patch_velocity(x,y,z,t)`` interface;
* checks the radial scale law and local solenoidality at three finite-
  difference resolutions; and
* calibrates the divergence detector by adding a fixed two-percent radial
  mutation to the public swirl field.

The object under test is only the local source reserved patch introduced by
Kokuno Agent 1 PR #300.  It is not a global leading field and this audit does
not evaluate the formal full-domain Navier--Stokes gate.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, localcontext
import argparse
import json
import math
import os
from pathlib import Path
from typing import Any, Callable

import numpy as np

from .kokuno_outer_reserved_patch_schedule import KokunoOuterReservedPatchSchedule


TASK_ID = "KOKUNO-A4-OUTER-RESERVED-PATCH-INDEPENDENT-011"
UPSTREAM_AGENT1_PR = 300
UPSTREAM_AGENT1_HEAD = "8aab23c07e0d239624d7c72bb9c3e6577674d984"
SEED = 9173061
RELATIVE_FD_STEPS = (8.0e-4, 4.0e-4, 2.0e-4)
VELOCITY_RELATIVE_TOLERANCE = 2.0e-10
SCALE_LOG_ABSOLUTE_TOLERANCE = 2.0e-11
SCALE_LAW_RELATIVE_TOLERANCE = 2.0e-10
Q_RELATION_RELATIVE_TOLERANCE = 2.0e-12
LOCAL_NORMALIZED_DIVERGENCE_TOLERANCE = 2.0e-5
MUTATION_MIN_NORMALIZED_DIVERGENCE = 5.0e-4
MUTATION_RADIAL_FRACTION = 0.02
FORMAL_MOMENTUM_GATE = 1.0e-3
FORMAL_DIVERGENCE_GATE = 1.0e-5

_PARAMETER_CASES = (
    {
        "name": "default",
        "M_d": 3.0,
        "log_p_star_margin": 1.0,
        "C": 2.0,
        "lambda_outer": 0.05,
        "h": 0.005,
        "repair_log_offset": -17.5,
    },
    {
        "name": "lower_scale",
        "M_d": 2.5,
        "log_p_star_margin": 0.8,
        "C": 1.7,
        "lambda_outer": 0.04,
        "h": 0.004,
        "repair_log_offset": -17.8,
    },
    {
        "name": "upper_scale",
        "M_d": 3.2,
        "log_p_star_margin": 1.2,
        "C": 2.4,
        "lambda_outer": 0.06,
        "h": 0.006,
        "repair_log_offset": -17.2,
    },
)

_SCALE_KEYS = (
    "T_d",
    "log_P_star",
    "T_w",
    "log_X_R",
    "log_X_w",
    "log_P1",
    "log_e_w",
    "log_c_patch",
    "log_X_star",
    "log_e_star",
)


def _d(value: float | int | str) -> Decimal:
    return Decimal(str(value))


def _reference_scales(schedule: KokunoOuterReservedPatchSchedule) -> dict[str, float]:
    """Rebuild source scales with high-precision arithmetic from raw parameters."""

    with localcontext() as context:
        context.prec = 70
        M_d = _d(schedule.M_d)
        margin = _d(schedule.log_p_star_margin)
        C = _d(schedule.C)
        lam = _d(schedule.lambda_outer)
        offset = _d(schedule.repair_log_offset)
        Td = M_d.exp() + _d(10)
        log_p_star = Td + margin
        T_w = _d(60) * (-lam.ln())
        log_X_R = _d(110).ln() + _d(10) * (C.ln() + log_p_star)
        log_X_w = log_X_R + Td + _d(2)
        log_P1 = log_p_star - _d("0.2")
        log_e_w = log_P1 - Td / _d(2) - _d("0.5") - lam / _d(2)
        log_c_patch = log_e_w + (_d("0.5") + lam) * log_X_w
        log_X_star = log_X_w + T_w + offset
        log_e_star = log_c_patch - (_d("0.5") + lam) * log_X_star
        base = log_X_w + T_w
        return {
            "T_d": float(Td),
            "log_P_star": float(log_p_star),
            "T_w": float(T_w),
            "log_X_R": float(log_X_R),
            "log_X_w": float(log_X_w),
            "log_P1": float(log_P1),
            "log_e_w": float(log_e_w),
            "log_c_patch": float(log_c_patch),
            "log_X_star": float(log_X_star),
            "log_e_star": float(log_e_star),
            "I1_low": float(base - _d(25)),
            "I1_high": float(base - _d(20)),
            "I2_low": float(base - _d(20)),
            "I2_high": float(base - _d(15)),
            "I3_low": float(base - _d(14)),
            "I3_high": float(base - _d(9)),
            "I4_low": float(base - _d(8)),
            "I4_high": float(base - _d(3)),
        }


def _q_bisection(z: float, t: float, h: float) -> float:
    """Independently invert q-z^2 q^(2h)=1-t by monotone bisection."""

    tau = 1.0 - float(t)
    if not (tau > 0.0):
        raise ValueError("independent Kokuno inversion requires t<1")
    D = 0.5 - float(h)
    q_star = abs(float(z)) ** (1.0 / D)
    z2 = float(z) * float(z)

    def residual(q: float) -> float:
        return q - z2 * q ** (2.0 * h) - tau

    low = q_star
    high = max(q_star + tau, tau)
    if residual(low) > 0.0:
        raise RuntimeError("independent lower bracket has the wrong sign")
    grow_count = 0
    while residual(high) <= 0.0:
        high = 2.0 * high + tau
        grow_count += 1
        if grow_count > 32 or not math.isfinite(high):
            raise RuntimeError("independent q bracket failed")

    for _ in range(120):
        middle = 0.5 * (low + high)
        if residual(middle) <= 0.0:
            low = middle
        else:
            high = middle
    return 0.5 * (low + high)


def _source_reference_velocity(
    x: float,
    y: float,
    z: float,
    t: float,
    schedule: KokunoOuterReservedPatchSchedule,
    scales: dict[str, float],
) -> tuple[np.ndarray, dict[str, float]]:
    """Evaluate the public source power law without production coordinate/profile helpers."""

    q = _q_bisection(z, t, schedule.h)
    D = 0.5 - schedule.h
    eta = z / q**D
    radius = math.hypot(x, y)
    X = radius * radius / (2.0 * q)
    if X <= 0.0 or not abs(eta) < 1.0:
        raise RuntimeError("independent source point left the patch coordinate domain")
    log_x_star = math.log(X) - scales["log_X_star"]
    if not scales["I2_low"] < math.log(X) < scales["I2_high"]:
        raise RuntimeError("independent source point left I2")

    e_star = math.exp(scales["log_e_star"])
    f_eta = 1.0 / (1.0 + eta * eta)
    E = e_star * f_eta * math.exp((-0.5 - schedule.lambda_outer) * log_x_star)
    u_theta = q ** (-0.5 - schedule.h) * E
    cos_phi = x / radius
    sin_phi = y / radius
    velocity = np.asarray((-sin_phi * u_theta, cos_phi * u_theta, 0.0), dtype=float)
    relation_residual = q - z * z * q ** (2.0 * schedule.h) - (1.0 - t)
    relation_scale = max(abs(q), abs(1.0 - t), np.finfo(float).tiny)
    metadata = {
        "q": q,
        "eta": eta,
        "X": X,
        "x_star": math.exp(log_x_star),
        "u_theta": u_theta,
        "q_relation_relative_residual": abs(relation_residual) / relation_scale,
    }
    return velocity, metadata


def _sample_point(
    rng: np.random.Generator,
    schedule: KokunoOuterReservedPatchSchedule,
    scales: dict[str, float],
    time_value: float,
) -> tuple[float, float, float, float]:
    eta = float(rng.uniform(-0.55, 0.55))
    x_star = float(rng.uniform(0.92, 1.45))
    phi = float(rng.uniform(0.15, 2.0 * math.pi - 0.15))
    q = (1.0 - time_value) / (1.0 - eta * eta)
    z = q ** (0.5 - schedule.h) * eta
    X = math.exp(scales["log_X_star"] + math.log(x_star))
    radius = math.sqrt(2.0 * q * X)
    return radius * math.cos(phi), radius * math.sin(phi), z, time_value


def _as_velocity(value: Any) -> np.ndarray:
    array = np.asarray(value, dtype=float)
    if array.shape != (3,) or not np.all(np.isfinite(array)):
        raise ValueError("velocity callable must return one finite Cartesian vector")
    return array


def _normalized_divergence(
    velocity: Callable[[float, float, float, float], np.ndarray],
    point: tuple[float, float, float, float],
    h_parameter: float,
    relative_step: float,
) -> float:
    x, y, z, t = point
    radius = math.hypot(x, y)
    q = _q_bisection(z, t, h_parameter)
    axial_scale = q ** (0.5 - h_parameter)
    dx = relative_step * radius
    dz = relative_step * axial_scale
    if dx <= 0.0 or dz <= 0.0:
        raise RuntimeError("independent finite-difference scale vanished")

    vx_plus = _as_velocity(velocity(x + dx, y, z, t))
    vx_minus = _as_velocity(velocity(x - dx, y, z, t))
    vy_plus = _as_velocity(velocity(x, y + dx, z, t))
    vy_minus = _as_velocity(velocity(x, y - dx, z, t))
    vz_plus = _as_velocity(velocity(x, y, z + dz, t))
    vz_minus = _as_velocity(velocity(x, y, z - dz, t))
    divergence = (
        (vx_plus[0] - vx_minus[0]) / (2.0 * dx)
        + (vy_plus[1] - vy_minus[1]) / (2.0 * dx)
        + (vz_plus[2] - vz_minus[2]) / (2.0 * dz)
    )
    center = _as_velocity(velocity(x, y, z, t))
    local_scale = max(float(np.linalg.norm(center)) / radius, np.finfo(float).tiny)
    return abs(float(divergence)) / local_scale


def _radial_mutation(
    schedule: KokunoOuterReservedPatchSchedule,
) -> Callable[[float, float, float, float], np.ndarray]:
    def mutated(x: float, y: float, z: float, t: float) -> np.ndarray:
        base = _as_velocity(schedule.patch_velocity(x, y, z, t))
        radius = math.hypot(x, y)
        e_r = np.asarray((x / radius, y / radius, 0.0), dtype=float)
        e_theta = np.asarray((-y / radius, x / radius, 0.0), dtype=float)
        u_theta = float(np.dot(base, e_theta))
        return base + MUTATION_RADIAL_FRACTION * u_theta * e_r

    return mutated


def _scale_law_error(
    schedule: KokunoOuterReservedPatchSchedule,
    scales: dict[str, float],
) -> float:
    t = 0.5
    eta = 0.23
    phi = 0.73
    q = (1.0 - t) / (1.0 - eta * eta)
    z = q ** (0.5 - schedule.h) * eta
    speeds: list[float] = []
    x_stars = (0.95, 1.35)
    for x_star in x_stars:
        X = math.exp(scales["log_X_star"] + math.log(x_star))
        radius = math.sqrt(2.0 * q * X)
        x = radius * math.cos(phi)
        y = radius * math.sin(phi)
        speeds.append(float(np.linalg.norm(_as_velocity(schedule.patch_velocity(x, y, z, t)))))
    measured = speeds[1] / speeds[0]
    expected = (x_stars[1] / x_stars[0]) ** (-0.5 - schedule.lambda_outer)
    return abs(measured - expected) / abs(expected)


def audit_schedule(
    schedule: KokunoOuterReservedPatchSchedule,
    *,
    seed: int,
    samples: int = 12,
) -> dict[str, Any]:
    if samples < 6:
        raise ValueError("independent audit requires at least six off-grid samples")
    scales = _reference_scales(schedule)
    public_scales = schedule.log_scale_report()
    log_differences = {
        key: abs(float(public_scales[key]) - scales[key]) for key in _SCALE_KEYS
    }
    intervals = schedule.reserved_log_intervals()
    interval_differences = []
    for index, name in enumerate(("I1", "I2", "I3", "I4"), start=1):
        low, high = intervals[name]
        interval_differences.extend(
            [
                abs(low - scales[f"I{index}_low"]),
                abs(high - scales[f"I{index}_high"]),
            ]
        )
    max_log_scale_error = max([*log_differences.values(), *interval_differences])

    rng = np.random.default_rng(seed)
    times = (0.375, 0.5, 0.625)
    points = [
        _sample_point(rng, schedule, scales, times[index % len(times)])
        for index in range(samples)
    ]

    velocity_relative_errors: list[float] = []
    q_relation_errors: list[float] = []
    public_speed: list[float] = []
    x_star_values: list[float] = []
    eta_values: list[float] = []
    for point in points:
        reference, metadata = _source_reference_velocity(*point, schedule, scales)
        public = _as_velocity(schedule.patch_velocity(*point))
        reference_norm = float(np.linalg.norm(reference))
        relative_error = float(np.linalg.norm(public - reference)) / max(
            reference_norm, np.finfo(float).tiny
        )
        velocity_relative_errors.append(relative_error)
        q_relation_errors.append(metadata["q_relation_relative_residual"])
        public_speed.append(float(np.linalg.norm(public)))
        x_star_values.append(metadata["x_star"])
        eta_values.append(metadata["eta"])

    def public_velocity(x: float, y: float, z: float, t: float) -> np.ndarray:
        return _as_velocity(schedule.patch_velocity(x, y, z, t))

    divergence_by_resolution: dict[str, dict[str, float]] = {}
    for relative_step in RELATIVE_FD_STEPS:
        values = np.asarray(
            [
                _normalized_divergence(
                    public_velocity, point, schedule.h, relative_step
                )
                for point in points
            ],
            dtype=float,
        )
        divergence_by_resolution[f"{relative_step:.1e}"] = {
            "rms": float(np.sqrt(np.mean(values * values))),
            "max": float(np.max(values)),
        }

    mutation = _radial_mutation(schedule)
    mutation_values = np.asarray(
        [
            _normalized_divergence(
                mutation, point, schedule.h, RELATIVE_FD_STEPS[-1]
            )
            for point in points
        ],
        dtype=float,
    )
    scale_law_error = _scale_law_error(schedule, scales)
    finest_key = f"{RELATIVE_FD_STEPS[-1]:.1e}"
    finest_divergence = divergence_by_resolution[finest_key]["max"]

    guards = {
        "independent_log_scales": max_log_scale_error <= SCALE_LOG_ABSOLUTE_TOLERANCE,
        "public_velocity_matches_source_reference": max(velocity_relative_errors)
        <= VELOCITY_RELATIVE_TOLERANCE,
        "independent_q_relation": max(q_relation_errors) <= Q_RELATION_RELATIVE_TOLERANCE,
        "source_radial_scale_law": scale_law_error <= SCALE_LAW_RELATIVE_TOLERANCE,
        "local_solenoidality_finest": finest_divergence
        <= LOCAL_NORMALIZED_DIVERGENCE_TOLERANCE,
        "mutation_detected": float(np.min(mutation_values))
        >= MUTATION_MIN_NORMALIZED_DIVERGENCE,
        "nontrivial_public_patch": min(public_speed) > 0.0,
    }
    return {
        "parameters": {
            "M_d": schedule.M_d,
            "log_p_star_margin": schedule.log_p_star_margin,
            "C": schedule.C,
            "lambda_outer": schedule.lambda_outer,
            "h": schedule.h,
            "repair_log_offset": schedule.repair_log_offset,
        },
        "independent_scales": scales,
        "max_public_vs_independent_log_scale_abs_error": max_log_scale_error,
        "max_public_vs_independent_velocity_relative_error": max(velocity_relative_errors),
        "rms_public_vs_independent_velocity_relative_error": float(
            np.sqrt(np.mean(np.square(velocity_relative_errors)))
        ),
        "max_independent_q_relation_relative_residual": max(q_relation_errors),
        "public_speed_min": min(public_speed),
        "public_speed_max": max(public_speed),
        "sample_x_star_range": [min(x_star_values), max(x_star_values)],
        "sample_eta_range": [min(eta_values), max(eta_values)],
        "radial_scale_law_relative_error": scale_law_error,
        "normalized_divergence_by_relative_step": divergence_by_resolution,
        "mutation_finest_normalized_divergence_rms": float(
            np.sqrt(np.mean(mutation_values * mutation_values))
        ),
        "mutation_finest_normalized_divergence_min": float(np.min(mutation_values)),
        "mutation_finest_normalized_divergence_max": float(np.max(mutation_values)),
        "guards": guards,
        "passed": all(guards.values()),
    }


def run_audit(*, seed: int = SEED, samples_per_case: int = 12) -> dict[str, Any]:
    cases = []
    for index, parameters in enumerate(_PARAMETER_CASES):
        name = str(parameters["name"])
        kwargs = {key: value for key, value in parameters.items() if key != "name"}
        schedule = KokunoOuterReservedPatchSchedule(**kwargs)
        case = audit_schedule(
            schedule,
            seed=seed + 1009 * index,
            samples=samples_per_case,
        )
        case["name"] = name
        cases.append(case)

    all_passed = all(case["passed"] for case in cases)
    return {
        "task_id": TASK_ID,
        "upstream": {
            "agent1_pr": UPSTREAM_AGENT1_PR,
            "agent1_exact_head": UPSTREAM_AGENT1_HEAD,
            "construction_path": "production fixed-point q solver + float source-scale formulas",
            "validation_path": "Decimal scale reconstruction + monotone bisection + black-box Cartesian FD",
        },
        "implementation_commit": os.environ.get("GITHUB_SHA"),
        "seed": seed,
        "samples_per_case": samples_per_case,
        "relative_fd_steps": list(RELATIVE_FD_STEPS),
        "parameter_case_count": len(cases),
        "cases": cases,
        "predeclared_local_guards": {
            "velocity_relative_tolerance": VELOCITY_RELATIVE_TOLERANCE,
            "log_scale_absolute_tolerance": SCALE_LOG_ABSOLUTE_TOLERANCE,
            "scale_law_relative_tolerance": SCALE_LAW_RELATIVE_TOLERANCE,
            "q_relation_relative_tolerance": Q_RELATION_RELATIVE_TOLERANCE,
            "local_normalized_divergence_tolerance": LOCAL_NORMALIZED_DIVERGENCE_TOLERANCE,
            "mutation_min_normalized_divergence": MUTATION_MIN_NORMALIZED_DIVERGENCE,
            "mutation_radial_fraction": MUTATION_RADIAL_FRACTION,
        },
        "formal_project_gate": {
            "normalized_full_momentum_threshold": FORMAL_MOMENTUM_GATE,
            "divergence_max_threshold": FORMAL_DIVERGENCE_GATE,
            "assessed": False,
            "reason": (
                "Agent 1 PR #300 is a local I2 reserved-patch component, not a globally "
                "matched leading+oscillatory+correction velocity/pressure/forcing candidate"
            ),
        },
        "st006_cross_route_boundary": {
            "directly_comparable": False,
            "retained_momentum_sampled_max": 0.1082289305112118,
            "retained_momentum_volume_l2": 0.10758432876230622,
            "reason": "this audit measures local source-map implementation consistency, not full-domain momentum",
        },
        "structural_preflight_passed": all_passed,
        "truth_boundary": {
            "training_loss_used": False,
            "candidate_internal_derivatives_used": False,
            "free_residual_defined_forcing_used": False,
            "global_leading_profile_reconstructed": False,
            "complete_kokuno_composite_velocity": False,
            "formal_full_domain_pde_gate_assessed": False,
            "normalized_ns_residual_le_1e3_claimed": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
    }


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument("--samples-per-case", type=int, default=12)
    return parser


def main() -> int:
    args = _parser().parse_args()
    report = run_audit(seed=args.seed, samples_per_case=args.samples_per_case)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if report["structural_preflight_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
