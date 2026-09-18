"""Independent value-only audit of Agent-2 autonomous signed rectangle geometry.

This audit is deliberately downstream of ``KokunoAutonomousSignedRectangleGeometry``.
It consumes only the public ``receipt()`` and ``instantiate(partition)`` values and
rebuilds the mathematical checks on an independent code path.  In particular it
never calls the Agent-2 rational-witness object, its torus-distance helpers, its
band-covering object, or its beta-center assignment helper.

The source does not publish the actual hidden rectangle centers, radius, partition
labels, or pulse integrals.  Therefore a PASS here certifies only the repository-
autonomous/source-compatible geometry handoff.  It is not a Navier--Stokes
residual, not actual-source recovery, and not permission to set ``pde_validated``.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_autonomous_signed_rectangle_geometry import (
    KokunoAutonomousSignedRectangleGeometry,
)

SEED = 9173211
J_G = ((3, 1), (1, 5))
T_G = 4.0 + math.sqrt(2.0)
B_G = math.sqrt(2.0) - 1.0

CASES = (
    {"h": 0.0025, "ell0": 12, "color_count": 5, "denominator": 19, "safety": 0.37},
    {"h": 0.0065, "ell0": 19, "color_count": 6, "denominator": 23, "safety": 0.41},
    {"h": 0.0090, "ell0": 37, "color_count": 5, "denominator": 29, "safety": 0.47},
)
ELL_OFFSETS = (0, 1, 3, 5)
PERMUTATION = (2, 0, 3, 1)

FROZEN_GUARDS = {
    "schedule_relative_max": 5.0e-13,
    "center_float_relative_max": 5.0e-15,
    "permutation_relative_max": 5.0e-15,
    "order_dependent_mutation_mismatch_fraction_min": 0.25,
    "wrong_pulse_formula_relative_rms_min": 0.25,
}


def _canonical_beta_bytes(label: tuple[int, tuple[int, int, int]]) -> bytes:
    # Independent implementation of the documented semantic-label convention.
    payload = [int(label[0]), [int(x) for x in label[1]]]
    return json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")


def _expected_center_indices(
    label: tuple[int, tuple[int, int, int]], color_count: int
) -> tuple[int, int]:
    plus = int.from_bytes(hashlib.sha256(_canonical_beta_bytes(label)).digest()[:8], "big") % color_count
    return plus, (plus + 1) % color_count


def _matmul_int(
    left: tuple[tuple[int, int], tuple[int, int]],
    right: tuple[tuple[int, int], tuple[int, int]],
) -> tuple[tuple[int, int], tuple[int, int]]:
    return (
        (
            left[0][0] * right[0][0] + left[0][1] * right[1][0],
            left[0][0] * right[0][1] + left[0][1] * right[1][1],
        ),
        (
            left[1][0] * right[0][0] + left[1][1] * right[1][0],
            left[1][0] * right[0][1] + left[1][1] * right[1][1],
        ),
    )


def _matrix_power_binary(power: int) -> tuple[tuple[int, int], tuple[int, int]]:
    if power < 0:
        raise ValueError("power must be nonnegative")
    result = ((1, 0), (0, 1))
    base = J_G
    exponent = int(power)
    while exponent:
        if exponent & 1:
            result = _matmul_int(result, base)
        base = _matmul_int(base, base)
        exponent >>= 1
    return result


def _mod1(value: Fraction) -> Fraction:
    return value - math.floor(value)


def _transform_center(
    center: tuple[Fraction, Fraction], power: int
) -> tuple[Fraction, Fraction]:
    matrix = _matrix_power_binary(power)
    x, y = center
    return (
        _mod1(matrix[0][0] * x + matrix[0][1] * y),
        _mod1(matrix[1][0] * x + matrix[1][1] * y),
    )


def _wrapped_abs(value: Fraction) -> Fraction:
    residue = _mod1(value)
    return min(residue, 1 - residue)


def _torus_distance_squared(
    left: tuple[Fraction, Fraction], right: tuple[Fraction, Fraction]
) -> Fraction:
    dx = _wrapped_abs(left[0] - right[0])
    dy = _wrapped_abs(left[1] - right[1])
    return dx * dx + dy * dy


def _parse_center(payload: list[dict[str, int]]) -> tuple[Fraction, Fraction]:
    if len(payload) != 2:
        raise ValueError("center payload must contain two coordinates")
    return tuple(
        Fraction(int(row["numerator"]), int(row["denominator"])) for row in payload
    )  # type: ignore[return-value]


def _minimum_separation_exact(
    centers: tuple[tuple[Fraction, Fraction], ...], delta_max: int
) -> Fraction:
    minimum: Fraction | None = None
    for mu, c_mu in enumerate(centers):
        for nu, c_nu in enumerate(centers):
            for delta in range(delta_max + 1):
                if delta == 0 and mu == nu:
                    continue
                value = _torus_distance_squared(c_mu, _transform_center(c_nu, delta))
                if value == 0:
                    return Fraction(0, 1)
                if minimum is None or value < minimum:
                    minimum = value
    if minimum is None:
        raise ValueError("at least two centers are required")
    return minimum


def _source_delta_max_independent(h: float, ell0: int) -> int:
    bound = 1.0 + (
        4.0 * (1.0 + h) * math.log(2.0) + 8.0 / float(ell0)
    ) / math.log(T_G)
    return int(math.ceil(bound))


def _band_schedule_independent(ell: int, h: float, r0: float) -> dict[str, float | int]:
    log_q = -float(ell) * math.log(2.0)
    q = math.exp(log_q)
    epsilon = math.exp(h * log_q)
    s_star = float(ell * ell)
    level_real = (-(1.0 + h) * log_q - math.log(s_star)) / math.log(T_G)
    level = math.floor(level_real)
    c_i = math.exp(level * math.log(T_G) + (1.0 + h) * log_q)
    return {
        "Q": q,
        "epsilon": epsilon,
        "covering_level": int(level),
        "c_i": c_i,
        "L_s": 2.0 * r0 / c_i,
    }


def _partition(case_index: int, ell0: int) -> dict[str, Any]:
    # Four-component unit-sphere partition with nonzero analytic r/z derivatives.
    rng = np.random.default_rng(SEED + 101 * case_index)
    points = rng.uniform(low=(-0.7, -0.6), high=(0.8, 0.7), size=(9, 2))
    eta_rows: list[np.ndarray] = []
    dr_rows: list[np.ndarray] = []
    dz_rows: list[np.ndarray] = []

    a_r, a_z = 0.17, 0.11
    b_r, b_z = -0.13, 0.09
    c_r, c_z = 0.07, -0.15
    for r, z in points:
        a = 0.62 + a_r * r + a_z * z
        b = 0.88 + b_r * r + b_z * z
        c = 1.14 + c_r * r + c_z * z
        sa, ca = math.sin(a), math.cos(a)
        sb, cb = math.sin(b), math.cos(b)
        sc, cc = math.sin(c), math.cos(c)

        eta = np.array((ca, sa * cb, sa * sb * cc, sa * sb * sc), dtype=float)
        d_da = np.array((-sa, ca * cb, ca * sb * cc, ca * sb * sc), dtype=float)
        d_db = np.array((0.0, -sa * sb, sa * cb * cc, sa * cb * sc), dtype=float)
        d_dc = np.array((0.0, 0.0, -sa * sb * sc, sa * sb * cc), dtype=float)
        eta_rows.append(eta)
        dr_rows.append(a_r * d_da + b_r * d_db + c_r * d_dc)
        dz_rows.append(a_z * d_da + b_z * d_db + c_z * d_dc)

    labels = tuple(
        (
            ell0 + ELL_OFFSETS[j],
            (2 * case_index + j, -3 * case_index + 2 * j, 5 - case_index - j),
        )
        for j in range(4)
    )
    return {
        "eta": np.stack(eta_rows, axis=0),
        "D_r_eta": np.stack(dr_rows, axis=0),
        "D_z_eta": np.stack(dz_rows, axis=0),
        "beta_labels": labels,
    }


def _relative_scalar(left: float, right: float) -> float:
    return abs(left - right) / max(1.0e-300, abs(right))


def _relative_rms(left: np.ndarray, right: np.ndarray) -> float:
    numerator = float(np.sqrt(np.mean((left - right) ** 2)))
    denominator = max(1.0e-300, float(np.sqrt(np.mean(right**2))))
    return numerator / denominator


def _same_exact_center(left: Any, right: Any) -> bool:
    return _parse_center(left) == _parse_center(right)


def _permutation_result(
    geometry: KokunoAutonomousSignedRectangleGeometry,
    partition: dict[str, Any],
) -> dict[str, Any]:
    perm = np.asarray(PERMUTATION, dtype=int)
    return geometry.instantiate(
        {
            "eta": np.asarray(partition["eta"], dtype=float)[..., perm],
            "D_r_eta": np.asarray(partition["D_r_eta"], dtype=float)[..., perm],
            "D_z_eta": np.asarray(partition["D_z_eta"], dtype=float)[..., perm],
            "beta_labels": tuple(partition["beta_labels"][int(i)] for i in perm),
        }
    )


def run_audit() -> dict[str, Any]:
    schedule_errors: list[float] = []
    center_float_errors: list[float] = []
    permutation_errors: list[float] = []
    wrong_pulse_errors: list[float] = []
    order_dependent_mutation_mismatches = 0
    order_dependent_mutation_total = 0
    covering_level_mismatches = 0
    center_index_mismatches = 0
    center_exact_mismatches = 0
    delta_max_mismatches = 0
    dc_exact_mismatches = 0
    strict_guard_violations = 0
    permutation_exact_mismatches = 0
    swapped_sign_assignment_mismatches = 0
    swapped_sign_assignment_total = 0
    forged_partition_rejections = 0
    case_receipts: list[dict[str, Any]] = []

    for case_index, cfg in enumerate(CASES):
        geometry = KokunoAutonomousSignedRectangleGeometry(**cfg)
        partition = _partition(case_index, int(cfg["ell0"]))
        public = geometry.instantiate(partition)
        receipt = geometry.receipt()

        expected_delta_max = _source_delta_max_independent(float(cfg["h"]), int(cfg["ell0"]))
        actual_delta_max = int(receipt["witness"]["delta_max"])
        delta_max_mismatches += int(actual_delta_max != expected_delta_max)

        centers = tuple(_parse_center(payload) for payload in receipt["witness"]["centers"])
        independent_dc2 = _minimum_separation_exact(centers, actual_delta_max)
        receipt_dc2 = Fraction(
            int(receipt["witness"]["d_c_squared"]["numerator"]),
            int(receipt["witness"]["d_c_squared"]["denominator"]),
        )
        dc_exact_mismatches += int(independent_dc2 != receipt_dc2)

        c_v = 2.0 * math.sqrt(1.0 + B_G * B_G)
        j_power = np.asarray(_matrix_power_binary(actual_delta_max), dtype=float)
        c_j = float(np.linalg.norm(j_power, ord=2))
        r0 = float(public["r0"])
        injectivity_margin = 1.0 - 4.0 * r0 * c_v
        separation_margin = math.sqrt(float(independent_dc2)) - 2.0 * r0 * c_v * (c_j + 1.0)
        strict_guard_violations += int(injectivity_margin <= 0.0 or separation_margin <= 0.0)

        labels = tuple(public["beta_labels"])
        for j, label in enumerate(labels):
            expected_pair = _expected_center_indices(label, int(cfg["color_count"]))
            actual_pair = tuple(int(x) for x in public["center_index_by_beta_sign"][j])
            center_index_mismatches += int(actual_pair != expected_pair)

            for sidx, index in enumerate(expected_pair):
                expected_center = centers[index]
                actual_exact = _parse_center(public["center_exact_by_beta_sign"][j][sidx])
                center_exact_mismatches += int(actual_exact != expected_center)
                actual_float = np.asarray(public["center_float_by_beta_sign"][j, sidx], dtype=float)
                expected_float = np.asarray(tuple(float(x) for x in expected_center), dtype=float)
                center_float_errors.append(
                    float(np.linalg.norm(actual_float - expected_float) / max(1.0, np.linalg.norm(expected_float)))
                )

            schedule = _band_schedule_independent(int(label[0]), float(cfg["h"]), r0)
            for key, public_key in (("Q", "Q_by_beta"), ("epsilon", "epsilon_by_beta"), ("L_s", "L_s_by_beta")):
                schedule_errors.append(
                    _relative_scalar(float(public[public_key][j]), float(schedule[key]))
                )
            covering_level_mismatches += int(
                int(public["covering_level_by_beta"][j]) != int(schedule["covering_level"])
            )

            wrong_l_s = 2.0 * r0 * float(schedule["c_i"])
            wrong_pulse_errors.append(
                _relative_scalar(wrong_l_s, float(schedule["L_s"]))
            )

            # Mutation 1: replace semantic SHA assignment by caller-axis order.
            wrong_pair = (j % int(cfg["color_count"]), (j + 1) % int(cfg["color_count"]))
            order_dependent_mutation_total += 1
            order_dependent_mutation_mismatches += int(wrong_pair != expected_pair)

            # Mutation 2: swap sigma+ and sigma- while keeping the sign labels fixed.
            swapped = (actual_pair[1], actual_pair[0])
            swapped_sign_assignment_total += 2
            swapped_sign_assignment_mismatches += int(swapped[0] != expected_pair[0])
            swapped_sign_assignment_mismatches += int(swapped[1] != expected_pair[1])

        permuted = _permutation_result(geometry, partition)
        original_lookup = {label: j for j, label in enumerate(public["beta_labels"])}
        permuted_lookup = {label: j for j, label in enumerate(permuted["beta_labels"])}
        for label in original_lookup:
            j0 = original_lookup[label]
            j1 = permuted_lookup[label]
            if tuple(public["center_index_by_beta_sign"][j0]) != tuple(permuted["center_index_by_beta_sign"][j1]):
                permutation_exact_mismatches += 1
            for sidx in range(2):
                if not _same_exact_center(
                    public["center_exact_by_beta_sign"][j0][sidx],
                    permuted["center_exact_by_beta_sign"][j1][sidx],
                ):
                    permutation_exact_mismatches += 1
            for key in ("Q_by_beta", "epsilon_by_beta", "L_s_by_beta"):
                permutation_errors.append(
                    _relative_scalar(float(public[key][j0]), float(permuted[key][j1]))
                )
            permutation_exact_mismatches += int(
                int(public["covering_level_by_beta"][j0]) != int(permuted["covering_level_by_beta"][j1])
            )

        forged = dict(partition)
        forged["eta"] = np.asarray(partition["eta"], dtype=float) * 1.001
        try:
            geometry.instantiate(forged)
        except ValueError:
            forged_partition_rejections += 1

        case_receipts.append(
            {
                "case_index": case_index,
                "parameters": dict(cfg),
                "delta_max": actual_delta_max,
                "independent_d_c_squared": {
                    "numerator": independent_dc2.numerator,
                    "denominator": independent_dc2.denominator,
                },
                "r0": r0,
                "injectivity_margin": injectivity_margin,
                "separation_margin": separation_margin,
                "beta_labels": [
                    [int(label[0]), [int(x) for x in label[1]]] for label in labels
                ],
            }
        )

    schedule_relative_max = max(schedule_errors, default=math.inf)
    center_float_relative_max = max(center_float_errors, default=math.inf)
    permutation_relative_max = max(permutation_errors, default=math.inf)
    wrong_pulse_formula_relative_rms = float(
        np.sqrt(np.mean(np.asarray(wrong_pulse_errors, dtype=float) ** 2))
    )
    order_dependent_mutation_mismatch_fraction = (
        order_dependent_mutation_mismatches / order_dependent_mutation_total
    )
    swapped_sign_mutation_mismatch_fraction = (
        swapped_sign_assignment_mismatches / swapped_sign_assignment_total
    )

    truth_flags_preserved = True
    for cfg in CASES:
        receipt = KokunoAutonomousSignedRectangleGeometry(**cfg).receipt()
        truth = receipt["truth_boundary"]
        truth_flags_preserved = truth_flags_preserved and (
            truth["source_actual_rectangle_labels_instantiated"] is False
            and truth["source_actual_partition_labels_instantiated"] is False
            and truth["actual_source_h_sigma_pulse_integrals_bound"] is False
            and truth["actual_positive_order_background_bound"] is False
            and truth["actual_auxiliary_torus_mode_family_bound"] is False
            and truth["formal_full_domain_pde_gate_assessed"] is False
            and truth["pde_validated"] is False
        )

    local_pass = (
        schedule_relative_max <= FROZEN_GUARDS["schedule_relative_max"]
        and center_float_relative_max <= FROZEN_GUARDS["center_float_relative_max"]
        and permutation_relative_max <= FROZEN_GUARDS["permutation_relative_max"]
        and covering_level_mismatches == 0
        and center_index_mismatches == 0
        and center_exact_mismatches == 0
        and delta_max_mismatches == 0
        and dc_exact_mismatches == 0
        and strict_guard_violations == 0
        and permutation_exact_mismatches == 0
        and forged_partition_rejections == len(CASES)
        and order_dependent_mutation_mismatch_fraction
        >= FROZEN_GUARDS["order_dependent_mutation_mismatch_fraction_min"]
        and wrong_pulse_formula_relative_rms
        >= FROZEN_GUARDS["wrong_pulse_formula_relative_rms_min"]
        and swapped_sign_mutation_mismatch_fraction == 1.0
        and truth_flags_preserved
    )

    return {
        "schema": "kokuno-agent4-autonomous-signed-rectangle-independent-audit-v1",
        "seed": SEED,
        "oracle": {
            "public_values_only": True,
            "agent2_rational_witness_helper_reused": False,
            "agent2_torus_distance_helper_reused": False,
            "agent2_band_covering_helper_reused": False,
            "agent2_center_assignment_helper_reused": False,
            "matrix_power_method": "binary integer exponentiation",
            "exact_torus_arithmetic": "fractions.Fraction",
            "partition": "independent four-component analytic unit-sphere partition",
        },
        "frozen_guards": dict(FROZEN_GUARDS),
        "metrics": {
            "schedule_relative_max": schedule_relative_max,
            "center_float_relative_max": center_float_relative_max,
            "covering_level_mismatches": covering_level_mismatches,
            "center_index_mismatches": center_index_mismatches,
            "center_exact_mismatches": center_exact_mismatches,
            "delta_max_mismatches": delta_max_mismatches,
            "d_c_squared_exact_mismatches": dc_exact_mismatches,
            "strict_guard_violations": strict_guard_violations,
            "permutation_relative_max": permutation_relative_max,
            "permutation_exact_mismatches": permutation_exact_mismatches,
            "forged_partition_rejections": forged_partition_rejections,
            "forged_partition_total": len(CASES),
            "order_dependent_mutation_mismatch_fraction": order_dependent_mutation_mismatch_fraction,
            "swapped_sign_mutation_mismatch_fraction": swapped_sign_mutation_mismatch_fraction,
            "wrong_pulse_formula_relative_rms": wrong_pulse_formula_relative_rms,
        },
        "cases": case_receipts,
        "local_structural_preflight_passed": bool(local_pass),
        "physical_coordinate_probe_applicable": False,
        "formal_full_domain_pde_gate_assessed": False,
        "heldout_ns_residual_assessed": False,
        "residual_reduction_claimed": False,
        "pde_validated": False,
        "truth_boundary": {
            "autonomous_geometry_only": True,
            "actual_source_rectangle_centers_recovered": False,
            "actual_source_rectangle_radius_recovered": False,
            "actual_source_partition_labels_bound": False,
            "actual_source_h_sigma_bound": False,
            "actual_positive_order_background_bound": False,
            "actual_auxiliary_torus_mode_family_bound": False,
            "public_source_bound_velocity_osc_materialized": False,
            "full_ns_momentum_gate_value": None,
            "full_ns_divergence_gate_value": None,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args()
    report = run_audit()
    text = json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if args.out is None:
        print(text, end="")
    else:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    main()
