"""Certify a robustness margin for one frozen Kokuno correction step.

Agent-3 PR #431 froze a scalar damping on calibration cells and transferred it
without retuning to disjoint theta/z compact-stress cells.  PR #440 then checked
that the same frozen step also respects the corrected-reader radial divergence

    (div T)_r = + d_z sigma_1.

This module does not choose a new damping.  Instead it asks a narrower finite-
correction-cycle question: how far can the *already frozen* scalar be perturbed
before any fresh-cell gain guard fails?

For a fixed physical update V, the retained covariance change is exactly

    Delta C(lambda) = lambda B(W,V) + lambda^2 C(V).

Consequently every squared theta/z stress residual and every squared radial
residual proxy is a quartic polynomial in lambda.  We form those quartics on the
same disjoint validation cells, solve their real sublevel intervals on [0,1],
and intersect the connected intervals containing the frozen #431 damping.

The interval is a repository correction-gain diagnostic, not a Kokuno theorem,
not validation-driven line search, and not a Navier--Stokes residual gate.  The
frozen damping, sign, coefficient budget, validation cells and all scientific
thresholds remain unchanged.
"""
from __future__ import annotations

import argparse
from dataclasses import replace
import json
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
from numpy.polynomial import Polynomial

from .kokuno_bounded_coordinate_covariance import _one_label_real_family
from .kokuno_bounded_family_coefficients import KokunoBoundedFamilyCoefficientCoordinates
from .kokuno_complete_curl import KokunoCompleteCurlCorrection
from .kokuno_leading_core_series import KokunoLeadingCoreSeriesCandidate
from .kokuno_mean_defect import PhaseMeanDefectContract
from .kokuno_missing_covariance_column_target import (
    DEFAULT_ANGULAR_COUNT,
    DEFAULT_PHASE_COUNT,
    DEFAULT_RADIAL_COUNT,
)
from .kokuno_quadratic_covariance_gain_guard import measure_bounded_coordinate_quadratic_change
from .kokuno_radial_aware_damping_holdout import (
    damped_axial_covariance_stress,
    evaluate_fixed_damping_radial_cell,
)
from .kokuno_radial_force_completeness_audit import (
    Z_DERIVATIVE_STEP_LADDER,
    centered_radial_force,
    sample_real_radial_mean_defect,
)
from .kokuno_signed_covariance_inverse import (
    ROUTED_OSCILLATORY_AMPLITUDE,
    ROUTED_OSCILLATORY_PHASE,
)
from .kokuno_spacetime_damping_holdout import (
    VALIDATION_TIMES,
    VALIDATION_Z,
    _build_cells,
    generate_actual_holdout_report,
)

TASK = "KOKUNO-A3-JOINT-GAIN-MARGIN-026"
SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_EDITION = "2026.09.09-consolidated"
SOURCE_DOI = "10.5281/zenodo.22678406"
SOURCE_FORMULAS = "R33-R34/R41 radial stress; signed covariance inverse with retained self-covariance"
AGENT3_DAMPING_HOLDOUT_PR = 431
AGENT3_RADIAL_AWARE_PR = 440
AGENT3_FAMILY_INVERSE_PR = 393

# Purely numerical root classification for degree <= 4 polynomials.  This is
# not a scientific/PDE tolerance and is recorded in every report.
ROOT_IMAG_TOLERANCE = 1.0e-9
ROOT_CLUSTER_TOLERANCE = 1.0e-10


def _squared_norm(values: np.ndarray) -> float:
    array = np.asarray(values, dtype=float)
    if array.size == 0 or not np.isfinite(array).all():
        raise ValueError("values must be nonempty and finite")
    return float(np.sum(array * array))


def quartic_relative_rms_sublevel_coefficients(
    base: np.ndarray,
    linear: np.ndarray,
    quadratic: np.ndarray,
    *,
    relative_bound: float,
) -> np.ndarray:
    """Return ascending quartic coefficients for an RMS sublevel condition.

    The residual is ``base + lambda*linear + lambda^2*quadratic``.  The
    returned polynomial p satisfies

        p(lambda) <= 0

    iff its Euclidean/RMS norm is no larger than ``relative_bound`` times the
    norm of ``base``.  Using sums instead of means is exact because the common
    number of entries cancels from the relative RMS comparison.
    """
    b = np.asarray(base, dtype=float)
    l = np.asarray(linear, dtype=float)
    q = np.asarray(quadratic, dtype=float)
    bound = float(relative_bound)
    if b.shape != l.shape or b.shape != q.shape or b.size == 0:
        raise ValueError("base, linear and quadratic must have one common nonempty shape")
    if not np.isfinite(b).all() or not np.isfinite(l).all() or not np.isfinite(q).all():
        raise ValueError("quartic inputs must be finite")
    if not np.isfinite(bound) or bound < 0.0:
        raise ValueError("relative_bound must be finite and nonnegative")
    b2 = _squared_norm(b)
    if b2 <= np.finfo(float).tiny:
        raise ValueError("base norm must be nonzero")
    dot_bl = float(np.sum(b * l))
    dot_bq = float(np.sum(b * q))
    dot_lq = float(np.sum(l * q))
    return np.asarray(
        [
            (1.0 - bound * bound) * b2,
            2.0 * dot_bl,
            _squared_norm(l) + 2.0 * dot_bq,
            2.0 * dot_lq,
            _squared_norm(q),
        ],
        dtype=float,
    )


def relative_rms(
    base: np.ndarray,
    linear: np.ndarray,
    quadratic: np.ndarray,
    damping: float,
) -> float:
    b = np.asarray(base, dtype=float)
    l = np.asarray(linear, dtype=float)
    q = np.asarray(quadratic, dtype=float)
    lam = float(damping)
    if b.shape != l.shape or b.shape != q.shape or b.size == 0:
        raise ValueError("base, linear and quadratic must share one nonempty shape")
    if not np.isfinite(lam):
        raise ValueError("damping must be finite")
    denominator = max(float(np.sqrt(np.mean(b * b))), np.finfo(float).tiny)
    residual = b + lam * l + lam * lam * q
    return float(np.sqrt(np.mean(residual * residual)) / denominator)


def _numeric_slack(coefficients: np.ndarray) -> float:
    coeff = np.asarray(coefficients, dtype=float)
    return float(1024.0 * np.finfo(float).eps * max(1.0, float(np.sum(np.abs(coeff)))))


def _cluster_real_roots(values: Sequence[float], *, lower: float, upper: float) -> list[float]:
    roots: list[float] = []
    for value in sorted(float(v) for v in values if lower <= float(v) <= upper):
        clipped = min(upper, max(lower, value))
        if not roots or abs(clipped - roots[-1]) > ROOT_CLUSTER_TOLERANCE:
            roots.append(clipped)
        else:
            roots[-1] = 0.5 * (roots[-1] + clipped)
    return roots


def connected_sublevel_interval(
    coefficients: Sequence[float],
    *,
    point: float,
    lower: float = 0.0,
    upper: float = 1.0,
) -> dict[str, Any]:
    """Find the connected p(lambda)<=0 component containing ``point``.

    Degree is at most four in this route.  Real roots are obtained from
    ``numpy.polynomial.Polynomial.roots``; tiny imaginary parts are treated as
    numerical root-solver noise under ``ROOT_IMAG_TOLERANCE``.  This numerical
    classification tolerance is not a scientific gain threshold.
    """
    coeff = np.asarray(tuple(float(v) for v in coefficients), dtype=float)
    x = float(point)
    lo = float(lower)
    hi = float(upper)
    if coeff.ndim != 1 or coeff.size < 1 or coeff.size > 5 or not np.isfinite(coeff).all():
        raise ValueError("coefficients must be finite with degree at most four")
    if not np.isfinite([x, lo, hi]).all() or not lo < hi or not lo <= x <= hi:
        raise ValueError("require finite lower <= point <= upper with lower < upper")

    polynomial = Polynomial(coeff)
    slack = _numeric_slack(coeff)
    point_value = float(polynomial(x))
    if point_value > slack:
        raise ValueError("frozen damping is outside this gain sublevel")

    real_roots: list[float] = []
    for root in polynomial.roots():
        real = float(np.real(root))
        imag = float(np.imag(root))
        if abs(imag) <= ROOT_IMAG_TOLERANCE * max(1.0, abs(real)):
            if lo - ROOT_CLUSTER_TOLERANCE <= real <= hi + ROOT_CLUSTER_TOLERANCE:
                real_roots.append(min(hi, max(lo, real)))
    roots = _cluster_real_roots(real_roots, lower=lo, upper=hi)
    cuts = _cluster_real_roots([lo, *roots, hi], lower=lo, upper=hi)
    if cuts[0] != lo:
        cuts.insert(0, lo)
    if cuts[-1] != hi:
        cuts.append(hi)

    passing: list[tuple[float, float]] = []
    for left, right in zip(cuts, cuts[1:]):
        if right - left <= ROOT_CLUSTER_TOLERANCE:
            continue
        mid = 0.5 * (left + right)
        if float(polynomial(mid)) <= slack:
            passing.append((left, right))

    # Handle the identically-zero / no-real-root case directly.
    if not passing and float(polynomial(0.5 * (lo + hi))) <= slack:
        passing = [(lo, hi)]

    containing = [
        interval
        for interval in passing
        if interval[0] - ROOT_CLUSTER_TOLERANCE <= x <= interval[1] + ROOT_CLUSTER_TOLERANCE
    ]
    if not containing:
        # ``point`` can sit exactly at a repeated root with a passing interval
        # on only one side.  Accept that side if the polynomial itself passes.
        candidates = [
            interval
            for interval in passing
            if min(abs(x - interval[0]), abs(x - interval[1])) <= ROOT_CLUSTER_TOLERANCE
        ]
        containing = candidates
    if not containing:
        raise RuntimeError("could not identify the connected sublevel component containing point")

    left = min(interval[0] for interval in containing)
    right = max(interval[1] for interval in containing)
    # Merge adjacent passing components separated only by a root at which the
    # polynomial is zero.  Repeating until fixed point handles multiple roots.
    changed = True
    while changed:
        changed = False
        for candidate_left, candidate_right in passing:
            if abs(candidate_right - left) <= ROOT_CLUSTER_TOLERANCE:
                left = candidate_left
                changed = True
            if abs(candidate_left - right) <= ROOT_CLUSTER_TOLERANCE:
                right = candidate_right
                changed = True

    return {
        "point": x,
        "point_polynomial_value": point_value,
        "sublevel_interval": [float(left), float(right)],
        "left_margin": float(x - left),
        "right_margin": float(right - x),
        "real_roots_in_domain": [float(v) for v in roots],
        "numerical_polynomial_slack": slack,
        "root_imag_tolerance": ROOT_IMAG_TOLERANCE,
        "root_cluster_tolerance": ROOT_CLUSTER_TOLERANCE,
    }


def intersect_gain_constraints(
    constraints: Iterable[Mapping[str, Any]],
    *,
    frozen_damping: float,
) -> dict[str, Any]:
    """Intersect all connected gain intervals containing the frozen damping."""
    lam = float(frozen_damping)
    rows: list[dict[str, Any]] = []
    for constraint in constraints:
        row = dict(constraint)
        name = str(row.get("name", ""))
        if not name:
            raise ValueError("every gain constraint requires a nonempty name")
        interval = connected_sublevel_interval(row["quartic_coefficients"], point=lam)
        row["gain_interval"] = interval
        rows.append(row)
    if not rows:
        raise ValueError("at least one gain constraint is required")

    left = max(float(row["gain_interval"]["sublevel_interval"][0]) for row in rows)
    right = min(float(row["gain_interval"]["sublevel_interval"][1]) for row in rows)
    if left > lam + ROOT_CLUSTER_TOLERANCE or right < lam - ROOT_CLUSTER_TOLERANCE:
        raise RuntimeError("gain-constraint intersection does not contain frozen damping")
    left = min(lam, left)
    right = max(lam, right)
    left_margin = lam - left
    right_margin = right - lam
    binding_left = [
        row["name"]
        for row in rows
        if abs(float(row["gain_interval"]["sublevel_interval"][0]) - left)
        <= 4.0 * ROOT_CLUSTER_TOLERANCE
    ]
    binding_right = [
        row["name"]
        for row in rows
        if abs(float(row["gain_interval"]["sublevel_interval"][1]) - right)
        <= 4.0 * ROOT_CLUSTER_TOLERANCE
    ]
    strict_interior = bool(left_margin > ROOT_CLUSTER_TOLERANCE and right_margin > ROOT_CLUSTER_TOLERANCE)
    return {
        "frozen_damping": lam,
        "constraint_count": len(rows),
        "joint_admissible_interval": [float(left), float(right)],
        "left_margin": float(left_margin),
        "right_margin": float(right_margin),
        "minimum_two_sided_margin": float(min(left_margin, right_margin)),
        "joint_interval_width": float(right - left),
        "binding_left_constraints": binding_left,
        "binding_right_constraints": binding_right,
        "frozen_damping_strictly_inside_joint_interval": strict_interior,
        "validation_used_for_damping_selection": False,
        "frozen_damping_changed": False,
        "constraints": rows,
    }


def _validate_axis(values: Iterable[float], *, name: str) -> tuple[float, ...]:
    result = tuple(float(value) for value in values)
    if not result or not np.isfinite(result).all():
        raise ValueError(f"{name} must contain finite values")
    if any(right <= left for left, right in zip(result, result[1:])):
        raise ValueError(f"{name} must be strictly increasing")
    return result


def _tangential_gain_constraints(
    rows: Iterable[Mapping[str, Any]],
    *,
    frozen_damping: float,
) -> list[dict[str, Any]]:
    constraints: list[dict[str, Any]] = []
    for source_row in rows:
        row = dict(source_row)
        target = np.asarray(row["target_stress"], dtype=float)
        active = np.asarray(row["active_target_mask"], dtype=bool)
        measurement = dict(row["measurement"])
        linear = np.asarray(measurement["linear_covariance_change_theta_axial"], dtype=float)
        quadratic = np.asarray(measurement["self_covariance_theta_axial"], dtype=float)
        if target.ndim != 2 or target.shape[1] != 2 or active.shape != (len(target),):
            raise ValueError("invalid tangential target shape")
        if linear.shape != target.shape or quadratic.shape != target.shape or not np.any(active):
            raise ValueError("invalid tangential covariance arrays or active mask")
        base = target[active]
        l = -linear[active]
        q = -quadratic[active]
        full_relative = relative_rms(base, l, q, 1.0)
        inherited_bound = min(1.0, full_relative)
        frozen_relative = relative_rms(base, l, q, frozen_damping)
        coefficients = quartic_relative_rms_sublevel_coefficients(
            base,
            l,
            q,
            relative_bound=inherited_bound,
        )
        constraints.append(
            {
                "name": f"tangential_t={float(row['time']):.10g}_z={float(row['z']):.10g}",
                "channel": "theta_z_compact_stress",
                "time": float(row["time"]),
                "z": float(row["z"]),
                "active_nodes": int(np.count_nonzero(active)),
                "inherited_relative_bound": float(inherited_bound),
                "full_step_relative_rms": float(full_relative),
                "frozen_relative_rms": float(frozen_relative),
                "quartic_coefficients": coefficients.tolist(),
            }
        )
    return constraints


def _radial_gain_constraint_for_cell(
    *,
    time: float,
    z: float,
    radii: np.ndarray,
    gated_radial_defect: np.ndarray,
    family: Any,
    coordinates: KokunoBoundedFamilyCoefficientCoordinates,
    delta_common: float,
    frozen_damping: float,
    angular_count: int,
    phase_count: int,
) -> tuple[dict[str, Any], dict[str, Any]]:
    active = np.ones(len(radii), dtype=bool)
    active[[0, -1]] = False
    linear_derivatives: list[np.ndarray] = []
    quadratic_derivatives: list[np.ndarray] = []
    combined_derivatives: list[np.ndarray] = []
    levels: list[dict[str, Any]] = []
    for step in Z_DERIVATIVE_STEP_LADDER:
        plus = measure_bounded_coordinate_quadratic_change(
            family,
            radii,
            coordinate_contract=coordinates,
            delta_common=delta_common,
            delta_band=0.0,
            time=float(time),
            z=float(z + step),
            angular_count=int(angular_count),
            phase_count=int(phase_count),
        )
        minus = measure_bounded_coordinate_quadratic_change(
            family,
            radii,
            coordinate_contract=coordinates,
            delta_common=delta_common,
            delta_band=0.0,
            time=float(time),
            z=float(z - step),
            angular_count=int(angular_count),
            phase_count=int(phase_count),
        )
        linear_plus = np.asarray(plus["linear_covariance_change_theta_axial"], dtype=float)[:, 1]
        linear_minus = np.asarray(minus["linear_covariance_change_theta_axial"], dtype=float)[:, 1]
        quadratic_plus = np.asarray(plus["self_covariance_theta_axial"], dtype=float)[:, 1]
        quadratic_minus = np.asarray(minus["self_covariance_theta_axial"], dtype=float)[:, 1]
        b_force = centered_radial_force(linear_plus, linear_minus, float(step))
        c_force = centered_radial_force(quadratic_plus, quadratic_minus, float(step))
        combined = frozen_damping * b_force + frozen_damping * frozen_damping * c_force
        # Independent algebra check against the same direct damped stress used by #440.
        direct = centered_radial_force(
            damped_axial_covariance_stress(plus, frozen_damping),
            damped_axial_covariance_stress(minus, frozen_damping),
            float(step),
        )
        scale = max(float(np.sqrt(np.mean(direct[active] * direct[active]))), np.finfo(float).tiny)
        identity = float(np.sqrt(np.mean((combined[active] - direct[active]) ** 2)) / scale)
        if identity > 4096.0 * np.finfo(float).eps:
            raise RuntimeError("radial linear/quadratic derivative decomposition failed")
        linear_derivatives.append(b_force)
        quadratic_derivatives.append(c_force)
        combined_derivatives.append(combined)
        levels.append(
            {
                "z_step": float(step),
                "linear_quadratic_to_direct_relative_rms": identity,
                "plus_quadratic_identity_passed": bool(plus["quadratic_identity_passed"]),
                "minus_quadratic_identity_passed": bool(minus["quadratic_identity_passed"]),
            }
        )

    frozen_radial = evaluate_fixed_damping_radial_cell(
        gated_radial_defect,
        combined_derivatives,
        z_steps=Z_DERIVATIVE_STEP_LADDER,
        active_mask=active,
    )
    b_finest = linear_derivatives[-1][active]
    c_finest = quadratic_derivatives[-1][active]
    base = np.asarray(gated_radial_defect, dtype=float)[active]
    coefficients = quartic_relative_rms_sublevel_coefficients(
        base,
        b_finest,
        c_finest,
        relative_bound=1.0,
    )
    frozen_relative = relative_rms(base, b_finest, c_finest, frozen_damping)
    if not np.isclose(
        frozen_relative,
        float(frozen_radial["radial_relative_residual_rms_after_source_sign"]),
        rtol=0.0,
        atol=2048.0 * np.finfo(float).eps,
    ):
        raise RuntimeError("quartic radial objective disagrees with #440 frozen-step objective")
    constraint = {
        "name": f"radial_t={float(time):.10g}_z={float(z):.10g}",
        "channel": "source_plus_d_z_sigma_1",
        "time": float(time),
        "z": float(z),
        "active_nodes": int(np.count_nonzero(active)),
        "inherited_relative_bound": 1.0,
        "frozen_relative_rms": float(frozen_relative),
        "quartic_coefficients": coefficients.tolist(),
    }
    audit = {
        "time": float(time),
        "z": float(z),
        "levels": levels,
        "frozen_radial_result": frozen_radial,
    }
    return constraint, audit


def generate_actual_gain_margin_report(
    *,
    output: str | Path = "artifacts/kokuno_agent3/joint_gain_margin_report.json",
    parent_output: str | Path = "artifacts/kokuno_agent3/joint_gain_margin_parent_431.json",
    validation_times: Iterable[float] = VALIDATION_TIMES,
    validation_z: Iterable[float] = VALIDATION_Z,
    radial_count: int = DEFAULT_RADIAL_COUNT,
    angular_count: int = DEFAULT_ANGULAR_COUNT,
    phase_count: int = DEFAULT_PHASE_COUNT,
) -> dict[str, Any]:
    fresh_times = _validate_axis(validation_times, name="validation_times")
    fresh_z = _validate_axis(validation_z, name="validation_z")
    parent = generate_actual_holdout_report(
        output=parent_output,
        validation_times=fresh_times,
        validation_z=fresh_z,
        radial_count=int(radial_count),
        angular_count=int(angular_count),
        phase_count=int(phase_count),
    )
    frozen_damping = float(parent["damping_calibration"]["shared_result"]["shared_selected_damping"])
    fractional_budget = float(parent["inputs"]["converted_fractional_budget"])
    delta_common = float(parent["inputs"]["frozen_delta_common"])
    coordinates = KokunoBoundedFamilyCoefficientCoordinates(max_l1_update=fractional_budget)
    family = _one_label_real_family()

    validation_rows, validation_receipts = _build_cells(
        times=fresh_times,
        z_values=fresh_z,
        family=family,
        coordinates=coordinates,
        delta_common=delta_common,
        radial_count=int(radial_count),
        angular_count=int(angular_count),
        phase_count=int(phase_count),
    )
    tangential_constraints = _tangential_gain_constraints(
        validation_rows,
        frozen_damping=frozen_damping,
    )

    leading = KokunoLeadingCoreSeriesCandidate()
    oscillatory = replace(
        KokunoCompleteCurlCorrection(),
        amplitude=ROUTED_OSCILLATORY_AMPLITUDE,
        phase=ROUTED_OSCILLATORY_PHASE,
    )
    defect_contract = PhaseMeanDefectContract(phase_count=int(phase_count))
    radial_constraints: list[dict[str, Any]] = []
    radial_audits: list[dict[str, Any]] = []
    receipt_by_pair = {
        (float(time), float(z)): (receipt, metadata)
        for time, z, receipt, metadata in validation_receipts
    }
    for z in fresh_z:
        for time in fresh_times:
            receipt, metadata = receipt_by_pair[(float(time), float(z))]
            radii = np.asarray(receipt.radii, dtype=float)
            _, gated_radial = sample_real_radial_mean_defect(
                leading.at_points,
                oscillatory,
                defect_contract,
                radii,
                time=float(time),
                z=float(z),
                angular_count=int(angular_count),
            )
            constraint, audit = _radial_gain_constraint_for_cell(
                time=float(time),
                z=float(z),
                radii=radii,
                gated_radial_defect=gated_radial,
                family=family,
                coordinates=coordinates,
                delta_common=delta_common,
                frozen_damping=frozen_damping,
                angular_count=int(angular_count),
                phase_count=int(phase_count),
            )
            constraint["surrogate_defect_used"] = False
            audit["leading_candidate_sha256"] = str(metadata["leading_candidate_sha256"])
            radial_constraints.append(constraint)
            radial_audits.append(audit)

    margin = intersect_gain_constraints(
        [*tangential_constraints, *radial_constraints],
        frozen_damping=frozen_damping,
    )
    parent_transfer = bool(parent["disjoint_validation"]["result"]["validation_transfer_preflight_passed"])
    radial_frozen_pass = all(
        bool(row["frozen_radial_result"]["radial_source_sign_not_worse"])
        and bool(
            row["frozen_radial_result"]["radial_force_derivative_audit"][
                "radial_force_derivative_stability_preflight_passed"
            ]
        )
        for row in radial_audits
    )
    joint_pass = bool(
        parent_transfer
        and radial_frozen_pass
        and margin["frozen_damping_strictly_inside_joint_interval"]
    )
    parent_rank_passed = bool(parent["routing"]["upstream_bounded_inverse_passed"])
    leading_shas = {str(row["leading_candidate_sha256"]) for row in radial_audits}
    if len(leading_shas) != 1:
        raise RuntimeError("leading candidate changed across gain-margin validation cells")

    report: dict[str, Any] = {
        "task": TASK,
        "source": {
            "repository": SOURCE_REPOSITORY,
            "commit": SOURCE_COMMIT,
            "path": SOURCE_PATH,
            "edition": SOURCE_EDITION,
            "doi": SOURCE_DOI,
            "formulas": SOURCE_FORMULAS,
            "scope": (
                "The source supplies the compact stress/radial-divergence structure and retains "
                "nonlinear covariance. Quartic gain-margin algebra and the connected-interval "
                "diagnostic are repository finite-step engineering."
            ),
        },
        "handoff": {
            "agent3_damping_holdout_pr": AGENT3_DAMPING_HOLDOUT_PR,
            "agent3_radial_aware_pr": AGENT3_RADIAL_AWARE_PR,
            "agent3_family_inverse_pr": AGENT3_FAMILY_INVERSE_PR,
            "agent1_latest_dependency": (
                "PR #443 instantiates the PA.10 T_sh numerical certificate but still leaves the "
                "source pressure datum/global leading completion fail-closed."
            ),
            "agent2_latest_dependency": (
                "PR #439 certifies rational rectangle separation but still does not materialize "
                "the actual positive-order/background multi-band public physical velocity family."
            ),
        },
        "inputs": {
            "validation_times": list(fresh_times),
            "validation_z": list(fresh_z),
            "radial_count": int(radial_count),
            "angular_count": int(angular_count),
            "phase_count": int(phase_count),
            "frozen_damping_from_pr431": frozen_damping,
            "frozen_delta_common": delta_common,
            "frozen_delta_band": 0.0,
            "converted_fractional_budget": fractional_budget,
            "budget_changed": False,
            "damping_changed": False,
            "validation_used_for_damping_selection": False,
            "surrogate_defect_used": False,
            "root_imag_tolerance_numerical_only": ROOT_IMAG_TOLERANCE,
            "root_cluster_tolerance_numerical_only": ROOT_CLUSTER_TOLERANCE,
        },
        "parent_disjoint_tangential_transfer": parent["disjoint_validation"]["result"],
        "radial_derivative_audits": radial_audits,
        "joint_gain_margin": margin,
        "routing": {
            "joint_gain_margin_preflight_passed": joint_pass,
            "parent_tangential_transfer_passed": parent_transfer,
            "radial_frozen_step_passed": radial_frozen_pass,
            "genuinely_independent_second_covariance_column_ready": False,
            "upstream_bounded_inverse_passed": parent_rank_passed,
            "public_velocity_correction_materialized": False,
            "finite_correction_cycle_rerun_allowed": False,
            "next_required": (
                "Agent 2 must still supply an actual source-motivated multi-band physical family "
                "with a genuine independent covariance direction. A future materialized correction "
                "must pass the existing unit/rank/budget/spacetime/quadratic/radial guards and this "
                "gain-margin diagnostic before a held-in/held-out NS correction cycle is run."
            ),
        },
        "truth_boundary": {
            "real_candidate_defect_targets_consumed": True,
            "surrogate_defect_used": False,
            "damping_calibration_validation_disjoint": True,
            "heldout_damping_reoptimized": False,
            "validation_objective_used_for_selection": False,
            "joint_gain_margin_is_diagnostic_only": True,
            "actual_source_multiband_family_consumed": False,
            "genuinely_independent_second_covariance_column_ready": False,
            "public_velocity_correction_materialized": False,
            "finite_correction_cycle_run": False,
            "residual_reduction_claimed": False,
            "formal_full_domain_normalized_pde_gate_assessed": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
    }
    destination = Path(output)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return report


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        default="artifacts/kokuno_agent3/joint_gain_margin_report.json",
    )
    parser.add_argument(
        "--parent-output",
        default="artifacts/kokuno_agent3/joint_gain_margin_parent_431.json",
    )
    parser.add_argument("--radial-count", type=int, default=DEFAULT_RADIAL_COUNT)
    parser.add_argument("--angular-count", type=int, default=DEFAULT_ANGULAR_COUNT)
    parser.add_argument("--phase-count", type=int, default=DEFAULT_PHASE_COUNT)
    return parser


def main() -> None:
    args = _parser().parse_args()
    report = generate_actual_gain_margin_report(
        output=args.output,
        parent_output=args.parent_output,
        radial_count=args.radial_count,
        angular_count=args.angular_count,
        phase_count=args.phase_count,
    )
    margin = report["joint_gain_margin"]
    print(
        json.dumps(
            {
                "task": report["task"],
                "frozen_damping": margin["frozen_damping"],
                "joint_admissible_interval": margin["joint_admissible_interval"],
                "left_margin": margin["left_margin"],
                "right_margin": margin["right_margin"],
                "minimum_two_sided_margin": margin["minimum_two_sided_margin"],
                "binding_left_constraints": margin["binding_left_constraints"],
                "binding_right_constraints": margin["binding_right_constraints"],
                "joint_gain_margin_preflight_passed": report["routing"][
                    "joint_gain_margin_preflight_passed"
                ],
                "finite_correction_cycle_rerun_allowed": report["routing"][
                    "finite_correction_cycle_rerun_allowed"
                ],
            },
            sort_keys=True,
            allow_nan=False,
        )
    )


if __name__ == "__main__":
    main()
