"""Independent Agent-4 audit of A1 #1179's precision-qualified total field.

A1 #1179 is a new representation identity: the late relative-swirl correction is
smaller than ordinary binary64 relative resolution, so it exposes one public
96-digit Decimal ``velocity(x,y,z,t)`` total that retains the nonzero correction.

This module treats that loaded public total as a black-box Cartesian evaluator.
Its scientific derivative path does not call A1 profile derivatives, split-profile
surfaces, or composition helpers.  It reconstructs the total Cartesian Jacobian
with an A4-owned five-point centered FD4 stencil on three preregistered local
relative physical resolutions.  Decimal subtraction/division is retained through
Jacobian formation so the precision-qualified total is not first collapsed back to
binary64.

The unchanged canonical project absolute FD ladder is checked separately for
binary64 representability.  The local relative ladder is scoped evidence only and
never replaces the canonical final-admission protocol.
"""

from __future__ import annotations

import copy
import inspect
import json
import math
import tempfile
from decimal import Decimal, ROUND_HALF_EVEN, localcontext
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .kokuno_pa16_current_cartesian_relative_swirl_decimal_composed import (
    KokunoPA16CurrentCartesianRelativeSwirlDecimalComposed,
)

TASK = "K4-VAL-122"
SCHEMA = "kokuno-a4-relative-swirl-decimal-total-independent-audit-v1"
UPSTREAM_PR = 1179
UPSTREAM_HEAD = "8c5c5b6d55a68de8285ed1dd0ebb55128f3078a4"
SEED = 9173931
TIMES = (0.31, 0.47, 0.63, 0.71)
RELATIVE_FD_LADDER = (2.0e-6, 1.0e-6, 5.0e-7)
CANONICAL_ABSOLUTE_FD_LADDER = (0.02, 0.01, 0.005)
AUDIT_DECIMAL_DIGITS = 112
EXPECTED_CANDIDATE_DECIMAL_DIGITS = 96

NORMALIZED_DIVERGENCE_GATE = 1.0e-5
ACTIVE_RATIO_FLOOR = 1.0e-45
REFINEMENT_FLOOR = 2.0e-8
REFINEMENT_WORSEN_FACTOR = 1.25
FINAL_PROJECT_MOMENTUM_GATE = 1.0e-3
FINAL_PROJECT_DIVERGENCE_GATE = 1.0e-5
FINAL_PROJECT_QUADRATURE = (24, 48, 96)

_TRUTH_BOUNDARY = {
    "precision_qualified_total_relative_swirl_independently_audited": True,
    "canonical_absolute_fd_transfer_ready": False,
    "binary64_total_relative_swirl_sum_is_resolved": False,
    "precision_qualified_relative_swirl_total_materialized": True,
    "leading_only_ns_residual_assessed": False,
    "leading_plus_oscillatory_ns_residual_assessed": False,
    "after_correction_ns_residual_assessed": False,
    "terminal_global_leading_materialized": False,
    "matching_self_contained_oscillatory_composite_materialized": False,
    "cartesian_correction_velocity_materialized": False,
    "matched_global_pressure_materialized": False,
    "restricted_forcing_materialized": False,
    "same_protocol_comparable_to_st006": False,
    "pde_validated": False,
}


def default_candidate() -> KokunoPA16CurrentCartesianRelativeSwirlDecimalComposed:
    return KokunoPA16CurrentCartesianRelativeSwirlDecimalComposed()


def _configure_decimal(ctx, precision: int = AUDIT_DECIMAL_DIGITS) -> None:
    ctx.prec = int(precision)
    ctx.rounding = ROUND_HALF_EVEN
    ctx.Emax = 999999999
    ctx.Emin = -999999999


def _as_decimal_vector(value: Any) -> tuple[Decimal, Decimal, Decimal]:
    arr = np.asarray(value, dtype=object).reshape(-1)
    if arr.size != 3:
        raise ValueError("public Decimal velocity must return a Cartesian 3-vector")
    out: list[Decimal] = []
    for raw in arr:
        if not isinstance(raw, Decimal):
            raise TypeError("scientific path requires Decimal-valued public velocity")
        if not raw.is_finite():
            raise ValueError("public Decimal velocity must be finite")
        out.append(raw)
    return out[0], out[1], out[2]


def _velocity(candidate: Any, point: Mapping[str, float | str]) -> tuple[Decimal, Decimal, Decimal]:
    return _as_decimal_vector(
        candidate.velocity(
            float(point["x"]),
            float(point["y"]),
            float(point["z"]),
            float(point["t"]),
        )
    )


def _public_similarity_coordinates(candidate: Any, x: float, y: float, z: float, t: float):
    # Probe construction may use the inherited public coordinate map.  Scientific
    # differentiation below still calls only candidate.velocity(...).
    return candidate.split.similarity_coordinates_logX(float(x), float(y), float(z), float(t))


def _q_at_z0(candidate: Any, t: float) -> float:
    coords = _public_similarity_coordinates(candidate, 1.0, 0.0, 0.0, float(t))
    q = float(np.asarray(coords["q"], dtype=float))
    if not (math.isfinite(q) and q > 0.0):
        raise RuntimeError("public similarity map returned invalid q at z=0")
    return q


def _physical_point(
    candidate: Any, *, s: float, theta: float, t: float, role: str
) -> dict[str, float | str]:
    q = _q_at_z0(candidate, t)
    log_X = float(candidate.log_X_flatten_end + s)
    log_r = 0.5 * (math.log(2.0) + math.log(q) + log_X)
    if log_r >= math.log(np.finfo(float).max):
        raise RuntimeError("held-out relative-swirl radius is not finite in binary64")
    r = math.exp(log_r)
    x = r * math.cos(theta)
    y = r * math.sin(theta)
    point: dict[str, float | str] = {
        "x": x,
        "y": y,
        "z": 0.0,
        "t": float(t),
        "r": r,
        "theta": float(theta),
        "q": q,
        "s_target": float(s),
        "log_X_target": log_X,
        "role": role,
    }
    coords = _public_similarity_coordinates(candidate, x, y, 0.0, float(t))
    got_log_X = float(np.asarray(coords["log_X"], dtype=float))
    got_eta = float(np.asarray(coords["eta"], dtype=float))
    point["log_X_reconstruction_error"] = abs(got_log_X - log_X)
    point["eta"] = got_eta
    if point["log_X_reconstruction_error"] > 5.0e-11 or abs(got_eta) > 5.0e-13:
        raise RuntimeError("held-out physical/source-coordinate reconstruction drift")
    return point


def _heldout_points(candidate: Any) -> tuple[list[dict[str, float | str]], list[dict[str, float | str]]]:
    rng = np.random.default_rng(SEED)
    width = float(candidate.compensator.width)
    y1 = float(candidate.compensator.y1)
    y2 = float(candidate.compensator.y2)
    active: list[dict[str, float | str]] = []
    controls: list[dict[str, float | str]] = []
    for t in TIMES:
        for label, center in (("bump1", y1), ("bump2", y2)):
            for k in range(2):
                offset = float(rng.uniform(-0.12 * width, 0.12 * width))
                theta = float(rng.uniform(0.31, 2.0 * math.pi - 0.31))
                active.append(
                    _physical_point(
                        candidate,
                        s=center + offset,
                        theta=theta,
                        t=t,
                        role=f"{label}_interior_{k}",
                    )
                )
        for role, s in (
            ("before_bump1", y1 - 2.0 * width),
            ("between_bumps", 0.5 * (y1 + y2)),
            ("after_bump2", y2 + 2.0 * width),
        ):
            controls.append(
                _physical_point(
                    candidate,
                    s=s,
                    theta=float(rng.uniform(0.37, 2.0 * math.pi - 0.37)),
                    t=t,
                    role=role,
                )
            )
    return active, controls


def _canonical_representability(points: list[Mapping[str, float | str]]) -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    all_ok = True
    for h in CANONICAL_ABSOLUTE_FD_LADDER:
        point_ok: list[bool] = []
        for p in points:
            axes_ok: list[bool] = []
            for key in ("x", "y", "z"):
                c = np.float64(float(p[key]))
                plus = c + np.float64(h)
                minus = c - np.float64(h)
                ok = bool(plus != c and minus != c and plus != minus)
                axes_ok.append(ok)
                all_ok = all_ok and ok
            point_ok.append(bool(all(axes_ok)))
        rows.append(
            {
                "absolute_step": h,
                "all_points_all_axes_representable": bool(all(point_ok)),
                "representable_point_fraction": float(np.mean(point_ok)),
            }
        )
    return {"ladder": rows, "all_representable": bool(all_ok)}


def _relative_steps(point: Mapping[str, float | str], rel_step: float) -> tuple[float, float, float]:
    r = float(point["r"])
    q = float(point["q"])
    h_xy = r * float(rel_step)
    h_z = max(math.sqrt(q), 1.0e-6) * float(rel_step)
    if not (math.isfinite(h_xy) and h_xy > 0.0 and math.isfinite(h_z) and h_z > 0.0):
        raise RuntimeError("invalid independent relative finite-difference step")
    return h_xy, h_xy, h_z


def _fd4_jacobian_decimal(
    candidate: Any, point: Mapping[str, float | str], rel_step: float
) -> tuple[tuple[Decimal, Decimal, Decimal], ...]:
    """A4-owned 5-point centered Cartesian FD4; total-field path calls velocity only."""
    steps = _relative_steps(point, rel_step)
    columns: list[tuple[Decimal, Decimal, Decimal]] = []
    with localcontext() as ctx:
        _configure_decimal(ctx)
        for key, h in zip(("x", "y", "z"), steps):
            c = float(point[key])
            p2 = dict(point)
            p1 = dict(point)
            m1 = dict(point)
            m2 = dict(point)
            p2[key] = c + 2.0 * h
            p1[key] = c + h
            m1[key] = c - h
            m2[key] = c - 2.0 * h
            coords = [float(m2[key]), float(m1[key]), c, float(p1[key]), float(p2[key])]
            if len(set(coords)) != 5:
                raise RuntimeError("scoped FD4 perturbation collapsed in binary64")
            fm2 = _velocity(candidate, m2)
            fm1 = _velocity(candidate, m1)
            fp1 = _velocity(candidate, p1)
            fp2 = _velocity(candidate, p2)
            denom = Decimal.from_float(12.0 * h)
            column: list[Decimal] = []
            for comp in range(3):
                num = fm2[comp] - Decimal(8) * fm1[comp] + Decimal(8) * fp1[comp] - fp2[comp]
                column.append(+(num / denom))
            columns.append((column[0], column[1], column[2]))
    # Return Jacobian by rows: J[i][j] = d_j u_i.
    return tuple(tuple(columns[j][i] for j in range(3)) for i in range(3))


def _decimal_norm(values: list[Decimal] | tuple[Decimal, ...]) -> Decimal:
    with localcontext() as ctx:
        _configure_decimal(ctx)
        total = sum((v * v for v in values), Decimal(0))
        return +total.sqrt()


def _decimal_ratio_to_float(num: Decimal, den: Decimal) -> float:
    with localcontext() as ctx:
        _configure_decimal(ctx)
        if den <= 0:
            raise ValueError("normalization denominator must be positive")
        ratio = +(abs(num) / den)
    out = float(ratio)
    if not math.isfinite(out):
        raise RuntimeError("normalized Decimal metric is not finite")
    return out


def _decimal_log10_abs(value: Decimal) -> float | None:
    if value == 0:
        return None
    with localcontext() as ctx:
        _configure_decimal(ctx)
        return float(abs(value).log10())


def _resolution_metrics(
    candidate: Any, points: list[Mapping[str, float | str]], rel_step: float
) -> dict[str, Any]:
    abs_div: list[Decimal] = []
    normalized: list[float] = []
    witnesses: list[dict[str, Any]] = []
    for p in points:
        jac = _fd4_jacobian_decimal(candidate, p, rel_step)
        with localcontext() as ctx:
            _configure_decimal(ctx)
            div = +(jac[0][0] + jac[1][1] + jac[2][2])
            flat = [jac[i][j] for i in range(3) for j in range(3)]
            grad_norm = _decimal_norm(flat)
            vel = _velocity(candidate, p)
            speed = _decimal_norm(vel)
            geom = +(speed / Decimal.from_float(max(float(p["r"]), 1.0e-300)))
            scale = max(grad_norm, geom, Decimal("1e-300"))
            nd = _decimal_ratio_to_float(div, scale)
        abs_div.append(abs(div))
        normalized.append(nd)
        witnesses.append(
            {
                "role": str(p["role"]),
                "t": float(p["t"]),
                "x": float(p["x"]),
                "y": float(p["y"]),
                "z": float(p["z"]),
                "r": float(p["r"]),
                "s_target": float(p["s_target"]),
                "abs_divergence_decimal": str(abs(div)),
                "abs_divergence_log10": _decimal_log10_abs(div),
                "gradient_frobenius_decimal": str(grad_norm),
                "speed_decimal": str(speed),
                "normalized_divergence": nd,
            }
        )
    worst = int(np.argmax(np.asarray(normalized, dtype=float)))
    with localcontext() as ctx:
        _configure_decimal(ctx)
        max_div = max(abs_div)
        rms_div = +(sum((v * v for v in abs_div), Decimal(0)) / Decimal(len(abs_div))).sqrt()
    values = np.asarray(normalized, dtype=float)
    return {
        "relative_step": float(rel_step),
        "sample_count": len(points),
        "divergence_sampled_max_decimal": str(max_div),
        "divergence_sampled_max_log10": _decimal_log10_abs(max_div),
        "divergence_sample_rms_decimal": str(rms_div),
        "divergence_sample_rms_log10": _decimal_log10_abs(rms_div),
        "normalized_divergence_sampled_max": float(np.max(values)),
        "normalized_divergence_sample_rms": float(np.sqrt(np.mean(values**2))),
        "worst_witness": witnesses[worst],
    }


def _independent_embed_float(value: float) -> Decimal:
    if not math.isfinite(float(value)):
        raise ValueError("independent parent split reference must be finite")
    with localcontext() as ctx:
        _configure_decimal(ctx, EXPECTED_CANDIDATE_DECIMAL_DIGITS)
        return +Decimal.from_float(float(value))


def _composition_checks(candidate: Any, active: list[Mapping[str, float | str]], controls: list[Mapping[str, float | str]]) -> dict[str, Any]:
    config = candidate.configuration()
    if int(config.get("decimal_digits", -1)) != EXPECTED_CANDIDATE_DECIMAL_DIGITS:
        raise ValueError("unexpected candidate Decimal precision/provenance")
    active_ratios: list[float] = []
    active_exact: list[bool] = []
    active_binary64_erased: list[bool] = []
    for p in active:
        total = _velocity(candidate, p)
        base_f, delta_f = candidate.split.velocity_split(
            float(p["x"]), float(p["y"]), float(p["z"]), float(p["t"])
        )
        base = np.asarray(base_f, dtype=float).reshape(3)
        delta = np.asarray(delta_f, dtype=float).reshape(3)
        active_binary64_erased.append(bool(np.array_equal(base + delta, base)))
        for i in range(3):
            base_d = _independent_embed_float(base[i])
            delta_d = _independent_embed_float(delta[i])
            with localcontext() as ctx:
                _configure_decimal(ctx)
                diff = +(total[i] - base_d)
            active_exact.append(bool(diff == delta_d))
        bnorm = float(np.linalg.norm(base))
        dnorm = float(np.linalg.norm(delta))
        if dnorm > 0.0:
            active_ratios.append(dnorm / max(bnorm, 1.0e-300))
    if not active_ratios:
        raise RuntimeError("active held-out set contains no nonzero relative-swirl correction")

    controls_equal: list[bool] = []
    for p in controls:
        total = _velocity(candidate, p)
        base_f, delta_f = candidate.split.velocity_split(
            float(p["x"]), float(p["y"]), float(p["z"]), float(p["t"])
        )
        base = np.asarray(base_f, dtype=float).reshape(3)
        delta = np.asarray(delta_f, dtype=float).reshape(3)
        controls_equal.append(bool(np.array_equal(delta, np.zeros(3))))
        for i in range(3):
            controls_equal.append(bool(total[i] == _independent_embed_float(base[i])))

    axis_transverse_zero: list[bool] = []
    for t in TIMES:
        axis = _as_decimal_vector(candidate.velocity(0.0, 0.0, 0.0, float(t)))
        axis_transverse_zero.extend((axis[0] == 0, axis[1] == 0))

    return {
        "active_total_minus_base_replays_delta_exact": bool(all(active_exact)),
        "active_binary64_base_plus_delta_erased": bool(all(active_binary64_erased)),
        "active_correction_ratio_min": float(min(active_ratios)),
        "active_correction_ratio_max": float(max(active_ratios)),
        "off_bump_total_equals_independent_base_exact": bool(all(controls_equal)),
        "exact_axis_transverse_velocity_zero": bool(all(axis_transverse_zero)),
    }


def _save_load_replay_exact(loaded: Any, original: Any, points: list[Mapping[str, float | str]]) -> bool:
    for p in points:
        if _velocity(loaded, p) != _velocity(original, p):
            return False
    return True


def _decimal_export_replay_exact(candidate: Any, points: list[Mapping[str, float | str]]) -> bool:
    if not points:
        return False
    chosen = points[: min(4, len(points))]
    x = np.asarray([float(p["x"]) for p in chosen])
    y = np.asarray([float(p["y"]) for p in chosen])
    z = np.asarray([float(p["z"]) for p in chosen])
    t = np.asarray([float(p["t"]) for p in chosen])
    expected = np.asarray(candidate.velocity(x, y, z, t), dtype=object)
    with tempfile.TemporaryDirectory() as td:
        path = Path(td) / "velocity_decimal.json"
        candidate.export_velocity_decimal_json(path, x, y, z, t)
        payload = json.loads(path.read_text())
    if payload.get("candidate_semantic_sha256") != candidate.semantic_sha256:
        return False
    if payload.get("shape") != list(expected.shape):
        return False
    replay = np.asarray([Decimal(v) for v in payload["values_row_major"]], dtype=object).reshape(expected.shape)
    return bool(all(replay[idx] == expected[idx] for idx in np.ndindex(expected.shape)))


def _assess(report: Mapping[str, Any]) -> list[str]:
    failures: list[str] = []
    fine = report["resolutions"][-1]
    medium = report["resolutions"][-2]
    comp = report["composition_checks"]
    if float(fine["normalized_divergence_sampled_max"]) > NORMALIZED_DIVERGENCE_GATE:
        failures.append("fine_total_normalized_divergence_max")
    if float(fine["normalized_divergence_sample_rms"]) > NORMALIZED_DIVERGENCE_GATE:
        failures.append("fine_total_normalized_divergence_rms")
    if not comp["active_total_minus_base_replays_delta_exact"]:
        failures.append("decimal_total_does_not_replay_parent_delta")
    if not comp["active_binary64_base_plus_delta_erased"]:
        failures.append("expected_binary64_erasure_not_observed")
    if float(comp["active_correction_ratio_max"]) < ACTIVE_RATIO_FLOOR:
        failures.append("relative_swirl_correction_collapsed")
    if float(comp["active_correction_ratio_max"]) >= float(np.finfo(float).eps):
        failures.append("correction_is_not_sub_binary64_epsilon")
    if not comp["off_bump_total_equals_independent_base_exact"]:
        failures.append("off_bump_total_base_mismatch")
    if not comp["exact_axis_transverse_velocity_zero"]:
        failures.append("axis_transverse_velocity_nonzero")
    if not report["save_load_replay_exact"]:
        failures.append("save_load_total_replay")
    if not report["decimal_string_export_replay_exact"]:
        failures.append("decimal_string_export_replay")
    f = float(fine["normalized_divergence_sampled_max"])
    m = float(medium["normalized_divergence_sampled_max"])
    if f > REFINEMENT_FLOOR and f > REFINEMENT_WORSEN_FACTOR * max(m, 1.0e-300):
        failures.append("medium_to_fine_total_divergence_worsened")
    return failures


def audit_loaded_decimal_total(loaded: Any, pre_serialization_reference: Any) -> dict[str, Any]:
    """Run fixed K4-VAL-122.  Public scientific API intentionally has no tuning knobs."""
    if not isinstance(loaded, KokunoPA16CurrentCartesianRelativeSwirlDecimalComposed):
        raise TypeError("loaded candidate must be exact #1179 Decimal-composed type")
    truth = dict(loaded.truth_boundary)
    required_false = (
        "binary64_total_relative_swirl_sum_is_resolved",
        "unified_global_cartesian_velocity_export_ready",
        "matched_global_pressure_materialized",
        "restricted_forcing_materialized",
        "heldout_ns_residual_assessed",
        "pde_validated",
    )
    for key in required_false:
        if bool(truth.get(key, False)):
            raise ValueError(f"upstream truth promotion is not admissible: {key}")
    if not bool(truth.get("precision_qualified_public_velocity_materialized", False)):
        raise ValueError("precision-qualified public total velocity is required")

    active, controls = _heldout_points(loaded)
    canonical = _canonical_representability(active)
    composition = _composition_checks(loaded, active, controls)
    resolutions = [_resolution_metrics(loaded, active, h) for h in RELATIVE_FD_LADDER]
    replay = _save_load_replay_exact(loaded, pre_serialization_reference, active + controls)
    export_replay = _decimal_export_replay_exact(loaded, active)

    report: dict[str, Any] = {
        "schema": SCHEMA,
        "task": TASK,
        "upstream_pr": UPSTREAM_PR,
        "upstream_head": UPSTREAM_HEAD,
        "candidate_semantic_sha256": loaded.semantic_sha256,
        "seed": SEED,
        "times": list(TIMES),
        "independent_operator": "centered_cartesian_five_point_fd4_on_public_decimal_total",
        "audit_decimal_digits": AUDIT_DECIMAL_DIGITS,
        "candidate_decimal_digits": EXPECTED_CANDIDATE_DECIMAL_DIGITS,
        "relative_fd_ladder": list(RELATIVE_FD_LADDER),
        "canonical_absolute_fd_ladder": list(CANONICAL_ABSOLUTE_FD_LADDER),
        "active_points": active,
        "control_points": controls,
        "resolutions": resolutions,
        "composition_checks": composition,
        "canonical_absolute_fd_representability": canonical,
        "save_load_replay_exact": replay,
        "decimal_string_export_replay_exact": export_replay,
        "gates": {
            "scoped_normalized_divergence_gate": NORMALIZED_DIVERGENCE_GATE,
            "active_ratio_floor": ACTIVE_RATIO_FLOOR,
            "active_ratio_ceiling_float64_epsilon": float(np.finfo(float).eps),
            "refinement_floor": REFINEMENT_FLOOR,
            "refinement_worsen_factor": REFINEMENT_WORSEN_FACTOR,
            "final_project_momentum_gate_unchanged": FINAL_PROJECT_MOMENTUM_GATE,
            "final_project_divergence_gate_unchanged": FINAL_PROJECT_DIVERGENCE_GATE,
            "final_project_quadrature_unchanged": list(FINAL_PROJECT_QUADRATURE),
        },
        "truth_boundary": copy.deepcopy(_TRUTH_BOUNDARY),
    }
    failures = _assess(report)
    report["failures"] = failures
    report["scoped_decimal_total_gate_passed"] = len(failures) == 0
    report["truth_boundary"]["canonical_absolute_fd_transfer_ready"] = bool(
        canonical["all_representable"] and len(failures) == 0
    )
    report["final_project_admission_ready"] = False
    blockers = [
        "terminal_global_leading_missing",
        "matching_self_contained_oscillatory_composite_missing",
        "cartesian_correction_velocity_missing",
        "matched_pressure_missing",
        "preregistered_restricted_forcing_missing",
        "complete_ns_residual_not_assessed",
    ]
    if not canonical["all_representable"]:
        blockers.insert(0, "canonical_absolute_fd_ladder_not_binary64_representable_here")
    report["final_project_blockers"] = blockers
    return report


def enforce_preregistered_gates(report: Mapping[str, Any]) -> None:
    failures = _assess(report)
    if failures:
        raise AssertionError("K4-VAL-122 scoped gate failure: " + ", ".join(failures))
    if report.get("truth_boundary", {}).get("pde_validated") is not False:
        raise AssertionError("pde_validated must remain false")
    if report.get("final_project_admission_ready") is not False:
        raise AssertionError("scoped Decimal-total audit cannot promote final project admission")


def public_api_has_no_scientific_tuning_knobs() -> bool:
    params = set(inspect.signature(audit_loaded_decimal_total).parameters)
    forbidden = {
        "step",
        "h",
        "threshold",
        "residual",
        "forcing",
        "force",
        "pressure",
        "viscosity",
        "precision",
        "decimal_digits",
        "gain",
        "optimizer",
        "tolerance",
    }
    return params == {"loaded", "pre_serialization_reference"} and not (params & forbidden)
