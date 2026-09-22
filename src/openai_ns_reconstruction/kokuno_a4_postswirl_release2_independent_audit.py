"""Independent A4 audit of A1 #1204's post-swirl ``l:-1 -> -h`` release.

Scientific derivatives use only the save/reloaded public Decimal
``velocity(x,y,z,t)``.  No A1 analytic/source-coordinate derivative or
construction tensor is used.  This is scoped leading-stage divergence evidence
only; it cannot promote complete Navier--Stokes or PDE validity.
"""
from __future__ import annotations

import copy
import inspect
import math
from decimal import Decimal, ROUND_HALF_EVEN, localcontext
from typing import Any, Mapping

import numpy as np

from .kokuno_pa16_current_cartesian_postswirl_release2 import (
    KokunoPA16CurrentCartesianPostSwirlRelease2,
)

TASK = "K4-VAL-125"
SCHEMA = "kokuno-a4-postswirl-release2-independent-audit-v1"
UPSTREAM_PR = 1204
UPSTREAM_HEAD = "18e75e6e2df43206147db34788c68ee6a7fa0f37"
SEED = 9173961
TIMES = (0.31, 0.47, 0.63, 0.71)
RELATIVE_FD_LADDER = (2.0e-6, 1.0e-6, 5.0e-7)
CANONICAL_ABSOLUTE_FD_LADDER = (0.02, 0.01, 0.005)
INTERIOR_S_BANDS = ((0.08, 0.24), (0.38, 0.58), (0.72, 0.88))
SEAM_S = (0.015, 0.035)
LATE_S = (0.92, 0.96, 0.985)
AUDIT_DECIMAL_DIGITS = 112
EXPECTED_CANDIDATE_DECIMAL_DIGITS = 96
NORMALIZED_DIVERGENCE_GATE = 1.0e-5
REFINEMENT_FLOOR = 2.0e-8
REFINEMENT_WORSEN_FACTOR = 1.25
FINAL_PROJECT_MOMENTUM_GATE = 1.0e-3
FINAL_PROJECT_DIVERGENCE_GATE = 1.0e-5
FINAL_PROJECT_QUADRATURE = (24, 48, 96)
SPEED_FLOOR = Decimal("1e-300")

_TRUTH_BOUNDARY = {
    "postswirl_release2_divergence_independently_audited": True,
    "canonical_absolute_fd_transfer_ready": False,
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


def default_candidate() -> KokunoPA16CurrentCartesianPostSwirlRelease2:
    return KokunoPA16CurrentCartesianPostSwirlRelease2()


def _ctx(ctx) -> None:
    ctx.prec = AUDIT_DECIMAL_DIGITS
    ctx.rounding = ROUND_HALF_EVEN
    ctx.Emax = 999999999
    ctx.Emin = -999999999


def _vec(value: Any) -> tuple[Decimal, Decimal, Decimal]:
    a = np.asarray(value, dtype=object).reshape(-1)
    if a.size != 3 or any(not isinstance(v, Decimal) or not v.is_finite() for v in a):
        raise TypeError("scientific path requires finite Decimal Cartesian velocity")
    return a[0], a[1], a[2]


def _velocity(candidate: Any, p: Mapping[str, Any]) -> tuple[Decimal, Decimal, Decimal]:
    return _vec(
        candidate.velocity(
            float(p["x"]), float(p["y"]), float(p["z"]), float(p["t"])
        )
    )


def _norm(values: tuple[Decimal, ...] | list[Decimal]) -> Decimal:
    with localcontext() as ctx:
        _ctx(ctx)
        return +sum((v * v for v in values), Decimal(0)).sqrt()


def _log10(value: Decimal) -> float | None:
    if value == 0:
        return None
    with localcontext() as ctx:
        _ctx(ctx)
        return float(abs(value).log10())


def _q_at_z0(candidate: Any, t: float) -> float:
    # Point construction may use the public coordinate map.  Scientific
    # derivatives below use only public velocity evaluations.
    c = candidate.similarity_coordinates_logX(1.0, 0.0, 0.0, float(t))
    q = float(np.asarray(c["q"], dtype=float))
    if not (math.isfinite(q) and q > 0.0):
        raise RuntimeError("invalid public q")
    return q


def _point(
    candidate: Any, s: float, theta: float, t: float, role: str
) -> dict[str, Any]:
    if not (0.0 < s < 1.0):
        raise ValueError("held-out release2 coordinate must be strict-interior")
    q = _q_at_z0(candidate, t)
    lx = float(candidate.log_X_release2_start + s)
    lr = 0.5 * (math.log(2.0) + math.log(q) + lx)
    if lr >= math.log(np.finfo(float).max):
        raise RuntimeError("held-out radius is not finite in binary64")
    r = math.exp(lr)
    x, y = r * math.cos(theta), r * math.sin(theta)
    c = candidate.similarity_coordinates_logX(x, y, 0.0, float(t))
    got_lx = float(np.asarray(c["log_X"], dtype=float))
    got_eta = float(np.asarray(c["eta"], dtype=float))
    if abs(got_lx - lx) > 8.0e-11 or abs(got_eta) > 5.0e-13:
        raise RuntimeError("held-out coordinate reconstruction drift")
    return {
        "x": x,
        "y": y,
        "z": 0.0,
        "t": float(t),
        "r": r,
        "q": q,
        "theta": float(theta),
        "release2_s_target": float(s),
        "role": role,
    }


def _heldout_points(candidate: Any) -> list[dict[str, Any]]:
    rng = np.random.default_rng(SEED)
    out: list[dict[str, Any]] = []
    for t in TIMES:
        for bi, (lo, hi) in enumerate(INTERIOR_S_BANDS, 1):
            for k in range(2):
                out.append(
                    _point(
                        candidate,
                        float(rng.uniform(lo, hi)),
                        float(rng.uniform(0.31, 2.0 * math.pi - 0.31)),
                        t,
                        f"interior_band{bi}_{k}",
                    )
                )
        for k, s in enumerate(SEAM_S):
            out.append(
                _point(
                    candidate,
                    s,
                    float(rng.uniform(0.33, 2.0 * math.pi - 0.33)),
                    t,
                    f"release2_entry_seam_{k}",
                )
            )
        for k, s in enumerate(LATE_S):
            out.append(
                _point(
                    candidate,
                    s,
                    float(rng.uniform(0.33, 2.0 * math.pi - 0.33)),
                    t,
                    f"late_release2_{k}",
                )
            )
    return out


def _canonical_representability(points: list[Mapping[str, Any]]) -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    all_ok = True
    for h in CANONICAL_ABSOLUTE_FD_LADDER:
        flags: list[bool] = []
        for p in points:
            axis_flags: list[bool] = []
            for key in ("x", "y", "z"):
                c = np.float64(float(p[key]))
                hp = np.float64(h)
                ok = bool(c + hp != c and c - hp != c and c + hp != c - hp)
                axis_flags.append(ok)
                all_ok = all_ok and ok
            flags.append(all(axis_flags))
        rows.append(
            {
                "absolute_step": h,
                "all_points_all_axes_representable": bool(all(flags)),
                "representable_point_fraction": float(np.mean(flags)),
            }
        )
    return {"ladder": rows, "all_representable": bool(all_ok)}


def _steps(p: Mapping[str, Any], rel: float) -> tuple[float, float, float]:
    hxy = float(p["r"]) * rel
    hz = max(math.sqrt(float(p["q"])), 1.0e-6) * rel
    if not (math.isfinite(hxy) and hxy > 0.0 and math.isfinite(hz) and hz > 0.0):
        raise RuntimeError("invalid scoped FD step")
    return hxy, hxy, hz


def _fd4(candidate: Any, p: Mapping[str, Any], rel: float):
    cols = []
    with localcontext() as ctx:
        _ctx(ctx)
        for key, h in zip(("x", "y", "z"), _steps(p, rel)):
            c = float(p[key])
            samples = []
            for m in (-2.0, -1.0, 1.0, 2.0):
                q = dict(p)
                q[key] = c + m * h
                samples.append(q)
            coords = [c - 2.0 * h, c - h, c, c + h, c + 2.0 * h]
            if len(set(coords)) != 5:
                raise RuntimeError("scoped FD4 perturbation collapsed")
            fm2, fm1, fp1, fp2 = (_velocity(candidate, q) for q in samples)
            den = Decimal.from_float(12.0 * h)
            col = []
            for i in range(3):
                num = fm2[i] - Decimal(8) * fm1[i] + Decimal(8) * fp1[i] - fp2[i]
                col.append(+(num / den))
            cols.append(tuple(col))
    return tuple(tuple(cols[j][i] for j in range(3)) for i in range(3))


def _family(role: str) -> str:
    if role.startswith("release2_entry"):
        return "seam"
    if role.startswith("late_release2"):
        return "late"
    for name in ("band1", "band2", "band3"):
        if name in role:
            return name
    return "other"


def _weighted_rms(values: np.ndarray, logw: np.ndarray) -> float:
    shift = float(np.max(logw))
    w = np.exp(logw - shift)
    return math.sqrt(float(np.sum(w * values * values)) / float(np.sum(w)))


def _resolution(
    candidate: Any, points: list[Mapping[str, Any]], rel: float
) -> dict[str, Any]:
    divs: list[Decimal] = []
    nds: list[float] = []
    logw: list[float] = []
    witness: list[dict[str, Any]] = []
    for p in points:
        J = _fd4(candidate, p, rel)
        with localcontext() as ctx:
            _ctx(ctx)
            div = +(J[0][0] + J[1][1] + J[2][2])
            grad = _norm([J[i][j] for i in range(3) for j in range(3)])
            speed = _norm(list(_velocity(candidate, p)))
            scale = max(Decimal(1), grad)
            nd = float(+(abs(div) / scale))
        divs.append(abs(div))
        nds.append(nd)
        logw.append(2.0 * math.log(float(p["r"])))
        witness.append(
            {
                "role": p["role"],
                "family": _family(str(p["role"])),
                "t": p["t"],
                "x": p["x"],
                "y": p["y"],
                "z": p["z"],
                "r": p["r"],
                "release2_s_target": p["release2_s_target"],
                "abs_divergence_decimal": str(abs(div)),
                "abs_divergence_log10": _log10(div),
                "normalized_divergence": nd,
                "gradient_frobenius_decimal": str(grad),
                "speed_decimal": str(speed),
            }
        )
    a = np.asarray(nds, dtype=float)
    lw = np.asarray(logw, dtype=float)
    fam: dict[str, float] = {}
    for name in ("band1", "band2", "band3", "seam", "late"):
        vals = [w["normalized_divergence"] for w in witness if w["family"] == name]
        fam[name] = float(max(vals)) if vals else 0.0
    with localcontext() as ctx:
        _ctx(ctx)
        rms = +(sum((v * v for v in divs), Decimal(0)) / Decimal(len(divs))).sqrt()
    return {
        "relative_step": rel,
        "sample_count": len(points),
        "divergence_sampled_max_decimal": str(max(divs)),
        "divergence_sampled_max_log10": _log10(max(divs)),
        "divergence_sample_rms_decimal": str(rms),
        "divergence_sample_rms_log10": _log10(rms),
        "normalized_divergence_sampled_max": float(np.max(a)),
        "normalized_divergence_sample_rms": float(np.sqrt(np.mean(a * a))),
        "logradial_measure_normalized_l2_proxy": _weighted_rms(a, lw),
        "family_normalized_divergence_max": fam,
        "worst_witness": witness[int(np.argmax(a))],
    }


def _seam(candidate: Any) -> dict[str, Any]:
    rng = np.random.default_rng(SEED + 17)
    rows: list[dict[str, Any]] = []
    exact: list[bool] = []
    for t in TIMES:
        q = _q_at_z0(candidate, t)
        r = math.exp(
            0.5 * (math.log(2.0) + math.log(q) + candidate.log_X_release2_start)
        )
        th = float(rng.uniform(0.29, 2.0 * math.pi - 0.29))
        x, y = r * math.cos(th), r * math.sin(th)
        c = _vec(candidate.velocity(x, y, 0.0, t))
        parent = _vec(candidate.parent.velocity(x, y, 0.0, t))
        same = c == parent
        exact.append(same)
        rows.append({"t": t, "r": r, "theta": th, "exact": same})
    return {"all_exact": bool(all(exact)), "rows": rows}


def _stage_axis(candidate: Any) -> dict[str, Any]:
    axis = []
    for t in TIMES:
        v = _vec(candidate.velocity(0.0, 0.0, 0.0, t))
        axis.append(v[0] == 0 and v[1] == 0)
    before = after = False
    try:
        candidate.profile_logX(candidate.log_X_release2_start - 1.0e-3, 0.0)
    except ValueError:
        before = True
    try:
        candidate.profile_logX(candidate.log_X_release2_end + 1.0e-3, 0.0)
    except ValueError:
        after = True
    return {
        "exact_axis_transverse_velocity_zero": bool(all(axis)),
        "before_release2_profile_fails_closed": before,
        "after_release2_profile_fails_closed": after,
    }


def _assess(report: Mapping[str, Any]) -> list[str]:
    fail: list[str] = []
    fine = report["resolutions"][-1]
    med = report["resolutions"][-2]
    for key in (
        "normalized_divergence_sampled_max",
        "normalized_divergence_sample_rms",
        "logradial_measure_normalized_l2_proxy",
    ):
        if float(fine[key]) > NORMALIZED_DIVERGENCE_GATE:
            fail.append("fine_" + key)
    for k, v in fine["family_normalized_divergence_max"].items():
        if float(v) > NORMALIZED_DIVERGENCE_GATE:
            fail.append(f"fine_{k}_divergence")
    if not report["save_load_replay_exact"]:
        fail.append("save_load_replay")
    if not report["release2_entry_seam_replay"]["all_exact"]:
        fail.append("release2_entry_seam")
    fw = report["axis_and_stage_firewall"]
    if not fw["exact_axis_transverse_velocity_zero"]:
        fail.append("axis_regular")
    if not fw["before_release2_profile_fails_closed"] or not fw["after_release2_profile_fails_closed"]:
        fail.append("stage_fail_close")
    if not report["nontriviality_passed"]:
        fail.append("velocity_collapsed")
    f = float(fine["normalized_divergence_sampled_max"])
    m = float(med["normalized_divergence_sampled_max"])
    if f > REFINEMENT_FLOOR and f > REFINEMENT_WORSEN_FACTOR * max(m, 1.0e-300):
        fail.append("medium_to_fine_worsened")
    return fail


def audit_loaded_release2(
    loaded: Any, pre_serialization_reference: Any
) -> dict[str, Any]:
    """Run frozen K4-VAL-125; no caller scientific tuning knobs."""
    if not isinstance(loaded, KokunoPA16CurrentCartesianPostSwirlRelease2):
        raise TypeError("loaded candidate must be exact #1204 release2 type")
    if not isinstance(
        pre_serialization_reference, KokunoPA16CurrentCartesianPostSwirlRelease2
    ):
        raise TypeError("reference must be exact #1204 release2 type")
    cfg = loaded.configuration()
    if int(cfg.get("decimal_digits", -1)) != EXPECTED_CANDIDATE_DECIMAL_DIGITS:
        raise ValueError("candidate precision provenance drift")
    if not bool(loaded.truth_boundary.get("source_l_minus1_to_minus_h_transition_materialized", False)):
        raise ValueError("exact #1204 -1 -> -h transition is required")
    if not bool(loaded.truth_boundary.get("current_cartesian_l_minus1_to_minus_h_composed", False)):
        raise ValueError("exact #1204 Cartesian release2 composition is required")
    for key in (
        "source_terminal_multiplier_materialized",
        "source_exterior_heat_replacement_materialized",
        "outer_global_leading_velocity_materialized",
        "unified_global_cartesian_velocity_export_ready",
        "matched_global_pressure_materialized",
        "restricted_forcing_materialized",
        "heldout_ns_residual_assessed",
        "pde_validated",
    ):
        if bool(loaded.truth_boundary.get(key, False)):
            raise ValueError("upstream truth promotion is inadmissible: " + key)

    points = _heldout_points(loaded)
    canonical = _canonical_representability(points)
    resolutions = [_resolution(loaded, points, h) for h in RELATIVE_FD_LADDER]
    replay = all(
        _velocity(loaded, p) == _velocity(pre_serialization_reference, p)
        for p in points
    )
    speeds = [_norm(list(_velocity(loaded, p))) for p in points]
    late_speeds = [
        speed for speed, p in zip(speeds, points) if str(p["role"]).startswith("late_release2")
    ]
    report: dict[str, Any] = {
        "schema": SCHEMA,
        "task": TASK,
        "upstream_pr": UPSTREAM_PR,
        "upstream_head": UPSTREAM_HEAD,
        "candidate_semantic_sha256": loaded.semantic_sha256,
        "seed": SEED,
        "times": list(TIMES),
        "release2_length": float(loaded.log_X_release2_end - loaded.log_X_release2_start),
        "h_current": float(loaded.h_value),
        "interior_s_bands": [list(x) for x in INTERIOR_S_BANDS],
        "seam_s": list(SEAM_S),
        "late_s": list(LATE_S),
        "heldout_points": points,
        "independent_operator": "centered_cartesian_five_point_fd4_on_loaded_public_decimal_velocity",
        "normalization": "abs(div_u)/max(1, frobenius_norm(Jac_u))",
        "scientific_derivative_uses_profile_derivatives": False,
        "relative_fd_ladder": list(RELATIVE_FD_LADDER),
        "canonical_absolute_fd_ladder": list(CANONICAL_ABSOLUTE_FD_LADDER),
        "resolutions": resolutions,
        "canonical_absolute_fd_representability": canonical,
        "save_load_replay_exact": bool(replay),
        "release2_entry_seam_replay": _seam(loaded),
        "axis_and_stage_firewall": _stage_axis(loaded),
        "nontriviality_passed": bool(all(s > SPEED_FLOOR for s in speeds)),
        "late_release2_velocity_finite_nonzero": bool(
            late_speeds and all(s > SPEED_FLOOR for s in late_speeds)
        ),
        "late_release2_speed_min_decimal": str(min(late_speeds)) if late_speeds else None,
        "gates": {
            "scoped_normalized_divergence": NORMALIZED_DIVERGENCE_GATE,
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
    report["scoped_postswirl_release2_divergence_gate_passed"] = len(failures) == 0
    report["truth_boundary"]["canonical_absolute_fd_transfer_ready"] = bool(
        canonical["all_representable"] and not failures
    )
    report["final_project_admission_ready"] = False
    report["final_project_blockers"] = [
        "terminal_multiplier_missing",
        "exterior_heat_replacement_missing",
        "global_leading_missing",
        "matching_self_contained_oscillatory_composite_missing",
        "cartesian_correction_velocity_missing",
        "matched_pressure_missing",
        "preregistered_restricted_forcing_missing",
        "complete_ns_residual_not_assessed",
    ]
    if not canonical["all_representable"]:
        report["final_project_blockers"].insert(
            0, "canonical_absolute_fd_ladder_not_binary64_representable_here"
        )
    return report


def enforce_preregistered_gates(report: Mapping[str, Any]) -> None:
    failures = _assess(report)
    if failures:
        raise AssertionError("K4-VAL-125 scoped gate failure: " + ", ".join(failures))
    if report.get("truth_boundary", {}).get("pde_validated") is not False:
        raise AssertionError("pde_validated must remain false")
    if report.get("final_project_admission_ready") is not False:
        raise AssertionError("scoped release2 audit cannot promote final admission")


def public_api_has_no_scientific_tuning_knobs() -> bool:
    params = set(inspect.signature(audit_loaded_release2).parameters)
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
