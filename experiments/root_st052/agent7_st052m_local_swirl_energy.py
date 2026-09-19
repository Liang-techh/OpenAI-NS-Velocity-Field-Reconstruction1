"""Screen one local swirl-energy compensation for the ST052-M tip taper.

Preregistered in issue #557 before evaluation.  Keep the #544 localized radial
Piola tip taper at fixed tau=.05 but remove its global energy-restoring scale.
Add exactly one off-midplane, axisymmetric swirl-only correction supported in
the frozen shoulder window.  Its coefficient beta is determined only by the
reference-energy equation E_child(t=.25)=E_control(t=.25) on order-64
quadrature.  No visual/path/PDE fit and no post-transform common scale are
allowed.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

import agent7_st052m_energy_neutral_taper as failed
import agent7_st052m_frozen_redistribution as base
import agent7_st052m_localized_taper as prior
import replay_st052

TASK_ID = "CR003-ST052M-LOCAL-SWIRL-ENERGY-082"
PREREG_ISSUE = 557
PARENT_ID = "ST052-M+frozen-redistribution-.05"
PARENT_HEAD = "779ffca71066e2864496d37de55a7aafc45d6f57"
SOURCE_TAPER_HEAD = "093c7171cd61c6bd439afa30b2da69598a02d182"
SOURCE_FAILED_COMP_HEAD = "62e8c170427d5d830d7f897ba31768e0fc4ce56a"
TAPER_TAU = 0.05
SHOULDER_WINDOW = (0.36, 0.49)
BETA_BRACKET = (0.0, 0.9)
ENERGY_ORDER = 64
ENERGY_ROOT_TOL = 1.0e-11
ENERGY_ROOT_MAXITER = 80
TIMES = base.TIMES
DIVERGENCE_FD_STEP = prior.DIVERGENCE_FD_STEP

CRITERIA = {
    "energy_relative_error_max": 1.0e-8,
    "central_velocity_identity_abs_max": 1.0e-11,
    "central_angular_retention_abs_error_max": 1.0e-8,
    "tip_radial_rms_reduction_min": 0.01,
    "central_radial_rms_abs_change_max": 0.0025,
    "global_axial_rms_loss_max": 0.01,
    "tip_absolute_enstrophy_ratio_min": 0.90,
    "collar_absolute_enstrophy_ratio_max": 1.15,
    "support_max_abs": 1.0e-12,
    "divergence_fd_max": 1.0e-5,
    "response_rank_min": 4,
    "response_condition_max": 8.0,
    "response_abs_cosine_max": 0.90,
}

TRUTH = {
    "canonical_velocity_changed": False,
    "saved_velocity_changed": False,
    "production_taper_selected": False,
    "production_compensation_selected": False,
    "production_candidate_selected": False,
    "pressure_or_force_changed": False,
    "held_out_pde_residual_evaluated": False,
    "parent_pde_receipt_transferred": False,
    "material_paths_integrated": False,
    "public_image_used": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "source_correspondence_verified": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def shoulder_profile(z):
    """Frozen symmetric C-infinity shoulder bump q_s(z)."""
    return failed._bump_and_derivative(z, SHOULDER_WINDOW)


def raw_taper_velocity(f, raw, points, time, redistribution_scale):
    """#544 Piola child before its global energy-restoring scale."""
    return prior.taper_velocity(
        f,
        raw,
        points,
        time,
        redistribution_scale=redistribution_scale,
        tau=TAPER_TAU,
        taper_scale=1.0,
    )


def apply_swirl_compensation(points, velocity, beta):
    """Apply u_theta -> (1-beta*q_s(z))*u_theta in Cartesian coordinates.

    The correction is axisymmetric and theta-directed only.  Thus, for the
    axisymmetric parent used here, it is divergence-free independently of the
    poloidal components.  The shoulder bump vanishes in the registered center
    and outside its compact z window.
    """
    points = np.asarray(points, dtype=float)
    velocity = np.asarray(velocity, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError("points must have shape (n,3)")
    if velocity.shape != points.shape:
        raise ValueError("velocity must match points shape")
    beta = float(beta)
    if not np.isfinite(beta) or not (BETA_BRACKET[0] <= beta <= BETA_BRACKET[1]):
        raise ValueError("beta outside frozen bracket")

    q, _ = shoulder_profile(points[:, 2])
    x = points[:, 0]
    y = points[:, 1]
    r = np.hypot(x, y)
    out = velocity.copy()
    mask = r > 1.0e-14
    if np.any(mask):
        xm = x[mask]
        ym = y[mask]
        rm = r[mask]
        u = velocity[mask]
        utheta = (-ym * u[:, 0] + xm * u[:, 1]) / rm
        dtheta = -beta * q[mask] * utheta
        out[mask, 0] += dtheta * (-ym / rm)
        out[mask, 1] += dtheta * (xm / rm)
    return out


def compensated_velocity(
    f,
    raw,
    points,
    time,
    *,
    redistribution_scale,
    beta,
):
    points = np.asarray(points, dtype=float)
    u = raw_taper_velocity(f, raw, points, time, redistribution_scale)
    return apply_swirl_compensation(points, u, beta)


def solve_energy_beta(f, raw, redistribution_scale):
    """Bisection on the frozen beta bracket using reference energy only."""
    control = lambda p, t: prior.control_velocity(f, raw, p, t, redistribution_scale)
    target = prior.axisymmetric_energy(control, order=ENERGY_ORDER)
    cache = {}

    def residual(beta):
        key = float(beta)
        if key not in cache:
            child = lambda p, t: compensated_velocity(
                f,
                raw,
                p,
                t,
                redistribution_scale=redistribution_scale,
                beta=key,
            )
            cache[key] = prior.axisymmetric_energy(child, order=ENERGY_ORDER)
        return cache[key] - target

    lo, hi = BETA_BRACKET
    flo = residual(lo)
    fhi = residual(hi)
    result = {
        "root_exists": False,
        "beta": None,
        "control_energy": target,
        "bracket": [lo, hi],
        "endpoint_energy_residuals": [flo, fhi],
        "iterations": 0,
        "post_transform_common_scale": 1.0,
    }
    if flo == 0.0:
        root = lo
    elif fhi == 0.0:
        root = hi
    elif flo * fhi > 0.0:
        return result
    else:
        a, b = lo, hi
        fa = flo
        root = 0.5 * (a + b)
        for it in range(1, ENERGY_ROOT_MAXITER + 1):
            root = 0.5 * (a + b)
            fm = residual(root)
            result["iterations"] = it
            if abs(fm) <= ENERGY_ROOT_TOL * max(1.0, target) or (b - a) <= 1.0e-12:
                break
            if fa * fm <= 0.0:
                b = root
            else:
                a, fa = root, fm
    child_energy = target + residual(root)
    result.update(
        {
            "root_exists": True,
            "beta": float(root),
            "child_energy": float(child_energy),
            "relative_energy_error": float(child_energy / target - 1.0),
            "evaluations": len(cache),
        }
    )
    return result


def morphology_rows(f, raw, redistribution_scale, beta):
    control = lambda p, t: prior.control_velocity(f, raw, p, t, redistribution_scale)
    child = lambda p, t: compensated_velocity(
        f,
        raw,
        p,
        t,
        redistribution_scale=redistribution_scale,
        beta=beta,
    )
    rows = []
    for t in TIMES:
        b = prior.morphology(control, t)
        c = prior.morphology(child, t)
        rows.append(
            {
                "time": float(t),
                "control": b,
                "child": c,
                "tip_radial_rms_relative": float(c["tip_radial_rms"] / b["tip_radial_rms"] - 1.0),
                "central_radial_rms_relative": float(c["central_radial_rms"] / b["central_radial_rms"] - 1.0),
                "global_axial_rms_relative": float(c["axial_rms"] / b["axial_rms"] - 1.0),
                "global_radial_rms_relative": float(c["radial_rms"] / b["radial_rms"] - 1.0),
                "tip_absolute_enstrophy_ratio": float(c["tip_absolute_enstrophy"] / b["tip_absolute_enstrophy"]),
                "central_absolute_enstrophy_ratio": float(c["central_absolute_enstrophy"] / b["central_absolute_enstrophy"]),
                "collar_absolute_ratio": float(c["collar_absolute_enstrophy"] / b["collar_absolute_enstrophy"]),
            }
        )
    return rows


def central_identity(f, raw, redistribution_scale, beta, seed=9175821):
    rng = np.random.default_rng(seed)
    n = 48
    r = rng.uniform(0.05, 1.45, size=n)
    th = rng.uniform(0.0, 2.0 * np.pi, size=n)
    z = rng.uniform(-0.68, 0.68, size=n)
    pts = np.column_stack((r * np.cos(th), r * np.sin(th), z))
    max_abs = 0.0
    max_rel = 0.0
    for t in TIMES:
        u0 = prior.control_velocity(f, raw, pts, t, redistribution_scale)
        u1 = compensated_velocity(
            f,
            raw,
            pts,
            t,
            redistribution_scale=redistribution_scale,
            beta=beta,
        )
        diff = np.abs(u1 - u0)
        max_abs = max(max_abs, float(np.max(diff)))
        denom = max(float(np.max(np.abs(u0))), 1.0e-15)
        max_rel = max(max_rel, float(np.max(diff) / denom))
    return {"max_abs": max_abs, "max_rel": max_rel}


def central_angular_retention(f, raw, redistribution_scale, beta):
    zs = np.array([-0.30, 0.0, 0.30])
    out = {}
    for r in (0.6, 0.9):
        pts = np.column_stack((np.full(len(zs), r), np.zeros(len(zs)), zs))
        u0 = prior.control_velocity(f, raw, pts, 0.50, redistribution_scale)
        u1 = compensated_velocity(
            f,
            raw,
            pts,
            0.50,
            redistribution_scale=redistribution_scale,
            beta=beta,
        )
        a0 = float(np.mean(np.abs(u0[:, 1]) / r))
        a1 = float(np.mean(np.abs(u1[:, 1]) / r))
        out[str(r)] = float(a1 / a0)
    return out


def structure_preflight(f, raw, redistribution_scale, beta, seed=9175822):
    child = lambda p, t: compensated_velocity(
        f,
        raw,
        p,
        t,
        redistribution_scale=redistribution_scale,
        beta=beta,
    )
    core = np.array(
        [[r, 0.0, z] for r in (0.10, 0.30, 0.60) for z in (-0.50, -0.20, 0.20, 0.50)],
        float,
    )
    uc = child(core, 0.50)
    signs = bool(
        np.all(uc[:, 0] < 0.0)
        and np.all(uc[:, 1] > 0.0)
        and np.all(core[:, 2] * uc[:, 2] > 0.0)
    )
    outside = np.array(
        [
            [2.05, 0, 0],
            [-2.05, 0, 0],
            [0, 2.05, 0],
            [0, -2.05, 0],
            [0, 0, 2.05],
            [0, 0, -2.05],
            [1.8, 0, 2.05],
        ],
        float,
    )
    support_max = max(float(np.max(np.abs(child(outside, t)))) for t in TIMES)

    rng = np.random.default_rng(seed)
    pts = rng.uniform(-1.35, 1.35, size=(36, 3))
    div = np.zeros(len(pts))
    for axis in range(3):
        d = np.zeros(3)
        d[axis] = DIVERGENCE_FD_STEP
        up = child(pts + d, 0.50)
        um = child(pts - d, 0.50)
        div += (up[:, axis] - um[:, axis]) / (2.0 * DIVERGENCE_FD_STEP)

    z = np.linspace(-2.0, 2.0, 16001)
    q, _ = shoulder_profile(z)
    factor = 1.0 - beta * q
    return {
        "core_signs_pass": signs,
        "support_max_abs": support_max,
        "divergence_fd_max": float(np.max(np.abs(div))),
        "swirl_factor_min": float(np.min(factor)),
        "swirl_factor_max": float(np.max(factor)),
        "central_compensation_identity_error": float(
            np.max(np.abs(factor[np.abs(z) / failed.AXIAL_SUPPORT <= 0.35] - 1.0))
        ),
    }


def response_diagnostics(f51, raw51, f52, raw52, redistribution_scale, seed=9175823):
    rng = np.random.default_rng(seed)
    pts = []
    ts = []
    for t in TIMES:
        for _ in range(12):
            r = rng.uniform(0.18, 1.45)
            th = rng.uniform(0.0, 2.0 * np.pi)
            z = rng.uniform(-1.35, 1.35)
            pts.append([r * np.cos(th), r * np.sin(th), z])
            ts.append(float(t))
    pts = np.asarray(pts, float)
    ts = np.asarray(ts, float)

    u51 = base.base_velocity_many(f51, raw51, pts, ts)
    u52 = base.base_velocity_many(f52, raw52, pts, ts)
    ured = base.transformed_velocity_many(
        f52, raw52, pts, ts, gain=base.GAIN, scale=redistribution_scale
    )

    eps = 1.0e-4
    plus = []
    minus = []
    taper0 = []
    comp_eps = []
    for t in TIMES:
        mask = ts == t
        p = pts[mask]
        plus.append(
            prior.taper_velocity(
                f52,
                raw52,
                p,
                t,
                redistribution_scale=redistribution_scale,
                tau=eps,
                taper_scale=1.0,
            ).ravel()
        )
        minus.append(
            prior.taper_velocity(
                f52,
                raw52,
                p,
                t,
                redistribution_scale=redistribution_scale,
                tau=-eps,
                taper_scale=1.0,
            ).ravel()
        )
        u0 = raw_taper_velocity(f52, raw52, p, t, redistribution_scale)
        taper0.append(u0.ravel())
        comp_eps.append(apply_swirl_compensation(p, u0, eps).ravel())
    taper_col = (np.concatenate(plus) - np.concatenate(minus)) / (2.0 * eps)
    compensation_col = (np.concatenate(comp_eps) - np.concatenate(taper0)) / eps

    cols = [
        (u52 - u51).ravel(),
        (ured - u52).ravel(),
        taper_col,
        compensation_col,
    ]
    norms = [float(np.linalg.norm(c)) for c in cols]
    if min(norms) <= 1.0e-14:
        return {
            "rank": 0,
            "condition_number": float("inf"),
            "max_pairwise_abs_cosine": 1.0,
            "column_norms": norms,
        }
    A = np.column_stack([c / n for c, n in zip(cols, norms)])
    s = np.linalg.svd(A, compute_uv=False)
    cosines = []
    for i in range(4):
        for j in range(i + 1, 4):
            cosines.append(float(np.dot(A[:, i], A[:, j])))
    return {
        "rank": int(np.linalg.matrix_rank(A, tol=1.0e-10)),
        "singular_values": [float(v) for v in s],
        "condition_number": float(s[0] / s[-1]),
        "pairwise_cosines": cosines,
        "max_pairwise_abs_cosine": float(max(abs(v) for v in cosines)),
        "column_norms": norms,
    }


def clean_rule(row):
    if not row["energy_solve"]["root_exists"]:
        return False
    c = CRITERIA
    return bool(
        abs(row["energy_solve"]["relative_energy_error"]) <= c["energy_relative_error_max"]
        and row["central_velocity_identity"]["max_abs"] <= c["central_velocity_identity_abs_max"]
        and max(abs(v - 1.0) for v in row["central_angular_retention"].values())
        <= c["central_angular_retention_abs_error_max"]
        and max(m["tip_radial_rms_relative"] for m in row["morphology"])
        <= -c["tip_radial_rms_reduction_min"]
        and max(abs(m["central_radial_rms_relative"]) for m in row["morphology"])
        <= c["central_radial_rms_abs_change_max"]
        and min(m["global_axial_rms_relative"] for m in row["morphology"])
        >= -c["global_axial_rms_loss_max"]
        and min(m["tip_absolute_enstrophy_ratio"] for m in row["morphology"])
        >= c["tip_absolute_enstrophy_ratio_min"]
        and max(m["collar_absolute_ratio"] for m in row["morphology"])
        <= c["collar_absolute_enstrophy_ratio_max"]
        and row["structure"]["core_signs_pass"]
        and row["structure"]["support_max_abs"] <= c["support_max_abs"]
        and row["structure"]["divergence_fd_max"] <= c["divergence_fd_max"]
        and row["response"]["rank"] >= c["response_rank_min"]
        and row["response"]["condition_number"] <= c["response_condition_max"]
        and row["response"]["max_pairwise_abs_cosine"] <= c["response_abs_cosine_max"]
    )


def run(out: Path):
    f52, raw52 = replay_st052.reconstruct()
    f51, raw51 = replay_st052.parent_field()
    redistribution_scale, parent_energy, red_unscaled_energy = base.child_scale(f52, raw52)
    energy_solve = solve_energy_beta(f52, raw52, redistribution_scale)

    row = {
        "redistribution_scale": redistribution_scale,
        "parent_energy": parent_energy,
        "redistribution_unscaled_energy": red_unscaled_energy,
        "energy_solve": energy_solve,
        "beta": energy_solve.get("beta"),
        "morphology": [],
        "central_velocity_identity": None,
        "central_angular_retention": None,
        "structure": None,
        "response": None,
        "clean_local_swirl_energy_capacity": False,
    }
    if energy_solve["root_exists"]:
        beta = float(energy_solve["beta"])
        row["morphology"] = morphology_rows(f52, raw52, redistribution_scale, beta)
        row["central_velocity_identity"] = central_identity(
            f52, raw52, redistribution_scale, beta
        )
        row["central_angular_retention"] = central_angular_retention(
            f52, raw52, redistribution_scale, beta
        )
        row["structure"] = structure_preflight(f52, raw52, redistribution_scale, beta)
        row["response"] = response_diagnostics(
            f51, raw51, f52, raw52, redistribution_scale
        )
        row["clean_local_swirl_energy_capacity"] = clean_rule(row)

    report = {
        "task_id": TASK_ID,
        "prereg_issue": PREREG_ISSUE,
        "parent_id": PARENT_ID,
        "parent_head": PARENT_HEAD,
        "source_taper_head": SOURCE_TAPER_HEAD,
        "source_failed_compensation_head": SOURCE_FAILED_COMP_HEAD,
        "frozen_transform": {
            "tau": TAPER_TAU,
            "tip_window_abs_z_over_2": list(prior.WINDOW),
            "shoulder_window_abs_z_over_2": list(SHOULDER_WINDOW),
            "beta_bracket": list(BETA_BRACKET),
            "redistribution_alpha": base.SOURCE_ALPHA,
            "redistribution_gain": base.GAIN,
            "redistribution_inner_window": list(base.INNER_WINDOW),
            "redistribution_outer_window": list(base.OUTER_WINDOW),
        },
        "compensation_representation": "axisymmetric_swirl_only_compact_shoulder",
        "parameter_scan_performed": False,
        "energy_only_dependent_root_solve": True,
        "post_transform_common_scale_allowed": False,
        "diagnostic_added_velocity_degrees": 1,
        "independently_tuned_new_coefficients": 0,
        "criteria": CRITERIA,
        "row": row,
        **TRUTH,
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    r = run(args.out)
    q = r["row"]
    e = q["energy_solve"]
    print("energy_root_exists=", e["root_exists"])
    print("endpoint_energy_residuals=", e["endpoint_energy_residuals"])
    print("beta=", q["beta"])
    print("clean_local_swirl_energy_capacity=", q["clean_local_swirl_energy_capacity"])
    if e["root_exists"]:
        print("energy_relative_error=", e["relative_energy_error"])
        print("central_velocity_identity=", q["central_velocity_identity"])
        for m in q["morphology"]:
            print(
                "t=", m["time"],
                "tip_rr=", m["tip_radial_rms_relative"],
                "central_rr=", m["central_radial_rms_relative"],
                "axial=", m["global_axial_rms_relative"],
                "tip_E=", m["tip_absolute_enstrophy_ratio"],
            )
        print("response=", q["response"])


if __name__ == "__main__":
    main()
