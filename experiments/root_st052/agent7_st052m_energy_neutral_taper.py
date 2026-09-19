"""Screen one energy-neutral localized Piola taper on ST052-M + frozen redistribution.

Preregistered in issue #548 before evaluation.  The #544 tip-localized Piola map
is retained, but one disjoint off-midplane shoulder bump is added to the same
radial map.  Its coefficient is *not* a visual/PDE/path fit: it is determined
only by E_child(t=.25) = E_control(t=.25) on the frozen order-64 quadrature.
No post-Piola common scale is permitted, so the central velocity can remain
pointwise identical where both bump profiles vanish.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

import agent7_st052m_frozen_redistribution as base
import agent7_st052m_localized_taper as prior
import replay_st052

TASK_ID = "CR003-ST052M-ENERGY-NEUTRAL-TAPER-081"
PREREG_ISSUE = 548
PARENT_ID = "ST052-M+frozen-redistribution-.05"
PARENT_HEAD = "779ffca71066e2864496d37de55a7aafc45d6f57"
SOURCE_TAPER_HEAD = "093c7171cd61c6bd439afa30b2da69598a02d182"
TAPER_TAU = 0.05
AXIAL_SUPPORT = 2.0
TIP_WINDOW = (0.50, 0.82)
SHOULDER_WINDOW = (0.36, 0.49)
ENERGY_ALPHA_BRACKET = (0.0, 4.0)
ENERGY_ORDER = 64
ENERGY_ROOT_TOL = 1.0e-11
ENERGY_ROOT_MAXITER = 80
TIMES = base.TIMES
MORPH_NR = prior.MORPH_NR
MORPH_NZ = prior.MORPH_NZ
VORTICITY_FD_STEP = prior.VORTICITY_FD_STEP
DIVERGENCE_FD_STEP = prior.DIVERGENCE_FD_STEP

CRITERIA = {
    "energy_relative_error_max": 1.0e-8,
    "tip_radial_rms_reduction_min": 0.01,
    "central_radial_rms_abs_change_max": 0.0025,
    "central_velocity_identity_abs_max": 1.0e-11,
    "global_axial_rms_loss_max": 0.01,
    "tip_absolute_enstrophy_ratio_min": 0.90,
    "shoulder_radial_rms_abs_change_max": 0.04,
    "shoulder_absolute_enstrophy_ratio_max": 1.25,
    "collar_absolute_enstrophy_ratio_max": 1.15,
    "central_angular_retention_abs_error_max": 1.0e-8,
    "support_max_abs": 1.0e-12,
    "divergence_fd_max": 1.0e-5,
    "determinant_min": 0.0,
    "response_rank_min": 3,
    "response_condition_max": 5.0,
    "response_abs_cosine_max": 0.90,
}

TRUTH = {
    "canonical_velocity_changed": False,
    "saved_velocity_changed": False,
    "production_taper_selected": False,
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


def _bump_and_derivative(z, window):
    """Symmetric C-infinity bump in window=(lo,hi) of |z|/2."""
    z = np.asarray(z, dtype=float)
    if not np.all(np.isfinite(z)):
        raise ValueError("z must be finite")
    lo, hi = map(float, window)
    if not (0.0 < lo < hi < 1.0):
        raise ValueError("window must lie strictly inside (0,1)")
    s = np.abs(z) / AXIAL_SUPPORT
    mid = 0.5 * (lo + hi)
    half = 0.5 * (hi - lo)
    x = (s - mid) / half
    q = np.zeros_like(s)
    dq = np.zeros_like(s)
    mask = np.abs(x) < 1.0
    if np.any(mask):
        xm = x[mask]
        den = 1.0 - xm * xm
        qm = np.exp(1.0 - 1.0 / den)
        q[mask] = qm
        dqdx = qm * (-2.0 * xm) / (den * den)
        dsdz = np.sign(z[mask]) / AXIAL_SUPPORT
        dq[mask] = dqdx * (1.0 / half) * dsdz
    return q, dq


def neutral_profile_and_derivative(z, alpha_e):
    alpha_e = float(alpha_e)
    if not np.isfinite(alpha_e) or not (ENERGY_ALPHA_BRACKET[0] <= alpha_e <= ENERGY_ALPHA_BRACKET[1]):
        raise ValueError("energy compensation alpha outside frozen bracket")
    qt, dqt = _bump_and_derivative(z, TIP_WINDOW)
    qs, dqs = _bump_and_derivative(z, SHOULDER_WINDOW)
    return qt - alpha_e * qs, dqt - alpha_e * dqs


def neutral_map(points, alpha_e, tau=TAPER_TAU):
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError("points must have shape (n,3)")
    tau = float(tau)
    if not np.isfinite(tau) or abs(tau) > 0.15:
        raise ValueError("tau outside audited local range")
    q, dq = neutral_profile_and_derivative(points[:, 2], alpha_e)
    a = 1.0 + tau * q
    ap = tau * dq
    if np.any(a <= 0.0):
        raise RuntimeError("energy-neutral taper lost orientation")
    mapped = points.copy()
    mapped[:, 0] *= a
    mapped[:, 1] *= a
    return mapped, a, ap


def control_velocity(f, raw, points, time, redistribution_scale):
    return base.transformed_velocity(
        f, raw, points, time, gain=base.GAIN, scale=redistribution_scale
    )


def neutral_velocity(
    f,
    raw,
    points,
    time,
    *,
    redistribution_scale,
    alpha_e,
    tau=TAPER_TAU,
):
    """Exact contravariant Piola pullback; intentionally no final common scale."""
    points = np.asarray(points, dtype=float)
    mapped, a, ap = neutral_map(points, alpha_e, tau=tau)
    u = control_velocity(f, raw, mapped, time, redistribution_scale)
    x = points[:, 0]
    y = points[:, 1]
    uz = u[:, 2].copy()
    out = np.empty_like(u)
    out[:, 0] = a * u[:, 0] - a * ap * x * uz
    out[:, 1] = a * u[:, 1] - a * ap * y * uz
    out[:, 2] = a * a * uz
    return out


def axisymmetric_energy(velocity_fn, order=ENERGY_ORDER):
    xr, wr = leggauss(int(order))
    xz, wz = leggauss(int(order))
    r = xr + 1.0
    z = 2.0 * xz
    R, Z = np.meshgrid(r, z, indexing="ij")
    pts = np.column_stack((R.ravel(), np.zeros(R.size), Z.ravel()))
    u = np.asarray(velocity_fn(pts, 0.25), float).reshape(R.shape + (3,))
    weight = wr[:, None] * (2.0 * wz[None, :]) * np.pi * R
    return float(np.sum(weight * np.sum(u * u, axis=-1)))


def solve_energy_alpha(f, raw, redistribution_scale):
    """Bisection on the frozen alpha bracket using energy equality only."""
    control = lambda p, t: control_velocity(f, raw, p, t, redistribution_scale)
    target = axisymmetric_energy(control, order=ENERGY_ORDER)

    cache = {}

    def residual(alpha):
        key = float(alpha)
        if key not in cache:
            child = lambda p, t: neutral_velocity(
                f, raw, p, t,
                redistribution_scale=redistribution_scale,
                alpha_e=key,
            )
            cache[key] = axisymmetric_energy(child, order=ENERGY_ORDER)
        return cache[key] - target

    lo, hi = ENERGY_ALPHA_BRACKET
    flo = residual(lo)
    fhi = residual(hi)
    result = {
        "root_exists": False,
        "alpha_e": None,
        "control_energy": target,
        "bracket": [lo, hi],
        "endpoint_energy_residuals": [flo, fhi],
        "iterations": 0,
        "post_piola_common_scale": 1.0,
    }
    if flo == 0.0:
        root = lo
    elif fhi == 0.0:
        root = hi
    elif flo * fhi > 0.0:
        return result
    else:
        a, b = lo, hi
        fa, fb = flo, fhi
        root = 0.5 * (a + b)
        for it in range(1, ENERGY_ROOT_MAXITER + 1):
            root = 0.5 * (a + b)
            fm = residual(root)
            result["iterations"] = it
            if abs(fm) <= ENERGY_ROOT_TOL * max(1.0, target) or (b - a) <= 1.0e-12:
                break
            if fa * fm <= 0.0:
                b, fb = root, fm
            else:
                a, fa = root, fm
    child_energy = target + residual(root)
    result.update({
        "root_exists": True,
        "alpha_e": float(root),
        "child_energy": float(child_energy),
        "relative_energy_error": float(child_energy / target - 1.0),
        "evaluations": len(cache),
    })
    return result


def morphology(velocity_fn, time):
    xr, wr = leggauss(MORPH_NR)
    xz, wz = leggauss(MORPH_NZ)
    r = xr + 1.0
    z = 2.0 * xz
    R, Z = np.meshgrid(r, z, indexing="ij")
    pts = np.column_stack((R.ravel(), np.zeros(R.size), Z.ravel()))
    omega = base.curl_fd(velocity_fn, pts, time, h=VORTICITY_FD_STEP)
    om2 = np.sum(omega * omega, axis=1).reshape(R.shape)
    weight = 2.0 * np.pi * (wr[:, None] * (2.0 * wz[None, :])) * R
    total = float(np.sum(weight * om2))
    if total <= np.finfo(float).tiny:
        raise RuntimeError("inactive vorticity")
    sabs = np.abs(Z) / AXIAL_SUPPORT
    tip = (sabs >= 0.50) & (sabs <= 0.80)
    central = sabs <= 0.35
    shoulder = (sabs >= SHOULDER_WINDOW[0]) & (sabs <= SHOULDER_WINDOW[1])
    collar = (R >= 1.50) | (np.abs(Z) >= 1.50)

    def band(mask):
        e = float(np.sum(weight * om2 * mask))
        rr = float(np.sqrt(np.sum(weight * R * R * om2 * mask) / e))
        return e, rr

    tip_e, tip_rr = band(tip)
    central_e, central_rr = band(central)
    shoulder_e, shoulder_rr = band(shoulder)
    return {
        "enstrophy": total,
        "axial_rms": float(np.sqrt(np.sum(weight * Z * Z * om2) / total)),
        "radial_rms": float(np.sqrt(np.sum(weight * R * R * om2) / total)),
        "tip_absolute_enstrophy": tip_e,
        "tip_radial_rms": tip_rr,
        "central_absolute_enstrophy": central_e,
        "central_radial_rms": central_rr,
        "shoulder_absolute_enstrophy": shoulder_e,
        "shoulder_radial_rms": shoulder_rr,
        "collar_absolute_enstrophy": float(np.sum(weight * om2 * collar)),
    }


def morphology_rows(f, raw, redistribution_scale, alpha_e):
    control = lambda p, t: control_velocity(f, raw, p, t, redistribution_scale)
    child = lambda p, t: neutral_velocity(
        f, raw, p, t,
        redistribution_scale=redistribution_scale,
        alpha_e=alpha_e,
    )
    rows = []
    for t in TIMES:
        b = morphology(control, t)
        c = morphology(child, t)
        rows.append({
            "time": float(t),
            "control": b,
            "child": c,
            "tip_radial_rms_relative": float(c["tip_radial_rms"] / b["tip_radial_rms"] - 1.0),
            "central_radial_rms_relative": float(c["central_radial_rms"] / b["central_radial_rms"] - 1.0),
            "shoulder_radial_rms_relative": float(c["shoulder_radial_rms"] / b["shoulder_radial_rms"] - 1.0),
            "global_axial_rms_relative": float(c["axial_rms"] / b["axial_rms"] - 1.0),
            "global_radial_rms_relative": float(c["radial_rms"] / b["radial_rms"] - 1.0),
            "tip_absolute_enstrophy_ratio": float(c["tip_absolute_enstrophy"] / b["tip_absolute_enstrophy"]),
            "central_absolute_enstrophy_ratio": float(c["central_absolute_enstrophy"] / b["central_absolute_enstrophy"]),
            "shoulder_absolute_enstrophy_ratio": float(c["shoulder_absolute_enstrophy"] / b["shoulder_absolute_enstrophy"]),
            "collar_absolute_ratio": float(c["collar_absolute_enstrophy"] / b["collar_absolute_enstrophy"]),
        })
    return rows


def central_identity(f, raw, redistribution_scale, alpha_e, seed=9175811):
    rng = np.random.default_rng(seed)
    n = 48
    r = rng.uniform(0.05, 1.45, size=n)
    th = rng.uniform(0.0, 2.0 * np.pi, size=n)
    z = rng.uniform(-0.68, 0.68, size=n)  # strictly inside |z|/2 <= .35
    pts = np.column_stack((r * np.cos(th), r * np.sin(th), z))
    max_abs = 0.0
    max_rel = 0.0
    for t in TIMES:
        u0 = control_velocity(f, raw, pts, t, redistribution_scale)
        u1 = neutral_velocity(
            f, raw, pts, t,
            redistribution_scale=redistribution_scale,
            alpha_e=alpha_e,
        )
        diff = np.abs(u1 - u0)
        max_abs = max(max_abs, float(np.max(diff)))
        denom = max(float(np.max(np.abs(u0))), 1.0e-15)
        max_rel = max(max_rel, float(np.max(diff) / denom))
    return {"max_abs": max_abs, "max_rel": max_rel}


def central_angular_retention(f, raw, redistribution_scale, alpha_e):
    zs = np.array([-0.30, 0.0, 0.30])
    out = {}
    for r in (0.6, 0.9):
        pts = np.column_stack((np.full(len(zs), r), np.zeros(len(zs)), zs))
        u0 = control_velocity(f, raw, pts, 0.50, redistribution_scale)
        u1 = neutral_velocity(
            f, raw, pts, 0.50,
            redistribution_scale=redistribution_scale,
            alpha_e=alpha_e,
        )
        a0 = float(np.mean(np.abs(u0[:, 1]) / r))
        a1 = float(np.mean(np.abs(u1[:, 1]) / r))
        out[str(r)] = float(a1 / a0)
    return out


def structure_preflight(f, raw, redistribution_scale, alpha_e, seed=9175812):
    child = lambda p, t: neutral_velocity(
        f, raw, p, t,
        redistribution_scale=redistribution_scale,
        alpha_e=alpha_e,
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
    outside = np.array([
        [2.05, 0, 0], [-2.05, 0, 0], [0, 2.05, 0], [0, -2.05, 0],
        [0, 0, 2.05], [0, 0, -2.05], [1.8, 0, 2.05],
    ], float)
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
    q, _ = neutral_profile_and_derivative(z, alpha_e)
    a = 1.0 + TAPER_TAU * q
    central = np.abs(z) / AXIAL_SUPPORT <= 0.35
    return {
        "core_signs_pass": signs,
        "support_max_abs": support_max,
        "divergence_fd_max": float(np.max(np.abs(div))),
        "radial_scale_min": float(np.min(a)),
        "radial_scale_max": float(np.max(a)),
        "determinant_min": float(np.min(a * a)),
        "determinant_max": float(np.max(a * a)),
        "central_map_identity_error": float(np.max(np.abs(a[central] - 1.0))),
        "support_endpoint_identity_error": float(max(abs(a[0] - 1.0), abs(a[-1] - 1.0))),
    }


def response_diagnostics(f51, raw51, f52, raw52, redistribution_scale, alpha_e, seed=9175813):
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
    uneutral = []
    for t in TIMES:
        mask = ts == t
        uneutral.append(neutral_velocity(
            f52, raw52, pts[mask], t,
            redistribution_scale=redistribution_scale,
            alpha_e=alpha_e,
        ))
    uneutral = np.vstack(uneutral)

    cols = [
        (u52 - u51).ravel(),
        (ured - u52).ravel(),
        (uneutral - ured).ravel(),
    ]
    norms = [float(np.linalg.norm(c)) for c in cols]
    if min(norms) <= 1.0e-14:
        return {"rank": 0, "condition_number": float("inf"), "max_pairwise_abs_cosine": 1.0}
    A = np.column_stack([c / n for c, n in zip(cols, norms)])
    s = np.linalg.svd(A, compute_uv=False)
    cosines = []
    for i in range(3):
        for j in range(i + 1, 3):
            cosines.append(float(np.dot(A[:, i], A[:, j])))
    return {
        "rank": int(np.linalg.matrix_rank(A, tol=1.0e-10)),
        "singular_values": [float(v) for v in s],
        "condition_number": float(s[0] / s[-1]),
        "pairwise_cosines": cosines,
        "max_pairwise_abs_cosine": float(max(abs(v) for v in cosines)),
        "column_norms": norms,
        "neutral_taper_column_is_finite_child_response": True,
    }


def clean_rule(row):
    if not row["energy_solve"]["root_exists"]:
        return False
    c = CRITERIA
    m = row["morphology"]
    return bool(
        abs(row["energy_solve"]["relative_energy_error"]) <= c["energy_relative_error_max"]
        and max(x["tip_radial_rms_relative"] for x in m) <= -c["tip_radial_rms_reduction_min"]
        and max(abs(x["central_radial_rms_relative"]) for x in m) <= c["central_radial_rms_abs_change_max"]
        and row["central_velocity_identity"]["max_abs"] <= c["central_velocity_identity_abs_max"]
        and min(x["global_axial_rms_relative"] for x in m) >= -c["global_axial_rms_loss_max"]
        and min(x["tip_absolute_enstrophy_ratio"] for x in m) >= c["tip_absolute_enstrophy_ratio_min"]
        and max(abs(x["shoulder_radial_rms_relative"]) for x in m) <= c["shoulder_radial_rms_abs_change_max"]
        and max(x["shoulder_absolute_enstrophy_ratio"] for x in m) <= c["shoulder_absolute_enstrophy_ratio_max"]
        and max(x["collar_absolute_ratio"] for x in m) <= c["collar_absolute_enstrophy_ratio_max"]
        and max(abs(v - 1.0) for v in row["central_angular_retention"].values()) <= c["central_angular_retention_abs_error_max"]
        and row["structure"]["core_signs_pass"]
        and row["structure"]["support_max_abs"] <= c["support_max_abs"]
        and row["structure"]["divergence_fd_max"] <= c["divergence_fd_max"]
        and row["structure"]["determinant_min"] > c["determinant_min"]
        and row["response"]["rank"] >= c["response_rank_min"]
        and row["response"]["condition_number"] <= c["response_condition_max"]
        and row["response"]["max_pairwise_abs_cosine"] <= c["response_abs_cosine_max"]
    )


def run(out: Path):
    f52, raw52 = replay_st052.reconstruct()
    f51, raw51 = replay_st052.parent_field()
    redistribution_scale, parent_energy, red_unscaled_energy = base.child_scale(f52, raw52)
    energy_solve = solve_energy_alpha(f52, raw52, redistribution_scale)

    if energy_solve["root_exists"]:
        alpha_e = energy_solve["alpha_e"]
        morph = morphology_rows(f52, raw52, redistribution_scale, alpha_e)
        identity = central_identity(f52, raw52, redistribution_scale, alpha_e)
        angular = central_angular_retention(f52, raw52, redistribution_scale, alpha_e)
        structure = structure_preflight(f52, raw52, redistribution_scale, alpha_e)
        response = response_diagnostics(f51, raw51, f52, raw52, redistribution_scale, alpha_e)
    else:
        alpha_e = None
        morph = []
        identity = {"max_abs": None, "max_rel": None}
        angular = {}
        structure = {}
        response = {}

    row = {
        "redistribution_scale": redistribution_scale,
        "parent_energy": parent_energy,
        "redistribution_unscaled_energy": red_unscaled_energy,
        "energy_solve": energy_solve,
        "alpha_e": alpha_e,
        "post_piola_common_scale": 1.0,
        "morphology": morph,
        "central_velocity_identity": identity,
        "central_angular_retention": angular,
        "structure": structure,
        "response": response,
    }
    row["clean_energy_neutral_taper_capacity"] = clean_rule(row)

    report = {
        "task_id": TASK_ID,
        "prereg_issue": PREREG_ISSUE,
        "parent_id": PARENT_ID,
        "parent_head": PARENT_HEAD,
        "source_taper_head": SOURCE_TAPER_HEAD,
        "frozen_transform": {
            "tau": TAPER_TAU,
            "tip_window_abs_z_over_2": list(TIP_WINDOW),
            "shoulder_window_abs_z_over_2": list(SHOULDER_WINDOW),
            "energy_alpha_bracket": list(ENERGY_ALPHA_BRACKET),
            "energy_order": ENERGY_ORDER,
            "redistribution_alpha": base.SOURCE_ALPHA,
            "redistribution_gain": base.GAIN,
            "redistribution_inner_window": list(base.INNER_WINDOW),
            "redistribution_outer_window": list(base.OUTER_WINDOW),
        },
        "parameter_scan_performed": False,
        "energy_only_dependent_root_solve": True,
        "post_piola_common_scale_allowed": False,
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
    print("root_exists=", q["energy_solve"]["root_exists"])
    print("alpha_e=", q["alpha_e"])
    print("clean_energy_neutral_taper_capacity=", q["clean_energy_neutral_taper_capacity"])
    if q["energy_solve"]["root_exists"]:
        print("energy_relative_error=", q["energy_solve"]["relative_energy_error"])
        print("central_velocity_identity=", q["central_velocity_identity"])
        for m in q["morphology"]:
            print(
                "t=", m["time"],
                "tip_rr=", m["tip_radial_rms_relative"],
                "central_rr=", m["central_radial_rms_relative"],
                "shoulder_rr=", m["shoulder_radial_rms_relative"],
                "axial=", m["global_axial_rms_relative"],
            )
        print("response=", q["response"])


if __name__ == "__main__":
    main()
