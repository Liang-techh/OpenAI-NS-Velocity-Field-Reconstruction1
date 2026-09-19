"""Screen one axially localized radial Piola taper on ST052-M + frozen swirl redistribution.

Preregistered in issue #542 before evaluation.  This Agent-7 increment adds
exactly one diagnostic geometry coordinate (tau=.05), with a fixed C-infinity
axial window.  There is no parameter scan, image fit, pressure/forcing refit,
material-path integration, or PDE receipt transfer.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

import agent7_st052m_frozen_redistribution as base
import replay_st052

TASK_ID = "CR003-ST052M-LOCALIZED-TAPER-080"
PREREG_ISSUE = 542
PARENT_ID = "ST052-M+frozen-redistribution-.05"
PARENT_HEAD = "779ffca71066e2864496d37de55a7aafc45d6f57"
TAPER_TAU = 0.05
AXIAL_SUPPORT = 2.0
WINDOW = (0.50, 0.82)  # in |z|/2
TIMES = base.TIMES
MORPH_NR = 40
MORPH_NZ = 56
VORTICITY_FD_STEP = 1.0e-3
DIVERGENCE_FD_STEP = 1.0e-5

CRITERIA = {
    "tip_radial_rms_reduction_min": 0.01,
    "central_radial_rms_abs_change_max": 0.0025,
    "global_axial_rms_loss_max": 0.01,
    "tip_absolute_enstrophy_ratio_min": 0.90,
    "collar_absolute_enstrophy_ratio_max": 1.15,
    "taper_scale_deviation_max": 0.02,
    "central_angular_retention_min": 0.99,
    "support_max_abs": 1.0e-12,
    "divergence_fd_max": 1.0e-5,
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


def localized_profile_and_derivative(z):
    """Return q(z), dq/dz for a symmetric C-infinity bump in .50<|z|/2<.82."""
    z = np.asarray(z, dtype=float)
    if not np.all(np.isfinite(z)):
        raise ValueError("z must be finite")
    s = np.abs(z) / AXIAL_SUPPORT
    lo, hi = WINDOW
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


def taper_map(points, tau=TAPER_TAU):
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError("points must have shape (n,3)")
    tau = float(tau)
    if not np.isfinite(tau) or abs(tau) > 0.15:
        raise ValueError("tau outside audited local range")
    q, dq = localized_profile_and_derivative(points[:, 2])
    a = 1.0 + tau * q
    ap = tau * dq
    if np.any(a <= 0.0):
        raise RuntimeError("localized taper lost orientation")
    mapped = points.copy()
    mapped[:, 0] *= a
    mapped[:, 1] *= a
    return mapped, a, ap


def control_velocity(f, raw, points, time, redistribution_scale):
    return base.transformed_velocity(
        f, raw, points, time, gain=base.GAIN, scale=redistribution_scale
    )


def taper_velocity(
    f,
    raw,
    points,
    time,
    *,
    redistribution_scale,
    tau=TAPER_TAU,
    taper_scale=1.0,
):
    """Contravariant Piola pullback for F=(a(z)x,a(z)y,z)."""
    points = np.asarray(points, dtype=float)
    mapped, a, ap = taper_map(points, tau)
    u = control_velocity(f, raw, mapped, time, redistribution_scale)
    x = points[:, 0]
    y = points[:, 1]
    uz = u[:, 2].copy()
    out = np.empty_like(u)
    out[:, 0] = a * u[:, 0] - a * ap * x * uz
    out[:, 1] = a * u[:, 1] - a * ap * y * uz
    out[:, 2] = a * a * uz
    return float(taper_scale) * out


def axisymmetric_energy(velocity_fn, order=64):
    xr, wr = leggauss(int(order))
    xz, wz = leggauss(int(order))
    r = xr + 1.0
    z = 2.0 * xz
    R, Z = np.meshgrid(r, z, indexing="ij")
    pts = np.column_stack((R.ravel(), np.zeros(R.size), Z.ravel()))
    u = np.asarray(velocity_fn(pts, 0.25), float).reshape(R.shape + (3,))
    weight = wr[:, None] * (2.0 * wz[None, :]) * np.pi * R
    return float(np.sum(weight * np.sum(u * u, axis=-1)))


def taper_energy_scale(f, raw, redistribution_scale):
    control = lambda p, t: control_velocity(f, raw, p, t, redistribution_scale)
    raw_child = lambda p, t: taper_velocity(
        f, raw, p, t, redistribution_scale=redistribution_scale, taper_scale=1.0
    )
    e0 = axisymmetric_energy(control, order=64)
    e1 = axisymmetric_energy(raw_child, order=64)
    scale = float(np.sqrt(e0 / e1))
    return scale, e0, e1


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
    collar = (R >= 1.50) | (np.abs(Z) >= 1.50)

    def band(mask):
        e = float(np.sum(weight * om2 * mask))
        rr = float(np.sqrt(np.sum(weight * R * R * om2 * mask) / e))
        return e, rr

    tip_e, tip_rr = band(tip)
    central_e, central_rr = band(central)
    return {
        "enstrophy": total,
        "axial_rms": float(np.sqrt(np.sum(weight * Z * Z * om2) / total)),
        "radial_rms": float(np.sqrt(np.sum(weight * R * R * om2) / total)),
        "tip_absolute_enstrophy": tip_e,
        "tip_radial_rms": tip_rr,
        "central_absolute_enstrophy": central_e,
        "central_radial_rms": central_rr,
        "collar_absolute_enstrophy": float(np.sum(weight * om2 * collar)),
    }


def morphology_rows(f, raw, redistribution_scale, taper_scale):
    control = lambda p, t: control_velocity(f, raw, p, t, redistribution_scale)
    child = lambda p, t: taper_velocity(
        f, raw, p, t,
        redistribution_scale=redistribution_scale,
        taper_scale=taper_scale,
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
            "global_axial_rms_relative": float(c["axial_rms"] / b["axial_rms"] - 1.0),
            "global_radial_rms_relative": float(c["radial_rms"] / b["radial_rms"] - 1.0),
            "tip_absolute_enstrophy_ratio": float(c["tip_absolute_enstrophy"] / b["tip_absolute_enstrophy"]),
            "central_absolute_enstrophy_ratio": float(c["central_absolute_enstrophy"] / b["central_absolute_enstrophy"]),
            "collar_absolute_ratio": float(c["collar_absolute_enstrophy"] / b["collar_absolute_enstrophy"]),
        })
    return rows


def central_angular_retention(f, raw, redistribution_scale, taper_scale):
    zs = np.array([-0.30, 0.0, 0.30])
    out = {}
    for r in (0.6, 0.9):
        pts = np.column_stack((np.full(len(zs), r), np.zeros(len(zs)), zs))
        u0 = control_velocity(f, raw, pts, 0.50, redistribution_scale)
        u1 = taper_velocity(
            f, raw, pts, 0.50,
            redistribution_scale=redistribution_scale,
            taper_scale=taper_scale,
        )
        a0 = float(np.mean(np.abs(u0[:, 1]) / r))
        a1 = float(np.mean(np.abs(u1[:, 1]) / r))
        out[str(r)] = float(a1 / a0)
    return out


def structure_preflight(f, raw, redistribution_scale, taper_scale, seed=9175801):
    child = lambda p, t: taper_velocity(
        f, raw, p, t,
        redistribution_scale=redistribution_scale,
        taper_scale=taper_scale,
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
    q, _ = localized_profile_and_derivative(z)
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


def response_diagnostics(f51, raw51, f52, raw52, redistribution_scale, seed=9175802):
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
    for t in TIMES:
        mask = ts == t
        plus.append(taper_velocity(
            f52, raw52, pts[mask], t,
            redistribution_scale=redistribution_scale,
            tau=eps,
            taper_scale=1.0,
        ).ravel())
        minus.append(taper_velocity(
            f52, raw52, pts[mask], t,
            redistribution_scale=redistribution_scale,
            tau=-eps,
            taper_scale=1.0,
        ).ravel())
    taper_col = (np.concatenate(plus) - np.concatenate(minus)) / (2.0 * eps)

    cols = [
        (u52 - u51).ravel(),
        (ured - u52).ravel(),
        taper_col,
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
    }


def clean_rule(row):
    c = CRITERIA
    return bool(
        max(m["tip_radial_rms_relative"] for m in row["morphology"]) <= -c["tip_radial_rms_reduction_min"]
        and max(abs(m["central_radial_rms_relative"]) for m in row["morphology"]) <= c["central_radial_rms_abs_change_max"]
        and min(m["global_axial_rms_relative"] for m in row["morphology"]) >= -c["global_axial_rms_loss_max"]
        and min(m["tip_absolute_enstrophy_ratio"] for m in row["morphology"]) >= c["tip_absolute_enstrophy_ratio_min"]
        and max(m["collar_absolute_ratio"] for m in row["morphology"]) <= c["collar_absolute_enstrophy_ratio_max"]
        and abs(row["taper_scale"] - 1.0) <= c["taper_scale_deviation_max"]
        and min(row["central_angular_retention"].values()) >= c["central_angular_retention_min"]
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
    taper_scale, control_energy, taper_raw_energy = taper_energy_scale(
        f52, raw52, redistribution_scale
    )
    morph = morphology_rows(f52, raw52, redistribution_scale, taper_scale)
    angular = central_angular_retention(f52, raw52, redistribution_scale, taper_scale)
    structure = structure_preflight(f52, raw52, redistribution_scale, taper_scale)
    response = response_diagnostics(f51, raw51, f52, raw52, redistribution_scale)

    row = {
        "redistribution_scale": redistribution_scale,
        "taper_scale": taper_scale,
        "parent_energy": parent_energy,
        "redistribution_unscaled_energy": red_unscaled_energy,
        "control_energy": control_energy,
        "taper_raw_energy": taper_raw_energy,
        "morphology": morph,
        "central_angular_retention": angular,
        "structure": structure,
        "response": response,
    }
    row["clean_localized_taper_capacity"] = clean_rule(row)

    report = {
        "task_id": TASK_ID,
        "prereg_issue": PREREG_ISSUE,
        "parent_id": PARENT_ID,
        "parent_head": PARENT_HEAD,
        "frozen_transform": {
            "tau": TAPER_TAU,
            "axial_window_abs_z_over_2": list(WINDOW),
            "redistribution_alpha": base.SOURCE_ALPHA,
            "redistribution_gain": base.GAIN,
            "redistribution_inner_window": list(base.INNER_WINDOW),
            "redistribution_outer_window": list(base.OUTER_WINDOW),
        },
        "parameter_scan_performed": False,
        "diagnostic_added_velocity_degrees": 1,
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
    print("clean_localized_taper_capacity=", q["clean_localized_taper_capacity"])
    print("taper_scale=", q["taper_scale"])
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
