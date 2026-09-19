"""Transfer one frozen swirl-redistribution degree onto ST052-M.

Preregistered in issue #527 before evaluation.  This is an Agent-7 expression-
capacity diagnostic stacked on the exact ST052-M head from PR #508.  The radial
profile and gain are frozen from earlier ST051-B work; there is no parameter
scan, pressure/forcing refit, image fitting, or candidate promotion here.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "experiments" / "root_st051"))
import replay_st052

TASK_ID = "CR003-ST052M-FROZEN-REDISTRIBUTION-078"
PREREG_ISSUE = 527
PARENT_ID = "ST052-M"
PARENT_HEAD = "b3b8bfdbe1077f9ec967d158602951997d81e17d"
UPSTREAM_PARENT_ID = "ST051-B"
SOURCE_ALPHA = 2.520520814687742
GAIN = 0.05
INNER_WINDOW = (0.30, 1.05)
OUTER_WINDOW = (0.95, 1.85)
TIMES = (0.25, 0.50, 0.75)
RADII = (0.6, 0.9, 1.2)
GRID_SIZE = 23
BOX_HALF_WIDTH = 1.85
COLLAR_FRACTION = 0.75

CRITERIA = {
    "inner_gain_floor": 0.04,
    "mid_gain_floor": 0.025,
    "outer_gain_ceiling": 0.0,
    "normalization_deviation_max": 0.005,
    "axial_rms_abs_change_max": 0.02,
    "radial_rms_growth_max": 0.03,
    "collar_ratio_max": 1.25,
    "support_max_abs": 1e-12,
    "divergence_fd_max": 1e-5,
    "response_rank_min": 2,
    "response_condition_max": 5.0,
    "response_abs_cosine_max": 0.90,
}


def compact_bump(r, lo, hi):
    r = np.asarray(r, dtype=float)
    mid = 0.5 * (lo + hi)
    half = 0.5 * (hi - lo)
    x = (r - mid) / half
    out = np.zeros_like(r)
    mask = np.abs(x) < 1.0
    xm = x[mask]
    out[mask] = np.exp(1.0 - 1.0 / (1.0 - xm * xm))
    return out


def h_profile(r):
    return compact_bump(r, *INNER_WINDOW) - SOURCE_ALPHA * compact_bump(r, *OUTER_WINDOW)


def base_velocity(f, raw, points, time):
    points = np.asarray(points, float)
    u, _ = f.fields(raw, points, np.full(len(points), float(time)))
    return np.asarray(u, float)


def base_velocity_many(f, raw, points, times):
    u, _ = f.fields(raw, np.asarray(points, float), np.asarray(times, float))
    return np.asarray(u, float)


def transformed_velocity(f, raw, points, time, gain=GAIN, scale=1.0):
    points = np.asarray(points, float)
    u = base_velocity(f, raw, points, time).copy()
    x, y = points[:, 0], points[:, 1]
    r = np.hypot(x, y)
    mask = r > 1e-14
    if np.any(mask):
        rx = x[mask] / r[mask]
        ry = y[mask] / r[mask]
        ux = u[mask, 0].copy()
        uy = u[mask, 1].copy()
        ur = rx * ux + ry * uy
        ut = -ry * ux + rx * uy
        ut *= 1.0 + float(gain) * h_profile(r[mask])
        u[mask, 0] = rx * ur - ry * ut
        u[mask, 1] = ry * ur + rx * ut
    return float(scale) * u


def transformed_velocity_many(f, raw, points, times, gain=GAIN, scale=1.0):
    points = np.asarray(points, float)
    times = np.asarray(times, float)
    u = base_velocity_many(f, raw, points, times).copy()
    x, y = points[:, 0], points[:, 1]
    r = np.hypot(x, y)
    mask = r > 1e-14
    if np.any(mask):
        rx = x[mask] / r[mask]
        ry = y[mask] / r[mask]
        ux = u[mask, 0].copy()
        uy = u[mask, 1].copy()
        ur = rx * ux + ry * uy
        ut = -ry * ux + rx * uy
        ut *= 1.0 + float(gain) * h_profile(r[mask])
        u[mask, 0] = rx * ur - ry * ut
        u[mask, 1] = ry * ur + rx * ut
    return float(scale) * u


def reference_energy(f, raw, gain=0.0, order=64):
    xr, wr = leggauss(order)
    xz, wz = leggauss(order)
    r = xr + 1.0
    z = 2.0 * xz
    R, Z = np.meshgrid(r, z, indexing="ij")
    pts = np.column_stack((R.ravel(), np.zeros(R.size), Z.ravel()))
    u = base_velocity(f, raw, pts, 0.25).reshape(R.shape + (3,)).copy()
    if gain:
        u[..., 1] *= 1.0 + float(gain) * h_profile(R)
    # pi*r is the axisymmetric 1/2 |u|^2 energy weight; the common factor
    # cancels in the normalization ratio.
    weight = wr[:, None] * (2.0 * wz[None, :]) * np.pi * R
    return float(np.sum(weight * np.sum(u * u, axis=-1)))


def child_scale(f, raw):
    e0 = reference_energy(f, raw, gain=0.0)
    e1 = reference_energy(f, raw, gain=GAIN)
    return float(np.sqrt(e0 / e1)), e0, e1


def angular_rate_changes(f, raw, scale):
    zs = np.array([-0.60, -0.30, 0.0, 0.30, 0.60])
    out = {}
    for r in RADII:
        pts = np.column_stack((np.full(len(zs), r), np.zeros(len(zs)), zs))
        parent = base_velocity(f, raw, pts, 0.50)
        child = transformed_velocity(f, raw, pts, 0.50, GAIN, scale)
        p = float(np.mean(np.abs(parent[:, 1]) / r))
        c = float(np.mean(np.abs(child[:, 1]) / r))
        out[str(r)] = {"parent": p, "child": c, "relative_change": c / p - 1.0}
    return out


def curl_fd(velocity_fn, points, time, h=1e-3):
    points = np.asarray(points, float)
    D = np.empty((len(points), 3, 3), float)  # velocity component, derivative axis
    for axis in range(3):
        d = np.zeros(3)
        d[axis] = h
        up = velocity_fn(points + d, time)
        um = velocity_fn(points - d, time)
        D[:, :, axis] = (up - um) / (2.0 * h)
    return np.column_stack((
        D[:, 2, 1] - D[:, 1, 2],
        D[:, 0, 2] - D[:, 2, 0],
        D[:, 1, 0] - D[:, 0, 1],
    ))


def morphology(f, raw, time, scale):
    axis = np.linspace(-BOX_HALF_WIDTH, BOX_HALF_WIDTH, GRID_SIZE)
    x, y, z = np.meshgrid(axis, axis, axis, indexing="ij")
    pts = np.column_stack((x.ravel(), y.ravel(), z.ravel()))
    r = np.hypot(pts[:, 0], pts[:, 1])
    absz = np.abs(pts[:, 2])
    base_fn = lambda q, tt: base_velocity(f, raw, q, tt)
    child_fn = lambda q, tt: transformed_velocity(f, raw, q, tt, GAIN, scale)
    wb = curl_fd(base_fn, pts, time)
    wc = curl_fd(child_fn, pts, time)

    def metrics(w):
        om2 = np.sum(w * w, axis=1)
        total = float(np.sum(om2))
        collar = (r >= 2.0 * COLLAR_FRACTION) | (absz >= 2.0 * COLLAR_FRACTION)
        return {
            "enstrophy_sum": total,
            "axial_rms": float(np.sqrt(np.sum(pts[:, 2] ** 2 * om2) / total)),
            "radial_rms": float(np.sqrt(np.sum(r * r * om2) / total)),
            "collar_absolute": float(np.sum(om2[collar])),
        }

    b = metrics(wb)
    c = metrics(wc)
    return {
        "time": float(time),
        "parent": b,
        "child": c,
        "axial_rms_relative_change": c["axial_rms"] / b["axial_rms"] - 1.0,
        "radial_rms_relative_change": c["radial_rms"] / b["radial_rms"] - 1.0,
        "collar_absolute_ratio": c["collar_absolute"] / b["collar_absolute"],
    }


def sign_support_preflight(f, raw, scale):
    core = np.array([[r, 0.0, z] for r in (0.1, 0.3, 0.6) for z in (-0.5, -0.2, 0.2, 0.5)], float)
    u = transformed_velocity(f, raw, core, 0.50, GAIN, scale)
    signs = bool(
        np.all(u[:, 0] < 0.0)
        and np.all(u[:, 1] > 0.0)
        and np.all(np.sign(core[:, 2]) * u[:, 2] > 0.0)
    )
    outside = np.array([
        [2.05, 0.0, 0.0], [-2.05, 0.0, 0.0], [0.0, 2.05, 0.0],
        [0.0, 0.0, 2.05], [0.0, 0.0, -2.05], [1.8, 0.0, 2.05],
    ])
    support_max = float(np.max(np.abs(transformed_velocity(f, raw, outside, 0.50, GAIN, scale))))
    return signs, support_max


def divergence_fd(f, raw, scale, seed=9175781, h=1e-5):
    rng = np.random.default_rng(seed)
    pts = rng.uniform(-1.20, 1.20, size=(24, 3))
    div = np.zeros(len(pts))
    for axis in range(3):
        d = np.zeros(3)
        d[axis] = h
        up = transformed_velocity(f, raw, pts + d, 0.50, GAIN, scale)
        um = transformed_velocity(f, raw, pts - d, 0.50, GAIN, scale)
        div += (up[:, axis] - um[:, axis]) / (2.0 * h)
    return float(np.max(np.abs(div)))


def response_independence(f51, raw51, f52, raw52, scale, seed=9175782):
    rng = np.random.default_rng(seed)
    points = []
    times = []
    for t in TIMES:
        count = 0
        while count < 12:
            r = rng.uniform(0.18, 1.45)
            th = rng.uniform(0.0, 2.0 * np.pi)
            z = rng.uniform(-1.15, 1.15)
            points.append([r * np.cos(th), r * np.sin(th), z])
            times.append(t)
            count += 1
    points = np.asarray(points, float)
    times = np.asarray(times, float)
    u51 = base_velocity_many(f51, raw51, points, times)
    u52 = base_velocity_many(f52, raw52, points, times)
    uc = transformed_velocity_many(f52, raw52, points, times, GAIN, scale)
    a = (u52 - u51).ravel()
    b = (uc - u52).ravel()
    na, nb = float(np.linalg.norm(a)), float(np.linalg.norm(b))
    if na <= 1e-14 or nb <= 1e-14:
        return {"rank": 0, "condition_number": float("inf"), "cosine": float("nan"), "norm_st052_displacement": na, "norm_redistribution_response": nb}
    A = np.column_stack((a / na, b / nb))
    s = np.linalg.svd(A, compute_uv=False)
    return {
        "rank": int(np.linalg.matrix_rank(A, tol=1e-10)),
        "singular_values": [float(v) for v in s],
        "condition_number": float(s[0] / s[-1]),
        "cosine": float(np.dot(a, b) / (na * nb)),
        "norm_st052_displacement": na,
        "norm_redistribution_response": nb,
    }


def clean_rule(row):
    c = CRITERIA
    return bool(
        row["angular_gain_r06"] >= c["inner_gain_floor"]
        and row["angular_gain_r09"] >= c["mid_gain_floor"]
        and row["angular_gain_r12"] <= c["outer_gain_ceiling"]
        and abs(row["normalization"] - 1.0) <= c["normalization_deviation_max"]
        and max(abs(v) for v in row["axial_rms_relative_changes"]) <= c["axial_rms_abs_change_max"]
        and max(row["radial_rms_relative_changes"]) <= c["radial_rms_growth_max"]
        and max(row["collar_ratios"]) <= c["collar_ratio_max"]
        and row["core_signs_pass"]
        and row["support_max_abs"] <= c["support_max_abs"]
        and row["divergence_fd_max"] <= c["divergence_fd_max"]
        and row["response_rank"] >= c["response_rank_min"]
        and row["response_condition_number"] <= c["response_condition_max"]
        and abs(row["response_cosine"]) <= c["response_abs_cosine_max"]
    )


def run(out: Path):
    f52, raw52 = replay_st052.reconstruct()
    f51, raw51 = replay_st052.parent_field()
    scale, parent_energy, unscaled_child_energy = child_scale(f52, raw52)
    angular = angular_rate_changes(f52, raw52, scale)
    morph = [morphology(f52, raw52, t, scale) for t in TIMES]
    signs, support_max = sign_support_preflight(f52, raw52, scale)
    divmax = divergence_fd(f52, raw52, scale)
    indep = response_independence(f51, raw51, f52, raw52, scale)

    row = {
        "gain": GAIN,
        "normalization": scale,
        "angular_gain_r06": angular["0.6"]["relative_change"],
        "angular_gain_r09": angular["0.9"]["relative_change"],
        "angular_gain_r12": angular["1.2"]["relative_change"],
        "axial_rms_relative_changes": [m["axial_rms_relative_change"] for m in morph],
        "radial_rms_relative_changes": [m["radial_rms_relative_change"] for m in morph],
        "collar_ratios": [m["collar_absolute_ratio"] for m in morph],
        "core_signs_pass": signs,
        "support_max_abs": support_max,
        "divergence_fd_max": divmax,
        "response_rank": indep["rank"],
        "response_condition_number": indep["condition_number"],
        "response_cosine": indep["cosine"],
    }
    row["clean_frozen_profile_transfer"] = clean_rule(row)

    report = {
        "task_id": TASK_ID,
        "prereg_issue": PREREG_ISSUE,
        "parent_id": PARENT_ID,
        "parent_head": PARENT_HEAD,
        "upstream_parent_id": UPSTREAM_PARENT_ID,
        "frozen_transform": {
            "alpha": SOURCE_ALPHA,
            "gain": GAIN,
            "inner_window": INNER_WINDOW,
            "outer_window": OUTER_WINDOW,
        },
        "criteria": CRITERIA,
        "parameter_scan_performed": False,
        "diagnostic_added_velocity_degrees": 1,
        "reference_energy": parent_energy,
        "unscaled_child_reference_energy": unscaled_child_energy,
        "angular_rate": angular,
        "morphology": morph,
        "response_independence": indep,
        "row": row,
        "canonical_velocity_changed": False,
        "candidate_artifact_changed": False,
        "pressure_or_force_changed": False,
        "held_out_pde_residual_evaluated": False,
        "material_paths_integrated": False,
        "public_image_used": False,
        "visualization_ready": False,
        "visual_correspondence_verified": False,
        "pde_validated": False,
        "source_correspondence_verified": False,
        "interpretation": (
            "A clean result means the already-frozen one-dimensional winding-control channel "
            "transfers to ST052-M with small morphology cost and remains locally distinct from "
            "the ST052 residual-fit displacement. It does not transfer ST052 PDE evidence to "
            "the externally transformed velocity child."
        ),
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(row, indent=2, sort_keys=True))
    print(json.dumps(indep, indent=2, sort_keys=True))
    return report


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    run(args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
