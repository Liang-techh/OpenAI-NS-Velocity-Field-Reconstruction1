"""Screen the frozen ST050R-C swirl redistribution on frozen ST051-B.

Preregistered in issue #468 before evaluation. This Agent-7 increment freezes the
radial redistribution profile found on ST050R-C (including alpha_C) and applies
exactly one gain to ST051-B. It is a target-free expression-capacity diagnostic:
no pressure/forcing fit, held-out PDE promotion, image fitting, or canonical
candidate mutation occurs here.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

import replay_st051
from diagnostics import jets

TASK_ID = "CR003-ST051B-FROZEN-REDISTRIBUTION-TRANSFER-072"
PREREG_ISSUE = 468
PARENT_ID = "ST051-B"
PARENT_HEAD = "4b784f1b8457af2ead49295631d834d4e882000b"
SOURCE_TRANSFER_PR = 458
SOURCE_PARENT_ID = "ST050R-C"
SOURCE_ALPHA = 2.520520814687742
GAIN = 0.025
TIMES = (0.25, 0.50, 0.75)
RADII = (0.6, 0.9, 1.2)
INNER_WINDOW = (0.30, 1.05)
OUTER_WINDOW = (0.95, 1.85)
BOX_HALF_WIDTH = 1.95
GRID_SIZE = 33
RADIAL_SUPPORT = 2.0
AXIAL_SUPPORT = 2.0
COLLAR_FRACTION = 0.75

CRITERIA = dict(
    inner_gain_floor=0.02,
    mid_gain_floor=0.01,
    outer_gain_ceiling=0.005,
    normalization_deviation_max=0.005,
    inward_speed_loss_max=0.005,
    axial_rms_abs_change_max=0.02,
    radial_rms_growth_max=0.03,
    collar_ratio_max=1.25,
    support_max_abs=1e-12,
    divergence_fd_max=1e-5,
)


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


def compact_bump_derivative(r, lo, hi):
    r = np.asarray(r, dtype=float)
    mid = 0.5 * (lo + hi)
    half = 0.5 * (hi - lo)
    x = (r - mid) / half
    out = np.zeros_like(r)
    mask = np.abs(x) < 1.0
    xm = x[mask]
    g = np.exp(1.0 - 1.0 / (1.0 - xm * xm))
    out[mask] = g * (-2.0 * xm / (1.0 - xm * xm) ** 2) / half
    return out


def h_profile(r):
    return compact_bump(r, *INNER_WINDOW) - SOURCE_ALPHA * compact_bump(r, *OUTER_WINDOW)


def h_derivative(r):
    return compact_bump_derivative(r, *INNER_WINDOW) - SOURCE_ALPHA * compact_bump_derivative(r, *OUTER_WINDOW)


def meridian_velocity(f, raw, r, z, t):
    r, z = np.broadcast_arrays(np.asarray(r, float), np.asarray(z, float))
    pts = np.column_stack((r.ravel(), np.zeros(r.size), z.ravel()))
    u, _ = f.fields(raw, pts, np.full(len(pts), float(t)))
    return np.asarray(u, float).reshape(r.shape + (3,))


def reference_energy_and_moment(f, raw, order=72):
    x, wx = leggauss(order)
    y, wz = leggauss(order)
    r = x + 1.0
    z = 2.0 * y
    R, Z = np.meshgrid(r, z, indexing="ij")
    U = meridian_velocity(f, raw, R, Z, 0.25)
    w = wx[:, None] * (2.0 * wz[None, :]) * np.pi * R
    energy = float(np.sum(w * np.sum(U * U, axis=-1)))
    swirl_moment = float(np.sum(w * U[..., 1] ** 2 * h_profile(R)))
    inner = float(np.sum(w * U[..., 1] ** 2 * compact_bump(R, *INNER_WINDOW)))
    outer = float(np.sum(w * U[..., 1] ** 2 * compact_bump(R, *OUTER_WINDOW)))
    return dict(
        parent_energy=energy,
        frozen_profile_swirl_moment=swirl_moment,
        kinetic_energy_first_derivative_at_zero=2.0 * swirl_moment,
        normalized_first_derivative=(2.0 * swirl_moment / energy),
        inner_moment=inner,
        outer_moment=outer,
        alpha_that_would_rebalance_parent=(inner / outer),
    )


def raw_child_energy(f, raw, gain, order=72):
    x, wx = leggauss(order)
    y, wz = leggauss(order)
    r = x + 1.0
    z = 2.0 * y
    R, Z = np.meshgrid(r, z, indexing="ij")
    U = meridian_velocity(f, raw, R, Z, 0.25).copy()
    U[..., 1] *= 1.0 + float(gain) * h_profile(R)
    w = wx[:, None] * (2.0 * wz[None, :]) * np.pi * R
    return float(np.sum(w * np.sum(U * U, axis=-1)))


def child_scale(f, raw, parent_energy, gain=GAIN, order=72):
    energy = raw_child_energy(f, raw, gain, order=order)
    return float(np.sqrt(parent_energy / energy)), energy


def transformed_velocity(f, raw, points, t, gain, scale):
    points = np.asarray(points, float)
    u, _ = f.fields(raw, points, np.full(len(points), float(t)))
    u = np.asarray(u, float).copy()
    x, y = points[:, 0], points[:, 1]
    r = np.hypot(x, y)
    mask = r > 1e-14
    if np.any(mask):
        rx, ry = x[mask] / r[mask], y[mask] / r[mask]
        ux, uy = u[mask, 0].copy(), u[mask, 1].copy()
        ur = rx * ux + ry * uy
        ut = -ry * ux + rx * uy
        ut *= 1.0 + gain * h_profile(r[mask])
        u[mask, 0] = rx * ur - ry * ut
        u[mask, 1] = ry * ur + rx * ut
    return scale * u


def quadrature_vorticity_metrics(f, raw, gain, scale, t, nr=48, nz=72):
    x, wx = leggauss(nr)
    y, wz = leggauss(nz)
    S, Z = np.meshgrid(2.0 * (x + 1.0), 2.0 * y, indexing="ij")
    s, z = S.ravel(), Z.ravel()
    r = np.sqrt(s)
    w = (4.0 * np.pi * np.outer(wx, wz)).ravel()
    j = jets(f, raw, r, z, float(t))
    ur, ut, uz = j[:, 0], j[:, 1], j[:, 2]
    wr, wt, wz0 = j[:, 3], j[:, 4], j[:, 5]
    mult = 1.0 + gain * h_profile(r)
    dmult = gain * h_derivative(r)
    wr1 = scale * mult * wr
    wt1 = scale * wt
    wz1 = scale * (mult * wz0 + dmult * ut)
    om2 = wr1 * wr1 + wt1 * wt1 + wz1 * wz1
    en = float(w @ om2)
    axial = float(np.sqrt((w * z * z) @ om2 / en))
    radial = float(np.sqrt((w * s) @ om2 / en))
    collar = (r > COLLAR_FRACTION * RADIAL_SUPPORT) | (np.abs(z) > COLLAR_FRACTION * AXIAL_SUPPORT)
    collar_abs = float(np.sum(w * om2 * collar))
    return dict(enstrophy=en, axial_rms=axial, radial_rms=radial, collar_absolute=collar_abs)


def weighted_quantile(values, weights, q):
    values = np.asarray(values, float).ravel()
    weights = np.asarray(weights, float).ravel()
    order = np.argsort(values)
    values, weights = values[order], weights[order]
    total = float(np.sum(weights))
    if total <= np.finfo(float).tiny:
        raise RuntimeError("inactive vorticity weights")
    idx = int(np.searchsorted(np.cumsum(weights), q * total, side="left"))
    return float(values[min(idx, len(values) - 1)])


def grid_morphology(f, raw, *, time, gain, scale, grid_size=GRID_SIZE):
    axis = np.linspace(-BOX_HALF_WIDTH, BOX_HALF_WIDTH, int(grid_size))
    h = float(axis[1] - axis[0])
    x, y, z = np.meshgrid(axis, axis, axis, indexing="ij")
    points = np.column_stack((x.ravel(), y.ravel(), z.ravel()))
    if gain == 0.0 and scale == 1.0:
        velocity, _ = f.fields(raw, points, np.full(len(points), float(time)))
        velocity = np.asarray(velocity, float)
    else:
        velocity = transformed_velocity(f, raw, points, time, gain, scale)
    V = velocity.reshape(grid_size, grid_size, grid_size, 3)
    u, v, w = V[..., 0], V[..., 1], V[..., 2]
    two_h = 2.0 * h
    ox = (w[1:-1, 2:, 1:-1] - w[1:-1, :-2, 1:-1]) / two_h - (v[1:-1, 1:-1, 2:] - v[1:-1, 1:-1, :-2]) / two_h
    oy = (u[1:-1, 1:-1, 2:] - u[1:-1, 1:-1, :-2]) / two_h - (w[2:, 1:-1, 1:-1] - w[:-2, 1:-1, 1:-1]) / two_h
    oz = (v[2:, 1:-1, 1:-1] - v[:-2, 1:-1, 1:-1]) / two_h - (u[1:-1, 2:, 1:-1] - u[1:-1, :-2, 1:-1]) / two_h
    om2 = ox * ox + oy * oy + oz * oz
    total = float(np.sum(om2))
    xi, yi, zi = x[1:-1, 1:-1, 1:-1], y[1:-1, 1:-1, 1:-1], z[1:-1, 1:-1, 1:-1]
    r2 = xi * xi + yi * yi
    absz = np.abs(zi)
    collar = (np.sqrt(r2) >= COLLAR_FRACTION * RADIAL_SUPPORT) | (absz >= COLLAR_FRACTION * AXIAL_SUPPORT)
    return dict(
        axial_rms=float(np.sqrt(np.sum(zi * zi * om2) / total)),
        radial_rms=float(np.sqrt(np.sum(r2 * om2) / total)),
        axial_q90_over_support=weighted_quantile(absz, om2, 0.90) / AXIAL_SUPPORT,
        axial_q99_over_support=weighted_quantile(absz, om2, 0.99) / AXIAL_SUPPORT,
        outer_065_enstrophy_fraction=float(np.sum(om2[absz >= 0.65 * AXIAL_SUPPORT]) / total),
        outer_075_enstrophy_fraction=float(np.sum(om2[absz >= 0.75 * AXIAL_SUPPORT]) / total),
        collar_fraction=float(np.sum(om2[collar]) / total),
        vorticity_rms=float(np.sqrt(np.mean(om2))),
        vorticity_max=float(np.sqrt(np.max(om2))),
    )


def sign_support_preflight(f, raw, gain, scale):
    core = np.array([[r, 0.0, z] for r in (0.1, 0.3, 0.6) for z in (-0.5, -0.2, 0.2, 0.5)], float)
    u = transformed_velocity(f, raw, core, 0.5, gain, scale)
    signs = bool(np.all(u[:, 0] < 0.0) and np.all(u[:, 1] > 0.0) and np.all(np.sign(core[:, 2]) * u[:, 2] > 0.0))
    outside = np.array([[2.05,0,0],[-2.05,0,0],[0,0,2.05],[0,0,-2.05],[1.8,0,2.05]], float)
    support_max = float(np.max(np.abs(transformed_velocity(f, raw, outside, 0.5, gain, scale))))
    return signs, support_max


def divergence_fd(f, raw, gain, scale, seed=9175721, h=1e-5):
    rng = np.random.default_rng(seed)
    pts = rng.uniform(-1.25, 1.25, size=(24, 3))
    div = np.zeros(len(pts))
    for axis in range(3):
        d = np.zeros(3)
        d[axis] = h
        up = transformed_velocity(f, raw, pts + d, 0.5, gain, scale)
        um = transformed_velocity(f, raw, pts - d, 0.5, gain, scale)
        div += (up[:, axis] - um[:, axis]) / (2.0 * h)
    return float(np.max(np.abs(div)))


def clean_rule(row):
    c = CRITERIA
    return bool(
        row["angular_gain_r06"] >= c["inner_gain_floor"]
        and row["angular_gain_r09"] >= c["mid_gain_floor"]
        and row["angular_gain_r12"] <= c["outer_gain_ceiling"]
        and abs(row["normalization"] - 1.0) <= c["normalization_deviation_max"]
        and row["inward_speed_loss"] <= c["inward_speed_loss_max"]
        and max(abs(x) for x in row["axial_rms_relative_changes"]) <= c["axial_rms_abs_change_max"]
        and max(row["radial_rms_relative_changes"]) <= c["radial_rms_growth_max"]
        and max(row["collar_ratios"]) <= c["collar_ratio_max"]
        and row["core_signs_pass"]
        and row["support_max_abs"] <= c["support_max_abs"]
        and row["divergence_fd_max"] <= c["divergence_fd_max"]
    )


def run(out: Path):
    f, raw = replay_st051.reconstruct(PARENT_ID)
    balance = reference_energy_and_moment(f, raw)
    scale, raw_energy = child_scale(f, raw, balance["parent_energy"])

    angular = {r: float(scale * (1.0 + GAIN * h_profile(np.array([r]))[0]) - 1.0) for r in RADII}
    base_quad = {t: quadrature_vorticity_metrics(f, raw, 0.0, 1.0, t) for t in TIMES}
    child_quad = {t: quadrature_vorticity_metrics(f, raw, GAIN, scale, t) for t in TIMES}
    axial = [child_quad[t]["axial_rms"] / base_quad[t]["axial_rms"] - 1.0 for t in TIMES]
    radial = [child_quad[t]["radial_rms"] / base_quad[t]["radial_rms"] - 1.0 for t in TIMES]
    collar = [child_quad[t]["collar_absolute"] / base_quad[t]["collar_absolute"] for t in TIMES]

    grid_rows = []
    for t in TIMES:
        base = grid_morphology(f, raw, time=t, gain=0.0, scale=1.0)
        child = grid_morphology(f, raw, time=t, gain=GAIN, scale=scale)
        grid_rows.append(dict(
            time=t,
            parent=base,
            child=child,
            child_minus_parent={k: float(child[k] - base[k]) for k in base},
            child_relative_to_parent={k: float(child[k] / max(abs(base[k]), 1e-300) - 1.0) for k in base},
        ))

    signs, support = sign_support_preflight(f, raw, GAIN, scale)
    div = divergence_fd(f, raw, GAIN, scale)
    row = dict(
        gain=GAIN,
        normalization=scale,
        raw_child_energy=raw_energy,
        inward_speed_loss=float(max(0.0, 1.0 - scale)),
        angular_gain_r06=angular[0.6],
        angular_gain_r09=angular[0.9],
        angular_gain_r12=angular[1.2],
        axial_rms_relative_changes=[float(x) for x in axial],
        radial_rms_relative_changes=[float(x) for x in radial],
        collar_ratios=[float(x) for x in collar],
        core_signs_pass=signs,
        support_max_abs=support,
        divergence_fd_max=div,
    )
    row["clean_frozen_profile_transfer"] = clean_rule(row)

    result = dict(
        task_id=TASK_ID,
        prereg_issue=PREREG_ISSUE,
        parent_id=PARENT_ID,
        parent_head=PARENT_HEAD,
        source_transfer_pr=SOURCE_TRANSFER_PR,
        source_parent_id=SOURCE_PARENT_ID,
        frozen_transform=dict(alpha=SOURCE_ALPHA, inner_window=list(INNER_WINDOW), outer_window=list(OUTER_WINDOW), gain=GAIN),
        diagnostic_added_velocity_degrees=1,
        parameter_scan_performed=False,
        balance_on_st051b=balance,
        criteria=CRITERIA,
        row=row,
        grid_contract=dict(grid_size=GRID_SIZE, box_half_width=BOX_HALF_WIDTH, times=list(TIMES), classification="target-free visualization fingerprint"),
        grid_rows=grid_rows,
        recommendation=(
            "Frozen ST050R-C redistribution remains clean on ST051-B; do not add another generic swirl basis for this local role. If a governed ST051-B visualization child is materialized, retain this frozen radial profile identity, then rebuild compatible pressure/restricted forcing and run fresh held-out full momentum plus material paths."
            if row["clean_frozen_profile_transfer"] else
            "Frozen ST050R-C redistribution does not transfer cleanly to ST051-B under the preregistered rule. Stop here; do not re-balance alpha or widen gain in this increment. Use the failed channel metrics to decide the next minimal representation change."
        ),
        canonical_velocity_changed=False,
        candidate_artifact_changed=False,
        pressure_or_force_changed=False,
        held_out_pde_residual_evaluated=False,
        material_paths_integrated=False,
        public_image_used=False,
        visualization_ready=False,
        visual_correspondence_verified=False,
        pde_validated=False,
        source_correspondence_verified=False,
        paper_exact=False,
        openai_field_identified=False,
        scope="One fixed cross-backbone expression-capacity transfer; no visual/PDE acceptance claim.",
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    run(args.out)
