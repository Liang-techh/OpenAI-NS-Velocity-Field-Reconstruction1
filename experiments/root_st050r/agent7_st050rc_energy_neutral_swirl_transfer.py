"""Transfer-screen one energy-neutral inner/mid swirl coordinate onto frozen ST050R-C.

Preregistered in #457 before evaluating this file.  This is an Agent-7 expression-
capacity diagnostic only: the stored ST050R-C velocity/pressure/force are not
changed, no PDE residual is inherited through the diagnostic transform, and no
OpenAI image value is fitted.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

import replay_recovery
from diagnostics import jets

TASK_ID = "CR003-ST050RC-ENERGY-NEUTRAL-SWIRL-TRANSFER-071"
PARENT_ID = "ST050R-C"
PARENT_HEAD = "dc1d5476ea9979566b774293add76196083288c5"
PREREG_ISSUE = 457
GAINS = (0.0, 0.015, 0.025, 0.035, 0.05)
TIMES = (0.25, 0.50, 0.75)
RADII = (0.6, 0.9, 1.2)
INNER_WINDOW = (0.30, 1.05)
OUTER_WINDOW = (0.95, 1.85)

CRITERIA = dict(
    inner_gain_floor=0.02,
    mid_gain_floor=0.01,
    outer_gain_ceiling=0.005,
    normalization_deviation_max=0.005,
    inward_speed_loss_max=0.005,
    axial_rms_abs_change_max=0.02,
    radial_rms_growth_max=0.03,
    collar_ratio_max=1.25,
)


def compact_bump(r, lo, hi):
    r = np.asarray(r, float)
    mid = 0.5 * (lo + hi); half = 0.5 * (hi - lo)
    x = (r - mid) / half
    out = np.zeros_like(r)
    m = np.abs(x) < 1.0
    xm = x[m]
    out[m] = np.exp(1.0 - 1.0 / (1.0 - xm * xm))
    return out


def compact_bump_derivative(r, lo, hi):
    r = np.asarray(r, float)
    mid = 0.5 * (lo + hi); half = 0.5 * (hi - lo)
    x = (r - mid) / half
    out = np.zeros_like(r)
    m = np.abs(x) < 1.0
    xm = x[m]
    g = np.exp(1.0 - 1.0 / (1.0 - xm * xm))
    out[m] = g * (-2.0 * xm / (1.0 - xm * xm) ** 2) / half
    return out


def h_profile(r, alpha):
    return compact_bump(r, *INNER_WINDOW) - alpha * compact_bump(r, *OUTER_WINDOW)


def h_derivative(r, alpha):
    return compact_bump_derivative(r, *INNER_WINDOW) - alpha * compact_bump_derivative(r, *OUTER_WINDOW)


def meridian_velocity(f, raw, r, z, t):
    r, z = np.broadcast_arrays(np.asarray(r, float), np.asarray(z, float))
    pts = np.column_stack((r.ravel(), np.zeros(r.size), z.ravel()))
    u, _ = f.fields(raw, pts, np.full(len(pts), float(t)))
    return np.asarray(u, float).reshape(r.shape + (3,))


def energy_and_moments(f, raw, order=72):
    x, wx = leggauss(order); y, wz = leggauss(order)
    r = x + 1.0; z = 2.0 * y
    R, Z = np.meshgrid(r, z, indexing="ij")
    U = meridian_velocity(f, raw, R, Z, 0.25)
    w = wx[:, None] * (2.0 * wz[None, :]) * np.pi * R
    u2 = np.sum(U * U, axis=-1)
    ut2 = U[..., 1] ** 2
    inner = float(np.sum(w * ut2 * compact_bump(R, *INNER_WINDOW)))
    outer = float(np.sum(w * ut2 * compact_bump(R, *OUTER_WINDOW)))
    if not (inner > 0.0 and outer > 0.0):
        raise RuntimeError("nonpositive swirl-energy moment")
    alpha = inner / outer
    defect = inner - alpha * outer
    parent_energy = float(np.sum(w * u2))
    return parent_energy, inner, outer, float(alpha), float(defect)


def raw_child_energy(f, raw, gain, alpha, order=72):
    x, wx = leggauss(order); y, wz = leggauss(order)
    r = x + 1.0; z = 2.0 * y
    R, Z = np.meshgrid(r, z, indexing="ij")
    U = meridian_velocity(f, raw, R, Z, 0.25)
    mult = 1.0 + float(gain) * h_profile(R, alpha)
    U = U.copy(); U[..., 1] *= mult
    w = wx[:, None] * (2.0 * wz[None, :]) * np.pi * R
    return float(np.sum(w * np.sum(U * U, axis=-1)))


def scale_for_gain(f, raw, gain, alpha, parent_energy, order=72):
    e = raw_child_energy(f, raw, gain, alpha, order=order)
    return float(np.sqrt(parent_energy / e)), e


def child_velocity(f, raw, points, t, gain, alpha, scale):
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
        ut *= 1.0 + gain * h_profile(r[mask], alpha)
        u[mask, 0] = rx * ur - ry * ut
        u[mask, 1] = ry * ur + rx * ut
    return scale * u


def vorticity_metrics(f, raw, gain, alpha, scale, t, nr=48, nz=72):
    x, wx = leggauss(nr); y, wz = leggauss(nz)
    S, Z = np.meshgrid(2.0 * (x + 1.0), 2.0 * y, indexing="ij")
    s, z = S.ravel(), Z.ravel(); r = np.sqrt(s)
    w = (4.0 * np.pi * np.outer(wx, wz)).ravel()
    j = jets(f, raw, r, z, float(t))
    ur, ut, uz = j[:, 0], j[:, 1], j[:, 2]
    wr, wt, wz0 = j[:, 3], j[:, 4], j[:, 5]
    m = 1.0 + gain * h_profile(r, alpha)
    mp = gain * h_derivative(r, alpha)
    wr1 = scale * m * wr
    wt1 = scale * wt
    wz1 = scale * (m * wz0 + mp * ut)
    om2 = wr1 * wr1 + wt1 * wt1 + wz1 * wz1
    en = float(w @ om2)
    axial = float(np.sqrt((w * z * z) @ om2 / en))
    radial = float(np.sqrt((w * s) @ om2 / en))
    collar_abs = float((w * om2 * ((r > 1.5) | (np.abs(z) > 1.5))) @ np.ones_like(om2))
    return dict(enstrophy=en, axial_rms=axial, radial_rms=radial, collar_absolute=collar_abs)


def sign_and_support_preflight(f, raw, gain, alpha, scale):
    core = []
    for r in (0.1, 0.3, 0.6):
        for z in (-0.5, -0.2, 0.2, 0.5):
            core.append((r, z))
    pts = np.array([[r, 0.0, z] for r, z in core], float)
    u = child_velocity(f, raw, pts, 0.5, gain, alpha, scale)
    bipolar = np.sign(pts[:, 2]) * u[:, 2] > 0
    signs = bool(np.all(u[:, 0] < 0) and np.all(u[:, 1] > 0) and np.all(bipolar))
    outside = np.array([[2.05,0,0],[-2.05,0,0],[0,0,2.05],[0,0,-2.05],[1.8,0,2.05]], float)
    uo = child_velocity(f, raw, outside, 0.5, gain, alpha, scale)
    support_max = float(np.max(np.abs(uo)))
    return signs, support_max


def divergence_fd(f, raw, gain, alpha, scale, seed=9175711, h=1e-5):
    rng = np.random.default_rng(seed)
    pts = rng.uniform(-1.25, 1.25, size=(24, 3))
    vals = np.zeros(len(pts))
    for axis in range(3):
        d = np.zeros(3); d[axis] = h
        up = child_velocity(f, raw, pts + d, 0.5, gain, alpha, scale)
        um = child_velocity(f, raw, pts - d, 0.5, gain, alpha, scale)
        vals += (up[:, axis] - um[:, axis]) / (2.0 * h)
    return float(np.max(np.abs(vals)))


def row_passes(row):
    c = CRITERIA
    return bool(
        row["angular_gain_r06"] >= c["inner_gain_floor"]
        and row["angular_gain_r09"] >= c["mid_gain_floor"]
        and row["angular_gain_r12"] <= c["outer_gain_ceiling"]
        and abs(row["normalization_relative_to_zero"] - 1.0) <= c["normalization_deviation_max"]
        and row["inward_speed_loss"] <= c["inward_speed_loss_max"]
        and max(abs(x) for x in row["axial_rms_relative_changes"]) <= c["axial_rms_abs_change_max"]
        and max(row["radial_rms_relative_changes"]) <= c["radial_rms_growth_max"]
        and max(row["collar_ratios"]) <= c["collar_ratio_max"]
        and row["core_signs_pass"] and row["support_max_abs"] <= 1e-12
        and row["divergence_fd_max"] <= 1e-5
    )


def run(out: Path):
    f, raw = replay_recovery.reconstruct(PARENT_ID)
    parent_energy, inner, outer, alpha, defect = energy_and_moments(f, raw)
    base_morph = {t: vorticity_metrics(f, raw, 0.0, alpha, 1.0, t) for t in TIMES}
    rows = []
    scale0, _ = scale_for_gain(f, raw, 0.0, alpha, parent_energy)
    for gain in GAINS:
        scale, raw_energy = scale_for_gain(f, raw, gain, alpha, parent_energy)
        angular = {r: float(scale / scale0 * (1.0 + gain * h_profile(np.array([r]), alpha)[0]) - 1.0) for r in RADII}
        morph = {t: vorticity_metrics(f, raw, gain, alpha, scale, t) for t in TIMES}
        axial = [morph[t]["axial_rms"] / base_morph[t]["axial_rms"] - 1.0 for t in TIMES]
        radial = [morph[t]["radial_rms"] / base_morph[t]["radial_rms"] - 1.0 for t in TIMES]
        collar = [morph[t]["collar_absolute"] / base_morph[t]["collar_absolute"] for t in TIMES]
        signs, support = sign_and_support_preflight(f, raw, gain, alpha, scale)
        div = divergence_fd(f, raw, gain, alpha, scale)
        row = dict(
            gain=float(gain), normalization=float(scale), normalization_relative_to_zero=float(scale/scale0),
            raw_energy=float(raw_energy), inward_speed_loss=float(max(0.0, 1.0-scale/scale0)),
            angular_gain_r06=angular[0.6], angular_gain_r09=angular[0.9], angular_gain_r12=angular[1.2],
            axial_rms_relative_changes=[float(x) for x in axial], radial_rms_relative_changes=[float(x) for x in radial],
            collar_ratios=[float(x) for x in collar], core_signs_pass=signs, support_max_abs=support,
            divergence_fd_max=div,
        )
        row["clean_transfer"] = row_passes(row) if gain > 0 else False
        rows.append(row)
    crossing = next((row for row in rows if row["clean_transfer"]), None)
    result = dict(
        task_id=TASK_ID, prereg_issue=PREREG_ISSUE, parent_id=PARENT_ID, parent_head=PARENT_HEAD,
        stored_parameter_count=2594, reduced_modifier_count=398, diagnostic_added_velocity_degrees=1,
        fixed_windows=dict(inner=list(INNER_WINDOW), outer=list(OUTER_WINDOW)), gains=list(GAINS), times=list(TIMES),
        balance=dict(parent_energy=parent_energy, inner_moment=inner, outer_moment=outer, alpha=alpha,
                     balanced_moment_defect=defect, kinetic_energy_first_derivative_at_zero=2.0*defect),
        criteria=CRITERIA, rows=rows,
        first_clean_crossing_gain=None if crossing is None else crossing["gain"],
        recommendation=("materialize exactly this one redistribution degree on an ST050R-C-derived governed child, then refit compatible pressure/restricted forcing and rerun fresh held-out full momentum" if crossing else "reject transfer; do not widen gain grid or add another basis in this increment"),
        canonical_velocity_changed=False, candidate_artifact_changed=False, pressure_or_force_changed=False,
        held_out_pde_residual_evaluated=False, pde_validated=False, source_correspondence_verified=False,
        paper_exact=False, openai_field_identified=False,
        scope="Target-free finite-grid expression-capacity transfer screen; no material-path integration and no visual-correspondence certificate.",
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"alpha":alpha,"first_clean_crossing_gain":result["first_clean_crossing_gain"],"rows":rows}, indent=2))
    return result


if __name__ == "__main__":
    p = argparse.ArgumentParser(); p.add_argument("--out", type=Path, required=True)
    a = p.parse_args(); run(a.out)
