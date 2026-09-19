"""Screen one smooth axial-conditioned radial Piola taper on callable ST051-B.

Preregistered in issue #519.  This Constrained Agent 7 increment adds exactly
one diagnostic geometry coordinate tau=.05 to the frozen callable ST051-B .025
swirl-redistribution child.  It does not scan coefficients, fit pressure or
forcing, or transfer any PDE receipt.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

import agent7_st051b_callable_piola_transfer as base
import replay_st051

TASK_ID = "CR003-ST051B-CALLABLE-RADIAL-TAPER-077"
PREREG_ISSUE = 519
TAPER_TAU = 0.05
TIMES = base.TIMES
ENERGY_ORDERS = base.ENERGY_ORDERS
VORTICITY_NR_NZ = base.VORTICITY_ORDER
VORTICITY_FD_STEP = base.VORTICITY_FD_STEP
DIVERGENCE_FD_STEP = base.DIVERGENCE_FD_STEP

CRITERIA = dict(
    tip_radial_rms_reduction_min=0.01,
    central_radial_rms_abs_change_max=0.01,
    global_axial_rms_loss_max=0.01,
    tip_absolute_enstrophy_ratio_min=0.90,
    support_collar_absolute_enstrophy_ratio_max=1.10,
    taper_scale_deviation_max=0.02,
    support_max_abs=1.0e-12,
    divergence_fd_max=1.0e-5,
    response_condition_max=5.0,
    response_abs_cosine_max=0.90,
)

TRUTH = dict(
    canonical_velocity_changed=False,
    saved_velocity_changed=False,
    production_taper_selected=False,
    production_beta_selected=False,
    production_candidate_selected=False,
    pressure_or_force_changed=False,
    held_out_pde_residual_evaluated=False,
    parent_pde_receipt_transferred=False,
    material_paths_integrated=False,
    public_image_used=False,
    visualization_ready=False,
    visual_correspondence_verified=False,
    pde_validated=False,
    source_correspondence_verified=False,
    paper_exact=False,
    openai_field_identified=False,
    blowup_proved=False,
)


def taper_profile_and_derivative(z):
    """Return q(z), dq/dz for q=(729/16)s^4(1-s^2)^4, s=z/2."""
    z = np.asarray(z, dtype=float)
    if not np.all(np.isfinite(z)):
        raise ValueError("z must be finite")
    s = z / base.AXIAL_SUPPORT
    inside = np.abs(s) < 1.0
    one_minus = np.maximum(1.0 - s * s, 0.0)
    k = 729.0 / 16.0
    q = np.where(inside, k * s**4 * one_minus**4, 0.0)
    # dq/dz = (1/2) dq/ds, with support half-width 2.
    dq = np.where(
        inside,
        2.0 * k * s**3 * one_minus**3 * (1.0 - 3.0 * s * s),
        0.0,
    )
    return q, dq


def taper_map(points, tau: float):
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError("points must have shape (n,3)")
    tau = float(tau)
    if not np.isfinite(tau) or abs(tau) > 0.15:
        raise ValueError("tau outside audited local taper range")
    q, dq = taper_profile_and_derivative(points[:, 2])
    a = 1.0 + tau * q
    ap = tau * dq
    if np.any(a <= 0.0) or not np.all(np.isfinite(a)) or not np.all(np.isfinite(ap)):
        raise RuntimeError("radial taper map lost orientation")
    mapped = np.array(points, copy=True)
    mapped[:, 0] *= a
    mapped[:, 1] *= a
    return mapped, a, ap


def taper_velocity(f, raw, points, time: float, *, redistribution_scale: float, tau: float, taper_scale: float):
    """Contravariant Piola pullback for F=(a(z)x,a(z)y,z)."""
    points = np.asarray(points, dtype=float)
    mapped, a, ap = taper_map(points, tau)
    u = base.redistributed_velocity(
        f,
        raw,
        mapped,
        time,
        gain=base.SOURCE_REDISTRIBUTION_GAIN,
        scale=redistribution_scale,
    )
    x = points[:, 0]
    y = points[:, 1]
    uz = u[:, 2].copy()
    out = np.empty_like(u)
    out[:, 0] = a * u[:, 0] - a * ap * x * uz
    out[:, 1] = a * u[:, 1] - a * ap * y * uz
    out[:, 2] = a * a * uz
    return float(taper_scale) * out


def axisymmetric_energy(velocity_fn, *, time: float, order: int):
    xr, wr = leggauss(int(order))
    xz, wz = leggauss(int(order))
    r = xr + 1.0
    z = 2.0 * xz
    wz = 2.0 * wz
    R, Z = np.meshgrid(r, z, indexing="ij")
    pts = np.column_stack((R.ravel(), np.zeros(R.size), Z.ravel()))
    u = np.asarray(velocity_fn(pts, float(time)), dtype=float).reshape(R.shape + (3,))
    return float(np.sum((wr[:, None] * wz[None, :]) * (np.pi * R * np.sum(u * u, axis=-1))))


def energy_rows(f, raw):
    rows = []
    for order in ENERGY_ORDERS:
        parent_e = axisymmetric_energy(
            lambda p, t: base.parent_velocity(f, raw, p, t), time=0.25, order=order
        )
        red_raw_e = axisymmetric_energy(
            lambda p, t: base.redistributed_velocity(
                f, raw, p, t, gain=base.SOURCE_REDISTRIBUTION_GAIN, scale=1.0
            ),
            time=0.25,
            order=order,
        )
        red_scale = float(np.sqrt(parent_e / red_raw_e))
        red_e = axisymmetric_energy(
            lambda p, t: base.redistributed_velocity(
                f, raw, p, t, gain=base.SOURCE_REDISTRIBUTION_GAIN, scale=red_scale
            ),
            time=0.25,
            order=order,
        )
        taper_raw_e = axisymmetric_energy(
            lambda p, t: taper_velocity(
                f, raw, p, t, redistribution_scale=red_scale, tau=TAPER_TAU, taper_scale=1.0
            ),
            time=0.25,
            order=order,
        )
        taper_scale = float(np.sqrt(red_e / taper_raw_e))
        final_e = axisymmetric_energy(
            lambda p, t: taper_velocity(
                f, raw, p, t,
                redistribution_scale=red_scale,
                tau=TAPER_TAU,
                taper_scale=taper_scale,
            ),
            time=0.25,
            order=order,
        )
        rows.append(dict(
            order=int(order),
            parent_energy=parent_e,
            redistribution_scale=red_scale,
            redistribution_energy=red_e,
            taper_raw_energy=taper_raw_e,
            taper_scale=taper_scale,
            final_energy=final_e,
            final_relative_to_redistribution=float(final_e / red_e - 1.0),
        ))
    return rows


def morphology(velocity_fn, *, time: float):
    xr, wr = leggauss(int(VORTICITY_NR_NZ[0]))
    xz, wz = leggauss(int(VORTICITY_NR_NZ[1]))
    r = xr + 1.0
    z = 2.0 * xz
    wz = 2.0 * wz
    R, Z = np.meshgrid(r, z, indexing="ij")
    pts = np.column_stack((R.ravel(), np.zeros(R.size), Z.ravel()))
    omega = base.cartesian_curl(velocity_fn, pts, time=time, h=VORTICITY_FD_STEP)
    om2 = np.sum(omega * omega, axis=1).reshape(R.shape)
    weight = 2.0 * np.pi * (wr[:, None] * wz[None, :]) * R
    total = float(np.sum(weight * om2))
    if total <= np.finfo(float).tiny:
        raise RuntimeError("inactive vorticity")

    tip = (np.abs(Z) >= 0.45 * base.AXIAL_SUPPORT) & (np.abs(Z) <= 0.80 * base.AXIAL_SUPPORT)
    central = np.abs(Z) <= 0.35 * base.AXIAL_SUPPORT
    collar = (R > 0.75 * base.RADIAL_SUPPORT) | (np.abs(Z) > 0.75 * base.AXIAL_SUPPORT)

    def band(mask):
        e = float(np.sum(weight * om2 * mask))
        if e <= np.finfo(float).tiny:
            raise RuntimeError("inactive morphology band")
        rr = float(np.sqrt(np.sum(weight * R * R * om2 * mask) / e))
        return e, rr

    tip_e, tip_rr = band(tip)
    central_e, central_rr = band(central)
    return dict(
        enstrophy=total,
        axial_rms=float(np.sqrt(np.sum(weight * Z * Z * om2) / total)),
        radial_rms=float(np.sqrt(np.sum(weight * R * R * om2) / total)),
        tip_absolute_enstrophy=tip_e,
        tip_radial_rms=tip_rr,
        central_absolute_enstrophy=central_e,
        central_radial_rms=central_rr,
        support_collar_absolute_enstrophy=float(np.sum(weight * om2 * collar)),
    )


def morphology_rows(f, raw, red_scale: float, taper_scale: float):
    control = lambda p, t: base.redistributed_velocity(
        f, raw, p, t, gain=base.SOURCE_REDISTRIBUTION_GAIN, scale=red_scale
    )
    child = lambda p, t: taper_velocity(
        f, raw, p, t, redistribution_scale=red_scale, tau=TAPER_TAU, taper_scale=taper_scale
    )
    rows = []
    for time in TIMES:
        b = morphology(control, time=time)
        c = morphology(child, time=time)
        rows.append(dict(
            time=float(time),
            control=b,
            taper=c,
            tip_radial_rms_relative=float(c["tip_radial_rms"] / b["tip_radial_rms"] - 1.0),
            central_radial_rms_relative=float(c["central_radial_rms"] / b["central_radial_rms"] - 1.0),
            global_axial_rms_relative=float(c["axial_rms"] / b["axial_rms"] - 1.0),
            global_radial_rms_relative=float(c["radial_rms"] / b["radial_rms"] - 1.0),
            tip_absolute_enstrophy_ratio=float(c["tip_absolute_enstrophy"] / b["tip_absolute_enstrophy"]),
            central_absolute_enstrophy_ratio=float(c["central_absolute_enstrophy"] / b["central_absolute_enstrophy"]),
            collar_absolute_ratio=float(c["support_collar_absolute_enstrophy"] / b["support_collar_absolute_enstrophy"]),
        ))
    return rows


def structure_preflight(f, raw, red_scale: float, taper_scale: float):
    child = lambda p, t: taper_velocity(
        f, raw, p, t, redistribution_scale=red_scale, tau=TAPER_TAU, taper_scale=taper_scale
    )
    core = np.array(
        [[r, 0.0, z] for r in (0.10, 0.30, 0.60) for z in (-0.50, -0.20, 0.20, 0.50)],
        dtype=float,
    )
    uc = child(core, 0.50)
    signs = bool(
        np.all(uc[:, 0] < 0.0)
        and np.all(uc[:, 1] > 0.0)
        and np.all(core[:, 2] * uc[:, 2] > 0.0)
    )
    outside = np.array(
        [[2.05,0,0],[-2.05,0,0],[0,2.05,0],[0,-2.05,0],[0,0,2.05],[0,0,-2.05],[1.8,0,2.05]],
        dtype=float,
    )
    support_max = max(float(np.max(np.abs(child(outside, t)))) for t in TIMES)

    rng = np.random.default_rng(9175771)
    points = rng.uniform(-1.15, 1.15, size=(36, 3))
    div = np.zeros(len(points), dtype=float)
    for axis in range(3):
        delta = np.zeros_like(points)
        delta[:, axis] = DIVERGENCE_FD_STEP
        up = child(points + delta, 0.50)
        um = child(points - delta, 0.50)
        div += (up[:, axis] - um[:, axis]) / (2.0 * DIVERGENCE_FD_STEP)

    z = np.linspace(-base.AXIAL_SUPPORT, base.AXIAL_SUPPORT, 8001)
    q, _ = taper_profile_and_derivative(z)
    a = 1.0 + TAPER_TAU * q
    return dict(
        core_signs_pass=signs,
        support_max_abs=support_max,
        divergence_fd_max=float(np.max(np.abs(div))),
        radial_scale_min=float(np.min(a)),
        radial_scale_max=float(np.max(a)),
        determinant_min=float(np.min(a * a)),
        determinant_max=float(np.max(a * a)),
        midplane_identity_error=float(abs(a[len(a)//2] - 1.0)),
        support_endpoint_identity_error=float(max(abs(a[0] - 1.0), abs(a[-1] - 1.0))),
    )


def response_diagnostics(f, raw, red_scale: float):
    theta = 0.37
    pts = []
    ts = []
    for time in TIMES:
        for r in (0.35, 0.75, 1.15):
            for z in (-0.65, -0.25, 0.25, 0.65):
                pts.append([r * np.cos(theta), r * np.sin(theta), z])
                ts.append(float(time))
    pts = np.asarray(pts, dtype=float)
    ts = np.asarray(ts, dtype=float)
    eps = 1.0e-4

    def stack(fn):
        chunks = []
        for time in TIMES:
            mask = ts == time
            chunks.append(np.asarray(fn(pts[mask], time), dtype=float).ravel())
        return np.concatenate(chunks)

    gain_plus = stack(lambda p, t: base.redistributed_velocity(
        f, raw, p, t, gain=base.SOURCE_REDISTRIBUTION_GAIN + eps, scale=red_scale
    ))
    gain_minus = stack(lambda p, t: base.redistributed_velocity(
        f, raw, p, t, gain=base.SOURCE_REDISTRIBUTION_GAIN - eps, scale=red_scale
    ))
    col_gain = (gain_plus - gain_minus) / (2.0 * eps)

    beta_plus = stack(lambda p, t: base.piola_velocity(
        f, raw, p, t, redistribution_scale=red_scale, beta=eps, piola_scale=1.0
    ))
    beta_minus = stack(lambda p, t: base.piola_velocity(
        f, raw, p, t, redistribution_scale=red_scale, beta=-eps, piola_scale=1.0
    ))
    col_beta = (beta_plus - beta_minus) / (2.0 * eps)

    taper_plus = stack(lambda p, t: taper_velocity(
        f, raw, p, t, redistribution_scale=red_scale, tau=eps, taper_scale=1.0
    ))
    taper_minus = stack(lambda p, t: taper_velocity(
        f, raw, p, t, redistribution_scale=red_scale, tau=-eps, taper_scale=1.0
    ))
    col_taper = (taper_plus - taper_minus) / (2.0 * eps)

    cols = [col_gain, col_beta, col_taper]
    norms = np.array([np.linalg.norm(c) for c in cols])
    if np.any(norms <= np.finfo(float).tiny):
        raise RuntimeError("degenerate response column")
    matrix = np.column_stack([c / n for c, n in zip(cols, norms)])
    singular = np.linalg.svd(matrix, compute_uv=False)
    gram = matrix.T @ matrix
    cosines = {
        "gain_beta": float(gram[0, 1]),
        "gain_taper": float(gram[0, 2]),
        "beta_taper": float(gram[1, 2]),
    }
    return dict(
        probe_count=int(len(pts)),
        normalized_rank=int(np.linalg.matrix_rank(matrix, tol=1.0e-10)),
        singular_values=singular.tolist(),
        condition_number=float(singular[0] / singular[-1]),
        pairwise_cosines=cosines,
        max_abs_pairwise_cosine=float(max(abs(v) for v in cosines.values())),
        raw_column_norms=norms.tolist(),
    )


def clean_rule(energy, morph, structure, response):
    c = CRITERIA
    finest = energy[-1]
    return bool(
        abs(finest["taper_scale"] - 1.0) <= c["taper_scale_deviation_max"]
        and abs(finest["final_relative_to_redistribution"]) <= 1.0e-3
        and all(row["tip_radial_rms_relative"] <= -c["tip_radial_rms_reduction_min"] for row in morph)
        and all(abs(row["central_radial_rms_relative"]) <= c["central_radial_rms_abs_change_max"] for row in morph)
        and all(row["global_axial_rms_relative"] >= -c["global_axial_rms_loss_max"] for row in morph)
        and all(row["tip_absolute_enstrophy_ratio"] >= c["tip_absolute_enstrophy_ratio_min"] for row in morph)
        and all(row["collar_absolute_ratio"] <= c["support_collar_absolute_enstrophy_ratio_max"] for row in morph)
        and structure["determinant_min"] > 0.0
        and structure["core_signs_pass"]
        and structure["support_max_abs"] <= c["support_max_abs"]
        and structure["divergence_fd_max"] <= c["divergence_fd_max"]
        and response["normalized_rank"] == 3
        and response["condition_number"] <= c["response_condition_max"]
        and response["max_abs_pairwise_cosine"] <= c["response_abs_cosine_max"]
    )


def run(out: Path):
    f, raw = replay_st051.reconstruct(base.PARENT_ID)
    energy = energy_rows(f, raw)
    red_scale = float(energy[-1]["redistribution_scale"])
    taper_scale = float(energy[-1]["taper_scale"])
    morph = morphology_rows(f, raw, red_scale, taper_scale)
    structure = structure_preflight(f, raw, red_scale, taper_scale)
    response = response_diagnostics(f, raw, red_scale)
    clean = clean_rule(energy, morph, structure, response)

    result = dict(
        task_id=TASK_ID,
        prereg_issue=PREREG_ISSUE,
        parent=dict(id=base.PARENT_ID, exact_head=base.PARENT_HEAD),
        frozen_redistribution=dict(
            alpha=base.SOURCE_REDISTRIBUTION_ALPHA,
            gain=base.SOURCE_REDISTRIBUTION_GAIN,
            inner_window=list(base.INNER_WINDOW),
            outer_window=list(base.OUTER_WINDOW),
        ),
        taper=dict(
            tau=TAPER_TAU,
            tau_scan_performed=False,
            added_geometry_degrees=1,
            callable_parent_used=True,
            formula="F_tau=(a(z)x,a(z)y,z), a=1+tau*(729/16)*(z/2)^4*(1-(z/2)^2)^4; contravariant Piola pullback",
        ),
        criteria=CRITERIA,
        energy_quadrature=energy,
        morphology=morph,
        structure=structure,
        response_diagnostics=response,
        clean_taper_capacity=clean,
        recommendation=(
            "The single fixed off-midplane radial Piola taper is useful independent morphology capacity. Keep tau=.05 diagnostic-only; do not grow another taper/ring family before a materialized/rendered child tests whether this channel is actually needed."
            if clean else
            "The single fixed off-midplane radial Piola taper does not satisfy the preregistered clean-capacity rule. Do not widen tau or add another taper family in this increment; retain the existing redistribution plus axial-Piola block and route from the failed component diagnostics."
        ),
        direct_contribution_to_final_velocity=("positive_expression_capacity_only" if clean else "no_positive_transfer"),
        truth=dict(TRUTH),
        scope="One fixed callable radial-taper Piola capacity test; no production, visual-correspondence, or PDE claim.",
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    run(args.out)


if __name__ == "__main__":
    main()
