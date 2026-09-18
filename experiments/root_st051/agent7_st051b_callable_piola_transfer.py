"""Replay one fixed axial Piola degree on the callable ST051-B .025 swirl child.

Preregistered in issue #509 before evaluation.  This is Constrained Agent 7's
representation-capacity lane.  It closes one numerical-representation seam in
PR #499 by applying the inherited beta=.075 Piola map directly to the callable
ST051-B field reconstructed from PR #460, rather than interpolating a frozen
33^3 visualization grid.

The experiment adds no spatial basis family and performs no coefficient scan.
It does not fit pressure/forcing or transfer any parent PDE receipt.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

import replay_st051

TASK_ID = "CR003-ST051B-CALLABLE-PIOLA-TRANSFER-076"
PREREG_ISSUE = 509
PARENT_ID = "ST051-B"
PARENT_HEAD = "4b784f1b8457af2ead49295631d834d4e882000b"
PARENT_RAW_SHA256 = "0071106ef10a5d77b620b942becc00b58c8a6765fd5fb7bc4dcd195ac65c970d"
SOURCE_REDISTRIBUTION_ALPHA = 2.520520814687742
SOURCE_REDISTRIBUTION_GAIN = 0.025
SOURCE_REDISTRIBUTION_SCALE = 1.0014791925672812
INNER_WINDOW = (0.30, 1.05)
OUTER_WINDOW = (0.95, 1.85)
PIOLA_BETA = 0.075
RADIAL_SUPPORT = 2.0
AXIAL_SUPPORT = 2.0
TIMES = (0.25, 0.50, 0.75)
ENERGY_ORDERS = (36, 52, 68)
VORTICITY_ORDER = (28, 40)
VORTICITY_FD_STEP = 5.0e-4
VORTICITY_FD_STEP_COARSE = 1.0e-3
DIVERGENCE_FD_STEP = 1.0e-5

CRITERIA = dict(
    redistribution_scale_replay_abs_max=5.0e-6,
    reference_energy_relative_error_max=1.0e-3,
    piola_scale_deviation_max=0.02,
    axial_vorticity_rms_gain_min=0.01,
    radial_vorticity_rms_growth_max=0.03,
    support_collar_absolute_enstrophy_ratio_max=1.10,
    support_max_abs=1.0e-12,
    divergence_fd_max=1.0e-5,
    response_condition_max=5.0,
    response_abs_cosine_max=0.90,
)

TRUTH = dict(
    canonical_velocity_changed=False,
    saved_velocity_changed=False,
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


def compact_bump(r, lo: float, hi: float):
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
    return compact_bump(r, *INNER_WINDOW) - SOURCE_REDISTRIBUTION_ALPHA * compact_bump(r, *OUTER_WINDOW)


def warp_z_and_jacobian(z, beta: float = PIOLA_BETA):
    z = np.asarray(z, dtype=float)
    beta = float(beta)
    if not np.isfinite(beta) or abs(beta) > 0.20:
        raise ValueError("beta outside audited local Piola range")
    if not np.all(np.isfinite(z)):
        raise ValueError("z must be finite")
    s = z / AXIAL_SUPPORT
    inside = np.abs(s) < 1.0
    one_minus = np.maximum(1.0 - s * s, 0.0)
    bump = one_minus**4
    mapped = np.where(inside, z * (1.0 - beta * bump), z)
    jac_inside = 1.0 - beta * bump + 8.0 * beta * s * s * one_minus**3
    jac = np.where(inside, jac_inside, 1.0)
    if not np.all(np.isfinite(mapped)) or not np.all(np.isfinite(jac)) or np.any(jac <= 0.0):
        raise RuntimeError("Piola coordinate map lost orientation")
    return mapped, jac


def parent_velocity(f, raw, points, time: float):
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError("points must have shape (n,3)")
    u, _ = f.fields(raw, points, np.full(len(points), float(time)))
    u = np.asarray(u, dtype=float)
    if u.shape != points.shape or not np.all(np.isfinite(u)):
        raise RuntimeError("unexpected callable ST051-B velocity")
    return u


def redistributed_velocity(f, raw, points, time: float, *, gain: float, scale: float):
    points = np.asarray(points, dtype=float)
    u = parent_velocity(f, raw, points, time).copy()
    x, y = points[:, 0], points[:, 1]
    r = np.hypot(x, y)
    mask = r > 1.0e-14
    if np.any(mask):
        rx, ry = x[mask] / r[mask], y[mask] / r[mask]
        ux, uy = u[mask, 0].copy(), u[mask, 1].copy()
        ur = rx * ux + ry * uy
        ut = -ry * ux + rx * uy
        ut *= 1.0 + float(gain) * h_profile(r[mask])
        u[mask, 0] = rx * ur - ry * ut
        u[mask, 1] = ry * ur + rx * ut
    return float(scale) * u


def piola_velocity(
    f,
    raw,
    points,
    time: float,
    *,
    redistribution_gain: float = SOURCE_REDISTRIBUTION_GAIN,
    redistribution_scale: float,
    beta: float,
    piola_scale: float,
):
    points = np.asarray(points, dtype=float)
    mapped_z, jac = warp_z_and_jacobian(points[:, 2], beta)
    mapped = np.array(points, copy=True)
    mapped[:, 2] = mapped_z
    u = redistributed_velocity(
        f,
        raw,
        mapped,
        time,
        gain=redistribution_gain,
        scale=redistribution_scale,
    ).copy()
    u[:, 0] *= jac
    u[:, 1] *= jac
    return float(piola_scale) * u


def axisymmetric_energy(velocity_fn, *, time: float, order: int):
    node_r, weight_r = leggauss(int(order))
    node_z, weight_z = leggauss(int(order))
    r = node_r + 1.0
    z = 2.0 * node_z
    wr = weight_r
    wz = 2.0 * weight_z
    R, Z = np.meshgrid(r, z, indexing="ij")
    pts = np.column_stack((R.ravel(), np.zeros(R.size), Z.ravel()))
    u = np.asarray(velocity_fn(pts, float(time)), dtype=float).reshape(R.shape + (3,))
    # E = (1/2) int |u|^2 dV = pi int |u|^2 r dr dz for axisymmetry.
    return float(np.sum((wr[:, None] * wz[None, :]) * (np.pi * R * np.sum(u * u, axis=-1))))


def energy_scales(f, raw):
    rows = []
    for order in ENERGY_ORDERS:
        parent_e = axisymmetric_energy(lambda p, t: parent_velocity(f, raw, p, t), time=0.25, order=order)
        red_raw_e = axisymmetric_energy(
            lambda p, t: redistributed_velocity(f, raw, p, t, gain=SOURCE_REDISTRIBUTION_GAIN, scale=1.0),
            time=0.25,
            order=order,
        )
        red_scale = float(np.sqrt(parent_e / red_raw_e))
        red_e = axisymmetric_energy(
            lambda p, t: redistributed_velocity(f, raw, p, t, gain=SOURCE_REDISTRIBUTION_GAIN, scale=red_scale),
            time=0.25,
            order=order,
        )
        piola_raw_e = axisymmetric_energy(
            lambda p, t: piola_velocity(
                f, raw, p, t,
                redistribution_scale=red_scale,
                beta=PIOLA_BETA,
                piola_scale=1.0,
            ),
            time=0.25,
            order=order,
        )
        piola_scale = float(np.sqrt(red_e / piola_raw_e))
        final_e = axisymmetric_energy(
            lambda p, t: piola_velocity(
                f, raw, p, t,
                redistribution_scale=red_scale,
                beta=PIOLA_BETA,
                piola_scale=piola_scale,
            ),
            time=0.25,
            order=order,
        )
        rows.append(dict(
            order=int(order),
            parent_energy=parent_e,
            redistribution_raw_energy=red_raw_e,
            redistribution_scale=red_scale,
            redistribution_energy=red_e,
            piola_raw_energy=piola_raw_e,
            piola_scale=piola_scale,
            final_energy=final_e,
            final_relative_to_redistribution=float(final_e / red_e - 1.0),
        ))
    return rows


def cartesian_curl(velocity_fn, points, *, time: float, h: float):
    points = np.asarray(points, dtype=float)
    deriv = np.empty((len(points), 3, 3), dtype=float)
    for axis in range(3):
        delta = np.zeros_like(points)
        delta[:, axis] = float(h)
        up = np.asarray(velocity_fn(points + delta, float(time)), dtype=float)
        um = np.asarray(velocity_fn(points - delta, float(time)), dtype=float)
        deriv[:, :, axis] = (up - um) / (2.0 * float(h))
    return np.column_stack((
        deriv[:, 2, 1] - deriv[:, 1, 2],
        deriv[:, 0, 2] - deriv[:, 2, 0],
        deriv[:, 1, 0] - deriv[:, 0, 1],
    ))


def vorticity_morphology(velocity_fn, *, time: float, nr: int, nz: int, h: float):
    xr, wr = leggauss(int(nr))
    xz, wz = leggauss(int(nz))
    r = xr + 1.0
    z = 2.0 * xz
    wz = 2.0 * wz
    R, Z = np.meshgrid(r, z, indexing="ij")
    pts = np.column_stack((R.ravel(), np.zeros(R.size), Z.ravel()))
    omega = cartesian_curl(velocity_fn, pts, time=time, h=h)
    om2 = np.sum(omega * omega, axis=1).reshape(R.shape)
    weight = 2.0 * np.pi * (wr[:, None] * wz[None, :]) * R
    total = float(np.sum(weight * om2))
    if not np.isfinite(total) or total <= np.finfo(float).tiny:
        raise RuntimeError("inactive vorticity quadrature")
    collar = (R > 0.75 * RADIAL_SUPPORT) | (np.abs(Z) > 0.75 * AXIAL_SUPPORT)
    return dict(
        enstrophy=total,
        axial_rms=float(np.sqrt(np.sum(weight * Z * Z * om2) / total)),
        radial_rms=float(np.sqrt(np.sum(weight * R * R * om2) / total)),
        support_collar_absolute_enstrophy=float(np.sum(weight * om2 * collar)),
        outer_065_fraction=float(np.sum(weight * om2 * (np.abs(Z) > 0.65 * AXIAL_SUPPORT)) / total),
        outer_075_fraction=float(np.sum(weight * om2 * (np.abs(Z) > 0.75 * AXIAL_SUPPORT)) / total),
    )


def morphology_block(f, raw, red_scale: float, piola_scale: float, *, h: float):
    base_fn = lambda p, t: redistributed_velocity(
        f, raw, p, t, gain=SOURCE_REDISTRIBUTION_GAIN, scale=red_scale
    )
    child_fn = lambda p, t: piola_velocity(
        f, raw, p, t,
        redistribution_scale=red_scale,
        beta=PIOLA_BETA,
        piola_scale=piola_scale,
    )
    rows = []
    for time in TIMES:
        base = vorticity_morphology(base_fn, time=time, nr=VORTICITY_ORDER[0], nz=VORTICITY_ORDER[1], h=h)
        child = vorticity_morphology(child_fn, time=time, nr=VORTICITY_ORDER[0], nz=VORTICITY_ORDER[1], h=h)
        rows.append(dict(
            time=float(time),
            base=base,
            child=child,
            axial_rms_relative=float(child["axial_rms"] / base["axial_rms"] - 1.0),
            radial_rms_relative=float(child["radial_rms"] / base["radial_rms"] - 1.0),
            collar_absolute_ratio=float(
                child["support_collar_absolute_enstrophy"] /
                max(base["support_collar_absolute_enstrophy"], 1.0e-300)
            ),
            outer_065_relative=float(child["outer_065_fraction"] / max(base["outer_065_fraction"], 1.0e-300) - 1.0),
            outer_075_relative=float(child["outer_075_fraction"] / max(base["outer_075_fraction"], 1.0e-300) - 1.0),
        ))
    return rows


def structure_preflight(f, raw, red_scale: float, piola_scale: float):
    child = lambda p, t: piola_velocity(
        f, raw, p, t,
        redistribution_scale=red_scale,
        beta=PIOLA_BETA,
        piola_scale=piola_scale,
    )
    core = np.array(
        [[r, 0.0, z] for r in (0.10, 0.30, 0.60) for z in (-0.50, -0.20, 0.20, 0.50)],
        dtype=float,
    )
    uc = child(core, 0.50)
    signs = bool(np.all(uc[:, 0] < 0.0) and np.all(uc[:, 1] > 0.0) and np.all(core[:, 2] * uc[:, 2] > 0.0))
    outside = np.array(
        [[2.05,0,0],[-2.05,0,0],[0,2.05,0],[0,-2.05,0],[0,0,2.05],[0,0,-2.05],[1.8,0,2.05]],
        dtype=float,
    )
    support_max = max(float(np.max(np.abs(child(outside, t)))) for t in TIMES)

    rng = np.random.default_rng(9175761)
    points = rng.uniform(-1.15, 1.15, size=(36, 3))
    div = np.zeros(len(points), dtype=float)
    for axis in range(3):
        delta = np.zeros_like(points)
        delta[:, axis] = DIVERGENCE_FD_STEP
        up = child(points + delta, 0.50)
        um = child(points - delta, 0.50)
        div += (up[:, axis] - um[:, axis]) / (2.0 * DIVERGENCE_FD_STEP)

    z = np.linspace(-AXIAL_SUPPORT, AXIAL_SUPPORT, 8001)
    mapped, jac = warp_z_and_jacobian(z, PIOLA_BETA)
    return dict(
        core_signs_pass=signs,
        support_max_abs=support_max,
        divergence_fd_max=float(np.max(np.abs(div))),
        coordinate_jacobian_min=float(np.min(jac)),
        coordinate_jacobian_max=float(np.max(jac)),
        coordinate_map_monotone=bool(np.all(np.diff(mapped) > 0.0)),
        support_endpoint_map_error=float(max(abs(mapped[0] + 2.0), abs(mapped[-1] - 2.0))),
    )


def response_diagnostics(f, raw, red_scale: float):
    theta = 0.37
    points = []
    times = []
    for time in TIMES:
        for r in (0.35, 0.75, 1.15):
            for z in (-0.65, -0.25, 0.25, 0.65):
                points.append([r * np.cos(theta), r * np.sin(theta), z])
                times.append(float(time))
    points = np.asarray(points, dtype=float)
    times = np.asarray(times, dtype=float)
    eps = 1.0e-4

    def stacked_gain(gain):
        chunks = []
        for t in TIMES:
            mask = times == t
            chunks.append(redistributed_velocity(f, raw, points[mask], t, gain=gain, scale=red_scale).ravel())
        return np.concatenate(chunks)

    col_gain = (stacked_gain(SOURCE_REDISTRIBUTION_GAIN + eps) - stacked_gain(SOURCE_REDISTRIBUTION_GAIN - eps)) / (2.0 * eps)

    chunks0 = []
    chunks1 = []
    for t in TIMES:
        mask = times == t
        p = points[mask]
        chunks0.append(redistributed_velocity(f, raw, p, t, gain=SOURCE_REDISTRIBUTION_GAIN, scale=red_scale).ravel())
        chunks1.append(piola_velocity(
            f, raw, p, t,
            redistribution_scale=red_scale,
            beta=eps,
            piola_scale=1.0,
        ).ravel())
    col_piola = (np.concatenate(chunks1) - np.concatenate(chunks0)) / eps

    norms = np.array([np.linalg.norm(col_gain), np.linalg.norm(col_piola)], dtype=float)
    if np.any(norms <= np.finfo(float).tiny):
        raise RuntimeError("degenerate sensitivity column")
    matrix = np.column_stack((col_gain / norms[0], col_piola / norms[1]))
    singular = np.linalg.svd(matrix, compute_uv=False)
    rank = int(np.linalg.matrix_rank(matrix, tol=1.0e-10))
    cosine = float(np.dot(matrix[:, 0], matrix[:, 1]))
    condition = float(singular[0] / singular[-1])
    return dict(
        probe_count=int(len(points)),
        spacetime_velocity_scalar_count=int(matrix.shape[0]),
        finite_difference_step=eps,
        normalization="each sensitivity column normalized to unit Euclidean norm; common energy scales held fixed",
        raw_column_norms=norms.tolist(),
        normalized_rank=rank,
        singular_values=singular.tolist(),
        condition_number=condition,
        column_cosine=cosine,
        piola_fraction_outside_redistribution_span=float(np.sqrt(max(0.0, 1.0 - cosine * cosine))),
    )


def clean_rule(energy_rows, morphology, structure, response):
    c = CRITERIA
    finest = energy_rows[-1]
    return bool(
        abs(finest["redistribution_scale"] - SOURCE_REDISTRIBUTION_SCALE) <= c["redistribution_scale_replay_abs_max"]
        and abs(finest["final_relative_to_redistribution"]) <= c["reference_energy_relative_error_max"]
        and abs(finest["piola_scale"] - 1.0) <= c["piola_scale_deviation_max"]
        and all(row["axial_rms_relative"] >= c["axial_vorticity_rms_gain_min"] for row in morphology)
        and all(row["radial_rms_relative"] <= c["radial_vorticity_rms_growth_max"] for row in morphology)
        and all(row["collar_absolute_ratio"] <= c["support_collar_absolute_enstrophy_ratio_max"] for row in morphology)
        and structure["coordinate_map_monotone"]
        and structure["coordinate_jacobian_min"] > 0.0
        and structure["core_signs_pass"]
        and structure["support_max_abs"] <= c["support_max_abs"]
        and structure["divergence_fd_max"] <= c["divergence_fd_max"]
        and response["normalized_rank"] == 2
        and response["condition_number"] <= c["response_condition_max"]
        and abs(response["column_cosine"]) <= c["response_abs_cosine_max"]
    )


def run(out: Path):
    f, raw = replay_st051.reconstruct(PARENT_ID)
    energy = energy_scales(f, raw)
    red_scale = float(energy[-1]["redistribution_scale"])
    piola_scale = float(energy[-1]["piola_scale"])

    morphology = morphology_block(f, raw, red_scale, piola_scale, h=VORTICITY_FD_STEP)
    morphology_coarse_fd = morphology_block(f, raw, red_scale, piola_scale, h=VORTICITY_FD_STEP_COARSE)
    structure = structure_preflight(f, raw, red_scale, piola_scale)
    response = response_diagnostics(f, raw, red_scale)
    clean = clean_rule(energy, morphology, structure, response)

    fd_refinement = []
    for fine, coarse in zip(morphology, morphology_coarse_fd):
        fd_refinement.append(dict(
            time=fine["time"],
            axial_gain_fine=fine["axial_rms_relative"],
            axial_gain_coarse=coarse["axial_rms_relative"],
            axial_gain_absolute_drift=float(abs(fine["axial_rms_relative"] - coarse["axial_rms_relative"])),
            radial_gain_absolute_drift=float(abs(fine["radial_rms_relative"] - coarse["radial_rms_relative"])),
        ))

    result = dict(
        task_id=TASK_ID,
        prereg_issue=PREREG_ISSUE,
        parent=dict(id=PARENT_ID, exact_head=PARENT_HEAD, reported_raw_candidate_sha256=PARENT_RAW_SHA256),
        frozen_redistribution=dict(
            alpha=SOURCE_REDISTRIBUTION_ALPHA,
            gain=SOURCE_REDISTRIBUTION_GAIN,
            inner_window=list(INNER_WINDOW),
            outer_window=list(OUTER_WINDOW),
            previously_audited_scale=SOURCE_REDISTRIBUTION_SCALE,
        ),
        piola=dict(
            beta=PIOLA_BETA,
            beta_scan_performed=False,
            added_geometry_degrees=1,
            callable_parent_used=True,
            sampled_grid_interpolation_used=False,
            formula="h_beta(z)=z*(1-beta*(1-(z/2)^2)^4); contravariant Piola pullback",
        ),
        criteria=CRITERIA,
        energy_quadrature=energy,
        morphology=morphology,
        morphology_fd_refinement=fd_refinement,
        structure=structure,
        response_diagnostics=response,
        clean_callable_transfer=clean,
        recommendation=(
            "The fixed beta=.075 axial Piola degree survives direct callable-parent replay on the frozen ST051-B .025 redistribution child. Treat it as one independent axial-geometry capacity coordinate; keep generic axial/swirl basis growth frozen and require a newly materialized child plus compatible pressure/restricted forcing and fresh held-out full momentum before any PDE promotion."
            if clean else
            "The fixed beta=.075 Piola degree does not satisfy the preregistered callable-parent transfer rule. Do not widen beta or add another axial basis in this increment; retain the callable ST051-B redistribution child and use the failed diagnostics to route the next minimal representation test."
        ),
        direct_contribution_to_final_velocity=(
            "positive_expression_capacity_only" if clean else "no_positive_transfer"
        ),
        truth=dict(TRUTH),
        scope="One fixed direct-callable Piola transfer on one frozen redistribution child; no production/PDE/visual-correspondence claim.",
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
