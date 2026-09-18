"""Screen one temporal coefficient on the frozen ST051-B swirl redistribution.

Preregistered in issue #488 before evaluation.  This Agent-7 increment adds no
new spatial basis.  It keeps the frozen redistribution profile and base gain
from #469/#480, then tests one affine-in-time coefficient

    k(t) = k0 + gamma * 2 * (t - .25),

so the t=.25 field and reference-energy normalization stay exactly frozen.
This is target-free expression-capacity evidence only: it is not a production
coefficient selection, image fit, material-path experiment, or PDE validation.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import agent7_st051b_frozen_redistribution_transfer as base
import replay_st051

TASK_ID = "CR003-ST051B-TEMPORAL-REDISTRIBUTION-074"
PREREG_ISSUE = 488
PARENT_ID = "ST051-B"
PARENT_HEAD = base.PARENT_HEAD
STACK_BASE_PR = 480
STACK_BASE_HEAD = "187157d377b14bae56415648ac100342f2b8bbb0"
SOURCE_PROFILE_PR = 458
TRANSFER_PR = 469
SOURCE_ALPHA = base.SOURCE_ALPHA
INNER_WINDOW = base.INNER_WINDOW
OUTER_WINDOW = base.OUTER_WINDOW
BASE_GAIN = 0.025
GAMMAS = (0.0, 0.010, 0.015, 0.020, 0.025)
TIMES = base.TIMES
RADII = base.RADII
GRID_SIZE = 33
SENSITIVITY_EPS = 1e-4

CRITERIA = dict(
    early_field_max_abs_delta=1e-12,
    late_additional_inner_gain_floor=0.02,
    late_additional_mid_gain_floor=0.01,
    late_additional_outer_gain_ceiling=0.0,
    axial_rms_abs_change_max=0.015,
    radial_rms_growth_max=0.02,
    collar_ratio_max=1.15,
    energy_min=0.1,
    energy_max=10.0,
    reference_energy_relative_error_max=1e-10,
    support_max_abs=1e-12,
    divergence_fd_max=1e-5,
    sensitivity_rank_min=2,
    sensitivity_condition_max=10.0,
    sensitivity_abs_cosine_max=0.95,
)


def time_shape(t: float) -> float:
    t = float(t)
    if t < TIMES[0] - 1e-14 or t > TIMES[-1] + 1e-14:
        raise ValueError("time outside registered interval")
    return 2.0 * (t - TIMES[0])


def gain_at_time(gamma: float, t: float, base_gain: float = BASE_GAIN) -> float:
    return float(base_gain) + float(gamma) * time_shape(float(t))


def kinetic_energy(f, raw, *, time: float, gain: float, scale: float, order: int = 72) -> float:
    x, wx = leggauss(order)
    y, wz = leggauss(order)
    r = x + 1.0
    z = 2.0 * y
    R, Z = np.meshgrid(r, z, indexing="ij")
    U = base.meridian_velocity(f, raw, R, Z, float(time)).copy()
    U[..., 1] *= 1.0 + float(gain) * base.h_profile(R)
    U *= float(scale)
    # 1/2 * 2*pi*r dr dz = pi*r dr dz, with z=2*y.
    w = wx[:, None] * (2.0 * wz[None, :]) * np.pi * R
    return float(np.sum(w * np.sum(U * U, axis=-1)))


def angular_summary(gamma: float, scale: float) -> dict:
    out = {}
    for t in TIMES:
        current_gain = gain_at_time(gamma, t)
        per_radius = {}
        for r in RADII:
            h = float(base.h_profile(np.array([r]))[0])
            static_multiplier = float(scale * (1.0 + BASE_GAIN * h))
            current_multiplier = float(scale * (1.0 + current_gain * h))
            per_radius[str(r)] = dict(
                absolute_gain_from_parent=current_multiplier - 1.0,
                additional_gain_from_static=current_multiplier / static_multiplier - 1.0,
            )
        out[str(t)] = per_radius
    return out


def sign_support_preflight_time(f, raw, *, time: float, gain: float, scale: float) -> tuple[bool, float]:
    core = np.array(
        [[r, 0.0, z] for r in (0.1, 0.3, 0.6) for z in (-0.5, -0.2, 0.2, 0.5)],
        dtype=float,
    )
    u = base.transformed_velocity(f, raw, core, float(time), float(gain), float(scale))
    signs = bool(
        np.all(u[:, 0] < 0.0)
        and np.all(u[:, 1] > 0.0)
        and np.all(np.sign(core[:, 2]) * u[:, 2] > 0.0)
    )
    outside = np.array(
        [[2.05, 0, 0], [-2.05, 0, 0], [0, 0, 2.05], [0, 0, -2.05], [1.8, 0, 2.05]],
        dtype=float,
    )
    support = float(
        np.max(np.abs(base.transformed_velocity(f, raw, outside, float(time), float(gain), float(scale))))
    )
    return signs, support


def divergence_fd_time(
    f,
    raw,
    *,
    time: float,
    gain: float,
    scale: float,
    seed: int,
    h: float = 1e-5,
) -> float:
    rng = np.random.default_rng(seed)
    pts = rng.uniform(-1.25, 1.25, size=(24, 3))
    div = np.zeros(len(pts), dtype=float)
    for axis in range(3):
        d = np.zeros(3, dtype=float)
        d[axis] = h
        up = base.transformed_velocity(f, raw, pts + d, float(time), float(gain), float(scale))
        um = base.transformed_velocity(f, raw, pts - d, float(time), float(gain), float(scale))
        div += (up[:, axis] - um[:, axis]) / (2.0 * h)
    return float(np.max(np.abs(div)))


def early_field_delta(f, raw, *, scale: float, gamma: float) -> float:
    points = np.array(
        [
            [0.31, 0.17, -0.51],
            [0.57, -0.24, -0.21],
            [0.82, 0.33, 0.18],
            [1.11, -0.27, 0.43],
            [-0.46, 0.62, -0.38],
            [-0.91, -0.34, 0.29],
        ],
        dtype=float,
    )
    reference = base.transformed_velocity(f, raw, points, TIMES[0], BASE_GAIN, scale)
    current = base.transformed_velocity(f, raw, points, TIMES[0], gain_at_time(gamma, TIMES[0]), scale)
    return float(np.max(np.abs(current - reference)))


def sensitivity_probes() -> list[tuple[float, np.ndarray]]:
    rows: list[tuple[float, np.ndarray]] = []
    for t in (0.30, 0.50, 0.70):
        points = []
        for r in (0.45, 0.85, 1.25):
            for z in (-0.55, 0.15):
                for phi in (0.37, 1.21):
                    points.append((r * np.cos(phi), r * np.sin(phi), z))
        rows.append((float(t), np.asarray(points, dtype=float)))
    return rows


def sensitivity_diagnostics(f, raw, parent_energy: float) -> dict:
    eps = SENSITIVITY_EPS
    columns = []

    # Static coordinate: shift k0 and recompute its reference-energy scale.
    vals_plus = []
    vals_minus = []
    for sign, sink in ((+1.0, vals_plus), (-1.0, vals_minus)):
        static_gain = BASE_GAIN + sign * eps
        scale, _ = base.child_scale(f, raw, parent_energy, gain=static_gain)
        chunks = []
        for t, points in sensitivity_probes():
            chunks.append(base.transformed_velocity(f, raw, points, t, static_gain, scale).ravel())
        sink.append(np.concatenate(chunks))
    columns.append((vals_plus[0] - vals_minus[0]) / (2.0 * eps))

    # Temporal coordinate: k0 and its reference-energy scale stay fixed.
    base_scale, _ = base.child_scale(f, raw, parent_energy, gain=BASE_GAIN)
    vals_plus = []
    vals_minus = []
    for sign, sink in ((+1.0, vals_plus), (-1.0, vals_minus)):
        chunks = []
        gamma = sign * eps
        for t, points in sensitivity_probes():
            chunks.append(
                base.transformed_velocity(
                    f,
                    raw,
                    points,
                    t,
                    gain_at_time(gamma, t),
                    base_scale,
                ).ravel()
            )
        sink.append(np.concatenate(chunks))
    columns.append((vals_plus[0] - vals_minus[0]) / (2.0 * eps))

    raw_matrix = np.column_stack(columns)
    norms = np.linalg.norm(raw_matrix, axis=0)
    if np.any(norms <= np.finfo(float).tiny):
        raise RuntimeError("degenerate sensitivity column")
    matrix = raw_matrix / norms[None, :]
    singular = np.linalg.svd(matrix, compute_uv=False)
    rank = int(np.sum(singular > 1e-8))
    condition = float(singular[0] / singular[-1])
    cosine = float(np.dot(matrix[:, 0], matrix[:, 1]))
    return dict(
        probe_count=sum(len(points) for _, points in sensitivity_probes()),
        velocity_component_rows=int(matrix.shape[0]),
        finite_difference_step=eps,
        column_order=["static_base_gain", "temporal_gamma"],
        raw_column_norms=[float(x) for x in norms],
        normalized_rank=rank,
        normalized_singular_values=[float(x) for x in singular],
        normalized_condition_number=condition,
        normalized_column_cosine=cosine,
    )


def clean_temporal_rule(row: dict, sensitivity: dict) -> bool:
    c = CRITERIA
    late = row["angular"][str(TIMES[-1])]
    return bool(
        row["gamma"] > 0.0
        and row["early_field_max_abs_delta"] <= c["early_field_max_abs_delta"]
        and late[str(0.6)]["additional_gain_from_static"] >= c["late_additional_inner_gain_floor"]
        and late[str(0.9)]["additional_gain_from_static"] >= c["late_additional_mid_gain_floor"]
        and late[str(1.2)]["additional_gain_from_static"] <= c["late_additional_outer_gain_ceiling"]
        and max(abs(x) for x in row["axial_rms_relative_changes_vs_static"]) <= c["axial_rms_abs_change_max"]
        and max(row["radial_rms_relative_changes_vs_static"]) <= c["radial_rms_growth_max"]
        and max(row["collar_ratios_vs_static"]) <= c["collar_ratio_max"]
        and min(row["energies"]) >= c["energy_min"]
        and max(row["energies"]) <= c["energy_max"]
        and row["reference_energy_relative_error"] <= c["reference_energy_relative_error_max"]
        and row["core_signs_pass"]
        and row["support_max_abs"] <= c["support_max_abs"]
        and row["divergence_fd_max"] <= c["divergence_fd_max"]
        and sensitivity["normalized_rank"] >= c["sensitivity_rank_min"]
        and sensitivity["normalized_condition_number"] <= c["sensitivity_condition_max"]
        and abs(sensitivity["normalized_column_cosine"]) <= c["sensitivity_abs_cosine_max"]
    )


def run(out: Path) -> dict:
    f, raw = replay_st051.reconstruct(PARENT_ID)
    balance = base.reference_energy_and_moment(f, raw)
    parent_energy = float(balance["parent_energy"])
    scale, static_raw_energy = base.child_scale(f, raw, parent_energy, gain=BASE_GAIN)

    static_quad = {
        t: base.quadrature_vorticity_metrics(f, raw, BASE_GAIN, scale, t)
        for t in TIMES
    }
    static_grid = {
        t: base.grid_morphology(f, raw, time=t, gain=BASE_GAIN, scale=scale, grid_size=GRID_SIZE)
        for t in TIMES
    }
    static_energy = {
        t: kinetic_energy(f, raw, time=t, gain=BASE_GAIN, scale=scale)
        for t in TIMES
    }
    sensitivity = sensitivity_diagnostics(f, raw, parent_energy)

    rows = []
    for gamma in GAMMAS:
        angular = angular_summary(gamma, scale)
        current_quad = {}
        current_grid = {}
        energies = []
        signs = []
        supports = []
        divergences = []
        for index, t in enumerate(TIMES):
            gain = gain_at_time(gamma, t)
            current_quad[t] = base.quadrature_vorticity_metrics(f, raw, gain, scale, t)
            current_grid[t] = base.grid_morphology(
                f, raw, time=t, gain=gain, scale=scale, grid_size=GRID_SIZE
            )
            energies.append(kinetic_energy(f, raw, time=t, gain=gain, scale=scale))
            sign_ok, support = sign_support_preflight_time(
                f, raw, time=t, gain=gain, scale=scale
            )
            signs.append(sign_ok)
            supports.append(support)
            divergences.append(
                divergence_fd_time(
                    f,
                    raw,
                    time=t,
                    gain=gain,
                    scale=scale,
                    seed=9175741 + index,
                )
            )

        axial_changes = [
            current_quad[t]["axial_rms"] / static_quad[t]["axial_rms"] - 1.0
            for t in TIMES
        ]
        radial_changes = [
            current_quad[t]["radial_rms"] / static_quad[t]["radial_rms"] - 1.0
            for t in TIMES
        ]
        collar_ratios = [
            current_quad[t]["collar_absolute"] / static_quad[t]["collar_absolute"]
            for t in TIMES
        ]
        fingerprint = []
        for t in TIMES:
            fingerprint.append(
                dict(
                    time=t,
                    gain=gain_at_time(gamma, t),
                    q90=current_grid[t]["axial_q90_over_support"],
                    q99=current_grid[t]["axial_q99_over_support"],
                    q90_delta_from_static=current_grid[t]["axial_q90_over_support"] - static_grid[t]["axial_q90_over_support"],
                    q99_delta_from_static=current_grid[t]["axial_q99_over_support"] - static_grid[t]["axial_q99_over_support"],
                    outer_065=current_grid[t]["outer_065_enstrophy_fraction"],
                    outer_075=current_grid[t]["outer_075_enstrophy_fraction"],
                    collar=current_grid[t]["collar_fraction"],
                    vorticity_rms=current_grid[t]["vorticity_rms"],
                    vorticity_max=current_grid[t]["vorticity_max"],
                )
            )

        row = dict(
            gamma=float(gamma),
            time_gains=[gain_at_time(gamma, t) for t in TIMES],
            normalization=float(scale),
            angular=angular,
            early_field_max_abs_delta=early_field_delta(f, raw, scale=scale, gamma=gamma),
            axial_rms_relative_changes_vs_static=[float(x) for x in axial_changes],
            radial_rms_relative_changes_vs_static=[float(x) for x in radial_changes],
            collar_ratios_vs_static=[float(x) for x in collar_ratios],
            energies=[float(x) for x in energies],
            energy_ratios_vs_static=[float(e / static_energy[t]) for e, t in zip(energies, TIMES)],
            reference_energy_relative_error=float(abs(energies[0] / static_energy[TIMES[0]] - 1.0)),
            grid33=fingerprint,
            core_signs_pass=bool(all(signs)),
            support_max_abs=float(max(supports)),
            divergence_fd_max=float(max(divergences)),
        )
        row["clean_temporal_capacity"] = clean_temporal_rule(row, sensitivity)
        rows.append(row)

    crossings = [row for row in rows if row["clean_temporal_capacity"]]
    first = crossings[0]["gamma"] if crossings else None
    result = dict(
        task_id=TASK_ID,
        prereg_issue=PREREG_ISSUE,
        parent_id=PARENT_ID,
        parent_head=PARENT_HEAD,
        stack_base_pr=STACK_BASE_PR,
        stack_base_head=STACK_BASE_HEAD,
        source_profile_pr=SOURCE_PROFILE_PR,
        transfer_pr=TRANSFER_PR,
        frozen_profile=dict(
            alpha=SOURCE_ALPHA,
            inner_window=list(INNER_WINDOW),
            outer_window=list(OUTER_WINDOW),
            base_gain=BASE_GAIN,
        ),
        temporal_shape="k(t)=0.025+gamma*2*(t-0.25)",
        preregistered_gammas=list(GAMMAS),
        gammas_outside_preregistered_grid_evaluated=False,
        criteria=CRITERIA,
        balance_on_st051b=balance,
        base_child_normalization=float(scale),
        base_child_raw_reference_energy=float(static_raw_energy),
        static_child_energies={str(t): float(static_energy[t]) for t in TIMES},
        sensitivity=sensitivity,
        rows=rows,
        smallest_preregistered_clean_temporal_capacity=first,
        recommendation=(
            "One temporal coefficient on the already-frozen redistribution channel supplies independent late-time winding control without another spatial basis. Keep generic swirl/ring/cap basis growth frozen; treat the crossing only as capacity and let governed 3-D render/material-path evidence decide whether this temporal degree is needed in a production child."
            if first is not None
            else
            "No preregistered temporal coefficient passes the frozen capacity rule. Do not widen gamma post hoc; wait for governed 3-D morphology evidence before choosing a different minimal basis channel."
        ),
        diagnostic_added_spatial_basis_degrees=0,
        diagnostic_temporal_degrees_screened=1,
        canonical_velocity_changed=False,
        candidate_artifact_changed=False,
        pressure_or_force_changed=False,
        held_out_pde_residual_evaluated=False,
        material_paths_integrated=False,
        public_image_used=False,
        production_temporal_coefficient_selected=False,
        visualization_ready=False,
        visual_correspondence_verified=False,
        pde_validated=False,
        source_correspondence_verified=False,
        paper_exact=False,
        openai_field_identified=False,
        scope="One preregistered temporal coefficient on an existing frozen redistribution coordinate; no new spatial basis or acceptance claim.",
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
