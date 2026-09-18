"""Screen bounded gain headroom of the frozen ST050R-C redistribution on ST051-B.

Preregistered in issue #479 before evaluation. This Agent-7 increment keeps the
radial profile identity from #458/#469 fixed and varies only its already-existing
one-dimensional gain on frozen ST051-B. It is a target-free expression-capacity
screen, not a production coefficient selection, PDE validation, image fit, or
material-path experiment.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import agent7_st051b_frozen_redistribution_transfer as base
import replay_st051

TASK_ID = "CR003-ST051B-FROZEN-REDISTRIBUTION-HEADROOM-073"
PREREG_ISSUE = 479
PARENT_ID = "ST051-B"
PARENT_HEAD = base.PARENT_HEAD
SOURCE_PROFILE_PR = 458
TRANSFER_PR = 469
SOURCE_ALPHA = base.SOURCE_ALPHA
INNER_WINDOW = base.INNER_WINDOW
OUTER_WINDOW = base.OUTER_WINDOW
GAINS = (0.0, 0.015, 0.025, 0.035, 0.05)
TIMES = base.TIMES
RADII = base.RADII
GRID_SIZES = (25, 33)

CRITERIA = dict(
    inner_gain_floor=0.04,
    mid_gain_floor=0.025,
    outer_gain_ceiling=0.0,
    normalization_deviation_max=0.005,
    inward_speed_loss_max=0.005,
    axial_rms_abs_change_max=0.02,
    radial_rms_growth_max=0.03,
    collar_ratio_max=1.25,
    support_max_abs=1e-12,
    divergence_fd_max=1e-5,
)


def clean_headroom_rule(row):
    c = CRITERIA
    return bool(
        row["gain"] > 0.025
        and row["angular_gain_r06"] >= c["inner_gain_floor"]
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


def angular_gains(gain, scale):
    out = {}
    for r in RADII:
        h = float(base.h_profile(np.array([r]))[0])
        out[r] = float(scale * (1.0 + float(gain) * h) - 1.0)
    return out


def grid_refinement_summary(grid_by_size):
    coarse = grid_by_size[str(GRID_SIZES[0])]
    fine = grid_by_size[str(GRID_SIZES[1])]
    rows = []
    for c, f in zip(coarse, fine):
        assert c["time"] == f["time"]
        keys = (
            "axial_rms",
            "radial_rms",
            "axial_q90_over_support",
            "axial_q99_over_support",
            "outer_065_enstrophy_fraction",
            "outer_075_enstrophy_fraction",
            "collar_fraction",
            "vorticity_rms",
            "vorticity_max",
        )
        rows.append(
            dict(
                time=c["time"],
                fine_minus_coarse={k: float(f["metrics"][k] - c["metrics"][k]) for k in keys},
                relative_change={
                    k: float(f["metrics"][k] / max(abs(c["metrics"][k]), 1e-300) - 1.0)
                    for k in keys
                },
            )
        )
    return rows


def run(out: Path):
    f, raw = replay_st051.reconstruct(PARENT_ID)
    balance = base.reference_energy_and_moment(f, raw)
    parent_energy = balance["parent_energy"]

    quad_parent = {t: base.quadrature_vorticity_metrics(f, raw, 0.0, 1.0, t) for t in TIMES}
    grid_parent = {
        str(n): [
            dict(time=t, metrics=base.grid_morphology(f, raw, time=t, gain=0.0, scale=1.0, grid_size=n))
            for t in TIMES
        ]
        for n in GRID_SIZES
    }

    rows = []
    for gain in GAINS:
        if gain == 0.0:
            scale = 1.0
            raw_energy = parent_energy
        else:
            scale, raw_energy = base.child_scale(f, raw, parent_energy, gain=gain)

        angular = angular_gains(gain, scale)
        quad = {t: base.quadrature_vorticity_metrics(f, raw, gain, scale, t) for t in TIMES}
        axial = [quad[t]["axial_rms"] / quad_parent[t]["axial_rms"] - 1.0 for t in TIMES]
        radial = [quad[t]["radial_rms"] / quad_parent[t]["radial_rms"] - 1.0 for t in TIMES]
        collar = [quad[t]["collar_absolute"] / quad_parent[t]["collar_absolute"] for t in TIMES]

        grids = {}
        for n in GRID_SIZES:
            vals = []
            for t in TIMES:
                if gain == 0.0:
                    metrics = next(x["metrics"] for x in grid_parent[str(n)] if x["time"] == t)
                else:
                    metrics = base.grid_morphology(f, raw, time=t, gain=gain, scale=scale, grid_size=n)
                vals.append(dict(time=t, metrics=metrics))
            grids[str(n)] = vals

        signs, support = base.sign_support_preflight(f, raw, gain, scale)
        div = base.divergence_fd(f, raw, gain, scale, seed=9175731)
        row = dict(
            gain=float(gain),
            normalization=float(scale),
            normalization_deviation=float(scale - 1.0),
            raw_child_energy=float(raw_energy),
            inward_speed_loss=float(max(0.0, 1.0 - scale)),
            angular_gain_r06=angular[0.6],
            angular_gain_r09=angular[0.9],
            angular_gain_r12=angular[1.2],
            axial_rms_relative_changes=[float(x) for x in axial],
            radial_rms_relative_changes=[float(x) for x in radial],
            collar_ratios=[float(x) for x in collar],
            q33=[
                dict(
                    time=x["time"],
                    q90=x["metrics"]["axial_q90_over_support"],
                    q99=x["metrics"]["axial_q99_over_support"],
                    outer_065=x["metrics"]["outer_065_enstrophy_fraction"],
                    outer_075=x["metrics"]["outer_075_enstrophy_fraction"],
                    collar=x["metrics"]["collar_fraction"],
                    vorticity_rms=x["metrics"]["vorticity_rms"],
                    vorticity_max=x["metrics"]["vorticity_max"],
                )
                for x in grids["33"]
            ],
            grid_refinement_25_to_33=grid_refinement_summary(grids),
            core_signs_pass=bool(signs),
            support_max_abs=float(support),
            divergence_fd_max=float(div),
        )
        row["clean_higher_gain_capacity"] = clean_headroom_rule(row)
        rows.append(row)

    crossings = [r for r in rows if r["clean_higher_gain_capacity"]]
    first = crossings[0]["gain"] if crossings else None

    parent_q33 = grid_parent["33"]
    result = dict(
        task_id=TASK_ID,
        prereg_issue=PREREG_ISSUE,
        parent_id=PARENT_ID,
        parent_head=PARENT_HEAD,
        source_profile_pr=SOURCE_PROFILE_PR,
        transfer_pr=TRANSFER_PR,
        frozen_profile=dict(alpha=SOURCE_ALPHA, inner_window=list(INNER_WINDOW), outer_window=list(OUTER_WINDOW)),
        preregistered_gains=list(GAINS),
        gains_outside_preregistered_grid_evaluated=False,
        criteria=CRITERIA,
        balance_on_st051b=balance,
        parent_grid33=parent_q33,
        rows=rows,
        smallest_preregistered_clean_higher_gain_capacity=first,
        recommendation=(
            "The already-frozen one-dimensional redistribution retains bounded higher-gain capacity on ST051-B. Keep generic swirl/ring basis dimension frozen; wait for the fixed-seed 3-D render to identify a different missing morphology channel before adding basis dimensions. Do not interpret the crossing as a production gain."
            if first is not None else
            "No preregistered higher-gain row passes the frozen capacity rule. Do not widen the gain grid or add generic swirl basis in this increment; use the failed channel to route the next minimal representation experiment after the fixed-seed 3-D render."
        ),
        diagnostic_added_velocity_degrees=0,
        existing_velocity_degree_screened=True,
        canonical_velocity_changed=False,
        candidate_artifact_changed=False,
        pressure_or_force_changed=False,
        held_out_pde_residual_evaluated=False,
        material_paths_integrated=False,
        public_image_used=False,
        production_redistribution_gain_selected=False,
        visualization_ready=False,
        visual_correspondence_verified=False,
        pde_validated=False,
        source_correspondence_verified=False,
        paper_exact=False,
        openai_field_identified=False,
        scope="One preregistered gain-headroom screen of an existing frozen redistribution coordinate; no new basis dimension or acceptance claim.",
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
