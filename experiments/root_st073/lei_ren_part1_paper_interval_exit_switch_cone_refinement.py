"""Directed refinement of the unresolved final stage-1 switch cell.

The frozen switch receipt stores cone summaries but not prefix states.  This
companion reconstructs the R = 100 prefix with the same directed cell update,
then subdivides only x in [.0046875, .005].  Stage-1 cone arithmetic keeps the
algebraic relations

    S_theta / F = -a,
    S_z / F = -I_z * epsilon_exit * blend / (F0 * Phi_bar),

before forming kappa.  This removes the avoidable interval dependency from
dividing by the tiny angular shear.  It is a local refinement diagnostic and
does not replace the frozen switch receipt or claim whole-switch closure.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import mpmath as mp

from lei_ren_part1_paper_candidate_pressure_function import ACCEPTED_SHA
from lei_ren_part1_paper_interval_comparison_jets import IntervalComparisonJets
from lei_ren_part1_paper_interval_comparison_enclosure import alpha_box
from lei_ren_part1_paper_interval_exit_continuation_enclosure import (
    _driver,
    _frozen_state,
    load_receipts,
    continue_exit,
)
from lei_ren_part1_paper_interval_exit_stress_enclosure import evaluate
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode


HERE = Path(__file__).resolve().parent
SWITCH_RECEIPT = HERE / "lei_ren_part1_paper_interval_exit_switch_enclosure.json"


def _pack_all(value):
    if hasattr(value, "coefficients") and hasattr(value, "order"):
        return [_pack_all(item) for item in value.coefficients]
    if isinstance(value, dict):
        return {key: _pack_all(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_pack_all(item) for item in value]
    return value


def _hash(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _validate_frozen_switch(calc):
    raw = json.loads(SWITCH_RECEIPT.read_text(encoding="utf-8"))
    if raw.get("state_sha256") != calc.state_hash:
        raise ValueError("Frozen switch receipt does not match the loaded core")
    if raw.get("cells_per_half") != 16:
        raise ValueError("Frozen switch receipt is not the 16-cell production receipt")
    if not raw.get("controlled_switch_cell_integrals", False):
        raise ValueError("Frozen switch cell enclosure is not marked controlled")
    rows = raw.get("switch_cell_cones", [])
    unresolved = [row for row in rows if row.get("index") == 15]
    if len(unresolved) != 1:
        raise ValueError("Frozen switch receipt has no unique final stage-1 cell")
    if not str(unresolved[0].get("cone", {}).get("status", "")).startswith(
        "not certified"
    ):
        raise ValueError("Frozen switch final stage-1 cell is no longer unresolved")
    for name, digest in raw.get("input_hashes", {}).items():
        dep = HERE / name
        if not dep.exists() or _hash(dep) != digest:
            raise ValueError("Frozen switch dependency hash changed: " + name)


def _nonnegative(c, value):
    lo, hi = endpoints(value)
    if hi <= 0:
        return c.mpf(0)
    return c.mpf([max(mp.mpf(0), lo), hi])


def _stage1_cone(calc, stress, a, blend, driver, eps_exit, bar_phi):
    """Cone box with stage-1 algebraic cancellation retained explicitly."""

    c = calc.ctx
    a0 = a[0] if hasattr(a, "coefficients") else a
    if endpoints(a0)[0] <= 0:
        return dict(
            prerequisites_certified=False,
            relaxed_cone_certified=False,
            admissible_cone_certified=False,
            status="not certified: stage-1 a interval is not positive",
        )
    if endpoints(calc.F0)[0] <= 0 or endpoints(bar_phi[0])[0] <= 0:
        return dict(
            prerequisites_certified=False,
            relaxed_cone_certified=False,
            admissible_cone_certified=False,
            status="not certified: amplitude denominator is not positive",
        )
    # B = -(phi0/bar_phi) I_z eps_exit sqrt(R/2) blend and
    # phi = phi0 exp(g), so the exp(g) and phi0 factors cancel exactly in
    # S_z/F.  I_z here is the already sqrt(epsilon)-scaled driver returned by
    # the validated continuation helper.
    sz = -(
        driver["I_z"][0] * eps_exit * blend
        / (calc.F0 * bar_phi[0])
    )
    st = -a0
    tt = stress["normalized_stress"]["T_theta_over_F"]
    tz = stress["normalized_stress"]["T_z_over_F"]
    kappa = a0 + (sz * sz) / a0
    dot = -tt * a0 + tz * sz
    cross = -tt * sz - tz * a0
    direction = endpoints(dot)[1] < 0
    if endpoints(kappa)[0] > 2:
        branch = "kappa>2"
        margin = 2 * (dot * dot) - (kappa - 2) * (cross * cross)
    elif endpoints(kappa)[1] <= 2:
        branch = "kappa<=2"
        # -dot/(-st) = tt - tz*sz/a, with positive a.  This is the
        # correlated form of the relaxed-cone margin.
        margin = tt - tz * sz / a0 - (2 - kappa)
    else:
        return dict(
            prerequisites_certified=True,
            kappa=kappa,
            dot=dot,
            cross=cross,
            direction_certified=direction,
            branch="not isolated",
            stage1_algebraic_cancellation=True,
            sz_over_F=sz,
            relaxed_cone_certified=False,
            admissible_cone_certified=False,
            status="not certified: correlated kappa interval crosses 2",
        )
    passed = direction and endpoints(margin)[0] > 0
    ruled_out = endpoints(dot)[0] >= 0 or endpoints(margin)[1] <= 0
    return dict(
        prerequisites_certified=True,
        kappa=kappa,
        dot=dot,
        cross=cross,
        direction_certified=direction,
        branch=branch,
        margin=margin,
        stage1_algebraic_cancellation=True,
        sz_over_F=sz,
        relaxed_cone_certified=passed,
        admissible_cone_certified=passed and branch == "kappa>2",
        family_cone_ruled_out=ruled_out,
        status="certified" if passed else (
            "ruled out over enclosed family" if ruled_out
            else "not certified: correlated margin/direction unresolved"
        ),
    )


def _advance_stage1(calc, continuation, state, x, width, *, inspect=True):
    """Advance one stage-1 cell and optionally return its cone packet."""

    c = calc.ctx
    phi0, u0 = state["phi0"], state["u0"]
    bar = continuation["endpoint_comparison_state"]
    se = continuation["endpoint_scaled_radius"]
    eps_exit = continuation["exit_shear_epsilon"]
    s = c.exp(x) * 100 / calc.eps
    ds = _nonnegative(c, s - se)
    ds2 = _nonnegative(c, s * s - se * se)
    comparison = _frozen_state(bar, bar["phi"], bar["U"], ds, ds2)
    driver = _driver(calc, bar["phi"], bar["U"], comparison, s)
    blend = alpha_box(c, x + c.mpf(".005"), c.mpf(".005"))
    a = driver["D"] * eps_exit
    A = -a / 2
    B = -(
        (phi0 / bar["phi"].truncate(1))
        * driver["I_z"]
        * eps_exit
        * c.sqrt(s * calc.eps / 2)
        * blend
    )
    partial = c.mpf([0, endpoints(width)[1]])
    gr = state["g"] + A * partial
    phi = phi0 * gr.exp()
    urhs = gr.exp() * B
    u = u0 + state["j"] + urhs * partial
    rhs = dict(
        theta=phi * (2 * s * s),
        z=u * s,
        theta_z=phi * u * (2 * s * s),
        p=phi * phi * s,
        u_squared=u * u * s,
        weighted_phi_squared=phi * phi * s * s,
    )
    range_m = {key: state["moments"][key] + rhs[key] * partial for key in rhs}
    next_state = dict(state)
    next_state["moments"] = {
        key: state["moments"][key] + rhs[key] * width for key in rhs
    }
    next_state["g"] = state["g"] + A * width
    next_state["j"] = state["j"] + urhs * width
    if not inspect:
        return next_state, None
    packet = dict(
        normalized_exit_state=dict(phi=phi, U=u, **range_m),
        physical_R=s * calc.eps,
        actual_terminal_g_y=A,
        actual_terminal_U_y=urhs,
    )
    stress = evaluate(calc, packet)
    cone = _stage1_cone(calc, stress, a, blend, driver, eps_exit, bar["phi"])
    return next_state, dict(
        x=x,
        width=width,
        scaled_radius=s,
        blend=blend,
        driver_D=driver["D"],
        driver_I_z=driver["I_z"],
        a=a,
        generic_cone=stress["cone"],
        cone=cone,
        normalized_stress=stress["normalized_stress"],
    )


def refine(calc, continuation, subdivisions=128):
    if not isinstance(subdivisions, int) or subdivisions < 1:
        raise ValueError("subdivisions must be a positive integer")
    c = calc.ctx
    with mp.workdps(calc.precision + 40):
        state0 = continuation["normalized_exit_state"]
        state = dict(
            phi0=state0["phi"],
            u0=state0["U"],
            g=state0["phi"] * 0,
            j=state0["phi"] * 0,
            moments={
                key: state0[key]
                for key in state0
                if key not in ("phi", "U")
            },
        )
        hb = c.mpf(".005")
        parent_step = hb / 16
        # Reconstruct the exact directed prefix through the first 15 stored
        # stage-1 cells; no core or source data is regenerated here.
        for index in range(15):
            left = parent_step * index
            right = parent_step * (index + 1)
            x = c.mpf([endpoints(left)[0], endpoints(right)[1]])
            state, _ = _advance_stage1(
                calc, continuation, state, x, parent_step, inspect=False
            )
        left = parent_step * 15
        right = hb
        refined_step = parent_step / subdivisions
        rows = []
        for index in range(subdivisions):
            x_left = left + refined_step * index
            x_right = left + refined_step * (index + 1)
            x = c.mpf([endpoints(x_left)[0], endpoints(x_right)[1]])
            state, row = _advance_stage1(
                calc, continuation, state, x, refined_step, inspect=True
            )
            row["index"] = index
            row["log_interval"] = [mp.nstr(endpoints(x)[0], 100), mp.nstr(endpoints(x)[1], 100)]
            rows.append(row)
        counts = {}
        for row in rows:
            status = row["cone"]["status"]
            counts[status] = counts.get(status, 0) + 1
        return dict(
            center_family=list(calc.center_family),
            state_sha256=calc.state_hash,
            accepted_schedule_sha256=ACCEPTED_SHA,
            source_cell_log_interval=[".0046875", ".005"],
            subdivisions=subdivisions,
            prefix_parent_cells=15,
            rows=rows,
            counts=counts,
            generic_last_parent_cell_status="not certified: kappa interval crosses 2",
            correlated_cone_used=True,
            stage1_algebraic_cancellation=True,
            no_point_sample_certification=True,
            no_core_regeneration=True,
            whole_switch_cone_certified=False,
            terminal_matching_complete=False,
            original_construction_parameter_errors_enclosed=False,
            temporal_recursion=False,
        )


def run(subdivisions=128):
    calc = IntervalComparisonJets(4)
    _validate_frozen_switch(calc)
    bridge, cells = load_receipts(calc)
    continuation = continue_exit(
        calc, bridge, cells["endpoint"], cells["initial_phi"], target_R="100"
    )
    with mp.workdps(calc.precision + 40):
        result = refine(calc, continuation, subdivisions=subdivisions)
    names = (
        "lei_ren_part1_paper_interval_exit_switch_enclosure.py",
        "lei_ren_part1_paper_interval_exit_switch_enclosure.json",
        "lei_ren_part1_paper_interval_exit_continuation_enclosure.py",
        "lei_ren_part1_paper_interval_exit_continuation_enclosure.json",
        "lei_ren_part1_paper_interval_exit_stress_enclosure.py",
        "lei_ren_part1_paper_interval_comparison_enclosure.py",
    )
    result["input_hashes"] = {name: _hash(HERE / name) for name in names}
    result["source_sha256"] = _hash(Path(__file__))
    output = HERE / "lei_ren_part1_paper_interval_exit_switch_cone_refinement.json"
    output.write_text(json.dumps(encode(_pack_all(result)), indent=2) + "\n", encoding="utf-8")
    print("Stage-1 final-cell refinement:", result["counts"], flush=True)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--subdivisions", type=int, default=128)
    args = parser.parse_args()
    run(args.subdivisions)


if __name__ == "__main__":
    main()
