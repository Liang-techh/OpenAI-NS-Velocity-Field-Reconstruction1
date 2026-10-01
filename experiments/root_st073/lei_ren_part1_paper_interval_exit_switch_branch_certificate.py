"""Piecewise branch certificate for the final stage-1 switch cell.

The frozen 128-way refinement leaves only
``x in [.00499755859375, .005]`` unresolved because the interval for
``kappa = a + sz**2/a`` crosses 2.  This companion makes a directed
three-piece split at two explicit endpoint distances.  It keeps the
correlated relations

    st = -a,
    sz = -I_z * epsilon_exit * blend / (F0 * Phi_bar)

and records the two certified outer branches separately.  The middle
transition strip is deliberately retained as ``not certified``; this is a
branch enclosure, not a point sample or a replacement for the whole switch
receipt.
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
from lei_ren_part1_paper_candidate_shared_inlet import normalized_inlet
from lei_ren_part1_paper_interval_exit_continuation_enclosure import (
    _driver,
    _frozen_state,
    continue_exit,
    load_receipts,
)
from lei_ren_part1_paper_interval_exit_switch_cone_refinement import (
    _advance_stage1,
    _stage1_cone,
    _validate_frozen_switch,
)
from lei_ren_part1_paper_interval_exit_stress_enclosure import evaluate
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode


HERE = Path(__file__).resolve().parent
SWITCH_RECEIPT = HERE / "lei_ren_part1_paper_interval_exit_switch_enclosure.json"
REFINEMENT_RECEIPT = (
    HERE / "lei_ren_part1_paper_interval_exit_switch_cone_refinement.json"
)


def _pack(value):
    """Convert IntervalTaylor-like objects to recursively serializable data."""

    if hasattr(value, "coefficients") and hasattr(value, "order"):
        return [_pack(item) for item in value.coefficients]
    if isinstance(value, dict):
        return {key: _pack(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_pack(item) for item in value]
    return value


def _hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _point_box(c, value):
    """Preserve a decimal cut as a directed point interval."""

    point = c.mpf(value)
    lo, hi = endpoints(point)
    return c.mpf([lo, hi])


def _prefix_state(calc, continuation):
    """Reconstruct the prefix before the unresolved 128-way subcell."""

    c = calc.ctx
    hb = c.mpf(".005")
    parent_step = hb / 16
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
    for index in range(15):
        left = parent_step * index
        right = parent_step * (index + 1)
        x = c.mpf([endpoints(left)[0], endpoints(right)[1]])
        state, _ = _advance_stage1(
            calc, continuation, state, x, parent_step, inspect=False
        )
    # The immutable refinement receipt resolves the first 127 subcells of
    # the final parent cell.  Rebuild that directed prefix so the companion
    # starts exactly at x = .00499755859375, rather than re-certifying the
    # already settled parent interval.
    refined_step = parent_step / 128
    left = parent_step * 15
    for index in range(127):
        x_left = left + refined_step * index
        x_right = left + refined_step * (index + 1)
        x = c.mpf([endpoints(x_left)[0], endpoints(x_right)[1]])
        state, _ = _advance_stage1(
            calc, continuation, state, x, refined_step, inspect=False
        )
    return state, left + refined_step * 127


def _small_branch_certificate(calc, row):
    """Certify the relaxed branch from signs, without dividing by tiny shear.

    If ``0 <= kappa <= 2``, ``a > 0``, ``tt > 2``, ``tz <= 0`` and ``sz >= 0``,
    then

        tt - tz*sz/a - (2-kappa) >= tt - 2 > 0,

    while ``dot = -tt*a + tz*sz < 0``.  This uses the exact stage-1
    cancellation and remains valid even when the generic angular interval is
    very wide.
    """

    c = calc.ctx
    stress = row["normalized_stress"]
    a = row["a"]
    a0 = a[0] if hasattr(a, "coefficients") else a
    tt = stress["T_theta_over_F"]
    tz = stress["T_z_over_F"]
    sz = stress["S_z_over_F"]
    if endpoints(a0)[0] <= 0:
        return dict(
            conditional_branch="kappa<=2",
            conditional_prerequisites_certified=False,
            conditional_relaxed_cone_certified=False,
            status="not certified: a is not strictly positive",
        )
    if endpoints(tt)[0] <= 2 or endpoints(tz)[1] > 0 or endpoints(sz)[0] < 0:
        return dict(
            conditional_branch="kappa<=2",
            conditional_prerequisites_certified=False,
            conditional_relaxed_cone_certified=False,
            status="not certified: sign or tt>2 prerequisite unresolved",
        )
    # This is a lower bound valid under the conditional kappa<=2 branch.
    margin_lower = tt - tz * sz / a0 - c.mpf(2)
    direction_upper = -tt * a0 + tz * sz
    passed = endpoints(margin_lower)[0] > 0 and endpoints(direction_upper)[1] < 0
    return dict(
        conditional_branch="kappa<=2",
        conditional_prerequisites_certified=True,
        conditional_relaxed_margin_lower=margin_lower,
        conditional_direction_upper=direction_upper,
        conditional_relaxed_cone_certified=passed,
        conditional_admissible_cone_certified=False,
        status="certified" if passed else "not certified: conditional margin unresolved",
    )


def _compact_row(calc, label, x, width, row):
    """Keep exact branch inputs and outputs without serializing prefix state."""

    cone = row["cone"]
    packet = dict(
        label=label,
        x=x,
        width=width,
        scaled_radius=row["scaled_radius"],
        blend=row["blend"],
        a=row["a"],
        driver_D=row["driver_D"],
        driver_I_z=row["driver_I_z"],
        normalized_stress=row["normalized_stress"],
        cone=cone,
    )
    if label == "small_blend_branch":
        packet["conditional_small_branch"] = _small_branch_certificate(calc, row)
    return packet


def _advance_full_cell_with_raw_tz(calc, continuation, state, x, width):
    """Advance one cell and expose the F0-cancelled physical Tz factors.

    The ordinary cone packet forms ``S_z/F`` and ``T_z/F`` separately, which
    leaves a huge F0-dependent interval when the blend reaches zero.  Here we
    retain the physical inlet value

        Tz = Iactual + Braw,
        Braw = epsilon_exit * (-Ibar) * phi / Phi_bar * blend,

    before dividing by F0.  Therefore ``b/v = Braw/(-Tz)`` and the strong
    branch sufficient condition can be tested without F0.
    """

    c = calc.ctx
    phi0, u0 = state["phi0"], state["u0"]
    bar = continuation["endpoint_comparison_state"]
    se = continuation["endpoint_scaled_radius"]
    eps_exit = continuation["exit_shear_epsilon"]
    s = c.exp(x) * 100 / calc.eps
    ds = c.mpf(
        [max(mp.mpf(0), endpoints(s - se)[0]), endpoints(s - se)[1]]
    )
    ds2 = c.mpf(
        [max(mp.mpf(0), endpoints(s * s - se * se)[0]), endpoints(s * s - se * se)[1]]
    )
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
    packet = dict(
        normalized_exit_state=dict(phi=phi, U=u, **range_m),
        physical_R=s * calc.eps,
        actual_terminal_g_y=A,
        actual_terminal_U_y=urhs,
    )
    stress = evaluate(calc, packet)
    cone = _stage1_cone(calc, stress, a, blend, driver, eps_exit, bar["phi"])

    # Recompute the physical inlet T_z from the same packet.  The boundary
    # slope contribution is Braw; subtracting it is an algebraic split, not
    # a subtraction of two nearly equal normalized stresses.
    inlet = normalized_inlet(
        [phi],
        [u],
        range_m,
        calc.S.truncate(1),
        calc.ell.truncate(1),
        calc.p0.truncate(1),
        calc.lam,
        s,
        calc.z.truncate(1),
        calc.delta,
        endpoints_override=dict(phi_exit=phi, u_exit=u, phi_s=phi * A / s, u_s=urhs / s),
    )
    Tz = inlet["tz"][0] * c.sqrt(calc.eps)
    Ibar = driver["I_z"][0]
    Braw = -Ibar * eps_exit * phi[0] / bar["phi"][0] * blend
    Braw_from_slope = (
        (urhs[0] / s) * c.sqrt(2 * s) * calc.lam * c.sqrt(calc.eps)
    )
    Iactual = Tz - Braw
    t = inlet["ttheta"][0] / phi[0]
    D = driver["D"][0]
    # b/v = Braw/(-Tz), and a = eps_exit*D.  This is the F0-cancelled beta.
    beta = eps_exit * (
        ((-Ibar) * phi[0] / bar["phi"][0] * blend / (-Iactual - Braw)) ** 2
        * t
        * t
        / D
    )
    return next_state, dict(
        x=x,
        width=width,
        scaled_radius=s,
        blend=blend,
        driver_D=driver["D"],
        driver_I_z=driver["I_z"],
        a=a,
        normalized_stress=stress["normalized_stress"],
        cone=cone,
        raw_physical_Tz=Tz,
        raw_Iactual=Iactual,
        raw_Braw=Braw,
        raw_Braw_from_slope=Braw_from_slope,
        raw_Braw_formula_overlap=(
            endpoints(Braw)[0] <= endpoints(Braw_from_slope)[1]
            and endpoints(Braw_from_slope)[0] <= endpoints(Braw)[1]
        ),
        raw_t=t,
        raw_beta=beta,
    )


def _both_branch_certificate(calc, raw_row):
    """Certify both conditional cone branches on one common full-cell box."""

    c = calc.ctx
    t = raw_row["raw_t"]
    Tz = raw_row["raw_physical_Tz"]
    Braw = raw_row["raw_Braw"]
    beta = raw_row["raw_beta"]
    a = raw_row["a"]
    a0 = a[0] if hasattr(a, "coefficients") else a
    # b >= 0 follows from Braw >= 0 and F0*phi > 0.  The sign of Tz gives
    # v=-Tz/(F0*phi)>0, while t is already T_theta/F.
    prerequisites = dict(
        a_positive=endpoints(a0)[0] > 0,
        a_less_than_two=endpoints(a0)[1] < 2,
        t_greater_than_two=endpoints(t)[0] > 2,
        physical_Tz_negative=endpoints(Tz)[1] < 0,
        Braw_nonnegative=endpoints(Braw)[0] >= 0,
    )
    strong_bound = c.mpf(2) - a0 - beta
    relaxed_bound = t - c.mpf(2)
    passed = all(prerequisites.values()) and endpoints(strong_bound)[0] > 0
    return dict(
        branch="both conditional branches",
        prerequisites=prerequisites,
        strong_beta=beta,
        strong_sufficient_margin=strong_bound,
        relaxed_sufficient_margin=relaxed_bound,
        both_conditional_branches_certified=passed,
        strong_cone_certified=False,
        relaxed_cone_certified=passed,
        status=(
            "certified: both kappa branches covered by common beta bound"
            if passed
            else "not certified: F0-cancelled beta prerequisites unresolved"
        ),
        expansion_identity=(
            "2dot^2-(kappa-2)cross^2 = "
            "[2a^2+(2-a)b^2-b^4/a]t^2 + "
            "2ab*kappa*t*v + (2-a)(a^2+b^2)v^2"
        ),
        sufficient_condition="beta=(b/v)^2*t^2/a < 2-a",
    )


def _validate_refinement_receipt(calc):
    raw = json.loads(REFINEMENT_RECEIPT.read_text(encoding="utf-8"))
    if raw.get("state_sha256") != calc.state_hash:
        raise ValueError("Existing refinement receipt does not match loaded core")
    if raw.get("subdivisions") != 128:
        raise ValueError("Expected immutable 128-way refinement receipt")
    if raw.get("counts", {}).get("certified") != 127:
        raise ValueError("Expected exactly 127 certified refinement cells")
    if raw.get("counts", {}).get("not certified: correlated kappa interval crosses 2") != 1:
        raise ValueError("Expected one unresolved transition cell")
    for name, digest in raw.get("input_hashes", {}).items():
        path = HERE / name
        if not path.exists() or _hash(path) != digest:
            raise ValueError("Refinement dependency hash changed: " + name)
    return raw


def certify(cuts=("1e-70", "1e-80")):
    """Build the directed three-piece certificate for the final parent cell."""

    calc = IntervalComparisonJets(4)
    _validate_frozen_switch(calc)
    refinement = _validate_refinement_receipt(calc)
    bridge, cells = load_receipts(calc)
    continuation = continue_exit(
        calc, bridge, cells["endpoint"], cells["initial_phi"], target_R="100"
    )
    c = calc.ctx
    with mp.workdps(calc.precision + 40):
        state, left = _prefix_state(calc, continuation)
        hb = c.mpf(".005")
        # Evaluate the complete unresolved subcell once with the raw physical
        # Tz split.  This is the branch-independent certificate candidate.
        full_x = c.mpf([endpoints(left)[0], endpoints(hb)[1]])
        _, full_raw_row = _advance_full_cell_with_raw_tz(
            calc, continuation, state, full_x, hb - left
        )
        full_branch = _both_branch_certificate(calc, full_raw_row)
        first_cut = hb - c.mpf(cuts[0])
        second_cut = hb - c.mpf(cuts[1])
        if not endpoints(left)[0] < endpoints(first_cut)[0] < endpoints(second_cut)[0] < endpoints(hb)[1]:
            raise ValueError("Branch cuts are not strictly ordered inside final cell")

        rows = []
        for label, right in (
            ("large_blend_branch", first_cut),
            ("transition_branch_gap", second_cut),
            ("small_blend_branch", hb),
        ):
            x = c.mpf([endpoints(left)[0], endpoints(right)[1]])
            width = right - left
            state, row = _advance_stage1(
                calc, continuation, state, x, width, inspect=True
            )
            rows.append(_compact_row(calc, label, x, width, row))
            left = right

        large = rows[0]["cone"]
        middle = rows[1]["cone"]
        small = rows[2]["cone"]
        small_conditional = rows[2]["conditional_small_branch"]
        result = dict(
            center_family=list(calc.center_family),
            state_sha256=calc.state_hash,
            accepted_schedule_sha256=ACCEPTED_SHA,
            source_final_cell_log_interval=[".00499755859375", ".005"],
            cuts_from_endpoint=list(cuts),
            rows=rows,
            full_unresolved_cell_raw=full_raw_row,
            full_cell_both_branch_certificate=full_branch,
            large_blend_branch_certified=(
                large.get("status") == "certified"
                and large.get("branch") == "kappa>2"
            ),
            transition_branch_status=middle.get("status"),
            small_blend_branch_certified=(
                small.get("status") == "certified"
                and small.get("branch") == "kappa<=2"
                and small_conditional.get("conditional_relaxed_cone_certified", False)
            ),
            small_blend_conditional=small_conditional,
            final_cell_fully_certified=full_branch["both_conditional_branches_certified"],
            final_cell_relaxed_cone_certified=full_branch[
                "both_conditional_branches_certified"
            ],
            final_cell_strong_cone_certified=False,
            transition_gap_isolated=True,
            correlated_stage1_relations_preserved=True,
            no_uniform_refinement=True,
            no_point_sample_certification=True,
            no_core_regeneration=True,
            whole_switch_cone_certified=False,
            terminal_matching_complete=False,
            original_construction_parameter_errors_enclosed=False,
            temporal_recursion=False,
        )
    dependency_names = (
        Path(__file__).name,
        SWITCH_RECEIPT.name,
        "lei_ren_part1_paper_interval_exit_switch_cone_refinement.py",
        REFINEMENT_RECEIPT.name,
        "lei_ren_part1_paper_interval_exit_continuation_enclosure.py",
        "lei_ren_part1_paper_interval_exit_continuation_enclosure.json",
        "lei_ren_part1_paper_interval_exit_stress_enclosure.py",
        "lei_ren_part1_paper_interval_comparison_enclosure.py",
        "lei_ren_part1_paper_interval_comparison_jets.py",
    )
    result["input_hashes"] = {name: _hash(HERE / name) for name in dependency_names}
    result["source_sha256"] = _hash(Path(__file__))
    output = HERE / "lei_ren_part1_paper_interval_exit_switch_branch_certificate.json"
    output.write_text(
        json.dumps(encode(_pack(result)), indent=2) + "\n", encoding="utf-8"
    )
    print(
        "Final stage-1 branch certificate:",
        result["large_blend_branch_certified"],
        result["transition_branch_status"],
        result["small_blend_branch_certified"],
        flush=True,
    )
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--large-cut", default="1e-70")
    parser.add_argument("--small-cut", default="1e-80")
    args = parser.parse_args()
    certify((args.large_cut, args.small_cut))


if __name__ == "__main__":
    main()
