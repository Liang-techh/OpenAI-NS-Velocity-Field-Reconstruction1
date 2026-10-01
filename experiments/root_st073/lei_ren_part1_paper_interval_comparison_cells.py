"""Cellwise directed boxes for the Section 9.25 comparison field.

The first half of the switch interval is represented directly by the full
finite-core plus analytic-tail box.  The second half uses the integrating
factor form of the comparison equations.  Each returned cell is a box over
the complete y-cell; its moment fields therefore include the partial-cell
range ``[0, step_upper]`` rather than only an endpoint value.

This module is limited to the stored center family ``[.49,.51]``.  It does
not solve the physical exit ODE or assert a terminal matching or cone
certificate.
"""

import hashlib
import json
import operator
from pathlib import Path

import mpmath as mp

from lei_ren_part1_paper_candidate_shared_inlet import normalized_inlet
from lei_ren_part1_paper_interval_comparison_enclosure import (
    alpha_box,
    core_box,
    tail_jet,
)
from lei_ren_part1_paper_interval_comparison_jets import IntervalComparisonJets
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode


STATE_KEYS = (
    "phi",
    "U",
    "theta",
    "z",
    "theta_z",
    "p",
    "u_squared",
    "weighted_phi_squared",
)
MOMENT_KEYS = STATE_KEYS[2:]
EXPECTED_CENTER_FAMILY = [".49", ".51"]


def _zero_like(value):
    """Return a Taylor zero with the same retained axial order."""

    return value * 0


def _interval_y(c, left, right):
    """Build a directed interval from two endpoint-compatible values."""

    return c.mpf([endpoints(left)[0], endpoints(right)[1]])


def _check_order_two(state, label):
    for key in STATE_KEYS:
        value = state[key]
        if not isinstance(value, IntervalTaylor) or value.order != 2:
            raise ValueError(f"{label} field {key} must be an order-2 Taylor jet")


def _check_positive_phi(state, label):
    lower = endpoints(state["phi"][0])[0]
    if lower <= 0:
        raise ValueError(f"{label} loses positive comparison phi")
    return lower


def _driver_packet(calc, state, scaled_radius, *, label):
    """Compute the normalized inlet driver for a complete state box.

    The endpoint override follows the established comparison enclosure: the
    propagated comparison field supplies ``phi_exit`` and ``u_exit`` while
    the radial derivative slots are zero.  D and the physical-scaled I_z are
    retained through axial order 1 for the bridge consumer.
    """

    _check_order_two(state, label)
    _check_positive_phi(state, label)
    zero = _zero_like(state["phi"])
    inlet = normalized_inlet(
        calc.phi,
        calc.u,
        state,
        calc.S.truncate(2),
        calc.ell.truncate(2),
        calc.p0.truncate(2),
        calc.lam,
        scaled_radius,
        calc.z.truncate(2),
        calc.delta,
        endpoints_override=dict(
            phi_exit=state["phi"],
            u_exit=state["U"],
            phi_s=zero,
            u_s=zero,
        ),
    )
    return dict(
        D=inlet["ratio"].truncate(1),
        I_z=(inlet["iz"] * calc.ctx.sqrt(calc.eps)).truncate(1),
        pressure=inlet["pressure"].truncate(2),
    )


def _cell_record(calc, *, index, phase, y, step, scaled_radius, state):
    packet = _driver_packet(calc, state, scaled_radius, label=f"cell {index}")
    return dict(
        index=index,
        phase=phase,
        y=y,
        step=step,
        scaled_radius=scaled_radius,
        state=state,
        D=packet["D"],
        I_z=packet["I_z"],
        pressure=packet["pressure"],
    )


def _rhs_ranges(calc, f, u, scaled_radius):
    """Moment right-hand-side ranges at one switched cell."""

    return dict(
        theta=f * (2 * scaled_radius * scaled_radius),
        z=u * scaled_radius,
        theta_z=f * u * (2 * scaled_radius * scaled_radius),
        p=f * f * scaled_radius,
        u_squared=u * u * scaled_radius,
        weighted_phi_squared=f * f * scaled_radius * scaled_radius,
    )


def _switch_cell(calc, *, index, y, step, start, G, J, prefix, ep, dp, du, hb):
    """Enclose one smooth-switch cell and return its endpoint prefix."""

    c = calc.ctx
    scaled_radius = 4 * c.exp(y)
    source_phi = calc.radial(calc.phi, scaled_radius).truncate(2) + ep
    source_lower = endpoints(source_phi[0])[0]
    if source_lower <= 0:
        raise ValueError(f"cell {index} source phi denominator crosses zero")
    alpha = alpha_box(c, y, hb)
    g = (calc.radial(calc.phi, scaled_radius, 1).truncate(2) + dp)
    g = g * scaled_radius * alpha / source_phi
    j = (calc.radial(calc.u, scaled_radius, 1).truncate(2) + du)
    j = j * scaled_radius * alpha

    # This is a range over every point of the cell, not an endpoint step.
    partial = c.mpf([0, endpoints(step)[1]])
    f = start["phi"] * (G + g * partial).exp()
    u = start["U"] + J + j * partial
    rhs = _rhs_ranges(calc, f, u, scaled_radius)
    state = dict(phi=f, U=u)
    state.update({key: prefix[key] + rhs[key] * partial for key in MOMENT_KEYS})
    _check_order_two(state, f"cell {index}")
    _check_positive_phi(state, f"cell {index}")

    record = _cell_record(
        calc,
        index=index,
        phase="switched",
        y=y,
        step=step,
        scaled_radius=scaled_radius,
        state=state,
    )

    # The prefix is advanced by the directed cell width, while ``state``
    # above retains the broader partial-cell range for all interior points.
    next_prefix = {
        key: prefix[key] + rhs[key] * step for key in MOMENT_KEYS
    }
    next_G = G + g * step
    next_J = J + j * step
    return record, next_G, next_J, next_prefix, source_lower


def _validate_calc(calc):
    family = list(calc.center_family)
    if family != EXPECTED_CENTER_FAMILY:
        raise ValueError(
            "Cellwise comparison is restricted to center family [.49,.51]; "
            f"received {family!r}"
        )
    if not calc.tail["target_met"]:
        raise ValueError("Cellwise comparison requires the analytic tail gate")
    if calc.degree != 124:
        raise ValueError("Cellwise comparison requires completed radial degree 124")


def pack(value):
    """Convert Taylor jets recursively before directed JSON encoding."""

    if isinstance(value, IntervalTaylor):
        return list(value.coefficients)
    if isinstance(value, dict):
        return {key: pack(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [pack(item) for item in value]
    return value


def build_cells(calc, cells_per_half=16):
    """Return directed comparison boxes over y in [0, 2 h_b].

    Parameters
    ----------
    calc:
        A production ``IntervalComparisonJets`` instance with center family
        ``[".49", ".51"]`` and a passing degree-124 tail receipt.
    cells_per_half:
        Positive integer number of cells in each half of the switch.

    Returns
    -------
    dict
        ``cells`` is chronological and contains complete state/driver boxes;
        ``endpoint`` is a separate exact switched-endpoint packet.  The
        top-level ``initial_phi`` is the full finite-plus-tail Phi box at
        scaled radius 4, retained for Fa normalization by the bridge.
    """

    cells_per_half = operator.index(cells_per_half)
    if cells_per_half < 1:
        raise ValueError("cells_per_half must be positive")
    _validate_calc(calc)
    c = calc.ctx
    with mp.workdps(calc.precision + 40):
        hb = c.mpf(".005")
        step = hb / cells_per_half
        zero_y = c.mpf(0)
        initial = core_box(calc, c.mpf(4))
        initial_phi = initial["phi"]
        _check_order_two(initial, "initial core")
        _check_positive_phi(initial, "initial core")

        cells = []
        # Direct full-core boxes cover [0, h_b].  Every y cell is retained as
        # an interval, so the tail and all axial coefficients remain directed.
        for n in range(cells_per_half):
            left = hb * c.mpf(n) / cells_per_half
            right = hb * c.mpf(n + 1) / cells_per_half
            y = _interval_y(c, left, right)
            scaled_radius = 4 * c.exp(y)
            state = core_box(calc, scaled_radius)
            cells.append(_cell_record(
                calc,
                index=n,
                phase="core",
                y=y,
                step=step,
                scaled_radius=scaled_radius,
                state=state,
            ))

        # Switch cells use the explicit integrating-factor solution.  Prefix
        # moments are endpoint enclosures; each cell's state uses partial as
        # [0, step_upper] and therefore covers all points in that cell.
        start = core_box(calc, 4 * c.exp(hb))
        _check_order_two(start, "switch start")
        _check_positive_phi(start, "switch start")
        zero = start["phi"] * 0
        G = zero
        J = zero
        prefix = {key: start[key] for key in MOMENT_KEYS}
        ep = tail_jet(calc, "Phi_tail")
        dp = tail_jet(calc, "Phi_tail", radial_order=1)
        du = tail_jet(calc, "Psi_tail", radial_order=1, scale=calc.eps)
        source_lower_min = None
        for n in range(cells_per_half):
            left = hb + hb * c.mpf(n) / cells_per_half
            right = hb + hb * c.mpf(n + 1) / cells_per_half
            y = _interval_y(c, left, right)
            record, G, J, prefix, source_lower = _switch_cell(
                calc,
                index=cells_per_half + n,
                y=y,
                step=step,
                start=start,
                G=G,
                J=J,
                prefix=prefix,
                ep=ep,
                dp=dp,
                du=du,
                hb=hb,
            )
            source_lower_min = source_lower if source_lower_min is None else min(
                source_lower_min, source_lower
            )
            cells.append(record)

        endpoint_scaled_radius = 4 * c.exp(2 * hb)
        endpoint_state = dict(
            phi=start["phi"] * G.exp(),
            U=start["U"] + J,
            **prefix,
        )
        _check_order_two(endpoint_state, "switched endpoint")
        _check_positive_phi(endpoint_state, "switched endpoint")
        endpoint_driver = _driver_packet(
            calc,
            endpoint_state,
            endpoint_scaled_radius,
            label="switched endpoint",
        )
        endpoint = dict(
            index=len(cells),
            phase="switched_endpoint",
            y=c.mpf([endpoints(2 * hb)[0], endpoints(2 * hb)[1]]),
            step=c.mpf(0),
            scaled_radius=endpoint_scaled_radius,
            state=endpoint_state,
            D=endpoint_driver["D"],
            I_z=endpoint_driver["I_z"],
            pressure=endpoint_driver["pressure"],
        )
        return dict(
            center_family=list(calc.center_family),
            state_sha256=calc.state_hash,
            completed_degree=calc.degree,
            cells_per_half=cells_per_half,
            hb=hb,
            y_interval=c.mpf([0, endpoints(2 * hb)[1]]),
            initial_phi=initial_phi,
            cells=cells,
            endpoint=endpoint,
            minimum_switch_source_phi_lower=source_lower_min,
            retained_axial_order=2,
            driver_axial_order=1,
            analytic_core_tail_propagated=True,
            cellwise_ranges_enclosed=True,
            full_core_first_half=True,
            integrating_factor_second_half=True,
            endpoint_is_separate=True,
            source_parameter_errors_enclosed=False,
            physical_exit_ODE_solved=False,
            terminal_five_moment_closure=False,
            stress_cone_certified=False,
        )


def run(cells_per_half=16):
    """Build one bounded production receipt for the stored center family."""

    calc = IntervalComparisonJets(4)
    result = build_cells(calc, cells_per_half)
    here = Path(__file__).resolve().parent
    result["input_hashes"] = {
        name: hashlib.sha256((here / name).read_bytes()).hexdigest()
        for name in (
            "lei_ren_part1_paper_interval_comparison_cells.py",
            "lei_ren_part1_paper_interval_comparison_enclosure.py",
            "lei_ren_part1_paper_interval_comparison_jets.py",
            "lei_ren_part1_paper_interval_core_inlet_bound.py",
            "lei_ren_part1_paper_candidate_shared_inlet.py",
            "lei_ren_part1_paper_interval_taylor.py",
        )
    }
    Path(__file__).with_suffix(".json").write_text(
        json.dumps(encode(pack(result)), indent=2) + "\n", encoding="utf-8"
    )
    print(
        "Directed comparison cells:",
        len(result["cells"]),
        "cells; endpoint separated; axial field order 2, driver order 1",
        flush=True,
    )
    return result


if __name__ == "__main__":
    run()
