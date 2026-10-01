"""Directed post-switch continuation of the actual Section 9.25 exit.

The comparison field is frozen after the stored switch endpoint.  Its
moments therefore have exact constant-field increments, while the actual
field receives a controlled small-epsilon correction.  The continuation is
an interval range over every radius in ``[R_e, target_R]`` and also returns a
separate target packet.  It does not assert outer matching, a stress cone, or
the temporal construction.
"""

import hashlib
import json
from pathlib import Path

import mpmath as mp

from lei_ren_part1_paper_candidate_pressure_function import ACCEPTED_SHA
from lei_ren_part1_paper_candidate_shared_inlet import normalized_inlet
from lei_ren_part1_paper_interval_comparison_jets import IntervalComparisonJets
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode


HERE = Path(__file__).resolve().parent
BRIDGE_RECEIPT = HERE / "lei_ren_part1_paper_interval_exit_bridge_enclosure.json"
CELLS_RECEIPT = HERE / "lei_ren_part1_paper_interval_comparison_cells.json"

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


def restore_value(ctx, value):
    """Restore exact interval/scalar packets recursively."""

    if isinstance(value, IntervalTaylor):
        return value
    if hasattr(value, "_mpi_"):
        return value if getattr(value, "ctx", None) is ctx else ctx.mpf(value)
    if hasattr(value, "_mpf_"):
        return ctx.mpf(value)
    if isinstance(value, dict):
        if "lower_exact_mpf_tuple" in value:
            lo = mp.make_mpf(tuple(value["lower_exact_mpf_tuple"]))
            hi = mp.make_mpf(tuple(value["upper_exact_mpf_tuple"]))
            return ctx.mpf([lo, hi])
        if "exact_mpf_tuple" in value:
            return ctx.mpf(mp.make_mpf(tuple(value["exact_mpf_tuple"])))
        return {key: restore_value(ctx, item) for key, item in value.items()}
    if isinstance(value, list):
        return [restore_value(ctx, item) for item in value]
    return value


def restore_jet(ctx, value, order=None):
    """Restore an encoded coefficient list as an ``IntervalTaylor`` jet."""

    if isinstance(value, IntervalTaylor):
        return value.truncate(value.order if order is None else order)
    if isinstance(value, (list, tuple)):
        values = [restore_value(ctx, item) for item in value]
        if not values:
            raise ValueError("Taylor packet cannot be empty")
        if order is not None:
            values = values[: order + 1]
        return IntervalTaylor(ctx, values)
    return IntervalTaylor(ctx, [restore_value(ctx, value)])


def restore_state(ctx, raw, order=None):
    if not isinstance(raw, dict):
        raise TypeError("state packet must be a mapping")
    missing = [key for key in STATE_KEYS if key not in raw]
    if missing:
        raise KeyError("state packet missing " + ", ".join(missing))
    return {key: restore_jet(ctx, raw[key], order) for key in STATE_KEYS}


def _pack(value):
    if isinstance(value, IntervalTaylor):
        return list(value.coefficients)
    if isinstance(value, dict):
        return {key: _pack(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_pack(item) for item in value]
    return value


def _load_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def _validate_hashes(receipt, label):
    for name, digest in receipt.get("input_hashes", {}).items():
        path = HERE / name
        if not path.exists():
            raise FileNotFoundError(f"{label} dependency missing: {name}")
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        if actual != digest:
            raise ValueError(f"{label} dependency hash changed: {name}")


def _overlap(left, right):
    return endpoints(left)[0] <= endpoints(right)[1] and endpoints(right)[0] <= endpoints(left)[1]


def _nonnegative(ctx, value, label):
    lo, hi = endpoints(value)
    if hi < 0:
        raise ValueError(f"{label} is negative")
    return ctx.mpf([max(mp.mpf(0), lo), max(mp.mpf(0), hi)])


def _max_abs(ctx, value):
    lo, hi = endpoints(value)
    return ctx.mpf(max(abs(lo), abs(hi)))


def _error_jet(ctx, delta_g):
    """Enclose exp(delta_g)-1 without subtracting two nearly equal intervals."""

    t0 = _max_abs(ctx, delta_g[0])
    t1 = _max_abs(ctx, delta_g[1]) if delta_g.order >= 1 else ctx.mpf(0)
    factor = ctx.exp(t0)
    q0 = factor * t0
    q1 = factor * t1
    return IntervalTaylor(ctx, [ctx.mpf([-q0, q0]), ctx.mpf([-q1, q1])])


def _state_from_endpoint(ctx, raw, order=1):
    return restore_state(ctx, raw, order)


def _truncate_state(state, order=1):
    return {key: value.truncate(order) for key, value in state.items()}


def _frozen_state(endpoint_state, phi, u, ds, ds2):
    """Extend normalized frozen fields and cumulative moments."""

    f = phi
    v = u
    state = {
        "phi": f,
        "U": v,
        "theta": endpoint_state["theta"] + f * ds2,
        "z": endpoint_state["z"] + v * ds,
        "theta_z": endpoint_state["theta_z"] + f * v * ds2,
        "p": endpoint_state["p"] + f * f * ds,
        "u_squared": endpoint_state["u_squared"] + v * v * ds,
        "weighted_phi_squared": endpoint_state["weighted_phi_squared"]
        + f * f * ds2 / 2,
    }
    return state


def _driver(calc, bar_phi, bar_u, state, scaled_radius):
    """Evaluate D and physical-scaled I_z for a singleton frozen field."""

    if endpoints(bar_phi[0])[0] <= 0:
        raise ValueError("Frozen comparison phi is not strictly positive")
    if any(state[key].order < 2 for key in MOMENT_KEYS):
        raise ValueError("Frozen driver requires C2 comparison moment jets")
    zero = bar_phi * 0
    inlet = normalized_inlet(
        [bar_phi],
        [bar_u],
        state,
        calc.S.truncate(2),
        calc.ell.truncate(2),
        calc.p0.truncate(2),
        calc.lam,
        scaled_radius,
        calc.z.truncate(2),
        calc.delta,
        endpoints_override=dict(
            phi_exit=bar_phi,
            u_exit=bar_u,
            phi_s=zero,
            u_s=zero,
        ),
    )
    return dict(
        D=inlet["ratio"].truncate(1),
        I_z=(inlet["iz"] * calc.ctx.sqrt(calc.eps)).truncate(1),
        pressure=inlet["pressure"].truncate(1),
    )


def _moment_errors(phi_e, u_e, delta_phi, delta_u, ds, ds2):
    return dict(
        theta=delta_phi * ds2,
        z=delta_u * ds,
        theta_z=(phi_e * delta_u + u_e * delta_phi + delta_phi * delta_u) * ds2,
        p=(2 * phi_e * delta_phi + delta_phi * delta_phi) * ds,
        u_squared=(2 * u_e * delta_u + delta_u * delta_u) * ds,
        weighted_phi_squared=(2 * phi_e * delta_phi + delta_phi * delta_phi)
        * ds2
        / 2,
    )


def _physical_packet(calc, state, scaled_radius):
    """Convert normalized C1 state to physical fields and moments."""

    c = calc.ctx
    amplitude = IntervalTaylor(c, [calc.F0, calc.F0 * calc.ell[0]])
    S_true = calc.S.truncate(1)
    P0_true = calc.p0.truncate(1)
    F = amplitude * state["phi"]
    physical_moments = dict(
        theta=amplitude * state["theta"] * calc.eps ** 2,
        z=state["z"] * calc.eps,
        theta_z=amplitude * state["theta_z"] * calc.eps ** 2,
        z_theta=state["u_squared"] * calc.eps
        - S_true * state["weighted_phi_squared"] * calc.eps ** 2,
        p=S_true * state["p"] * calc.eps,
    )
    P = P0_true + physical_moments["p"]
    R = scaled_radius * calc.eps
    z = calc.z[0]
    Mz = physical_moments["z"]
    Ur = (
        2 * z * R * state["U"][0]
        - (1 - calc.delta) * z * Mz[0]
        - (1 - z * z) * Mz[1]
    ) / ((1 - calc.delta * z * z) * c.sqrt(2 * R))
    return dict(
        physical_R=R,
        F=F,
        Uz=state["U"],
        P=P,
        Ur_value=Ur,
        physical_moments=physical_moments,
    )


def _validate_receipts(calc, bridge, cells):
    if bridge.get("state_sha256") != calc.state_hash:
        raise ValueError("Bridge receipt state hash does not match calc")
    if cells.get("state_sha256") != calc.state_hash:
        raise ValueError("Comparison-cell receipt state hash does not match calc")
    if bridge.get("accepted_schedule_sha256") != ACCEPTED_SHA:
        raise ValueError("Bridge accepted pressure datum changed")
    if cells.get("center_family") != list(calc.center_family):
        raise ValueError("Comparison center family changed")
    if bridge.get("center_family") != list(calc.center_family):
        raise ValueError("Bridge center family changed")
    if bridge.get("completed_core_degree") != calc.degree or cells.get("completed_degree") != calc.degree:
        raise ValueError("Receipt radial degree does not match calc")
    if bridge.get("cells_per_half") != cells.get("cells_per_half"):
        raise ValueError("Bridge and comparison cell partitions differ")
    _validate_hashes(bridge, "bridge")
    _validate_hashes(cells, "comparison")


def load_receipts(calc, bridge_path=BRIDGE_RECEIPT, cells_path=CELLS_RECEIPT):
    """Load, validate, and restore the two frozen production receipts."""

    bridge_raw = _load_json(bridge_path)
    cells_raw = _load_json(cells_path)
    _validate_receipts(calc, bridge_raw, cells_raw)
    c = calc.ctx
    bridge = restore_value(c, bridge_raw)
    cells = restore_value(c, cells_raw)
    bridge["normalized_exit_state"] = _state_from_endpoint(
        c, bridge_raw["normalized_exit_state"], order=1
    )
    cells["initial_phi"] = restore_jet(c, cells_raw["initial_phi"], order=2)
    cells["endpoint"]["state"] = _state_from_endpoint(
        c, cells_raw["endpoint"]["state"], order=2
    )
    cells["endpoint"]["D"] = restore_jet(c, cells_raw["endpoint"]["D"], order=1)
    cells["endpoint"]["I_z"] = restore_jet(c, cells_raw["endpoint"]["I_z"], order=1)
    return bridge, cells


def continue_exit(calc, bridge_data, comparison_endpoint, initial_phi, target_R="100"):
    """Continue the actual exit from the switched endpoint to ``target_R``.

    ``bridge_data`` and ``comparison_endpoint`` may be raw JSON mappings or
    mappings already restored by :func:`load_receipts`.  All returned field,
    moment, and driver entries are C1 ``IntervalTaylor`` objects unless
    explicitly marked as scalar intervals.
    """

    c = calc.ctx
    bridge = restore_value(c, bridge_data)
    endpoint_raw = comparison_endpoint
    if isinstance(endpoint_raw, dict) and "state" in endpoint_raw:
        endpoint = restore_value(c, endpoint_raw)
    else:
        raise TypeError("comparison_endpoint must contain a state packet")
    phi_initial = restore_jet(c, initial_phi, order=2)
    if endpoints(phi_initial[0])[0] <= 0:
        raise ValueError("Initial Phi(4) is not strictly positive")

    with mp.workdps(calc.precision + 40):
        phi_e = restore_jet(c, bridge["normalized_exit_state"]["phi"], order=1)
        u_e = restore_jet(c, bridge["normalized_exit_state"]["U"], order=1)
        actual_e = {
            key: restore_jet(c, bridge["normalized_exit_state"][key], order=1)
            for key in STATE_KEYS
        }
        if endpoints(phi_e[0])[0] <= 0:
            raise ValueError("Actual switched-endpoint Phi is not strictly positive")
        # Keep the comparison state through C2: normalized_inlet differentiates
        # its moments once, so C2 moments are required to retain C1 D/I_z.
        bar_e = _state_from_endpoint(c, endpoint["state"], order=2)
        bar_phi = bar_e["phi"]
        bar_u = bar_e["U"]
        se = restore_value(c, endpoint["scaled_radius"])
        eps_exit = restore_value(c, bridge["exit_shear_epsilon"])
        if endpoints(eps_exit)[0] <= 0 or endpoints(eps_exit)[1] >= 1:
            raise ValueError("Exit shear epsilon is outside (0,1)")
        if endpoints(se)[0] <= 0:
            raise ValueError("Comparison endpoint scaled radius is nonpositive")

        # The bridge and comparison receipts must describe the same physical
        # endpoint and the same radial scaling.
        bridge_R = restore_value(c, bridge["physical_R"])
        if not _overlap(bridge_R, se * calc.eps):
            raise ValueError("Bridge and comparison endpoint radii disagree")
        receipt_eps = restore_value(c, bridge["radial_scale_epsilon"])
        if not _overlap(receipt_eps, calc.eps):
            raise ValueError("Bridge radial epsilon disagrees with calc")

        target_physical_R = c.mpf(target_R)
        if endpoints(target_physical_R)[0] <= endpoints(bridge_R)[1]:
            raise ValueError("target_R must be beyond the switched endpoint")
        st = target_physical_R / calc.eps
        s_range = c.mpf([endpoints(se)[0], endpoints(st)[1]])
        ds_target = _nonnegative(c, st - se, "target scaled-radius increment")
        ds2_target = _nonnegative(c, st * st - se * se, "target squared-radius increment")
        ds_range = c.mpf([0, endpoints(ds_target)[1]])
        ds2_range = c.mpf([0, endpoints(ds2_target)[1]])
        dy_target = _nonnegative(c, c.ln(st / se), "target logarithmic increment")
        partial_dy = c.mpf([0, endpoints(dy_target)[1]])

        # Frozen comparison moments, both over the complete radius range and
        # at the requested target.  No old center tensor is extrapolated.
        bar_target = _frozen_state(bar_e, bar_phi, bar_u, ds_target, ds2_target)
        bar_range = _frozen_state(bar_e, bar_phi, bar_u, ds_range, ds2_range)
        driver_target = _driver(calc, bar_phi, bar_u, bar_target, st)
        driver_range = _driver(calc, bar_phi, bar_u, bar_range, s_range)

        # g_delta is the path integral of D over the whole continuation
        # interval.  The target driver is reserved for the terminal slope;
        # the range driver encloses D(s) at every intermediate radius.  Use a
        # symmetric explicit enclosure for exp(g_delta)-1 to avoid
        # subtracting nearly equal intervals around one.
        # Keep the Taylor jet on the left of interval scalars: ivmpf's
        # multiplication does not dispatch to the custom Taylor class.
        delta_g_target = -(driver_range["D"] * eps_exit * dy_target / 2)
        delta_g_range = -(driver_range["D"] * eps_exit * partial_dy / 2)
        error_target = _error_jet(c, delta_g_target)
        error_range = _error_jet(c, delta_g_range)
        one_target = IntervalTaylor(c, [c.mpf(1), c.mpf(0)]) + error_target
        one_range = IntervalTaylor(c, [c.mpf(1), c.mpf(0)]) + error_range
        delta_phi_target = phi_e * error_target
        delta_phi_range = phi_e * error_range

        ratio = phi_e / bar_phi
        B_range = -(
            ratio
            * eps_exit
            * driver_range["I_z"]
            * c.sqrt(s_range * calc.eps / 2)
        )
        # U_delta is also a path integral.  Use the full-range I_z and
        # sqrt(R) box; the endpoint driver is used only for terminal U_y.
        delta_u_target = B_range * one_target * dy_target
        delta_u_range = B_range * one_range * partial_dy

        base_target = _frozen_state(actual_e, phi_e, u_e, ds_target, ds2_target)
        base_range = _frozen_state(actual_e, phi_e, u_e, ds_range, ds2_range)
        errors_target = _moment_errors(
            phi_e, u_e, delta_phi_target, delta_u_target, ds_target, ds2_target
        )
        errors_range = _moment_errors(
            phi_e, u_e, delta_phi_range, delta_u_range, ds_range, ds2_range
        )
        target_state = dict(phi=phi_e + delta_phi_target, U=u_e + delta_u_target)
        target_state.update({key: base_target[key] + errors_target[key] for key in MOMENT_KEYS})
        range_state = dict(phi=phi_e + delta_phi_range, U=u_e + delta_u_range)
        range_state.update({key: base_range[key] + errors_range[key] for key in MOMENT_KEYS})

        target_physical = _physical_packet(calc, target_state, st)
        range_physical = _physical_packet(calc, range_state, s_range)
        terminal_gy = -(driver_target["D"] * eps_exit / 2)
        terminal_uy = -(
            (target_state["phi"] / bar_phi)
            * eps_exit
            * driver_target["I_z"]
            * c.sqrt(st * calc.eps / 2)
        )
        terminal_Fy = target_physical["F"] * terminal_gy
        terminal_R = target_physical_R

        return dict(
            center_family=list(calc.center_family),
            state_sha256=calc.state_hash,
            accepted_schedule_sha256=ACCEPTED_SHA,
            endpoint_scaled_radius=se,
            endpoint_physical_R=bridge_R,
            target_scaled_radius=st,
            target_physical_R=terminal_R,
            target_R=str(target_R),
            radial_range=dict(s=s_range, dy=partial_dy),
            exit_shear_epsilon=eps_exit,
            initial_phi=phi_initial,
            endpoint_actual_state=actual_e,
            endpoint_comparison_state=bar_e,
            comparison_target_state=bar_target,
            comparison_range_state=bar_range,
            comparison_target_driver=driver_target,
            comparison_range_driver=driver_range,
            target_g_integral_uses_range_driver=True,
            target_u_integral_uses_range_driver=True,
            terminal_derivatives_use_target_driver=True,
            delta_g_target=delta_g_target,
            delta_g_range=delta_g_range,
            delta_phi_target=delta_phi_target,
            delta_phi_range=delta_phi_range,
            delta_u_target=delta_u_target,
            delta_u_range=delta_u_range,
            moment_error_target=errors_target,
            moment_error_range=errors_range,
            normalized_exit_state=target_state,
            normalized_range_state=range_state,
            F=target_physical["F"],
            Uz=target_physical["Uz"],
            Ur_value=target_physical["Ur_value"],
            P=target_physical["P"],
            physical_moments=target_physical["physical_moments"],
            physical_range=range_physical,
            actual_terminal_g_y=terminal_gy,
            actual_terminal_U_y=terminal_uy,
            actual_F_y=terminal_Fy,
            actual_F_R=terminal_Fy / terminal_R,
            actual_Uz_R=terminal_uy / terminal_R,
            retained_axial_order=1,
            comparison_fields_frozen_after_switch=True,
            constant_field_moment_increments_exact=True,
            controlled_post_collar_range=True,
            physical_exit_ODE_error_enclosed=True,
            Ur_Z_available=False,
            construction_parameter_errors_enclosed=False,
            switch_100_110_enclosed=False,
            terminal_matching_complete=False,
            stress_cone_certified=False,
            temporal_recursion=False,
            whole_axis_transition=False,
            terminal_five_moment_closure=False,
            midpoint_projection_used=False,
        )


def run(target_R="100"):
    calc = IntervalComparisonJets(4)
    bridge, cells = load_receipts(calc)
    result = continue_exit(
        calc,
        bridge,
        cells["endpoint"],
        cells["initial_phi"],
        target_R=target_R,
    )
    result["input_hashes"] = {
        name: hashlib.sha256((HERE / name).read_bytes()).hexdigest()
        for name in (
            "lei_ren_part1_paper_interval_exit_continuation_enclosure.py",
            "lei_ren_part1_paper_interval_exit_bridge_enclosure.json",
            "lei_ren_part1_paper_interval_comparison_cells.json",
        )
    }
    Path(__file__).with_suffix(".json").write_text(
        json.dumps(encode(_pack(result)), indent=2) + "\n", encoding="utf-8"
    )
    print(
        "Controlled post-switch continuation to physical R",
        target_R,
        "complete; C1 target/range boxes retained",
        flush=True,
    )
    return result


if __name__ == "__main__":
    run()
