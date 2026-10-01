"""Cellwise radial cone atlas for the controlled post-collar continuation.

This module reuses the frozen continuation receipt.  It does not reload or
recompute a core at each radius: the stored global ``delta_phi_range`` and
``delta_u_range`` boxes are valid at every radius in the continuation range.
Each logarithmic radial cell gets a fresh frozen comparison state and a fresh
directed inlet-driver evaluation.  An unresolved cell is retained as
``not certified``; it is never interpreted as a physical cone failure.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import time
from pathlib import Path

import mpmath as mp

from lei_ren_part1_paper_candidate_pressure_function import ACCEPTED_SHA
from lei_ren_part1_paper_interval_comparison_jets import IntervalComparisonJets
from lei_ren_part1_paper_interval_exit_continuation_enclosure import (
    MOMENT_KEYS,
    _driver,
    _frozen_state,
    _moment_errors,
    restore_jet,
    restore_state,
    restore_value,
)
from lei_ren_part1_paper_interval_exit_stress_enclosure import evaluate as evaluate_stress
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode


HERE = Path(__file__).resolve().parent
CONTINUATION_RECEIPT = HERE / "lei_ren_part1_paper_interval_exit_continuation_enclosure.json"


def _pack(value):
    """Turn IntervalTaylor objects into recursively encodable coefficient lists."""

    if hasattr(value, "coefficients") and hasattr(value, "order"):
        return [_pack(item) for item in value.coefficients]
    if isinstance(value, dict):
        return {key: _pack(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_pack(item) for item in value]
    return value


def _hash(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _nonnegative(ctx, value):
    """Clamp a directed difference at zero without losing its upper endpoint."""

    lo, hi = endpoints(value)
    if hi <= 0:
        return ctx.mpf(0)
    return ctx.mpf([max(mp.mpf(0), lo), hi])


def _log_text(value):
    return mp.nstr(mp.mpf(value), 120)


def _exp_box(ctx, left, right, lower, upper):
    """Directed exp box, clipped to the intended global radial interval."""

    raw = ctx.exp(ctx.mpf([left, right]))
    lo = max(mp.mpf(lower), endpoints(raw)[0])
    hi = min(mp.mpf(upper), endpoints(raw)[1])
    if lo > hi:
        raise ValueError("logarithmic cell does not intersect radial interval")
    return ctx.mpf([lo, hi])


def _restore_continuation(calc, path=CONTINUATION_RECEIPT):
    """Load and validate the immutable continuation packet once."""

    raw = json.loads(Path(path).read_text(encoding="utf-8"))
    if raw.get("state_sha256") != calc.state_hash:
        raise ValueError("Continuation state hash does not match the loaded core")
    if raw.get("accepted_schedule_sha256") != ACCEPTED_SHA:
        raise ValueError("Continuation accepted pressure datum changed")
    if raw.get("target_R") != "100":
        raise ValueError("Atlas is restricted to the stored target R = 100")
    for name, digest in raw.get("input_hashes", {}).items():
        dep = HERE / name
        if not dep.exists() or _hash(dep) != digest:
            raise ValueError("Continuation dependency hash changed: " + name)
    if not raw.get("controlled_post_collar_range", False):
        raise ValueError("Continuation range is not marked controlled")
    if not raw.get("physical_exit_ODE_error_enclosed", False):
        raise ValueError("Continuation does not enclose the post-collar ODE error")

    c = calc.ctx
    data = dict(
        raw=raw,
        state_sha256=raw["state_sha256"],
        accepted_schedule_sha256=raw["accepted_schedule_sha256"],
        endpoint_scaled_radius=restore_value(c, raw["endpoint_scaled_radius"]),
        target_scaled_radius=restore_value(c, raw["target_scaled_radius"]),
        target_physical_R=restore_value(c, raw["target_physical_R"]),
        exit_shear_epsilon=restore_value(c, raw["exit_shear_epsilon"]),
        endpoint_actual_state=restore_state(c, raw["endpoint_actual_state"], order=1),
        endpoint_comparison_state=restore_state(c, raw["endpoint_comparison_state"], order=2),
        delta_phi_range=restore_jet(c, raw["delta_phi_range"], order=1),
        delta_u_range=restore_jet(c, raw["delta_u_range"], order=1),
    )
    return data


def _cell_packet(calc, continuation, scaled_radius):
    """Build one actual/comparison state pair and its local stress packet."""

    c = calc.ctx
    se = continuation["endpoint_scaled_radius"]
    actual_e = continuation["endpoint_actual_state"]
    comparison_e = continuation["endpoint_comparison_state"]
    phi_e, u_e = actual_e["phi"], actual_e["U"]
    bar_phi, bar_u = comparison_e["phi"], comparison_e["U"]
    ds = _nonnegative(c, scaled_radius - se)
    ds2 = _nonnegative(c, scaled_radius * scaled_radius - se * se)

    # These are the stored whole-path field-error boxes.  They enclose the
    # correction at every intermediate radius, so moment errors below remain
    # cumulative from the switched endpoint and need no new core evaluation.
    delta_phi = continuation["delta_phi_range"]
    delta_u = continuation["delta_u_range"]
    actual_base = _frozen_state(actual_e, phi_e, u_e, ds, ds2)
    actual_errors = _moment_errors(phi_e, u_e, delta_phi, delta_u, ds, ds2)
    actual = dict(phi=phi_e + delta_phi, U=u_e + delta_u)
    actual.update(
        {
            key: actual_base[key] + actual_errors[key]
            for key in MOMENT_KEYS
        }
    )

    comparison = _frozen_state(comparison_e, bar_phi, bar_u, ds, ds2)
    driver = _driver(calc, bar_phi, bar_u, comparison, scaled_radius)
    eps_exit = continuation["exit_shear_epsilon"]
    gy = -(driver["D"] * eps_exit / 2)
    uy = -(
        (actual["phi"] / bar_phi)
        * eps_exit
        * driver["I_z"]
        * c.sqrt(scaled_radius * calc.eps / 2)
    )
    packet = dict(
        normalized_exit_state=actual,
        physical_R=scaled_radius * calc.eps,
        actual_terminal_g_y=gy,
        actual_terminal_U_y=uy,
    )
    stress = evaluate_stress(calc, packet)
    return dict(
        stress=stress,
        packet=packet,
        actual_state=actual,
        comparison_state=comparison,
        driver=driver,
        ds=ds,
        ds2=ds2,
    )


def _cell_result(calc, continuation, cell, started, deadline):
    """Evaluate a cell or return an explicit full-coverage timeout record."""

    if time.monotonic() >= deadline:
        return dict(
            status="not evaluated: walltime budget",
            path=cell["path"],
            depth=cell["depth"],
            log_interval=[_log_text(cell["y0"]), _log_text(cell["y1"])],
            scaled_radius=_pack(cell["scaled_radius"]),
            physical_R=_pack(cell["scaled_radius"] * calc.eps),
        )
    try:
        result = _cell_packet(calc, continuation, cell["scaled_radius"])
        stress = result["stress"]
        cone = stress["cone"]
        status = cone.get("status", "not certified")
        return dict(
            status=status,
            path=cell["path"],
            depth=cell["depth"],
            log_interval=[_log_text(cell["y0"]), _log_text(cell["y1"])],
            scaled_radius=_pack(cell["scaled_radius"]),
            physical_R=_pack(result["packet"]["physical_R"]),
            driver_D=_pack(result["driver"]["D"]),
            driver_I_z=_pack(result["driver"]["I_z"]),
            phi=_pack(result["actual_state"]["phi"]),
            normalized_stress=_pack(stress["normalized_stress"]),
            cone=_pack(cone),
            relaxed_cone_certified=bool(cone.get("relaxed_cone_certified", False)),
            admissible_cone_certified=bool(cone.get("admissible_cone_certified", False)),
            evaluated=True,
            elapsed_seconds=time.monotonic() - started,
        )
    except Exception as exc:  # retain the cell in the atlas, never certify it
        return dict(
            status="evaluation error: " + str(exc),
            path=cell["path"],
            depth=cell["depth"],
            log_interval=[_log_text(cell["y0"]), _log_text(cell["y1"])],
            scaled_radius=_pack(cell["scaled_radius"]),
            physical_R=_pack(cell["scaled_radius"] * calc.eps),
            evaluated=False,
        )


def _needs_refinement(record):
    return record.get("evaluated", False) and not record.get(
        "admissible_cone_certified", False
    )


def _child_cell(calc, parent, left, right, suffix, radial_lower, radial_upper):
    return dict(
        y0=left,
        y1=right,
        scaled_radius=_exp_box(calc.ctx, left, right, radial_lower, radial_upper),
        depth=parent["depth"] + 1,
        path=parent["path"] + suffix,
    )


def _refine(calc, continuation, cell, started, deadline, max_depth, radial_lower, radial_upper):
    record = _cell_result(calc, continuation, cell, started, deadline)
    if (
        _needs_refinement(record)
        and cell["depth"] < max_depth
        and time.monotonic() < deadline
    ):
        mid = (cell["y0"] + cell["y1"]) / 2
        left = _child_cell(calc, cell, cell["y0"], mid, "L", radial_lower, radial_upper)
        right = _child_cell(calc, cell, mid, cell["y1"], "R", radial_lower, radial_upper)
        return (
            _refine(calc, continuation, left, started, deadline, max_depth, radial_lower, radial_upper)
            + _refine(calc, continuation, right, started, deadline, max_depth, radial_lower, radial_upper)
        )
    return [record]


def _coverage(records, y0, y1):
    """Check that retained leaf log intervals cover the requested range."""

    if not records:
        return False
    ordered = sorted(records, key=lambda item: mp.mpf(item["log_interval"][0]))
    tolerance = mp.mpf("1e-100")
    previous = None
    for record in ordered:
        left = mp.mpf(record["log_interval"][0])
        right = mp.mpf(record["log_interval"][1])
        if right < left:
            return False
        if previous is not None and left > previous + tolerance:
            return False
        previous = right if previous is None else max(previous, right)
    return (
        mp.mpf(ordered[0]["log_interval"][0]) <= y0 + tolerance
        and mp.mpf(ordered[-1]["log_interval"][1]) >= y1 - tolerance
    )


def _summary(calc, records, started, deadline, cells, max_depth, y0, y1):
    counts = dict(total=len(records), strong=0, relaxed=0, unresolved=0, errors=0)
    unresolved_locations = []
    min_margin = None
    max_shear_upper = None
    for record in records:
        if record.get("admissible_cone_certified"):
            counts["strong"] += 1
        if record.get("relaxed_cone_certified"):
            counts["relaxed"] += 1
        status = record.get("status", "")
        if status.startswith("evaluation error"):
            counts["errors"] += 1
        if not record.get("admissible_cone_certified", False):
            counts["unresolved"] += 1
            unresolved_locations.append(
                dict(
                    path=record.get("path"),
                    depth=record.get("depth"),
                    log_interval=record.get("log_interval"),
                    physical_R=record.get("physical_R"),
                    status=status,
                )
            )
        cone = record.get("cone", {})
        margin = cone.get("margin")
        if hasattr(margin, "_mpi_"):
            lo, hi = endpoints(margin)
            min_margin = lo if min_margin is None else min(min_margin, lo)
        stress = record.get("normalized_stress", {})
        shear = stress.get("S_theta_over_F")
        if hasattr(shear, "_mpi_"):
            hi = endpoints(shear)[1]
            max_shear_upper = hi if max_shear_upper is None else max(max_shear_upper, hi)
        elif isinstance(shear, list) and shear and hasattr(shear[0], "_mpi_"):
            hi = endpoints(shear[0])[1]
            max_shear_upper = hi if max_shear_upper is None else max(max_shear_upper, hi)
    return dict(
        counts=counts,
        full_coverage_partition=_coverage(records, y0, y1),
        initial_log_cells=cells,
        max_refinement_depth=max_depth,
        unresolved_locations=unresolved_locations,
        minimum_margin_lower=None if min_margin is None else mp.nstr(min_margin, 40),
        maximum_shear_upper=None if max_shear_upper is None else mp.nstr(max_shear_upper, 40),
        walltime_seconds=time.monotonic() - started,
        walltime_expired=time.monotonic() >= deadline,
    )


def _run_at_precision(calc, seconds=60, initial_cells=256, max_depth=4):
    started = time.monotonic()
    if initial_cells < 1 or max_depth < 0 or seconds <= 0:
        raise ValueError("seconds, initial_cells, and max_depth must be positive")
    continuation = _restore_continuation(calc)
    c = calc.ctx
    se = continuation["endpoint_scaled_radius"]
    target_s = continuation["target_scaled_radius"]
    radial_lower = endpoints(se)[0]
    radial_upper = endpoints(target_s)[1]
    y0 = endpoints(c.ln(se))[0]
    y1 = endpoints(c.ln(target_s))[1]
    if not y1 > y0:
        raise ValueError("Invalid continuation logarithmic interval")
    deadline = time.monotonic() + float(seconds)

    cells = []
    for index in range(initial_cells):
        left = y0 + (y1 - y0) * index / initial_cells
        right = y0 + (y1 - y0) * (index + 1) / initial_cells
        cells.append(
            dict(
                y0=left,
                y1=right,
                scaled_radius=_exp_box(c, left, right, radial_lower, radial_upper),
                depth=0,
                path=str(index),
            )
        )

    records = []
    next_cell = 0
    for index, cell in enumerate(cells):
        records.extend(
            _refine(
                calc,
                continuation,
                cell,
                started,
                deadline,
                max_depth,
                radial_lower,
                radial_upper,
            )
        )
        next_cell = index + 1
        if time.monotonic() >= deadline:
            break
    # Preserve complete radial coverage when a bounded run stops early.
    if next_cell < len(cells):
        for cell in cells[next_cell:]:
            records.append(
                dict(
                    status="not evaluated: walltime budget",
                    path=cell["path"],
                    depth=cell["depth"],
                    log_interval=[_log_text(cell["y0"]), _log_text(cell["y1"])],
                    scaled_radius=_pack(cell["scaled_radius"]),
                    physical_R=_pack(cell["scaled_radius"] * calc.eps),
                )
            )

    records.sort(key=lambda item: (mp.mpf(item["log_interval"][0]), item["path"]))
    result = dict(
        state_sha256=continuation["state_sha256"],
        accepted_schedule_sha256=continuation["accepted_schedule_sha256"],
        center_family=list(calc.center_family),
        Lambda=mp.nstr(calc.lam, 50),
        epsilon=mp.nstr(calc.eps, 50),
        endpoint_scaled_radius=_pack(se),
        target_scaled_radius=_pack(target_s),
        target_physical_R=_pack(continuation["target_physical_R"]),
        log_interval=[mp.nstr(y0, 60), mp.nstr(y1, 60)],
        cells=records,
        summary=_summary(calc, records, started, deadline, initial_cells, max_depth, y0, y1),
        core_loaded_once=True,
        global_field_error_boxes_reused=True,
        global_error_sources=["delta_phi_range", "delta_u_range"],
        cumulative_moment_errors_from_endpoint=True,
        point_sample_certification=False,
        whole_axis_or_whole_annulus_cone_certified=False,
        original_construction_parameter_errors_enclosed=False,
        temporal_recursion=False,
        terminal_matching_complete=False,
        stress_cone_certified=False,
        source_sha256=_hash(Path(__file__)),
        input_hashes={
            "lei_ren_part1_paper_interval_exit_continuation_enclosure.py": _hash(
                HERE / "lei_ren_part1_paper_interval_exit_continuation_enclosure.py"
            ),
            "lei_ren_part1_paper_interval_exit_continuation_enclosure.json": _hash(
                HERE / "lei_ren_part1_paper_interval_exit_continuation_enclosure.json"
            ),
            "lei_ren_part1_paper_interval_exit_stress_enclosure.py": _hash(
                HERE / "lei_ren_part1_paper_interval_exit_stress_enclosure.py"
            ),
            "lei_ren_part1_paper_interval_comparison_jets.py": _hash(
                HERE / "lei_ren_part1_paper_interval_comparison_jets.py"
            ),
            "lei_ren_part1_paper_candidate_shared_inlet.py": _hash(
                HERE / "lei_ren_part1_paper_candidate_shared_inlet.py"
            ),
            "lei_ren_part1_paper_candidate_pressure_function.py": _hash(
                HERE / "lei_ren_part1_paper_candidate_pressure_function.py"
            ),
            "lei_ren_part1_paper_schedule_endpoint_enclosures.py": _hash(
                HERE / "lei_ren_part1_paper_schedule_endpoint_enclosures.py"
            ),
            "lei_ren_part1_paper_schedule_endpoint_enclosures_check.py": _hash(
                HERE / "lei_ren_part1_paper_schedule_endpoint_enclosures_check.py"
            ),
        },
    )
    output = HERE / "lei_ren_part1_paper_interval_exit_radial_cone_atlas.json"
    output.write_text(json.dumps(encode(_pack(result)), indent=2) + "\n", encoding="utf-8")
    print(
        "Radial cone atlas:",
        result["summary"]["counts"],
        "walltime",
        mp.nstr(result["summary"]["walltime_seconds"], 8),
        "seconds",
        flush=True,
    )
    return result


def run(seconds=60, initial_cells=256, max_depth=4):
    calc = IntervalComparisonJets(4)
    with mp.workdps(calc.precision + 40):
        return _run_at_precision(
            calc,
            seconds=seconds,
            initial_cells=initial_cells,
            max_depth=max_depth,
        )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seconds", type=float, default=60.0)
    parser.add_argument("--cells", type=int, default=256)
    parser.add_argument("--max-depth", type=int, default=4)
    args = parser.parse_args()
    run(seconds=args.seconds, initial_cells=args.cells, max_depth=args.max_depth)


if __name__ == "__main__":
    main()
