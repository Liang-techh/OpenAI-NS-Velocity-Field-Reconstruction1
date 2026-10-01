"""Separate resumable Lambda120 core using the unchanged gauge recurrence.

The exact-amplitude primitive supplies the anchored F0 interval. The accepted
physical pressure datum is unchanged. Existing Lambda48 artifacts are never
opened for writing. This finite axial-center calculation is not matching or
temporal coefficient recursion.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

import lei_ren_part1_paper_candidate_gauge_core as core
from lei_ren_part1_paper_uniform_axis_jets import uniform_axis_jets
from lei_ren_part1_paper_factored_core_positivity import squared_axis_rows

HERE = Path(__file__).resolve().parent
LAMBDA = "1e120"
PRESSURE = "lei_ren_part1_paper_candidate_pressure_axis_jets_refined.json"
AMPLITUDE_SOURCE = "lei_ren_part1_paper_candidate_exact_amplitude.py"
STATE_PATH = Path(__file__).with_name(Path(__file__).stem + "_state.json")
RECEIPT_PATH = Path(__file__).with_suffix(".json")


def target(pressure_identity):
    result = core._target(pressure_identity)
    result.update(Lambda=LAMBDA, driver_file=Path(__file__).name,
                  amplitude_source=AMPLITUDE_SOURCE,
                  amplitude_method="anchored rational primitive with enclosed complex poles")
    return result


def source_hashes(pressure_path):
    result = core._source_hashes(HERE, pressure_path)
    result[Path(core.__file__).name] = core._sha256(Path(core.__file__))
    result["driver"] = core._sha256(Path(__file__))
    result[AMPLITUDE_SOURCE] = core._sha256(HERE / AMPLITUDE_SOURCE)
    return result


def build_fixed(ctx, pressure):
    from lei_ren_part1_paper_candidate_exact_amplitude import exact_amplitude

    amplitude = exact_amplitude(ctx, Lambda=LAMBDA, Z=core.Z_CENTER)
    f0 = amplitude["F0_interval"]
    axis = uniform_axis_jets(ctx, radius=1, j=core.J, Lambda=LAMBDA,
                             logC=core.LOG_C, delta=core.DELTA,
                             length=core.AXIS_LENGTH,
                             axial_interval=ctx.mpf(core.Z_CENTER))
    if not core.endpoints(axis["F0_interval"])[0] <= core.endpoints(f0)[0]:
        raise AssertionError("anchored amplitude below the coarse amplitude enclosure")
    if not core.endpoints(f0)[1] <= core.endpoints(axis["F0_interval"])[1]:
        raise AssertionError("anchored amplitude above the coarse amplitude enclosure")
    gradient = list(axis["gradient_coefficients"])
    lam = ctx.mpf(LAMBDA)
    return dict(ell_Z_taylor=[-lam*v for v in gradient],
                S_Z_taylor=squared_axis_rows(ctx, gradient, lam, f0, core.AXIS_LENGTH),
                U0_Z_taylor=list(axis["U0"]), P0_Z_taylor=list(pressure),
                gradient_coefficients=gradient, F0_interval=f0)


def validate(state, hashes, pressure_sha, accepted_sha, identity, nested):
    expected = dict(target=target(identity), source_hashes=hashes,
                    input_pressure_receipt_sha256=pressure_sha,
                    accepted_schedule_sha256=accepted_sha,
                    nested_pressure_input_hashes=nested,
                    initial_count=core.AXIS_LENGTH)
    for key, value in expected.items():
        if state.get(key) != value:
            raise ValueError("Lambda120 resume dependency changed: " + key)
    degree = state["completed_radial_order"]
    if not isinstance(degree, int) or not 0 <= degree <= core.RADIAL_DEGREE:
        raise ValueError("invalid completed radial order")
    if state["next_radial_order"] != degree:
        raise ValueError("inconsistent next radial order")
    for name in ("A_rows", "Uz_rows", "P_rows"):
        if len(state[name]) != degree + 1:
            raise ValueError("inconsistent row count: " + name)
        for n, row in enumerate(state[name]):
            if len(row) != core.AXIS_LENGTH - n:
                raise ValueError("inconsistent axial depth: " + name)
    if not state["self_tests"]["all_checks_passed"]:
        raise ValueError("recurrence self-check evidence missing")


def run(seconds=180.0, pressure_file=PRESSURE):
    if seconds < 0:
        raise ValueError("seconds must be nonnegative")
    pressure_path = Path(pressure_file)
    if not pressure_path.is_absolute():
        pressure_path = HERE / pressure_path
    identity = core._pressure_identity(HERE, pressure_path)
    ctx = MPIntervalContext()
    ctx.dps = core.PRECISION
    with mp.workdps(ctx.dps + 40):
        pressure, accepted_sha, pressure_sha, nested = core._load_pressure_axis(ctx, pressure_path)
        hashes = source_hashes(pressure_path)
        if STATE_PATH.exists():
            state = json.loads(STATE_PATH.read_text(encoding="utf-8"))
            validate(state, hashes, pressure_sha, accepted_sha, identity, nested)
            runtime = core._runtime_from_state(ctx, state)
        else:
            fixed = build_fixed(ctx, pressure)
            rows = core._initial_rows(ctx, fixed)
            checks = core.run_self_tests()
            state = core._new_state(fixed, rows, hashes, pressure_sha, accepted_sha,
                                    checks, identity, nested)
            state["target"] = target(identity)
            runtime = dict(fixed=fixed, rows=rows, completed=0)
            core._atomic_write_json(STATE_PATH, state)
        print("Lambda120 starting at radial order", runtime["completed"], flush=True)
        start = time.monotonic()
        while runtime["completed"] < core.RADIAL_DEGREE and time.monotonic()-start < seconds:
            n = runtime["completed"]
            step_start = time.monotonic()
            core._advance_one(ctx, runtime["fixed"], runtime["rows"], n)
            runtime["completed"] += 1
            elapsed = time.monotonic()-step_start
            state["step_timings_seconds"].append(elapsed)
            state["elapsed_seconds"] += elapsed
            core._state_from_runtime(state, runtime)
            core._atomic_write_json(STATE_PATH, state)
            core._write_receipt(RECEIPT_PATH, STATE_PATH, state, state["self_tests"])
            print("Lambda120 radial order", runtime["completed"], "of", core.RADIAL_DEGREE,
                  "step seconds", round(elapsed, 3), flush=True)
        core._state_from_runtime(state, runtime)
        core._atomic_write_json(STATE_PATH, state)
        core._write_receipt(RECEIPT_PATH, STATE_PATH, state, state["self_tests"])
        return state


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--seconds", type=float, default=180.0)
    parser.add_argument("--pressure-file", default=PRESSURE)
    args = parser.parse_args()
    run(args.seconds, args.pressure_file)
