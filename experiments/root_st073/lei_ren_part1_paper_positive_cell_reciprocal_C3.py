"""Resumable reciprocal mixed C3 bounds on the accepted positive cells.

The input cells are an immutable snapshot of the accepted finite-core
positivity receipt.  This driver does not retest positivity.  On each cell it
replays the exact amplitude-factored degree-18 gauge recurrence, converts the
radial polynomial to Bernstein form on ``x=R/R_a`` in ``[0,1]``, and computes
the reciprocal multivariate Taylor coefficients through total order three.

Only the accepted snapshot domain is processed.  The finite polynomial,
whole-core, and infinite radial remainder limitations remain explicit in the
receipt.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
import time
from typing import Any

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext


HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from lei_ren_part1_paper_amplitude_factored_core import (  # noqa: E402
    factored_core_coefficients,
)
from lei_ren_part1_paper_factored_core_positivity import (  # noqa: E402
    bernstein_enclosure,
    squared_axis_rows,
)
from lei_ren_part1_paper_global_finite_core_bounds import (  # noqa: E402
    accepted_pressure_axis_rows,
)
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints  # noqa: E402
from lei_ren_part1_paper_uniform_axis_jets import uniform_axis_jets  # noqa: E402


PRECISION = 260
RADIAL_DEGREE = 18
JET_LENGTH = RADIAL_DEGREE + 4
LAMBDA = "1e36"
J = "1e-14"
LOG_C = "5e151"
DELTA = "1e-200"
SNAPSHOT_NAME = "lei_ren_part1_paper_reciprocal_input_snapshot.json"


def pack(value: Any) -> list[int]:
    """Store an exact MP scalar tuple in JSON."""

    return list(value._mpf_)


def unpack(value: list[int]) -> mp.mpf:
    return mp.make_mpf(tuple(value))


def interval_pack(value: Any) -> dict[str, list[int]]:
    lo, hi = endpoints(value)
    return {"lower": pack(lo), "upper": pack(hi)}


def scalar_upper(value: Any) -> mp.mpf:
    lo, hi = endpoints(value)
    return max(abs(lo), abs(hi))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dependency_hashes(snapshot: dict[str, Any], base: Path) -> dict[str, str]:
    expected = dict(snapshot["input_hashes"])
    actual = {}
    for name, expected_hash in expected.items():
        path = base / name
        if not path.exists():
            raise FileNotFoundError(path)
        actual[name] = sha256(path)
        if actual[name] != expected_hash:
            raise ValueError(
                f"Dependency hash mismatch for {name}: "
                f"expected {expected_hash}, got {actual[name]}"
            )
    return actual


def multi_indices(max_total: int = 3) -> list[tuple[int, int]]:
    return [
        (i, k)
        for total in range(max_total + 1)
        for i in range(total + 1)
        for k in (total - i,)
    ]


def _axis_derivative_ranges(
    ctx: MPIntervalContext,
    core: dict[str, Any],
    radius: Any,
    *,
    saved_lower: Any,
) -> dict[tuple[int, int], Any]:
    """Bernstein-enclose normalized A Taylor coefficients through C3."""

    ranges: dict[tuple[int, int], Any] = {}
    for radial_order, axial_order in multi_indices(3):
        coefficients = []
        for m in range(RADIAL_DEGREE - radial_order + 1):
            source_row = core["A"][m + radial_order]
            coefficient = source_row[axial_order]
            coefficient *= math.comb(m + radial_order, radial_order)
            coefficient *= radius ** (m + radial_order)
            coefficients.append(coefficient)
        ranges[(radial_order, axial_order)] = bernstein_enclosure(
            ctx, coefficients
        )["range"]

    # The stored lower endpoint is the already accepted positivity proof.
    # The upper endpoint is freshly obtained from the current recurrence;
    # do not reuse the old upper endpoint for the reciprocal calculation.
    new_a00 = ranges[(0, 0)]
    _, new_upper = endpoints(new_a00)
    lower = unpack(saved_lower)
    if lower <= 0 or new_upper < lower:
        raise AssertionError("Saved positive lower is outside new A00 upper")
    ranges[(0, 0)] = ctx.mpf([lower, new_upper])
    return ranges


def _reciprocal_coefficients(
    ctx: MPIntervalContext,
    a_ranges: dict[tuple[int, int], Any],
) -> dict[tuple[int, int], Any]:
    """Return B coefficients for the reciprocal B=1/A through total C3."""

    b_ranges: dict[tuple[int, int], Any] = {}
    a00 = a_ranges[(0, 0)]
    b_ranges[(0, 0)] = 1 / a00
    for alpha in multi_indices(3)[1:]:
        i, k = alpha
        convolution = ctx.mpf(0)
        for beta_i in range(i + 1):
            for beta_k in range(k + 1):
                if beta_i == 0 and beta_k == 0:
                    continue
                convolution += a_ranges[(beta_i, beta_k)] * b_ranges[
                    (i - beta_i, k - beta_k)
                ]
        b_ranges[alpha] = -convolution / a00
    return b_ranges


def _inverse_axis_coefficients(
    ctx: MPIntervalContext, gradient: list[Any], Lambda: Any
) -> list[Any]:
    """Taylor coefficients for the normalized reciprocal axis factor."""

    lam = ctx.mpf(Lambda)
    inverse = [ctx.mpf(1)]
    for n in range(3):
        inverse.append(
            lam
            * sum(
                (gradient[k] * inverse[n - k] for k in range(n + 1)),
                ctx.mpf(0),
            )
            / (n + 1)
        )
    return inverse


def _reciprocal_mixed_derivatives(
    ctx: MPIntervalContext,
    b_ranges: dict[tuple[int, int], Any],
    inverse_axis: list[Any],
) -> dict[tuple[int, int], Any]:
    """Convolve B with the inverse-axis factor and restore factorials."""

    result = {}
    for radial_order, axial_order in multi_indices(3):
        coefficient = ctx.mpf(0)
        for axis_order in range(axial_order + 1):
            coefficient += inverse_axis[axis_order] * b_ranges[
                (radial_order, axial_order - axis_order)
            ]
        result[(radial_order, axial_order)] = (
            math.factorial(radial_order)
            * math.factorial(axial_order)
            * coefficient
        )
    return result


def _record_cell(
    ctx: MPIntervalContext,
    cell: dict[str, Any],
    *,
    pressure_rows: list[Any],
    radius: Any,
    degree: int,
    precision: int,
    pressure_schedule_sha256: str,
) -> dict[str, Any]:
    left = unpack(cell["left"])
    right = unpack(cell["right"])
    z = ctx.mpf([left, right])
    axis = uniform_axis_jets(
        ctx,
        radius=1,
        j=J,
        Lambda=LAMBDA,
        logC=LOG_C,
        delta=DELTA,
        length=JET_LENGTH,
        axial_interval=z,
    )
    gradient = axis["gradient_coefficients"]
    core = factored_core_coefficients(
        z,
        DELTA,
        ell_Z_taylor=[-ctx.mpf(LAMBDA) * value for value in gradient],
        S_Z_taylor=squared_axis_rows(
            ctx,
            gradient,
            ctx.mpf(LAMBDA),
            axis["F0_interval"],
            JET_LENGTH,
        ),
        U0_Z_taylor=axis["U0"],
        P0_Z_taylor=pressure_rows,
        radial_degree=degree,
        precision=precision,
        scalar_converter=ctx.mpf,
    )
    a_ranges = _axis_derivative_ranges(
        ctx,
        core,
        radius,
        saved_lower=cell["lower"],
    )
    b_ranges = _reciprocal_coefficients(ctx, a_ranges)
    inverse_axis = _inverse_axis_coefficients(
        ctx, gradient, ctx.mpf(LAMBDA)
    )
    mixed = _reciprocal_mixed_derivatives(ctx, b_ranges, inverse_axis)

    ratio_sum = ctx.mpf(0)
    for value in mixed.values():
        ratio_sum += ctx.mpf([0, scalar_upper(value)])
    inverse_axis_log_prefactor = -ctx.log(axis["F0_interval"])
    upper_log = inverse_axis_log_prefactor + ctx.log(ratio_sum)

    def keyed(values: dict[tuple[int, int], Any]) -> dict[str, Any]:
        return {
            f"x{i}_Z{k}": interval_pack(values[(i, k)])
            for i, k in multi_indices(3)
        }

    return {
        "left": cell["left"],
        "right": cell["right"],
        "depth": cell["depth"],
        "accepted_lower": cell["lower"],
        "A_normalized_mixed": keyed(a_ranges),
        "B_normalized_mixed": keyed(b_ranges),
        "inverse_F_normalized_mixed": keyed(mixed),
        "inverse_axis_log_prefactor": interval_pack(
            inverse_axis_log_prefactor
        ),
        "inverse_F_C3_log_upper": interval_pack(upper_log),
        "upper_log": interval_pack(upper_log),
        "ratio_sum": interval_pack(ratio_sum),
        "accepted_schedule_sha256": pressure_schedule_sha256,
        "finite_radial_polynomial_degree": degree,
        "normalized_radial_domain": ["0", "1"],
    }


def _new_state(
    snapshot: dict[str, Any],
    snapshot_sha256: str,
    dependencies: dict[str, str],
    driver_sha256: str,
) -> dict[str, Any]:
    return {
        "input_snapshot": SNAPSHOT_NAME,
        "input_snapshot_sha256": snapshot_sha256,
        "dependency_hashes": dependencies,
        "driver_sha256": driver_sha256,
        "precision": PRECISION,
        "radial_degree": RADIAL_DEGREE,
        "axis_domain": ["-1", "1"],
        "normalized_radial_domain": ["0", "1"],
        "pending": [dict(cell) for cell in snapshot["accepted"]],
        "completed": [],
        "evaluations": 0,
        "accepted_schedule_sha256": snapshot["accepted_schedule_sha256"],
        "accepted_cell_count": len(snapshot["accepted"]),
        "accepted_axial_length_fraction": snapshot[
            "accepted_axial_length_fraction"
        ],
        "accepted_domain_only": True,
        "whole_axis_reciprocal_C3_certified": False,
        "full_K_certified": False,
        "infinite_radial_remainder_enclosed": False,
        "temporal_recursion": False,
        "recomputation_is_not_positivity_testing": True,
        "aggregate_max_upper_log": None,
    }


def _validate_state(
    state: dict[str, Any],
    snapshot: dict[str, Any],
    snapshot_sha256: str,
    dependencies: dict[str, str],
    driver_sha256: str,
) -> None:
    if state.get("input_snapshot_sha256") != snapshot_sha256:
        raise ValueError("Input snapshot SHA changed; refusing resume")
    if state.get("dependency_hashes") != dependencies:
        raise ValueError("Dependency hashes changed; refusing resume")
    if state.get("driver_sha256") != driver_sha256:
        raise ValueError("Driver hash changed; refusing resume")
    if state.get("precision") != PRECISION or state.get("radial_degree") != RADIAL_DEGREE:
        raise ValueError("Precision or radial degree changed; refusing resume")
    if state.get("accepted_schedule_sha256") != snapshot["accepted_schedule_sha256"]:
        raise ValueError("Accepted pressure schedule changed; refusing resume")


def _save_state(path: Path, state: dict[str, Any]) -> None:
    values = [
        unpack(record["upper_log"]["upper"])
        for record in state["completed"]
    ]
    state["completed_count"] = len(state["completed"])
    state["pending_count"] = len(state["pending"])
    state["all_snapshot_cells_processed"] = not state["pending"]
    state["aggregate_max_upper_log"] = (
        pack(max(values)) if values else None
    )
    path.write_text(json.dumps(state, indent=2) + "\n")


def run(seconds: float = 180.0) -> dict[str, Any]:
    base = HERE
    snapshot_path = base / SNAPSHOT_NAME
    output_path = Path(__file__).with_suffix(".json")
    if not snapshot_path.exists():
        raise FileNotFoundError(snapshot_path)
    snapshot = json.loads(snapshot_path.read_text())
    snapshot_sha256 = sha256(snapshot_path)
    dependencies = dependency_hashes(snapshot, base)
    driver_sha256 = sha256(Path(__file__))
    if snapshot.get("precision") != PRECISION or snapshot.get("radial_degree") != RADIAL_DEGREE:
        raise ValueError("Snapshot precision or radial degree does not match driver")

    if output_path.exists():
        state = json.loads(output_path.read_text())
        _validate_state(
            state, snapshot, snapshot_sha256, dependencies, driver_sha256
        )
    else:
        state = _new_state(
            snapshot, snapshot_sha256, dependencies, driver_sha256
        )

    ctx = MPIntervalContext()
    ctx.dps = PRECISION
    error_path = base / "lei_ren_part1_paper_global_pressure_high_derivatives.json"
    with mp.workdps(PRECISION + 40):
        pressure_rows, alignment = accepted_pressure_axis_rows(
            ctx,
            ctx.mpf([-1, 1]),
            JET_LENGTH,
            error_path,
        )
        pressure_rows = [value.value for value in pressure_rows]
        if alignment["accepted_schedule"]["sha256"] != state[
            "accepted_schedule_sha256"
        ]:
            raise ValueError("Accepted schedule SHA changed")

        start = time.monotonic()
        while state["pending"] and time.monotonic() - start < seconds:
            cell = state["pending"].pop()
            record = _record_cell(
                ctx,
                cell,
                pressure_rows=pressure_rows,
                radius=ctx.mpf(4) / ctx.mpf(LAMBDA),
                degree=RADIAL_DEGREE,
                precision=PRECISION,
                pressure_schedule_sha256=state["accepted_schedule_sha256"],
            )
            state["completed"].append(record)
            state["evaluations"] += 1
            _save_state(output_path, state)
            print(
                "processed",
                state["evaluations"],
                "completed",
                len(state["completed"]),
                "pending",
                len(state["pending"]),
                "aggregate_upper_log",
                mp.nstr(
                    unpack(state["aggregate_max_upper_log"]), 18
                ),
                flush=True,
            )
        _save_state(output_path, state)
        print(
            "saved completed",
            len(state["completed"]),
            "pending",
            len(state["pending"]),
            "all processed",
            state["all_snapshot_cells_processed"],
            flush=True,
        )
    return state


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--seconds", type=float, default=180.0)
    args = parser.parse_args()
    run(args.seconds)
