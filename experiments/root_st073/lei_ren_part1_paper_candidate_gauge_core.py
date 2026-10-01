"""Resumable candidate gauge-core radial recurrence at ``Z = 0.3``.

The fixed axis jets come from the accepted pressure-axis receipt and from the
directed ``uniform_axis_jets`` construction with the candidate ``Lambda``.
Each radial order is advanced with the exact algebra of
``factored_core_coefficients`` and saved atomically, so a later invocation can
resume without recomputing completed rows.  The state file stores every
interval endpoint as an exact MP tuple.

This is a centered finite radial jet at one axial point.  It has no infinite
radial remainder, whole-axis, temporal-recursion, or matching claim.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
import os
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
from lei_ren_part1_paper_core_recursion import core_coefficients  # noqa: E402
from lei_ren_part1_paper_factored_core_positivity import (  # noqa: E402
    squared_axis_rows,
)
from lei_ren_part1_paper_schedule_endpoint_enclosures import (  # noqa: E402
    endpoints,
)
from lei_ren_part1_paper_uniform_axis_jets import (  # noqa: E402
    uniform_axis_jets,
)


PRECISION = 260
RADIAL_DEGREE = 124
AXIS_LENGTH = 128
LAMBDA = "1e48"
J = "1e-14"
DELTA = "1e-200"
LOG_C = "5e151"
Z_CENTER = ".3"
PRESSURE_RECEIPT_NAME = "lei_ren_part1_paper_candidate_pressure_axis_jets.json"
STATE_NAME = Path(__file__).stem + "_state.json"
ACCEPTED_SCHEDULE_SHA256 = (
    "736bbadbde99bc2f3d098d279d61ef4cb64418368263a4aba7b275e7f8892de4"
)

SOURCE_NAMES = (
    "lei_ren_part1_paper_amplitude_factored_core.py",
    "lei_ren_part1_paper_core_recursion.py",
    "lei_ren_part1_paper_factored_core_positivity.py",
    "lei_ren_part1_paper_schedule_endpoint_enclosures.py",
    "lei_ren_part1_paper_uniform_axis_jets.py",
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _target(pressure_identity: str = PRESSURE_RECEIPT_NAME) -> dict[str, Any]:
    return {
        "Lambda": LAMBDA,
        "j": J,
        "delta": DELTA,
        "logC": LOG_C,
        "Z": Z_CENTER,
        "precision": PRECISION,
        "radial_degree": RADIAL_DEGREE,
        "axis_length": AXIS_LENGTH,
        "pressure_file": pressure_identity,
    }


def _source_hashes(base: Path, pressure_path: Path) -> dict[str, str]:
    paths = {name: base / name for name in SOURCE_NAMES}
    paths["pressure_input"] = pressure_path
    paths["driver"] = Path(__file__)
    return {name: _sha256(path) for name, path in paths.items()}


def _pressure_identity(base: Path, pressure_path: Path) -> str:
    try:
        return str(pressure_path.resolve().relative_to(base.resolve()))
    except ValueError:
        return str(pressure_path.resolve())


def _artifact_paths(pressure_path: Path) -> tuple[Path, Path]:
    """Keep alternate pressure inputs in separate state and receipt files."""

    if pressure_path.name == PRESSURE_RECEIPT_NAME:
        return base_state_path(), Path(__file__).with_suffix(".json")
    stem = Path(__file__).stem + "_" + pressure_path.stem
    return HERE / (stem + "_state.json"), HERE / (stem + ".json")


def base_state_path() -> Path:
    return HERE / STATE_NAME


def _pack_interval(value: Any) -> dict[str, list[int]]:
    if hasattr(value, "_mpi_"):
        lo, hi = endpoints(value)
    else:
        scalar = mp.mpf(value)
        lo = scalar
        hi = scalar
    return {
        "lower_exact_mpf_tuple": list(lo._mpf_),
        "upper_exact_mpf_tuple": list(hi._mpf_),
    }


def _unpack_interval(ctx: MPIntervalContext, value: dict[str, Any]):
    return ctx.mpf(
        [
            mp.make_mpf(tuple(value["lower_exact_mpf_tuple"])),
            mp.make_mpf(tuple(value["upper_exact_mpf_tuple"])),
        ]
    )


def _pack_vector(values: list[Any]) -> list[dict[str, list[int]]]:
    return [_pack_interval(value) for value in values]


def _unpack_vector(ctx: MPIntervalContext, values: list[Any]) -> list[Any]:
    return [_unpack_interval(ctx, value) for value in values]


def _pack_rows(rows: list[list[Any]]) -> list[list[dict[str, list[int]]]]:
    return [_pack_vector(row) for row in rows]


def _unpack_rows(ctx: MPIntervalContext, rows: list[Any]) -> list[list[Any]]:
    return [_unpack_vector(ctx, row) for row in rows]


def _pack_fixed(fixed: dict[str, Any]) -> dict[str, Any]:
    return {
        "ell_Z_taylor": _pack_vector(fixed["ell_Z_taylor"]),
        "S_Z_taylor": _pack_vector(fixed["S_Z_taylor"]),
        "U0_Z_taylor": _pack_vector(fixed["U0_Z_taylor"]),
        "P0_Z_taylor": _pack_vector(fixed["P0_Z_taylor"]),
        "gradient_coefficients": _pack_vector(fixed["gradient_coefficients"]),
        "F0_interval": _pack_interval(fixed["F0_interval"]),
    }


def _unpack_fixed(ctx: MPIntervalContext, fixed: dict[str, Any]) -> dict[str, Any]:
    return {
        "ell_Z_taylor": _unpack_vector(ctx, fixed["ell_Z_taylor"]),
        "S_Z_taylor": _unpack_vector(ctx, fixed["S_Z_taylor"]),
        "U0_Z_taylor": _unpack_vector(ctx, fixed["U0_Z_taylor"]),
        "P0_Z_taylor": _unpack_vector(ctx, fixed["P0_Z_taylor"]),
        "gradient_coefficients": _unpack_vector(
            ctx, fixed["gradient_coefficients"]
        ),
        "F0_interval": _unpack_interval(ctx, fixed["F0_interval"]),
    }


def _atomic_write_json(path: Path, value: dict[str, Any]) -> None:
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(
        json.dumps(value, indent=2, separators=(",", ": ")) + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, path)


def _interval_midpoint(value: Any) -> mp.mpf:
    lo, hi = endpoints(value)
    return (lo + hi) / 2


def _constant_vector(ctx: MPIntervalContext, value: Any, length: int):
    return [ctx.mpf(value)] + [ctx.mpf(0)] * (length - 1)


def _load_pressure_axis(ctx: MPIntervalContext, path: Path):
    receipt = json.loads(path.read_text(encoding="utf-8"))
    if int(receipt.get("Taylor_length", 0)) < AXIS_LENGTH:
        raise ValueError("pressure-axis receipt is shorter than the fixed axis length")
    if int(receipt.get("stage_count_total", 0)) != 14:
        raise ValueError("pressure-axis receipt does not contain all 14 stages")
    if receipt.get("all_14_true_pressure_stages_included") is not True:
        raise ValueError("pressure-axis receipt does not mark all 14 stages included")
    if receipt.get("pressure_units") != "physical_P":
        raise ValueError("pressure-axis receipt is not in physical_P units")
    if receipt.get("finite_pressure_quadrature_approximation_used"):
        raise ValueError("pressure-axis receipt used finite pressure quadrature")
    if receipt.get("pressure_parameter_order_truncated"):
        raise ValueError("pressure-axis receipt is pressure-parameter truncated")
    accepted_sha = receipt.get("accepted_schedule_sha256")
    if accepted_sha != ACCEPTED_SCHEDULE_SHA256:
        raise ValueError("pressure-axis receipt has an unexpected schedule hash")
    center = receipt.get("center_Z")
    if not isinstance(center, dict):
        raise ValueError("pressure-axis receipt lacks an exact center interval")
    center_interval = _unpack_interval(ctx, center)
    center_lo, center_hi = endpoints(center_interval)
    z = ctx.mpf(Z_CENTER)
    if not center_lo <= z <= center_hi:
        raise ValueError("pressure-axis center interval does not enclose Z = .3")
    rows = receipt["Taylor_rows"][:AXIS_LENGTH]
    if any(int(row.get("order", -1)) != index for index, row in enumerate(rows)):
        raise ValueError("pressure-axis Taylor rows are not ordered")
    nested_hashes = receipt.get("input_hashes")
    if not isinstance(nested_hashes, dict) or not nested_hashes:
        raise ValueError("pressure-axis receipt lacks nested input hashes")
    resolved_nested_hashes = {}
    for name, expected in nested_hashes.items():
        dependency = Path(name)
        if not dependency.is_absolute():
            dependency = path.parent / dependency
        dependency = dependency.resolve()
        if not dependency.exists():
            raise FileNotFoundError(dependency)
        actual = _sha256(dependency)
        if actual != expected:
            raise ValueError(
                f"nested pressure input hash mismatch for {name}: "
                f"expected {expected}, got {actual}"
            )
        resolved_nested_hashes[name] = actual
    pressure = [
        _unpack_interval(ctx, row["physical_pressure_coefficient"])
        for row in rows
    ]
    return pressure, accepted_sha, _sha256(path), resolved_nested_hashes


def _build_fixed_jets(ctx: MPIntervalContext, pressure: list[Any]) -> dict[str, Any]:
    z = ctx.mpf(Z_CENTER)
    lam = ctx.mpf(LAMBDA)
    axis = uniform_axis_jets(
        ctx,
        radius=1,
        j=J,
        Lambda=LAMBDA,
        logC=LOG_C,
        delta=DELTA,
        length=AXIS_LENGTH,
        axial_interval=z,
    )
    gradient = list(axis["gradient_coefficients"])
    return {
        "ell_Z_taylor": [-lam * value for value in gradient],
        "S_Z_taylor": squared_axis_rows(
            ctx,
            gradient,
            lam,
            axis["F0_interval"],
            AXIS_LENGTH,
        ),
        "U0_Z_taylor": list(axis["U0"]),
        "P0_Z_taylor": list(pressure),
        "gradient_coefficients": gradient,
        "F0_interval": axis["F0_interval"],
    }


def _advance_one(
    ctx: MPIntervalContext,
    fixed: dict[str, Any],
    rows: dict[str, list[list[Any]]],
    n: int,
    *,
    delta_value: Any = DELTA,
) -> None:
    """Append radial order ``n + 1`` using the source recurrence verbatim."""

    # The production target has 128 axis coefficients; the independent
    # degree-4 fixture deliberately uses its shorter seven-coefficient jet.
    initial_count = len(fixed["ell_Z_taylor"])
    K = initial_count - n - 1
    if K < 1:
        raise ValueError("no remaining axial coefficients for this radial order")
    convert = ctx.mpf
    zero = convert(0)

    def pad(values):
        values = list(values)
        if len(values) >= K:
            return values[:K]
        return values + [zero] * (K - len(values))

    def const(value):
        return pad([value])

    def add(*args):
        return [sum((arg[k] for arg in args), zero) for k in range(K)]

    def scale(values, scalar):
        return [value * scalar for value in values]

    def mul(left, right):
        return [
            sum((left[i] * right[k - i] for i in range(k + 1)), zero)
            for k in range(K)
        ]

    def diff(values):
        return [
            (k + 1) * values[k + 1] if k + 1 < len(values) else zero
            for k in range(K)
        ]

    def inverse(values):
        quotient = [convert(1) / values[0]]
        for k in range(1, K):
            quotient.append(
                -sum(
                    (values[i] * quotient[k - i] for i in range(1, k + 1)),
                    zero,
                )
                / values[0]
            )
        return quotient

    z = pad([convert(Z_CENTER), convert(1)])
    delta = convert(delta_value)
    one = const(1)
    d = add(one, scale(mul(z, z), -1))
    L = add(one, scale(mul(z, z), -delta))
    invL = inverse(L)
    ell = pad(fixed["ell_Z_taylor"])
    S = pad(fixed["S_Z_taylor"])
    a = rows["A"]
    u = rows["Uz"]
    pressure = rows["P"]

    W = []
    H = []
    for i in range(n + 1):
        wi = scale(
            add(
                scale(mul(z, u[i]), 1 - delta),
                mul(d, diff(u[i])),
            ),
            -convert(1) / (i + 1),
        )
        if i == 0:
            wi = add(one, wi)
        W.append(wi)
        hi = mul(d, u[i])
        if i == 0:
            hi = add(hi, scale(z, (1 - delta) / 2))
        H.append(hi)

    aa_u = const(0)
    uu = const(0)
    rhs_a = const(0)
    rhs_u = const(0)
    for i in range(n + 1):
        aa_u = add(aa_u, mul(u[i], a[n - i]))
        uu = add(uu, mul(u[i], u[n - i]))
        a_derivative = add(diff(a[n - i]), mul(ell, a[n - i]))
        rhs_a = add(
            rhs_a,
            scale(mul(W[i], a[n - i]), n - i + 1),
            mul(H[i], a_derivative),
        )
        rhs_u = add(
            rhs_u,
            scale(mul(W[i], u[n - i]), n - i),
            mul(H[i], diff(u[n - i])),
        )

    rhs_a = add(
        rhs_a,
        scale(add(a[n], scale(mul(z, aa_u), -2)), delta / 2),
    )
    rhs_u = add(
        rhs_u,
        scale(add(u[n], scale(mul(z, uu), -2)), (1 + delta) / 2),
        mul(d, diff(pressure[n])),
        scale(mul(z, pressure[n]), -2 * (1 + delta)),
    )

    if n:
        aa_previous = const(0)
        for i in range(n):
            aa_previous = add(
                aa_previous,
                mul(a[i], a[n - 1 - i]),
            )
        rhs_u = add(
            rhs_u,
            scale(mul(z, mul(S, aa_previous)), -2),
        )

    a_next = scale(
        mul(invL, rhs_a),
        convert(1) / (2 * (n + 1) * (n + 2)),
    )
    u_next = scale(
        mul(invL, rhs_u),
        convert(1) / (2 * (n + 1) ** 2),
    )
    aa_current = const(0)
    for i in range(n + 1):
        aa_current = add(aa_current, mul(a[i], a[n - i]))
    p_next = scale(mul(S, aa_current), convert(1) / (n + 1))

    rows["A"].append(a_next)
    rows["Uz"].append(u_next)
    rows["P"].append(p_next)


def _initial_rows(ctx: MPIntervalContext, fixed: dict[str, Any]):
    return {
        "A": [_constant_vector(ctx, 1, AXIS_LENGTH)],
        "Uz": [list(fixed["U0_Z_taylor"])],
        "P": [list(fixed["P0_Z_taylor"])],
    }


def _new_state(
    fixed: dict[str, Any],
    rows: dict[str, list[list[Any]]],
    source_hashes: dict[str, str],
    pressure_sha: str,
    accepted_sha: str,
    self_tests: dict[str, Any],
    pressure_identity: str,
    nested_hashes: dict[str, str],
) -> dict[str, Any]:
    return {
        "schema_version": 1,
        "target": _target(pressure_identity),
        "input_pressure_receipt": pressure_identity,
        "input_pressure_receipt_sha256": pressure_sha,
        "nested_pressure_input_hashes": nested_hashes,
        "source_hashes": source_hashes,
        "accepted_schedule_sha256": accepted_sha,
        "fixed_jets": _pack_fixed(fixed),
        "initial_count": AXIS_LENGTH,
        "completed_radial_order": 0,
        "next_radial_order": 0,
        "A_rows": _pack_rows(rows["A"]),
        "Uz_rows": _pack_rows(rows["Uz"]),
        "P_rows": _pack_rows(rows["P"]),
        "step_timings_seconds": [],
        "elapsed_seconds": 0.0,
        "self_tests": self_tests,
        "candidate_center_only": True,
        "finite_radial_jet_only": True,
        "infinite_radial_remainder_enclosed": False,
        "whole_axis_certified": False,
        "global_matching_claim": False,
        "temporal_recursion": False,
    }


def _validate_state(
    state: dict[str, Any],
    source_hashes: dict[str, str],
    pressure_sha: str,
    accepted_sha: str,
    pressure_identity: str,
    nested_hashes: dict[str, str],
) -> None:
    if state.get("target") != _target(pressure_identity):
        raise ValueError("candidate target/settings changed; refusing resume")
    if state.get("input_pressure_receipt_sha256") != pressure_sha:
        raise ValueError("pressure-axis receipt changed; refusing resume")
    if state.get("nested_pressure_input_hashes") != nested_hashes:
        raise ValueError("nested pressure input hashes changed; refusing resume")
    if state.get("source_hashes") != source_hashes:
        raise ValueError("input/source hashes changed; refusing resume")
    if state.get("accepted_schedule_sha256") != accepted_sha:
        raise ValueError("accepted schedule changed; refusing resume")
    if state.get("initial_count") != AXIS_LENGTH:
        raise ValueError("axis length changed; refusing resume")
    completed = int(state.get("completed_radial_order", -1))
    if state.get("next_radial_order") != completed:
        raise ValueError("inconsistent completed radial order")
    for key in ("A_rows", "Uz_rows", "P_rows"):
        if len(state.get(key, [])) != completed + 1:
            raise ValueError(f"inconsistent {key} row count")
    if completed < 0 or completed > RADIAL_DEGREE:
        raise ValueError("completed radial order outside target")


def _runtime_from_state(ctx: MPIntervalContext, state: dict[str, Any]):
    return {
        "fixed": _unpack_fixed(ctx, state["fixed_jets"]),
        "rows": {
            "A": _unpack_rows(ctx, state["A_rows"]),
            "Uz": _unpack_rows(ctx, state["Uz_rows"]),
            "P": _unpack_rows(ctx, state["P_rows"]),
        },
        "completed": int(state["completed_radial_order"]),
    }


def _state_from_runtime(state: dict[str, Any], runtime: dict[str, Any]) -> None:
    state["A_rows"] = _pack_rows(runtime["rows"]["A"])
    state["Uz_rows"] = _pack_rows(runtime["rows"]["Uz"])
    state["P_rows"] = _pack_rows(runtime["rows"]["P"])
    state["completed_radial_order"] = runtime["completed"]
    state["next_radial_order"] = runtime["completed"]


def _write_receipt(
    receipt_path: Path,
    state_path: Path,
    state: dict[str, Any],
    self_tests: dict[str, Any],
) -> None:
    completed = int(state["completed_radial_order"])
    payload = {
        "schema_version": 1,
        "target": state["target"],
        "candidate_center_Z": Z_CENTER,
        "state_file": state_path.name,
        "state_sha256": _sha256(state_path),
        "input_pressure_receipt": state["input_pressure_receipt"],
        "input_pressure_receipt_sha256": state[
            "input_pressure_receipt_sha256"
        ],
        "nested_pressure_input_hashes": state[
            "nested_pressure_input_hashes"
        ],
        "source_hashes": state["source_hashes"],
        "accepted_schedule_sha256": state["accepted_schedule_sha256"],
        "completed_radial_order": completed,
        "target_radial_degree": RADIAL_DEGREE,
        "next_radial_order": completed,
        "step_count": len(state["step_timings_seconds"]),
        "step_timings_seconds": state["step_timings_seconds"],
        "elapsed_seconds": state["elapsed_seconds"],
        "exact_interval_storage": {
            "file": state_path.name,
            "fixed_jets": [
                "fixed_jets.ell_Z_taylor",
                "fixed_jets.S_Z_taylor",
                "fixed_jets.U0_Z_taylor",
                "fixed_jets.P0_Z_taylor",
                "fixed_jets.gradient_coefficients",
                "fixed_jets.F0_interval",
            ],
            "rows": ["A_rows", "Uz_rows", "P_rows"],
            "format": "each interval stores lower_exact_mpf_tuple and upper_exact_mpf_tuple",
        },
        "self_tests": self_tests,
        "candidate_center_only": True,
        "finite_radial_jet_only": True,
        "infinite_radial_remainder_enclosed": False,
        "whole_axis_certified": False,
        "global_matching_claim": False,
        "temporal_recursion": False,
    }
    _atomic_write_json(receipt_path, payload)


def _taylor(function, center: mp.mpf, length: int):
    return [
        mp.diff(function, center, order) / mp.factorial(order)
        for order in range(length)
    ]


def _series_product(left: list[Any], right: list[Any], length: int):
    return [
        sum((left[i] * right[k - i] for i in range(k + 1)), mp.mpf(0))
        for k in range(length)
    ]


def _scaled_error(value: mp.mpf, reference: mp.mpf) -> mp.mpf:
    return abs(value - reference) / max(
        mp.mpf(1), abs(value), abs(reference)
    )


def _interval_contains(value: Any, reference: mp.mpf, tolerance: mp.mpf) -> bool:
    lo, hi = endpoints(value)
    return lo - tolerance <= reference <= hi + tolerance


def _run_small_recurrence(
    ctx: MPIntervalContext,
    fixed: dict[str, Any],
    degree: int,
    *,
    delta_value: Any,
):
    rows = _initial_rows(ctx, fixed)
    for n in range(degree):
        _advance_one(ctx, fixed, rows, n, delta_value=delta_value)
    return rows


def run_self_tests() -> dict[str, Any]:
    """Check degree 4 against the original full recurrence and JSON resume."""

    degree = 4
    depth = 3
    length = degree + depth + 1
    with mp.workdps(240):
        ctx = MPIntervalContext()
        ctx.dps = 220
        center = mp.mpf(".3")
        f0 = _taylor(lambda x: mp.exp(-2 * x * x), center, length)
        ell = _taylor(lambda x: -4 * x, center, length)
        squared = _taylor(lambda x: mp.exp(-4 * x * x), center, length)
        u0 = _taylor(lambda x: 4 * x + mp.mpf(".01"), center, length)
        p0 = _taylor(
            lambda x: -3 / (1 + x * x) ** 2,
            center,
            length,
        )
        fixed = {
            "ell_Z_taylor": [ctx.mpf(value) for value in ell],
            "S_Z_taylor": [ctx.mpf(value) for value in squared],
            "U0_Z_taylor": [ctx.mpf(value) for value in u0],
            "P0_Z_taylor": [ctx.mpf(value) for value in p0],
            "gradient_coefficients": [ctx.mpf(0)] * length,
            "F0_interval": ctx.mpf(f0[0]),
        }
        full_rows = _run_small_recurrence(
            ctx, fixed, degree, delta_value=".01"
        )
        original = core_coefficients(
            center,
            ".01",
            F0_Z_taylor=f0,
            U0_Z_taylor=u0,
            P0_Z_taylor=p0,
            radial_degree=degree,
            precision=200,
        )
        factored = factored_core_coefficients(
            center,
            ".01",
            ell_Z_taylor=ell,
            S_Z_taylor=squared,
            U0_Z_taylor=u0,
            P0_Z_taylor=p0,
            radial_degree=degree,
            precision=200,
        )
        tolerance = mp.mpf("1e-165")
        maximum_reconstruction_error = mp.mpf(0)
        maximum_direct_error = mp.mpf(0)
        containment_checks = 0
        for n in range(degree + 1):
            reconstructed = _series_product(
                f0,
                [
                    _interval_midpoint(value)
                    for value in full_rows["A"][n]
                ],
                len(original["F"][n]),
            )
            for value, reference in zip(reconstructed, original["F"][n]):
                maximum_reconstruction_error = max(
                    maximum_reconstruction_error,
                    _scaled_error(value, reference),
                )
            for name in ("Uz", "P"):
                for value, reference in zip(
                    full_rows[name][n], original[name][n]
                ):
                    midpoint = _interval_midpoint(value)
                    maximum_direct_error = max(
                        maximum_direct_error,
                        _scaled_error(midpoint, reference),
                    )
                    if not _interval_contains(value, reference, tolerance):
                        raise AssertionError(
                            f"degree-4 {name} row {n} misses original coefficient"
                        )
                    containment_checks += 1
            for value, reference in zip(full_rows["A"][n], factored["A"][n]):
                if not _interval_contains(value, reference, tolerance):
                    raise AssertionError(
                        f"degree-4 A row {n} misses factored coefficient"
                    )
                containment_checks += 1

        split_rows = _run_small_recurrence(
            ctx, fixed, 2, delta_value=".01"
        )
        serialized = {
            "A_rows": _pack_rows(split_rows["A"]),
            "Uz_rows": _pack_rows(split_rows["Uz"]),
            "P_rows": _pack_rows(split_rows["P"]),
        }
        resumed_rows = {
            "A": _unpack_rows(ctx, serialized["A_rows"]),
            "Uz": _unpack_rows(ctx, serialized["Uz_rows"]),
            "P": _unpack_rows(ctx, serialized["P_rows"]),
        }
        for n in range(2, degree):
            _advance_one(ctx, fixed, resumed_rows, n, delta_value=".01")
        resume_max_error = mp.mpf(0)
        resume_checks = 0
        for name in ("A", "Uz", "P"):
            for left_row, right_row in zip(
                full_rows[name], resumed_rows[name]
            ):
                for left, right in zip(left_row, right_row):
                    resume_max_error = max(
                        resume_max_error,
                        _scaled_error(
                            _interval_midpoint(left),
                            _interval_midpoint(right),
                        ),
                    )
                    resume_checks += 1

        if maximum_reconstruction_error >= tolerance:
            raise AssertionError(
                f"degree-4 F reconstruction error {maximum_reconstruction_error}"
            )
        if maximum_direct_error >= tolerance:
            raise AssertionError(
                f"degree-4 direct recurrence error {maximum_direct_error}"
            )
        if resume_max_error != 0:
            raise AssertionError(f"resume mismatch {resume_max_error}")
        return {
            "all_checks_passed": True,
            "radial_degree": degree,
            "Z_derivative_depth": depth,
            "jet_length": length,
            "maximum_reconstruction_scaled_error": mp.nstr(
                maximum_reconstruction_error, 40
            ),
            "maximum_direct_scaled_error": mp.nstr(
                maximum_direct_error, 40
            ),
            "resume_maximum_scaled_error": mp.nstr(resume_max_error, 40),
            "containment_checks": containment_checks,
            "resume_checks": resume_checks,
            "comparison": "Reconstructed F = F0 * A and direct Uz/P rows against original core_coefficients; A also checked against factored_core_coefficients.",
            "scope": "Resolved degree-4 algebra and JSON row serialization only; no candidate source or global claim.",
        }


def run(
    seconds: float = 180.0,
    pressure_file: str = PRESSURE_RECEIPT_NAME,
) -> dict[str, Any]:
    if seconds < 0:
        raise ValueError("seconds must be nonnegative")
    base = HERE
    pressure_path = Path(pressure_file)
    if not pressure_path.is_absolute():
        pressure_path = base / pressure_path
    pressure_path = pressure_path.resolve()
    if not pressure_path.exists():
        raise FileNotFoundError(pressure_path)
    pressure_identity = _pressure_identity(base, pressure_path)
    state_path, receipt_path = _artifact_paths(pressure_path)

    self_tests = run_self_tests()
    ctx = MPIntervalContext()
    ctx.dps = PRECISION
    with mp.workdps(PRECISION + 40):
        pressure, accepted_sha, pressure_sha, nested_hashes = _load_pressure_axis(
            ctx, pressure_path
        )
        source_hashes = _source_hashes(base, pressure_path)
        if state_path.exists():
            state = json.loads(state_path.read_text(encoding="utf-8"))
            _validate_state(
                state,
                source_hashes,
                pressure_sha,
                accepted_sha,
                pressure_identity,
                nested_hashes,
            )
            runtime = _runtime_from_state(ctx, state)
            fixed = runtime["fixed"]
            print(
                "resuming candidate gauge core at radial order",
                runtime["completed"],
                flush=True,
            )
        else:
            fixed = _build_fixed_jets(ctx, pressure)
            rows = _initial_rows(ctx, fixed)
            runtime = {"fixed": fixed, "rows": rows, "completed": 0}
            state = _new_state(
                fixed,
                rows,
                source_hashes,
                pressure_sha,
                accepted_sha,
                self_tests,
                pressure_identity,
                nested_hashes,
            )
            _atomic_write_json(state_path, state)
            _write_receipt(receipt_path, state_path, state, self_tests)
            print("initialized candidate gauge core state", flush=True)

        if state.get("self_tests") != self_tests:
            state["self_tests"] = self_tests

        start = time.monotonic()
        while (
            runtime["completed"] < RADIAL_DEGREE
            and time.monotonic() - start < seconds
        ):
            n = runtime["completed"]
            step_start = time.monotonic()
            _advance_one(ctx, fixed, runtime["rows"], n)
            runtime["completed"] += 1
            step_seconds = time.monotonic() - step_start
            state["step_timings_seconds"].append(step_seconds)
            state["elapsed_seconds"] += step_seconds
            _state_from_runtime(state, runtime)
            _atomic_write_json(state_path, state)
            _write_receipt(receipt_path, state_path, state, self_tests)
            print(
                "completed radial order",
                runtime["completed"],
                "of",
                RADIAL_DEGREE,
                "step_seconds",
                round(step_seconds, 3),
                "elapsed_seconds",
                round(state["elapsed_seconds"], 3),
                flush=True,
            )

        _state_from_runtime(state, runtime)
        _atomic_write_json(state_path, state)
        _write_receipt(receipt_path, state_path, state, self_tests)
        print(
            "saved candidate gauge core radial order",
            runtime["completed"],
            "of",
            RADIAL_DEGREE,
            "pending",
            RADIAL_DEGREE - runtime["completed"],
            flush=True,
        )
        return state


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--seconds", type=float, default=180.0)
    parser.add_argument(
        "--pressure-file",
        default=PRESSURE_RECEIPT_NAME,
        help="accepted pressure-axis receipt; alternate files get separate artifacts",
    )
    parser.add_argument("--self-test-only", action="store_true")
    args = parser.parse_args()
    if args.self_test_only:
        print(json.dumps(run_self_tests(), indent=2, default=str))
    else:
        run(args.seconds, pressure_file=args.pressure_file)
