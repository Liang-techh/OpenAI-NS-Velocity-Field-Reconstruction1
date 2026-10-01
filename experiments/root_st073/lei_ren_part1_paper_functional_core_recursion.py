"""Functional-center directed radial core recurrence.

This module is a small, side-effect-free generalization of the directed
``_advance_one`` recurrence in ``lei_ren_part1_paper_candidate_gauge_core``.
The axial center is an explicit scalar or real interval, and all input jets
are consumed with one axial coefficient removed per radial step.  Thus an
input length of ``degree + required_depth + 1`` leaves ``required_depth + 1``
coefficients in the final row.

The result is an interval enclosure over the supplied center interval and
input jet intervals.  With a non-degenerate center interval it is not a
single Taylor series about one point.  The recurrence retains the physical
``S = F0**2`` swirl and pressure factors and all pressure/cross terms from
the Section 8.2 equations.  It provides no infinite-radial, whole-axis, or
terminal matching claim.
"""

from __future__ import annotations

import hashlib
import json
import operator
from pathlib import Path
import sys
from typing import Any

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext


HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from lei_ren_part1_paper_amplitude_factored_core import (  # noqa: E402
    factored_core_coefficients,
)
from lei_ren_part1_paper_candidate_gauge_core import (  # noqa: E402
    _advance_one as _legacy_advance_one,
    _interval_contains,
    _pack_rows,
    _unpack_rows,
)
from lei_ren_part1_paper_core_recursion import core_coefficients  # noqa: E402
from lei_ren_part1_paper_schedule_endpoint_enclosures import (  # noqa: E402
    endpoints,
)


PRECISION = 220
DEFAULT_REQUIRED_DEPTH = 3
SOURCE_NAMES = (
    "lei_ren_part1_paper_amplitude_factored_core.py",
    "lei_ren_part1_paper_candidate_gauge_core.py",
    "lei_ren_part1_paper_core_recursion.py",
    "lei_ren_part1_paper_schedule_endpoint_enclosures.py",
    "lei_ren_part1_paper_functional_core_recursion.py",
)


def _source_hashes() -> dict[str, str]:
    """Return hashes for the copied recurrence and its algebra sources."""

    result: dict[str, str] = {}
    for name in SOURCE_NAMES:
        path = Path(name) if Path(name).is_absolute() else HERE / name
        result[name] = hashlib.sha256(path.read_bytes()).hexdigest()
    return result


def _as_interval(ctx: MPIntervalContext, value: Any):
    """Convert a scalar or preserve an ``MPIntervalContext`` interval."""

    if hasattr(value, "_mpi_"):
        return value
    return ctx.mpf(value)


def _coerce_vector(ctx: MPIntervalContext, values: Any) -> list[Any]:
    return [_as_interval(ctx, value) for value in values]


def _validate_fixed(fixed: dict[str, Any], degree: int, required_depth: int) -> int:
    names = ("ell_Z_taylor", "S_Z_taylor", "U0_Z_taylor", "P0_Z_taylor")
    missing = [name for name in names if name not in fixed]
    if missing:
        raise KeyError("fixed jets missing: " + ", ".join(missing))
    lengths = [len(fixed[name]) for name in names]
    initial_count = min(lengths)
    needed = int(degree) + int(required_depth) + 1
    if initial_count < needed:
        raise ValueError(
            "axis jets must have degree + required_depth + 1 coefficients "
            f"({needed}); got {initial_count}"
        )
    return initial_count


def coupled_rows(
    ctx: MPIntervalContext,
    fixed: dict[str, Any],
    center: Any,
    degree: int,
    delta: Any,
    *,
    required_depth: int = DEFAULT_REQUIRED_DEPTH,
) -> dict[str, Any]:
    """Compute directed interval radial rows at a real axial center.

    ``fixed`` contains ordinary Taylor coefficients of ``ell = F0'/F0``,
    ``S = F0**2``, ``U0`` and physical ``P0``.  Every coefficient may be a
    scalar or an interval.  For radial degree ``degree`` and retained axial
    depth ``required_depth``, each input jet must contain at least
    ``degree + required_depth + 1`` coefficients.  One axial coefficient is
    consumed at each radial step; no requested final coefficient is silently
    truncated.

    A non-degenerate ``center`` is an interval enclosure over centers.  The
    returned coefficient intervals therefore represent a family of local
    Taylor rows, rather than coefficients of a single Taylor series.
    """

    try:
        degree = operator.index(degree)
        required_depth = operator.index(required_depth)
    except TypeError as exc:
        raise ValueError("degree and required_depth must be integers") from exc
    if degree < 0:
        raise ValueError("degree must be a nonnegative integer")
    if required_depth < 0:
        raise ValueError("required_depth must be nonnegative")
    initial_count = _validate_fixed(fixed, degree, required_depth)

    convert = ctx.mpf
    z_center = _as_interval(ctx, center)
    dt = _as_interval(ctx, delta)
    zero = convert(0)

    def pad(values: Any, length: int):
        values = list(values)
        if len(values) >= length:
            return [_as_interval(ctx, value) for value in values[:length]]
        return [_as_interval(ctx, value) for value in values] + [zero] * (
            length - len(values)
        )

    def const(value: Any, length: int):
        return pad([value], length)

    def add(*args: list[Any]):
        return [sum((arg[k] for arg in args), zero) for k in range(len(args[0]))]

    def scale(values: list[Any], scalar: Any):
        return [value * scalar for value in values]

    def mul(left: list[Any], right: list[Any]):
        length = len(left)
        return [
            sum((left[i] * right[k - i] for i in range(k + 1)), zero)
            for k in range(length)
        ]

    def diff(values: list[Any]):
        length = len(values)
        return [
            (k + 1) * values[k + 1] if k + 1 < length else zero
            for k in range(length)
        ]

    def inverse(values: list[Any]):
        quotient = [convert(1) / values[0]]
        for k in range(1, len(values)):
            quotient.append(
                -sum(
                    (values[i] * quotient[k - i] for i in range(1, k + 1)),
                    zero,
                )
                / values[0]
            )
        return quotient

    ell = _coerce_vector(ctx, fixed["ell_Z_taylor"])[:initial_count]
    squared = _coerce_vector(ctx, fixed["S_Z_taylor"])[:initial_count]
    u0 = _coerce_vector(ctx, fixed["U0_Z_taylor"])[:initial_count]
    p0 = _coerce_vector(ctx, fixed["P0_Z_taylor"])[:initial_count]
    z = pad([z_center, convert(1)], initial_count)
    one = const(1, initial_count)
    d = add(one, scale(mul(z, z), -1))
    L = add(one, scale(mul(z, z), -dt))
    invL = inverse(L)

    # Rows intentionally shrink by one coefficient per radial order.  Since
    # initial_count was checked against degree + required_depth + 1, the final
    # row retains exactly (or more than) required_depth + 1 coefficients.
    a: list[list[Any]] = [const(1, initial_count)]
    u: list[list[Any]] = [u0[:initial_count]]
    pressure: list[list[Any]] = [p0[:initial_count]]

    for n in range(degree):
        K = initial_count - n - 1
        # Keep the full input rows while producing only K output entries.
        # In particular, ndiff needs the (K + 1)-st input coefficient for
        # the last retained derivative.  Slicing inputs to K here would drop
        # that term and silently corrupt higher axial coefficients.
        z_n = z
        d_n = d
        invL_n = invL
        ell_n = ell
        squared_n = squared

        def nadd(*args: list[Any]):
            return [sum((arg[k] for arg in args), zero) for k in range(K)]

        def nscale(values: list[Any], scalar: Any):
            return [value * scalar for value in values]

        def nmul(left: list[Any], right: list[Any]):
            return [
                sum((left[i] * right[k - i] for i in range(k + 1)), zero)
                for k in range(K)
            ]

        def ndiff(values: list[Any]):
            return [
                (k + 1) * values[k + 1] if k + 1 < len(values) else zero
                for k in range(K)
            ]

        def nconst(value: Any):
            return [value] + [zero] * (K - 1)

        W: list[list[Any]] = []
        H: list[list[Any]] = []
        for i in range(n + 1):
            wi = nscale(
                nadd(
                    nscale(nmul(z_n, u[i]), 1 - dt),
                    nmul(d_n, ndiff(u[i])),
                ),
                -convert(1) / (i + 1),
            )
            if i == 0:
                wi = nadd(nconst(1), wi)
            W.append(wi)
            hi = nmul(d_n, u[i])
            if i == 0:
                hi = nadd(hi, nscale(z_n, (1 - dt) / 2))
            H.append(hi)

        aa_u = nconst(0)
        uu = nconst(0)
        rhs_a = nconst(0)
        rhs_u = nconst(0)
        for i in range(n + 1):
            ai = a[n - i]
            ui = u[n - i]
            aa_u = nadd(aa_u, nmul(u[i], ai))
            uu = nadd(uu, nmul(u[i], ui))
            a_derivative = nadd(ndiff(ai), nmul(ell_n, ai))
            rhs_a = nadd(
                rhs_a,
                nscale(nmul(W[i], ai), n - i + 1),
                nmul(H[i], a_derivative),
            )
            rhs_u = nadd(
                rhs_u,
                nscale(nmul(W[i], ui), n - i),
                nmul(H[i], ndiff(ui)),
            )

        rhs_a = nadd(
            rhs_a,
            nscale(nadd(a[n], nscale(nmul(z_n, aa_u), -2)), dt / 2),
        )
        rhs_u = nadd(
            rhs_u,
            nscale(nadd(u[n], nscale(nmul(z_n, uu), -2)), (1 + dt) / 2),
            nmul(d_n, ndiff(pressure[n])),
            nscale(nmul(z_n, pressure[n]), -2 * (1 + dt)),
        )

        if n:
            aa_previous = nconst(0)
            for i in range(n):
                aa_previous = nadd(
                    aa_previous,
                    nmul(a[i], a[n - 1 - i]),
                )
            rhs_u = nadd(
                rhs_u,
                nscale(nmul(z_n, nmul(squared_n, aa_previous)), -2),
            )

        a_next = nscale(
            nmul(invL_n, rhs_a),
            convert(1) / (2 * (n + 1) * (n + 2)),
        )
        u_next = nscale(
            nmul(invL_n, rhs_u),
            convert(1) / (2 * (n + 1) ** 2),
        )
        aa_current = nconst(0)
        for i in range(n + 1):
            aa_current = nadd(
                aa_current,
                nmul(a[i], a[n - i]),
            )
        p_next = nscale(nmul(squared_n, aa_current), convert(1) / (n + 1))

        a.append(a_next)
        u.append(u_next)
        pressure.append(p_next)

    return {
        "A": a,
        "Uz": u,
        "P": pressure,
        "center": z_center,
        "delta": dt,
        "radial_degree": degree,
        "required_depth": required_depth,
        "initial_Z_degree": initial_count - 1,
        "axis_jet_length": initial_count,
        "gauge_factored": True,
        "axis_amplitude_normalized": True,
        "physical_pressure_rows": True,
        "center_interval_coefficient_enclosure": True,
        "single_point_taylor_series": False,
        "equation_provenance": {
            "paper": "Lei--Ren Part I, Section 8.2, equations (8.1)--(8.2)",
            "paper_text_lines": "work_paper_cache/lei_ren_part1.txt:12275-12332",
            "local_recurrence": "lei_ren_part1_paper_candidate_gauge_core.py:_advance_one",
            "factored_identity": "F=F0*A, F_Z=F0*(A_Z+(F0'/F0)*A)",
            "pressure_identity": "P[0]=P0; P[n+1]=S*sum(A[i]*A[n-i])/(n+1)",
            "swirl_identity": "F^2=S*A^2 in the axial cross term",
        },
        "scope": "Finite directed interval radial jet only; no infinite radial remainder, whole-axis certification, or terminal matching.",
    }


def _taylor(function, center: mp.mpf, length: int):
    return [
        mp.diff(function, center, order) / mp.factorial(order)
        for order in range(length)
    ]


def _interval_coefficients_for_center(
    ctx: MPIntervalContext, center: Any, length: int
) -> dict[str, list[Any]]:
    """Analytic test jets enclosing the fixture data over a center interval."""

    c = _as_interval(ctx, center)
    zero = ctx.mpf(0)

    def series_add(left, right):
        return [left[k] + right[k] for k in range(length)]

    def series_mul(left, right):
        return [
            sum((left[i] * right[k - i] for i in range(k + 1)), zero)
            for k in range(length)
        ]

    def series_inverse(values):
        result = [1 / values[0]]
        for k in range(1, length):
            result.append(
                -sum(
                    values[i] * result[k - i] for i in range(1, k + 1)
                )
                / values[0]
            )
        return result

    f0 = [ctx.exp(-2 * c * c)]
    f0.extend([zero] * (length - 1))
    for n in range(length - 1):
        previous = f0[n - 1] if n else zero
        f0[n + 1] = (-4 * c * f0[n] - 4 * previous) / (n + 1)

    squared = [ctx.exp(-4 * c * c)]
    squared.extend([zero] * (length - 1))
    for n in range(length - 1):
        previous = squared[n - 1] if n else zero
        squared[n + 1] = (-8 * c * squared[n] - 8 * previous) / (n + 1)

    ell = [-4 * c, ctx.mpf(-4)] + [zero] * (length - 2)
    u0 = [4 * c + ctx.mpf(".01"), ctx.mpf(4)] + [zero] * (length - 2)

    q = [1 + c * c, 2 * c, ctx.mpf(1)] + [zero] * (length - 3)
    q_inverse = series_inverse(q)
    p0 = [-3 * value for value in series_mul(q_inverse, q_inverse)]
    return {
        "F0_Z_taylor": f0,
        "ell_Z_taylor": ell,
        "S_Z_taylor": squared,
        "U0_Z_taylor": u0,
        "P0_Z_taylor": p0,
    }


def _point_fixture(center: mp.mpf, length: int) -> dict[str, list[mp.mpf]]:
    """Independent point Taylor data for the original recurrence checks."""

    return {
        "F0_Z_taylor": _taylor(lambda x: mp.exp(-2 * x * x), center, length),
        "ell_Z_taylor": _taylor(lambda x: -4 * x, center, length),
        "S_Z_taylor": _taylor(lambda x: mp.exp(-4 * x * x), center, length),
        "U0_Z_taylor": _taylor(
            lambda x: 4 * x + mp.mpf(".01"), center, length
        ),
        "P0_Z_taylor": _taylor(
            lambda x: -3 / (1 + x * x) ** 2, center, length
        ),
    }


def _contains_all(interval_rows, point_rows, tolerance: mp.mpf = mp.mpf("1e-170")):
    checks = 0
    maximum_gap = mp.mpf(0)
    for interval_row, point_row in zip(interval_rows, point_rows):
        for interval_value, point_value in zip(interval_row, point_row):
            lo, hi = _bounds(interval_value)
            maximum_gap = max(maximum_gap, lo - point_value, point_value - hi)
            if not (lo - tolerance <= point_value <= hi + tolerance):
                return False, checks + 1, maximum_gap
            checks += 1
    return True, checks, maximum_gap


def _bounds(value: Any):
    if hasattr(value, "_mpi_"):
        return endpoints(value)
    scalar = mp.mpf(value)
    return scalar, scalar


def run_self_tests() -> dict[str, Any]:
    """Run degree-4 point and non-degenerate center-interval tests."""

    degree = 4
    required_depth = 3
    length = degree + required_depth + 1
    centers = ("-.4", "0")
    records: list[dict[str, Any]] = []
    with mp.workdps(PRECISION + 20):
        ctx = MPIntervalContext()
        ctx.dps = PRECISION
        for center_text in centers:
            center = mp.mpf(center_text)
            width = mp.mpf("1e-40")
            center_interval = ctx.mpf(
                [mp.mpf(center - width), mp.mpf(center + width)]
            )
            interval_fixed = _interval_coefficients_for_center(
                ctx, center_interval, length
            )
            fixed = {
                key: list(values)
                for key, values in interval_fixed.items()
                if key != "F0_Z_taylor"
            }
            rows = coupled_rows(
                ctx,
                fixed,
                center_interval,
                degree,
                ".01",
                required_depth=required_depth,
            )

            point_fixed = _point_fixture(center, length)
            original = core_coefficients(
                center,
                ".01",
                F0_Z_taylor=point_fixed["F0_Z_taylor"],
                U0_Z_taylor=point_fixed["U0_Z_taylor"],
                P0_Z_taylor=point_fixed["P0_Z_taylor"],
                radial_degree=degree,
                precision=PRECISION,
            )
            factored = factored_core_coefficients(
                center,
                ".01",
                ell_Z_taylor=point_fixed["ell_Z_taylor"],
                S_Z_taylor=point_fixed["S_Z_taylor"],
                U0_Z_taylor=point_fixed["U0_Z_taylor"],
                P0_Z_taylor=point_fixed["P0_Z_taylor"],
                radial_degree=degree,
                precision=PRECISION,
            )

            checks = 0
            maximum_scaled_error = mp.mpf(0)
            maximum_reconstruction_error = mp.mpf(0)
            # Compare A with the independent factored recurrence and the
            # physical Uz/P rows with the independent unfactored recurrence.
            for name, reference_rows in (
                ("A", factored["A"]),
                ("Uz", original["Uz"]),
                ("P", original["P"]),
            ):
                ok, count, gap = _contains_all(rows[name], reference_rows)
                if not ok:
                    raise AssertionError(
                        f"center {center_text} {name} interval misses point result; gap {gap}"
                    )
                checks += count
                for interval_row, reference_row in zip(rows[name], reference_rows):
                    for interval_value, reference in zip(interval_row, reference_row):
                        lo, hi = _bounds(interval_value)
                        midpoint = (lo + hi) / 2
                        maximum_scaled_error = max(
                            maximum_scaled_error,
                            abs(midpoint - reference)
                            / max(mp.mpf(1), abs(reference), abs(midpoint)),
                        )

            # Reconstruct physical F from the interval A and point F0.
            for n, (a_row, f_row) in enumerate(zip(rows["A"], original["F"])):
                f0 = [
                    ctx.mpf(mp.nstr(value, 190))
                    for value in point_fixed["F0_Z_taylor"]
                ]
                reconstructed = [
                    sum(
                        f0[i] * a_row[k - i]
                        for i in range(k + 1)
                    )
                    for k in range(len(f_row))
                ]
                ok, count, gap = _contains_all(
                    [reconstructed], [f_row]
                )
                if not ok:
                    raise AssertionError(
                        f"center {center_text} reconstructed F misses point result; gap {gap}"
                    )
                checks += count
                for values, reference_row in zip([reconstructed], [f_row]):
                    for interval_value, reference in zip(values, reference_row):
                        lo, hi = _bounds(interval_value)
                        maximum_reconstruction_error = max(
                            maximum_reconstruction_error,
                            abs((lo + hi) / 2 - reference)
                            / max(mp.mpf(1), abs(reference)),
                        )

            records.append(
                {
                    "center": center_text,
                    "center_interval": [
                        mp.nstr(center - width, 60),
                        mp.nstr(center + width, 60),
                    ],
                    "radial_degree": degree,
                    "required_depth": required_depth,
                    "jet_length": length,
                    "row_lengths": {
                        name: [len(row) for row in rows[name]]
                        for name in ("A", "Uz", "P")
                    },
                    "containment_checks": checks,
                    "maximum_midpoint_scaled_error": mp.nstr(
                        maximum_scaled_error, 30
                    ),
                    "maximum_reconstructed_F_scaled_error": mp.nstr(
                        maximum_reconstruction_error, 30
                    ),
                    "all_checks_passed": True,
                }
            )

        # Also verify the exact scalar-center .3 recurrence against the
        # unchanged directed production algebra at degree 4 by comparing it
        # with both independent recurrences.
        center = mp.mpf(".3")
        point_fixed = _point_fixture(center, length)
        scalar_interval_fixed = _interval_coefficients_for_center(
            ctx, ctx.mpf(".3"), length
        )
        scalar_fixed = {
            key: list(values)
            for key, values in scalar_interval_fixed.items()
            if key != "F0_Z_taylor"
        }
        scalar_rows = coupled_rows(
            ctx,
            scalar_fixed,
            ctx.mpf(".3"),
            degree,
            ".01",
            required_depth=required_depth,
        )
        legacy_rows = {
            "A": [[ctx.mpf(1)] + [ctx.mpf(0)] * (length - 1)],
            "Uz": [list(scalar_fixed["U0_Z_taylor"])],
            "P": [list(scalar_fixed["P0_Z_taylor"])],
        }
        for n in range(degree):
            _legacy_advance_one(
                ctx,
                scalar_fixed,
                legacy_rows,
                n,
                delta_value=".01",
            )
        exact_interval_count = 0
        exact_endpoint_count = 0
        for name in ("A", "Uz", "P"):
            for new_row, old_row in zip(scalar_rows[name], legacy_rows[name]):
                if len(new_row) != len(old_row):
                    raise AssertionError(
                        f"scalar .3 {name} row length differs from legacy recurrence"
                    )
                for new_value, old_value in zip(new_row, old_row):
                    if not hasattr(new_value, "_mpi_") or not hasattr(
                        old_value, "_mpi_"
                    ):
                        raise AssertionError(
                            "legacy comparison requires interval coefficients"
                        )
                    if new_value._mpi_ != old_value._mpi_:
                        raise AssertionError(
                            f"scalar .3 exact legacy mismatch in {name}"
                        )
                    exact_interval_count += 1
                    exact_endpoint_count += 2
        original = core_coefficients(
            center,
            ".01",
            F0_Z_taylor=point_fixed["F0_Z_taylor"],
            U0_Z_taylor=point_fixed["U0_Z_taylor"],
            P0_Z_taylor=point_fixed["P0_Z_taylor"],
            radial_degree=degree,
            precision=PRECISION,
        )
        factored = factored_core_coefficients(
            center,
            ".01",
            ell_Z_taylor=point_fixed["ell_Z_taylor"],
            S_Z_taylor=point_fixed["S_Z_taylor"],
            U0_Z_taylor=point_fixed["U0_Z_taylor"],
            P0_Z_taylor=point_fixed["P0_Z_taylor"],
            radial_degree=degree,
            precision=PRECISION,
        )
        point_checks = 0
        for name, reference_rows in (
            ("A", factored["A"]),
            ("Uz", original["Uz"]),
            ("P", original["P"]),
        ):
            ok, count, gap = _contains_all(scalar_rows[name], reference_rows)
            if not ok:
                raise AssertionError(f"scalar .3 {name} misses point result; gap {gap}")
            point_checks += count

        # Exercise the existing pack/unpack helpers without changing any
        # production state or file.
        serialized = _pack_rows(scalar_rows["A"])
        unpacked = _unpack_rows(ctx, serialized)
        if len(unpacked) != len(scalar_rows["A"]):
            raise AssertionError("pack/unpack changed radial row count")

    return {
        "schema_version": 1,
        "all_checks_passed": True,
        "api": "coupled_rows(ctx, fixed, center, degree, delta, required_depth=3)",
        "precision": PRECISION,
        "scalar_center_point_checks": point_checks,
        "exact_interval_comparison_count": exact_interval_count,
        "exact_endpoint_comparison_count": exact_endpoint_count,
        "nondegenerate_centers": records,
        "source_hashes": _source_hashes(),
        "equation_provenance": {
            "paper": "Lei--Ren Part I, Section 8.2, equations (8.1)--(8.2)",
            "paper_text_lines": "work_paper_cache/lei_ren_part1.txt:12275-12332",
            "production_recurrence": "candidate_gauge_core.py:_advance_one",
            "factored_recurrence": "amplitude_factored_core.py:factored_core_coefficients",
            "unfactored_recurrence": "core_recursion.py:core_coefficients",
        },
        "scope": "Degree-4 directed algebra and interval-center containment only; no actual candidate 124-order run, whole-axis core, or terminal moment closure.",
    }


def main() -> None:
    result = run_self_tests()
    output = HERE / "lei_ren_part1_paper_functional_core_recursion.json"
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
