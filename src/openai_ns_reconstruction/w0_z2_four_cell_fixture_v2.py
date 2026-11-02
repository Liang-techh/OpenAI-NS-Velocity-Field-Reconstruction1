"""Fresh synthetic audit for :mod:`w0_z2_four_cell_v2`.

The fixture uses only ``run_fixture(leaf.build, c)``'s public build callable.
Its source intervals and identifiers are independent of the v1 fixture.  The
reference integral is over ``u`` in ``[2, X]`` and checks the stated kernel
bound without introducing a physical or Cartesian measure.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Any


FIELDS = ("F1", "Uz1", "K1")
ROWS = {"F1": "5", "Uz1": "6", "K1": "6"}
ATOMS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")


def _rational(c: Any, value: Fraction) -> Any:
    return c.mpf(value.numerator) / value.denominator


def _interval(c: Any, lower: Fraction, upper: Fraction) -> Any:
    lo = _rational(c, lower)
    hi = _rational(c, upper)
    return c.make_mpf((lo._mpi_[0], hi._mpi_[1]))


def _same_interval(left: Any, right: Any) -> bool:
    try:
        return left._mpi_ == right._mpi_
    except AttributeError:
        return left == right


def _term(
    c: Any,
    lower: Fraction,
    upper: Fraction,
    source_id: str,
    source_path: str,
    *,
    powers: list[int] | None = None,
    uncertainty: bool = False,
    error_classes: list[str] | None = None,
) -> dict[str, Any]:
    classes = list(error_classes or ["synthetic_v2_source"])
    refs = [{"source_id": source_id, "source_path": source_path}]
    return {
        "powers": list(powers or [0, 0, 0, 0]),
        "coefficient": _interval(c, lower, upper),
        "source_ids": [source_id],
        "source_path": source_path,
        "source_paths": [source_path],
        "source_node_refs": refs,
        "uncertainty": uncertainty,
        "error_classes": classes,
        "origin": {
            "source_ids": [source_id],
            "source_paths": [source_path],
            "source_node_refs": refs,
            "uncertainty": uncertainty,
            "error_classes": classes,
            "operation": "synthetic_v2_source",
        },
    }


def _core(c: Any) -> dict[str, list[dict[str, Any]]]:
    out: dict[str, list[dict[str, Any]]] = {}
    for index, field in enumerate(FIELDS):
        out[field] = [_term(
            c,
            Fraction(index + 2, 29),
            Fraction(index + 3, 29),
            f"v2:core:{field}",
            f"synthetic-v2/core/{field}",
            powers=[1 + index % 2, -index, 2 + index, -2],
            uncertainty=index == 1,
            error_classes=["synthetic_v2_core", f"v2_core_{field}"],
        )]
    return out


def _sources(c: Any) -> tuple[dict[str, list[dict[str, Any]]], dict[str, Fraction]]:
    rows = {
        "5": [
            (-Fraction(7, 4), -Fraction(3, 5), "v2:g5:negative", -Fraction(1, 1)),
            (-Fraction(3, 5), Fraction(2, 5), "v2:g5:cross", Fraction(1, 7)),
        ],
        "6": [
            (Fraction(2, 5), Fraction(7, 5), "v2:g6:positive", Fraction(4, 5)),
            (-Fraction(1, 2), Fraction(3, 10), "v2:g6:cross", -Fraction(1, 5)),
        ],
    }
    result: dict[str, list[dict[str, Any]]] = {"5": [], "6": []}
    selected: dict[str, Fraction] = {}
    for row, entries in rows.items():
        for lower, upper, source_id, scalar in entries:
            path = f"synthetic-v2/g/{row}/{source_id}"
            result[row].append(_term(
                c,
                lower,
                upper,
                source_id,
                path,
                uncertainty=row == "6",
                error_classes=["synthetic_v2_g", f"v2_row_{row}"],
            ))
            selected[source_id] = scalar
    return result, selected


def _formal_value(c: Any, terms: Any, hb: Any) -> Any:
    bases = (hb, c.mpf(1), c.mpf(1), c.mpf(1))
    total = c.mpf(0)
    for item in terms:
        value = item["coefficient"]
        for base, power in zip(bases, item["powers"]):
            value *= base ** int(power)
        total += value
    return total


def _kernel(mp: Any, x: Any, u: Any, field: str) -> Any:
    ratio = u / x
    f1 = u / 2 * (1 - ratio ** 2)
    uz1 = u * mp.log(x / u)
    return {"F1": f1, "Uz1": uz1, "K1": -f1}[field]


def _direct_integral(mp: Any, hb: Any, s: Any, field: str, source: Any) -> Any:
    x = 2 * mp.exp(hb * s / 2)
    return mp.quad(
        lambda u: source * _kernel(mp, x, u, field),
        [mp.mpf(2), x],
    )


def run_fixture(build: Any, c: Any) -> dict[str, Any]:
    """Run fresh synthetic checks against a callable ``leaf.build``."""

    import mpmath as mp

    mp.mp.dps = 100
    hb_mp = mp.mpf(2) / 5
    sc_mp = mp.mpf(1) / 4096
    sc = c.mpf(1) / 4096
    core = _core(c)
    g, selected = _sources(c)
    result = build(c, sc=sc, core=core, g=g)
    checks: dict[str, bool] = {}

    checks["fresh_v2_source_ids"] = all(
        term["origin"]["source_ids"][0].startswith("v2:")
        for row in g.values()
        for term in row
    )
    checks["ordinary_Z2"] = result.get("ordinary_Z_order") == 2
    checks["ordinary_no_factorial"] = result.get("ordinary_Z_factorial_divided") is False
    checks["formal_basis"] = tuple(result.get("source_power_basis", ())) == ATOMS
    checks["source_marker"] = result.get("source_z2_is_ordinary_derivative") is True
    checks["X_half_exponent"] = result.get("X_formula") == "X=2*exp(hb*s/2)"
    checks["hb_stays_formal"] = (
        result.get("hb_materialized") is False
        and result.get("formal_hb2_materialized") is False
    )
    checks["no_physical_cartesian_measure"] = (
        result.get("no_physical_cartesian_measure") is True
        and result.get("physical_or_cartesian_measure_included") is False
    )
    checks["source_rows_nonempty"] = (
        all(isinstance(core[field], (list, tuple)) and bool(core[field]) for field in FIELDS)
        and all(isinstance(g[row], (list, tuple)) and bool(g[row]) for row in ("5", "6"))
    )

    core_out = result.get("core_transport")
    checks["core_fields"] = isinstance(core_out, dict) and set(core_out) == set(FIELDS)
    if checks["core_fields"]:
        for field in FIELDS:
            checks[f"core_preserved_{field}"] = (
                len(core_out[field]) == len(core[field])
                and all(
                    tuple(a["powers"]) == tuple(b["powers"])
                    and _same_interval(a["coefficient"], b["coefficient"])
                    and a.get("origin") == b.get("origin")
                    for a, b in zip(core_out[field], core[field])
                )
            )

    q_cells = result.get("q_cells")
    local = result.get("local")
    checks["four_q_cells"] = isinstance(q_cells, (list, tuple)) and len(q_cells) == 4
    checks["four_local_cells"] = isinstance(local, dict) and set(local) == {0, 1, 2, 3}
    if checks["four_q_cells"]:
        for cell, data in enumerate(q_cells):
            q_lower = c.mpf(cell) / 4
            q_upper = c.mpf(cell + 1) / 4
            s_lower = sc * q_lower / 4
            s_upper = sc * q_upper / 4
            checks[f"coordinates_{cell}"] = (
                _same_interval(data["q_lower"], q_lower)
                and _same_interval(data["q_upper"], q_upper)
                and _same_interval(data["s_lower"], s_lower)
                and _same_interval(data["s_upper"], s_upper)
                and _same_interval(data["q_hi"], q_upper)
            )
            checks[f"S_q_formula_{cell}"] = _same_interval(
                data["S_q"], (sc / 4) * q_upper * c.exp(sc / 16)
            )
            bound = data["kernel_bound_without_formal_hb2"]
            checks[f"positive_bound_{cell}"] = bool(bound.a >= 0 and bound.b >= bound.a)
            checks[f"bound_without_hb2_{cell}"] = (
                data.get("kernel_bound_is_without_formal_hb2") is True
                and data.get("formal_hb2_materialized") is False
                and data.get("physical_or_cartesian_measure_included") is False
            )

    if checks["four_local_cells"]:
        for cell in range(4):
            cell_local = local[cell]
            checks[f"complete_local_{cell}"] = set(cell_local) == set(FIELDS)
            if not checks[f"complete_local_{cell}"]:
                continue
            q_mid = mp.mpf(2 * cell + 1) / 8
            s = sc_mp * q_mid / 4
            s_q_mp = (sc_mp / 4) * mp.mpf(cell + 1) / 4 * mp.exp(sc_mp / 16)
            physical_upper = hb_mp ** 2 * s_q_mp ** 2 / 2
            for field in FIELDS:
                row = ROWS[field]
                terms = cell_local[field]
                checks[f"nonempty_{cell}_{field}"] = bool(terms)
                checks[f"row_{cell}_{field}"] = all(
                    item["origin"].get("source_row") == row for item in terms
                )
                checks[f"hb_power_{cell}_{field}"] = all(
                    int(item["powers"][0]) == 2 for item in terms
                )
                checks[f"provenance_{cell}_{field}"] = all(
                    item.get("source_node_refs")
                    and item["origin"].get("source_ids")
                    and item["origin"].get("error_classes")
                    for item in terms
                )
                for item in terms:
                    source_id = item["origin"]["source_ids"][0]
                    reference = _direct_integral(
                        mp,
                        hb_mp,
                        s,
                        field,
                        selected[source_id],
                    )
                    expected = c.mpf(mp.nstr(reference, 78))
                    got = _formal_value(c, [item], c.mpf(2) / 5)
                    checks[f"direct_{cell}_{field}_{source_id}"] = bool(
                        got.a <= expected.a and got.b >= expected.b
                    )
                    checks[f"physical_bound_{cell}_{field}_{source_id}"] = bool(
                        abs(reference) <= abs(selected[source_id]) * physical_upper
                    )
                    if field == "K1" and source_id == "v2:g6:positive":
                        checks[f"negative_K_{cell}"] = bool(got.b <= 0 and got.a < 0)

    checks["negative_source_interval"] = bool(g["5"][0]["coefficient"].a < 0)
    checks["crossing_source_interval"] = bool(
        g["5"][1]["coefficient"].a < 0 < g["5"][1]["coefficient"].b
    )
    checks["core_not_duplicated"] = sum(len(core_out[field]) for field in FIELDS) == sum(
        len(core[field]) for field in FIELDS
    )

    missing_g6 = dict(g)
    missing_g6.pop("6")
    try:
        build(c, sc=sc, core=core, g=missing_g6)
    except (KeyError, TypeError, ValueError):
        checks["missing_g6_rejected"] = True
    else:
        checks["missing_g6_rejected"] = False

    missing_core = dict(core)
    missing_core.pop("K1")
    try:
        build(c, sc=sc, core=missing_core, g=g)
    except (KeyError, TypeError, ValueError):
        checks["missing_core_rejected"] = True
    else:
        checks["missing_core_rejected"] = False

    zero = c.mpf(0)
    zero_s = (sc / 4) * zero
    zero_sq = (sc / 4) * zero * c.exp(sc / 16)
    checks["q_zero_limit"] = (
        _same_interval(q_cells[0]["q_lower"], zero)
        and _same_interval(q_cells[0]["s_lower"], zero_s)
        and _same_interval(zero_sq, zero)
    )
    checks["scope_guards"] = (
        result.get("allword_tail_added") is False
        and result.get("m15_computed") is False
        and result.get("old_W1_consumed") is False
        and result.get("same_source_cancellation_claimed") is False
    )

    return {
        "passed": all(checks.values()),
        "check_count": len(checks),
        "checks": checks,
        "synthetic_only": True,
        "reference": "fresh-v2 high-precision direct u-integral of the stated kernels",
        "hb": "2/5",
        "sc": "1/4096",
        "cell_count": 4,
    }


def run(leaf: Any, c: Any) -> dict[str, Any]:
    """Compatibility alias accepting either a leaf object or its build."""

    return run_fixture(leaf if callable(leaf) else leaf.build, c)


__all__ = ["ATOMS", "FIELDS", "ROWS", "run", "run_fixture"]
