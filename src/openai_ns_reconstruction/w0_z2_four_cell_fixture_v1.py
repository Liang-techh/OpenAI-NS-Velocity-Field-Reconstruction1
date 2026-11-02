"""Independent synthetic checks for :mod:`w0_z2_four_cell_v1`.

The fixture uses only the public ``build(c, *, sc, core, g)`` contract.  Its
source intervals include negative and sign-crossing values, and its physical
reference is an unrelated high-precision direct integral of the three stated
Green kernels.  It is synthetic bookkeeping, not a physical receipt.
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
    return c.make_mpf((_rational(c, lower)._mpi_[0], _rational(c, upper)._mpi_[1]))


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
    classes = list(error_classes or ["synthetic_source"])
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
            "operation": "synthetic_source",
        },
    }


def _core(c: Any) -> dict[str, list[dict[str, Any]]]:
    out: dict[str, list[dict[str, Any]]] = {}
    for index, field in enumerate(FIELDS):
        out[field] = [_term(
            c,
            Fraction(index + 1, 19),
            Fraction(index + 2, 19),
            f"core:{field}",
            f"synthetic/core/{field}",
            powers=[index % 2, 1 + index, -index, -1],
            uncertainty=index == 2,
            error_classes=["synthetic_core", f"core_{field}"],
        )]
    return out


def _sources(c: Any) -> tuple[dict[str, list[dict[str, Any]]], dict[str, Fraction]]:
    # g5 has a negative and a sign-crossing interval; g6 has a positive and a
    # sign-crossing interval.  The selected scalar is only for the independent
    # synthetic integral and is never passed to the leaf.
    rows = {
        "5": [
            (-Fraction(3, 2), -Fraction(1, 2), "g5:negative", -Fraction(1, 1)),
            (-Fraction(1, 2), Fraction(1, 3), "g5:cross", Fraction(1, 5)),
        ],
        "6": [
            (Fraction(1, 2), Fraction(3, 2), "g6:positive", Fraction(1, 1)),
            (-Fraction(2, 3), Fraction(1, 4), "g6:cross", -Fraction(1, 4)),
        ],
    }
    out: dict[str, list[dict[str, Any]]] = {"5": [], "6": []}
    selected: dict[str, Fraction] = {}
    for row, entries in rows.items():
        for lower, upper, source_id, scalar in entries:
            path = f"synthetic/g/{row}/{source_id}"
            out[row].append(_term(
                c,
                lower,
                upper,
                source_id,
                path,
                uncertainty=row == "6",
                error_classes=["synthetic_g", f"row_{row}"],
            ))
            selected[source_id] = scalar
    return out, selected


def _lookup(mapping: Any, key: Any) -> Any:
    if isinstance(mapping, dict):
        if key in mapping:
            return mapping[key]
        text = str(key)
        if text in mapping:
            return mapping[text]
    if isinstance(mapping, (list, tuple)) and isinstance(key, int):
        return mapping[key]
    raise KeyError(key)


def _formal_value(c: Any, terms: Any, hb: Any) -> Any:
    bases = (hb, c.mpf(1), c.mpf(1), c.mpf(1))
    total = c.mpf(0)
    for item in terms:
        value = item["coefficient"]
        for base, power in zip(bases, item["powers"]):
            value *= base ** int(power)
        total += value
    return total


def _direct_kernel(mp: Any, x: Any, u: Any, field: str) -> Any:
    """Kernel from the stated inequalities, as a function of ``u``."""

    ratio = u / x
    f1 = u / 2 * (1 - ratio ** 2)
    uz1 = u * mp.log(x / u)
    return {"F1": f1, "Uz1": uz1, "K1": -f1}[field]


def _direct_integral(mp: Any, hb: Any, s: Any, field: str, source: Any) -> Any:
    x = 2 * mp.exp(hb * s / 2)

    def integrand(u: Any) -> Any:
        return source * _direct_kernel(mp, x, u, field)

    return mp.quad(integrand, [mp.mpf(2), x])


def run(leaf: Any, c: Any) -> dict[str, Any]:
    """Run the synthetic four-cell and independent kernel-bound audit."""

    import mpmath as mp

    mp.mp.dps = 90
    hb_mp = mp.mpf(1) / 3
    sc_mp = mp.mpf(1) / 2048
    sc = c.mpf(1) / 2048
    core = _core(c)
    g, selected = _sources(c)
    build = leaf if callable(leaf) else leaf.build
    result = build(c, sc=sc, core=core, g=g)
    checks: dict[str, bool] = {}

    checks["fixed_fields"] = tuple(getattr(leaf, "FIELDS", FIELDS)) == FIELDS
    checks["fixed_rows"] = getattr(leaf, "ROWS", ROWS) == ROWS
    checks["ordinary_Z2"] = result.get("ordinary_Z_order") == 2
    checks["ordinary_no_factorial"] = result.get("ordinary_Z_factorial_divided") is False
    checks["ordinary_source_marker"] = result.get("source_z2_is_ordinary_derivative") is True
    checks["formal_basis"] = tuple(result.get("source_power_basis", ())) == ATOMS
    checks["hb_stays_formal"] = result.get("hb_materialized") is False
    checks["formal_hb2_stays_formal"] = result.get("formal_hb2_materialized") is False
    checks["no_physical_cartesian_measure"] = (
        result.get("no_physical_cartesian_measure") is True
        and result.get("physical_or_cartesian_measure_included") is False
    )
    checks["source_rows_nonempty"] = all(
        isinstance(core[field], (list, tuple)) and bool(core[field]) for field in FIELDS
    ) and all(isinstance(g[row], (list, tuple)) and bool(g[row]) for row in ("5", "6"))

    core_out = result.get("core_transport")
    checks["core_fields"] = isinstance(core_out, dict) and set(core_out) == set(FIELDS)
    if checks["core_fields"]:
        for field in FIELDS:
            saved = core_out[field]
            original = core[field]
            checks[f"core_preserved_{field}"] = (
                len(saved) == len(original)
                and all(
                    tuple(a["powers"]) == tuple(b["powers"])
                    and _same_interval(a["coefficient"], b["coefficient"])
                    and a.get("origin") == b.get("origin")
                    for a, b in zip(saved, original)
                )
            )

    q_cells = result.get("q_cells")
    local = result.get("local")
    checks["four_q_cells"] = isinstance(q_cells, (list, tuple)) and len(q_cells) == 4
    checks["four_local_cells"] = isinstance(local, dict) and set(local) == {0, 1, 2, 3}
    if checks["four_q_cells"]:
        for index, cell in enumerate(q_cells):
            q_lower = c.mpf(index) / 4
            q_upper = c.mpf(index + 1) / 4
            s_lower = sc * q_lower / 4
            s_upper = sc * q_upper / 4
            checks[f"coordinates_{index}"] = (
                _same_interval(cell["q_lower"], q_lower)
                and _same_interval(cell["q_upper"], q_upper)
                and _same_interval(cell["s_lower"], s_lower)
                and _same_interval(cell["s_upper"], s_upper)
            )
            checks[f"endpoint_formula_{index}"] = _same_interval(
                cell["S_q"], (sc / 4) * q_upper * c.exp(sc / 16)
            )
            bound = cell["kernel_bound_without_formal_hb2"]
            checks[f"positive_bound_{index}"] = bool(bound.a >= 0 and bound.b >= bound.a)
            checks[f"bound_excludes_hb2_{index}"] = (
                cell.get("kernel_bound_is_without_formal_hb2") is True
                and cell.get("formal_hb2_materialized") is False
                and cell.get("physical_or_cartesian_measure_included") is False
            )
            checks[f"core_ref_{index}"] = set(cell["core_ref"]) == set(FIELDS)

    if checks["four_local_cells"]:
        for index in range(4):
            cell_out = _lookup(local, index)
            checks[f"complete_local_{index}"] = isinstance(cell_out, dict) and set(cell_out) == set(FIELDS)
            if not checks[f"complete_local_{index}"]:
                continue
            q_mid = mp.mpf(2 * index + 1) / 8
            s = sc_mp * q_mid / 4
            s_q_mp = (sc_mp / 4) * mp.mpf(index + 1) / 4 * mp.exp(sc_mp / 16)
            physical_upper = hb_mp ** 2 * s_q_mp ** 2 / 2
            for field in FIELDS:
                terms = cell_out[field]
                row = ROWS[field]
                checks[f"nonempty_{index}_{field}"] = bool(terms)
                checks[f"row_{index}_{field}"] = all(
                    item["origin"].get("source_row") == row for item in terms
                )
                checks[f"hb_power_{index}_{field}"] = all(
                    int(item["powers"][0]) == 2 for item in terms
                )
                checks[f"provenance_{index}_{field}"] = all(
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
                    expected = c.mpf(mp.nstr(reference, 72))
                    got = _formal_value(c, [item], c.mpf(1) / 3)
                    checks[f"direct_{index}_{field}_{source_id}"] = bool(
                        got.a <= expected.a and got.b >= expected.b
                    )
                    if field == "K1" and source_id == "g6:positive":
                        checks[f"negative_K_{index}"] = bool(got.b <= 0 and got.a < 0)
                    # The independent scalar integral also checks the stated
                    # positive physical bound, with the source magnitude kept
                    # outside the geometric factor.
                    checks[f"physical_bound_{index}_{field}_{source_id}"] = bool(
                        abs(reference) <= abs(selected[source_id]) * physical_upper
                    )

    # The source intervals themselves include a negative interval and an
    # interval crossing zero; this guards against accidental endpoint/value
    # selection in the new leaf.
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
        and q_cells[0]["x_minus_2_bound"]["interval_without_hb"].a >= 0
    )
    checks["no_tail_or_m15"] = (
        result.get("allword_tail_added") is False
        and result.get("m15_computed") is False
    )
    checks["no_old_W1"] = result.get("old_W1_consumed") is False

    return {
        "passed": all(checks.values()),
        "check_count": len(checks),
        "checks": checks,
        "synthetic_only": True,
        "reference": "independent high-precision direct integration of the F1, Uz1, and K1 kernels",
        "hb": "1/3",
        "sc": "1/2048",
        "cell_count": 4,
    }


def run_fixture(build: Any, c: Any) -> dict[str, Any]:
    """Runner-compatible alias accepting ``leaf.build`` directly."""

    return run(build, c)


__all__ = ["ATOMS", "FIELDS", "ROWS", "run", "run_fixture"]
