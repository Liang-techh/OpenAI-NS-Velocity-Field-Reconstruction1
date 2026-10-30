"""Independent synthetic checks for the four-cell W0 collar leaf.

The fixture deliberately knows only the public ``build`` contract.  It makes
constant signed source rows, evaluates the six Green kernels by an unrelated
high-precision quadrature, and checks the corresponding output cell.  The
formal source atoms are never materialised: only the formal ``hb`` exponent
introduced by a kernel is inspected/evaluated.

This is a bookkeeping fixture, not a physical source or a global collar
certificate.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Any


FIELDS = ("F1", "Uz1", "K1", "P1", "d_x_F1", "d_x_Uz1")
ROWS = ("5", "6", "6", "4", "5", "6")
ATOMS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")


def _rational(c: Any, value: Fraction) -> Any:
    """Make an exact directed point interval from a rational."""

    return c.mpf(value.numerator) / value.denominator


def _lookup(mapping: Any, key: Any) -> Any:
    """Read either integer- or string-keyed mappings used by leaf receipts."""

    if isinstance(mapping, dict):
        if key in mapping:
            return mapping[key]
        text = str(key)
        if text in mapping:
            return mapping[text]
        if isinstance(key, str) and key.isdigit() and int(key) in mapping:
            return mapping[int(key)]
    if isinstance(mapping, (list, tuple)) and isinstance(key, int):
        return mapping[key]
    raise KeyError(key)


def _same_interval(left: Any, right: Any) -> bool:
    """Compare directed intervals without converting them to floats."""

    try:
        return left._mpi_ == right._mpi_
    except AttributeError:
        return left == right


def _term(c: Any, value: Fraction, source_id: str, *, uncertainty: bool = False,
          error_class: str = "synthetic_source") -> dict[str, Any]:
    return {
        "powers": [0, 0, 0, 0],
        "coefficient": _rational(c, value),
        "origin": {
            "source_ids": [source_id],
            "uncertainty": uncertainty,
            "error_classes": [error_class],
        },
    }


def _core_inputs(c: Any) -> dict[int, dict[str, list[dict[str, Any]]]]:
    """Create core terms with nontrivial formal powers and signed receipts."""

    out: dict[int, dict[str, list[dict[str, Any]]]] = {}
    for order in (0, 1):
        out[order] = {}
        for index, field in enumerate(FIELDS):
            # These atoms are intentionally left formal.  The leaf must copy
            # them to core_transport rather than evaluating or dropping them.
            powers = [0, (order + index) % 3, 1 + (index % 2), -1 - order]
            sign = -1 if (order + index) % 2 else 1
            value = Fraction(sign * (order + 2) * (index + 1), 23 + index)
            out[order][field] = [{
                "powers": powers,
                "coefficient": _rational(c, value),
                "origin": {
                    "source_ids": [f"core:{order}:{field}"],
                    "uncertainty": order == 1,
                    "error_classes": ["synthetic_core"],
                },
            }]
    return out


def _forcing_inputs(c: Any) -> tuple[dict[int, dict[str, list[dict[str, Any]]]],
                    dict[tuple[int, str], Fraction]]:
    """Create constant, known-signed ``g0`` and ``g1`` rows."""

    values: dict[int, dict[str, Fraction]] = {
        0: {"4": Fraction(-2, 7), "5": Fraction(-3, 11), "6": Fraction(5, 13)},
        1: {"4": Fraction(7, 17), "5": Fraction(-11, 19), "6": Fraction(13, 23)},
    }
    out: dict[int, dict[str, list[dict[str, Any]]]] = {}
    for order in (0, 1):
        out[order] = {}
        for row, value in values[order].items():
            out[order][row] = [_term(
                c, value, f"g{order}:{row}", uncertainty=order == 1,
                error_class=f"synthetic_g{order}",
            )]
    return out, {(order, row): value for order in values for row, value in values[order].items()}


def _formal_value(c: Any, terms: Any, hb: Any) -> Any:
    """Evaluate only the synthetic formal point ``(hb,1,1,1)`` in intervals."""

    bases = (hb, c.mpf(1), c.mpf(1), c.mpf(1))
    total = c.mpf(0)
    for item in terms:
        value = item["coefficient"]
        for base, power in zip(bases, item["powers"]):
            value *= base ** int(power)
        total += value
    return total


def _direct_kernel(mp: Any, hb: Any, s: Any, v: Any, field_index: int) -> Any:
    """Return one of the six requested kernels, excluding ``dx``."""

    u = 2 * mp.exp(hb * v / 2)
    d = hb * (s - v)
    if d == 0:
        exprel = mp.mpf(1)
    else:
        exprel = (1 - mp.exp(-d)) / d
    first = u / 2 * hb * (s - v) * exprel
    return (
        first,
        u / 2 * hb * (s - v),
        -first,
        mp.mpf(1),
        mp.exp(-3 * d / 2),
        mp.exp(-d / 2),
    )[field_index]


def run(leaf: Any, c: Any) -> dict[str, Any]:
    """Run an independent synthetic direct-integral audit of ``leaf.build``."""

    import mpmath as mp

    # Keep the quadrature independent of the interval context used by the
    # leaf.  The four representatives are q=(1/8,3/8,5/8,7/8), with
    # s=sc*q/4 in the first collar.
    mp.mp.dps = 90
    hb_mp = mp.mpf(1) / 3
    sc_mp = mp.mpf(1) / 2048
    sc = c.mpf(1) / 2048
    core = _core_inputs(c)
    g, g_values = _forcing_inputs(c)
    source_atoms = {name: f"synthetic:{name}" for name in ATOMS}
    certificate = {
        "hb_strictly_between_zero_and_half": True,
        "sc_phase_Z_independent": True,
        "hb_epsilon_Z_independent": True,
    }

    result = leaf.build(
        c,
        sc=sc,
        core=core,
        g=g,
        source_family="synthetic_w0_collar",
        source_atoms=source_atoms,
        certificate=certificate,
    )

    checks: dict[str, bool] = {}

    # The leaf has six fixed field/forcing-row pairings and must retain the
    # separate formal core rather than replacing it by a numerical value.
    checks["fixed_fields"] = tuple(getattr(leaf, "FIELDS", ())) == FIELDS
    checks["fixed_rows"] = tuple(getattr(leaf, "ROWS", ())) == ROWS
    checks["source_family"] = result.get("source_family") == "synthetic_w0_collar"
    checks["source_atoms"] = result.get("source_atoms") == source_atoms
    checks["formal_power_basis"] = tuple(result.get("source_power_basis", ())) == ATOMS

    core_out = result.get("core_transport")
    checks["core_orders"] = isinstance(core_out, dict) and set(core_out) in ({0, 1}, {"0", "1"})
    if checks["core_orders"]:
        for order in (0, 1):
            order_out = _lookup(core_out, order)
            checks[f"core_fields_{order}"] = isinstance(order_out, dict) and set(order_out) == set(FIELDS)
            if checks[f"core_fields_{order}"]:
                for field in FIELDS:
                    original = core[order][field]
                    saved = _lookup(order_out, field)
                    checks[f"core_preserved_{order}_{field}"] = (
                        len(saved) == len(original)
                        and all(
                            tuple(a["powers"]) == tuple(b["powers"])
                            and _same_interval(a["coefficient"], b["coefficient"])
                            and a.get("origin", {}).get("source_ids")
                            == b.get("origin", {}).get("source_ids")
                            for a, b in zip(saved, original)
                        )
                    )

    local = result.get("local")
    checks["local_orders"] = isinstance(local, dict) and set(local) in ({0, 1}, {"0", "1"})
    if checks["local_orders"]:
        for order in (0, 1):
            order_local = _lookup(local, order)
            checks[f"local_cells_{order}"] = isinstance(order_local, dict) and set(order_local) in (
                {0, 1, 2, 3}, {"0", "1", "2", "3"}
            )
            if not checks[f"local_cells_{order}"]:
                continue
            for cell in range(4):
                cell_out = _lookup(order_local, cell)
                complete = isinstance(cell_out, dict) and set(cell_out) == set(FIELDS)
                checks[f"complete_cell_{order}_{cell}"] = complete
                if not complete:
                    continue

                q_mid = mp.mpf(2 * cell + 1) / 8
                s = sc_mp * q_mid / 4
                for index, field in enumerate(FIELDS):
                    terms = _lookup(cell_out, field)
                    nonempty = isinstance(terms, (list, tuple)) and bool(terms)
                    checks[f"nonempty_{order}_{cell}_{field}"] = nonempty
                    if not nonempty:
                        continue

                    row = ROWS[index]
                    source_id = f"g{order}:{row}"
                    expected_hb_power = 2 if index < 3 else 1
                    checks[f"hb_power_{order}_{cell}_{field}"] = all(
                        int(item["powers"][0]) == expected_hb_power
                        and tuple(item["powers"][1:]) == (0, 0, 0)
                        for item in terms
                    )
                    checks[f"provenance_{order}_{cell}_{field}"] = all(
                        source_id in item.get("origin", {}).get("source_ids", [])
                        for item in terms
                    )

                    # The integral is evaluated at an interior point of this
                    # q-cell and must lie in the output interval for the
                    # corresponding formal sector.
                    source_value = mp.mpf(g_values[(order, row)].numerator) / g_values[(order, row)].denominator
                    reference = mp.quad(
                        lambda v: source_value
                        * _direct_kernel(mp, hb_mp, s, v, index)
                        * (hb_mp * (2 * mp.exp(hb_mp * v / 2)) / 2),
                        [0, s],
                    )
                    expected = c.mpf(mp.nstr(reference, 72))
                    got = _formal_value(c, terms, c.mpf(1) / 3)
                    checks[f"direct_{order}_{cell}_{field}"] = bool(
                        got.a <= expected.a and got.b >= expected.b
                    )
                    if index == 2:
                        # g[order]['6'] is positive in both synthetic rows;
                        # K1 therefore retains the signed negative kernel.
                        checks[f"negative_K_{order}_{cell}"] = bool(got.b <= 0 and got.a < 0)

    q_cells = result.get("q_cells")
    checks["four_q_cells"] = isinstance(q_cells, (list, tuple)) and len(q_cells) == 4
    if checks["four_q_cells"]:
        for cell, q_data in enumerate(q_cells):
            checks[f"q_cell_{cell}_complete"] = all(
                key in q_data for key in ("q_lower", "q_upper", "s_lower", "s_upper")
            )
            if checks[f"q_cell_{cell}_complete"]:
                q_lower = c.mpf(cell) / 4
                q_upper = c.mpf(cell + 1) / 4
                s_lower = sc * cell / 16
                s_upper = sc * (cell + 1) / 16
                checks[f"q_cell_{cell}_coordinates"] = (
                    q_data["q_lower"]._mpi_ == q_lower._mpi_
                    and q_data["q_upper"]._mpi_ == q_upper._mpi_
                    and q_data["s_lower"]._mpi_ == s_lower._mpi_
                    and q_data["s_upper"]._mpi_ == s_upper._mpi_
                )

    # These are explicit scope guards: no tail, history reset, hb
    # materialisation, or same-source cancellation may be inferred by the
    # synthetic fixture.
    checks["core_history_preserved"] = result.get("core_history_reset") is False
    checks["local_added_once"] = result.get("local_added_once_in_definition") is True
    checks["no_tail"] = result.get("allword_tail_added") is False
    checks["no_m15"] = result.get("m15_computed") is False
    checks["hb_remains_formal"] = result.get("hb_materialized") is False
    checks["no_same_source_cancellation"] = result.get("same_source_cancellation_claimed") is False

    return {
        "passed": all(checks.values()),
        "check_count": len(checks),
        "checks": checks,
        "synthetic_only": True,
        "reference": "independent high-precision direct integration at four interior q-cell representatives",
        "hb": "1/3",
        "sc": "1/2048",
        "cell_count": 4,
    }


__all__ = ["FIELDS", "ROWS", "ATOMS", "run"]
