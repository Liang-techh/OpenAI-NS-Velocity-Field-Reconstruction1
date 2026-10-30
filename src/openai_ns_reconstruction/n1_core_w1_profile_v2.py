"""Cellwise target-radius Green enclosures for the core W1 word.

The old core leaf used one interval value such as ``K(t)/4`` for a source
cell.  That is a useful coarse bound, but it is not the Green mass of the
cell.  This leaf keeps the same formal source-sector convention and replaces
that value by the directed integral of each principal Green kernel over the
source cell.  The independent ``x=2`` exit is assembled separately; the last
target-cell interval is never used as an exit proxy.

The normalized radius is ``t=x/2``.  Eight cells are ``[i/8,(i+1)/8]``.
For a target cell ``i`` only source cells ``j <= i`` are used.  The diagonal
cell is triangular: its source endpoint is the moving target radius.  At the
axis the diagonal masses have polynomial continuations, so no interval
division or logarithm at ``x=s=0`` is attempted.

This module is a formal enclosure adapter.  It does not run an ancestor,
choose midpoint values, materialize ``F0``/``Lambda``, or claim a physical
field.  Coefficients are caller-supplied directed intervals.
"""
from __future__ import annotations

from collections.abc import Mapping
from typing import Any


FIELDS = (
    ("F1", "5"),
    ("Uz1", "6"),
    ("K1", "6"),
    ("P1", "4"),
    ("d_x_F1", "5"),
    ("d_x_Uz1", "6"),
)
FIELD_ROW = dict(FIELDS)
EXPECTED_C0 = {
    "4,1", "5,1", "5,2", "5,3", "5,5", "6,1", "6,2", "6,3", "6,4", "6,6"
}
EXPECTED_B1 = {"5,1", "5,2", "5,3", "6,2", "6,3", "6,4"}
ORDERS = ("0", "1", "2")
MATRIX_ORDERS = ("0", "1")

# The profile computes W1 ordinary Z0 and Z1.  W0 Z2 is an input because the
# B1*W0_ZZ term is needed for W1_Z1; no matrix order 2 or W0 Z3 is touched.
PRODUCT_RULE_Z0 = (
    ("C0", "0", "0", 1, "C0 W0"),
    ("B1", "0", "1", 1, "B1 W0_Z"),
)
PRODUCT_RULE_Z1 = (
    ("C0", "1", "0", 1, "C0_Z W0"),
    ("C0", "0", "1", 1, "C0 W0_Z"),
    ("B1", "1", "1", 1, "B1_Z W0_Z"),
    ("B1", "0", "2", 1, "B1 W0_ZZ"),
)
PRODUCT_RULE = {"Z0": PRODUCT_RULE_Z0, "Z1": PRODUCT_RULE_Z1}


def _lookup(mapping: Mapping[Any, Any], key: Any) -> Any:
    """Read an integer/string keyed mapping without changing its shape."""

    if key in mapping:
        return mapping[key]
    if isinstance(key, str) and key.isdigit():
        alternate = int(key)
    else:
        alternate = str(key)
    if alternate in mapping:
        return mapping[alternate]
    raise KeyError(key)


def _endpoint_bounds(c: Any, value: Any) -> tuple[Any, Any]:
    """Return comparable zero-width intervals at both raw endpoints."""

    value = _interval(c, value)
    lo, hi = value._mpi_
    return c.make_mpf((lo, lo)), c.make_mpf((hi, hi))


def _interval(c: Any, value: Any) -> Any:
    """Coerce a scalar or interval into the caller's interval context."""

    if hasattr(value, "_mpi_"):
        return value
    if hasattr(value, "a") and hasattr(value, "b"):
        return c.mpf([value.a, value.b])
    return c.mpf(value)


def _exact_cell(c: Any, index: int, cells: int = 8) -> Any:
    left = c.mpf(index) / cells
    right = c.mpf(index + 1) / cells
    return c.make_mpf((left._mpi_[0], right._mpi_[1]))


def _same_interval(left: Any, right: Any) -> bool:
    try:
        return left._mpi_ == right._mpi_
    except AttributeError:
        return left == right


def _degenerate_endpoint(c: Any, value: Any, side: int) -> Any:
    """Return a zero-width interval at one raw directed endpoint."""

    value = _interval(c, value)
    raw = value._mpi_[side]
    return c.make_mpf((raw, raw))


def _nonnegative(c: Any, value: Any, label: str) -> Any:
    """Intersect a directed mass enclosure with its analytic nonnegative range."""

    value = _interval(c, value)
    raw_lo, raw_hi = value._mpi_
    lo = c.make_mpf((raw_lo, raw_lo))
    hi = c.make_mpf((raw_hi, raw_hi))
    if hi < 0:
        raise ValueError(f"{label} has a negative upper bound")
    if lo < 0:
        raw_lo = c.mpf(0)._mpi_[0]
    return c.make_mpf((raw_lo, raw_hi))


def _positive(c: Any, value: Any, label: str) -> None:
    lo, hi = _endpoint_bounds(c, value)
    if hi <= 0:
        raise ValueError(f"{label} must be positive")
    if lo <= 0:
        raise ValueError(f"{label} interval contains the axis")


def _primitive_uz(c: Any, target_t: Any, source_t: Any) -> Any:
    """Primitive of ``2*s^2 log(T/s)+s^2`` in normalized radius.

    The value at ``s=0`` is exactly zero.  Keeping this branch explicit is
    important: interval division ``T/0`` and ``log(T/0)`` are not valid
    enclosures of the integrable axis endpoint.
    """

    source_t = _interval(c, source_t)
    target_t = _interval(c, target_t)
    a, _ = _endpoint_bounds(c, source_t)
    if a == 0:
        return c.mpf(0)
    _positive(c, source_t, "positive source endpoint")
    _positive(c, target_t, "positive target radius")
    return 2 * source_t * source_t * c.log(target_t / source_t) + source_t * source_t


def cell_kernel_masses(
    c: Any,
    target_t: Any,
    source_left_t: Any,
    source_right_t: Any,
    *,
    triangular: bool = False,
) -> dict[str, Any]:
    """Return exact directed Green masses for one normalized source cell.

    ``target_t`` is the target radius interval ``x/2``.  The source interval
    is ``[source_left_t, source_right_t]``.  When ``triangular`` is true, the
    source right endpoint is the moving target radius; this is the diagonal
    source/target cell.  Returned masses correspond to the K0 kernels
    ``F1: s(1-s^2/x^2)/2``, ``Uz1: s log(x/s)``, ``K1=-F1``, ``P1: 1``,
    ``d_x_F1: s^3/x^3``, and ``d_x_Uz1: s/x``.

    The first triangular cell uses the analytic axis continuation
    ``(T^2/2,T^2,-T^2/2,2T,T/2,T)``.
    """

    T = _interval(c, target_t)
    a = _interval(c, source_left_t)
    b = T if triangular else _interval(c, source_right_t)
    Tlo, Thi = _endpoint_bounds(c, T)
    alo, ahi = _endpoint_bounds(c, a)
    blo, bhi = _endpoint_bounds(c, b)
    if alo < 0 or ahi < alo or blo < alo:
        raise ValueError("invalid nonnegative source interval")
    if not triangular and bhi > Thi:
        raise ValueError("source cell extends beyond target radius")
    if triangular and not _same_interval(b, T):
        raise ValueError("triangular source endpoint must be the target interval")

    # At x=s=0, evaluate removable limits as polynomials.  This branch also
    # covers the complete first target cell [0,1/8], not only x=0.
    if triangular and alo == 0 and Tlo == 0:
        f1 = T * T / 2
        uz1 = T * T
        p1 = 2 * T
        dxf1 = T / 2
        dxuz1 = T
    else:
        _positive(c, T, "target radius")
        T2 = T * T
        a2 = a * a
        b2 = b * b
        a4 = a2 * a2
        b4 = b2 * b2
        f1 = (b2 - a2) - (b4 - a4) / (2 * T2)
        uz1 = _primitive_uz(c, T, b) - _primitive_uz(c, T, a)
        p1 = 2 * (b - a)
        dxf1 = (b4 - a4) / (2 * T * T2)
        dxuz1 = (b2 - a2) / T

    f1_nonnegative = _nonnegative(c, f1, "F1 Green mass")
    return {
        "F1": f1_nonnegative,
        "Uz1": _nonnegative(c, uz1, "Uz1 Green mass"),
        "K1": -f1_nonnegative,
        "P1": _nonnegative(c, p1, "P1 Green mass"),
        "d_x_F1": _nonnegative(c, dxf1, "d_x_F1 Green mass"),
        "d_x_Uz1": _nonnegative(c, dxuz1, "d_x_Uz1 Green mass"),
    }


def target_radius_kernel_masses(
    c: Any,
    target_x: Any,
    source_left_x: Any,
    source_right_x: Any,
    *,
    triangular: bool = False,
) -> dict[str, Any]:
    """Physical-radius spelling of :func:`cell_kernel_masses`."""

    return cell_kernel_masses(
        c,
        _interval(c, target_x) / 2,
        _interval(c, source_left_x) / 2,
        _interval(c, source_right_x) / 2,
        triangular=triangular,
    )


green_cell_kernel_masses = cell_kernel_masses
kernel_masses = cell_kernel_masses


def _cell_t(c: Any, index: int, cell: Mapping[str, Any] | None = None) -> Any:
    expected = _exact_cell(c, index)
    if cell is None or "t" not in cell:
        return expected
    actual = _interval(c, cell["t"])
    if not _same_interval(actual, expected):
        raise ValueError(f"cell {index} has a noncanonical t interval")
    if "rho" in cell and not _same_interval(_interval(c, cell["rho"]), 4 * actual * actual):
        raise ValueError(f"cell {index} has inconsistent rho=4*t^2")
    return actual


def _formal(term: Mapping[str, Any]) -> tuple[int, int, Any]:
    try:
        return int(term["F0_power"]), int(term["Lambda_power"]), term["coefficient"]
    except (KeyError, TypeError) as exc:
        raise ValueError("every source term needs F0_power, Lambda_power and coefficient") from exc


def _scaled(term: Mapping[str, Any], factor: Any, **extra: Any) -> dict[str, Any]:
    f0, lam, coefficient = _formal(term)
    result = {"F0_power": f0, "Lambda_power": lam, "coefficient": coefficient * factor}
    result.update(extra)
    return result


def _matrix_sector(cell: Mapping[str, Any], name: str, order: str, key: str) -> list[Mapping[str, Any]]:
    sector = _lookup(cell, name)
    values = _lookup(_lookup(sector, order), key)
    if isinstance(values, Mapping):
        return [values]
    if not isinstance(values, (list, tuple)):
        raise TypeError("matrix source sectors must be a list")
    return list(values)


def _validate_matrix(cell: Mapping[str, Any], label: str) -> None:
    expected = EXPECTED_C0 if label == "C0" else EXPECTED_B1
    try:
        sector = _lookup(cell, label)
    except (KeyError, TypeError) as exc:
        raise ValueError(f"matrix cell is missing {label}") from exc
    for order in MATRIX_ORDERS:
        entries = _lookup(sector, order)
        if set(str(key) for key in entries) != expected:
            raise ValueError(f"{label} support changed at order {order}")
        for key in expected:
            for term in _matrix_sector(cell, label, order, key):
                _formal(term)


def _product(
    matrix_term: Mapping[str, Any],
    w0_term: Mapping[str, Any],
    multiplier: int,
    origin: str,
    *,
    source_index: int,
    target_index: int,
    matrix_name: str,
    matrix_order: str,
    matrix_key: str,
    matrix_term_index: int,
    w0_order: str,
    w0_term_index: int,
) -> dict[str, Any]:
    mf0, mlam, mcoefficient = _formal(matrix_term)
    wf0, wlam, wcoefficient = _formal(w0_term)
    matrix_id = matrix_term.get("source_cell", matrix_term.get("source_id"))
    w0_id = w0_term.get("source_cell", w0_term.get("W0_source_cell"))
    return {
        "F0_power": mf0 + wf0,
        "Lambda_power": mlam + wlam,
        "coefficient": multiplier * mcoefficient * wcoefficient,
        "origin": origin,
        "leibniz_multiplier": multiplier,
        "row": int(matrix_key.split(",")[0]),
        "column": int(matrix_key.split(",")[1]),
        "matrix_name": matrix_name,
        "matrix_order": matrix_order,
        "matrix_key": matrix_key,
        "matrix_term_index": matrix_term_index,
        "matrix_source_cell": matrix_id,
        "matrix_source_cell_index": source_index,
        "matrix_source_id": matrix_term.get("source_id", matrix_id),
        "matrix_sector_identity": f"matrix/{source_index}/{matrix_name}/{matrix_order}/{matrix_key}/{matrix_term_index}",
        "W0_order": w0_order,
        "W0_term_index": w0_term_index,
        "W0_column": int(matrix_key.split(",")[1]),
        "W0_column_index": int(matrix_key.split(",")[1]),
        "W0_source_cell": w0_id,
        "W0_source_cell_index": w0_term.get("W0_source_cell_index", source_index),
        "W0_saved_cell_index": w0_term.get("W0_saved_cell_index", source_index),
        "W0_provenance": "saved_outer_cells",
        "W0_source_id": w0_term.get("source_id", w0_id),
        "W0_sector_identity": f"W0/{source_index}/{w0_order}/{int(matrix_key.split(',')[1])}/{w0_term_index}",
        "source_cell": source_index,
        "source_cell_id": matrix_id,
        "target_cell": target_index,
    }


def _feedback(
    matrix_cell: Mapping[str, Any],
    w0: Mapping[str, Any],
    *,
    rules: tuple[tuple[str, str, str, int, str], ...],
    source_index: int,
    target_index: int,
) -> dict[str, list[dict[str, Any]]]:
    source = {"4": [], "5": [], "6": []}
    for matrix_name, matrix_order, w0_order, multiplier, origin in rules:
        expected = EXPECTED_C0 if matrix_name == "C0" else EXPECTED_B1
        for key in expected:
            row, column = (int(value) for value in key.split(","))
            matrix_terms = _matrix_sector(matrix_cell, matrix_name, matrix_order, key)
            w0_terms = _lookup(_lookup(w0, w0_order), str(column))
            for matrix_term_index, matrix_term in enumerate(matrix_terms):
                for w0_term_index, w0_value in enumerate(w0_terms):
                    source[str(row)].append(
                        _product(
                            matrix_term,
                            w0_value,
                            multiplier,
                            origin,
                            source_index=source_index,
                            target_index=target_index,
                            matrix_name=matrix_name,
                            matrix_order=matrix_order,
                            matrix_key=key,
                            matrix_term_index=matrix_term_index,
                            w0_order=w0_order,
                            w0_term_index=w0_term_index,
                        )
                    )
    return source


def _apply_kernel(
    source: Mapping[str, list[Mapping[str, Any]]],
    masses: Mapping[str, Any],
    *,
    target_index: int,
    source_index: int,
    target_t: Any,
    source_t: Any,
) -> dict[str, list[dict[str, Any]]]:
    output = {name: [] for name, _ in FIELDS}
    for field, row in FIELDS:
        for term_index, term in enumerate(source[row]):
            output[field].append(
                _scaled(
                    term,
                    masses[field],
                    target_cell=target_index,
                    source_cell=source_index,
                    outer_cell=source_index,
                    target_radius_cell=target_index,
                    source_cell_id=term.get("source_cell_id", source_index),
                    source_term_index=term_index,
                    target_t=target_t,
                    source_t=source_t,
                    kernel_mass=masses[field],
                    kernel_mass_field=field,
                    matrix_source_cell=term.get("matrix_source_cell"),
                    matrix_source_cell_index=term.get("matrix_source_cell_index", source_index),
                    matrix_source_id=term.get("matrix_source_id"),
                    matrix_sector_identity=term.get("matrix_sector_identity"),
                    W0_source_cell=term.get("W0_source_cell"),
                    W0_source_cell_index=term.get("W0_source_cell_index"),
                    W0_saved_cell_index=term.get("W0_saved_cell_index", source_index),
                    W0_provenance=term.get("W0_provenance", "saved_outer_cells"),
                    W0_source_id=term.get("W0_source_id"),
                    W0_sector_identity=term.get("W0_sector_identity"),
                    W0_column=term.get("W0_column"),
                    W0_column_index=term.get("W0_column_index"),
                    origin=term.get("origin"),
                    leibniz_multiplier=term.get("leibniz_multiplier", 1),
                    matrix_name=term.get("matrix_name"),
                    matrix_order=term.get("matrix_order"),
                    matrix_key=term.get("matrix_key"),
                    matrix_term_index=term.get("matrix_term_index"),
                    W0_order=term.get("W0_order"),
                    W0_term_index=term.get("W0_term_index"),
                )
            )
    return output


def _w0_terms(cell: Mapping[str, Any], order: str, column: int) -> list[Mapping[str, Any]]:
    try:
        w0 = _lookup(cell, "W0_ordinary_Z")
        terms = _lookup(_lookup(w0, order), str(column))
    except (KeyError, TypeError) as exc:
        raise ValueError(f"saved outer cell is missing W0 order {order}, column {column}") from exc
    if isinstance(terms, Mapping):
        terms = [terms]
    if not isinstance(terms, (list, tuple)):
        raise TypeError("saved W0 source sectors must be a list")
    return list(terms)


def _validate_inputs(
    outer_cells: list[Mapping[str, Any]],
    matrix_cells: list[Mapping[str, Any]],
    c: Any,
) -> None:
    if len(outer_cells) != 8 or len(matrix_cells) != 8:
        raise ValueError("eight saved outer cells and eight matrix cells required")
    for index, outer_cell in enumerate(outer_cells):
        _cell_t(c, index, outer_cell)
        for order in ORDERS:
            for column in range(1, 7):
                for term in _w0_terms(outer_cell, order, column):
                    _formal(term)
    for index, matrix_cell in enumerate(matrix_cells):
        _cell_t(c, index, matrix_cell)
        _validate_matrix(matrix_cell, "C0")
        _validate_matrix(matrix_cell, "B1")


def build(
    c: Any,
    outer_cells: list[Mapping[str, Any]],
    matrix_cells: list[Mapping[str, Any]],
) -> dict[str, Any]:
    """Build dynamic target-cell and fixed ``x=2`` W1 Z0/Z1 enclosures.

    ``outer_cells`` are the already saved eight cells; their
    ``W0_ordinary_Z`` entries supply ordinary Z0/Z1/Z2 sectors directly.
    ``matrix_cells`` supplies source-cell C0/B1 sectors at orders 0 and 1.
    No W0 source or ancestor is reconstructed here.
    Every product remains unmerged and carries matrix, W0, target, and source
    cell provenance.  The two fixed-exit outputs are assembled independently
    at target x=2.
    """

    _validate_inputs(outer_cells, matrix_cells, c)
    target_cells = []
    for target_index in range(8):
        target_t = _cell_t(c, target_index, matrix_cells[target_index])
        by_source = {"Z0": {}, "Z1": {}}
        aggregate = {"Z0": {name: [] for name, _ in FIELDS},
                     "Z1": {name: [] for name, _ in FIELDS}}
        for source_index in range(target_index + 1):
            source_t = _cell_t(c, source_index, matrix_cells[source_index])
            source_left_t = _degenerate_endpoint(c, source_t, 0)
            source_right_t = (
                target_t
                if source_index == target_index
                else _degenerate_endpoint(c, source_t, 1)
            )
            masses = cell_kernel_masses(
                c,
                target_t,
                source_left_t,
                source_right_t,
                triangular=source_index == target_index,
            )
            for order_name, rules in (("Z0", PRODUCT_RULE_Z0), ("Z1", PRODUCT_RULE_Z1)):
                feedback = _feedback(
                    matrix_cells[source_index],
                    outer_cells[source_index]["W0_ordinary_Z"],
                    rules=rules,
                    source_index=source_index,
                    target_index=target_index,
                )
                terms = _apply_kernel(
                    feedback,
                    masses,
                    target_index=target_index,
                    source_index=source_index,
                    target_t=target_t,
                    source_t=source_t,
                )
                by_source[order_name][str(source_index)] = terms
                for field, values in terms.items():
                    aggregate[order_name][field].extend(values)
        target_cells.append(
            {
                "index": target_index,
                "t": target_t,
                "x": 2 * target_t,
                "rho": 4 * target_t * target_t,
                "source_indices": list(range(target_index + 1)),
                "source_j_le_target_i": True,
                "triangular_diagonal_source": target_index,
                "axis_endpoint_continuation": target_index == 0,
                "W0_ordinary_Z": outer_cells[target_index]["W0_ordinary_Z"],
                "W1_Z0": aggregate["Z0"],
                "W1_Z1": aggregate["Z1"],
                "W1_by_order": aggregate,
                "W1_by_order_by_source_cell": by_source,
            }
        )

    # Independent fixed exit: target x=2 means normalized T=1 exactly.  It
    # is not obtained by reading target cell 7's interval as a point.
    exit_target_t = c.mpf(1)
    exit_by_source = {"Z0": {}, "Z1": {}}
    exit_aggregate = {"Z0": {name: [] for name, _ in FIELDS},
                      "Z1": {name: [] for name, _ in FIELDS}}
    for source_index in range(8):
        source_t = _cell_t(c, source_index, matrix_cells[source_index])
        source_left_t = _degenerate_endpoint(c, source_t, 0)
        source_right_t = _degenerate_endpoint(c, source_t, 1)
        masses = cell_kernel_masses(
            c, exit_target_t, source_left_t, source_right_t, triangular=False
        )
        for order_name, rules in (("Z0", PRODUCT_RULE_Z0), ("Z1", PRODUCT_RULE_Z1)):
            feedback = _feedback(
                matrix_cells[source_index],
                outer_cells[source_index]["W0_ordinary_Z"],
                rules=rules,
                source_index=source_index,
                target_index=-1,
            )
            terms = _apply_kernel(
                feedback,
                masses,
                target_index=-1,
                source_index=source_index,
                target_t=exit_target_t,
                source_t=source_t,
            )
            exit_by_source[order_name][str(source_index)] = terms
            for field, values in terms.items():
                exit_aggregate[order_name][field].extend(values)

    return {
        "outer_cells": target_cells,
        "target_cells": target_cells,
        "saved_outer_cells": outer_cells,
        "W0_ordinary_Z_by_source_cell": outer_cells,
        "saved_outer_cells_reused": True,
        "W0_reconstructed": False,
        "W0_ancestor_rerun": False,
        "W0_reweighted": False,
        "no_source_cells_collapsed": True,
        "W1_Z0_target_cells": {
            str(cell["index"]): cell["W1_Z0"] for cell in target_cells
        },
        "W1_Z1_target_cells": {
            str(cell["index"]): cell["W1_Z1"] for cell in target_cells
        },
        "W1_target_cells_by_order_by_source": {
            str(cell["index"]): cell["W1_by_order_by_source_cell"]
            for cell in target_cells
        },
        "W1_Z0_core_exit": exit_aggregate["Z0"],
        "W1_Z1_core_exit": exit_aggregate["Z1"],
        "W1_Z0_exit_x2": exit_aggregate["Z0"],
        "W1_Z1_exit_x2": exit_aggregate["Z1"],
        "W1_exit_x2_by_source_cell": exit_by_source,
        "exit_x2_source_enclosures": exit_by_source,
        "exit_target": {
            "x": c.mpf(2),
            "t": exit_target_t,
            "rho": c.mpf(4),
            "fixed": True,
            "not_last_target_cell_interval": True,
        },
        "ordinary_Z_orders": [0, 1],
        "input_W0_ordinary_Z_orders": [0, 1, 2],
        "output_Z_order": 1,
        "coordinate": "rho=4 target exit and rho=4*t^2 target cells, real Z in [-1,1]",
        "normalized_radius": "t=x/2",
        "green_kernel_mass_definition": "exact directed source-cell integral; diagonal source cell is triangular",
        "green_kernel_masses_exact": True,
        "target_radius_kernel_not_last_interval_proxy": True,
        "source_j_le_target_i": True,
        "axis_endpoint_continuation": True,
        "axis_log_division_avoided": True,
        "source_ids_retained": True,
        "formal_powers_retained": True,
        "signed_coefficients_retained": True,
        "ordinary_Z0_product_rule": PRODUCT_RULE_Z0,
        "ordinary_Z1_product_rule": PRODUCT_RULE_Z1,
        "principal_kernel": "K0 only; no Kchi or higher swirl words",
        "source_exponentials_materialized": False,
        "first_collar_added": False,
        "allword_tail_added": False,
        "m15_computed": False,
        "physical_field_claimed": False,
    }


__all__ = [
    "FIELDS",
    "EXPECTED_C0",
    "EXPECTED_B1",
    "PRODUCT_RULE_Z0",
    "PRODUCT_RULE_Z1",
    "PRODUCT_RULE",
    "cell_kernel_masses",
    "target_radius_kernel_masses",
    "green_cell_kernel_masses",
    "kernel_masses",
    "build",
]
