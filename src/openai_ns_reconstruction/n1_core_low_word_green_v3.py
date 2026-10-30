"""Green assembly for the low first-word core exit.

The eight ``outer_cells`` are saved, already integrated ``W0`` enclosures.
This module only multiplies those enclosures by the caller-supplied zeroth
and first ``C0``/``B1`` matrix sectors and applies the exact six leaf kernels.
It deliberately does not reconstruct or reweight ``W0`` and it does not form
the second ``Z`` derivative of ``W1``.
"""
from __future__ import annotations


FIELDS = (
    ("F1", "5"),
    ("Uz1", "6"),
    ("K1", "6"),
    ("P1", "4"),
    ("d_x_F1", "5"),
    ("d_x_Uz1", "6"),
)

EXPECTED_C0 = {
    "4,1",
    "5,1",
    "5,2",
    "5,3",
    "5,5",
    "6,1",
    "6,2",
    "6,3",
    "6,4",
    "6,6",
}
EXPECTED_B1 = {"5,1", "5,2", "5,3", "6,2", "6,3", "6,4"}
_ORDERS = ("0", "1", "2")
_MATRIX_ORDERS = ("0", "1")


def _lookup(mapping, key):
    """Read a required string-keyed order or column without changing input."""
    if key in mapping:
        return mapping[key]
    alternate = int(key) if isinstance(key, str) and key.isdigit() else str(key)
    if alternate in mapping:
        return mapping[alternate]
    raise KeyError(key)


def _scale(term, factor, **extra):
    """Scale one formal term while retaining only its formal powers and data."""
    result = {
        "F0_power": term["F0_power"],
        "Lambda_power": term["Lambda_power"],
        "coefficient": term["coefficient"] * factor,
    }
    result.update(extra)
    return result


def _product(
    matrix_term,
    w0_term,
    multiplier,
    origin,
    row,
    column,
    outer_cell,
    matrix_key,
    matrix_order,
    matrix_term_index,
    w0_term_index,
    w0_order,
):
    """Form one unmerged Leibniz product and retain both source identities."""
    result = {
        "F0_power": matrix_term["F0_power"] + w0_term["F0_power"],
        "Lambda_power": matrix_term["Lambda_power"] + w0_term["Lambda_power"],
        "coefficient": (
            multiplier
            * matrix_term["coefficient"]
            * w0_term["coefficient"]
        ),
        "origin": origin,
        "row": row,
        "column": column,
        "matrix_row": row,
        "matrix_column": column,
        "matrix_key": matrix_key,
        "matrix_order": matrix_order,
        "matrix_term_index": matrix_term_index,
        "W0_term_index": w0_term_index,
        "W0_order": w0_order,
        "outer_cell": outer_cell,
        "matrix_source_cell": matrix_term.get("source_cell", outer_cell),
        "W0_source_cell": w0_term.get("source_cell"),
    }
    # Keep a matrix-sector source_cell verbatim when one was supplied.  The
    # explicit matrix_source_cell field above still distinguishes it from the
    # saved W0 source cell.
    if "source_cell" in matrix_term:
        result["source_cell"] = matrix_term["source_cell"]
    return result


def _validate_interval_cell(c, index, cell, label):
    """Validate one directed outer cell using the old leaf convention."""
    try:
        t = cell["t"]
        rho = cell["rho"]
        left = c.mpf(index) / 8
        right = c.mpf(index + 1) / 8
        expected_rho = 4 * t * t
        if t.a != left.a or t.b != right.b or rho._mpi_ != expected_rho._mpi_:
            raise ValueError(f"{label} is not the directed rho=4*t^2 outer cell")
    except (AttributeError, KeyError, TypeError) as exc:
        raise ValueError(f"{label} is missing a valid t/rho interval") from exc
    return t, rho


def _validate_w0(w0, label):
    if not isinstance(w0, dict):
        raise ValueError(f"{label} W0_ordinary_Z must be a mapping")
    for order in _ORDERS:
        try:
            columns = _lookup(w0, order)
        except KeyError as exc:
            raise ValueError(f"{label} is missing W0 order {order}") from exc
        if not isinstance(columns, dict):
            raise ValueError(f"{label} W0 order {order} must be a mapping")
        if {str(key) for key in columns} != {str(col) for col in range(1, 7)}:
            raise ValueError(f"{label} W0 order {order} must contain columns 1 through 6")


def _validate_matrix(cell, name):
    try:
        sector = cell[name]
    except (KeyError, TypeError) as exc:
        raise ValueError(f"matrix cell is missing {name}") from exc
    if not isinstance(sector, dict):
        raise ValueError(f"matrix cell {name} must be a mapping")
    expected = EXPECTED_C0 if name == "C0" else EXPECTED_B1
    for order in _MATRIX_ORDERS:
        try:
            entries = _lookup(sector, order)
        except KeyError as exc:
            raise ValueError(f"matrix cell {name} is missing order {order}") from exc
        if not isinstance(entries, dict):
            raise ValueError(f"matrix cell {name} order {order} must be a mapping")
        if {str(key) for key in entries} != expected:
            raise ValueError(f"raw dynamic {name} support changed")


def _same_interval(left, right):
    try:
        return left._mpi_ == right._mpi_
    except AttributeError:
        return left == right


def _append_products(target, matrix_cell, w0, specs, outer_index):
    """Append all products for one requested W1 Z order, without merging."""
    for matrix_name, matrix_order, w0_order, multiplier, origin in specs:
        entries = _lookup(matrix_cell[matrix_name], matrix_order)
        for key, matrix_terms in entries.items():
            row, column = (int(value) for value in str(key).split(","))
            if row not in (4, 5, 6) or column not in range(1, 7):
                raise ValueError("invalid raw matrix coordinate")
            w0_terms = _lookup(_lookup(w0, w0_order), str(column))
            for matrix_term_index, matrix_term in enumerate(matrix_terms):
                for w0_term_index, w0_term in enumerate(w0_terms):
                    target[str(row)].append(
                        _product(
                            matrix_term,
                            w0_term,
                            multiplier,
                            origin,
                            row,
                            column,
                            outer_index,
                            str(key),
                            matrix_order,
                            matrix_term_index,
                            w0_term_index,
                            w0_order,
                        )
                    )


def build(c, outer_cells, matrix_cells):
    """Build the local first-word ``W1_Z0`` and ``W1_Z1`` Green terms.

    ``outer_cells`` must be the eight saved outer cells from the ordinary
    ``W0_ordinary_Z`` construction.  Their order and column keys are read as
    strings when present; every saved term is used as-is, including its
    integrated coefficient, formal powers, and ``source_cell``.  Each matching
    entry in ``matrix_cells`` supplies only the ``C0``/``B1`` matrix sectors of
    orders 0 and 1.
    """
    if len(outer_cells) != 8 or len(matrix_cells) != 8:
        raise ValueError("eight saved outer cells and eight matrix cells required")

    output_z0 = {name: [] for name, _ in FIELDS}
    output_z1 = {name: [] for name, _ in FIELDS}

    z0_specs = (
        ("C0", "0", "0", 1, "C0 W0"),
        ("B1", "0", "1", 1, "B1 W0_Z"),
    )
    z1_specs = (
        ("C0", "1", "0", 1, "C0_Z W0"),
        ("C0", "0", "1", 1, "C0 W0_Z"),
        ("B1", "1", "1", 1, "B1_Z W0_Z"),
        ("B1", "0", "2", 1, "B1 W0_ZZ"),
    )

    for index, (outer, matrix_cell) in enumerate(zip(outer_cells, matrix_cells)):
        outer_t, outer_rho = _validate_interval_cell(c, index, outer, "saved outer cell")
        matrix_t, matrix_rho = _validate_interval_cell(c, index, matrix_cell, "matrix cell")
        if not _same_interval(outer_t, matrix_t) or not _same_interval(outer_rho, matrix_rho):
            raise ValueError("saved outer cell and matrix cell coordinates differ")

        try:
            w0 = outer["W0_ordinary_Z"]
        except (KeyError, TypeError) as exc:
            raise ValueError("saved outer cell is missing W0_ordinary_Z") from exc
        _validate_w0(w0, "saved outer cell")
        _validate_matrix(matrix_cell, "C0")
        _validate_matrix(matrix_cell, "B1")

        source_z0 = {"4": [], "5": [], "6": []}
        source_z1 = {"4": [], "5": [], "6": []}
        _append_products(source_z0, matrix_cell, w0, z0_specs, index)
        _append_products(source_z1, matrix_cell, w0, z1_specs, index)

        # These are the exact directed leaf kernels used by the prior Green
        # assembly.  The final factor 1/4 is part of the leaf convention.
        tlog = c.mpf([0, "0.5"]) if index == 0 else outer_t * c.log(1 / outer_t)
        K = outer_t * (1 - outer_t * outer_t)
        kernels = {
            "F1": K,
            "Uz1": 2 * tlog,
            "K1": -K,
            "P1": c.mpf(1),
            "d_x_F1": outer_t**3,
            "d_x_Uz1": outer_t,
        }
        for field, row in FIELDS:
            for term in source_z0[row]:
                output_z0[field].append(
                    _scale(
                        term,
                        kernels[field] / 4,
                        outer_cell=index,
                        source_row=row,
                        row=term["row"],
                        origin=term["origin"],
                        column=term["column"],
                        matrix_row=term["matrix_row"],
                        matrix_column=term["matrix_column"],
                        matrix_key=term["matrix_key"],
                        matrix_order=term["matrix_order"],
                        matrix_term_index=term["matrix_term_index"],
                        W0_term_index=term["W0_term_index"],
                        W0_order=term["W0_order"],
                        matrix_source_cell=term["matrix_source_cell"],
                        W0_source_cell=term["W0_source_cell"],
                        source_cell=term.get("source_cell"),
                    )
                )
            for term in source_z1[row]:
                output_z1[field].append(
                    _scale(
                        term,
                        kernels[field] / 4,
                        outer_cell=index,
                        source_row=row,
                        row=term["row"],
                        origin=term["origin"],
                        column=term["column"],
                        matrix_row=term["matrix_row"],
                        matrix_column=term["matrix_column"],
                        matrix_key=term["matrix_key"],
                        matrix_order=term["matrix_order"],
                        matrix_term_index=term["matrix_term_index"],
                        W0_term_index=term["W0_term_index"],
                        W0_order=term["W0_order"],
                        matrix_source_cell=term["matrix_source_cell"],
                        W0_source_cell=term["W0_source_cell"],
                        source_cell=term.get("source_cell"),
                    )
                )

    return {
        "outer_cells": outer_cells,
        "W1_Z0_core_exit": output_z0,
        "W1_Z1_core_exit": output_z1,
        "ordinary_Z_orders": [0, 1],
        "input_W0_ordinary_Z_orders": [0, 1, 2],
        "output_Z_order": 1,
        "coordinate": "rho=4, core exit, real Z in [-1,1]",
        "W0_reweighted": False,
        "no_source_cells_collapsed": True,
        "axis_endpoint_continuation": True,
        "source_exponentials_materialized": False,
        "first_collar_added": False,
        "allword_tail_added": False,
        "m15_computed": False,
        "allword_included": False,
        "m15_included": False,
        "W1_Z2_computed": False,
    }
