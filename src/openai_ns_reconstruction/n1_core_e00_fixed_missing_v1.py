"""Supplementary fixed-endpoint E00 Green sectors.

This leaf fills only the four E00 components which are not supplied by the
saved fixed ``x=2`` ``F1``/``Uz1`` packet.  It consumes caller-owned saved
``W0_ordinary_Z`` source cells and the corresponding ``C0``/``B1`` matrix
cells.  No ancestor is rebuilt and no source coefficient is sampled or
merged.

The normalized source coordinate is ``t=x/2``.  For source cell ``j`` the
physical source measure is applied exactly once through ``dx=2*dt``.  The
four masses are the polynomial fixed-target moments

``K1: -2 integral(t*(1-t**2) dt)``
``P1:  2 integral(1 dt)``
``d_x_F1: 2 integral(t**3 dt)``
``d_x_Uz1: 2 integral(t dt)``.

Only the ordinary product rules ``C0 W0 + B1 W0_Z`` and
``C0_Z W0 + C0 W0_Z + B1_Z W0_Z + B1 W0_ZZ`` are formed.  Products remain
unmerged and carry both matrix and saved-W0 identities and paths.  ``P1`` is
the row-4 state component; it is not a pressure/m15 contribution.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any


SCHEMA = "n1_core_e00_fixed_missing_v1"
SOURCE_SCHEMA = "n1_core_low_word_v3"
CELL_COUNT = 8
ORDERS = (0, 1)
INPUT_W0_ORDERS = (0, 1, 2)
ROWS = ("4", "5", "6")
FIELDS = ("K1", "P1", "d_x_F1", "d_x_Uz1")
FIELD_ROWS = {
    "K1": "6",
    "P1": "4",
    "d_x_F1": "5",
    "d_x_Uz1": "6",
}

C0_KEYS = (
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
)
B1_KEYS = ("5,1", "5,2", "5,3", "6,2", "6,3", "6,4")
EXPECTED_C0 = set(C0_KEYS)
EXPECTED_B1 = set(B1_KEYS)

PRODUCT_RULE_Z0 = (
    ("C0", "0", "0", "C0 W0"),
    ("B1", "0", "1", "B1 W0_Z"),
)
PRODUCT_RULE_Z1 = (
    ("C0", "1", "0", "C0_Z W0"),
    ("C0", "0", "1", "C0 W0_Z"),
    ("B1", "1", "1", "B1_Z W0_Z"),
    ("B1", "0", "2", "B1 W0_ZZ"),
)
PRODUCT_RULE = {"Z0": PRODUCT_RULE_Z0, "Z1": PRODUCT_RULE_Z1}

__all__ = [
    "SCHEMA",
    "SOURCE_SCHEMA",
    "CELL_COUNT",
    "ORDERS",
    "INPUT_W0_ORDERS",
    "ROWS",
    "FIELDS",
    "FIELD_ROWS",
    "C0_KEYS",
    "B1_KEYS",
    "EXPECTED_C0",
    "EXPECTED_B1",
    "PRODUCT_RULE_Z0",
    "PRODUCT_RULE_Z1",
    "PRODUCT_RULE",
    "fixed_cell_kernel_masses",
    "cell_kernel_masses",
    "endpoint_kernel_masses",
    "build",
    "build_fixed_exit",
    "build_x2",
]


def _lookup(mapping: Mapping[Any, Any], key: Any) -> Any:
    """Read a string/integer keyed mapping without changing its shape."""

    if key in mapping:
        return mapping[key]
    if isinstance(key, str) and key.isdigit():
        alternate: Any = int(key)
    else:
        alternate = str(key)
    if alternate in mapping:
        return mapping[alternate]
    raise KeyError(key)


def _interval(c: Any, value: Any, label: str = "interval") -> Any:
    """Coerce a scalar or directed interval without taking a midpoint."""

    if hasattr(value, "_mpi_") and hasattr(value, "a") and hasattr(value, "b"):
        return value
    if hasattr(value, "a") and hasattr(value, "b"):
        lo = value.a._mpi_[0] if hasattr(value.a, "_mpi_") else value.a
        hi = value.b._mpi_[1] if hasattr(value.b, "_mpi_") else value.b
        return c.make_mpf((lo, hi))
    if isinstance(value, Mapping):
        if "lower_exact_mpf_tuple" in value and "upper_exact_mpf_tuple" in value:
            return c.make_mpf(
                (
                    tuple(value["lower_exact_mpf_tuple"]),
                    tuple(value["upper_exact_mpf_tuple"]),
                )
            )
        if "lower" in value and "upper" in value:
            return c.mpf([value["lower"], value["upper"]])
        if "lo" in value and "hi" in value:
            return c.mpf([value["lo"], value["hi"]])
        raise TypeError(f"unsupported {label} mapping")
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        if len(value) == 2:
            return c.mpf([value[0], value[1]])
    try:
        return c.mpf(value)
    except (TypeError, ValueError) as exc:
        raise TypeError(f"unsupported {label}") from exc


def _same_interval(left: Any, right: Any) -> bool:
    try:
        return left._mpi_ == right._mpi_
    except AttributeError:
        return left == right


def _canonical_t(c: Any, index: int) -> Any:
    left = c.mpf(index) / CELL_COUNT
    right = c.mpf(index + 1) / CELL_COUNT
    return c.make_mpf((left._mpi_[0], right._mpi_[1]))


def _endpoint(c: Any, value: Any, upper: bool) -> Any:
    value = _interval(c, value, "cell t")
    raw = value._mpi_[1 if upper else 0]
    return c.make_mpf((raw, raw))


def _cell_t(c: Any, cell: Mapping[str, Any], index: int, label: str) -> Any:
    expected = _canonical_t(c, index)
    if not isinstance(cell, Mapping):
        raise TypeError(f"{label} {index} must be a mapping")
    candidate = cell.get("t", cell.get("normalized_t", cell.get("radius_t")))
    if candidate is None:
        raise ValueError(f"{label} {index} is missing t")
    actual = _interval(c, candidate, f"{label} t")
    if not _same_interval(actual, expected):
        raise ValueError(f"{label} {index} is not canonical t=[i/8,(i+1)/8]")
    if "rho" in cell and not _same_interval(
        _interval(c, cell["rho"], f"{label} rho"), 4 * actual * actual
    ):
        raise ValueError(f"{label} {index} has inconsistent rho=4*t^2")
    if "x" in cell and not _same_interval(
        _interval(c, cell["x"], f"{label} x"), 2 * actual
    ):
        raise ValueError(f"{label} {index} has inconsistent x=2*t")
    return actual


def _as_cells(value: Any, label: str) -> list[Mapping[str, Any]]:
    if isinstance(value, Mapping):
        for wrapper in ("cells", "source_cells", "outer_cells", "matrix_cells"):
            if wrapper in value:
                value = value[wrapper]
                break
    if isinstance(value, Mapping):
        try:
            cells = [_lookup(value, index) for index in range(CELL_COUNT)]
        except KeyError as exc:
            raise TypeError(f"{label} must contain eight canonical cells") from exc
    elif isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        cells = list(value)
    else:
        raise TypeError(f"{label} must contain eight canonical cells")
    if len(cells) != CELL_COUNT:
        raise ValueError(f"eight {label}s are required")
    if not all(isinstance(cell, Mapping) for cell in cells):
        raise TypeError(f"every {label} must be a mapping")
    return cells


def _formal(term: Mapping[str, Any], label: str) -> tuple[int, int, Any]:
    if not isinstance(term, Mapping):
        raise TypeError(f"{label} must be a formal term mapping")
    f0 = term.get("F0_power")
    lam = term.get("Lambda_power")
    if isinstance(f0, bool) or not isinstance(f0, int):
        raise ValueError(f"{label} is missing integer F0_power")
    if isinstance(lam, bool) or not isinstance(lam, int):
        raise ValueError(f"{label} is missing integer Lambda_power")
    if "coefficient" not in term:
        raise ValueError(f"{label} is missing coefficient")
    return f0, lam, term["coefficient"]


def _terms(raw: Any, label: str) -> list[Mapping[str, Any]]:
    if isinstance(raw, Mapping):
        values = [raw]
    elif isinstance(raw, Sequence) and not isinstance(raw, (str, bytes, bytearray)):
        values = list(raw)
    else:
        raise TypeError(f"{label} must be a term list")
    result: list[Mapping[str, Any]] = []
    for index, value in enumerate(values):
        if not isinstance(value, Mapping):
            raise TypeError(f"{label}[{index}] must be a term mapping")
        result.append(value)
    return result


def _w0_terms(cell: Mapping[str, Any], order: int, column: int) -> list[Mapping[str, Any]]:
    try:
        grid = _lookup(cell, "W0_ordinary_Z")
        order_grid = _lookup(grid, order)
        raw = _lookup(order_grid, str(column))
    except (KeyError, TypeError) as exc:
        raise ValueError(f"saved outer cell is missing W0 order {order}, column {column}") from exc
    values = _terms(raw, f"W0 order {order}, column {column}")
    if not values:
        raise ValueError(f"W0 order {order}, column {column} is empty")
    return values


def _matrix_terms(
    cell: Mapping[str, Any], name: str, order: str, key: str
) -> list[Mapping[str, Any]]:
    try:
        sector = _lookup(cell, name)
        order_sector = _lookup(sector, order)
        raw = _lookup(order_sector, key)
    except (KeyError, TypeError) as exc:
        raise ValueError(f"matrix cell is missing {name} {order} {key}") from exc
    values = _terms(raw, f"matrix {name} {order} {key}")
    if not values:
        raise ValueError(f"matrix {name} {order} {key} is empty")
    return values


def _validate_inputs(
    c: Any,
    outer_cells: list[Mapping[str, Any]],
    matrix_cells: list[Mapping[str, Any]],
) -> None:
    if len(outer_cells) != CELL_COUNT or len(matrix_cells) != CELL_COUNT:
        raise ValueError("eight saved outer cells and eight matrix cells required")
    for index, cell in enumerate(outer_cells):
        _cell_t(c, cell, index, "outer cell")
        for order in INPUT_W0_ORDERS:
            for column in range(1, 7):
                for term_index, term in enumerate(_w0_terms(cell, order, column)):
                    _formal(term, f"outer cell {index} W0 Z{order} column {column} term {term_index}")
    for index, cell in enumerate(matrix_cells):
        _cell_t(c, cell, index, "matrix cell")
        for name, expected in (("C0", EXPECTED_C0), ("B1", EXPECTED_B1)):
            try:
                sector = _lookup(cell, name)
            except (KeyError, TypeError) as exc:
                raise ValueError(f"matrix cell {index} is missing {name}") from exc
            for order in ("0", "1"):
                try:
                    order_sector = _lookup(sector, order)
                except (KeyError, TypeError) as exc:
                    raise ValueError(f"matrix cell {index} is missing {name} order {order}") from exc
                if set(str(key) for key in order_sector) != expected:
                    raise ValueError(f"{name} support changed at order {order} in matrix cell {index}")
                for key in expected:
                    for term_index, term in enumerate(_matrix_terms(cell, name, order, key)):
                        _formal(term, f"matrix cell {index} {name} {order} {key} term {term_index}")


def fixed_cell_kernel_masses(c: Any, source_index: int) -> dict[str, Any]:
    """Return exact directed polynomial moments for one fixed source cell."""

    if not isinstance(source_index, int) or isinstance(source_index, bool):
        raise TypeError("source_index must be an integer")
    if source_index < 0 or source_index >= CELL_COUNT:
        raise ValueError("source_index must be in range(8)")
    t = _canonical_t(c, source_index)
    a = _endpoint(c, t, False)
    b = _endpoint(c, t, True)
    a2 = a * a
    b2 = b * b
    a4 = a2 * a2
    b4 = b2 * b2
    masses = {
        "K1": -(b2 - a2) + (b4 - a4) / 2,
        "P1": 2 * (b - a),
        "d_x_F1": (b4 - a4) / 2,
        "d_x_Uz1": b2 - a2,
    }
    masses.update(
        {
            "target_t": c.mpf(1),
            "target_x": c.mpf(2),
            "target_rho": c.mpf(4),
            "source_left_t": a,
            "source_right_t": b,
            "source_left_x": 2 * a,
            "source_right_x": 2 * b,
            "physical_source_measure": "dx=2*dt",
            "physical_measure_applied_once": True,
            "principal_kernel": "G0 fixed-target E00 supplement",
            "exact_source_cell_mass": True,
            "axis_endpoint_continuation": source_index == 0,
            "negative_kernel_K1": True,
        }
    )
    return masses


cell_kernel_masses = fixed_cell_kernel_masses
endpoint_kernel_masses = fixed_cell_kernel_masses


def _sequence_values(value: Any) -> list[Any]:
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return list(value)
    return [value]


def _unique(values: Sequence[Any]) -> list[Any]:
    result: list[Any] = []
    markers: set[Any] = set()
    for value in values:
        try:
            marker = (type(value).__name__, value)
            hash(marker)
        except TypeError:
            marker = (type(value).__name__, repr(value))
        if marker not in markers:
            markers.add(marker)
            result.append(value)
    return result


def _term_ids(term: Mapping[str, Any]) -> list[Any]:
    values: list[Any] = []
    # ``source_cell`` is an index/owner label, not an identity.  Keep it in
    # the dedicated cell fields below; only explicit IDs (or origin IDs) enter
    # the combined ownership list.
    for key in ("source_ids", "source_id", "source_cell_id", "W0_source_id"):
        if key in term and term[key] is not None:
            values.extend(_sequence_values(term[key]))
    origin = term.get("origin")
    if isinstance(origin, Mapping):
        for key in ("source_ids", "source_id"):
            if key in origin and origin[key] is not None:
                values.extend(_sequence_values(origin[key]))
    return _unique(values)


def _term_paths(term: Mapping[str, Any]) -> list[Any]:
    values: list[Any] = []
    for key in ("source_paths", "source_path", "path"):
        if key not in term or term[key] is None:
            continue
        raw = term[key]
        if key == "source_paths" and isinstance(raw, Sequence) and not isinstance(raw, (str, bytes, bytearray)):
            entries = list(raw)
            if entries and all(
                isinstance(entry, Sequence) and not isinstance(entry, (str, bytes, bytearray))
                for entry in entries
            ):
                values.extend(entries)
            else:
                values.append(raw)
        else:
            values.append(raw)
    origin = term.get("origin")
    if isinstance(origin, Mapping):
        for key in ("source_paths", "source_path", "path"):
            if key not in origin or origin[key] is None:
                continue
            raw = origin[key]
            if key == "source_paths" and isinstance(raw, Sequence) and not isinstance(raw, (str, bytes, bytearray)):
                entries = list(raw)
                if entries and all(
                    isinstance(entry, Sequence) and not isinstance(entry, (str, bytes, bytearray))
                    for entry in entries
                ):
                    values.extend(entries)
                else:
                    values.append(raw)
            else:
                values.append(raw)
    return _unique(values)


def _path_id(path: Any) -> str:
    if isinstance(path, Sequence) and not isinstance(path, (str, bytes, bytearray)):
        return "path:" + "/".join(str(value) for value in path)
    return "path:" + str(path)


def _product(
    c: Any,
    matrix_term: Mapping[str, Any],
    w0_term: Mapping[str, Any],
    *,
    matrix_name: str,
    matrix_order: str,
    matrix_key: str,
    matrix_term_index: int,
    w0_order: str,
    w0_term_index: int,
    source_index: int,
    output_order: int,
    origin: str,
    outer_cell: Mapping[str, Any],
    matrix_cell: Mapping[str, Any],
) -> dict[str, Any]:
    mf0, mlam, mcoefficient = _formal(matrix_term, "matrix term")
    wf0, wlam, wcoefficient = _formal(w0_term, "W0 term")
    matrix_ids = _term_ids(matrix_term)
    w0_ids = _term_ids(w0_term)
    matrix_paths = _term_paths(matrix_term)
    w0_paths = _term_paths(w0_term)
    matrix_path = matrix_paths[0] if matrix_paths else [
        "matrix_cells", source_index, matrix_name, matrix_order, matrix_key, matrix_term_index
    ]
    w0_path = w0_paths[0] if w0_paths else [
        "outer_cells", source_index, "W0_ordinary_Z", w0_order,
        matrix_key.split(",")[1], w0_term_index,
    ]
    if not matrix_ids:
        matrix_ids = [_path_id(matrix_path)]
    if not w0_ids:
        w0_ids = [_path_id(w0_path)]
    source_ids = _unique(matrix_ids + w0_ids)
    matrix_source_id = matrix_ids[0]
    w0_source_id = w0_ids[0]
    matrix_cell_id = matrix_term.get(
        "source_cell", matrix_term.get("source_cell_id", matrix_cell.get("cell_id", source_index))
    )
    w0_cell_id = w0_term.get(
        "source_cell", w0_term.get("W0_source_cell", outer_cell.get("cell_id", source_index))
    )
    result = {
        "F0_power": mf0 + wf0,
        "Lambda_power": mlam + wlam,
        "coefficient": _interval(c, mcoefficient, "matrix coefficient")
        * _interval(c, wcoefficient, "W0 coefficient"),
        "origin": origin,
        "input_origins": {
            "matrix": matrix_term.get("origin"),
            "W0": w0_term.get("origin"),
        },
        "source_ids": source_ids,
        "source_id": source_ids[0],
        "source_paths": [matrix_path, w0_path],
        "matrix_source_id": matrix_source_id,
        "W0_source_id": w0_source_id,
        "matrix_source_path": matrix_path,
        "W0_source_path": w0_path,
        "matrix_name": matrix_name,
        "matrix_order": matrix_order,
        "matrix_key": matrix_key,
        "matrix_term_index": matrix_term_index,
        "matrix_source_cell": matrix_cell_id,
        "matrix_source_cell_index": source_index,
        "matrix_sector_identity": f"matrix/{source_index}/{matrix_name}/{matrix_order}/{matrix_key}/{matrix_term_index}",
        "W0_order": w0_order,
        "W0_term_index": w0_term_index,
        "W0_column": int(matrix_key.split(",")[1]),
        "W0_source_cell": w0_cell_id,
        "W0_source_cell_index": source_index,
        "W0_saved_cell_index": source_index,
        "W0_sector_identity": f"W0/{source_index}/{w0_order}/{matrix_key.split(',')[1]}/{w0_term_index}",
        "W0_provenance": "saved_outer_cells",
        "source_cell": source_index,
        "source_cell_index": source_index,
        "source_cell_id": w0_cell_id,
        "outer_cell": source_index,
        "outer_cell_index": source_index,
        "target_cell": -1,
        "target_cell_index": -1,
        "ordinary_Z_order": output_order,
        "output_Z_order": output_order,
        "leibniz_multiplier": 1,
        "signed_source_term": True,
        "formal_powers_retained": True,
    }
    return result


def _feedback(
    c: Any,
    outer_cell: Mapping[str, Any],
    matrix_cell: Mapping[str, Any],
    *,
    field: str,
    output_order: int,
    source_index: int,
) -> list[dict[str, Any]]:
    row = FIELD_ROWS[field]
    rules = PRODUCT_RULE_Z0 if output_order == 0 else PRODUCT_RULE_Z1
    result: list[dict[str, Any]] = []
    for matrix_name, matrix_order, w0_order, origin in rules:
        keys = C0_KEYS if matrix_name == "C0" else B1_KEYS
        for matrix_key in keys:
            if matrix_key.split(",", 1)[0] != row:
                continue
            matrix_values = _matrix_terms(matrix_cell, matrix_name, matrix_order, matrix_key)
            column = int(matrix_key.split(",", 1)[1])
            w0_values = _w0_terms(outer_cell, int(w0_order), column)
            for matrix_term_index, matrix_term in enumerate(matrix_values):
                for w0_term_index, w0_term in enumerate(w0_values):
                    result.append(
                        _product(
                            c,
                            matrix_term,
                            w0_term,
                            matrix_name=matrix_name,
                            matrix_order=matrix_order,
                            matrix_key=matrix_key,
                            matrix_term_index=matrix_term_index,
                            w0_order=w0_order,
                            w0_term_index=w0_term_index,
                            source_index=source_index,
                            output_order=output_order,
                            origin=origin,
                            outer_cell=outer_cell,
                            matrix_cell=matrix_cell,
                        )
                    )
    return result


def _empty_fields() -> dict[str, list[dict[str, Any]]]:
    return {field: [] for field in FIELDS}


def _scale(
    term: Mapping[str, Any],
    mass: Any,
    *,
    field: str,
    source_index: int,
    source_term_index: int,
) -> dict[str, Any]:
    result = dict(term)
    result["coefficient"] = term["coefficient"] * mass
    result.update(
        {
            "field": field,
            "row": int(FIELD_ROWS[field]),
            "field_row": FIELD_ROWS[field],
            "source_row": FIELD_ROWS[field],
            "source_term_index": source_term_index,
            "kernel_mass": mass,
            "kernel_mass_bare": mass,
            "kernel_mass_field": field,
            "principal_kernel": "G0 fixed-target E00 supplement",
            "negative_kernel_K1": field == "K1",
            "exact_source_cell_mass": True,
            "physical_source_measure": "dx=2*dt",
            "physical_measure_applied_once": True,
            "source_cell": source_index,
            "source_cell_index": source_index,
            "forcing_source_cell": source_index,
            "target_cell": -1,
            "target_cell_index": -1,
            "target_radius_cell": -1,
        }
    )
    return result


def build(
    c: Any,
    outer_cells: Any = None,
    matrix_cells: Any = None,
    source_family: Any = None,
    source_logs: Any = None,
) -> dict[str, Any]:
    """Build the four missing fixed-target E00 fields through ordinary Z1."""

    if not hasattr(c, "mpf") or not hasattr(c, "make_mpf"):
        raise TypeError("c must provide directed interval arithmetic")
    if outer_cells is None or matrix_cells is None:
        raise ValueError("outer_cells and matrix_cells are required")
    outer = _as_cells(outer_cells, "outer cell")
    matrices = _as_cells(matrix_cells, "matrix cell")
    _validate_inputs(c, outer, matrices)

    aggregate = {f"Z{order}": _empty_fields() for order in ORDERS}
    by_source: dict[str, Any] = {}
    for source_index, (outer_cell, matrix_cell) in enumerate(zip(outer, matrices)):
        source_t = _cell_t(c, outer_cell, source_index, "outer cell")
        matrix_t = _cell_t(c, matrix_cell, source_index, "matrix cell")
        masses = fixed_cell_kernel_masses(c, source_index)
        cell_output = {f"Z{order}": _empty_fields() for order in ORDERS}
        for output_order in ORDERS:
            for field in FIELDS:
                products = _feedback(
                    c,
                    outer_cell,
                    matrix_cell,
                    field=field,
                    output_order=output_order,
                    source_index=source_index,
                )
                for term_index, product in enumerate(products):
                    term = _scale(
                        product,
                        masses[field],
                        field=field,
                        source_index=source_index,
                        source_term_index=term_index,
                    )
                    term["source_t"] = source_t
                    term["matrix_t"] = matrix_t
                    term["target_t"] = masses["target_t"]
                    term["target_x"] = masses["target_x"]
                    term["target_rho"] = masses["target_rho"]
                    cell_output[f"Z{output_order}"][field].append(term)
                    aggregate[f"Z{output_order}"][field].append(term)
        by_source[str(source_index)] = {
            "source_index": source_index,
            "source_cell": source_index,
            "source_t": source_t,
            "source_x": 2 * source_t,
            "target_t": masses["target_t"],
            "target_x": masses["target_x"],
            "target_rho": masses["target_rho"],
            "fixed_target": True,
            "kernel_masses": masses,
            "E00": cell_output,
            "E00_by_order": cell_output,
            "E00_Z0": cell_output["Z0"],
            "E00_Z1": cell_output["Z1"],
            "fields": list(FIELDS),
            "ordinary_Z_orders": [0, 1],
            "source_ids_retained": True,
            "formal_powers_retained": True,
            "signed_coefficients_retained": True,
        }

    fixed_metadata = {
        "fields": list(FIELDS),
        "ordinary_Z_orders": [0, 1],
        "target_t": c.mpf(1),
        "target_x": c.mpf(2),
        "target_rho": c.mpf(4),
        "fixed": True,
    }
    return {
        "schema": SCHEMA,
        "source_schema": SOURCE_SCHEMA,
        "coordinate": "rho=4, x=2, t=1, real Z in [-1,1]",
        "source_family": source_family,
        "source_logs": source_logs,
        "source_cell_count": CELL_COUNT,
        "source_cells": list(range(CELL_COUNT)),
        "fields": list(FIELDS),
        "field_rows": dict(FIELD_ROWS),
        "ordinary_Z_orders": [0, 1],
        "input_W0_ordinary_Z_orders": [0, 1, 2],
        "input_matrix_orders": [0, 1],
        "E00_fixed": fixed_metadata,
        "E00_fixed_by_source_cell": by_source,
        "E00_fixed_Z0": aggregate["Z0"],
        "E00_fixed_Z1": aggregate["Z1"],
        "fixed_target": {
            "x": c.mpf(2),
            "t": c.mpf(1),
            "rho": c.mpf(4),
            "fixed": True,
            "independent": True,
            "not_last_source_cell_proxy": True,
        },
        "saved_outer_cells_reused": True,
        "W0_reconstructed": False,
        "W0_reweighted": False,
        "ancestor_rerun": False,
        "ancestor_helpers_recomputed": False,
        "prior_fixed_F1_Uz1_reused": True,
        "prior_fixed_F1_Uz1_recomputed": False,
        "F1_Uz1_included": False,
        "only_missing_E00_fields": True,
        "source_ids_retained": True,
        "source_paths_retained": True,
        "formal_powers_retained": True,
        "signed_coefficients_retained": True,
        "ordinary_Z1_product_rule": PRODUCT_RULE_Z1,
        "physical_source_measure": "dx=2*dt",
        "physical_measure_applied_once": True,
        "axis_endpoint_continuation": True,
        "axis_log_division_avoided": True,
        "green_kernel_masses_exact": True,
        "P1_is_state_component": True,
        "P1_is_m15": False,
        "P1_is_pressure_contribution": False,
        "pressure_g4_added": False,
        "N1p_pressure_added": False,
        "m15_includes_P1": False,
        "QP_added": False,
        "Q_P_added": False,
        "Q_branch_added": False,
        "P_branch_added": False,
        "allword_tail_added": False,
        "feedback_tail_added": False,
        "m15_computed": False,
        "ledger_appended": False,
        "complete_field_claimed": False,
    }


build_fixed_exit = build
build_x2 = build

