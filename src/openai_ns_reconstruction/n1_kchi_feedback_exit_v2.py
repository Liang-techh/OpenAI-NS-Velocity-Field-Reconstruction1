"""Independent v2 raw-Q fixed-core feedback exit at the physical point ``x=2``.

This is a new endpoint leaf.  It consumes saved data and does not call an
older core-exit builder or reconstruct an ancestor field.  For each of the
eight canonical source cells it forms the raw matrix product

``Q = C0 W0 + B1 W0_Z``

through ordinary ``Z1`` (the ``B1 W0_ZZ`` term is retained), and exposes the
three new branches

``E01 = G0 Q DeltaW0``
``E10 = D Q W0``
``E11 = D Q DeltaW0``

at the fixed target ``x=2``.  Here ``D=(-2 chi) T K0`` is only the first
physical correction.  The potential/conjugation compensation ``P=+2 chi
E51`` is deliberately not combined with these raw-Q sectors.

The source packets are caller-owned.  Saved ``W0`` packets have ordinary
``Z0/Z1/Z2`` and six columns; ``DeltaW0`` may be sparse and uses columns 1
and 5 for ``F1`` and ``d_x_F1``.  Matrix cells contain the required ``C0``
and ``B1`` supports through ``Z1``.  Every product keeps source IDs,
uncertainty classes, signed coefficients, and the four formal powers
``(hb, F0_anchor, Pstar_squared, Lambda)``.

The endpoint is independent of any target-cell profile.  Source-cell masses
are exact directed endpoint integrals over ``[i/8,(i+1)/8]`` in normalized
radius ``t=x/2``; the last source cell is never used as an endpoint proxy.
Caller-supplied kernel-tail packets are evaluated in separate symmetric
branches and are never merged into a signed branch.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any


SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
CELL_COUNT = 8
ORDERS = (0, 1, 2)
OUTPUT_ORDERS = (0, 1)
MATRIX_ORDERS = (0, 1)
FIELDS = ("F1", "Uz1", "K1", "P1", "d_x_F1", "d_x_Uz1")
DELTA_FIELDS = ("F1", "d_x_F1")
ROWS = ("4", "5", "6")

# Q rows are the rows of C0/B1.  The arrows below are the complete G0
# row-to-field contract at this endpoint.  D intentionally has only row 5.
FIELD_ROWS = {
    "F1": "5",
    "Uz1": "6",
    "K1": "6",
    "P1": "4",
    "d_x_F1": "5",
    "d_x_Uz1": "6",
}
D_ROW = "5"
ZERO_D_FIELDS = ("Uz1", "K1", "P1", "d_x_Uz1")

EXPECTED_C0 = frozenset(
    ("4,1", "5,1", "5,2", "5,3", "5,5",
     "6,1", "6,2", "6,3", "6,4", "6,6")
)
EXPECTED_B1 = frozenset(("5,1", "5,2", "5,3", "6,2", "6,3", "6,4"))

PRODUCT_RULE_Z0 = (
    ("C0", 0, 0, 1, "C0 W0"),
    ("B1", 0, 1, 1, "B1 W0_Z"),
)
PRODUCT_RULE_Z1 = (
    ("C0", 1, 0, 1, "C0_Z W0"),
    ("C0", 0, 1, 1, "C0 W0_Z"),
    ("B1", 1, 1, 1, "B1_Z W0_Z"),
    ("B1", 0, 2, 1, "B1 W0_ZZ"),
)
PRODUCT_RULE = {"Z0": PRODUCT_RULE_Z0, "Z1": PRODUCT_RULE_Z1}
CHI_PRODUCT_RULE_Z0 = ((0, 0, 1, "chi Q_Z0"),)
CHI_PRODUCT_RULE_Z1 = ((1, 0, 1, "chi_Z Q_Z0"),
                       (0, 1, 1, "chi Q_Z1"))
CHI_PRODUCT_RULE = {"Z0": CHI_PRODUCT_RULE_Z0,
                     "Z1": CHI_PRODUCT_RULE_Z1}


__all__ = [
    "SOURCE_POWER_BASIS", "CELL_COUNT", "ORDERS", "OUTPUT_ORDERS",
    "MATRIX_ORDERS", "FIELDS", "DELTA_FIELDS", "ROWS", "FIELD_ROWS",
    "D_ROW", "ZERO_D_FIELDS", "EXPECTED_C0", "EXPECTED_B1",
    "PRODUCT_RULE_Z0", "PRODUCT_RULE_Z1", "PRODUCT_RULE", "CHI_PRODUCT_RULE_Z0",
    "CHI_PRODUCT_RULE", "cell_kernel_masses", "delta_kernel_masses",
    "exit_kernel_masses", "green_cell_kernel_masses", "kernel_masses",
    "target_radius_kernel_masses", "correction_kernel_masses", "ordinary_product",
    "build_feedback_exit", "build", "build_core_feedback_exit", "build_endpoint",
    "build_fixed_exit", "build_x2", "build_delta_w0", "build_signed_delta_g5_exit",
    "signed_delta_g5_exit", "feedback_exit", "exit_x2",
]


def _lookup(mapping: Mapping[Any, Any], key: Any) -> Any:
    if key in mapping:
        return mapping[key]
    alternate = str(key) if not isinstance(key, str) else (
        int(key) if key.isdigit() else key
    )
    if alternate in mapping:
        return mapping[alternate]
    raise KeyError(key)


def _is_interval(value: Any) -> bool:
    return hasattr(value, "_mpi_") and hasattr(value, "a") and hasattr(value, "b")


def _interval(c: Any, value: Any) -> Any:
    if _is_interval(value):
        return value
    if hasattr(value, "a") and hasattr(value, "b"):
        lower = value.a._mpi_[0] if hasattr(value.a, "_mpi_") else value.a
        upper = value.b._mpi_[1] if hasattr(value.b, "_mpi_") else value.b
        return c.make_mpf((lower, upper))
    if isinstance(value, Mapping):
        if "lower" in value and "upper" in value:
            return c.mpf([value["lower"], value["upper"]])
        if "lo" in value and "hi" in value:
            return c.mpf([value["lo"], value["hi"]])
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        if len(value) == 2:
            return c.mpf([value[0], value[1]])
    return c.mpf(value)


def _bounds(c: Any, value: Any) -> tuple[Any, Any]:
    value = _interval(c, value)
    lo, hi = value._mpi_
    return c.make_mpf((lo, lo)), c.make_mpf((hi, hi))


def _endpoint(c: Any, value: Any, upper: bool) -> Any:
    value = _interval(c, value)
    raw = value._mpi_[1 if upper else 0]
    return c.make_mpf((raw, raw))


def _same_interval(left: Any, right: Any) -> bool:
    try:
        return left._mpi_ == right._mpi_
    except AttributeError:
        return left == right


def _exact_t(c: Any, index: int) -> Any:
    left = c.mpf(index) / CELL_COUNT
    right = c.mpf(index + 1) / CELL_COUNT
    return c.make_mpf((left._mpi_[0], right._mpi_[1]))


def _cell_t(c: Any, cell: Mapping[str, Any], index: int, label: str) -> Any:
    expected = _exact_t(c, index)
    candidate = None
    for key in ("t", "normalized_t", "radius_t"):
        if key in cell:
            candidate = cell[key]
            break
    if candidate is None:
        for key in ("x", "target_x", "radius"):
            if key in cell:
                candidate = _interval(c, cell[key]) / 2
                break
    if candidate is None and "rho" in cell:
        candidate = c.sqrt(_interval(c, cell["rho"])) / 2
    if candidate is None:
        return expected
    candidate = _interval(c, candidate)
    if not _same_interval(candidate, expected):
        raise ValueError(f"{label} {index} is not canonical t=[i/8,(i+1)/8]")
    return candidate


def _as_items(value: Any, *, keep_none: bool = False) -> list[Any]:
    if value is None:
        return [None] if keep_none else []
    if isinstance(value, (str, bytes, bytearray)):
        return [value]
    if isinstance(value, Sequence):
        return list(value)
    return [value]


def _unique(values: Sequence[Any]) -> list[Any]:
    result: list[Any] = []
    seen: set[Any] = set()
    for value in values:
        try:
            key = (type(value).__name__, value)
            hash(key)
        except TypeError:
            key = (type(value).__name__, repr(value))
        if key not in seen:
            seen.add(key)
            result.append(value)
    return result


def _origin(raw: Any) -> Mapping[str, Any]:
    if not isinstance(raw, Mapping):
        return {}
    value = raw.get("origin")
    return value if isinstance(value, Mapping) else {}


def _metadata_children(raw: Any) -> list[Any]:
    if not isinstance(raw, Mapping):
        return []
    result: list[Any] = []
    for key in ("origin", "metadata", "parents"):
        value = raw.get(key)
        if value is None or value is raw:
            continue
        if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
            result.extend(value)
        else:
            result.append(value)
    return result


def _declared(raw: Mapping[str, Any], names: Sequence[str]) -> tuple[list[Any], bool]:
    values: list[Any] = []
    present = False
    for name in names:
        if name in raw:
            present = True
            values.extend(_as_items(raw[name], keep_none=True))
    return values, present


def _source_ids(raw: Any, *, _seen: set[int] | None = None) -> list[Any]:
    if _seen is None:
        _seen = set()
    if not isinstance(raw, Mapping):
        if isinstance(raw, Sequence) and not isinstance(raw, (str, bytes, bytearray)):
            values: list[Any] = []
            for item in raw:
                values.extend(_source_ids(item, _seen=_seen))
            return _unique(values)
        return []
    marker = id(raw)
    if marker in _seen:
        return []
    _seen.add(marker)
    values, present = _declared(
        raw, ("source_ids", "sources", "source_id", "source",
              "source_cell_id", "input_id", "source_cell"),
    )
    if not present:
        # An explicit ``source_id=None`` is present above and must not be
        # replaced by a synthetic ID.  Nested metadata is only a fallback.
        for child in _metadata_children(raw):
            values.extend(_source_ids(child, _seen=_seen))
    else:
        for child in _metadata_children(raw):
            values.extend(_source_ids(child, _seen=_seen))
    return _unique(values)


def _error_classes(raw: Any, *, _seen: set[int] | None = None) -> list[Any]:
    if _seen is None:
        _seen = set()
    if not isinstance(raw, Mapping):
        return []
    marker = id(raw)
    if marker in _seen:
        return []
    _seen.add(marker)
    values: list[Any] = []
    values.extend(_as_items(raw.get("error_classes"), keep_none=False))
    if "error_classes" not in raw:
        values.extend(_as_items(raw.get("error_class"), keep_none=False))
        values.extend(_as_items(raw.get("uncertainty_class"), keep_none=False))
    for child in _metadata_children(raw):
        values.extend(_error_classes(child, _seen=_seen))
    return _unique(values)


def _uncertain(raw: Any, *, _seen: set[int] | None = None) -> bool:
    if _seen is None:
        _seen = set()
    if not isinstance(raw, Mapping):
        return False
    marker = id(raw)
    if marker in _seen:
        return False
    _seen.add(marker)
    operation = str(raw.get("operation", "")).lower()
    kind = str(raw.get("error_kind", raw.get("kind", ""))).lower()
    result = bool(raw.get("uncertainty", raw.get("is_uncertainty", False))
                  or raw.get("symmetric_error", False)
                  or any(token in operation or token in kind
                         for token in ("tail", "symmetric", "uncertainty", "cauchy")))
    return result or any(_uncertain(child, _seen=_seen)
                        for child in _metadata_children(raw))


def _powers(raw: Any, label: str) -> tuple[int, int, int, int]:
    if not isinstance(raw, Mapping):
        return (0, 0, 0, 0)
    value = raw.get("powers")
    if value is None:
        value = (
            raw.get("hb_power", 0),
            raw.get("F0_power", raw.get("F0_anchor_power", 0)),
            raw.get("Pstar_squared_power", raw.get("Pstar_power", 0)),
            raw.get("Lambda_power", raw.get("lambda_power", 0)),
        )
    if (not isinstance(value, Sequence)
            or isinstance(value, (str, bytes, bytearray))
            or len(value) != 4
            or any(isinstance(item, bool) or not isinstance(item, int)
                   for item in value)):
        raise ValueError(f"{label} must carry four integer formal powers")
    return tuple(int(item) for item in value)


def _term(c: Any, raw: Any, *, label: str, default_id: Any,
          error_only: bool = False, provenance: str = "saved_W0") -> dict[str, Any]:
    if isinstance(raw, Mapping):
        if "coefficient" in raw:
            coefficient = raw["coefficient"]
        elif "bound" in raw:
            coefficient = raw["bound"]
        elif "radius" in raw:
            coefficient = raw["radius"]
        else:
            raise ValueError(f"{label} is missing coefficient")
        result = dict(raw)
    else:
        coefficient = raw
        result = {}
    coefficient = _interval(c, coefficient)
    powers = _powers(raw, label)
    ids = _source_ids(raw)
    if not ids and isinstance(raw, Mapping):
        # No declared key means a synthetic fallback is appropriate.  An
        # explicit source_id=None was retained by _source_ids as [None].
        ids = [default_id]
    elif not ids:
        ids = [default_id]
    classes = _error_classes(raw)
    uncertain = bool(error_only or _uncertain(raw))
    origin = dict(_origin(raw))
    if error_only:
        classes = _unique(list(classes) + ["symmetric_W0_kernel_tail"])
    origin.update({
        "source_ids": list(ids),
        "error_classes": list(classes),
        "uncertainty": uncertain,
    })
    result.update({
        "powers": powers,
        "hb_power": powers[0],
        "F0_power": powers[1],
        "Pstar_squared_power": powers[2],
        "Lambda_power": powers[3],
        "coefficient": coefficient,
        "source_ids": list(ids),
        "error_classes": list(classes),
        "uncertainty": uncertain,
        "origin": origin,
        "W0_provenance": provenance,
    })
    if isinstance(raw, Mapping) and "source_id" in raw:
        result["source_id"] = raw["source_id"]
    elif ids:
        result["source_id"] = ids[0]
    if error_only:
        result.update({"error_only": True, "symmetric_error": True})
    return result


def _term_list(c: Any, raw: Any, *, label: str, default_id: Any,
               error_only: bool = False, provenance: str = "saved_W0") -> list[dict[str, Any]]:
    if raw is None:
        return []
    if isinstance(raw, Mapping) and any(
            key in raw for key in ("coefficient", "powers", "bound", "radius")):
        values = [raw]
    elif isinstance(raw, Sequence) and not isinstance(raw, (str, bytes, bytearray)):
        values = list(raw)
    else:
        values = [raw]
    return [
        _term(c, value, label=f"{label}[{index}]", default_id=default_id,
              error_only=error_only, provenance=provenance)
        for index, value in enumerate(values)
    ]


def _order_key(key: Any, label: str) -> int:
    if isinstance(key, bool):
        raise ValueError(f"boolean is not an order in {label}")
    if isinstance(key, int):
        return key
    if isinstance(key, str):
        text = key.strip().upper()
        if text.startswith("Z"):
            text = text[1:]
        try:
            return int(text)
        except ValueError as exc:
            raise ValueError(f"invalid {label} order {key!r}") from exc
    raise ValueError(f"invalid {label} order {key!r}")


def _order_map(raw: Any, required: Sequence[int], label: str) -> dict[int, Any]:
    if isinstance(raw, Mapping):
        result: dict[int, Any] = {}
        for key, value in raw.items():
            order = _order_key(key, label)
            if order in result:
                raise ValueError(f"duplicate {label} order {order}")
            result[order] = value
    elif isinstance(raw, Sequence) and not isinstance(raw, (str, bytes, bytearray)):
        result = {index: value for index, value in enumerate(raw)}
    else:
        raise TypeError(f"{label} must be an order mapping")
    if set(result) != set(required):
        raise ValueError(f"{label} must provide orders {sorted(required)}")
    return result


def _packet_mapping(cell: Mapping[str, Any]) -> Any:
    for key in (
        "W0_ordinary_Z", "ordinary_Z", "W0", "source_rows", "source_jets",
        "delta_W0", "delta_w0", "DeltaW0", "R", "error", "kernel_tail",
        "tail", "W0_kernel_tail", "W0_separate_symmetric_kernel_tail",
        "symmetric_error",
    ):
        if key in cell:
            return cell[key]
    return cell


def _direct_sparse_packet(cell: Mapping[str, Any], label: str) -> Any:
    if not ("F1" in cell or "d_x_F1" in cell):
        return None
    converted: dict[int, dict[int, Any]] = {}
    for order in ORDERS:
        f1 = cell.get("F1", {})
        dxf1 = cell.get("d_x_F1", {})
        f1_rows = _order_map(f1, ORDERS, f"{label} F1") if isinstance(f1, Mapping) else {}
        dxf1_rows = (_order_map(dxf1, ORDERS, f"{label} d_x_F1")
                     if isinstance(dxf1, Mapping) else {})
        converted[order] = {
            1: f1_rows.get(order, []),
            5: dxf1_rows.get(order, []),
            2: [], 3: [], 4: [], 6: [],
        }
    return converted


def _normalise_packet(c: Any, cell: Mapping[str, Any], index: int, *, label: str,
                      sparse: bool, error_only: bool,
                      provenance: str) -> dict[int, dict[int, list[dict[str, Any]]]]:
    raw = _packet_mapping(cell)
    if sparse and isinstance(cell, Mapping):
        direct = _direct_sparse_packet(cell, label)
        if direct is not None and raw is cell:
            raw = direct
    orders = _order_map(raw, ORDERS, f"{label} W0")
    packet: dict[int, dict[int, list[dict[str, Any]]]] = {}
    for order in ORDERS:
        row = orders[order]
        if not isinstance(row, Mapping):
            raise TypeError(f"{label} W0 Z{order} must be a column mapping")
        packet[order] = {}
        for column in range(1, 7):
            if str(column) in row:
                value = row[str(column)]
            elif column in row:
                value = row[column]
            elif sparse:
                value = []
            else:
                raise ValueError(f"{label} W0 Z{order} is missing column {column}")
            packet[order][column] = _term_list(
                c, value, label=f"{label} W0 Z{order} column {column}",
                default_id=f"W0:{index}:Z{order}:column{column}",
                error_only=error_only, provenance=provenance,
            )
            for term_index, term in enumerate(packet[order][column]):
                term.setdefault("W0_source_cell", index)
                term.setdefault("W0_source_cell_index", index)
                term.setdefault("W0_source_id", term.get("source_id", index))
                term.setdefault("source_cell", index)
                term.setdefault("source_cell_index", index)
                term.setdefault("source_term_index", term_index)
                term.setdefault("W0_order", order)
                term.setdefault("W0_column", column)
                term.setdefault("W0_column_index", column)
    return packet


def _normalise_cells(c: Any, cells: Any, *, label: str, sparse: bool,
                     error_only: bool = False,
                     provenance: str = "saved_W0") -> tuple[list[Any], list[dict[int, dict[int, list[dict[str, Any]]]]]]:
    if isinstance(cells, Mapping):
        raw_cells = [_lookup(cells, index) for index in range(CELL_COUNT)]
    elif isinstance(cells, Sequence) and not isinstance(cells, (str, bytes, bytearray)):
        raw_cells = list(cells)
    else:
        raise TypeError(f"{label}s must contain eight cells")
    if len(raw_cells) != CELL_COUNT:
        raise ValueError(f"eight {label}s are required")
    radii: list[Any] = []
    packets: list[dict[int, dict[int, list[dict[str, Any]]]]] = []
    for index, cell in enumerate(raw_cells):
        if not isinstance(cell, Mapping):
            raise TypeError(f"{label} {index} must be a mapping")
        radii.append(_cell_t(c, cell, index, label))
        packets.append(_normalise_packet(
            c, cell, index, label=label, sparse=sparse, error_only=error_only,
            provenance=provenance,
        ))
    return radii, packets


def _matrix_terms(c: Any, cell: Mapping[str, Any], name: str, order: int,
                  key: str, index: int) -> list[dict[str, Any]]:
    try:
        sector = cell[name]
        rows = _order_map(sector, MATRIX_ORDERS, f"matrix {name}")
        row = rows[order]
        if not isinstance(row, Mapping):
            raise TypeError
        if key in row:
            raw = row[key]
        elif any(str(candidate) == key for candidate in row):
            raw = next(row[candidate] for candidate in row if str(candidate) == key)
        else:
            raise KeyError(key)
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"matrix {index} is missing {name} Z{order} {key}") from exc
    terms = _term_list(
        c, raw, label=f"matrix:{index}:{name}:Z{order}:{key}",
        default_id=f"matrix:{index}:{name}:Z{order}:{key}",
        provenance="saved_C0_B1",
    )
    for term_index, term in enumerate(terms):
        term.setdefault("matrix_source_cell", index)
        term.setdefault("matrix_source_cell_index", index)
        term.setdefault("matrix_source_id", term.get("source_id", index))
        term.setdefault("matrix_name", name)
        term.setdefault("matrix_order", order)
        term.setdefault("matrix_key", key)
        term.setdefault("matrix_term_index", term_index)
    return terms


def _normalise_matrix(c: Any, cells: Any) -> list[dict[str, Any]]:
    if isinstance(cells, Mapping):
        raw_cells = [_lookup(cells, index) for index in range(CELL_COUNT)]
    elif isinstance(cells, Sequence) and not isinstance(cells, (str, bytes, bytearray)):
        raw_cells = list(cells)
    else:
        raise TypeError("matrix_cells must contain eight cells")
    if len(raw_cells) != CELL_COUNT:
        raise ValueError("eight saved matrix cells are required")
    output: list[dict[str, Any]] = []
    for index, raw in enumerate(raw_cells):
        if not isinstance(raw, Mapping):
            raise TypeError(f"matrix cell {index} must be a mapping")
        _cell_t(c, raw, index, "matrix cell")
        for name, expected in (("C0", EXPECTED_C0), ("B1", EXPECTED_B1)):
            if name not in raw or not isinstance(raw[name], Mapping):
                raise ValueError(f"matrix cell {index} is missing {name}")
            order_map = _order_map(raw[name], MATRIX_ORDERS, f"matrix {name}")
            for order in MATRIX_ORDERS:
                row = order_map[order]
                if not isinstance(row, Mapping) or {str(key) for key in row} != set(expected):
                    raise ValueError(f"matrix cell {index} {name} support changed at Z{order}")
                for key in expected:
                    _matrix_terms(c, raw, name, order, key, index)
        output.append(dict(raw))
    return output


def _combine_metadata(*terms: Mapping[str, Any], extra_ids: Sequence[Any] = (),
                      extra_classes: Sequence[Any] = (),
                      force_uncertainty: bool = False) -> tuple[list[Any], list[Any], bool]:
    ids: list[Any] = []
    classes: list[Any] = []
    uncertain = force_uncertainty
    for term in terms:
        ids.extend(_source_ids(term))
        classes.extend(_error_classes(term))
        uncertain = uncertain or _uncertain(term)
    ids.extend(extra_ids)
    classes.extend(extra_classes)
    return _unique(ids), _unique(classes), uncertain


def _q_product(matrix_term: Mapping[str, Any], w0_term: Mapping[str, Any], *,
               source_index: int, matrix_name: str, matrix_order: int,
               matrix_key: str, matrix_term_index: int, w0_order: int,
               w0_term_index: int, branch: str, origin: str) -> dict[str, Any]:
    powers = tuple(a + b for a, b in zip(matrix_term["powers"], w0_term["powers"]))
    ids, classes, uncertain = _combine_metadata(
        matrix_term, w0_term, extra_ids=(f"source_cell:{source_index}",)
    )
    row, column = (int(value) for value in matrix_key.split(","))
    result = {
        "powers": powers,
        "hb_power": powers[0], "F0_power": powers[1],
        "Pstar_squared_power": powers[2], "Lambda_power": powers[3],
        "coefficient": matrix_term["coefficient"] * w0_term["coefficient"],
        "source_ids": ids,
        "error_classes": classes,
        "uncertainty": uncertain,
        "origin": origin,
        "origin_record": {
            "operation": "raw_Q_product",
            "source_ids": ids,
            "error_classes": classes,
            "uncertainty": uncertain,
            "symbolic_cancellation_claimed": False,
        },
        "branch": branch,
        "source_cell": source_index,
        "source_cell_index": source_index,
        "source_cell_id": w0_term.get("source_cell_id", source_index),
        "row": row,
        "column": column,
        "matrix_row": row,
        "matrix_column": column,
        "matrix_name": matrix_name,
        "matrix_key": matrix_key,
        "matrix_order": matrix_order,
        "matrix_term_index": matrix_term_index,
        "matrix_source_cell": matrix_term.get("matrix_source_cell", source_index),
        "matrix_source_cell_index": matrix_term.get("matrix_source_cell_index", source_index),
        "matrix_source_id": matrix_term.get("matrix_source_id", matrix_term.get("source_id")),
        "W0_order": w0_order,
        "W0_term_index": w0_term_index,
        "W0_column": column,
        "W0_column_index": column,
        "W0_source_cell": w0_term.get("W0_source_cell", source_index),
        "W0_source_cell_index": w0_term.get("W0_source_cell_index", source_index),
        "W0_source_id": w0_term.get("W0_source_id", w0_term.get("source_id")),
        "W0_provenance": w0_term.get("W0_provenance", "saved_W0"),
        "leibniz_multiplier": 1,
        "signed_source_term": not bool(w0_term.get("error_only", False)),
    }
    if matrix_term.get("error_only") or w0_term.get("error_only"):
        result.update({"error_only": True, "symmetric_error": True,
                       "signed_source_term": False})
    return result


def _build_q(c: Any, matrix_cell: Mapping[str, Any],
             packet: Mapping[int, Mapping[int, Sequence[Mapping[str, Any]]]], *,
             source_index: int, branch: str) -> dict[str, dict[str, list[dict[str, Any]]]]:
    result: dict[str, dict[str, list[dict[str, Any]]]] = {
        "Z0": {row: [] for row in ROWS},
        "Z1": {row: [] for row in ROWS},
    }
    for order_name, rules in PRODUCT_RULE.items():
        for matrix_name, matrix_order, w0_order, multiplier, origin in rules:
            expected = EXPECTED_C0 if matrix_name == "C0" else EXPECTED_B1
            for matrix_key in sorted(expected):
                row = matrix_key.split(",")[0]
                column = int(matrix_key.split(",")[1])
                matrix_terms = _matrix_terms(
                    c, matrix_cell, matrix_name, matrix_order, matrix_key, source_index
                )
                for matrix_term_index, matrix_term in enumerate(matrix_terms):
                    for w0_term_index, w0_term in enumerate(packet[w0_order][column]):
                        product = _q_product(
                            matrix_term, w0_term,
                            source_index=source_index,
                            matrix_name=matrix_name,
                            matrix_order=matrix_order,
                            matrix_key=matrix_key,
                            matrix_term_index=matrix_term_index,
                            w0_order=w0_order,
                            w0_term_index=w0_term_index,
                            branch=branch,
                            origin=origin,
                        )
                        if multiplier != 1:
                            product["coefficient"] = multiplier * product["coefficient"]
                            product["leibniz_multiplier"] = multiplier
                        result[order_name][row].append(product)
    return result


def _nonnegative(c: Any, value: Any, label: str) -> Any:
    value = _interval(c, value)
    lo, hi = value._mpi_
    if c.make_mpf((hi, hi)) < 0:
        raise ValueError(f"{label} has a negative upper bound")
    if c.make_mpf((lo, lo)) < 0:
        lo = c.mpf(0)._mpi_[0]
    return c.make_mpf((lo, hi))


def _primitive_uz(c: Any, target: Any, source: Any) -> Any:
    source = _interval(c, source)
    if _bounds(c, source)[0] == 0:
        return c.mpf(0)
    return 2 * source * source * c.log(target / source) + source * source


def cell_kernel_masses(c: Any, target_t: Any, source_left_t: Any,
                       source_right_t: Any, *, triangular: bool = False) -> dict[str, Any]:
    """Return exact directed six-component ``G0`` masses in normalized radius."""

    target = _interval(c, target_t)
    left = _interval(c, source_left_t)
    right = target if triangular else _interval(c, source_right_t)
    tlo, thi = _bounds(c, target)
    alo, ahi = _bounds(c, left)
    blo, bhi = _bounds(c, right)
    if alo < 0 or ahi < alo or blo < alo or bhi < blo:
        raise ValueError("source endpoints must be ordered and nonnegative")
    if not triangular and bhi > thi:
        raise ValueError("source cell extends beyond target radius")
    if triangular and not _same_interval(right, target):
        raise ValueError("triangular endpoint must equal target interval")
    if triangular and alo == 0 and tlo == 0:
        f1 = target * target / 2
        uz1 = target * target
        p1 = 2 * target
        dxf1 = target / 2
        dxuz1 = target
    else:
        if thi <= 0:
            raise ValueError("positive target radius required")
        target2 = target * target
        left2 = left * left
        right2 = right * right
        f1 = (right2 - left2) - (right2 * right2 - left2 * left2) / (2 * target2)
        uz1 = _primitive_uz(c, target, right) - _primitive_uz(c, target, left)
        p1 = 2 * (right - left)
        dxf1 = (right2 * right2 - left2 * left2) / (2 * target * target2)
        dxuz1 = (right2 - left2) / target
    f1 = _nonnegative(c, f1, "G0 F1 mass")
    uz1 = _nonnegative(c, uz1, "G0 Uz1 mass")
    p1 = _nonnegative(c, p1, "G0 P1 mass")
    dxf1 = _nonnegative(c, dxf1, "G0 d_x_F1 mass")
    dxuz1 = _nonnegative(c, dxuz1, "G0 d_x_Uz1 mass")
    return {
        "F1": f1, "Uz1": uz1, "K1": -f1, "P1": p1,
        "d_x_F1": dxf1, "d_x_Uz1": dxuz1,
        "target_t": target, "source_left_t": left, "source_right_t": right,
        "target_x": 2 * target, "source_left_x": 2 * left,
        "source_right_x": 2 * right,
        "principal_kernel": "G0", "negative_kernel_K1": True,
        "triangular": bool(triangular),
        "axis_continuation": bool(triangular and alo == 0 and tlo == 0),
        "axis_endpoint_continuation": bool(alo == 0),
        "axis_log_division_avoided": bool(alo == 0),
        "exact_source_cell_mass": True,
    }


def _k1_primitive(c: Any, x: Any, source: Any) -> Any:
    source = _interval(c, source)
    if _bounds(c, source)[0] == 0:
        return c.mpf(0)
    if _bounds(c, x)[1] <= 0:
        raise ValueError("positive physical exit radius required")
    return (x * x * source * source / 32
            - source**4 * c.log(x / source) / 16
            - source**4 / 64
            - source**6 / (96 * x * x))


def _k1_dx_primitive(c: Any, x: Any, source: Any) -> Any:
    source = _interval(c, source)
    if _bounds(c, source)[0] == 0:
        return c.mpf(0)
    if _bounds(c, x)[1] <= 0:
        raise ValueError("positive physical exit radius required")
    return x * source * source / 16 - source**4 / (16 * x) + source**6 / (48 * x**3)


def delta_kernel_masses(c: Any, target_x: Any, source_left_x: Any,
                        source_right_x: Any) -> dict[str, Any]:
    """Return exact bare physical ``T K0`` masses before ``-2 chi``."""

    target = _interval(c, target_x)
    left = _interval(c, source_left_x)
    right = _interval(c, source_right_x)
    tlo, thi = _bounds(c, target)
    alo, ahi = _bounds(c, left)
    blo, bhi = _bounds(c, right)
    if alo < 0 or ahi < alo or blo < alo or bhi < blo:
        raise ValueError("source physical endpoints must be ordered and nonnegative")
    if bhi > thi:
        raise ValueError("source physical cell extends beyond exit radius")
    if thi <= 0:
        raise ValueError("positive physical exit radius required")
    f1 = _k1_primitive(c, target, right) - _k1_primitive(c, target, left)
    dxf1 = _k1_dx_primitive(c, target, right) - _k1_dx_primitive(c, target, left)
    return {
        "F1": _nonnegative(c, f1, "D F1 mass"),
        "d_x_F1": _nonnegative(c, dxf1, "D d_x_F1 mass"),
        "target_x": target, "source_left_x": left, "source_right_x": right,
        "target_t": target / 2, "source_left_t": left / 2,
        "source_right_t": right / 2,
        "principal_kernel": "T K0 first correction",
        "negative_kernel_K1": True, "kernel_order": 1,
        "axis_log_division_avoided": bool(alo == 0),
        "exact_source_cell_mass": True, "physical_units": True,
    }


exit_kernel_masses = cell_kernel_masses
green_cell_kernel_masses = cell_kernel_masses
kernel_masses = cell_kernel_masses
target_radius_kernel_masses = cell_kernel_masses


def correction_kernel_masses(c: Any, target_x: Any, source_left_x: Any,
                             source_right_x: Any, *, chi: Any) -> dict[str, Any]:
    """Return signed ``(-2 chi) T K0`` masses without forming a Q product."""

    bare = delta_kernel_masses(c, target_x, source_left_x, source_right_x)
    factor = -2 * _interval(c, chi)
    return {
        "F1": factor * bare["F1"],
        "d_x_F1": factor * bare["d_x_F1"],
        "bare_F1": bare["F1"], "bare_d_x_F1": bare["d_x_F1"],
        "correction_factor": factor,
        "principal_kernel": "(-2 chi) T K0 first correction",
        "exact_source_cell_mass": True,
    }


def ordinary_product(left: Sequence[Any], right: Sequence[Any]) -> list[Any]:
    """Multiply ordinary rows through ``Z2`` with the ordinary factor 2."""

    if len(left) != 3 or len(right) != 3:
        raise ValueError("ordinary rows must contain Z0, Z1 and Z2")
    return [left[0] * right[0],
            left[1] * right[0] + left[0] * right[1],
            left[2] * right[0] + 2 * left[1] * right[1] + left[0] * right[2]]


def _empty_fields() -> dict[str, list[dict[str, Any]]]:
    return {field: [] for field in FIELDS}


def _scale_g0(term: Mapping[str, Any], mass: Any, *, field: str,
              source_index: int, branch: str) -> dict[str, Any]:
    result = dict(term)
    result.update({
        "coefficient": term["coefficient"] * mass,
        "field": field,
        "source_cell": source_index,
        "source_cell_index": source_index,
        "branch": branch,
        "kernel_mass": mass,
        "kernel_mass_bare": mass,
        "kernel_mass_field": field,
        "principal_kernel": "G0",
        "negative_kernel_K1": field == "K1",
        "exact_source_cell_mass": True,
        "signed_source_term": not bool(term.get("error_only", False)),
    })
    return result


def _g0_result(q: Mapping[str, Mapping[str, Sequence[Mapping[str, Any]]]],
               masses: Mapping[str, Any], *, source_index: int,
               branch: str) -> dict[str, dict[str, list[dict[str, Any]]]]:
    output = {f"Z{order}": _empty_fields() for order in OUTPUT_ORDERS}
    for order in OUTPUT_ORDERS:
        for field in FIELDS:
            row = FIELD_ROWS[field]
            output[f"Z{order}"][field] = [
                _scale_g0(term, masses[field], field=field,
                          source_index=source_index, branch=branch)
                for term in q[f"Z{order}"][row]
            ]
    return output


def _chi_rows(c: Any, value: Any) -> dict[int, list[dict[str, Any]]]:
    if value is None:
        raise ValueError("chi Z0/Z1 sectors are required for D branches")
    if isinstance(value, Mapping) and any(key in value for key in ("coefficient", "powers")):
        value = {0: value}
    elif isinstance(value, Mapping):
        for wrapper in ("ordinary_Z", "chi_jets", "jets", "rows"):
            if wrapper in value:
                value = value[wrapper]
                break
    if isinstance(value, Mapping):
        parsed: dict[int, Any] = {}
        for key, raw in value.items():
            order = _order_key(key, "chi")
            if order not in OUTPUT_ORDERS:
                raise ValueError(f"chi supports only Z0 and Z1, got {key!r}")
            if order in parsed:
                raise ValueError(f"duplicate chi order {order}")
            parsed[order] = raw
    elif isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        if (value and all(isinstance(item, Mapping)
                          and any(key in item for key in ("coefficient", "powers"))
                          for item in value)):
            parsed = {0: value}
        elif len(value) == 2:
            parsed = {0: value[0], 1: value[1]}
        else:
            parsed = {0: value}
    else:
        parsed = {0: value}
    if 0 not in parsed:
        raise ValueError("chi must provide Z0")
    if 1 not in parsed:
        parsed[1] = {"coefficient": c.mpf(0), "powers": (0, 0, 0, 0),
                     "source_id": "chi:Z1:structural_zero"}
    result: dict[int, list[dict[str, Any]]] = {}
    for order in OUTPUT_ORDERS:
        result[order] = _term_list(
            c, parsed[order], label=f"chi Z{order}",
            default_id=f"chi:Z{order}", provenance="saved_chi",
        )
        for index, term in enumerate(result[order]):
            term.setdefault("chi_order", order)
            term.setdefault("chi_term_index", index)
            term.setdefault("chi_source_id", term.get("source_id", f"chi:Z{order}"))
    return result


def _d_result(q: Mapping[str, Mapping[str, Sequence[Mapping[str, Any]]]],
              masses: Mapping[str, Any], chi: Mapping[int, Sequence[Mapping[str, Any]]], *,
              source_index: int, branch: str) -> dict[str, dict[str, list[dict[str, Any]]]]:
    """Apply D to the row-5 Q sector, retaining explicit structural zeros."""

    output = {f"Z{order}": _empty_fields() for order in OUTPUT_ORDERS}
    rules_by_order = {0: CHI_PRODUCT_RULE_Z0, 1: CHI_PRODUCT_RULE_Z1}
    for outer_order in OUTPUT_ORDERS:
        for chi_order, q_order, multiplier, origin in rules_by_order[outer_order]:
            for chi_index, chi_term in enumerate(chi[chi_order]):
                for q_index, q_term in enumerate(q[f"Z{q_order}"][D_ROW]):
                    powers = tuple(a + b for a, b in zip(
                        chi_term["powers"], q_term["powers"]
                    ))
                    ids, classes, uncertain = _combine_metadata(
                        chi_term, q_term,
                        extra_ids=(f"source_cell:{source_index}",),
                    )
                    coefficient = -2 * multiplier * chi_term["coefficient"] * q_term["coefficient"]
                    error_only = bool(chi_term.get("error_only") or q_term.get("error_only"))
                    base: dict[str, Any] = {
                        "powers": powers,
                        "hb_power": powers[0], "F0_power": powers[1],
                        "Pstar_squared_power": powers[2], "Lambda_power": powers[3],
                        "coefficient": coefficient,
                        "source_ids": ids,
                        "error_classes": classes,
                        "uncertainty": uncertain,
                        "origin": origin,
                        "origin_record": {
                            "operation": "raw_D_Q_product",
                            "source_ids": ids,
                            "error_classes": classes,
                            "uncertainty": uncertain,
                            "chi_order": chi_order,
                            "Q_order": q_order,
                            "outer_order": outer_order,
                            "symbolic_cancellation_claimed": False,
                        },
                        "branch": branch,
                        "source_cell": source_index,
                        "source_cell_index": source_index,
                        "source_cell_id": q_term.get("source_cell_id", source_index),
                        "chi_order": chi_order,
                        "chi_term_index": chi_index,
                        "chi_source_id": chi_term.get("chi_source_id", "chi"),
                        "Q_order": q_order,
                        "Q_term_index": q_index,
                        "Q_origin": q_term.get("origin"),
                        "matrix_name": q_term.get("matrix_name"),
                        "matrix_key": q_term.get("matrix_key"),
                        "matrix_order": q_term.get("matrix_order"),
                        "matrix_source_id": q_term.get("matrix_source_id"),
                        "W0_order": q_term.get("W0_order"),
                        "W0_column": q_term.get("W0_column"),
                        "W0_source_id": q_term.get("W0_source_id"),
                        "W0_provenance": q_term.get("W0_provenance"),
                        "outer_leibniz_multiplier": multiplier,
                        "correction_factor": -2 * chi_term["coefficient"],
                        "principal_kernel": "D=(-2 chi) T K0",
                        "outer_source": "D Q",
                        "signed_source_term": not error_only,
                        "finite_branch": not error_only,
                    }
                    if error_only:
                        base.update({"error_only": True, "symmetric_error": True})
                    for field in DELTA_FIELDS:
                        item = dict(base)
                        item.update({
                            "field": field,
                            "coefficient": coefficient * masses[field],
                            "kernel_mass": masses[field],
                            "kernel_mass_bare": masses[field],
                            "kernel_mass_field": field,
                            "exact_source_cell_mass": True,
                        })
                        output[f"Z{outer_order}"][field].append(item)
    return output


def _branch_summary(source_results: Mapping[str, Any], *, q_by_source: Mapping[str, Any],
                    masses_by_source: Mapping[str, Any], branch: str,
                    formula: str, operator: str, complete_six: bool,
                    error_only: bool = False) -> dict[str, Any]:
    aggregate = {f"Z{order}": _empty_fields() for order in OUTPUT_ORDERS}
    for value in source_results.values():
        fields = value["fields"]
        for order in OUTPUT_ORDERS:
            for field in FIELDS:
                aggregate[f"Z{order}"][field].extend(fields[f"Z{order}"][field])
    result = {
        "branch": branch,
        "formula": formula,
        "operator": operator,
        "Q_by_source_cell": dict(q_by_source),
        "kernel_masses_by_source_cell": dict(masses_by_source),
        "exit_by_source_cell": dict(source_results),
        "exit_x2_by_source_cell": dict(source_results),
        "aggregate": aggregate,
        "fields": aggregate,
        "complete_six_component": complete_six,
        "nonzero_fields": list(FIELDS if complete_six else DELTA_FIELDS),
        "structural_zero_fields": [] if complete_six else list(ZERO_D_FIELDS),
        "source_cells_collapsed": False,
        "exact_source_cell_masses": True,
        "fixed_exit": True,
        "signed_coefficients_retained": not error_only,
        "source_ids_retained": True,
        "formal_powers_retained": True,
        "error_only": error_only,
        "symmetric_error": error_only,
        "finite_branch": not error_only,
        "potential_P_compensation_included": False,
    }
    # Keep the operation-specific view names explicit for consumers that do
    # not want to infer whether a branch is G0 or D from ``operator``.
    if operator == "G0":
        result.update({
            "G0_exit_by_source_cell": dict(source_results),
            "G0_Z0_core_exit": aggregate["Z0"],
            "G0_Z1_core_exit": aggregate["Z1"],
            "DeltaG1_exit_by_source_cell": {},
            "DeltaG1_Z0_core_exit": _empty_fields(),
            "DeltaG1_Z1_core_exit": _empty_fields(),
        })
    else:
        result.update({
            "G0_exit_by_source_cell": {},
            "G0_Z0_core_exit": _empty_fields(),
            "G0_Z1_core_exit": _empty_fields(),
            "DeltaG1_exit_by_source_cell": dict(source_results),
            "DeltaG1_Z0_core_exit": aggregate["Z0"],
            "DeltaG1_Z1_core_exit": aggregate["Z1"],
        })
    return result


def _q_only_summary(q_by_source: Mapping[str, Any], *, branch: str,
                    provenance: str) -> dict[str, Any]:
    return {
        "branch": branch,
        "formula": "Q = C0 W0 + B1 W0_Z",
        "Q_by_source_cell": dict(q_by_source),
        "source_cells_collapsed": False,
        "saved_W0_reused": True,
        "W0_reconstructed": False,
        "E00_G0_recomputed": False,
        "signed_coefficients_retained": provenance != "caller_W0_kernel_tail",
        "source_ids_retained": True,
        "formal_powers_retained": True,
        "provenance": provenance,
    }


def _resolve(cells: Any, *aliases: Any) -> Any:
    if cells is not None:
        return cells
    for value in aliases:
        if value is not None:
            return value
    return None


def _embedded_tail(cells: Any) -> Any:
    if isinstance(cells, Mapping):
        raw = [_lookup(cells, index) for index in range(CELL_COUNT)]
    elif isinstance(cells, Sequence) and not isinstance(cells, (str, bytes, bytearray)):
        raw = list(cells)
    else:
        return None
    if len(raw) != CELL_COUNT or not all(isinstance(cell, Mapping) for cell in raw):
        return None
    keys = ("W0_kernel_tail", "kernel_tail", "tail",
            "W0_separate_symmetric_kernel_tail", "symmetric_error")
    if not any(any(key in cell for key in keys) for cell in raw):
        return None
    return [{**{key: cell[key] for key in ("t", "x", "target_x", "radius", "rho")
               if key in cell},
             "W0_ordinary_Z": next(cell[key] for key in keys if key in cell)}
            for cell in raw]


def _source_result(index: int, radius: Any, target_t: Any, masses: Mapping[str, Any],
                   fields: Mapping[str, Any], *, operator: str, branch: str,
                   q: Mapping[str, Any], error_only: bool) -> dict[str, Any]:
    result = {
        "source_index": index,
        "source_cell": index,
        "source_t": radius,
        "source_radius": radius,
        "target_t": target_t,
        "target_x": 2 * target_t,
        "target_radius": 4 * target_t * target_t,
        "fixed_exit": True,
        "kernel_masses": masses,
        "Q": q,
        "fields": fields,
        "operator": operator,
        "branch": branch,
        "error_only": error_only,
        "structural_zero_fields": list(ZERO_D_FIELDS) if operator == "D" else [],
        "nonzero_fields": list(DELTA_FIELDS) if operator == "D" else list(FIELDS),
    }
    result["G0"] = fields if operator == "G0" else {}
    result["DeltaG1"] = fields if operator == "D" else {}
    return result


def build_feedback_exit(
    c: Any,
    saved_w0_cells: Any = None,
    matrix_cells: Any = None,
    chi: Any = None,
    *,
    outer_cells: Any = None,
    saved_outer_cells: Any = None,
    saved_W0_cells: Any = None,
    w0_cells: Any = None,
    W0_cells: Any = None,
    C0_B1_cells: Any = None,
    delta_w0_cells: Any = None,
    delta_W0_cells: Any = None,
    delta_w0: Any = None,
    delta_W0: Any = None,
    chi_jets: Any = None,
    W0_kernel_tail_cells: Any = None,
    w0_kernel_tail_cells: Any = None,
    w0_tail_cells: Any = None,
    saved_w0_kernel_tail_cells: Any = None,
    symmetric_W0_kernel_tail_cells: Any = None,
    kernel_tail_cells: Any = None,
    symmetric_error_cells: Any = None,
    symmetric_R_cells: Any = None,
    R_cells: Any = None,
    remainder_cells: Any = None,
    branch: str | None = None,
) -> dict[str, Any]:
    """Build independent raw-Q fixed ``x=2`` feedback branches.

    Positional arguments are ``(c, saved_W0_cells, matrix_cells, chi)``.
    Keyword aliases retain the saved-packet vocabulary used by earlier
    leaves.  The returned signed branches are never merged with the optional
    caller-supplied W0 kernel-tail branches.
    """

    if not hasattr(c, "mpf") or not hasattr(c, "make_mpf"):
        raise TypeError("c must be a directed MPIntervalContext")
    saved_w0_cells = _resolve(
        saved_w0_cells, outer_cells, saved_outer_cells, saved_W0_cells,
        w0_cells, W0_cells
    )
    matrix_cells = _resolve(matrix_cells, C0_B1_cells)
    if saved_w0_cells is None:
        raise ValueError("caller-owned saved W0 cells are required")
    if matrix_cells is None:
        raise ValueError("caller-owned C0/B1 matrix cells are required")
    delta_input = _resolve(delta_w0_cells, delta_W0_cells, delta_w0, delta_W0)
    tail_input = _resolve(
        W0_kernel_tail_cells, w0_kernel_tail_cells, w0_tail_cells,
        saved_w0_kernel_tail_cells, symmetric_W0_kernel_tail_cells,
        kernel_tail_cells, symmetric_error_cells, symmetric_R_cells,
        R_cells, remainder_cells,
    )
    if tail_input is None:
        tail_input = _embedded_tail(saved_w0_cells)
    chi_value = _resolve(chi, chi_jets)
    matrices = _normalise_matrix(c, matrix_cells)
    radii, saved_packets = _normalise_cells(
        c, saved_w0_cells, label="saved W0 cell", sparse=False,
        provenance="saved_W0",
    )
    saved_q: dict[str, Any] = {}
    g0_masses: dict[str, Any] = {}
    d_masses: dict[str, Any] = {}
    target_t = c.mpf(1)
    target_x = c.mpf(2)
    for index in range(CELL_COUNT):
        saved_q[str(index)] = _build_q(
            c, matrices[index], saved_packets[index], source_index=index,
            branch="saved_W0_Q",
        )
        left = _endpoint(c, radii[index], False)
        right = _endpoint(c, radii[index], True)
        g0_masses[str(index)] = cell_kernel_masses(
            c, target_t, left, right, triangular=False,
        )
        d_masses[str(index)] = delta_kernel_masses(
            c, target_x, 2 * left, 2 * right,
        )

    delta_radii: list[Any] | None = None
    delta_q: dict[str, Any] = {}
    if delta_input is not None:
        delta_radii, delta_packets = _normalise_cells(
            c, delta_input, label="Delta W0 cell", sparse=True,
            provenance="saved_DeltaW0",
        )
        for index in range(CELL_COUNT):
            if not _same_interval(delta_radii[index], radii[index]):
                raise ValueError(f"Delta W0 cell {index} coordinate differs from saved W0")
            delta_q[str(index)] = _build_q(
                c, matrices[index], delta_packets[index], source_index=index,
                branch="DeltaW0_Q",
            )

    chi_rows: dict[int, list[dict[str, Any]]] | None = None
    if chi_value is not None:
        chi_rows = _chi_rows(c, chi_value)

    mode = None if branch is None else str(branch).strip().lower().replace("-", "_")
    aliases = {
        "e01": "E01", "e01_g0_q_delta": "E01", "g0_q_delta": "E01",
        "e10": "E10", "e10_d_q_w0": "E10", "d_q_w0": "E10",
        "e11": "E11", "e11_d_q_delta": "E11", "d_q_delta": "E11",
    }
    selected: set[str]
    if mode in (None, "all", "both", "raw_q", "rawq"):
        selected = {"E01", "E10", "E11"}
    elif mode in aliases:
        selected = {aliases[mode]}
    elif mode in ("q", "baseline_q", "saved_w0_q"):
        selected = set()
    else:
        raise ValueError(f"unknown raw-Q feedback branch {branch!r}")
    if "E01" in selected and delta_input is None:
        raise ValueError("E01 requires caller-owned Delta W0 cells")
    if ("E10" in selected or "E11" in selected) and chi_rows is None:
        raise ValueError("E10/E11 require chi Z0/Z1 sectors")
    if "E11" in selected and delta_input is None:
        raise ValueError("E11 requires caller-owned Delta W0 cells")

    branches: dict[str, Any] = {}
    if "E01" in selected:
        source_results: dict[str, Any] = {}
        for index in range(CELL_COUNT):
            fields = _g0_result(
                delta_q[str(index)], g0_masses[str(index)],
                source_index=index, branch="E01_G0_Q_Delta",
            )
            source_results[str(index)] = _source_result(
                index, (delta_radii or radii)[index], target_t, g0_masses[str(index)],
                fields, operator="G0", branch="E01_G0_Q_Delta",
                q=delta_q[str(index)], error_only=False,
            )
        branches["E01"] = _branch_summary(
            source_results, q_by_source=delta_q, masses_by_source=g0_masses,
            branch="E01_G0_Q_Delta", formula="G0 Q DeltaW0", operator="G0",
            complete_six=True,
        )
    if "E10" in selected:
        source_results = {}
        for index in range(CELL_COUNT):
            fields = _d_result(
                saved_q[str(index)], d_masses[str(index)], chi_rows,
                source_index=index, branch="E10_D_Q_W0",
            )
            source_results[str(index)] = _source_result(
                index, radii[index], target_t, d_masses[str(index)], fields,
                operator="D", branch="E10_D_Q_W0", q=saved_q[str(index)],
                error_only=False,
            )
        branches["E10"] = _branch_summary(
            source_results, q_by_source=saved_q, masses_by_source=d_masses,
            branch="E10_D_Q_W0", formula="D Q W0", operator="D",
            complete_six=False,
        )
    if "E11" in selected:
        source_results = {}
        for index in range(CELL_COUNT):
            fields = _d_result(
                delta_q[str(index)], d_masses[str(index)], chi_rows,
                source_index=index, branch="E11_D_Q_Delta",
            )
            source_results[str(index)] = _source_result(
                index, (delta_radii or radii)[index], target_t,
                d_masses[str(index)], fields, operator="D",
                branch="E11_D_Q_Delta", q=delta_q[str(index)], error_only=False,
            )
        branches["E11"] = _branch_summary(
            source_results, q_by_source=delta_q, masses_by_source=d_masses,
            branch="E11_D_Q_Delta", formula="D Q DeltaW0", operator="D",
            complete_six=False,
        )

    tail_branches: dict[str, Any] = {}
    tail_q: dict[str, Any] = {}
    tail_radii: list[Any] | None = None
    if tail_input is not None:
        tail_radii, tail_packets = _normalise_cells(
            c, tail_input, label="W0 kernel-tail cell", sparse=True,
            error_only=True, provenance="caller_W0_kernel_tail",
        )
        for index in range(CELL_COUNT):
            if not _same_interval(tail_radii[index], radii[index]):
                raise ValueError(f"W0 kernel-tail cell {index} coordinate differs from saved W0")
            tail_q[str(index)] = _build_q(
                c, matrices[index], tail_packets[index], source_index=index,
                branch="W0_kernel_tail_Q",
            )
        tail_g0_sources: dict[str, Any] = {}
        for index in range(CELL_COUNT):
            fields = _g0_result(
                tail_q[str(index)], g0_masses[str(index)], source_index=index,
                branch="W0_kernel_tail_G0_Q",
            )
            tail_g0_sources[str(index)] = _source_result(
                index, tail_radii[index], target_t, g0_masses[str(index)], fields,
                operator="G0", branch="W0_kernel_tail_G0_Q", q=tail_q[str(index)],
                error_only=True,
            )
        tail_branches["G0"] = _branch_summary(
            tail_g0_sources, q_by_source=tail_q, masses_by_source=g0_masses,
            branch="W0_kernel_tail_G0_Q", formula="G0 Q R_W0", operator="G0",
            complete_six=True, error_only=True,
        )
        if chi_rows is not None:
            tail_d_sources: dict[str, Any] = {}
            for index in range(CELL_COUNT):
                fields = _d_result(
                    tail_q[str(index)], d_masses[str(index)], chi_rows,
                    source_index=index, branch="W0_kernel_tail_D_Q",
                )
                tail_d_sources[str(index)] = _source_result(
                    index, tail_radii[index], target_t, d_masses[str(index)], fields,
                    operator="D", branch="W0_kernel_tail_D_Q", q=tail_q[str(index)],
                    error_only=True,
                )
            tail_branches["D"] = _branch_summary(
                tail_d_sources, q_by_source=tail_q, masses_by_source=d_masses,
                branch="W0_kernel_tail_D_Q", formula="D Q R_W0", operator="D",
                complete_six=False, error_only=True,
            )
        branches["W0_kernel_tail_G0"] = tail_branches["G0"]
        if "D" in tail_branches:
            branches["W0_kernel_tail_D"] = tail_branches["D"]

    full_g0 = cell_kernel_masses(c, target_t, c.mpf(0), target_t, triangular=True)
    full_d = delta_kernel_masses(c, target_x, c.mpf(0), target_x)
    baseline_q = _q_only_summary(saved_q, branch="baseline_Q_only",
                                 provenance="saved_W0")
    baseline = dict(baseline_q)
    baseline.update({
        "G0_exit_by_source_cell": {},
        "DeltaG1_exit_by_source_cell": {},
        "G0_Z0_core_exit": None,
        "G0_Z1_core_exit": None,
        "E00_G0_recomputed": False,
        "complete_six_component": False,
        "raw_Q_only": True,
    })
    tail_summary = {
        "branches": tail_branches,
        "error_only": True,
        "symmetric_error": True,
        "caller_supplied": tail_input is not None,
        "W0_kernel_tail_uncertainty": True,
        "branch": tail_branches.get("D", tail_branches.get("G0")),
    } if tail_input is not None else None
    result = {
        "schema": "n1_kchi_feedback_exit_v2",
        "coordinate": "rho=4, x=2, t=1, real Z in [-1,1]",
        "exit_target": {
            "x": target_x, "t": target_t, "rho": c.mpf(4),
            "fixed": True, "independent": True,
            "not_last_target_cell_proxy": True,
            "source_cells": "canonical t=[i/8,(i+1)/8]",
        },
        "branches": branches,
        "raw_Q_branches": {key: value for key, value in branches.items()
                           if key in ("E01", "E10", "E11")},
        "E01_G0_Q_Delta": branches.get("E01"),
        "E10_D_Q_W0": branches.get("E10"),
        "E11_D_Q_Delta": branches.get("E11"),
        "E01": branches.get("E01"),
        "E10": branches.get("E10"),
        "E11": branches.get("E11"),
        "baseline_Q": baseline_q,
        "baseline": baseline,
        "Q_W0_by_source_cell": saved_q,
        "Q_DeltaW0_by_source_cell": delta_q,
        "Q_W0_kernel_tail_by_source_cell": tail_q,
        "symmetric_error": tail_summary,
        "symmetric_R_error": tail_summary,
        "error_branch_separate": True,
        "W0_kernel_tail_uncertainty": tail_summary,
        "symmetric_R_branches": tail_branches,
        "chi_jets": chi_rows,
        "source_power_basis": list(SOURCE_POWER_BASIS),
        "fields": list(FIELDS),
        "delta_fields": list(DELTA_FIELDS),
        "rows": dict(FIELD_ROWS),
        "D_row": D_ROW,
        "structural_zero_fields_D": list(ZERO_D_FIELDS),
        "matrix_orders": [0, 1],
        "ordinary_Z_orders": [0, 1],
        "input_W0_ordinary_Z_orders": [0, 1, 2],
        "ordinary_Z1_product_rule": PRODUCT_RULE_Z1,
        "ordinary_Z1_B1_W0_ZZ_included": True,
        "ordinary_Z1_chi_Z_included": True,
        "negative_kernel_K1": True,
        "green_kernel_masses_exact": True,
        "delta_kernel_masses_exact": True,
        "full_green_interval_kernel_masses": full_g0,
        "full_physical_interval_kernel_masses": full_d,
        "full_interval_exterior_masses": {
            "F1": "-chi/6", "d_x_F1": "-chi/3",
        },
        "delta_endpoint_source": "independent fixed x=2, not cell 7 endpoint",
        "saved_W0_reused": True,
        "W0_reconstructed": False,
        "W0_computed": False,
        "E00_G0_recomputed": False,
        "baseline_E00_recomputed": False,
        "ancestor_helpers_recomputed": False,
        "old_builder_called": False,
        "source_ids_retained": True,
        "formal_powers_retained": True,
        "signed_coefficients_retained": True,
        "baseline_and_tail_separate": True,
        "tail_input_supplied": tail_input is not None,
        "delta_W0_input_supplied": delta_input is not None,
        "raw_Q_only": True,
        "raw_Q_vs_potential_P": "separate",
        "potential_P_compensation_included": False,
        "potential_P_compensation_merged": False,
        "potential_compensation_status": "not included; raw-Q sectors only",
        "complete_six_component_E01": "E01" in branches,
        "complete_six_component_E10": False,
        "complete_six_component_E11": False,
        "full_W1_claimed": False,
        "complete_feedback": False,
        "Kchi_higher_words_included": False,
        "allword_tail_added": False,
        "m15_computed": False,
        "branch_selection": mode or "all",
        "new_branches_only": True,
    }
    return result


build = build_feedback_exit
build_core_feedback_exit = build_feedback_exit
build_endpoint = build_feedback_exit
build_fixed_exit = build_feedback_exit
build_x2 = build_feedback_exit
feedback_exit = build_feedback_exit
exit_x2 = build_feedback_exit
build_delta_w0 = build_feedback_exit
build_signed_delta_g5_exit = build_feedback_exit
signed_delta_g5_exit = build_feedback_exit


