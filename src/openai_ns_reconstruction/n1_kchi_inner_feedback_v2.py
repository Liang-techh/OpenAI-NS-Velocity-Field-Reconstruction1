"""Finite K0 feedback of the first inner Kchi correction.

The caller supplies two already enclosed, cellwise packets.  The first packet
is the signed ``Delta W0`` correction (only the F1 and ``d_x_F1`` columns are
used); the optional second packet is a symmetric uncertainty for the omitted
Kchi kernel words.  This leaf only forms the finite subterm

    G0 ( C0 DeltaW0 + B1 partial_Z DeltaW0 ).

No old W0 producer, matrix producer, Green builder, or Green helper is called.
The matrix and DeltaW0 packets are read as saved data and every product is
kept as an unmerged signed interval sector.  The uncertainty packet is
returned separately and is never presented as a signed field value.

The normalized radius is ``t=x/2`` and the eight canonical cells are
``[i/8,(i+1)/8]``.  A target cell i uses source cells j <= i; the diagonal
source cell is triangular.  The first triangular cell has the removable axis
continuation, so no ``log(0)`` or interval division by zero is formed.

Each output sector carries four formal powers in the order
``(hb, F0_anchor, Pstar_squared, Lambda)``.  Source IDs, error classes,
uncertainty flags, and explicit ``None`` source IDs are retained without
inventing provenance.
"""
from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any


SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
CELL_COUNT = 8
ORDERS = (0, 1, 2)
MATRIX_ORDERS = (0, 1)
OUTPUT_FIELDS = ("F1", "Uz1")
DELTA_FIELDS = ("F1", "d_x_F1")
DELTA_COLUMNS = {"F1": 1, "d_x_F1": 5}
FIELDS = (("F1", "5"), ("Uz1", "6"))
SOURCE_FIELDS = (("F1", "1"), ("d_x_F1", "5"))

EXPECTED_C0 = {
    "4,1", "5,1", "5,2", "5,3", "5,5",
    "6,1", "6,2", "6,3", "6,4", "6,6",
}
EXPECTED_B1 = {"5,1", "5,2", "5,3", "6,2", "6,3", "6,4"}
MATRIX_OUTPUT_ROWS = {"F1": 5, "Uz1": 6}

# These are ordinary derivatives.  In particular, Z1 has no Taylor-factorial
# conversion and the four terms below are the complete finite product rule.
PRODUCT_RULE_Z0 = (
    ("C0", "0", "0", 1, "C0 DeltaW0"),
    ("B1", "0", "1", 1, "B1 DeltaW0_Z"),
)
PRODUCT_RULE_Z1 = (
    ("C0", "1", "0", 1, "C0_Z DeltaW0"),
    ("C0", "0", "1", 1, "C0 DeltaW0_Z"),
    ("B1", "1", "1", 1, "B1_Z DeltaW0_Z"),
    ("B1", "0", "2", 1, "B1 DeltaW0_ZZ"),
)
PRODUCT_RULE = {"Z0": PRODUCT_RULE_Z0, "Z1": PRODUCT_RULE_Z1}

__all__ = [
    "SOURCE_POWER_BASIS", "CELL_COUNT", "ORDERS", "MATRIX_ORDERS",
    "OUTPUT_FIELDS", "FIELDS", "SOURCE_FIELDS", "DELTA_FIELDS", "DELTA_COLUMNS", "EXPECTED_C0",
    "EXPECTED_B1", "PRODUCT_RULE_Z0", "PRODUCT_RULE_Z1", "PRODUCT_RULE",
    "cell_kernel_masses", "target_radius_kernel_masses",
    "green_cell_kernel_masses", "kernel_masses", "build", "build_signed",
    "build_uncertainty", "build_from_cells", "build_g0_feedback",
    "build_kchi_inner_feedback",
]


def _is_interval(value: Any) -> bool:
    return hasattr(value, "_mpi_")


def _interval(c: Any, value: Any) -> Any:
    """Coerce a scalar/interval-like value in the caller's context."""

    if _is_interval(value):
        return value
    if hasattr(value, "a") and hasattr(value, "b"):
        return c.mpf([value.a, value.b])
    if isinstance(value, Mapping) and "lower" in value and "upper" in value:
        return c.mpf([value["lower"], value["upper"]])
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        if len(value) == 2:
            return c.mpf([value[0], value[1]])
    return c.mpf(value)


def _endpoint_bounds(c: Any, value: Any) -> tuple[Any, Any]:
    value = _interval(c, value)
    lo, hi = value._mpi_
    return c.make_mpf((lo, lo)), c.make_mpf((hi, hi))


def _same_interval(left: Any, right: Any) -> bool:
    try:
        return left._mpi_ == right._mpi_
    except AttributeError:
        return left == right


def _exact_cell(c: Any, index: int, cells: int = CELL_COUNT) -> Any:
    left = c.mpf(index) / cells
    right = c.mpf(index + 1) / cells
    return c.make_mpf((left._mpi_[0], right._mpi_[1]))


def _degenerate_endpoint(c: Any, value: Any, side: int) -> Any:
    value = _interval(c, value)
    raw = value._mpi_[side]
    return c.make_mpf((raw, raw))


def _nonnegative(c: Any, value: Any, label: str) -> Any:
    value = _interval(c, value)
    lo, hi = value._mpi_
    lo_i = c.make_mpf((lo, lo))
    hi_i = c.make_mpf((hi, hi))
    if hi_i < 0:
        raise ValueError(f"{label} has a negative upper bound")
    if lo_i < 0:
        lo = c.mpf(0)._mpi_[0]
    return c.make_mpf((lo, hi))


def _positive(c: Any, value: Any, label: str) -> None:
    lo, hi = _endpoint_bounds(c, value)
    if hi <= 0:
        raise ValueError(f"{label} must be positive")
    if lo <= 0:
        raise ValueError(f"{label} interval contains the axis")


def _primitive_uz(c: Any, target_t: Any, source_t: Any) -> Any:
    """Primitive of ``2*s^2 log(T/s)+s^2`` in normalized radius."""

    source_t = _interval(c, source_t)
    target_t = _interval(c, target_t)
    source_lo, _ = _endpoint_bounds(c, source_t)
    if source_lo == 0:
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
    """Return directed K0 source-cell masses in normalized radius.

    The six returned names retain the standard raw-state convention.  This
    leaf consumes only ``F1`` and ``Uz1``; the other four are exposed so a
    caller can inspect the independent kernel contract without invoking an
    older Green helper.
    """

    target = _interval(c, target_t)
    left = _interval(c, source_left_t)
    right = target if triangular else _interval(c, source_right_t)
    tlo, thi = _endpoint_bounds(c, target)
    alo, ahi = _endpoint_bounds(c, left)
    blo, bhi = _endpoint_bounds(c, right)
    if alo < 0 or ahi < alo or blo < alo or bhi < blo:
        raise ValueError("source endpoints must be ordered and nonnegative")
    if not triangular and bhi > thi:
        raise ValueError("source cell extends beyond target radius")
    if triangular and not _same_interval(right, target):
        raise ValueError("triangular source endpoint must be the target interval")

    # The whole first diagonal cell has a removable axis continuation.
    if triangular and alo == 0 and tlo == 0:
        f1 = target * target / 2
        uz1 = target * target
        p1 = 2 * target
        dxf1 = target / 2
        dxuz1 = target
    else:
        _positive(c, target, "target radius")
        t2 = target * target
        left2 = left * left
        right2 = right * right
        left4 = left2 * left2
        right4 = right2 * right2
        f1 = (right2 - left2) - (right4 - left4) / (2 * t2)
        uz1 = _primitive_uz(c, target, right) - _primitive_uz(c, target, left)
        p1 = 2 * (right - left)
        dxf1 = (right4 - left4) / (2 * target * t2)
        dxuz1 = (right2 - left2) / target

    f1 = _nonnegative(c, f1, "F1 Green mass")
    uz1 = _nonnegative(c, uz1, "Uz1 Green mass")
    p1 = _nonnegative(c, p1, "P1 Green mass")
    dxf1 = _nonnegative(c, dxf1, "d_x_F1 Green mass")
    dxuz1 = _nonnegative(c, dxuz1, "d_x_Uz1 Green mass")
    return {
        "F1": f1, "Uz1": uz1, "K1": -f1, "P1": p1,
        "d_x_F1": dxf1, "d_x_Uz1": dxuz1,
        "triangular": triangular,
        "axis_continuation": bool(triangular and alo == 0 and tlo == 0),
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
        c, _interval(c, target_x) / 2, _interval(c, source_left_x) / 2,
        _interval(c, source_right_x) / 2, triangular=triangular,
    )


green_cell_kernel_masses = cell_kernel_masses
kernel_masses = cell_kernel_masses


def _freeze(value: Any) -> Any:
    try:
        hash(value)
    except TypeError:
        return (type(value).__name__, repr(value))
    return (type(value).__name__, value)


def _unique(values: Sequence[Any]) -> list[Any]:
    result: list[Any] = []
    seen: set[Any] = set()
    for value in values:
        key = _freeze(value)
        if key not in seen:
            seen.add(key)
            result.append(value)
    return result


def _as_items(value: Any, *, keep_none: bool = False) -> list[Any]:
    if value is None:
        return [None] if keep_none else []
    if isinstance(value, (str, bytes, bytearray)):
        return [value]
    if isinstance(value, Sequence):
        return list(value)
    return [value]


def _origin_info(value: Any, *, _seen: set[int] | None = None) -> dict[str, Any]:
    """Collect provenance while retaining an explicit ``source_id=None``."""

    if _seen is None:
        _seen = set()
    marker = id(value)
    if marker in _seen:
        return {"source_ids": [], "error_classes": [], "uncertainty": False}
    if isinstance(value, Mapping):
        _seen.add(marker)
        ids: list[Any] = []
        if "source_ids" in value:
            ids.extend(_as_items(value["source_ids"], keep_none=True))
        elif "sources" in value:
            ids.extend(_as_items(value["sources"], keep_none=True))
        elif "source_id" in value:
            ids.extend(_as_items(value["source_id"], keep_none=True))
        elif "source" in value:
            ids.extend(_as_items(value["source"], keep_none=True))
        elif "source_cell" in value:
            ids.extend(_as_items(value["source_cell"], keep_none=True))
        elif "input_id" in value:
            ids.extend(_as_items(value["input_id"], keep_none=True))
        classes = list(_as_items(value.get("error_classes"), keep_none=False))
        if "error_classes" not in value:
            classes.extend(_as_items(value.get("uncertainty_class"), keep_none=False))
            classes.extend(_as_items(value.get("error_class"), keep_none=False))
        operation = str(value.get("operation", "")).lower()
        kind = str(value.get("error_kind", value.get("kind", ""))).lower()
        uncertain = bool(value.get("uncertainty", value.get("is_uncertainty", False)))
        uncertain = uncertain or any(
            token in operation or token in kind
            for token in ("symmetric", "cauchy", "uncertainty", "tail")
        )
        for key in ("origin", "metadata", "parents"):
            child = value.get(key)
            if child is None or child is value:
                continue
            if isinstance(child, Sequence) and not isinstance(child, (str, bytes, bytearray)):
                nested = [_origin_info(item, _seen=_seen) for item in child]
            else:
                nested = [_origin_info(child, _seen=_seen)]
            for info in nested:
                ids.extend(info["source_ids"])
                classes.extend(info["error_classes"])
                uncertain = uncertain or info["uncertainty"]
        return {
            "source_ids": _unique(ids),
            "error_classes": _unique(classes),
            "uncertainty": bool(uncertain),
        }
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        ids: list[Any] = []
        classes: list[Any] = []
        uncertain = False
        for item in value:
            info = _origin_info(item, _seen=_seen)
            ids.extend(info["source_ids"])
            classes.extend(info["error_classes"])
            uncertain = uncertain or info["uncertainty"]
        return {"source_ids": _unique(ids), "error_classes": _unique(classes),
                "uncertainty": uncertain}
    return {"source_ids": [], "error_classes": [], "uncertainty": False}


def _declared_source_id(value: Any) -> Any:
    """Return the first declared source ID; never replace explicit ``None``."""

    if not isinstance(value, Mapping):
        return None
    for key in ("source_id", "source_ids", "sources", "source", "input_id", "source_cell"):
        if key not in value:
            continue
        raw = value[key]
        if key == "source_ids" and isinstance(raw, Sequence) and not isinstance(raw, (str, bytes, bytearray)):
            return raw[0] if raw else None
        if key in {"sources"} and isinstance(raw, Sequence) and not isinstance(raw, (str, bytes, bytearray)):
            return raw[0] if raw else None
        return raw
    origin = value.get("origin")
    if isinstance(origin, Mapping):
        return _declared_source_id(origin)
    return None


def _powers(value: Mapping[str, Any], label: str) -> tuple[int, int, int, int]:
    if "powers" in value:
        raw = value["powers"]
    else:
        raw = (
            value.get("hb_power", 0),
            value.get("F0_power", 0),
            value.get("Pstar_squared_power", value.get("Pstar_power", 0)),
            value.get("Lambda_power", 0),
        )
    if (
        not isinstance(raw, Sequence) or isinstance(raw, (str, bytes, bytearray))
        or len(raw) != 4
        or any(not isinstance(item, int) or isinstance(item, bool) for item in raw)
    ):
        raise ValueError(f"{label} must contain four integer formal powers")
    return tuple(int(item) for item in raw)


def _term(c: Any, raw: Any, *, label: str, uncertainty: bool = False) -> dict[str, Any]:
    if isinstance(raw, Mapping):
        if "coefficient" not in raw:
            # A compact bound spelling is useful for a separate tail packet.
            if "bound" in raw:
                raw = dict(raw, coefficient=raw["bound"])
            elif "radius" in raw:
                raw = dict(raw, coefficient=raw["radius"])
            else:
                raise ValueError(f"{label} requires coefficient")
        coefficient = raw["coefficient"]
        powers = _powers(raw, label)
        extras = dict(raw)
    else:
        coefficient = raw
        powers = (0, 0, 0, 0)
        extras = {}
    coefficient = _interval(c, coefficient)
    info = _origin_info(raw if isinstance(raw, Mapping) else {})
    origin_value = raw.get("origin") if isinstance(raw, Mapping) else None
    origin = dict(origin_value) if isinstance(origin_value, Mapping) else {}
    uncertain = bool(uncertainty or info["uncertainty"])
    origin.update({
        "source_ids": list(info["source_ids"]),
        "error_classes": list(info["error_classes"]),
        "uncertainty": uncertain,
    })
    result = dict(extras)
    result.update({
        "powers": powers,
        "hb_power": powers[0], "F0_power": powers[1],
        "Pstar_squared_power": powers[2], "Lambda_power": powers[3],
        "coefficient": coefficient,
        "origin": origin,
        "source_ids": list(info["source_ids"]),
        "error_classes": list(info["error_classes"]),
        "uncertainty": uncertain,
    })
    return result


def _term_list(c: Any, value: Any, *, label: str, uncertainty: bool = False) -> list[dict[str, Any]]:
    if value is None:
        return []
    if isinstance(value, Mapping) and ("coefficient" in value or "powers" in value or "bound" in value or "radius" in value):
        values = [value]
    elif isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        values = list(value)
    else:
        values = [value]
    return [_term(c, item, label=f"{label}[{index}]", uncertainty=uncertainty)
            for index, item in enumerate(values)]


def _order_key(key: Any, label: str) -> int:
    if isinstance(key, bool):
        raise ValueError(f"boolean is not a Z order in {label}")
    if isinstance(key, int):
        return key
    if isinstance(key, str):
        text = key.upper().removeprefix("Z")
        try:
            return int(text)
        except ValueError as exc:
            raise ValueError(f"invalid Z order {key!r} in {label}") from exc
    raise ValueError(f"invalid Z order {key!r} in {label}")


def _lookup(mapping: Mapping[Any, Any], key: Any) -> Any:
    if key in mapping:
        return mapping[key]
    other = str(key) if not isinstance(key, str) else (int(key) if key.isdigit() else key)
    if other in mapping:
        return mapping[other]
    raise KeyError(key)


def _ordered(c: Any, value: Any, *, label: str, uncertainty: bool) -> dict[int, list[dict[str, Any]]]:
    """Normalize one F1/d_x_F1 value through ordinary Z2."""

    if isinstance(value, Mapping):
        for wrapper in ("ordinary_Z", "source_jets", "source_rows", "jets", "rows"):
            if wrapper in value:
                value = value[wrapper]
                break
    if isinstance(value, Mapping):
        parsed: dict[int, Any] = {}
        for key, item in value.items():
            try:
                order = _order_key(key, label)
            except ValueError:
                continue
            if order in parsed:
                raise ValueError(f"duplicate Z order {order} in {label}")
            parsed[order] = item
    elif isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        # A term mapping is a one-term list, not an order vector.
        if len(value) == 3 and not all(isinstance(item, Mapping) and "coefficient" in item for item in value):
            parsed = {index: item for index, item in enumerate(value)}
        else:
            parsed = {0: value}
    else:
        parsed = {0: value}
    if set(parsed) != set(ORDERS):
        raise ValueError(f"{label} must provide exactly Z0, Z1, and Z2")
    return {
        order: _term_list(c, parsed[order], label=f"{label}[Z{order}]", uncertainty=uncertainty)
        for order in ORDERS
    }


def _empty_rows() -> dict[int, dict[int, list[dict[str, Any]]]]:
    return {column: {order: [] for order in ORDERS} for column in DELTA_COLUMNS.values()}


def _assign_field_rows(c: Any, rows: dict[int, dict[int, list[dict[str, Any]]]], field: str,
                       value: Any, *, label: str, uncertainty: bool) -> None:
    column = DELTA_COLUMNS[field]
    ordered = _ordered(c, value, label=label, uncertainty=uncertainty)
    for order in ORDERS:
        rows[column][order].extend(ordered[order])


def _collect_delta(c: Any, candidate: Any, rows: dict[int, dict[int, list[dict[str, Any]]]],
                   *, label: str, uncertainty: bool) -> None:
    """Read sparse field/column/order forms without deriving absent rows."""

    if candidate is None:
        return
    if isinstance(candidate, Mapping):
        # Cell-level wrappers used by both the first-correction leaf and saved
        # error packets.  The first matching wrapper is structural only.
        for wrapper in ("signed", "finite", "signed_first_correction", "delta_W0",
                        "DeltaW0", "W0_correction", "first_correction",
                        "W0_by_order", "W0", "source_jets", "source_rows",
                        "ordinary_Z"):
            if wrapper in candidate and candidate[wrapper] is not candidate:
                _collect_delta(c, candidate[wrapper], rows, label=f"{label}.{wrapper}", uncertainty=uncertainty)
                return
        for field in DELTA_FIELDS:
            if field in candidate:
                _assign_field_rows(c, rows, field, candidate[field], label=f"{label}.{field}", uncertainty=uncertainty)
        # Explicit order wrappers (W0_Z0/W0_Z1/W0_Z2) contain fields or
        # six-column mappings.  This form is emitted by first-correction cells.
        found_prefixed = False
        for order in ORDERS:
            for key in (f"W0_Z{order}", f"DeltaW0_Z{order}", f"delta_W0_Z{order}"):
                if key not in candidate:
                    continue
                found_prefixed = True
                value = candidate[key]
                if isinstance(value, Mapping):
                    for field in DELTA_FIELDS:
                        if field in value:
                            rows[DELTA_COLUMNS[field]][order].extend(
                                _term_list(c, value[field], label=f"{label}.{key}.{field}", uncertainty=uncertainty)
                            )
                    for column in DELTA_COLUMNS.values():
                        try:
                            terms = _lookup(value, str(column))
                        except KeyError:
                            continue
                        rows[column][order].extend(_term_list(c, terms, label=f"{label}.{key}.{column}", uncertainty=uncertainty))
                else:
                    raise TypeError(f"{label}.{key} must be a field/column mapping")
        if found_prefixed:
            return
        # Ordinary order-first mapping: order -> field or source column.
        parsed_orders: dict[int, Any] = {}
        for key, value in candidate.items():
            try:
                parsed = _order_key(key, label)
            except ValueError:
                continue
            if parsed in ORDERS:
                parsed_orders[parsed] = value
        if parsed_orders:
            if set(parsed_orders) != set(ORDERS):
                raise ValueError(f"{label} must provide exactly Z0, Z1, and Z2")
            for order, value in parsed_orders.items():
                if not isinstance(value, Mapping):
                    raise TypeError(f"{label}[Z{order}] must be a field/column mapping")
                for field in DELTA_FIELDS:
                    if field in value:
                        rows[DELTA_COLUMNS[field]][order].extend(
                            _term_list(c, value[field], label=f"{label}[Z{order}].{field}", uncertainty=uncertainty)
                        )
                for column in DELTA_COLUMNS.values():
                    try:
                        terms = _lookup(value, str(column))
                    except KeyError:
                        continue
                    rows[column][order].extend(_term_list(c, terms, label=f"{label}[Z{order}].{column}", uncertainty=uncertainty))
            return
        # A single field value may be wrapped as ``{"value": ...}``.
        if "value" in candidate:
            _collect_delta(c, candidate["value"], rows, label=f"{label}.value", uncertainty=uncertainty)
            return
        return
    if isinstance(candidate, Sequence) and not isinstance(candidate, (str, bytes, bytearray)):
        if len(candidate) == CELL_COUNT:
            raise TypeError(f"{label} needs one cell, not the complete cell sequence")
        raise TypeError(f"{label} must be a field/column mapping")
    raise TypeError(f"{label} must be a field/column mapping")


def _cell_wrapper(cell: Mapping[str, Any], *, uncertainty: bool) -> Any:
    if uncertainty:
        for key in ("uncertainty", "uncertainty_W0", "W0_uncertainty",
                    "kernel_tail", "kernel_tail_uncertainty", "tail",
                    "error", "errors", "W0_error", "W0_tail"):
            if key in cell:
                return cell[key]
        # A separate uncertainty packet may use exactly the same sparse
        # F1/d_x_F1 or order-first shape as the signed packet, with the
        # uncertainty meaning carried by this call rather than a wrapper.
        if any(key in cell for key in (*DELTA_FIELDS, "W0_Z0", "W0_Z1", "W0_Z2", 0, 1, 2, "0", "1", "2")):
            return cell
        return None
    for key in ("signed", "finite", "signed_first_correction", "delta_W0",
                "DeltaW0", "W0_correction", "first_correction", "W0_by_order",
                "W0_ordinary_Z", "W0", "source_jets", "source_rows", "source_terms"):
        if key in cell:
            return cell[key]
    return cell


def _delta_rows(c: Any, cell: Any, index: int, *, uncertainty: bool) -> dict[int, dict[int, list[dict[str, Any]]]]:
    if not isinstance(cell, Mapping):
        raise TypeError(f"DeltaW0 cell {index} must be a mapping")
    rows = _empty_rows()
    candidate = _cell_wrapper(cell, uncertainty=uncertainty)
    if candidate is None:
        return rows
    _collect_delta(c, candidate, rows, label=f"DeltaW0[{index}]", uncertainty=uncertainty)
    return rows


def _cell_t(c: Any, index: int, cell: Mapping[str, Any] | None, *, label: str) -> Any:
    expected = _exact_cell(c, index)
    if cell is None:
        return expected
    candidate = None
    coordinate = "t"
    if "t" in cell:
        candidate = cell["t"]
    elif "x" in cell:
        candidate = _interval(c, cell["x"]) / 2
        coordinate = "x"
    elif "target_x" in cell:
        candidate = _interval(c, cell["target_x"]) / 2
        coordinate = "target_x"
    elif "radius" in cell:
        candidate = _interval(c, cell["radius"])
    if candidate is None:
        return expected
    actual = _interval(c, candidate)
    if not _same_interval(actual, expected):
        raise ValueError(f"{label} cell {index} has a noncanonical {coordinate} interval")
    if "rho" in cell and not _same_interval(_interval(c, cell["rho"]), 4 * actual * actual):
        raise ValueError(f"{label} cell {index} has inconsistent rho=4*t^2")
    return actual


def _as_cells(value: Any, *, label: str) -> list[Mapping[str, Any]]:
    if isinstance(value, Mapping):
        for key in ("target_cells", "outer_cells", "source_cells", "matrix_cells",
                    "saved_matrix_cells", "cells"):
            if key in value:
                value = value[key]
                break
        else:
            try:
                return [_lookup(value, index) for index in range(CELL_COUNT)]
            except KeyError as exc:
                raise ValueError(f"{label} must contain eight canonical cells") from exc
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes, bytearray)):
        raise TypeError(f"{label} must be an eight-cell sequence or mapping")
    if len(value) != CELL_COUNT:
        raise ValueError(f"{label} must contain eight canonical cells")
    if any(not isinstance(item, Mapping) for item in value):
        raise TypeError(f"{label} cells must be mappings")
    return list(value)


def _matrix_sector(cell: Mapping[str, Any], name: str, order: int, key: str) -> list[Any]:
    sector = _lookup(cell, name)
    values = _lookup(_lookup(sector, str(order)), key)
    if isinstance(values, Mapping):
        return [values]
    if not isinstance(values, (list, tuple)):
        raise TypeError("matrix source sectors must be a list")
    return list(values)


def _validate_matrix(c: Any, cell: Mapping[str, Any], index: int) -> None:
    _cell_t(c, index, cell, label="matrix")
    for name, expected in (("C0", EXPECTED_C0), ("B1", EXPECTED_B1)):
        try:
            sector = _lookup(cell, name)
        except (KeyError, TypeError) as exc:
            raise ValueError(f"matrix cell {index} is missing {name}") from exc
        for order in MATRIX_ORDERS:
            try:
                entries = _lookup(sector, str(order))
            except (KeyError, TypeError) as exc:
                raise ValueError(f"matrix cell {index} is missing {name} order {order}") from exc
            if {str(key) for key in entries} != expected:
                raise ValueError(f"{name} support changed at cell {index}, order {order}")
            for key in expected:
                for term_index, raw in enumerate(_matrix_sector(cell, name, order, key)):
                    if not isinstance(raw, Mapping):
                        raise TypeError(f"{name}[{order}][{key}][{term_index}] must be a term mapping")
                    _powers(raw, f"{name}[{order}][{key}][{term_index}]")
                    if "coefficient" not in raw:
                        raise ValueError(f"{name}[{order}][{key}][{term_index}] requires coefficient")


def _normalise_matrix_term(c: Any, raw: Mapping[str, Any], *, name: str, order: int,
                           key: str, term_index: int) -> dict[str, Any]:
    result = _term(c, raw, label=f"{name}[{order}][{key}][{term_index}]")
    result["matrix_name"] = name
    result["matrix_order"] = str(order)
    result["matrix_key"] = key
    result["matrix_term_index"] = term_index
    result["matrix_source_id"] = _declared_source_id(raw)
    result["matrix_source_cell"] = raw.get("source_cell")
    return result


def _symmetrise(c: Any, value: Any) -> Any:
    value = _interval(c, value)
    # abs() is interval-safe in the caller context; use its upper endpoint.
    radius = abs(value)
    radius_hi = radius._mpi_[1]
    radius_point = c.make_mpf((radius_hi, radius_hi))
    negative = -radius_point
    return c.make_mpf((negative._mpi_[0], radius_hi))


def _product_term(
    c: Any,
    matrix_term: Mapping[str, Any],
    delta_term: Mapping[str, Any],
    multiplier: int,
    origin_label: str,
    *,
    field: str,
    target_index: int,
    source_index: int,
    source_column: int,
    delta_order: int,
    delta_term_index: int,
    mass: Any,
    uncertainty: bool,
) -> dict[str, Any]:
    matrix_powers = matrix_term["powers"]
    delta_powers = delta_term["powers"]
    powers = tuple(a + b for a, b in zip(matrix_powers, delta_powers))
    info = _origin_info((matrix_term, delta_term))
    coefficient = multiplier * matrix_term["coefficient"] * delta_term["coefficient"] * mass
    if uncertainty:
        coefficient = _symmetrise(c, coefficient)
    source_id = _declared_source_id(delta_term)
    matrix_source_id = matrix_term.get("matrix_source_id")
    origin = {
        "operation": "kchi_inner_feedback_uncertainty" if uncertainty else "kchi_inner_feedback_product",
        "product": origin_label,
        "origin_label": origin_label,
        "source_ids": list(info["source_ids"]),
        "error_classes": list(info["error_classes"]),
        "uncertainty": bool(uncertainty or info["uncertainty"]),
        "symbolic_cancellation_claimed": False,
        "parents": (
            {"source_ids": list(_origin_info(matrix_term)["source_ids"])},
            {"source_ids": list(_origin_info(delta_term)["source_ids"])},
        ),
        "finite_subterm": True,
        "kernel_tail_not_a_field": bool(uncertainty),
    }
    return {
        "powers": powers,
        "hb_power": powers[0], "F0_power": powers[1],
        "Pstar_squared_power": powers[2], "Lambda_power": powers[3],
        "coefficient": coefficient,
        "origin": origin,
        "source_ids": list(info["source_ids"]),
        "error_classes": list(info["error_classes"]),
        "uncertainty": bool(uncertainty or info["uncertainty"]),
        "signed_interval": not uncertainty,
        "symmetric_uncertainty": bool(uncertainty),
        "finite_subterm": True,
        "kernel_tail_not_a_field": bool(uncertainty),
        "target_cell": target_index,
        "source_cell": source_index,
        "source_index": source_index,
        "source_column": source_column,
        "source_cell_id": delta_term.get("source_cell"),
        "source_term_id": source_id,
        "source_id": source_id,
        "delta_W0_source_id": source_id,
        "delta_W0_source_cell": delta_term.get("source_cell"),
        "delta_W0_order": str(delta_order),
        "delta_W0_term_index": delta_term_index,
        "delta_W0_provenance": delta_term.get("origin", {}).get("operation") if isinstance(delta_term.get("origin"), Mapping) else None,
        "matrix_source_id": matrix_source_id,
        "matrix_source_cell": matrix_term.get("matrix_source_cell"),
        "matrix_name": matrix_term.get("matrix_name"),
        "matrix_order": matrix_term.get("matrix_order"),
        "matrix_key": matrix_term.get("matrix_key"),
        "matrix_term_index": matrix_term.get("matrix_term_index"),
        "W0_source_id": source_id,
        "W0_source_cell": delta_term.get("source_cell"),
        "W0_order": str(delta_order),
        "W0_term_index": delta_term_index,
        "row": int(str(matrix_term["matrix_key"]).split(",")[0]),
        "column": int(str(matrix_term["matrix_key"]).split(",")[1]),
        "origin_label": origin_label,
        "leibniz_multiplier": multiplier,
        "kernel_mass": mass,
        "kernel_mass_field": field,
        "principal_kernel": "K0",
    }


def _feedback_terms(
    c: Any,
    matrix_cell: Mapping[str, Any],
    delta_rows: Mapping[int, Mapping[int, Sequence[Mapping[str, Any]]]],
    masses: Mapping[str, Any],
    *,
    target_index: int,
    source_index: int,
    uncertainty: bool,
) -> dict[str, dict[str, list[dict[str, Any]]]]:
    output = {"Z0": {field: [] for field in OUTPUT_FIELDS},
              "Z1": {field: [] for field in OUTPUT_FIELDS}}
    rules = {"Z0": PRODUCT_RULE_Z0, "Z1": PRODUCT_RULE_Z1}
    for order_name, rule_set in rules.items():
        for matrix_name, matrix_order_text, delta_order_text, multiplier, origin_label in rule_set:
            matrix_order = int(matrix_order_text)
            delta_order = int(delta_order_text)
            for key in sorted(EXPECTED_C0 if matrix_name == "C0" else EXPECTED_B1):
                row_text, column_text = key.split(",")
                row = int(row_text)
                column = int(column_text)
                if row not in MATRIX_OUTPUT_ROWS.values() or column not in delta_rows:
                    continue
                matrix_values = _matrix_sector(matrix_cell, matrix_name, matrix_order, key)
                matrix_terms = [
                    _normalise_matrix_term(c, raw, name=matrix_name,
                                           order=matrix_order, key=key,
                                           term_index=index)
                    for index, raw in enumerate(matrix_values)
                ]
                delta_terms = delta_rows[column][delta_order]
                field = "F1" if row == 5 else "Uz1"
                for matrix_term in matrix_terms:
                    for delta_term_index, delta_term in enumerate(delta_terms):
                        output[order_name][field].append(
                            _product_term(
                                c, matrix_term, delta_term, multiplier, origin_label,
                                field=field, target_index=target_index,
                                source_index=source_index, source_column=column,
                                delta_order=delta_order,
                                delta_term_index=delta_term_index,
                                mass=masses[field], uncertainty=uncertainty,
                            )
                        )
    return output


def _empty_feedback() -> dict[str, dict[str, list[dict[str, Any]]]]:
    return {"Z0": {field: [] for field in OUTPUT_FIELDS},
            "Z1": {field: [] for field in OUTPUT_FIELDS}}


def _target_packet(
    c: Any,
    matrix_cells: Sequence[Mapping[str, Any]],
    delta_cells: Sequence[Mapping[str, Any]],
    *,
    uncertainty: bool,
) -> list[dict[str, Any]]:
    targets: list[dict[str, Any]] = []
    for target_index in range(CELL_COUNT):
        target_t = _cell_t(c, target_index, matrix_cells[target_index], label="matrix")
        aggregate = _empty_feedback()
        by_source: dict[str, Any] = {}
        for source_index in range(target_index + 1):
            source_t = _cell_t(c, source_index, matrix_cells[source_index], label="matrix")
            source_left = _degenerate_endpoint(c, source_t, 0)
            source_right = (
                target_t if source_index == target_index
                else _degenerate_endpoint(c, source_t, 1)
            )
            masses = cell_kernel_masses(
                c, target_t, source_left, source_right,
                triangular=source_index == target_index,
            )
            rows = _delta_rows(c, delta_cells[source_index], source_index,
                               uncertainty=uncertainty)
            feedback = _feedback_terms(
                c, matrix_cells[source_index], rows, masses,
                target_index=target_index, source_index=source_index,
                uncertainty=uncertainty,
            )
            by_source[str(source_index)] = {
                "source_index": source_index,
                "target_index": target_index,
                "target_t": target_t,
                "source_t": source_t,
                "triangular": source_index == target_index,
                "kernel_masses": masses,
                "feedback": feedback,
                "G0_feedback": feedback,
                "finite_subterm": True,
                "kernel_tail_not_a_field": bool(uncertainty),
            }
            for order_name in ("Z0", "Z1"):
                for field in OUTPUT_FIELDS:
                    aggregate[order_name][field].extend(feedback[order_name][field])
        target = {
            "index": target_index,
            "target_index": target_index,
            "t": target_t,
            "x": 2 * target_t,
            "rho": 4 * target_t * target_t,
            "source_indices": list(range(target_index + 1)),
            "source_j_le_target_i": True,
            "triangular_diagonal_source": target_index,
            "axis_endpoint_continuation": target_index == 0,
            "finite_subterm": True,
            "kernel_tail_not_a_field": bool(uncertainty),
            "signed_interval_enclosure": not uncertainty,
            "symmetric_uncertainty_enclosure": bool(uncertainty),
            "G0_feedback": aggregate,
            "G0_by_order": aggregate,
            "G0_feedback_by_order": aggregate,
            "G0_feedback_by_source_cell": by_source,
            "G0_by_order_by_source_cell": by_source,
            "feedback": aggregate,
            "feedback_by_source_cell": by_source,
            "W1_by_order": aggregate,
            "W1_by_order_by_source_cell": by_source,
            "F1": {order: aggregate[order]["F1"] for order in ("Z0", "Z1")},
            "Uz1": {order: aggregate[order]["Uz1"] for order in ("Z0", "Z1")},
            "F1_by_order": {order: aggregate[order]["F1"] for order in ("Z0", "Z1")},
            "Uz1_by_order": {order: aggregate[order]["Uz1"] for order in ("Z0", "Z1")},
            "kernel_masses_by_source_cell": {
                key: value["kernel_masses"] for key, value in by_source.items()
            },
        }
        target["G0_Z0"] = aggregate["Z0"]
        target["G0_Z1"] = aggregate["Z1"]
        target["feedback_Z0"] = aggregate["Z0"]
        target["feedback_Z1"] = aggregate["Z1"]
        target["inner_feedback_Z0"] = aggregate["Z0"]
        target["inner_feedback_Z1"] = aggregate["Z1"]
        target["W1_Z0"] = aggregate["Z0"]
        target["W1_Z1"] = aggregate["Z1"]
        targets.append(target)
    return targets


def _empty_cells(c: Any, matrix_cells: Sequence[Mapping[str, Any]], *, uncertainty: bool) -> list[dict[str, Any]]:
    empty = [{"x": _cell_t(c, index, matrix_cells[index], label="matrix") * 2,
              "t": _cell_t(c, index, matrix_cells[index], label="matrix"),
              "rho": 4 * _cell_t(c, index, matrix_cells[index], label="matrix") ** 2,
              "signed_interval_enclosure": not uncertainty,
              "symmetric_uncertainty_enclosure": bool(uncertainty),
              "finite_subterm": True,
              "kernel_tail_not_a_field": bool(uncertainty)} for index in range(CELL_COUNT)]
    for index, target in enumerate(empty):
        aggregate = _empty_feedback()
        target.update({
            "index": index, "target_index": index,
            "source_indices": list(range(index + 1)),
            "source_j_le_target_i": True,
            "triangular_diagonal_source": index,
            "axis_endpoint_continuation": index == 0,
            "G0_feedback": aggregate, "G0_by_order": aggregate,
            "feedback": aggregate,
            "F1": {order: aggregate[order]["F1"] for order in ("Z0", "Z1")},
            "Uz1": {order: aggregate[order]["Uz1"] for order in ("Z0", "Z1")},
        })
    return empty


def _split_packets(value: Any, uncertainty_value: Any) -> tuple[Any, Any]:
    signed = value
    tail = uncertainty_value
    if isinstance(value, Mapping):
        for key in ("signed", "finite", "signed_first_correction", "delta_W0", "DeltaW0"):
            if key in value:
                signed = value[key]
                break
        if tail is None:
            for key in ("uncertainty", "kernel_tail", "kernel_tail_uncertainty", "tail", "W0_error"):
                if key in value:
                    tail = value[key]
                    break
    return signed, tail


def build(
    c: Any,
    delta_w0_cells: Any = None,
    matrix_cells: Any = None,
    *,
    w0_correction_cells: Any = None,
    signed_delta_w0: Any = None,
    uncertainty_cells: Any = None,
    kernel_tail_cells: Any = None,
    kernel_tail: Any = None,
    w0_uncertainty_cells: Any = None,
    uncertainty: Any = None,
) -> dict[str, Any]:
    """Build signed and separately reported uncertainty K0 feedback sectors.

    ``delta_w0_cells`` (or ``w0_correction_cells``) is an eight-cell sequence
    of saved first-correction packets.  A cell may use the first-correction
    form ``F1/d_x_F1 -> Z0/Z1/Z2`` or an order-first six-column form where
    only columns 1 and 5 are populated.  ``matrix_cells`` is the saved
    eight-cell C0/B1 packet with orders 0 and 1.  ``uncertainty_cells`` and
    its aliases have the same sparse shape; each of their coefficients is
    widened to a symmetric interval and remains outside the signed branch.
    """

    if not hasattr(c, "mpf") or not hasattr(c, "make_mpf"):
        raise TypeError("c must be the caller's MPIntervalContext")
    if matrix_cells is None:
        raise ValueError("saved matrix_cells are required")
    if delta_w0_cells is None:
        delta_w0_cells = w0_correction_cells if w0_correction_cells is not None else signed_delta_w0
    signed_packet, split_tail = _split_packets(delta_w0_cells, uncertainty_cells)
    if split_tail is not None and kernel_tail_cells is None and kernel_tail is None and w0_uncertainty_cells is None:
        kernel_tail_cells = split_tail
    if kernel_tail_cells is None:
        kernel_tail_cells = kernel_tail if kernel_tail is not None else w0_uncertainty_cells
    if kernel_tail_cells is None:
        kernel_tail_cells = uncertainty

    matrices = _as_cells(matrix_cells, label="matrix_cells")
    signed_cells = _as_cells(signed_packet, label="delta_w0_cells") if signed_packet is not None else [dict() for _ in range(CELL_COUNT)]
    tail_cells = _as_cells(kernel_tail_cells, label="uncertainty_cells") if kernel_tail_cells is not None else [dict() for _ in range(CELL_COUNT)]
    for index, cell in enumerate(matrices):
        _validate_matrix(c, cell, index)
    for index, cell in enumerate(signed_cells):
        _cell_t(c, index, cell, label="DeltaW0")
    for index, cell in enumerate(tail_cells):
        _cell_t(c, index, cell, label="uncertainty")

    signed_targets = _target_packet(c, matrices, signed_cells, uncertainty=False)
    uncertainty_present = kernel_tail_cells is not None
    uncertainty_targets = (
        _target_packet(c, matrices, tail_cells, uncertainty=True)
        if uncertainty_present else _empty_cells(c, matrices, uncertainty=True)
    )

    signed_by_order = {
        "Z0": {str(target["index"]): target["G0_Z0"] for target in signed_targets},
        "Z1": {str(target["index"]): target["G0_Z1"] for target in signed_targets},
    }
    uncertainty_by_order = {
        "Z0": {str(target["index"]): target["G0_Z0"] for target in uncertainty_targets},
        "Z1": {str(target["index"]): target["G0_Z1"] for target in uncertainty_targets},
    }
    return {
        "schema": "n1_kchi_inner_feedback_v2",
        "scope": "finite K0 feedback of signed first Kchi DeltaW0; separate symmetric kernel uncertainty",
        "target_cells": signed_targets,
        "signed_target_cells": signed_targets,
        "G0_feedback_target_cells": signed_targets,
        "inner_feedback_target_cells": signed_targets,
        "G0_feedback_Z0": signed_by_order["Z0"],
        "G0_feedback_Z1": signed_by_order["Z1"],
        "inner_feedback_Z0": signed_by_order["Z0"],
        "inner_feedback_Z1": signed_by_order["Z1"],
        "W1_Z0_target_cells": signed_by_order["Z0"],
        "W1_Z1_target_cells": signed_by_order["Z1"],
        "signed_G0_feedback_Z0": signed_by_order["Z0"],
        "signed_G0_feedback_Z1": signed_by_order["Z1"],
        "uncertainty_target_cells": uncertainty_targets,
        "kernel_tail_target_cells": uncertainty_targets,
        "uncertainty_feedback_target_cells": uncertainty_targets,
        "kernel_tail_feedback_target_cells": uncertainty_targets,
        "uncertainty_G0_feedback_Z0": uncertainty_by_order["Z0"],
        "uncertainty_G0_feedback_Z1": uncertainty_by_order["Z1"],
        "kernel_tail_G0_feedback_Z0": uncertainty_by_order["Z0"],
        "kernel_tail_G0_feedback_Z1": uncertainty_by_order["Z1"],
        "saved_delta_W0_cells": signed_cells,
        "saved_uncertainty_cells": tail_cells if uncertainty_present else None,
        "saved_matrix_cells": matrices,
        "saved_delta_W0_reused": True,
        "saved_matrix_reused": True,
        "ordinary_Z_orders": [0, 1],
        "input_delta_W0_ordinary_Z_orders": [0, 1, 2],
        "matrix_ordinary_Z_orders": [0, 1],
        "output_fields": list(OUTPUT_FIELDS),
        "delta_W0_nonzero_fields": list(DELTA_FIELDS),
        "delta_W0_columns": dict(DELTA_COLUMNS),
        "source_power_basis": list(SOURCE_POWER_BASIS),
        "source_j_le_target_i": True,
        "triangular_diagonal_source": True,
        "axis_endpoint_continuation": True,
        "axis_log_division_avoided": True,
        "green_kernel_mass_definition": "exact directed K0 source-cell integral; diagonal source cell is triangular",
        "green_kernel_masses_exact": True,
        "ordinary_Z0_product_rule": PRODUCT_RULE_Z0,
        "ordinary_Z1_product_rule": PRODUCT_RULE_Z1,
        "signed_first_kchi_input_only": True,
        "signed_coefficients_retained": True,
        "symmetric_kernel_uncertainty_input_separate": uncertainty_present,
        "uncertainty_output_is_not_a_signed_field": True,
        "uncertainty_is_not_function_value": True,
        "uncertainty_branch_returned_separately": True,
        "finite_subterm": True,
        "finite_subterm_only": True,
        "kernel_tail_merged": False,
        "kernel_tail_added_to_signed_field": False,
        "tail_added": False,
        "tail_included": False,
        "allword_tail_added": False,
        "allword_tail_included": False,
        "Kchi_higher_words_included": False,
        "W0_reconstructed": False,
        "W0_ancestor_rerun": False,
        "W0_reweighted": False,
        "matrix_reconstructed": False,
        "matrix_ancestor_rerun": False,
        "old_builder_called": False,
        "old_green_helper_called": False,
        "no_source_cells_collapsed": True,
        "formal_powers_retained": True,
        "source_ids_retained": True,
        "null_provenance_retained": True,
        "error_classes_retained": True,
        "uncertainty_retained": True,
        "signed_field_claimed": False,
        "physical_field_claimed": False,
    }


def build_signed(c: Any, delta_w0_cells: Any, matrix_cells: Any, **kwargs: Any) -> dict[str, Any]:
    """Return the full schema with only the finite signed branch populated."""

    kwargs.pop("uncertainty_cells", None)
    kwargs.pop("kernel_tail_cells", None)
    kwargs.pop("kernel_tail", None)
    kwargs.pop("w0_uncertainty_cells", None)
    return build(c, delta_w0_cells, matrix_cells, **kwargs)


def build_uncertainty(c: Any, uncertainty_cells: Any, matrix_cells: Any, **kwargs: Any) -> dict[str, Any]:
    """Return the full schema with the symmetric uncertainty branch populated."""

    kwargs["uncertainty_cells"] = uncertainty_cells
    result = build(c, None, matrix_cells, **kwargs)
    return result


build_from_cells = build
build_g0_feedback = build
build_kchi_inner_feedback = build
