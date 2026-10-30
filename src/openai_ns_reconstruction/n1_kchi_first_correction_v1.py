"""Cellwise first ``Kchi`` correction for the original n1 radial solve.

This leaf is deliberately narrower than the original regular-series owner.  It
only applies the first term

``(-2*chi) * T K0``

to caller-supplied source cells.  Here ``K0(x,s) = s*(1-s**2/x**2)/2`` is the
swirl Green kernel and ``T`` is the same radial inverse used by the original
regular series.  No higher ``Kchi`` word, feedback word, or all-word tail is
constructed.  Source coefficients are already enclosed by the caller; this
module never replays an old source producer or recomputes an old word.

The public cell-mass function returns directed masses of ``T K0`` for
``F1`` and ``d_x_F1``.  :func:`build` multiplies those masses by the signed
factor ``-2*chi`` and combines the caller's ordinary source jets through
orders ``Z0/Z1/Z2``.  The second ordinary derivative uses the Leibniz factor
``2`` in ``2*chi_Z*source_Z``.  Four integer formal source powers
``(hb, F0_anchor, Pstar_squared, Lambda)`` and source identities are retained
on every output term.

For a source interval ``[a,b]`` and target radius ``x`` the two bare kernels
are

``K1(x,s) = s*x**2/16 - s**3*log(x/s)/4 - s**5/(16*x**2)``

and ``d_x K1(x,s) = s*(x**2-s**2)**2/(8*x**3)``.  Their cell primitives are
evaluated directly.  On a diagonal cell ``b=x`` and on the first triangular
cell ``a=0`` polynomial continuations are used, so no ``log(0)`` or division
by zero is formed.  All intervals remain directed and signed; no midpoint or
absolute-value replacement is made.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from typing import Any


SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
# Eight cells of width 1/4 cover the physical radius x in [0, 2].
CELL_COUNT = 8
ORDERS = (0, 1, 2)
FIELDS = ("F1", "d_x_F1")

# Ordinary derivatives, not Taylor coefficients.  The middle Z2 product has
# the usual binomial factor two.
PRODUCT_RULE_Z0 = ((0, 0, 1),)
PRODUCT_RULE_Z1 = ((1, 0, 1), (0, 1, 1))
PRODUCT_RULE_Z2 = ((2, 0, 1), (1, 1, 2), (0, 2, 1))
PRODUCT_RULE = {"Z0": PRODUCT_RULE_Z0, "Z1": PRODUCT_RULE_Z1,
                "Z2": PRODUCT_RULE_Z2}


__all__ = [
    "SOURCE_POWER_BASIS", "CELL_COUNT", "ORDERS", "FIELDS",
    "PRODUCT_RULE_Z0", "PRODUCT_RULE_Z1", "PRODUCT_RULE_Z2", "PRODUCT_RULE",
    "cell_kernel_masses", "target_radius_kernel_masses",
    "green_cell_kernel_masses", "kernel_masses", "correction_kernel_masses",
    "first_correction_cell_kernel_masses", "kchi_cell_kernel_masses",
    "first_correction_kernel_masses", "source_cell_kernel_masses",
    "ordinary_product", "build", "build_w0_correction",
    "build_from_global_low_source",
]


def _is_interval(value: Any) -> bool:
    return hasattr(value, "_mpi_") and hasattr(value, "a") and hasattr(value, "b")


def _raw_endpoint(value: Any, upper: bool) -> Any:
    raw = getattr(value, "_mpi_", None)
    if raw is not None:
        return raw[1 if upper else 0]
    return value


def _interval(c: Any, value: Any) -> Any:
    if _is_interval(value):
        return value
    if hasattr(value, "a") and hasattr(value, "b"):
        return c.make_mpf((_raw_endpoint(value.a, False), _raw_endpoint(value.b, True)))
    return c.mpf(value)


def _same_interval(left: Any, right: Any) -> bool:
    try:
        return left._mpi_ == right._mpi_
    except AttributeError:
        return left == right


def _endpoint_bounds(c: Any, value: Any) -> tuple[Any, Any]:
    value = _interval(c, value)
    lo, hi = value._mpi_
    return c.make_mpf((lo, lo)), c.make_mpf((hi, hi))


def _exact_cell(c: Any, index: int) -> Any:
    left = c.mpf(index) / 4
    right = c.mpf(index + 1) / 4
    return c.make_mpf((left._mpi_[0], right._mpi_[1]))


def _degenerate_endpoint(c: Any, value: Any, side: int) -> Any:
    value = _interval(c, value)
    raw = value._mpi_[side]
    return c.make_mpf((raw, raw))


def _is_zero(c: Any, value: Any) -> bool:
    lo, hi = _endpoint_bounds(c, value)
    return lo == 0 and hi == 0


def _positive(c: Any, value: Any, label: str) -> None:
    lo, hi = _endpoint_bounds(c, value)
    if hi <= 0:
        raise ValueError(f"{label} must be positive")
    if lo <= 0:
        raise ValueError(f"{label} interval contains the axis")


def _validate_geometry(c: Any, target_x: Any, source_left_x: Any,
                       source_right_x: Any, triangular: bool) -> tuple[Any, Any, Any]:
    """Validate ``0 <= source <= target`` without evaluating old fields."""

    x = _interval(c, target_x)
    a = _interval(c, source_left_x)
    b = x if triangular else _interval(c, source_right_x)
    xlo, xhi = _endpoint_bounds(c, x)
    alo, ahi = _endpoint_bounds(c, a)
    blo, bhi = _endpoint_bounds(c, b)
    if alo < 0 or ahi < alo or blo < alo or bhi < blo:
        raise ValueError("source endpoints must be ordered and nonnegative")
    if not triangular and bhi > xhi:
        raise ValueError("source cell extends beyond target radius")
    if triangular and not _same_interval(b, x):
        raise ValueError("triangular source endpoint must be the target interval")
    return x, a, b


def _k1_primitive(c: Any, x: Any, s: Any) -> Any:
    """Primitive in ``s`` of ``T K0``'s F1 kernel."""

    s = _interval(c, s)
    if _is_zero(c, s):
        return c.mpf(0)
    _positive(c, x, "target radius")
    _positive(c, s, "source endpoint")
    return (x * x * s * s / 32
            - s**4 * c.log(x / s) / 16
            - s**4 / 64
            - s**6 / (96 * x * x))


def _k1_dx_primitive(c: Any, x: Any, s: Any) -> Any:
    """Primitive in ``s`` of ``d_x(T K0)``."""

    s = _interval(c, s)
    if _is_zero(c, s):
        return c.mpf(0)
    _positive(c, x, "target radius")
    _positive(c, s, "source endpoint")
    return x * s * s / 16 - s**4 / (16 * x) + s**6 / (48 * x**3)


def _nonnegative(c: Any, value: Any, label: str) -> Any:
    """Intersect an analytically nonnegative mass with zero, outwardly."""

    value = _interval(c, value)
    lo, hi = value._mpi_
    hi_i = c.make_mpf((hi, hi))
    if hi_i < 0:
        raise ValueError(f"{label} has a negative upper bound")
    lo_i = c.make_mpf((lo, lo))
    if lo_i < 0:
        lo = c.mpf(0)._mpi_[0]
    return c.make_mpf((lo, hi))


def cell_kernel_masses(
    c: Any,
    target_x: Any,
    source_left_x: Any,
    source_right_x: Any = None,
    *,
    triangular: bool = False,
) -> dict[str, Any]:
    """Return exact directed bare ``T K0`` masses on one source cell.

    ``target_x`` and source endpoints are physical ``x`` radii.  The eight
    canonical cells are ``[i/4,(i+1)/4]`` for ``0 <= i < 8``; the diagonal
    source cell is triangular when ``triangular=True`` and its right endpoint
    is then the moving target interval.  The first triangular cell has the
    removable limits ``F1=x**4/192`` and ``d_x_F1=x**3/48``.
    """

    if source_right_x is None and not triangular:
        raise ValueError("source_right_x is required for a nontriangular cell")
    x, a, b = _validate_geometry(c, target_x, source_left_x,
                                  source_right_x if source_right_x is not None else target_x,
                                  triangular)
    xlo, _ = _endpoint_bounds(c, x)
    alo, _ = _endpoint_bounds(c, a)

    # This branch covers the complete first target cell [0, 1/4], not just a
    # point at x=0.  It is the analytic axis continuation of the triangular
    # integral, not a zero padding or a midpoint evaluation.
    _, xhi = _endpoint_bounds(c, x)
    if _is_zero(c, b):
        # A degenerate source at the axis has zero mass even when a caller
        # asks for a nontriangular zero-radius limit.
        f1 = c.mpf(0)
        dxf1 = c.mpf(0)
    elif triangular and alo == 0 and xlo == 0:
        f1 = x**4 / 192
        dxf1 = x**3 / 48
    elif xlo == 0 and xhi == 0:
        raise ValueError("a nonzero source cannot be integrated at x=0")
    else:
        _positive(c, x, "target radius")
        # On a diagonal cell H(x,x) and H_x(x,x) are polynomial values.  This
        # avoids evaluating log(x/x) on a nondegenerate interval, which would
        # introduce artificial dependency width.
        if triangular:
            f1_hi = x**4 / 192
            dxf1_hi = x**3 / 48
        else:
            f1_hi = _k1_primitive(c, x, b)
            dxf1_hi = _k1_dx_primitive(c, x, b)
        if _is_zero(c, a):
            f1_lo = c.mpf(0)
            dxf1_lo = c.mpf(0)
        else:
            f1_lo = _k1_primitive(c, x, a)
            dxf1_lo = _k1_dx_primitive(c, x, a)
        f1 = f1_hi - f1_lo
        dxf1 = dxf1_hi - dxf1_lo

    f1 = _nonnegative(c, f1, "T K0 F1 mass")
    dxf1 = _nonnegative(c, dxf1, "T K0 d_x_F1 mass")
    return {
        "F1": f1,
        "d_x_F1": dxf1,
        "T_K0_F1": f1,
        "T_K0_d_x_F1": dxf1,
        "kernel_order": 1,
        "principal_kernel": "T K0",
        "triangular": bool(triangular),
        "axis_continuation": bool(triangular and alo == 0 and xlo == 0),
    }


def target_radius_kernel_masses(
    c: Any,
    target_x: Any,
    source_left_x: Any,
    source_right_x: Any = None,
    *,
    triangular: bool = False,
) -> dict[str, Any]:
    """Alias with explicit physical-radius naming."""

    return cell_kernel_masses(c, target_x, source_left_x, source_right_x,
                              triangular=triangular)


green_cell_kernel_masses = cell_kernel_masses
kernel_masses = cell_kernel_masses


def correction_kernel_masses(
    c: Any,
    target_x: Any,
    source_left_x: Any,
    source_right_x: Any = None,
    *,
    chi: Any,
    triangular: bool = False,
) -> dict[str, Any]:
    """Return the signed first correction masses ``(-2*chi) T K0``."""

    bare = cell_kernel_masses(c, target_x, source_left_x, source_right_x,
                              triangular=triangular)
    q = -2 * _interval(c, chi)
    return {
        "F1": q * bare["F1"],
        "d_x_F1": q * bare["d_x_F1"],
        "bare_F1": bare["F1"],
        "bare_d_x_F1": bare["d_x_F1"],
        "correction_factor": q,
        "kernel_order": 1,
        "principal_kernel": "(-2*chi) T K0",
        "triangular": bare["triangular"],
        "axis_continuation": bare["axis_continuation"],
    }


first_correction_cell_kernel_masses = correction_kernel_masses
kchi_cell_kernel_masses = correction_kernel_masses
first_correction_kernel_masses = correction_kernel_masses
source_cell_kernel_masses = cell_kernel_masses


def ordinary_product(left: Sequence[Any], right: Sequence[Any]) -> list[Any]:
    """Multiply ordinary rows through ``Z2`` (derivatives, not Taylor data)."""

    if (not isinstance(left, Sequence) or isinstance(left, (str, bytes))
            or len(left) != 3):
        raise ValueError("left ordinary row must contain Z0, Z1, Z2")
    if (not isinstance(right, Sequence) or isinstance(right, (str, bytes))
            or len(right) != 3):
        raise ValueError("right ordinary row must contain Z0, Z1, Z2")
    a0, a1, a2 = left
    b0, b1, b2 = right
    return [a0 * b0, a1 * b0 + a0 * b1,
            a2 * b0 + 2 * a1 * b1 + a0 * b2]


def _as_items(value: Any) -> list[Any]:
    if value is None:
        return []
    if isinstance(value, (str, bytes, bytearray)):
        return [value]
    if isinstance(value, Sequence):
        return list(value)
    return [value]


def _freeze(value: Any) -> Any:
    try:
        hash(value)
    except TypeError:
        return (type(value).__name__, repr(value))
    return (type(value).__name__, value)


def _unique(values: Iterable[Any]) -> list[Any]:
    out: list[Any] = []
    seen: set[Any] = set()
    for value in values:
        key = _freeze(value)
        if key not in seen:
            seen.add(key)
            out.append(value)
    return out


def _origin_info(value: Any) -> dict[str, Any]:
    """Collect shallow source IDs/classes while retaining caller metadata."""

    if not isinstance(value, Mapping):
        return {"source_ids": [], "error_classes": [], "uncertainty": False}
    ids = list(_as_items(value.get("source_ids", value.get("source_id"))))
    classes = list(_as_items(value.get("error_classes", value.get("error_class"))))
    uncertain = bool(value.get("uncertainty", value.get("is_uncertainty", False)))
    for child_key in ("origin", "metadata"):
        child = value.get(child_key)
        if isinstance(child, Mapping):
            nested = _origin_info(child)
            ids.extend(nested["source_ids"])
            classes.extend(nested["error_classes"])
            uncertain = uncertain or nested["uncertainty"]
    return {"source_ids": _unique(ids), "error_classes": _unique(classes),
            "uncertainty": uncertain}


def _new_origin(operation: str, *inputs: Any, **metadata: Any) -> dict[str, Any]:
    ids: list[Any] = []
    classes: list[Any] = []
    uncertain = bool(metadata.get("uncertainty", False))
    for item in inputs:
        info = _origin_info(item)
        ids.extend(info["source_ids"])
        classes.extend(info["error_classes"])
        uncertain = uncertain or info["uncertainty"]
    ids.extend(_as_items(metadata.get("source_ids")))
    classes.extend(_as_items(metadata.get("error_classes")))
    out = {
        "operation": operation,
        "source_ids": _unique(ids),
        "error_classes": _unique(classes),
        "uncertainty": bool(uncertain),
        "symbolic_cancellation_claimed": False,
    }
    if inputs:
        out["parents"] = tuple({"source_ids": _origin_info(item)["source_ids"]}
                                for item in inputs)
    for key, value in metadata.items():
        if key not in {"source_ids", "error_classes", "uncertainty"}:
            out[key] = value
    return out


def _powers(value: Mapping[str, Any], label: str) -> tuple[int, int, int, int]:
    if "powers" in value:
        raw = value["powers"]
    else:
        raw = (value.get("hb_power", 0), value.get("F0_power", 0),
               value.get("Pstar_squared_power", value.get("Pstar_power", 0)),
               value.get("Lambda_power", 0))
    if (not isinstance(raw, Sequence) or isinstance(raw, (str, bytes, bytearray))
            or len(raw) != 4
            or any(not isinstance(item, int) or isinstance(item, bool) for item in raw)):
        raise ValueError(f"{label} must contain four integer formal powers")
    return tuple(int(item) for item in raw)


def _term(c: Any, raw: Any, *, label: str, default_source_id: Any) -> dict[str, Any]:
    if isinstance(raw, Mapping):
        coefficient = raw.get("coefficient")
        powers = _powers(raw, label)
        source_ids = list(_as_items(raw.get("source_ids", raw.get("source_id"))))
        origin = dict(raw.get("origin", {})) if isinstance(raw.get("origin"), Mapping) else {}
        info = _origin_info(raw)
        source_ids.extend(info["source_ids"])
        source_ids = _unique(source_ids)
        if not source_ids:
            source_ids = [raw.get("source_cell", default_source_id)]
        classes = _unique(info["error_classes"])
        uncertain = bool(info["uncertainty"])
        extras = dict(raw)
    else:
        coefficient = raw
        powers = (0, 0, 0, 0)
        source_ids = [default_source_id]
        classes = []
        uncertain = False
        origin = {}
        extras = {}
    if coefficient is None:
        raise ValueError(f"{label} requires coefficient")
    coefficient = _interval(c, coefficient)
    origin.update(source_ids=source_ids, error_classes=classes, uncertainty=uncertain)
    result = dict(extras)
    result.update({
        "powers": powers,
        "hb_power": powers[0],
        "F0_power": powers[1],
        "Pstar_squared_power": powers[2],
        "Lambda_power": powers[3],
        "coefficient": coefficient,
        "origin": origin,
        "source_ids": source_ids,
        "error_classes": classes,
        "uncertainty": uncertain,
    })
    return result


def _term_list(c: Any, value: Any, *, label: str, default_source_id: Any,
               allow_empty: bool = True) -> list[dict[str, Any]]:
    if isinstance(value, Mapping) and ("coefficient" in value or "powers" in value):
        values = [value]
    elif isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        values = list(value)
    else:
        values = [value]
    if not values and not allow_empty:
        raise ValueError(f"{label} must not be empty")
    return [_term(c, item, label=f"{label}[{index}]",
                  default_source_id=default_source_id) for index, item in enumerate(values)]


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
    raise ValueError(f"integer Z order required in {label}")


def _rows(c: Any, value: Any, *, label: str, default_source_id: Any) -> dict[int, list[dict[str, Any]]]:
    """Normalize a complete source/chi ordinary row through Z2."""

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
        parsed = {index: item for index, item in enumerate(value)}
    else:
        raise TypeError(f"{label} must provide Z0/Z1/Z2 rows")
    if set(parsed) != set(ORDERS):
        raise ValueError(f"{label} must provide exactly Z0, Z1, and Z2")
    return {order: _term_list(c, parsed[order], label=f"{label}[Z{order}]",
                              default_source_id=f"{default_source_id}:Z{order}")
            for order in ORDERS}


def _chi_rows(c: Any, value: Any) -> dict[int, list[dict[str, Any]]]:
    if isinstance(value, Mapping):
        for wrapper in ("ordinary_Z", "chi_jets", "jets", "rows"):
            if wrapper in value:
                value = value[wrapper]
                break
    if not isinstance(value, Mapping) and not isinstance(value, Sequence):
        # A scalar chi is a legitimate Z-independent input; its zero first
        # and second derivatives are explicit, rather than inferred source
        # rows.  This convenience does not affect source rows, which remain
        # mandatory through Z2.
        z = _term(c, value, label="chi[Z0]", default_source_id="chi")
        zero = c.mpf(0)
        return {0: [z], 1: [_term(c, zero, label="chi[Z1]", default_source_id="chi_Z1")],
                2: [_term(c, zero, label="chi[Z2]", default_source_id="chi_Z2")]}
    return _rows(c, value, label="chi", default_source_id="chi")


def _lookup(mapping: Mapping[Any, Any], key: Any) -> Any:
    if key in mapping:
        return mapping[key]
    other = str(key) if not isinstance(key, str) else (int(key) if key.isdigit() else key)
    if other in mapping:
        return mapping[other]
    raise KeyError(key)


def _cell_radius(c: Any, cell: Mapping[str, Any] | None, index: int) -> Any:
    expected = _exact_cell(c, index)
    if cell is None:
        return expected
    candidate = None
    for key in ("x", "target_x", "radius", "t"):
        if key in cell:
            candidate = cell[key]
            break
    if candidate is None:
        return expected
    actual = _interval(c, candidate)
    if not _same_interval(actual, expected):
        raise ValueError(f"cell {index} must be the canonical x=[i/4,(i+1)/4] interval")
    return actual


def _source_rows_from_cell(c: Any, cell: Mapping[str, Any], index: int) -> dict[int, list[dict[str, Any]]]:
    for key in ("source_jets", "source_rows", "source_terms", "W0_ordinary_Z", "W0", "source"):
        if key in cell:
            candidate = cell[key]
            # A six-column W0 packet is accepted only when its F1/row-5
            # column is explicitly selected; no other old field is rebuilt.
            if key in {"W0_ordinary_Z", "W0"} and isinstance(candidate, Mapping):
                try:
                    candidate = {order: _lookup(_lookup(candidate, order), 5)
                                 for order in ORDERS}
                except (KeyError, TypeError):
                    try:
                        candidate = {order: _lookup(_lookup(candidate, order), "F1")
                                     for order in ORDERS}
                    except (KeyError, TypeError):
                        pass
            return _rows(c, candidate, label=f"source_cell[{index}]",
                         default_source_id=f"source_cell:{index}")
    raise ValueError(f"source cell {index} is missing complete source Z0/Z1/Z2 rows")


def _global_source_orientation(source_rows: Any, *, row_key: Any) -> Any:
    """Accept the original ``global_low_source[Z][row][cell]`` orientation.

    The leaf still consumes exactly the caller's source terms.  This adapter
    only changes a bookkeeping axis: ordinary order first becomes source cell
    first.  It is intentionally structural and never derives a missing jet.
    """

    order_values: dict[int, Any] = {}
    if isinstance(source_rows, Mapping):
        for order in ORDERS:
            try:
                value = _lookup(source_rows, order)
            except KeyError:
                return source_rows
            if isinstance(value, Mapping):
                try:
                    value = _lookup(value, row_key)
                except KeyError:
                    try:
                        value = _lookup(value, str(row_key))
                    except KeyError:
                        return source_rows
            order_values[order] = value
    elif isinstance(source_rows, Sequence) and not isinstance(source_rows, (str, bytes, bytearray)):
        if len(source_rows) != 3:
            return source_rows
        order_values = {order: source_rows[order] for order in ORDERS}
    else:
        return source_rows
    if any(not isinstance(order_values[order], Sequence)
           or isinstance(order_values[order], (str, bytes, bytearray))
           or len(order_values[order]) != CELL_COUNT for order in ORDERS):
        return source_rows
    return [{order: order_values[order][index] for order in ORDERS}
            for index in range(CELL_COUNT)]


def _normalize_source_cells(c: Any, source_cells: Any, source_rows: Any,
                            *, source_row_key: Any = 5) -> tuple[list[Any], list[dict[int, list[dict[str, Any]]]]]:
    if source_cells is None:
        cells = [None] * CELL_COUNT
    elif isinstance(source_cells, Mapping):
        cells = [_lookup(source_cells, index) for index in range(CELL_COUNT)]
    elif isinstance(source_cells, Sequence) and not isinstance(source_cells, (str, bytes, bytearray)):
        if len(source_cells) != CELL_COUNT:
            raise ValueError("source_cells must contain eight canonical cells")
        cells = list(source_cells)
    else:
        raise TypeError("source_cells must be an eight-cell sequence or mapping")
    radii = [_cell_radius(c, item if isinstance(item, Mapping) else None, index)
             for index, item in enumerate(cells)]
    if source_rows is not None:
        source_rows = _global_source_orientation(source_rows, row_key=source_row_key)
        if isinstance(source_rows, Mapping):
            raw_rows = [_lookup(source_rows, index) for index in range(CELL_COUNT)]
        elif isinstance(source_rows, Sequence) and not isinstance(source_rows, (str, bytes, bytearray)):
            if len(source_rows) != CELL_COUNT:
                raise ValueError("source_rows must contain eight source-cell rows")
            raw_rows = list(source_rows)
        else:
            raise TypeError("source_rows must be an eight-cell sequence or mapping")
        rows = [_rows(c, raw, label=f"source_rows[{index}]",
                      default_source_id=f"source_cell:{index}")
                for index, raw in enumerate(raw_rows)]
    else:
        if any(not isinstance(item, Mapping) for item in cells):
            raise ValueError("source_rows or source-cell mappings are required")
        rows = [_source_rows_from_cell(c, item, index)
                for index, item in enumerate(cells)]
    return radii, rows


def _correction_terms(
    c: Any,
    rows: Mapping[int, list[Mapping[str, Any]]],
    chi_rows: Mapping[int, list[Mapping[str, Any]]],
    masses: Mapping[str, Any],
    *,
    target_index: int,
    source_index: int,
    target_x: Any,
    source_x: Any,
) -> dict[str, dict[int, list[dict[str, Any]]]]:
    output = {field: {order: [] for order in ORDERS} for field in FIELDS}
    rules = {0: PRODUCT_RULE_Z0, 1: PRODUCT_RULE_Z1, 2: PRODUCT_RULE_Z2}
    for order in ORDERS:
        for chi_order, source_order, leibniz in rules[order]:
            for chi_index, chi_term in enumerate(chi_rows[chi_order]):
                for source_term_index, source_term in enumerate(rows[source_order]):
                    cp = tuple(a + b for a, b in zip(chi_term["powers"], source_term["powers"]))
                    origin = _new_origin(
                        "kchi_first_correction_product", chi_term["origin"], source_term["origin"],
                        source_ids=[f"target_cell:{target_index}", f"source_cell:{source_index}"],
                        target_cell=target_index, source_cell=source_index,
                        chi_order=chi_order, source_order=source_order,
                        leibniz_multiplier=leibniz,
                        kernel_order=1, kernel="(-2*chi) T K0",
                        source_cell_id=source_term.get("source_cell", source_term.get("source_id", source_index)),
                        chi_term_index=chi_index, source_term_index=source_term_index,
                    )
                    for field in FIELDS:
                        mass = masses[field]
                        output[field][order].append({
                            "powers": cp,
                            "hb_power": cp[0], "F0_power": cp[1],
                            "Pstar_squared_power": cp[2], "Lambda_power": cp[3],
                            "coefficient": -2 * leibniz * chi_term["coefficient"] * source_term["coefficient"] * mass,
                            "origin": origin,
                            "source_ids": list(origin["source_ids"]),
                            "error_classes": list(origin["error_classes"]),
                            "uncertainty": bool(origin["uncertainty"]),
                            "target_cell": target_index,
                            "source_cell": source_index,
                            "source_cell_id": source_term.get("source_cell", source_term.get("source_id", source_index)),
                            "source_term_id": source_term.get("source_id", source_term.get("source_cell", source_index)),
                            "source_index": source_index,
                            "source_term_index": source_term_index,
                            "chi_order": chi_order,
                            "source_order": source_order,
                            "leibniz_multiplier": leibniz,
                            "kernel_mass": mass,
                            "kernel_mass_field": field,
                            "kernel_mass_bare": mass,
                            "correction_factor": -2 * chi_term["coefficient"],
                            "target_x": target_x,
                            "source_x": source_x,
                            "principal_kernel": "(-2*chi) T K0",
                        })
    return output


def _aggregate_by_order(values: Mapping[int, Mapping[str, Sequence[Mapping[str, Any]]]]) -> dict[str, dict[int, list[dict[str, Any]]]]:
    return {field: {order: list(values[field][order]) for order in ORDERS} for field in FIELDS}


def build(
    c: Any,
    source_cells: Any = None,
    chi: Any = None,
    *,
    source_rows: Any = None,
    chi_jets: Any = None,
    cells: Any = None,
    source_family: Any = None,
    source_atoms: Any = None,
    certificate: Any = None,
    source_row_key: Any = 5,
    global_low_source: Any = None,
) -> dict[str, Any]:
    """Assemble the first ``(-2chi) T K0`` correction by target/source cell.

    ``source_cells`` may contain a ``source_jets``/``source_rows`` mapping
    through Z2.  The eight canonical cells are ``[i/4,(i+1)/4]`` and cover
    ``0 <= x <= 2``.  Alternatively pass ``source_rows`` explicitly.  Each
    source row is a list of already-enclosed formal terms; no old source term
    is recomputed.  ``chi`` or ``chi_jets`` supplies the three ordinary chi rows.
    A scalar chi is treated as explicitly Z-independent; source rows are
    never padded.
    """

    if not hasattr(c, "mpf") or not hasattr(c, "make_mpf"):
        raise TypeError("c must be the caller's MPIntervalContext")
    if source_cells is None and cells is not None:
        source_cells = cells
    if source_rows is None and global_low_source is not None:
        source_rows = global_low_source
    if chi is None:
        chi = chi_jets
    if chi is None:
        raise ValueError("chi or chi_jets is required")
    radii, rows = _normalize_source_cells(c, source_cells, source_rows,
                                           source_row_key=source_row_key)
    chi_rows = _chi_rows(c, chi)
    target_cells: list[dict[str, Any]] = []
    for target_index, target_x in enumerate(radii):
        by_source: dict[str, Any] = {}
        aggregate = {field: {order: [] for order in ORDERS} for field in FIELDS}
        mass_by_source: dict[str, Any] = {}
        for source_index in range(target_index + 1):
            source_x = radii[source_index]
            source_left = _degenerate_endpoint(c, source_x, 0)
            source_right = target_x if source_index == target_index else _degenerate_endpoint(c, source_x, 1)
            masses = cell_kernel_masses(c, target_x, source_left, source_right,
                                        triangular=source_index == target_index)
            correction = _correction_terms(
                c, rows[source_index], chi_rows, masses,
                target_index=target_index, source_index=source_index,
                target_x=target_x, source_x=source_x,
            )
            order_view = {field: {f"Z{order}": correction[field][order]
                                  for order in ORDERS} for field in FIELDS}
            by_source[str(source_index)] = {
                "source_index": source_index,
                "target_index": target_index,
                "target_x": target_x,
                "source_x": source_x,
                "triangular": source_index == target_index,
                "kernel_masses": masses,
                "correction": order_view,
                "W0_correction": order_view,
                "F1": order_view["F1"],
                "d_x_F1": order_view["d_x_F1"],
            }
            mass_by_source[str(source_index)] = masses
            for field in FIELDS:
                for order in ORDERS:
                    aggregate[field][order].extend(correction[field][order])
        aggregate_view = {field: {f"Z{order}": aggregate[field][order]
                                  for order in ORDERS} for field in FIELDS}
        target_cells.append({
            "index": target_index,
            "target_index": target_index,
            "x": target_x,
            "target_x": target_x,
            "source_indices": list(range(target_index + 1)),
            "source_j_le_target_i": True,
            "triangular_diagonal_source": target_index,
            "axis_endpoint_continuation": target_index == 0,
            "W0_correction": aggregate_view,
            "first_correction": aggregate_view,
            "W0_by_order": aggregate_view,
            "F1": aggregate_view["F1"],
            "d_x_F1": aggregate_view["d_x_F1"],
            "W0_Z0": {field: aggregate[field][0] for field in FIELDS},
            "W0_Z1": {field: aggregate[field][1] for field in FIELDS},
            "W0_Z2": {field: aggregate[field][2] for field in FIELDS},
            "W0_correction_by_source_cell": by_source,
            "W0_by_order_by_source_cell": by_source,
            "kernel_masses_by_source_cell": mass_by_source,
        })

    atom_names = list(SOURCE_POWER_BASIS) if source_atoms is None else source_atoms
    return {
        "source_family": source_family,
        "source_atoms": atom_names,
        "source_power_basis": list(SOURCE_POWER_BASIS),
        "cell_count": CELL_COUNT,
        "target_cells": target_cells,
        "outer_cells": target_cells,
        "source_cells": source_cells,
        "chi_jets": chi_rows,
        "source_rows_preserved_through_Z2": True,
        "chi_rows_preserved_through_Z2": True,
        "ordinary_Z_orders": list(ORDERS),
        "ordinary_product_rule": PRODUCT_RULE,
        "ordinary_Z0_product_rule": PRODUCT_RULE_Z0,
        "ordinary_Z1_product_rule": PRODUCT_RULE_Z1,
        "ordinary_Z2_product_rule": PRODUCT_RULE_Z2,
        "source_j_le_target_i": True,
        "triangular_diagonal_source": True,
        "axis_endpoint_continuation": True,
        "axis_log_division_avoided": True,
        "green_kernel_mass_definition": "directed source-cell integral of T K0",
        "green_kernel_masses_exact": True,
        "signed_first_correction": "(-2*chi) T K0",
        "principal_kernel": "K0 plus first correction only",
        "Kchi_first_correction_only": True,
        "Kchi_higher_words_included": False,
        "allword_tail_added": False,
        "source_to_W0_only": True,
        "formal_powers_retained": True,
        "source_indices_retained": True,
        "source_ids_retained": True,
        "signed_coefficients_retained": True,
        "source_terms_recomputed": False,
        "old_words_recomputed": False,
        "certificate_passthrough": certificate,
        "geometry_certificate_passthrough": certificate is not None,
        "source_row_key": source_row_key,
        "global_low_source_orientation_adapter": global_low_source is not None,
    }


build_w0_correction = build


def build_from_global_low_source(
    c: Any,
    global_low_source: Any,
    chi: Any,
    **kwargs: Any,
) -> dict[str, Any]:
    """Convenience entry using original ``global_low_source[Z][row][cell]``.

    ``row`` defaults to 5 (the swirl source row) and can be changed with the
    forwarded ``source_row_key`` keyword.  The source terms themselves are
    passed through unchanged.
    """

    return build(c, chi=chi, global_low_source=global_low_source, **kwargs)

