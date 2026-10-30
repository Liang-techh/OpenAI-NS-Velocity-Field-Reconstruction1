"""Cellwise outer first ``Kchi`` feedback for saved ``W0`` sectors.

This leaf is intentionally smaller than a ``W1`` builder.  It consumes two
already saved objects only:

* the zeroth/first ordinary ``C0`` and ``B1`` matrix sectors in each source
  cell, and
* a six-column ``W0_ordinary_Z`` packet through ordinary ``Z2``.

For every target cell ``i`` and source cell ``j <= i`` it forms

``(-2*chi) T K0 (C0 W0 + B1 d_Z W0)``

as unmerged signed terms.  The output contains the ``F1`` and ``d_x_F1``
source-cell masses of ``T K0``.  The diagonal source cell is triangular, and
the first triangular cell uses the removable axis continuations
``F1=x**4/192`` and ``d_x_F1=x**3/48``.  Matrix terms, W0 terms and chi terms
retain their four formal powers and source IDs; no coefficient is made
positive and no old producer is imported or replayed.

For reference, the physical kernels integrated over each directed source
cell are
``K_F1(x,s)=s*x**2/16-s**3*log(x/s)/4-s**5/(16*x**2)`` and
``d_x K_F1(x,s)=s*(x**2-s**2)**2/(8*x**3)``.  Their endpoint primitives are
used directly, with the diagonal upper endpoint replaced by the polynomial
values above.

``build`` can also receive an independent ``delta_w0_cells`` packet.  That
packet is evaluated by the same leaf and returned in a separate
``delta_W0`` branch.  It is useful for the cross correction
``Delta G1 Q Delta W0``: columns 1 (``F1``) and 5 (``d_x_F1``) may be signed
terms while columns 2, 3, 4 and 6 are explicit empty sequences.  Empty
sequences are preserved as empty; no values are inferred for them.

The canonical coordinate is physical ``x`` with eight cells
``[i/4,(i+1)/4]``.  Saved core packets commonly expose the equivalent
normalized ``t=x/2`` cells ``[i/8,(i+1)/8]``; both spellings are accepted.
The implementation is independent of the older Kchi and W1 builders.

For downstream remainder propagation, each branch exposes
``feedback_q5_products_by_source_cell`` (the raw C0/B1 q5 products before
chi) and ``unweighted_outer_delta_source_by_source_cell`` (the same q5
products after ``-2*chi``/``chi_Z`` but before the radial ``T K0`` mass).
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any


SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
CELL_COUNT = 8
ORDERS = (0, 1, 2)
OUTPUT_ORDERS = (0, 1)
MATRIX_ORDERS = ("0", "1")
FIELDS = ("F1", "d_x_F1")

EXPECTED_C0 = frozenset(
    ("4,1", "5,1", "5,2", "5,3", "5,5",
     "6,1", "6,2", "6,3", "6,4", "6,6")
)
EXPECTED_B1 = frozenset(("5,1", "5,2", "5,3", "6,2", "6,3", "6,4"))

# Ordinary derivatives are not Taylor coefficients.  These are the source
# products before multiplication by the outer chi jet.
FEEDBACK_PRODUCT_RULE_Z0 = (
    ("C0", "0", "0", 1, "C0 W0"),
    ("B1", "0", "1", 1, "B1 W0_Z"),
)
FEEDBACK_PRODUCT_RULE_Z1 = (
    ("C0", "1", "0", 1, "C0_Z W0"),
    ("C0", "0", "1", 1, "C0 W0_Z"),
    ("B1", "1", "1", 1, "B1_Z W0_Z"),
    ("B1", "0", "2", 1, "B1 W0_ZZ"),
)
FEEDBACK_PRODUCT_RULE = {"Z0": FEEDBACK_PRODUCT_RULE_Z0,
                         "Z1": FEEDBACK_PRODUCT_RULE_Z1}
# Compatibility spelling for callers that use the core-profile convention;
# these are still only the source-feedback products, not a W1 result.
PRODUCT_RULE_Z0 = FEEDBACK_PRODUCT_RULE_Z0
PRODUCT_RULE_Z1 = FEEDBACK_PRODUCT_RULE_Z1
PRODUCT_RULE = FEEDBACK_PRODUCT_RULE

# The outer kernel depends on Z through chi.  In particular, Z1 is
# chi_Z * source_Z0 + chi * source_Z1; the chi_Z terms must not be dropped.
OUTER_PRODUCT_RULE_Z0 = ((0, 0, 1, "chi * feedback_Z0"),)
OUTER_PRODUCT_RULE_Z1 = (
    (1, 0, 1, "chi_Z * feedback_Z0"),
    (0, 1, 1, "chi * feedback_Z1"),
)
OUTER_PRODUCT_RULE = {"Z0": OUTER_PRODUCT_RULE_Z0,
                      "Z1": OUTER_PRODUCT_RULE_Z1}


__all__ = [
    "SOURCE_POWER_BASIS", "CELL_COUNT", "ORDERS", "OUTPUT_ORDERS",
    "MATRIX_ORDERS", "FIELDS", "EXPECTED_C0", "EXPECTED_B1",
    "FEEDBACK_PRODUCT_RULE_Z0", "FEEDBACK_PRODUCT_RULE_Z1",
    "FEEDBACK_PRODUCT_RULE", "OUTER_PRODUCT_RULE_Z0",
    "OUTER_PRODUCT_RULE_Z1", "OUTER_PRODUCT_RULE", "cell_kernel_masses",
    "target_radius_kernel_masses", "green_cell_kernel_masses",
    "kernel_masses", "correction_kernel_masses", "ordinary_product",
    "PRODUCT_RULE_Z0", "PRODUCT_RULE_Z1", "PRODUCT_RULE", "build",
    "build_outer_feedback", "build_delta_w0",
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
        return c.make_mpf((_raw_endpoint(value.a, False),
                           _raw_endpoint(value.b, True)))
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


def _nonnegative(c: Any, value: Any, label: str) -> Any:
    """Intersect an analytically nonnegative kernel mass with zero."""

    value = _interval(c, value)
    raw_lo, raw_hi = value._mpi_
    hi = c.make_mpf((raw_hi, raw_hi))
    if hi < 0:
        raise ValueError(f"{label} has a negative upper bound")
    lo = c.make_mpf((raw_lo, raw_lo))
    if lo < 0:
        raw_lo = c.mpf(0)._mpi_[0]
    return c.make_mpf((raw_lo, raw_hi))


def _exact_x_cell(c: Any, index: int) -> Any:
    left = c.mpf(index) / 4
    right = c.mpf(index + 1) / 4
    return c.make_mpf((left._mpi_[0], right._mpi_[1]))


def _exact_t_cell(c: Any, index: int) -> Any:
    left = c.mpf(index) / 8
    right = c.mpf(index + 1) / 8
    return c.make_mpf((left._mpi_[0], right._mpi_[1]))


def _validate_geometry(c: Any, target_x: Any, source_left_x: Any,
                       source_right_x: Any, triangular: bool) -> tuple[Any, Any, Any]:
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
    """Primitive in ``s`` of the ``F1`` kernel of ``T K0``."""

    s = _interval(c, s)
    if _is_zero(c, s):
        return c.mpf(0)
    _positive(c, x, "target radius")
    _positive(c, s, "positive source endpoint")
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
    _positive(c, s, "positive source endpoint")
    return x * s * s / 16 - s**4 / (16 * x) + s**6 / (48 * x**3)


def cell_kernel_masses(
    c: Any,
    target_x: Any,
    source_left_x: Any,
    source_right_x: Any = None,
    *,
    triangular: bool = False,
) -> dict[str, Any]:
    """Return exact directed ``T K0`` masses on one physical source cell.

    ``target_x`` and source endpoints are physical radii.  For a triangular
    cell the right source endpoint is the moving target interval.  The
    polynomial upper primitive on that cell avoids evaluating ``log(x/x)``
    with dependent intervals.  The first triangular cell additionally avoids
    every division by an interval containing zero.
    """

    if source_right_x is None and not triangular:
        raise ValueError("source_right_x is required for a nontriangular cell")
    x, a, b = _validate_geometry(
        c, target_x, source_left_x,
        source_right_x if source_right_x is not None else target_x,
        triangular,
    )
    xlo, xhi = _endpoint_bounds(c, x)
    alo, _ = _endpoint_bounds(c, a)
    if _is_zero(c, b):
        f1 = c.mpf(0)
        dxf1 = c.mpf(0)
    elif triangular and alo == 0 and xlo == 0:
        # Exact removable limits for the complete first target cell.
        f1 = x**4 / 192
        dxf1 = x**3 / 48
    elif xlo == 0 and xhi == 0:
        raise ValueError("a nonzero source cannot be integrated at x=0")
    else:
        _positive(c, x, "target radius")
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


target_radius_kernel_masses = cell_kernel_masses
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
    """Return the signed bare-mass factor ``(-2*chi) T K0``."""

    bare = cell_kernel_masses(c, target_x, source_left_x, source_right_x,
                              triangular=triangular)
    factor = -2 * _interval(c, chi)
    return {
        "F1": factor * bare["F1"],
        "d_x_F1": factor * bare["d_x_F1"],
        "bare_F1": bare["F1"],
        "bare_d_x_F1": bare["d_x_F1"],
        "correction_factor": factor,
        "kernel_order": 1,
        "principal_kernel": "(-2*chi) T K0",
        "triangular": bare["triangular"],
        "axis_continuation": bare["axis_continuation"],
    }


def ordinary_product(left: Sequence[Any], right: Sequence[Any]) -> list[Any]:
    """Multiply ordinary rows through Z2 using derivative conventions."""

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


def _lookup(mapping: Mapping[Any, Any], key: Any) -> Any:
    if key in mapping:
        return mapping[key]
    if isinstance(key, str) and key.isdigit():
        alternates = (int(key), f"Z{key}")
    elif isinstance(key, int) and not isinstance(key, bool):
        alternates = (str(key), f"Z{key}")
    else:
        alternates = (str(key),)
    for alternate in alternates:
        if alternate in mapping:
            return mapping[alternate]
    raise KeyError(key)


def _as_items(value: Any) -> list[Any]:
    if value is None:
        return []
    if isinstance(value, (str, bytes, bytearray)):
        return [value]
    if isinstance(value, Mapping):
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


def _unique(values: Sequence[Any]) -> list[Any]:
    result: list[Any] = []
    seen: set[Any] = set()
    for value in values:
        key = _freeze(value)
        if key not in seen:
            seen.add(key)
            result.append(value)
    return result


def _origin_info(value: Any) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        return {"source_ids": [], "error_classes": [], "uncertainty": False}
    ids = list(_as_items(value.get("source_ids", value.get("source_id"))))
    classes = list(_as_items(value.get(
        "error_classes",
        value.get("error_class", value.get(
            "source_classes",
            value.get("source_class", value.get("uncertainty_class")))),
    )))
    uncertain = bool(value.get("uncertainty", value.get("is_uncertainty", False)))
    # source_cell is an identity when a producer does not also expose a
    # source_id.  It is intentionally retained, not replaced by a cell index.
    if not ids and value.get("source_cell") is not None:
        ids.append(value["source_cell"])
    for child_key in ("origin", "metadata", "parents"):
        child = value.get(child_key)
        if isinstance(child, Sequence) and not isinstance(child, (str, bytes, bytearray)):
            children = child
        else:
            children = (child,)
        for item in children:
            if isinstance(item, Mapping):
                nested = _origin_info(item)
                ids.extend(nested["source_ids"])
                classes.extend(nested["error_classes"])
                uncertain = uncertain or nested["uncertainty"]
    return {"source_ids": _unique(ids), "error_classes": _unique(classes),
            "uncertainty": uncertain}


def _power(value: Any, label: str) -> int:
    if not isinstance(value, int) or isinstance(value, bool):
        raise ValueError(f"{label} must be an integer formal power")
    return int(value)


def _powers(raw: Mapping[str, Any], label: str) -> tuple[int, int, int, int]:
    if "powers" in raw:
        values = raw["powers"]
        if (not isinstance(values, Sequence)
                or isinstance(values, (str, bytes, bytearray))
                or len(values) != 4):
            raise ValueError(f"{label} must contain four integer formal powers")
        return tuple(_power(value, label) for value in values)  # type: ignore[return-value]
    return (
        _power(raw.get("hb_power", 0), label),
        _power(raw.get("F0_power", raw.get("F0_anchor_power", 0)), label),
        _power(raw.get("Pstar_squared_power", raw.get("Pstar_power", 0)), label),
        _power(raw.get("Lambda_power", 0), label),
    )


def _term(c: Any, raw: Any, *, label: str, default_source_id: Any) -> dict[str, Any]:
    if isinstance(raw, Mapping):
        if "coefficient" not in raw:
            raise ValueError(f"{label} requires coefficient")
        coefficient = raw["coefficient"]
        powers = _powers(raw, label)
        info = _origin_info(raw)
        source_ids = list(info["source_ids"])
        if not source_ids:
            source_ids = [raw.get("source_cell", raw.get("source_id", default_source_id))]
        classes = list(info["error_classes"])
        uncertain = bool(info["uncertainty"])
        source_cell = raw.get("source_cell", default_source_id)
        source_id = raw.get("source_id", source_cell)
    else:
        coefficient = raw
        powers = (0, 0, 0, 0)
        source_ids = [default_source_id]
        classes = []
        uncertain = False
        source_cell = default_source_id
        source_id = default_source_id
    if coefficient is None:
        raise ValueError(f"{label} requires coefficient")
    return {
        "powers": powers,
        "hb_power": powers[0],
        "F0_power": powers[1],
        "Pstar_squared_power": powers[2],
        "Lambda_power": powers[3],
        "coefficient": _interval(c, coefficient),
        "source_ids": _unique(source_ids),
        "error_classes": _unique(classes),
        "source_classes": _unique(classes),
        "uncertainty_class": _unique(classes),
        "uncertainty": uncertain,
        "source_cell": source_cell,
        "source_id": source_id,
        "input_origin": raw.get("origin") if isinstance(raw, Mapping) else None,
    }


def _term_list(c: Any, value: Any, *, label: str,
               default_source_id: Any, allow_empty: bool = True) -> list[dict[str, Any]]:
    if value is None:
        values: list[Any] = []
    elif isinstance(value, Mapping) and ("coefficient" in value or "powers" in value):
        values = [value]
    elif isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        values = list(value)
    else:
        values = [value]
    if not values and not allow_empty:
        raise ValueError(f"{label} must not be empty")
    return [_term(c, item, label=f"{label}[{index}]",
                  default_source_id=default_source_id)
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
    raise ValueError(f"integer Z order required in {label}")


def _normalise_cells(value: Any, label: str) -> list[Mapping[str, Any]]:
    if isinstance(value, Mapping):
        cells = [_lookup(value, index) for index in range(CELL_COUNT)]
    elif isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        if len(value) != CELL_COUNT:
            raise ValueError(f"{label} must contain eight canonical cells")
        cells = list(value)
    else:
        raise TypeError(f"{label} must be an eight-cell sequence or mapping")
    if any(not isinstance(cell, Mapping) for cell in cells):
        raise TypeError(f"{label} entries must be mappings")
    return list(cells)


def _cell_x(c: Any, cell: Mapping[str, Any] | None, index: int, label: str) -> Any:
    expected_x = _exact_x_cell(c, index)
    if cell is None:
        return expected_x
    actual_x: Any = None
    if "x" in cell:
        actual_x = _interval(c, cell["x"])
        if not _same_interval(actual_x, expected_x):
            raise ValueError(f"{label} {index} has a noncanonical x interval")
    elif "target_x" in cell:
        actual_x = _interval(c, cell["target_x"])
        if not _same_interval(actual_x, expected_x):
            raise ValueError(f"{label} {index} has a noncanonical target_x interval")
    elif "radius" in cell:
        actual_x = _interval(c, cell["radius"])
        if not _same_interval(actual_x, expected_x):
            raise ValueError(f"{label} {index} has a noncanonical radius interval")

    if "t" in cell:
        actual_t = _interval(c, cell["t"])
        if not _same_interval(actual_t, _exact_t_cell(c, index)):
            raise ValueError(f"{label} {index} has a noncanonical t interval")
        t_x = 2 * actual_t
        if actual_x is not None and not _same_interval(actual_x, t_x):
            raise ValueError(f"{label} {index} has inconsistent x=2*t")
        actual_x = t_x
        if "rho" in cell and not _same_interval(_interval(c, cell["rho"]), 4 * actual_t * actual_t):
            raise ValueError(f"{label} {index} has inconsistent rho=4*t^2")
    return expected_x if actual_x is None else actual_x


def _cell_source_packet(cell: Mapping[str, Any], branch: str) -> Any:
    if branch == "K0_W0":
        keys = ("W0_ordinary_Z", "W0_K0_ordinary_Z", "W0_K0",
                "W0", "saved_W0", "saved_W0_K0", "w0")
    else:
        keys = ("delta_W0_ordinary_Z", "DeltaW0_ordinary_Z",
                "delta_W0", "DeltaW0", "delta_w0",
                "correction_W0", "W0_delta", "Delta_W0",
                "W0_ordinary_Z", "W0_K0_ordinary_Z", "W0")
    for key in keys:
        if key in cell:
            return cell[key]
    # A packet with no wrapper is useful for a compact delta input.  Do not
    # mistake the coordinate keys of a saved cell for ordinary orders.
    if any(str(key).upper().removeprefix("Z") in {"0", "1", "2"}
           for key in cell):
        return cell
    raise ValueError(f"saved {branch} cell is missing a six-column W0 packet")


def _packet_rows(c: Any, packet: Any, *, label: str, branch: str) -> dict[int, dict[int, list[dict[str, Any]]]]:
    if isinstance(packet, Mapping):
        for wrapper in ("ordinary_Z", "source_jets", "source_rows", "jets", "rows"):
            if wrapper in packet:
                packet = packet[wrapper]
                break
    if not isinstance(packet, Mapping):
        raise TypeError(f"{label} must be a Z0/Z1/Z2 six-column mapping")
    rows: dict[int, dict[int, list[dict[str, Any]]]] = {}
    for order in ORDERS:
        try:
            raw_row = _lookup(packet, order)
        except KeyError as exc:
            raise ValueError(f"{label} is missing W0 order Z{order}") from exc
        if not isinstance(raw_row, Mapping):
            raise TypeError(f"{label} order Z{order} must be a column mapping")
        row: dict[int, list[dict[str, Any]]] = {}
        for column in range(1, 7):
            try:
                value = _lookup(raw_row, str(column))
            except KeyError:
                # Named field spelling is convenient for a delta packet.  It
                # maps only the two meaningful correction columns.
                if column == 1 and "F1" in raw_row:
                    value = raw_row["F1"]
                elif column == 5 and "d_x_F1" in raw_row:
                    value = raw_row["d_x_F1"]
                elif branch == "delta_W0":
                    value = []
                else:
                    raise ValueError(f"{label} order Z{order} is missing column {column}") from exc
            row[column] = _term_list(
                c, value, label=f"{label}[Z{order}][{column}]",
                default_source_id=f"{branch}:Z{order}:column:{column}",
                allow_empty=True,
            )
        rows[order] = row
    return rows


def _matrix_terms(c: Any, cell: Mapping[str, Any], name: str,
                  order: str, key: str, *, index: int) -> list[dict[str, Any]]:
    try:
        sector = cell[name]
        row = _lookup(sector, order)
        value = _lookup(row, key)
    except (KeyError, TypeError) as exc:
        raise ValueError(f"matrix cell {index} is missing {name} {order} {key}") from exc
    return _term_list(c, value, label=f"matrix[{index}].{name}[{order}][{key}]",
                      default_source_id=f"matrix:{index}:{name}:{order}:{key}",
                      allow_empty=True)


def _validate_matrix(c: Any, cell: Mapping[str, Any], name: str, *, index: int) -> None:
    expected = EXPECTED_C0 if name == "C0" else EXPECTED_B1
    try:
        sector = cell[name]
    except (KeyError, TypeError) as exc:
        raise ValueError(f"matrix cell {index} is missing {name}") from exc
    if not isinstance(sector, Mapping):
        raise TypeError(f"matrix cell {index} {name} must be a mapping")
    for order in MATRIX_ORDERS:
        try:
            row = _lookup(sector, order)
        except KeyError as exc:
            raise ValueError(f"matrix cell {index} {name} is missing order {order}") from exc
        if not isinstance(row, Mapping):
            raise TypeError(f"matrix cell {index} {name} order {order} must be a mapping")
        if {str(key) for key in row} != set(expected):
            raise ValueError(f"matrix cell {index} {name} support changed at order {order}")
        for key in expected:
            _matrix_terms(c, cell, name, order, key, index=index)


def _chi_rows(c: Any, value: Any) -> dict[int, list[dict[str, Any]]]:
    if value is None:
        raise ValueError("chi or chi_jets is required")
    if isinstance(value, Mapping):
        if "coefficient" in value or "powers" in value:
            value = {0: value}
        else:
            for wrapper in ("ordinary_Z", "chi_jets", "jets", "rows"):
                if wrapper in value:
                    value = value[wrapper]
                    break
    if isinstance(value, Mapping):
        parsed: dict[int, Any] = {}
        for key, item in value.items():
            try:
                order = _order_key(key, "chi")
            except ValueError:
                continue
            parsed[order] = item
        if not parsed or 0 not in parsed:
            raise ValueError("chi must provide a Z0 row")
    elif isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        # A two/three-row sequence is the ordinary jet spelling.  A sequence
        # of scalar values is still accepted as one Z0 term list.  A list of
        # coefficient-bearing term mappings is also one Z0 row, even when it
        # happens to contain two or three terms.
        term_sequence = all(
            isinstance(item, Mapping)
            and ("coefficient" in item or "powers" in item)
            for item in value
        )
        if not term_sequence and len(value) in (2, 3):
            parsed = {index: item for index, item in enumerate(value)}
        else:
            parsed = {0: value}
    else:
        parsed = {0: value}
    result: dict[int, list[dict[str, Any]]] = {}
    for order in (0, 1, 2):
        if order in parsed:
            result[order] = _term_list(c, parsed[order], label=f"chi[Z{order}]",
                                       default_source_id=f"chi:Z{order}",
                                       allow_empty=True)
        elif order == 1:
            result[order] = [_term(c, {"coefficient": c.mpf(0),
                                       "powers": (0, 0, 0, 0),
                                       "source_id": "chi_Z"},
                                  label="chi[Z1]", default_source_id="chi_Z")]
    return result


def _source_ids(*terms: Any, extra: Sequence[Any] = ()) -> list[Any]:
    ids: list[Any] = []
    for term in terms:
        if isinstance(term, Mapping):
            ids.extend(_as_items(term.get("source_ids", term.get("source_id"))))
            if not ids or term.get("source_cell") is not None:
                ids.extend(_as_items(term.get("source_cell")))
            info = _origin_info(term)
            ids.extend(info["source_ids"])
    ids.extend(extra)
    return _unique(ids)


def _feedback_product(
    matrix_term: Mapping[str, Any], w0_term: Mapping[str, Any], multiplier: int,
    origin: str, *, source_index: int, matrix_name: str, matrix_order: str,
    matrix_key: str, matrix_term_index: int, w0_order: int,
    w0_term_index: int, branch: str,
) -> dict[str, Any]:
    row, column = (int(value) for value in matrix_key.split(","))
    powers = tuple(a + b for a, b in zip(matrix_term["powers"], w0_term["powers"]))
    ids = _source_ids(matrix_term, w0_term,
                      extra=(f"source_cell:{source_index}",))
    classes = _unique(list(matrix_term.get("error_classes", ()))
                      + list(w0_term.get("error_classes", ())))
    return {
        "powers": powers,
        "hb_power": powers[0],
        "F0_power": powers[1],
        "Pstar_squared_power": powers[2],
        "Lambda_power": powers[3],
        "coefficient": multiplier * matrix_term["coefficient"] * w0_term["coefficient"],
        "source_ids": ids,
        "error_classes": classes,
        "source_classes": classes,
        "uncertainty_class": classes,
        "uncertainty": bool(matrix_term.get("uncertainty", False)
                            or w0_term.get("uncertainty", False)),
        "origin": origin,
        "origin_record": {
            "operation": "outer_feedback_source_product",
            "source_ids": ids,
            "error_classes": classes,
            "source_classes": classes,
            "uncertainty_class": classes,
            "symbolic_cancellation_claimed": False,
        },
        "row": row,
        "column": column,
        "matrix_row": row,
        "matrix_column": column,
        "matrix_name": matrix_name,
        "matrix_key": matrix_key,
        "matrix_order": matrix_order,
        "matrix_term_index": matrix_term_index,
        "matrix_source_cell": matrix_term.get("source_cell", source_index),
        "matrix_source_cell_index": source_index,
        "matrix_source_id": matrix_term.get("source_id", matrix_term.get("source_cell", source_index)),
        "W0_order": w0_order,
        "W0_term_index": w0_term_index,
        "W0_column": column,
        "W0_column_index": column,
        "W0_source_cell": w0_term.get("source_cell", source_index),
        "W0_source_cell_index": w0_term.get("source_cell_index", source_index),
        "W0_source_id": w0_term.get("source_id", w0_term.get("source_cell", source_index)),
        "W0_provenance": "saved_outer_cells" if branch == "K0_W0" else "saved_delta_W0_cells",
        "source_cell": source_index,
        "source_cell_index": source_index,
        "source_cell_id": w0_term.get("source_cell", w0_term.get("source_id", source_index)),
        "branch": branch,
        "leibniz_multiplier": multiplier,
    }


def _feedback_terms(
    c: Any, matrix_cell: Mapping[str, Any], packet: Mapping[int, Mapping[int, Sequence[Mapping[str, Any]]]],
    *, source_index: int, branch: str, order: int,
) -> dict[str, list[dict[str, Any]]]:
    output = {"4": [], "5": [], "6": []}
    rules = FEEDBACK_PRODUCT_RULE_Z0 if order == 0 else FEEDBACK_PRODUCT_RULE_Z1
    for matrix_name, matrix_order, w0_order_text, multiplier, origin in rules:
        w0_order = int(w0_order_text)
        expected = EXPECTED_C0 if matrix_name == "C0" else EXPECTED_B1
        for matrix_key in sorted(expected):
            row, column = (int(value) for value in matrix_key.split(","))
            matrix_terms = _matrix_terms(c, matrix_cell, matrix_name,
                                          matrix_order, matrix_key,
                                          index=source_index)
            for matrix_term_index, matrix_term in enumerate(matrix_terms):
                for w0_term_index, w0_term in enumerate(packet[w0_order][column]):
                    output[str(row)].append(
                        _feedback_product(
                            matrix_term, w0_term, multiplier, origin,
                            source_index=source_index,
                            matrix_name=matrix_name,
                            matrix_order=matrix_order,
                            matrix_key=matrix_key,
                            matrix_term_index=matrix_term_index,
                            w0_order=w0_order,
                            w0_term_index=w0_term_index,
                            branch=branch,
                        )
                    )
    return output


def _outer_term(
    chi_term: Mapping[str, Any], feedback_term: Mapping[str, Any], mass: Any,
    *, chi_order: int, feedback_order: int, outer_order: int,
    chi_term_index: int, target_index: int, source_index: int,
    target_x: Any, source_x: Any, field: str,
) -> dict[str, Any]:
    powers = tuple(a + b for a, b in zip(chi_term["powers"], feedback_term["powers"]))
    location_ids = [f"source_cell:{source_index}"]
    if target_index >= 0:
        location_ids.insert(0, f"target_cell:{target_index}")
    ids = _source_ids(chi_term, feedback_term,
                      extra=tuple(location_ids))
    classes = _unique(list(chi_term.get("error_classes", ()))
                      + list(feedback_term.get("error_classes", ())))
    chi_label = "chi_Z" if chi_order == 1 else "chi"
    origin = f"{chi_label} * {feedback_term['origin']}"
    coefficient = -2 * chi_term["coefficient"] * feedback_term["coefficient"] * mass
    origin_record = {
        "operation": "outer_delta_G1_product",
        "source_ids": ids,
        "error_classes": classes,
        "source_classes": classes,
        "uncertainty_class": classes,
        "uncertainty": bool(chi_term.get("uncertainty", False)
                             or feedback_term.get("uncertainty", False)),
        "symbolic_cancellation_claimed": False,
        "chi_order": chi_order,
        "feedback_order": feedback_order,
        "outer_order": outer_order,
    }
    result = {
        "powers": powers,
        "hb_power": powers[0],
        "F0_power": powers[1],
        "Pstar_squared_power": powers[2],
        "Lambda_power": powers[3],
        # This is deliberately signed.  In particular, no abs() or interval
        # majorant is substituted for either source coefficient.
        "coefficient": coefficient,
        "source_ids": ids,
        "error_classes": classes,
        "source_classes": classes,
        "uncertainty_class": classes,
        "uncertainty": origin_record["uncertainty"],
        "origin": origin,
        "origin_record": origin_record,
        "target_cell": target_index,
        "target_cell_index": target_index,
        "source_cell": source_index,
        "source_cell_index": source_index,
        "source_cell_id": feedback_term.get("source_cell_id", source_index),
        "source_term_id": feedback_term.get("W0_source_id", source_index),
        "target_x": target_x,
        "source_x": source_x,
        "target_radius": target_x,
        "source_radius": source_x,
        "outer_order": outer_order,
        "chi_order": chi_order,
        "chi_term_index": chi_term_index,
        "chi_source_id": chi_term.get("source_id", "chi"),
        "feedback_order": feedback_order,
        "feedback_origin": feedback_term.get("origin"),
        "feedback_origin_record": feedback_term.get("origin_record"),
        "matrix_name": feedback_term.get("matrix_name"),
        "matrix_key": feedback_term.get("matrix_key"),
        "matrix_order": feedback_term.get("matrix_order"),
        "matrix_term_index": feedback_term.get("matrix_term_index"),
        "matrix_source_cell": feedback_term.get("matrix_source_cell"),
        "matrix_source_cell_index": feedback_term.get("matrix_source_cell_index", source_index),
        "matrix_source_id": feedback_term.get("matrix_source_id"),
        "row": feedback_term.get("row"),
        "column": feedback_term.get("column"),
        "W0_order": feedback_term.get("W0_order"),
        "W0_term_index": feedback_term.get("W0_term_index"),
        "W0_column": feedback_term.get("W0_column"),
        "W0_column_index": feedback_term.get("W0_column_index"),
        "W0_source_cell": feedback_term.get("W0_source_cell"),
        "W0_source_cell_index": feedback_term.get("W0_source_cell_index", source_index),
        "W0_source_id": feedback_term.get("W0_source_id"),
        "W0_provenance": feedback_term.get("W0_provenance"),
        "branch": feedback_term.get("branch"),
        "leibniz_multiplier": feedback_term.get("leibniz_multiplier", 1),
        "outer_leibniz_multiplier": 1,
        "kernel_mass": mass,
        "kernel_mass_bare": mass,
        "kernel_mass_field": field,
        "correction_factor": -2 * chi_term["coefficient"],
        "principal_kernel": "(-2*chi) T K0",
        "outer_source": "DeltaG1(C0 W0+B1 d_Z W0)",
        "signed_source_term": True,
    }
    return result


def _order_view(values: Mapping[str, Mapping[int, Sequence[Mapping[str, Any]]]]) -> dict[str, dict[str, list[dict[str, Any]]]]:
    return {field: {f"Z{order}": list(values[field][order])
                    for order in OUTPUT_ORDERS}
            for field in values}


def _build_branch(
    c: Any, outer_cells: Any, matrix_cells: Any, chi_rows: Mapping[int, Sequence[Mapping[str, Any]]],
    *, branch: str, include_d_x_F1: bool,
) -> dict[str, Any]:
    outer = _normalise_cells(outer_cells, f"{branch} outer cells")
    matrices = _normalise_cells(matrix_cells, "matrix cells")
    outer_x = [_cell_x(c, cell, index, f"{branch} outer cell")
               for index, cell in enumerate(outer)]
    matrix_x = [_cell_x(c, cell, index, "matrix cell")
                for index, cell in enumerate(matrices)]
    for index, (left, right) in enumerate(zip(outer_x, matrix_x)):
        if not _same_interval(left, right):
            raise ValueError(f"saved {branch} cell and matrix cell coordinates differ at {index}")

    packets: list[dict[int, dict[int, list[dict[str, Any]]]]] = []
    for index, cell in enumerate(outer):
        packet = _cell_source_packet(cell, branch)
        packets.append(_packet_rows(c, packet, label=f"{branch} cell {index}", branch=branch))
    for index, cell in enumerate(matrices):
        _validate_matrix(c, cell, "C0", index=index)
        _validate_matrix(c, cell, "B1", index=index)

    fields = ("F1", "d_x_F1") if include_d_x_F1 else ("F1",)
    # The source products do not depend on the target radius.  Build them once
    # per source cell so the q5 ledger can be consumed without repeating a
    # large term packet for every later target cell.
    source_feedback: list[dict[int, dict[str, list[dict[str, Any]]]]] = []
    source_feedback_q5: list[dict[str, list[dict[str, Any]]]] = []
    source_unweighted_q5: list[dict[str, list[dict[str, Any]]]] = []
    unit_mass = c.mpf(1)
    for source_index in range(CELL_COUNT):
        feedback = {
            0: _feedback_terms(c, matrices[source_index], packets[source_index],
                               source_index=source_index, branch=branch, order=0),
            1: _feedback_terms(c, matrices[source_index], packets[source_index],
                               source_index=source_index, branch=branch, order=1),
        }
        source_feedback.append(feedback)
        source_feedback_q5.append({
            f"Z{order}": feedback[order]["5"] for order in OUTPUT_ORDERS
        })
        q5: dict[str, list[dict[str, Any]]] = {"Z0": [], "Z1": []}
        source_x = outer_x[source_index]
        for outer_order in OUTPUT_ORDERS:
            outer_rules = (OUTER_PRODUCT_RULE_Z0 if outer_order == 0
                           else OUTER_PRODUCT_RULE_Z1)
            for chi_order, feedback_order, _, _ in outer_rules:
                for chi_term_index, chi_term in enumerate(chi_rows[chi_order]):
                    for feedback_term in feedback[feedback_order]["5"]:
                        q5[f"Z{outer_order}"].append(
                            _outer_term(
                                chi_term, feedback_term, unit_mass,
                                chi_order=chi_order,
                                feedback_order=feedback_order,
                                outer_order=outer_order,
                                chi_term_index=chi_term_index,
                                target_index=-1,
                                source_index=source_index,
                                target_x=None,
                                source_x=source_x,
                                field="unweighted",
                            )
                        )
        for values in q5.values():
            for term in values:
                term["kernel_mass"] = None
                term["kernel_mass_bare"] = None
                term["kernel_mass_field"] = None
                term["unweighted_outer_source"] = True
        source_unweighted_q5.append(q5)

    target_cells: list[dict[str, Any]] = []
    for target_index, target_x in enumerate(outer_x):
        aggregate: dict[str, dict[int, list[dict[str, Any]]]] = {
            field: {order: [] for order in OUTPUT_ORDERS} for field in fields
        }
        by_source: dict[str, Any] = {}
        for source_index in range(target_index + 1):
            source_x = outer_x[source_index]
            source_left = _degenerate_endpoint(c, source_x, 0)
            diagonal = source_index == target_index
            source_right = target_x if diagonal else _degenerate_endpoint(c, source_x, 1)
            masses = cell_kernel_masses(c, target_x, source_left, source_right,
                                        triangular=diagonal)
            feedback = source_feedback[source_index]
            correction: dict[str, dict[int, list[dict[str, Any]]]] = {
                field: {order: [] for order in OUTPUT_ORDERS} for field in fields
            }
            for outer_order in OUTPUT_ORDERS:
                outer_rules = (OUTER_PRODUCT_RULE_Z0 if outer_order == 0
                               else OUTER_PRODUCT_RULE_Z1)
                for chi_order, feedback_order, outer_multiplier, _ in outer_rules:
                    for chi_term_index, chi_term in enumerate(chi_rows[chi_order]):
                        # The outer F1 leaf consumes the q5/source row only.
                        # Rows q4 and q6 remain available in the unweighted
                        # feedback ledger but are not fed through this radial
                        # kernel.
                        for feedback_term in feedback[feedback_order]["5"]:
                            for field in fields:
                                correction[field][outer_order].append(
                                    _outer_term(
                                        chi_term, feedback_term,
                                        masses[field],
                                        chi_order=chi_order,
                                        feedback_order=feedback_order,
                                        outer_order=outer_order,
                                        chi_term_index=chi_term_index,
                                        target_index=target_index,
                                        source_index=source_index,
                                        target_x=target_x,
                                        source_x=source_x,
                                        field=field,
                                    )
                                )
            order_view = _order_view(correction)
            by_source[str(source_index)] = {
                "source_index": source_index,
                "target_index": target_index,
                "target_x": target_x,
                "source_x": source_x,
                "target_radius": target_x,
                "source_radius": source_x,
                "triangular": diagonal,
                "axis_continuation": masses["axis_continuation"],
                "kernel_masses": masses,
                "q5_source_cell_index": source_index,
                "outer_delta_G1": order_view,
                "F1": order_view["F1"],
                "d_x_F1": order_view.get("d_x_F1", {}),
            }
            for field in fields:
                for order in OUTPUT_ORDERS:
                    aggregate[field][order].extend(correction[field][order])

        aggregate_view = _order_view(aggregate)
        target_cells.append({
            "index": target_index,
            "target_index": target_index,
            "x": target_x,
            "target_x": target_x,
            "target_radius": target_x,
            "t": target_x / 2,
            "source_indices": list(range(target_index + 1)),
            "source_j_le_target_i": True,
            "triangular_diagonal_source": target_index,
            "axis_endpoint_continuation": target_index == 0,
            "branch": branch,
            "outer_delta_G1": aggregate_view,
            "F1": aggregate_view["F1"],
            "d_x_F1": aggregate_view.get("d_x_F1", {}),
            "outer_delta_G1_by_source_cell": by_source,
        })

    return {
        "branch": branch,
        "target_cells": target_cells,
        "source_j_le_target_i": True,
        "triangular_diagonal_source": True,
        "axis_endpoint_continuation": True,
        "saved_outer_cells_reused": True,
        "saved_W0_reused": True,
        "W0_reconstructed": False,
        "W0_computed": False,
        "matrix_orders": ["0", "1"],
        "input_W0_ordinary_Z_orders": [0, 1, 2],
        "ordinary_Z_orders": [0, 1],
        "output_Z_orders": [0, 1],
        "fields": list(fields),
        "source_power_basis": list(SOURCE_POWER_BASIS),
        "source_ids_retained": True,
        "formal_powers_retained": True,
        "signed_coefficients_retained": True,
        "empty_delta_columns_preserved": branch == "delta_W0",
        "delta_empty_columns": [2, 3, 4, 6] if branch == "delta_W0" else [],
        "delta_signed_columns": [1, 5] if branch == "delta_W0" else [],
        "green_kernel_mass_definition": "directed source-cell integral of T K0",
        "green_kernel_masses_exact": True,
        "principal_kernel": "(-2*chi) T K0",
        "outer_operator": "DeltaG1(C0 W0+B1 d_Z W0)",
        "outer_feedback_only": True,
        "full_W1_claimed": False,
        "W1_claimed": False,
        "complete_W1": False,
        "Kchi_first_correction_only": True,
        "Kchi_higher_words_included": False,
        "allword_tail_added": False,
        "source_terms_recomputed": False,
        "old_builder_called": False,
        "axis_log_division_avoided": True,
        "include_d_x_F1": include_d_x_F1,
        "chi_jets": {f"Z{order}": list(values) for order, values in chi_rows.items()},
        "feedback_q5_products_by_source_cell": {
            str(index): source_feedback_q5[index]
            for index in range(CELL_COUNT)
        },
        "unweighted_outer_delta_source_by_source_cell": {
            str(index): source_unweighted_q5[index]
            for index in range(CELL_COUNT)
        },
    }


def build(
    c: Any,
    outer_cells: Any = None,
    matrix_cells: Any = None,
    chi: Any = None,
    *,
    chi_jets: Any = None,
    source_cells: Any = None,
    saved_outer_cells: Any = None,
    w0_cells: Any = None,
    delta_w0_cells: Any = None,
    delta_W0_cells: Any = None,
    delta_w0: Any = None,
    delta_W0: Any = None,
    branch: str | None = None,
    include_d_x_F1: bool = True,
    include_dx: bool | None = None,
) -> dict[str, Any]:
    """Build saved-W0 outer ``Delta G1`` branches by target/source cell.

    The positional API follows :mod:`n1_core_w1_profile_v2`: ``outer_cells``
    are saved W0 cells and ``matrix_cells`` contain C0/B1.  ``chi`` may be a
    scalar or ordinary Z0/Z1(/Z2) rows; a scalar has an explicit zero chi_Z.
    ``delta_w0_cells`` is optional and has the same six-column shape.  When it
    is supplied, both ``K0_W0`` and ``delta_W0`` are returned under
    ``branches``; each branch remains unmerged.  Set ``branch='delta_W0'`` to
    evaluate only an independent delta packet (using ``outer_cells`` as that
    packet when no delta packet is supplied).
    """

    if not hasattr(c, "mpf") or not hasattr(c, "make_mpf"):
        raise TypeError("c must be the caller's MPIntervalContext")
    if outer_cells is None:
        outer_cells = source_cells if source_cells is not None else saved_outer_cells
    if outer_cells is None:
        outer_cells = w0_cells
    if matrix_cells is None:
        raise ValueError("matrix_cells is required")
    if chi is None:
        chi = chi_jets
    if include_dx is not None:
        include_d_x_F1 = bool(include_dx)
    delta_input = delta_w0_cells
    if delta_input is None:
        delta_input = delta_W0_cells
    if delta_input is None:
        delta_input = delta_w0
    if delta_input is None:
        delta_input = delta_W0

    mode = None if branch is None else str(branch).lower().replace("-", "_")
    if mode in {"k0", "w0", "k0_w0", "ordinary", "ordinary_w0"}:
        selected = ("K0_W0",)
    elif mode in {"delta", "delta_w0", "deltaw0", "correction", "delta_w0_only"}:
        selected = ("delta_W0",)
    elif mode is None or mode in {"both", "all"}:
        selected = ("K0_W0", "delta_W0") if delta_input is not None else ("K0_W0",)
    else:
        raise ValueError(f"unknown outer feedback branch {branch!r}")

    chi_rows = _chi_rows(c, chi)
    branches: dict[str, dict[str, Any]] = {}
    for name in selected:
        if name == "K0_W0":
            if outer_cells is None:
                raise ValueError("saved outer W0 cells are required for K0_W0")
            branch_cells = outer_cells
        else:
            branch_cells = delta_input if delta_input is not None else outer_cells
            if branch_cells is None:
                raise ValueError("delta_W0 cells or outer cells are required for delta_W0")
        branches[name] = _build_branch(
            c, branch_cells, matrix_cells, chi_rows,
            branch=name, include_d_x_F1=include_d_x_F1,
        )

    primary_name = "K0_W0" if "K0_W0" in branches else "delta_W0"
    primary = dict(branches[primary_name])
    primary.update({
        "branches": branches,
        "branch_names": list(branches),
        "primary_branch": primary_name,
        "delta_W0_input_supplied": delta_input is not None,
        "same_leaf_dual_branch_capable": True,
        "full_W1_claimed": False,
        "W1_claimed": False,
        "complete_W1": False,
        "physical_field_claimed": False,
    })
    return primary


build_outer_feedback = build


def build_delta_w0(c: Any, outer_cells: Any, matrix_cells: Any, chi: Any = None,
                   **kwargs: Any) -> dict[str, Any]:
    """Convenience wrapper returning the independent ``delta_W0`` branch."""

    kwargs["branch"] = "delta_W0"
    return build(c, outer_cells, matrix_cells, chi, **kwargs)
