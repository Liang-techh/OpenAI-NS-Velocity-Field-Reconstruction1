"""Independent fixed-endpoint principal Green leaf for ``G0 g``.

The low-word source uses three signed forcing rows (``4``, ``5`` and ``6``)
and the six natural fields ``(F1, Uz1, K1, P1, d_x_F1, d_x_Uz1)``.  This
module consumes those rows as caller-owned ordinary ``Z0``/``Z1`` sectors on
the eight canonical source cells ``t=[i/8,(i+1)/8]``.  It applies the exact
fixed-target ``x=2`` source-cell moments of the principal ``G0`` kernel.  The
endpoint is an independent target ``t=1``; the last source cell is never
treated as an endpoint value or midpoint sample.

The accepted low-word spelling is the order-major packet used by
``n1_core_word1_z2_v2``::

    {"0": {"4": [terms], "5": [terms], "6": [terms]},
     "1": {"4": [terms], "5": [terms], "6": [terms]}}

An equivalent ``{"Z0": ..., "Z1": ...}`` spelling and row-major spelling
are accepted.  A source cell may wrap the packet in ``g``, ``forcing``,
``source_rows`` or ``ordinary_Z``.  Terms are copied without collapsing their
formal powers or source identities.  No producer, receipt, old builder, or
source value is imported or called here.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any


SCHEMA = "n1_g0_six_exit_v2"
SOURCE_SCHEMA = "n1_core_word1_z2_v2"
CELL_COUNT = 8
ORDERS = (0, 1)
INPUT_ORDERS = (0, 1, 2)
ROWS = ("4", "5", "6")
FORCING_ROWS = ROWS
FIELDS = ("F1", "Uz1", "K1", "P1", "d_x_F1", "d_x_Uz1")
FIELD_ROWS = {
    "F1": "5",
    "Uz1": "6",
    "K1": "6",
    "P1": "4",
    "d_x_F1": "5",
    "d_x_Uz1": "6",
}
SOURCE_POWER_KEYS = ("F0_power", "Lambda_power")
SOURCE_POWER_BASIS = ("F0", "Lambda")


__all__ = [
    "SCHEMA",
    "SOURCE_SCHEMA",
    "CELL_COUNT",
    "ORDERS",
    "INPUT_ORDERS",
    "ROWS",
    "FORCING_ROWS",
    "FIELDS",
    "FIELD_ROWS",
    "SOURCE_POWER_KEYS",
    "SOURCE_POWER_BASIS",
    "g0_cell_kernel_masses",
    "cell_kernel_masses",
    "endpoint_kernel_masses",
    "build_g0g",
    "build",
    "build_core_exit",
    "build_endpoint",
    "build_fixed_exit",
    "build_x2",
    "g0g_exit",
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


def _order_key(value: Any, label: str) -> int:
    if isinstance(value, bool):
        raise ValueError(f"boolean is not an ordinary-Z order in {label}")
    if isinstance(value, int):
        return value
    if isinstance(value, str):
        text = value.strip().upper()
        if text.startswith("Z"):
            text = text[1:]
        try:
            return int(text)
        except ValueError as exc:
            raise ValueError(f"invalid ordinary-Z order {value!r} in {label}") from exc
    raise ValueError(f"invalid ordinary-Z order {value!r} in {label}")


def _interval(c: Any, value: Any) -> Any:
    """Coerce a value while retaining directed endpoints only."""

    if hasattr(value, "_mpi_") and hasattr(value, "a") and hasattr(value, "b"):
        return value
    if hasattr(value, "a") and hasattr(value, "b"):
        lo = value.a._mpi_[0] if hasattr(value.a, "_mpi_") else value.a
        hi = value.b._mpi_[1] if hasattr(value.b, "_mpi_") else value.b
        return c.make_mpf((lo, hi))
    if isinstance(value, Mapping):
        if "lower" in value and "upper" in value:
            return c.mpf([value["lower"], value["upper"]])
        if "lo" in value and "hi" in value:
            return c.mpf([value["lo"], value["hi"]])
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        if len(value) == 2:
            return c.mpf([value[0], value[1]])
    return c.mpf(value)


def _endpoint(c: Any, value: Any, upper: bool) -> Any:
    value = _interval(c, value)
    raw = value._mpi_[1 if upper else 0]
    return c.make_mpf((raw, raw))


def _same_interval(left: Any, right: Any) -> bool:
    try:
        return left._mpi_ == right._mpi_
    except AttributeError:
        return left == right


def _canonical_t(c: Any, index: int) -> Any:
    left = c.mpf(index) / CELL_COUNT
    right = c.mpf(index + 1) / CELL_COUNT
    return c.make_mpf((left._mpi_[0], right._mpi_[1]))


def _cell_t(c: Any, cell: Any, index: int, label: str) -> Any:
    expected = _canonical_t(c, index)
    if not isinstance(cell, Mapping):
        return expected
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
    if candidate is None:
        actual = expected
    else:
        actual = _interval(c, candidate)
        if not _same_interval(actual, expected):
            raise ValueError(f"{label} {index} is not canonical t=[i/8,(i+1)/8]")
    if "rho" in cell and not _same_interval(_interval(c, cell["rho"]), 4 * actual * actual):
        raise ValueError(f"{label} {index} has inconsistent rho=4*t^2")
    if "x" in cell and not _same_interval(_interval(c, cell["x"]), 2 * actual):
        raise ValueError(f"{label} {index} has inconsistent x=2*t")
    return actual


def _order_map(raw: Any, required: Sequence[int], label: str) -> dict[int, Any]:
    if not isinstance(raw, Mapping):
        raise TypeError(f"{label} must be an ordinary-Z order mapping")
    result: dict[int, Any] = {}
    for key, value in raw.items():
        order = _order_key(key, label)
        if order in result:
            raise ValueError(f"duplicate {label} order {order}")
        result[order] = value
    if not set(required).issubset(result):
        raise ValueError(f"{label} must provide ordinary orders {list(required)}")
    return result


def _row_map(raw: Any, label: str) -> dict[str, Any]:
    if not isinstance(raw, Mapping):
        raise TypeError(f"{label} must provide forcing rows 4, 5, 6")
    result: dict[str, Any] = {}
    for key, value in raw.items():
        try:
            row = str(int(key))
        except (TypeError, ValueError) as exc:
            raise ValueError(f"invalid forcing row {key!r} in {label}") from exc
        if row in result:
            raise ValueError(f"duplicate forcing row {row} in {label}")
        result[row] = value
    if set(result) != set(ROWS):
        raise ValueError(f"{label} must provide exactly forcing rows 4, 5, 6")
    return result


def _term_list(raw: Any, label: str) -> list[Mapping[str, Any]]:
    if isinstance(raw, Mapping) and ("coefficient" in raw or "powers" in raw):
        values = [raw]
    elif isinstance(raw, Sequence) and not isinstance(raw, (str, bytes, bytearray)):
        values = list(raw)
    else:
        values = [raw]
    result: list[Mapping[str, Any]] = []
    for index, value in enumerate(values):
        if not isinstance(value, Mapping):
            raise TypeError(f"{label}[{index}] must be a formal source term")
        if "coefficient" not in value:
            raise ValueError(f"{label}[{index}] is missing coefficient")
        if "powers" in value:
            powers = value["powers"]
            if (
                not isinstance(powers, Sequence)
                or isinstance(powers, (str, bytes, bytearray))
                or not all(isinstance(power, int) and not isinstance(power, bool) for power in powers)
            ):
                raise ValueError(f"{label}[{index}].powers must contain integer formal powers")
        else:
            f0 = value.get("F0_power", value.get("f0_power"))
            lam = value.get("Lambda_power", value.get("lambda_power"))
            if not isinstance(f0, int) or isinstance(f0, bool):
                raise ValueError(f"{label}[{index}] is missing integer F0_power")
            if not isinstance(lam, int) or isinstance(lam, bool):
                raise ValueError(f"{label}[{index}] is missing integer Lambda_power")
        result.append(value)
    return result


def _packet_wrapper(cell: Any) -> Any:
    if not isinstance(cell, Mapping):
        return cell
    for key in (
        "g",
        "forcing",
        "source_rows",
        "source_jets",
        "ordinary_Z",
        "low_source",
        "global_low_source",
    ):
        if key in cell:
            return cell[key]
    return cell


def _normalise_packet(raw: Any, label: str) -> dict[int, dict[str, list[Mapping[str, Any]]]]:
    """Normalize order-major or row-major ``g`` sectors without Z3 leakage."""

    raw = _packet_wrapper(raw)
    if not isinstance(raw, Mapping):
        raise TypeError(f"{label} must provide ordinary-Z forcing rows")

    # The n1_core_word1_z2_v2 spelling is order -> row -> terms.  Explicit
    # order keys are required; a bare row map from the separate Z3 packet is
    # intentionally rejected rather than silently promoted to Z0.
    parsed_orders: dict[int, Any] = {}
    order_like = True
    for key, value in raw.items():
        try:
            order = _order_key(key, label)
        except ValueError:
            order_like = False
            break
        if order in parsed_orders:
            raise ValueError(f"duplicate {label} order {order}")
        parsed_orders[order] = value

    result: dict[int, dict[str, list[Mapping[str, Any]]]] = {}
    if order_like and 0 in parsed_orders and 1 in parsed_orders:
        for order in INPUT_ORDERS:
            if order not in parsed_orders:
                # Z2 is optional for this endpoint; only Z0 and Z1 are used.
                continue
            rows = _row_map(parsed_orders[order], f"{label} Z{order}")
            result[order] = {
                row: _term_list(rows[row], f"{label} Z{order} row {row}")
                for row in ROWS
            }
        return result

    # Equivalent row -> order -> terms spelling.  It must still expose both
    # ordinary Z orders, so a Z3-only source cannot enter this leaf.
    rows = _row_map(raw, label)
    for row in ROWS:
        row_orders = _order_map(rows[row], ORDERS, f"{label} row {row}")
        for order in ORDERS:
            result.setdefault(order, {})[row] = _term_list(
                row_orders[order], f"{label} row {row} Z{order}"
            )
    return result


def _looks_like_global_packet(value: Any) -> bool:
    if not isinstance(value, Mapping):
        return False
    keys = set()
    for key in value:
        try:
            keys.add(_order_key(key, "source"))
        except ValueError:
            return False
    return 0 in keys and 1 in keys


def _as_cells(value: Any, label: str) -> list[Any]:
    if isinstance(value, Mapping):
        for wrapper in ("cells", "source_cells", "forcing_cells", "g_cells"):
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
    return cells


def _bounds(c: Any, value: Any) -> tuple[Any, Any]:
    value = _interval(c, value)
    lo, hi = value._mpi_
    return c.make_mpf((lo, lo)), c.make_mpf((hi, hi))


def _nonnegative(c: Any, value: Any, label: str) -> Any:
    value = _interval(c, value)
    lo, hi = value._mpi_
    if c.make_mpf((hi, hi)) < 0:
        raise ValueError(f"{label} has a negative upper bound")
    if c.make_mpf((lo, lo)) < 0:
        lo = c.mpf(0)._mpi_[0]
    return c.make_mpf((lo, hi))


def _primitive_uz(c: Any, target_t: Any, source_t: Any) -> Any:
    source_t = _interval(c, source_t)
    if _bounds(c, source_t)[0] == 0:
        return c.mpf(0)
    if _bounds(c, target_t)[1] <= 0:
        raise ValueError("positive fixed target t is required")
    # d/ds [2 s^2 log(T/s)+s^2] = 4 s log(T/s), including the
    # integrable logarithmic endpoint by the explicit s=0 branch above.
    return 2 * source_t * source_t * c.log(_interval(c, target_t) / source_t) + source_t * source_t


def g0_cell_kernel_masses(
    c: Any,
    target_t: Any,
    source_left_t: Any,
    source_right_t: Any,
) -> dict[str, Any]:
    """Return exact directed fixed-target ``G0`` moments for one source cell.

    The physical source measure is ``dx=2 dt``.  Consequently the pressure
    row has the weighted mass ``P1=2*(b-a)``; it is not a constant copied from
    a pointwise target kernel.  All other rows are endpoint primitives in the
    normalized radius ``t=x/2``.
    """

    T = _interval(c, target_t)
    a = _interval(c, source_left_t)
    b = _interval(c, source_right_t)
    Tlo, Thi = _bounds(c, T)
    alo, ahi = _bounds(c, a)
    blo, bhi = _bounds(c, b)
    if alo < 0 or ahi < alo or blo < alo or bhi < blo:
        raise ValueError("source endpoints must be ordered and nonnegative")
    if bhi > Thi:
        raise ValueError("source cell extends beyond fixed target")
    if Thi <= 0:
        raise ValueError("positive fixed target t is required")

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
    f1 = _nonnegative(c, f1, "G0 F1 mass")
    uz1 = _nonnegative(c, uz1, "G0 Uz1 mass")
    p1 = _nonnegative(c, p1, "G0 P1 mass")
    dxf1 = _nonnegative(c, dxf1, "G0 d_x_F1 mass")
    dxuz1 = _nonnegative(c, dxuz1, "G0 d_x_Uz1 mass")
    return {
        "F1": f1,
        "Uz1": uz1,
        "K1": -f1,
        "P1": p1,
        "d_x_F1": dxf1,
        "d_x_Uz1": dxuz1,
        "target_t": T,
        "source_left_t": a,
        "source_right_t": b,
        "target_x": 2 * T,
        "source_left_x": 2 * a,
        "source_right_x": 2 * b,
        "principal_kernel": "G0",
        "negative_kernel_K1": True,
        "physical_source_measure": "dx=2*dt",
        "P1_weighted_mass": "2*(source_right_t-source_left_t)",
        "exact_source_cell_mass": True,
    }


cell_kernel_masses = g0_cell_kernel_masses
endpoint_kernel_masses = g0_cell_kernel_masses


def _source_ids(term: Mapping[str, Any]) -> list[Any]:
    values: list[Any] = []
    for key in ("source_ids", "source_id", "source_cell_id", "source_cell"):
        value = term.get(key)
        if value is None:
            continue
        if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
            values.extend(value)
        else:
            values.append(value)
    origin = term.get("origin")
    if isinstance(origin, Mapping):
        value = origin.get("source_ids", origin.get("source_id"))
        if value is not None:
            if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
                values.extend(value)
            else:
                values.append(value)
    unique: list[Any] = []
    seen: set[Any] = set()
    for value in values:
        try:
            marker = (type(value).__name__, value)
            hash(marker)
        except TypeError:
            marker = (type(value).__name__, repr(value))
        if marker not in seen:
            seen.add(marker)
            unique.append(value)
    return unique


def _scale_term(
    term: Mapping[str, Any],
    coefficient: Any,
    *,
    field: str,
    row: str,
    order: int,
    source_index: int,
    source_term_index: int,
    mass: Any,
) -> dict[str, Any]:
    result = dict(term)
    result["coefficient"] = coefficient * mass
    ids = _source_ids(term)
    if not ids:
        ids = [f"g0g:{source_index}:Z{order}:row{row}:term{source_term_index}"]
    if ids and "source_ids" not in result:
        result["source_ids"] = list(ids)
    if ids and "source_id" not in result:
        result["source_id"] = ids[0]
    result.update({
        "field": field,
        "row": row,
        "source_row": row,
        "ordinary_Z_order": order,
        "source_cell_index": source_index,
        "forcing_source_cell": source_index,
        "source_term_index": source_term_index,
        "kernel_mass": mass,
        "kernel_mass_bare": mass,
        "kernel_mass_field": field,
        "principal_kernel": "G0",
        "negative_kernel_K1": field == "K1",
        "exact_source_cell_mass": True,
        "signed_source_term": True,
    })
    return result


def _empty_fields() -> dict[str, list[dict[str, Any]]]:
    return {field: [] for field in FIELDS}


def _resolve_source_cells(
    source_cells: Any,
    low_source: Any,
    *,
    forcing_cells: Any,
    forcing_source: Any,
    global_low_source: Any,
    source: Any,
    cells: Any,
    outer_cells: Any,
) -> tuple[list[Any], Any, str]:
    # Handle the two-positional spelling first.  In
    # ``build(c, global_low_source, source_cells)`` the second positional
    # argument is the cells, not a second global packet.
    if (
        source_cells is not None
        and _looks_like_global_packet(source_cells)
        and low_source is not None
        and not _looks_like_global_packet(low_source)
        and not any(value is not None for value in (forcing_cells, source, cells, outer_cells))
    ):
        source_cells, low_source = low_source, source_cells
    aliases = [
        value for value in (forcing_cells, source, cells, outer_cells)
        if value is not None
    ]
    if len(aliases) > 1:
        raise TypeError("source_cells aliases were supplied more than once")
    if source_cells is None:
        source_cells = aliases[0] if aliases else None
    elif aliases:
        raise TypeError("source_cells and a source_cells alias were supplied twice")

    global_aliases = [
        value for value in (low_source, forcing_source, global_low_source)
        if value is not None
    ]
    if len(global_aliases) > 1:
        raise TypeError("low_source aliases were supplied more than once")
    global_packet = global_aliases[0] if global_aliases else None

    # Support the natural two-positional spelling build(c, low_source,
    # source_cells) while keeping build(c, source_cells) concise.
    if source_cells is not None and _looks_like_global_packet(source_cells):
        if global_packet is not None:
            raise TypeError("two global low-source packets were supplied")
        if low_source is None or not isinstance(low_source, (Mapping, Sequence)):
            raise TypeError("a global low source requires eight source cells")
        global_packet, source_cells = source_cells, low_source
    if source_cells is None and global_packet is None:
        raise ValueError("eight caller-owned source cells or a global low source are required")

    if source_cells is None:
        source_cells = [None for _ in range(CELL_COUNT)]
    raw_cells = _as_cells(source_cells, "source cell")
    scope = "per_cell"
    if global_packet is not None:
        scope = "global_low_source_reused_per_cell"
    return raw_cells, global_packet, scope


def build_g0g(
    c: Any,
    source_cells: Any = None,
    low_source: Any = None,
    *,
    forcing_cells: Any = None,
    forcing_source: Any = None,
    global_low_source: Any = None,
    source: Any = None,
    cells: Any = None,
    outer_cells: Any = None,
) -> dict[str, Any]:
    """Build the fixed ``x=2`` baseline ``G0g`` through ordinary ``Z1``.

    ``source_cells`` is either eight cells each carrying an explicit ``g``
    packet or eight coordinate cells paired with ``low_source``.  The latter
    is the whole-core ``global_low_source`` spelling from the low-word source;
    it is copied as a formal enclosure for each cell and is never evaluated at
    a point.  The direct source-cell form is preferred when the caller owns
    per-cell source sectors.
    """

    if not hasattr(c, "mpf") or not hasattr(c, "make_mpf"):
        raise TypeError("c must provide directed interval arithmetic")
    raw_cells, global_packet, source_scope = _resolve_source_cells(
        source_cells,
        low_source,
        forcing_cells=forcing_cells,
        forcing_source=forcing_source,
        global_low_source=global_low_source,
        source=source,
        cells=cells,
        outer_cells=outer_cells,
    )
    normalized_global = (
        _normalise_packet(global_packet, "global_low_source")
        if global_packet is not None else None
    )

    aggregate = {f"Z{order}": _empty_fields() for order in ORDERS}
    by_source: dict[str, Any] = {}
    source_t_by_cell: dict[str, Any] = {}
    kernel_masses_by_source_cell: dict[str, Any] = {}
    for index, raw_cell in enumerate(raw_cells):
        t = _cell_t(c, raw_cell, index, "source cell")
        source_t_by_cell[str(index)] = t
        packet = normalized_global
        if packet is None:
            packet_raw = _packet_wrapper(raw_cell)
            if isinstance(raw_cell, Mapping):
                packet_raw = _packet_wrapper(raw_cell)
            packet = _normalise_packet(packet_raw, f"source cell {index}")
        # A cell packet may include Z2 for compatibility with the low-word
        # source schema; only ordinary Z0 and Z1 are consumed here.
        for order in ORDERS:
            if order not in packet:
                raise ValueError(f"source cell {index} is missing ordinary Z{order}")
        left = _endpoint(c, t, False)
        right = _endpoint(c, t, True)
        masses = g0_cell_kernel_masses(c, c.mpf(1), left, right)
        kernel_masses_by_source_cell[str(index)] = masses
        cell_output = {f"Z{order}": _empty_fields() for order in ORDERS}
        for order in ORDERS:
            for field in FIELDS:
                row = FIELD_ROWS[field]
                for term_index, term in enumerate(packet[order][row]):
                    coefficient = _interval(c, term["coefficient"])
                    item = _scale_term(
                        term,
                        coefficient,
                        field=field,
                        row=row,
                        order=order,
                        source_index=index,
                        source_term_index=term_index,
                        mass=masses[field],
                    )
                    cell_output[f"Z{order}"][field].append(item)
                    aggregate[f"Z{order}"][field].append(item)
        cell_record = {
            "source_index": index,
            "source_t": t,
            "source_x": 2 * t,
            "target_t": c.mpf(1),
            "target_x": c.mpf(2),
            "kernel_masses": masses,
            "G0g": cell_output,
            "G0": cell_output,
            "fields": cell_output,
            "forcing_rows": list(ROWS),
            "source_scope": source_scope,
        }
        by_source[str(index)] = cell_record

    baseline = {
        "schema": "n1_g0_six_exit_baseline_G0g_v1",
        "branch": "baseline_G0g",
        "G0g_exit_by_source_cell": by_source,
        "G0_exit_by_source_cell": by_source,
        "exit_by_source_cell": by_source,
        "kernel_masses_by_source_cell": kernel_masses_by_source_cell,
        "source_t_by_cell": source_t_by_cell,
        "aggregate": aggregate,
        "G0g_Z0": aggregate["Z0"],
        "G0g_Z1": aggregate["Z1"],
        "fields": list(FIELDS),
        "forcing_rows": list(ROWS),
        "ordinary_Z_orders": [0, 1],
        "complete_six": True,
        "source_cells": CELL_COUNT,
        "source_ids_retained": True,
        "formal_powers_retained": True,
        "signed_coefficients_retained": True,
        "green_kernel_masses_exact": True,
    }
    result = {
        "schema": SCHEMA,
        "source_schema": SOURCE_SCHEMA,
        "coordinate": "rho=4, x=2, t=1, real Z in [-1,1]",
        "exit_target": {
            "x": c.mpf(2),
            "t": c.mpf(1),
            "rho": c.mpf(4),
            "fixed": True,
            "independent": True,
            "not_last_source_cell_proxy": True,
        },
        "baseline_G0g": baseline,
        "baseline": baseline,
        "G0g": baseline,
        "G0": baseline,
        "G0g_exit_by_source_cell": by_source,
        "G0_exit_by_source_cell": by_source,
        "exit_by_source_cell": by_source,
        "G0g_Z0_core_exit": aggregate["Z0"],
        "G0g_Z1_core_exit": aggregate["Z1"],
        "G0_Z0_core_exit": aggregate["Z0"],
        "G0_Z1_core_exit": aggregate["Z1"],
        "ordinary_Z_orders": [0, 1],
        "fields": list(FIELDS),
        "forcing_rows": list(ROWS),
        "source_cell_count": CELL_COUNT,
        "source_scope": source_scope,
        "source_power_basis": list(SOURCE_POWER_BASIS),
        "source_power_keys": list(SOURCE_POWER_KEYS),
        "source_ids_retained": True,
        "formal_powers_retained": True,
        "signed_coefficients_retained": True,
        "green_kernel_masses_exact": True,
        "fixed_x2_endpoint": True,
        "W0_reconstructed": False,
        "ancestor_helpers_recomputed": False,
        "old_builder_called": False,
        "source_values_selected": False,
        "interval_midpoint_used": False,
        "source_exponentials_materialized": False,
        "first_collar_added": False,
        "allword_tail_added": False,
        "m15_computed": False,
    }
    return result


build = build_g0g
build_core_exit = build_g0g
build_endpoint = build_g0g
build_fixed_exit = build_g0g
build_x2 = build_g0g
g0g_exit = build_g0g
