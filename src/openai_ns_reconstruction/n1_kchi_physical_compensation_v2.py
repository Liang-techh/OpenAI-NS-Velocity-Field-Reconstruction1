"""Finite physical ``+2 chi E51`` compensation sectors for the n1 inlet.

This module is intentionally a small, source-bound compensation leaf.  It
consumes two caller-owned eight-cell packets containing only the saved
``F1`` ordinary rows ``Z0``/``Z1``: the ordinary ``W0`` packet and a signed
``Delta W0`` packet.  It forms

``P = +2 chi E51``

with the ordinary Leibniz rows
``P_Z0 = 2 chi V`` and ``P_Z1 = 2(chi_Z V + chi V_Z)``.  It then evaluates
four finite branches independently:

``G0 P W0``, ``G0 P DeltaW0``, ``D P W0``, and ``D P DeltaW0``.

Here ``G0`` is the exact six-cell radial ``F1`` Green mass and ``D`` is the
first physical correction ``(-2 chi) T K0``.  The public result includes
target-cell profiles (including axis and diagonal continuations) and a
separate fixed ``x=2`` exit.  Each branch carries only ``F1`` and, when
requested, ``d_x_F1``; no other matrix row, residual C0, all-word tail, or
complete Kchi feedback is inferred.

The four formal powers ``(hb, F0_anchor, Pstar_squared, Lambda)`` and every
caller source identity are retained.  The two chi factors are recorded at
separate levels in the D branches.  Symmetric/uncertainty metadata is copied
without adding it to a signed finite term.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any


SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
CELL_COUNT = 8
ORDERS = (0, 1)
FIELDS = ("F1", "d_x_F1")


__all__ = [
    "SOURCE_POWER_BASIS", "CELL_COUNT", "ORDERS", "FIELDS",
    "g0_cell_kernel_masses", "d_cell_kernel_masses",
    "build_compensation", "build", "physical_compensation",
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
    if isinstance(value, (list, tuple)) and len(value) == 2:
        return c.mpf([value[0], value[1]])
    return c.mpf(value)


def _bounds(c: Any, value: Any) -> tuple[Any, Any]:
    value = _interval(c, value)
    raw = value._mpi_
    return c.make_mpf((raw[0], raw[0])), c.make_mpf((raw[1], raw[1]))


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
    if candidate is None:
        return expected
    candidate = _interval(c, candidate)
    if not _same_interval(candidate, expected):
        raise ValueError(f"{label} {index} is not canonical t=[i/8,(i+1)/8]")
    return candidate


def _as_items(value: Any) -> list[Any]:
    if value is None:
        return []
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


def _ids(raw: Any) -> list[Any]:
    if not isinstance(raw, Mapping):
        return []
    result = list(_as_items(raw.get("source_ids", raw.get("source_id"))))
    nested = _origin(raw)
    result.extend(_as_items(nested.get("source_ids", nested.get("source_id"))))
    return _unique(result)


def _classes(raw: Any) -> list[Any]:
    if not isinstance(raw, Mapping):
        return []
    result = list(_as_items(raw.get("error_classes", raw.get("error_class"))))
    nested = _origin(raw)
    result.extend(_as_items(nested.get("error_classes", nested.get("error_class"))))
    return _unique(result)


def _uncertainty(raw: Any) -> bool:
    if not isinstance(raw, Mapping):
        return False
    nested = _origin(raw)
    return bool(raw.get("uncertainty", raw.get("is_uncertainty", False))
                or raw.get("symmetric_error", False)
                or nested.get("uncertainty", nested.get("is_uncertainty", False)))


def _powers(raw: Any, label: str) -> tuple[int, int, int, int]:
    if not isinstance(raw, Mapping):
        return (0, 0, 0, 0)
    powers = raw.get("powers")
    if powers is None:
        powers = (
            raw.get("hb_power", 0),
            raw.get("F0_power", raw.get("F0_anchor_power", 0)),
            raw.get("Pstar_squared_power", raw.get("Pstar_power", 0)),
            raw.get("Lambda_power", raw.get("lambda_power", 0)),
        )
    if (not isinstance(powers, Sequence)
            or isinstance(powers, (str, bytes, bytearray))
            or len(powers) != 4
            or any(isinstance(item, bool) or not isinstance(item, int)
                   for item in powers)):
        raise ValueError(f"{label} must carry four integer formal powers")
    return tuple(int(item) for item in powers)


def _term(c: Any, raw: Any, *, label: str, default_id: Any,
          error_only: bool = False) -> dict[str, Any]:
    if isinstance(raw, Mapping):
        if "coefficient" not in raw:
            raise ValueError(f"{label} is missing coefficient")
        result = dict(raw)
        coefficient = raw["coefficient"]
    else:
        result = {}
        coefficient = raw
    coefficient = _interval(c, coefficient)
    powers = _powers(raw, label)
    ids = _ids(raw) or [default_id]
    classes = _classes(raw)
    uncertain = _uncertainty(raw) or error_only
    result.update({
        "powers": powers,
        "hb_power": powers[0],
        "F0_power": powers[1],
        "Pstar_squared_power": powers[2],
        "Lambda_power": powers[3],
        "coefficient": coefficient,
        "source_ids": _unique(ids),
        "error_classes": list(classes),
        "uncertainty": uncertain,
        "origin": dict(_origin(raw)),
    })
    result["origin"].update({
        "source_ids": list(result["source_ids"]),
        "error_classes": list(classes),
        "uncertainty": uncertain,
    })
    if error_only:
        result["error_only"] = True
        result["symmetric_error"] = True
        result["error_classes"] = _unique(list(classes) + ["symmetric_R"])
        result["origin"]["error_classes"] = list(result["error_classes"])
    result.setdefault("source_id", result["source_ids"][0])
    return result


def _terms(c: Any, raw: Any, *, label: str, default_id: Any,
           error_only: bool = False) -> list[dict[str, Any]]:
    if isinstance(raw, Mapping) and ("coefficient" in raw or "powers" in raw):
        values = [raw]
    elif isinstance(raw, Sequence) and not isinstance(raw, (str, bytes, bytearray)):
        values = list(raw)
    else:
        values = [raw]
    return [_term(c, value, label=f"{label}[{index}]",
                  default_id=default_id, error_only=error_only)
            for index, value in enumerate(values)]


def _order_map(raw: Any, required: Sequence[int], label: str) -> dict[int, Any]:
    if isinstance(raw, Mapping):
        result: dict[int, Any] = {}
        for key, value in raw.items():
            text = key.upper().removeprefix("Z") if isinstance(key, str) else key
            try:
                order = int(text)
            except (TypeError, ValueError) as exc:
                raise ValueError(f"invalid {label} order {key!r}") from exc
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


def _packet_wrapper(cell: Mapping[str, Any]) -> Any:
    for key in ("W0_ordinary_Z", "ordinary_Z", "W0", "source_rows",
                "source_jets", "delta_W0", "delta_w0", "R", "error"):
        if key in cell:
            return cell[key]
    return cell


def _normalise_source_cells(c: Any, cells: Any, *, label: str,
                            source_column: int = 1,
                            error_only: bool = False) -> tuple[list[Any], list[dict[int, list[dict[str, Any]]]]]:
    if isinstance(cells, Mapping):
        raw_cells = [_lookup(cells, index) for index in range(CELL_COUNT)]
    elif isinstance(cells, Sequence) and not isinstance(cells, (str, bytes, bytearray)):
        raw_cells = list(cells)
    else:
        raise TypeError(f"{label} must contain eight cells")
    if len(raw_cells) != CELL_COUNT:
        raise ValueError(f"eight {label}s are required")
    radii: list[Any] = []
    rows: list[dict[int, list[dict[str, Any]]]] = []
    for index, cell in enumerate(raw_cells):
        if not isinstance(cell, Mapping):
            raise TypeError(f"{label} {index} must be a mapping")
        radii.append(_cell_t(c, cell, index, label))
        direct = None
        for key in ("F1", "g5", "G5", "row5", "row_5"):
            if key in cell and not (key == "F1" and "W0_ordinary_Z" in cell):
                direct = cell[key]
                break
        if direct is not None:
            order_map = _order_map(direct, ORDERS, f"{label} F1")
            rows.append({
                order: _terms(c, order_map[order],
                              label=f"{label}[{index}] Z{order}",
                              default_id=f"{label}:{index}:Z{order}",
                              error_only=error_only)
                for order in ORDERS
            })
            continue
        packet = _packet_wrapper(cell)
        order_map = _order_map(packet, ORDERS, f"{label} W0")
        result: dict[int, list[dict[str, Any]]] = {}
        for order in ORDERS:
            columns = order_map[order]
            if not isinstance(columns, Mapping):
                raise TypeError(f"{label} Z{order} must be a column mapping")
            value = columns.get(source_column, columns.get(str(source_column)))
            if value is None:
                # A sparse packet is allowed only for the requested source
                # column; all other columns are intentionally not inferred.
                value = []
            result[order] = _terms(
                c, value, label=f"{label}[{index}] Z{order} column {source_column}",
                default_id=f"{label}:{index}:Z{order}:column{source_column}",
                error_only=error_only,
            )
        rows.append(result)
    return radii, rows


def _chi_rows(c: Any, chi: Any) -> dict[int, list[dict[str, Any]]]:
    if chi is None:
        raise ValueError("chi Z0/Z1 is required")
    if isinstance(chi, Mapping) and ("coefficient" in chi or "powers" in chi):
        chi = {0: chi}
    elif isinstance(chi, Mapping):
        for key in ("ordinary_Z", "chi_jets", "jets", "rows"):
            if key in chi:
                chi = chi[key]
                break
    if isinstance(chi, Mapping):
        parsed: dict[int, Any] = {}
        for key, value in chi.items():
            text = key.upper().removeprefix("Z") if isinstance(key, str) else key
            try:
                parsed[int(text)] = value
            except (TypeError, ValueError) as exc:
                raise ValueError(f"invalid chi order {key!r}") from exc
    elif isinstance(chi, Sequence) and not isinstance(chi, (str, bytes, bytearray)):
        parsed = {index: value for index, value in enumerate(chi)}
    else:
        parsed = {0: chi}
    if 0 not in parsed:
        raise ValueError("chi Z0 is required")
    if 1 not in parsed:
        parsed[1] = {"coefficient": c.mpf(0),
                     "powers": (0, 0, 0, 0), "source_id": "chi:Z1"}
    if set(parsed) - {0, 1}:
        raise ValueError("chi supports only Z0 and Z1")
    result: dict[int, list[dict[str, Any]]] = {}
    for order in ORDERS:
        result[order] = _terms(c, parsed[order], label=f"chi Z{order}",
                                default_id=f"chi:Z{order}")
        for index, term in enumerate(result[order]):
            term.setdefault("chi_order", order)
            term.setdefault("chi_term_index", index)
            term.setdefault("chi_source_id", term.get("source_id", f"chi:Z{order}"))
    return result


def _nonnegative(c: Any, value: Any, label: str) -> Any:
    value = _interval(c, value)
    lo, hi = value._mpi_
    if c.make_mpf((hi, hi)) < 0:
        raise ValueError(f"{label} has negative upper bound")
    if c.make_mpf((lo, lo)) < 0:
        lo = c.mpf(0)._mpi_[0]
    return c.make_mpf((lo, hi))


def g0_cell_kernel_masses(c: Any, target_t: Any, source_left_t: Any,
                          source_right_t: Any, *, triangular: bool = False) -> dict[str, Any]:
    """Exact normalized ``G0`` F1 and d_x F1 masses on one source cell."""

    T = _interval(c, target_t)
    a = _interval(c, source_left_t)
    b = T if triangular else _interval(c, source_right_t)
    Tlo, Thi = _bounds(c, T)
    alo, ahi = _bounds(c, a)
    blo, bhi = _bounds(c, b)
    if alo < 0 or ahi < alo or blo < alo or bhi < blo:
        raise ValueError("invalid nonnegative G0 source cell")
    if not triangular and bhi > Thi:
        raise ValueError("G0 source cell extends beyond target")
    if triangular and not _same_interval(b, T):
        raise ValueError("triangular G0 endpoint must equal target")
    if triangular and alo == 0 and Tlo == 0:
        f1 = T * T / 2
        dxf1 = T / 2
    else:
        if Thi <= 0:
            raise ValueError("positive G0 target required")
        T2 = T * T
        a2 = a * a
        b2 = b * b
        if triangular:
            f1_hi = T2 / 2
            dxf1_hi = T / 2
        else:
            f1_hi = b2 - b2 * b2 / (2 * T2)
            dxf1_hi = b2 * b2 / (2 * T * T2)
        f1_lo = c.mpf(0) if alo == 0 else a2 - a2 * a2 / (2 * T2)
        dxf1_lo = c.mpf(0) if alo == 0 else a2 * a2 / (2 * T * T2)
        f1 = f1_hi - f1_lo
        dxf1 = dxf1_hi - dxf1_lo
    return {
        "F1": _nonnegative(c, f1, "G0 F1 mass"),
        "d_x_F1": _nonnegative(c, dxf1, "G0 d_x_F1 mass"),
        "target_t": T,
        "source_left_t": a,
        "source_right_t": b,
        "target_x": 2 * T,
        "source_left_x": 2 * a,
        "source_right_x": 2 * b,
        "operator": "G0",
        "triangular": bool(triangular),
        "axis_continuation": bool(triangular and alo == 0 and Tlo == 0),
        "exact_source_cell_mass": True,
    }


def _k1_primitive(c: Any, x: Any, s: Any) -> Any:
    s = _interval(c, s)
    if _bounds(c, s)[0] == 0:
        return c.mpf(0)
    return (x * x * s * s / 32
            - s**4 * c.log(x / s) / 16
            - s**4 / 64
            - s**6 / (96 * x * x))


def _k1_dx_primitive(c: Any, x: Any, s: Any) -> Any:
    s = _interval(c, s)
    if _bounds(c, s)[0] == 0:
        return c.mpf(0)
    return x * s * s / 16 - s**4 / (16 * x) + s**6 / (48 * x**3)


def d_cell_kernel_masses(c: Any, target_t: Any, source_left_t: Any,
                         source_right_t: Any, *, triangular: bool = False) -> dict[str, Any]:
    """Exact physical first-correction masses for ``D=(-2chi)TK0``."""

    T = _interval(c, target_t)
    a = _interval(c, source_left_t)
    b = T if triangular else _interval(c, source_right_t)
    target_x = 2 * T
    left_x = 2 * a
    right_x = 2 * b
    Tlo, Thi = _bounds(c, T)
    alo, ahi = _bounds(c, a)
    blo, bhi = _bounds(c, b)
    if alo < 0 or ahi < alo or blo < alo or bhi < blo:
        raise ValueError("invalid nonnegative D source cell")
    if not triangular and bhi > Thi:
        raise ValueError("D source cell extends beyond target")
    if triangular and not _same_interval(b, T):
        raise ValueError("triangular D endpoint must equal target")
    if triangular and alo == 0 and Tlo == 0:
        f1 = target_x**4 / 192
        dxf1 = target_x**3 / 48
    else:
        if Thi <= 0:
            raise ValueError("positive D target required")
        if triangular:
            f1_hi = target_x**4 / 192
            dxf1_hi = target_x**3 / 48
        else:
            f1_hi = _k1_primitive(c, target_x, right_x)
            dxf1_hi = _k1_dx_primitive(c, target_x, right_x)
        f1_lo = c.mpf(0) if alo == 0 else _k1_primitive(c, target_x, left_x)
        dxf1_lo = c.mpf(0) if alo == 0 else _k1_dx_primitive(c, target_x, left_x)
        f1 = f1_hi - f1_lo
        dxf1 = dxf1_hi - dxf1_lo
    return {
        "F1": _nonnegative(c, f1, "D F1 bare mass"),
        "d_x_F1": _nonnegative(c, dxf1, "D d_x_F1 bare mass"),
        "target_t": T,
        "source_left_t": a,
        "source_right_t": b,
        "target_x": target_x,
        "source_left_x": left_x,
        "source_right_x": right_x,
        "operator": "T K0 first correction; exterior factor -2 chi",
        "triangular": bool(triangular),
        "axis_continuation": bool(triangular and alo == 0 and Tlo == 0),
        "exact_source_cell_mass": True,
        "physical_units": True,
    }


def _combine(c: Any, left: Mapping[str, Any], right: Mapping[str, Any],
             factor: int, *, operation: str, metadata: Mapping[str, Any]) -> dict[str, Any]:
    powers = tuple(a + b for a, b in zip(left["powers"], right["powers"]))
    ids = _unique(_ids(left) + _ids(right))
    classes = _unique(_classes(left) + _classes(right))
    uncertain = _uncertainty(left) or _uncertainty(right)
    result = {
        "powers": powers,
        "hb_power": powers[0], "F0_power": powers[1],
        "Pstar_squared_power": powers[2], "Lambda_power": powers[3],
        "coefficient": factor * left["coefficient"] * right["coefficient"],
        "source_ids": ids, "error_classes": classes,
        "uncertainty": uncertain,
        "origin": operation,
        "origin_record": {
            "operation": operation, "source_ids": ids,
            "error_classes": classes, "uncertainty": uncertain,
            "symbolic_cancellation_claimed": False,
            **dict(metadata),
        },
        **dict(metadata),
    }
    if left.get("error_only") or right.get("error_only"):
        result["error_only"] = True
        result["symmetric_error"] = True
    result["source_id"] = ids[0] if ids else operation
    return result


def _build_p(c: Any, rows: Sequence[Mapping[int, Sequence[Mapping[str, Any]]]],
             chi: Mapping[int, Sequence[Mapping[str, Any]]], *, source_kind: str,
             source_index: int) -> dict[int, list[dict[str, Any]]]:
    """Form P=+2chi E51 rows while retaining both chi Leibniz levels."""

    values = rows[source_index]
    out: dict[int, list[dict[str, Any]]] = {0: [], 1: []}
    for chi_index, chi_term in enumerate(chi[0]):
        for value_index, value in enumerate(values[0]):
            out[0].append(_combine(
                c, chi_term, value, 2,
                operation="P=+2chi E51; chi*V",
                metadata={
                    "source_kind": source_kind, "source_cell": source_index,
                    "source_cell_index": source_index,
                    "chi_level": "P", "chi_order": 0,
                    "chi_term_index": chi_index, "V_order": 0,
                    "V_term_index": value_index,
                    "leibniz_multiplier": 1,
                    "P_formula": "+2 chi E51 V",
                },
            ))
    for chi_index, chi_term in enumerate(chi[1]):
        for value_index, value in enumerate(values[0]):
            out[1].append(_combine(
                c, chi_term, value, 2,
                operation="P=+2chi E51; chi_Z*V",
                metadata={
                    "source_kind": source_kind, "source_cell": source_index,
                    "source_cell_index": source_index,
                    "chi_level": "P", "chi_order": 1,
                    "chi_term_index": chi_index, "V_order": 0,
                    "V_term_index": value_index,
                    "leibniz_multiplier": 1,
                    "P_formula": "+2(chi_Z V + chi V_Z)",
                },
            ))
    for chi_index, chi_term in enumerate(chi[0]):
        for value_index, value in enumerate(values[1]):
            out[1].append(_combine(
                c, chi_term, value, 2,
                operation="P=+2chi E51; chi*V_Z",
                metadata={
                    "source_kind": source_kind, "source_cell": source_index,
                    "source_cell_index": source_index,
                    "chi_level": "P", "chi_order": 0,
                    "chi_term_index": chi_index, "V_order": 1,
                    "V_term_index": value_index,
                    "leibniz_multiplier": 1,
                    "P_formula": "+2(chi_Z V + chi V_Z)",
                },
            ))
    return out


def _apply_g0(c: Any, p: Mapping[int, Sequence[Mapping[str, Any]]], mass: Mapping[str, Any],
              *, source_index: int, target_index: int, branch: str,
              include_d_x_f1: bool) -> dict[str, dict[str, list[dict[str, Any]]]]:
    result = {"Z0": {field: [] for field in FIELDS},
              "Z1": {field: [] for field in FIELDS}}
    for order in ORDERS:
        for term in p[order]:
            for field in (FIELDS if include_d_x_f1 else ("F1",)):
                item = dict(term)
                item["coefficient"] = term["coefficient"] * mass[field]
                item.update({
                    "branch": branch, "operator": "G0",
                    "field": field, "target_cell": target_index,
                    "source_cell": source_index, "kernel_mass": mass[field],
                    "kernel_mass_field": field,
                    "radial_stage": "G0 after P",
                    "exact_source_cell_mass": True,
                })
                result[f"Z{order}"][field].append(item)
    if not include_d_x_f1:
        result["Z0"]["d_x_F1"] = []
        result["Z1"]["d_x_F1"] = []
    return result


def _apply_d(c: Any, p: Mapping[int, Sequence[Mapping[str, Any]]], mass: Mapping[str, Any],
             chi: Mapping[int, Sequence[Mapping[str, Any]]], *, source_index: int,
             target_index: int, branch: str, include_d_x_f1: bool) -> dict[str, dict[str, list[dict[str, Any]]]]:
    result = {"Z0": {field: [] for field in FIELDS},
              "Z1": {field: [] for field in FIELDS}}
    routes = {
        0: ((0, 0, 1, "D chi*P"),),
        1: ((1, 0, 1, "D chi_Z*P"), (0, 1, 1, "D chi*P_Z")),
    }
    for outer_order in ORDERS:
        for chi_order, p_order, multiplier, operation in routes[outer_order]:
            for chi_index, chi_term in enumerate(chi[chi_order]):
                for p_index, p_term in enumerate(p[p_order]):
                    ids = _unique(_ids(chi_term) + _ids(p_term))
                    classes = _unique(_classes(chi_term) + _classes(p_term))
                    uncertain = _uncertainty(chi_term) or _uncertainty(p_term)
                    powers = tuple(a + b for a, b in zip(chi_term["powers"], p_term["powers"]))
                    base = {
                        "powers": powers,
                        "hb_power": powers[0], "F0_power": powers[1],
                        "Pstar_squared_power": powers[2], "Lambda_power": powers[3],
                        "coefficient": -2 * multiplier * chi_term["coefficient"] * p_term["coefficient"],
                        "source_ids": ids, "error_classes": classes,
                        "uncertainty": uncertain,
                        "origin": operation,
                        "origin_record": {
                            "operation": operation, "source_ids": ids,
                            "error_classes": classes, "uncertainty": uncertain,
                            "symbolic_cancellation_claimed": False,
                        },
                        "branch": branch, "operator": "D=(-2chi)TK0",
                        "field": "F1", "target_cell": target_index,
                        "source_cell": source_index,
                        "chi_level": "D", "outer_chi_order": chi_order,
                        "outer_chi_term_index": chi_index,
                        "P_order": p_order, "P_term_index": p_index,
                        "outer_leibniz_multiplier": multiplier,
                        "chi_levels": [
                            {"level": "P", "recorded": True},
                            {"level": "D", "order": chi_order,
                             "leibniz_multiplier": multiplier},
                        ],
                        "radial_stage": "D after P",
                        "principal_kernel": "(-2chi) T K0",
                        "signed_source_term": True,
                    }
                    if chi_term.get("error_only") or p_term.get("error_only"):
                        base["error_only"] = True
                        base["symmetric_error"] = True
                    for field in (FIELDS if include_d_x_f1 else ("F1",)):
                        item = dict(base)
                        item["field"] = field
                        item["coefficient"] = base["coefficient"] * mass[field]
                        item["kernel_mass"] = mass[field]
                        item["kernel_mass_field"] = field
                        item["exact_source_cell_mass"] = True
                        result[f"Z{outer_order}"][field].append(item)
    if not include_d_x_f1:
        result["Z0"]["d_x_F1"] = []
        result["Z1"]["d_x_F1"] = []
    return result


def _branch(c: Any, radii: Sequence[Any], source_rows: Sequence[Mapping[int, Sequence[Mapping[str, Any]]]],
            chi: Mapping[int, Sequence[Mapping[str, Any]]], *, source_kind: str,
            operator: str, name: str, include_d_x_f1: bool) -> dict[str, Any]:
    p_rows = [_build_p(c, source_rows, chi, source_kind=source_kind,
                       source_index=index) for index in range(CELL_COUNT)]
    target_cells: list[dict[str, Any]] = []
    for target_index in range(CELL_COUNT):
        target_t = radii[target_index]
        aggregate = {f"Z{order}": {field: [] for field in FIELDS} for order in ORDERS}
        by_source: dict[str, Any] = {}
        for source_index in range(target_index + 1):
            source_t = radii[source_index]
            triangular = source_index == target_index
            left = _endpoint(c, source_t, False)
            right = target_t if triangular else _endpoint(c, source_t, True)
            mass = (g0_cell_kernel_masses(c, target_t, left, right,
                                           triangular=triangular)
                    if operator == "G0" else
                    d_cell_kernel_masses(c, target_t, left, right,
                                         triangular=triangular))
            applied = (_apply_g0(c, p_rows[source_index], mass,
                                 source_index=source_index, target_index=target_index,
                                 branch=name, include_d_x_f1=include_d_x_f1)
                       if operator == "G0" else
                       _apply_d(c, p_rows[source_index], mass, chi,
                                source_index=source_index, target_index=target_index,
                                branch=name, include_d_x_f1=include_d_x_f1))
            by_source[str(source_index)] = {
                "source_index": source_index, "target_index": target_index,
                "source_t": source_t, "target_t": target_t,
                "triangular": triangular,
                "kernel_masses": mass, "P_q5": p_rows[source_index],
                "fields": applied,
            }
            for order in ORDERS:
                for field in FIELDS:
                    aggregate[f"Z{order}"][field].extend(applied[f"Z{order}"][field])
        target_cells.append({
            "index": target_index, "target_index": target_index,
            "t": target_t, "x": 2 * target_t, "rho": 4 * target_t * target_t,
            "source_indices": list(range(target_index + 1)),
            "triangular_diagonal_source": target_index,
            "axis_endpoint_continuation": target_index == 0,
            "source_kind": source_kind, "operator": operator, "branch": name,
            "P_q5_by_source_cell": {str(i): p_rows[i] for i in range(target_index + 1)},
            "fields": aggregate, "by_source": by_source,
        })

    exit_t = c.mpf(1)
    exit_by_source: dict[str, Any] = {}
    exit_aggregate = {f"Z{order}": {field: [] for field in FIELDS} for order in ORDERS}
    exit_masses: dict[str, Any] = {}
    for source_index in range(CELL_COUNT):
        source_t = radii[source_index]
        left = _endpoint(c, source_t, False)
        right = _endpoint(c, source_t, True)
        mass = (g0_cell_kernel_masses(c, exit_t, left, right)
                if operator == "G0" else
                d_cell_kernel_masses(c, exit_t, left, right))
        applied = (_apply_g0(c, p_rows[source_index], mass,
                             source_index=source_index, target_index=-1,
                             branch=name, include_d_x_f1=include_d_x_f1)
                   if operator == "G0" else
                   _apply_d(c, p_rows[source_index], mass, chi,
                            source_index=source_index, target_index=-1,
                            branch=name, include_d_x_f1=include_d_x_f1))
        exit_masses[str(source_index)] = mass
        exit_by_source[str(source_index)] = {
            "source_index": source_index, "source_t": source_t,
            "target_t": exit_t, "target_x": c.mpf(2),
            "kernel_masses": mass, "P_q5": p_rows[source_index],
            "fields": applied,
        }
        for order in ORDERS:
            for field in FIELDS:
                exit_aggregate[f"Z{order}"][field].extend(applied[f"Z{order}"][field])
    return {
        "branch": name, "source_kind": source_kind, "operator": operator,
        "formula": ("G0 P " if operator == "G0" else "D P ") + source_kind,
        "kernel_units": "normalized t=x/2" if operator == "G0" else "physical x",
        "P_q5_by_source_cell": {str(i): p_rows[i] for i in range(CELL_COUNT)},
        "target_cells": target_cells,
        "exit_target": {"x": c.mpf(2), "t": exit_t, "rho": c.mpf(4),
                         "fixed": True, "independent": True,
                         "not_last_target_cell_proxy": True},
        "exit_by_source_cell": exit_by_source,
        "exit_x2_by_source_cell": exit_by_source,
        "exit_aggregate": exit_aggregate,
        "exit_kernel_masses_by_source_cell": exit_masses,
        "fields": list(FIELDS if include_d_x_f1 else ("F1",)),
        "ordinary_Z_orders": [0, 1],
        "source_ids_retained": True, "formal_powers_retained": True,
        "signed_coefficients_retained": True,
        "uncertainty_separate": True,
        "source_cells_collapsed": False,
        "complete_feedback": False,
    }


def build_compensation(
    c: Any,
    W0_cells: Any = None,
    delta_W0_cells: Any = None,
    chi: Any = None,
    *,
    w0_cells: Any = None,
    signed_delta_cells: Any = None,
    chi_jets: Any = None,
    source_column: int = 1,
    include_d_x_F1: bool = True,
) -> dict[str, Any]:
    """Build the four finite ``E51`` compensation branches.

    ``W0_cells`` and ``delta_W0_cells`` each contain eight canonical cells
    and ordinary F1 rows Z0/Z1.  ``chi`` supplies chi Z0/Z1.  A direct
    ``F1`` wrapper is preferred; a sparse ``W0_ordinary_Z`` wrapper selects
    ``source_column`` (default 1).  No raw C0/B1 or W0 reconstruction is
    accepted or required.
    """

    if not hasattr(c, "mpf") or not hasattr(c, "make_mpf"):
        raise TypeError("c must be a directed MPIntervalContext")
    if W0_cells is None:
        W0_cells = w0_cells
    if delta_W0_cells is None:
        delta_W0_cells = signed_delta_cells
    if W0_cells is None or delta_W0_cells is None:
        raise ValueError("both caller-owned W0 and signed Delta W0 F1 cells are required")
    if chi is None:
        chi = chi_jets
    chi_rows = _chi_rows(c, chi)
    w0_radii, w0_rows = _normalise_source_cells(
        c, W0_cells, label="W0 F1 source", source_column=source_column,
    )
    delta_radii, delta_rows = _normalise_source_cells(
        c, delta_W0_cells, label="Delta W0 F1 source", source_column=source_column,
    )
    for index in range(CELL_COUNT):
        if not _same_interval(w0_radii[index], delta_radii[index]):
            raise ValueError(f"W0 and Delta W0 cell {index} coordinates differ")
    branches = {
        "G0_P_W0": _branch(c, w0_radii, w0_rows, chi_rows,
                            source_kind="W0", operator="G0",
                            name="G0_P_W0", include_d_x_f1=include_d_x_F1),
        "G0_P_DeltaW0": _branch(c, w0_radii, delta_rows, chi_rows,
                                source_kind="DeltaW0", operator="G0",
                                name="G0_P_DeltaW0", include_d_x_f1=include_d_x_F1),
        "D_P_W0": _branch(c, w0_radii, w0_rows, chi_rows,
                           source_kind="W0", operator="D",
                           name="D_P_W0", include_d_x_f1=include_d_x_F1),
        "D_P_DeltaW0": _branch(c, w0_radii, delta_rows, chi_rows,
                                source_kind="DeltaW0", operator="D",
                                name="D_P_DeltaW0", include_d_x_f1=include_d_x_F1),
    }
    return {
        "schema": "n1_kchi_physical_compensation_v2",
        "branches": branches,
        "target_cells_by_branch": {
            name: branch["target_cells"] for name, branch in branches.items()
        },
        "exit_x2_by_branch": {
            name: branch["exit_by_source_cell"] for name, branch in branches.items()
        },
        "G0_P_W0_exit_x2_by_source_cell": branches["G0_P_W0"]["exit_by_source_cell"],
        "G0_P_DeltaW0_exit_x2_by_source_cell": branches["G0_P_DeltaW0"]["exit_by_source_cell"],
        "D_P_W0_exit_x2_by_source_cell": branches["D_P_W0"]["exit_by_source_cell"],
        "D_P_DeltaW0_exit_x2_by_source_cell": branches["D_P_DeltaW0"]["exit_by_source_cell"],
        "chi_jets": chi_rows,
        "source_power_basis": list(SOURCE_POWER_BASIS),
        "fields": list(FIELDS if include_d_x_F1 else ("F1",)),
        "ordinary_Z_orders": [0, 1],
        "source_column": source_column,
        "P_formula": "P=+2 chi E51; P_Z1=2(chi_Z V+chi V_Z)",
        "D_formula": "D=(-2 chi) T K0 applied after P",
        "chi_leibniz_levels_retained": ["P", "D"],
        "uncertainty_separate": True,
        "source_ids_retained": True,
        "formal_powers_retained": True,
        "signed_coefficients_retained": True,
        "axis_and_diagonal_target_cells": True,
        "fixed_x2_exit_included": True,
        "fixed_x2_bare_D_full_interval_masses": d_cell_kernel_masses(
            c, c.mpf(1), c.mpf(0), c.mpf(1), triangular=True,
        ),
        "raw_C0_consumed": False,
        "raw_B1_consumed": False,
        "W0_reconstructed": False,
        "complete_feedback": False,
        "complete_Kchi_W1_inlet": False,
        "finite_compensation_only": True,
    }


build = build_compensation
physical_compensation = build_compensation
