"""Independent v2 signed-Delta-W0 g5 endpoint.

This file deliberately does not import or call the frozen v1 source.  It is
the narrow endpoint leaf only: caller-owned signed row-5 source cells are
multiplied by ``P = -2 chi`` and integrated with the exact physical
first-correction kernels at fixed ``x=2``.  The order parser accepts integer
keys, ``0``/``1`` strings, and ``Z0``/``Z1`` spellings for chi rows.

The result is not a complete Kchi W1 inlet.  Raw C0/B1, residual
conjugation compensation, and ancestor producers are intentionally absent.
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
    "delta_kernel_masses", "build_signed_delta_g5_exit", "build",
    "build_core_exit", "build_endpoint", "signed_delta_g5_exit",
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


def _order_key(key: Any, label: str) -> int:
    """Normalize int, ``0``/``1``, and ``Z0``/``Z1`` order spellings."""

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


def _is_interval(value: Any) -> bool:
    return hasattr(value, "_mpi_") and hasattr(value, "a") and hasattr(value, "b")


def _interval(c: Any, value: Any) -> Any:
    if _is_interval(value):
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


def _items(value: Any) -> list[Any]:
    if value is None:
        return []
    if isinstance(value, (str, bytes, bytearray)):
        return [value]
    if isinstance(value, Sequence):
        return list(value)
    return [value]


def _unique(values: Sequence[Any]) -> list[Any]:
    out: list[Any] = []
    seen: set[Any] = set()
    for value in values:
        try:
            key = (type(value).__name__, value)
            hash(key)
        except TypeError:
            key = (type(value).__name__, repr(value))
        if key not in seen:
            seen.add(key)
            out.append(value)
    return out


def _nested(raw: Any) -> Mapping[str, Any]:
    if not isinstance(raw, Mapping):
        return {}
    value = raw.get("origin")
    return value if isinstance(value, Mapping) else {}


def _ids(raw: Any) -> list[Any]:
    if not isinstance(raw, Mapping):
        return []
    values = list(_items(raw.get("source_ids", raw.get("source_id"))))
    nested = _nested(raw)
    values.extend(_items(nested.get("source_ids", nested.get("source_id"))))
    return _unique(values)


def _classes(raw: Any) -> list[Any]:
    if not isinstance(raw, Mapping):
        return []
    values = list(_items(raw.get("error_classes", raw.get("error_class"))))
    nested = _nested(raw)
    values.extend(_items(nested.get("error_classes", nested.get("error_class"))))
    return _unique(values)


def _uncertain(raw: Any) -> bool:
    if not isinstance(raw, Mapping):
        return False
    nested = _nested(raw)
    return bool(raw.get("uncertainty", raw.get("is_uncertainty", False))
                or raw.get("symmetric_error", False)
                or nested.get("uncertainty", nested.get("is_uncertainty", False)))


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
            or any(isinstance(item, bool) or not isinstance(item, int) for item in value)):
        raise ValueError(f"{label} must carry four integer formal powers")
    return tuple(int(item) for item in value)


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
    uncertain = _uncertain(raw) or error_only
    result.update({
        "powers": powers,
        "hb_power": powers[0], "F0_power": powers[1],
        "Pstar_squared_power": powers[2], "Lambda_power": powers[3],
        "coefficient": coefficient, "source_ids": _unique(ids),
        "error_classes": list(classes), "uncertainty": uncertain,
        "origin": dict(_nested(raw)),
    })
    result["origin"].update({
        "source_ids": list(result["source_ids"]),
        "error_classes": list(classes), "uncertainty": uncertain,
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


def _packet_wrapper(cell: Mapping[str, Any]) -> Any:
    for key in ("W0_ordinary_Z", "ordinary_Z", "W0", "source_rows",
                "source_jets", "delta_W0", "delta_w0", "R", "error"):
        if key in cell:
            return cell[key]
    return cell


def _source_cells(c: Any, cells: Any, *, label: str,
                  source_column: int = 5,
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
            order_map = _order_map(direct, (0, 1), f"{label} F1")
            rows.append({
                order: _terms(c, order_map[order],
                              label=f"{label}[{index}] Z{order}",
                              default_id=f"{label}:{index}:Z{order}",
                              error_only=error_only)
                for order in ORDERS
            })
            continue
        packet = _packet_wrapper(cell)
        order_map = _order_map(packet, (0, 1, 2), f"{label} W0")
        result: dict[int, list[dict[str, Any]]] = {}
        for order in ORDERS:
            columns = order_map[order]
            if not isinstance(columns, Mapping):
                raise TypeError(f"{label} Z{order} must be a column mapping")
            value = columns.get(source_column, columns.get(str(source_column), []))
            result[order] = _terms(
                c, value, label=f"{label}[{index}] Z{order} column {source_column}",
                default_id=f"{label}:{index}:Z{order}:column{source_column}",
                error_only=error_only,
            )
        rows.append(result)
    return radii, rows


def _chi_rows(c: Any, chi: Any) -> dict[int, list[dict[str, Any]]]:
    """Parse chi rows; importantly, ``Z0`` is normalized before int()."""

    if chi is None:
        raise ValueError("chi Z0/Z1 is required")
    if isinstance(chi, Mapping) and ("coefficient" in chi or "powers" in chi):
        chi = {0: chi}
    elif isinstance(chi, Mapping):
        for wrapper in ("ordinary_Z", "chi_jets", "jets", "rows"):
            if wrapper in chi:
                chi = chi[wrapper]
                break
    if isinstance(chi, Mapping):
        parsed: dict[int, Any] = {}
        for key, value in chi.items():
            order = _order_key(key, "chi")
            if order not in ORDERS:
                raise ValueError(f"chi supports only Z0 and Z1, got {key!r}")
            if order in parsed:
                raise ValueError(f"duplicate chi order {order}")
            parsed[order] = value
    elif isinstance(chi, Sequence) and not isinstance(chi, (str, bytes, bytearray)):
        parsed = {index: value for index, value in enumerate(chi)}
    else:
        parsed = {0: chi}
    if 0 not in parsed:
        raise ValueError("chi Z0 is required")
    if 1 not in parsed:
        parsed[1] = {"coefficient": c.mpf(0),
                     "powers": (0, 0, 0, 0), "source_id": "chi:Z1"}
    result: dict[int, list[dict[str, Any]]] = {}
    for order in ORDERS:
        result[order] = _terms(c, parsed[order], label=f"chi Z{order}",
                                default_id=f"chi:Z{order}")
        for index, term in enumerate(result[order]):
            term["chi_order"] = order
            term["chi_term_index"] = index
            term["chi_source_id"] = term.get("source_id", f"chi:Z{order}")
    return result


def _nonnegative(c: Any, value: Any, label: str) -> Any:
    value = _interval(c, value)
    lo, hi = value._mpi_
    if c.make_mpf((hi, hi)) < 0:
        raise ValueError(f"{label} has negative upper bound")
    if c.make_mpf((lo, lo)) < 0:
        lo = c.mpf(0)._mpi_[0]
    return c.make_mpf((lo, hi))


def _primitive(c: Any, x: Any, s: Any) -> Any:
    s = _interval(c, s)
    if _bounds(c, s)[0] == 0:
        return c.mpf(0)
    return (x * x * s * s / 32 - s**4 * c.log(x / s) / 16
            - s**4 / 64 - s**6 / (96 * x * x))


def _dx_primitive(c: Any, x: Any, s: Any) -> Any:
    s = _interval(c, s)
    if _bounds(c, s)[0] == 0:
        return c.mpf(0)
    return x * s * s / 16 - s**4 / (16 * x) + s**6 / (48 * x**3)


def delta_kernel_masses(c: Any, target_t: Any, source_left_t: Any,
                        source_right_t: Any, *, triangular: bool = False) -> dict[str, Any]:
    """Bare physical K1 masses before the exterior ``-2 chi`` factor."""

    T = _interval(c, target_t)
    a = _interval(c, source_left_t)
    b = T if triangular else _interval(c, source_right_t)
    Tlo, Thi = _bounds(c, T)
    alo, ahi = _bounds(c, a)
    blo, bhi = _bounds(c, b)
    if alo < 0 or ahi < alo or blo < alo or bhi < blo:
        raise ValueError("invalid source cell")
    if not triangular and bhi > Thi:
        raise ValueError("source extends past target")
    if triangular and not _same_interval(b, T):
        raise ValueError("triangular source endpoint mismatch")
    x = 2 * T
    ax = 2 * a
    bx = 2 * b
    if triangular and alo == 0 and Tlo == 0:
        f1 = x**4 / 192
        dxf1 = x**3 / 48
    else:
        if Thi <= 0:
            raise ValueError("positive target required")
        f_hi = x**4 / 192 if triangular else _primitive(c, x, bx)
        dx_hi = x**3 / 48 if triangular else _dx_primitive(c, x, bx)
        f_lo = c.mpf(0) if alo == 0 else _primitive(c, x, ax)
        dx_lo = c.mpf(0) if alo == 0 else _dx_primitive(c, x, ax)
        f1 = f_hi - f_lo
        dxf1 = dx_hi - dx_lo
    return {
        "F1": _nonnegative(c, f1, "K1 F1 mass"),
        "d_x_F1": _nonnegative(c, dxf1, "K1 d_x_F1 mass"),
        "target_t": T, "source_left_t": a, "source_right_t": b,
        "target_x": x, "source_left_x": ax, "source_right_x": bx,
        "operator": "T K0 first correction", "physical_units": True,
        "triangular": bool(triangular),
        "axis_continuation": bool(triangular and alo == 0 and Tlo == 0),
        "exact_source_cell_mass": True,
    }


def _combine(c: Any, left: Mapping[str, Any], right: Mapping[str, Any],
             factor: int, *, operation: str, metadata: Mapping[str, Any]) -> dict[str, Any]:
    powers = tuple(a + b for a, b in zip(left["powers"], right["powers"]))
    ids = _unique(_ids(left) + _ids(right))
    classes = _unique(_classes(left) + _classes(right))
    uncertain = _uncertain(left) or _uncertain(right)
    result = {
        "powers": powers, "hb_power": powers[0], "F0_power": powers[1],
        "Pstar_squared_power": powers[2], "Lambda_power": powers[3],
        "coefficient": factor * left["coefficient"] * right["coefficient"],
        "source_ids": ids, "error_classes": classes,
        "uncertainty": uncertain, "origin": operation,
        "origin_record": {"operation": operation, "source_ids": ids,
                           "error_classes": classes, "uncertainty": uncertain,
                           "symbolic_cancellation_claimed": False,
                           **dict(metadata)},
        **dict(metadata),
    }
    if left.get("error_only") or right.get("error_only"):
        result["error_only"] = True
        result["symmetric_error"] = True
    result["source_id"] = ids[0] if ids else operation
    return result


def _apply_d(c: Any, rows: Mapping[int, Sequence[Mapping[str, Any]]],
             chi: Mapping[int, Sequence[Mapping[str, Any]]], mass: Mapping[str, Any],
             *, source_index: int, branch: str,
             include_d_x_f1: bool) -> dict[str, dict[str, list[dict[str, Any]]]]:
    out = {f"Z{order}": {field: [] for field in FIELDS} for order in ORDERS}
    routes = {0: ((0, 0, "chi*V"),),
              1: ((1, 0, "chi_Z*V"), (0, 1, "chi*V_Z"))}
    for outer_order in ORDERS:
        for chi_order, source_order, operation in routes[outer_order]:
            for chi_index, chi_term in enumerate(chi[chi_order]):
                for source_index_term, source_term in enumerate(rows[source_order]):
                    ids = _unique(_ids(chi_term) + _ids(source_term))
                    classes = _unique(_classes(chi_term) + _classes(source_term))
                    uncertain = _uncertain(chi_term) or _uncertain(source_term)
                    powers = tuple(a + b for a, b in zip(chi_term["powers"], source_term["powers"]))
                    base = {
                        "powers": powers, "hb_power": powers[0],
                        "F0_power": powers[1], "Pstar_squared_power": powers[2],
                        "Lambda_power": powers[3],
                        "coefficient": -2 * chi_term["coefficient"] * source_term["coefficient"],
                        "source_ids": ids, "error_classes": classes,
                        "uncertainty": uncertain, "origin": operation,
                        "origin_record": {"operation": operation,
                                           "source_ids": ids,
                                           "error_classes": classes,
                                           "uncertainty": uncertain,
                                           "symbolic_cancellation_claimed": False},
                        "branch": branch, "source_cell": source_index,
                        "chi_order": chi_order, "chi_term_index": chi_index,
                        "source_order": source_order,
                        "source_term_index": source_index_term,
                        "outer_leibniz_multiplier": 1,
                        "principal_kernel": "(-2chi) T K0",
                        "chi_level": "D", "signed_source_term": True,
                    }
                    if chi_term.get("error_only") or source_term.get("error_only"):
                        base["error_only"] = True
                        base["symmetric_error"] = True
                    fields = FIELDS if include_d_x_f1 else ("F1",)
                    for field in fields:
                        item = dict(base)
                        item["field"] = field
                        item["coefficient"] = base["coefficient"] * mass[field]
                        item["kernel_mass"] = mass[field]
                        item["kernel_mass_field"] = field
                        item["exact_source_cell_mass"] = True
                        out[f"Z{outer_order}"][field].append(item)
    if not include_d_x_f1:
        out["Z0"]["d_x_F1"] = []
        out["Z1"]["d_x_F1"] = []
    return out


def _evaluate(c: Any, radii: Sequence[Any], rows: Sequence[Mapping[int, Sequence[Mapping[str, Any]]]],
              chi: Mapping[int, Sequence[Mapping[str, Any]]], *, branch: str,
              include_d_x_f1: bool, error_only: bool = False) -> dict[str, Any]:
    masses: dict[str, Any] = {}
    exits: dict[str, Any] = {}
    aggregate = {f"Z{order}": {field: [] for field in FIELDS} for order in ORDERS}
    for index in range(CELL_COUNT):
        source_t = radii[index]
        mass = delta_kernel_masses(
            c, c.mpf(1), _endpoint(c, source_t, False),
            _endpoint(c, source_t, True), triangular=False,
        )
        masses[str(index)] = mass
        fields = _apply_d(
            c, rows[index], chi, mass, source_index=index,
            branch=branch, include_d_x_f1=include_d_x_f1,
        )
        exits[str(index)] = {
            "source_index": index, "source_t": source_t,
            "target_t": c.mpf(1), "target_x": c.mpf(2),
            "kernel_masses": mass, "fields": fields,
        }
        for order in ORDERS:
            for field in FIELDS:
                aggregate[f"Z{order}"][field].extend(fields[f"Z{order}"][field])
    return {
        "branch": branch, "exit_by_source_cell": exits,
        "exit_x2_by_source_cell": exits, "exit_aggregate": aggregate,
        "kernel_masses_by_source_cell": masses,
        "target_cells": [], "fields": list(FIELDS if include_d_x_f1 else ("F1",)),
        "error_only": error_only, "symmetric_error": error_only,
        "source_ids_retained": True, "formal_powers_retained": True,
        "signed_coefficients_retained": True, "complete_feedback": False,
    }


def build_signed_delta_g5_exit(
    c: Any,
    g5_cells: Any = None,
    chi: Any = None,
    *,
    signed_delta_g5_cells: Any = None,
    delta_w0_cells: Any = None,
    symmetric_error_cells: Any = None,
    symmetric_R_cells: Any = None,
    R_cells: Any = None,
    chi_jets: Any = None,
    source_column: int = 5,
    include_d_x_f1: bool = True,
) -> dict[str, Any]:
    """Evaluate the independent signed g5 endpoint at fixed ``x=2``."""

    if g5_cells is None:
        g5_cells = signed_delta_g5_cells if signed_delta_g5_cells is not None else delta_w0_cells
    if g5_cells is None:
        raise ValueError("caller-owned signed Delta W0 g5 cells are required")
    if chi is None:
        chi = chi_jets
    chi_rows = _chi_rows(c, chi)
    radii, rows = _source_cells(c, g5_cells, label="signed Delta W0 g5",
                                source_column=source_column)
    branch = _evaluate(c, radii, rows, chi_rows,
                       branch="signed_DeltaW0_g5_v2",
                       include_d_x_f1=include_d_x_f1)
    result: dict[str, Any] = {
        "schema": "n1_kchi_core_exit_v2",
        "signed_delta_g5_exit": branch,
        "branch": branch,
        "exit_by_source_cell": branch["exit_by_source_cell"],
        "exit_x2_by_source_cell": branch["exit_by_source_cell"],
        "aggregate": branch["exit_aggregate"],
        "chi_jets": chi_rows,
        "exit_target": {"x": c.mpf(2), "t": c.mpf(1), "rho": c.mpf(4),
                         "fixed": True, "independent": True,
                         "not_last_target_cell_proxy": True},
        "axis_scale_identity": {
            "physical_x": c.mpf(2), "normalized_t": c.mpf(1),
            "rho": c.mpf(4), "identity": "t=x/2; rho=4*t^2=x^2",
            "proven": True,
        },
        "full_physical_interval_kernel_masses": delta_kernel_masses(
            c, c.mpf(1), c.mpf(0), c.mpf(1), triangular=True,
        ),
        "full_interval_exterior_masses": {"F1": "-chi/6", "d_x_F1": "-chi/3"},
        "kernel_mass_definition": "exact physical T K0 first-correction source-cell mass",
        "source_power_basis": list(SOURCE_POWER_BASIS),
        "fields": list(FIELDS if include_d_x_f1 else ("F1",)),
        "ordinary_Z_orders": [0, 1], "source_column": source_column,
        "source_ids_retained": True, "formal_powers_retained": True,
        "signed_coefficients_retained": True, "green_kernel_masses_exact": True,
        "delta_kernel_masses_exact": True, "negative_kernel_K1": True,
        "saved_W0_reused": True, "W0_reconstructed": False,
        "ancestor_helpers_recomputed": False,
        "baseline_E00_recomputed": False,
        "normalization_conjugation_compensation_included": False,
        "normalization_conjugation_compensation_status": "unproved; +2 chi E51 omitted",
        "complete_Kchi_W1_inlet": False, "complete_feedback": False,
        "finite_compensation_only": True,
    }
    error_input = symmetric_error_cells
    if error_input is None:
        error_input = symmetric_R_cells
    if error_input is None:
        error_input = R_cells
    if error_input is not None:
        error_radii, error_rows = _source_cells(
            c, error_input, label="symmetric R g5", source_column=source_column,
            error_only=True,
        )
        for index in range(CELL_COUNT):
            if not _same_interval(error_radii[index], radii[index]):
                raise ValueError(f"symmetric R cell {index} coordinate differs")
        error_branch = _evaluate(
            c, error_radii, error_rows, chi_rows,
            branch="signed_DeltaW0_g5_symmetric_R_v2",
            include_d_x_f1=include_d_x_f1, error_only=True,
        )
        result["symmetric_R_error"] = {
            "branch": error_branch, "error_only": True,
            "symmetric_error": True,
        }
        result["error_branch_separate"] = True
    else:
        result["symmetric_R_error"] = None
        result["error_branch_separate"] = True
    return result


build = build_signed_delta_g5_exit
build_core_exit = build_signed_delta_g5_exit
build_endpoint = build_signed_delta_g5_exit
signed_delta_g5_exit = build_signed_delta_g5_exit
