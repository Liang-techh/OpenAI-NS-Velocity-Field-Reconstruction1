"""Source-bound physical ``m15`` cell assembly for the E00 packet.

This leaf is intentionally a small consumer boundary.  The caller supplies
directed interval enclosures for the leading ``U``/``F`` fields and for the
already-built K0, J1, and J2 field sectors.  The sectors are indexed by the
whole q-cell that they enclose; no point in a cell is selected here and no
ancestor leaf is imported or called.

For one cell the two physical measures used below are

``+2 Lambda**-1 d rho = hb * Lambda**-1 * sc/2 * rho dq`` and
``-2 Lambda**-2 rho d rho = -hb * Lambda**-2 * sc/2 * rho**2 dq``.

The formal ``hb`` power is kept in the returned power tuple.  Only the
caller-certified cutoff mass and directed ``rho``/``rho**2`` intervals enter
the coefficient.  J1 and J2 are therefore consumed as fields, never applied
again as operators.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from typing import Any

SCHEMA = "n1_e00_physical_m15_cells_v2"
SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
ORDINARY_Z_ORDERS = (0, 1)
PARTS = ("K0", "J1", "J2")
CHANNELS = ("axial", "swirl")
ALL_CELLS = (0, 1, 2, 3)
ACTIVE_CELLS = (0, 1, 2)
ETA_CELL_MASS = ((1, 4), (1, 4), (1, 8), (0, 1))
PHYSICAL_EXTERNAL_POWERS = {"axial": (1, 0, 0, -1), "swirl": (1, 0, 0, -2)}
PART_CHANNELS = {"K0": ("axial", "swirl"), "J1": ("swirl",), "J2": ("swirl",)}

_PART_ALIASES = {"k0": "K0", "k0_u_f": "K0", "j1": "J1", "j2": "J2"}
_FIELD_ALIASES = {
    "K0": {
        "axial": ("U", "u", "Uz1", "uz1", "axial"),
        "swirl": ("F", "f", "F1", "f1", "swirl"),
    },
    "J1": {"swirl": ("F", "f", "F1", "f1", "swirl")},
    "J2": {"swirl": ("F", "f", "F1", "f1", "swirl")},
}

__all__ = [
    "SCHEMA",
    "SOURCE_POWER_BASIS",
    "ORDINARY_Z_ORDERS",
    "PARTS",
    "CHANNELS",
    "ALL_CELLS",
    "ACTIVE_CELLS",
    "ETA_CELL_MASS",
    "build",
    "build_physical_m15",
    "evaluate",
]


def _is_interval(value: Any) -> bool:
    return hasattr(value, "_mpi_") and hasattr(value, "a") and hasattr(value, "b")


def _raw_endpoint(value: Any, upper: bool) -> Any:
    raw = getattr(value, "_mpi_", None)
    return raw[1 if upper else 0] if raw is not None else value


def _make_interval(c: Any, lower: Any, upper: Any) -> Any:
    """Construct a directed interval without taking a midpoint."""

    return c.make_mpf((_raw_endpoint(lower, False), _raw_endpoint(upper, True)))


def _scalar_bounds(value: Any) -> tuple[Any, Any]:
    raw = getattr(value, "_mpi_", None)
    if raw is None:
        return value, value
    # Keep singleton endpoints in the caller's directed interval context.
    return value.a, value.b


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
    result: list[Any] = []
    seen: set[Any] = set()
    for value in values:
        key = _freeze(value)
        if key not in seen:
            seen.add(key)
            result.append(value)
    return result


def _origin_info(value: Any, *, _seen: set[int] | None = None) -> dict[str, Any]:
    """Read provenance without dropping nested source/error annotations."""

    if _seen is None:
        _seen = set()
    marker = id(value)
    if marker in _seen:
        return {
            "source_ids": [],
            "error_classes": [],
            "uncertainty": False,
            "error_only": False,
            "symmetric_error": False,
        }
    if isinstance(value, Mapping):
        _seen.add(marker)
        ids = list(_as_items(value.get("source_ids", value.get("sources"))))
        ids.extend(_as_items(value.get("source_id")))
        classes = list(_as_items(value.get("error_classes", value.get("error_class"))))
        uncertain = bool(value.get("uncertainty", value.get("is_uncertainty", False)))
        error_only = bool(value.get("error_only", value.get("is_error_only", False)))
        symmetric = bool(value.get("symmetric_error", value.get("is_symmetric_error", False)))
        operation = str(value.get("operation", "")).lower()
        kind = str(value.get("error_kind", value.get("kind", ""))).lower()
        uncertain = uncertain or any(token in operation or token in kind for token in ("error", "uncertainty", "symmetric"))
        error_only = error_only or "error" in operation or "error" in kind
        symmetric = symmetric or "symmetric" in operation or "symmetric" in kind
        nested: list[dict[str, Any]] = []
        for key in ("origin", "metadata", "parents"):
            child = value.get(key)
            if child is None or child is value:
                continue
            if isinstance(child, Sequence) and not isinstance(child, (str, bytes, bytearray)):
                nested.extend(_origin_info(item, _seen=_seen) for item in child)
            else:
                nested.append(_origin_info(child, _seen=_seen))
        for info in nested:
            ids.extend(info["source_ids"])
            classes.extend(info["error_classes"])
            uncertain = uncertain or info["uncertainty"]
            error_only = error_only or info["error_only"]
            symmetric = symmetric or info["symmetric_error"]
        return {
            "source_ids": _unique(ids),
            "error_classes": _unique(classes),
            "uncertainty": bool(uncertain),
            "error_only": bool(error_only),
            "symmetric_error": bool(symmetric),
        }
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        infos = [_origin_info(item, _seen=_seen) for item in value]
        ids: list[Any] = []
        classes: list[Any] = []
        uncertain = error_only = symmetric = False
        for info in infos:
            ids.extend(info["source_ids"])
            classes.extend(info["error_classes"])
            uncertain = uncertain or info["uncertainty"]
            error_only = error_only or info["error_only"]
            symmetric = symmetric or info["symmetric_error"]
        return {
            "source_ids": _unique(ids),
            "error_classes": _unique(classes),
            "uncertainty": uncertain,
            "error_only": error_only,
            "symmetric_error": symmetric,
        }
    return {
        "source_ids": [],
        "error_classes": [],
        "uncertainty": False,
        "error_only": False,
        "symmetric_error": False,
    }


def _source_paths(value: Any, *, _seen: set[int] | None = None) -> list[Any]:
    """Collect explicit source paths while leaving their structure untouched."""

    if _seen is None:
        _seen = set()
    marker = id(value)
    if marker in _seen:
        return []
    if isinstance(value, Mapping):
        _seen.add(marker)
        paths: list[Any] = []
        for key in ("source_paths", "source_path", "path", "paths"):
            if key in value and value[key] is not None:
                raw = value[key]
                if key in ("source_paths", "paths") and isinstance(raw, Sequence) and not isinstance(raw, (str, bytes, bytearray)):
                    paths.extend(raw)
                else:
                    paths.append(raw)
        for key in ("origin", "metadata", "parents"):
            child = value.get(key)
            if child is None or child is value:
                continue
            if isinstance(child, Sequence) and not isinstance(child, (str, bytes, bytearray)):
                for item in child:
                    paths.extend(_source_paths(item, _seen=_seen))
            else:
                paths.extend(_source_paths(child, _seen=_seen))
        return _unique(paths)
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        paths: list[Any] = []
        for item in value:
            paths.extend(_source_paths(item, _seen=_seen))
        return _unique(paths)
    return []


def _check_powers(value: Any, label: str) -> tuple[int, int, int, int]:
    if (
        not isinstance(value, Sequence)
        or isinstance(value, (str, bytes, bytearray))
        or len(value) != 4
        or any(not isinstance(item, int) or isinstance(item, bool) for item in value)
    ):
        raise ValueError(f"{label} must contain four integer formal powers")
    return tuple(int(item) for item in value)  # type: ignore[return-value]


def _check_term(value: Any, label: str, *, allow_missing_source: bool = False) -> dict[str, Any]:
    if not isinstance(value, Mapping) or "powers" not in value or "coefficient" not in value:
        raise ValueError(f"{label} requires powers and coefficient")
    if not _is_interval(value["coefficient"]):
        raise TypeError(f"{label}.coefficient must be a directed interval")
    powers = _check_powers(value["powers"], f"{label}.powers")
    info = _origin_info(value)
    if not info["source_ids"] and not allow_missing_source:
        raise ValueError(f"{label} must retain source_ids")
    result = dict(value)
    result.update(
        {
            "powers": powers,
            "coefficient": value["coefficient"],
            "source_ids": list(info["source_ids"]),
            "source_paths": _source_paths(value),
            "error_classes": list(info["error_classes"]),
            "uncertainty": bool(info["uncertainty"]),
            "error_only": bool(info["error_only"]),
            "symmetric_error": bool(info["symmetric_error"]),
            "origin": value.get("origin", value),
        }
    )
    return result


def _order_key(value: Any, label: str) -> int:
    if isinstance(value, bool):
        raise ValueError(f"boolean is not an ordinary-Z order in {label}")
    if isinstance(value, int):
        order = value
    else:
        text = str(value).strip().lower()
        if text.startswith("z"):
            text = text[1:]
        try:
            order = int(text)
        except ValueError as exc:
            raise ValueError(f"invalid ordinary-Z order {value!r} in {label}") from exc
    if order not in ORDINARY_Z_ORDERS:
        raise ValueError(f"ordinary-Z order must be 0 or 1 in {label}")
    return order


def _cell_key(value: Any, label: str) -> int:
    if isinstance(value, bool):
        raise ValueError(f"boolean is not a q-cell index in {label}")
    if isinstance(value, int):
        cell = value
    else:
        text = str(value).strip().lower()
        if text.startswith("cell"):
            text = text[4:]
        try:
            cell = int(text)
        except ValueError as exc:
            raise ValueError(f"invalid q-cell {value!r} in {label}") from exc
    if cell not in ALL_CELLS:
        raise ValueError(f"q-cell must be 0, 1, 2, or 3 in {label}")
    return cell


def _term_list(value: Any, label: str, *, allow_empty: bool) -> list[dict[str, Any]]:
    if isinstance(value, Mapping) and "powers" in value:
        values = [value]
    elif isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        values = list(value)
    else:
        raise TypeError(f"{label} must be a sector sequence")
    if not values and not allow_empty:
        raise ValueError(f"{label} must not be empty")
    return [_check_term(item, f"{label}[{index}]") for index, item in enumerate(values)]


def _unwrap_cells(value: Any) -> Any:
    if not isinstance(value, Mapping):
        return value
    for key in ("cells", "by_cell", "q_cells", "cell_sectors", "sectors_by_cell"):
        child = value.get(key)
        if isinstance(child, Mapping):
            return child
    return value


def _normalize_cell_orders(value: Any, label: str) -> dict[int, dict[int, list[dict[str, Any]]]]:
    """Accept either ``cell -> Z -> terms`` or ``Z -> cell -> terms``."""

    value = _unwrap_cells(value)
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        if len(value) not in (3, 4):
            raise ValueError(f"{label} sequence must contain three or four q-cells")
        value = {index: item for index, item in enumerate(value)}
    if not isinstance(value, Mapping):
        raise TypeError(f"{label} must be a q-cell/order mapping")
    keys = list(value)
    # A direct Z0/Z1 mapping may carry a cell index on every sector.  This is
    # accepted only when the index is explicit; a single point value is never
    # silently copied to all cells.
    try:
        order_keys = [_order_key(item, label) for item in keys]
    except ValueError:
        order_keys = []
    if len(keys) == 2 and set(order_keys) == set(ORDINARY_Z_ORDERS):
        rows = {order: value[key] for key, order in zip(keys, order_keys)}
        buckets: dict[int, dict[int, list[dict[str, Any]]]] = {
            cell: {0: [], 1: []} for cell in ALL_CELLS
        }
        for order in ORDINARY_Z_ORDERS:
            raw_terms = rows[order]
            if isinstance(raw_terms, Mapping) and "powers" not in raw_terms:
                # More specific Z -> cell form.
                for raw_cell, cell_terms in raw_terms.items():
                    cell = _cell_key(raw_cell, f"{label}.Z{order}")
                    buckets[cell][order] = _term_list(
                        cell_terms,
                        f"{label}.Z{order}.cell{cell}",
                        allow_empty=cell == 3,
                    )
            else:
                terms = _term_list(raw_terms, f"{label}.Z{order}", allow_empty=False)
                for index, term in enumerate(terms):
                    raw_cell = term.get("cell", term.get("q_cell", term.get("phase_cell")))
                    if raw_cell is None:
                        raise ValueError(f"{label}.Z{order}[{index}] needs an explicit cell")
                    cell = _cell_key(raw_cell, f"{label}.Z{order}[{index}]")
                    buckets[cell][order].append(term)
        return buckets

    # Cell -> order form.  Parse every key as a cell and require explicit
    # ordinary-Z rows in each nonempty cell.
    cells: dict[int, dict[int, list[dict[str, Any]]]] = {}
    for raw_cell, raw_orders in value.items():
        cell = _cell_key(raw_cell, label)
        if not isinstance(raw_orders, Mapping):
            raise TypeError(f"{label}.cell{cell} must provide Z0 and Z1")
        rows: dict[int, Any] = {}
        for raw_order, row in raw_orders.items():
            order = _order_key(raw_order, f"{label}.cell{cell}")
            if order in rows:
                raise ValueError(f"duplicate Z{order} in {label}.cell{cell}")
            rows[order] = row
        if set(rows) != set(ORDINARY_Z_ORDERS):
            raise ValueError(f"{label}.cell{cell} must provide exactly Z0 and Z1")
        cells[cell] = {
            order: _term_list(rows[order], f"{label}.cell{cell}.Z{order}", allow_empty=cell == 3)
            for order in ORDINARY_Z_ORDERS
        }
    return cells


def _check_active_cells(value: Mapping[int, Mapping[int, Sequence[Any]]], label: str) -> None:
    missing = [cell for cell in ACTIVE_CELLS if cell not in value]
    if missing:
        raise ValueError(f"{label} must provide active q-cells {missing}")
    for cell in ACTIVE_CELLS:
        for order in ORDINARY_Z_ORDERS:
            if not value[cell].get(order):
                raise ValueError(f"{label}.cell{cell}.Z{order} must not be empty")


def _normalize_leading(value: Any, label: str) -> dict[int, dict[int, list[dict[str, Any]]]]:
    result = _normalize_cell_orders(value, label)
    _check_active_cells(result, label)
    for cell in ALL_CELLS:
        result.setdefault(cell, {0: [], 1: []})
    return result


def _canonical_part(value: Any, part: str) -> Any:
    if not isinstance(value, Mapping):
        return value
    if "fields" in value and isinstance(value["fields"], Mapping):
        value = value["fields"]
    names = _FIELD_ALIASES[part]
    channels: dict[str, Any] = {}
    for channel, aliases in names.items():
        present = [alias for alias in aliases if alias in value]
        if len(present) > 1:
            raise ValueError(
                f"field_sectors[{part}] uses multiple aliases for {channel}: {present}"
            )
        if present:
            channels[channel] = value[present[0]]
    if channels:
        return channels
    # J1/J2 are F-only and are often passed directly as a cell/order packet.
    if part in ("J1", "J2"):
        return {"swirl": value}
    return value


def _normalize_fields(value: Any) -> dict[str, dict[str, dict[int, dict[int, list[dict[str, Any]]]]]]:
    if not isinstance(value, Mapping):
        raise TypeError("field_sectors must be a mapping")
    if "parts" in value and isinstance(value["parts"], Mapping):
        value = value["parts"]
    raw_parts: dict[str, Any] = {}
    for raw_part, item in value.items():
        key = _PART_ALIASES.get(str(raw_part).strip().lower(), str(raw_part).strip())
        if key in PARTS:
            raw_parts[key] = item
    missing = [part for part in PARTS if part not in raw_parts]
    if missing:
        raise ValueError("field_sectors must provide K0, J1, and J2: missing " + ", ".join(missing))
    result: dict[str, dict[str, dict[int, dict[int, list[dict[str, Any]]]]]] = {}
    for part in PARTS:
        packet = _canonical_part(raw_parts[part], part)
        if not isinstance(packet, Mapping):
            raise TypeError(f"field_sectors[{part}] must provide field channels")
        result[part] = {}
        for channel in PART_CHANNELS[part]:
            raw = packet.get(channel)
            if raw is None:
                aliases = _FIELD_ALIASES[part][channel]
                present = [alias for alias in aliases if alias in packet]
                if len(present) > 1:
                    raise ValueError(
                        f"field_sectors[{part}] uses multiple aliases for {channel}: {present}"
                    )
                if present:
                    raw = packet[present[0]]
            if raw is None:
                raise ValueError(f"field_sectors[{part}] is missing {channel} field sectors")
            cells = _normalize_cell_orders(raw, f"field_sectors[{part}].{channel}")
            _check_active_cells(cells, f"field_sectors[{part}].{channel}")
            for cell in ALL_CELLS:
                cells.setdefault(cell, {0: [], 1: []})
            result[part][channel] = cells
    return result


def _validate_source_atoms(source_atoms: Any) -> None:
    if isinstance(source_atoms, Mapping):
        names = set(source_atoms)
    elif isinstance(source_atoms, Sequence) and not isinstance(source_atoms, (str, bytes, bytearray)):
        names = {
            item.get("name", item.get("atom", item.get("power"))) if isinstance(item, Mapping) else item
            for item in source_atoms
        }
    else:
        raise ValueError("source_atoms must name the four formal source atoms")
    if names != set(SOURCE_POWER_BASIS):
        raise ValueError("source_atoms must contain exactly hb, F0_anchor, Pstar_squared, Lambda")


def _certified(certificate: Mapping[str, Any], *names: str) -> bool:
    return any(certificate.get(name) is True for name in names)


def _validate_certificate(certificate: Any) -> None:
    if not isinstance(certificate, Mapping):
        raise ValueError("caller geometry and cutoff certificates are required")
    requirements = (
        ("hb_strictly_between_zero_and_half", "strict_hb_admitted"),
        ("sc_phase_Z_independent", "hb_sc_Z_independent"),
        ("hb_epsilon_Z_independent", "rho_radius_Z_independent"),
        ("original_cutoff_mass_5_over_8", "cutoff_mass_5_over_8_exact", "actual_cutoff_complete_support_measure"),
        ("original_cutoff_mass_by_cell", "cutoff_mass_by_cell_exact", "eta_cell_masses_exact"),
    )
    missing = [names[0] for names in requirements if not _certified(certificate, *names)]
    if missing:
        raise ValueError("missing caller certificate: " + ", ".join(missing))


def _symmetric(c: Any, value: Any) -> Any:
    lower, upper = _scalar_bounds(value)
    radius = max(abs(lower), abs(upper))
    return _make_interval(c, -radius, radius)


def _cell_bounds(c: Any, *, four_cells: bool, q_cells: Any = None) -> dict[int, tuple[Any, Any]]:
    if four_cells:
        fractions = ((0, 4, 1, 4), (1, 4, 1, 2), (1, 2, 3, 4), (3, 4, 1, 1))
    else:
        # The third supplied field enclosure may cover [1/2,1].  The exact
        # cutoff mass is still 1/8 because the [3/4,1] tail has zero mass.
        fractions = ((0, 4, 1, 4), (1, 4, 1, 2), (1, 2, 1, 1), (3, 4, 1, 1))
    result = {
        cell: (c.mpf(lo_num) / c.mpf(lo_den), c.mpf(hi_num) / c.mpf(hi_den))
        for cell, (lo_num, lo_den, hi_num, hi_den) in enumerate(fractions)
    }
    if q_cells is None:
        return result
    if isinstance(q_cells, Mapping) and "cells" in q_cells and isinstance(q_cells["cells"], Mapping):
        q_cells = q_cells["cells"]
    if not isinstance(q_cells, Mapping):
        raise TypeError("q_cells must be a mapping when supplied")
    if any(_cell_key(raw_cell, "q_cells") == 3 for raw_cell in q_cells):
        raise ValueError("q_cells must use exactly the three active intervals; cell 3 has zero mass")
    for raw_cell, raw in q_cells.items():
        cell = _cell_key(raw_cell, "q_cells")
        if not isinstance(raw, Mapping):
            raise TypeError(f"q_cells.cell{cell} must be a mapping")
        lower = raw.get("q_lower", raw.get("lower"))
        upper = raw.get("q_upper", raw.get("upper"))
        if lower is None or upper is None:
            raise ValueError(f"q_cells.cell{cell} requires q_lower and q_upper")
        result[cell] = (c.mpf(lower), c.mpf(upper))
    expected = ((0, 4, 1, 4), (1, 4, 1, 2), (1, 2, 1, 1))
    for cell, (lo_num, lo_den, hi_num, hi_den) in enumerate(expected):
        lo, hi = result[cell]
        if lo != c.mpf(lo_num) / c.mpf(lo_den) or hi != c.mpf(hi_num) / c.mpf(hi_den):
            raise ValueError("q_cells must strictly cover [0,1/4], [1/4,1/2], [1/2,1] without overlap")
    return result


def _has_explicit_cell3(value: Any) -> bool:
    """Detect an explicitly supplied fourth q-cell before empty defaults are added."""

    value = _unwrap_cells(value)
    if isinstance(value, Mapping) and ("powers" in value or "coefficient" in value):
        return False
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        if value and all(isinstance(item, Mapping) and ("powers" in item or "coefficient" in item) for item in value):
            return False
        return len(value) == 4
    if not isinstance(value, Mapping):
        return False
    for raw_key, raw_value in value.items():
        try:
            if _cell_key(raw_key, "q-cell input") == 3:
                return True
        except ValueError:
            pass
        if isinstance(raw_value, (Mapping, Sequence)) and not isinstance(raw_value, (str, bytes, bytearray)):
            if _has_explicit_cell3(raw_value):
                return True
    return False


def _geometry(c: Any, sc: Any, cell_bounds: Mapping[int, tuple[Any, Any]], cell: int) -> dict[str, Any]:
    q_lower, q_upper = cell_bounds[cell]
    q_interval = _make_interval(c, q_lower, q_upper)
    mass_num, mass_den = ETA_CELL_MASS[cell]
    eta_mass = c.mpf(mass_num) / c.mpf(mass_den)
    sc_iv = c.mpf(sc)
    _, sc_upper = _scalar_bounds(sc_iv)
    exponent_upper = sc_upper * q_upper / 8
    exp_upper = c.exp(_make_interval(c, c.mpf(0), exponent_upper))
    rho_interval = _make_interval(c, c.mpf(4), c.mpf(4) * exp_upper)
    rho_squared_interval = rho_interval * rho_interval
    axial_factor = c.mpf(1) / c.mpf(2) * sc_iv * eta_mass * rho_interval
    swirl_factor = -c.mpf(1) / c.mpf(2) * sc_iv * eta_mass * rho_squared_interval
    return {
        "cell": cell,
        "q_lower": q_lower,
        "q_upper": q_upper,
        "q_interval": q_interval,
        "eta_mass": eta_mass,
        "eta_mass_exact": True,
        "rho_interval": rho_interval,
        "rho_squared_interval": rho_squared_interval,
        "rho_formula": "4*exp(hb*sc*q/4)",
        "d_rho_formula": "hb*(sc/4)*rho*dq",
        "hb_directed_range": "0 < hb < 1/2; upper bound 1/2 used for rho enclosure",
        "factors": {"axial": axial_factor, "swirl": swirl_factor},
        "factor_formulas": {
            "axial": "+sc*eta_mass*rho/2, with formal hb and Lambda^-1",
            "swirl": "-sc*eta_mass*rho^2/2, with formal hb and Lambda^-2",
        },
    }


def _product(
    leading: Mapping[str, Any],
    field: Mapping[str, Any],
    factor: Any,
    *,
    c: Any,
    part: str,
    channel: str,
    cell: int,
    order: int,
    piece: str,
    geometry: Mapping[str, Any],
) -> dict[str, Any]:
    lead_info = _origin_info(leading)
    field_info = _origin_info(field)
    field_error = bool(part == "J2" or field_info["error_only"])
    error_only = bool(field_error or lead_info["error_only"])
    symmetric_error = bool(part == "J2" or field_info["symmetric_error"] or lead_info["symmetric_error"])
    classes = _unique(list(lead_info["error_classes"]) + list(field_info["error_classes"]))
    if symmetric_error and "symmetric_error" not in classes:
        classes.append("symmetric_error")
    uncertainty = bool(error_only or lead_info["uncertainty"] or field_info["uncertainty"])
    source_ids = _unique(list(lead_info["source_ids"]) + list(field_info["source_ids"]))
    source_paths = _unique(list(_source_paths(leading)) + list(_source_paths(field)))
    short_piece = piece
    legacy_piece = {
        "lead0*field0": "leading0*field0",
        "lead1*field0": "leading1*field0",
        "lead0*field1": "leading0*field1",
    }.get(piece, piece)
    powers0 = tuple(a + b for a, b in zip(leading["powers"], field["powers"]))
    external = PHYSICAL_EXTERNAL_POWERS[channel]
    powers = tuple(powers0[index] + external[index] for index in range(4))
    coefficient = leading["coefficient"] * field["coefficient"] * factor
    if symmetric_error:
        coefficient = _symmetric(c, coefficient)
    path = {
        "part": part,
        "channel": channel,
        "field": "U" if channel == "axial" else "F",
        "cell": cell,
        "phase_cell": cell,
        "ordinary_Z_order": order,
        "ordinary_Z_piece": short_piece,
        "ordinary_Z_piece_legacy": legacy_piece,
        "leading_source_ids": list(lead_info["source_ids"]),
        "field_source_ids": list(field_info["source_ids"]),
        "leading_source_paths": list(_source_paths(leading)),
        "field_source_paths": list(_source_paths(field)),
    }
    origin = {
        "operation": "e00_physical_m15_product",
        "part": part,
        "channel": channel,
        "field": "U" if channel == "axial" else "F",
        "cell": cell,
        "phase_cell": cell,
        "ordinary_Z_order": order,
        "ordinary_Z_piece": short_piece,
        "ordinary_Z_piece_legacy": legacy_piece,
        "path": path,
        "paths": [path],
        "source_ids": source_ids,
        "source_paths": source_paths,
        "leading_source_ids": list(lead_info["source_ids"]),
        "field_source_ids": list(field_info["source_ids"]),
        "leading_source_paths": list(_source_paths(leading)),
        "field_source_paths": list(_source_paths(field)),
        "error_classes": classes,
        "uncertainty": uncertainty,
        "error_only": error_only,
        "symmetric_error": symmetric_error,
        "eta_mass": geometry["eta_mass"],
        "q_interval": geometry["q_interval"],
        "rho_interval": geometry["rho_interval"],
        "rho_squared_interval": geometry["rho_squared_interval"],
        "physical_factor": factor,
        "formal_external_powers": external,
        "coefficient_contains_formal_hb": False,
        "field_is_already_source_bound": True,
        "j1_j2_reapplied": False,
        "symbolic_cancellation_claimed": False,
    }
    return {
        "powers": powers,
        "coefficient": coefficient,
        "origin": origin,
        "source_ids": source_ids,
        "source_paths": source_paths,
        "leading_source_ids": list(lead_info["source_ids"]),
        "field_source_ids": list(field_info["source_ids"]),
        "leading_source_paths": list(_source_paths(leading)),
        "field_source_paths": list(_source_paths(field)),
        "error_classes": classes,
        "uncertainty": uncertainty,
        "error_only": error_only,
        "symmetric_error": symmetric_error,
        "part": part,
        "channel": channel,
        "field": "U" if channel == "axial" else "F",
        "cell": cell,
        "phase_cell": cell,
        "ordinary_Z_order": order,
        "ordinary_Z_piece": piece,
        "ordinary_Z_piece_legacy": legacy_piece,
        "physical_factor": factor,
        "formal_external_powers": external,
        "source_term_is_whole_q_cell_enclosure": True,
        "leading_term_is_whole_q_cell_enclosure": True,
        "field_is_already_source_bound": True,
        "j1_j2_reapplied": False,
    }


def _assemble_part(
    c: Any,
    part: str,
    leading: Mapping[str, Mapping[int, Mapping[int, Sequence[Mapping[str, Any]]]]],
    fields: Mapping[str, Mapping[int, Mapping[int, Sequence[Mapping[str, Any]]]]],
    geometries: Mapping[int, Mapping[str, Any]],
) -> dict[str, Any]:
    channels: dict[str, dict[int, dict[int, list[dict[str, Any]]]]] = {}
    for channel in PART_CHANNELS[part]:
        channel_result = {order: {cell: [] for cell in ALL_CELLS} for order in ORDINARY_Z_ORDERS}
        channel_leading = leading[channel]
        channel_fields = fields[channel]
        for cell in ACTIVE_CELLS:
            geometry = geometries[cell]
            factor = geometry["factors"][channel]
            for order in ORDINARY_Z_ORDERS:
                if order == 0:
                    pieces = (("lead0*field0", 0, 0),)
                else:
                    pieces = (("lead1*field0", 1, 0), ("lead0*field1", 0, 1))
                for piece, leading_order, field_order in pieces:
                    for lead in channel_leading[cell][leading_order]:
                        for field in channel_fields[cell][field_order]:
                            channel_result[order][cell].append(
                                _product(
                                    lead,
                                    field,
                                    factor,
                                    c=c,
                                    part=part,
                                    channel=channel,
                                    cell=cell,
                                    order=order,
                                    piece=piece,
                                    geometry=geometry,
                                )
                            )
        channels[channel] = channel_result
    counts = {
        channel: {
            f"Z{order}": {str(cell): len(channels[channel][order][cell]) for cell in ALL_CELLS}
            for order in ORDINARY_Z_ORDERS
        }
        for channel in PART_CHANNELS[part]
    }
    return {
        "part": part,
        "channels": channels,
        "paths": channels,
        "swirl": channels["swirl"],
        "active_channels": list(PART_CHANNELS[part]),
        "axial_not_added": part != "K0",
        "input_field_channel": "U/F" if part == "K0" else "F",
        "input_sector_counts": counts,
        "ordinary_Z0_rule": "leading0*field0",
        "ordinary_Z1_rule": "leading1*field0 + leading0*field1",
        "all_terms_retained": True,
        "source_paths_retained": True,
        "error_only_retained": part == "J2",
        "symmetric_error_retained": part == "J2",
        "j1_j2_reapplied": False,
    }


def _raw_part_aliases(field_sectors: Any, k0: Any, j1: Any, j2: Any) -> Any:
    if field_sectors is not None and any(item is not None for item in (k0, j1, j2)):
        raise ValueError("field_sectors cannot be mixed with k0, j1, or j2 aliases")
    if field_sectors is not None:
        return field_sectors
    if k0 is None or j1 is None or j2 is None:
        raise ValueError("field_sectors or all k0, j1, and j2 inputs are required")
    return {"K0": k0, "J1": j1, "J2": j2}


def build(
    c: Any,
    *,
    sc: Any,
    leading_U: Any,
    leading_F: Any,
    field_sectors: Any = None,
    source_family: Any = None,
    source_atoms: Any = None,
    certificate: Any = None,
    q_cells: Any = None,
    k0: Any = None,
    j1: Any = None,
    j2: Any = None,
) -> dict[str, Any]:
    """Combine caller-owned leading and field q-cell enclosures through Z1.

    ``leading_U`` and ``leading_F`` are ``cell -> {Z0, Z1} -> sectors``.
    ``field_sectors`` is ``{"K0": {"U": ..., "F": ...},
    "J1": {"F": ...}, "J2": {"F": ...}}`` with the same cell layout.
    Three-cell packets use the non-overlapping q intervals ``[0,1/4]``,
    ``[1/4,1/2]``, and ``[1/2,1]``.  A separate ``zero_mass_tail`` record
    documents ``[3/4,1]`` without making it a fourth overlapping input cell.
    The active certified masses are exactly ``(1/4, 1/4, 1/8)``.
    """

    if not hasattr(c, "mpf") or not hasattr(c, "make_mpf") or not hasattr(c, "exp"):
        raise TypeError("c must provide directed interval arithmetic")
    if source_family is None or (isinstance(source_family, str) and not source_family.strip()):
        raise ValueError("source_family is required")
    if source_atoms is None:
        raise ValueError("source_atoms is required")
    _validate_source_atoms(source_atoms)
    _validate_certificate(certificate)

    sc_iv = c.mpf(sc)
    if not _is_interval(sc_iv):
        raise TypeError("sc must be a directed interval")
    sc_lower, sc_upper = _scalar_bounds(sc_iv)
    limit_lower, _ = _scalar_bounds(c.mpf(1) / c.mpf(1024))
    if sc_lower <= 0 or sc_upper >= limit_lower:
        raise ValueError("sc must satisfy 0 < sc < 1/1024")

    leading_u = _normalize_leading(leading_U, "leading_U")
    leading_f = _normalize_leading(leading_F, "leading_F")
    raw_fields = _raw_part_aliases(field_sectors, k0, j1, j2)
    if _has_explicit_cell3(leading_U) or _has_explicit_cell3(leading_F) or _has_explicit_cell3(raw_fields):
        raise ValueError("input q cells must be exactly the three active intervals; cell 3 is a zero-mass tail")
    fields = _normalize_fields(raw_fields)
    # The caller's third field enclosure is the single [1/2,1] packet.  A
    # separate [3/4,1] packet is represented only by the exact zero mass in
    # q_cells[3], never as another input sector.
    cell_bounds = _cell_bounds(c, four_cells=False, q_cells=q_cells)
    geometries = {cell: _geometry(c, sc_iv, cell_bounds, cell) for cell in ALL_CELLS}

    parts = {
        "K0": _assemble_part(c, "K0", {"axial": leading_u, "swirl": leading_f}, fields["K0"], geometries),
        "J1": _assemble_part(c, "J1", {"swirl": leading_f}, fields["J1"], geometries),
        "J2": _assemble_part(c, "J2", {"swirl": leading_f}, fields["J2"], geometries),
    }
    eta_masses = [c.mpf(num) / c.mpf(den) for num, den in ETA_CELL_MASS]
    q_rows = [
        {
            "cell": cell,
            "q_lower": geometries[cell]["q_lower"],
            "q_upper": geometries[cell]["q_upper"],
            "q_interval": geometries[cell]["q_interval"],
            "eta_mass": eta_masses[cell],
            "eta_mass_exact": True,
            "velocity_cutoff_active": cell in ACTIVE_CELLS,
            "field_and_leading_are_whole_q_cell_enclosures": cell in ACTIVE_CELLS,
        }
        for cell in ACTIVE_CELLS
    ]
    zero_mass_tail = {
        "q_lower": c.mpf(3) / c.mpf(4),
        "q_upper": c.mpf(1),
        "q_interval": _make_interval(c, c.mpf(3) / c.mpf(4), c.mpf(1)),
        "eta_mass": c.mpf(0),
        "eta_mass_exact": True,
        "velocity_cutoff_active": False,
        "zero_mass": True,
    }
    active_eta_masses = eta_masses[:3]
    return {
        "schema": SCHEMA,
        "status": "source_bound_physical_m15_cell_projection",
        "source_family": source_family,
        "source_atoms": source_atoms,
        "source_power_basis": list(SOURCE_POWER_BASIS),
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "parts": parts,
        "K0": parts["K0"],
        "J1": parts["J1"],
        "J2": parts["J2"],
        "q_cells": q_rows,
        "zero_mass_tail": zero_mass_tail,
        "radial_phase_enclosures": geometries,
        "cell_layout": "three_active_enclosures_plus_zero_mass_tail",
        "active_cells": list(ACTIVE_CELLS),
        "zero_mass_cell": 3,
        "eta_profile": "1-sigma(4q-2)",
        "eta_cell_masses": active_eta_masses,
        "cutoff_mass_by_cell": active_eta_masses,
        "eta_cell_masses_with_zero_tail": eta_masses,
        "cutoff_mass_by_cell_with_zero_tail": eta_masses,
        "eta_cell_masses_exact": True,
        "eta_integral": c.mpf(5) / c.mpf(8),
        "cutoff_mass": c.mpf(5) / c.mpf(8),
        "positive_cutoff_mass_exact": True,
        "certificate": dict(certificate),
        "cutoff_certificate": {
            "mass_by_cell": [(1, 4), (1, 4), (1, 8), (0, 1)],
            "total_mass": (5, 8),
            "mass_is_caller_certified": True,
            "tail_3_over_4_to_1_zero": True,
        },
        "rho_formula": "4*exp(hb*sc*q/4)",
        "d_rho_formula": "hb*(sc/4)*rho*dq",
        "axial_formula": "+2*Lambda^-1*integral(eta*Ulead*U1*d_rho)",
        "swirl_formula": "-2*Lambda^-2*integral(eta*rho*Flead*F1*d_rho)",
        "physical_external_powers": dict(PHYSICAL_EXTERNAL_POWERS),
        "formal_hb_measure_retained": True,
        "formal_powers_not_materialized": True,
        "hb_directed_bound_used_for_geometry": True,
        "hb_point_value_selected": False,
        "q_point_value_selected": False,
        "geometry_Z_independent": True,
        "field_q_cell_enclosures_required": True,
        "leading_q_cell_enclosures_required": True,
        "source_ids_retained": True,
        "source_paths_retained": True,
        "error_only_retained": True,
        "symmetric_error_retained": True,
        "ordinary_Z1_product_rule": "leading1*field0 + leading0*field1",
        "ordinary_Z1_leibniz_terms_complete": True,
        "j1_j2_reapplied": False,
        "chi_reapplied": False,
        "cutoff_reapplied": False,
        "P1_is_m15": False,
        "prefix_mutual_exclusivity_claimed": False,
        "old_prefix_disjointness_not_proved": True,
        "ledger_appended": False,
        "appended_to_ledger": False,
        "ledger_written": False,
        "allword_tail_added": False,
        "complete_m15_computed": False,
        "complete_m15": False,
        "scope": "standalone source-bound E00 physical K0/J1/J2 cell projection only",
    }


build_physical_m15 = build
evaluate = build
