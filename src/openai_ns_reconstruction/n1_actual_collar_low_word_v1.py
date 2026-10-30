"""Local first-word ``Z0``/``Z1`` Green transport on the first collar.

This module is a source-only leaf.  It is deliberately based on the compact
signed Green-sector implementation, but it does not import that (or any
other) producer.  The public :func:`enclose_low` function accepts all caller
receipts needed to validate the complete low-order input contract:

* ``core_W0_Z`` contains the four core homogeneous orders ``0..3``;
* ``forcing_Z`` contains actual ``g_Z0..g_Z3`` rows ``4,5,6``;
* ``matrix_Z`` contains actual ``C0_Z0..Z2`` and ``B1_Z0..Z2`` entries.

Only the local first-word contributions ``W1_Z0`` and ``W1_Z1`` are formed.
The missing core ``W1`` inlet is not an argument, is not replaced by zero,
and is not represented as a completed word.  Likewise, no all-word tail,
endpoint cancellation, or second ordinary-Z result is claimed.  The four
source powers remain formal, signed kernels retain their orientation, and the
``hb`` Green measure remains formal in every sector.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any


SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
FIELDS = ("F1", "Uz1", "K1", "P1", "d_x_F1", "d_x_Uz1")
ROWS = ("5", "6", "6", "4", "5", "6")
C0_KEYS = frozenset(("4,1", "5,1", "5,2", "5,3", "5,5", "6,1", "6,2", "6,3", "6,4", "6,6"))
B1_KEYS = frozenset(("5,1", "5,2", "5,3", "6,2", "6,3", "6,4"))


def _is_interval(value: Any) -> bool:
    return hasattr(value, "_mpi_") and hasattr(value, "a") and hasattr(value, "b")


def _raw_endpoint(value: Any, upper: bool) -> Any:
    raw = getattr(value, "_mpi_", None)
    if raw is not None:
        return raw[1 if upper else 0]
    return value


def _make_interval(c: Any, lower: Any, upper: Any) -> Any:
    return c.make_mpf((_raw_endpoint(lower, False), _raw_endpoint(upper, True)))


def _as_ids(value: Any) -> list[Any]:
    if value is None:
        return []
    if isinstance(value, (str, bytes, bytearray)):
        return [value]
    if isinstance(value, Mapping):
        return [value]
    if isinstance(value, Sequence):
        return list(value)
    return [value]


def _as_labels(value: Any) -> list[Any]:
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


def _unique(values: Sequence[Any]) -> list[Any]:
    result: list[Any] = []
    seen: set[Any] = set()
    for value in values:
        key = _freeze(value)
        if key not in seen:
            seen.add(key)
            result.append(value)
    return result


def _positive_int(value: Any, default: int = 1) -> int:
    try:
        result = int(value)
    except (TypeError, ValueError):
        result = default
    return max(1, result)


def _origin_info(value: Any) -> dict[str, Any]:
    """Read shallow/deep provenance without discarding source/error unions."""

    if isinstance(value, Mapping):
        ids = list(_as_ids(value.get("source_ids", value.get("sources"))))
        classes = list(_as_labels(value.get(
            "error_classes", value.get("uncertainty_class", value.get("error_class")))))
        operation = str(value.get("operation", "")).lower()
        error_kind = str(value.get("error_kind", value.get("kind", ""))).lower()
        uncertain = bool(value.get("uncertainty", value.get("is_uncertainty", False)))
        uncertain = uncertain or any(token in operation or token in error_kind
                                     for token in ("cauchy", "symmetric", "history"))
        nested: list[dict[str, Any]] = []
        for key in ("origin", "metadata", "parents"):
            child = value.get(key)
            if child is None or child is value:
                continue
            if isinstance(child, Sequence) and not isinstance(child, (str, bytes, bytearray)):
                nested.extend(_origin_info(item) for item in child)
            else:
                nested.append(_origin_info(child))
        for info in nested:
            ids.extend(info.get("source_ids", ()))
            classes.extend(info.get("error_classes", ()))
            uncertain = uncertain or bool(info.get("uncertainty", False))
        raw_path = value.get("path_count", value.get("aggregate_path_count"))
        if raw_path is None:
            raw_path = max((info.get("path_count", 1) for info in nested), default=1)
        raw_source = value.get("source_count", value.get("aggregate_source_count"))
        if raw_source is None:
            raw_source = max((info.get("source_count", info.get("path_count", 1))
                              for info in nested), default=raw_path)
        return {
            "source_ids": _unique(ids),
            "uncertainty": bool(uncertain),
            "error_classes": _unique(classes),
            "path_count": _positive_int(raw_path),
            "source_count": _positive_int(raw_source),
        }
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        nested = [_origin_info(item) for item in value
                  if isinstance(item, (Mapping, Sequence))
                  and not isinstance(item, (str, bytes, bytearray))]
        if nested:
            ids: list[Any] = []
            classes: list[Any] = []
            uncertain = False
            path_count = 1
            source_count = 1
            for info in nested:
                ids.extend(info["source_ids"])
                classes.extend(info["error_classes"])
                uncertain = uncertain or info["uncertainty"]
                path_count = max(path_count, info["path_count"])
                source_count = max(source_count, info["source_count"])
            return {"source_ids": _unique(ids), "uncertainty": uncertain,
                    "error_classes": _unique(classes), "path_count": path_count,
                    "source_count": source_count}
    return {"source_ids": [], "uncertainty": False, "error_classes": [],
            "path_count": 1, "source_count": 1}


def _provenance(term: Any) -> dict[str, Any]:
    return _origin_info(term)


def _new_origin(operation: str, *inputs: Any, **metadata: Any) -> dict[str, Any]:
    infos = [_provenance(item) for item in inputs]
    ids: list[Any] = []
    classes: list[Any] = []
    uncertain = bool(metadata.get("uncertainty", False))
    if operation in {"product", "multiply"}:
        path_count = 1
        source_count = 1
        for info in infos:
            path_count *= int(info.get("path_count", 1))
            source_count *= int(info.get("source_count", 1))
    elif operation in {"aggregate", "add", "sum"}:
        path_count = 0
        source_count = 0
        for info in infos:
            path_count += int(info.get("path_count", 1))
            source_count += int(info.get("source_count", 1))
    else:
        path_count = int(infos[0].get("path_count", 1)) if infos else 1
        source_count = int(infos[0].get("source_count", 1)) if infos else 1
    for info in infos:
        ids.extend(info.get("source_ids", ()))
        classes.extend(info.get("error_classes", ()))
        uncertain = uncertain or bool(info.get("uncertainty", False))
    if "source_ids" in metadata:
        ids.extend(_as_ids(metadata["source_ids"]))
    if "error_classes" in metadata:
        classes.extend(_as_labels(metadata["error_classes"]))
    if metadata.get("cauchy_error") or metadata.get("symmetric_error"):
        uncertain = True
    origin: dict[str, Any] = {
        "operation": operation,
        "source_ids": _unique(ids),
        "uncertainty": bool(uncertain),
        "error_classes": _unique(classes),
        "path_count": max(1, int(path_count)),
        "aggregate_path_count": max(1, int(path_count)),
        "source_count": max(1, int(source_count)),
        "aggregate_source_count": max(1, int(source_count)),
        "source_id_count": len(_unique(ids)),
    }
    for key, value in metadata.items():
        if key not in {"source_ids", "uncertainty", "error_classes", "cauchy_error", "symmetric_error"}:
            origin[key] = value
    return origin


def _check_sector(term: Any) -> None:
    if not isinstance(term, Mapping):
        raise TypeError("each sector must be a mapping")
    if set(("powers", "coefficient")) - set(term):
        raise ValueError("sector metadata missing")
    powers = term["powers"]
    if (not isinstance(powers, Sequence) or isinstance(powers, (str, bytes, bytearray))
            or len(powers) != 4
            or any(not isinstance(power, int) or isinstance(power, bool) for power in powers)):
        raise ValueError("four integer original source powers required")
    if not _is_interval(term["coefficient"]):
        raise TypeError("directed interval coefficient required")


def _term_list(value: Any, label: str) -> list[Mapping[str, Any]]:
    if isinstance(value, Mapping) and "powers" in value:
        values = [value]
    elif isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        values = list(value)
    else:
        raise TypeError(f"{label} must be a sector sequence")
    for term in values:
        _check_sector(term)
    return values


def _scale(terms: Sequence[Mapping[str, Any]], factor: Any, hb_power: int,
           origin: Any) -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    for term in _term_list(terms, "scale input"):
        powers = list(term["powers"])
        powers[0] += hb_power
        result.append({
            "powers": tuple(powers),
            "coefficient": term["coefficient"] * factor,
            "origin": _new_origin("scale", term, operation_label=origin,
                                   hb_power=hb_power),
        })
    return result


def _product(left: Sequence[Mapping[str, Any]], right: Sequence[Mapping[str, Any]],
             multiplier: Any, origin: Any) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    for left_term in _term_list(left, "left product input"):
        for right_term in _term_list(right, "right product input"):
            out.append({
                "powers": tuple(a + b for a, b in zip(left_term["powers"], right_term["powers"])),
                "coefficient": left_term["coefficient"] * right_term["coefficient"] * multiplier,
                "origin": _new_origin("product", left_term, right_term,
                                       operation_label=origin),
            })
    return out


def _group_key(term: Mapping[str, Any]) -> tuple[Any, bool, tuple[Any, ...]]:
    info = _provenance(term)
    classes = tuple(sorted((_freeze(value) for value in info["error_classes"]), key=repr))
    return tuple(term["powers"]), bool(info["uncertainty"]), classes


def _aggregate(terms: Sequence[Mapping[str, Any]], operation: str = "aggregate") -> list[dict[str, Any]]:
    """Add directed intervals outward while retaining provenance unions."""

    groups: dict[tuple[Any, bool, tuple[Any, ...]], dict[str, Any]] = {}
    for term in _term_list(terms, "aggregation input"):
        key = _group_key(term)
        info = _provenance(term)
        group = groups.get(key)
        if group is None:
            group = {"powers": tuple(term["powers"]), "coefficient": term["coefficient"],
                     "ids": {}, "classes": {}, "uncertainty": False,
                     "path_count": 0, "source_count": 0}
            groups[key] = group
        else:
            group["coefficient"] = group["coefficient"] + term["coefficient"]
        for value in info["source_ids"]:
            group["ids"][_freeze(value)] = value
        for value in info["error_classes"]:
            group["classes"][_freeze(value)] = value
        group["uncertainty"] = group["uncertainty"] or bool(info["uncertainty"])
        group["path_count"] += int(info["path_count"])
        group["source_count"] += int(info["source_count"])
    return [{
        "powers": group["powers"],
        "coefficient": group["coefficient"],
        "origin": {
            "operation": operation,
            "source_ids": list(group["ids"].values()),
            "error_classes": list(group["classes"].values()),
            "uncertainty": group["uncertainty"],
            "path_count": group["path_count"],
            "source_count": group["source_count"],
            "source_id_count": len(group["ids"]),
        },
    } for group in groups.values()]


def _unit(c: Any, lower: Any) -> Any:
    if lower.a < 0:
        raise ValueError("microscopic exponential bound left positive range")
    return _make_interval(c, lower, c.mpf(1))


def _factors(c: Any, phase: Any) -> tuple[Any, Any, Any, Any]:
    # hb < 1/2, so exp(-alpha*hb*phase) >= 1-alpha*phase/2.
    return (
        _unit(c, 1 - phase / 2),
        _unit(c, 1 - 3 * phase / 4),
        _unit(c, 1 - phase / 4),
        _unit(c, 1 - phase / 4),
    )


def _components(value: Any, label: str) -> tuple[Sequence[Mapping[str, Any]], ...]:
    if isinstance(value, Mapping):
        if set(FIELDS).issubset(value):
            values = [value[name] for name in FIELDS]
        elif set(range(6)).issubset(value):
            values = [value[index] for index in range(6)]
        elif all(str(index) in value for index in range(6)):
            values = [value[str(index)] for index in range(6)]
        else:
            raise ValueError(f"{label} must provide six core components")
    elif isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        values = list(value)
    else:
        raise TypeError(f"{label} must provide six core components")
    if len(values) != 6:
        raise ValueError(f"{label} must provide six core components")
    return tuple(_term_list(component, f"{label}[{index}]")
                 for index, component in enumerate(values))


def _transition(c: Any, phase: Any, initial: Any) -> tuple[list[dict[str, Any]], ...]:
    f, u, k, p, fx, ux = _components(initial, "core initial")
    a, b, e, exprel = _factors(c, phase)
    return (
        _aggregate(list(f) + _scale(fx, phase * exprel, 1, "F_transport_cross"), "core_transport"),
        _aggregate(list(u) + _scale(ux, phase, 1, "U_transport_cross"), "core_transport"),
        _aggregate(_scale(k, a, 0, "K_transport_diagonal")
                   + _scale(ux, -phase * exprel, 1, "K_transport_cross"), "core_transport"),
        _aggregate(list(p), "core_transport"),
        _aggregate(_scale(fx, b, 0, "Fx_transport_diagonal"), "core_transport"),
        _aggregate(_scale(ux, e, 0, "Ux_transport_diagonal"), "core_transport"),
    )


def _partition(c: Any, count: int) -> list[tuple[Any, Any]]:
    if not isinstance(count, int) or isinstance(count, bool) or count < 1:
        raise ValueError("positive cell count required")
    return [(c.mpf(index) / count, c.mpf(index + 1) / count)
            for index in range(count)]


def _order_key(key: Any, label: str) -> int:
    if isinstance(key, bool):
        raise ValueError(f"boolean is not an integer {label}")
    if isinstance(key, int):
        return key
    if isinstance(key, str):
        try:
            return int(key)
        except ValueError as exc:
            raise ValueError(f"invalid integer {label}") from exc
    raise ValueError(f"integer {label} required")


def _order_mapping(raw: Any, orders: Sequence[int], label: str) -> dict[int, Any]:
    if isinstance(raw, Mapping):
        result: dict[int, Any] = {}
        for key, value in raw.items():
            order = _order_key(key, label)
            if order in result:
                raise ValueError(f"duplicate integer {label} {order}")
            result[order] = value
    elif isinstance(raw, Sequence) and not isinstance(raw, (str, bytes, bytearray)):
        result = {index: value for index, value in enumerate(raw)}
    else:
        raise TypeError(f"{label} must be an integer-order mapping")
    expected = set(orders)
    if set(result) != expected:
        raise ValueError(f"{label} orders {sorted(expected)} required")
    return result


def _row_mapping(raw: Any, label: str) -> dict[str, list[Mapping[str, Any]]]:
    if not isinstance(raw, Mapping):
        raise TypeError(f"{label} must provide rows 4, 5, 6")
    result: dict[str, list[Mapping[str, Any]]] = {}
    for key, value in raw.items():
        try:
            row = str(int(key))
        except (TypeError, ValueError) as exc:
            raise ValueError(f"invalid forcing row {key!r}") from exc
        if row in result:
            raise ValueError(f"duplicate forcing row {row}")
        result[row] = _term_list(value, f"{label}[{row}]")
    if set(result) != {"4", "5", "6"}:
        raise ValueError("signed forcing rows 4, 5, 6 required")
    return result


def _forcing_mapping(raw: Any) -> dict[int, dict[str, list[Mapping[str, Any]]]]:
    """Normalize and require all actual ``g_Z0..g_Z3`` receipts."""

    if isinstance(raw, Mapping) and set(raw).issubset({"g_Z0", "g_Z1", "g_Z2", "g_Z3"}):
        raw = {int(str(key)[3:]): value for key, value in raw.items()}
    orders = _order_mapping(raw, range(4), "forcing")
    return {order: _row_mapping(value, f"forcing[{order}]")
            for order, value in orders.items()}


def _matrix_mapping(raw: Any) -> dict[str, dict[int, dict[str, list[Mapping[str, Any]]]]]:
    """Normalize and require all ``C0/B1`` orders ``Z0..Z2``."""

    if isinstance(raw, Mapping) and set(raw).issubset(
            {f"{name}_Z{order}" for name in ("C0", "B1") for order in range(3)}):
        converted: dict[str, dict[int, Any]] = {"C0": {}, "B1": {}}
        for key, value in raw.items():
            name, order_text = str(key).split("_Z")
            converted[name][int(order_text)] = value
        raw = converted
    if not isinstance(raw, Mapping) or set(raw) != {"C0", "B1"}:
        raise ValueError("C0/B1 required")
    result: dict[str, dict[int, dict[str, list[Mapping[str, Any]]]]] = {}
    for name, expected_keys in (("C0", C0_KEYS), ("B1", B1_KEYS)):
        order_map = _order_mapping(raw[name], range(3), f"{name} matrix")
        result[name] = {}
        for order, entries in order_map.items():
            if not isinstance(entries, Mapping) or set(entries) != set(expected_keys):
                raise ValueError(f"raw {name} Z0..Z2 support required")
            result[name][order] = {
                str(key): _term_list(value, f"{name}[{order}][{key}]")
                for key, value in entries.items()
            }
    return result


def _callback(value: Any, label: str):
    if callable(value):
        return value
    if isinstance(value, Mapping):
        return lambda _phase: value
    raise TypeError(f"{label} callback or order mapping required")


def _green(c: Any, phase: Any, callback: Any, count: int) -> dict[str, list[dict[str, Any]]]:
    """Integrate six signed kernels with the formal ``hb`` measure."""

    if phase.b == 0:
        return {name: [] for name in FIELDS}
    callback = _callback(callback, "Green source")
    fields: dict[str, list[Mapping[str, Any]]] = {name: [] for name in FIELDS}
    for index, (left, right) in enumerate(_partition(c, count)):
        t = _make_interval(c, left, right)
        v = phase * t
        p = phase * (1 - t)
        # u/2 = exp(hb*v/2) <= 1/(1-v/4); exp(-hb*p) lies in [1-p/4, 1].
        upper = (1 / (1 - v.b / 4)).b
        uhalf = _make_interval(c, c.mpf(1), upper)
        a, b, e, exprel = _factors(c, p)
        kernels = (
            (uhalf * p * exprel, 1),
            (uhalf * p, 1),
            (-uhalf * p * exprel, 1),
            (c.mpf(1), 0),
            (b, 0),
            (e, 0),
        )
        source = _row_mapping(callback(v), f"forcing cell {index}")
        measure = phase * (right - left) * uhalf
        for field, row, (kernel, kernel_power) in zip(FIELDS, ROWS, kernels):
            fields[field].extend(
                _scale(source[row], measure * kernel, 1 + kernel_power,
                       ("local_green", index, row, field))
            )
    return {field: _aggregate(values, "local_green")
            for field, values in fields.items()}


def _validate_source_atoms(source_atoms: Any) -> None:
    if isinstance(source_atoms, Mapping):
        names = set(source_atoms)
    elif isinstance(source_atoms, Sequence) and not isinstance(source_atoms, (str, bytes, bytearray)):
        names = set()
        for atom in source_atoms:
            if isinstance(atom, Mapping):
                name = atom.get("name", atom.get("atom", atom.get("power")))
                if name is not None:
                    names.add(name)
            else:
                try:
                    names.add(atom)
                except TypeError:
                    continue
    else:
        raise ValueError("four same-source atom records required")
    if names != set(SOURCE_POWER_BASIS):
        raise ValueError("four same-source atom records required")


def _phase_interval(c: Any, value: Any) -> Any:
    return c.mpf(value)


def enclose_low(c: Any, *, sc: Any, phase: Any, core_W0_Z: Any,
                forcing_Z: Any, matrix_Z: Any, inner_cells: int = 1,
                outer_cells: int = 1, source_family: Any = None,
                source_atoms: Any = None,
                geometry_certificate: Any = None) -> dict[str, Any]:
    """Return only local first-word sectors through ordinary ``Z1``.

    The feedback source is assembled as

    ``Z0 = C0 W0 + B1 W0_Z``

    and

    ``Z1 = C0_Z W0 + C0 W0_Z + B1_Z W0_Z + B1 W0_ZZ``.

    Here ``W0`` means the actual core ``W0_Z0..Z2`` homogeneous transport
    plus its actual ``g_Z0..Z2`` Green contribution.  The function emits the
    local Green result for orders 0 and 1 only.  The ``B1 W0_ZZ`` term is
    represented by the supplied order-2 ``W0`` and is therefore fully
    included in the local Z1 feedback; this is not a computation of a Z2
    output word.
    """

    if not hasattr(c, "mpf") or not hasattr(c, "make_mpf"):
        raise TypeError("c must be the caller's MPIntervalContext")
    if forcing_Z is None or matrix_Z is None:
        raise TypeError("actual g_Z0..g_Z3 and C0/B1_Z0..Z2 callbacks are required")
    if source_family is None:
        raise ValueError("source family label required")
    _validate_source_atoms(source_atoms)
    if not isinstance(geometry_certificate, Mapping):
        raise ValueError("caller geometry certificates required")
    required_geometry = (
        "hb_strictly_between_zero_and_half",
        "sc_phase_Z_independent",
        "hb_epsilon_Z_independent",
    )
    if any(geometry_certificate.get(key) is not True for key in required_geometry):
        raise ValueError("caller geometry certificates required")

    sc_iv, phase_iv = c.mpf(sc), _phase_interval(c, phase)
    if sc_iv.a <= 0 or sc_iv.b >= c.mpf(1) / 1024:
        raise ValueError("microscopic sc required")
    if sc_iv.a != sc_iv.b:
        raise ValueError("exact sc required; nondegenerate sc needs a joint phase certificate")
    if phase_iv.a < 0 or phase_iv.b > (sc_iv / 4).a:
        raise ValueError("phase outside actual first collar")

    core_orders = _order_mapping(core_W0_Z, range(4), "core W0")
    core_w0 = {order: _components(value, f"core_W0_Z[{order}]")
               for order, value in core_orders.items()}
    _partition(c, inner_cells)
    _partition(c, outer_cells)
    # Mapping receipts are validated even for a zero-length phase.  A
    # callable receipt remains lazy so a caller can provide a zero-phase
    # callback that is intentionally never evaluated.
    if not callable(forcing_Z):
        _forcing_mapping(forcing_Z)
    if not callable(matrix_Z):
        _matrix_mapping(matrix_Z)
    forcing_callback = _callback(forcing_Z, "forcing")
    matrix_callback = _callback(matrix_Z, "matrix")

    def word0(v: Any) -> dict[int, tuple[list[dict[str, Any]], ...]]:
        out: dict[int, tuple[list[dict[str, Any]], ...]] = {}
        # Low feedback through ordinary Z1 needs W0_Z0, W0_Z1, and W0_Z2.
        for order in range(3):
            initial = _transition(c, v, core_w0[order])

            def body(w: Any, order: int = order) -> dict[str, list[Mapping[str, Any]]]:
                data = _forcing_mapping(forcing_callback(w))
                return data[order]

            local = _green(c, v, body, inner_cells)
            out[order] = tuple(
                _aggregate(initial[index] + local[field], "W0_core_plus_local")
                for index, field in enumerate(FIELDS)
            )
        return out

    def feedback(v: Any, order: int) -> dict[str, list[Mapping[str, Any]]]:
        word = word0(v)
        matrix = _matrix_mapping(matrix_callback(v))
        result: dict[str, list[Mapping[str, Any]]] = {str(row): [] for row in (4, 5, 6)}
        for name, shift in (("C0", 0), ("B1", 1)):
            for matrix_order in range(order + 1):
                word_order = order - matrix_order + shift
                for key, terms in matrix[name][matrix_order].items():
                    row, column = key.split(",")
                    result[row].extend(
                        _product(terms, word[word_order][int(column) - 1], 1,
                                 (name, matrix_order, word_order, key, order))
                    )
        return {row: _aggregate(terms, "feedback_Z" + str(order))
                for row, terms in result.items()}

    local_signed_sectors: dict[int, dict[str, list[dict[str, Any]]]] = {}
    for order in (0, 1):
        local_signed_sectors[order] = _green(
            c, phase_iv, lambda v, order=order: feedback(v, order), outer_cells)

    return {
        "source_family": source_family,
        "source_atoms": source_atoms,
        "caller_geometry_certificate": geometry_certificate,
        "geometry_certificate_independently_verified": False,
        "source_power_basis": list(SOURCE_POWER_BASIS),
        "local_signed_sectors": local_signed_sectors,
        "local_word_orders": [0, 1],
        "core_W1_inlet_included": False,
        "core_W1_inlet_missing": True,
        "core_W1_inlet_status": "not_supplied_and_not_replaced_by_zero",
        "allword_tail_added": False,
        "allword_tail_missing": True,
        "allword_tail_status": "not_added",
        "z2_computed": False,
        "m15_computed": False,
        "core_history_reset": False,
        "source_callback_identity_proved": False,
        "ordinary_Z_fixed_phase_and_R": True,
        "error_classes_preserved": True,
        "uncertainty_preserved": True,
        "aggregation_key": "(powers, uncertainty, error_classes)",
        "aggregation_is_outward_interval_addition_only": True,
        "source_ids_union_retained": True,
        "source_count_retained": True,
        "same_source_cancellation_claimed": False,
        "endpoint_cancellation_claimed": False,
        "hb_materialized": False,
        "formal_hb_measure": True,
        "green_kernel_count": 6,
    }


__all__ = [
    "SOURCE_POWER_BASIS", "FIELDS", "ROWS", "C0_KEYS", "B1_KEYS", "enclose_low",
]
