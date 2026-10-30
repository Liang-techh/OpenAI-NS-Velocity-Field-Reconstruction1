"""Formal first-collar Green transport for signed source sectors.

The public entry point is :func:`enclose`.  It accepts complete, caller-bound
core-exit sectors ``core_W0_Z`` (orders ``0`` through ``3``) and
``core_W1_Z2`` together with callbacks for the actual uniform-collar forcing
``g_Z0..g_Z3`` and the actual ``C0/B1_Z0..Z2`` matrix.  A sector is a mapping
with four integer powers, a directed interval coefficient, and provenance in
``origin`` (or in the equivalent direct metadata fields).

The four powers remain formal: ``(hb, F0_anchor, Pstar_squared, Lambda)`` are
never evaluated.  Products and Green transport therefore preserve the
Laurent/source bookkeeping instead of replacing it with a numerical norm.
Terms are compacted only when powers, ``uncertainty``, and the same set of
``error_classes`` agree.  Compaction is outward directed-interval addition;
it does not assert a common source or a cancellation.  Every compact group
keeps the union of ``source_ids`` and both a path/source count and a unique
source-id count.  Cauchy and historical errors are consequently retained in
the group metadata.

Mathematical limits are intentional.  The caller certifies source identity,
the actual callback enclosures, the geometry of the first collar, and all
tail/history/Cauchy estimates.  This leaf proves none of those facts, adds no
tail, and does not materialize ``hb`` or the other formal atoms.  It only
performs the six signed Green kernels on ``0 <= phase <= sc/4`` and adds the
finite core transport to the separately reported local contribution.
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
    """Return whether *value* looks like a directed interval mpf."""

    return hasattr(value, "_mpi_") and hasattr(value, "a") and hasattr(value, "b")


def _raw_endpoint(value: Any, upper: bool) -> Any:
    """Return one raw endpoint accepted by ``MPIntervalContext.make_mpf``."""

    raw = getattr(value, "_mpi_", None)
    if raw is not None:
        return raw[1 if upper else 0]
    return value


def _make_interval(c: Any, lower: Any, upper: Any) -> Any:
    return c.make_mpf((_raw_endpoint(lower, False), _raw_endpoint(upper, True)))


def _bounds(value: Any) -> tuple[Any, Any]:
    if not _is_interval(value):
        raise TypeError("directed interval required")
    return value.a, value.b


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
    """Make an opaque metadata value usable in a grouping key."""

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


def _combine_infos(infos: Sequence[Mapping[str, Any]], *, path_count: int | None = None,
                   source_count: int | None = None) -> dict[str, Any]:
    ids: list[Any] = []
    classes: list[Any] = []
    uncertain = False
    if path_count is None:
        path_count = 1
    if source_count is None:
        source_count = 1
    for info in infos:
        ids.extend(info.get("source_ids", ()))
        classes.extend(info.get("error_classes", ()))
        uncertain = uncertain or bool(info.get("uncertainty", False))
    return {
        "source_ids": _unique(ids),
        "uncertainty": bool(uncertain),
        "error_classes": _unique(classes),
        "path_count": max(1, int(path_count)),
        "source_count": max(1, int(source_count)),
    }


def _origin_info(value: Any) -> dict[str, Any]:
    """Read shallow provenance from a term/origin without retaining a tree."""

    if isinstance(value, Mapping):
        ids = list(_as_ids(value.get("source_ids", value.get("sources"))))
        classes = list(_as_labels(value.get("error_classes", value.get("uncertainty_class", value.get("error_class")))))
        operation = str(value.get("operation", "")).lower()
        error_kind = str(value.get("error_kind", value.get("kind", ""))).lower()
        uncertain = bool(value.get("uncertainty", value.get("is_uncertainty", False)))
        # A caller may identify a Cauchy/symmetric/history operation rather
        # than repeat the boolean.  Preserve that declaration as uncertainty.
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
            raw_source = max((info.get("source_count", info.get("path_count", 1)) for info in nested), default=raw_path)
        return {
            "source_ids": _unique(ids),
            "uncertainty": bool(uncertain),
            "error_classes": _unique(classes),
            "path_count": _positive_int(raw_path),
            "source_count": _positive_int(raw_source),
            "operation": value.get("operation"),
        }
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        nested = [_origin_info(item) for item in value
                  if isinstance(item, (Mapping, Sequence)) and not isinstance(item, (str, bytes, bytearray))]
        if not nested:
            return {"source_ids": [], "uncertainty": False, "error_classes": [],
                    "path_count": 1, "source_count": 1, "operation": None}
        return _combine_infos(nested,
                              path_count=max(info.get("path_count", 1) for info in nested),
                              source_count=max(info.get("source_count", 1) for info in nested)) | {
            "operation": None
        }
    return {"source_ids": [], "uncertainty": False, "error_classes": [],
            "path_count": 1, "source_count": 1, "operation": None}


def _provenance(term: Any) -> dict[str, Any]:
    """Read term-level and origin-level metadata together."""

    if not isinstance(term, Mapping):
        return _origin_info(term)
    # _origin_info handles direct metadata and the nested ``origin`` field.
    return _origin_info(term)


def _origin_ref(value: Any) -> dict[str, Any]:
    info = _provenance(value)
    return {
        "source_ids": list(info["source_ids"]),
        "uncertainty": bool(info["uncertainty"]),
        "error_classes": list(info["error_classes"]),
        "path_count": int(info["path_count"]),
        "source_count": int(info["source_count"]),
    }


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
    # ``origin`` is the preferred compact provenance container, but direct
    # ``source_ids``/``error_classes``/``uncertainty`` fields are accepted as
    # well so a caller can pass a small sector receipt without wrapping it.
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


def _scale(c: Any, terms: Sequence[Mapping[str, Any]], factor: Any, hb_power: int,
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


def _product(c: Any, left: Sequence[Mapping[str, Any]], right: Sequence[Mapping[str, Any]],
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


def _aggregate(c: Any, terms: Sequence[Mapping[str, Any]], operation: str = "aggregate") -> list[dict[str, Any]]:
    """Outward interval addition; keep a compact, complete input-ID union."""
    groups: dict[tuple[Any, bool, tuple[Any, ...]], dict[str, Any]] = {}
    for term in _term_list(terms, "aggregation input"):
        key = _group_key(term)
        info = _provenance(term)
        group = groups.get(key)
        if group is None:
            group = {"powers": tuple(term["powers"]), "coefficient": term["coefficient"],
                     "ids": {}, "classes": {}, "uncertainty": False, "path_count": 0, "source_count": 0}
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
    return [{"powers": group["powers"], "coefficient": group["coefficient"],
             "origin": {"operation": operation, "source_ids": list(group["ids"].values()),
                        "error_classes": list(group["classes"].values()),
                        "uncertainty": group["uncertainty"], "path_count": group["path_count"],
                        "source_count": group["source_count"], "source_id_count": len(group["ids"])}}
            for group in groups.values()]


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
    return tuple(_term_list(component, f"{label}[{index}]") for index, component in enumerate(values))


def _transition(c: Any, phase: Any, initial: Any) -> tuple[list[dict[str, Any]], ...]:
    f, u, k, p, fx, ux = _components(initial, "core initial")
    a, b, e, exprel = _factors(c, phase)
    return (
        _aggregate(c, list(f) + _scale(c, fx, phase * exprel, 1, "F_transport_cross"), "core_transport"),
        _aggregate(c, list(u) + _scale(c, ux, phase, 1, "U_transport_cross"), "core_transport"),
        _aggregate(c, _scale(c, k, a, 0, "K_transport_diagonal")
                   + _scale(c, ux, -phase * exprel, 1, "K_transport_cross"), "core_transport"),
        _aggregate(c, list(p), "core_transport"),
        _aggregate(c, _scale(c, fx, b, 0, "Fx_transport_diagonal"), "core_transport"),
        _aggregate(c, _scale(c, ux, e, 0, "Ux_transport_diagonal"), "core_transport"),
    )


def _partition(c: Any, count: int) -> list[tuple[Any, Any]]:
    if not isinstance(count, int) or isinstance(count, bool) or count < 1:
        raise ValueError("positive cell count required")
    return [(c.mpf(index) / count, c.mpf(index + 1) / count) for index in range(count)]


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
    """Normalize g_Z0..g_Z3, including the named callback form."""

    if isinstance(raw, Mapping) and set(raw).issubset({"g_Z0", "g_Z1", "g_Z2", "g_Z3"}):
        raw = {int(str(key)[3:]): value for key, value in raw.items()}
    orders = _order_mapping(raw, range(4), "forcing")
    return {order: _row_mapping(value, f"forcing[{order}]") for order, value in orders.items()}


def _matrix_mapping(raw: Any) -> dict[str, dict[int, dict[str, list[Mapping[str, Any]]]]]:
    """Normalize C0/B1 order maps, including C0_Z0/B1_Z0 named records."""

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
    """Integrate the six signed kernels with the formal ``hb`` measure."""

    if phase.b == 0:
        return {name: [] for name in FIELDS}
    callback = _callback(callback, "Green source")
    fields: dict[str, list[Mapping[str, Any]]] = {name: [] for name in FIELDS}
    for index, (left, right) in enumerate(_partition(c, count)):
        t = _make_interval(c, left, right)
        v = phase * t
        p = phase * (1 - t)
        # u/2 = exp(hb*v/2) <= 1/(1-v/4), and exp(-hb*p) is enclosed by
        # [1-p/4, 1].  The signs of K and -K are kept in the kernels.
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
                _scale(c, source[row], measure * kernel, 1 + kernel_power,
                       ("local_green", index, row, field))
            )
    return {field: _aggregate(c, values, "local_green") for field, values in fields.items()}


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


def enclose(c: Any, *, sc: Any, phase: Any, core_W0_Z: Any, core_W1_Z2: Any,
            forcing_Z: Any = None, matrix_Z: Any = None, inner_cells: int = 1,
            outer_cells: int = 1, source_family: Any = None, source_atoms: Any = None,
            geometry_certificate: Any = None, g_Z: Any = None,
            C0_B1_Z: Any = None) -> dict[str, Any]:
    """Return compact formal core/local Green sectors.

    ``forcing_Z`` may be a callback ``phase -> {order: {row: sectors}}`` or a
    constant order map; ``matrix_Z`` has the analogous
    ``{'C0': {order: {row,col: sectors}}, 'B1': ...}`` form.  Integer order
    keys and their string spellings are accepted.  ``g_Z`` and ``C0_B1_Z``
    are aliases for callers that name the physical receipts explicitly.
    """

    if not hasattr(c, "mpf") or not hasattr(c, "make_mpf"):
        raise TypeError("c must be the caller's MPIntervalContext")
    if forcing_Z is None:
        forcing_Z = g_Z
    elif g_Z is not None and forcing_Z is not g_Z:
        raise ValueError("forcing_Z and g_Z disagree")
    if matrix_Z is None:
        matrix_Z = C0_B1_Z
    elif C0_B1_Z is not None and matrix_Z is not C0_B1_Z:
        raise ValueError("matrix_Z and C0_B1_Z disagree")
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
    core_w0 = {order: _components(value, f"core_W0_Z[{order}]") for order, value in core_orders.items()}
    core_w1 = _components(core_W1_Z2, "core_W1_Z2")
    _partition(c, inner_cells)
    _partition(c, outer_cells)

    forcing_callback = _callback(forcing_Z, "forcing")
    matrix_callback = _callback(matrix_Z, "matrix")

    def word0(v: Any) -> dict[int, tuple[list[dict[str, Any]], ...]]:
        out: dict[int, tuple[list[dict[str, Any]], ...]] = {}
        for order in range(4):
            initial = _transition(c, v, core_w0[order])

            def body(w: Any, order: int = order) -> dict[str, list[Mapping[str, Any]]]:
                data = _forcing_mapping(forcing_callback(w))
                return data[order]

            local = _green(c, v, body, inner_cells)
            out[order] = tuple(
                _aggregate(c, initial[index] + local[field], "W0_core_plus_local")
                for index, field in enumerate(FIELDS)
            )
        return out

    def feedback(v: Any) -> dict[str, list[Mapping[str, Any]]]:
        word = word0(v)
        matrix = _matrix_mapping(matrix_callback(v))
        result: dict[str, list[Mapping[str, Any]]] = {str(row): [] for row in (4, 5, 6)}
        specs = (
            ("C0", 2, 0, 1), ("C0", 1, 1, 2), ("C0", 0, 2, 1),
            ("B1", 2, 1, 1), ("B1", 1, 2, 2), ("B1", 0, 3, 1),
        )
        for name, matrix_order, word_order, factor in specs:
            for key, terms in matrix[name][matrix_order].items():
                row, column = key.split(",")
                result[row].extend(
                    _product(c, terms, word[word_order][int(column) - 1], factor,
                             (name, matrix_order, word_order, key))
                )
        return {row: _aggregate(c, terms, "feedback_row") for row, terms in result.items()}

    core = _transition(c, phase_iv, core_w1)
    local = _green(c, phase_iv, feedback, outer_cells)
    finite = {
        field: _aggregate(c, core[index] + local[field], "finite_W1_Z2")
        for index, field in enumerate(FIELDS)
    }
    return {
        "source_family": source_family,
        "source_atoms": source_atoms,
        "caller_geometry_certificate": geometry_certificate,
        "geometry_certificate_independently_verified": False,
        "source_power_basis": list(SOURCE_POWER_BASIS),
        "core_transport": dict(zip(FIELDS, core)),
        "core_W0_transport": {
            order: dict(zip(FIELDS, _transition(c, phase_iv, core_w0[order])))
            for order in range(4)
        },
        "local_signed_sectors": local,
        "finite_W1_Z2": finite,
        "source_callback_identity_proved": False,
        "ordinary_Z_fixed_phase_and_R": True,
        "allword_tail_added": False,
        "m15_computed": False,
        "core_history_reset": False,
        "error_classes_preserved": True,
        "uncertainty_preserved": True,
        "aggregation_key": "(powers, uncertainty, error_classes)",
        "aggregation_is_outward_interval_addition_only": True,
        "source_ids_union_retained": True,
        "source_count_retained": True,
        "same_source_cancellation_claimed": False,
        "hb_materialized": False,
        "formal_hb_measure": True,
        "green_kernel_count": 6,
    }


__all__ = [
    "SOURCE_POWER_BASIS", "FIELDS", "ROWS", "C0_KEYS", "B1_KEYS", "enclose",
]
