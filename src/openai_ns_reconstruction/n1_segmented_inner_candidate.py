"""Conditional segmented inner finite candidate.

This is a small, source-only consumer.  The caller supplies two already
authenticated, already physically weighted packets: the core and the first
collar.  Each packet contains ordinary ``Z0``/``Z1`` sector records and all
coefficients are directed intervals created by one caller-owned
``MPIntervalContext``.  The leaf only groups and adds those records.  It does
not evaluate a source, attach a physical factor, add a pressure endpoint,
construct a Q/P packet, select a point value, or write a ledger.

``segment_groups`` groups records separately by the exact key
``(powers, uncertainty, ordered error_classes)``.  ``formal_sum`` then adds
groups having the same key across the two segments.  Every group and every
formal-sum item retains original record and provenance references; a shared
key is a bookkeeping key, not a source-equality or endpoint-cancellation
claim.  Full original core ``W1`` and the ``Kchi-K0``, feedback, and
all-word-tail pieces remain explicitly incomplete, while a separately saved
``K0``/``W1`` subterm may be supplied as its own word.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any


SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
ORDINARY_Z_ORDERS = (0, 1)
SEGMENTS = ("core", "collar")
SCHEMA = "n1_segmented_inner_candidate_v1"
REQUIRED_RECORD_FIELDS = frozenset(
    (
        "powers",
        "coefficient",
        "segment",
        "word",
        "channel",
        "cell",
        "receipt",
        "path",
        "origin",
    )
)
GROUPING_KEY = "(powers, uncertainty, ordered error_classes)"

__all__ = [
    "SOURCE_POWER_BASIS",
    "ORDINARY_Z_ORDERS",
    "SEGMENTS",
    "SCHEMA",
    "GROUPING_KEY",
    "build",
    "assemble",
    "combine",
    "construct",
    "enclose",
]


def _is_sequence(value: Any) -> bool:
    return isinstance(value, Sequence) and not isinstance(
        value, (str, bytes, bytearray)
    )


def _reject_float(value: Any, label: str, *, allow_bool: bool = False) -> None:
    """Reject host floats in structural/provenance input.

    Coefficients are checked separately and are intentionally not traversed:
    an interval object is an opaque value owned by the caller's context.
    ``uncertainty`` is the one provenance boolean accepted by this leaf.
    """

    if isinstance(value, float):
        raise TypeError(f"{label} must not contain floating-point values")
    if isinstance(value, bool) and not allow_bool:
        raise TypeError(f"{label} must not contain boolean values")
    if isinstance(value, Mapping):
        for key, item in value.items():
            _reject_float(key, f"{label}.key", allow_bool=allow_bool)
            _reject_float(item, label, allow_bool=allow_bool)
    elif _is_sequence(value):
        for index, item in enumerate(value):
            _reject_float(item, f"{label}[{index}]", allow_bool=allow_bool)


def _typed_key(value: Any, label: str) -> tuple[type[Any], Any]:
    """Return a type-sensitive hash key without coercing provenance."""

    if isinstance(value, (bool, float)):
        raise TypeError(f"{label} must not be bool or float")
    try:
        hash(value)
    except TypeError as exc:
        raise TypeError(f"{label} must be hashable") from exc
    return type(value), value


def _is_directed_interval(value: Any) -> bool:
    """Recognize an existing mpmath interval without importing mpmath."""

    return (
        not isinstance(value, (bool, int, float, complex, str, bytes, bytearray))
        and hasattr(value, "_mpi_")
        and hasattr(value, "a")
        and hasattr(value, "b")
    )


def _check_context(context: Any, coefficient: Any, label: str) -> None:
    """Ensure an interval can operate with the caller's context."""

    if not _is_directed_interval(coefficient):
        raise TypeError(f"{label} must be an existing directed interval")
    try:
        # This compatibility check does not coerce the coefficient or inspect
        # a midpoint.  mpmath exposes no stable public owner attribute on all
        # interval versions.
        _ = coefficient.a + context.mpf(0)
    except (AttributeError, TypeError, ValueError, ArithmeticError) as exc:
        raise TypeError(
            f"{label} is not compatible with the caller's interval context"
        ) from exc


def _powers(value: Any, label: str) -> tuple[int, int, int, int]:
    if not _is_sequence(value) or len(value) != 4:
        raise ValueError(f"{label} must contain exactly four integer powers")
    result: list[int] = []
    for index, power in enumerate(value):
        if isinstance(power, bool) or not isinstance(power, int):
            raise TypeError(f"{label}[{index}] must be an integer, not bool or float")
        result.append(int(power))
    return tuple(result)  # type: ignore[return-value]


def _order(raw: Any, label: str) -> int:
    if isinstance(raw, bool) or isinstance(raw, float):
        raise TypeError(f"{label} must be ordinary-Z integer 0 or 1")
    if isinstance(raw, int):
        result = raw
    elif isinstance(raw, str):
        text = raw.strip()
        if text.lower() in {"z0", "z1"}:
            result = int(text[-1])
        elif text in {"0", "1"}:
            result = int(text)
        else:
            raise ValueError(f"{label} must be ordinary-Z order 0 or 1")
    else:
        raise TypeError(f"{label} must be ordinary-Z integer 0 or 1")
    if result not in ORDINARY_Z_ORDERS:
        raise ValueError(f"{label} must be ordinary-Z order 0 or 1")
    return result


def _rows(value: Any, label: str) -> dict[int, list[Any]]:
    """Normalize a packet to two order-major sector lists."""

    if isinstance(value, Mapping):
        normalized: dict[int, Any] = {}
        for raw_order, row in value.items():
            parsed = _order(raw_order, f"{label} order")
            if parsed in normalized:
                raise ValueError(f"{label} contains duplicate ordinary-Z order {parsed}")
            normalized[parsed] = row
        if set(normalized) != set(ORDINARY_Z_ORDERS):
            raise ValueError(f"{label} must provide exactly ordinary-Z orders 0 and 1")
    elif _is_sequence(value):
        if len(value) != 2:
            raise ValueError(f"{label} must contain exactly Z0 and Z1 rows")
        normalized = {order: value[order] for order in ORDINARY_Z_ORDERS}
    else:
        raise TypeError(f"{label} must be an order mapping or a two-row sequence")

    result: dict[int, list[Any]] = {}
    for order in ORDINARY_Z_ORDERS:
        row = normalized[order]
        if isinstance(row, Mapping) and "powers" in row:
            row = [row]
        if not _is_sequence(row):
            raise TypeError(f"{label}[Z{order}] must be a sector record list")
        result[order] = list(row)
    return result


def _prepared_record(
    raw: Any,
    *,
    segment: str,
    word: str,
    channel: str,
    cell: Any,
    receipt: str,
    path: Any,
    origin: Any = None,
) -> Any:
    """Add the small leaf envelope without mutating a caller record."""

    if not isinstance(raw, Mapping):
        raise TypeError("sector record must be a mapping")
    if REQUIRED_RECORD_FIELDS.issubset(raw):
        return raw
    view = dict(raw)
    view.setdefault("segment", segment)
    view.setdefault("word", word)
    view.setdefault("channel", channel)
    view.setdefault("cell", cell)
    view.setdefault("receipt", receipt)
    view.setdefault("path", path)
    if "origin" not in view and origin is not None:
        view["origin"] = origin
    # ``_original_record`` is internal only; output members expose the raw
    # object through ``record`` and never claim that this view is a source.
    view["_original_record"] = raw
    return view


def _original_record(value: Mapping[str, Any]) -> Any:
    return value.get("_original_record", value)


def _core_packet_rows(value: Any) -> tuple[dict[int, list[Any]], Any, Any]:
    """Flatten the saved core hierarchy or accept a direct order packet."""

    if not isinstance(value, Mapping) or not any(
        key in value for key in ("axial", "swirl", "pressure", "nested_w1")
    ):
        return _rows(value, "core"), None, None

    result: dict[int, list[Any]] = {0: [], 1: []}
    for channel in ("axial", "swirl"):
        channel_rows = value.get(channel, {})
        if not isinstance(channel_rows, Mapping):
            raise TypeError(f"core.{channel} must be an order/cell mapping")
        for raw_order, cells in channel_rows.items():
            order = _order(raw_order, f"core.{channel} order")
            if not isinstance(cells, Mapping):
                raise TypeError(f"core.{channel}[Z{order}] must be a cell mapping")
            for raw_cell, terms in cells.items():
                cell = raw_cell
                if isinstance(cell, (bool, float)):
                    raise TypeError("core cell must not be bool or float")
                if not _is_sequence(terms):
                    raise TypeError(f"core.{channel}[Z{order}][{cell!r}] must be a sector list")
                for index, raw in enumerate(terms):
                    result[order].append(
                        _prepared_record(
                            raw,
                            segment="core",
                            word="W0",
                            channel=channel,
                            cell=cell,
                            receipt="core",
                            path=["core", channel, order, cell, index],
                        )
                    )

    pressure_rows = value.get("pressure", {})
    if not isinstance(pressure_rows, Mapping):
        raise TypeError("core.pressure must be an order mapping")
    for raw_order, terms in pressure_rows.items():
        order = _order(raw_order, "core.pressure order")
        if not _is_sequence(terms):
            raise TypeError(f"core.pressure[Z{order}] must be a sector list")
        for index, raw in enumerate(terms):
            result[order].append(
                _prepared_record(
                    raw,
                    segment="core",
                    word="W0",
                    channel="pressure",
                    cell=None,
                    receipt="core",
                    path=["core", "pressure", order, index],
                )
            )

    nested = value.get("nested_w1", {})
    if nested is not None:
        if not isinstance(nested, Mapping):
            raise TypeError("core.nested_w1 must be a channel mapping")
        for channel in ("axial", "swirl"):
            channel_rows = nested.get(channel, {})
            if not isinstance(channel_rows, Mapping):
                raise TypeError(f"core.nested_w1.{channel} must be an order/cell mapping")
            for raw_order, cells in channel_rows.items():
                order = _order(raw_order, f"core.nested_w1.{channel} order")
                if not isinstance(cells, Mapping):
                    raise TypeError(f"core.nested_w1.{channel}[Z{order}] must be a cell mapping")
                for raw_cell, terms in cells.items():
                    if isinstance(raw_cell, (bool, float)):
                        raise TypeError("core nested W1 cell must not be bool or float")
                    if not _is_sequence(terms):
                        raise TypeError("core nested W1 cell must be a sector list")
                    for index, raw in enumerate(terms):
                        result[order].append(
                            _prepared_record(
                                raw,
                                segment="core",
                                word="W1-K0",
                                channel=channel,
                                cell=raw_cell,
                                receipt="core_W1_K0",
                                path=["core", "nested_w1", channel, order, raw_cell, index],
                            )
                        )
    return result, value.get("cells"), value.get("source_index", value.get("input_index", {}))


def _collar_packet_rows(value: Any) -> tuple[dict[int, list[Any]], Any, Any]:
    """Flatten prior collar groups or accept a direct order packet."""

    if not isinstance(value, Mapping) or "groups" not in value:
        return _rows(value, "collar"), None, None
    groups = value["groups"]
    if not isinstance(groups, Mapping):
        raise TypeError("collar.groups must be an order mapping")
    result: dict[int, list[Any]] = {0: [], 1: []}
    for raw_order, row in groups.items():
        order = _order(raw_order, "collar.groups order")
        if not _is_sequence(row):
            raise TypeError(f"collar.groups[Z{order}] must be a group list")
        for group_index, group in enumerate(row):
            if not isinstance(group, Mapping):
                raise TypeError("collar group must be a mapping")
            members = group.get("members", ())
            if not _is_sequence(members):
                raise TypeError("collar group members must be a list")
            for member_index, member in enumerate(members):
                if not isinstance(member, Mapping):
                    raise TypeError("collar group member must be a mapping")
                raw = member.get("record", member)
                if not isinstance(raw, Mapping):
                    raise TypeError("collar member record must be a mapping")
                result[order].append(
                    _prepared_record(
                        raw,
                        segment="collar",
                        word=str(member.get("word", raw.get("word", "collar-prior"))),
                        channel=str(member.get("channel", raw.get("channel", "unknown"))),
                        cell=member.get("cell", raw.get("cell")),
                        receipt=str(member.get("receipt", raw.get("receipt", "collar"))),
                        path=member.get(
                            "path",
                            raw.get(
                                "path",
                                ["collar", "groups", order, group_index, "members", member_index, "record"],
                            ),
                        ),
                        origin=member.get("provenance", group.get("origin")),
                    )
                )
    return result, groups, value.get("source_index", value.get("input_index", {}))


def _source_index_from_rows(rows: Mapping[int, Sequence[Any]]) -> dict[Any, dict[str, Any]]:
    result: dict[Any, dict[str, Any]] = {}
    for order in ORDINARY_Z_ORDERS:
        for index, record in enumerate(rows[order]):
            origin = record.get("origin", {})
            source_ids = origin.get("source_ids", ()) if isinstance(origin, Mapping) else ()
            for source_id in source_ids:
                result.setdefault(source_id, {"source_id": source_id, "occurrences": []})["occurrences"].append({
                    "source_id": source_id,
                    "receipt": record.get("receipt"),
                    "path": record.get("path"),
                    "ordinary_Z_order": order,
                    "record_index": index,
                })
    return result


def _labels(value: Any, label: str) -> None:
    if not isinstance(value, str) or not value:
        raise TypeError(f"{label} must be a non-empty string")
    _reject_float(value, label)


def _origin_data(origin: Any, label: str) -> tuple[bool, list[Any], list[Any]]:
    if not isinstance(origin, Mapping):
        raise TypeError(f"{label} must be a provenance mapping")
    missing = {"uncertainty", "error_classes"}.difference(origin)
    if missing:
        raise ValueError(f"{label} is missing {', '.join(sorted(missing))}")

    uncertainty = origin["uncertainty"]
    if not isinstance(uncertainty, bool):
        raise TypeError(f"{label}.uncertainty must be bool")
    classes = origin["error_classes"]
    if not _is_sequence(classes):
        raise TypeError(f"{label}.error_classes must be an ordered sequence")
    classes_copy = list(classes)
    for index, item in enumerate(classes_copy):
        _typed_key(item, f"{label}.error_classes[{index}]")

    ids_raw = origin.get("source_ids", ())
    if isinstance(ids_raw, (str, bytes, bytearray)) or not _is_sequence(ids_raw):
        raise TypeError(f"{label}.source_ids must be a sequence")
    ids_copy = list(ids_raw)
    if not ids_copy or any(not isinstance(item, str) or not item for item in ids_copy):
        raise ValueError(f"{label}.source_ids must contain nonempty strings")
    for index, item in enumerate(ids_copy):
        _typed_key(item, f"{label}.source_ids[{index}]")

    # Preserve arbitrary provenance fields, including boolean status flags;
    # host floating-point values remain inadmissible.
    for key, item in origin.items():
        if key == "uncertainty":
            continue
        _reject_float(item, f"{label}.{key}", allow_bool=True)
    return bool(uncertainty), classes_copy, ids_copy


def _reject_forbidden_labels(record: Mapping[str, Any], label: str) -> None:
    for key in ("word", "channel"):
        value = record[key]
        if isinstance(value, str) and value.strip().lower() == "p1":
            raise ValueError(f"{label}.{key} P1 is not an m15 sector")
        # A direct Q/P packet is outside this leaf.  Path/receipt text remains
        # opaque provenance and is deliberately not inspected.
        if isinstance(value, str) and value.strip().lower() in {"q", "p", "q/p", "qp"}:
            raise ValueError(f"{label}.{key} independent Q/P input is not accepted")


def _record(
    context: Any,
    raw: Any,
    *,
    segment: str,
    order: int,
    index: int,
) -> tuple[tuple[int, int, int, int], bool, list[Any], list[Any]]:
    label = f"{segment}[Z{order}][{index}]"
    if not isinstance(raw, Mapping):
        raise TypeError(f"{label} must be a sector record mapping")
    missing = REQUIRED_RECORD_FIELDS.difference(raw)
    if missing:
        raise ValueError(f"{label} is missing {', '.join(sorted(missing))}")

    _labels(raw["segment"], f"{label}.segment")
    if raw["segment"] != segment:
        raise ValueError(f"{label}.segment must be {segment!r}")
    _labels(raw["word"], f"{label}.word")
    _labels(raw["channel"], f"{label}.channel")
    _labels(raw["receipt"], f"{label}.receipt")

    cell = raw["cell"]
    if isinstance(cell, (bool, float)):
        raise TypeError(f"{label}.cell must not be bool or float")
    _reject_float(cell, f"{label}.cell")
    path = raw["path"]
    if not isinstance(path, (str, bytes, bytearray)) and not _is_sequence(path):
        raise TypeError(f"{label}.path must be a string or sequence")
    _reject_float(path, f"{label}.path")

    powers = _powers(raw["powers"], f"{label}.powers")
    _check_context(context, raw["coefficient"], f"{label}.coefficient")
    uncertainty, classes, source_ids = _origin_data(raw["origin"], f"{label}.origin")

    for optional_order in ("ordinary_Z_order", "order"):
        if optional_order in raw and _order(raw[optional_order], f"{label}.{optional_order}") != order:
            raise ValueError(f"{label}.{optional_order} disagrees with its row")
    _reject_forbidden_labels(raw, label)
    word = raw["word"].strip().lower()
    channel = raw["channel"].strip().lower()
    # A separately saved core K0/W1 packet is allowed, but W1 must never
    # carry a pressure record in this candidate.
    if segment == "core" and "w1" in word and ("pressure" in channel or "n1p" in word):
        raise ValueError(f"{label} core W1 pressure is not part of the candidate")
    return powers, uncertainty, classes, source_ids


def _key(
    powers: tuple[int, int, int, int],
    uncertainty: bool,
    classes: Sequence[Any],
) -> tuple[Any, bool, tuple[tuple[type[Any], Any], ...]]:
    return (
        powers,
        uncertainty,
        tuple(_typed_key(item, "error class") for item in classes),
    )


def _group_segment(
    context: Any,
    segment: str,
    rows: Mapping[int, Sequence[Any]],
) -> dict[int, list[dict[str, Any]]]:
    """Group each ordinary order while retaining original member objects."""

    grouped: dict[int, list[dict[str, Any]]] = {0: [], 1: []}
    lookup: dict[int, dict[Any, dict[str, Any]]] = {0: {}, 1: {}}

    for order in ORDINARY_Z_ORDERS:
        for index, record in enumerate(rows[order]):
            powers, uncertainty, classes, source_ids = _record(
                context,
                record,
                segment=segment,
                order=order,
                index=index,
            )
            key = _key(powers, uncertainty, classes)
            group = lookup[order].get(key)
            if group is None:
                group = {
                    "segment": segment,
                    "ordinary_Z_order": order,
                    "powers": powers,
                    "coefficient": record["coefficient"],
                    "origin": {
                        "operation": "within_segment_outward_interval_group",
                        "segment": segment,
                        "ordinary_Z_order": order,
                        "source_ids": [],
                        "uncertainty": uncertainty,
                        "error_classes": list(classes),
                        "member_count": 0,
                        "source_ids_by_member": [],
                        "error_classes_by_member": [],
                        "provenance_retained": True,
                        "source_equality_proved": False,
                    },
                    "members": [],
                    "original_records": [],
                    "original_origins": [],
                    "group_key": {
                        "powers": powers,
                        "uncertainty": uncertainty,
                        "error_classes": list(classes),
                    },
                }
                lookup[order][key] = group
                grouped[order].append(group)
            else:
                # The only arithmetic in this leaf is directed interval
                # addition of records sharing the exact grouping key.
                group["coefficient"] = group["coefficient"] + record["coefficient"]

            origin = record["origin"]
            original = _original_record(record)
            member = {
                "segment": record["segment"],
                "word": record["word"],
                "channel": record["channel"],
                "cell": record["cell"],
                "receipt": record["receipt"],
                "path": record["path"],
                "ordinary_Z_order": order,
                "index": index,
                "record_index": index,
                "record": original,
                "origin": origin,
                "provenance": origin,
            }
            group["members"].append(member)
            group["original_records"].append(original)
            group["original_origins"].append(origin)
            group["origin"]["member_count"] += 1
            group["origin"]["source_ids_by_member"].append(list(source_ids))
            group["origin"]["error_classes_by_member"].append(list(classes))
            group["origin"]["source_ids"].extend(source_ids)
    return grouped


def _group_key(group: Mapping[str, Any]) -> tuple[Any, bool, tuple[tuple[type[Any], Any], ...]]:
    marker = group.get("group_key")
    if not isinstance(marker, Mapping):
        raise ValueError("internal group key missing")
    return _key(
        tuple(marker["powers"]),
        bool(marker["uncertainty"]),
        marker["error_classes"],
    )


def _formal_sum(
    grouped: Mapping[str, Mapping[int, Sequence[Mapping[str, Any]]]],
) -> dict[int, list[dict[str, Any]]]:
    """Add matching core/collar groups by key and retain both group refs."""

    result: dict[int, list[dict[str, Any]]] = {0: [], 1: []}
    for order in ORDINARY_Z_ORDERS:
        source_groups: dict[Any, dict[str, Mapping[str, Any] | None]] = {}
        ordered_keys: list[Any] = []
        for segment in SEGMENTS:
            for group in grouped[segment][order]:
                key = _group_key(group)
                if key not in source_groups:
                    source_groups[key] = {"core": None, "collar": None}
                    ordered_keys.append(key)
                source_groups[key][segment] = group

        for key in ordered_keys:
            refs = source_groups[key]
            groups = [
                refs[segment]
                for segment in SEGMENTS
                if refs[segment] is not None
            ]
            first = groups[0]
            coefficient = first["coefficient"]
            if len(groups) == 2:
                coefficient = coefficient + groups[1]["coefficient"]

            source_ids_by_segment = {
                segment: list(refs[segment]["origin"].get("source_ids", ()))
                if refs[segment] is not None
                else []
                for segment in SEGMENTS
            }
            classes_by_segment = {
                segment: list(refs[segment]["origin"].get("error_classes", ()))
                if refs[segment] is not None
                else []
                for segment in SEGMENTS
            }
            uncertainty_by_segment = {
                segment: refs[segment]["origin"].get("uncertainty")
                if refs[segment] is not None
                else None
                for segment in SEGMENTS
            }
            all_members = [
                member
                for segment in SEGMENTS
                if refs[segment] is not None
                for member in refs[segment]["members"]
            ]
            parent_origins = {
                segment: refs[segment]["origin"]
                for segment in SEGMENTS
                if refs[segment] is not None
            }
            present_segments = [
                segment for segment in SEGMENTS if refs[segment] is not None
            ]

            term = {
                "segment": "segmented",
                "segments": present_segments,
                "ordinary_Z_order": order,
                "powers": tuple(first["powers"]),
                "coefficient": coefficient,
                "word_by_segment": {
                    segment: [member["word"] for member in refs[segment]["members"]]
                    if refs[segment] is not None
                    else []
                    for segment in SEGMENTS
                },
                "channel_by_segment": {
                    segment: [member["channel"] for member in refs[segment]["members"]]
                    if refs[segment] is not None
                    else []
                    for segment in SEGMENTS
                },
                "cell_by_segment": {
                    segment: [member["cell"] for member in refs[segment]["members"]]
                    if refs[segment] is not None
                    else []
                    for segment in SEGMENTS
                },
                "receipt_by_segment": {
                    segment: [member["receipt"] for member in refs[segment]["members"]]
                    if refs[segment] is not None
                    else []
                    for segment in SEGMENTS
                },
                "path_by_segment": {
                    segment: [member["path"] for member in refs[segment]["members"]]
                    if refs[segment] is not None
                    else []
                    for segment in SEGMENTS
                },
                # These are references to the actual grouped dictionaries;
                # they are not copied summaries and are never used to infer
                # source equality.
                "group_refs": {segment: refs[segment] for segment in SEGMENTS},
                "source_group_refs": {segment: refs[segment] for segment in SEGMENTS},
                "source_groups": tuple(groups),
                "members": tuple(all_members),
                "original_records": tuple(member["record"] for member in all_members),
                "original_origins": tuple(member["origin"] for member in all_members),
                "origin": {
                    "operation": "cross_segment_outward_interval_sum",
                    "segments": present_segments,
                    "ordinary_Z_order": order,
                    "source_ids": [
                        source_id
                        for segment in SEGMENTS
                        for source_id in source_ids_by_segment[segment]
                    ],
                    "source_ids_by_segment": source_ids_by_segment,
                    "uncertainty": bool(first["origin"]["uncertainty"]),
                    "uncertainty_by_segment": uncertainty_by_segment,
                    "error_classes": list(first["origin"]["error_classes"]),
                    "error_classes_by_segment": classes_by_segment,
                    "parent_origins": parent_origins,
                    "source_group_count": len(groups),
                    "member_count": len(all_members),
                    "provenance_retained": True,
                    "source_equality_proved": False,
                    "endpoint_cancellation_claimed": False,
                },
                "group_key": {
                    "powers": tuple(first["powers"]),
                    "uncertainty": bool(first["origin"]["uncertainty"]),
                    "error_classes": list(first["origin"]["error_classes"]),
                },
            }
            result[order].append(term)
    return result


def _word_flags(rows: Mapping[int, Sequence[Any]]) -> dict[str, bool]:
    words = [
        record.get("word", "")
        for order in ORDINARY_Z_ORDERS
        for record in rows[order]
        if isinstance(record, Mapping)
    ]
    normalized = [word.strip().lower() for word in words if isinstance(word, str)]
    has_w0 = any("w0" in word for word in normalized)
    has_w1 = any("w1" in word for word in normalized)
    has_k0_w1 = any("w1" in word and "k0" in word for word in normalized)
    has_n1p = any("n1p" in word for word in normalized) or any(
        isinstance(record, Mapping) and record.get("channel") == "pressure"
        for order in ORDINARY_Z_ORDERS
        for record in rows[order]
    )
    return {
        "finite_W0_present": has_w0,
        "W1_present": has_w1,
        "K0_W1_present": has_k0_w1,
        "N1p_pressure_present": has_n1p,
    }


def _scope(
    core_rows: Mapping[int, Sequence[Any]],
    collar_rows: Mapping[int, Sequence[Any]],
) -> dict[str, Any]:
    core = _word_flags(core_rows)
    collar = _word_flags(collar_rows)
    return {
        "core": {
            "finite_W0": core["finite_W0_present"],
            "finite_W0_present": core["finite_W0_present"],
            "N1p_pressure": core["N1p_pressure_present"],
            "N1p_pressure_present": core["N1p_pressure_present"],
            "K0_W1_separately_saved": core["K0_W1_present"],
            "K0_W1_in_candidate": core["K0_W1_present"],
            "K0_W1_not_yet_consumed": not core["K0_W1_present"],
            "full_original_W1_missing": True,
            "full_original_core_W1_missing": True,
            "core_W1_radial_m15_complete": False,
            "Kchi_minus_K0_missing": True,
            "feedback_missing": True,
            "allword_tail_missing": True,
            "tail_missing": True,
        },
        "collar": {
            "finite_W0": None,
            "finite_W0_present": None,
            "W1": None,
            "W1_present": None,
            "N1p_pressure": None,
            "N1p_pressure_present": None,
            "prior_group_word_contents_require_caller_provenance": True,
            "pressure_endpoint_added": False,
            "independent_QP_added": False,
        },
    }


def build(
    context: Any,
    core: Any = None,
    collar: Any = None,
    *,
    segments: Any = None,
    source_family: Any = None,
    source_atoms: Any = None,
    source_logs: Any = None,
    geometry_certificate: Any = None,
    owner_certificate: Any = None,
    scale_certificate: Any = None,
) -> dict[str, Any]:
    """Build the conditional segmented candidate from two caller packets.

    ``core`` and ``collar`` are either ``{0: records, 1: records}`` (string
    order spellings are accepted) or ``[Z0_records, Z1_records]``.  Each
    record must contain the nine fields listed in
    :data:`REQUIRED_RECORD_FIELDS`.  ``segments={"core": ..., "collar": ...}``
    is accepted as a keyword alias.  No input is copied or numerically
    re-evaluated; member records and origins remain the original caller
    objects.
    """

    if not hasattr(context, "mpf") or not hasattr(context, "make_mpf"):
        raise TypeError("context must be the caller's MPIntervalContext")

    if segments is not None:
        if not isinstance(segments, Mapping) or set(segments) != set(SEGMENTS):
            raise ValueError("segments must provide exactly core and collar")
        if core is not None and core is not segments["core"]:
            raise ValueError("core and segments['core'] disagree")
        if collar is not None and collar is not segments["collar"]:
            raise ValueError("collar and segments['collar'] disagree")
        core, collar = segments["core"], segments["collar"]
    if core is None or collar is None:
        raise TypeError("core and collar ordinary-Z sector packets are required")

    # Optional caller metadata is retained as supplied.  Floats are rejected,
    # while boolean certificate flags remain legal metadata.
    for value, label in (
        (source_family, "source_family"),
        (source_atoms, "source_atoms"),
        (source_logs, "source_logs"),
        (geometry_certificate, "geometry_certificate"),
        (owner_certificate, "owner_certificate"),
        (scale_certificate, "scale_certificate"),
    ):
        if value is not None:
            _reject_float(value, label, allow_bool=True)

    core_rows, core_cells, core_source_index = _core_packet_rows(core)
    collar_rows, collar_groups, collar_source_index = _collar_packet_rows(collar)
    grouped = {
        "core": _group_segment(context, "core", core_rows),
        "collar": _group_segment(context, "collar", collar_rows),
    }
    formal_sum = _formal_sum(grouped)

    source_index = {
        "core": core_source_index if isinstance(core_source_index, Mapping) and core_source_index else _source_index_from_rows(core_rows),
        "collar": collar_source_index if isinstance(collar_source_index, Mapping) and collar_source_index else _source_index_from_rows(collar_rows),
    }

    segments_output = {
        "core": grouped["core"],
        "collar": grouped["collar"],
        "C0": formal_sum[0],
        "Z1": formal_sum[1],
    }

    # The aliases point at the same in-memory records/lists.  This lets a
    # caller use either the segment-oriented or sum-oriented spelling without
    # creating a second arithmetic result.
    return {
        "schema": SCHEMA,
        "source_family": source_family,
        "source_logs": source_logs,
        "source_atoms": source_atoms,
        "source_power_basis": list(SOURCE_POWER_BASIS),
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "segments": segments_output,
        "segment_names": list(SEGMENTS),
        "grouping_key": GROUPING_KEY,
        "segment_groups": grouped,
        "groups_by_segment": grouped,
        "groups": grouped,
        "formal_sum": formal_sum,
        "segmented_candidate": formal_sum,
        "candidate": formal_sum,
        "scope": _scope(core_rows, collar_rows),
        "core_cells": core_cells if core_cells is not None else core_rows,
        "collar_groups": collar_groups if collar_groups is not None else collar_rows,
        "source_index": source_index,
        "input_index": source_index,
        "core_input": core,
        "collar_input": collar,
        "geometry_certificate": geometry_certificate,
        "owner_certificate": owner_certificate,
        "scale_certificate": scale_certificate,
        "records_already_physically_weighted": None,
        "physical_weighting_is_caller_obligation": True,
        "additional_physical_factor_applied": False,
        "physical_factors_applied": False,
        "pressure_endpoint_added": False,
        "independent_QP_added": False,
        "P1_mixed": None,
        "all_member_error_classes_retained": True,
        "provenance_retained": True,
        "same_interval_context_validated": False,
        "same_interval_context_is_caller_obligation": True,
        "source_equality_proved": False,
        "source_equality_from_intervals": False,
        "endpoint_cancellation_claimed": False,
        "same_source_cancellation_claimed": False,
        "point_values_selected": False,
        "real_point_values_computed": False,
        "midpoint_values_used": False,
        "ledger_appended": False,
        "writes_ledger": False,
        "no_ledger": True,
        "conditional_owner_scale": True,
        "owner_identity_proved": False,
        "physical_scale_identity_proved": False,
        "same_measure_proved": False,
        "candidate_is_conditional": True,
        "finite_W0_plus_N1p_core": None,
        "finite_W0_W1_plus_N1p_collar": None,
        "nested_K0_W1_saved_not_full_W1": True,
        "rho_seam": (
            geometry_certificate.get("rho_seam", False)
            if isinstance(geometry_certificate, Mapping)
            else False
        ),
        "geometry_seam": (
            geometry_certificate.get("rho_seam", False)
            if isinstance(geometry_certificate, Mapping)
            else False
        ),
        "P1_is_moment": False,
        "pressure_is_moment": False,
        "physical_scale_reapplied": False,
        "endpoint_cancellation_proved": False,
        "terminal_closure_proved": False,
        "ledger_212_changed": False,
        "ledger_212_consulted": False,
    }


def assemble(context: Any, core: Any = None, collar: Any = None, **kwargs: Any) -> dict[str, Any]:
    """Alias for :func:`build`."""

    return build(context, core, collar, **kwargs)


def combine(context: Any, core: Any = None, collar: Any = None, **kwargs: Any) -> dict[str, Any]:
    """Alias for :func:`build`."""

    return build(context, core, collar, **kwargs)


def construct(context: Any, core: Any = None, collar: Any = None, **kwargs: Any) -> dict[str, Any]:
    """Alias for :func:`build`."""

    return build(context, core, collar, **kwargs)


def enclose(context: Any, core: Any = None, collar: Any = None, **kwargs: Any) -> dict[str, Any]:
    """Alias for :func:`build`."""

    return build(context, core, collar, **kwargs)
