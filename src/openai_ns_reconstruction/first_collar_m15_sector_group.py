"""Common first-collar finite-prefix sector grouping.

This leaf combines already weighted axial, swirl, and pressure sectors.  It
does not evaluate a source, apply a physical factor, materialize a formal
atom, or infer a cancellation.  Each input record is retained as a member of
the output group together with its original provenance.

The public :func:`combine` API accepts one ordinary-Z ``0``/``1`` row for
each of ``axial``, ``swirl``, and ``pressure``.  A row may be a two-element
sequence ``[Z0, Z1]`` or a mapping with exactly the keys ``0`` and ``1``
(the JSON spellings ``"0"`` and ``"1"`` are also accepted).  Every sector
coefficient must already be a caller-owned directed interval.  Coefficients
are combined only when the exact ``(powers, uncertainty, error_classes)``
key agrees.

The module deliberately has no dependency on repository producers or stored
receipts.  In particular, it is suitable for a caller that has already
hydrated the saved sector lists without requiring an unsaved function.
The caller must authenticate common source semantics and hydrate all intervals
with one outward-rounded context.  Atom names alone do not prove identical
atom definitions, weights, domain, or owner.  Member records are retained by
reference; callers must not mutate them before serializing the result.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any


SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
ORDINARY_Z_ORDERS = (0, 1)
CHANNELS = ("axial", "swirl", "pressure")

__all__ = [
    "SOURCE_POWER_BASIS",
    "ORDINARY_Z_ORDERS",
    "CHANNELS",
    "combine",
    "build",
]


def _is_sequence(value: Any) -> bool:
    return isinstance(value, Sequence) and not isinstance(
        value, (str, bytes, bytearray)
    )


def _reject_float(value: Any, label: str) -> None:
    """Reject Python floating-point values in structural input."""

    if isinstance(value, float):
        raise TypeError(f"{label} must not contain floating-point values")


def _typed_key(value: Any, label: str) -> tuple[type[Any], Any]:
    _reject_float(value, label)
    try:
        hash(value)
    except TypeError as exc:
        raise TypeError(f"{label} must contain hashable values") from exc
    return type(value), value


def _unique_preserving(values: Sequence[Any], label: str) -> list[Any]:
    result: list[Any] = []
    seen: set[tuple[type[Any], Any]] = set()
    for index, value in enumerate(values):
        key = _typed_key(value, f"{label}[{index}]")
        if key not in seen:
            seen.add(key)
            result.append(value)
    return result


def _validate_source_atoms(source_atoms: Any) -> tuple[str, ...]:
    """Validate the canonical source-atom names and their order."""

    if isinstance(source_atoms, Mapping):
        names = tuple(source_atoms.keys())
        _reject_float(names, "source_atoms")
        if set(names) != set(SOURCE_POWER_BASIS):
            raise ValueError(
                "source_atoms must contain exactly hb, F0_anchor, "
                "Pstar_squared, and Lambda"
            )
        # Mappings carry named atoms; their key order is not semantic.  The
        # returned basis remains the one used by every powers four-tuple.
        return SOURCE_POWER_BASIS

    if not _is_sequence(source_atoms):
        raise TypeError("source_atoms must be the canonical atom sequence or mapping")
    names = tuple(source_atoms)
    if names != SOURCE_POWER_BASIS:
        raise ValueError(
            "source_atoms must be exactly "
            "('hb', 'F0_anchor', 'Pstar_squared', 'Lambda')"
        )
    return SOURCE_POWER_BASIS


def _validate_order_key(raw: Any, label: str) -> int:
    if isinstance(raw, bool) or isinstance(raw, float):
        raise TypeError(f"{label} must use ordinary-Z integer keys 0 and 1")
    if isinstance(raw, int):
        order = raw
    elif isinstance(raw, str) and raw in {"0", "1"}:
        order = int(raw)
    else:
        raise TypeError(f"{label} must use ordinary-Z integer keys 0 and 1")
    if order not in ORDINARY_Z_ORDERS:
        raise ValueError(f"{label} contains unsupported ordinary-Z order {order}")
    return order


def _normalize_channel(value: Any, label: str) -> dict[int, list[Any]]:
    """Normalize exactly two ordinary-Z rows while preserving row order."""

    if isinstance(value, Mapping):
        rows: dict[int, Any] = {}
        for raw_order, row in value.items():
            order = _validate_order_key(raw_order, f"{label} order")
            if order in rows:
                raise ValueError(f"{label} contains duplicate ordinary-Z order {order}")
            rows[order] = row
        if set(rows) != set(ORDINARY_Z_ORDERS):
            raise ValueError(f"{label} must provide exactly ordinary-Z orders 0 and 1")
    elif _is_sequence(value):
        if len(value) != len(ORDINARY_Z_ORDERS):
            raise ValueError(f"{label} must contain exactly two rows in order Z0, Z1")
        rows = {order: value[order] for order in ORDINARY_Z_ORDERS}
    else:
        raise TypeError(f"{label} must be a Z0/Z1 sequence or mapping")

    result: dict[int, list[Any]] = {}
    for order in ORDINARY_Z_ORDERS:
        row = rows[order]
        if not _is_sequence(row):
            raise TypeError(f"{label}[Z{order}] must be a sector list")
        if not row:
            raise ValueError(f"{label}[Z{order}] must not be empty")
        result[order] = list(row)
    return result


def _is_directed_interval(value: Any) -> bool:
    """Recognize an existing interval without importing an arithmetic package."""

    return (
        not isinstance(value, (bool, int, float, complex, str, bytes, bytearray))
        and hasattr(value, "_mpi_")
        and hasattr(value, "a")
        and hasattr(value, "b")
    )


def _validate_powers(value: Any, label: str) -> tuple[int, int, int, int]:
    if not _is_sequence(value) or len(value) != 4:
        raise ValueError(f"{label} must be a four-integer source-power tuple")
    powers: list[int] = []
    for index, power in enumerate(value):
        if isinstance(power, bool) or not isinstance(power, int):
            raise TypeError(f"{label}[{index}] must be an integer, not a float or bool")
        powers.append(int(power))
    return tuple(powers)  # type: ignore[return-value]


def _validate_origin(origin: Any, label: str) -> tuple[list[Any], bool, tuple[Any, ...]]:
    if not isinstance(origin, Mapping):
        raise TypeError(f"{label} must be a provenance mapping")
    required = {"source_ids", "uncertainty", "error_classes"}
    missing = required.difference(origin)
    if missing:
        raise ValueError(f"{label} is missing {', '.join(sorted(missing))}")

    source_ids = origin["source_ids"]
    if not _is_sequence(source_ids) or not source_ids:
        raise ValueError(f"{label}.source_ids must be a non-empty sequence")
    source_ids_copy = list(source_ids)
    # Validate each ID and retain the original order in the member record.
    _unique_preserving(source_ids_copy, f"{label}.source_ids")

    uncertainty = origin["uncertainty"]
    if not isinstance(uncertainty, bool):
        raise TypeError(f"{label}.uncertainty must be a bool")

    error_classes = origin["error_classes"]
    if not _is_sequence(error_classes):
        raise TypeError(f"{label}.error_classes must be a sequence")
    classes = tuple(error_classes)
    for index, item in enumerate(classes):
        _typed_key(item, f"{label}.error_classes[{index}]")
    return source_ids_copy, uncertainty, classes


def _validate_term(
    term: Any,
    label: str,
) -> tuple[
    tuple[int, int, int, int],
    Any,
    dict[str, Any],
    tuple[tuple[int, int, int, int], bool, tuple[tuple[type[Any], Any], ...]],
]:
    if not isinstance(term, Mapping):
        raise TypeError(f"{label} must be a sector mapping")
    required = {"powers", "coefficient", "origin"}
    missing = required.difference(term)
    if missing:
        raise ValueError(f"{label} is missing {', '.join(sorted(missing))}")

    powers = _validate_powers(term["powers"], f"{label}.powers")
    coefficient = term["coefficient"]
    if not _is_directed_interval(coefficient):
        raise TypeError(f"{label}.coefficient must be an existing directed interval")
    origin = term["origin"]
    source_ids, uncertainty, error_classes = _validate_origin(origin, f"{label}.origin")

    # If a caller carries a basis marker on a sector, validate it rather than
    # silently combining sectors written in a different atom convention.
    for key in ("source_power_basis", "source_atoms"):
        if key in term and term[key] is not None:
            _validate_source_atoms(term[key])
    if "source_power_basis" in origin and origin["source_power_basis"] is not None:
        _validate_source_atoms(origin["source_power_basis"])

    # Include the value type in the structural key so labels such as ``1``
    # and ``True`` cannot collapse through Python's ordinary equality rules.
    classes_key = tuple(
        _typed_key(item, f"{label}.origin.error_classes[{index}]")
        for index, item in enumerate(error_classes)
    )
    key = (powers, uncertainty, classes_key)
    return powers, coefficient, dict(origin), key


def _append_member(
    group: dict[str, Any],
    *,
    channel: str,
    order: int,
    index: int,
    record: Any,
    origin: Mapping[str, Any],
    source_ids: Sequence[Any],
) -> None:
    group["members"].append(
        {
            "channel": channel,
            "ordinary_Z_order": order,
            "index": index,
            "record": record,
            "provenance": dict(origin),
        }
    )
    group["origin"]["source_ids"] = _unique_preserving(
        [*group["origin"]["source_ids"], *source_ids],
        "group.origin.source_ids",
    )


def combine(
    axial: Any,
    swirl: Any,
    pressure: Any,
    *,
    source_atoms: Any,
) -> dict[str, Any]:
    """Group three already-weighted sector channels through ordinary Z1.

    Parameters
    ----------
    axial, swirl, pressure:
        Each is either ``[Z0_terms, Z1_terms]`` or a mapping with exactly
        ordinary-Z keys ``0`` and ``1``.  Terms must contain ``powers``, an
        existing directed-interval ``coefficient``, and an ``origin`` with
        ``source_ids``, boolean ``uncertainty``, and sequence-valued
        ``error_classes``.
    source_atoms:
        The canonical four atoms, either in exact sequence order or as a
        mapping with exactly those names.

    Returns
    -------
    dict
        ``groups[0]`` and ``groups[1]`` contain the first-seen grouped
        sectors.  Every group has ``members`` entries containing each
        original record and copied provenance.  Coefficients are combined by
        interval addition only; no physical multiplication or cancellation is
        performed.
    """

    basis = _validate_source_atoms(source_atoms)
    channels = {
        "axial": _normalize_channel(axial, "axial"),
        "swirl": _normalize_channel(swirl, "swirl"),
        "pressure": _normalize_channel(pressure, "pressure"),
    }

    grouped: dict[int, list[dict[str, Any]]] = {0: [], 1: []}
    by_key: dict[int, dict[tuple[Any, bool, tuple[Any, ...]], dict[str, Any]]] = {
        0: {},
        1: {},
    }

    for order in ORDINARY_Z_ORDERS:
        for channel in CHANNELS:
            for index, record in enumerate(channels[channel][order]):
                powers, coefficient, origin, key = _validate_term(
                    record,
                    f"{channel}[Z{order}][{index}]",
                )
                group = by_key[order].get(key)
                if group is None:
                    group = {
                        "powers": powers,
                        "coefficient": coefficient,
                        "origin": {
                            "operation": "outward_interval_sector_group",
                            "source_ids": [],
                            "uncertainty": key[1],
                            "error_classes": list(origin["error_classes"]),
                            "member_count": 0,
                        },
                        "members": [],
                    }
                    by_key[order][key] = group
                    grouped[order].append(group)
                else:
                    # Existing interval arithmetic supplies outward rounding;
                    # no scalar conversion or new physical factor is used.
                    group["coefficient"] = group["coefficient"] + coefficient

                _append_member(
                    group,
                    channel=channel,
                    order=order,
                    index=index,
                    record=record,
                    origin=origin,
                    source_ids=origin["source_ids"],
                )
                group["origin"]["member_count"] += 1

    return {
        "source_power_basis": list(basis),
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "channels": list(CHANNELS),
        "groups": grouped,
        "group_key": "(powers, uncertainty, error_classes)",
        "member_records_retained": True,
        "provenance_retained": True,
        "source_ids_retained": True,
        "outward_interval_addition_only": True,
        "physical_factors_applied": False,
        "cancellation_performed": False,
    }


def build(
    axial: Any,
    swirl: Any,
    pressure: Any,
    *,
    source_atoms: Any,
) -> dict[str, Any]:
    """Alias for :func:`combine` with the same argument order."""

    return combine(axial, swirl, pressure, source_atoms=source_atoms)

