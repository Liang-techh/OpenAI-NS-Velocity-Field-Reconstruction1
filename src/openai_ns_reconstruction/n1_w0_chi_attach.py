"""Attach a saved finite ``W0_chi`` correction to a segmented candidate.

This module is intentionally a leaf.  ``parent`` is an already constructed
segmented candidate ``formal_sum`` and the two other packets contain already
physically integrated records.  The leaf does not call a producer, evaluate a
source, integrate a radial weight, or apply geometry, ``P1`` or a tail factor.

The only arithmetic performed here is outward interval addition between
signed correction records with exactly equal formal powers.  Parent records
remain separate objects, and error-only records remain in ``error_stream``;
neither stream is silently converted into the other.  Source paths,
references, and complete input origins are retained by reference in the
returned expression and in every aggregate's ``input_terms``/``input_origins``
fields.  Enclosures are opaque directed intervals, never point values.

Contract
========

``build(context, parent, signed_correction, error_only)`` requires a caller
owned interval context exposing ``mpf`` and ``make_mpf``.  Each packet is a
mapping with exactly ordinary-Z keys ``0`` and ``1`` (``"Z0"``/``"Z1"`` are
accepted as spelling aliases) and non-empty ``list`` rows.  Every record is a
mapping with a four-integer ``powers`` sequence, an existing directed
``coefficient`` interval, and a non-empty provenance ``origin`` mapping with
``source_ids`` (or an equivalent singular source-id field).  No record is
copied or numerically coerced.

Signed records must carry ``signed_field is True`` either on the record or in
its origin (when present in both places both must be true).  Error records must
carry ``origin.error_only is True`` and may not carry a true signed marker.
The explicit correction branch is ``W0_chi``; ``E01``, ``E10`` and ``E11``
labels are rejected.  Same-family or matching source labels are retained as
metadata only and never establish a live owner.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any


SCHEMA = "n1_w0_chi_attach_v1"
SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
ORDINARY_Z_ORDERS = (0, 1)
BRANCH = "W0_chi"
FORBIDDEN_BRANCHES = frozenset(("E01", "E10", "E11"))
GROUPING_KEY = "powers"

__all__ = [
    "SCHEMA",
    "SOURCE_POWER_BASIS",
    "ORDINARY_Z_ORDERS",
    "BRANCH",
    "GROUPING_KEY",
    "build",
    "attach",
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
    """Reject host floats in metadata without traversing interval objects."""

    if isinstance(value, float):
        raise TypeError(f"{label} must not contain floating-point values")
    if isinstance(value, bool) and not allow_bool:
        raise TypeError(f"{label} must not contain boolean values")
    if isinstance(value, Mapping):
        for key, item in value.items():
            _reject_float(key, f"{label}.key", allow_bool=allow_bool)
            _reject_float(item, f"{label}[{key!r}]", allow_bool=allow_bool)
    elif _is_sequence(value):
        for index, item in enumerate(value):
            _reject_float(item, f"{label}[{index}]", allow_bool=allow_bool)


def _is_interval(value: Any) -> bool:
    """Recognize a caller-created mpmath interval without importing mpmath."""

    return (
        not isinstance(value, (bool, int, float, complex, str, bytes, bytearray))
        and hasattr(value, "_mpi_")
        and hasattr(value, "a")
        and hasattr(value, "b")
    )


def _check_context(context: Any, coefficient: Any, label: str) -> None:
    """Check context compatibility without converting or midpointing a value."""

    if not _is_interval(coefficient):
        raise TypeError(f"{label} must be an existing directed interval")
    try:
        # There is no stable public interval-owner attribute.  This operation
        # exercises the caller context against an endpoint and does not read a
        # midpoint, bound, or numerical value from the enclosure.
        zero = context.mpf(0)
        _ = coefficient.a + zero
    except (AttributeError, TypeError, ValueError, ArithmeticError) as exc:
        raise TypeError(
            f"{label} is not compatible with the caller's interval context"
        ) from exc


def _typed_key(value: Any, label: str) -> tuple[type[Any], Any]:
    if isinstance(value, (bool, float)):
        raise TypeError(f"{label} must not be bool or float")
    try:
        hash(value)
    except TypeError as exc:
        raise TypeError(f"{label} must be hashable") from exc
    return type(value), value


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
    if isinstance(raw, bool):
        raise TypeError(f"{label} must be ordinary-Z order 0 or 1")
    if isinstance(raw, int):
        value = raw
    elif isinstance(raw, str):
        text = raw.strip().lower()
        if text in {"z0", "0"}:
            value = 0
        elif text in {"z1", "1"}:
            value = 1
        else:
            raise ValueError(f"{label} must be ordinary-Z order 0 or 1")
    else:
        raise TypeError(f"{label} must be ordinary-Z order 0 or 1")
    if value not in ORDINARY_Z_ORDERS:
        raise ValueError(f"{label} must be ordinary-Z order 0 or 1")
    return value


def _rows(value: Any, label: str, *, allow_formal_wrapper: bool = False) -> tuple[dict[int, list[Any]], Any]:
    """Validate an order packet and return rows plus the original packet."""

    original = value
    has_order_key = False
    if isinstance(value, Mapping):
        has_order_key = any(
            isinstance(key, int)
            and not isinstance(key, bool)
            and key in ORDINARY_Z_ORDERS
            or isinstance(key, str)
            and key.strip().lower() in {"0", "1", "z0", "z1"}
            for key in value
        )
    if (
        allow_formal_wrapper
        and isinstance(value, Mapping)
        and "formal_sum" in value
        and not has_order_key
    ):
        value = value["formal_sum"]
    if not isinstance(value, Mapping):
        raise TypeError(f"{label} must be a mapping with Z0 and Z1 list rows")

    normalized: dict[int, Any] = {}
    for raw_order, row in value.items():
        order = _order(raw_order, f"{label} order")
        if order in normalized:
            raise ValueError(f"{label} contains duplicate ordinary-Z order {order}")
        normalized[order] = row
    if set(normalized) != set(ORDINARY_Z_ORDERS):
        raise ValueError(f"{label} must provide exactly ordinary-Z orders 0 and 1")

    result: dict[int, list[Any]] = {}
    for order in ORDINARY_Z_ORDERS:
        row = normalized[order]
        if not isinstance(row, list):
            raise TypeError(f"{label}[Z{order}] must be a list")
        if not row:
            raise ValueError(f"{label}[Z{order}] must not be empty")
        result[order] = row
    return result, original


def _branch_rows(value: Any, label: str) -> tuple[dict[int, list[Any]], dict[int, list[Any]], Any]:
    """Extract saved ``swirl`` and ``swirl_error`` rows without re-expanding cells."""

    if not isinstance(value, Mapping):
        raise TypeError(f"{label} must be a W0_chi branch mapping")
    signed: dict[int, list[Any]] = {}
    errors: dict[int, list[Any]] = {}
    seen_orders: set[int] = set()
    for raw_order, payload in value.items():
        # Branch metadata (if present) is retained by the caller but is not an
        # ordinary-Z row.  A direct order packet is handled by _rows instead.
        try:
            order = _order(raw_order, f"{label} order")
        except (TypeError, ValueError):
            continue
        if order in seen_orders:
            raise ValueError(f"{label} contains duplicate ordinary-Z order {order}")
        seen_orders.add(order)
        if not isinstance(payload, Mapping):
            raise TypeError(f"{label}[Z{order}] must be a channel mapping")
        for channel in ("axial", "axial_error"):
            if channel in payload and (not isinstance(payload[channel], list) or payload[channel]):
                raise ValueError(f"{label}[Z{order}].{channel} must be an empty list")
        signed_row = payload.get("swirl")
        error_row = payload.get("swirl_error")
        if not isinstance(signed_row, list) or not signed_row:
            raise ValueError(f"{label}[Z{order}].swirl must be a non-empty list")
        if not isinstance(error_row, list) or not error_row:
            raise ValueError(
                f"{label}[Z{order}].swirl_error must be a non-empty list"
            )
        signed[order] = signed_row
        errors[order] = error_row
    if set(seen_orders) != set(ORDINARY_Z_ORDERS):
        raise ValueError(f"{label} must provide exactly ordinary-Z orders 0 and 1")
    return signed, errors, value


def _extract_correction(
    value: Any,
    label: str,
    *,
    branch: Any = None,
) -> tuple[dict[int, list[Any]], dict[int, list[Any]], Any, Any]:
    """Accept either direct rows or a saved branch wrapper.

    ``correction`` is metadata/receipt-shaped in the saved result, while
    ``w0_chi`` is the already selected branch.  The branch itself is the only
    payload consumed; no cells or producer records are traversed.
    """

    selected = branch if branch is not None else value
    if branch is not None and value is not None and value is not branch:
        # A correction wrapper and an explicitly selected branch are both
        # allowed, but the wrapper must contain the same object by reference.
        if isinstance(value, Mapping):
            branches = value.get("branches")
            if isinstance(branches, Mapping) and branches.get(BRANCH) is branch:
                pass
            elif value is not branch:
                raise ValueError("correction and w0_chi disagree")
    if isinstance(selected, Mapping) and "branches" in selected:
        branches = selected["branches"]
        if not isinstance(branches, Mapping) or BRANCH not in branches:
            raise ValueError(f"{label}.branches must retain {BRANCH}")
        selected = branches[BRANCH]
    if isinstance(selected, Mapping):
        # A direct order packet is a signed-only spelling; errors must still
        # be supplied separately by the caller in that form.
        order_keys = []
        for raw_order in selected:
            try:
                order_keys.append(_order(raw_order, f"{label} order"))
            except (TypeError, ValueError):
                pass
        if set(order_keys) == set(ORDINARY_Z_ORDERS) and all(
            isinstance(selected.get(raw_order), list)
            for raw_order in selected
            if raw_order in (0, 1, "0", "1", "Z0", "z0", "Z1", "z1")
        ):
            rows, original = _rows(selected, label)
            return rows, {0: [], 1: []}, original, None
        signed, errors, original = _branch_rows(selected, label)
        return signed, errors, original, original
    rows, original = _rows(selected, label)
    return rows, {0: [], 1: []}, original, None


def _source_items(value: Any) -> list[Any]:
    if value is None:
        return []
    if isinstance(value, (str, bytes, bytearray)):
        return [value]
    if _is_sequence(value):
        return list(value)
    return [value]


def _origin_sources(origin: Mapping[str, Any]) -> list[Any]:
    values: list[Any] = []
    for key in ("source_ids", "source_id", "sources", "source", "input_id"):
        if key in origin:
            values.extend(_source_items(origin[key]))
    return values


def _term_sources(term: Mapping[str, Any], origin: Mapping[str, Any]) -> list[Any]:
    values = _origin_sources(origin)
    for key in ("source_ids", "source_id", "sources", "source", "input_id"):
        if key in term:
            values.extend(_source_items(term[key]))
    return values


def _marker(value: Any, label: str) -> bool | None:
    if value is None:
        return None
    if not isinstance(value, bool):
        raise TypeError(f"{label} must be boolean")
    return value


def _marker_value(term: Mapping[str, Any], origin: Mapping[str, Any], name: str) -> bool | None:
    top = _marker(term.get(name), f"term.{name}") if name in term else None
    nested = _marker(origin.get(name), f"origin.{name}") if name in origin else None
    if top is not None and nested is not None and top != nested:
        raise ValueError(f"term.{name} and origin.{name} disagree")
    return top if top is not None else nested


def _text_labels(term: Mapping[str, Any], origin: Mapping[str, Any]) -> set[str]:
    labels: set[str] = set()
    for mapping in (term, origin):
        for key in (
            "branch",
            "component",
            "route",
            "word",
            "field",
            "channel",
            "source_family",
            "family",
        ):
            value = mapping.get(key)
            if isinstance(value, str):
                labels.add(value.strip().lower())
    return labels


def _reject_forbidden_branch(term: Mapping[str, Any], origin: Mapping[str, Any], label: str) -> None:
    labels = _text_labels(term, origin)
    for forbidden in FORBIDDEN_BRANCHES:
        if forbidden.lower() in labels:
            raise ValueError(f"{label} contains forbidden correction branch {forbidden}")
    if "p1" in labels:
        raise ValueError(f"{label} P1 is not part of this m15 attachment")


def _validate_origin(origin: Any, label: str) -> Mapping[str, Any]:
    if not isinstance(origin, Mapping):
        raise TypeError(f"{label} must be a provenance mapping")
    sources = _origin_sources(origin)
    if not sources:
        raise ValueError(f"{label} must retain at least one source id/reference")
    for index, source in enumerate(sources):
        _typed_key(source, f"{label}.source[{index}]")
    for key in ("error_only", "signed_field"):
        if key in origin:
            _marker(origin[key], f"{label}.{key}")
    _reject_float(origin, label, allow_bool=True)
    return origin


def _validate_term(
    context: Any,
    raw: Any,
    label: str,
    *,
    stream: str,
) -> tuple[tuple[int, int, int, int], Mapping[str, Any], bool | None, bool | None]:
    if not isinstance(raw, Mapping):
        raise TypeError(f"{label} must be a term mapping")
    if "powers" not in raw:
        raise ValueError(f"{label} requires powers")
    if "coefficient" not in raw:
        raise ValueError(f"{label} requires coefficient")
    powers = _powers(raw["powers"], f"{label}.powers")
    _check_context(context, raw["coefficient"], f"{label}.coefficient")
    origin = _validate_origin(raw.get("origin"), f"{label}.origin")
    signed = _marker_value(raw, origin, "signed_field")
    error_only = _marker_value(raw, origin, "error_only")
    _reject_forbidden_branch(raw, origin, label)
    _reject_float(
        {key: value for key, value in raw.items() if key != "coefficient"},
        label,
        allow_bool=True,
    )

    if stream == "signed":
        if signed is not True:
            raise ValueError(f"{label} must have signed_field=True")
        if error_only is True:
            raise ValueError(f"{label} may not be error_only")
    elif stream == "error":
        # This marker is intentionally required in the nested origin.  A
        # top-level alias alone is insufficient: it would hide a malformed
        # saved tail origin and blur the signed/error distinction.
        if origin.get("error_only") is not True:
            raise ValueError(f"{label}.origin.error_only must be True")
        if signed is True:
            raise ValueError(f"{label} error-only term may not be signed_field=True")
    elif stream == "parent":
        # Parent terms are already authenticated by their producer.  Validate
        # their shape and preserve their markers without assigning a new one.
        pass
    else:  # pragma: no cover - internal call contract
        raise ValueError(f"unknown stream {stream!r}")
    return powers, origin, signed, error_only


def _source_paths(term: Mapping[str, Any], origin: Mapping[str, Any]) -> tuple[Any, ...]:
    values: list[Any] = []
    for mapping in (term, origin):
        # A singular ``path`` is itself commonly a list of path components;
        # retain that list object rather than flattening it into components.
        for key in ("source_path", "path", "input_path"):
            if key in mapping:
                values.append(mapping[key])
        for key in ("source_paths", "paths", "input_paths"):
            if key in mapping:
                raw = mapping[key]
                if _is_sequence(raw):
                    values.extend(raw)
                else:
                    values.append(raw)
    return tuple(values)


def _source_refs(term: Mapping[str, Any], origin: Mapping[str, Any]) -> tuple[Any, ...]:
    values: list[Any] = []
    for mapping in (term, origin):
        # Tuples are valid opaque references (for example ``(order, index)``)
        # and must not be mistaken for a sequence of independent references.
        for key in ("source_ref", "source_reference", "reference", "ref"):
            if key in mapping:
                values.append(mapping[key])
        for key in ("source_refs", "source_references", "references", "refs"):
            if key in mapping:
                raw = mapping[key]
                if _is_sequence(raw):
                    values.extend(raw)
                else:
                    values.append(raw)
    return tuple(values)


def _aggregate_signed(
    context: Any,
    rows: Mapping[int, Sequence[Any]],
) -> dict[int, list[dict[str, Any]]]:
    """Add signed correction intervals only within equal formal powers."""

    del context  # interval operators remain owned by the input coefficients
    output: dict[int, list[dict[str, Any]]] = {0: [], 1: []}
    for order in ORDINARY_Z_ORDERS:
        groups: dict[tuple[int, int, int, int], dict[str, Any]] = {}
        order_terms = rows[order]
        for index, raw in enumerate(order_terms):
            powers = tuple(raw["powers"])
            group = groups.get(powers)
            if group is None:
                group = {
                    "powers": powers,
                    "coefficient": raw["coefficient"],
                    "input_terms": [raw],
                    "input_origins": [raw["origin"]],
                    "source_paths": [_source_paths(raw, raw["origin"])],
                    "source_refs": [_source_refs(raw, raw["origin"])],
                    "source_ids": [_term_sources(raw, raw["origin"])],
                    "indices": [index],
                }
                groups[powers] = group
            else:
                # This is the sole numerical operation of the leaf.  It is an
                # interval sum, not a midpoint or a point-value computation.
                group["coefficient"] = group["coefficient"] + raw["coefficient"]
                group["input_terms"].append(raw)
                group["input_origins"].append(raw["origin"])
                group["source_paths"].append(_source_paths(raw, raw["origin"]))
                group["source_refs"].append(_source_refs(raw, raw["origin"]))
                group["source_ids"].append(_term_sources(raw, raw["origin"]))
                group["indices"].append(index)

        for group in groups.values():
            input_terms = tuple(group["input_terms"])
            input_origins = tuple(group["input_origins"])
            source_ids = tuple(group["source_ids"])
            source_paths = tuple(group["source_paths"])
            source_refs = tuple(group["source_refs"])
            origin = {
                "operation": "W0_chi_signed_correction_net",
                "branch": BRANCH,
                "ordinary_Z_order": order,
                "powers": tuple(group["powers"]),
                "signed_field": True,
                "error_only": False,
                "uncertainty": any(
                    bool(item.get("uncertainty", False))
                    for item in input_origins
                    if isinstance(item, Mapping)
                ),
                "source_ids": [item for row in group["source_ids"] for item in row],
                "source_ids_by_term": source_ids,
                "source_paths_by_term": source_paths,
                "source_refs_by_term": source_refs,
                "input_origins": input_origins,
                "input_term_count": len(input_terms),
                "input_indices": tuple(group["indices"]),
                "provenance_retained": True,
                "source_equality_proved": False,
                "symbolic_cancellation_claimed": False,
                "physical_weight_reapplied": False,
                "geometry_reapplied": False,
                "P1_added": False,
                "tail_added": False,
            }
            output[order].append(
                {
                    "kind": "signed_correction",
                    "source": BRANCH,
                    "source_label": BRANCH,
                    "branch": BRANCH,
                    "ordinary_Z_order": order,
                    "powers": tuple(group["powers"]),
                    "coefficient": group["coefficient"],
                    "signed_field": True,
                    "error_only": False,
                    "input_terms": input_terms,
                    "original_terms": input_terms,
                    "input_origins": input_origins,
                    "original_origins": input_origins,
                    "source_ids_by_term": source_ids,
                    "source_paths_by_term": source_paths,
                    "source_refs_by_term": source_refs,
                    "source_ids": tuple(
                        item for row in group["source_ids"] for item in row
                    ),
                    "source_paths": source_paths,
                    "source_refs": source_refs,
                    "source_path": source_paths,
                    "source_ref": source_refs,
                    "paths": source_paths,
                    "refs": source_refs,
                    "origin": origin,
                    "formal_power_source": "input terms unchanged",
                    "coefficient_is_interval_sum": len(input_terms) > 1,
                    "parent_term_retained_separately": True,
                    "error_stream_separate": True,
                }
            )
    return output


def _formal_envelopes(
    order: int,
    parent_terms: Sequence[Any],
    correction_terms: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    for index, raw in enumerate(parent_terms):
        result.append(
            {
                "kind": "prior",
                "source": "segmented_parent",
                "source_label": "segmented_parent",
                "ordinary_Z_order": order,
                "index": index,
                "term": raw,
                "powers": tuple(raw["powers"]),
                "coefficient": raw["coefficient"],
                "origin": raw["origin"],
                "source_paths": _source_paths(raw, raw["origin"]),
                "source_refs": _source_refs(raw, raw["origin"]),
                "signed_field": _marker_value(raw, raw["origin"], "signed_field"),
                "error_only": _marker_value(raw, raw["origin"], "error_only"),
                "by_reference": True,
            }
        )
    for index, raw in enumerate(correction_terms):
        result.append(
            {
                "kind": "signed_correction",
                "source": BRANCH,
                "source_label": BRANCH,
                "ordinary_Z_order": order,
                "index": index,
                "term": raw,
                "powers": tuple(raw["powers"]),
                "coefficient": raw["coefficient"],
                "origin": raw["origin"],
                "source_paths": raw["source_paths"],
                "source_refs": raw["source_refs"],
                "signed_field": True,
                "error_only": False,
                "by_reference": True,
            }
        )
    return result


def _net_by_powers(net: Mapping[int, Sequence[Mapping[str, Any]]]) -> dict[int, dict[tuple[int, int, int, int], Any]]:
    return {
        order: {tuple(term["powers"]): term["coefficient"] for term in net[order]}
        for order in ORDINARY_Z_ORDERS
    }


def build(
    context: Any,
    parent: Any = None,
    signed_correction: Any = None,
    error_only: Any = None,
    *,
    parent_formal_sum: Any = None,
    prior_formal_sum: Any = None,
    correction: Any = None,
    w0_chi: Any = None,
    error_tail: Any = None,
    source_family: Any = None,
    source_logs: Any = None,
    source_atoms: Any = None,
    source_index: Any = None,
    prior_source_index: Any = None,
    correction_source_index: Any = None,
    owner_certificate: Any = None,
    scale_certificate: Any = None,
    certificate: Any = None,
    source_ref: Any = None,
) -> dict[str, Any]:
    """Attach saved W0-chi terms to a saved segmented formal sum.

    ``parent_formal_sum``, ``correction`` and ``error_tail`` are keyword
    aliases for the three positional packets.  Supplying an alias alongside
    its positional counterpart is rejected unless both objects are identical;
    this prevents an accidental second parent or correction stream.
    """

    if not hasattr(context, "mpf") or not hasattr(context, "make_mpf"):
        raise TypeError("context must be the caller's MPIntervalContext")
    if parent_formal_sum is not None:
        if parent is not None and parent is not parent_formal_sum:
            raise ValueError("parent and parent_formal_sum disagree")
        parent = parent_formal_sum
    if prior_formal_sum is not None:
        if parent is not None and parent is not prior_formal_sum:
            raise ValueError("parent and prior_formal_sum disagree")
        parent = prior_formal_sum
    if correction is not None:
        if signed_correction is not None and signed_correction is not correction:
            raise ValueError("signed_correction and correction disagree")
        if signed_correction is None:
            signed_correction = correction
    if w0_chi is not None and signed_correction is None:
        signed_correction = w0_chi
    if error_tail is not None:
        if error_only is not None and error_only is not error_tail:
            raise ValueError("error_only and error_tail disagree")
        error_only = error_tail
    if parent is None or signed_correction is None:
        raise TypeError("parent and signed_correction are required")

    parent_rows, parent_input = _rows(
        parent, "parent", allow_formal_wrapper=True
    )
    correction_rows, derived_errors, correction_input, selected_branch = _extract_correction(
        signed_correction,
        "signed_correction",
        branch=w0_chi,
    )
    if error_only is None:
        if selected_branch is None:
            raise TypeError("error_only is required when correction has no branch tail")
        error_rows = derived_errors
        error_input = selected_branch
        if any(not error_rows[order] for order in ORDINARY_Z_ORDERS):
            raise ValueError("selected W0_chi branch must retain non-empty swirl_error rows")
    else:
        error_rows, error_input = _rows(error_only, "error_only")
        if selected_branch is not None:
            if error_only is not selected_branch:
                # The selected saved branch carries its error stream; an
                # explicit packet may replace it only when it is identical.
                raise ValueError("explicit error_only cannot replace selected branch tail")

    for order in ORDINARY_Z_ORDERS:
        for index, raw in enumerate(parent_rows[order]):
            _validate_term(
                context, raw, f"parent[Z{order}][{index}]", stream="parent"
            )
        for index, raw in enumerate(correction_rows[order]):
            _validate_term(
                context,
                raw,
                f"signed_correction[Z{order}][{index}]",
                stream="signed",
            )
        for index, raw in enumerate(error_rows[order]):
            _validate_term(
                context,
                raw,
                f"error_only[Z{order}][{index}]",
                stream="error",
            )

    for value, label in (
        (source_family, "source_family"),
        (source_logs, "source_logs"),
        (source_atoms, "source_atoms"),
        (source_index, "source_index"),
        (prior_source_index, "prior_source_index"),
        (correction_source_index, "correction_source_index"),
        (owner_certificate, "owner_certificate"),
        (scale_certificate, "scale_certificate"),
        (certificate, "certificate"),
        (source_ref, "source_ref"),
    ):
        if value is not None:
            _reject_float(value, label, allow_bool=True)

    signed_net = _aggregate_signed(context, correction_rows)
    formal_sum: dict[int, list[Any]] = {
        order: list(parent_rows[order]) + list(signed_net[order])
        for order in ORDINARY_Z_ORDERS
    }
    expression = {
        order: _formal_envelopes(order, parent_rows[order], signed_net[order])
        for order in ORDINARY_Z_ORDERS
    }
    net_by_powers = _net_by_powers(signed_net)

    parent_error_stream: dict[int, list[Any]] = {0: [], 1: []}
    for order in ORDINARY_Z_ORDERS:
        for raw in parent_rows[order]:
            origin = raw["origin"]
            if origin.get("error_only") is True or raw.get("error_only") is True:
                parent_error_stream[order].append(raw)

    return {
        "schema": SCHEMA,
        "source_power_basis": list(SOURCE_POWER_BASIS),
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "branch": BRANCH,
        "source_family": source_family,
        "source_logs": source_logs,
        "source_atoms": source_atoms,
        "source_ref": source_ref,
        "parent_input": parent_input,
        "parent_formal_sum_input": parent_input,
        "prior_formal_sum": parent_rows,
        "signed_correction_input": signed_correction,
        "error_only_input": error_input,
        "correction_input": signed_correction,
        "selected_correction_input": correction_input,
        "w0_chi_input": selected_branch if selected_branch is not None else signed_correction,
        "source_index": source_index,
        "prior_source_index": prior_source_index,
        "correction_source_index": correction_source_index,
        "owner_certificate": owner_certificate,
        "scale_certificate": scale_certificate,
        "certificate": certificate,
        "prior": parent_rows,
        "parent": parent_rows,
        "prior_terms": parent_rows,
        "signed_correction": correction_rows,
        "correction": correction_rows,
        "signed_rows": signed_net,
        "signed": {
            0: list(signed_net[0]),
            1: list(signed_net[1]),
            "signed": {0: list(signed_net[0]), 1: list(signed_net[1])},
        },
        "W0_chi": {
            "signed": {0: list(signed_net[0]), 1: list(signed_net[1])},
            "error_only": {0: list(error_rows[0]), 1: list(error_rows[1])},
        },
        "signed_correction_net": signed_net,
        "signed_correction_sum": signed_net,
        "signed_net": signed_net,
        "correction_net": signed_net,
        "signed_correction_net_by_order": signed_net,
        "signed_correction_net_by_powers": net_by_powers,
        "error_only_terms": error_rows,
        "error_only_packet": {"error_only": {0: list(error_rows[0]), 1: list(error_rows[1])}},
        "error_stream": error_rows,
        "error_only": {
            0: list(error_rows[0]),
            1: list(error_rows[1]),
            "error_only": {0: list(error_rows[0]), 1: list(error_rows[1])},
        },
        "error_only_rows": error_rows,
        "error_tail": error_rows,
        "parent_error_stream": parent_error_stream,
        "formal_sum": formal_sum,
        "segmented_candidate": formal_sum,
        "candidate": formal_sum,
        "prior_plus_correction": formal_sum,
        "incremental": formal_sum,
        "formal_expression": expression,
        "incremental_formal_expression": expression,
        "formal_expression_by_order": {
            order: {
                "prior": list(parent_rows[order]),
                "correction": list(signed_net[order]),
                "signed": list(formal_sum[order]),
                "error": list(error_rows[order]),
                "terms": expression[order],
            }
            for order in ORDINARY_Z_ORDERS
        },
        "grouping_key": GROUPING_KEY,
        "signed_correction_added": True,
        "signed_correction_same_powers_only": True,
        "signed_correction_interval_sum_only": True,
        "parent_terms_retained_separately": True,
        "new_terms_retained_separately": True,
        "input_terms_retained_by_reference": True,
        "source_paths_retained": True,
        "source_refs_retained": True,
        "complete_origins_retained": True,
        "error_stream_retained_separately": True,
        "error_only_added_to_signed": False,
        "error_only_as_signed_field": False,
        "error_stream_summed": False,
        "error_stream_numeric_sum": False,
        "error_tail_is_independent_reference": False,
        "error_tail_kept_as_separate_reference": True,
        "error_tail_independence_not_proved": True,
        "physical_weight_reapplied": False,
        "radial_weight_reapplied": False,
        "geometry_reapplied": False,
        "P1_added": False,
        "tail_added": False,
        "pressure_endpoint_added": False,
        "independent_QP_added": False,
        "E01_added": False,
        "E10_added": False,
        "E11_added": False,
        "ledger_appended": False,
        "writes_ledger": False,
        "no_ledger": True,
        "ledger_212_changed": False,
        "same_family_is_live_owner": False,
        "owner_identity_proved": False,
        "source_equality_proved": False,
        "endpoint_cancellation_claimed": False,
        "terminal_closure_proved": False,
        "point_values_selected": False,
        "numeric_values_computed": False,
        "midpoint_values_used": False,
        "enclosures_are_numeric_values": False,
        "same_interval_context_validated": False,
        "same_interval_context_is_caller_obligation": True,
        "parent_recomputed": False,
        "parent_formal_sum_recomputed": False,
        "cells_resummed": False,
        "reaggregated_cells": False,
        "W0_chi_already_physically_integrated": True,
        "scope": (
            "attach already physically integrated finite W0_chi signed Z0/Z1 "
            "terms to a saved segmented parent; independent error-only tail "
            "is retained outside the signed formal sum"
        ),
    }


def attach(context: Any, parent: Any = None, signed_correction: Any = None, error_only: Any = None, **kwargs: Any) -> dict[str, Any]:
    """Alias for :func:`build`."""

    return build(context, parent, signed_correction, error_only, **kwargs)


def assemble(context: Any, parent: Any = None, signed_correction: Any = None, error_only: Any = None, **kwargs: Any) -> dict[str, Any]:
    """Alias for :func:`build`."""

    return build(context, parent, signed_correction, error_only, **kwargs)


def combine(context: Any, parent: Any = None, signed_correction: Any = None, error_only: Any = None, **kwargs: Any) -> dict[str, Any]:
    """Alias for :func:`build`."""

    return build(context, parent, signed_correction, error_only, **kwargs)


def construct(context: Any, parent: Any = None, signed_correction: Any = None, error_only: Any = None, **kwargs: Any) -> dict[str, Any]:
    """Alias for :func:`build`."""

    return build(context, parent, signed_correction, error_only, **kwargs)


def enclose(context: Any, parent: Any = None, signed_correction: Any = None, error_only: Any = None, **kwargs: Any) -> dict[str, Any]:
    """Alias for :func:`build`."""

    return build(context, parent, signed_correction, error_only, **kwargs)
