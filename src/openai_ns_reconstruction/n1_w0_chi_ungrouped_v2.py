"""Direct formal-expression leaf for the ungrouped ``W0_chi`` attachment.

The three inputs to :func:`build` are authenticated packets.  This module
does not produce an interval, add coefficients, replay a parent builder, or
choose a point value from an interval.  It only records the three streams as
separate formal segments, retaining every input term by reference.

``parent`` is the prior ordinary-Z formal packet.  ``signed`` contains the
raw ``W0_chi`` signed rows and ``errors`` contains the direct ``j >= 2``
kernel-tail rows.  The latter is an error-only stream; it is deliberately not
inserted into the signed expression and is not represented as a total error.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any


SCHEMA = "n1_w0_chi_ungrouped_v2"
BRANCH = "W0_chi"
CORRECTION_MARKER = BRANCH
ATTACHMENT_MARKER = SCHEMA
ORDINARY_Z_ORDERS = (0, 1)
SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
FORBIDDEN_LABELS = frozenset(("e01", "e10", "e11", "p1"))

__all__ = [
    "SCHEMA",
    "BRANCH",
    "CORRECTION_MARKER",
    "ATTACHMENT_MARKER",
    "ORDINARY_Z_ORDERS",
    "SOURCE_POWER_BASIS",
    "build",
]


def _is_sequence(value: Any) -> bool:
    return isinstance(value, Sequence) and not isinstance(
        value, (str, bytes, bytearray)
    )


def _order_key(value: Any) -> int | None:
    if isinstance(value, bool):
        return None
    if isinstance(value, int):
        return value if value in ORDINARY_Z_ORDERS else None
    if isinstance(value, str):
        text = value.strip().lower()
        if text in {"0", "z0"}:
            return 0
        if text in {"1", "z1"}:
            return 1
    return None


def _rows(packet: Any, label: str) -> tuple[dict[int, Sequence[Any]], Any]:
    """Return ordinary-Z rows without copying the caller's row objects."""

    original = packet
    if isinstance(packet, Mapping):
        found: dict[int, Sequence[Any]] = {}
        for raw_order, row in packet.items():
            order = _order_key(raw_order)
            if order is None:
                continue
            if order in found:
                raise ValueError(f"{label} contains duplicate ordinary-Z order {order}")
            if not _is_sequence(row):
                raise TypeError(f"{label}[Z{order}] must be a term sequence")
            if not row:
                raise ValueError(f"{label}[Z{order}] must not be empty")
            found[order] = row
        if set(found) == set(ORDINARY_Z_ORDERS):
            return {order: found[order] for order in ORDINARY_Z_ORDERS}, original

        # Saved receipts sometimes put the direct packet under a channel
        # name.  Unwrap only a direct order packet; no producer object is
        # called and no term is expanded.
        for key in (
            "formal_sum",
            "prior_formal_sum",
            "prior",
            "rows",
            "signed",
            "signed_rows",
            "error_only",
            "errors",
            "error_rows",
            "tail",
        ):
            if key in packet and packet[key] is not packet:
                try:
                    return _rows(packet[key], f"{label}.{key}")
                except (TypeError, ValueError):
                    continue
        raise ValueError(f"{label} must provide ordinary-Z rows 0 and 1")

    if _is_sequence(packet) and len(packet) == len(ORDINARY_Z_ORDERS):
        result: dict[int, Sequence[Any]] = {}
        for order, row in zip(ORDINARY_Z_ORDERS, packet):
            if not _is_sequence(row):
                raise TypeError(f"{label}[Z{order}] must be a term sequence")
            if not row:
                raise ValueError(f"{label}[Z{order}] must not be empty")
            result[order] = row
        return result, original

    raise TypeError(f"{label} must be an ordinary-Z packet")


def _label_tokens(value: Any) -> set[str]:
    """Extract conservative marker tokens from provenance text."""

    if not isinstance(value, str):
        return set()
    text = value.strip().lower()
    for separator in ("-", ":", "/", "\\", ".", " "):
        text = text.replace(separator, "_")
    return {token for token in text.split("_") if token}


def _has_forbidden_label(value: Any, seen: set[int] | None = None) -> bool:
    if seen is None:
        seen = set()
    if isinstance(value, str):
        return bool(_label_tokens(value) & FORBIDDEN_LABELS)
    if isinstance(value, Mapping):
        identity = id(value)
        if identity in seen:
            return False
        seen.add(identity)
        return any(
            _has_forbidden_label(key, seen) or _has_forbidden_label(item, seen)
            for key, item in value.items()
        )
    if _is_sequence(value):
        identity = id(value)
        if identity in seen:
            return False
        seen.add(identity)
        return any(_has_forbidden_label(item, seen) for item in value)
    return False


def _is_marker_value(value: Any) -> bool:
    if isinstance(value, str):
        return value.strip().lower() in {
            CORRECTION_MARKER.lower(),
            ATTACHMENT_MARKER.lower(),
        }
    if isinstance(value, bool):
        return value
    if isinstance(value, Mapping):
        return any(_is_marker_value(item) for item in value.values())
    if _is_sequence(value):
        return any(_is_marker_value(item) for item in value)
    return False


def _already_attached(value: Any, seen: set[int] | None = None) -> bool:
    """Find this leaf's marker in a parent without inspecting coefficients."""

    if seen is None:
        seen = set()
    if isinstance(value, str):
        return value.strip().lower() in {
            CORRECTION_MARKER.lower(),
            ATTACHMENT_MARKER.lower(),
        }
    if isinstance(value, Mapping):
        identity = id(value)
        if identity in seen:
            return False
        seen.add(identity)
        for key, item in value.items():
            key_text = str(key).strip().lower()
            if key_text in {
                CORRECTION_MARKER.lower(),
                ATTACHMENT_MARKER.lower(),
            }:
                return True
            if key_text in {
                "correction_marker",
                "attachment_marker",
                "attached_correction",
                "attached_corrections",
                "correction_attached",
                "w0_chi_attached",
            } and _is_marker_value(item):
                return True
            if _already_attached(item, seen):
                return True
        return False
    if _is_sequence(value):
        identity = id(value)
        if identity in seen:
            return False
        seen.add(identity)
        return any(_already_attached(item, seen) for item in value)
    return False


def _require_term(raw: Any, label: str, stream: str) -> Mapping[str, Any]:
    if not isinstance(raw, Mapping):
        raise TypeError(f"{label} must be a term mapping")
    for field in ("powers", "coefficient", "origin"):
        if field not in raw:
            raise ValueError(f"{label} requires {field}")
    powers = raw["powers"]
    if not _is_sequence(powers) or len(powers) != 4 or any(type(power) is not int for power in powers):
        raise ValueError(f"{label}.powers must contain four integers")
    if not hasattr(raw["coefficient"], "_mpi_"):
        raise TypeError(f"{label}.coefficient must be a caller-supplied interval")
    origin = raw["origin"]
    if not isinstance(origin, Mapping):
        raise TypeError(f"{label}.origin must be a mapping")
    if stream != "parent":
        for field in ("branch", "word", "channel"):
            value = raw.get(field)
            if isinstance(value, str) and value.strip().lower() in FORBIDDEN_LABELS:
                raise ValueError(f"{label} contains a forbidden branch label")

    if stream == "signed":
        if raw.get("signed_field") is not True:
            raise ValueError(f"{label} must have signed_field=True")
        if raw.get("error_only") is True or origin.get("error_only") is True:
            raise ValueError(f"{label} may not be error_only")
        if raw.get("symmetric_error") is True or origin.get("symmetric_error") is True:
            raise ValueError(f"{label} may not be a symmetric error")
        if "signed_field" in origin and origin["signed_field"] is not True:
            raise ValueError(f"{label}.origin.signed_field must be True")
    elif stream == "error":
        if origin.get("error_only") is not True:
            raise ValueError(f"{label}.origin.error_only must be True")
        if "error_only" in raw and raw["error_only"] is not True:
            raise ValueError(f"{label}.error_only must be True when supplied")
        if raw.get("signed_field") is True or origin.get("signed_field") is True:
            raise ValueError(f"{label} may not have signed_field=True")
        signed_markers = []
        if "signed_field" in raw:
            signed_markers.append(raw["signed_field"])
        if "signed_field" in origin:
            signed_markers.append(origin["signed_field"])
        if not signed_markers or not all(marker is False for marker in signed_markers):
            raise ValueError(f"{label} must have signed_field=False")
        if not (raw.get("symmetric_error") is True or origin.get("symmetric_error") is True or origin.get("bound_type") == "symmetric_error"):
            raise ValueError(f"{label} must retain the symmetric-error marker")
    elif stream != "parent":
        raise ValueError(f"unknown stream {stream!r}")
    return raw


def _field(raw: Mapping[str, Any], origin: Mapping[str, Any], names: Sequence[str]) -> Any:
    for mapping in (raw, origin):
        for name in names:
            if name in mapping:
                return mapping[name]
    return None


def _term_reference(
    raw: Mapping[str, Any],
    order: int,
    index: int,
    kind: str,
) -> dict[str, Any]:
    """Describe a raw term while retaining every payload field by reference."""

    origin = raw["origin"]
    path = _field(
        raw,
        origin,
        ("path", "source_path", "input_path", "source_paths", "paths", "input_paths"),
    )
    source_ids = _field(
        raw,
        origin,
        ("source_ids", "source_id", "sources", "source", "input_id"),
    )
    error_classes = _field(
        raw,
        origin,
        ("error_classes", "error_class", "source_error_classes", "error_category", "error_categories"),
    )
    return {
        "kind": kind,
        "segment": kind,
        "ordinary_Z_order": order,
        "index": index,
        "term": raw,
        "raw_term": raw,
        "input_term": raw,
        "powers": raw["powers"],
        "coefficient": raw["coefficient"],
        "origin": origin,
        "path": path,
        "source_ids": source_ids,
        "error_classes": error_classes,
        "by_reference": True,
    }


def _segment(
    rows: Mapping[int, Sequence[Any]],
    order: int,
    kind: str,
) -> dict[str, Any]:
    refs = [
        _term_reference(raw, order, index, kind)
        for index, raw in enumerate(rows[order])
    ]
    return {
        "kind": kind,
        "ordinary_Z_order": order,
        "terms": rows[order],
        "raw_terms": rows[order],
        "term_references": refs,
        "term_refs": refs,
        "count": len(rows[order]),
        "direct": True,
        "by_reference": True,
    }


def build(parent: Any, signed: Any, errors: Any) -> dict[str, Any]:
    """Return a direct, ungrouped formal expression.

    ``parent``, ``signed`` and ``errors`` are caller-authenticated interval
    arrays.  Their coefficient objects are never inspected or operated on;
    only structural markers are checked and all raw terms are retained.
    """

    if _already_attached(parent):
        raise ValueError("parent already contains the W0_chi correction marker")

    parent_rows, parent_input = _rows(parent, "parent")
    signed_rows, signed_input = _rows(signed, "signed")
    error_rows, errors_input = _rows(errors, "errors")

    for order in ORDINARY_Z_ORDERS:
        for index, raw in enumerate(parent_rows[order]):
            _require_term(raw, f"parent[Z{order}][{index}]", "parent")
        for index, raw in enumerate(signed_rows[order]):
            _require_term(raw, f"signed[Z{order}][{index}]", "signed")
        for index, raw in enumerate(error_rows[order]):
            _require_term(raw, f"errors[Z{order}][{index}]", "error")

    prior_segments = {
        order: _segment(parent_rows, order, "prior")
        for order in ORDINARY_Z_ORDERS
    }
    signed_segments = {
        order: _segment(signed_rows, order, "signed")
        for order in ORDINARY_Z_ORDERS
    }
    error_segments = {
        order: _segment(error_rows, order, "error_only")
        for order in ORDINARY_Z_ORDERS
    }

    # The expression itself contains raw rows in separate channels.  The
    # detailed segment records are supplied alongside it so consumers can
    # inspect powers, coefficients, origins and paths without reconstituting
    # or regrouping any term.
    formal_expression = {
        order: {
            "prior": parent_rows[order],
            "signed": signed_rows[order],
            "prior_terms": parent_rows[order],
            "signed_terms": signed_rows[order],
            "segments": {
                "prior": prior_segments[order],
                "signed": signed_segments[order],
            },
        }
        for order in ORDINARY_Z_ORDERS
    }
    formal_segments = {
        order: {
            "prior": prior_segments[order],
            "signed": signed_segments[order],
            "error_only": error_segments[order],
        }
        for order in ORDINARY_Z_ORDERS
    }

    return {
        "schema": SCHEMA,
        "branch": BRANCH,
        "correction_marker": CORRECTION_MARKER,
        "attachment_marker": ATTACHMENT_MARKER,
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "source_power_basis": list(SOURCE_POWER_BASIS),
        "parent_input": parent_input,
        "signed_input": signed_input,
        "errors_input": errors_input,
        "prior_formal_sum": parent_rows,
        "prior": parent_rows,
        "parent": parent_rows,
        "signed": signed_rows,
        "signed_rows": signed_rows,
        "W0_chi": {
            "branch": BRANCH,
            "signed": signed_rows,
            "signed_rows": signed_rows,
        },
        "errors": error_rows,
        "error_only": error_rows,
        "error_rows": error_rows,
        "error_stream": error_rows,
        "error_tail": error_rows,
        "direct_j_ge_2_kernel_tail": {
            "kind": "direct_j_ge_2_kernel_tail",
            "rows": error_rows,
            "terms": error_rows,
            "error_only": True,
            "is_total_error": False,
        },
        "formal_expression": formal_expression,
        "formal_expression_by_order": formal_expression,
        "formal_segments": formal_segments,
        "prior_segments": prior_segments,
        "signed_segments": signed_segments,
        "error_segments": error_segments,
        "prior_counts": {
            0: len(parent_rows[0]),
            1: len(parent_rows[1]),
        },
        "signed_counts": {
            0: len(signed_rows[0]),
            1: len(signed_rows[1]),
        },
        "error_counts": {
            0: len(error_rows[0]),
            1: len(error_rows[1]),
        },
        "correction_attached": True,
        "duplicate_marker_rejected": True,
        "source_level_duplicate_not_proved": True,
        "branch_binding_is_caller_obligation": True,
        "input_terms_retained_by_reference": True,
        "prior_and_signed_segments_separate": True,
        "error_terms_in_signed_expression": False,
        "error_terms_are_total_error": False,
        "error_terms_are_independent_references": False,
        "error_terms_are_separate_references": True,
        "error_independence_not_proved": True,
        "coefficients_operated_on": False,
        "numeric_addition_performed": False,
        "midpoint_selected": False,
        "point_value_selected": False,
        "physical_weight_reapplied": False,
    }

