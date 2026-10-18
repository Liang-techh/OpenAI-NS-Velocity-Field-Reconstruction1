"""Source-only consumer for a corrected reserved ``m15`` endpoint graph.

The consumer accepts already materialised graph rows and endpoint handles.  It
does not load a saved report, import an ancestor producer, evaluate a source
function, or choose a numerical enclosure.  A corrected endpoint patch is
written in the namespace immediately after ``base_nodes``; ``construct``
relocates that patch after ``full_nodes`` and appends the small signed source
view needed by the reserved/full-pulse comparison.

Only graph rows are transferred.  The ordinary-C1 source formula and the
positive-unit receipt are returned as metadata, so a caller can bind them to a
later physical consumer without treating this source view as a numerical
claim.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from copy import deepcopy
import hashlib
import json
from typing import Any


_SOURCE_PAIR_NAMES = ("inner", "lower", "preRp", "pulse", "full")
_PATCH_OPERATIONS = {
    "sum": ("arguments",),
    "product": ("arguments",),
    "negative": ("argument",),
    "positive_quotient": ("numerator", "denominator"),
}


def canonical_digest(value: Any) -> str:
    """Return the canonical JSON SHA-256 digest used by this source view."""

    encoded = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _sequence(value: Any, name: str) -> Sequence[Any]:
    if isinstance(value, (str, bytes, bytearray)) or not isinstance(value, Sequence):
        raise TypeError(f"{name} must be a sequence")
    return value


def _nodes(value: Any, name: str) -> list[dict[str, Any]]:
    rows = _sequence(value, name)
    result: list[dict[str, Any]] = []
    for index, row in enumerate(rows):
        if not isinstance(row, Mapping):
            raise TypeError(f"{name}[{index}] must be a mapping")
        result.append(deepcopy(dict(row)))
    return result


def _handle(value: Any, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{name} must be an integer graph handle")
    if value < 0:
        raise ValueError(f"{name} must be nonnegative")
    return value


def _pair(value: Any, name: str, limit: int) -> tuple[int, int]:
    values = _sequence(value, name)
    if len(values) != 2:
        raise ValueError(f"{name} must contain ordinary Z0 and Z1")
    result = (_handle(values[0], f"{name}[0]"), _handle(values[1], f"{name}[1]"))
    if any(item >= limit for item in result):
        raise ValueError(f"{name} contains a handle outside the graph")
    return result


def _patch_references(row: Mapping[str, Any], index: int) -> list[tuple[str, int]]:
    operation = row.get("operation")
    if operation not in _PATCH_OPERATIONS:
        raise ValueError(
            f"patch row {index} uses unsupported operation {operation!r}; "
            "only sum, product, negative, and positive_quotient are allowed"
        )
    references: list[tuple[str, int]] = []
    for field in _PATCH_OPERATIONS[operation]:
        if field not in row:
            raise ValueError(f"patch row {index} ({operation}) is missing {field}")
        value = row[field]
        if field == "arguments":
            values = _sequence(value, f"patch row {index} arguments")
            if not values:
                raise ValueError(f"patch row {index} arguments cannot be empty")
        else:
            values = (value,)
        for offset, child in enumerate(values):
            references.append((f"{field}[{offset}]" if field == "arguments" else field,
                               _handle(child, f"patch row {index} {field}")))
    return references


def _relocate_patch(
    rows: Sequence[Mapping[str, Any]],
    *,
    old_base: int,
    new_base: int,
) -> tuple[list[dict[str, Any]], dict[int, int]]:
    """Copy and rebase a patch while rejecting forward/out-of-range edges."""

    suffix_count = len(rows)
    relocated: list[dict[str, Any]] = []
    mapping: dict[int, int] = {}
    for offset, source_row in enumerate(rows):
        row = deepcopy(dict(source_row))
        references = _patch_references(row, offset)
        mapped_by_field: dict[str, list[int]] = {}
        for field_name, child in references:
            if child < old_base:
                # The prefix check in construct has established that this
                # handle has exactly the same row in both namespaces.
                mapped = child
            else:
                if child >= old_base + suffix_count:
                    raise ValueError(
                        f"patch row {offset} {field_name} points outside appended range"
                    )
                if child >= old_base + offset:
                    raise ValueError(
                        f"patch row {offset} {field_name} points forward or to itself"
                    )
                mapped = new_base + child - old_base
            mapped_by_field.setdefault(field_name.split("[", 1)[0], []).append(mapped)

        for field in _PATCH_OPERATIONS[row["operation"]]:
            if field == "arguments":
                row[field] = mapped_by_field[field]
            else:
                row[field] = mapped_by_field[field][0]

        relocated.append(row)
        mapping[old_base + offset] = new_base + offset

    return relocated, mapping


def _append_delta_row(
    delta: list[dict[str, Any]],
    row: Mapping[str, Any],
    *,
    full_node_count: int,
) -> int:
    """Append one local row, checking that every graph edge is backward."""

    copied = deepcopy(dict(row))
    index = full_node_count + len(delta)
    references = _patch_references(copied, len(delta))
    for field_name, child in references:
        if child >= index:
            raise ValueError(f"derived row {len(delta)} {field_name} points forward")
    delta.append(copied)
    return index


def _append_pair(
    delta: list[dict[str, Any]],
    rows: Sequence[Mapping[str, Any]],
    *,
    full_node_count: int,
) -> list[int]:
    return [
        _append_delta_row(delta, row, full_node_count=full_node_count)
        for row in rows
    ]


def _verify_full_roots(
    full_nodes: Sequence[Mapping[str, Any]],
    source_pairs: Mapping[str, tuple[int, int]],
) -> None:
    """Check only the four physical ports at each original full root.

    Saved full roots may use a short sum tree rather than one flat row.  The
    walk is deliberately bounded and sum-only: it stops at the four supplied
    physical ports, rejects any other non-sum node, and never performs a
    general graph reduction.
    """

    components = ("inner", "lower", "preRp", "pulse")
    for order, root in enumerate(source_pairs["full"]):
        row = full_nodes[root]
        if row.get("operation") != "sum":
            raise ValueError(f"original full root at ordinary Z{order} is not a sum")
        expected = [source_pairs[name][order] for name in components]
        expected_set = set(expected)

        def flatten(index: int, depth: int, path: frozenset[int]) -> list[int]:
            if index in path:
                raise ValueError(f"cycle in original full root at ordinary Z{order}")
            # A component port is an allowed leaf even when its own row is a
            # sum.  This is the stopping rule for the bounded source walk.
            if index in expected_set:
                return [index]
            if depth >= 8:
                raise ValueError(
                    f"original full root at ordinary Z{order} exceeds bounded sum depth"
                )
            current = full_nodes[index]
            if current.get("operation") != "sum":
                raise ValueError(
                    f"non-sum node before a physical port at ordinary Z{order}"
                )
            arguments = current.get("arguments")
            if not isinstance(arguments, (list, tuple)) or not arguments:
                raise ValueError(
                    f"sum node in original full root at ordinary Z{order} has no arguments"
                )
            children = [
                _handle(value, f"full root Z{order} argument") for value in arguments
            ]
            if any(child >= len(full_nodes) for child in children):
                raise ValueError(
                    f"original full root at ordinary Z{order} points outside full_nodes"
                )
            if len(children) > 4:
                raise ValueError(
                    f"original full root at ordinary Z{order} has too many leaves"
                )
            next_path = path | {index}
            leaves: list[int] = []
            for child in children:
                leaves.extend(flatten(child, depth + 1, next_path))
                if len(leaves) > 4:
                    raise ValueError(
                        f"original full root at ordinary Z{order} has too many leaves"
                    )
            return leaves

        handles = flatten(root, 0, frozenset())
        if len(handles) != 4 or sorted(handles) != sorted(expected):
            raise ValueError(
                f"original full root at ordinary Z{order} does not contain the four ports"
            )


def _verify_unit_receipt(unit_receipt: Any) -> dict[str, Any]:
    if not isinstance(unit_receipt, Mapping):
        raise TypeError("unit_receipt must be a mapping")
    for key in ("positive", "Z_independent", "ordinary_Z_derivative"):
        if unit_receipt.get(key) is not True:
            raise ValueError(f"unit_receipt must prove {key}=True")
    return deepcopy(dict(unit_receipt))


def construct(
    base_nodes: Sequence[Mapping[str, Any]],
    full_nodes: Sequence[Mapping[str, Any]],
    patch: Mapping[str, Any],
    source_pairs: Mapping[str, Sequence[int]],
    *,
    expected_base_digest: str,
    expected_full_digest: str,
    unit_receipt: Mapping[str, Any],
) -> dict[str, Any]:
    """Consume one corrected endpoint graph in a fresh namespace.

    ``base_nodes`` and ``full_nodes`` are caller-owned snapshots.  The patch
    rows use the namespace whose first appended handle is ``len(base_nodes)``;
    the returned ``delta_nodes`` use the namespace whose first appended handle
    is ``len(full_nodes)``.  Existing rows are never returned as part of the
    delta and are never mutated.
    """

    # Materialise private copies before doing any append or field rewrite.
    base = _nodes(base_nodes, "base_nodes")
    full = _nodes(full_nodes, "full_nodes")
    base_count = len(base)
    full_count = len(full)
    if full_count < base_count:
        raise ValueError("full_nodes is shorter than base_nodes")

    actual_base_digest = canonical_digest(base)
    actual_full_digest = canonical_digest(full)
    if expected_base_digest != actual_base_digest:
        raise ValueError("base_nodes canonical digest does not match expected_base_digest")
    if expected_full_digest != actual_full_digest:
        raise ValueError("full_nodes canonical digest does not match expected_full_digest")

    # This exact comparison is deliberately done before any namespace join.
    if full[:base_count] != base:
        raise ValueError("full_nodes does not retain base_nodes as an exact prefix")

    if not isinstance(patch, Mapping):
        raise TypeError("patch must be a mapping")
    if patch.get("base_graph_node_count") != base_count:
        raise ValueError("patch base_graph_node_count does not equal len(base_nodes)")
    patch_rows_value = patch.get("appended_nodes")
    patch_rows = _sequence(patch_rows_value, "patch appended_nodes")
    if any(not isinstance(row, Mapping) for row in patch_rows):
        raise TypeError("patch appended_nodes must contain mappings")

    if not isinstance(source_pairs, Mapping):
        raise TypeError("source_pairs must be a mapping")
    if set(source_pairs) != set(_SOURCE_PAIR_NAMES):
        raise ValueError(
            "source_pairs must contain exactly inner, lower, preRp, pulse, and full"
        )
    pairs = {
        name: _pair(source_pairs[name], f"source_pairs.{name}", full_count)
        for name in _SOURCE_PAIR_NAMES
    }

    original_pairs = patch.get("original_pairs")
    corrected_pairs = patch.get("corrected_pairs")
    if not isinstance(original_pairs, Mapping) or not isinstance(corrected_pairs, Mapping):
        raise TypeError("patch original_pairs and corrected_pairs must be mappings")
    original_pre_rp = _pair(
        original_pairs.get("preRp_m15"),
        "patch.original_pairs.preRp_m15",
        full_count,
    )
    if original_pre_rp != pairs["preRp"]:
        raise ValueError("patch original_pairs.preRp_m15 differs from source_pairs.preRp")

    corrected_pre_rp_old = _pair(
        corrected_pairs.get("preRp_m15"),
        "patch.corrected_pairs.preRp_m15",
        base_count + len(patch_rows),
    )
    if any(
        value < base_count or value >= base_count + len(patch_rows)
        for value in corrected_pre_rp_old
    ):
        raise ValueError("patch corrected_pairs.preRp_m15 must point into appended_nodes")

    _verify_full_roots(full, pairs)
    receipt = _verify_unit_receipt(unit_receipt)

    relocated_patch, mapping = _relocate_patch(
        patch_rows,
        old_base=base_count,
        new_base=full_count,
    )
    corrected_pre_rp = tuple(mapping[value] for value in corrected_pre_rp_old)

    delta: list[dict[str, Any]] = list(relocated_patch)
    # Relocated rows were validated against their old namespace; check their
    # new namespace edges as a second guard before appending derived rows.
    for offset, row in enumerate(delta):
        for field_name, child in _patch_references(row, offset):
            if child >= full_count + offset:
                raise ValueError(
                    f"relocated row {offset} {field_name} points forward after rebasing"
                )

    def append_pair(rows: Sequence[Mapping[str, Any]]) -> list[int]:
        return _append_pair(delta, rows, full_node_count=full_count)

    corrected_full = append_pair(
        [
            {
                "operation": "sum",
                "arguments": [
                    pairs["inner"][order],
                    pairs["lower"][order],
                    corrected_pre_rp[order],
                    pairs["pulse"][order],
                ],
            }
            for order in range(2)
        ]
    )
    corrected_reserved = append_pair(
        [
            {"operation": "negative", "argument": corrected_full[order]}
            for order in range(2)
        ]
    )
    independent_target = append_pair(
        [
            {"operation": "negative", "argument": pairs["pulse"][order]}
            for order in range(2)
        ]
    )
    correlated_sum = append_pair(
        [
            {
                "operation": "sum",
                "arguments": [
                    pairs["inner"][order],
                    pairs["lower"][order],
                    corrected_pre_rp[order],
                ],
            }
            for order in range(2)
        ]
    )
    correlated_discrepancy = append_pair(
        [
            {"operation": "negative", "argument": correlated_sum[order]}
            for order in range(2)
        ]
    )
    correction_negative = append_pair(
        [
            {"operation": "negative", "argument": corrected_pre_rp[order]}
            for order in range(2)
        ]
    )
    correction = append_pair(
        [
            {
                "operation": "sum",
                "arguments": [pairs["preRp"][order], correction_negative[order]],
            }
            for order in range(2)
        ]
    )

    physical_c1_new_roots = {
        "corrected_preRp_m15": list(corrected_pre_rp),
        "corrected_full_m15": corrected_full,
        "corrected_reserved_m15": corrected_reserved,
        "independent_target_m15": independent_target,
        "correlated_discrepancy_m15": correlated_discrepancy,
        "correction_m15": correction,
    }
    # Build the final namespace only to certify append references and its
    # resulting digest.  It is intentionally not returned.
    joined = full + deepcopy(delta)
    final_digest = canonical_digest(joined)
    source_formula = {
        "name": "Dcorr",
        "formula": "Dcorr=Dold+(oldpreRp-correctedpreRp)/U",
        "sign": "additive",
        "ordinary_Z_orders": [0, 1],
        "old_preRp_m15": list(pairs["preRp"]),
        "corrected_preRp_m15": list(corrected_pre_rp),
        "correction_m15": correction,
        "unit_symbol": "U",
    }
    graph_patch = {
        "base_graph_node_count": full_count,
        "appended_nodes": deepcopy(delta),
        "appended_node_count": len(delta),
        "full_node_count": len(joined),
        "full_canonical_nodes_sha256": final_digest,
    }
    return {
        "passed": True,
        "source_only": True,
        "source_executed": False,
        "accepted_ancestors_replayed": False,
        "delta_nodes": deepcopy(delta),
        "delta_node_count": len(delta),
        "graph_patch": graph_patch,
        "base_prefix": {"node_count": base_count, "sha256": actual_base_digest},
        "full_prefix": {"node_count": full_count, "sha256": actual_full_digest},
        "full_prefix_node_count": full_count,
        "full_prefix_sha256": actual_full_digest,
        "mapping": dict(mapping),
        "relocated_pairs": {"preRp_m15": list(corrected_pre_rp)},
        "physical_C1_new_roots": physical_c1_new_roots,
        "source_formula": source_formula,
        "unit_receipt": receipt,
        "numerical_bound": False,
        "global_closure": False,
        "existing_ledger_edited": False,
        "own_swirl_added": False,
    }


__all__ = ["canonical_digest", "construct"]
