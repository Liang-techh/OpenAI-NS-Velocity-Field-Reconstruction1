"""Conditional cell-local sums for already integrated finite ``P`` sectors.

The caller supplies the four physical ``+2 chi E51`` finite-swirl branches
after their physical scale and measure have already been applied.  The input
shape is ``branch -> ordinary-Z order -> target cell -> term list``.  This
leaf does not call a producer, re-integrate a kernel, or apply another scale.

Only terms with the same ``(powers, uncertainty, ordered error_classes,
lead_order, word_order)`` are added, and the addition is never allowed to
cross an ordinary-Z row, target cell, or branch.  Every input list and term is
retained by reference in the result.  Error-only records are retained in a
separate stream and are never treated as signed finite terms.

Signedness is intentionally conditional.  A missing signed marker is
accepted only because the exact finite-P branch roster is authenticated by
the caller; an explicit caller obligation may additionally state that source
authentication and one common directed-interval context were used.  This
module does not prove those obligations, same-source or same-gauge identity,
physical completeness, an actual field, or independent error certification.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any


SCHEMA = "n1_p_compensation_cell_sum_v2"
SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
BRANCHES = (
    "G0_P_W0",
    "G0_P_DeltaW0",
    "D_P_W0",
    "D_P_DeltaW0",
)
ORDINARY_Z_ORDERS = (0, 1)
TARGET_CELLS = tuple(range(8))
GROUPING_KEY = (
    "powers",
    "uncertainty",
    "ordered_error_classes",
    "lead_order",
    "word_order",
)

__all__ = [
    "SCHEMA",
    "SOURCE_POWER_BASIS",
    "BRANCHES",
    "ORDINARY_Z_ORDERS",
    "TARGET_CELLS",
    "GROUPING_KEY",
    "build",
    "sum_compensation_cells",
    "cell_sum",
    "aggregate",
    "combine",
]


def _is_sequence(value: Any) -> bool:
    return isinstance(value, Sequence) and not isinstance(
        value, (str, bytes, bytearray)
    )


def _is_context(value: Any) -> bool:
    return hasattr(value, "mpf") and hasattr(value, "make_mpf")


def _is_interval(value: Any) -> bool:
    return (
        not isinstance(value, (bool, int, float, complex, str, bytes, bytearray))
        and hasattr(value, "_mpi_")
        and hasattr(value, "a")
        and hasattr(value, "b")
    )


def _typed(value: Any, label: str) -> tuple[type[Any], Any]:
    """Return a type-sensitive hash key for provenance labels."""

    try:
        hash(value)
    except TypeError as exc:
        raise TypeError(f"{label} must be hashable") from exc
    return type(value), value


def _unique(values: Sequence[Any], label: str) -> list[Any]:
    result: list[Any] = []
    seen: set[tuple[type[Any], Any]] = set()
    for index, value in enumerate(values):
        key = _typed(value, f"{label}[{index}]")
        if key not in seen:
            seen.add(key)
            result.append(value)
    return result


def _order(raw: Any, label: str) -> int:
    if isinstance(raw, bool):
        raise TypeError(f"{label} must be ordinary-Z order 0 or 1")
    if isinstance(raw, int):
        value = raw
    elif isinstance(raw, str) and raw.strip().upper() in {"0", "Z0", "1", "Z1"}:
        value = int(raw.strip().upper().removeprefix("Z"))
    else:
        raise TypeError(f"{label} must be ordinary-Z order 0 or 1")
    if value not in ORDINARY_Z_ORDERS:
        raise ValueError(f"{label} must be ordinary-Z order 0 or 1")
    return value


def _cell(raw: Any, label: str) -> int:
    if isinstance(raw, bool):
        raise TypeError(f"{label} must be target cell 0 through 7")
    if isinstance(raw, int):
        value = raw
    elif isinstance(raw, str) and raw.strip().isdigit():
        value = int(raw.strip())
    else:
        raise TypeError(f"{label} must be target cell 0 through 7")
    if value not in TARGET_CELLS:
        raise ValueError(f"{label} must be target cell 0 through 7")
    return value


def _items(value: Any, label: str, *, allow_empty: bool = True) -> list[Any] | Sequence[Any]:
    """Return a term sequence without copying the caller's sequence object."""

    if not _is_sequence(value):
        raise TypeError(f"{label} must be a term list")
    if not allow_empty and not value:
        raise ValueError(f"{label} must not be empty")
    return value


def _unwrap(mapping: Mapping[Any, Any], names: Sequence[str], label: str) -> Any:
    found = [mapping[name] for name in names if name in mapping]
    if len(found) > 1 and any(item is not found[0] for item in found[1:]):
        raise ValueError(f"{label} contains conflicting aliases")
    return found[0] if found else None


def _mapping_orders(value: Any, label: str) -> dict[int, Any]:
    if isinstance(value, Mapping):
        rows: dict[int, Any] = {}
        for raw_order, row in value.items():
            order = _order(raw_order, f"{label} key")
            if order in rows:
                raise ValueError(f"{label} contains duplicate ordinary-Z order {order}")
            rows[order] = row
        if set(rows) != set(ORDINARY_Z_ORDERS):
            raise ValueError(f"{label} must provide ordinary-Z orders 0 and 1")
        return rows
    if _is_sequence(value) and len(value) == len(ORDINARY_Z_ORDERS):
        return {order: value[order] for order in ORDINARY_Z_ORDERS}
    raise TypeError(f"{label} must provide ordinary-Z orders 0 and 1")


def _mapping_cells(value: Any, label: str) -> dict[int, Any]:
    if isinstance(value, Mapping):
        cells: dict[int, Any] = {}
        for raw_cell, terms in value.items():
            cell = _cell(raw_cell, f"{label} key")
            if cell in cells:
                raise ValueError(f"{label} contains duplicate target cell {cell}")
            cells[cell] = terms
        if set(cells) != set(TARGET_CELLS):
            raise ValueError(f"{label} must provide target cells 0 through 7")
        return cells
    if _is_sequence(value) and len(value) == len(TARGET_CELLS):
        return {cell: value[cell] for cell in TARGET_CELLS}
    raise TypeError(f"{label} must provide target cells 0 through 7")


def _branch_rows(branch: Any, label: str) -> dict[int, dict[int, Sequence[Any]]]:
    if not isinstance(branch, Mapping):
        raise TypeError(f"{label} must be an ordinary-Z/cell mapping")
    nested = _unwrap(branch, ("ordinary_Z", "ordinary_z", "orders", "by_order"), label)
    if nested is not None:
        branch = nested
    orders = _mapping_orders(branch, label)
    return {
        order: {
            cell: _items(terms, f"{label}.Z{order}.cell{cell}")
            for cell, terms in _mapping_cells(orders[order], f"{label}.Z{order}").items()
        }
        for order in ORDINARY_Z_ORDERS
    }


def _branches(packet: Any) -> tuple[Mapping[Any, Any], dict[str, dict[int, dict[int, Sequence[Any]]]]]:
    if not isinstance(packet, Mapping):
        raise TypeError("branches must be a branch mapping")
    original = packet
    nested = packet.get("branches")
    if isinstance(nested, Mapping):
        packet = nested
    if set(packet) != set(BRANCHES):
        raise ValueError("branches must contain exactly the four finite P branches")
    rows = {
        branch: _branch_rows(packet[branch], f"branches[{branch!r}]")
        for branch in BRANCHES
    }
    return original, rows


def _sequence_field(value: Any, label: str, *, nonempty: bool = False) -> list[Any]:
    if not _is_sequence(value):
        raise TypeError(f"{label} must be a sequence")
    if nonempty and not value:
        raise ValueError(f"{label} must not be empty")
    # This is only metadata normalization; the complete origin mapping stays
    # available by reference on each member.
    return list(value)


def _field(mapping: Mapping[Any, Any], names: Sequence[str]) -> Any:
    for name in names:
        if name in mapping:
            return mapping[name]
    return None


def _obligation(value: Any, branch_names: Sequence[str]) -> dict[str, Any]:
    """Validate an optional caller assertion without proving it here."""

    if value is None:
        return {
            "provided": False,
            "authenticated_branches": (),
            "upstream_source_authenticated": None,
            "same_interval_context": None,
            "condition": (
                "missing signed markers are accepted only under the exact "
                "caller-authenticated finite-P branch roster"
            ),
        }
    if not isinstance(value, Mapping):
        raise TypeError("caller_obligation must be a mapping when supplied")

    raw_branches = _field(
        value,
        (
            "authenticated_branches",
            "authenticated_finite_P_branches",
            "finite_P_branches",
            "branch_roster",
        ),
    )
    if isinstance(raw_branches, Mapping):
        selected = [name for name in branch_names if raw_branches.get(name) is True]
    elif raw_branches is None:
        selected = list(branch_names)
    elif _is_sequence(raw_branches):
        selected = list(raw_branches)
    else:
        raise TypeError("caller_obligation authenticated branches must be a sequence or mapping")
    if set(selected) != set(branch_names):
        raise ValueError("caller_obligation must authenticate all four finite P branches")

    auth = _field(
        value,
        (
            "upstream_source_authenticated",
            "source_authenticated",
            "authenticated",
        ),
    )
    if auth is not None and auth is not True:
        raise ValueError("caller_obligation source authentication must be true")
    same = _field(
        value,
        (
            "same_interval_context",
            "same_mpmath_interval_context",
            "interval_context_same",
        ),
    )
    if same is not None and same is not True:
        raise ValueError("caller_obligation must state one common interval context")
    return {
        "provided": True,
        "authenticated_branches": tuple(selected),
        "upstream_source_authenticated": auth,
        "same_interval_context": same,
        "condition": value.get(
            "condition",
            "caller asserts upstream source authentication and one hydrated interval context",
        ),
        "record": value,
    }


def _origin(term: Mapping[str, Any], label: str, branch: str, order: int, cell: int) -> dict[str, Any]:
    raw = term.get("origin")
    if not isinstance(raw, Mapping):
        raise TypeError(f"{label}.origin must be a provenance mapping")
    required = {"source_ids", "error_classes", "uncertainty", "target_cell", "lead_order", "word_order"}
    missing = required.difference(raw)
    if missing:
        raise ValueError(f"{label}.origin is missing {', '.join(sorted(missing))}")

    source_ids = _sequence_field(raw["source_ids"], f"{label}.origin.source_ids", nonempty=True)
    _unique(source_ids, f"{label}.origin.source_ids")
    error_classes = _sequence_field(raw["error_classes"], f"{label}.origin.error_classes")
    # Error class order is semantic for this leaf; do not sort or deduplicate.
    for index, item in enumerate(error_classes):
        _typed(item, f"{label}.origin.error_classes[{index}]")
    uncertainty = raw["uncertainty"]
    if not isinstance(uncertainty, bool):
        raise TypeError(f"{label}.origin.uncertainty must be boolean")

    target = _cell(raw["target_cell"], f"{label}.origin.target_cell")
    if target != cell:
        raise ValueError(
            f"{label}.origin.target_cell {target} does not match target cell {cell}"
        )
    lead = raw["lead_order"]
    word = raw["word_order"]
    if type(lead) is not int or type(word) is not int or lead not in (0, 1) or word not in (0, 1) or lead + word != order:
        raise ValueError(f"{label} has incompatible ordinary-Z product orders")

    if "ordinary_Z_order" in raw and _order(raw["ordinary_Z_order"], f"{label}.origin.ordinary_Z_order") != order:
        raise ValueError(f"{label}.origin ordinary-Z order does not match its input row")
    for mapping, field in ((term, "branch"), (raw, "branch")):
        if field in mapping and mapping[field] != branch:
            raise ValueError(f"{label}.{field} does not match its input branch")
    if "target_cell" in term and _cell(term["target_cell"], f"{label}.target_cell") != cell:
        raise ValueError(f"{label}.target_cell does not match its input cell")

    error_only = term.get("error_only") is True or raw.get("error_only") is True
    symmetric = term.get("symmetric_error") is True or raw.get("symmetric_error") is True
    is_error = error_only or symmetric
    # A missing signed marker is covered only by the caller's source contract.
    for mapping, prefix in ((term, label), (raw, f"{label}.origin")):
        if "signed_field" in mapping and mapping["signed_field"] is not (not is_error):
            raise ValueError(f"{prefix}.signed_field disagrees with its stream")
        if "signed" in mapping and mapping["signed"] is not (not is_error):
            raise ValueError(f"{prefix}.signed disagrees with its stream")
        if "error_only" in mapping and mapping["error_only"] is not True and not isinstance(mapping["error_only"], bool):
            raise TypeError(f"{prefix}.error_only must be boolean")

    if is_error:
        return {
            "raw": raw,
            "source_ids": source_ids,
            "error_classes": error_classes,
            "uncertainty": uncertainty,
            "target_cell": target,
            "lead_order": lead,
            "word_order": word,
            "error_only": True,
            "source_path": _field(
                term, ("source_path", "source_paths", "input_path", "input_paths", "path", "paths")
            )
            if _field(term, ("source_path", "source_paths", "input_path", "input_paths", "path", "paths")) is not None
            else _field(raw, ("source_path", "source_paths", "input_path", "input_paths", "path", "paths")),
        }
    return {
        "raw": raw,
        "source_ids": source_ids,
        "error_classes": error_classes,
        "uncertainty": uncertainty,
        "target_cell": target,
        "lead_order": lead,
        "word_order": word,
        "error_only": False,
        "source_path": _field(
            term, ("source_path", "source_paths", "input_path", "input_paths", "path", "paths")
        )
        if _field(term, ("source_path", "source_paths", "input_path", "input_paths", "path", "paths")) is not None
        else _field(raw, ("source_path", "source_paths", "input_path", "input_paths", "path", "paths")),
    }


def _validate_term(
    context: Any,
    term: Any,
    label: str,
    branch: str,
    order: int,
    cell: int,
) -> tuple[dict[str, Any], bool]:
    if not isinstance(term, Mapping):
        raise TypeError(f"{label} must be a term mapping")
    if "powers" not in term:
        raise ValueError(f"{label} requires powers")
    powers = term["powers"]
    if not _is_sequence(powers) or len(powers) != 4:
        raise ValueError(f"{label}.powers must contain four integers")
    for index, power in enumerate(powers):
        if isinstance(power, bool) or not isinstance(power, int):
            raise TypeError(f"{label}.powers[{index}] must be an integer")
    coefficient = term.get("coefficient", term.get("coefficient_interval"))
    if coefficient is None or not _is_interval(coefficient):
        raise TypeError(f"{label}.coefficient must be a caller-supplied directed interval")
    # The context is intentionally not used to convert or rehydrate a value.
    # Its presence and the caller obligation are the only context contract.
    del context
    info = _origin(term, label, branch, order, cell)
    info["powers"] = tuple(powers)
    info["coefficient"] = coefficient
    return info, bool(info["error_only"])


def _member(
    term: Mapping[str, Any],
    origin: Mapping[str, Any],
    *,
    branch: str,
    order: int,
    cell: int,
    index: int,
) -> dict[str, Any]:
    return {
        "branch": branch,
        "ordinary_Z_order": order,
        "target_cell": cell,
        "lead_order": origin["lead_order"],
        "word_order": origin["word_order"],
        "input_index": index,
        "index": index,
        "record": term,
        "raw_record": term,
        "input_term": term,
        "source": term,
        "origin": term["origin"],
        "provenance": term["origin"],
        "source_path": origin["source_path"],
        "source_paths": origin["source_path"],
        "source_ids": origin["source_ids"],
        "error_classes": origin["error_classes"],
        "by_reference": True,
    }


def _aggregate(
    terms: Sequence[tuple[Mapping[str, Any], dict[str, Any], int]],
    *,
    branch: str,
    order: int,
    cell: int,
) -> list[dict[str, Any]]:
    groups: dict[Any, dict[str, Any]] = {}
    for term, info, index in terms:
        classes = tuple(info["error_classes"])
        key = (
            tuple(info["powers"]),
            bool(info["uncertainty"]),
            tuple(_typed(value, "error class") for value in classes),
            _typed(info["lead_order"], "lead_order"),
            _typed(info["word_order"], "word_order"),
        )
        group = groups.get(key)
        if group is None:
            group = {
                "powers": tuple(info["powers"]),
                "coefficient": info["coefficient"],
                "branch": branch,
                "ordinary_Z_order": order,
                "target_cell": cell,
                "lead_order": info["lead_order"],
                "word_order": info["word_order"],
                "origin": {
                    "operation": "conditional_directed_cell_sum",
                    "branch": branch,
                    "ordinary_Z_order": order,
                    "target_cell": cell,
                    "lead_order": info["lead_order"],
                    "word_order": info["word_order"],
                    "source_ids": [],
                    "source_paths": [],
                    "uncertainty": bool(info["uncertainty"]),
                    "error_classes": list(classes),
                    "input_indices": [],
                    "member_count": 0,
                    "symbolic_cancellation_claimed": False,
                },
                "members": [],
                "records": [],
                "input_indices": [],
                "group_key": (
                    tuple(info["powers"]),
                    bool(info["uncertainty"]),
                    classes,
                    info["lead_order"],
                    info["word_order"],
                ),
                "by_reference": True,
            }
            groups[key] = group
        else:
            # This is the only arithmetic in this leaf: directed addition of
            # equal-key terms that already belong to this one cell and row.
            group["coefficient"] = group["coefficient"] + info["coefficient"]
        member = _member(
            term,
            info,
            branch=branch,
            order=order,
            cell=cell,
            index=index,
        )
        group["members"].append(member)
        group["records"].append(term)
        group["input_indices"].append(index)
        metadata = group["origin"]
        metadata["input_indices"].append(index)
        metadata["member_count"] += 1
        metadata["source_ids"] = _unique(
            [*metadata["source_ids"], *info["source_ids"]],
            "group.origin.source_ids",
        )
        if info["source_path"] is not None:
            metadata["source_paths"].append(info["source_path"])
    return list(groups.values())


def build(
    context: Any,
    branches: Any = None,
    *,
    caller_obligation: Any = None,
    obligation: Any = None,
    caller_condition: Any = None,
    source_certificate: Any = None,
) -> dict[str, Any]:
    """Conditionally sum already integrated finite-P terms by branch/Z/cell.

    ``context`` is the caller's mpmath interval context.  ``branches`` is a
    mapping with exactly the four names in :data:`BRANCHES`; each branch must
    provide ordinary-Z rows ``0`` and ``1``, each with cells ``0`` through
    ``7``.  A cell's list is never copied into the output.  ``obligation`` is
    an alias for ``caller_obligation`` and, when supplied, must authenticate
    all four branches and one common interval context.
    """

    # Accept the natural reversed spelling as a small compatibility aid while
    # keeping the canonical API context-first.
    if not _is_context(context) and _is_context(branches):
        context, branches = branches, context
    if not _is_context(context):
        raise TypeError("context must be the caller's directed mpmath interval context")
    if branches is None:
        raise TypeError("branches are required")
    for alias, value in (
        ("obligation", obligation),
        ("caller_condition", caller_condition),
        ("source_certificate", source_certificate),
    ):
        if value is None:
            continue
        if caller_obligation is not None and caller_obligation is not value:
            raise ValueError(f"caller_obligation and {alias} disagree")
        caller_obligation = value

    input_packet, rows = _branches(branches)
    obligation_record = _obligation(caller_obligation, BRANCHES)

    branch_outputs: dict[str, Any] = {}
    all_groups: dict[str, dict[int, dict[int, list[dict[str, Any]]]]] = {}
    all_records: dict[str, dict[int, dict[int, Sequence[Any]]]] = {}
    all_errors: dict[str, dict[int, dict[int, list[Mapping[str, Any]]]]] = {}

    for branch in BRANCHES:
        order_outputs: dict[int, Any] = {}
        branch_groups: dict[int, dict[int, list[dict[str, Any]]]] = {}
        branch_records: dict[int, dict[int, Sequence[Any]]] = {}
        branch_errors: dict[int, dict[int, list[Mapping[str, Any]]]] = {}
        for order in ORDINARY_Z_ORDERS:
            cell_outputs: dict[int, Any] = {}
            cell_groups: dict[int, list[dict[str, Any]]] = {}
            cell_records: dict[int, Sequence[Any]] = {}
            cell_errors: dict[int, list[Mapping[str, Any]]] = {}
            for cell in TARGET_CELLS:
                raw_terms = rows[branch][order][cell]
                signed: list[tuple[Mapping[str, Any], dict[str, Any], int]] = []
                errors: list[Mapping[str, Any]] = []
                for index, raw_term in enumerate(raw_terms):
                    info, error_only = _validate_term(
                        context,
                        raw_term,
                        f"{branch}.Z{order}.cell{cell}[{index}]",
                        branch,
                        order,
                        cell,
                    )
                    if error_only:
                        errors.append(raw_term)
                    else:
                        signed.append((raw_term, info, index))
                groups = _aggregate(signed, branch=branch, order=order, cell=cell)
                cell_groups[cell] = groups
                cell_records[cell] = raw_terms
                cell_errors[cell] = errors
                cell_outputs[cell] = {
                    "branch": branch,
                    "ordinary_Z_order": order,
                    "target_cell": cell,
                    "records": raw_terms,
                    "raw_records": raw_terms,
                    "input_records": raw_terms,
                    "groups": groups,
                    "aggregates": groups,
                    "error_only_records": errors,
                    "signed_record_count": sum(len(group["members"]) for group in groups),
                    "error_only_count": len(errors),
                    "input_record_count": len(raw_terms),
                    "records_retained_by_reference": True,
                }
            order_outputs[order] = {
                "ordinary_Z_order": order,
                "cells": cell_outputs,
                "target_cells": cell_outputs,
                "groups": cell_groups,
                "aggregates": cell_groups,
                "records": cell_records,
                "raw_records": cell_records,
                "error_only_records": cell_errors,
            }
            branch_groups[order] = cell_groups
            branch_records[order] = cell_records
            branch_errors[order] = cell_errors
        branch_outputs[branch] = {
            "branch": branch,
            "ordinary_Z": order_outputs,
            "orders": order_outputs,
            "by_order": order_outputs,
            "groups": branch_groups,
            "aggregates": branch_groups,
            "records": branch_records,
            "raw_records": branch_records,
            "error_only_records": branch_errors,
            "branch_input": input_packet.get("branches", input_packet)[branch]
            if isinstance(input_packet.get("branches", input_packet), Mapping)
            else None,
        }
        all_groups[branch] = branch_groups
        all_records[branch] = branch_records
        all_errors[branch] = branch_errors

    return {
        "schema": SCHEMA,
        "source_power_basis": list(SOURCE_POWER_BASIS),
        "branches": branch_outputs,
        "branch_names": list(BRANCHES),
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "target_cells": list(TARGET_CELLS),
        "groups": all_groups,
        "aggregates": all_groups,
        "records": all_records,
        "raw_records": all_records,
        "error_only_records": all_errors,
        "input": input_packet,
        "caller_obligation": obligation_record.get("record", caller_obligation),
        "caller_obligation_record": obligation_record,
        "authenticated_branches_condition": obligation_record["authenticated_branches"],
        "upstream_source_authenticated": obligation_record["upstream_source_authenticated"],
        "same_interval_context": obligation_record["same_interval_context"],
        "grouping_key": GROUPING_KEY,
        "grouping_key_description": "(powers, uncertainty, ordered error_classes, lead_order, word_order)",
        "input_records_retained_by_reference": True,
        "member_records_retained_by_reference": True,
        "source_paths_retained": True,
        "provenance_retained": True,
        "input_indices_retained": True,
        "outward_interval_addition_only": True,
        "physical_scale_reapplied": False,
        "physical_measure_reapplied": False,
        "radial_weight_reapplied": False,
        "pressure_added": False,
        "P1_added": False,
        "outer_Q_added": False,
        "tail_added": False,
        "ledger_appended": False,
        "cross_cell_grouping": False,
        "cross_branch_grouping": False,
        "error_only_as_signed": False,
        "error_only_records_separate": True,
        "signed_marker_missing_condition": obligation_record["condition"],
        "upstream_source_authentication_is_caller_obligation": True,
        "same_interval_context_is_caller_obligation": True,
        "same_interval_context_validated": False,
        "same_source_claimed": False,
        "same_gauge_claimed": False,
        "physical_field_claimed": False,
        "physical_completeness_claimed": False,
        "independent_error_certification": False,
        "finite_P_branch_context_only": True,
        "no_signed_field_label_added": True,
    }


sum_compensation_cells = build
cell_sum = build
aggregate = build
combine = build

