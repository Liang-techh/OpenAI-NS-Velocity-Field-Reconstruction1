"""Bind the saved first ``Kchi`` core cells without evaluating them.

The one-entry ``run_n1_core_kchi_first_v1_once.py`` output already contains
the finite ``j=1`` source-to-``W0`` correction.  This module is deliberately
only a structural adapter for that output.  It copies the saved ``F1`` and
``d_x_F1`` terms, their source/target cell ownership, source identifiers,
ordinary ``Z0/Z1/Z2`` rows, four formal source powers, and exact binary
endpoint records.  It does not import the producer, select a centre, perform
interval arithmetic, sum terms, or evaluate a field.

The ``j>=2`` packet is kept in a different branch.  Its coefficients are
required to be symmetric error intervals and are marked ``error_only``;
they are never exposed as signed finite fields and are never merged into the
finite correction.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from copy import deepcopy
from typing import Any


SCHEMA = "original_kchi_first_saved_cells_v2"
CELL_COUNT = 8
ORDERS = ("Z0", "Z1", "Z2")
ORDER_NUMBERS = (0, 1, 2)
FIELDS = ("F1", "d_x_F1")
SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
FAMILY_FIELDS = (
    "actual_five_defect_family_sha256",
    "implicit_source_sha256",
    "datum_enclosure_sha256",
)
TAIL_ERROR_CLASS = "Kchi_j_ge_2_kernel_tail"

# Ordinary derivatives, not Taylor coefficients.  This is the exact product
# rule used by the saved producer for ``chi * source``.
LEIBNIZ_RULE = {
    0: {(0, 0): 1},
    1: {(1, 0): 1, (0, 1): 1},
    2: {(2, 0): 1, (1, 1): 2, (0, 2): 1},
}
LEIBNIZ_RULE_JSON = {
    "Z0": [[0, 0, 1]],
    "Z1": [[1, 0, 1], [0, 1, 1]],
    "Z2": [[2, 0, 1], [1, 1, 2], [0, 2, 1]],
}


class SavedKchiInputError(ValueError):
    """Raised when a saved Kchi packet is structurally incomplete."""


def _copy(value: Any) -> Any:
    """Copy caller data without asking it for a numerical representation."""

    return deepcopy(value)


def _is_int(value: Any) -> bool:
    return type(value) is int


def _mapping(value: Any, label: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise SavedKchiInputError(label + " must be a mapping")
    return value


def _sequence(value: Any, label: str, *, nonempty: bool = False) -> Sequence[Any]:
    if isinstance(value, (str, bytes, bytearray)) or not isinstance(value, Sequence):
        raise SavedKchiInputError(label + " must be a sequence")
    if nonempty and not value:
        raise SavedKchiInputError(label + " must not be empty")
    return value


def _binary_tuple(value: Any, label: str) -> tuple[int, int, int, int]:
    if (
        not isinstance(value, list)
        or len(value) != 4
        or any(type(item) is not int for item in value)
    ):
        raise SavedKchiInputError(label + " must be a four-integer binary tuple")
    sign, mantissa, exponent, count = value
    if sign not in (0, 1) or mantissa < 0 or count < 0:
        raise SavedKchiInputError(label + " has invalid sign, mantissa, or count")
    if mantissa == 0:
        if sign != 0 or exponent != 0 or count != 0:
            raise SavedKchiInputError(label + " has an invalid zero tuple")
    elif count != mantissa.bit_length():
        raise SavedKchiInputError(label + " has an invalid mantissa count")
    return sign, mantissa, exponent, count


def _compare_binary_tuples(left: tuple[int, int, int, int], right: tuple[int, int, int, int]) -> int:
    """Compare exact dyadic tuples without converting them to a numeric value."""

    lsign, lmantissa, _lexponent, _lcount = left
    rsign, rmantissa, _rexponent, _rcount = right
    if lmantissa == 0 or rmantissa == 0:
        if lmantissa == rmantissa == 0:
            return 0
        if lmantissa == 0:
            return 1 if rsign == 1 else -1
        return -1 if lsign == 1 else 1
    if lsign != rsign:
        return -1 if lsign == 1 else 1
    # A normalized positive magnitude has binary top bit exponent+count.
    ltop = left[2] + left[3]
    rtop = right[2] + right[3]
    if ltop != rtop:
        magnitude = -1 if ltop < rtop else 1
    else:
        lscaled = lmantissa << right[3]
        rscaled = rmantissa << left[3]
        magnitude = -1 if lscaled < rscaled else (1 if lscaled > rscaled else 0)
    return -magnitude if lsign else magnitude


def _endpoint_record(value: Any, label: str) -> Mapping[str, Any]:
    """Validate an exact binary interval record without decoding its value."""

    record = _mapping(value, label)
    _binary_tuple(record.get("lower_exact_mpf_tuple"), label + ".lower")
    _binary_tuple(record.get("upper_exact_mpf_tuple"), label + ".upper")
    if record.get("directed") is not True:
        raise SavedKchiInputError(label + " must retain directed=True")
    _require_endpoint_order(record, label)
    return record


def _require_endpoint_order(value: Mapping[str, Any], label: str) -> None:
    lower = _binary_tuple(value["lower_exact_mpf_tuple"], label + ".lower")
    upper = _binary_tuple(value["upper_exact_mpf_tuple"], label + ".upper")
    if _compare_binary_tuples(lower, upper) > 0:
        raise SavedKchiInputError(label + " has reversed binary endpoints")


def _powers(value: Any, label: str) -> list[int]:
    raw = _sequence(value, label)
    if len(raw) != 4 or any(not _is_int(item) for item in raw):
        raise SavedKchiInputError(label + " must contain four integer source powers")
    return [int(item) for item in raw]


def _ids(value: Any, label: str) -> list[Any]:
    raw = _sequence(value, label, nonempty=True)
    result = [_copy(item) for item in raw]
    if any(item is None or item == "" for item in result):
        raise SavedKchiInputError(label + " must contain nonempty term IDs")
    return result


def _require_bool(row: Mapping[str, Any], key: str, value: bool, label: str) -> None:
    if row.get(key) is not value:
        raise SavedKchiInputError(label + " must have " + key + "=" + str(value))


def _require_scope_flags(row: Mapping[str, Any], label: str) -> None:
    for key in (
        "Kchi_higher_words_included",
        "allword_tail_added",
        "source_terms_recomputed",
        "old_words_recomputed",
        "appended_to_cumulative_m15",
        "complete_W0_computed",
        "complete_m15_computed",
        "terminal_closure_proved",
    ):
        if key in row:
            _require_bool(row, key, False, label)


def _order_value(row: Mapping[str, Any], order: str, label: str) -> Any:
    """Read the runner's string order keys, accepting integer aliases only."""

    number = int(order.removeprefix("Z"))
    candidates = (order, str(number), number)
    found = [candidate for candidate in candidates if candidate in row]
    if not found:
        raise SavedKchiInputError(label + " is missing " + order)
    if len(found) > 1:
        first = row[found[0]]
        if any(row[candidate] != first for candidate in found[1:]):
            raise SavedKchiInputError(label + " has duplicate order keys for " + order)
    return row[found[0]]


def _validate_source_family(value: Any, label: str) -> dict[str, Any]:
    family = _mapping(value, label)
    for key in FAMILY_FIELDS:
        digest = family.get(key)
        if not isinstance(digest, str) or len(digest) != 64:
            raise SavedKchiInputError(label + " is missing " + key)
    return _copy(dict(family))


def _validate_source_logs(value: Any, label: str) -> dict[str, Any]:
    """Require four encoded same-source logs, allowing legacy no-directed records."""

    logs = _mapping(value, label)
    if set(logs) != set(SOURCE_POWER_BASIS):
        raise SavedKchiInputError(label + " must contain exactly the four source logs")
    for name in SOURCE_POWER_BASIS:
        record = _mapping(logs.get(name), label + "." + name)
        _binary_tuple(record.get("lower_exact_mpf_tuple"), label + "." + name + ".lower")
        _binary_tuple(record.get("upper_exact_mpf_tuple"), label + "." + name + ".upper")
        if "directed" in record and record.get("directed") is not True:
            raise SavedKchiInputError(label + "." + name + " must retain directed=True when present")
        _require_endpoint_order(record, label + "." + name)
    return _copy(dict(logs))


def _reject_finite_error_markers(term: Mapping[str, Any], label: str) -> None:
    """A saved error branch may never be relabelled as a signed finite term."""

    if term.get("error_only") is True or term.get("symmetric_error") is True:
        raise SavedKchiInputError(label + " is an error-only term, not finite j=1")
    if term.get("signed_field") is False or term.get("finite_branch") is False:
        raise SavedKchiInputError(label + " explicitly rejects signed finite use")


def _validate_origin(term: Mapping[str, Any], label: str) -> Mapping[str, Any]:
    origin = _mapping(term.get("origin"), label + ".origin")
    _ids(origin.get("source_ids"), label + ".origin.source_ids")
    classes = _sequence(origin.get("error_classes"), label + ".origin.error_classes", nonempty=True)
    if any(not isinstance(item, str) or not item for item in classes):
        raise SavedKchiInputError(label + ".origin.error_classes must be nonempty strings")
    if type(origin.get("uncertainty")) is not bool:
        raise SavedKchiInputError(label + ".origin.uncertainty must be boolean")
    if origin.get("kernel_order") != 1:
        raise SavedKchiInputError(label + ".origin is not the finite j=1 term")
    return origin


def _validate_finite_term(
    term: Any,
    *,
    field: str,
    order: str,
    source_index: int,
    target_index: int,
    position: int,
) -> dict[str, Any]:
    label = (
        "first_correction.target_cells[{}].source[{}].{}.{}[{}]".format(
            target_index, source_index, field, order, position
        )
    )
    raw = _mapping(term, label)
    _reject_finite_error_markers(raw, label)
    _endpoint_record(raw.get("coefficient"), label + ".coefficient")
    powers = _powers(raw.get("powers"), label + ".powers")

    if raw.get("source_cell") != source_index or raw.get("source_index") != source_index:
        raise SavedKchiInputError(label + " source-cell index differs from its container")
    if raw.get("target_cell") != target_index:
        raise SavedKchiInputError(label + " target-cell index differs from its container")
    for key in ("source_cell_id", "source_term_id"):
        value = raw.get(key)
        if value is None or value == "":
            raise SavedKchiInputError(label + " is missing " + key)
    if not _is_int(raw.get("source_term_index")) or raw["source_term_index"] < 0:
        raise SavedKchiInputError(label + " has no source-term index")
    _ids(raw.get("source_ids"), label + ".source_ids")
    _validate_origin(raw, label)

    if raw.get("principal_kernel") != "(-2*chi) T K0":
        raise SavedKchiInputError(label + " has the wrong principal kernel")
    for key in ("chi_order", "source_order", "leibniz_multiplier"):
        if not _is_int(raw.get(key)) or raw[key] < 0:
            raise SavedKchiInputError(label + " is missing " + key)
    z_order = int(order.removeprefix("Z"))
    pair = (raw["chi_order"], raw["source_order"])
    if pair not in LEIBNIZ_RULE[z_order] or raw["chi_order"] + raw["source_order"] != z_order:
        raise SavedKchiInputError(label + " has a Leibniz pair from the wrong Z row")
    if raw["leibniz_multiplier"] != LEIBNIZ_RULE[z_order][pair]:
        raise SavedKchiInputError(label + " has the wrong ordinary Leibniz multiplier")
    if raw.get("kernel_mass_field") != field:
        raise SavedKchiInputError(label + " kernel mass field differs from its saved field")
    for key in ("target_x", "source_x", "kernel_mass", "kernel_mass_bare", "correction_factor"):
        _endpoint_record(raw.get(key), label + "." + key)

    # Keep the producer's record intact and add only nonnumeric binder aliases.
    selected = _copy(dict(raw))
    selected["formal_source_powers"] = _copy(powers)
    selected["original_term"] = _copy(dict(raw))
    selected["term_role"] = "finite_j1_signed"
    selected["error_only"] = False
    selected["symmetric_error"] = False
    return selected


def _validate_source_rows(
    rows: Any,
    *,
    source_index: int,
    target_index: int,
) -> dict[str, dict[str, list[dict[str, Any]]]]:
    label = "first_correction.target_cells[{}].source[{}].correction".format(
        target_index, source_index
    )
    mapping = _mapping(rows, label)
    output: dict[str, dict[str, list[dict[str, Any]]]] = {}
    # The saved producer stores ``correction`` as field-first
    # ``{F1: {Z0: ...}, d_x_F1: {Z0: ...}}``.  An order-first packet is not an
    # equivalent source record and is rejected rather than silently adapted.
    if any(order in mapping for order in ORDERS):
        raise SavedKchiInputError(label + " must retain field-first producer rows")
    if any(field not in mapping for field in FIELDS):
        raise SavedKchiInputError(label + " is missing a field-first producer row")
    for order in ORDERS:
        output[order] = {}
        for field in FIELDS:
            field_rows = _mapping(mapping.get(field), label + "." + field)
            terms = _sequence(
                _order_value(field_rows, order, label + "." + field),
                label + "." + field + "." + order,
                nonempty=True,
            )
            output[order][field] = [
                _validate_finite_term(
                    term,
                    field=field,
                    order=order,
                    source_index=source_index,
                    target_index=target_index,
                    position=position,
                )
                for position, term in enumerate(terms)
            ]
    return output


def _validate_target_cell(raw: Any, index: int) -> dict[str, Any]:
    label = "first_correction.target_cells[{}]".format(index)
    cell = _mapping(raw, label)
    if cell.get("index") != index or cell.get("target_index") != index:
        raise SavedKchiInputError(label + " index differs from its position")
    expected_sources = list(range(index + 1))
    if cell.get("source_indices") != expected_sources:
        raise SavedKchiInputError(label + " source_indices differ from the saved triangular cells")
    if cell.get("triangular_diagonal_source") != index:
        raise SavedKchiInputError(label + " diagonal source differs from target cell")
    if cell.get("source_j_le_target_i") is not True:
        raise SavedKchiInputError(label + " must retain source_j_le_target_i=True")
    _endpoint_record(cell.get("x"), label + ".x")
    _endpoint_record(cell.get("target_x"), label + ".target_x")

    by_source = _mapping(cell.get("W0_correction_by_source_cell"), label + ".W0_correction_by_source_cell")
    selected_sources: list[dict[str, Any]] = []
    for source_index in expected_sources:
        key = str(source_index)
        if key not in by_source and source_index not in by_source:
            raise SavedKchiInputError(label + " is missing source cell " + key)
        source = by_source[key] if key in by_source else by_source[source_index]
        source_label = label + ".source[{}]".format(source_index)
        source_mapping = _mapping(source, source_label)
        if source_mapping.get("source_index") != source_index or source_mapping.get("target_index") != index:
            raise SavedKchiInputError(source_label + " source/target indices differ")
        if source_mapping.get("triangular") is not (source_index == index):
            raise SavedKchiInputError(source_label + " triangular marker differs")
        _endpoint_record(source_mapping.get("source_x"), source_label + ".source_x")
        _endpoint_record(source_mapping.get("target_x"), source_label + ".target_x")
        ordinary = _validate_source_rows(
            source_mapping.get("correction"),
            source_index=source_index,
            target_index=index,
        )
        selected_sources.append(
            {
                "source_index": source_index,
                "target_index": index,
                "triangular": source_index == index,
                "source_x": _copy(source_mapping["source_x"]),
                "target_x": _copy(source_mapping["target_x"]),
                "ordinary_Z": ordinary,
                "original_source_cell": _copy(dict(source_mapping)),
            }
        )

    return {
        "index": index,
        "target_cell": index,
        "target_index": index,
        "source_indices": _copy(expected_sources),
        "triangular_diagonal_source": index,
        "axis_endpoint_continuation": _copy(cell.get("axis_endpoint_continuation")),
        "x": _copy(cell["x"]),
        "target_x": _copy(cell["target_x"]),
        "source_cells": selected_sources,
        "source_cells_reused_without_aggregation": True,
        "original_target_cell": _copy(dict(cell)),
    }


def _symmetric_endpoint_record(value: Any, label: str) -> Mapping[str, Any]:
    record = _endpoint_record(value, label)
    _require_endpoint_order(record, label)
    lower = record["lower_exact_mpf_tuple"]
    upper = record["upper_exact_mpf_tuple"]
    zero = lower[1] == 0 and upper[1] == 0
    if zero and (lower != [0, 0, 0, 0] or upper != [0, 0, 0, 0]):
        raise SavedKchiInputError(label + " has a malformed zero symmetric interval")
    if not zero and (lower[0], upper[0]) != (1, 0):
        raise SavedKchiInputError(label + " is not a symmetric error interval")
    if not zero and lower[1:] != upper[1:]:
        raise SavedKchiInputError(label + " does not retain symmetric binary endpoints")
    return record


def _validate_tail_term(term: Any, *, field: str, order: str, position: int) -> dict[str, Any]:
    label = "W0_separate_symmetric_kernel_tail.{}.{}[{}]".format(field, order, position)
    raw = _mapping(term, label)
    _symmetric_endpoint_record(raw.get("coefficient"), label + ".coefficient")
    _powers(raw.get("powers"), label + ".powers")
    _ids(raw.get("source_ids"), label + ".source_ids")
    classes = _sequence(raw.get("error_classes"), label + ".error_classes", nonempty=True)
    if TAIL_ERROR_CLASS not in classes:
        raise SavedKchiInputError(label + " is missing the j>=2 tail error class")
    if raw.get("uncertainty") is not True:
        raise SavedKchiInputError(label + " must retain uncertainty=True")
    source_order = raw.get("source_order")
    kernel_z_order = raw.get("kernel_Z_order")
    if not _is_int(source_order) or not _is_int(kernel_z_order):
        raise SavedKchiInputError(label + " is missing tail order provenance")
    if source_order < 0 or kernel_z_order < 0:
        raise SavedKchiInputError(label + " has negative tail order provenance")
    ordinary_order = int(order.removeprefix("Z"))
    if source_order + kernel_z_order != ordinary_order:
        raise SavedKchiInputError(label + " tail orders do not sum to ordinary " + order)
    if raw.get("physical_source_width_bound") != 2:
        raise SavedKchiInputError(label + " is missing the saved source-width bound")
    if "error_only" in raw and raw.get("error_only") is not True:
        raise SavedKchiInputError(label + " explicitly claims signed use")
    if "symmetric_error" in raw and raw.get("symmetric_error") is not True:
        raise SavedKchiInputError(label + " explicitly claims nonsymmetric use")
    if "signed_field" in raw and raw.get("signed_field") is not False:
        raise SavedKchiInputError(label + " explicitly claims a signed field")
    selected = _copy(dict(raw))
    selected["ordinary_Z"] = order
    selected["field"] = field
    selected["term_role"] = "j_ge_2_symmetric_error"
    selected["error_only"] = True
    selected["symmetric_error"] = True
    selected["signed_field"] = False
    return selected


def _validate_tail_bounds(value: Any) -> dict[str, Any]:
    bounds = _mapping(value, "kernel_tail_bounds")
    if bounds.get("ordinary_Z_orders") != [0, 1, 2]:
        raise SavedKchiInputError("kernel_tail_bounds ordinary orders are incomplete")
    if bounds.get("first_omitted_index") != 2:
        raise SavedKchiInputError("kernel_tail_bounds must begin at j>=2")
    if bounds.get("bounds_are_not_field_values") is not True:
        raise SavedKchiInputError("kernel_tail_bounds must remain error-only")
    if bounds.get("source_or_feedback_included") is not False:
        raise SavedKchiInputError("kernel_tail_bounds source scope differs")
    for field in FIELDS:
        rows = _sequence(bounds.get(field), "kernel_tail_bounds." + field)
        if len(rows) != 3:
            raise SavedKchiInputError("kernel_tail_bounds." + field + " must retain Z0/Z1/Z2")
        for order, interval in enumerate(rows):
            _endpoint_record(interval, "kernel_tail_bounds.{}.Z{}".format(field, order))
    return _copy(dict(bounds))


def _validate_top_scope(saved: Mapping[str, Any]) -> None:
    for key, expected in (
        ("first_correction_signed", True),
        ("Kchi_tail_is_uncertainty_not_signed_field", True),
        ("tail_includes_feedback", False),
        ("feedback_cross_terms_added", False),
        ("appended_to_cumulative_m15", False),
        ("complete_W0_computed", False),
        ("complete_m15_computed", False),
        ("terminal_closure_proved", False),
        ("source_cells_reuse_one_whole_core_enclosure", True),
        ("source_boxes_not_point_values", True),
    ):
        _require_bool(saved, key, expected, "saved output")


def bind(saved_output: Any) -> dict[str, Any]:
    """Copy the saved finite ``j=1`` cells and separate ``j>=2`` error.

    ``saved_output`` is the decoded ``result.output`` object returned by the
    one-entry runner.  The function is intentionally source-bound: it accepts
    only already saved binary records and does not call any producer or
    numerical evaluator.
    """

    saved = _mapping(saved_output, "saved_output")
    _validate_top_scope(saved)
    family = _validate_source_family(saved.get("source_family"), "saved_output.source_family")
    logs = _validate_source_logs(saved.get("source_logs"), "saved_output.source_logs")
    domain = saved.get("domain")
    if not isinstance(domain, str) or not domain:
        raise SavedKchiInputError("saved_output.domain must retain the producer domain")
    input_index = _mapping(saved.get("input_index"), "saved_output.input_index")

    finite = _mapping(saved.get("first_correction"), "saved_output.first_correction")
    if finite.get("source_family") != family:
        raise SavedKchiInputError("first_correction source_family differs")
    if finite.get("source_atoms") != logs:
        raise SavedKchiInputError("first_correction source_atoms differ from source_logs")
    if finite.get("source_power_basis") != list(SOURCE_POWER_BASIS):
        raise SavedKchiInputError("first_correction four-source basis is missing")
    if finite.get("ordinary_Z_orders") != [0, 1, 2]:
        raise SavedKchiInputError("first_correction must retain ordinary Z0/Z1/Z2")
    if finite.get("ordinary_product_rule") != LEIBNIZ_RULE_JSON:
        raise SavedKchiInputError("first_correction ordinary product rule differs")
    for order in ORDERS:
        if finite.get("ordinary_{}_product_rule".format(order)) != LEIBNIZ_RULE_JSON[order]:
            raise SavedKchiInputError("first_correction {} product rule differs".format(order))
    if finite.get("cell_count") != CELL_COUNT:
        raise SavedKchiInputError("first_correction must retain eight cells")
    if finite.get("signed_first_correction") != "(-2*chi) T K0":
        raise SavedKchiInputError("saved finite correction is not the j=1 Kchi term")
    for key, expected in (
        ("Kchi_first_correction_only", True),
        ("Kchi_higher_words_included", False),
        ("allword_tail_added", False),
        ("source_to_W0_only", True),
        ("source_rows_preserved_through_Z2", True),
        ("source_ids_retained", True),
        ("source_indices_retained", True),
        ("signed_coefficients_retained", True),
        ("formal_powers_retained", True),
        ("source_terms_recomputed", False),
        ("old_words_recomputed", False),
    ):
        _require_bool(finite, key, expected, "first_correction")
    _require_scope_flags(finite, "first_correction")

    raw_cells = _sequence(finite.get("target_cells"), "first_correction.target_cells")
    if len(raw_cells) != CELL_COUNT:
        raise SavedKchiInputError("first_correction.target_cells must contain eight cells")
    cells = [_validate_target_cell(raw, index) for index, raw in enumerate(raw_cells)]

    raw_tail = _mapping(
        saved.get("W0_separate_symmetric_kernel_tail"),
        "W0_separate_symmetric_kernel_tail",
    )
    tail: dict[str, dict[str, list[dict[str, Any]]]] = {}
    tail_counts: dict[str, dict[str, int]] = {}
    for field in FIELDS:
        field_rows = _mapping(raw_tail.get(field), "W0_separate_symmetric_kernel_tail." + field)
        tail[field] = {}
        tail_counts[field] = {}
        for order in ORDERS:
            number = str(int(order.removeprefix("Z")))
            terms = _order_value(field_rows, order, "W0_separate_symmetric_kernel_tail." + field)
            terms = _sequence(terms, "W0_separate_symmetric_kernel_tail.{}.{}".format(field, number), nonempty=True)
            selected = [
                _validate_tail_term(term, field=field, order=order, position=position)
                for position, term in enumerate(terms)
            ]
            tail[field][order] = selected
            tail_counts[field][order] = len(selected)

    bounds = _validate_tail_bounds(saved.get("kernel_tail_bounds"))
    return {
        "schema": SCHEMA,
        "status": "saved_kchi_first_cells",
        "source_family": family,
        "source_logs": _copy(dict(logs)),
        "domain": _copy(domain),
        "input_index": _copy(dict(input_index)),
        "formal_source_basis": list(SOURCE_POWER_BASIS),
        "ordinary_Z_orders": list(ORDERS),
        "fields": list(FIELDS),
        "cells": cells,
        "finite_j1_correction": True,
        "finite_correction_is_signed": True,
        "finite_correction_kernel": "(-2*chi) T K0",
        "symmetric_kernel_tail": tail,
        "symmetric_kernel_tail_counts": tail_counts,
        "kernel_tail_bounds": bounds,
        "tail_is_error_only": True,
        "tail_is_not_a_signed_field": True,
        "terms_were_not_summed": True,
        "ledger_was_not_extended": True,
        "point_value_selected": False,
        "source_values_computed": False,
        "numeric_operator_called": False,
        "complete_W0_computed": False,
        "complete_response": False,
        "physical_m15_computed": False,
        "appended_to_cumulative_m15": False,
        "endpoint_closure_proved": False,
        "terminal_closure_proved": False,
        "source_identity_verified": False,
        "source_identity_status": "source identity unverified",
        "coverage_recertified": False,
    }


__all__ = [
    "SCHEMA",
    "CELL_COUNT",
    "ORDERS",
    "FIELDS",
    "SOURCE_POWER_BASIS",
    "LEIBNIZ_RULE",
    "LEIBNIZ_RULE_JSON",
    "SavedKchiInputError",
    "bind",
]

