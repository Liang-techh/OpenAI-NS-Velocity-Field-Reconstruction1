"""Conditional grouping consumer for the saved K0 symbolic m15 terms.

This leaf consumes the already materialized output shape of
``original_k0_symbolic_m15_cells_v1``.  It does not import or call that
producer, and it does not perform any background multiplication, integration,
source lookup, or ledger operation.  Its only numeric operation is the
caller-supplied bounded ``sum_intervals`` function applied to each group's
raw binary coefficient intervals.

The grouping key is deliberately exact and includes the cell index, channel,
ordinary order, symbolic piece, signed factor, six formal source powers,
background name and record, rho record, product-factor sequence, and local
factor record.  Source IDs are retained as provenance; they are never used as
a same-source proof.  Every input term is copied into exactly one member
record with its original term, coefficient, source IDs, cell/order location,
and term index.
"""

from __future__ import annotations

import operator
from collections.abc import Mapping, Sequence
from copy import deepcopy
from typing import Any, Callable


SCHEMA = "original_k0_grouped_dyadic_v1"
SYMBOLIC_SCHEMA = "original_k0_symbolic_m15_cells_v1"
DEFAULT_PRECISION_BITS = 2048
MAX_PRECISION_BITS = 4096
MAX_EXPONENT_BITS = 256
MAX_ENDPOINT_MANTISSA_BITS = 65_536
MAX_TERMS = 100_000
MAX_GROUPS = 100_000

ORDERS = ("Z0", "Z1")
CHANNELS = ("axial", "swirl")
EXPECTED_TERM_COUNTS = {
    "axial": {"Z0": 540, "Z1": 1620},
    "swirl": {"Z0": 396, "Z1": 1584},
}
EXPECTED_TOTAL_TERMS = 4140
REQUIRED_UNKNOWN_NAMES = ("Kchi", "N1p", "allword_tail")

GROUP_KEY_FIELDS = (
    "cell_index",
    "channel",
    "ordinary_order",
    "piece",
    "signed_factor",
    "formal_source_powers",
    "background_factor_name",
    "background_factor",
    "rho_factor",
    "product_factors",
    "local_factor",
)


class GroupedDyadicInputError(ValueError):
    """Raised when the symbolic receipt is not safe to consume."""


class GroupedDyadicResourceError(GroupedDyadicInputError):
    """Raised when explicit grouping or precision limits would be exceeded."""


def _index(value: Any, label: str) -> int:
    if isinstance(value, bool):
        raise GroupedDyadicInputError(f"{label} must be an integer, not bool")
    try:
        result = operator.index(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise GroupedDyadicInputError(f"{label} must support operator.index") from exc
    if isinstance(result, bool) or type(result) is not int:
        raise GroupedDyadicInputError(f"{label} must produce a plain integer")
    return result


def _precision(value: Any) -> int:
    result = _index(value, "precision")
    if result <= 0:
        raise GroupedDyadicInputError("precision must be positive")
    if result > MAX_PRECISION_BITS:
        raise GroupedDyadicResourceError("precision exceeds MAX_PRECISION_BITS")
    return result


def _copy(value: Any) -> Any:
    return deepcopy(value)


def _freeze(value: Any, label: str = "value") -> tuple[Any, ...]:
    """Validate JSON-shaped data and make an order-independent key form."""

    if value is None:
        return ("none",)
    if isinstance(value, bool):
        return ("bool", value)
    if type(value) is int:
        return ("int", value)
    if isinstance(value, str):
        return ("str", value)
    if isinstance(value, Mapping):
        pairs: list[tuple[str, tuple[Any, ...]]] = []
        for key, item in value.items():
            if not isinstance(key, str):
                raise GroupedDyadicInputError(f"{label} mapping keys must be strings")
            pairs.append((key, _freeze(item, f"{label}.{key}")))
        pairs.sort(key=lambda pair: pair[0])
        return ("mapping", tuple(pairs))
    if isinstance(value, tuple):
        return ("tuple", tuple(_freeze(item, f"{label}[]") for item in value))
    if isinstance(value, list):
        return ("list", tuple(_freeze(item, f"{label}[]") for item in value))
    raise GroupedDyadicInputError(
        f"{label} contains a non-JSON value of type {type(value).__name__}"
    )


def _safe_copy(value: Any, label: str) -> Any:
    _freeze(value, label)
    return _copy(value)


def _safe_exponent(value: int, label: str) -> int:
    if abs(value).bit_length() > MAX_EXPONENT_BITS:
        raise GroupedDyadicResourceError(
            f"{label} exceeds the {MAX_EXPONENT_BITS}-bit exponent limit"
        )
    return value


def _endpoint_tuple(value: Any, label: str, *, allow_list: bool) -> tuple[int, int, int, int]:
    if allow_list:
        valid_container = isinstance(value, (tuple, list))
    else:
        valid_container = isinstance(value, tuple)
    if not valid_container or len(value) != 4:
        raise GroupedDyadicInputError(
            f"{label} must contain exactly four integer fields"
        )
    fields = tuple(_index(item, f"{label}[{index}]") for index, item in enumerate(value))
    sign, mantissa, exponent, bitcount = fields
    _safe_exponent(exponent, f"{label}.exponent")
    if sign not in (0, 1):
        raise GroupedDyadicInputError(f"{label}.sign must be 0 or 1")
    if mantissa < 0:
        raise GroupedDyadicInputError(f"{label}.mantissa must be nonnegative")
    if bitcount < 0:
        raise GroupedDyadicInputError(f"{label}.bitcount is a special sentinel")
    if mantissa == 0:
        if exponent != 0 or bitcount != 0:
            raise GroupedDyadicInputError(
                f"{label} is a nonstandard zero or special endpoint"
            )
        return (0, 0, 0, 0)
    if bitcount == 0 or bitcount != mantissa.bit_length():
        raise GroupedDyadicInputError(
            f"{label}.bitcount must equal mantissa.bit_length()"
        )
    if bitcount > MAX_ENDPOINT_MANTISSA_BITS:
        raise GroupedDyadicResourceError(
            f"{label}.mantissa exceeds MAX_ENDPOINT_MANTISSA_BITS"
        )
    return fields


def _compare_abs(left: tuple[int, int, int, int], right: tuple[int, int, int, int]) -> int:
    if left[1] == 0 or right[1] == 0:
        if left[1] == 0 and right[1] == 0:
            return 0
        return -1 if left[1] == 0 else 1
    left_top = left[2] + left[3]
    right_top = right[2] + right[3]
    if left_top != right_top:
        return 1 if left_top > right_top else -1
    shift = left[2] - right[2]
    if abs(shift) > max(left[3], right[3]):
        raise GroupedDyadicResourceError("comparison shift exceeded mantissa bounds")
    if shift >= 0:
        left_mantissa = left[1] << shift
        right_mantissa = right[1]
    else:
        left_mantissa = left[1]
        right_mantissa = right[1] << (-shift)
    return (left_mantissa > right_mantissa) - (left_mantissa < right_mantissa)


def _compare_signed(left: tuple[int, int, int, int], right: tuple[int, int, int, int]) -> int:
    left_zero = left[1] == 0
    right_zero = right[1] == 0
    if left_zero or right_zero:
        if left_zero and right_zero:
            return 0
        if left_zero:
            return 1 if right[0] else -1
        return -1 if right[0] else 1
    if left[0] != right[0]:
        return -1 if left[0] else 1
    magnitude = _compare_abs(left, right)
    return magnitude if left[0] == 0 else -magnitude


def _coefficient_interval(value: Any, label: str) -> tuple[tuple[int, int, int, int], tuple[int, int, int, int]]:
    if not isinstance(value, Mapping):
        raise GroupedDyadicInputError(f"{label} must be a raw coefficient mapping")
    for key in ("lower_exact_mpf_tuple", "upper_exact_mpf_tuple"):
        if key not in value:
            raise GroupedDyadicInputError(f"{label} is missing {key}")
    if value.get("directed") is not True:
        raise GroupedDyadicInputError(f"{label}.directed must be true")
    _safe_copy(value, label)
    lower = _endpoint_tuple(
        value["lower_exact_mpf_tuple"],
        label + ".lower_exact_mpf_tuple",
        allow_list=True,
    )
    upper = _endpoint_tuple(
        value["upper_exact_mpf_tuple"],
        label + ".upper_exact_mpf_tuple",
        allow_list=True,
    )
    if _compare_signed(lower, upper) > 0:
        raise GroupedDyadicInputError(f"{label} endpoints are reversed")
    return (lower, upper)


def _required_sequence(value: Any, label: str, *, permit_empty: bool) -> Sequence[Any]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes, bytearray)):
        raise GroupedDyadicInputError(f"{label} must be a sequence")
    if not permit_empty and len(value) == 0:
        raise GroupedDyadicInputError(f"{label} must not be empty")
    return value


def _validate_unknowns(saved_symbolic: Mapping[str, Any]) -> list[Any]:
    unknowns = saved_symbolic.get("unknown_terms")
    if not isinstance(unknowns, Sequence) or isinstance(unknowns, (str, bytes, bytearray)):
        raise GroupedDyadicInputError("unknown_terms must be a sequence")
    found: set[str] = set()
    for index, item in enumerate(unknowns):
        if not isinstance(item, Mapping):
            raise GroupedDyadicInputError(f"unknown_terms[{index}] must be a mapping")
        name = item.get("name")
        if not isinstance(name, str) or item.get("status") != "unknown":
            raise GroupedDyadicInputError("unknown terms must retain name/status unknown")
        found.add(name)
        _safe_copy(item, f"unknown_terms[{index}]")
    missing = [name for name in REQUIRED_UNKNOWN_NAMES if name not in found]
    if missing:
        raise GroupedDyadicInputError("required unknown terms are missing: " + ",".join(missing))
    if saved_symbolic.get("unknown_terms_are_not_zero") is not True:
        raise GroupedDyadicInputError("unknown_terms_are_not_zero must remain true")
    return _safe_copy(list(unknowns), "unknown_terms")


def _validate_symbolic(saved_symbolic: Any) -> tuple[Mapping[str, Any], list[Mapping[str, Any]], dict[str, dict[str, int]], int, bool]:
    if not isinstance(saved_symbolic, Mapping):
        raise GroupedDyadicInputError("saved_symbolic mapping required")
    if saved_symbolic.get("schema") != SYMBOLIC_SCHEMA:
        raise GroupedDyadicInputError("saved_symbolic schema is not K0 symbolic m15")
    _safe_copy(saved_symbolic.get("source_family"), "source_family")
    _safe_copy(saved_symbolic.get("source_logs"), "source_logs")
    scope = saved_symbolic.get("source_kernel_scope")
    if not isinstance(scope, str) or not scope:
        raise GroupedDyadicInputError("source_kernel_scope is required")
    if saved_symbolic.get("source_identity_verified") is not False:
        raise GroupedDyadicInputError("source identity must remain unverified")
    f0_identity = saved_symbolic.get("F0_scale_identity")
    if not isinstance(f0_identity, Mapping) or f0_identity.get("identified") is not False:
        raise GroupedDyadicInputError("F0_anchor/F0_local must remain independent")
    _validate_unknowns(saved_symbolic)
    for key in ("N1p_included", "Kchi_included", "allword_tail_included", "appended_to_212_ledger"):
        if saved_symbolic.get(key) is not False:
            raise GroupedDyadicInputError(f"{key} must remain false")

    raw_cells = saved_symbolic.get("cells")
    if not isinstance(raw_cells, Sequence) or isinstance(raw_cells, (str, bytes, bytearray)):
        raise GroupedDyadicInputError("saved_symbolic.cells must be a sequence")
    if len(raw_cells) != 8:
        raise GroupedDyadicInputError("saved_symbolic must contain exactly eight cells")

    synthetic = saved_symbolic.get("synthetic_only") is True
    counts = {channel: {order: 0 for order in ORDERS} for channel in CHANNELS}
    cells: list[Mapping[str, Any]] = []
    total = 0
    for index, cell in enumerate(raw_cells):
        if not isinstance(cell, Mapping) or cell.get("index") != index:
            raise GroupedDyadicInputError("cell order/index is not exact")
        terms = cell.get("terms")
        if not isinstance(terms, Mapping):
            raise GroupedDyadicInputError(f"cells[{index}].terms must be a mapping")
        for channel in CHANNELS:
            rows = terms.get(channel)
            if not isinstance(rows, Mapping):
                raise GroupedDyadicInputError(f"cells[{index}].terms.{channel} is missing")
            for order in ORDERS:
                sequence = _required_sequence(
                    rows.get(order),
                    f"cells[{index}].terms.{channel}.{order}",
                    permit_empty=synthetic,
                )
                counts[channel][order] += len(sequence)
                total += len(sequence)
        _safe_copy(cell, f"cells[{index}]")
        cells.append(cell)

    declared_counts = saved_symbolic.get("term_counts")
    if not isinstance(declared_counts, Mapping):
        raise GroupedDyadicInputError("saved_symbolic.term_counts is required")
    for channel in CHANNELS:
        declared_row = declared_counts.get(channel)
        if not isinstance(declared_row, Mapping):
            raise GroupedDyadicInputError("term_counts channel row is missing")
        for order in ORDERS:
            if _index(declared_row.get(order), f"term_counts.{channel}.{order}") != counts[channel][order]:
                raise GroupedDyadicInputError("declared term_counts differs from input terms")

    if total <= 0 or total > MAX_TERMS:
        raise GroupedDyadicResourceError("input term count is outside the bounded range")
    if not synthetic and (total != EXPECTED_TOTAL_TERMS or counts != EXPECTED_TERM_COUNTS):
        raise GroupedDyadicInputError(
            "non-synthetic K0 receipt must contain the authenticated 4140-term counts"
        )
    return saved_symbolic, cells, counts, total, synthetic


def _validate_term(
    term: Any,
    label: str,
    cell_index: int,
    channel: str,
    ordinary_order: str,
    term_index: int,
) -> tuple[dict[str, Any], tuple[tuple[int, int, int, int], tuple[int, int, int, int]], dict[str, Any]]:
    if not isinstance(term, Mapping):
        raise GroupedDyadicInputError(f"{label} must be a term mapping")
    term_copy = _safe_copy(term, label)
    if term.get("channel") != channel:
        raise GroupedDyadicInputError(f"{label}.channel differs from its container")
    if term.get("ordinary_Z") != ordinary_order:
        raise GroupedDyadicInputError(f"{label}.ordinary_Z differs from its container")
    if _index(term.get("ordinary_Z_order"), label + ".ordinary_Z_order") != ORDERS.index(ordinary_order):
        raise GroupedDyadicInputError(f"{label}.ordinary_Z_order is inconsistent")
    piece = term.get("piece")
    if not isinstance(piece, str) or not piece:
        raise GroupedDyadicInputError(f"{label}.piece must be nonempty")
    signed_factor = _index(term.get("signed_factor"), label + ".signed_factor")
    expected_sign = 1 if signed_factor > 0 else -1
    if _index(term.get("sign"), label + ".sign") != expected_sign:
        raise GroupedDyadicInputError(f"{label}.sign is inconsistent with signed_factor")
    if term.get("coefficient_is_unmodified_saved_binary_record") is not True:
        raise GroupedDyadicInputError(
            f"{label} must retain its unmodified saved coefficient marker"
        )
    if term.get("source_identity_status") not in (None, "source identity unverified"):
        raise GroupedDyadicInputError(f"{label} must not assert source identity")
    powers = term.get("formal_source_powers")
    if not isinstance(powers, Sequence) or isinstance(powers, (str, bytes, bytearray)) or len(powers) != 6:
        raise GroupedDyadicInputError(f"{label}.formal_source_powers must contain six integers")
    if any(type(power) is not int for power in powers):
        raise GroupedDyadicInputError(f"{label}.formal_source_powers must contain plain integers")

    background_name = term.get("background_factor_name")
    if not isinstance(background_name, str) or not background_name:
        raise GroupedDyadicInputError(f"{label}.background_factor_name must be nonempty")
    background_factor = term.get("background_factor")
    if not isinstance(background_factor, Mapping):
        raise GroupedDyadicInputError(f"{label}.background_factor must be a mapping")
    _safe_copy(background_factor, label + ".background_factor")

    rho_factor = term.get("rho_factor")
    if rho_factor is not None:
        _safe_copy(rho_factor, label + ".rho_factor")
    product_factors = term.get("product_factors")
    if not isinstance(product_factors, Sequence) or isinstance(product_factors, (str, bytes, bytearray)):
        raise GroupedDyadicInputError(f"{label}.product_factors must be a sequence")
    if any(not isinstance(item, str) for item in product_factors):
        raise GroupedDyadicInputError(f"{label}.product_factors must contain strings")
    local_factor = term.get("local_factor")
    if local_factor is not None and not isinstance(local_factor, Mapping):
        raise GroupedDyadicInputError(f"{label}.local_factor must be a mapping when present")
    if local_factor is not None:
        _safe_copy(local_factor, label + ".local_factor")

    source_cell = _index(term.get("source_cell"), label + ".source_cell")
    target_cell = _index(term.get("target_cell"), label + ".target_cell")
    if source_cell < 0 or source_cell > cell_index or target_cell != cell_index:
        raise GroupedDyadicInputError(f"{label} source/target cell provenance is invalid")
    source_ids = term.get("source_ids")
    if not isinstance(source_ids, Sequence) or isinstance(source_ids, (str, bytes, bytearray)):
        raise GroupedDyadicInputError(f"{label}.source_ids must be a sequence")
    _safe_copy(source_ids, label + ".source_ids")

    original_term = term.get("original_term")
    if original_term is not None and not isinstance(original_term, Mapping):
        raise GroupedDyadicInputError(f"{label}.original_term must be a mapping when present")
    if original_term is not None:
        _safe_copy(original_term, label + ".original_term")

    coefficient = term.get("coefficient")
    coefficient_interval = _coefficient_interval(coefficient, label + ".coefficient")
    key_payload = {
        "cell_index": cell_index,
        "channel": channel,
        "ordinary_order": ordinary_order,
        "piece": piece,
        "signed_factor": signed_factor,
        "formal_source_powers": _copy(list(powers)),
        "background_factor_name": background_name,
        "background_factor": _copy(background_factor),
        "rho_factor": _copy(rho_factor),
        "product_factors": _copy(list(product_factors)),
        "local_factor": _copy(local_factor),
    }
    _freeze(key_payload, label + ".group_key")
    member = {
        "cell_index": cell_index,
        "channel": channel,
        "ordinary_order": ordinary_order,
        "ordinary_order_index": ORDERS.index(ordinary_order),
        "term_index": term_index,
        "source_cell": source_cell,
        "target_cell": target_cell,
        "source_ids": _copy(source_ids),
        "coefficient": _copy(coefficient),
        "original_term": _copy(original_term),
        "term": term_copy,
    }
    return key_payload, coefficient_interval, member


def _resolve_sum_intervals(value: Any) -> Callable[..., Any]:
    if callable(value):
        return value
    candidate = value.get("sum_intervals") if isinstance(value, Mapping) else getattr(value, "sum_intervals", None)
    if not callable(candidate):
        raise GroupedDyadicInputError("sum_intervals callable is required")
    return candidate


def _validate_result_interval(value: Any, label: str) -> tuple[tuple[int, int, int, int], tuple[int, int, int, int]]:
    if not isinstance(value, tuple) or len(value) != 2:
        raise GroupedDyadicInputError(f"{label} must return a two-endpoint tuple")
    lower = _endpoint_tuple(value[0], label + ".lower", allow_list=False)
    upper = _endpoint_tuple(value[1], label + ".upper", allow_list=False)
    if _compare_signed(lower, upper) > 0:
        raise GroupedDyadicInputError(f"{label} returned reversed endpoints")
    return (lower, upper)


def group(
    saved_symbolic: Any,
    sum_intervals: Any,
    precision: int = DEFAULT_PRECISION_BITS,
) -> dict[str, Any]:
    """Group the saved K0 terms and sum each raw coefficient interval once."""

    working_precision = _precision(precision)
    sum_operation = _resolve_sum_intervals(sum_intervals)
    saved, cells, counts, total, synthetic = _validate_symbolic(saved_symbolic)

    grouped: dict[tuple[Any, ...], dict[str, Any]] = {}
    for cell in cells:
        cell_index = _index(cell.get("index"), "cell.index")
        terms = cell["terms"]
        for channel in CHANNELS:
            rows = terms[channel]
            for ordinary_order in ORDERS:
                sequence = rows[ordinary_order]
                for term_index, term in enumerate(sequence):
                    key_payload, coefficient_interval, member = _validate_term(
                        term,
                        f"cells[{cell_index}].terms.{channel}.{ordinary_order}[{term_index}]",
                        cell_index,
                        channel,
                        ordinary_order,
                        term_index,
                    )
                    key = _freeze(key_payload, "group_key")
                    if key not in grouped:
                        if len(grouped) >= MAX_GROUPS:
                            raise GroupedDyadicResourceError("group count exceeds MAX_GROUPS")
                        grouped[key] = {
                            "group_key": key_payload,
                            "coefficient_inputs": [],
                            "members": [],
                        }
                    grouped[key]["coefficient_inputs"].append(coefficient_interval)
                    grouped[key]["members"].append(member)

    member_total = sum(len(item["members"]) for item in grouped.values())
    if member_total != total:
        raise GroupedDyadicInputError("not every input term entered exactly one group")

    output_groups: list[dict[str, Any]] = []
    for item in grouped.values():
        inputs = tuple(item["coefficient_inputs"])
        try:
            summed = sum_operation(inputs, precision=working_precision)
        except Exception as exc:
            raise GroupedDyadicInputError("sum_intervals failed for a coefficient group") from exc
        result_interval = _validate_result_interval(summed, "sum_intervals result")
        output_groups.append(
            {
                "group_key": _safe_copy(item["group_key"], "group_key"),
                "coefficient_interval": result_interval,
                "member_count": len(item["members"]),
                "members": _safe_copy(item["members"], "group members"),
            }
        )

    unknown_terms = _safe_copy(saved["unknown_terms"], "unknown_terms")
    source_family = _safe_copy(saved.get("source_family"), "source_family")
    source_logs = _safe_copy(saved.get("source_logs"), "source_logs")
    source_receipt = getattr(sum_intervals, "ROUNDING_RECEIPT", None)
    if source_receipt is None and isinstance(sum_intervals, Mapping):
        source_receipt = sum_intervals.get("ROUNDING_RECEIPT")
    if source_receipt is not None:
        source_receipt = _safe_copy(source_receipt, "sum_intervals.ROUNDING_RECEIPT")

    return {
        "schema": SCHEMA,
        "status": "grouped_conditional_dyadic_sum",
        "input_schema": SYMBOLIC_SCHEMA,
        "precision_bits": working_precision,
        "synthetic_only": synthetic,
        "cell_count": 8,
        "expected_total_terms": EXPECTED_TOTAL_TERMS,
        "input_term_count": total,
        "member_count": total,
        "group_count": len(output_groups),
        "input_term_counts": counts,
        "groups": output_groups,
        "source_family": source_family,
        "source_logs": source_logs,
        "source_kernel_scope": saved["source_kernel_scope"],
        "unknown_terms": unknown_terms,
        "unknown_terms_are_not_zero": True,
        "source_identity_status": "source identity unverified",
        "source_identity_verified": False,
        "F0_scale_identity": _safe_copy(saved.get("F0_scale_identity"), "F0_scale_identity"),
        "L_local_binding": "formal L_local retained; no source binding performed",
        "sum_intervals_receipt": source_receipt,
        "members_consumed_exactly_once": True,
        "source_ids_are_provenance_only": True,
        "background_products_computed": False,
        "integration_performed": False,
        "numeric_operator_called": False,
        "point_value_selected": False,
        "N1p_included": False,
        "Kchi_included": False,
        "allword_tail_included": False,
        "appended_to_212_ledger": False,
        "ledger_consulted": False,
        "complete_W1": False,
        "complete_m15_computed": False,
    }


__all__ = [
    "SCHEMA",
    "SYMBOLIC_SCHEMA",
    "DEFAULT_PRECISION_BITS",
    "MAX_PRECISION_BITS",
    "MAX_EXPONENT_BITS",
    "MAX_ENDPOINT_MANTISSA_BITS",
    "MAX_TERMS",
    "MAX_GROUPS",
    "EXPECTED_TERM_COUNTS",
    "EXPECTED_TOTAL_TERMS",
    "GroupedDyadicInputError",
    "GroupedDyadicResourceError",
    "group",
]
