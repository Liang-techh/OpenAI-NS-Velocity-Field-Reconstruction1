"""Synthetic raw-endpoint checks for ``n1_p_sector_signs_v3``.

The packet has all four branches, both ordinary-Z rows, and all eight cells.
It includes strict positive/negative, zero-containing, exact-zero,
error-only, and two distinct ordinary-Z1 derivative records.  All values are
exact endpoint tokens; no coefficient operation is performed here.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any, Callable


SCHEMA = "n1_p_sector_signs_fixture_v3"
BRANCHES = ("G0_P_W0", "G0_P_DeltaW0", "D_P_W0", "D_P_DeltaW0")
ORDERS = (0, 1)
CELLS = tuple(range(8))
POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")


def _endpoint(sign: int, mantissa: int, exponent: int) -> list[int]:
    return [0, 0, 0, 0] if mantissa == 0 else [sign, mantissa, exponent, mantissa.bit_length()]


def _interval(lower: list[int], upper: list[int]) -> dict[str, Any]:
    return {"lower_exact_mpf_tuple": list(lower), "upper_exact_mpf_tuple": list(upper), "directed": True}


POSITIVE = _endpoint(0, 1, 0)
POSITIVE_WIDE = _endpoint(0, 1, 1)
NEGATIVE = _endpoint(1, 1, 0)
NEGATIVE_WIDE = _endpoint(1, 1, 1)
ZERO = _endpoint(0, 0, 0)


def _term(
    source_id: str,
    lower: list[int],
    upper: list[int],
    *,
    branch: str,
    order: int,
    cell: int,
    lead_order: int,
    word_order: int,
    powers: tuple[int, int, int, int],
    uncertainty: bool = False,
    error_only: bool = False,
    error_classes: tuple[str, ...] = ("fixture:signed",),
    index: int = 0,
) -> dict[str, Any]:
    path = ["fixture", "result", "output", "potential_compensation_m15", branch, str(order), str(cell), index]
    classes = list(error_classes)
    origin: dict[str, Any] = {
        "source_ids": [source_id],
        "error_classes": classes,
        "uncertainty": uncertainty,
        "target_cell": cell,
        "lead_order": lead_order,
        "word_order": word_order,
        "source_path": path,
        "receipt": "fixture://n1_p_sector_signs_v3",
    }
    if error_only:
        origin["error_only"] = True
    term: dict[str, Any] = {
        "powers": list(powers),
        "coefficient": _interval(lower, upper),
        "origin": origin,
        "source_ids": [source_id],
        "error_classes": classes,
        "uncertainty": uncertainty,
        "target_cell": cell,
        "lead_order": lead_order,
        "word_order": word_order,
        "receipt": "fixture://n1_p_sector_signs_v3",
        "source_path": path,
        "input_index": index,
    }
    if error_only:
        term.update({"error_only": True, "symmetric_error": True})
    return term


def _packet() -> tuple[dict[str, Any], list[str]]:
    packet: dict[str, Any] = {}
    expected: list[str] = []
    for branch in BRANCHES:
        packet[branch] = {str(order): {} for order in ORDERS}
        for order in ORDERS:
            for cell in CELLS:
                base = f"fixture:{branch}:Z{order}:cell{cell}:base"
                rows = [_term(
                    base, POSITIVE, POSITIVE_WIDE, branch=branch, order=order, cell=cell,
                    lead_order=0, word_order=order, powers=(1, 0, 0, -1),
                    uncertainty=(branch == BRANCHES[0] and order == 0 and cell == 0), index=0,
                )]
                expected.append(base)
                if branch == BRANCHES[0] and order == 0 and cell == 0:
                    rows.extend([
                        _term("fixture:strict-negative", NEGATIVE_WIDE, NEGATIVE, branch=branch, order=order, cell=cell, lead_order=0, word_order=0, powers=(0, 1, 0, -2), index=1),
                        _term("fixture:contains-zero", NEGATIVE, POSITIVE, branch=branch, order=order, cell=cell, lead_order=0, word_order=0, powers=(0, 0, 1, -1), uncertainty=True, error_classes=("fixture:contains-zero",), index=2),
                        _term("fixture:exact-zero", ZERO, ZERO, branch=branch, order=order, cell=cell, lead_order=0, word_order=0, powers=(0, 0, 0, 0), error_classes=("fixture:exact-zero",), index=3),
                        _term("fixture:error-only", NEGATIVE, POSITIVE, branch=branch, order=order, cell=cell, lead_order=0, word_order=0, powers=(0, 0, 0, -3), uncertainty=True, error_only=True, error_classes=("fixture:error-only",), index=4),
                    ])
                    expected.extend(["fixture:strict-negative", "fixture:contains-zero", "fixture:exact-zero", "fixture:error-only"])
                if branch == BRANCHES[0] and order == 1 and cell == 0:
                    rows.extend([
                        _term("fixture:Z1-lead-derivative", NEGATIVE, POSITIVE_WIDE, branch=branch, order=order, cell=cell, lead_order=1, word_order=0, powers=(1, 0, 0, -1), error_classes=("fixture:Z1-lead",), index=1),
                        _term("fixture:Z1-word-derivative", POSITIVE, POSITIVE_WIDE, branch=branch, order=order, cell=cell, lead_order=0, word_order=1, powers=(1, 0, 0, -1), error_classes=("fixture:Z1-word",), index=2),
                    ])
                    expected.extend(["fixture:Z1-lead-derivative", "fixture:Z1-word-derivative"])
                packet[branch][str(order)][str(cell)] = rows
    return packet, expected


def _records(result: Mapping[str, Any]) -> tuple[list[Mapping[str, Any]], bool]:
    branches = result.get("branches")
    if not isinstance(branches, Mapping):
        return [], False
    records: list[Mapping[str, Any]] = []
    complete = True
    for branch in BRANCHES:
        branch_value = branches.get(branch)
        orders = branch_value.get("ordinary_Z", branch_value.get("orders")) if isinstance(branch_value, Mapping) else None
        if not isinstance(orders, Mapping):
            complete = False
            continue
        for order in ORDERS:
            row = orders.get(order, orders.get(str(order)))
            cells = row.get("cells", row.get("target_cells")) if isinstance(row, Mapping) else None
            if not isinstance(cells, Mapping):
                complete = False
                continue
            for cell in CELLS:
                cell_value = cells.get(cell, cells.get(str(cell)))
                rows = cell_value.get("records", cell_value.get("raw_records")) if isinstance(cell_value, Mapping) else None
                if not isinstance(rows, list):
                    complete = False
                    continue
                for record in rows:
                    if isinstance(record, Mapping):
                        records.append(record)
                    else:
                        complete = False
    return records, complete


def _ids(record: Mapping[str, Any]) -> list[str]:
    direct = record.get("source_id")
    if isinstance(direct, str):
        return [direct]
    direct_many = record.get("source_ids")
    if isinstance(direct_many, Sequence) and not isinstance(direct_many, (str, bytes, bytearray)):
        return [str(value) for value in direct_many]
    origin = record.get("origin")
    if isinstance(origin, Mapping):
        source_ids = origin.get("source_ids")
        if isinstance(source_ids, Sequence) and not isinstance(source_ids, (str, bytes, bytearray)):
            return [str(value) for value in source_ids]
    return []


def run_fixture(leaf: Any) -> dict[str, Any]:
    packet, expected_ids = _packet()
    checks = {
        "callable_available": False,
        "call_completed": False,
        "result_mapping": False,
        "complete_branch_Z_cell_roster": False,
        "strict_positive_retained": False,
        "strict_negative_retained": False,
        "contains_zero_retained": False,
        "exact_zero_retained": False,
        "error_only_retained_separately": False,
        "two_Z1_derivatives_distinct": False,
        "each_record_once": False,
        "no_coefficient_arithmetic": False,
        "no_midpoint_or_function_value": False,
        "no_physical_weight_reapplication": False,
    }
    try:
        call = leaf if callable(leaf) else getattr(leaf, "build")
        checks["callable_available"] = callable(call)
        result = call(packet)
        checks["call_completed"] = True
    except BaseException as exc:
        return {"schema": SCHEMA, "passed": False, "checks": checks, "check_count": len(checks), "synthetic_only": True, "error": f"{type(exc).__name__}: {exc}"}
    checks["result_mapping"] = isinstance(result, Mapping)
    if not checks["result_mapping"]:
        return {"schema": SCHEMA, "passed": False, "checks": checks, "check_count": len(checks), "synthetic_only": True}
    records, roster = _records(result)
    seen = [source_id for record in records for source_id in _ids(record)]
    checks["complete_branch_Z_cell_roster"] = roster and all(source_id in seen for source_id in expected_ids)
    checks["each_record_once"] = len(seen) == len(expected_ids) and set(seen) == set(expected_ids)

    def find(source_id: str) -> list[Mapping[str, Any]]:
        return [record for record in records if source_id in _ids(record)]

    def has_class(source_id: str, name: str) -> bool:
        return any(record.get("classification") == name or record.get("sign_class") == name for record in find(source_id))

    checks["strict_positive_retained"] = has_class("fixture:G0_P_W0:Z0:cell0:base", "strictly_positive")
    checks["strict_negative_retained"] = has_class("fixture:strict-negative", "strictly_negative")
    checks["contains_zero_retained"] = has_class("fixture:contains-zero", "contains_zero")
    checks["exact_zero_retained"] = has_class("fixture:exact-zero", "exact_zero")
    errors = find("fixture:error-only")
    checks["error_only_retained_separately"] = bool(errors) and all(record.get("error_only") is True and record.get("excluded_from_signed_counts") is True for record in errors)
    lead = find("fixture:Z1-lead-derivative")
    word = find("fixture:Z1-word-derivative")
    checks["two_Z1_derivatives_distinct"] = bool(lead and word) and all(record.get("lead_order") == 1 and record.get("word_order") == 0 for record in lead) and all(record.get("lead_order") == 0 and record.get("word_order") == 1 for record in word)
    checks["no_coefficient_arithmetic"] = (
        result.get("coefficient_arithmetic_performed") is False
        and result.get("interval_addition_performed") is False
        and result.get("interval_subtraction_performed") is False
        and result.get("interval_multiplication_performed") is False
        and result.get("integer_count_aggregation_only") is True
        and result.get("raw_endpoint_tuples_used") is True
    )
    checks["no_midpoint_or_function_value"] = result.get("midpoint_used", False) is False and result.get("function_value_evaluated", False) is False and result.get("coefficient_selected", False) is False and result.get("exact_binary_endpoint_comparison") is True
    checks["no_physical_weight_reapplication"] = result.get("physical_weight_reapplied", False) is False and result.get("physical_scale_reapplied", False) is False and result.get("physical_measure_reapplied", False) is False and result.get("physical_field_claimed", False) is False
    return {
        "schema": SCHEMA,
        "passed": all(checks.values()),
        "checks": checks,
        "check_count": len(checks),
        "synthetic_only": True,
        "raw_endpoint_tokens": True,
        "expected_record_count": len(expected_ids),
        "error_only_is_not_signed": True,
        "ordinary_Z1_derivatives_checked": True,
        "branch_and_cell_roster_checked": True,
        "records_exactly_once_checked": True,
    }


def run(leaf: Any) -> dict[str, Any]:
    return run_fixture(leaf)

