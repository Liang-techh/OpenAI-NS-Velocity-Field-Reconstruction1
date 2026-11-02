"""Synthetic checks for :mod:`row5_g5_binding_v1`.

The fixture uses exact JSON interval dictionaries rather than decimal values.
It creates sixteen source-only ``g5`` occurrences in each ordinary-Z row and
checks both the successful 32-entry binding and fail-closed mutations: a
duplicate, a missing occurrence, a wrong row, an identity collision with the
same interval, scalar/power changes, endpoint changes, and baseline emission.
"""

from __future__ import annotations

from typing import Any, Mapping


SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
OUTPUT_ROWS = ("h5", "h5_Z1")
EXPECTED_G5_PER_ROW = 16

__all__ = ["SOURCE_POWER_BASIS", "OUTPUT_ROWS", "EXPECTED_G5_PER_ROW", "run_fixture", "run"]


def _clone(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {key: _clone(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_clone(item) for item in value]
    return value


def _mpf_tuple(value: int) -> list[int]:
    sign = 1 if value < 0 else 0
    magnitude = abs(int(value))
    if magnitude == 0:
        return [0, 0, 0, 0]
    exponent = 0
    while magnitude % 2 == 0:
        magnitude //= 2
        exponent += 1
    return [sign, magnitude, exponent, magnitude.bit_length()]


def _interval(index: int, row_order: int) -> dict[str, Any]:
    lower = -(1000 + row_order * 100 + index * 3)
    upper = lower + 2
    return {
        "lower_exact_mpf_tuple": _mpf_tuple(lower),
        "upper_exact_mpf_tuple": _mpf_tuple(upper),
        "directed": True,
        "function_value": False,
    }


def _term(row: str, index: int) -> dict[str, Any]:
    row_order = 0 if row == "h5" else 1
    source_id = f"fixture:g5:{row}:{index}"
    source_path = ["fixture", "affine", row, index]
    powers = [row_order, index % 5 - 2, (index + row_order) % 4, -index % 3]
    interval = _interval(index, row_order)
    role = "source_only_g5" if row == "h5" else "source_only_g5_Z1"
    return {
        "coefficient": _clone(interval),
        "coefficient_interval": _clone(interval),
        "formal_powers": list(powers),
        "powers": list(powers),
        "power_basis": list(SOURCE_POWER_BASIS),
        "scalar": 1,
        "state_reference": None,
        "state_ref": None,
        "state": None,
        "reference": None,
        "source_only": True,
        "source_only_g5": True,
        "ordinary_Z_order": row_order,
        "source_ordinary_Z_order": row_order,
        "source_family": "g5",
        "source_key": "g5",
        "matrix_name": None,
        "role": role,
        "term_role": role,
        "source_ids": [source_id],
        "source_paths": [_clone(source_path)],
        "source_path": _clone(source_path),
        "uncertainty": index % 2 == 0,
        "error_only": False,
        "symmetric_error": False,
        "signed_field": True,
        "error_classes": [f"fixture:error:{row}:{index}"],
        "source_error_classes": [f"fixture:source-error:{row}:{index}"],
        "midpoint_values_used": False,
    }


def _expression() -> dict[str, Any]:
    return {
        "schema": "row5_affine_forcing_v1",
        "source_power_basis": list(SOURCE_POWER_BASIS),
        "rows": {row: [_term(row, index) for index in range(EXPECTED_G5_PER_ROW)] for row in OUTPUT_ROWS},
    }


def _firstword(expression: Mapping[str, Any]) -> dict[str, Any]:
    source_nodes: dict[str, dict[str, Any]] = {}
    rows: dict[str, list[dict[str, Any]]] = {row: [] for row in OUTPUT_ROWS}
    for row in OUTPUT_ROWS:
        for index, term in enumerate(expression["rows"][row]):
            node_id = f"expression:{row}:{index}"
            interval = term["coefficient"]
            node = _clone(term)
            node.update(
                {
                    "node_id": node_id,
                    "node_kind": "source_expression",
                    "coefficient_endpoints": {
                        "lower": list(interval["lower_exact_mpf_tuple"]),
                        "upper": list(interval["upper_exact_mpf_tuple"]),
                    },
                    "provenance_ids": [term["source_ids"][0]],
                    "parent_input_paths": [_clone(term["source_path"])],
                    "source_ids": [term["source_ids"][0]],
                    "source_paths": [_clone(term["source_path"])],
                    "original_node": _clone(term),
                    "function_value": False,
                }
            )
            source_nodes[node_id] = node
            rows[row].append(
                {
                    "product_id": f"product:{row}:{index}:source-only",
                    "row": row,
                    "source_node_id": node_id,
                    "state_node_id": None,
                    "source_only": True,
                    "state_reference": None,
                    "state_index": None,
                    "source_index": index,
                    "scalar": 1,
                    "formal_powers": list(term["formal_powers"]),
                    "powers": list(term["formal_powers"]),
                    "coefficient_enclosure": _clone(interval),
                    "coefficient_is_function_value": False,
                    "scalar_applied_to_coefficient": False,
                    "term_semantics": "scalar * coefficient_enclosure * formal_monomial; scalar deferred",
                    "directed": True,
                    "function_value": False,
                    "provenance": {"source_node_id": node_id, "state_node_id": None, "cancellation_performed": False},
                }
            )
    return {
        "schema": "row5_firstword_product_v4",
        "source_schema": "row5_affine_forcing_v1",
        "source_power_basis": list(SOURCE_POWER_BASIS),
        "rows": rows,
        "source_nodes": source_nodes,
        "state_nodes": {},
        "unknown_residual": {"symbol": "R", "preserved": True, "included_in_products": False, "added_to_ledger": False},
    }


def _baseline(expression: Mapping[str, Any]) -> dict[str, Any]:
    source_nodes: dict[str, dict[str, Any]] = {}
    skipped: list[dict[str, Any]] = []
    rows = {
        row: {"cell-left": {"core": [], "local": []}, "cell-right": {"core": [], "local": []}}
        for row in OUTPUT_ROWS
    }
    for row in OUTPUT_ROWS:
        for index, term in enumerate(expression["rows"][row]):
            node_id = f"expression:{row}:{index}"
            node = _clone(term)
            node.update(
                {
                    "node_id": node_id,
                    "node_kind": "row5_affine_source",
                    "source_interval": _clone(term["coefficient"]),
                    "source_interval_endpoints": _clone(term["coefficient"]),
                    "coefficient_endpoints": _clone(term["coefficient"]),
                    "provenance_ids": [term["source_ids"][0]],
                    "source_ids": [term["source_ids"][0]],
                    "source_id": term["source_ids"][0],
                    "parent_input_paths": [_clone(term["source_path"])],
                    "source_paths": [_clone(term["source_path"])],
                    "original_term": _clone(term),
                    "function_value": False,
                }
            )
            source_nodes[node_id] = node
            skipped.append(
                {
                    "source_node_id": node_id,
                    "row": row,
                    "source_index": index,
                    "source_family": "g5",
                    "retained_by": "row5_firstword_product_v4",
                    "product_emitted": False,
                }
            )
    return {
        "schema": "row5_baseline_product_v4",
        "source_schema": "row5_affine_forcing_v1",
        "source_power_basis": list(SOURCE_POWER_BASIS),
        "rows": rows,
        "source_nodes": source_nodes,
        "skipped_g5": skipped,
        "cells": [{"index": 0, "cell": "cell-left"}, {"index": 1, "cell": "cell-right"}],
        "R_hi": {
            "symbol": "R_hi",
            "definition": "W - (W0_baseline + W_firstword)",
            "preserved": True,
            "known": False,
            "included_in_products": False,
            "added_to_ledger": False,
        },
        "unknown_residual": {
            "symbol": "R_hi",
            "preserved": True,
            "known": False,
            "included_in_products": False,
            "added_to_ledger": False,
        },
    }


def _callable(value: Any) -> Any:
    candidate = getattr(value, "build", value)
    if not callable(candidate):
        raise TypeError("build must be callable or expose build")
    return candidate


def _rejected(call: Any) -> bool:
    try:
        call()
    except Exception:
        return True
    return False


def run_fixture(build: Any, c: Any = None) -> dict[str, Any]:
    build_fn = _callable(build)
    expression = _expression()
    firstword = _firstword(expression)
    baseline = _baseline(expression)
    affine_wrapper = {"receipt_path": "fixture/affine_receipt.json.gz", "expression": expression}
    firstword_wrapper = {"receipt_path": "fixture/firstword_receipt.json.gz", "result": firstword}
    baseline_wrapper = {"receipt_path": "fixture/baseline_receipt.json.gz", "products": baseline}
    try:
        result = build_fn(c, affine_wrapper, firstword_wrapper, baseline_wrapper)
    except Exception as exc:
        return {
            "schema": "row5_g5_binding_fixture_v1",
            "passed": False,
            "checks": {"build_completed": False},
            "check_count": 1,
            "error": f"{type(exc).__name__}: {exc}",
            "synthetic_only": True,
            "implementation_imported": False,
        }

    checks: dict[str, bool] = {}
    checks["build_completed"] = isinstance(result, Mapping)
    checks["schema"] = result.get("schema") == "row5_g5_binding_v1"
    checks["mapping_count"] = result.get("mapping_count") == 32 and len(result.get("mappings", [])) == 32
    checks["row_counts"] = all(
        len(result.get("rows", {}).get(row, [])) == EXPECTED_G5_PER_ROW for row in OUTPUT_ROWS
    )
    checks["source_node_refs"] = all(
        mapping.get("source_node_references", {}).get("firstword", {}).get("node_id")
        == mapping.get("firstword_source_node_id")
        and mapping.get("source_node_references", {}).get("baseline", {}).get("node_id")
        == mapping.get("baseline_source_node_id")
        for mapping in result.get("mappings", [])
    )
    checks["exact_interval_payloads"] = all(
        mapping.get("directed_exact_interval_endpoints", {}).get("directed") is True
        and len(mapping.get("directed_exact_interval_endpoints", {}).get("lower_exact_mpf_tuple", [])) == 4
        and len(mapping.get("directed_exact_interval_endpoints", {}).get("upper_exact_mpf_tuple", [])) == 4
        for mapping in result.get("mappings", [])
    )
    checks["cell_reference_plan"] = all(
        len(item.get("references", [])) == 6
        and sum(ref.get("branch") == "firstword" for ref in item["references"]) == 2
        and sum(ref.get("branch") == "baseline" for ref in item["references"]) == 4
        and item.get("R_hi", {}).get("known") is False
        for item in result.get("cell_combinations", [])
    ) and len(result.get("cell_combinations", [])) == 2
    checks["unknowns_retained"] = (
        result.get("unknowns", {}).get("firstword_existing_W", {}).get("preserved") is True
        and result.get("unknowns", {}).get("R_hi", {}).get("known") is False
    )
    checks["ownership_not_checked"] = result.get("conditions", {}).get("physical_ownership_validated") is False

    duplicate = _clone(firstword_wrapper)
    duplicate["result"]["rows"]["h5"].append(_clone(duplicate["result"]["rows"]["h5"][0]))
    checks["duplicate_rejected"] = _rejected(lambda: build_fn(c, affine_wrapper, duplicate, baseline_wrapper))

    g5_feedback = _clone(firstword_wrapper)
    rogue = _clone(g5_feedback["result"]["rows"]["h5"][0])
    rogue["source_only"] = False
    rogue["state_node_id"] = "fixture:unexpected-state"
    g5_feedback["result"]["rows"]["h5"].append(rogue)
    checks["g5_feedback_reuse_rejected"] = _rejected(lambda: build_fn(c, affine_wrapper, g5_feedback, baseline_wrapper))

    wrong_product_id = _clone(firstword_wrapper)
    wrong_product_id["result"]["rows"]["h5"][0]["product_id"] = "fixture:wrong-product"
    checks["product_id_rejected"] = _rejected(lambda: build_fn(c, affine_wrapper, wrong_product_id, baseline_wrapper))
    wrong_source_index = _clone(firstword_wrapper)
    rogue = _clone(wrong_source_index["result"]["rows"]["h5"][0])
    rogue["source_index"] = 91
    wrong_source_index["result"]["rows"]["h5"].append(rogue)
    checks["extra_wrong_index_rejected"] = _rejected(lambda: build_fn(c, affine_wrapper, wrong_source_index, baseline_wrapper))
    bad_original_metadata = _clone(firstword_wrapper)
    bad_original_metadata["result"]["source_nodes"]["expression:h5:0"]["original_node"]["uncertainty"] = False
    checks["original_metadata_rejected"] = _rejected(lambda: build_fn(c, affine_wrapper, bad_original_metadata, baseline_wrapper))

    missing = _clone(firstword_wrapper)
    missing["result"]["rows"]["h5"].pop(0)
    checks["missing_rejected"] = _rejected(lambda: build_fn(c, affine_wrapper, missing, baseline_wrapper))

    wrong_row = _clone(firstword_wrapper)
    wrong_row["result"]["rows"]["h5"][0]["row"] = "h5_Z1"
    checks["wrong_row_rejected"] = _rejected(lambda: build_fn(c, affine_wrapper, wrong_row, baseline_wrapper))

    same_interval_different_source = _clone(firstword_wrapper)
    same_interval_different_source["result"]["source_nodes"]["expression:h5:0"]["source_ids"] = ["fixture:wrong-source"]
    checks["same_interval_different_source_rejected"] = _rejected(
        lambda: build_fn(c, affine_wrapper, same_interval_different_source, baseline_wrapper)
    )

    reversed_interval = _clone(affine_wrapper)
    bad_term = reversed_interval["expression"]["rows"]["h5"][0]
    for alias in ("coefficient", "coefficient_interval"):
        interval = bad_term[alias]
        interval["lower_exact_mpf_tuple"], interval["upper_exact_mpf_tuple"] = interval["upper_exact_mpf_tuple"], interval["lower_exact_mpf_tuple"]
    checks["reversed_interval_rejected"] = _rejected(lambda: build_fn(c, reversed_interval, firstword_wrapper, baseline_wrapper))

    bad_scalar = _clone(firstword_wrapper)
    bad_scalar["result"]["rows"]["h5"][0]["scalar"] = 2
    checks["scalar_rejected"] = _rejected(lambda: build_fn(c, affine_wrapper, bad_scalar, baseline_wrapper))

    bad_powers = _clone(baseline_wrapper)
    bad_powers["products"]["source_nodes"]["expression:h5:0"]["formal_powers"][0] += 1
    bad_powers["products"]["source_nodes"]["expression:h5:0"]["powers"][0] += 1
    checks["powers_rejected"] = _rejected(lambda: build_fn(c, affine_wrapper, firstword_wrapper, bad_powers))

    bad_endpoint = _clone(firstword_wrapper)
    bad_endpoint["result"]["rows"]["h5"][0]["coefficient_enclosure"]["upper_exact_mpf_tuple"][1] += 1
    checks["endpoint_rejected"] = _rejected(lambda: build_fn(c, affine_wrapper, bad_endpoint, baseline_wrapper))

    baseline_duplicate = _clone(baseline_wrapper)
    baseline_duplicate["products"]["skipped_g5"].append(_clone(baseline_duplicate["products"]["skipped_g5"][0]))
    checks["baseline_duplicate_rejected"] = _rejected(
        lambda: build_fn(c, affine_wrapper, firstword_wrapper, baseline_duplicate)
    )

    baseline_missing = _clone(baseline_wrapper)
    baseline_missing["products"]["skipped_g5"].pop()
    checks["baseline_missing_rejected"] = _rejected(
        lambda: build_fn(c, affine_wrapper, firstword_wrapper, baseline_missing)
    )

    baseline_wrong_row = _clone(baseline_wrapper)
    baseline_wrong_row["products"]["skipped_g5"][0]["row"] = "h5_Z1"
    checks["baseline_wrong_row_rejected"] = _rejected(
        lambda: build_fn(c, affine_wrapper, firstword_wrapper, baseline_wrong_row)
    )

    baseline_product = _clone(baseline_wrapper)
    baseline_product["products"]["rows"]["h5"]["cell-left"]["core"].append(
        {"source_node_id": "expression:h5:0"}
    )
    checks["baseline_product_rejected"] = _rejected(
        lambda: build_fn(c, affine_wrapper, firstword_wrapper, baseline_product)
    )

    return {
        "schema": "row5_g5_binding_fixture_v1",
        "passed": all(checks.values()),
        "check_count": len(checks),
        "checks": checks,
        "synthetic_only": True,
        "fresh_values_and_ids": True,
        "exact_json_intervals": True,
        "duplicate_missing_wrong_row_identity_scalar_power_endpoint_guards": True,
        "cell_reference_plan_checked": True,
        "ownership_not_checked": True,
        "implementation_imported": False,
    }


def run(build: Any, c: Any = None) -> dict[str, Any]:
    return run_fixture(build, c)
