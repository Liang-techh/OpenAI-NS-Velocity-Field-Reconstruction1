"""Small independent fixture for the K0 grouped dyadic consumer.

The fixture builds three JSON-shaped symbolic terms locally and supplies the
dyadic summation callable separately.  Two terms have exactly the same group
key and one differs only in formal source powers, so the expected result is
two groups with member counts two and one.  No producer, source receipt,
background arithmetic, integration, or archive is accessed.
"""

from __future__ import annotations

from copy import deepcopy
from typing import Any


SCHEMA = "original_k0_grouped_dyadic_fixture_v1"
SYMBOLIC_SCHEMA = "original_k0_symbolic_m15_cells_v1"


def _member(value: Any, name: str, default: Any = None) -> Any:
    if isinstance(value, dict):
        return value.get(name, default)
    return getattr(value, name, default)


def _raw_interval(lower: int, upper: int) -> dict[str, Any]:
    return {
        "lower_exact_mpf_tuple": [0, lower, 0, lower.bit_length()],
        "upper_exact_mpf_tuple": [0, upper, 0, upper.bit_length()],
        "directed": True,
    }


def _term(index: int, coefficient: dict[str, Any], powers: list[int]) -> dict[str, Any]:
    original = {
        "coefficient": deepcopy(coefficient),
        "F0_power": powers[1],
        "Lambda_power": powers[3],
        "source_cell": 0,
        "target_cell": 0,
        "matrix_sector_identity": f"fixture-matrix-{index}",
        "W0_sector_identity": f"fixture-w0-{index}",
        "source_id": f"fixture-original-{index}",
    }
    return {
        "channel": "axial",
        "ordinary_Z": "Z0",
        "ordinary_Z_order": 0,
        "piece": "U0*Uz1_0",
        "expression": "2*U0*Uz1_0",
        "signed_factor": 2,
        "sign": 1,
        "coefficient": deepcopy(coefficient),
        "coefficient_is_unmodified_saved_binary_record": True,
        "original_term": original,
        "saved_field": "Uz1",
        "saved_Z": "Z0",
        "saved_formal_source_powers": list(powers[:4]),
        "formal_source_powers": list(powers),
        "background_factor_name": "U_lead_Z0",
        "background_factor": {
            "lower_exact_mpf_tuple": [0, 7, 0, 3],
            "upper_exact_mpf_tuple": [0, 8, 0, 4],
            "directed": True,
            "source_id": "fixture-background-U-Z0",
        },
        "target_cell": 0,
        "source_cell": 0,
        "matrix_sector_identity": original["matrix_sector_identity"],
        "W0_sector_identity": original["W0_sector_identity"],
        "source_indices": [0],
        "source_cell_provenance": {
            "target_cell": 0,
            "source_cell": 0,
            "source_indices": [0],
            "triangular_diagonal_source": 0,
        },
        "source_ids": [f"fixture-source-{index}"],
        "source_identity_status": "source identity unverified",
        "product_factors": ["U_lead_Z0", "Uz1_0"],
        "unknown": False,
    }


def _saved_symbolic() -> dict[str, Any]:
    first = _term(0, _raw_interval(1, 2), [0, 1, 0, -2, 0, 0])
    second = _term(1, _raw_interval(3, 4), [0, 1, 0, -2, 0, 0])
    different_power = _term(2, _raw_interval(5, 5), [0, 2, 0, -2, 0, 0])
    cells: list[dict[str, Any]] = []
    for index in range(8):
        terms = {
            "axial": {"Z0": [], "Z1": []},
            "swirl": {"Z0": [], "Z1": []},
        }
        if index == 0:
            terms["axial"]["Z0"] = [first, second, different_power]
        cells.append({"index": index, "terms": terms})
    return {
        "schema": SYMBOLIC_SCHEMA,
        "synthetic_only": True,
        "source_family": {"fixture": "synthetic source family"},
        "source_logs": {"fixture": "synthetic source logs"},
        "source_kernel_scope": "synthetic K0 only",
        "term_counts": {
            "axial": {"Z0": 3, "Z1": 0},
            "swirl": {"Z0": 0, "Z1": 0},
        },
        "unknown_terms": [
            {"name": "N1p", "status": "unknown", "reason": "fixture"},
            {"name": "Kchi", "status": "unknown", "reason": "fixture"},
            {"name": "allword_tail", "status": "unknown", "reason": "fixture"},
        ],
        "unknown_terms_are_not_zero": True,
        "source_identity_status": "source identity unverified",
        "source_identity_verified": False,
        "F0_scale_identity": {"F0_anchor": "independent", "F0_local": "independent", "identified": False},
        "N1p_included": False,
        "Kchi_included": False,
        "allword_tail_included": False,
        "appended_to_212_ledger": False,
        "cells": cells,
    }


def _sum_operation(dyad: Any) -> Any:
    if callable(dyad):
        return dyad
    return _member(dyad, "sum_intervals")


def _group_operation(leaf: Any) -> Any:
    return _member(leaf, "group")


def _reject(group_operation: Any, saved: Any, sum_operation: Any) -> bool:
    try:
        group_operation(saved, sum_operation, precision=8)
    except Exception:
        return True
    return False


def run_fixture(leaf: Any, dyad: Any) -> dict[str, Any]:
    """Run three-term grouping, provenance and fail-closed checks."""

    group_operation = _group_operation(leaf)
    sum_operation = _sum_operation(dyad)
    checks: dict[str, bool] = {
        "group_api_callable": callable(group_operation),
        "sum_api_callable": callable(sum_operation),
    }
    if not all(checks.values()):
        return {
            "schema": SCHEMA,
            "passed": False,
            "checks": checks,
            "test_count": len(checks),
            "synthetic_only": True,
            "source_provider_called": False,
            "source_graph_walked": False,
        }

    saved = _saved_symbolic()
    before = deepcopy(saved)
    result = group_operation(saved, sum_operation, precision=8)
    groups = result.get("groups", [])
    sizes = sorted(group.get("member_count") for group in groups)
    two_member = next((item for item in groups if item.get("member_count") == 2), None)
    one_member = next((item for item in groups if item.get("member_count") == 1), None)

    checks.update(
        {
            "result_schema": result.get("schema") == "original_k0_grouped_dyadic_v1",
            "two_exact_groups": result.get("group_count") == 2 and sizes == [1, 2],
            "three_members_consumed": result.get("member_count") == 3
            and result.get("input_term_count") == 3,
            "channel_order_count_retained": result.get("input_term_counts")
            == {"axial": {"Z0": 3, "Z1": 0}, "swirl": {"Z0": 0, "Z1": 0}},
            "same_factor_coefficients_summed": two_member is not None
            and two_member.get("coefficient_interval")
            == ((0, 1, 2, 1), (0, 3, 1, 2)),
            "different_formal_power_not_merged": one_member is not None
            and one_member.get("coefficient_interval")
            == ((0, 5, 0, 3), (0, 5, 0, 3))
            and one_member.get("group_key", {}).get("formal_source_powers")
            == [0, 2, 0, -2, 0, 0],
            "term_indices_and_full_terms_retained": two_member is not None
            and [item.get("term_index") for item in two_member.get("members", [])] == [0, 1]
            and all(item.get("term", {}).get("channel") == "axial" for item in two_member.get("members", [])),
            "unknown_terms_and_scope_retained": result.get("unknown_terms_are_not_zero") is True
            and {item.get("name") for item in result.get("unknown_terms", [])}
            >= {"N1p", "Kchi", "allword_tail"}
            and result.get("N1p_included") is False
            and result.get("Kchi_included") is False
            and result.get("allword_tail_included") is False,
            "source_identity_unverified": result.get("source_identity_verified") is False
            and result.get("source_ids_are_provenance_only") is True,
            "no_background_product_or_integration": result.get("background_products_computed") is False
            and result.get("integration_performed") is False
            and result.get("appended_to_212_ledger") is False,
            "saved_input_unmodified": saved == before,
        }
    )

    bad = deepcopy(saved)
    bad["cells"][0]["terms"]["axial"]["Z0"][0]["coefficient"]["lower_exact_mpf_tuple"] = [0, 1, 0, 2]
    checks["malformed_coefficient_rejected"] = _reject(group_operation, bad, sum_operation)
    bad = deepcopy(saved)
    bad["source_identity_verified"] = True
    checks["verified_source_claim_rejected"] = _reject(group_operation, bad, sum_operation)

    return {
        "schema": SCHEMA,
        "passed": bool(checks) and all(checks.values()),
        "checks": checks,
        "test_count": len(checks),
        "synthetic_only": True,
        "source_provider_called": False,
        "source_graph_walked": False,
        "integration_performed": False,
        "source_identity_verified": False,
    }


__all__ = ["SCHEMA", "run_fixture"]
