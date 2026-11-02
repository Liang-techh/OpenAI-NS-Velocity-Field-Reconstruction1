"""Fresh synthetic checks for :mod:`core_z2_seam_bridge_v3`.

The fixture is source-only.  It constructs complete result-shaped mappings in
memory and never imports a producer, reads a receipt, invokes old mathematics,
or performs a numerical rerun.  There are twenty-four old row-5 unknowns
(eight source cells times B1 columns 1, 2, and 3), two profile state terms for
each seam state field, and two output fields (``F1`` and ``d_x_F1``).
"""

from __future__ import annotations

import traceback
from collections.abc import Mapping, Sequence
from typing import Any


SCHEMA = "core_z2_seam_bridge_fixture_v3"
OLD_SCHEMA = "core_secondword_partial_v2"
PROFILE_SCHEMA = "core_firstword_z2_profile_v2"
SEAM_SCHEMA = PROFILE_SCHEMA + ":missing_seam_v1"
SEAM_FIELDS = ("F1", "d_x_F1")
STATE_FIELDS = ("F1", "Uz1", "K1")
STATE_FIELD_BY_COLUMN = {1: "F1", 2: "Uz1", 3: "K1"}
CELL_COUNT = 8

__all__ = [
    "SCHEMA",
    "OLD_SCHEMA",
    "PROFILE_SCHEMA",
    "SEAM_SCHEMA",
    "synthetic_inputs",
    "synthetic_low_output",
    "run_fixture",
    "run",
]


def _is_sequence(value: Any) -> bool:
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray))


def _clone(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {key: _clone(item) for key, item in value.items()}
    if _is_sequence(value):
        return [_clone(item) for item in value]
    return value


def _mpf(value: int, exponent: int = 0) -> list[int]:
    """Canonical finite exact MPF tuple for an integer times ``2**exponent``."""

    if value == 0:
        return [0, 0, 0, 0]
    sign = 1 if value < 0 else 0
    mantissa = abs(value)
    while mantissa % 2 == 0:
        mantissa //= 2
        exponent += 1
    return [sign, mantissa, exponent, mantissa.bit_length()]


def _iv(lower: int, upper: int, *, exponent: int = 0) -> dict[str, Any]:
    return {
        "lower_exact_mpf_tuple": _mpf(lower, exponent),
        "upper_exact_mpf_tuple": _mpf(upper, exponent),
        "directed": True,
    }


def _cell_t(index: int) -> dict[str, Any]:
    # [index/8, (index+1)/8] has exact binary endpoints.
    return {
        "lower_exact_mpf_tuple": _mpf(index, -3),
        "upper_exact_mpf_tuple": _mpf(index + 1, -3),
        "directed": True,
    }


def _path(source_cell: int, matrix_key: str, matrix_term_index: int = 0) -> list[Any]:
    return [
        "experiments/root_st073/n1_core_low_word_v3_receipt.json.gz",
        "result",
        "output",
        "matrix_cells",
        source_cell,
        "B1",
        "0",
        matrix_key,
        matrix_term_index,
    ]


def _state_term(source_cell: int, field: str, state_term_index: int) -> dict[str, Any]:
    term_id = f"fixture_v3:profile:cell:{source_cell}:Z2:{field}:{state_term_index}"
    return {
        "term_id": term_id,
        "input_id": term_id + ":input",
        "source_id": term_id,
        "source_ids": [term_id, f"fixture_v3:profile:cell:{source_cell}"],
        "F0_power": source_cell + state_term_index - 2,
        "Lambda_power": -1 - state_term_index,
        "coefficient": _iv(50 + source_cell + state_term_index, 51 + source_cell + state_term_index),
        "error_classes": [f"fixture_v3:profile:error:{source_cell}:{field}:{state_term_index}"],
        "source_error_classes": [f"fixture_v3:profile:error:{source_cell}:{field}:{state_term_index}"],
        "uncertainty": True,
        "signed_field": True,
        "formal_atoms_materialized": False,
        "target_cell": source_cell,
        "source_cell": max(0, source_cell - 1) if state_term_index == 0 else source_cell,
        "source_cell_index": max(0, source_cell - 1) if state_term_index == 0 else source_cell,
    }


def _matrix_term(source_cell: int, matrix_key: str) -> dict[str, Any]:
    term_id = f"fixture_v3:matrix:cell:{source_cell}:B1:0:{matrix_key}:0"
    path = _path(source_cell, matrix_key)
    return {
        "term_id": term_id,
        "input_id": term_id,
        "source_id": term_id,
        "source_ids": [term_id],
        "F0_power": source_cell - 1,
        "Lambda_power": -2,
        "coefficient": _iv(60 + source_cell, 61 + source_cell),
        "error_classes": [f"fixture_v3:matrix:error:{source_cell}:{matrix_key}"],
        "source_error_classes": [f"fixture_v3:matrix:error:{source_cell}:{matrix_key}"],
        "uncertainty": True,
        "signed_field": True,
        "formal_atoms_materialized": False,
        "source_cell": source_cell,
        "source_cell_index": source_cell,
        "receipt_input_path": path,
    }


def synthetic_inputs() -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    """Return ``(old_result, profile_result, seam_result)`` fresh each call."""

    old_unknowns: list[dict[str, Any]] = []
    profile_cells: list[dict[str, Any]] = []
    seam_terms = {field: [] for field in SEAM_FIELDS}
    seam_kernels: dict[str, dict[str, Any]] = {}
    old_kernels: dict[str, dict[str, Any]] = {}
    seam_source_records: list[dict[str, Any]] = []
    old_source_records: list[dict[str, Any]] = []
    target_x = _iv(2, 2)
    target_t = _iv(1, 1)

    for source_cell in range(CELL_COUNT):
        t = _cell_t(source_cell)
        masses = {
            "F1": _iv(1000 + source_cell, 1001 + source_cell),
            "Fx": _iv(2000 + source_cell, 2001 + source_cell),
            "d_x_F1": _iv(2000 + source_cell, 2001 + source_cell),
        }
        old_kernels[str(source_cell)] = {
            "source_cell": source_cell,
            "source_cell_index": source_cell,
            "t": _clone(t),
            "F1": _clone(masses["F1"]),
            "d_x_F1": _clone(masses["d_x_F1"]),
            "exact": True,
            "positive": True,
        }
        seam_kernels[str(source_cell)] = {
            "source_cell": source_cell,
            "source_cell_index": source_cell,
            "t": _clone(t),
            "F1": _clone(masses["F1"]),
            "Fx": _clone(masses["Fx"]),
            "d_x_F1": _clone(masses["d_x_F1"]),
            "exact": True,
            "positive": True,
        }

        state_family: dict[str, list[dict[str, Any]]] = {}
        for field in STATE_FIELDS:
            state_family[field] = [_state_term(source_cell, field, state_term_index) for state_term_index in range(2)]
        profile_cells.append(
            {
                "index": source_cell,
                "t": _clone(t),
                "x": _iv(source_cell, source_cell + 1, exponent=-2),
                "rho": _iv(1 + source_cell, 2 + source_cell),
                "W1_Z2": state_family,
                "W1_by_order": {"Z2": {field: _clone(terms) for field, terms in state_family.items()}},
                "finite_K0_condition_package": True,
            }
        )

        for column in (1, 2, 3):
            matrix_key = f"5,{column}"
            matrix_term = _matrix_term(source_cell, matrix_key)
            path = _path(source_cell, matrix_key)
            matrix_ref = {
                "role": "matrix_source_node",
                "node_id": matrix_term["term_id"],
                "ids": [matrix_term["term_id"]],
                "source_ids": [matrix_term["term_id"]],
                "source_paths": [_clone(path)],
                "error_classes": list(matrix_term["error_classes"]),
                "source_cell": source_cell,
                "source_cell_index": source_cell,
                "term_index": 0,
            }
            unknown_id = f"fixture_v3:old:cell:{source_cell}:unknown:{column}"
            old_unknowns.append(
                {
                    "unknown_id": unknown_id,
                    "unknown_kind": "missing_B1_W1_Z2",
                    "reason": "fixture omitted B1[0]*W[2]",
                    "origin": "B1[0]*W[2] (unknown)",
                    "ordinary_Z_order": 1,
                    "matrix_name": "B1",
                    "matrix_order": "0",
                    "matrix_key": matrix_key,
                    "row": 5,
                    "column": column,
                    "state_field": STATE_FIELD_BY_COLUMN[column],
                    "state_order": "2",
                    "source_cell": source_cell,
                    "source_cell_index": source_cell,
                    "source_t": _clone(t),
                    "target_x": _clone(target_x),
                    "target_t": _clone(target_t),
                    "matrix_term_index": 0,
                    "matrix_node_ref": matrix_ref,
                    "source_node_ref": _clone(matrix_ref),
                    "source_ids": [matrix_term["term_id"]],
                    "source_paths": [_clone(path)],
                    "error_classes": ["fixture_v3:old:unknown", *matrix_term["error_classes"]],
                    "source_error_classes": list(matrix_term["error_classes"]),
                    "state_error_classes": ["unknown:firstword_Z2_not_saved"],
                    "uncertainty": True,
                    "error_only": True,
                    "included": False,
                    "materialized": False,
                    "coefficient": None,
                    "formal_powers": None,
                    "kernel_mass": None,
                    "no_cancellation": True,
                }
            )

            state_field = STATE_FIELD_BY_COLUMN[column]
            for state_term_index, state_term in enumerate(state_family[state_field]):
                for target_field in SEAM_FIELDS:
                    product_id = f"fixture_v3:seam:cell:{source_cell}:{matrix_key}:0:{state_field}:{state_term_index}:{target_field}"
                    product = {
                        "term_id": product_id,
                        "input_id": product_id,
                        "input_ids": [matrix_term["term_id"], state_term["term_id"]],
                        "source_ids": [matrix_term["term_id"], state_term["term_id"]],
                        "F0_power": matrix_term["F0_power"] + state_term["F0_power"],
                        "Lambda_power": matrix_term["Lambda_power"] + state_term["Lambda_power"],
                        "coefficient": _iv(300 + source_cell, 301 + source_cell),
                        "input_matrix_coefficient": _clone(matrix_term["coefficient"]),
                        "input_state_coefficient": _clone(state_term["coefficient"]),
                        "coefficient_is_directed": True,
                        "ordinary_Z_order": 1,
                        "origin": "B1[0]*W1_Z2 (missing Z1 branch)",
                        "missing_product_rule": "B1[0]*W1_Z2",
                        "matrix_name": "B1",
                        "matrix_order": "0",
                        "matrix_key": matrix_key,
                        "row": 5,
                        "column": column,
                        "state_field": state_field,
                        "state_order": "2",
                        "target_field": target_field,
                        "matrix_term_index": 0,
                        "state_term_index": state_term_index,
                        "source_cell": source_cell,
                        "source_cell_index": source_cell,
                        "source_t": _clone(t),
                        "target_x": _clone(target_x),
                        "target_t": _clone(target_t),
                        "kernel_mass": _clone(masses[target_field]),
                        "kernel_mass_field": target_field,
                        "kernel_mass_exact": True,
                        "matrix_input_id": matrix_term["term_id"],
                        "state_input_id": state_term["term_id"],
                        "state_profile_term_id": state_term["term_id"],
                        "matrix_receipt_input_path": _clone(path),
                        "state_receipt_input_path": ["profile", "target_cells", source_cell, "W1_Z2", state_field, state_term_index],
                        "matrix_source_cell": source_cell,
                        "matrix_source_cell_index": source_cell,
                        "state_source_cell": state_term["source_cell"],
                        "state_source_cell_index": state_term["source_cell_index"],
                        "error_classes": ["fixture_v3:seam:error", *matrix_term["error_classes"], *state_term["error_classes"]],
                        "source_error_classes": list(matrix_term["error_classes"]),
                        "state_error_classes": list(state_term["error_classes"]),
                        "uncertainty": True,
                        "signed": True,
                        "signed_field": True,
                        "formal_atoms_materialized": False,
                        "F0_materialized": False,
                        "Lambda_materialized": False,
                        "no_cancellation": True,
                        "cancellation_performed": False,
                    }
                    seam_terms[target_field].append(product)

        seam_source_records.append(
            {
                "index": source_cell,
                "source_cell": source_cell,
                "source_cell_index": source_cell,
                "t": _clone(t),
                "columns": [1, 2, 3],
            }
        )
        old_source_records.append(
            {
                "index": source_cell,
                "source_cell": source_cell,
                "source_cell_index": source_cell,
                "t": _clone(t),
                "unknown_ids": [item["unknown_id"] for item in old_unknowns if item["source_cell"] == source_cell],
            }
        )

    old_result = {
        "schema": OLD_SCHEMA,
        "source_cells": old_source_records,
        "target_seam": {"x": _clone(target_x), "t": _clone(target_t), "rho": _iv(4, 4), "fixed": True},
        "kernels": old_kernels,
        "unknowns": old_unknowns,
        "unknown_ids_by_kind": {"missing_B1_W1_Z2": [item["unknown_id"] for item in old_unknowns]},
        "missing_product_rule": "B1[0]*W[2]",
        "ordinary_Z_orders": [0, 1],
        "raw_matrix_columns": {"C0": [1, 2, 3, 5], "B1": [1, 2, 3]},
        "error_classes_retained": True,
        "signed_coefficients_retained": True,
        "finite_K0_condition_package": True,
        "complete_secondword": False,
        "physical_field_claimed": False,
        "no_cancellation": True,
        "cancellation_performed": False,
    }
    profile_result = {
        "schema": PROFILE_SCHEMA,
        "target_cells": profile_cells,
        "finite_K0_condition_package": True,
        "complete_firstword_Z2": False,
        "physical_field_claimed": False,
        "ordinary_Z_orders": [2],
        "principal_kernel": "K0 only; no Kchi or higher kernel words",
    }
    seam_result = {
        "schema": SEAM_SCHEMA,
        "input_schema": PROFILE_SCHEMA,
        "source_cells": seam_source_records,
        "target_seam": {"x": _clone(target_x), "t": _clone(target_t), "rho": _iv(4, 4), "fixed": True},
        "kernels": seam_kernels,
        "terms_by_order": {"Z1": seam_terms},
        "missing_product_rule": "B1[0]*W1_Z2",
        "missing_branch_included": True,
        "known_Z1_records_untouched": True,
        "old_records_merged": False,
        "ordinary_Z_orders": [1],
        "state_Z_order": 2,
        "kernel_mass_exact": True,
        "F_and_Fx_exact": True,
        "finite_K0_condition_package": True,
        "principal_kernel": "K0 only; missing seam branch only",
        "complete_secondword": False,
        "physical_field_claimed": False,
        "m15_computed": False,
    }
    return old_result, profile_result, seam_result


def synthetic_low_output() -> dict[str, Any]:
    """Return the saved LOW-shaped matrix packet used for source validation."""

    cells: list[dict[str, Any]] = []
    for source_cell in range(CELL_COUNT):
        cells.append(
            {
                "index": source_cell,
                "t": _cell_t(source_cell),
                "rho": _iv(1 + source_cell, 2 + source_cell),
                "B1": {
                    "0": {
                        f"5,{column}": [_matrix_term(source_cell, f"5,{column}")]
                        for column in (1, 2, 3)
                    }
                },
            }
        )
    return {
        "schema": "n1_core_low_word_v3",
        "matrix_cells": cells,
        "finite_K0_condition_package": True,
        "complete_secondword": False,
    }


def _rejected(call: Any) -> bool:
    try:
        call()
    except (TypeError, ValueError, KeyError, IndexError):
        return True
    return False


def run_fixture(build: Any, c: Any = None) -> dict[str, Any]:
    """Run success plus path-mismatch, duplicate, and omission checks.

    ``c`` is accepted solely for compatibility with the other independent
    fixtures; it is intentionally unused because this bridge performs no
    interval arithmetic.
    """

    del c
    build_fn = getattr(build, "build", build)
    if not callable(build_fn):
        raise TypeError("build must be callable or expose build")
    old, profile, seam = synthetic_inputs()
    low_output = synthetic_low_output()
    old_before = _clone(old)
    try:
        result = build_fn(old, profile, seam, low_output)
    except Exception as exc:
        return {
            "schema": SCHEMA,
            "passed": False,
            "checks": {"build_completed": False},
            "check_count": 1,
            "error": f"{type(exc).__name__}: {exc}",
            "traceback": traceback.format_exc(),
            "synthetic_only": True,
            "implementation_imported": False,
            "numeric_rerun": False,
        }

    checks: dict[str, bool] = {}
    checks["build_completed"] = isinstance(result, Mapping)
    checks["schema"] = result.get("schema") == "core_z2_seam_bridge_v3"
    checks["twenty_four_old_unknowns"] = result.get("old_unknown_count") == 24
    checks["ninety_six_products"] = result.get("product_count") == 96
    checks["two_fields_per_unknown"] = all(
        set(item.get("product_ids_by_field", {})) == set(SEAM_FIELDS)
        and all(len(item["product_ids_by_field"][field]) == 2 for field in SEAM_FIELDS)
        for item in result.get("bindings", [])
    )
    checks["finite_conditional_scope"] = (
        result.get("finite_K0_condition_package") is True
        and result.get("physical_owner_claimed") is False
        and result.get("physical_field_claimed") is False
        and result.get("complete_secondword") is False
    )
    checks["exact_checks_reported"] = (
        result.get("source_t_exact_checked") is True
        and result.get("target_scales_exact_checked") is True
        and result.get("F_and_Fx_kernel_masses_exact_checked") is True
        and result.get("midpoint_used") is False
    )
    checks["raw_retention_reported"] = result.get("old_error_classes_retained") is True and all(
        binding.get("raw_directed_coefficients_retained") is True and binding.get("raw_ids_retained") is True
        for binding in result.get("bindings", [])
    )
    checks["old_input_untouched"] = old == old_before

    mismatch = _clone(seam)
    mismatch["terms_by_order"]["Z1"]["F1"][0]["matrix_receipt_input_path"] = ["wrong", "LOW", 0]
    checks["source_path_mismatch_rejected"] = _rejected(lambda: build_fn(old, profile, mismatch, low_output))

    duplicate = _clone(seam)
    duplicate["terms_by_order"]["Z1"]["F1"].append(_clone(duplicate["terms_by_order"]["Z1"]["F1"][0]))
    checks["duplicate_rejected"] = _rejected(lambda: build_fn(old, profile, duplicate, low_output))

    omission = _clone(seam)
    omission["terms_by_order"]["Z1"]["d_x_F1"].pop()
    checks["omission_rejected"] = _rejected(lambda: build_fn(old, profile, omission, low_output))

    wrong_state = _clone(seam)
    wrong_state["terms_by_order"]["Z1"]["F1"][0]["state_profile_term_id"] = "fixture_v3:missing-state-term"
    checks["missing_profile_state_rejected"] = _rejected(lambda: build_fn(old, profile, wrong_state, low_output))

    wrong_inner_source = _clone(seam)
    wrong_inner_source["terms_by_order"]["Z1"]["F1"][0]["state_source_cell"] = 7
    checks["wrong_inner_source_rejected"] = _rejected(lambda: build_fn(old, profile, wrong_inner_source, low_output))

    wrong_power = _clone(seam)
    wrong_power["terms_by_order"]["Z1"]["F1"][0]["F0_power"] += 1
    checks["formal_power_mismatch_rejected"] = _rejected(lambda: build_fn(old, profile, wrong_power, low_output))

    wrong_input_coefficient = _clone(seam)
    wrong_input_coefficient["terms_by_order"]["Z1"]["F1"][0]["input_matrix_coefficient"] = _iv(999, 1000)
    checks["input_coefficient_mismatch_rejected"] = _rejected(lambda: build_fn(old, profile, wrong_input_coefficient, low_output))

    return {
        "schema": SCHEMA,
        "passed": all(checks.values()),
        "checks": checks,
        "check_count": len(checks),
        "synthetic_only": True,
        "implementation_imported": False,
        "numeric_rerun": False,
        "ancestor_called": False,
        "path_identity_checked": True,
        "duplicate_and_omission_checked": True,
    }


def run(build: Any, c: Any = None) -> dict[str, Any]:
    return run_fixture(build, c)

