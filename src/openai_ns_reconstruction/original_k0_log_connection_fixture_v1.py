"""Fresh synthetic records for ``original_k0_log_connection_v1``.

The fixture is source-only.  It creates decoded JSON-shaped records in
memory, calls a supplied ``bind`` API, and reports structural observations.
It does not import a producer, read a receipt, decode gzip, or evaluate an
interval.
"""

from __future__ import annotations

from copy import deepcopy
from typing import Any


SCHEMA = "original_k0_log_connection_fixture_v1"
_BASIS = ["hb", "F0_anchor", "Pstar_squared", "Lambda", "F0_local", "L_local"]
_HASH = "d" * 64


def _interval(seed: int) -> dict[str, Any]:
    lower, upper = seed + 1, seed + 2
    return {
        "lower_exact_mpf_tuple": [0, lower, -8, lower.bit_length()],
        "upper_exact_mpf_tuple": [0, upper, -8, upper.bit_length()],
        "directed": True,
    }


def _family() -> dict[str, str]:
    return {
        "actual_five_defect_family_sha256": "a" * 64,
        "implicit_source_sha256": "b" * 64,
        "datum_enclosure_sha256": "c" * 64,
    }


def _logs() -> dict[str, str]:
    return {
        "hb": "synthetic hb source log",
        "F0_anchor": "synthetic F0 source log",
        "Pstar_squared": "synthetic pressure source log",
        "Lambda": "synthetic Lambda source log",
    }


def _original_term(index: int, term_name: str, powers: list[int]) -> dict[str, Any]:
    coefficient = _interval(100 + index * 10 + (1 if term_name == "local" else 2))
    return {
        "coefficient": deepcopy(coefficient),
        "F0_power": powers[1],
        "Lambda_power": powers[3],
        "source_cell": index,
        "target_cell": index,
        "matrix_sector_identity": "fixture-matrix-" + term_name + "-" + str(index),
        "W0_sector_identity": "fixture-w0-" + term_name + "-" + str(index),
        "source_id": "fixture-source-" + term_name + "-" + str(index),
        "error_metadata": {
            "uncertainty": True,
            "error_classes": ["fixture_saved_error", "inherited_interval_error"],
        },
    }


def _term(index: int, term_name: str, powers: list[int], *, local_factor: bool) -> dict[str, Any]:
    input_powers = list(powers)
    if local_factor:
        input_powers[3] += 2
    original = _original_term(index, term_name, input_powers)
    coefficient = deepcopy(original["coefficient"])
    product_factors = ["rho", "Phi_Z0"]
    if local_factor:
        product_factors.append("L_local")
    product_factors.append("F1_0")
    result = {
        "channel": "swirl",
        "ordinary_Z": "Z1",
        "piece": "L_local*Phi0*F1_0" if local_factor else "fixture_nonlocal",
        "coefficient": coefficient,
        "coefficient_is_unmodified_saved_binary_record": True,
        "original_term": original,
        "saved_formal_source_powers": input_powers[:4],
        "formal_source_powers": list(powers),
        "signed_factor": -2,
        "source_cell": index,
        "target_cell": index,
        "matrix_sector_identity": original["matrix_sector_identity"],
        "W0_sector_identity": original["W0_sector_identity"],
        "source_ids": [original["source_id"], "fixture-background-Phi-Z0"],
        "source_cell_provenance": {
            "source_cell": index,
            "target_cell": index,
            "matrix_sector_identity": original["matrix_sector_identity"],
            "W0_sector_identity": original["W0_sector_identity"],
        },
        "uncertainty": True,
        "error_classes": ["fixture_saved_error", "inherited_interval_error"],
        "error_metadata": deepcopy(original["error_metadata"]),
        "source_identity_status": "source identity unverified",
        "product_factors": product_factors,
        "unknown": False,
    }
    if local_factor:
        result["local_factor"] = {"symbol": "L_local", "formal_power": 1}
    return result


def _symbolic() -> dict[str, Any]:
    cells = []
    for index in range(8):
        affected = _term(index, "local", [1, 2, 3, -7, 1, 1], local_factor=True)
        unaffected = _term(index, "ordinary", [1, 2, 3, -7, 1, 0], local_factor=False)
        axial = _term(index, "axial", [1, 1, 3, -6, 0, 0], local_factor=False)
        cells.append(
            {
                "index": index,
                "target_cell": index,
                "source_indices": list(range(index + 1)),
                "triangular_diagonal_source": index,
                "axis_endpoint_continuation": index == 0,
                "t": _interval(index),
                "x": _interval(index + 20),
                "rho": _interval(index + 40),
                "terms": {
                    "axial": {"Z0": [deepcopy(axial)], "Z1": [deepcopy(axial)]},
                    "swirl": {
                        "Z0": [deepcopy(unaffected)],
                        "Z1": [affected, unaffected],
                    },
                },
            }
        )
    return {
        "schema": "original_k0_symbolic_m15_cells_v1",
        "status": "symbolic_unintegrated",
        "source_family": _family(),
        "source_logs": _logs(),
        "source_kernel_scope": "synthetic K0-only W1 symbolic products",
        "source_power_basis": list(_BASIS),
        "ordinary_Z_orders": ["Z0", "Z1"],
        "channels": ["axial", "swirl"],
        "cells": cells,
        "unknown_terms": [
            {"name": "N1p", "status": "unknown", "reason": "pressure outside this leaf"},
            {"name": "Kchi", "status": "unknown", "reason": "Kchi correction is not supplied"},
            {"name": "allword_tail", "status": "unknown", "reason": "tail is not supplied"},
        ],
        "unknown_terms_are_not_zero": True,
        "source_identity_status": "source identity unverified",
        "source_identity_verified": False,
        "F0_scale_identity": {
            "F0_anchor": "independent formal atom",
            "F0_local": "independent formal atom",
            "identified": False,
        },
        "source_values_computed": False,
        "point_value_selected": False,
        "interval_arithmetic_performed": False,
        "integration_performed": False,
        "numeric_operator_called": False,
        "N1p_included": False,
        "Kchi_included": False,
        "allword_tail_included": False,
        "appended_to_212_ledger": False,
        "complete_W1": False,
        "complete_m15_computed": False,
    }


def _connection() -> dict[str, Any]:
    return {
        "source_family": _family(),
        "source_logs": _logs(),
        "ell_over_lambda": _interval(900),
        "saved_core_Z_domain": {"lower_exact_mpf_tuple": [1, 1, 0, 1],
                                "upper_exact_mpf_tuple": [0, 1, 0, 1]},
        "authenticated_raw_receipt_hash": _HASH,
        "provenance_correction": {
            "declared_by_runner": True,
            "receipt_sha256": _HASH,
            "correction_id": "fixture-provenance-correction",
        },
        "same_owner_identity_declared": False,
        "analytic_recipe_declared": False,
    }


def _call(leaf: Any, symbolic: Any, connection: Any) -> Any:
    function = getattr(leaf, "bind", leaf)
    if not callable(function):
        raise TypeError("leaf must expose bind")
    return function(symbolic, connection)


def _reject(leaf: Any, symbolic: Any, connection: Any) -> bool:
    try:
        _call(leaf, symbolic, connection)
    except (TypeError, ValueError, KeyError):
        return True
    return False


def run_fixture(leaf: Any) -> dict[str, Any]:
    """Apply structural checks to a caller-supplied bridge."""

    symbolic = _symbolic()
    connection = _connection()
    before_symbolic = deepcopy(symbolic)
    before_connection = deepcopy(connection)
    result = _call(leaf, symbolic, connection)
    affected = [term for cell in result["cells"] for term in cell["substituted_terms"]]
    expected_ratio = connection["ell_over_lambda"]
    checks: dict[str, bool] = {
        "schema": result.get("schema") == "original_k0_log_connection_v1",
        "eight_cells": len(result.get("cells", ())) == 8,
        "one_affected_term_per_cell": len(affected) == 8 and all(cell["substituted_term_count"] == 1 for cell in result["cells"]),
        "six_axis_transform": all(term["formal_source_powers"] == [1, 2, 3, -6, 1, 0] for term in affected),
        "original_six_powers_preserved": all(term["original_six_powers"] == [1, 2, 3, -7, 1, 1] for term in affected),
        "ratio_copied_byte_exact": all(term["saved_ell_over_lambda"] == expected_ratio for term in affected)
        and result["saved_ell_over_lambda"] == expected_ratio,
        "coefficient_endpoints_unchanged": all(
            term["coefficient"] == symbolic["cells"][index]["terms"]["swirl"]["Z1"][0]["coefficient"]
            for index, term in enumerate(affected)
        ),
        "source_target_matrix_w0_ids_retained": all(
            term["source_cell"] == symbolic["cells"][index]["terms"]["swirl"]["Z1"][0]["source_cell"]
            and term["target_cell"] == symbolic["cells"][index]["terms"]["swirl"]["Z1"][0]["target_cell"]
            and term["matrix_sector_identity"]
            == symbolic["cells"][index]["terms"]["swirl"]["Z1"][0]["matrix_sector_identity"]
            and term["W0_sector_identity"]
            == symbolic["cells"][index]["terms"]["swirl"]["Z1"][0]["W0_sector_identity"]
            for index, term in enumerate(affected)
        ),
        "error_metadata_retained": all(term.get("uncertainty") is True and term.get("error_classes") for term in affected),
        "no_double_L": all(
            term["formal_source_powers"][5] == 0
            and term["product_factors"].count("L_local") == 0
            and term["product_factors"].count("saved_ell_over_lambda") == 1
            for term in affected
        ),
        "scales_remain_distinct": result["F0_scale_identity"]["identified"] is False,
        "unknowns_preserved": result["unknown_terms"] == symbolic["unknown_terms"]
        and result["unknown_terms_are_not_zero"] is True
        and result["Kchi_included"] is False
        and result["N1p_included"] is False
        and result["allword_tail_included"] is False,
        "conditional_not_function_identity": result["conditional_identity"]["status"] == "conditional"
        and result["conditional_identity"]["callable_identity_proved"] is False
        and result["conditional_identity"]["endpoint_identity_proved"] is False,
        "inputs_unmodified": symbolic == before_symbolic and connection == before_connection,
        "no_ledger_append": result["appended_to_212_ledger"] is False,
    }

    bad_family = deepcopy(connection)
    bad_family["source_family"]["implicit_source_sha256"] = "e" * 64
    checks["source_family_mismatch_rejected"] = _reject(leaf, symbolic, bad_family)

    bad_logs = deepcopy(connection)
    bad_logs["source_logs"]["Lambda"] = "different source log"
    checks["source_logs_mismatch_rejected"] = _reject(leaf, symbolic, bad_logs)

    missing_ratio = deepcopy(connection)
    del missing_ratio["ell_over_lambda"]
    checks["missing_ratio_rejected"] = _reject(leaf, symbolic, missing_ratio)

    malformed_ratio = deepcopy(connection)
    del malformed_ratio["ell_over_lambda"]["lower_exact_mpf_tuple"]
    checks["missing_binary_tuple_rejected"] = _reject(leaf, symbolic, malformed_ratio)

    double_l = deepcopy(symbolic)
    double_l["cells"][0]["terms"]["swirl"]["Z1"][0]["formal_source_powers"][5] = 2
    checks["double_L_rejected"] = _reject(leaf, double_l, connection)

    wrong_basis = deepcopy(symbolic)
    wrong_basis["source_power_basis"][5] = "unknown_scale"
    checks["wrong_basis_rejected"] = _reject(leaf, wrong_basis, connection)

    wrong_branch = deepcopy(symbolic)
    wrong_branch["cells"][0]["terms"]["swirl"]["Z1"][0]["piece"] = "unrelated_branch"
    checks["wrong_derivative_branch_rejected"] = _reject(leaf, wrong_branch, connection)

    return {
        "schema": SCHEMA,
        "passed": all(checks.values()),
        "checks": checks,
        "check_count": len(checks),
        "synthetic_only": True,
        "source_values_computed": False,
        "interval_arithmetic_performed": False,
        "producer_imported": False,
        "old_source_called": False,
    }


__all__ = ["SCHEMA", "run_fixture"]
