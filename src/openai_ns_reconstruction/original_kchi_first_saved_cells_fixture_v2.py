"""Synthetic positive and negative checks for the saved Kchi-cell binder.

The fixture uses only small structural records.  It does not import a
producer, evaluate an interval, or represent an original physical source.
"""

from __future__ import annotations

from copy import deepcopy
from typing import Any


def _endpoint(sign: int = 0, mantissa: int = 1, exponent: int = 0, count: int = 1) -> list[int]:
    return [sign, mantissa, exponent, count]


def _interval(*, symmetric: bool = False) -> dict[str, Any]:
    if symmetric:
        return {
            "lower_exact_mpf_tuple": _endpoint(1, 2, 0, 2),
            "upper_exact_mpf_tuple": _endpoint(0, 2, 0, 2),
            "directed": True,
        }
    return {
        "lower_exact_mpf_tuple": _endpoint(0, 2, 0, 2),
        "upper_exact_mpf_tuple": _endpoint(0, 4, 0, 3),
        "directed": True,
    }


def _legacy_log_interval() -> dict[str, Any]:
    """Model the historical logs that predate the directed marker."""

    return {
        "lower": _endpoint(0, 2, 0, 2),
        "upper": _endpoint(0, 4, 0, 3),
        "lower_exact_mpf_tuple": _endpoint(0, 2, 0, 2),
        "upper_exact_mpf_tuple": _endpoint(0, 4, 0, 3),
    }


def _finite_term(
    field: str,
    order: str,
    source: int,
    target: int,
    chi_order: int,
    source_order: int,
    leibniz_multiplier: int,
    term_index: int,
) -> dict[str, Any]:
    source_id = "synthetic:v2:g7:{}:{}:{}:{}:{}".format(
        source, field, order, chi_order, source_order, term_index
    )
    cell_id = "synthetic-v2-cell:{}".format(source)
    term = {
        "coefficient": _interval(),
        "powers": [1, 2, 0, -2],
        "hb_power": 1,
        "F0_power": 2,
        "Pstar_squared_power": 0,
        "Lambda_power": -2,
        "origin": {
            "operation": "kchi_first_correction_product",
            "source_ids": [source_id],
            "error_classes": ["synthetic_v2_source_interval"],
            "uncertainty": True,
            "target_cell": target,
            "source_cell": source,
            "chi_order": chi_order,
            "source_order": source_order,
            "leibniz_multiplier": leibniz_multiplier,
            "kernel_order": 1,
            "kernel": "(-2*chi) T K0",
            "source_cell_id": cell_id,
            "chi_term_index": 0,
            "source_term_index": term_index,
        },
        "source_ids": [source_id],
        "error_classes": ["synthetic_v2_source_interval"],
        "uncertainty": True,
        "target_cell": target,
        "source_cell": source,
        "source_cell_id": cell_id,
        "source_term_id": source_id,
        "source_index": source,
        "source_term_index": term_index,
        "chi_order": chi_order,
        "source_order": source_order,
        "leibniz_multiplier": leibniz_multiplier,
        "principal_kernel": "(-2*chi) T K0",
        "kernel_mass": _interval(),
        "kernel_mass_bare": _interval(),
        "correction_factor": _interval(),
        "target_x": _interval(),
        "source_x": _interval(),
        "kernel_mass_field": field,
    }
    return term


def _source_packet(source: int, target: int) -> dict[str, Any]:
    # This is the producer's field-first orientation:
    # {F1: {Z0: [...], Z1: [...]}, d_x_F1: {...}}.
    pairs = {
        "Z0": ((0, 0, 1),),
        "Z1": ((1, 0, 1), (0, 1, 1)),
        "Z2": ((2, 0, 1), (1, 1, 2), (0, 2, 1)),
    }
    correction: dict[str, dict[str, list[dict[str, Any]]]] = {
        field: {
            order: [
                _finite_term(
                    field,
                    order,
                    source,
                    target,
                    chi_order,
                    source_order,
                    multiplier,
                    term_index,
                )
                for term_index, (chi_order, source_order, multiplier) in enumerate(pairs[order])
            ]
            for order in ("Z0", "Z1", "Z2")
        }
        for field in ("F1", "d_x_F1")
    }
    return {
        "source_index": source,
        "target_index": target,
        "source_x": _interval(),
        "target_x": _interval(),
        "triangular": source == target,
        "correction": correction,
    }


def _tail_term(field: str, order: int, source_order: int, kernel_z_order: int) -> dict[str, Any]:
    source_id = "synthetic:v2:tail:{}:{}:{}:{}".format(
        field, order, source_order, kernel_z_order
    )
    return {
        "coefficient": _interval(symmetric=True),
        "powers": [1, 2, 0, -2],
        "source_ids": [source_id],
        "error_classes": ["Kchi_j_ge_2_kernel_tail", "synthetic_v2_tail"],
        "uncertainty": True,
        "source_order": source_order,
        "kernel_Z_order": kernel_z_order,
        "physical_source_width_bound": 2,
    }


def _synthetic() -> dict[str, Any]:
    family = {
        "actual_five_defect_family_sha256": "d" * 64,
        "implicit_source_sha256": "e" * 64,
        "datum_enclosure_sha256": "f" * 64,
    }
    logs = {
        "hb": _legacy_log_interval(),
        "F0_anchor": _legacy_log_interval(),
        "Pstar_squared": _legacy_log_interval(),
        "Lambda": _legacy_log_interval(),
    }
    cells = []
    for target in range(8):
        cells.append(
            {
                "index": target,
                "target_index": target,
                "x": _interval(),
                "target_x": _interval(),
                "source_indices": list(range(target + 1)),
                "source_j_le_target_i": True,
                "triangular_diagonal_source": target,
                "axis_endpoint_continuation": target == 0,
                "W0_correction_by_source_cell": {
                    str(source): _source_packet(source, target)
                    for source in range(target + 1)
                },
            }
        )

    tail = {
        field: {
            str(order): [
                _tail_term(field, order, source_order, order - source_order)
                for source_order in range(order + 1)
            ]
            for order in range(3)
        }
        for field in ("F1", "d_x_F1")
    }
    bounds = {
        "F1": [_interval(), _interval(), _interval()],
        "d_x_F1": [_interval(), _interval(), _interval()],
        "ordinary_Z_orders": [0, 1, 2],
        "scope": "pointwise potential-kernel remainder j>=2 only",
        "source_or_feedback_included": False,
        "first_omitted_index": 2,
        "bounds_are_not_field_values": True,
    }
    finite = {
        "source_family": family,
        "source_atoms": logs,
        "source_power_basis": ["hb", "F0_anchor", "Pstar_squared", "Lambda"],
        "ordinary_Z_orders": [0, 1, 2],
        "ordinary_product_rule": {
            "Z0": [[0, 0, 1]],
            "Z1": [[1, 0, 1], [0, 1, 1]],
            "Z2": [[2, 0, 1], [1, 1, 2], [0, 2, 1]],
        },
        "ordinary_Z0_product_rule": [[0, 0, 1]],
        "ordinary_Z1_product_rule": [[1, 0, 1], [0, 1, 1]],
        "ordinary_Z2_product_rule": [[2, 0, 1], [1, 1, 2], [0, 2, 1]],
        "cell_count": 8,
        "target_cells": cells,
        "signed_first_correction": "(-2*chi) T K0",
        "Kchi_first_correction_only": True,
        "Kchi_higher_words_included": False,
        "allword_tail_added": False,
        "source_to_W0_only": True,
        "source_rows_preserved_through_Z2": True,
        "source_ids_retained": True,
        "source_indices_retained": True,
        "signed_coefficients_retained": True,
        "formal_powers_retained": True,
        "source_terms_recomputed": False,
        "old_words_recomputed": False,
    }
    input_index = {}
    for order in range(3):
        input_index["synthetic:v2:g7:0:F1:Z{}:0".format(order)] = {
            "receipt": "synthetic-v2-g7-receipt",
            "path": ["global_low_source", str(order), "5", 0],
        }
        input_index["chi:Z{}".format(order)] = {
            "inputs": ["synthetic-v2-chi-source"],
            "definition": "synthetic v2 chi jet",
            "ordinary_Z_order": order,
        }
    return {
        "source_family": family,
        "source_logs": logs,
        "domain": "x in [0,2], rho in [0,4], real Z in [-1,1]",
        "input_index": input_index,
        "first_correction": finite,
        "kernel_tail_bounds": bounds,
        "W0_separate_symmetric_kernel_tail": tail,
        "first_correction_signed": True,
        "Kchi_tail_is_uncertainty_not_signed_field": True,
        "tail_includes_feedback": False,
        "feedback_cross_terms_added": False,
        "appended_to_cumulative_m15": False,
        "complete_W0_computed": False,
        "complete_m15_computed": False,
        "terminal_closure_proved": False,
        "source_cells_reuse_one_whole_core_enclosure": True,
        "source_boxes_not_point_values": True,
    }


def _rejects(leaf: Any, value: Any) -> bool:
    try:
        leaf.bind(value)
    except (TypeError, ValueError, KeyError):
        return True
    return False


def run_fixture(leaf: Any) -> dict[str, Any]:
    """Return independent positive and negative structural checks."""

    data = _synthetic()
    result = leaf.bind(data)
    finite_term = result["cells"][4]["source_cells"][2]["ordinary_Z"]["Z2"]["F1"][0]
    tail_term = result["symmetric_kernel_tail"]["d_x_F1"]["Z2"][0]
    checks = {
        "schema": result["schema"] == "original_kchi_first_saved_cells_v2",
        "eight_target_cells": len(result["cells"]) == 8,
        "source_target_cell_ids_preserved": finite_term["source_cell"] == 2 and finite_term["target_cell"] == 4,
        "four_powers_preserved": finite_term["formal_source_powers"] == [1, 2, 0, -2],
        "binary_endpoints_preserved": finite_term["coefficient"] == _interval(),
        "z2_preserved": "Z2" in result["cells"][7]["source_cells"][7]["ordinary_Z"],
        "finite_term_is_signed": finite_term["term_role"] == "finite_j1_signed" and finite_term["error_only"] is False,
        "tail_is_separate_error": tail_term["error_only"] is True and tail_term["symmetric_error"] is True and tail_term["signed_field"] is False,
        "tail_endpoints_symmetric": tail_term["coefficient"] == _interval(symmetric=True),
        "tail_order_decomposition": tail_term["source_order"] + tail_term["kernel_Z_order"] == 2,
        "family_and_logs_preserved": result["source_family"] == data["source_family"] and result["source_logs"] == data["source_logs"],
        "legacy_logs_without_directed_preserved": all(
            "directed" not in result["source_logs"][name]
            and "lower" in result["source_logs"][name]
            and "upper" in result["source_logs"][name]
            for name in ("hb", "F0_anchor", "Pstar_squared", "Lambda")
        ),
        "no_sum_or_ledger": result["terms_were_not_summed"] is True and result["ledger_was_not_extended"] is True,
        "no_full_claim": result["complete_response"] is False and result["physical_m15_computed"] is False,
    }

    bad = deepcopy(data)
    del bad["first_correction"]["target_cells"][0]["W0_correction_by_source_cell"]["0"]["correction"]["F1"]["Z2"]
    checks["missing_z2_rejected"] = _rejects(leaf, bad)

    bad = deepcopy(data)
    del bad["first_correction"]["target_cells"][1]["W0_correction_by_source_cell"]["0"]["correction"]["F1"]["Z0"][0]["source_term_id"]
    checks["missing_term_id_rejected"] = _rejects(leaf, bad)

    bad = deepcopy(data)
    bad["first_correction"]["target_cells"][0]["W0_correction_by_source_cell"]["0"]["correction"]["F1"]["Z0"][0]["chi_order"] = 1
    checks["wrong_leibniz_pair_rejected"] = _rejects(leaf, bad)

    bad = deepcopy(data)
    bad["first_correction"]["target_cells"][0]["W0_correction_by_source_cell"]["0"]["correction"]["F1"]["Z0"][0]["error_only"] = True
    checks["finite_error_marker_rejected"] = _rejects(leaf, bad)

    bad = deepcopy(data)
    bad["W0_separate_symmetric_kernel_tail"]["F1"]["1"][0]["coefficient"] = _interval()
    checks["nonsymmetric_tail_rejected"] = _rejects(leaf, bad)

    bad = deepcopy(data)
    bad["appended_to_cumulative_m15"] = True
    checks["ledger_claim_rejected"] = _rejects(leaf, bad)

    bad = deepcopy(data)
    del bad["source_logs"]
    checks["missing_source_logs_rejected"] = _rejects(leaf, bad)

    bad = deepcopy(data)
    bad["source_logs"]["hb"]["directed"] = False
    checks["legacy_directed_false_rejected"] = _rejects(leaf, bad)

    return {
        "schema": "original_kchi_first_saved_cells_fixture_v2",
        "passed": all(checks.values()),
        "checks": checks,
        "check_count": len(checks),
        "synthetic_only": True,
        "source_values_computed": False,
        "method": "structural copy with independent missing-field and tail-role rejections",
    }


run = run_fixture


__all__ = ["run_fixture", "run"]

