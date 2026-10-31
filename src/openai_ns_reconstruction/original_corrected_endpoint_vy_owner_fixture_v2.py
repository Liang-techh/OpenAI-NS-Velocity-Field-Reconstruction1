"""Independent in-memory fixture for the corrected endpoint ``V_y`` owner.

This fixture owns synthetic reports and source rows.  It does not import the
owner, read a receipt/report file, evaluate a graph, or use a numerical
midpoint.  A caller supplies the consumer callable to :func:`run_fixture`.
"""

from __future__ import annotations

from collections.abc import Mapping
from copy import deepcopy
from fractions import Fraction
from typing import Any, Callable


SCHEMA = "original_corrected_endpoint_vy_owner_fixture_v2"
SOURCE_FAMILY = "synthetic_original_preRc_shared_family_v2_fresh"
CELLS = (("-1", "-.5"), ("-.5", "0"), ("0", ".5"), (".5", "1"))

SOURCE_BINDINGS = {
    "lower": "physical.source_rows(first_micro, ('sc_half',))",
    "lower_V_y": "actual_background_source.V_y[0:4]",
    "lower_U_ZZ": "factorial source rows V[2:4]",
    "lower_U_yZ": "factorial source rows V_y[1:4]",
    "common_Pstar": "outer.owner.owner(ends).Pstar",
}
RECEIPT_CHECKS = {
    "exact_four_Z_cell_cover": True,
    "lower_actual_V_y_rows_and_U_ZZ_checked": True,
    "physical_common_Pstar_once": True,
}


def _row(cell: tuple[str, str], order: int, coefficient: int) -> dict[str, Any]:
    return {
        "formal_positive_scale": {
            "basis": "fixture_v2_source_scale",
            "powers": [1, -2, 0, order],
            "expression": "B*R^3*Pstar^-1*Z^order",
        },
        "coefficient_interval": [coefficient - 1, coefficient + 1],
        "log_absolute_upper": "fixture-v2-cap-not-a-value",
        "encloses_original_source_function": True,
        "point_value_selected": False,
        "ordinary_Z_order": order,
        "source_family": SOURCE_FAMILY,
        "provenance": {
            "source_path": "fixture_v2.lower.background_V_y_rows",
            "cell": list(cell),
            "row_order": order,
            "same_source_id": "fixture-v2-lower-sc-half-source",
            "factorial_conversion": "already_once_v2" if order else "value_row_v2",
        },
    }


def _scale_row(name: str, cell: tuple[str, str]) -> dict[str, Any]:
    return {
        "formal_positive_scale": {"basis": f"fixture_v2_{name}", "expression": name},
        "coefficient_interval": [8, 10] if name == "R" else [12, 14],
        "log_absolute_upper": f"fixture-v2-{name}-cap-not-a-value",
        "encloses_original_source_function": True,
        "point_value_selected": False,
        "source_family": SOURCE_FAMILY,
        "provenance": {
            "source_path": f"fixture_v2.lower.source_rows.{name}",
            "cell": list(cell),
            "same_source_id": "fixture-v2-lower-sc-half-source",
        },
    }


def build_synthetic_receipt() -> dict[str, Any]:
    return {
        "all_passed": True,
        "source_family": SOURCE_FAMILY,
        "cells": [list(cell) for cell in CELLS],
        "checks": dict(RECEIPT_CHECKS),
    }


def build_synthetic_endpoint_report() -> dict[str, Any]:
    samples: dict[str, Any] = {}
    for cell_index, cell in enumerate(CELLS):
        coefficients = (11 + 2 * cell_index, -(13 + 2 * cell_index), 17 + 2 * cell_index, -(19 + 2 * cell_index))
        rows = [_row(cell, order, coefficient) for order, coefficient in enumerate(coefficients)]
        lower_proof = {
            "endpoint": "first_micro/sc_half singleton",
            "source_family": SOURCE_FAMILY,
            "actual_V_y_rows_0_through_3_factorial_converted_once": True,
            "actual_U_ZZ_from_V_rows_2_and_3": True,
            "actual_U_yZ_from_V_y_rows_1_through_3": True,
            "common_Pstar_applied_once": True,
            "same_source": True,
            "source_quantity_paths": {
                "V_y": "actual_background_source.V_y",
                "R": "original_generic_source.original_physical_radius",
                "Pstar": SOURCE_BINDINGS["common_Pstar"],
            },
        }
        lower_rows = {"V_y": rows, "R": _scale_row("R", cell), "Pstar": _scale_row("Pstar", cell)}
        samples[str(cell)] = {
            "source_family": SOURCE_FAMILY,
            "exact_Z_cell": list(cell),
            "lower": {"source_rows": lower_rows, "proof": lower_proof, "caps": {"not_consumed": True}},
            "upper": {"proof": {"source_family": SOURCE_FAMILY}},
            "proof": {
                "exact_Z_cell": list(cell),
                "source_family": SOURCE_FAMILY,
                "physical_common_Pstar_applied_once": True,
            },
        }
    return {
        "candidate_preRc_endpoint_C1_constructed": True,
        "source_family": SOURCE_FAMILY,
        "cells": [list(cell) for cell in CELLS],
        "samples": samples,
        "source_quantity_scope": "actual lower sc_half singleton and O3_power right endpoint",
        "source_bindings": dict(SOURCE_BINDINGS),
        "endpoint_receipt": build_synthetic_receipt(),
    }


def build_synthetic_corrected_reserved_report(endpoint_report: Mapping[str, Any]) -> dict[str, Any]:
    """Build the terminal-source wrapper used by the corrected report check."""

    selected = endpoint_report["samples"][str(("0", ".5"))]["lower"]
    return {"corrected_terminal_source_inputs": {
        "numeric_endpoint_cell": ["0", ".5"],
        "physical_endpoint_source_rows": deepcopy(selected["source_rows"]),
        "physical_endpoint_source_proof": deepcopy(selected["proof"]),
    }}


def _evaluate_formula_node(node: Mapping[str, Any], values: Mapping[str, Fraction]) -> Fraction:
    """Tiny independent evaluator for the emitted synthetic formula AST."""

    operation = node.get("operation")
    if operation == "literal":
        return Fraction(node["value"])
    if operation == "factor":
        name = node.get("name")
        if name not in values:
            raise ValueError("fixture evaluator encountered an unbound symbolic factor")
        return values[name]
    if operation == "product":
        result = Fraction(1)
        for factor in node.get("factors", ()):
            result *= _evaluate_formula_node(factor, values)
        return result
    if operation == "sum":
        return sum((_evaluate_formula_node(term, values) for term in node.get("terms", ())), Fraction(0))
    if operation == "power":
        return _evaluate_formula_node(node["base"], values) ** int(node["exponent"])
    if operation == "quotient":
        denominator = _evaluate_formula_node(node["denominator"], values)
        if denominator == 0:
            raise ValueError("fixture formula denominator is zero")
        return _evaluate_formula_node(node["numerator"], values) / denominator
    raise ValueError("fixture evaluator encountered an unsupported formula operation")


def independent_rational_formula_check(formula: Mapping[str, Any] | None = None) -> dict[str, Any]:
    """Evaluate the leaf F/Fprime AST against fresh rational data.

    At Z=1/3 and delta=1/4, the independent quotient-rule derivative is
    compared with the value obtained by recursively evaluating the emitted AST.
    """

    z = Fraction(1, 3)
    delta = Fraction(1, 4)
    L = 1 - delta * z * z
    F = Fraction(2) * z / L
    quotient_rule = Fraction(2, 1) / L + Fraction(4, 1) * delta * z * z / (L * L)
    single_denominator = (Fraction(2, 1) + Fraction(2, 1) * delta * z * z) / (L * L)
    emitted_F = emitted_Fprime = None
    emitted_ok = False
    if isinstance(formula, Mapping):
        definitions = formula.get("definitions")
        if isinstance(definitions, Mapping):
            values = {"Z": z, "delta": delta}
            emitted_F = _evaluate_formula_node(definitions["F"], values)
            emitted_Fprime = _evaluate_formula_node(definitions["Fprime"], values)
            emitted_ok = emitted_F == F and emitted_Fprime == quotient_rule
    return {
        "Z": [1, 3],
        "delta": [1, 4],
        "F": F,
        "emitted_F": emitted_F,
        "quotient_rule_Fprime": quotient_rule,
        "single_denominator_Fprime": single_denominator,
        "emitted_Fprime": emitted_Fprime,
        "passed": quotient_rule == single_denominator and emitted_ok,
        "emitted_formula_evaluated_only_on_synthetic_input": emitted_ok,
    }


def mutate_duplicate_cell(report: Mapping[str, Any]) -> dict[str, Any]:
    result = deepcopy(dict(report))
    result["cells"][1] = list(result["cells"][0])
    return result


def mutate_missing_cell(report: Mapping[str, Any]) -> dict[str, Any]:
    result = deepcopy(dict(report))
    result["cells"] = result["cells"][:-1]
    result["samples"].pop(str(CELLS[-1]))
    return result


def mutate_source_corruption(report: Mapping[str, Any]) -> dict[str, Any]:
    result = deepcopy(dict(report))
    result["samples"][str(CELLS[2])]["lower"]["source_rows"]["V_y"][1]["provenance"]["same_source_id"] = "fixture-v2-broken-source-identifier"
    return result


def _call(consumer: Callable[..., Any], report: Mapping[str, Any], corrected: Mapping[str, Any]) -> Any:
    return consumer(report, corrected_reserved_report=corrected, frozen=True, entered=True)


def run_fixture(consumer: Callable[..., Any]) -> dict[str, Any]:
    """Run structural acceptance and fail-closed mutation checks."""

    if not callable(consumer):
        raise TypeError("consumer must be a callable source-layer owner")
    report = build_synthetic_endpoint_report()
    corrected = build_synthetic_corrected_reserved_report(report)
    accepted = _call(consumer, report, corrected)
    checks: dict[str, bool] = {
        "accepted": isinstance(accepted, Mapping),
        "four_cells": isinstance(accepted, Mapping) and accepted.get("cells") == [list(cell) for cell in CELLS],
        "v2_nonzero_C0_C1_rows": False,
        "v2_formal_positive_scale_retained": False,
        "v2_whole_rows_and_provenance_retained": False,
        "v2_ordinary_Z1_formula_is_present": False,
        "v2_structured_repair_ast": False,
        "v2_independent_rational_Fprime": False,
        "v2_duplicate_rejected": False,
        "v2_missing_rejected": False,
        "v2_source_corruption_rejected": False,
        "v2_unattested_comparison_rejected": False,
    }
    if isinstance(accepted, Mapping):
        cells = accepted.get("cell_sources", {})
        try:
            entries = [cells[str(cell)] for cell in CELLS]
            rational = independent_rational_formula_check(entries[0]["formula"])
            checks["v2_independent_rational_Fprime"] = rational["passed"] is True
            checks["v2_nonzero_C0_C1_rows"] = all(
                interval[0] != 0 and interval[1] != 0
                for entry in entries for order in (0, 1)
                for interval in (entry["lower"]["source_rows"]["V_y_C0_C1"][order]
                                 ["row"]["coefficient_interval"],)
            )
            checks["v2_formal_positive_scale_retained"] = all(
                entry["lower"]["source_rows"]["V_y_C0_C1"][order]["symbolic_scaling"]["source_key"]
                == "formal_positive_scale"
                for entry in entries for order in (0, 1)
            )
            checks["v2_whole_rows_and_provenance_retained"] = all(
                entry["lower"]["source_rows"]["V_y_C0_C1"][order]["raw_row_preserved"]
                and entry["lower"]["source_rows"]["V_y_C0_C1"][order]["provenance"]["same_source_id"]
                == "fixture-v2-lower-sc-half-source"
                for entry in entries for order in (0, 1)
            )
            checks["v2_ordinary_Z1_formula_is_present"] = all(
                "Y1 is the ordinary derivative" in entry["formula"]["ordinary_Z1_rule"]
                for entry in entries
            )
            checks["v2_structured_repair_ast"] = all(
                entry["formula"]["repair0"]["operation"] == "product"
                and entry["formula"]["repair1"]["operation"] == "product"
                and entry["formula"]["definitions"]["F"]["operation"] == "quotient"
                and entry["formula"]["definitions"]["Fprime"]["operation"] == "quotient"
                and entry["formula"]["repair0"]["factors"][2]["path"].endswith(".formula.definitions.F")
                and entry["formula"]["repair1"]["factors"][2]["terms"][0]["factors"][0]["path"].endswith(".formula.definitions.F")
                and entry["formula"]["repair1"]["factors"][2]["terms"][1]["factors"][0]["path"].endswith(".formula.definitions.Fprime")
                for entry in entries
            )
        except (KeyError, TypeError):
            pass

    for name, mutated in (
        ("v2_duplicate_rejected", mutate_duplicate_cell(report)),
        ("v2_missing_rejected", mutate_missing_cell(report)),
        ("v2_source_corruption_rejected", mutate_source_corruption(report)),
    ):
        try:
            _call(consumer, mutated, corrected)
        except Exception:
            checks[name] = True

    try:
        consumer(report, corrected_reserved_report=corrected, frozen=False, entered=False)
    except Exception:
        checks["v2_unattested_comparison_rejected"] = True

    return {
        "schema": SCHEMA,
        "passed": all(checks.values()),
        "checks": checks,
        "synthetic_only": True,
        "nonzero_derivatives_supplied": True,
        "source_rows_independent": True,
        "receipt_read": False,
        "report_file_read": False,
        "ancestor_imported": False,
        "saved_source_graph_evaluated": False,
        "synthetic_formula_tree_evaluated": True,
        "midpoint_values_used": False,
        "log_absolute_upper_used_as_value": False,
        "old_error_removed": False,
    }


__all__ = [
    "SCHEMA", "SOURCE_FAMILY", "CELLS", "SOURCE_BINDINGS",
    "build_synthetic_receipt", "build_synthetic_endpoint_report",
    "build_synthetic_corrected_reserved_report", "independent_rational_formula_check",
    "mutate_duplicate_cell",
    "mutate_missing_cell", "mutate_source_corruption", "run_fixture",
]


