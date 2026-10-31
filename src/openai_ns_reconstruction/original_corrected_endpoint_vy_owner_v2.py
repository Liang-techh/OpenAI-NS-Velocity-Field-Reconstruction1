"""Data-only owner for corrected reserved lower-endpoint ``V_y`` rows.

The caller supplies already-decoded Python mappings.  This leaf never opens a
report, decodes a receipt, imports the endpoint producer, walks a graph, or
turns a logarithmic bound into a value.  Its output is a source-layer packet,
not a numerical endpoint evaluation.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from copy import deepcopy
from typing import Any


SCHEMA = "original_corrected_endpoint_vy_owner_v2"
COMPONENT = "corrected_reserved_lower_endpoint_actual_physical_V_y_v2"
CELLS = (("-1", "-.5"), ("-.5", "0"), ("0", ".5"), (".5", "1"))
ORDERS = (0, 1)
ALL_SOURCE_ORDERS = (0, 1, 2, 3)

SOURCE_BINDINGS = {
    "lower": "physical.source_rows(first_micro, ('sc_half',))",
    "lower_V_y": "actual_background_source.V_y[0:4]",
    "lower_U_ZZ": "factorial source rows V[2:4]",
    "lower_U_yZ": "factorial source rows V_y[1:4]",
    "common_Pstar": "outer.owner.owner(ends).Pstar",
}
LOWER_PROOF_VY_PATH = "actual_background_source.V_y"

RECEIPT_CHECKS = (
    "exact_four_Z_cell_cover",
    "lower_actual_V_y_rows_and_U_ZZ_checked",
    "physical_common_Pstar_once",
)


class EndpointVyOwnerError(ValueError):
    """Raised when a saved endpoint report cannot be admitted as a source."""


EndpointVYOwnerError = EndpointVyOwnerError


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise EndpointVyOwnerError(message)


def _mapping(value: Any, label: str) -> Mapping[str, Any]:
    _require(isinstance(value, Mapping), f"{label} must be a mapping")
    return value


def _sequence(value: Any, label: str) -> Sequence[Any]:
    _require(isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)),
             f"{label} must be a sequence")
    return value


def _cell(value: Any, label: str) -> tuple[str, str]:
    values = _sequence(value, label)
    _require(len(values) == 2 and all(isinstance(item, str) for item in values),
             f"{label} must contain two exact serialized endpoints")
    result = (values[0], values[1])
    _require(result in CELLS, f"{label} is not one of the four exact Z cells")
    return result


def _canonical_cells(value: Any, label: str) -> tuple[tuple[str, str], ...]:
    values = _sequence(value, label)
    result = tuple(_cell(item, f"{label}[{index}]") for index, item in enumerate(values))
    _require(len(result) == len(CELLS), f"{label} must contain exactly four cells")
    _require(len(set(result)) == len(CELLS), f"{label} contains a duplicate or missing Z cell")
    _require(set(result) == set(CELLS), f"{label} does not exactly cover the four Z cells")
    return result


def _copy(value: Any) -> Any:
    return deepcopy(value)


def _family(report_family: Any, value: Any, label: str) -> None:
    _require(value == report_family, f"{label} changed source_family")


def _unwrap_report(endpoint_report: Any) -> tuple[Mapping[str, Any], Any]:
    outer = _mapping(endpoint_report, "endpoint_report")
    report: Mapping[str, Any] = outer
    embedded_receipt: Any = None
    for key in ("endpoint_report", "report"):
        candidate = outer.get(key)
        if isinstance(candidate, Mapping) and (
            key == "endpoint_report"
            or any(name in outer for name in ("endpoint_receipt", "receipt", "checker", "report_checker"))
        ):
            report = candidate
            break
    for key in ("endpoint_receipt", "receipt", "checker", "report_checker", "check"):
        candidate = outer.get(key)
        if isinstance(candidate, Mapping):
            embedded_receipt = candidate
            break
    return report, embedded_receipt


def _find_receipt(report: Mapping[str, Any], explicit: Any, embedded: Any) -> Mapping[str, Any]:
    if explicit is not None:
        return _mapping(explicit, "endpoint_receipt")
    if embedded is not None:
        return _mapping(embedded, "embedded endpoint_receipt")
    for key in ("endpoint_receipt", "receipt", "checker", "report_checker", "check"):
        candidate = report.get(key)
        if isinstance(candidate, Mapping):
            return candidate
    if report.get("all_passed") is True:
        return report
    raise EndpointVyOwnerError("decoded endpoint receipt/checker is required")


def _receipt_check(receipt: Mapping[str, Any], name: str) -> bool:
    checks = receipt.get("checks")
    if isinstance(checks, Mapping) and checks.get(name) is True:
        return True
    return receipt.get(name) is True


def _validate_receipt(receipt: Mapping[str, Any], report_family: Any) -> dict[str, Any]:
    _require(receipt.get("all_passed") is True, "endpoint receipt/checker all_passed is not true")
    if "source_family" in receipt:
        _family(report_family, receipt.get("source_family"), "endpoint receipt")
    if "cells" in receipt:
        _canonical_cells(receipt.get("cells"), "endpoint receipt.cells")
    for name in RECEIPT_CHECKS:
        _require(_receipt_check(receipt, name), f"endpoint receipt checker missing {name}")
    return _copy(dict(receipt))


def _validate_binding(bindings: Mapping[str, Any], name: str, expected: str) -> None:
    _require(name in bindings and bindings.get(name) == expected,
             f"endpoint source_bindings.{name} changed or is missing")


def _row_order(row: Mapping[str, Any], order: int, label: str) -> None:
    for key in ("ordinary_Z_order", "ordinary_Z_derivative_order", "Z_order"):
        if key in row:
            _require(row.get(key) == order, f"{label}.{key} is not ordinary Z{order}")
    if "ordinary_Z" in row:
        _require(row.get("ordinary_Z") in (order, f"Z{order}"), f"{label}.ordinary_Z changed")


def _row_provenance(row: Mapping[str, Any]) -> Any:
    for key in ("provenance", "source_provenance", "origin", "source_origin"):
        if key in row:
            return _copy(row.get(key))
    return None


def _metadata_value(row: Mapping[str, Any], keys: Sequence[str], path: str) -> dict[str, Any]:
    for key in keys:
        if key in row:
            return {"source_key": key, "value": _copy(row.get(key)),
                    "materialized_as_function_value": False}
    return {"source_key": None, "symbolic_ref": path, "value_materialized": False,
            "source_field_absent": True}


def _source_record(row: Any, path: str, *, order: int | None = None,
                   require_enclosure_flags: bool = False) -> dict[str, Any]:
    original = _mapping(row, path)
    if order is not None:
        _row_order(original, order, path)
    if require_enclosure_flags:
        _require(original.get("encloses_original_source_function") is True,
                 f"{path} is not marked as an original source-function enclosure")
        _require(original.get("point_value_selected") is False,
                 f"{path} selected a point value")
        _require("coefficient_interval" in original,
                 f"{path} lost its original signed coefficient interval")
        _require("formal_positive_scale" in original,
                 f"{path} lost its formal_positive_scale provenance")
    return {
        "row": _copy(dict(original)),
        "provenance": _row_provenance(original),
        "source_path": path,
        "raw_row_preserved": True,
        "signed_coefficient": _metadata_value(
            original,
            ("signed_coefficient", "signed_coefficient_interval", "coefficient_interval", "coefficient"),
            path + ".signed_coefficient",
        ),
        "symbolic_scaling": _metadata_value(
            original,
            ("symbolic_scaling", "formal_scaling", "formal_positive_scale",
             "source_scaling", "scaling"),
            path + ".symbolic_scaling",
        ),
        "log_absolute_upper_not_used_as_value": True,
        "reported_sign_not_used_as_value": True,
        "midpoint_not_used_as_value": True,
    }


def _validate_family_in_row(row: Mapping[str, Any], family: Any, label: str) -> None:
    for key in ("source_family", "family"):
        if key in row:
            _family(family, row.get(key), label + "." + key)


def _same_source_provenance(rows: Sequence[Mapping[str, Any]], label: str) -> None:
    """Reject an explicitly split source while accepting legacy rows w/o ids."""

    ids: list[Any] = []
    for row in rows:
        provenance = _row_provenance(row)
        if isinstance(provenance, Mapping) and "same_source_id" in provenance:
            ids.append(provenance.get("same_source_id"))
    if ids:
        _require(all(item == ids[0] for item in ids), f"{label} does not retain one same-source provenance")


def _validate_vy_provenance(row: Mapping[str, Any], cell: tuple[str, str], order: int, label: str) -> None:
    provenance = _row_provenance(row)
    if not isinstance(provenance, Mapping):
        return
    if "source_path" in provenance:
        _require(
            isinstance(provenance.get("source_path"), str)
            and bool(provenance.get("source_path")),
            f"{label}.provenance source path missing",
        )
    if "row_order" in provenance:
        _require(provenance.get("row_order") == order, f"{label}.provenance row order changed")
    if "cell" in provenance:
        _require(provenance.get("cell") == list(cell), f"{label}.provenance cell changed")


def _validate_lower(sample: Mapping[str, Any], cell: tuple[str, str], family: Any,
                    sample_label: str) -> tuple[dict[str, Any], dict[str, Any]]:
    if "source_family" in sample:
        _family(family, sample.get("source_family"), sample_label)
    _require(sample.get("exact_Z_cell") == list(cell), f"{sample_label}.exact_Z_cell changed")
    sample_proof = _mapping(sample.get("proof"), sample_label + ".proof")
    _family(family, sample_proof.get("source_family"), sample_label + ".proof")
    _require(sample_proof.get("physical_common_Pstar_applied_once") is True,
             f"{sample_label}.proof lacks physical_common_Pstar_applied_once")

    lower = _mapping(sample.get("lower"), sample_label + ".lower")
    upper = sample.get("upper")
    if isinstance(upper, Mapping):
        if "source_family" in upper:
            _family(family, upper.get("source_family"), sample_label + ".upper")
        upper_proof = upper.get("proof")
        if isinstance(upper_proof, Mapping) and "source_family" in upper_proof:
            _family(family, upper_proof.get("source_family"), sample_label + ".upper.proof")
    lower_proof = _mapping(lower.get("proof"), sample_label + ".lower.proof")
    _family(family, lower_proof.get("source_family"), sample_label + ".lower.proof")
    _require(lower_proof.get("endpoint") == "first_micro/sc_half singleton",
             f"{sample_label}.lower is not the actual lower sc_half endpoint")
    _require(lower_proof.get("actual_V_y_rows_0_through_3_factorial_converted_once") is True,
             f"{sample_label}.lower proof lacks four-row factorial conversion guard")
    _require(lower_proof.get("common_Pstar_applied_once") is True,
             f"{sample_label}.lower proof lacks common_Pstar_applied_once")
    if "same_source" in lower_proof:
        _require(lower_proof.get("same_source") is True,
                 f"{sample_label}.lower proof does not attest same-source ownership")
    paths = lower_proof.get("source_quantity_paths")
    if isinstance(paths, Mapping):
        _require(paths.get("V_y") == LOWER_PROOF_VY_PATH,
                 f"{sample_label}.lower V_y source path changed")
        _require(paths.get("Pstar") == SOURCE_BINDINGS["common_Pstar"],
                 f"{sample_label}.lower Pstar source path changed")

    source_rows = _mapping(lower.get("source_rows"), sample_label + ".lower.source_rows")
    vy = _sequence(source_rows.get("V_y"), sample_label + ".lower.source_rows.V_y")
    _require(len(vy) >= 4, f"{sample_label}.lower.source_rows.V_y needs rows 0 through 3")
    raw_vy = [_mapping(row, f"{sample_label}.lower.source_rows.V_y[{order}]")
              for order, row in enumerate(vy[:4])]
    _same_source_provenance(raw_vy, sample_label + ".lower.source_rows.V_y")
    all_rows: list[dict[str, Any]] = []
    selected_rows: list[dict[str, Any]] = []
    for order, row_mapping in enumerate(raw_vy):
        _validate_family_in_row(row_mapping, family, f"{sample_label}.lower.source_rows.V_y[{order}]")
        _validate_vy_provenance(row_mapping, cell, order,
                                f"{sample_label}.lower.source_rows.V_y[{order}]")
        retained = _source_record(
            row_mapping,
            f"{sample_label}.lower.source_rows.V_y[{order}]",
            order=order,
            require_enclosure_flags=True,
        )
        all_rows.append(retained)
        if order in ORDERS:
            selected_rows.append(retained)
    _require(len(selected_rows) == 2, f"{sample_label} did not retain V_y C0/C1")

    R = _mapping(source_rows.get("R"), sample_label + ".lower.source_rows.R")
    Pstar = _mapping(source_rows.get("Pstar"), sample_label + ".lower.source_rows.Pstar")
    _validate_family_in_row(R, family, sample_label + ".lower.source_rows.R")
    _validate_family_in_row(Pstar, family, sample_label + ".lower.source_rows.Pstar")
    retained_rows = {
        "V_y": all_rows,
        "V_y_C0_C1": selected_rows,
        "R": _source_record(R, sample_label + ".lower.source_rows.R"),
        "Pstar": _source_record(Pstar, sample_label + ".lower.source_rows.Pstar"),
    }
    return (
        {"source_rows": retained_rows,
         "original_source_rows": _copy(dict(source_rows)),
         "proof": _copy(dict(lower_proof)),
         "original_lower": _copy(dict(lower))},
        _copy(dict(sample_proof)),
    )


def _factor(name: str, path: str) -> dict[str, Any]:
    return {"kind": "factor", "operation": "factor", "name": name, "path": path}


def _literal(value: int) -> dict[str, Any]:
    return {"kind": "literal", "operation": "literal", "value": value}


def _product(*factors: dict[str, Any]) -> dict[str, Any]:
    return {"kind": "product", "operation": "product", "factors": list(factors)}


def _sum(*terms: dict[str, Any]) -> dict[str, Any]:
    return {"kind": "sum", "operation": "sum", "terms": list(terms)}


def _difference(left: dict[str, Any], right: dict[str, Any]) -> dict[str, Any]:
    # Keep the evaluator vocabulary to factor/product/sum, with a semantic
    # annotation for signed subtraction rather than a hidden reduction.
    return {"kind": "sum", "operation": "sum",
            "terms": [left, _product(_literal(-1), right)],
            "semantic": "difference"}


def _quotient(numerator: dict[str, Any], denominator: dict[str, Any]) -> dict[str, Any]:
    return {"kind": "quotient", "operation": "quotient", "numerator": numerator,
            "denominator": denominator}


def _power(base: dict[str, Any], exponent: int) -> dict[str, Any]:
    return {"kind": "power", "operation": "power", "base": base, "exponent": exponent}


def _repair_formula(cell_label: str) -> dict[str, Any]:
    """Return the unevaluated shared-source repair AST for one cell."""

    R = _factor("R", f"{cell_label}.lower.source_rows.R")
    Pstar = _factor("Pstar", f"{cell_label}.lower.source_rows.Pstar")
    Y0 = _factor("Y0", f"{cell_label}.lower.source_rows.V_y[0]")
    Y1 = _factor("Y1", f"{cell_label}.lower.source_rows.V_y[1]")
    Z = _factor("Z", f"{cell_label}.exact_Z_cell")
    delta = _factor("delta", f"{cell_label}.symbolic_delta")
    denominator = _difference(_literal(1), _product(delta, _power(Z, 2)))
    F = _quotient(_product(_literal(2), Z), denominator)
    # d/dZ [2Z/(1-delta Z^2)] = (2+2*delta Z^2)/(1-delta Z^2)^2.
    Fprime = _quotient(_sum(_literal(2), _product(_literal(2), delta, _power(Z, 2))),
                       _power(denominator, 2))
    pstar_minus_one = _difference(Pstar, _literal(1))
    return {
        "scope": "shared-source only; no old error term is removed or inferred",
        "symbols": {"R": R, "Pstar": Pstar, "Y0": Y0, "Y1": Y1, "Z": Z, "delta": delta},
        "definitions": {"F": F, "Fprime": Fprime},
        "repair0": _product(R, pstar_minus_one,
                             _factor("F", f"{cell_label}.formula.definitions.F"), Y0),
        "repair1": _product(
            R, pstar_minus_one,
            _sum(
                _product(_factor("F", f"{cell_label}.formula.definitions.F"), Y1),
                _product(_factor("Fprime", f"{cell_label}.formula.definitions.Fprime"), Y0),
            ),
        ),
        "ordinary_Z1_rule": "R*(Pstar-1)*(F*Y1+Fprime*Y0); Y1 is the ordinary derivative",
        "F_formula": "2Z/(1-delta*Z^2)",
        "Fprime_formula": "(2+2*delta*Z^2)/(1-delta*Z^2)^2",
        "evaluated": False,
        "point_value_selected": False,
        "log_absolute_upper_used_as_value": False,
        "reported_sign_used_as_value": False,
        "midpoint_used_as_value": False,
    }


def _terminal_inputs(corrected_report: Any) -> Mapping[str, Any]:
    report, _ = _unwrap_report(corrected_report)
    candidate = report.get("corrected_terminal_source_inputs")
    if isinstance(candidate, Mapping):
        return candidate
    if all(key in report for key in ("physical_endpoint_source_rows", "physical_endpoint_source_proof",
                                     "numeric_endpoint_cell")):
        return report
    raise EndpointVyOwnerError("corrected report lacks corrected_terminal_source_inputs")


def _compare_corrected_terminal_source(
    frozen_report: Mapping[str, Any], corrected_report: Any, *, caller_attested: bool
) -> dict[str, Any]:
    _require(caller_attested is True,
             "caller must attest frozen=True and entered=True before terminal comparison")
    terminal = _terminal_inputs(corrected_report)
    if "source_family" in terminal:
        _family(frozen_report.get("source_family"), terminal.get("source_family"),
                "corrected terminal source")
    _require(_cell(terminal.get("numeric_endpoint_cell"), "corrected numeric endpoint cell") == ("0", ".5"),
             "corrected numeric endpoint cell is not [0,.5]")
    samples = _mapping(frozen_report.get("samples"), "frozen endpoint_report.samples")
    endpoint_sample = _mapping(samples.get(str(("0", ".5"))), "frozen lower endpoint sample")
    lower = _mapping(endpoint_sample.get("lower"), "frozen lower endpoint")
    _require(terminal.get("physical_endpoint_source_rows") == lower.get("source_rows"),
             "corrected terminal physical_endpoint_source_rows differs from saved lower.source_rows")
    _require(terminal.get("physical_endpoint_source_proof") == lower.get("proof"),
             "corrected terminal physical_endpoint_source_proof differs from saved lower.proof")
    return {
        "checked": True,
        "numeric_endpoint_cell": ["0", ".5"],
        "physical_endpoint_source_rows_match": True,
        "physical_endpoint_source_proof_match": True,
        "whole_saved_lower_rows_and_proof_compared": True,
        "comparison_after_caller_attested_freeze_and_enter": True,
        "caller_freeze_enter_attestation_only": True,
        "old_error_terms_removed": False,
    }


class EndpointVyOwner:
    """Source-layer owner built from an already-decoded report.

    The caller may attest that the input was frozen and entered; the class
    records that attestation but cannot establish it independently.
    """

    def __init__(self, endpoint_report: Any, endpoint_receipt: Any = None, *,
                 expected_source_family: Any = None, corrected_reserved_report: Any = None,
                 frozen: bool = False, entered: bool = False) -> None:
        _require(type(frozen) is bool and type(entered) is bool,
                 "frozen and entered must be explicit booleans")
        report, embedded_receipt = _unwrap_report(endpoint_report)
        self._report = _copy(dict(report))
        # The copy is local ownership; it does not prove that the caller froze
        # or entered the source.  Those are explicit runner attestations.
        self._frozen_attested = frozen
        self._entered_attested = entered
        source_receipt = _find_receipt(report, endpoint_receipt, embedded_receipt)
        self._receipt = _copy(dict(source_receipt))
        self._result = self._build(expected_source_family)
        if corrected_reserved_report is not None:
            self._result["corrected_terminal_source_comparison"] = _compare_corrected_terminal_source(
                self._report, corrected_reserved_report,
                caller_attested=self._frozen_attested and self._entered_attested)
        else:
            self._result["corrected_terminal_source_comparison"] = {
                "checked": False,
                "comparison_requires_corrected_reserved_report": True,
                "old_error_terms_removed": False,
            }

    def _build(self, expected_source_family: Any) -> dict[str, Any]:
        report = self._report
        family = report.get("source_family")
        _require(family is not None, "endpoint report.source_family is missing")
        if expected_source_family is not None:
            _family(expected_source_family, family, "endpoint report")
        _require(report.get("candidate_preRc_endpoint_C1_constructed") is True,
                 "endpoint report is not a constructed pre-Rc C1 report")
        _canonical_cells(report.get("cells"), "endpoint report.cells")
        bindings = _mapping(report.get("source_bindings"), "endpoint report.source_bindings")
        for name, expected in SOURCE_BINDINGS.items():
            _validate_binding(bindings, name, expected)
        receipt = _validate_receipt(self._receipt, family)
        samples = _mapping(report.get("samples"), "endpoint report.samples")
        _require(set(samples) == {str(cell) for cell in CELLS},
                 "endpoint report.samples has a duplicate, missing, or extra cell key")

        cell_sources: dict[str, Any] = {}
        for cell in CELLS:
            key = str(cell)
            label = f"samples[{key!r}]"
            sample = _mapping(samples.get(key), label)
            lower, sample_proof = _validate_lower(sample, cell, family, label)
            cell_sources[key] = {
                "cell": list(cell),
                "exact_Z_cell": list(cell),
                "source_family": _copy(family),
                "lower": lower,
                "source_enclosures": lower["source_rows"]["V_y_C0_C1"],
                "sample_proof": sample_proof,
                "formula": _repair_formula(f"cells[{key!r}]"),
                "same_source": {
                    "lower_endpoint": SOURCE_BINDINGS["lower"],
                    "source_family_checked": True,
                    "shared_Pstar_checked": True,
                    "Y_already_has_Pstar_once": True,
                    "no_extra_Pstar_applied_by_owner": True,
                },
            }
            cell_sources[key]["repair_formula"] = _copy(cell_sources[key]["formula"])
            cell_sources[key]["repair0"] = _copy(cell_sources[key]["formula"]["repair0"])
            cell_sources[key]["repair1"] = _copy(cell_sources[key]["formula"]["repair1"])
        result = {
            "schema": SCHEMA,
            "component": COMPONENT,
            "source_family": _copy(family),
            "cells": [list(cell) for cell in CELLS],
            "source_bindings": _copy(dict(bindings)),
            "receipt": receipt,
            "cell_sources": cell_sources,
            "ordinary_Z_orders": [0, 1],
            "all_source_orders": [0, 1, 2, 3],
            "ordinary_Z1_is_fixed_coordinate_derivative": True,
            "ordinary_Z1_is_not_a_second_factorial_conversion": True,
            "V_y_rows_0_through_3_retained": True,
            "V_y_C0_C1_rows_are_source_enclosures": True,
            "physical_common_Pstar_applied_once": True,
            "Y_already_includes_Pstar_once": True,
            "additional_Pstar_multiplication": False,
            "frozen_input_attested_by_caller": self._frozen_attested,
            "entered_input_attested_by_caller": self._entered_attested,
            "input_report_read": False,
            "report_file_opened": False,
            "receipt_decoded_here": False,
            "ancestor_imported": False,
            "graph_walked": False,
            "numeric_source_values_selected": False,
            "full_radial_field_claimed": False,
            "limits": {
                "source_function_identity_beyond_saved_report": False,
                "numeric_endpoint_function_values": False,
                "log_absolute_upper_as_function_value": False,
                "reported_sign_or_midpoint_as_function_value": False,
                "old_error_removed": False,
                "global_corrected_reserved_closure": False,
            },
        }
        result["repair_formulas"] = {
            key: _copy(value["formula"]) for key, value in cell_sources.items()
        }
        result["guarantees"] = {
            "four_exact_Z_cells": True,
            "lower_actual_sc_half_source": True,
            "same_source_family": True,
            "same_source_provenance_checked_when_serialized": True,
            "signed_rows_and_symbolic_scaling_retained": True,
            "R_Pstar_Y0_Y1_rows_retained_whole": True,
            "ordinary_Z1_not_reinterpreted_as_Taylor_coefficient": True,
            "formula_is_structural_and_unevaluated": True,
            "Y_Pstar_not_reapplied": True,
        }
        result["report_api"] = {
            "input": "already-decoded endpoint report and optional decoded checker receipt",
            "entrypoint": "build_owner(endpoint_report, endpoint_receipt=None, corrected_reserved_report=None, frozen=True, entered=True) [v2]",
            "cell_source_key": "cell_sources[str(tuple(cell))]",
            "formula_key": "cell_sources[str(tuple(cell))].formula",
            "terminal_compare": "corrected_terminal_source_inputs physical_endpoint_source_rows/proof",
        }
        return result

    @property
    def result(self) -> dict[str, Any]:
        return _copy(self._result)

    def source_layer(self) -> dict[str, Any]:
        return self.result

    def compare_corrected_terminal_source(self, corrected_reserved_report: Any) -> dict[str, Any]:
        comparison = _compare_corrected_terminal_source(
            self._report, corrected_reserved_report,
            caller_attested=self._frozen_attested and self._entered_attested)
        self._result["corrected_terminal_source_comparison"] = comparison
        return _copy(comparison)


def build_owner(endpoint_report: Any, endpoint_receipt: Any = None, *,
                expected_source_family: Any = None, corrected_reserved_report: Any = None,
                frozen: bool = False, entered: bool = False) -> dict[str, Any]:
    return EndpointVyOwner(
        endpoint_report, endpoint_receipt,
        expected_source_family=expected_source_family,
        corrected_reserved_report=corrected_reserved_report,
        frozen=frozen, entered=entered,
    ).result


def extract_actual_physical_v_y_owner(endpoint_report: Any, endpoint_receipt: Any = None, *,
                                      expected_source_family: Any = None,
                                      corrected_reserved_report: Any = None,
                                      frozen: bool = False, entered: bool = False) -> dict[str, Any]:
    return build_owner(endpoint_report, endpoint_receipt,
                       expected_source_family=expected_source_family,
                       corrected_reserved_report=corrected_reserved_report,
                       frozen=frozen, entered=entered)


def consume(endpoint_report: Any, endpoint_receipt: Any = None, *,
            expected_source_family: Any = None, corrected_reserved_report: Any = None,
            frozen: bool = False, entered: bool = False) -> dict[str, Any]:
    return build_owner(endpoint_report, endpoint_receipt,
                       expected_source_family=expected_source_family,
                       corrected_reserved_report=corrected_reserved_report,
                       frozen=frozen, entered=entered)


def build(endpoint_report: Any, endpoint_receipt: Any = None, *,
          expected_source_family: Any = None, corrected_reserved_report: Any = None,
          frozen: bool = False, entered: bool = False) -> dict[str, Any]:
    """Generic builder alias for one-shot runners."""

    return build_owner(
        endpoint_report, endpoint_receipt,
        expected_source_family=expected_source_family,
        corrected_reserved_report=corrected_reserved_report,
        frozen=frozen, entered=entered,
    )


def compare_corrected_terminal_source(
    endpoint_report: Any,
    corrected_reserved_report: Any,
    endpoint_receipt: Any = None,
    *,
    expected_source_family: Any = None,
    frozen: bool = False,
    entered: bool = False,
) -> dict[str, Any]:
    """Compare terminal source rows only under caller freeze/enter attestation."""

    owner = EndpointVyOwner(
        endpoint_report,
        endpoint_receipt,
        expected_source_family=expected_source_family,
        corrected_reserved_report=corrected_reserved_report,
        frozen=frozen,
        entered=entered,
    )
    return owner.result["corrected_terminal_source_comparison"]


__all__ = [
    "SCHEMA", "COMPONENT", "CELLS", "ORDERS", "ALL_SOURCE_ORDERS",
    "SOURCE_BINDINGS", "RECEIPT_CHECKS", "EndpointVyOwnerError",
    "EndpointVYOwnerError", "EndpointVyOwner", "build", "build_owner",
    "extract_actual_physical_v_y_owner", "consume", "compare_corrected_terminal_source",
]


