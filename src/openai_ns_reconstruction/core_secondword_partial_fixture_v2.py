"""Independent synthetic checks for :mod:`core_secondword_partial_v2`.

The fixture builds a saved-profile-shaped input directly.  It does not import
the firstword producer, read a receipt, or construct an ancestor.  The values
are deliberately signed and nonzero so that a missing Leibniz branch or an
accidental absolute-value/cancellation operation is observable.
"""

from __future__ import annotations

import traceback
from collections.abc import Mapping, Sequence
from fractions import Fraction
from typing import Any


SCHEMA = "core_secondword_partial_fixture_v2"
INPUT_SCHEMA = "n1_core_w1_profile_v2"
OUTPUT_FIELDS = ("F1", "d_x_F1")
STATE_FIELDS = ("F1", "Uz1", "K1", "P1", "d_x_F1", "d_x_Uz1")
C0_KEYS = ("4,1", "5,1", "5,2", "5,3", "5,5", "6,1", "6,2", "6,3", "6,4", "6,6")
B1_KEYS = ("5,1", "5,2", "5,3", "6,2", "6,3", "6,4")

__all__ = ["SCHEMA", "INPUT_SCHEMA", "C0_KEYS", "B1_KEYS", "run_fixture", "run"]


def _q(c: Any, value: Fraction | int) -> Any:
    rational = Fraction(value)
    return c.mpf(rational.numerator) / rational.denominator


def _point(c: Any, value: Fraction | int) -> Any:
    return c.mpf(Fraction(value).numerator) / Fraction(value).denominator


def _term(c: Any, coefficient: Fraction | int, f0: int, lam: int, node_id: str, *, errors: Sequence[str] = (), origin: str | None = None, source_cell: int | None = None) -> dict[str, Any]:
    result = {
        "F0_power": f0,
        "Lambda_power": lam,
        "coefficient": _point(c, coefficient),
        "node_id": node_id,
        "source_id": node_id,
        "source_ids": [node_id],
        "source_paths": [["fixture_v2", "source", node_id]],
        "error_classes": list(errors),
        "source_error_classes": list(errors),
        "uncertainty": bool(errors),
        "error_only": False,
        "symmetric_error": False,
        "signed_field": True,
    }
    if origin is not None:
        result["origin"] = origin
    if source_cell is not None:
        result["source_cell"] = source_cell
        result["source_cell_index"] = source_cell
    return result


def _cell_t(c: Any, index: int) -> Any:
    return c.make_mpf((_point(c, Fraction(index, 8))._mpi_[0], _point(c, Fraction(index + 1, 8))._mpi_[1]))


def _profile(c: Any) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    targets: list[dict[str, Any]] = []
    matrices: list[dict[str, Any]] = []
    for index in range(8):
        t = _cell_t(c, index)
        w1: dict[str, dict[str, list[dict[str, Any]]]] = {"Z0": {}, "Z1": {}}
        for order in ("Z0", "Z1"):
            order_index = int(order[-1])
            for field_index, field in enumerate(STATE_FIELDS, start=1):
                base = (-1 if (2 * index + 3 * order_index + field_index) % 5 == 0 else 1) * (3 * index + 2 * field_index + order_index + 5)
                if field == "Uz1":
                    base *= -2
                if field == "K1":
                    base *= 3
                if field == "d_x_F1":
                    base *= -3
                powers = ((2 * index + field_index) % 5 - 2, -2 * order_index - ((field_index + 1) % 3))
                errors = (f"fixture_v2:firstword:cell:{index}",) if (index + 2 * field_index + order_index) % 3 == 0 else ()
                w1[order][field] = [
                    _term(
                        c,
                        base,
                        powers[0],
                        powers[1],
                        f"state:W1:cell:{index}:{order}:{field}",
                        errors=errors,
                        origin=f"fixture W1 {order} {field}",
                        source_cell=index,
                    )
                ]
                if field in ("F1", "Uz1") and order == "Z1":
                    # A second signed state term makes Leibniz/product ordering
                    # visible without ever summing or cancelling terms.
                    w1[order][field].append(
                        _term(
                            c,
                            -base // 3 if base % 3 == 0 else -(base + index + 1),
                            powers[0] + 1,
                            powers[1] - 1,
                            f"state:W1:cell:{index}:{order}:{field}:extra",
                            errors=(f"fixture_v2:firstword:extra:{index}",),
                            origin=f"fixture W1 {order} {field} extra",
                            source_cell=index,
                        )
                    )
        targets.append({
            "index": index,
            "t": t,
            "rho": 4 * t * t,
            "W1_by_order": w1,
        })

        c0: dict[str, dict[str, list[dict[str, Any]]]] = {"0": {}, "1": {}}
        b1: dict[str, dict[str, list[dict[str, Any]]]] = {"0": {}, "1": {}}
        for order in ("0", "1"):
            for key_index, key in enumerate(C0_KEYS):
                row, column = (int(value) for value in key.split(","))
                sign = -1 if (2 * index + 3 * int(order) + key_index + 1) % 5 == 0 else 1
                c0[order][key] = [
                    _term(c, sign * (2 * index + 3 * key_index + 5), (2 * index + key_index + 1) % 4 - 1, -2 * int(order) - 2, f"matrix:C0:{index}:{order}:{key}", errors=(f"fixture_v2:matrix:C0:{index}",) if (index + 2 * key_index) % 4 == 0 else (), source_cell=index)
                ]
            for key_index, key in enumerate(B1_KEYS):
                sign = -1 if (3 * index + int(order) + 2 * key_index + 1) % 4 == 0 else 1
                b1[order][key] = [
                    _term(c, sign * (2 * index + 2 * key_index + 7), (index + 2 * int(order) + 2) % 4, -2 * int(order) - 3, f"matrix:B1:{index}:{order}:{key}", errors=(f"fixture_v2:matrix:B1:{index}",) if (index + key_index + int(order)) % 3 == 1 else (), source_cell=index)
                ]
        matrices.append({"index": index, "t": t, "rho": 4 * t * t, "C0": c0, "B1": b1})
    return {"schema": INPUT_SCHEMA, "target_cells": targets}, matrices


def _clone(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {key: _clone(item) for key, item in value.items()}
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return [_clone(item) for item in value]
    return value


def _rejected(call: Any) -> bool:
    try:
        call()
    except (TypeError, ValueError, KeyError, ZeroDivisionError):
        return True
    return False


def _inside(value: Any, expected: Fraction) -> bool:
    point = _q(_inside.context, expected)
    return bool(value.a <= point.a and point.b <= value.b)


def run_fixture(build: Any, c: Any) -> dict[str, Any]:
    """Run independent shape, sign, exact-kernel, and rejection checks."""

    _inside.context = c
    build_fn = getattr(build, "build", build)
    if not callable(build_fn):
        raise TypeError("build must be callable or expose build")
    profile, matrices = _profile(c)
    try:
        result = build_fn(c, profile, matrices)
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
        }

    checks: dict[str, bool] = {}
    checks["build_completed"] = isinstance(result, Mapping)
    checks["schema"] = result.get("schema") == "core_secondword_partial_v2"
    checks["eight_source_cells"] = len(result.get("source_cells", [])) == 8
    checks["finite_subset"] = result.get("complete_secondword") is False and result.get("full_secondword_claimed") is False and result.get("physical_field_claimed") is False
    checks["no_false_Z2"] = result.get("Z2_materialized") is False and 2 not in result.get("ordinary_Z_orders", []) and "Z2" not in result.get("terms_by_order", {})
    terms = result.get("terms_by_order", {})
    f1_terms = list(terms.get("Z0", {}).get("F1", [])) + list(terms.get("Z1", {}).get("F1", []))
    fx_terms = list(terms.get("Z0", {}).get("d_x_F1", [])) + list(terms.get("Z1", {}).get("d_x_F1", []))
    all_f1_terms = f1_terms
    checks["nonzero_F1"] = bool(terms.get("Z0", {}).get("F1")) and bool(terms.get("Z1", {}).get("F1"))
    checks["nonzero_Fx"] = bool(terms.get("Z0", {}).get("d_x_F1")) and bool(terms.get("Z1", {}).get("d_x_F1"))
    checks["signed_coefficients"] = any(term["coefficient"].a < 0 for term in f1_terms) and any(term["coefficient"].b > 0 for term in f1_terms)
    checks["Z1_Leibniz_origins"] = {term.get("origin") for term in result.get("terms_by_order", {}).get("Z1", {}).get("F1", [])} >= {"C0[1]*W[0]", "C0[0]*W[1]", "B1[1]*W[1]"}
    z1_terms = result.get("terms_by_order", {}).get("Z1", {}).get("F1", [])
    checks["Z1_state_order_provenance"] = all(term.get("state_order") == ("0" if term.get("origin") == "C0[1]*W[0]" else "1") for term in z1_terms)
    missing_w2_ids = set(result.get("unknown_ids_by_kind", {}).get("missing_B1_W1_Z2", []))
    missing_w2 = [item for item in result.get("unknowns", []) if item.get("unknown_id") in missing_w2_ids]
    checks["missing_B1_W2_explicit"] = bool(missing_w2) and {item.get("column") for item in missing_w2} == {1, 2, 3} and all(item.get("included") is False and item.get("materialized") is False and item.get("state_order") == "2" for item in missing_w2)
    checks["known_C0_columns_1_2_3_5"] = {term.get("column") for term in all_f1_terms if term.get("matrix_name") == "C0"} == {1, 2, 3, 5}
    checks["known_B1_columns_1_2_3"] = {term.get("column") for term in all_f1_terms if term.get("matrix_name") == "B1"} == {1, 2, 3}
    checks["no_fabricated_B1_column5"] = all(term.get("matrix_key") != "5,5" for term in f1_terms if term.get("matrix_name") == "B1")
    checks["kernel_flags"] = result.get("kernel_mass_exact") is True and result.get("kernel_mass_positive") is True and result.get("kernel_mass_independent_polynomial") is True
    checks["formal_atoms_not_materialized"] = result.get("formal_atoms_materialized") is False and result.get("F0_materialized") is False and result.get("Lambda_materialized") is False
    checks["no_cancellation"] = result.get("no_cancellation") is True and result.get("cancellation_performed") is False and all(term.get("no_cancellation") is True for term in f1_terms)
    checks["source_state_refs"] = all(term.get("source_node_ref") and term.get("state_node_ref") and term.get("source_cell_index") == term.get("source_cell") for term in f1_terms)
    checks["error_classes_retained"] = any(term.get("error_classes") for term in f1_terms) and all("error_classes" in item for item in result.get("unknowns", []))
    checks["formal_powers_combined"] = all(isinstance(term.get("F0_power"), int) and isinstance(term.get("Lambda_power"), int) for term in f1_terms)

    expected_kernels = {}
    for index in range(8):
        a = Fraction(index, 8)
        b = Fraction(index + 1, 8)
        expected_kernels[str(index)] = {
            "F1": (b * b - a * a) - (b**4 - a**4) / 2,
            "d_x_F1": (b**4 - a**4) / 2,
        }
    checks["F_kernel_exact_polynomial"] = all(_inside(result["kernels"][key]["F1"], value["F1"]) for key, value in expected_kernels.items())
    checks["Fx_kernel_exact_polynomial"] = all(_inside(result["kernels"][key]["d_x_F1"], value["d_x_F1"]) for key, value in expected_kernels.items())
    checks["kernel_independent_of_profile"] = result["kernels"]["0"]["F1"]._mpi_ == result["kernels"][result["source_cells"][0]["kernel_id"]]["F1"]._mpi_

    wrong_cells = _clone(profile)
    wrong_cells["target_cells"][3]["t"] = _cell_t(c, 4)
    checks["wrong_cells_rejected"] = _rejected(lambda: build_fn(c, wrong_cells, matrices))
    missing_matrix = _clone(matrices)
    del missing_matrix[2]["C0"]["0"]["5,1"]
    checks["missing_matrix_rejected"] = _rejected(lambda: build_fn(c, profile, missing_matrix))
    missing_state = _clone(profile)
    del missing_state["target_cells"][1]["W1_by_order"]["Z1"]["F1"]
    checks["missing_state_rejected"] = _rejected(lambda: build_fn(c, missing_state, matrices))
    malformed_power = _clone(profile)
    del malformed_power["target_cells"][0]["W1_by_order"]["Z0"]["F1"][0]["F0_power"]
    checks["missing_power_rejected"] = _rejected(lambda: build_fn(c, malformed_power, matrices))

    return {
        "schema": SCHEMA,
        "passed": all(checks.values()),
        "checks": checks,
        "check_count": len(checks),
        "synthetic_only": True,
        "implementation_imported": False,
        "ordinary_Z1_Leibniz_checked": True,
        "exact_F_and_Fx_kernel_checked": True,
        "unknown_B1_W2_and_columns_checked": True,
    }


def run(build: Any, c: Any) -> dict[str, Any]:
    return run_fixture(build, c)



