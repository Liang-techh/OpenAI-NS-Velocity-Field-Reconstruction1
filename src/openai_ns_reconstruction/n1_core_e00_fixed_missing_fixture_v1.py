"""Independent synthetic checks for the fixed E00 missing-field leaf.

The fixture owns fresh signed, nonconstant W0 and C0/B1 packets.  It checks
the polynomial fixed-target masses independently, including the axis cell,
the four ordinary-Z1 Leibniz products, formal-power addition, and unmerged
matrix/W0 provenance.  It does not read a receipt or import an ancestor.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Any


FIELDS = ("K1", "P1", "d_x_F1", "d_x_Uz1")
FIELD_ROWS = {"K1": "6", "P1": "4", "d_x_F1": "5", "d_x_Uz1": "6"}
C0_KEYS = (
    "4,1",
    "5,1",
    "5,2",
    "5,3",
    "5,5",
    "6,1",
    "6,2",
    "6,3",
    "6,4",
    "6,6",
)
B1_KEYS = ("5,1", "5,2", "5,3", "6,2", "6,3", "6,4")
Z0_ORIGINS = ("C0 W0", "B1 W0_Z")
Z1_ORIGINS = ("C0_Z W0", "C0 W0_Z", "B1_Z W0_Z", "B1 W0_ZZ")


def _q(c: Any, value: Fraction) -> Any:
    return c.mpf(value.numerator) / value.denominator


def _box(c: Any, lower: Fraction, upper: Fraction) -> Any:
    lo = _q(c, lower)
    hi = _q(c, upper)
    return c.make_mpf((lo._mpi_[0], hi._mpi_[1]))


def _t(c: Any, index: int) -> Any:
    left = _q(c, Fraction(index, 8))
    right = _q(c, Fraction(index + 1, 8))
    return c.make_mpf((left._mpi_[0], right._mpi_[1]))


def _signed_value(c: Any, cell: int, order: int, key_index: int, role: int) -> Any:
    sign = -1 if (3 * cell + 2 * order + key_index + role) % 5 in (0, 3) else 1
    magnitude = 2 + ((cell + 1) * (order + 2) * (key_index + 1) + role) % 13
    center = Fraction(sign * magnitude, 19)
    width = Fraction(1, 20000 + 101 * cell + 17 * order + key_index + role)
    return _box(c, center - width, center + width)


def _w0_term(c: Any, cell: int, order: int, column: int, term_index: int) -> dict[str, Any]:
    source_id = f"e00-fixture:W0:{cell}:Z{order}:col{column}:term{term_index}"
    return {
        "F0_power": (cell + order + column + term_index) % 3,
        "Lambda_power": -1 - ((2 * cell + order + column + term_index) % 4),
        "coefficient": _signed_value(c, cell, order, column, term_index),
        "source_id": source_id,
        "source_ids": [source_id],
        "source_path": ["outer_cells", cell, "W0_ordinary_Z", str(order), str(column), term_index],
        "source_cell": f"fixture-W0-cell-{cell}",
        "origin": {"source_ids": [source_id], "source_path": ["fixture", source_id]},
        "uncertainty": True,
        "error_classes": ["fixture_signed_source"],
    }


def _matrix_term(
    c: Any,
    cell: int,
    name: str,
    order: int,
    key: str,
    key_index: int,
    term_index: int,
) -> dict[str, Any]:
    source_id = f"e00-fixture:{name}:{cell}:Z{order}:{key}:term{term_index}"
    return {
        "F0_power": (cell + order + key_index + term_index + 1) % 3,
        "Lambda_power": -1 - ((cell + 2 * order + key_index + term_index) % 4),
        "coefficient": _signed_value(c, cell, order, key_index, term_index + 3),
        "source_id": source_id,
        "source_ids": [source_id],
        "source_path": ["matrix_cells", cell, name, str(order), key, term_index],
        "source_cell": f"fixture-{name}-cell-{cell}",
        "origin": {"source_ids": [source_id], "source_path": ["fixture", source_id]},
        "uncertainty": True,
        "error_classes": ["fixture_signed_matrix"],
    }


def _inputs(c: Any) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    outer: list[dict[str, Any]] = []
    matrices: list[dict[str, Any]] = []
    for cell in range(8):
        t = _t(c, cell)
        w0: dict[str, dict[str, list[dict[str, Any]]]] = {}
        for order in range(3):
            w0[str(order)] = {}
            for column in range(1, 7):
                count = 2 if (cell + order + column) % 3 == 0 else 1
                w0[str(order)][str(column)] = [
                    _w0_term(c, cell, order, column, term_index)
                    for term_index in range(count)
                ]
        outer.append(
            {
                "t": t,
                "rho": 4 * t * t,
                "W0_ordinary_Z": w0,
                "cell_id": f"fixture-W0-cell-{cell}",
            }
        )

        c0: dict[str, dict[str, list[dict[str, Any]]]] = {"0": {}, "1": {}}
        b1: dict[str, dict[str, list[dict[str, Any]]]] = {"0": {}, "1": {}}
        for order in range(2):
            for key_index, key in enumerate(C0_KEYS):
                count = 2 if (cell + order + key_index) % 5 == 0 else 1
                c0[str(order)][key] = [
                    _matrix_term(c, cell, "C0", order, key, key_index, term_index)
                    for term_index in range(count)
                ]
            for key_index, key in enumerate(B1_KEYS):
                count = 2 if (cell + 2 * order + key_index) % 4 == 0 else 1
                b1[str(order)][key] = [
                    _matrix_term(c, cell, "B1", order, key, key_index + 11, term_index)
                    for term_index in range(count)
                ]
        matrices.append(
            {
                "t": t,
                "rho": 4 * t * t,
                "C0": c0,
                "B1": b1,
                "cell_id": f"fixture-matrix-cell-{cell}",
            }
        )
    return outer, matrices


def _analytic_masses(c: Any, index: int) -> dict[str, Any]:
    """Independent endpoint primitives for target ``t=1``."""

    a = _q(c, Fraction(index, 8))
    b = _q(c, Fraction(index + 1, 8))
    a2 = a * a
    b2 = b * b
    a4 = a2 * a2
    b4 = b2 * b2
    return {
        "K1": -(b2 - a2) + (b4 - a4) / 2,
        "P1": 2 * (b - a),
        "d_x_F1": (b4 - a4) / 2,
        "d_x_Uz1": b2 - a2,
    }


def _contains(actual: Any, expected: Any) -> bool:
    return actual.a <= expected.a and actual.b >= expected.b


def _expected_count(
    outer_cell: dict[str, Any],
    matrix_cell: dict[str, Any],
    field: str,
    output_order: int,
) -> int:
    row = FIELD_ROWS[field]
    rules = (
        (("C0", "0", "0"), ("B1", "0", "1"))
        if output_order == 0
        else (("C0", "1", "0"), ("C0", "0", "1"), ("B1", "1", "1"), ("B1", "0", "2"))
    )
    count = 0
    for name, matrix_order, w0_order in rules:
        keys = C0_KEYS if name == "C0" else B1_KEYS
        for key in keys:
            if key.split(",", 1)[0] != row:
                continue
            matrix_terms = matrix_cell[name][matrix_order][key]
            column = key.split(",", 1)[1]
            w0_terms = outer_cell["W0_ordinary_Z"][w0_order][column]
            count += len(matrix_terms) * len(w0_terms)
    return count


def run_fixture(leaf: Any, c: Any) -> dict[str, Any]:
    """Run fresh signed mass, Z1, and provenance checks."""

    outer, matrices = _inputs(c)
    family = {"fixture_family": "e00-fixed-missing", "generation": 1}
    logs = {"F0": "fixture-formal-F0", "Lambda": "fixture-formal-Lambda"}
    builder = getattr(leaf, "build", leaf)
    result = builder(
        c,
        outer_cells=outer,
        matrix_cells=matrices,
        source_family=family,
        source_logs=logs,
    )
    checks: dict[str, bool] = {}
    checks["schema"] = result.get("schema") == "n1_core_e00_fixed_missing_v1"
    checks["source_identity"] = result.get("source_family") == family and result.get("source_logs") == logs
    checks["four_fields_only"] = result.get("fields") == list(FIELDS)
    checks["old_F1_Uz1_not_recomputed"] = (
        result.get("prior_fixed_F1_Uz1_recomputed") is False
        and result.get("F1_Uz1_included") is False
        and "F1" not in result.get("fields", [])
        and "Uz1" not in result.get("fields", [])
    )
    checks["fixed_target"] = (
        result.get("fixed_target", {}).get("fixed") is True
        and result["fixed_target"]["x"].a == 2
        and result["fixed_target"]["x"].b == 2
        and result["fixed_target"]["t"].a == 1
        and result["fixed_target"]["t"].b == 1
    )
    output = result.get("E00_fixed_by_source_cell", {})
    checks["eight_source_cells"] = set(output) == {str(index) for index in range(8)}
    checks["ordinary_orders"] = result.get("ordinary_Z_orders") == [0, 1]
    checks["input_orders"] = result.get("input_W0_ordinary_Z_orders") == [0, 1, 2]
    checks["axis_marked"] = bool(output) and output["0"]["kernel_masses"].get("axis_endpoint_continuation") is True

    mass_checks = True
    sign_checks = True
    count_checks = True
    metadata_checks = True
    z1_checks = True
    for index in range(8):
        cell = output.get(str(index), {})
        actual_masses = cell.get("kernel_masses", {})
        expected_masses = _analytic_masses(c, index)
        for field in FIELDS:
            actual = actual_masses.get(field)
            expected = expected_masses[field]
            mass_checks = mass_checks and actual is not None and _contains(actual, expected)
            if field == "K1":
                sign_checks = sign_checks and actual is not None and actual.b <= 0
            else:
                sign_checks = sign_checks and actual is not None and actual.a >= 0
            for order in (0, 1):
                terms = cell.get("E00", {}).get(f"Z{order}", {}).get(field, [])
                count_checks = count_checks and len(terms) == _expected_count(
                    outer[index], matrices[index], field, order
                )
                for term in terms:
                    metadata_checks = metadata_checks and (
                        term.get("outer_cell") == index
                        and term.get("source_cell_index") == index
                        and term.get("matrix_source_cell_index") == index
                        and term.get("matrix_key")
                        and term.get("matrix_order") in ("0", "1")
                        and term.get("matrix_term_index") is not None
                        and term.get("W0_order") in ("0", "1", "2")
                        and term.get("W0_term_index") is not None
                        and len(term.get("source_ids", [])) >= 2
                        and all(str(source_id).startswith("e00-fixture:") for source_id in term.get("source_ids", []))
                        and len(term.get("source_paths", [])) == 2
                        and term.get("formal_powers_retained", True) is True
                    )
                if order == 1:
                    expected_origins = set(Z1_ORIGINS)
                    if FIELD_ROWS[field] == "4":
                        expected_origins = {"C0_Z W0", "C0 W0_Z"}
                    origins = {term.get("origin") for term in terms}
                    z1_checks = z1_checks and expected_origins.issubset(origins)
        if index == 0:
            first_z0 = cell.get("E00", {}).get("Z0", {})
            first_z1 = cell.get("E00", {}).get("Z1", {})
            for field in FIELDS:
                z0_terms = first_z0.get(field, [])
                z1_terms = first_z1.get(field, [])
                if z0_terms and z1_terms:
                    metadata_checks = metadata_checks and (
                        z0_terms[0]["kernel_mass"]._mpi_ == z1_terms[0]["kernel_mass"]._mpi_
                    )

    checks["analytic_polynomial_masses"] = mass_checks
    checks["signed_kernel_masses"] = sign_checks
    checks["all_unmerged_term_counts"] = count_checks
    checks["ordinary_Z1_product_rule"] = z1_checks
    checks["provenance_and_key_metadata"] = metadata_checks
    checks["axis_polynomial_values"] = (
        _contains(output["0"]["kernel_masses"]["P1"], _q(c, Fraction(1, 4)))
        and _contains(output["0"]["kernel_masses"]["d_x_F1"], _q(c, Fraction(1, 8192)))
        and _contains(output["0"]["kernel_masses"]["d_x_Uz1"], _q(c, Fraction(1, 64)))
        and _contains(output["0"]["kernel_masses"]["K1"], -_q(c, Fraction(127, 8192)))
    )
    checks["source_signs_are_fresh"] = (
        any(term["coefficient"].a < 0 for cell in outer for order in cell["W0_ordinary_Z"].values() for terms in order.values() for term in terms)
        and any(term["coefficient"].a > 0 for cell in outer for order in cell["W0_ordinary_Z"].values() for terms in order.values() for term in terms)
    )
    checks["no_QP_pressure_tail_ledger"] = all(
        result.get(key) is False
        for key in (
            "QP_added",
            "Q_branch_added",
            "P_branch_added",
            "allword_tail_added",
            "feedback_tail_added",
            "m15_computed",
            "ledger_appended",
        )
    ) and result.get("P1_is_m15") is False
    checks["saved_W0_not_reconstructed"] = (
        result.get("saved_outer_cells_reused") is True
        and result.get("W0_reconstructed") is False
        and result.get("ancestor_rerun") is False
    )

    return {
        "passed": all(checks.values()),
        "checks": checks,
        "check_count": len(checks),
        "synthetic_only": True,
        "fixture_is_not_actual_field": True,
        "independent_source_is_nonconstant": True,
        "analytic_polynomial_masses": True,
        "method": "fresh signed W0/C0/B1 cells with fixed-target polynomial masses and ordinary Z1 products",
        "scope": "bookkeeping fixture only; no receipt or physical-field certification",
    }


run = run_fixture


__all__ = ["run_fixture", "run"]

