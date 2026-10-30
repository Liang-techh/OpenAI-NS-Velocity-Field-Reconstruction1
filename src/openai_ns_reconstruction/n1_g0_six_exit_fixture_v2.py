"""Independent synthetic checks for :mod:`n1_g0_six_exit_v2`.

Eight source cells carry deliberately nonconstant, signed formal sectors for
all three forcing rows and both ordinary derivative orders.  The expected
coefficients use independent endpoint primitives for the six fixed-target
Green moments.  No producer, old builder, receipt, midpoint, or source value
is imported or selected.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Any


FIELDS = ("F1", "Uz1", "K1", "P1", "d_x_F1", "d_x_Uz1")
FIELD_ROWS = {
    "F1": "5",
    "Uz1": "6",
    "K1": "6",
    "P1": "4",
    "d_x_F1": "5",
    "d_x_Uz1": "6",
}
ROWS = ("4", "5", "6")


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


def _signed_coefficient(c: Any, cell: int, order: int, row: str, term: int) -> Any:
    row_index = int(row) - 4
    sign = -1 if (3 * cell + 2 * order + row_index + term) % 5 in (0, 3) else 1
    magnitude = 2 + ((cell + 1) * (order + 2) * (row_index + 1) + term) % 11
    center = Fraction(sign * magnitude, 17)
    width = Fraction(1, 17000 + 100 * cell + 10 * order + term)
    return _box(c, center - width, center + width)


def _term(c: Any, cell: int, order: int, row: str, term: int) -> dict[str, Any]:
    f0_power = (cell + order + int(row) + term) % 3
    lambda_power = -1 - ((2 * cell + order + int(row) + term) % 4)
    source = f"synthetic_g:{cell}:Z{order}:row{row}:term{term}"
    return {
        "F0_power": f0_power,
        "Lambda_power": lambda_power,
        "coefficient": _signed_coefficient(c, cell, order, row, term),
        "source_id": source,
        "source_ids": [source],
        "source_cell": f"synthetic-cell-{cell}",
        "origin": {
            "source_ids": [source],
            "error_classes": ["synthetic_signed_interval"],
            "uncertainty": True,
        },
        "error_classes": ["synthetic_signed_interval"],
        "uncertainty": True,
    }


def _inputs(c: Any) -> list[dict[str, Any]]:
    cells: list[dict[str, Any]] = []
    for cell in range(8):
        t = _t(c, cell)
        ordinary: dict[str, dict[str, list[dict[str, Any]]]] = {}
        for order in (0, 1):
            ordinary[str(order)] = {}
            for row in ROWS:
                count = 2 if (cell + order + int(row)) % 3 == 0 else 1
                ordinary[str(order)][row] = [
                    _term(c, cell, order, row, term) for term in range(count)
                ]
        cells.append({
            "t": t,
            "x": 2 * t,
            "rho": 4 * t * t,
            "g": ordinary,
            "cell_id": f"synthetic-cell-{cell}",
        })
    return cells


def _analytic_masses(c: Any, index: int) -> dict[str, Any]:
    """Independent endpoint primitives at target ``T=1``."""

    a = _q(c, Fraction(index, 8))
    b = _q(c, Fraction(index + 1, 8))
    a2 = a * a
    b2 = b * b
    a4 = a2 * a2
    b4 = b2 * b2
    f1 = (b2 - a2) - (b4 - a4) / 2
    if index == 0:
        primitive_a = c.mpf(0)
    else:
        primitive_a = 2 * a2 * c.log(1 / a) + a2
    primitive_b = 2 * b2 * c.log(1 / b) + b2
    uz1 = primitive_b - primitive_a
    p1 = 2 * (b - a)
    dxf1 = (b4 - a4) / 2
    dxuz1 = b2 - a2
    return {
        "F1": f1,
        "Uz1": uz1,
        "K1": -f1,
        "P1": p1,
        "d_x_F1": dxf1,
        "d_x_Uz1": dxuz1,
    }


def _endpoint(value: Any) -> Any:
    return value.a, value.b


def _contains(actual: Any, expected: Any) -> bool:
    actual_lo, actual_hi = _endpoint(actual)
    expected_lo, expected_hi = _endpoint(expected)
    return actual_lo <= expected_lo and actual_hi >= expected_hi


def _term_key(term: Any) -> tuple[Any, ...]:
    return (
        term.get("source_id"),
        term.get("field"),
        term.get("ordinary_Z_order"),
        term.get("source_cell_index"),
    )


def run_fixture(leaf: Any, c: Any) -> dict[str, Any]:
    """Run exact source-cell, sign, provenance, and mass checks."""

    cells = _inputs(c)
    builder = getattr(leaf, "build", leaf)
    result = builder(c, cells)
    checks: dict[str, bool] = {}

    checks["schema"] = result.get("schema") == "n1_g0_six_exit_v2"
    checks["fixed_x2"] = (
        result.get("exit_target", {}).get("fixed") is True
        and result["exit_target"]["x"].a == 2
        and result["exit_target"]["x"].b == 2
        and result["exit_target"].get("independent") is True
    )
    checks["eight_source_cells"] = (
        result.get("source_cell_count") == 8
        and set(result.get("G0g_exit_by_source_cell", {})) == {str(i) for i in range(8)}
    )
    checks["three_forcing_rows"] = result.get("forcing_rows") == ["4", "5", "6"]
    checks["six_fields"] = result.get("fields") == list(FIELDS)
    checks["ordinary_orders"] = result.get("ordinary_Z_orders") == [0, 1]
    checks["exact_first_F1_mass"] = (
        checks["eight_source_cells"]
        and _contains(
            result["G0g_exit_by_source_cell"]["0"]["kernel_masses"]["F1"],
            _q(c, Fraction(127, 8192)),
        )
    )
    checks["exact_first_P1_weighted_mass"] = (
        checks["eight_source_cells"]
        and _contains(
            result["G0g_exit_by_source_cell"]["0"]["kernel_masses"]["P1"],
            _q(c, Fraction(1, 4)),
        )
    )
    checks["nonconstant_cell_masses"] = (
        result["G0g_exit_by_source_cell"]["0"]["kernel_masses"]["F1"]._mpi_
        != result["G0g_exit_by_source_cell"]["7"]["kernel_masses"]["F1"]._mpi_
    )

    signs = [
        term["coefficient"]
        for cell in cells
        for order in cell["g"].values()
        for terms in order.values()
        for term in terms
    ]
    checks["signed_nonconstant_sources"] = any(value.b < 0 for value in signs) and any(value.a > 0 for value in signs)

    expected_counts = 0
    actual_counts = 0
    provenance = True
    coefficient_inclusion = True
    output = result["G0g_exit_by_source_cell"]
    for cell_index, cell in enumerate(cells):
        masses = _analytic_masses(c, cell_index)
        actual_masses = output[str(cell_index)]["kernel_masses"]
        for field in FIELDS:
            provenance = provenance and _contains(actual_masses[field], masses[field])
            provenance = provenance and _contains(masses[field], actual_masses[field])
        for order in (0, 1):
            for field in FIELDS:
                row = FIELD_ROWS[field]
                expected_terms = cell["g"][str(order)][row]
                actual_terms = output[str(cell_index)]["G0g"][f"Z{order}"][field]
                expected_counts += len(expected_terms)
                actual_counts += len(actual_terms)
                provenance = provenance and all(
                    term.get("source_cell_index") == cell_index
                    and term.get("source_row") == row
                    and term.get("ordinary_Z_order") == order
                    and term.get("principal_kernel") == "G0"
                    and term.get("exact_source_cell_mass") is True
                    for term in actual_terms
                )
                expected_by_id = {term["source_id"]: term for term in expected_terms}
                actual_by_id = {term["source_id"]: term for term in actual_terms}
                provenance = provenance and set(expected_by_id) == set(actual_by_id)
                for source_id, expected in expected_by_id.items():
                    got = actual_by_id.get(source_id)
                    if got is None:
                        coefficient_inclusion = False
                        continue
                    coefficient_inclusion = coefficient_inclusion and got["F0_power"] == expected["F0_power"]
                    coefficient_inclusion = coefficient_inclusion and got["Lambda_power"] == expected["Lambda_power"]
                    coefficient_inclusion = coefficient_inclusion and _contains(
                        got["coefficient"], expected["coefficient"] * masses[field]
                    )

    checks["all_term_counts"] = actual_counts == expected_counts
    checks["source_ids_and_powers"] = provenance
    checks["coefficient_inclusion"] = coefficient_inclusion
    checks["Z1_nonzero"] = any(
        term["coefficient"].a != 0 or term["coefficient"].b != 0
        for term in result["G0g_Z1_core_exit"]["F1"]
    )
    checks["K1_negative_kernel"] = all(
        output[str(index)]["kernel_masses"]["K1"].b <= 0
        for index in range(8)
    )
    checks["no_midpoint_or_reconstruction"] = (
        result.get("interval_midpoint_used") is False
        and result.get("W0_reconstructed") is False
        and result.get("old_builder_called") is False
    )

    return {
        "passed": all(checks.values()),
        "checks": checks,
        "check_count": len(checks),
        "synthetic_only": True,
        "fixture_is_not_actual_field": True,
        "independent_source_is_nonconstant": True,
        "meaningful_analytic_integrals": True,
        "method": "independent signed cell sources with exact fixed-x=2 endpoint Green primitives",
        "scope": "bookkeeping fixture only; no receipt or physical-field certification",
    }


run = run_fixture


__all__ = ["run_fixture", "run"]
