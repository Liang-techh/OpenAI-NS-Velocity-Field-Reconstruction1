"""Independent synthetic fixture for ``core_firstword_z2_profile_v1``.

The fixture creates saved-cell-shaped ``feedback_Z2`` inputs directly.  It
does not import the producer or read a receipt.  Coefficients are signed
directed intervals, every source cell has all six already-formed Leibniz
branch labels, and the checks exercise two-source causality, the axis,
triangular diagonal, all six fields, provenance/error retention, and the
distinct derivative scales.
"""

from __future__ import annotations

import traceback
from collections.abc import Mapping, Sequence
from fractions import Fraction
from typing import Any


SCHEMA = "core_firstword_z2_profile_fixture_v1"
INPUT_SCHEMA = "n1_core_word1_z2_v2"
FIELDS = ("F1", "Uz1", "K1", "P1", "d_x_F1", "d_x_Uz1")
SEAM_FIELDS = ("F1", "d_x_F1")
ROWS = {"F1": "5", "Uz1": "6", "K1": "6", "P1": "4", "d_x_F1": "5", "d_x_Uz1": "6"}
B1_SEAM_KEYS = ("5,1", "5,2", "5,3")
BRANCHES = (
    ("C0_ZZ W0", 1),
    ("2 C0_Z W0_Z", 2),
    ("C0 W0_ZZ", 1),
    ("B1_ZZ W0_Z", 1),
    ("2 B1_Z W0_ZZ", 2),
    ("B1 W0_ZZZ", 1),
)

__all__ = [
    "SCHEMA",
    "INPUT_SCHEMA",
    "FIELDS",
    "SEAM_FIELDS",
    "ROWS",
    "B1_SEAM_KEYS",
    "BRANCHES",
    "synthetic_cells",
    "synthetic_matrix_cells",
    "run_fixture",
    "run",
]


def _point(c: Any, value: int | Fraction) -> Any:
    rational = Fraction(value)
    return c.mpf(rational.numerator) / rational.denominator


def _cell_t(c: Any, index: int) -> Any:
    left = _point(c, Fraction(index, 8))
    right = _point(c, Fraction(index + 1, 8))
    return c.make_mpf((left._mpi_[0], right._mpi_[1]))


def _term(
    c: Any,
    coefficient: int,
    source_cell: int,
    row: str,
    branch_index: int,
    origin: str,
    multiplier: int,
) -> dict[str, Any]:
    input_id = f"fixture:feedback_Z2:source:{source_cell}:row:{row}:branch:{branch_index}"
    # Use a genuinely directed signed interval.  The sign alternates by
    # source/row/branch, while the interval width is never collapsed.
    if coefficient < 0:
        directed = c.mpf([coefficient - 1, coefficient])
    else:
        directed = c.mpf([coefficient, coefficient + 1])
    errors = [f"fixture:upstream-error:source:{source_cell}"]
    if branch_index % 2 == 0:
        errors.append("fixture:nonlinear-tail")
    return {
        "F0_power": (source_cell + branch_index) % 3 - 1,
        "Lambda_power": -1 - ((row == "6") + branch_index % 2),
        "coefficient": directed,
        "input_id": input_id,
        "source_id": input_id,
        "source_ids": [input_id, f"fixture:source-cell:{source_cell}"],
        "provenance": {
            "input_id": input_id,
            "source_cell": source_cell,
            "row": row,
            "branch": origin,
            "leibniz_multiplier": multiplier,
            "error_classes": list(errors),
        },
        "source_cell": source_cell,
        "source_cell_index": source_cell,
        "origin": origin,
        "leibniz_multiplier": multiplier,
        "error_classes": list(errors),
        "source_error_classes": list(errors),
        "uncertainty": True,
        "signed_field": True,
        "formal_atoms_materialized": False,
    }


def synthetic_cells(c: Any) -> list[dict[str, Any]]:
    """Return eight independent saved-cell-shaped source packets."""

    cells: list[dict[str, Any]] = []
    for source_cell in range(8):
        t = _cell_t(c, source_cell)
        rows: dict[str, list[dict[str, Any]]] = {"4": [], "5": [], "6": []}
        for row in rows:
            for branch_index, (origin, multiplier) in enumerate(BRANCHES):
                sign = -1 if (source_cell + int(row) + branch_index) % 2 else 1
                magnitude = (source_cell + 2) * (branch_index + 2) + int(row)
                rows[row].append(
                    _term(
                        c,
                        sign * magnitude,
                        source_cell,
                        row,
                        branch_index,
                        origin,
                        multiplier,
                    )
                )
        cells.append(
            {
                "index": source_cell,
                "t": t,
                "rho": 4 * t * t,
                "feedback_Z2": rows,
                "source_provenance": {
                    "input_id": f"fixture:source-cell:{source_cell}",
                    "error_classes": [f"fixture:cell-error:{source_cell}"],
                },
            }
        )
    return cells


def synthetic_matrix_cells(c: Any) -> list[dict[str, Any]]:
    """Return signed raw B1[0] seam cells for the optional missing branch."""

    cells: list[dict[str, Any]] = []
    for source_cell in range(8):
        t = _cell_t(c, source_cell)
        b1: dict[str, dict[str, list[dict[str, Any]]]] = {"0": {}}
        for key_index, key in enumerate(B1_SEAM_KEYS):
            sign = -1 if (source_cell + key_index) % 2 else 1
            value = sign * (source_cell + key_index + 7)
            input_id = f"fixture:matrix:B1:cell:{source_cell}:0:{key}:0"
            b1["0"][key] = [
                {
                    "F0_power": (source_cell + key_index) % 2,
                    "Lambda_power": -1 - key_index,
                    "coefficient": c.mpf([value, value + 1])
                    if value > 0
                    else c.mpf([value - 1, value]),
                    "input_id": input_id,
                    "source_id": input_id,
                    "source_ids": [input_id],
                    "origin": "fixture raw B1[0]",
                    "source_cell": source_cell,
                    "source_cell_index": source_cell,
                    "error_classes": [f"fixture:matrix-error:{source_cell}"],
                    "source_error_classes": [f"fixture:matrix-error:{source_cell}"],
                    "uncertainty": True,
                    "signed_field": True,
                    "formal_atoms_materialized": False,
                }
            ]
        cells.append(
            {
                "index": source_cell,
                "t": t,
                "rho": 4 * t * t,
                "B1": b1,
            }
        )
    return cells


def _inside(interval: Any, rational: Fraction) -> bool:
    numerator = rational.numerator
    denominator = rational.denominator
    point = _inside.context.mpf(numerator) / denominator
    return bool(interval.a <= point.a and point.b <= interval.b)


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


def run_fixture(build: Any, c: Any) -> dict[str, Any]:
    """Run independent shape, causality, kernel and provenance checks."""

    _inside.context = c
    build_fn = getattr(build, "build", build)
    seam_fn = getattr(build, "build_missing_seam", None)
    if not callable(build_fn):
        raise TypeError("build must be callable or expose build")
    cells = synthetic_cells(c)
    try:
        result = build_fn(c, cells)
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
    target_cells = result.get("target_cells", [])
    checks["eight_target_cells"] = len(target_cells) == 8
    checks["finite_K0_package"] = (
        result.get("finite_K0_condition_package") is True
        and result.get("complete_firstword_Z2") is False
        and result.get("physical_field_claimed") is False
    )
    checks["six_fields"] = all(
        set(cell.get("W1_Z2", {})) == set(FIELDS) for cell in target_cells
    )
    checks["strict_causality"] = all(
        cell.get("source_indices") == list(range(index + 1))
        and cell.get("strict_causality") is True
        and cell.get("source_j_gt_target_i_excluded") is True
        for index, cell in enumerate(target_cells)
    )
    checks["two_source_target"] = (
        set(target_cells[1].get("W1_Z2", {}).get("F1", [])[0].get("source_cell", -1) for _ in [0])
        <= {0, 1}
        and {term.get("source_cell") for term in target_cells[1]["W1_Z2"]["F1"]} == {0, 1}
    ) if len(target_cells) > 1 else False
    checks["axis_continuation"] = (
        len(target_cells) == 8
        and target_cells[0].get("axis_endpoint_continuation") is True
        and result.get("axis_log_division_avoided") is True
    )
    checks["triangular_diagonal"] = all(
        cell.get("triangular_diagonal_source") == index
        and cell.get("kernel_masses", {}).get(str(index)) is not None
        for index, cell in enumerate(target_cells)
    )
    terms = [
        term
        for cell in target_cells
        for field in FIELDS
        for term in cell.get("W1_Z2", {}).get(field, [])
    ]
    checks["signed_coefficients"] = any(term["coefficient"].a < 0 for term in terms) and any(
        term["coefficient"].b > 0 for term in terms
    )
    checks["provenance_and_ids"] = all(
        term.get("input_id")
        and term.get("source_term_id")
        and term.get("provenance")
        and isinstance(term.get("F0_power"), int)
        and isinstance(term.get("Lambda_power"), int)
        for term in terms
    )
    checks["errors_retained"] = any(term.get("error_classes") for term in terms) and all(
        "error_classes" in term for term in terms
    )
    checks["no_cancellation"] = (
        result.get("no_cancellation") is True
        and result.get("cancellation_performed") is False
        and all(term.get("no_cancellation") is True for term in terms)
    )
    checks["formal_atoms_not_materialized"] = (
        result.get("formal_atoms_materialized") is False
        and result.get("F0_materialized") is False
        and result.get("Lambda_materialized") is False
    )
    checks["six_input_branches_retained"] = {
        term.get("origin")
        for term in target_cells[1]["W1_Z2"]["F1"]
    } >= {branch[0] for branch in BRANCHES} if len(target_cells) > 1 else False

    # Independent rational checks for all six exact kernel scales at the
    # fixed non-axis source cell [1/8, 2/8] and target T=2/8.
    a = Fraction(1, 8)
    b = Fraction(2, 8)
    T = b
    expected = {
        "F1": (b * b - a * a) - (b**4 - a**4) / (2 * T * T),
        "Uz1": None,  # logarithmic scale is checked by positivity below.
        "K1": -((b * b - a * a) - (b**4 - a**4) / (2 * T * T)),
        "P1": 2 * (b - a),
        "d_x_F1": (b**4 - a**4) / (2 * T * T * T),
        "d_x_Uz1": (b * b - a * a) / T,
    }
    masses = target_cells[1].get("kernel_masses", {}).get("1", {}) if len(target_cells) > 1 else {}
    checks["F1_kernel_exact"] = _inside(masses.get("F1"), expected["F1"])
    checks["K1_kernel_exact"] = _inside(masses.get("K1"), expected["K1"])
    checks["P1_kernel_exact"] = _inside(masses.get("P1"), expected["P1"])
    checks["d_x_F1_kernel_exact"] = _inside(masses.get("d_x_F1"), expected["d_x_F1"])
    checks["d_x_Uz1_kernel_exact"] = _inside(masses.get("d_x_Uz1"), expected["d_x_Uz1"])
    checks["Uz1_positive"] = masses.get("Uz1").b > 0
    checks["derivative_scales_distinct"] = (
        masses.get("d_x_F1")._mpi_ != masses.get("F1")._mpi_
        and masses.get("d_x_Uz1")._mpi_ != masses.get("Uz1")._mpi_
    )

    wrong_cells = _clone(cells)
    wrong_cells[3]["t"] = _cell_t(c, 4)
    checks["wrong_cell_rejected"] = _rejected(lambda: build_fn(c, wrong_cells))
    missing_feedback = _clone(cells)
    del missing_feedback[2]["feedback_Z2"]["5"]
    checks["missing_row_rejected"] = _rejected(lambda: build_fn(c, missing_feedback))

    seam_summary: dict[str, Any] | None = None
    if callable(seam_fn):
        seam_result = seam_fn(c, result, synthetic_matrix_cells(c))
        seam_terms = seam_result.get("terms_by_order", {}).get("Z1", {})
        seam_flat = [term for field in SEAM_FIELDS for term in seam_terms.get(field, [])]
        seam_kernels = seam_result.get("kernels", {})
        checks["missing_seam_completed"] = isinstance(seam_result, Mapping)
        checks["missing_seam_only_F_and_Fx"] = (
            set(seam_terms) == set(SEAM_FIELDS)
            and set(seam_result.get("output_fields", [])) == set(SEAM_FIELDS)
            and not (set(seam_terms) - set(SEAM_FIELDS))
        )
        checks["missing_seam_marked"] = (
            seam_result.get("missing_branch_included") is True
            and seam_result.get("missing_product_rule") == "B1[0]*W1_Z2"
            and seam_result.get("known_Z1_records_untouched") is True
            and seam_result.get("old_records_merged") is False
        )
        checks["missing_seam_three_columns"] = seam_result.get("B1_zero_columns") == [1, 2, 3] and all(
            term.get("matrix_key") in B1_SEAM_KEYS and term.get("target_field") in SEAM_FIELDS
            for term in seam_flat
        )
        checks["missing_seam_sources"] = {term.get("source_cell") for term in seam_flat} >= {0, 1}
        checks["missing_seam_unique_term_ids"] = len({term.get("term_id") for term in seam_flat}) == len(seam_flat)
        checks["missing_seam_refs_and_errors"] = all(
            term.get("matrix_input_id")
            and term.get("state_input_id")
            and term.get("state_profile_term_id")
            and term.get("source_ids")
            and term.get("error_classes")
            and term.get("state_order") == "2"
            and term.get("ordinary_Z_order") == 1
            for term in seam_flat
        )
        checks["missing_seam_no_cancellation"] = (
            seam_result.get("no_cancellation") is True
            and seam_result.get("cancellation_performed") is False
            and all(term.get("no_cancellation") is True for term in seam_flat)
        )
        checks["missing_seam_finite_K0"] = (
            seam_result.get("finite_K0_condition_package") is True
            and seam_result.get("complete_secondword") is False
            and seam_result.get("physical_field_claimed") is False
        )
        checks["missing_seam_F_Fx_flags"] = (
            seam_result.get("F_and_Fx_exact") is True
            and all(item.get("exact") is True for item in seam_kernels.values())
            and all(item.get("F1") is not None and item.get("Fx") is not None for item in seam_kernels.values())
        )
        seam_summary = {
            "schema": seam_result.get("schema"),
            "term_count": len(seam_flat),
            "columns": seam_result.get("B1_zero_columns"),
            "missing_product_rule": seam_result.get("missing_product_rule"),
        }

    return {
        "schema": SCHEMA,
        "passed": all(checks.values()),
        "checks": checks,
        "check_count": len(checks),
        "synthetic_only": True,
        "implementation_imported": False,
        "seam": seam_summary,
        "strict_j_less_than_i_checked": True,
        "axis_and_triangular_checked": True,
        "derivative_scales_checked": True,
    }


def run(build: Any, c: Any) -> dict[str, Any]:
    return run_fixture(build, c)

