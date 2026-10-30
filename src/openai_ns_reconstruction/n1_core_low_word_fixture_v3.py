"""Independent synthetic checks for the low core first-word Green leaf.

The fixture supplies eight deliberately different outer cells.  Their W0
columns and the C0/B1 matrix entries are formal source sectors with signed
directed coefficients.  The expected source is made by an independent
polynomial-coefficient convolution: ``C0*W0 + B1*d_Z(W0)``.  The check then
applies the six Green kernels to ordinary Z0 and Z1 and compares every
unmerged product, not only a merged total.

This is a synthetic bookkeeping test.  It does not import a producer, read a
saved receipt, select a point from an interval, or certify a physical field.
"""
from __future__ import annotations

from collections import defaultdict
from fractions import Fraction
from math import factorial


FIELDS = (
    ("F1", "5"),
    ("Uz1", "6"),
    ("K1", "6"),
    ("P1", "4"),
    ("d_x_F1", "5"),
    ("d_x_Uz1", "6"),
)

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


def _rational(c, value):
    """Create a zero-width directed interval from an exact Fraction."""

    return c.mpf(value.numerator) / value.denominator


def _box(c, lower, upper):
    """Create a directed interval with exact rational endpoints."""

    left = _rational(c, lower)
    right = _rational(c, upper)
    return c.make_mpf((left._mpi_[0], right._mpi_[1]))


def _signed_box(c, numerator):
    """Return a narrow nonzero box whose sign follows ``numerator``."""

    center = Fraction(numerator, 19)
    width = Fraction(1, 1900)
    return _box(c, center - width, center + width)


def _outer_t(c, index):
    left = _rational(c, Fraction(index, 8))
    right = _rational(c, Fraction(index + 1, 8))
    return c.make_mpf((left._mpi_[0], right._mpi_[1]))


def _w0_terms(c, cell, order, column):
    """Make one or two signed formal sectors for one saved W0 column."""

    sign = -1 if (5 * cell + 2 * order + column) % 4 in (0, 3) else 1
    magnitude = 2 + ((cell + 1) * (order + 2) * (column + 1)) % 8
    f0_power = (cell + 2 * order + column) % 3
    lambda_power = -1 - ((2 * cell + order + column) % 3)
    terms = [
        {
            "F0_power": f0_power,
            "Lambda_power": lambda_power,
            "coefficient": _signed_box(c, sign * magnitude),
            "source_cell": 100000 + 1000 * cell + 100 * order + 2 * column,
        }
    ]
    # A few columns carry a second sector, so the convolution has genuine
    # cross-products rather than one-term scalar multiplication everywhere.
    if (cell + order + column) % 5 == 0:
        second_sign = -sign
        terms.append(
            {
                "F0_power": (f0_power + 1) % 4,
                "Lambda_power": -1 - ((lambda_power * -1 + cell + column) % 4),
                "coefficient": _signed_box(c, second_sign * (1 + (column % 4))),
                "source_cell": 100000 + 1000 * cell + 100 * order + 2 * column + 1,
            }
        )
    return terms


def _matrix_terms(c, cell, name, order, key_index):
    """Make signed C0/B1 terms with independent source-cell identities."""

    family_shift = 0 if name == "C0" else 1
    sign = -1 if (cell + 3 * order + key_index + family_shift) % 3 == 0 else 1
    magnitude = 1 + ((2 * cell + order + 2 * key_index + family_shift) % 9)
    f0_power = (cell + order + key_index + family_shift) % 3
    lambda_power = -1 - ((cell + 2 * order + key_index) % 3)
    source_base = (
        200000
        + (0 if name == "C0" else 100000)
        + 10000 * cell
        + 1000 * order
        + 10 * key_index
    )
    terms = [
        {
            "F0_power": f0_power,
            "Lambda_power": lambda_power,
            "coefficient": _signed_box(c, sign * magnitude),
            "source_cell": source_base,
        }
    ]
    if (cell + order + key_index + family_shift) % 6 == 0:
        terms.append(
            {
                "F0_power": (f0_power + 2) % 4,
                "Lambda_power": -2 - ((cell + key_index) % 3),
                "coefficient": _signed_box(c, -sign * (2 + (key_index % 3))),
                "source_cell": source_base + 1,
            }
        )
    return terms


def _inputs(c):
    """Build eight coordinate-matched outer and matrix cells."""

    outer_cells = []
    matrix_cells = []
    for cell in range(8):
        t = _outer_t(c, cell)
        rho = 4 * t * t
        w0 = {
            str(order): {
                str(column): _w0_terms(c, cell, order, column)
                for column in range(1, 7)
            }
            for order in range(3)
        }
        outer_cells.append({"t": t, "rho": rho, "W0_ordinary_Z": w0})
        matrix_cells.append(
            {
                "t": t,
                "rho": rho,
                "C0": {
                    str(order): {
                        key: _matrix_terms(c, cell, "C0", order, key_index)
                        for key_index, key in enumerate(C0_KEYS)
                    }
                    for order in range(2)
                },
                "B1": {
                    str(order): {
                        key: _matrix_terms(c, cell, "B1", order, key_index)
                        for key_index, key in enumerate(B1_KEYS)
                    }
                    for order in range(2)
                },
            }
        )
    return outer_cells, matrix_cells


def _w0_order(order):
    if order == 0:
        return "W0"
    if order == 1:
        return "W0_Z"
    return "W0_ZZ"


def _matrix_order(name, order):
    return name if order == 0 else f"{name}_Z"


def _expected_sources(outer, matrix, order):
    """Convolve polynomial coefficients without using the leaf's specs."""

    result = {"4": [], "5": [], "6": []}
    # The source polynomial is C0(Z) W0(Z) + B1(Z) dW0(Z).  Expanding this
    # with coefficient arrays gives the requested ordinary Z0/Z1 rows.
    for name, keys in (("C0", C0_KEYS), ("B1", B1_KEYS)):
        for key in keys:
            row, column = (int(value) for value in key.split(","))
            key_index = keys.index(key)
            for matrix_order in range(order + 1):
                matrix_terms = matrix[name][str(matrix_order)][key]
                if name == "C0":
                    w0_order = order - matrix_order
                    derivative_factor = 1
                else:
                    # d_Z W0 has coefficient (q+1) W0[q+1].
                    w0_order = order - matrix_order + 1
                    derivative_factor = w0_order
                if w0_order > 2:
                    continue
                multiplier = Fraction(factorial(order)*derivative_factor, factorial(matrix_order)*factorial(w0_order))
                origin = f"{_matrix_order(name, matrix_order)} {_w0_order(w0_order)}"
                for matrix_term in matrix_terms:
                    for w0_term in outer["W0_ordinary_Z"][str(w0_order)][str(column)]:
                        result[str(row)].append(
                            {
                                "F0_power": int(matrix_term["F0_power"])
                                + int(w0_term["F0_power"]),
                                "Lambda_power": int(matrix_term["Lambda_power"])
                                + int(w0_term["Lambda_power"]),
                                "base_coefficient": (
                                    multiplier.numerator
                                    * matrix_term["coefficient"]
                                    * w0_term["coefficient"]
                                    / multiplier.denominator
                                ),
                                "origin": origin,
                                "row": row,
                                "column": column,
                                "matrix_row": row,
                                "matrix_column": column,
                                "outer_cell": None,
                                "matrix_source_cell": matrix_term["source_cell"],
                                "W0_source_cell": w0_term["source_cell"],
                                "source_cell": matrix_term["source_cell"],
                                "matrix_key_index": key_index,
                            }
                        )
    return result


def _kernel(c, index, t, field):
    """Return one signed whole-cell Green kernel including the leaf factor."""

    tlog = c.mpf([0, "0.5"]) if index == 0 else t * c.log(1 / t)
    K = t * (1 - t * t)
    raw = {
        "F1": K,
        "Uz1": 2 * tlog,
        "K1": -K,
        "P1": c.mpf(1),
        "d_x_F1": t**3,
        "d_x_Uz1": t,
    }
    return raw[field] / 4


def _integral_kernel(c, index, field):
    """Exact integral of a synthetic cell's polynomial kernel over its t-cell."""

    left = _rational(c, Fraction(index, 8))
    right = _rational(c, Fraction(index + 1, 8))

    def f1_antiderivative(x):
        return x * x / 8 - x**4 / 16

    def uz1_antiderivative(x, zero=False):
        if zero:
            return c.mpf(0)
        return x * x * c.log(1 / x) / 4 + x * x / 8

    if field == "F1":
        return f1_antiderivative(right) - f1_antiderivative(left)
    if field == "K1":
        return -(f1_antiderivative(right) - f1_antiderivative(left))
    if field == "Uz1":
        return uz1_antiderivative(right) - uz1_antiderivative(left, index == 0)
    if field == "P1":
        return (right - left) / 4
    if field == "d_x_F1":
        return (right**4 - left**4) / 16
    if field == "d_x_Uz1":
        return (right * right - left * left) / 8
    raise KeyError(field)


def _signature(term):
    return (
        term.get("outer_cell"),
        term.get("origin"),
        term.get("row"),
        term.get("column"),
        term.get("matrix_row"),
        term.get("matrix_column"),
        term.get("matrix_source_cell"),
        term.get("W0_source_cell"),
        term.get("source_cell"),
        int(term.get("F0_power")),
        int(term.get("Lambda_power")),
    )


def _contains(container, value):
    """Check interval inclusion without selecting a midpoint."""

    return container.a <= value.a and container.b >= value.b


def _actual_index(terms):
    indexed = defaultdict(list)
    for term in terms:
        indexed[_signature(term)].append(term)
    return indexed


def _expected_terms(c, outer_cells, matrix_cells, order, field, row):
    expected = []
    for index, (outer, matrix) in enumerate(zip(outer_cells, matrix_cells)):
        sources = _expected_sources(outer, matrix, order)
        factor = _kernel(c, index, outer["t"], field)
        for source in sources[row]:
            term = dict(source)
            term["outer_cell"] = index
            term["coefficient"] = source["base_coefficient"] * factor
            del term["base_coefficient"]
            expected.append(term)
    return expected


def _continuous_total(c, outer_cells, matrix_cells, order, field, row):
    total = c.mpf(0)
    for index, (outer, matrix) in enumerate(zip(outer_cells, matrix_cells)):
        sources = _expected_sources(outer, matrix, order)
        factor = _integral_kernel(c, index, field)
        for source in sources[row]:
            total += source["base_coefficient"] * factor
    return total


def _check_family(c, result, outer_cells, matrix_cells, order, field, row, checks):
    key = f"W1_Z{order}_core_exit"
    actual = result[key][field]
    expected = _expected_terms(c, outer_cells, matrix_cells, order, field, row)
    actual_index = _actual_index(actual)
    expected_index = _actual_index(expected)
    label = f"Z{order}_{field}"

    checks[f"{label}_term_count"] = len(actual) == len(expected)
    checks[f"{label}_eight_cells"] = {
        item.get("outer_cell") for item in actual
    } == set(range(8))
    checks[f"{label}_no_duplicate_sectors"] = (
        len(actual_index) == len(actual)
        and sorted(actual_index) == sorted(expected_index)
    )
    checks[f"{label}_source_powers"] = set(actual_index) == set(expected_index)

    coefficient_ok = True
    for signature, expected_items in expected_index.items():
        actual_items = actual_index.get(signature, [])
        if len(actual_items) != len(expected_items):
            coefficient_ok = False
            continue
        for got, want in zip(actual_items, expected_items):
            coefficient_ok = coefficient_ok and _contains(
                got["coefficient"], want["coefficient"]
            )
    checks[f"{label}_term_coefficients"] = coefficient_ok

    # The source cell has width 1/8 in t and dx=2 dt.  The leaf's
    # range kernel times 1/4 therefore encloses 2 times the t integral.
    # Check every unmerged sector individually, avoiding fragile comparisons
    # of interval sums made in different rounding orders.
    integrated = []
    for index, (outer, matrix) in enumerate(zip(outer_cells, matrix_cells)):
        for source in _expected_sources(outer, matrix, order)[row]:
            term = dict(source)
            term["outer_cell"] = index
            term["coefficient"] = source["base_coefficient"] * (8 * _integral_kernel(c, index, field))
            del term["base_coefficient"]
            integrated.append(term)
    integrated_index = _actual_index(integrated)
    valid = set(actual_index) == set(integrated_index)
    for signature, expected_items in integrated_index.items():
        actual_items = actual_index.get(signature, [])
        valid = valid and len(actual_items) == len(expected_items)
        for got, want in zip(actual_items, expected_items):
            valid = valid and _contains(got["coefficient"], want["coefficient"])
    checks[f"{label}_polynomial_integral"] = valid



def run(leaf, c):
    """Run the independent synthetic low-word checks against ``leaf``."""

    outer_cells, matrix_cells = _inputs(c)
    result = leaf.build(c, outer_cells, matrix_cells)
    checks = {}

    checks["eight_outer_cells"] = len(outer_cells) == 8
    checks["eight_matrix_cells"] = len(matrix_cells) == 8
    checks["coordinate_cells_match"] = all(
        outer["t"]._mpi_ == matrix["t"]._mpi_
        and outer["rho"]._mpi_ == matrix["rho"]._mpi_
        for outer, matrix in zip(outer_cells, matrix_cells)
    )
    checks["synthetic_W0_has_both_signs"] = any(
        term["coefficient"].b < 0
        for outer in outer_cells
        for order in outer["W0_ordinary_Z"].values()
        for terms in order.values()
        for term in terms
    ) and any(
        term["coefficient"].a > 0
        for outer in outer_cells
        for order in outer["W0_ordinary_Z"].values()
        for terms in order.values()
        for term in terms
    )
    checks["synthetic_matrix_has_both_signs"] = any(
        term["coefficient"].b < 0
        for matrix in matrix_cells
        for name in ("C0", "B1")
        for order in matrix[name].values()
        for terms in order.values()
        for term in terms
    ) and any(
        term["coefficient"].a > 0
        for matrix in matrix_cells
        for name in ("C0", "B1")
        for order in matrix[name].values()
        for terms in order.values()
        for term in terms
    )

    for order in (0, 1):
        for field, row in FIELDS:
            _check_family(c, result, outer_cells, matrix_cells, order, field, row, checks)

    checks["signed_K_kernel"] = all(
        _kernel(c, index, outer_cells[index]["t"], "K1")
        == -_kernel(c, index, outer_cells[index]["t"], "F1")
        for index in range(8)
    )
    checks["ordinary_input_order"] = result.get("input_W0_ordinary_Z_orders") == [0, 1, 2]
    checks["ordinary_output_order"] = result.get("output_Z_order") == 1
    checks["W0_not_reweighted"] = result.get("W0_reweighted") is False
    checks["source_cells_not_collapsed"] = result.get("no_source_cells_collapsed") is True
    checks["formal_exponentials_kept"] = result.get("source_exponentials_materialized") is False

    return {
        "passed": all(checks.values()),
        "checks": checks,
        "check_count": len(checks),
        "synthetic_only": True,
        "fixture_is_not_actual_field": True,
        "method": "independent ordinary-Z polynomial convolution and directed whole-cell six-kernel inclusion",
        "certification_scope": "synthetic bookkeeping only; no real-point or physical-field certification",
    }


__all__ = ["run"]
