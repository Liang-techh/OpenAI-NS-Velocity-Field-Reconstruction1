"""Self-contained synthetic checks for ``n1_kchi_inner_feedback_v1``.

The fixture supplies eight variable source cells, sparse DeltaW0 columns 1 and
5, complete saved C0/B1 support at matrix orders 0 and 1, and a separate
symmetric uncertainty packet.  It checks source ordering, the triangular
diagonal and axis continuation, ordinary Z1 Leibniz terms, four formal powers,
source/error provenance (including an explicit ``None`` source ID), and the
fact that uncertainty sectors stay outside the signed finite field.

This file only defines synthetic data and ``run_fixture``.  It does not import
an old producer, read an archive, or execute a fixture at import time.
"""
from __future__ import annotations

from fractions import Fraction
from typing import Any


C0_KEYS = (
    "4,1", "5,1", "5,2", "5,3", "5,5",
    "6,1", "6,2", "6,3", "6,4", "6,6",
)
B1_KEYS = ("5,1", "5,2", "5,3", "6,2", "6,3", "6,4")


def _q(c: Any, value: Fraction) -> Any:
    return c.mpf(value.numerator) / value.denominator


def _box(c: Any, lower: Fraction, upper: Fraction) -> Any:
    lo = _q(c, lower)
    hi = _q(c, upper)
    return c.make_mpf((lo._mpi_[0], hi._mpi_[1]))


def _cell_x(c: Any, index: int) -> Any:
    return _box(c, Fraction(index, 4), Fraction(index + 1, 4))


def _cell_t(c: Any, index: int) -> Any:
    return _box(c, Fraction(index, 8), Fraction(index + 1, 8))


def _term(c: Any, coefficient: Any, source_id: Any, *, powers=(1, 2, 0, -1),
          error_class="synthetic_signed") -> dict[str, Any]:
    return {
        "powers": list(powers),
        "coefficient": coefficient if hasattr(coefficient, "_mpi_") else c.mpf(coefficient),
        "source_id": source_id,
        "origin": {
            "source_ids": [source_id],
            "error_classes": [error_class],
            "uncertainty": False,
        },
    }


def _tail_term(c: Any, width: Fraction, source_id: Any, *, order: int) -> dict[str, Any]:
    return {
        "powers": [1, 2, 0, -1],
        "coefficient": _box(c, -width, width),
        "source_id": source_id,
        "origin": {
            "source_ids": [source_id],
            "error_classes": ["synthetic_kernel_tail"],
            "uncertainty": True,
            "operation": "symmetric_kernel_tail",
            "order": order,
        },
        "uncertainty": True,
    }


def _delta_cells(c: Any) -> list[dict[str, Any]]:
    cells: list[dict[str, Any]] = []
    for index in range(8):
        # The source varies with the cell and with ordinary Z order; this is
        # deliberately not a constant source that could hide a Leibniz error.
        lo = Fraction(1 + index, 10)
        hi = Fraction(2 + index, 10)
        f1 = {
            "Z0": [_term(c, _box(c, lo, hi), None if index == 0 else f"dF:{index}:0")],
            "Z1": [_term(c, _box(c, lo + Fraction(1, 10), hi + Fraction(1, 10)), f"dF:{index}:1")],
            "Z2": [_term(c, _box(c, lo + Fraction(2, 10), hi + Fraction(2, 10)), f"dF:{index}:2")],
        }
        dxf1 = {
            "Z0": [_term(c, _box(c, lo + Fraction(1, 20), hi + Fraction(1, 20)), f"dFx:{index}:0")],
            "Z1": [_term(c, _box(c, lo + Fraction(2, 20), hi + Fraction(2, 20)), f"dFx:{index}:1")],
            "Z2": [_term(c, _box(c, lo + Fraction(3, 20), hi + Fraction(3, 20)), f"dFx:{index}:2")],
        }
        x = _cell_x(c, index)
        cells.append({
            "x": x,
            "rho": x * x,
            "W0_correction": {"F1": f1, "d_x_F1": dxf1},
            # These are intentionally empty unsupported columns.  The leaf
            # consumes only DeltaW0 columns 1 and 5.
            "unused_columns": {"2": [], "3": [], "4": [], "6": []},
        })
    return cells


def _matrix_term(c: Any, index: int, name: str, order: int, key: str) -> list[dict[str, Any]]:
    row, column = (int(value) for value in key.split(","))
    key_index = (C0_KEYS if name == "C0" else B1_KEYS).index(key)
    # Leave most rows empty while preserving the complete saved support map.
    if column not in (1, 5) or row not in (5, 6):
        return []
    sign = -1 if (index + order + key_index + (name == "B1")) % 3 == 0 else 1
    magnitude = Fraction(1 + (index + key_index + order) % 4, 20)
    coefficient = _box(c, magnitude - Fraction(1, 1000), magnitude + Fraction(1, 1000))
    if sign < 0:
        coefficient = -coefficient
    return [{
        "powers": [0, order, 0, -1 - order],
        "coefficient": coefficient,
        "source_id": f"{name}:{index}:{order}:{key}",
        "origin": {
            "source_ids": [f"{name}:{index}:{order}:{key}"],
            "error_classes": ["synthetic_matrix"],
            "uncertainty": False,
        },
    }]


def _matrix_cells(c: Any) -> list[dict[str, Any]]:
    cells: list[dict[str, Any]] = []
    for index in range(8):
        cells.append({
            "t": _cell_t(c, index),
            "rho": 4 * _cell_t(c, index) * _cell_t(c, index),
            "C0": {
                str(order): {
                    key: _matrix_term(c, index, "C0", order, key)
                    for key in C0_KEYS
                }
                for order in (0, 1)
            },
            "B1": {
                str(order): {
                    key: _matrix_term(c, index, "B1", order, key)
                    for key in B1_KEYS
                }
                for order in (0, 1)
            },
        })
    return cells


def _uncertainty_cells(c: Any) -> list[dict[str, Any]]:
    cells: list[dict[str, Any]] = []
    for index in range(8):
        width = Fraction(1 + index, 1000)
        cells.append({
            "x": _cell_x(c, index),
            "rho": _cell_x(c, index) * _cell_x(c, index),
            "uncertainty": {
                "F1": {
                    "Z0": [_tail_term(c, width, None, order=0)],
                    "Z1": [_tail_term(c, width + Fraction(1, 1000), f"tailF:{index}:1", order=1)],
                    "Z2": [_tail_term(c, width + Fraction(2, 1000), f"tailF:{index}:2", order=2)],
                },
                "d_x_F1": {
                    "Z0": [_tail_term(c, width, f"tailFx:{index}:0", order=0)],
                    "Z1": [_tail_term(c, width + Fraction(1, 1000), f"tailFx:{index}:1", order=1)],
                    "Z2": [_tail_term(c, width + Fraction(2, 1000), f"tailFx:{index}:2", order=2)],
                },
            },
        })
    return cells


def _endpoint(bits: Any) -> Fraction:
    sign, mantissa, exponent, _ = bits
    value = Fraction(int(mantissa), 1)
    value = value * (1 << exponent) if exponent >= 0 else value / (1 << (-exponent))
    return -value if sign else value


def _contains(interval: Any, expected: Fraction) -> bool:
    return _endpoint(interval._mpi_[0]) <= expected <= _endpoint(interval._mpi_[1])


def run_fixture(build: Any, c: Any) -> dict[str, Any]:
    """Run synthetic bookkeeping checks against a supplied build callable."""

    delta = _delta_cells(c)
    matrices = _matrix_cells(c)
    uncertainty = _uncertainty_cells(c)
    result = build(c, delta, matrices, uncertainty_cells=uncertainty)
    targets = result.get("target_cells", [])
    checks: dict[str, bool] = {
        "eight_targets": len(targets) == 8,
        "source_j_le_target_i": all(
            target.get("source_indices") == list(range(index + 1))
            for index, target in enumerate(targets)
        ),
        "triangular_diagonal_source": all(
            target.get("triangular_diagonal_source") == index
            for index, target in enumerate(targets)
        ),
        "axis_continuation_marked": bool(targets)
        and targets[0].get("axis_endpoint_continuation") is True,
        "signed_branch_separate": result.get("signed_target_cells") == targets,
        "uncertainty_branch_separate": len(result.get("uncertainty_target_cells", [])) == 8,
        "uncertainty_not_field": result.get("uncertainty_output_is_not_a_signed_field") is True
        and result.get("kernel_tail_merged") is False,
        "four_powers": result.get("source_power_basis") == ["hb", "F0_anchor", "Pstar_squared", "Lambda"],
        "z1_rule": tuple(item[0:3] for item in result.get("ordinary_Z1_product_rule", ())) == (
            ("C0", "1", "0"), ("C0", "0", "1"),
            ("B1", "1", "1"), ("B1", "0", "2"),
        ),
    }
    checks["axis_f1_mass"] = bool(targets) and bool(targets[0]["G0_Z0"]["F1"])
    checks["axis_uz1_mass"] = bool(targets) and bool(targets[0]["G0_Z0"]["Uz1"])

    for target in targets:
        for order in ("Z0", "Z1"):
            for field in ("F1", "Uz1"):
                terms = target["G0_feedback"][order][field]
                checks[f"{order}_{field}_directed"] = bool(terms) and all(
                    hasattr(term["coefficient"], "_mpi_") for term in terms
                )
                checks[f"{order}_{field}_powers"] = bool(terms) and all(
                    len(term["powers"]) == 4 and all(isinstance(item, int) for item in term["powers"])
                    for term in terms
                )
                checks[f"{order}_{field}_sparse_columns"] = all(
                    term["source_column"] in (1, 5) for term in terms
                )
                if order == "Z1":
                    checks["Z1_leibniz_middle"] = checks.get("Z1_leibniz_middle", False) or any(
                        term["origin_label"] == "C0 DeltaW0"
                        or term["origin_label"] == "B1 DeltaW0_ZZ"
                        for term in terms
                    )
    checks["z1_has_all_four_origins"] = bool(targets) and {
        term["origin_label"]
        for term in targets[1]["G0_Z1"]["F1"]
    } >= {"C0_Z DeltaW0", "C0 DeltaW0_Z", "B1_Z DeltaW0_Z", "B1 DeltaW0_ZZ"}
    checks["null_source_preserved"] = bool(targets) and any(
        None in term.get("source_ids", [])
        for term in targets[0]["G0_Z0"]["F1"]
    )
    uncertainty_terms = [
        term
        for target in result.get("uncertainty_target_cells", [])
        for order in ("Z0", "Z1")
        for field in ("F1", "Uz1")
        for term in target.get("G0_feedback", {}).get(order, {}).get(field, [])
    ]
    checks["uncertainty_terms_symmetric"] = bool(uncertainty_terms) and all(
        term.get("symmetric_uncertainty") is True
        and term.get("kernel_tail_not_a_field") is True
        and term.get("uncertainty") is True
        for term in uncertainty_terms
    )
    checks["uncertainty_classes_retained"] = bool(uncertainty_terms) and all(
        "synthetic_kernel_tail" in term.get("error_classes", [])
        for term in uncertainty_terms
    )
    checks["no_future_source"] = all(
        max(target.get("source_indices", [-1])) == index
        for index, target in enumerate(targets)
    )
    checks["finite_only"] = result.get("finite_subterm_only") is True and result.get("allword_tail_added") is False
    return {
        "passed": all(checks.values()),
        "checks": checks,
        "check_count": len(checks),
        "synthetic_only": True,
        "fixture_is_not_original_source": True,
        "method": "independent variable-cell source, sparse DeltaW0 columns, K0 axis continuation, and separate symmetric tail",
    }


run = run_fixture


__all__ = ["run_fixture", "run"]
