"""Self-contained synthetic checks for the outer Kchi feedback leaf.

The fixture creates no receipt and does not import an old producer.  Its
saved packets have nonconstant signed source coefficients, ordinary Z0/Z1/Z2
W0 sectors, an explicitly nonzero chi_Z, and a second delta-W0 branch with
only columns 1 and 5 populated.  Checks cover Z1 chi differentiation,
triangular and axis masses, source ordering, unweighted q5 products, formal
powers, source IDs and uncertainty classes.
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


def _t_cell(c: Any, index: int) -> Any:
    left = _q(c, Fraction(index, 8))
    right = _q(c, Fraction(index + 1, 8))
    return c.make_mpf((left._mpi_[0], right._mpi_[1]))


def _interval(c: Any, lo: Fraction, hi: Fraction) -> Any:
    lower = _q(c, lo)
    upper = _q(c, hi)
    return c.make_mpf((lower._mpi_[0], upper._mpi_[1]))


def _term(c: Any, coefficient: Any, source_id: str, *, powers: tuple[int, int, int, int],
          error_class: str = "synthetic_source", uncertainty: bool = False) -> dict[str, Any]:
    return {
        "powers": list(powers),
        "coefficient": coefficient if hasattr(coefficient, "_mpi_") else c.mpf(coefficient),
        "source_id": source_id,
        "source_cell": source_id,
        "origin": {
            "source_ids": [source_id],
            "error_classes": [error_class],
            "uncertainty": uncertainty,
        },
    }


def _saved_w0(c: Any, index: int, *, delta: bool = False) -> dict[str, Any]:
    packet: dict[str, Any] = {}
    for order in range(3):
        row: dict[str, Any] = {}
        for column in range(1, 7):
            if delta and column not in (1, 5):
                row[str(column)] = []
                continue
            # The coefficient changes with cell, ordinary order and column;
            # the bounds retain a genuinely nonconstant signed source.
            base = Fraction(index + order + column, 8)
            signed = -1 if (index + order + column) % 3 == 0 else 1
            width = Fraction(1, 32)
            if signed < 0:
                coefficient = _interval(c, -(base + width), -(base - width))
            else:
                coefficient = _interval(c, base - width, base + width)
            row[str(column)] = [_term(
                c, coefficient,
                f"{'dW0' if delta else 'W0'}:{index}:Z{order}:col{column}",
                powers=(0, 1 + order, column % 2, -order),
                error_class="synthetic_delta" if delta else "synthetic_saved_W0",
                uncertainty=(column == 5 and delta),
            )]
        packet[str(order)] = row
    return packet


def _inputs(c: Any) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    outer: list[dict[str, Any]] = []
    matrices: list[dict[str, Any]] = []
    delta: list[dict[str, Any]] = []
    for index in range(8):
        t = _t_cell(c, index)
        outer.append({
            "t": t,
            "rho": 4 * t * t,
            "W0_ordinary_Z": _saved_w0(c, index),
        })
        delta.append({
            "t": t,
            "rho": 4 * t * t,
            "delta_W0_ordinary_Z": _saved_w0(c, index, delta=True),
        })
        c0: dict[str, Any] = {}
        b1: dict[str, Any] = {}
        for order in range(2):
            c0_row: dict[str, Any] = {}
            b1_row: dict[str, Any] = {}
            for key_index, key in enumerate(C0_KEYS):
                sign = -1 if (index + order + key_index) % 4 == 0 else 1
                coefficient = c.mpf(sign * (index + order + key_index + 1)) / 16
                c0_row[key] = [_term(
                    c, coefficient, f"C0:{index}:Z{order}:{key}",
                    powers=(1, -order, key_index % 2, order),
                    error_class="synthetic_matrix",
                    uncertainty=(key_index == 2),
                )]
            for key_index, key in enumerate(B1_KEYS):
                sign = -1 if (index + 2 * order + key_index) % 3 == 0 else 1
                coefficient = c.mpf(sign * (index + key_index + 2)) / 24
                b1_row[key] = [_term(
                    c, coefficient, f"B1:{index}:Z{order}:{key}",
                    powers=(0, -1 - order, (key_index + 1) % 2, order),
                    error_class="synthetic_matrix",
                    uncertainty=(key_index == 1),
                )]
            c0[str(order)] = c0_row
            b1[str(order)] = b1_row
        matrices.append({"t": t, "rho": 4 * t * t, "C0": c0, "B1": b1})
    return outer, matrices, delta


def _endpoint(value: Any) -> Fraction:
    sign, mantissa, exponent, _ = value
    result = Fraction(int(mantissa), 1)
    result = result * (1 << exponent) if exponent >= 0 else result / (1 << (-exponent))
    return -result if sign else result


def _contains(interval: Any, expected: Fraction) -> bool:
    return _endpoint(interval._mpi_[0]) <= expected <= _endpoint(interval._mpi_[1])


def run_fixture(build: Any, c: Any) -> dict[str, Any]:
    """Run synthetic structural, sign and source-provenance checks."""

    outer, matrices, delta = _inputs(c)
    chi = {
        0: [_term(c, c.mpf(1) / 3, "chi:Z0",
                  powers=(0, 0, 0, 0), error_class="synthetic_chi")],
        1: [_term(c, -c.mpf(1) / 5, "chi:Z1",
                  powers=(0, 0, 0, 0), error_class="synthetic_chi")],
    }
    result = build(c, outer, matrices, chi, delta_w0_cells=delta)
    branches = result.get("branches", {})
    targets = result.get("target_cells", [])
    delta_targets = branches.get("delta_W0", {}).get("target_cells", [])
    checks: dict[str, bool] = {
        "eight_targets": len(targets) == 8,
        "eight_delta_targets": len(delta_targets) == 8,
        "dual_branches": set(branches) == {"K0_W0", "delta_W0"},
        "saved_w0_reused": result.get("saved_outer_cells_reused") is True,
        "w0_not_reconstructed": result.get("W0_reconstructed") is False,
        "source_j_le_target_i": all(
            cell.get("source_indices") == list(range(index + 1))
            for index, cell in enumerate(targets)
        ),
        "delta_source_j_le_target_i": all(
            cell.get("source_indices") == list(range(index + 1))
            for index, cell in enumerate(delta_targets)
        ),
        "triangular_source_marked": all(
            cell.get("triangular_diagonal_source") == index
            for index, cell in enumerate(targets)
        ),
        "axis_continuation_marked": bool(targets)
        and targets[0].get("axis_endpoint_continuation") is True,
        "formal_powers_retained": bool(targets)
        and all(len(term.get("powers", ())) == 4
                for term in targets[1]["outer_delta_G1"]["F1"]["Z1"]),
        "source_ids_retained": bool(targets)
        and all(term.get("source_ids") for term in targets[1]["F1"]["Z0"]),
        "outer_uses_q5_only": all(
            term.get("row") == 5
            for cell in targets
            for field in ("F1", "d_x_F1")
            for order in ("Z0", "Z1")
            for term in cell[field][order]
        ),
        "uncertainty_retained": any(
            term.get("uncertainty") for cell in targets
            for term in cell["F1"]["Z0"]
        ),
        "no_full_w1_claim": result.get("full_W1_claimed") is False,
    }

    axis_entry = targets[0]["outer_delta_G1_by_source_cell"]["0"]
    checks["axis_mass"] = _contains(
        axis_entry["kernel_masses"]["F1"], Fraction(1, 192 * 4**4)
    )
    checks["axis_dx_mass"] = _contains(
        axis_entry["kernel_masses"]["d_x_F1"], Fraction(1, 48 * 4**3)
    )
    checks["diagonal_triangular_mass"] = (
        targets[1]["outer_delta_G1_by_source_cell"]["1"]["triangular"]
        is True
    )
    checks["z1_contains_chi_z"] = any(
        "chi_Z" in term.get("origin", "")
        for term in targets[1]["F1"]["Z1"]
    )
    checks["z1_contains_source_derivative"] = any(
        any(label in term.get("origin", "")
            for label in ("C0_Z", "C0 W0_Z", "B1_Z", "B1 W0_ZZ"))
        for term in targets[1]["F1"]["Z1"]
    )
    k0_branch = branches["K0_W0"]
    checks["unweighted_q5_kept"] = bool(
        k0_branch["feedback_q5_products_by_source_cell"]["0"]["Z0"]
    )
    feedback_q5 = k0_branch["feedback_q5_products_by_source_cell"]["0"]["Z1"]
    checks["feedback_q5_is_before_chi"] = bool(feedback_q5) and all(
        "chi" not in term.get("origin", "")
        and "correction_factor" not in term
        and "kernel_mass" not in term
        for term in feedback_q5
    )
    q5_unweighted = k0_branch["unweighted_outer_delta_source_by_source_cell"]["0"]["Z1"]
    checks["unweighted_q5_has_chi_z"] = any(
        "chi_Z" in term.get("origin", "") for term in q5_unweighted
    )
    checks["unweighted_q5_has_no_green_mass"] = bool(q5_unweighted) and all(
        term.get("unweighted_outer_source") is True
        and term.get("kernel_mass") is None
        for term in q5_unweighted
    )

    delta_terms = [
        term
        for cell in delta_targets
        for term in cell["F1"]["Z0"]
    ]
    checks["delta_nonempty_columns_only"] = bool(delta_terms) and all(
        term.get("W0_column") in (1, 5) for term in delta_terms
    )
    checks["delta_signed_terms"] = any(
        term["coefficient"].a < 0 or term["coefficient"].b < 0
        for term in delta_terms
    )
    checks["delta_branch_separate"] = (
        delta_targets is not targets
        and branches["delta_W0"].get("branch") == "delta_W0"
    )
    checks["delta_uncertainty_retained"] = any(
        term.get("uncertainty") for cell in delta_targets
        for term in cell["F1"]["Z0"]
    )

    return {
        "passed": all(checks.values()),
        "checks": checks,
        "check_count": len(checks),
        "synthetic_only": True,
        "fixture_is_not_original_source": True,
        "method": "independent signed polynomial source, exact axis/triangular TK0 masses and dual-branch bookkeeping",
    }


run = run_fixture


__all__ = ["run_fixture", "run"]
