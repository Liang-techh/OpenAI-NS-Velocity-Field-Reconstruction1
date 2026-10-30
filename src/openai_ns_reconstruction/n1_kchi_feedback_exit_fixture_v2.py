"""Independent v2 synthetic fixture for the raw-Q fixed core feedback exit.

This fixture constructs nonconstant signed source cells, ordinary Z1/Z2 rows,
all required C0/B1 supports, and a separate caller-owned symmetric W0
kernel-tail packet.  It does not import the endpoint module or read an
archive.  ``run_fixture`` is a caller-owned checker and is not run at import.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Any


C0_KEYS = (
    "4,1", "5,1", "5,2", "5,3", "5,5",
    "6,1", "6,2", "6,3", "6,4", "6,6",
)
B1_KEYS = ("5,1", "5,2", "5,3", "6,2", "6,3", "6,4")
FIELDS = ("F1", "Uz1", "K1", "P1", "d_x_F1", "d_x_Uz1")
DELTA_FIELDS = ("F1", "d_x_F1")
ZERO_D_FIELDS = ("Uz1", "K1", "P1", "d_x_Uz1")


def _q(c: Any, value: Fraction) -> Any:
    return c.mpf(value.numerator) / value.denominator


def _t(c: Any, index: int) -> Any:
    left = _q(c, Fraction(index, 8))
    right = _q(c, Fraction(index + 1, 8))
    return c.make_mpf((left._mpi_[0], right._mpi_[1]))


def _term(c: Any, value: int, powers: tuple[int, int, int, int], source: str,
          *, error: bool = False) -> dict[str, Any]:
    classes = ["synthetic_W0_kernel_tail" if error else "synthetic_signed"]
    return {
        "powers": powers,
        "coefficient": c.mpf(value),
        "source_id": source,
        "source_ids": [source],
        "error_classes": classes,
        "uncertainty": error,
        "origin": {"source_ids": [source], "error_classes": classes,
                    "uncertainty": error},
        "error_only": error,
        "symmetric_error": error,
    }


def _symmetric_term(c: Any, radius: Fraction, powers: tuple[int, int, int, int],
                    source: str) -> dict[str, Any]:
    bound = _q(c, radius)
    term = _term(c, 0, powers, source, error=True)
    term["coefficient"] = c.mpf([-bound.b, bound.b])
    term["symmetric_error"] = True
    return term


def make_inputs(c: Any) -> tuple[list[dict[str, Any]], list[dict[str, Any]],
                                list[dict[str, Any]], list[dict[str, Any]],
                                dict[str, list[dict[str, Any]]]]:
    """Return saved W0, matrices, DeltaW0, W0-tail, and chi sectors."""

    saved: list[dict[str, Any]] = []
    matrices: list[dict[str, Any]] = []
    delta: list[dict[str, Any]] = []
    tail: list[dict[str, Any]] = []
    for cell in range(8):
        t = _t(c, cell)
        w0: dict[str, dict[str, list[dict[str, Any]]]] = {}
        d0: dict[str, dict[str, list[dict[str, Any]]]] = {}
        r0: dict[str, dict[str, list[dict[str, Any]]]] = {}
        for order in range(3):
            w0[str(order)] = {}
            d0[str(order)] = {}
            r0[str(order)] = {}
            for column in range(1, 7):
                sign = -1 if (cell + 2 * order + column) % 4 == 0 else 1
                powers = (cell % 2, 1 + (column % 2), order % 2, -1 - order)
                w0[str(order)][str(column)] = [_term(
                    c, sign * (1 + cell + order + column), powers,
                    f"W0:{cell}:Z{order}:C{column}")]
            # Delta W0 is intentionally sparse: only columns 1 and 5.
            for column in (1, 5):
                sign = -1 if (cell + order) % 2 else 1
                d0[str(order)][str(column)] = [_term(
                    c, sign * (1 + cell + order),
                    (0, 0, order % 2, -2 - order),
                    f"DeltaW0:{cell}:Z{order}:C{column}")]
                r0[str(order)][str(column)] = [_symmetric_term(
                    c, Fraction(1 + cell + order, 10000),
                    (1, 0, 0, -2 - order),
                    f"W0-tail:{cell}:Z{order}:C{column}")]
        saved.append({"t": t, "x": 2 * t, "rho": 4 * t * t,
                      "W0_ordinary_Z": w0})
        delta.append({"t": t, "x": 2 * t, "W0_ordinary_Z": d0})
        tail.append({"t": t, "x": 2 * t, "W0_ordinary_Z": r0})

        c0: dict[str, dict[str, list[dict[str, Any]]]] = {"0": {}, "1": {}}
        b1: dict[str, dict[str, list[dict[str, Any]]]] = {"0": {}, "1": {}}
        for order in range(2):
            for key_index, key in enumerate(C0_KEYS):
                sign = -1 if (cell + order + key_index) % 5 == 0 else 1
                c0[str(order)][key] = [_term(
                    c, sign * (2 + cell + order),
                    (0, 1 + (key_index % 2), order % 2, -1 - order),
                    f"C0:{cell}:Z{order}:{key}")]
            for key_index, key in enumerate(B1_KEYS):
                sign = -1 if (cell + 2 * order + key_index) % 3 == 0 else 1
                b1[str(order)][key] = [_term(
                    c, sign * (3 + cell + key_index),
                    (1, 0, (key_index + order) % 2, -1 - order),
                    f"B1:{cell}:Z{order}:{key}")]
        matrices.append({"t": t, "x": 2 * t, "rho": 4 * t * t,
                         "C0": c0, "B1": b1})
    chi = {
        "Z0": [_term(c, 2, (0, 0, 1, 0), "chi:Z0")],
        "Z1": [_term(c, -1, (0, 0, 0, 1), "chi:Z1")],
    }
    return saved, matrices, delta, tail, chi


def _endpoint(bits: Any) -> Fraction:
    sign, mantissa, exponent, _ = bits
    value = Fraction(int(mantissa), 1)
    value = value * (1 << exponent) if exponent >= 0 else value / (1 << (-exponent))
    return -value if sign else value


def _contains(interval: Any, expected: Fraction) -> bool:
    return (_endpoint(interval._mpi_[0]) <= expected <=
            _endpoint(interval._mpi_[1]))


def _all_terms(value: Any):
    for source in value.values():
        fields = source.get("fields", {})
        for order in ("Z0", "Z1"):
            for terms in fields.get(order, {}).values():
                yield from terms


def run_fixture(build: Any, c: Any) -> dict[str, Any]:
    """Run independent fixed-exit, sign, Z1, mass, and tail checks."""

    saved, matrices, delta, tail, chi = make_inputs(c)
    result = build(c, saved, matrices, chi, delta_w0_cells=delta,
                   W0_kernel_tail_cells=tail)
    checks: dict[str, bool] = {}
    target = result.get("exit_target", {})
    checks["fixed_x2_not_last_cell"] = (
        target.get("fixed") is True
        and target.get("independent") is True
        and target.get("not_last_target_cell_proxy") is True
        and target.get("x").a == 2 and target.get("x").b == 2
    )
    checks["three_raw_q_branches"] = all(
        result.get(name, {}).get("formula") == formula
        for name, formula in (("E01", "G0 Q DeltaW0"),
                              ("E10", "D Q W0"),
                              ("E11", "D Q DeltaW0"))
    )
    e01 = result.get("E01", {})
    e10 = result.get("E10", {})
    e11 = result.get("E11", {})
    e01_exits = e01.get("exit_by_source_cell", {})
    e10_exits = e10.get("exit_by_source_cell", {})
    e11_exits = e11.get("exit_by_source_cell", {})
    checks["eight_source_cells"] = all(
        set(branch) == {str(index) for index in range(8)}
        for branch in (e01_exits, e10_exits, e11_exits)
    )
    checks["E01_six_fields"] = (
        checks["eight_source_cells"]
        and all(set(e01_exits[str(index)]["fields"]["Z0"]) == set(FIELDS)
                and set(e01_exits[str(index)]["fields"]["Z1"]) == set(FIELDS)
                for index in range(8))
    )
    checks["D_structural_zeros_explicit"] = (
        e10.get("nonzero_fields") == list(DELTA_FIELDS)
        and e10.get("structural_zero_fields") == list(ZERO_D_FIELDS)
        and all(not e10_exits[str(index)]["fields"]["Z0"][field]
                and not e10_exits[str(index)]["fields"]["Z1"][field]
                for index in range(8) for field in ZERO_D_FIELDS)
    )
    first_mass = e01_exits.get("0", {}).get("kernel_masses", {}).get("F1")
    checks["exact_first_G0_mass"] = (
        first_mass is not None and _contains(first_mass, Fraction(127, 8192))
    )
    full_g0 = result.get("full_green_interval_kernel_masses", {})
    full_d = result.get("full_physical_interval_kernel_masses", {})
    checks["exact_full_masses"] = (
        _contains(full_g0.get("F1"), Fraction(1, 2))
        and _contains(full_d.get("F1"), Fraction(1, 12))
        and _contains(full_d.get("d_x_F1"), Fraction(1, 6))
    )
    checks["nonconstant_source_masses"] = (
        e01_exits.get("0", {}).get("kernel_masses", {}).get("F1")._mpi_
        != e01_exits.get("7", {}).get("kernel_masses", {}).get("F1")._mpi_
    )
    # Q rows are intentionally preserved by row number (4/5/6), not
    # collapsed into one list.  Flatten the row values before inspecting
    # term metadata; iterating the dict itself yields string row labels.
    z1_e10_rows = e10.get("Q_by_source_cell", {}).get("0", {}).get("Z1", {})
    z1_e10 = [term for row in z1_e10_rows.values() for term in row]
    checks["Z1_B1_W0_ZZ"] = any(
        term.get("matrix_name") == "B1"
        and term.get("matrix_order") == 0
        and term.get("W0_order") == 2
        for term in z1_e10
    )
    z1_d_terms = e10_exits.get("0", {}).get("fields", {}).get("Z1", {}).get("F1", [])
    checks["Z1_chi_Z"] = any(term.get("chi_order") == 1 for term in z1_d_terms)
    checks["signed_terms_and_provenance"] = any(
        term.get("coefficient").a < 0 and term.get("source_ids")
        and len(term.get("powers", ())) == 4
        for term in _all_terms(e01_exits)
    )
    tail_result = result.get("W0_kernel_tail_uncertainty") or {}
    tail_branches = tail_result.get("branches", {})
    checks["tail_separate"] = (
        tail_result.get("caller_supplied") is True
        and tail_result.get("error_only") is True
        and tail_result.get("symmetric_error") is True
        and tail_branches.get("G0", {}).get("error_only") is True
        and tail_branches.get("D", {}).get("error_only") is True
    )
    tail_terms = _all_terms(tail_branches.get("G0", {}).get(
        "exit_by_source_cell", {}))
    checks["tail_terms_marked"] = any(
        term.get("error_only") is True
        and term.get("symmetric_error") is True
        and term.get("uncertainty") is True
        for term in tail_terms
    )
    checks["raw_Q_not_P"] = (
        result.get("raw_Q_only") is True
        and result.get("potential_P_compensation_included") is False
        and result.get("potential_P_compensation_merged") is False
    )
    checks["no_E00_or_reconstruction"] = (
        result.get("E00_G0_recomputed") is False
        and result.get("W0_reconstructed") is False
        and result.get("old_builder_called") is False
    )
    return {
        "passed": all(checks.values()),
        "checks": checks,
        "check_count": len(checks),
        "synthetic_only": True,
        "independent_source_is_nonconstant": True,
        "ordinary_Z1_and_Z2_present": True,
        "signed_and_W0_kernel_tail_separate": True,
        "method": "saved raw-Q products with exact fixed-x=2 G0 and D masses",
        "limitations": {
            "potential_compensation": "not included; raw-Q only",
            "complete_Kchi_W1_inlet": False,
            "ancestor_or_receipt_read": False,
        },
    }


run = run_fixture


__all__ = ["C0_KEYS", "B1_KEYS", "FIELDS", "DELTA_FIELDS", "make_inputs",
           "run_fixture", "run"]


