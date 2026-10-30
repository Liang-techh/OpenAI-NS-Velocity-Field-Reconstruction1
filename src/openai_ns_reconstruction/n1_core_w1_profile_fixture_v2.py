"""Independent synthetic checks for the core W1 target-radius profile.

This fixture passes saved ``W0_ordinary_Z`` outer cells directly to the leaf.
It does not construct an ancestor source, read a receipt, or import the
producer.  Checks cover the two ordinary product rules, source ordering,
axis/triangular masses, and the separate fixed ``x=2`` exit.
"""
from __future__ import annotations

from fractions import Fraction
from typing import Any


FIELDS = (("F1", "5"), ("Uz1", "6"), ("K1", "6"), ("P1", "4"), ("d_x_F1", "5"), ("d_x_Uz1", "6"))
C0_KEYS = ("4,1", "5,1", "5,2", "5,3", "5,5", "6,1", "6,2", "6,3", "6,4", "6,6")
B1_KEYS = ("5,1", "5,2", "5,3", "6,2", "6,3", "6,4")


def _q(c: Any, value: Fraction) -> Any:
    return c.mpf(value.numerator) / value.denominator


def _cell_t(c: Any, index: int) -> Any:
    left = _q(c, Fraction(index, 8))
    right = _q(c, Fraction(index + 1, 8))
    return c.make_mpf((left._mpi_[0], right._mpi_[1]))


def _term(c: Any, f0: int, lam: int, value: int, source: str) -> dict[str, Any]:
    return {"F0_power": f0, "Lambda_power": lam, "coefficient": c.mpf(value), "source_cell": source}


def _endpoint(bits: Any) -> Fraction:
    sign, mantissa, exponent, _ = bits
    value = Fraction(int(mantissa), 1)
    value = value * (1 << exponent) if exponent >= 0 else value / (1 << (-exponent))
    return -value if sign else value


def _contains(interval: Any, expected: Fraction) -> bool:
    return _endpoint(interval._mpi_[0]) <= expected <= _endpoint(interval._mpi_[1])


def _inputs(c: Any) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    outer = []
    matrices = []
    for cell in range(8):
        t = _cell_t(c, cell)
        w0 = {}
        for order in range(3):
            w0[str(order)] = {}
            for column in range(1, 7):
                sign = -1 if (cell + order + column) % 3 == 0 else 1
                w0[str(order)][str(column)] = [_term(c, (cell + column) % 3, -1 - order, sign * (order + 1), f"W0-{cell}-{order}-{column}")]
        outer.append({"t": t, "rho": 4 * t * t, "W0_ordinary_Z": w0})

        c0 = {str(order): {} for order in range(2)}
        b1 = {str(order): {} for order in range(2)}
        for order in range(2):
            for key_index, key in enumerate(C0_KEYS):
                sign = -1 if (cell + order + key_index) % 4 == 0 else 1
                c0[str(order)][key] = [_term(c, (cell + key_index) % 3, -1 - order, sign * (1 + order), f"C0-{cell}-{order}-{key}")]
            for key_index, key in enumerate(B1_KEYS):
                sign = -1 if (cell + 2 * order + key_index) % 4 == 0 else 1
                b1[str(order)][key] = [_term(c, (cell + key_index + 1) % 3, -1 - order, sign * (2 + order), f"B1-{cell}-{order}-{key}")]
        matrices.append({"t": t, "rho": 4 * t * t, "C0": c0, "B1": b1})
    return outer, matrices


def run_fixture(build: Any, c: Any) -> dict[str, Any]:
    """Run synthetic inclusion and provenance checks against ``build``."""

    outer, matrices = _inputs(c)
    result = build(c, outer, matrices)
    checks: dict[str, bool] = {}
    targets = result.get("target_cells", [])
    checks["eight_targets"] = len(targets) == 8
    checks["saved_w0_reused"] = result.get("saved_outer_cells_reused") is True and result.get("W0_reconstructed") is False
    checks["source_j_le_target_i"] = all(cell.get("source_indices") == list(range(i + 1)) for i, cell in enumerate(targets))
    checks["triangular_source_marked"] = all(cell.get("triangular_diagonal_source") == i for i, cell in enumerate(targets))
    checks["axis_continuation_marked"] = bool(targets) and targets[0].get("axis_endpoint_continuation") is True

    for order in ("Z0", "Z1"):
        for cell in targets:
            family = cell.get("W1_by_order", {}).get(order, {})
            checks[f"{order}_fields"] = checks.get(f"{order}_fields", True) and set(family) == {name for name, _ in FIELDS}
            for field in ("F1", "Uz1"):
                for term in family.get(field, []):
                    checks[f"{order}_directed"] = checks.get(f"{order}_directed", True) and hasattr(term["coefficient"], "_mpi_")
                    checks[f"{order}_metadata"] = checks.get(f"{order}_metadata", True) and term.get("matrix_order") in ("0", "1") and term.get("matrix_key") and term.get("matrix_term_index") == 0 and term.get("W0_order") in ("0", "1", "2") and term.get("W0_term_index") == 0 and term.get("W0_column") in range(1, 7)
                    checks[f"{order}_source_ids"] = checks.get(f"{order}_source_ids", True) and term.get("matrix_source_cell") is not None and term.get("W0_source_cell") is not None and term.get("W0_provenance") == "saved_outer_cells"
    checks["z1_origins"] = set(term.get("origin") for term in targets[1]["W1_by_order"]["Z1"]["F1"]) >= {"C0_Z W0", "C0 W0_Z", "B1_Z W0_Z", "B1 W0_ZZ"}
    checks["no_matrix_order2"] = all(term.get("matrix_order") in ("0", "1") for cell in targets for term in cell["W1_by_order"]["Z1"]["F1"])
    checks["no_w0_order3"] = all(term.get("W0_order") in ("0", "1", "2") for cell in targets for term in cell["W1_by_order"]["Z1"]["F1"])

    axis_f1 = targets[0]["W1_by_order"]["Z0"]["F1"]
    axis_uz = targets[0]["W1_by_order"]["Z0"]["Uz1"]
    checks["axis_f1_mass"] = bool(axis_f1) and _contains(axis_f1[0]["kernel_mass"], Fraction(1, 128))
    checks["axis_uz_mass"] = bool(axis_uz) and _contains(axis_uz[0]["kernel_mass"], Fraction(1, 64))

    exits = result.get("W1_exit_x2_by_source_cell", {})
    checks["exit_orders"] = set(exits) == {"Z0", "Z1"}
    checks["exit_all_sources"] = all(set(exits.get(order, {})) == {str(i) for i in range(8)} for order in ("Z0", "Z1"))
    checks["exit_fixed_x2"] = result.get("exit_target", {}).get("fixed") is True and result["exit_target"]["x"].a == 2 and result["exit_target"]["x"].b == 2
    e_f1 = exits["Z0"]["0"]["F1"] if checks["exit_all_sources"] else []
    e_p1 = exits["Z0"]["0"]["P1"] if checks["exit_all_sources"] else []
    e_uz_outer = exits["Z0"]["1"]["Uz1"] if checks["exit_all_sources"] else []
    checks["exit_f1_exact_mass"] = bool(e_f1) and _contains(e_f1[0]["kernel_mass"], Fraction(127, 8192))
    checks["exit_p1_exact_mass"] = bool(e_p1) and _contains(e_p1[0]["kernel_mass"], Fraction(1, 4))
    checks["exit_nonzero_outer_Uz1_log_mass"] = bool(e_uz_outer) and e_uz_outer[0]["kernel_mass"].a > 0
    last_target = targets[7]["W1_by_order"]["Z0"]["F1"]
    last_exit = exits["Z0"]["7"]["F1"] if checks["exit_all_sources"] else []
    checks["last_target_not_exit_proxy"] = bool(last_target and last_exit) and last_target[-1]["kernel_mass"]._mpi_ != last_exit[-1]["kernel_mass"]._mpi_

    return {
        "passed": all(checks.values()),
        "checks": checks,
        "check_count": len(checks),
        "synthetic_only": True,
        "fixture_is_not_actual_field": True,
        "method": "independent saved-W0 cell bookkeeping and exact K0 source-cell mass checks",
        "model_metadata": {"principal_kernel": "K0 only", "Kchi_higher_words_included": False, "full_W1_claimed": False, "fixed_exit": "x=2 exact cell masses", "axis_log_singularity_divided_directly": False},
    }


__all__ = ["run_fixture"]
