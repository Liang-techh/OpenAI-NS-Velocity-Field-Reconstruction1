"""Independent v2 fixture for chi order-key normalization and g5 exit."""

from __future__ import annotations

from fractions import Fraction
from typing import Any


def _q(c: Any, value: Fraction) -> Any:
    return c.mpf(value.numerator) / value.denominator


def _t(c: Any, index: int) -> Any:
    left = _q(c, Fraction(index, 8))
    right = _q(c, Fraction(index + 1, 8))
    return c.make_mpf((left._mpi_[0], right._mpi_[1]))


def _term(c: Any, coefficient: int, source: str,
          *, uncertainty: bool = False) -> dict[str, Any]:
    return {
        "powers": (coefficient % 2, 1, 0, -1),
        "coefficient": c.mpf(coefficient),
        "source_id": source,
        "source_ids": [source],
        "error_classes": (["synthetic_R"] if uncertainty else ["synthetic_signed"]),
        "uncertainty": uncertainty,
        "origin": {"source_ids": [source],
                    "error_classes": (["synthetic_R"] if uncertainty else ["synthetic_signed"]),
                    "uncertainty": uncertainty},
    }


def _cells(c: Any, *, error: bool = False) -> list[dict[str, Any]]:
    cells: list[dict[str, Any]] = []
    for index in range(8):
        cells.append({
            "t": _t(c, index),
            "F1": {
                "Z0": [_term(c, (-1 if index % 2 else 1) * (index + 1),
                              f"g5:{index}:Z0", uncertainty=error)],
                "Z1": [_term(c, (-1 if index % 3 == 0 else 1) * (index + 2),
                              f"g5:{index}:Z1", uncertainty=error)],
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
    """Check all supported chi key spellings without importing v1."""

    cells = _cells(c)
    error = _cells(c, error=True)
    chi_terms = {
        "Z0": [_term(c, 2, "chi:Z0")],
        "Z1": [_term(c, -1, "chi:Z1")],
    }
    key_forms = (
        {"Z0": chi_terms["Z0"], "Z1": chi_terms["Z1"]},
        {"0": chi_terms["Z0"], "1": chi_terms["Z1"]},
        {0: chi_terms["Z0"], 1: chi_terms["Z1"]},
        {"Z0": chi_terms["Z0"], 1: chi_terms["Z1"]},
    )
    outputs = [build(c, cells, chi, symmetric_error_cells=error)
               for chi in key_forms]
    result = outputs[0]
    checks: dict[str, bool] = {}
    checks["all_key_forms_accept"] = all(
        output.get("exit_target", {}).get("fixed") is True
        and set(output.get("chi_jets", {})) == {0, 1}
        for output in outputs
    )
    exits = result.get("exit_by_source_cell", {})
    checks["eight_sources"] = set(exits) == {str(index) for index in range(8)}
    checks["z0_z1_present"] = (
        checks["eight_sources"]
        and set(exits["0"]["fields"]) == {"Z0", "Z1"}
    )
    checks["signed_source"] = any(
        term.get("coefficient").a < 0
        for source in exits.values()
        for term in source.get("fields", {}).get("Z0", {}).get("F1", [])
    ) if exits else False
    checks["first_source_mass_positive"] = (
        checks["eight_sources"]
        and exits["0"]["kernel_masses"]["F1"].a > 0
    )
    full = result.get("full_physical_interval_kernel_masses", {})
    checks["full_mass"] = (
        _contains(full.get("F1"), Fraction(1, 12))
        and _contains(full.get("d_x_F1"), Fraction(1, 6))
    ) if full else False
    first = exits.get("0", {}).get("fields", {}).get("Z1", {}).get("F1", [])
    checks["chi_order_metadata"] = bool(first) and all(
        term.get("chi_order") in (0, 1)
        and term.get("source_order") in (0, 1)
        and len(term.get("powers", ())) == 4
        and term.get("source_ids")
        for term in first
    )
    error_result = result.get("symmetric_R_error") or {}
    checks["uncertainty_separate"] = (
        result.get("error_branch_separate") is True
        and error_result.get("error_only") is True
        and error_result.get("branch") is not result.get("branch")
    )
    checks["no_v1_or_complete_claim"] = (
        result.get("complete_Kchi_W1_inlet") is False
        and result.get("complete_feedback") is False
        and result.get("normalization_conjugation_compensation_included") is False
    )
    return {
        "passed": all(checks.values()),
        "checks": checks,
        "check_count": len(checks),
        "synthetic_only": True,
        "key_forms": ["int", "0/1", "Z0/Z1"],
        "no_import_or_call_v1": True,
        "limitations": {"complete_Kchi_W1_inlet": False,
                         "normalization_compensation": "not claimed"},
    }


run = run_fixture


__all__ = ["run_fixture", "run"]
