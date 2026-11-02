"""Small independent checks for the target-cell P-tail error leaf.

The fixture creates only eight tiny synthetic cells.  It does not import a
producer or any old kernel helper, and it does not resemble the saved physical
source packet.  Its purpose is to check route bookkeeping, symmetric error
intervals, saved-mass use, the conservative G0 bounds, and the one-time outer
R width.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Any


FIELDS = ("F1", "d_x_F1")
BRANCHES = {
    "G0_P_Rg",
    "D_P_Rg",
    "R_P_G0g",
    "R_P_Dg",
    "R_P_Rg",
}


def _q(c: Any, value: Fraction) -> Any:
    return c.mpf(value.numerator) / value.denominator


def _interval(c: Any, low: Fraction, high: Fraction) -> Any:
    left = _q(c, low)
    right = _q(c, high)
    return c.make_mpf((left._mpi_[0], right._mpi_[1]))


def _record(c: Any, low: Fraction, high: Fraction) -> dict[str, Any]:
    value = _interval(c, low, high)
    return {
        "lower_exact_mpf_tuple": list(value._mpi_[0]),
        "upper_exact_mpf_tuple": list(value._mpi_[1]),
        "directed": True,
    }


def _symmetric_record(c: Any, radius: Fraction) -> dict[str, Any]:
    return _record(c, -radius, radius)


def _term(
    c: Any,
    coefficient: tuple[Fraction, Fraction],
    source: str,
    *,
    powers: tuple[int, int, int, int] = (0, 1, 0, -1),
    error_only: bool = False,
) -> dict[str, Any]:
    classes = (["synthetic_tail", "Kchi_j_ge_2_kernel_tail"] if error_only
               else ["synthetic_source"])
    return {
        "coefficient": _record(c, coefficient[0], coefficient[1]),
        "powers": list(powers),
        "source_id": source,
        "source_ids": [source],
        "error_classes": classes,
        "uncertainty": error_only,
        "error_only": error_only,
        "origin": {
            "source_ids": [source],
            "error_classes": classes,
            "uncertainty": error_only,
        },
    }


def _cell_x(c: Any, index: int) -> dict[str, Any]:
    return {
        "x": _record(c, Fraction(index, 4), Fraction(index + 1, 4)),
        "t": _record(c, Fraction(index, 8), Fraction(index + 1, 8)),
    }


def _synthetic_inputs(c: Any) -> tuple[dict[str, Any], dict[str, Any]]:
    family = {
        "actual_five_defect_family_sha256": "a" * 64,
        "implicit_source_sha256": "b" * 64,
        "datum_enclosure_sha256": "c" * 64,
    }
    logs = {
        name: _record(c, Fraction(1), Fraction(1))
        for name in ("hb", "F0_anchor", "Pstar_squared", "Lambda")
    }

    low_cells: list[dict[str, Any]] = []
    for index in range(8):
        z0 = Fraction(0) if index == 1 else Fraction(-index - 1 if index % 2 == 0 else index + 1)
        z1 = Fraction(-2) if index % 2 == 0 else Fraction(index + 2)
        low_cells.append({
            **_cell_x(c, index),
            "W0_ordinary_Z": {
                "0": {"1": [_term(c, (z0, z0), f"W0:{index}:Z0")]},
                "1": {"1": [_term(c, (z1, z1), f"W0:{index}:Z1")]},
            },
        })

    kchi_cells: list[dict[str, Any]] = []
    for target in range(8):
        d0 = Fraction(-target - 1 if target % 2 == 0 else target + 1)
        d1 = Fraction(0) if target == 2 else Fraction(target + 3)
        by_source: dict[str, Any] = {}
        for source in range(target + 1):
            # These are saved positive TK0 masses, not recomputed by the leaf.
            f1 = Fraction(source + 1, 64 + target)
            dxf1 = Fraction(source + 2, 128 + target)
            by_source[str(source)] = {
                "source_index": source,
                "target_index": target,
                "kernel_masses": {
                    "F1": _record(c, f1, f1),
                    "d_x_F1": _record(c, dxf1, dxf1),
                },
            }
        kchi_cells.append({
            **_cell_x(c, target),
            "index": target,
            "target_index": target,
            "F1": {
                "Z0": [_term(c, (d0, d0), f"Dg:{target}:Z0")],
                "Z1": [_term(c, (d1, d1), f"Dg:{target}:Z1")],
            },
            "W0_correction_by_source_cell": by_source,
        })

    tail = {
        "F1": {
            "0": [_term(c, (Fraction(-1), Fraction(1)), "Rg:F1:Z0", error_only=True)
                  | {"physical_source_width_bound": 2}],
            "1": [
                _term(c, (Fraction(-2), Fraction(2)), "Rg:F1:Z1:main", error_only=True)
                | {"physical_source_width_bound": 2},
                _term(c, (Fraction(0), Fraction(0)), "Rg:F1:Z1:zero", error_only=True)
                | {"physical_source_width_bound": 2},
            ],
        },
        "d_x_F1": {
            "0": [_term(c, (Fraction(-10), Fraction(10)), "Rg:d_x:Z0", error_only=True)],
            "1": [_term(c, (Fraction(-20), Fraction(20)), "Rg:d_x:Z1", error_only=True)],
        },
    }
    bounds = {
        "F1": [
            _record(c, Fraction(1, 4), Fraction(1, 4)),
            _record(c, Fraction(1, 2), Fraction(1, 2)),
            _record(c, Fraction(1), Fraction(1)),
        ],
        "d_x_F1": [
            _record(c, Fraction(1, 8), Fraction(1, 8)),
            _record(c, Fraction(1, 4), Fraction(1, 4)),
            _record(c, Fraction(1, 2), Fraction(1, 2)),
        ],
        "ordinary_Z_orders": [0, 1, 2],
        "scope": "pointwise potential-kernel remainder j>=2 only",
        "source_or_feedback_included": False,
        "bounds_are_not_field_values": True,
    }

    low = {
        "source_family": family,
        "source_logs": logs,
        "outer_cells": low_cells,
    }
    kchi = {
        "source_family": family,
        "source_logs": logs,
        "chi": {
            "chi": [
                _term(c, (Fraction(1, 2), Fraction(1, 2)), "chi:Z0", powers=(0, 0, 0, 0)),
                _term(c, (Fraction(-1, 3), Fraction(-1, 3)), "chi:Z1", powers=(0, 0, 0, 0)),
            ],
        },
        "first_correction": {"target_cells": kchi_cells},
        "W0_separate_symmetric_kernel_tail": tail,
        "kernel_tail_bounds": bounds,
    }
    return low, kchi


def _fraction_endpoint(bits: Any) -> Fraction:
    sign, mantissa, exponent, _count = bits
    value = Fraction(int(mantissa), 1)
    if exponent >= 0:
        value *= 1 << int(exponent)
    else:
        value /= 1 << int(-exponent)
    return -value if sign else value


def _lower(value: Any) -> Fraction:
    return _fraction_endpoint(value._mpi_[0])


def _upper(value: Any) -> Fraction:
    return _fraction_endpoint(value._mpi_[1])


def _all_terms(result: dict[str, Any]) -> list[dict[str, Any]]:
    values: list[dict[str, Any]] = []
    for branch in result["branches"].values():
        for target in branch["target_cells"]:
            for source in target["by_source_cell"].values():
                for order in ("Z0", "Z1"):
                    for field in FIELDS:
                        values.extend(source["fields"][order][field])
    return values


def run_fixture(build: Any, c: Any) -> dict[str, Any]:
    """Run independent small-rational target-cell checks."""

    low, kchi = _synthetic_inputs(c)
    result = build(c, low, kchi)
    checks: dict[str, bool] = {}
    checks["five_branches"] = set(result.get("branches", {})) == BRANCHES
    checks["eight_target_cells"] = all(
        len(result["branches"][name].get("target_cells", [])) == 8
        for name in BRANCHES
    ) if checks["five_branches"] else False
    checks["triangular_source_order"] = all(
        target["source_indices"] == list(range(index + 1))
        and target["not_fixed_exit"] is True
        for index, target in enumerate(result["branches"]["G0_P_Rg"]["target_cells"])
    ) if checks["eight_target_cells"] else False
    checks["target_radius_endpoints"] = all(
        _lower(target["target_radius_segment"][0]) == Fraction(index, 4)
        and _upper(target["target_radius_segment"][1]) == Fraction(index + 1, 4)
        and _upper(target["target_radius_upper_bound"]) == Fraction(index + 1, 4)
        for index, target in enumerate(result["branches"]["G0_P_Rg"]["target_cells"])
    ) if checks["eight_target_cells"] else False

    terms = _all_terms(result)
    checks["nonempty_error_stream"] = bool(terms)
    checks["all_error_only"] = bool(terms) and all(
        term.get("error_only") is True
        and term.get("uncertainty") is True
        and term.get("symmetric_error") is True
        and term.get("signed_field") is False
        for term in terms
    )
    checks["symmetric_endpoints"] = bool(terms) and all(
        _lower(term["coefficient"]) == -_upper(term["coefficient"])
        for term in terms
    )
    checks["four_powers_and_parent_ids"] = bool(terms) and all(
        len(term.get("powers", ())) == 4
        and term.get("source_ids")
        and term.get("parent_source_ids")
        and term.get("parent_source_ids_by_factor")
        and all(term.get("parent_powers_by_factor", ()))
        and all(len(powers) == 4 for powers in term.get("parent_powers_by_factor", ()))
        and term.get("target_cell") is not None
        and term.get("source_cell") is not None
        for term in terms
    )
    checks["P_Z1_both_routes"] = (
        any("chi_Z*V_Z0" in term.get("parent_routes", [])
            for term in result["branches"]["G0_P_Rg"]["target_cells"][3]
            ["by_source_cell"]["0"]["fields"]["Z1"]["F1"])
        and any("chi*V_Z1" in term.get("parent_routes", [])
                for term in result["branches"]["G0_P_Rg"]["target_cells"][3]
                ["by_source_cell"]["0"]["fields"]["Z1"]["F1"])
    )
    checks["D_Z1_both_routes"] = (
        any(term.get("route") == "chi_Z*P_Z0"
            for term in result["branches"]["D_P_Rg"]["target_cells"][3]
            ["by_source_cell"]["0"]["fields"]["Z1"]["F1"])
        and any(term.get("route") == "chi*P_Z1"
                for term in result["branches"]["D_P_Rg"]["target_cells"][3]
                ["by_source_cell"]["0"]["fields"]["Z1"]["F1"])
    )
    checks["R_Z1_derivative_routes"] = (
        any(term.get("route") == "R_Z1*P_Z0"
            for term in result["branches"]["R_P_Rg"]["target_cells"][3]
            ["by_source_cell"]["0"]["fields"]["Z1"]["F1"])
        and any(term.get("route") == "R_Z0*P_Z1"
                for term in result["branches"]["R_P_Rg"]["target_cells"][3]
                ["by_source_cell"]["0"]["fields"]["Z1"]["F1"])
    )

    g0_term = result["branches"]["G0_P_Rg"]["target_cells"][0]
    g0_term = g0_term["by_source_cell"]["0"]["fields"]["Z0"]["F1"][0]
    checks["G0_full_triangle_upper"] = (
        _lower(g0_term["kernel_mass_upper_bound"]) <= Fraction(1, 128)
        <= _upper(g0_term["kernel_mass_upper_bound"])
        and g0_term.get("conservative_G0_mass_upper") is True
    )
    checks["independent_z0_rational_bounds"] = all(
        _lower(result["branches"][branch]["target_cells"][0]["by_source_cell"]["0"]["fields"]["Z0"]["F1"][0]["coefficient"]) == -radius
        and _upper(result["branches"][branch]["target_cells"][0]["by_source_cell"]["0"]["fields"]["Z0"]["F1"][0]["coefficient"]) == radius
        for branch, radius in (("G0_P_Rg", Fraction(1, 128)),
                               ("D_P_Rg", Fraction(1, 64)),
                               ("R_P_G0g", Fraction(1, 16)),
                               ("R_P_Dg", Fraction(1, 16)),
                               ("R_P_Rg", Fraction(1, 16)))
    )
    checks["E51_uses_F1_input_for_dx_output"] = not any(
        "Rg:d_x:" in str(source_id)
        for term in terms for source_id in term.get("source_ids", [])
    )
    d_terms = result["branches"]["D_P_Rg"]["target_cells"][2]
    d_terms = d_terms["by_source_cell"]["1"]["fields"]["Z0"]["F1"]
    checks["D_saved_mass_only"] = bool(d_terms) and all(
        term.get("saved_TK0_mass_used") is True
        and term.get("saved_TK0_mass_reintegrated") is False
        and term.get("kernel_mass_source") == "saved_kchi_TK0_mass"
        for term in d_terms
    )
    r_term = result["branches"]["R_P_G0g"]["target_cells"][2]
    r_term = r_term["by_source_cell"]["1"]["fields"]["Z1"]["F1"][0]
    checks["R_width_once"] = (
        _lower(r_term["outer_source_width_bound"]) == Fraction(1, 4)
        and _upper(r_term["outer_source_width_bound"]) == Fraction(1, 4)
        and r_term.get("outer_source_width_applied_once") is True
        and r_term.get("inner_Rg_width_readded") is False
    )
    checks["zero_and_negative_inputs_retained"] = (
        any(_upper(term["coefficient"]) == 0 for term in result["branches"]["R_P_G0g"]
            ["target_cells"][1]["by_source_cell"]["1"]["fields"]["Z0"]["F1"])
        and any(term.get("parent_input_may_be_negative") is True for term in terms)
    )
    checks["scope_flags"] = (
        result.get("fixed_x2_exit_computed") is False
        and result.get("feedback_word_tail_included") is False
        and result.get("collar_included") is False
        and result.get("pressure_included") is False
        and result.get("appended_to_ledger") is False
        and result.get("complete_W1") is False
    )
    checks["directed_formal_policy"] = (
        result.get("formal_powers_not_expanded") is True
        and result.get("coefficients_are_field_values") is False
        and "c.make_mpf" in result.get("directed_endpoint_policy", "")
    )
    checks["conditional_source_binding"] = (
        result.get("source_binding", {}).get("family_equal_conditionally") is True
        and result.get("source_binding", {}).get("logs_equal_conditionally") is True
        and result.get("source_binding", {}).get("live_source_owner_proved") is False
    )
    checks["no_aggregate_or_cancel"] = all(
        "aggregate" not in result["branches"][name]
        and all(source.get("no_aggregation") is True
                and source.get("no_cancellation") is True
                for target in result["branches"][name]["target_cells"]
                for source in target["by_source_cell"].values())
        for name in BRANCHES
    )
    return {
        "passed": all(checks.values()),
        "checks": checks,
        "check_count": len(checks),
        "synthetic_only": True,
        "small_rational_inputs": True,
        "nonzero_chi_Z1": True,
        "no_old_module_called": True,
        "limitations": {
            "fixed_x2_exit": False,
            "feedback_word_tail": False,
            "collar": False,
            "pressure": False,
            "complete_W1": False,
        },
    }


run = run_fixture


__all__ = ["run_fixture", "run"]
