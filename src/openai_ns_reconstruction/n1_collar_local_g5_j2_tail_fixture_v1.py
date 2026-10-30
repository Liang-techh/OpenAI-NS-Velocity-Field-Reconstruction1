"""Independent synthetic checks for the local ``g5`` Kchi tail.

The data below are deliberately distinct from the entered first-word and
inlet-tail fixtures.  The fixture checks the whole-collar interface,
source-specific signs and formal powers, both ordinary-Z1 insertion routes,
the separate symmetric remainder, and the zero-chi/nonzero-chi_Z case.  It
does not import a producer, read a receipt, call an older helper, or perform a
quadrature calculation.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Any


SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
BRANCHES = ("F_from_g5", "Fx_from_g5")
FORMAL_HB_POWER_ADDED = {"F_from_g5": 6, "Fx_from_g5": 5}


def _q(c: Any, value: Fraction | int) -> Any:
    if isinstance(value, Fraction):
        return c.mpf(value.numerator) / value.denominator
    return c.mpf(value)


def _box(c: Any, lower: Fraction | int, upper: Fraction | int) -> Any:
    lo, hi = _q(c, lower), _q(c, upper)
    return c.make_mpf((lo._mpi_[0], hi._mpi_[1]))


def _term(
    c: Any,
    coefficient: Fraction | int | Any,
    source_id: str,
    powers: tuple[int, int, int, int],
    *,
    origin: Any,
    error_only: bool = False,
    error_class: str | None = None,
) -> dict[str, Any]:
    class_name = error_class or ("fixture_local_g5_error" if error_only else "fixture_local_g5_signed")
    result = {
        "coefficient": coefficient if hasattr(coefficient, "_mpi_") else _q(c, coefficient),
        "powers": list(powers),
        "source_ids": [source_id],
        "error_classes": [class_name],
        "origin": origin,
        "uncertainty": error_only,
    }
    if error_only:
        result.update({"error_only": True, "symmetric_error": True})
    return result


def _bounds(c: Any, value: Any) -> tuple[Any, Any]:
    lo, hi = value._mpi_
    return c.make_mpf((lo, lo)), c.make_mpf((hi, hi))


def _upper(c: Any, value: Any) -> Any:
    return _bounds(c, value)[1]


def _nonzero(c: Any, value: Any) -> bool:
    lo, hi = _bounds(c, value)
    return lo != 0 or hi != 0


def _contains_zero(c: Any, value: Any) -> bool:
    lo, hi = _bounds(c, value)
    return lo <= 0 <= hi


def _source_powers(source: Any) -> dict[str, tuple[int, int, int, int]]:
    result: dict[str, tuple[int, int, int, int]] = {}
    for order in ("Z0", "Z1"):
        for item in source[order]:
            result[item["source_ids"][0]] = tuple(item["powers"])
    return result


def _all_terms(result: Any, branch: str) -> list[dict[str, Any]]:
    return list(result["branches"][branch]["Z0"]) + list(result["branches"][branch]["Z1"])


def run_fixture(build: Any, c: Any) -> dict[str, Any]:
    """Run independent structural and directed-interval checks."""

    if hasattr(build, "build"):
        build = build.build
    if not callable(build):
        raise TypeError("build must be callable or expose build")

    # New rationals and powers: these are not the v3 first-word or inlet-tail
    # fixture values.  Both signed and inherited symmetric source sectors are
    # retained as separate whole-collar inputs.
    hb = _box(c, 0, Fraction(1, 3))
    sc = _q(c, Fraction(1, 96))
    source = {
        "Z0": [
            _term(c, Fraction(7, 19), "v1:g5:positive", (-2, 1, -1, 2), origin="v1:g5:positive"),
            _term(c, Fraction(-5, 23), "v1:g5:negative", (1, 0, 2, -1), origin={"tag": "v1:g5:negative"}),
        ],
        "Z1": [
            _term(c, Fraction(-4, 29), "v1:g5:z1_signed", (-3, 2, 0, 1), origin="v1:g5:z1_signed"),
            _term(c, _box(c, Fraction(-1, 31), Fraction(1, 31)), "v1:g5:z1_error", (-1, 0, 1, -2), origin={"tag": "v1:g5:z1_error"}, error_only=True),
        ],
    }
    chi = {
        "Z0": [_term(c, Fraction(-3, 8), "v1:chi:Z0", (0, 0, 0, 0), origin="v1:chi:Z0", error_class="fixture_chi_z0")],
        "Z1": [_term(c, Fraction(5, 13), "v1:chi:Z1", (0, 0, 0, 0), origin="v1:chi:Z1", error_class="fixture_chi_z1")],
    }
    certificate = {
        "hb_sc_Z_independent": True,
        "rho_radius_Z_independent": True,
        "strict_hb_admitted": True,
        "fixture_geometry_owner": "independent_local_g5_v1",
    }
    atoms = {name: "v1:" + name for name in SOURCE_POWER_BASIS}
    kwargs = {
        "hb": hb,
        "sc": sc,
        "g5": source,
        "chi": chi,
        "source_family": "synthetic_local_g5_j2_v1_distinct",
        "source_atoms": atoms,
        "geometry_certificate": certificate,
    }
    result = build(c, **kwargs)
    source_power_map = _source_powers(source)

    checks: dict[str, bool] = {
        "schema": result.get("schema") == "n1_collar_local_g5_j2_tail_v1",
        "fields": tuple(result.get("fields", ())) == ("F", "Fx"),
        "legacy_aliases": tuple(result.get("legacy_fields", ())) == ("F1", "d_x_F1"),
        "formal_basis": tuple(result.get("source_power_basis", ())) == SOURCE_POWER_BASIS,
        "source_family": result.get("source_family") == "synthetic_local_g5_j2_v1_distinct",
        "whole_collar_input": len(result.get("g5_jets", {}).get("Z0", [])) == 2 and len(result.get("g5_jets", {}).get("Z1", [])) == 2,
        "nonzero_source_sectors": all(_nonzero(c, item["coefficient"]) for items in source.values() for item in items),
        "groups_are_not_merged": set(result.get("branches", {})) == set(BRANCHES),
        "geometry_certificate": result.get("geometry_certificate", {}).get("strict_hb_admitted") is True,
        "closed_hb_geometry": result.get("geometry", {}).get("closed_hb_enclosure_allowed") is True,
        "formal_hb_not_materialized": result.get("geometry", {}).get("formal_hb_not_materialized") is True,
        "no_midpoint": result.get("geometry", {}).get("midpoint_values_used") is False,
        "initial_K0_once": result.get("initial_K0_attached_once") is True and result.get("geometry", {}).get("initial_K0_physical_dx_count") == 1,
        "no_double_integration": result.get("physical_dx_double_counted") is False and all(
            branch.get("double_integration") is False for branch in result.get("branches", {}).values()
        ),
        "j2_only": result.get("j2_kernel_remainder_included") is True and result.get("j1_included") is False,
        "operator_scope": all(branch.get("operator") == "(Kchi-K0-D)g5" for branch in result.get("branches", {}).values()),
        "chi_majorant_both_rows_once": result.get("ordinary_Z1_chi_majorant_uses_both_chi_rows") is True and result.get("ordinary_Z1_chi_factor_count") == 1,
        "no_pressure_or_m15": result.get("P1_included") is False and result.get("m15_computed") is False,
        "separate_symmetric_outputs": all(
            term.get("error_only") is True
            and term.get("symmetric_error") is True
            and _contains_zero(c, term["coefficient"])
            for branch in BRANCHES for term in _all_terms(result, branch)
        ),
        "signed_source_provenance": all(
            term.get("signed_source_term") is False
            and term.get("signed_tail_value") is False
            and term.get("g5_source_ids")
            for branch in BRANCHES for term in _all_terms(result, branch)
        ),
    }

    for branch in BRANCHES:
        addition = FORMAL_HB_POWER_ADDED[branch]
        terms = _all_terms(result, branch)
        checks[branch + ":nonempty"] = bool(terms)
        checks[branch + ":formal_increment"] = all(
            tuple(term["powers"]) == (
                source_power_map[term["g5_source_ids"][0]][0] + addition,
                source_power_map[term["g5_source_ids"][0]][1],
                source_power_map[term["g5_source_ids"][0]][2],
                source_power_map[term["g5_source_ids"][0]][3],
            )
            for term in terms
        )
        checks[branch + ":chi_ids"] = all(
            "v1:chi:Z0" in term.get("chi_source_ids", [])
            or "v1:chi:Z1" in term.get("chi_source_ids", [])
            for term in terms
        )
        checks[branch + ":tail_metadata"] = all(
            term.get("tail_starts_at_j") == 2
            and term.get("kernel_tail_order") == "j>=2"
            and term.get("operator_word_sign") == "(-1)^j"
            and term.get("physical_measure_attached_once") is True
            for term in terms
        )

    for field, branch in (("F", "F_from_g5"), ("Fx", "Fx_from_g5")):
        z1 = result["response"][field]["Z1"]
        checks[field + ":Z1_source_route"] = any(
            term.get("ordinary_Z_piece") == "chi*g5_Z1" and term.get("source_order") == 1 and term.get("chi_order") == 0
            for term in z1
        )
        checks[field + ":Z1_chi_route"] = any(
            term.get("ordinary_Z_piece") == "chi_Z*g5_Z0" and term.get("source_order") == 0 and term.get("chi_order") == 1
            for term in z1
        )
        checks[field + ":Z1_chi_both_provenance"] = any(
            term.get("ordinary_Z_piece") == "chi_Z*g5_Z0"
            and set(term.get("chi_source_ids", ())) == {"v1:chi:Z0", "v1:chi:Z1"}
            and set(term.get("chi_source_classes", ())) == {"fixture_chi_z0", "fixture_chi_z1"}
            and set(term.get("chi_origins", ())) == {"v1:chi:Z0", "v1:chi:Z1"}
            for term in z1
        )
        checks[field + ":Z1_chi_no_power_double"] = all(
            term.get("chi_formal_powers_multiplied") == [[0, 0, 0, 0]]
            and term.get("chi_factor_count") == 1
            and term.get("chi_provenance_not_multiplied") is True
            for term in z1 if term.get("ordinary_Z_piece") == "chi_Z*g5_Z0"
        )

    checks["formal_counts"] = (
        result["branches"]["F_from_g5"]["formal_hb_power_added"] == 6
        and result["branches"]["Fx_from_g5"]["formal_hb_power_added"] == 5
        and result["kernel_bounds"]["m0_formal_hb_power"] == 2
    )
    checks["whole_collar_formulas"] = (
        result["kernel_bounds"]["formulas"]["m0"] == "2*c**2*E**2"
        and result["component_error_radii"]["F"]["R_F_Z0"] == "S0*hb**6*m0*k0**2*T/2"
        and result["component_error_radii"]["Fx"]["R_Fx_Z0"] == "S0*hb**5*m0*kx*k0*T/2"
    )

    # Independent route guard: with chi_Z0 = 0 but chi_Z1 nonzero, the Z0
    # tail vanishes while the chi_Z insertion in Z1 remains available.
    zero_chi = {
        "Z0": [_term(c, 0, "v1:zero:chi:Z0", (0, 0, 0, 0), origin="v1:zero:chi:Z0")],
        "Z1": chi["Z1"],
    }
    zero_result = build(c, **dict(kwargs, chi=zero_chi))
    checks["zero_chi_nonzero_chiZ"] = (
        all(_upper(c, term["coefficient"]) == 0 for term in zero_result["response"]["F"]["Z0"])
        and any(
            term.get("ordinary_Z_piece") == "chi_Z*g5_Z0"
            for term in zero_result["response"]["F"]["Z1"]
        )
    )

    duplicate_chi = {
        "Z0": chi["Z0"] + [_term(c, Fraction(1, 37), "v1:chi:extra", (0, 0, 0, 0), origin="v1:chi:extra")],
        "Z1": chi["Z1"],
    }
    try:
        build(c, **dict(kwargs, chi=duplicate_chi))
    except (TypeError, ValueError):
        checks["single_chi_sector_guard"] = True
    else:
        checks["single_chi_sector_guard"] = False

    return {
        "passed": all(checks.values()),
        "check_count": len(checks),
        "checks": checks,
        "synthetic_only": True,
        "distinct_fixture_inputs": True,
        "independent_quadrature": False,
        "old_helpers_called": False,
        "receipt_reads": False,
        "mathematical_scope": "symmetric j>=2 local g5 tail only; no P1/m15/feedback/closure",
    }


run = run_fixture

__all__ = [
    "SOURCE_POWER_BASIS", "BRANCHES", "FORMAL_HB_POWER_ADDED",
    "run_fixture", "run",
]

