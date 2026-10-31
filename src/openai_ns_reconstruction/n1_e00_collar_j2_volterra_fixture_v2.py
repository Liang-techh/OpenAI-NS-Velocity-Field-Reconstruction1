"""Independent synthetic checks for the fixed-E00 Volterra tail leaf.

The fixture owns all terms, labels, powers and q cells.  It does not import an
ancestor producer, read a receipt, or run on import.  ``run_fixture`` receives
the leaf and the caller-owned directed interval context explicitly.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Any, Mapping


SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
BRANCHES = ("F_from_Fc", "F_from_Fxc", "Fx_from_Fc", "Fx_from_Fxc")
FORMAL_HB_POWER_ADDED = {
    "F_from_Fc": 4,
    "F_from_Fxc": 5,
    "Fx_from_Fc": 3,
    "Fx_from_Fxc": 4,
}
SENTINEL_CLASS = "fixture:v2:E00:j2:sentinel_error"
SENTINEL_ID = "fixture:v2:E00:Fx:Z1:sentinel"
Q = Fraction

__all__ = ["SOURCE_POWER_BASIS", "BRANCHES", "SENTINEL_CLASS", "SENTINEL_ID", "run_fixture", "run"]


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
    path: str,
    error_class: str,
    error_only: bool = False,
) -> dict[str, Any]:
    value = coefficient if hasattr(coefficient, "_mpi_") else _q(c, coefficient)
    result = {
        "coefficient_interval": value,
        "powers": list(powers),
        "source_ids": [source_id],
        "source_error_classes": [error_class],
        "error_classes": [error_class],
        "uncertainty": error_only,
        "origin": {"operation": "independent_fixture_v2", "term_path": path, "source_ids": [source_id]},
        "term_path": path,
    }
    if error_only:
        result.update({"error_only": True, "symmetric_error": True})
    return result


def _rows(c: Any, prefix: str, values: tuple[Any, Any], powers: tuple[int, int, int, int], *, error_z1: bool = False) -> dict[str, list[dict[str, Any]]]:
    return {
        "Z0": [_term(c, values[0], f"{prefix}:Z0", powers, path=f"{prefix}.Z0", error_class=f"{prefix}:signed")],
        "Z1": [
            _term(
                c,
                values[1],
                SENTINEL_ID if error_z1 else f"{prefix}:Z1",
                powers,
                path=f"{prefix}.Z1",
                error_class=SENTINEL_CLASS if error_z1 else f"{prefix}:signed",
                error_only=error_z1,
            )
        ],
    }


def _inlet(c: Any) -> dict[str, dict[str, list[dict[str, Any]]]]:
    return {
        "F1": _rows(c, "fixture:v2:E00:F1", (Q(-7, 31), Q(11, 37)), (4, -3, 0, 2)),
        "Fx": _rows(c, "fixture:v2:E00:Fx", (Q(13, 41), Q(-17, 43)), (-5, 4, 2, -1), error_z1=True),
    }


def _chi(c: Any, *, zero: bool = False) -> tuple[dict[str, Any], dict[str, Any]]:
    chi0_value = _box(c, 0, 0) if zero else _box(c, Q(-3, 19), Q(5, 23))
    chi0 = _term(c, chi0_value, "fixture:v2:E00:chi:Z0", (0, 0, 0, 0), path="fixture.chi.Z0", error_class="fixture:v2:chi")
    chi1 = _term(c, _box(c, Q(3, 29), Q(4, 17)), "fixture:v2:E00:chi:Z1", (0, 0, 0, 0), path="fixture.chi.Z1", error_class="fixture:v2:chi_Z")
    return chi0, chi1


def _zero_interval(value: Any) -> bool:
    raw = value._mpi_
    return raw[0][1] == 0 and raw[1][1] == 0


def _symmetric_interval(value: Any) -> bool:
    return value.a == -value.b


def _terms(result: Mapping[str, Any]):
    for branch in BRANCHES:
        yield from result["branches"][branch]["Z0"]
        yield from result["branches"][branch]["Z1"]


def _guard_checks(build: Any, c: Any, inlet: Any, chi0: Any, chi1: Any, common: Mapping[str, Any]) -> dict[str, bool]:
    checks: dict[str, bool] = {}
    try:
        build(c, inlet=inlet, chi=chi0, hb=common["hb"], sc=common["sc"], q=common["q"])
    except (TypeError, ValueError):
        checks["missing_chi_Z_rejected"] = True
    else:
        checks["missing_chi_Z_rejected"] = False

    try:
        build(c, inlet=inlet, chi=chi0, chi_Z=chi1, hb=_q(c, Q(6, 10)), sc=common["sc"], q=common["q"])
    except (TypeError, ValueError):
        checks["hb_above_half_rejected"] = True
    else:
        checks["hb_above_half_rejected"] = False

    try:
        build(c, inlet=inlet, chi=chi0, chi_Z=chi1, hb=common["hb"], sc=common["sc"], q=_box(c, -1, 2))
    except (TypeError, ValueError):
        checks["negative_q_rejected"] = True
    else:
        checks["negative_q_rejected"] = False

    try:
        build(c, inlet=inlet, chi=chi0, chi_Z=chi1, hb=common["hb"], sc=common["sc"], q=_box(c, Q(9, 10), Q(11, 10)))
    except (TypeError, ValueError):
        checks["q_above_one_rejected"] = True
    else:
        checks["q_above_one_rejected"] = False

    try:
        build(c, inlet=inlet, chi=dict(chi0, powers=[1, 0, 0, 0]), chi_Z=chi1, hb=common["hb"], sc=common["sc"], q=common["q"])
    except (TypeError, ValueError):
        checks["nonzero_chi_formal_power_rejected"] = True
    else:
        checks["nonzero_chi_formal_power_rejected"] = False

    no_id = {name: {order: [dict(row[0], source_ids=[])] for order, row in rows.items()} for name, rows in inlet.items()}
    try:
        build(c, inlet=no_id, chi=chi0, chi_Z=chi1, hb=common["hb"], sc=common["sc"], q=common["q"])
    except (TypeError, ValueError):
        checks["missing_input_source_id_rejected"] = True
    else:
        checks["missing_input_source_id_rejected"] = False

    bad_power = {name: {order: list(row) for order, row in rows.items()} for name, rows in inlet.items()}
    bad_power["F1"]["Z0"] = [dict(bad_power["F1"]["Z0"][0], powers=[1, 2, 3])]
    try:
        build(c, inlet=bad_power, chi=chi0, chi_Z=chi1, hb=common["hb"], sc=common["sc"], q=common["q"])
    except (TypeError, ValueError):
        checks["bad_powers_rejected"] = True
    else:
        checks["bad_powers_rejected"] = False

    try:
        build(c, inlet=inlet, chi=chi0, chi_Z=chi1, hb=common["hb"], sc=_q(c, 0), q=common["q"])
    except (TypeError, ValueError):
        checks["zero_sc_rejected"] = True
    else:
        checks["zero_sc_rejected"] = False
    return checks


def run_fixture(leaf: Any, c: Any) -> dict[str, Any]:
    """Run independent interval, lineage, power and guard checks."""
    build = leaf.build if hasattr(leaf, "build") else leaf
    if not callable(build):
        raise TypeError("leaf must be a build callable or expose build")
    inlet = _inlet(c)
    chi0, chi1 = _chi(c)
    common = {"hb": _q(c, Q(4, 17)), "sc": _box(c, Q(5, 7), Q(7, 4)), "q": _box(c, Q(3, 8), Q(1, 2))}
    atoms = {name: f"fixture:v2:E00:{name}" for name in SOURCE_POWER_BASIS}
    result = build(
        c,
        inlet=inlet,
        chi=chi0,
        chi_Z=chi1,
        hb=common["hb"],
        sc=common["sc"],
        q=common["q"],
        source_family="synthetic_fixed_E00_j2_v2",
        source_atoms=atoms,
    )
    checks: dict[str, bool] = {
        "schema": result.get("schema") == "n1_e00_collar_j2_volterra_v2",
        "four_branches": tuple(result.get("branch_names", ())) == BRANCHES,
        "explicit_orders": tuple(result.get("ordinary_Z_orders", ())) == (0, 1),
        "source_power_basis": tuple(result.get("source_power_basis", ())) == SOURCE_POWER_BASIS,
        "no_ledger": result.get("no_ledger") is True and "source_term_ledger" not in result,
        "all_terms_symmetric": all(_symmetric_interval(term["coefficient"]) for term in _terms(result)),
        "all_terms_error_only": all(term.get("error_only") is True and term.get("symmetric_error") is True for term in _terms(result)),
        "all_terms_tail_order": all(term.get("tail_starts_at_j") == 2 for term in _terms(result)),
        "chi_derivative_has_two_ids": all(
            len(term.get("chi_source_ids", ())) == 2
            for branch in BRANCHES
            for term in result["branches"][branch]["Z1"]
            if term.get("chi_order") == 1
        ),
        "chi_derivative_piece_separate": all(
            term.get("leibniz_piece") == "chi_Z*inlet_Z0"
            for branch in BRANCHES
            for term in result["branches"][branch]["Z1"]
            if term.get("chi_order") == 1
        ),
        "sentinel_id_preserved": any(SENTINEL_ID in term.get("input_source_ids", ()) for term in _terms(result)),
        "sentinel_class_preserved": any(SENTINEL_CLASS in term.get("source_error_classes", ()) for term in _terms(result)),
        "formal_hb_powers": all(
            result["branches"][branch]["Z0"][0]["powers"][0]
            == inlet["F1" if branch.endswith("Fc") else "Fx"]["Z0"][0]["powers"][0] + FORMAL_HB_POWER_ADDED[branch]
            for branch in BRANCHES
        ),
        "Fxc_bound_attached": all(
            result["branches"][branch]["Z0"][0].get("h_fxc_bound")
            for branch in ("F_from_Fxc", "Fx_from_Fxc")
        ),
        "geometry_upper_only": result["kernel"].get("formal_hb_not_materialized") is True and result["kernel"].get("hb_upper_for_exponent") == _q(c, Q(1, 2)),
    }

    zero_chi0, zero_chi1 = _chi(c, zero=True)
    zero_m = build(c, inlet=inlet, chi=zero_chi0, chi_Z=zero_chi1, hb=common["hb"], sc=common["sc"], q=common["q"])
    checks["M_zero_Mz_nonzero_is_zero"] = all(_zero_interval(term["coefficient"]) for term in _terms(zero_m))
    checks["M_zero_keeps_chi_Z_sources"] = any(
        term.get("chi_order") == 1 and "fixture:v2:E00:chi:Z1" in term.get("chi_source_ids", ())
        for term in _terms(zero_m)
    )

    seam = build(c, inlet=inlet, chi=chi0, chi_Z=chi1, hb=common["hb"], sc=common["sc"], q=_q(c, 0))
    checks["zero_seam_exact"] = seam["seam"]["componentwise_zero"] is True and all(_zero_interval(term["coefficient"]) for term in _terms(seam))
    outer = build(c, inlet=inlet, chi=chi0, chi_Z=chi1, hb=common["hb"], sc=common["sc"], q_cells=((0, Q(1, 4)), (Q(1, 4), Q(1, 2)), (Q(1, 2), 1)))
    checks["outer_q_cells"] = outer.get("cell_count") == 3 and len(outer.get("q_cells", ())) == 3
    checks["outer_cell_contains_seam"] = outer["q_cells"][0]["q_interval"].a == 0 and outer["q_cells"][0]["seam"]["componentwise_zero"] is False
    checks.update(_guard_checks(build, c, inlet, chi0, chi1, common))
    return {
        "passed": all(checks.values()),
        "checks": checks,
        "check_count": len(checks),
        "synthetic_only": True,
        "independent_exact_fraction_inputs": True,
        "directed_interval_arithmetic": True,
        "old_fixture_imported": False,
        "project_imports": False,
        "receipt_reads": False,
        "midpoint_values_used": False,
        "quadrature_used": False,
        "limits": {
            "tail": "only homogeneous Kchi words j>=2",
            "geometry": "L<=hb*B with hb<=1/2 used only for exponent and formal hb powers",
            "Z1": "chi_Z*input_Z0 is separate and retains both chi source IDs",
        },
    }


run = run_fixture

