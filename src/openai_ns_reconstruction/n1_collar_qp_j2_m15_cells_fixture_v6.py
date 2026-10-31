"""Independent synthetic checks for the Q/P Kchi-tail physical m15 leaf.

The fixture owns all sectors and does not import a producer, receipt, prior
transport, or helper module.  It checks the four-cell physical factors,
formal powers, ordinary-Z routes, chi provenance, symmetric error semantics,
and independent high-precision positive ``j=2`` target integrals.
"""

from __future__ import annotations

import json
from fractions import Fraction
from typing import Any


SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
GROUP_NAMES = (
    "Q_E01", "Q_E10", "Q_E11",
    "P_G0_P_W0", "P_G0_P_DeltaW0", "P_D_P_W0", "P_D_P_DeltaW0",
    "Q_inner_R_G0", "Q_inner_R_D",
    "Q_outer_R_Q_W0", "Q_outer_R_Q_Delta", "Q_outer_R_Q_Rg",
    "P_tail_G0_P_Rg", "P_tail_D_P_Rg", "P_tail_R_P_W0",
    "P_tail_R_P_Delta", "P_tail_R_P_Rg",
)
ERROR_GROUPS = frozenset(GROUP_NAMES[7:])
BRANCHES = ("F_from_Fc", "F_from_Fxc")


def _q(c: Any, value: Fraction | int) -> Any:
    if isinstance(value, Fraction):
        return c.mpf(value.numerator) / c.mpf(value.denominator)
    return c.mpf(value)


def _box(c: Any, lower: Fraction | int, upper: Fraction | int) -> Any:
    lo, hi = _q(c, lower), _q(c, upper)
    return c.make_mpf((lo._mpi_[0], hi._mpi_[1]))


def _term(
    c: Any, value: Fraction | int | Any, source_id: str, *,
    powers: tuple[int, int, int, int] = (0, 0, 0, 0),
    error_only: bool = False, error_class: str = "fixture_qp_j2_signed",
) -> dict[str, Any]:
    coefficient = value if hasattr(value, "_mpi_") else _q(c, value)
    classes = [error_class] + (["symmetric_error"] if error_only else [])
    origin = {
        "source_ids": [source_id], "error_classes": classes,
        "uncertainty": error_only, "error_only": error_only,
        "symmetric_error": error_only,
    }
    return {
        "powers": list(powers), "coefficient": coefficient,
        "source_ids": [source_id], "error_classes": classes,
        "uncertainty": error_only, "error_only": error_only,
        "symmetric_error": error_only, "origin": origin,
    }


def _packet(c: Any, prefix: str, index: int, *, error_only: bool = False, negative_powers: bool = False) -> dict[str, Any]:
    if error_only:
        fc0 = _box(c, Fraction(-1, 1000 + index), Fraction(1, 1000 + index))
        fc1 = _box(c, Fraction(-1, 1100 + index), Fraction(1, 1100 + index))
        fxc0 = _box(c, Fraction(-1, 1300 + index), Fraction(1, 1300 + index))
        fxc1 = _box(c, Fraction(-1, 1400 + index), Fraction(1, 1400 + index))
    else:
        fc0 = Fraction(23 + index, 41 + index)
        fc1 = Fraction(-19 - index, 59 + index)
        fxc0 = Fraction(-17 - index, 101 + index)
        fxc1 = Fraction(21 + index, 107 + index)
    fxc_power = (-2, 1, -1, 0) if negative_powers else (0, 0, 0, 0)
    suffix = "error" if error_only else "signed"
    return {
        "Fc": {
            "Z0": [_term(c, fc0, f"{prefix}:Fc:Z0", error_only=error_only, error_class=f"{prefix}:{suffix}")],
            "Z1": [_term(c, fc1, f"{prefix}:Fc:Z1", error_only=error_only, error_class=f"{prefix}:{suffix}")],
        },
        "Fxc": {
            "Z0": [_term(c, fxc0, f"{prefix}:Fxc:Z0", powers=fxc_power, error_only=error_only, error_class=f"{prefix}:{suffix}")],
            "Z1": [_term(c, fxc1, f"{prefix}:Fxc:Z1", powers=fxc_power, error_only=error_only, error_class=f"{prefix}:{suffix}")],
        },
    }


def _term_paths(term: Any) -> list[dict[str, Any]]:
    origin = term.get("origin", {})
    if origin.get("path") is not None:
        return [origin["path"]]
    return [entry.get("path", {}) for entry in origin.get("paths", ()) if isinstance(entry, dict)]


def _bounds_mp(c: Any, value: Any, mp: Any) -> tuple[Any, Any]:
    raw = value._mpi_
    return mp.make_mpf(raw[0]), mp.make_mpf(raw[1])


def _formal_interval(c: Any, term: Any, bases: tuple[Any, Any, Any, Any], mp: Any) -> tuple[Any, Any]:
    lo, hi = _bounds_mp(c, term["coefficient"], mp)
    scale = mp.mpf(1)
    for base, power in zip(bases, term["powers"]):
        scale *= base ** int(power)
    if scale >= 0:
        return lo * scale, hi * scale
    return hi * scale, lo * scale


def _find_path_terms(result: Any, group: str, order: int, cell: int, branch: str, piece: str) -> list[dict[str, Any]]:
    values = result["group_results"][group]["paths"]["swirl"][order][cell][branch]
    return [term for term in values if any(path.get("ordinary_Z_piece") == piece for path in _term_paths(term))]


def _gauss8(mp: Any, fn: Any, lower: Any, upper: Any) -> Any:
    """Bounded 8-node Gauss-Legendre rule generated at run time."""
    nodes, weights = mp.gauss_quadrature(8, "legendre")
    midpoint = (lower + upper) / 2
    half_width = (upper - lower) / 2
    total = mp.mpf(0)
    for node, weight in zip(nodes, weights):
        total += weight * fn(midpoint + half_width * node)
    return half_width * total


def _j2_outer(mp: Any, *, hb: Any, sc: Any, target: Any, branch: str = "F_from_Fc") -> Any:
    """Closed-form nested positive K0 word at target q.

    With u=hb*sc*q/4, the physical K0 measure is
    exp(u)*(1-exp(u-r))*du.  The expression below is the exact nested
    j=2 value; it is used only as a high-precision synthetic cross-check.
    """
    r = hb * sc * target / 4
    if branch == "F_from_Fxc":
        return (
            mp.exp(2 * r) / 12
            + (mp.mpf(3) / 4 - r / 2) * mp.exp(r)
            - r / 2
            - mp.mpf(3) / 4
            - mp.exp(-r) / 12
        )
    if branch != "F_from_Fc":
        raise ValueError(f"unknown branch {branch!r}")
    return (
        mp.exp(2 * r) / 12
        - mp.exp(r) / 2
        + mp.exp(-r) / 6
        + mp.mpf(1) / 4
        + r / 2
    )


def _j2_outer_numeric(mp: Any, *, hb: Any, sc: Any, target: Any, branch: str) -> Any:
    """Independent dynamic 8-node target check for the closed-form masses."""
    phase = sc / 8

    def x(q: Any) -> Any:
        return 2 * mp.exp(hb * phase * q)

    def dx(q: Any) -> Any:
        return hb * phase * x(q)

    def k0(X: Any, t: Any) -> Any:
        return t * (1 - (t / X) ** 2) / 2

    X = x(target)

    def inner(t: Any) -> Any:
        T = x(t)

        def source_shape(s: Any) -> Any:
            if branch == "F_from_Fc":
                return mp.mpf(1)
            return 1 - (2 / x(s)) ** 2

        return _gauss8(
            mp,
            lambda s: k0(T, x(s)) * source_shape(s) * dx(s),
            mp.mpf(0), t,
        )

    return _gauss8(
        mp,
        lambda t: k0(X, x(t)) * dx(t) * inner(t),
        mp.mpf(0), target,
    )


def run_fixture(leaf: Any, c: Any) -> dict[str, Any]:
    """Run independent structural and numerical containment checks."""

    if hasattr(leaf, "build"):
        leaf = leaf.build
    if not callable(leaf):
        raise TypeError("leaf must be callable or expose build")

    sc = _q(c, Fraction(1, 65536))
    leading_F = {
        "Z0": [_term(c, Fraction(25, 48), "fixture_v6:leading:F:Z0")],
        "Z1": [_term(c, Fraction(-23, 103), "fixture_v6:leading:F:Z1")],
    }
    chi = {
        "Z0": [_term(c, Fraction(29, 83), "fixture_v6:chi:Z0", error_class="fixture_chi_Z0")],
        "Z1": [_term(c, Fraction(-31, 97), "fixture_v6:chi:Z1", error_class="fixture_chi_Z1")],
    }
    groups = {
        group: _packet(
            c, "fixture_v6:" + group, index,
            error_only=group in ERROR_GROUPS,
            negative_powers=group == "Q_E10",
        )
        for index, group in enumerate(GROUP_NAMES)
    }
    atoms = {name: "fixture_v6:" + name for name in SOURCE_POWER_BASIS}
    certificate = {
        "hb_strictly_between_zero_and_half": True,
        "sc_phase_Z_independent": True,
        "hb_epsilon_Z_independent": True,
        "original_cutoff_mass_by_cell": True,
    }
    result = leaf(
        c, sc=sc, inlet_groups=groups, leading_F=leading_F, chi=chi,
        source_family="synthetic_qp_j2_m15_cells_v6_distinct",
        source_atoms=atoms, certificate=certificate,
    )

    checks: dict[str, bool] = {
        "schema": result.get("schema") == "n1_collar_qp_j2_m15_cells_v6",
        "formal_basis": tuple(result.get("source_power_basis", ())) == SOURCE_POWER_BASIS,
        "all_17_groups": tuple(result.get("inlet_groups", ())) == GROUP_NAMES and set(result.get("group_results", {})) == set(GROUP_NAMES),
        "four_cells": tuple(item["cell"] for item in result.get("q_cells", ())) == (0, 1, 2, 3),
        "exact_eta_masses": result.get("eta_cell_masses_exact") is True and all(
            result["q_cells"][index]["eta_mass"]._mpi_ == _q(c, Fraction(*mass))._mpi_
            for index, mass in enumerate(((1, 4), (1, 4), (1, 8), (0, 1)))
        ),
        "formal_not_materialized": result.get("formal_powers_not_materialized") is True,
        "tail_only_flags": result.get("first_potential_added") is False and result.get("feedback_tail_included") is False and result.get("allword_tail_included") is False,
        "no_ledger": result.get("appended_to_212_ledger") is False and result.get("terminal_closure_proved") is False,
        "z1_rule": result.get("ordinary_Z1_product_rule") == "leading_Z*tail_Z0 + leading*tail_inlet_Z1 + leading*tail_chi_Z1",
        "external_physical_powers": tuple(result.get("formal_external_powers", ())) == (1, 0, 0, -2),
        "no_signed_values": result.get("signed_values_returned") is False,
        "zero_cell_velocity_only": True,
        "all_outputs_symmetric_error": True,
        "source_paths_retained": True,
        "distinct_z1_routes": True,
        "chi_z_both_ids": True,
        "independent_quadrature": False,
    }

    for group in GROUP_NAMES:
        for order in (0, 1):
            for cell in (0, 1, 2):
                for branch in BRANCHES:
                    terms = result["group_results"][group]["paths"]["swirl"][order][cell][branch]
                    checks["all_outputs_symmetric_error"] = checks["all_outputs_symmetric_error"] and all(
                        term.get("error_only") is True
                        and term.get("symmetric_error") is True
                        and term["coefficient"].a == -term["coefficient"].b
                        and term["coefficient"].a <= 0 <= term["coefficient"].b
                        for term in terms
                    )
                    checks["source_paths_retained"] = checks["source_paths_retained"] and all(
                        term.get("source_ids") and term.get("incoming_source_ids")
                        and term.get("leading_source_ids") and term.get("chi_source_ids")
                        and _term_paths(term)
                        for term in terms
                    )
            for branch in BRANCHES:
                checks["zero_cell_velocity_only"] = checks["zero_cell_velocity_only"] and not result["group_results"][group]["paths"]["swirl"][order][3][branch]

    for group in ERROR_GROUPS:
        checks[f"{group}:error_only"] = all(
            term.get("error_only") is True and term.get("symmetric_error") is True
            for order in (0, 1) for cell in (0, 1, 2) for branch in BRANCHES
            for term in result["group_results"][group]["paths"]["swirl"][order][cell][branch]
        )

    finite = result["group_results"]["Q_E01"]
    route_pieces = {
        path.get("ordinary_Z_piece")
        for branch_terms in finite["paths"]["swirl"][1][0].values()
        for term in branch_terms
        for path in _term_paths(term)
    }
    checks["distinct_z1_routes"] = {
        "leading_Z*tail_Z0", "leading*tail_inlet_Z1", "leading*tail_chi_Z1",
    }.issubset(route_pieces)
    chi_route_terms = [
        term for branch in BRANCHES
        for term in _find_path_terms(result, "Q_E01", 1, 0, branch, "leading*tail_chi_Z1")
    ]
    checks["chi_z_both_ids"] = bool(chi_route_terms) and all(
        set(term.get("chi_source_ids", ())) == {"fixture_v6:chi:Z0", "fixture_v6:chi:Z1"}
        and set(term.get("chi_source_classes", ())) == {"fixture_chi_Z0", "fixture_chi_Z1"}
        for term in chi_route_terms
    )

    incoming_power = {}
    for group_packet in groups.values():
        for component in ("Fc", "Fxc"):
            for order in ("Z0", "Z1"):
                for term in group_packet[component][order]:
                    for source_id in term["source_ids"]:
                        incoming_power[source_id] = term["powers"][0]
    leading_power = {source_id: term["powers"][0] for order in ("Z0", "Z1") for term in leading_F[order] for source_id in term["source_ids"]}
    checks["formal_hb_increments"] = all(
        term["powers"][0] - incoming_power.get(term["incoming_source_ids"][0], 0) - leading_power.get(term["leading_source_ids"][0], 0)
        in (5, 6)
        for group in GROUP_NAMES for order in (0, 1) for cell in (0, 1, 2) for branch in BRANCHES
        for term in result["group_results"][group]["paths"]["swirl"][order][cell][branch]
    )
    checks["negative_input_power_preserved"] = any(
        term["powers"][1] != 0 or term["powers"][2] != 0 or term["powers"][3] < -2
        for group in GROUP_NAMES for cell in (0, 1, 2) for branch in BRANCHES
        for term in result["group_results"][group]["paths"]["swirl"][0][cell][branch]
    )

    # Independent positive j=2 nested integral on cell 0.  The closed-form
    # nested mass already carries its hb^4 behavior; only the physical hb and
    # Lambda scale is materialised in this synthetic check.
    import mpmath as mp

    mp.mp.dps = 100
    hb_mp = mp.mpf(5) / 17
    sc_mp = mp.mpf(1) / 65536
    chi0_mp = mp.mpf(29) / 83
    chi1_mp = mp.mpf(-31) / 97
    fc_mp = mp.mpf(23) / 41
    flead_mp = mp.mpf(25) / 48
    lambda_mp = mp.mpf(17) / 9

    q_probe = mp.mpf(5) / 32
    closed_fc = _j2_outer(mp, hb=hb_mp, sc=sc_mp, target=q_probe, branch="F_from_Fc")
    numeric_fc = _j2_outer_numeric(mp, hb=hb_mp, sc=sc_mp, target=q_probe, branch="F_from_Fc")
    closed_fxc = _j2_outer(mp, hb=hb_mp, sc=sc_mp, target=q_probe, branch="F_from_Fxc")
    numeric_fxc = _j2_outer_numeric(mp, hb=hb_mp, sc=sc_mp, target=q_probe, branch="F_from_Fxc")
    fc_relative_difference = abs(closed_fc - numeric_fc) / max(abs(closed_fc), abs(numeric_fc), mp.mpf("1e-90"))
    fxc_relative_difference = abs(closed_fxc - numeric_fxc) / max(abs(closed_fxc), abs(numeric_fxc), mp.mpf("1e-90"))
    checks["closed_form_fc_matches_gauss8"] = fc_relative_difference <= mp.mpf("1e-30")
    checks["closed_form_fxc_matches_gauss8"] = fxc_relative_difference <= mp.mpf("1e-30")
    quadrature_crosscheck = {
        "fc": {
            "closed": mp.nstr(closed_fc, 80),
            "numeric": mp.nstr(numeric_fc, 80),
            "relative_difference": mp.nstr(fc_relative_difference, 80),
        },
        "fxc": {
            "closed": mp.nstr(closed_fxc, 80),
            "numeric": mp.nstr(numeric_fxc, 80),
            "relative_difference": mp.nstr(fxc_relative_difference, 80),
        },
        "rule": "mp.gauss_quadrature(8, \"legendre\")",
        "relative_tolerance": "1e-30",
    }

    def rho(q: Any) -> Any:
        return 4 * mp.exp(hb_mp * sc_mp * q / 4)

    def tail_fc_at(q: Any, chi_value: Any, source_value: Any) -> Any:
        return 4 * chi_value * chi_value * source_value * _j2_outer(
            mp, hb=hb_mp, sc=sc_mp, target=q,
        )

    def tail_fxc_at(q: Any, chi_value: Any, source_value: Any) -> Any:
        return 4 * chi_value * chi_value * source_value * _j2_outer(
            mp, hb=hb_mp, sc=sc_mp, target=q, branch="F_from_Fxc",
        )

    fc1_mp = mp.mpf(-19) / 59
    fxc0_mp = mp.mpf(-17) / 101
    flead_z_mp = mp.mpf(-23) / 103
    physical_scale = hb_mp * lambda_mp**-2 * 8 * sc_mp
    target_z0 = physical_scale * flead_mp * mp.quad(
        lambda q: rho(q) ** 2 / 16 * tail_fc_at(q, chi0_mp, fc_mp), [0, mp.mpf(1) / 4],
    )
    target_fxc = physical_scale * flead_mp * mp.quad(
        lambda q: rho(q) ** 2 / 16 * tail_fxc_at(q, chi0_mp, fxc0_mp), [0, mp.mpf(1) / 4],
    )
    target_z1_leading = physical_scale * flead_z_mp * mp.quad(
        lambda q: rho(q) ** 2 / 16 * tail_fc_at(q, chi0_mp, fc_mp), [0, mp.mpf(1) / 4],
    )
    target_z1_inlet = physical_scale * flead_mp * mp.quad(
        lambda q: rho(q) ** 2 / 16 * tail_fc_at(q, chi0_mp, fc1_mp), [0, mp.mpf(1) / 4],
    )
    z0_terms = _find_path_terms(result, "Q_E01", 0, 0, "F_from_Fc", "leading*tail_Z0")
    z0_bounds = [
        _formal_interval(c, term, (hb_mp, mp.mpf(1), mp.mpf(1), lambda_mp), mp)
        for term in z0_terms
    ]
    lower_z0 = sum((item[0] for item in z0_bounds), mp.mpf(0))
    upper_z0 = sum((item[1] for item in z0_bounds), mp.mpf(0))
    checks["j2_z0_target_enclosed"] = bool(z0_terms) and lower_z0 <= target_z0 <= upper_z0

    def tail_z_at(q: Any) -> Any:
        return 8 * chi0_mp * chi1_mp * fc_mp * _j2_outer(
            mp, hb=hb_mp, sc=sc_mp, target=q,
        )

    target_z1_chi = physical_scale * flead_mp * mp.quad(
        lambda q: rho(q) ** 2 / 16 * tail_z_at(q), [0, mp.mpf(1) / 4],
    )
    z1_terms = _find_path_terms(result, "Q_E01", 1, 0, "F_from_Fc", "leading*tail_chi_Z1")
    z1_bounds = [
        _formal_interval(c, term, (hb_mp, mp.mpf(1), mp.mpf(1), lambda_mp), mp)
        for term in z1_terms
    ]
    lower_z1 = sum((item[0] for item in z1_bounds), mp.mpf(0))
    upper_z1 = sum((item[1] for item in z1_bounds), mp.mpf(0))
    checks["j2_chi_z1_target_enclosed"] = bool(z1_terms) and lower_z1 <= target_z1_chi <= upper_z1

    fxc_terms = _find_path_terms(result, "Q_E01", 0, 0, "F_from_Fxc", "leading*tail_Z0")
    fxc_bounds = [
        _formal_interval(c, term, (hb_mp, mp.mpf(1), mp.mpf(1), lambda_mp), mp)
        for term in fxc_terms
    ]
    lower_fxc = sum((item[0] for item in fxc_bounds), mp.mpf(0))
    upper_fxc = sum((item[1] for item in fxc_bounds), mp.mpf(0))
    checks["j2_fxc_target_enclosed"] = bool(fxc_terms) and lower_fxc <= target_fxc <= upper_fxc

    leading_z_terms = _find_path_terms(result, "Q_E01", 1, 0, "F_from_Fc", "leading_Z*tail_Z0")
    leading_z_bounds = [
        _formal_interval(c, term, (hb_mp, mp.mpf(1), mp.mpf(1), lambda_mp), mp)
        for term in leading_z_terms
    ]
    lower_leading_z = sum((item[0] for item in leading_z_bounds), mp.mpf(0))
    upper_leading_z = sum((item[1] for item in leading_z_bounds), mp.mpf(0))
    checks["j2_leading_z1_target_enclosed"] = bool(leading_z_terms) and lower_leading_z <= target_z1_leading <= upper_leading_z

    inlet_z_terms = _find_path_terms(result, "Q_E01", 1, 0, "F_from_Fc", "leading*tail_inlet_Z1")
    inlet_z_bounds = [
        _formal_interval(c, term, (hb_mp, mp.mpf(1), mp.mpf(1), lambda_mp), mp)
        for term in inlet_z_terms
    ]
    lower_inlet_z = sum((item[0] for item in inlet_z_bounds), mp.mpf(0))
    upper_inlet_z = sum((item[1] for item in inlet_z_bounds), mp.mpf(0))
    checks["j2_inlet_z1_target_enclosed"] = bool(inlet_z_terms) and lower_inlet_z <= target_z1_inlet <= upper_inlet_z
    checks["independent_quadrature"] = all(
        checks[name] for name in (
            "j2_z0_target_enclosed", "j2_chi_z1_target_enclosed",
            "j2_fxc_target_enclosed", "j2_leading_z1_target_enclosed",
            "j2_inlet_z1_target_enclosed",
        )
    )

    return {
        "passed": all(checks.values()), "checks": checks, "check_count": len(checks),
        "quadrature_crosscheck_json": json.dumps(quadrature_crosscheck, sort_keys=True),
        "synthetic_only": True, "independent_nested_j2": True,
        "independent_chi_z1_reference": True, "true_reference_containment": True,
        "source_receipts_read": False, "prior_helpers_imported": False,
        "formal_scales_materialized": False, "automatic_execution": False,
    }


def run(leaf: Any, c: Any) -> dict[str, Any]:
    return run_fixture(leaf, c)


__all__ = ["SOURCE_POWER_BASIS", "GROUP_NAMES", "ERROR_GROUPS", "BRANCHES", "run_fixture", "run"]












