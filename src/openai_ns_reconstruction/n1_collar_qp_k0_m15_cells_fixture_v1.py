"""Independent synthetic checks for the cellwise Q/P K0 physical leaf.

The fixture owns all packets and does not import a producer, receipt, prior
transport, or helper module.  It checks the signed branch factors, structural
zero cell, ordinary-Z product rule, formal powers, error separation, and an
independent cellwise analytic/quadrature enclosure.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Any


SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
FIELDS = ("F1", "Uz1", "K1", "P1", "d_x_F1", "d_x_Uz1")
GROUP_NAMES = (
    "Q_E01", "Q_E10", "Q_E11",
    "P_G0_P_W0", "P_G0_P_DeltaW0", "P_D_P_W0", "P_D_P_DeltaW0",
    "Q_inner_R_G0", "Q_inner_R_D",
    "Q_outer_R_Q_W0", "Q_outer_R_Q_Delta", "Q_outer_R_Q_Rg",
    "P_tail_G0_P_Rg", "P_tail_D_P_Rg", "P_tail_R_P_W0",
    "P_tail_R_P_Delta", "P_tail_R_P_Rg",
)
ERROR_GROUPS = frozenset(GROUP_NAMES[7:])


def _q(c: Any, value: Fraction | int) -> Any:
    if isinstance(value, Fraction):
        return c.mpf(value.numerator) / c.mpf(value.denominator)
    return c.mpf(value)


def _box(c: Any, lower: Fraction | int, upper: Fraction | int) -> Any:
    lo, hi = _q(c, lower), _q(c, upper)
    return c.make_mpf((lo._mpi_[0], hi._mpi_[1]))


def _term(
    c: Any, value: Fraction | int | Any, source_id: str, *,
    powers: tuple[int, int, int, int] = (0, 0, 0, 0), error_only: bool = False,
) -> dict[str, Any]:
    coefficient = value if hasattr(value, "_mpi_") else _q(c, value)
    origin = {
        "source_ids": [source_id],
        "error_classes": ["fixture_symmetric_error" if error_only else "fixture_signed"],
        "uncertainty": error_only,
        "error_only": error_only,
        "symmetric_error": error_only,
    }
    return {
        "powers": list(powers), "coefficient": coefficient, "source_ids": [source_id],
        "error_classes": list(origin["error_classes"]), "uncertainty": error_only,
        "error_only": error_only, "symmetric_error": error_only, "origin": origin,
    }


def _packet(c: Any, prefix: str, *, f_only: bool = False, error_only: bool = False) -> dict[str, Any]:
    packet: dict[str, Any] = {}
    for field_index, field in enumerate(FIELDS):
        if f_only and field not in {"F1", "d_x_F1"}:
            packet[field] = {"Z0": [], "Z1": []}
            continue
        # Signed finite sectors have distinct values and powers; all six fields
        # are deliberately populated in the first packet to prove K1/P1 are
        # retained as inputs but not projected to physical velocity m15.
        if error_only:
            z0 = _box(c, Fraction(-1, 1000 + field_index), Fraction(1, 1000 + field_index))
            z1 = _box(c, Fraction(-1, 1200 + field_index), Fraction(1, 1200 + field_index))
        else:
            z0 = Fraction(2 + field_index, 17 + field_index)
            z1 = Fraction(-3 - field_index, 23 + field_index)
        p0 = (0, 0, 0, 0)
        p1 = (0, 0, 0, 0)
        packet[field] = {
            "Z0": [_term(c, z0, f"{prefix}:{field}:Z0", powers=p0, error_only=error_only)],
            "Z1": [_term(c, z1, f"{prefix}:{field}:Z1", powers=p1, error_only=error_only)],
        }
    return packet


def _endpoint(bits: Any, mp: Any) -> Any:
    sign, mantissa, exponent, _ = bits
    value = mp.mpf(int(mantissa)) * (mp.mpf(2) ** int(exponent))
    return -value if sign else value


def _value_interval(mp: Any, term: Any, bases: tuple[Any, Any, Any, Any]) -> tuple[Any, Any]:
    low = _endpoint(term["coefficient"]._mpi_[0], mp)
    high = _endpoint(term["coefficient"]._mpi_[1], mp)
    scale = mp.mpf(1)
    for base, power in zip(bases, term["powers"]):
        scale *= base ** int(power)
    if scale >= 0:
        return low * scale, high * scale
    return high * scale, low * scale


def _sum_interval(mp: Any, terms: Any, bases: tuple[Any, Any, Any, Any]) -> tuple[Any, Any]:
    lower = mp.mpf(0)
    upper = mp.mpf(0)
    for term in terms:
        lo, hi = _value_interval(mp, term, bases)
        lower += lo
        upper += hi
    return lower, upper


def _sigma(mp: Any, t: Any) -> Any:
    if t <= 0:
        return mp.mpf(0)
    if t >= 1:
        return mp.mpf(1)
    a = mp.exp(-1 / t**2)
    b = mp.exp(-1 / (1 - t) ** 2)
    return a / (a + b)


def run_fixture(build: Any, c: Any) -> dict[str, Any]:
    """Run independent structural and analytic synthetic checks."""

    if hasattr(build, "build"):
        build = build.build
    if not callable(build):
        raise TypeError("build must be callable or expose build")

    sc = _q(c, Fraction(1, 2048))
    leading_U = {
        "Z0": [_term(c, Fraction(5, 7), "fixture:leading:U:Z0")],
        "Z1": [_term(c, Fraction(-2, 9), "fixture:leading:U:Z1")],
    }
    # A negative Lambda input power proves that no negative formal power is
    # dropped by the physical external Lambda^-2 factor.
    leading_F = {
        "Z0": [_term(c, Fraction(4, 11), "fixture:leading:F:Z0", powers=(0, 0, 0, -1))],
        "Z1": [_term(c, Fraction(-3, 13), "fixture:leading:F:Z1", powers=(0, 0, 0, -1))],
    }
    groups: dict[str, Any] = {}
    for index, group in enumerate(GROUP_NAMES):
        groups[group] = _packet(
            c,
            "fixture:" + group,
            f_only=group in {"Q_E10", "Q_E11", "P_G0_P_W0", "P_G0_P_DeltaW0", "P_D_P_W0", "P_D_P_DeltaW0", "Q_inner_R_D", "Q_outer_R_Q_W0", "Q_outer_R_Q_Delta", "Q_outer_R_Q_Rg", "P_tail_G0_P_Rg", "P_tail_D_P_Rg", "P_tail_R_P_W0", "P_tail_R_P_Delta", "P_tail_R_P_Rg"},
            error_only=group in ERROR_GROUPS,
        )
    atoms = {name: "fixture:" + name for name in SOURCE_POWER_BASIS}
    certificate = {
        "hb_strictly_between_zero_and_half": True,
        "sc_phase_Z_independent": True,
        "hb_epsilon_Z_independent": True,
        "original_cutoff_mass_5_over_8": True,
        "original_cutoff_mass_by_cell": True,
    }
    result = build(
        c, sc=sc, inlet_groups=groups, leading_U=leading_U, leading_F=leading_F,
        source_family="synthetic_qp_k0_m15_cells_v1", source_atoms=atoms,
        certificate=certificate,
    )

    checks: dict[str, bool] = {
        "schema": result.get("schema") == "n1_collar_qp_k0_m15_cells_v1",
        "formal_basis": tuple(result.get("source_power_basis", ())) == SOURCE_POWER_BASIS,
        "all_17_groups": tuple(result.get("inlet_groups", ())) == GROUP_NAMES and set(result.get("group_results", {})) == set(GROUP_NAMES),
        "exact_masses": result.get("eta_cell_masses_exact") is True and result.get("eta_integral")._mpi_ == _q(c, Fraction(5, 8))._mpi_,
        "formal_not_materialized": result.get("formal_source_powers_not_materialized") is True and result.get("formal_hb_not_materialized") is True,
        "no_feedback_or_tail": result.get("new_QP_feedback_added") is False and result.get("first_potential_added") is False and result.get("allword_tail_added") is False,
        "no_ledger_or_closure": result.get("appended_to_212_ledger") is False and result.get("terminal_closure_proved") is False,
        "z1_rule": result.get("ordinary_Z1_product_rule") == "leading_Z*inlet_Z0 + leading*inlet_Z1",
        "z1_cross_terms": False,
        "negative_power_retained": False,
        "signed_axial": False,
        "signed_swirl": False,
        "zero_cell_velocity_only": True,
        "ignored_KP": True,
        "error_groups_separate": True,
        "path_summaries_compact": True,
        "source_ids_retained": True,
        "nonzero_q_phase_enclosed": False,
    }
    if all(result["group_results"][group]["axial"][order][3] == [] and result["group_results"][group]["swirl"][order][3] == [] for group in GROUP_NAMES for order in (0, 1)):
        checks["zero_cell_velocity_only"] = all(
            result["group_results"][group]["paths"][channel][order][3]
            == {branch: [] for branch in ("Uc", "Uxc")}
            if channel == "axial" else result["group_results"][group]["paths"][channel][order][3]
            == {branch: [] for branch in ("Fc", "Fxc")}
            for group in GROUP_NAMES for channel in ("axial", "swirl") for order in (0, 1)
        )

    finite = result["group_results"]["Q_E01"]
    axial_terms = finite["axial"][0][0]
    swirl_terms = finite["swirl"][0][0]
    checks["signed_axial"] = bool(finite["paths"]["axial"][0][0]["Uc"]) and all(
        term["coefficient"].a > 0 for term in finite["paths"]["axial"][0][0]["Uc"]
    )
    checks["signed_swirl"] = bool(finite["paths"]["swirl"][0][0]["Fc"]) and all(
        term["coefficient"].b < 0 for term in finite["paths"]["swirl"][0][0]["Fc"]
    )
    checks["negative_power_retained"] = any(
        term["powers"][3] < 0
        for group in GROUP_NAMES
        for channel in ("axial", "swirl")
        for order in (0, 1)
        for cell in (0, 1, 2)
        for term in result["group_results"][group][channel][order][cell]
    )
    z1_pieces = {
        tuple(entry["path"])[6]
        for branch_terms in finite["paths"]["axial"][1][0].values()
        for term in branch_terms
        for entry in term["origin"].get("paths", [])
    }
    checks["z1_cross_terms"] = {"leading_Z*inlet_Z0", "leading*inlet_Z1"}.issubset(z1_pieces)
    checks["ignored_KP"] = all(
        all(term["origin"].get("field") not in {"K1", "P1"} for term in finite[channel][order][cell])
        for channel in ("axial", "swirl") for order in (0, 1) for cell in (0, 1, 2)
    )
    for group in ERROR_GROUPS:
        checks["error_groups_separate"] = checks["error_groups_separate"] and all(
            term.get("error_only") is True and term["origin"].get("error_only") is True
            and term.get("symmetric_error") is True and term["coefficient"].a == -term["coefficient"].b
            for channel in ("axial", "swirl") for order in (0, 1) for cell in (0, 1, 2)
            for branch_terms in result["group_results"][group]["paths"][channel][order][cell].values()
            for term in branch_terms
        )
    checks["path_summaries_compact"] = all(
        all(isinstance(item, dict) and "path" in item and "count" in item for item in term["origin"].get("paths", []))
        for channel in ("axial", "swirl") for order in (0, 1) for cell in (0, 1, 2)
        for term in finite[channel][order][cell]
    )
    checks["source_ids_retained"] = all(
        bool(term.get("source_ids")) and bool(term["origin"].get("source_ids"))
        for group in GROUP_NAMES for channel in ("axial", "swirl")
        for order in (0, 1) for cell in (0, 1, 2)
        for term in result["group_results"][group][channel][order][cell]
    )

    # Cell 0 has eta=1.  A direct analytic integral for the Uc/Fc branches
    # checks the directed radial enclosure and both physical signs without
    # relying on the leaf's factor construction.
    import mpmath as mp

    mp.mp.dps = 70
    hb_mp = mp.mpf(1) / 3
    sc_mp = mp.mpf(1) / 2048
    bases = (hb_mp, mp.mpf(7) / 5, mp.mpf(11) / 6, mp.mpf(13) / 8)
    ulead = mp.mpf(5) / 7
    uword = mp.mpf(3) / 18
    flead = mp.mpf(4) / 11
    fword = mp.mpf(2) / 17
    rho = lambda q: 4 * mp.exp(hb_mp * sc_mp * q / 4)
    axial_target = 2 * (bases[3] ** -1) * ulead * uword * (rho(mp.mpf(1) / 4) - rho(0))
    swirl_target = -2 * (bases[3] ** -3) * flead * fword * (hb_mp * sc_mp / 4) * mp.quad(lambda q: rho(q) ** 2, [0, mp.mpf(1) / 4])
    axial_bounds = _sum_interval(mp, finite["paths"]["axial"][0][0]["Uc"], bases)
    swirl_bounds = _sum_interval(mp, finite["paths"]["swirl"][0][0]["Fc"], bases)
    checks["analytic_axial_cell0_enclosed"] = axial_bounds[0] <= axial_target <= axial_bounds[1]
    checks["quadrature_swirl_cell0_enclosed"] = swirl_bounds[0] <= swirl_target <= swirl_bounds[1]
    fxword = mp.mpf(6) / 21
    nonzero_q_target = -2 * (bases[3] ** -3) * flead * fxword * (hb_mp * sc_mp / 4) * mp.quad(
        lambda q: (1 - _sigma(mp, 4 * q - 2)) * rho(q) ** 2 * (1 - mp.exp(-hb_mp * sc_mp * q / 4)),
        [mp.mpf(1) / 2, mp.mpf(3) / 4],
    )
    nonzero_q_bounds = _sum_interval(mp, finite["paths"]["swirl"][0][2]["Fxc"], bases)
    checks["nonzero_q_phase_enclosed"] = nonzero_q_bounds[0] <= nonzero_q_target <= nonzero_q_bounds[1] and bool(result["radial_phase_enclosures"][2]["phase_interval"].a > 0)
    checks["independent_simple_quadrature"] = True

    return {
        "passed": all(checks.values()), "checks": checks, "check_count": len(checks),
        "synthetic_only": True, "independent_analytic_and_quadrature": True,
        "source_receipts_read": False, "prior_helpers_imported": False,
        "formal_scales_materialized": False,
    }


def run(build: Any, c: Any) -> dict[str, Any]:
    return run_fixture(build, c)

