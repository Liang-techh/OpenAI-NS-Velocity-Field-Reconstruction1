"""Independent synthetic checks for the Q/P first-potential m15 cell leaf.

The fixture owns every synthetic sector and does not import a producer,
receipt, prior leaf, or helper module.  It checks the variable-target mass
bounds against direct analytic and high-precision quadrature values on a
nonzero-q cell, all three ordinary-Z1 product pieces, formal powers, source
provenance, ignored fields, and symmetric error separation.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Any


SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
FIELDS = ("F1", "Uz1", "K1", "P1", "d_x_F1", "d_x_Uz1")
BRANCHES = ("F_from_Fc", "F_from_Fxc")
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
    c: Any,
    value: Fraction | int | Any,
    source_id: str,
    *,
    powers: tuple[int, int, int, int] = (0, 0, 0, 0),
    error_only: bool = False,
    uncertainty: bool = False,
    error_class: str = "fixture_signed",
) -> dict[str, Any]:
    coefficient = value if hasattr(value, "_mpi_") else _q(c, value)
    classes = [error_class]
    if error_only:
        classes.append("symmetric_error")
    origin = {
        "source_ids": [source_id],
        "error_classes": classes,
        "uncertainty": bool(uncertainty or error_only),
        "error_only": error_only,
        "symmetric_error": error_only,
    }
    return {
        "powers": list(powers),
        "coefficient": coefficient,
        "source_ids": [source_id],
        "error_classes": classes,
        "uncertainty": bool(uncertainty or error_only),
        "error_only": error_only,
        "symmetric_error": error_only,
        "origin": origin,
    }


def _packet(c: Any, prefix: str, *, error_only: bool = False) -> dict[str, Any]:
    packet: dict[str, Any] = {}
    for index, field in enumerate(FIELDS):
        if error_only:
            z0 = _box(c, Fraction(-1, 1000 + index), Fraction(1, 1000 + index))
            z1 = _box(c, Fraction(-1, 1300 + index), Fraction(1, 1300 + index))
        else:
            z0 = Fraction(2 + index, 17 + index)
            z1 = Fraction(-3 - index, 23 + index)
            if field == "d_x_F1":
                z0 = -z0
        packet[field] = {
            "Z0": [_term(c, z0, f"{prefix}:{field}:Z0", error_only=error_only, error_class="fixture_error" if error_only else "fixture_signed")],
            "Z1": [_term(c, z1, f"{prefix}:{field}:Z1", error_only=error_only, error_class="fixture_error" if error_only else "fixture_signed")],
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
    return (low * scale, high * scale) if scale >= 0 else (high * scale, low * scale)


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


def _eta(mp: Any, q: Any) -> Any:
    return 1 - _sigma(mp, 4 * q - 2)


def _mass(mp: Any, hb: Any, sc: Any, q: Any, branch: str) -> Any:
    y = hb * sc * q / 8
    if branch == "F_from_Fc":
        return mp.cosh(2 * y) - 1
    return mp.sinh(2 * y) - 2 * y


def _target(
    mp: Any,
    *,
    hb: Any,
    sc: Any,
    lam: Any,
    q_lower: Any,
    q_upper: Any,
    leading: Any,
    chi: Any,
    inlet: Any,
    branch: str,
) -> Any:
    integral = mp.quad(
        lambda q: _eta(mp, q) * (4 * mp.exp(hb * sc * q / 4)) ** 2 * _mass(mp, hb, sc, q, branch),
        [q_lower, q_upper],
    )
    return hb * sc * lam ** -2 * leading * chi * inlet * integral


def run_fixture(build: Any, c: Any) -> dict[str, Any]:
    """Run independent structural, analytic, and quadrature checks."""

    if hasattr(build, "build"):
        build = build.build
    if not callable(build):
        raise TypeError("build must be callable or expose build")

    sc = _q(c, Fraction(1, 2048))
    leading_f = {
        "Z0": [_term(c, Fraction(3, 7), "fixture:leading_F:Z0", powers=(0, 0, 0, -1))],
        "Z1": [_term(c, Fraction(-5, 11), "fixture:leading_F:Z1", powers=(0, 0, 0, -1))],
    }
    chi = {
        "Z0": [_term(c, Fraction(2, 5), "fixture:chi:Z0", uncertainty=True, error_class="fixture_chi")],
        "Z1": [_term(c, Fraction(-1, 6), "fixture:chi:Z1", uncertainty=True, error_class="fixture_chi")],
    }
    groups = {
        group: _packet(c, "fixture:" + group, error_only=group in ERROR_GROUPS)
        for group in GROUP_NAMES
    }
    atoms = {name: "fixture:" + name for name in SOURCE_POWER_BASIS}
    certificate = {
        "hb_strictly_between_zero_and_half": True,
        "strict_hb_admitted": True,
        "sc_phase_Z_independent": True,
        "hb_sc_Z_independent": True,
        "hb_epsilon_Z_independent": True,
        "rho_radius_Z_independent": True,
        "original_cutoff_mass_5_over_8": True,
        "original_cutoff_mass_by_cell": True,
    }
    result = build(
        c,
        sc=sc,
        inlet_groups=groups,
        leading_F=leading_f,
        chi=chi,
        source_family="synthetic_qp_j1_m15_cells_v1",
        source_atoms=atoms,
        certificate=certificate,
    )

    checks: dict[str, bool] = {
        "schema": result.get("schema") == "n1_collar_qp_j1_m15_cells_v1",
        "formal_basis": tuple(result.get("source_power_basis", ())) == SOURCE_POWER_BASIS,
        "all_17_groups": tuple(result.get("inlet_groups", ())) == GROUP_NAMES and set(result.get("group_results", {})) == set(GROUP_NAMES),
        "exact_cell_masses": all(
            result["eta_cell_masses"][index]._mpi_ == _q(c, Fraction(*pair))._mpi_
            for index, pair in enumerate(((1, 4), (1, 4), (1, 8), (0, 1)))
        ),
        "formal_not_materialized": result.get("formal_source_powers_not_materialized") is True and result.get("formal_hb_not_materialized") is True,
        "variable_target": result.get("target_q_is_variable") is True and result.get("endpoint_profile_reused") is False,
        "only_F_swirl": result.get("only_F_swirl_added") is True and result.get("axial_added") is False and result.get("pressure_added") is False,
        "physical_measure_once": result.get("physical_measure_integrated_once") is True,
        "ordinary_z1_rule": "leading_Z*chi_Z*inlet_Z0" not in result.get("ordinary_Z1_product_rule", "") and "leading_Z*chi*inlet_Z0" in result.get("ordinary_Z1_product_rule", ""),
        "positive_mass_bounds": True,
        "mass_h0f_factor": "1-4/t^2" in result.get("H0F_formula", "") and "1-4/t^2" in result["radial_phase_enclosures"][2]["mass_formulas"]["F_from_Fxc"],
        "zero_cell_velocity_only": True,
        "ignored_fields": result.get("ignored_fields") == ["Uz1", "K1", "P1", "d_x_Uz1"],
        "signed_factors": False,
        "z1_pieces": False,
        "both_chi_ids": False,
        "error_flags": True,
        "source_ids_retained": True,
        "analytic_fc_cell2": False,
        "quadrature_fxc_cell2": False,
        "analytic_z1_derivative": False,
    }

    for cell in (0, 1, 2):
        for branch in BRANCHES:
            mass = result["radial_phase_enclosures"][cell]["positive_kernel_masses"][branch]
            checks["positive_mass_bounds"] = checks["positive_mass_bounds"] and mass.a >= 0 and mass.b >= mass.a
    checks["zero_cell_velocity_only"] = all(
        result["group_results"][group]["paths"]["swirl"][order][3][branch] == []
        for group in GROUP_NAMES for order in (0, 1) for branch in BRANCHES
    )

    finite = result["group_results"]["Q_E01"]
    fc_z0 = finite["paths"]["swirl"][0][2]["F_from_Fc"]
    fxc_z0 = finite["paths"]["swirl"][0][2]["F_from_Fxc"]
    checks["signed_factors"] = bool(fc_z0) and bool(fxc_z0) and fc_z0[0]["coefficient"].a > 0 and fxc_z0[0]["coefficient"].b < 0
    z1_terms = finite["paths"]["swirl"][1][2]["F_from_Fc"]
    pieces = {
        tuple(entry["path"])[6]
        for term in z1_terms
        for entry in term["origin"].get("paths", [])
    }
    checks["z1_pieces"] = {
        "leading_Z*chi*inlet_Z0", "leading*chi_Z*inlet_Z0", "leading*chi*inlet_Z1",
    }.issubset(pieces)
    checks["both_chi_ids"] = all(
        sid in z1_terms[0]["origin"].get("chi_source_ids", [])
        for sid in ("fixture:chi:Z0", "fixture:chi:Z1")
    ) if z1_terms else False
    for group in ERROR_GROUPS:
        for order in (0, 1):
            for cell in (0, 1, 2):
                for branch in BRANCHES:
                    for term in result["group_results"][group]["paths"]["swirl"][order][cell][branch]:
                        checks["error_flags"] = checks["error_flags"] and term.get("error_only") is True and term.get("symmetric_error") is True and term["coefficient"].a == -term["coefficient"].b
    checks["source_ids_retained"] = all(
        term.get("source_ids") and term["origin"].get("source_ids") and term["origin"].get("chi_source_ids")
        for group in GROUP_NAMES for order in (0, 1) for cell in (0, 1, 2)
        for branch in BRANCHES for term in result["group_results"][group]["paths"]["swirl"][order][cell][branch]
    )

    import mpmath as mp

    mp.mp.dps = 80
    hb_mp = mp.mpf(1) / 3
    sc_mp = mp.mpf(1) / 2048
    lam_mp = mp.mpf(13) / 8
    bases = (hb_mp, mp.mpf(7) / 5, mp.mpf(11) / 6, lam_mp)
    fc_bounds = _sum_interval(mp, fc_z0, bases)
    fxc_bounds = _sum_interval(mp, fxc_z0, bases)
    fc_target = _target(
        mp, hb=hb_mp, sc=sc_mp, lam=lam_mp,
        q_lower=mp.mpf(1) / 2, q_upper=mp.mpf(3) / 4,
        leading=mp.mpf(3) / 7, chi=mp.mpf(2) / 5,
        inlet=mp.mpf(2) / 17, branch="F_from_Fc",
    )
    fxc_target = _target(
        mp, hb=hb_mp, sc=sc_mp, lam=lam_mp,
        q_lower=mp.mpf(1) / 2, q_upper=mp.mpf(3) / 4,
        leading=mp.mpf(3) / 7, chi=mp.mpf(2) / 5,
        inlet=mp.mpf(-2) / 7, branch="F_from_Fxc",
    )
    checks["analytic_fc_cell2"] = fc_bounds[0] <= fc_target <= fc_bounds[1]
    checks["quadrature_fxc_cell2"] = fxc_bounds[0] <= fxc_target <= fxc_bounds[1]

    z1_bounds = _sum_interval(mp, z1_terms, bases)
    z1_amplitude = (
        mp.mpf(-5) / 11 * mp.mpf(2) / 5 * mp.mpf(2) / 17
        + mp.mpf(3) / 7 * mp.mpf(-1) / 6 * mp.mpf(2) / 17
        + mp.mpf(3) / 7 * mp.mpf(2) / 5 * mp.mpf(-3) / 23
    )
    z1_target = _target(
        mp, hb=hb_mp, sc=sc_mp, lam=lam_mp,
        q_lower=mp.mpf(1) / 2, q_upper=mp.mpf(3) / 4,
        leading=z1_amplitude, chi=mp.mpf(1), inlet=mp.mpf(1),
        branch="F_from_Fc",
    )
    checks["analytic_z1_derivative"] = z1_bounds[0] <= z1_target <= z1_bounds[1]

    return {
        "passed": all(checks.values()),
        "checks": checks,
        "check_count": len(checks),
        "synthetic_only": True,
        "independent_analytic_and_quadrature": True,
        "source_receipts_read": False,
        "prior_helpers_imported": False,
        "formal_scales_materialized": False,
    }


def run(build: Any, c: Any) -> dict[str, Any]:
    return run_fixture(build, c)
