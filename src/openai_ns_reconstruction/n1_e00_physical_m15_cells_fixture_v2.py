"""Independent synthetic checks for :mod:`n1_e00_physical_m15_cells_v2`.

The fixture creates fresh signed, nonconstant interval sectors for all three
active q-cell enclosures.  It does not import a project module, read a receipt,
or run on import.  ``run_fixture`` receives the leaf explicitly and checks the
physical signs, ordinary-Z1 product rule, formal powers, symmetric J2 tail,
zero-mass tail, source identities, and one independent constant-cell integral.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Any

SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
PARTS = ("K0", "J1", "J2")
ACTIVE_CELLS = (0, 1, 2)


def _q(c: Any, value: Fraction | int) -> Any:
    if isinstance(value, Fraction):
        return c.mpf(value.numerator) / c.mpf(value.denominator)
    return c.mpf(value)


def _box(c: Any, lower: Fraction | int, upper: Fraction | int) -> Any:
    lo, hi = _q(c, lower), _q(c, upper)
    return c.make_mpf((lo._mpi_[0], hi._mpi_[1]))


def _term(
    c: Any,
    lower: Fraction | int,
    upper: Fraction | int,
    source_id: str,
    *,
    powers: tuple[int, int, int, int] = (0, 0, 0, 0),
    path: tuple[Any, ...],
    error_only: bool = False,
) -> dict[str, Any]:
    classes = ["fixture_signed_source"]
    if error_only:
        classes.append("fixture_symmetric_tail")
    origin = {
        "source_ids": [source_id],
        "source_paths": [list(path)],
        "error_classes": classes,
        "uncertainty": bool(error_only),
        "error_only": bool(error_only),
        "symmetric_error": bool(error_only),
    }
    return {
        "powers": list(powers),
        "coefficient": _box(c, lower, upper),
        "source_ids": [source_id],
        "source_paths": [list(path)],
        "source_path": list(path),
        "error_classes": classes,
        "uncertainty": bool(error_only),
        "error_only": bool(error_only),
        "symmetric_error": bool(error_only),
        "origin": origin,
    }


def _cell_packet(c: Any, prefix: str, *, tail: bool = False) -> dict[int, dict[int, list[dict[str, Any]]]]:
    result: dict[int, dict[int, list[dict[str, Any]]]] = {}
    for cell in ACTIVE_CELLS:
        if tail:
            width0 = Fraction(1, 1300 + cell)
            width1 = Fraction(1, 1500 + cell)
            z0 = _term(
                c,
                -width0,
                width0,
                f"fixture:v2:{prefix}:cell{cell}:Z0",
                path=(prefix, "cell", cell, "Z0"),
                error_only=True,
            )
            z1 = _term(
                c,
                -width1,
                width1,
                f"fixture:v2:{prefix}:cell{cell}:Z1",
                path=(prefix, "cell", cell, "Z1"),
                error_only=True,
            )
        else:
            base = Fraction(5 + cell, 2)
            z0 = _term(
                c,
                base,
                base + Fraction(1, 100),
                f"fixture:v2:{prefix}:cell{cell}:Z0",
                powers=(0, 0, 0, 0),
                path=(prefix, "cell", cell, "Z0"),
            )
            z1 = _term(
                c,
                -Fraction(2, 11 + cell),
                -Fraction(2, 13 + cell),
                f"fixture:v2:{prefix}:cell{cell}:Z1",
                powers=(0, 1, 0, 0),
                path=(prefix, "cell", cell, "Z1"),
            )
        result[cell] = {0: [z0], 1: [z1]}
    return result


def _leading_packet(c: Any, prefix: str, *, swirl: bool) -> dict[int, dict[int, list[dict[str, Any]]]]:
    result: dict[int, dict[int, list[dict[str, Any]]]] = {}
    for cell in ACTIVE_CELLS:
        if swirl:
            z0_lo, z0_hi = Fraction(3 + cell, 1), Fraction(3 + cell, 1) + Fraction(1, 100)
            z1_lo, z1_hi = Fraction(-2, 17 + cell), Fraction(-2, 19 + cell)
            powers0, powers1 = (0, 0, 0, 0), (0, 0, 1, 0)
        else:
            z0_lo, z0_hi = Fraction(3 + cell, 1), Fraction(3 + cell, 1) + Fraction(1, 100)
            z1_lo, z1_hi = Fraction(2, 11 + cell), Fraction(2, 9 + cell)
            powers0, powers1 = (0, 0, 0, 0), (0, 0, 0, 0)
        result[cell] = {
            0: [
                _term(c, z0_lo, z0_hi, f"fixture:v2:{prefix}:cell{cell}:Z0", powers=powers0, path=(prefix, "cell", cell, "Z0"))
            ],
            1: [
                _term(c, z1_lo, z1_hi, f"fixture:v2:{prefix}:cell{cell}:Z1", powers=powers1, path=(prefix, "cell", cell, "Z1"))
            ],
        }
    return result


def _inputs(c: Any) -> tuple[Any, Any, Any]:
    leading_u = _leading_packet(c, "leading_U", swirl=False)
    leading_f = _leading_packet(c, "leading_F", swirl=True)
    fields = {
        "K0": {
            "U": _cell_packet(c, "K0_U"),
            "F": _cell_packet(c, "K0_F"),
        },
        "J1": {"F": _cell_packet(c, "J1_F")},
        "J2": {"F": _cell_packet(c, "J2_F", tail=True)},
    }
    return leading_u, leading_f, fields


def _endpoint(value: Any, upper: bool) -> Any:
    raw = value._mpi_
    return raw[1 if upper else 0]


def _formal_interval(mp: Any, term: Any, bases: tuple[Any, Any, Any, Any]) -> tuple[Any, Any]:
    lo = mp.make_mpf(_endpoint(term["coefficient"], False))
    hi = mp.make_mpf(_endpoint(term["coefficient"], True))
    lower, upper = lo, hi
    for base, power in zip(bases, term["powers"]):
        scale = base**power
        if scale >= 0:
            lower *= scale
            upper *= scale
        else:
            lower, upper = upper * scale, lower * scale
    return lower, upper


def _sum_interval(mp: Any, terms: Any, bases: tuple[Any, Any, Any, Any]) -> tuple[Any, Any]:
    lower = mp.mpf(0)
    upper = mp.mpf(0)
    for term in terms:
        lo, hi = _formal_interval(mp, term, bases)
        lower += lo
        upper += hi
    return lower, upper


def _all_terms(result: Any, part: str, channel: str, order: int, cell: int) -> list[Any]:
    channels = result["parts"][part]["channels"]
    if channel not in channels:
        return []
    return channels[channel][order][cell]


def run_fixture(leaf: Any, c: Any) -> dict[str, Any]:
    """Run independent signed, scale, Z1, and error-tail checks."""

    builder = getattr(leaf, "build", leaf)
    if not callable(builder):
        raise TypeError("leaf must be callable or expose build")
    leading_u, leading_f, fields = _inputs(c)
    certificate = {
        "hb_strictly_between_zero_and_half": True,
        "sc_phase_Z_independent": True,
        "hb_epsilon_Z_independent": True,
        "original_cutoff_mass_5_over_8": True,
        "original_cutoff_mass_by_cell": True,
    }
    atoms = {name: f"fixture:v2:{name}" for name in SOURCE_POWER_BASIS}
    sc = c.mpf(1) / c.mpf(4096)
    result = builder(
        c,
        sc=sc,
        leading_U=leading_u,
        leading_F=leading_f,
        field_sectors=fields,
        source_family={"fixture": "e00-physical-m15", "generation": 1},
        source_atoms=atoms,
        certificate=certificate,
    )
    checks: dict[str, bool] = {}
    checks["schema"] = result.get("schema") == "n1_e00_physical_m15_cells_v2"
    checks["parts_separated"] = set(result.get("parts", {})) == {"K0", "J1", "J2"}
    checks["three_q_cells"] = result.get("active_cells") == [0, 1, 2]
    checks["exact_cutoff_masses"] = [
        (row["eta_mass"].a, row["eta_mass"].b) for row in result.get("q_cells", [])
    ] == [
        (c.mpf(1) / c.mpf(4), c.mpf(1) / c.mpf(4)),
        (c.mpf(1) / c.mpf(4), c.mpf(1) / c.mpf(4)),
        (c.mpf(1) / c.mpf(8), c.mpf(1) / c.mpf(8)),
    ] and result.get("zero_mass_tail", {}).get("zero_mass") is True
    checks["no_ledger"] = all(result.get(key) is False for key in ("ledger_appended", "appended_to_ledger", "ledger_written"))
    checks["no_operator_reapplication"] = result.get("j1_j2_reapplied") is False and result.get("chi_reapplied") is False and result.get("cutoff_reapplied") is False
    checks["no_P1_m15"] = result.get("P1_is_m15") is False
    checks["no_prefix_claim"] = result.get("prefix_mutual_exclusivity_claimed") is False

    # Every active product has exactly one Z0 term and both Leibniz terms at
    # Z1.  The J1/J2 parts expose no axial channel.
    for part in PARTS:
        channels = ("axial", "swirl") if part == "K0" else ("swirl",)
        for channel in channels:
            for cell in ACTIVE_CELLS:
                z0 = _all_terms(result, part, channel, 0, cell)
                z1 = _all_terms(result, part, channel, 1, cell)
                checks[f"{part}_{channel}_z0_count_cell{cell}"] = len(z0) == 1
                checks[f"{part}_{channel}_z1_count_cell{cell}"] = len(z1) == 2
                checks[f"{part}_{channel}_z1_pieces_cell{cell}"] = {
                    term["ordinary_Z_piece"] for term in z1
                } == {"lead1*field0", "lead0*field1"}
                checks[f"{part}_{channel}_ids_cell{cell}"] = all(
                    term.get("source_ids")
                    and term["origin"].get("leading_source_ids")
                    and term["origin"].get("field_source_ids")
                    and term["origin"].get("source_paths")
                    for term in z0 + z1
                )
                checks[f"{part}_{channel}_powers_cell{cell}"] = all(
                    term["powers"][0] == 1
                    and term["powers"][3] == (-1 if channel == "axial" else -2)
                    for term in z0 + z1
                )
        checks[f"{part}_no_axial_for_non_k0"] = part == "K0" or all(
            not _all_terms(result, part, "axial", order, cell)
            for order in (0, 1)
            for cell in range(4)
        )

    k0_axial = _all_terms(result, "K0", "axial", 0, 0)
    k0_swirl = _all_terms(result, "K0", "swirl", 0, 0)
    j1_swirl = _all_terms(result, "J1", "swirl", 0, 0)
    checks["K0_axial_uses_Ulead"] = all("fixture:v2:leading_U:cell0:Z0" in term["leading_source_ids"] for term in k0_axial)
    checks["K0_swirl_uses_Flead"] = all("fixture:v2:leading_F:cell0:Z0" in term["leading_source_ids"] for term in k0_swirl)
    checks["K0_swirl_Z1_uses_Flead_Z"] = any("fixture:v2:leading_F:cell0:Z1" in term["leading_source_ids"] for term in _all_terms(result, "K0", "swirl", 1, 0))
    checks["positive_axial_sign"] = bool(k0_axial) and k0_axial[0]["coefficient"].a > 0
    checks["negative_swirl_sign"] = bool(k0_swirl) and k0_swirl[0]["coefficient"].b < 0
    checks["negative_j1_swirl_sign"] = bool(j1_swirl) and j1_swirl[0]["coefficient"].b < 0

    j2_terms = [
        term
        for cell in ACTIVE_CELLS
        for order in (0, 1)
        for term in _all_terms(result, "J2", "swirl", order, cell)
    ]
    checks["j2_error_only"] = bool(j2_terms) and all(term.get("error_only") is True for term in j2_terms)
    checks["j2_symmetric"] = bool(j2_terms) and all(
        term.get("symmetric_error") is True and term["coefficient"].a == -term["coefficient"].b
        for term in j2_terms
    )
    checks["zero_mass_tail"] = all(
        not _all_terms(result, part, channel, order, 3)
        for part in PARTS
        for channel in ("axial", "swirl")
        for order in (0, 1)
    )

    # Independent nondegenerate constant-cell integrals.  The source
    # intervals enclose the chosen constants and the returned rho/rho^2
    # factors use only the directed hb range, so both exact references must
    # lie in the output intervals after formal-power evaluation.
    import mpmath as mp

    mp.mp.dps = 70
    hb = mp.mpf(2) / 5
    lam = mp.mpf(3)
    sc_mp = mp.mpf(1) / 4096
    rho = lambda q: 4 * mp.exp(hb * sc_mp * q / 4)
    axial_ref = 2 * lam**-1 * mp.mpf("3.005") * mp.mpf("2.505") * (rho(mp.mpf(1) / 4) - rho(0))
    swirl_ref = -2 * lam**-2 * mp.mpf("3.005") * mp.mpf("2.505") * (hb * sc_mp / 4) * mp.quad(
        lambda q: rho(q) ** 2, [0, mp.mpf(1) / 4]
    )
    bases = (hb, mp.mpf(1), mp.mpf(1), lam)
    axial_bounds = _sum_interval(mp, k0_axial, bases)
    swirl_bounds = _sum_interval(mp, k0_swirl, bases)
    checks["analytic_axial_cell0"] = axial_bounds[0] <= axial_ref <= axial_bounds[1]
    checks["analytic_swirl_cell0"] = swirl_bounds[0] <= swirl_ref <= swirl_bounds[1]
    checks["physical_scale_distinct"] = all(
        term["powers"][3] == -1 for term in k0_axial
    ) and all(term["powers"][3] == -2 for term in k0_swirl)

    return {
        "passed": all(checks.values()),
        "checks": checks,
        "check_count": len(checks),
        "synthetic_only": True,
        "independent_constant_cell_integrals": True,
        "source_receipts_read": False,
        "prior_helpers_imported": False,
        "formal_scales_materialized": False,
        "ledger_touched": False,
    }


def run(leaf: Any, c: Any) -> dict[str, Any]:
    return run_fixture(leaf, c)


__all__ = ["SOURCE_POWER_BASIS", "PARTS", "run_fixture", "run"]
