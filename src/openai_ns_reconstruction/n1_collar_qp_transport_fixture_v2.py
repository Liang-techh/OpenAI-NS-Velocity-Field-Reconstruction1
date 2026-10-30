"""Independent synthetic checks for :mod:`n1_collar_qp_transport_v2`.

The fixture owns all of its rational data.  It deliberately does not import a
producer, a previous transport, a checker, or a receipt.  The high-precision
quadrature is only a local cross-check of the first potential word; it is not
used as a source certificate.  Importing this module does not execute the
fixture.
"""

from __future__ import annotations

from collections.abc import Mapping
from fractions import Fraction
from typing import Any


SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
FIELDS = ("F1", "Uz1", "K1", "P1", "d_x_F1", "d_x_Uz1")
BRANCHES = ("F_from_Fc", "F_from_Fxc", "Fx_from_Fc", "Fx_from_Fxc")
FIRST_HB_POWER_ADDED = {
    "F_from_Fc": 2,
    "F_from_Fxc": 3,
    "Fx_from_Fc": 1,
    "Fx_from_Fxc": 2,
}
J2_HB_POWER_ADDED = {
    "F_from_Fc": 4,
    "F_from_Fxc": 5,
    "Fx_from_Fc": 3,
    "Fx_from_Fxc": 4,
}


def _q(c: Any, value: Fraction | int) -> Any:
    if isinstance(value, Fraction):
        return c.mpf(value.numerator) / value.denominator
    return c.mpf(value)


def _box(c: Any, lower: Fraction | int, upper: Fraction | int) -> Any:
    lo, hi = _q(c, lower), _q(c, upper)
    return c.make_mpf((lo._mpi_[0], hi._mpi_[1]))


def _term(
    c: Any,
    value: Fraction | int | Any,
    source_id: str,
    powers: tuple[int, int, int, int],
    *,
    origin: Any,
    error_only: bool = False,
    error_class: str = "v2_fixture_qp_signed",
) -> dict[str, Any]:
    result = {
        "powers": list(powers),
        "coefficient": value if hasattr(value, "_mpi_") else _q(c, value),
        "source_ids": [source_id],
        "error_classes": [error_class],
        "uncertainty": error_only,
        "origin": origin,
    }
    if error_only:
        result.update({"error_only": True, "symmetric_error": True})
    return result


def _component_packet(
    c: Any,
    prefix: str,
    values: Mapping[str, tuple[Fraction, Fraction]],
    powers: Mapping[str, tuple[int, int, int, int]],
    *,
    error_only: bool = False,
    structural_zero: tuple[str, ...] = (),
) -> dict[str, Any]:
    packet: dict[str, Any] = {}
    for field in FIELDS:
        if field in structural_zero:
            packet[field] = {"Z0": [], "Z1": []}
            continue
        canonical = {
            "F1": "F",
            "Uz1": "U",
            "K1": "K",
            "P1": "P",
            "d_x_F1": "Fx",
            "d_x_Uz1": "Ux",
        }[field]
        z0, z1 = values[canonical]
        source_power = powers[canonical]
        z0_id = f"{prefix}:{field}:Z0"
        z1_id = f"{prefix}:{field}:Z1"
        z0_origin: Any = f"{prefix}:text:{field}:Z0"
        z1_origin: Any = {"tag": f"{prefix}:map:{field}:Z1", "source_ids": [z1_id]}
        packet[field] = {
            "Z0": [_term(c, z0, z0_id, source_power, origin=z0_origin, error_only=error_only, error_class=prefix + (":error" if error_only else ":signed"))],
            "Z1": [_term(c, z1, z1_id, source_power, origin=z1_origin, error_only=error_only, error_class=prefix + (":error" if error_only else ":signed"))],
        }
    return packet


def _radius(c: Any, value: Any) -> Any:
    lo, hi = value._mpi_
    lo, hi = c.make_mpf((lo, lo)), c.make_mpf((hi, hi))
    return max(-lo, hi)


def _sector(mp: Any, terms: Any, bases: tuple[Any, Any, Any, Any]) -> tuple[Any, Any]:
    lower, upper = mp.mpf(0), mp.mpf(0)
    for term in terms:
        coeff = term["coefficient"]
        lo_bits, hi_bits = coeff._mpi_
        def endpoint(bits: Any) -> Any:
            sign, mantissa, exponent, _ = bits
            value = mp.mpf(int(mantissa)) * (mp.mpf(2) ** int(exponent))
            return -value if sign else value
        lo, hi = endpoint(lo_bits), endpoint(hi_bits)
        scale = mp.mpf(1)
        for base, power in zip(bases, term["powers"]):
            scale *= base ** int(power)
        if scale >= 0:
            lower, upper = lower + lo * scale, upper + hi * scale
        else:
            lower, upper = lower + hi * scale, upper + lo * scale
    return lower, upper


def _first_integral(mp: Any, branch: str, hb: Any, sc: Any, chi: Any, inlet: Any) -> Any:
    """Independent quadrature for one first-potential endpoint branch."""

    X = 2 * mp.exp(hb * sc / 8)
    kernel = (lambda t: t * (1 - (t / X) ** 2) / 2) if branch.startswith("F_") else (lambda t: (t / X) ** 3)
    basis = (lambda t: mp.mpf(1)) if branch.endswith("Fc") else (lambda t: 1 - 4 / t**2)
    return -2 * chi * inlet * mp.quad(lambda t: kernel(t) * basis(t), [2, X])


def _formal(mp: Any, coefficient: Any, powers: tuple[int, int, int, int], bases: tuple[Any, Any, Any, Any]) -> Any:
    value = mp.mpf(coefficient)
    for base, power in zip(bases, powers):
        value *= base ** int(power)
    return value


def run_fixture(build: Any, c: Any) -> dict[str, Any]:
    """Run independent structural, enclosure, and quadrature checks."""

    if hasattr(build, "build"):
        build = build.build
    if not callable(build):
        raise TypeError("build must be callable or expose build")

    hb = _box(c, 0, Fraction(1, 2))
    sc = _q(c, Fraction(3, 91))
    signed_powers = {
        "F": (-5, 2, 0, 0),
        "U": (2, -2, 0, 0),
        "K": (0, 0, -2, 0),
        "P": (3, 0, 2, -2),
        "Fx": (-4, 0, 0, 2),
        "Ux": (2, 0, -2, 0),
    }
    signed_values = {
        "F": (Fraction(7, 23), Fraction(-11, 31)),
        "U": (Fraction(-9, 37), Fraction(8, 41)),
        "K": (Fraction(5, 43), Fraction(-11, 67)),
        "P": (Fraction(9, 59), Fraction(13, 61)),
        "Fx": (Fraction(-11, 67), Fraction(17, 71)),
        "Ux": (Fraction(7, 73), Fraction(-5, 79)),
    }
    error_powers = {
        "F": (-6, 0, 0, 0),
        "U": (0, -2, 0, 0),
        "K": (0, 0, -2, 0),
        "P": (2, 0, 2, -2),
        "Fx": (-5, 0, 0, 2),
        "Ux": (0, 0, -2, 0),
    }
    error_values = {
        "F": (_box(c, Fraction(-1, 811), Fraction(1, 811)), _box(c, Fraction(-1, 787), Fraction(1, 787))),
        "U": (_box(c, Fraction(-1, 769), Fraction(1, 769)), _box(c, Fraction(-1, 751), Fraction(1, 751))),
        "K": (_box(c, Fraction(-1, 829), Fraction(1, 829)), _box(c, Fraction(-1, 827), Fraction(1, 827))),
        "P": (_box(c, Fraction(-1, 839), Fraction(1, 839)), _box(c, Fraction(-1, 853), Fraction(1, 853))),
        "Fx": (_box(c, Fraction(-1, 877), Fraction(1, 877)), _box(c, Fraction(-1, 881), Fraction(1, 881))),
        "Ux": (_box(c, Fraction(-1, 613), Fraction(1, 613)), _box(c, Fraction(-1, 601), Fraction(1, 601))),
    }
    e10_powers = {
        key: tuple(power - (1 if key in {"F", "Fx"} else 0) for power in value)
        for key, value in signed_powers.items()
    }
    groups = {
        "Q01_v2": _component_packet(c, "v2fixture:Q01_v2", signed_values, signed_powers),
        "Q10_v2": _component_packet(c, "v2fixture:Q10_v2", signed_values, e10_powers),
        "Q11_v2": _component_packet(c, "v2fixture:Q11_v2", signed_values, signed_powers, structural_zero=("Uz1", "K1", "P1", "d_x_Uz1")),
        "P_v2": _component_packet(c, "v2fixture:P_v2", signed_values, signed_powers),
        "error_only_v2": _component_packet(c, "v2fixture:error_only_v2", error_values, error_powers, error_only=True),
    }
    # Keep one row genuinely multi-sector: the transport must expand this row
    # linearly, without adding a synthetic sector or cross-multiplying sources.
    groups["Q01_v2"]["F1"]["Z0"].append(
        _term(c, Fraction(5, 113), "v2fixture:Q01_v2:F1:Z0:extra", signed_powers["F"], origin="v2fixture:Q01_v2:extra")
    )
    # The synthetic Q11_v2 packet intentionally exercises explicit structural
    # zero rows.  legacy_A00_v2 is omitted: a caller that already entered it must not
    # silently cause this adapter to recompute its K0 branch.
    chi = {
        "Z0": [_term(c, Fraction(5, 19), "v2fixture:chi:Z0", (0, 0, 0, 0), origin="v2fixture:chi:Z0", error_class="v2_fixture_chi")],
        "Z1": [_term(c, Fraction(-7, 23), "v2fixture:chi:Z1", (0, 0, 0, 0), origin={"tag": "v2fixture:chi:Z1"}, error_class="v2_fixture_chi")],
    }
    certificate = {
        "hb_sc_Z_independent": True,
        "rho_radius_Z_independent": True,
        "strict_hb_admitted": True,
        "fixture_geometry_owner": "independent_qp_transport_v2",
    }
    atoms = {name: "v2fixture:" + name for name in SOURCE_POWER_BASIS}
    result = build(
        c,
        incoming={"groups": groups},
        chi=chi,
        hb=hb,
        sc=sc,
        source_family="synthetic_qp_transport_v2",
        source_atoms=atoms,
        geometry_certificate=certificate,
    )

    checks: dict[str, bool] = {
        "schema": result.get("schema") == "n1_collar_qp_transport_v2",
        "formal_basis": tuple(result.get("source_power_basis", ())) == SOURCE_POWER_BASIS,
        "fields": tuple(result.get("fields", ())) == ("F", "U", "K", "P", "Fx", "Ux"),
        "groups_preserved": tuple(result.get("inlet_groups", ())) == tuple(groups),
        "groups_not_merged": result.get("grouped_inlet_packets_merged") is False and set(result.get("group_results", {})) == set(groups),
        "closed_hb_certified": result.get("geometry", {}).get("closed_hb_enclosure_allowed") is True,
        "formal_not_materialized": result.get("formal_source_powers_not_materialized") is True and result.get("geometry", {}).get("formal_hb_not_materialized") is True,
        "fixed_geometry": result.get("geometry", {}).get("ordinary_Z_geometry_fixed") is True and result.get("geometry", {}).get("ordinary_Z_geometry_differentiated") is False,
        "stages_disjoint": result.get("K0_j1_j2_disjoint") is True and result.get("baseline_K0_reentered") is False,
        "no_QP_claim": result.get("complete_QP_feedback") is False and result.get("m15_computed") is False,
        "multi_sector_row_preserved": len(result["group_results"]["Q01_v2"]["first_potential"]["branches"]["F_from_Fc"]["Z0"]) == 2,
    }

    # K0 retains all six state routes for the ordinary groups, while explicit
    # structural zero rows stay empty rather than becoming fake zero sectors.
    for group in groups:
        stages = result["group_results"][group]
        source_power_table = (
            error_powers
            if group == "error_only_v2"
            else e10_powers
            if group == "Q10_v2"
            else signed_powers
        )
        k0_fields = stages["k0"]["fields"]
        first_fields = stages["first_potential"]["fields"]
        tail_fields = stages["j2_tail"]["fields"]
        checks[group + ":K0_keys"] = set(k0_fields) == {"F", "U", "K", "P", "Fx", "Ux"}
        checks[group + ":structural_zero_lists"] = all(
            first_fields[field] == {"Z0": [], "Z1": []} and tail_fields[field] == {"Z0": [], "Z1": []}
            for field in ("U", "K", "P", "Ux")
        )
        checks[group + ":first_four_branches"] = tuple(stages["first_potential"]["branches"]) == BRANCHES
        checks[group + ":tail_four_branches"] = tuple(stages["j2_tail"]["branches"]) == BRANCHES
        checks[group + ":Z1_product_rule"] = stages["first_potential"]["ordinary_Z1_product_rule"] == "chi_Z*inlet_Z0 + chi*inlet_Z1"
        if group == "Q11_v2":
            checks[group + ":input_structural_zero"] = all(
                not stages["k0"]["fields"][field]["Z0"] and not stages["k0"]["fields"][field]["Z1"]
                for field in ("U", "K", "P", "Ux")
            )
        else:
            checks[group + ":K0_nonempty"] = any(stages["k0"]["fields"][field]["Z0"] for field in ("F", "U", "K", "P", "Fx", "Ux"))
        for branch in BRANCHES:
            first_terms = stages["first_potential"]["branches"][branch]["Z0"] + stages["first_potential"]["branches"][branch]["Z1"]
            tail_terms = stages["j2_tail"]["branches"][branch]["Z0"] + stages["j2_tail"]["branches"][branch]["Z1"]
            checks[f"{group}:{branch}:first_terms"] = bool(first_terms) if group != "Q11_v2" or branch.endswith("Fc") else True
            checks[f"{group}:{branch}:tail_terms"] = bool(tail_terms) if group != "Q11_v2" or branch.endswith("Fc") else True
            checks[f"{group}:{branch}:first_power"] = all(
                term["powers"][0] - source_power_table["F" if branch.endswith("Fc") else "Fx"][0] == FIRST_HB_POWER_ADDED[branch]
                for term in first_terms
            )
            checks[f"{group}:{branch}:tail_power"] = all(
                term["powers"][0] - source_power_table["F" if branch.endswith("Fc") else "Fx"][0] == J2_HB_POWER_ADDED[branch]
                for term in tail_terms
            )
            checks[f"{group}:{branch}:provenance"] = all(term.get("source_ids") and term.get("incoming_source_ids") and term.get("chi_source_ids") for term in first_terms + tail_terms)
            checks[f"{group}:{branch}:tail_symmetric"] = all(
                term.get("error_only") is True
                and term.get("symmetric_error") is True
                and term["coefficient"].a <= 0 <= term["coefficient"].b
                and (-term["coefficient"].a) == term["coefficient"].b
                for term in tail_terms
            )

    checks["negative_formal_powers_preserved"] = any(
        term["powers"][0] < 0
        for field in ("F", "U", "K", "P", "Fx", "Ux")
        for order in ("Z0", "Z1")
        for term in result["group_results"]["Q10_v2"]["inlet"][field][order]
    )
    checks["error_group_separate"] = result["group_results"]["error_only_v2"]["input_error_only"] is True and all(
        term.get("error_only") is True and term.get("symmetric_error") is True
        for branch in BRANCHES
        for stage in ("first_potential", "j2_tail")
        for order in ("Z0", "Z1")
        for term in result["group_results"]["error_only_v2"][stage]["branches"][branch][order]
    )

    # Independent high-precision first-word containment.  The returned mass is
    # evaluated against a concrete hb and source basis only here.
    import mpmath as mp

    mp.mp.dps = 83
    hb_mp, sc_mp = mp.mpf(3) / 10, mp.mpf(3) / 91
    bases = (hb_mp, mp.mpf(9) / 8, mp.mpf(17) / 5, mp.mpf(19) / 7)
    chi_mp = mp.mpf(5) / 19
    incoming_mp = _formal(mp, mp.mpf(7) / 23, signed_powers["F"], bases)
    incoming_mp += _formal(mp, mp.mpf(5) / 113, signed_powers["F"], bases)
    target = _first_integral(mp, "F_from_Fc", hb_mp, sc_mp, chi_mp, incoming_mp)
    terms = result["group_results"]["Q01_v2"]["first_potential"]["branches"]["F_from_Fc"]["Z0"]
    got = _sector(mp, terms, bases)
    checks["analytic_first_integral_enclosed"] = got[0] <= target <= got[1]
    checks["analytic_quadrature_independent"] = True

    # The strict-certificate guard is active even when the numerical hb box is
    # closed: removing it must reject the same input.
    weak = dict(certificate)
    weak.pop("strict_hb_admitted", None)
    try:
        build(c, incoming={"groups": groups}, chi=chi, hb=hb, sc=sc, source_family="synthetic_qp_transport_v2", source_atoms=atoms, geometry_certificate=weak)
    except (TypeError, ValueError):
        checks["missing_strict_certificate_rejected"] = True
    else:
        checks["missing_strict_certificate_rejected"] = False

    # Direct aliases are intentionally exclusive.  Test every six-state pair
    # so a future API change cannot silently prefer one spelling and discard
    # the other.
    duplicate_alias_cases = {
        "F/Fc": ("F", "Fc", "F1"),
        "U/Uc": ("U", "Uc", "Uz1"),
        "K/Kc": ("K", "Kc", "K1"),
        "P/Pc": ("P", "Pc", "P1"),
        "Fx/Fxc": ("Fx", "Fxc", "d_x_F1"),
        "Ux/Uxc": ("Ux", "Uxc", "d_x_Uz1"),
    }
    for label, (canonical, alias, field) in duplicate_alias_cases.items():
        try:
            build(
                c,
                chi=chi,
                hb=hb,
                sc=sc,
                geometry_certificate=certificate,
                **{canonical: groups["Q01_v2"][field], alias: groups["Q01_v2"][field]},
            )
        except (TypeError, ValueError) as exc:
            checks["duplicate_direct_" + label.replace("/", "_") + "_rejected"] = "supplied twice" in str(exc)
        else:
            checks["duplicate_direct_" + label.replace("/", "_") + "_rejected"] = False

    try:
        build(
            c,
            incoming={"groups": groups},
            chi=chi,
            hb=hb,
            sc=sc,
            geometry_certificate=certificate,
            source_atoms=atoms,
            F=groups["Q01_v2"]["F1"],
        )
    except (TypeError, ValueError) as exc:
        checks["direct_grouped_input_mix_rejected"] = "cannot be mixed" in str(exc)
    else:
        checks["direct_grouped_input_mix_rejected"] = False

    return {
        "passed": all(checks.values()),
        "check_count": len(checks),
        "checks": checks,
        "synthetic_only": True,
        "independent_high_precision_quadrature": True,
        "negative_formal_powers_checked": True,
        "input_error_separation_checked": True,
        "structural_zero_fields_checked": True,
        "receipt_reads": False,
        "old_builder_called": False,
        "automatic_code_execution": False,
    }


run = run_fixture

__all__ = [
    "SOURCE_POWER_BASIS",
    "FIELDS",
    "BRANCHES",
    "run_fixture",
    "run",
]


