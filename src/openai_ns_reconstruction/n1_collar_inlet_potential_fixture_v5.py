"""Independent v5 fixture for the finite first ``Kchi`` inlet leaf.

The fixture owns its rational data and does not import a previous leaf or
fixture.  It exercises string and mapping ``origin`` metadata, signed and
symmetric error-only terms, the three caller-owned inlet groups, a certified
closed ``hb`` enclosure, all four relative formal ``hb`` increments, and
independent high-precision physical integral enclosures.  The integral
comparison is performed only by :func:`run_fixture`; importing this module
does not execute it.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Any

SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
BRANCHES = ("F_from_Fc", "F_from_Fxc", "Fx_from_Fc", "Fx_from_Fxc")
SOURCE_POWERS = {"Fc": (4, -3, 1, 0), "Fxc": (-4, 2, -1, 2)}
HB_ADDED = {"F_from_Fc": 2, "F_from_Fxc": 3, "Fx_from_Fc": 1, "Fx_from_Fxc": 2}


def _q(c: Any, value: Fraction | int) -> Any:
    return c.mpf(value.numerator) / value.denominator if isinstance(value, Fraction) else c.mpf(value)


def _box(c: Any, left: Fraction | int, right: Fraction | int) -> Any:
    lo, hi = _q(c, left), _q(c, right)
    return c.make_mpf((lo._mpi_[0], hi._mpi_[1]))


def _term(
    c: Any,
    value: Any,
    source_id: str,
    powers: tuple[int, int, int, int],
    *,
    origin: Any,
    error_only: bool = False,
    error_class: str = "v5_signed",
) -> dict[str, Any]:
    result = {
        "powers": list(powers),
        "coefficient": value if hasattr(value, "_mpi_") else c.mpf(value),
        "source_ids": [source_id],
        "error_classes": [error_class],
        "uncertainty": error_only,
        "origin": origin,
    }
    if error_only:
        result.update({"error_only": True, "symmetric_error": True})
    return result


def _packet(c: Any, component: str, z0: Fraction, z1: Fraction, *, error: Any = None) -> dict[str, Any]:
    powers = SOURCE_POWERS[component]
    z0_id = f"v5:{component}:Z0"
    z1_id = f"v5:{component}:Z1"
    if component == "Fc":
        z0_origin: Any = "v5-origin-text:Fc:Z0"
        z1_origin: Any = {"tag": "v5-origin-map:Fc:Z1", "source_ids": [z1_id], "error_classes": ["v5_fc_map"]}
    else:
        z0_origin = {"tag": "v5-origin-map:Fxc:Z0", "source_ids": [z0_id], "error_classes": ["v5_fxc_map"]}
        z1_origin = "v5-origin-text:Fxc:Z1"
    z0_terms = [_term(c, _q(c, z0), z0_id, powers, origin=z0_origin)]
    if error is not None:
        z0_terms.append(
            _term(
                c,
                error,
                f"v5:{component}:Z0:error",
                powers,
                origin=f"v5-origin-text:{component}:Z0:error",
                error_only=True,
                error_class="v5_symmetric_tail",
            )
        )
    return {"Z0": z0_terms, "Z1": [_term(c, _q(c, z1), z1_id, powers, origin=z1_origin)]}


def _mp_endpoint(bits: Any, mp: Any) -> Any:
    sign, mantissa, exponent, _ = bits
    value = mp.mpf(int(mantissa)) * (mp.mpf(2) ** int(exponent))
    return -value if sign else value


def _sector(mp: Any, terms: Any, bases: tuple[Any, Any, Any, Any]) -> tuple[Any, Any]:
    lower, upper = mp.mpf(0), mp.mpf(0)
    for term in terms:
        coeff = term["coefficient"]
        lo, hi = _mp_endpoint(coeff._mpi_[0], mp), _mp_endpoint(coeff._mpi_[1], mp)
        scale = mp.mpf(1)
        for base, power in zip(bases, term["powers"]):
            scale *= base ** int(power)
        if scale >= 0:
            lower, upper = lower + lo * scale, upper + hi * scale
        else:
            lower, upper = lower + hi * scale, upper + lo * scale
    return lower, upper


def _formal(mp: Any, coefficient: Any, powers: tuple[int, int, int, int], bases: tuple[Any, Any, Any, Any]) -> Any:
    value = mp.mpf(coefficient)
    for base, power in zip(bases, powers):
        value *= base ** int(power)
    return value


def _integral(mp: Any, branch: str, hb: Any, sc: Any, chi: Any, incoming: Any) -> Any:
    X = 2 * mp.exp(hb * sc / 8)
    kernel = (lambda t: t * (1 - (t / X) ** 2) / 2) if branch.startswith("F_") else (lambda t: (t / X) ** 3)
    basis = (lambda t: mp.mpf(1)) if branch.endswith("Fc") else (lambda t: 1 - 4 / t**2)
    return -2 * chi * incoming * mp.quad(lambda t: kernel(t) * basis(t), [2, X])


def run_fixture(build: Any, c: Any) -> dict[str, Any]:
    """Run independent v5 checks against a supplied build callable."""

    if hasattr(build, "build"):
        build = build.build
    if not callable(build):
        raise TypeError("build must be callable or expose build")

    hb = _box(c, Fraction(1, 9), Fraction(7, 16))
    sc = _q(c, Fraction(1, 96))
    source = {
        "Fc": _packet(c, "Fc", Fraction(5, 19), Fraction(-7, 26)),
        "Fxc": _packet(c, "Fxc", Fraction(-11, 32), Fraction(9, 35), error=_box(c, Fraction(-1, 733), Fraction(1, 733))),
    }
    chi = {
        "Z0": [_term(c, _q(c, Fraction(-5, 12)), "v5:chi:Z0", (0, 0, 0, 0), origin="v5-origin-text:chi:Z0", error_class="v5_chi")],
        "Z1": [_term(c, _q(c, Fraction(9, 16)), "v5:chi:Z1", (0, 0, 0, 0), origin={"tag": "v5-origin-map:chi:Z1", "source_ids": ["v5:chi:Z1"], "error_classes": ["v5_chi_map"]}, error_class="v5_chi")],
    }
    certificate = {"hb_sc_Z_independent": True, "rho_radius_Z_independent": True, "strict_hb_admitted": True, "fixture_geometry_owner": "v5_independent"}
    atoms = {name: "v5:" + name for name in SOURCE_POWER_BASIS}
    kwargs = {
        "incoming": source,
        "chi": chi,
        "hb": hb,
        "sc": sc,
        "source_family": "synthetic_v5_inlet",
        "source_atoms": atoms,
        "geometry_certificate": certificate,
    }
    result = build(c, **kwargs)

    import mpmath as mp

    mp.mp.dps = 100
    hb_mp, sc_mp = mp.mpf(1) / 4, mp.mpf(1) / 96
    bases = (hb_mp, mp.mpf(7) / 5, mp.mpf(13) / 6, mp.mpf(11) / 4)
    chi_values = (mp.mpf(-5) / 12, mp.mpf(9) / 16)
    incoming_values = {
        "Fc": (_formal(mp, mp.mpf(5) / 19, SOURCE_POWERS["Fc"], bases), _formal(mp, -mp.mpf(7) / 26, SOURCE_POWERS["Fc"], bases)),
        "Fxc": (_formal(mp, -mp.mpf(11) / 32, SOURCE_POWERS["Fxc"], bases), _formal(mp, mp.mpf(9) / 35, SOURCE_POWERS["Fxc"], bases)),
    }
    checks: dict[str, bool] = {
        "schema": result.get("schema") == "n1_collar_inlet_potential_v5",
        "formal_atoms": tuple(result.get("source_power_basis", ())) == SOURCE_POWER_BASIS,
        "four_branches": tuple(result.get("branch_names", ())) == BRANCHES,
        "geometry_fixed": result.get("geometry", {}).get("ordinary_Z_geometry_fixed") is True,
        "geometry_not_differentiated": result.get("geometry", {}).get("ordinary_Z_geometry_differentiated") is False,
        "geometry_certificate": result.get("geometry_certificate", {}).get("strict_hb_admitted") is True,
        "chi_Z_required": result.get("chi_Z_required") is True,
        "no_omissions": result.get("baseline_K0_reentered") is False and result.get("allword_tail_added") is False and result.get("m15_computed") is False and result.get("j2_kernel_remainder_included") is False,
    }
    branches = result["branches"]
    for branch in BRANCHES:
        component = "Fc" if branch.endswith("Fc") else "Fxc"
        source_hb = SOURCE_POWERS[component][0]
        terms0, terms1 = branches[branch]["Z0"], branches[branch]["Z1"]
        checks[branch + ":nonempty"] = bool(terms0) and len(terms1) >= 2
        checks[branch + ":hb_increment_relative"] = all(term["powers"][0] - source_hb == HB_ADDED[branch] for term in terms0 + terms1)
        checks[branch + ":signed_Z0_id"] = all(term.get("error_only") is True or ("v5:chi:Z0" in term["source_ids"] and f"v5:{component}:Z0" in term["source_ids"]) for term in terms0)
        checks[branch + ":exact_Z1_route_ids"] = all(f"v5:chi:Z{term['chi_order']}" in term["source_ids"] and (f"v5:{component}:Z{term['incoming_order']}:error" if term.get("error_only") else f"v5:{component}:Z{term['incoming_order']}") in term["source_ids"] for term in terms1)
        checks[branch + ":both_Z1_routes"] = any(term["chi_order"] == 1 and term["incoming_order"] == 0 for term in terms1) and any(term["chi_order"] == 0 and term["incoming_order"] == 1 for term in terms1)
        checks[branch + ":error_classes"] = all("v5_chi" in term["error_classes"] for term in terms0 + terms1)
        nominal0, nominal1 = incoming_values[component]
        target0 = _integral(mp, branch, hb_mp, sc_mp, chi_values[0], nominal0)
        target1 = _integral(mp, branch, hb_mp, sc_mp, chi_values[1], nominal0) + _integral(mp, branch, hb_mp, sc_mp, chi_values[0], nominal1)
        got0, got1 = _sector(mp, terms0, bases), _sector(mp, terms1, bases)
        checks[branch + ":analytic_Z0"] = got0[0] <= target0 <= got0[1]
        checks[branch + ":analytic_Z1"] = got1[0] <= target1 <= got1[1]

    incoming_jets = result.get("incoming_jets", {})
    fc0 = incoming_jets.get("Fc", {}).get("Z0", [])[0]
    fc1 = incoming_jets.get("Fc", {}).get("Z1", [])[0]
    fxc0 = incoming_jets.get("Fxc", {}).get("Z0", [])[0]
    fxc1 = incoming_jets.get("Fxc", {}).get("Z1", [])[0]
    checks["origin_string_preserved"] = fc0.get("origin") == "v5-origin-text:Fc:Z0" and fxc1.get("origin") == "v5-origin-text:Fxc:Z1"
    checks["origin_mapping_preserved"] = isinstance(fc1.get("origin"), dict) and fc1["origin"].get("tag") == "v5-origin-map:Fc:Z1" and isinstance(fxc0.get("origin"), dict) and fxc0["origin"].get("tag") == "v5-origin-map:Fxc:Z0"
    checks["origin_source_metadata"] = "v5:Fc:Z1" in fc1["origin"].get("source_ids", []) and "v5:Fxc:Z0" in fxc0["origin"].get("source_ids", [])
    checks["output_origin_string_metadata"] = any("v5-origin-text:chi:Z0" in term.get("origin_texts", []) and term.get("origin", {}).get("chi_origin") == "v5-origin-text:chi:Z0" for term in branches["F_from_Fc"]["Z0"])
    checks["output_origin_mapping_metadata"] = any(isinstance(term.get("origin", {}).get("incoming_origin"), dict) and term["origin"]["incoming_origin"].get("tag") == "v5-origin-map:Fc:Z1" for term in branches["F_from_Fc"]["Z1"])
    checks["error_only_symmetric"] = any(term.get("error_only") is True and term.get("symmetric_error") is True and term["coefficient"].a <= 0 <= term["coefficient"].b and "v5_symmetric_tail" in term["error_classes"] for term in branches["F_from_Fxc"]["Z0"])

    closure_kwargs = dict(kwargs)
    closure_kwargs["hb"] = _box(c, 0, Fraction(1, 2))
    try:
        closure = build(c, **closure_kwargs)
    except (TypeError, ValueError):
        checks["certified_closed_hb"] = False
        checks["closed_hb_formal"] = False
    else:
        checks["certified_closed_hb"] = True
        checks["closed_hb_formal"] = closure["geometry"]["formal_hb_not_materialized"] is True and closure["geometry"]["closed_hb_enclosure_allowed"] is True
    weak = dict(certificate)
    weak.pop("strict_hb_admitted", None)
    try:
        build(c, **dict(kwargs, geometry_certificate=weak))
    except (TypeError, ValueError):
        checks["missing_strict_certificate_rejected"] = True
    else:
        checks["missing_strict_certificate_rejected"] = False

    grouped = build(c, **dict(kwargs, incoming=None, inlet_groups={"G0g": source, "core_Dg": source, "core_symmetric_tail": source}))
    names = tuple(grouped.get("inlet_groups", ()))
    checks["three_groups"] = names == ("G0g", "core_Dg", "core_symmetric_tail")
    checks["groups_not_merged"] = set(grouped.get("group_results", {})) == set(names) and all(f"{name}:F_from_Fc" in grouped["all_branches"] for name in names)
    checks["group_error_only_symmetric"] = any(term.get("error_only") is True and term.get("symmetric_error") is True and term["coefficient"].a <= 0 <= term["coefficient"].b for term in grouped["group_results"]["core_symmetric_tail"]["F_from_Fxc"]["Z0"])

    return {"passed": all(checks.values()), "checks": checks, "check_count": len(checks), "synthetic_only": True, "independent_high_precision_integrals": True, "relative_hb_deltas_checked": True, "source_id_Z1_routes_checked": True, "origin_string_and_mapping_checked": True, "closed_hb_certificate_checked": True, "j2_remainder_future": True, "receipt_reads": False, "old_builder_called": False, "automatic_execution": False}


run = run_fixture

__all__ = ["SOURCE_POWER_BASIS", "BRANCHES", "SOURCE_POWERS", "run_fixture", "run"]
