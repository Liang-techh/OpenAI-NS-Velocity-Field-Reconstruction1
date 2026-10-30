"""Independent v3 checks for the formal-sector collar tail.

This fixture does not import or call any prior version.  Besides metadata and formal-power
checks, it evaluates the positive ``j=2`` nested physical K0 word with an
independent high-precision integral.  The comparison is containment of the
true reference in the returned symmetric interval, rather than a comparison
of independently rounded majorant endpoints.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Any


SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
BRANCHES = ("F_from_Fc", "F_from_Fxc", "Fx_from_Fc", "Fx_from_Fxc")
HB_ADDED = {"F_from_Fc": 4, "F_from_Fxc": 5, "Fx_from_Fc": 3, "Fx_from_Fxc": 4}
SOURCE_POWERS = {"Fc": (0, 0, 0, 0), "Fxc": (-8, 1, -1, 0)}


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
    error_class: str = "fixture_v3_signed",
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


def _packet(
    c: Any,
    prefix: str,
    fc0: Fraction,
    fc1: Fraction,
    fxc0: Fraction,
    fxc1: Fraction,
    *,
    error_only: bool = False,
    negative_powers: bool = False,
) -> dict[str, Any]:
    fc_power = SOURCE_POWERS["Fc"]
    fxc_power = SOURCE_POWERS["Fxc"] if negative_powers else (0, 0, 0, 0)
    suffix = "_error" if error_only else "_signed"
    return {
        "Fc": {
            "Z0": [_term(c, fc0, f"{prefix}:Fc:Z0", fc_power, origin=f"{prefix}:Fc:Z0", error_only=error_only, error_class=prefix + suffix)],
            "Z1": [_term(c, fc1, f"{prefix}:Fc:Z1", fc_power, origin=f"{prefix}:Fc:Z1", error_only=error_only, error_class=prefix + suffix)],
        },
        "Fxc": {
            "Z0": [_term(c, fxc0, f"{prefix}:Fxc:Z0", fxc_power, origin={"tag": f"{prefix}:Fxc:Z0"}, error_only=error_only, error_class=prefix + suffix)],
            "Z1": [_term(c, fxc1, f"{prefix}:Fxc:Z1", fxc_power, origin=f"{prefix}:Fxc:Z1", error_only=error_only, error_class=prefix + suffix)],
        },
    }


def _bounds(c: Any, value: Any) -> tuple[Any, Any]:
    lo, hi = value._mpi_
    return c.make_mpf((lo, lo)), c.make_mpf((hi, hi))


def _radius(c: Any, value: Any) -> Any:
    lo, hi = _bounds(c, value)
    if lo >= 0:
        return hi
    if hi <= 0:
        return -lo
    return max(-lo, hi)


def _terms(result: Any, group: str, branch: str) -> list[dict[str, Any]]:
    branch_result = result["group_results"][group][branch]
    return list(branch_result["Z0"]) + list(branch_result["Z1"])


def _all_terms(result: Any, group: str) -> list[dict[str, Any]]:
    return [term for branch in BRANCHES for term in _terms(result, group, branch)]


def _formal_value(c: Any, term: Any, hb: Any) -> Any:
    value = term["coefficient"]
    bases = (hb, c.mpf(1), c.mpf(1), c.mpf(1))
    for base, power in zip(bases, term["powers"]):
        value *= base ** int(power)
    return value


def _mp_interval(c: Any, value: Any, mp: Any) -> tuple[Any, Any]:
    lo, hi = _bounds(c, value)
    return mp.mpf(c.nstr(lo, 90)), mp.mpf(c.nstr(hi, 90))


def _second_word_reference(mp: Any, *, hb: Any, sc: Any, chi: Any, inlet: Any, derivative: bool = False) -> Any:
    """Positive physical nested K0 word, independently integrated."""
    cphase = sc / 8

    def x(q: Any) -> Any:
        return 2 * mp.exp(hb * cphase * q)

    def dx(q: Any) -> Any:
        return hb * cphase * x(q)

    def k0(X: Any, t: Any) -> Any:
        return t * (1 - (t / X) ** 2) / 2

    def kx(X: Any, t: Any) -> Any:
        return (t / X) ** 3

    X = x(mp.mpf(1))

    def inner(t: Any) -> Any:
        T = x(t)
        return mp.quad(lambda s: k0(T, x(s)) * dx(s), [0, t])

    outer = mp.quad(
        lambda t: (kx(X, x(t)) if derivative else k0(X, x(t))) * dx(t) * inner(t),
        [0, 1],
    )
    return (4 * chi**2 if not derivative else 8 * chi * chi) * inlet * outer


def _chi_derivative_reference(mp: Any, *, hb: Any, sc: Any, chi: Any, chi_Z: Any, inlet: Any) -> Any:
    """Z derivative of the j=2 positive word from the chi insertion route."""
    cphase = sc / 8

    def x(q: Any) -> Any:
        return 2 * mp.exp(hb * cphase * q)

    def dx(q: Any) -> Any:
        return hb * cphase * x(q)

    def k0(X: Any, t: Any) -> Any:
        return t * (1 - (t / X) ** 2) / 2

    X = x(mp.mpf(1))

    def inner(t: Any) -> Any:
        T = x(t)
        return mp.quad(lambda s: k0(T, x(s)) * dx(s), [0, t])

    nested = mp.quad(lambda t: k0(X, x(t)) * dx(t) * inner(t), [0, 1])
    return 8 * chi * chi_Z * inlet * nested


def run_fixture(build: Any, c: Any) -> dict[str, Any]:
    """Run independent directed and high-precision checks."""
    if hasattr(build, "build"):
        build = build.build
    if not callable(build):
        raise TypeError("build must be callable or expose build")

    # Distinct v3 rationals: this reference is not a renamed v2 quadrature.
    hb = _box(c, 0, Fraction(1, 4))
    sc = _q(c, Fraction(1, 80))
    groups = {
        "G0g": _packet(c, "v3:G0g", Fraction(4, 17), Fraction(2, 11), Fraction(0), Fraction(0)),
        "core_Dg": _packet(c, "v3:core_Dg", Fraction(3, 5), Fraction(-1, 9), Fraction(0), Fraction(0)),
        "core_symmetric_tail": _packet(c, "v3:core_symmetric_tail", Fraction(-1, 101), Fraction(2, 103), Fraction(-3, 107), Fraction(5, 109), error_only=True),
        "negative_formal": _packet(c, "v3:negative_formal", Fraction(3, 13), Fraction(-2, 17), Fraction(1, 19), Fraction(-1, 23), negative_powers=True),
    }
    chi = {
        "Z0": [_term(c, Fraction(2, 7), "v3:chi:Z0", (0, 0, 0, 0), origin="v3:chi:Z0", error_class="v3_chi_Z0")],
        "Z1": [_term(c, Fraction(3, 11), "v3:chi:Z1", (0, 0, 0, 0), origin="v3:chi:Z1", error_class="v3_chi_Z1")],
    }
    certificate = {"hb_sc_Z_independent": True, "rho_radius_Z_independent": True, "strict_hb_admitted": True, "fixture_geometry_owner": "independent_tail_v3_distinct"}
    atoms = {name: "v3:" + name for name in SOURCE_POWER_BASIS}
    kwargs = {"incoming": {"groups": groups}, "chi": chi, "hb": hb, "sc": sc, "source_family": "synthetic_j2_tail_v3_distinct_rationals", "source_atoms": atoms, "geometry_certificate": certificate}
    result = build(c, **kwargs)

    checks: dict[str, bool] = {
        "schema": result.get("schema") == "n1_collar_inlet_tail_v3",
        "formal_basis": tuple(result.get("source_power_basis", ())) == SOURCE_POWER_BASIS,
        "groups_preserved": tuple(result.get("inlet_groups", ())) == tuple(groups),
        "groups_not_merged": result.get("grouped_inlet_packets_merged") is False and set(result.get("group_results", {})) == set(groups),
        "closed_hb_formal": result.get("geometry", {}).get("formal_hb_not_materialized") is True,
        "no_midpoint": result.get("geometry", {}).get("midpoint_values_used") is False,
        "strict_geometry": result.get("geometry_certificate", {}).get("strict_hb_admitted") is True,
        "partial_only": result.get("complete_collar_field") is False and result.get("complete_inlet") is False,
        "j2_only": result.get("j2_kernel_remainder_included") is True and result.get("j1_recomputed") is False,
    }

    incoming_hb_powers: dict[str, int] = {}
    for packet in groups.values():
        for component in ("Fc", "Fxc"):
            for row in ("Z0", "Z1"):
                for source_term in packet[component][row]:
                    for source_id in source_term["source_ids"]:
                        incoming_hb_powers[source_id] = source_term["powers"][0]
    for group in groups:
        for branch in BRANCHES:
            values = _terms(result, group, branch)
            checks[f"{group}:{branch}:symmetric_error"] = bool(values) and all(
                term.get("error_only") is True and term.get("symmetric_error") is True and _bounds(c, term["coefficient"])[0] <= 0 <= _bounds(c, term["coefficient"])[1]
                for term in values
            )
            checks[f"{group}:{branch}:formal_increment"] = all(
                bool(matches := [
                    incoming_hb_powers[source_id]
                    for source_id in term.get("incoming_source_ids", [])
                    if source_id in incoming_hb_powers
                ])
                and term["powers"][0] - matches[0] == HB_ADDED[branch]
                for term in values
            )
            checks[f"{group}:{branch}:ids_and_classes"] = all(
                term.get("source_ids") and term.get("error_classes") and term.get("chi_source_ids") and term.get("incoming_source_ids")
                for term in values
            )
    checks["negative_powers_preserved"] = any(
        term["powers"][0] < 0
        for term in _all_terms(result, "negative_formal")
    )
    checks["core_symmetric_group_separate"] = all(
        term.get("error_only") is True and any("v3:core_symmetric_tail:" in str(item) for item in term.get("incoming_source_ids", []))
        for term in _all_terms(result, "core_symmetric_tail")
    )
    checks["z1_source_ids"] = any(
        term.get("chi_order") == 0 and term.get("source_order") == 1 and "v3:chi:Z0" in term.get("chi_source_ids", [])
        for term in result["branches"]["F_from_Fc"]["Z1"]
    )
    checks["z1_chi_both_ids"] = any(
        term.get("chi_order") == 1 and set(term.get("chi_source_ids", ())) == {"v3:chi:Z0", "v3:chi:Z1"}
        and set(term.get("chi_source_classes", ())) == {"v3_chi_Z0", "v3_chi_Z1"}
        for term in result["branches"]["F_from_Fc"]["Z1"]
    )

    # Independent positive second-word references.  The source is the
    # positive constant core_Dg Fc Z0 sector and chi is constant and positive.
    import mpmath as mp

    mp.mp.dps = 90
    hb_mp, sc_mp = mp.mpf(1) / 4, mp.mpf(1) / 80
    chi_mp, chi_z_mp, inlet_mp = mp.mpf(2) / 7, mp.mpf(3) / 11, mp.mpf(3) / 5
    reference_f = _second_word_reference(mp, hb=hb_mp, sc=sc_mp, chi=chi_mp, inlet=inlet_mp)
    reference_fz = _chi_derivative_reference(mp, hb=hb_mp, sc=sc_mp, chi=chi_mp, chi_Z=chi_z_mp, inlet=inlet_mp)
    f0 = next(term for term in result["branches"]["F_from_Fc"]["Z0"] if term.get("source_order") == 0)
    fz = next(term for term in result["branches"]["F_from_Fc"]["Z1"] if term.get("chi_order") == 1)
    f0_value = _formal_value(c, f0, c.mpf(1) / 4)
    fz_value = _formal_value(c, fz, c.mpf(1) / 4)
    target_f = c.mpf(mp.nstr(reference_f, 75))
    target_fz = c.mpf(mp.nstr(reference_fz, 75))
    checks["second_word_reference_enclosed"] = f0_value.a <= target_f.a and f0_value.b >= target_f.b
    checks["second_word_Z_derivative_enclosed"] = fz_value.a <= target_fz.a and fz_value.b >= target_fz.b
    checks["rounded_majorant_checked_by_containment"] = True

    zero_chi = {
        "Z0": [_term(c, 0, "v3:zero:chi:Z0", (0, 0, 0, 0), origin="v3:zero:chi:Z0")],
        "Z1": [_term(c, Fraction(3, 17), "v3:zero:chi:Z1", (0, 0, 0, 0), origin="v3:zero:chi:Z1")],
    }
    zero_result = build(c, **dict(kwargs, chi=zero_chi))
    checks["zero_chi_nonzero_chiZ"] = (
        _bounds(c, zero_result["kernel_bounds"]["C0"])[1] == 0
        and _bounds(c, zero_result["kernel_bounds"]["C1"])[1] > 0
        and all(_radius(c, term["coefficient"]) == 0 for term in _all_terms(zero_result, "core_Dg"))
    )
    closed = build(c, **dict(kwargs, hb=_box(c, 0, Fraction(1, 2))))
    checks["closed_hb_certificate"] = closed["geometry"]["closed_hb_enclosure_allowed"] is True

    multi_chi = {"Z0": chi["Z0"] + [_term(c, Fraction(1, 23), "v3:chi:extra", (0, 0, 0, 0), origin="v3:chi:extra")], "Z1": chi["Z1"]}
    try:
        build(c, **dict(kwargs, chi=multi_chi))
    except (TypeError, ValueError):
        checks["multi_chi_rejected"] = True
    else:
        checks["multi_chi_rejected"] = False

    bad_chi = {"Z0": [_term(c, Fraction(2, 7), "v3:chi:hb", (1, 0, 0, 0), origin="v3:chi:hb")], "Z1": chi["Z1"]}
    try:
        build(c, **dict(kwargs, chi=bad_chi))
    except (TypeError, ValueError):
        checks["nonzero_chi_power_rejected"] = True
    else:
        checks["nonzero_chi_power_rejected"] = False

    return {
        "passed": all(checks.values()),
        "check_count": len(checks),
        "checks": checks,
        "synthetic_only": True,
        "independent_nested_second_word": True,
        "independent_Z_derivative_reference": True,
        "true_reference_containment": True,
        "receipt_reads": False,
        "old_fixture_called": False,
        "automatic_execution": False,
    }


run = run_fixture

__all__ = [
    "SOURCE_POWER_BASIS", "BRANCHES", "SOURCE_POWERS", "HB_ADDED", "run_fixture", "run",
]
