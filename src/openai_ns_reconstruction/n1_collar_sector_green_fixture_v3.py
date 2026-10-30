"""Independent v3 signed-rational fixture for the first-collar Dg leaf.

This fixture uses a different signed source/potential packet from the earlier
synthetic checks.  The source carries negative ``hb`` power, nontrivial
``F0``/``Lambda`` ownership, and an explicit ordinary Z1.  The potential has
opposite-sign Z0/Z1 values so the two signed product rows can be inspected
separately.  The reference is an independent high-precision nested phase
integral; it does not import or invoke another fixture or source producer.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Any


FIELDS = ("F1", "d_x_F1")
SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
SOURCE_POWERS = (-1, 2, 0, 1)


def _iv(c: Any, value: Fraction | int | str) -> Any:
    if isinstance(value, Fraction):
        return c.mpf(value.numerator) / value.denominator
    return c.mpf(value)


def _term(c: Any, coefficient: Fraction, source_id: str,
          powers: tuple[int, int, int, int]) -> dict[str, Any]:
    return {
        "powers": list(powers),
        "coefficient": _iv(c, coefficient),
        "origin": {
            "source_ids": [source_id],
            "uncertainty": True,
            "error_classes": ["v3_signed_rational"] ,
        },
    }


def _formal_value(c: Any, terms: Any, bases: tuple[Any, Any, Any, Any]) -> Any:
    total = c.mpf(0)
    for term in terms:
        value = term["coefficient"]
        for base, power in zip(bases, term["powers"]):
            value *= base ** int(power)
        total += value
    return total


def _phi(mp: Any, z: Any) -> Any:
    return mp.mpf(1) if z == 0 else -mp.expm1(-z) / z


def _nested_reference(mp: Any, *, hb: Any, sc: Any,
                      source0: Any, source1: Any,
                      chi0: Any, chi1: Any) -> dict[str, Any]:
    cphase = sc / 8

    def xx(q: Any) -> Any:
        return 2 * mp.exp(hb * cphase * q)

    def k0(q: Any, r: Any) -> Any:
        return hb**2 * cphase**2 * xx(r)**2 * (q - r) * _phi(
            mp, 2 * hb * cphase * (q - r)
        )

    def kx(q: Any, r: Any) -> Any:
        return hb * cphase * xx(r) * mp.exp(-3 * hb * cphase * (q - r))

    mass_f = mp.quad(
        lambda r: k0(1, r) * mp.quad(lambda p: k0(r, p), [0, r]),
        [0, 1],
    )
    mass_x = mp.quad(
        lambda r: kx(1, r) * mp.quad(lambda p: k0(r, p), [0, r]),
        [0, 1],
    )
    return {
        "Z0": {
            "F1": -2 * chi0 * source0 * mass_f,
            "d_x_F1": -2 * chi0 * source0 * mass_x,
        },
        "Z1": {
            "F1": -2 * (chi1 * source0 + chi0 * source1) * mass_f,
            "d_x_F1": -2 * (chi1 * source0 + chi0 * source1) * mass_x,
        },
    }


def run(leaf: Any, c: Any) -> dict[str, Any]:
    """Run independent v3 signed, geometry, derivative, and guard checks."""

    import mpmath as mp

    hb_interval = c.mpf([c.mpf(1) / 5, c.mpf(2) / 5])
    sc = _iv(c, Fraction(1, 1024))
    cells = [(c.mpf(0), c.mpf(1))]
    source = {
        "Z0": [_term(c, Fraction(5, 13), "v3:g5:Z0", SOURCE_POWERS)],
        "Z1": [_term(c, Fraction(-4, 17), "v3:g5:Z1", SOURCE_POWERS)],
    }
    chi = {
        "Z0": [_term(c, Fraction(-2, 7), "v3:chi:Z0", (0, 0, 0, 0))],
        "Z1": [_term(c, Fraction(3, 11), "v3:chi:Z1", (0, 0, 0, 0))],
    }
    certificate = {
        "hb_sc_Z_independent": True,
        "rho_radius_Z_independent": True,
        "alternate_v3_geometry_owner": "synthetic_v3",
    }
    atoms = {name: f"v3:{name}" for name in SOURCE_POWER_BASIS}
    kwargs = {
        "hb": hb_interval,
        "sc": sc,
        "g5": source,
        "chi": chi,
        "cells": cells,
        "q_star": 1,
        "source_family": "synthetic_v3_signed_rational",
        "source_atoms": atoms,
        "certificate": certificate,
    }
    result = leaf.build(c, **kwargs)

    mp.mp.dps = 100
    hb_mp = mp.mpf(3) / 10
    sc_mp = mp.mpf(1) / 1024
    f0_mp = mp.mpf(7) / 3
    pstar_squared_mp = mp.mpf(5) / 2
    lambda_mp = mp.mpf(11) / 4
    source0 = mp.mpf(5) / 13 * hb_mp**SOURCE_POWERS[0] * f0_mp**SOURCE_POWERS[1] * pstar_squared_mp**SOURCE_POWERS[2] * lambda_mp**SOURCE_POWERS[3]
    source1 = mp.mpf(-4) / 17 * hb_mp**SOURCE_POWERS[0] * f0_mp**SOURCE_POWERS[1] * pstar_squared_mp**SOURCE_POWERS[2] * lambda_mp**SOURCE_POWERS[3]
    chi0 = -mp.mpf(2) / 7
    chi1 = mp.mpf(3) / 11
    reference = _nested_reference(
        mp,
        hb=hb_mp,
        sc=sc_mp,
        source0=source0,
        source1=source1,
        chi0=chi0,
        chi1=chi1,
    )
    bases = (
        c.mpf(3) / 10,
        c.mpf(7) / 3,
        c.mpf(5) / 2,
        c.mpf(11) / 4,
    )

    checks: dict[str, bool] = {
        "schema": result.get("schema") == "n1_collar_sector_green_v3",
        "whole_domain_box": result.get("whole_domain_source_box") is True,
        "source_partition_not_required": result.get("source_partition_required") is False,
        "geometry_certificate": (
            result.get("geometry_certificate", {}).get("hb_sc_Z_independent") is True
            and result.get("geometry_certificate", {}).get("rho_radius_Z_independent") is True
        ),
        "fixed_Z_geometry": result.get("geometry", {}).get("ordinary_Z_geometry_fixed") is True,
        "exact_geometry_formula": result.get("geometry", {}).get("dx_dq_formula") == "hb*sc*x(q)/8",
        "source_row_5": result.get("source_row") == "5",
        "formal_atoms": tuple(result.get("source_power_basis", ())) == SOURCE_POWER_BASIS,
        "baseline_not_reentered": result.get("baseline_K0_reentered") is False,
        "no_tail": result.get("allword_tail_added") is False,
    }

    response = result.get("response", {})
    for field in FIELDS:
        hb_increment = 4 if field == "F1" else 3
        expected_powers = (
            SOURCE_POWERS[0] + hb_increment,
            SOURCE_POWERS[1], SOURCE_POWERS[2], SOURCE_POWERS[3],
        )
        for order in (0, 1):
            terms = response.get(field, {}).get(f"Z{order}", [])
            checks[f"nonempty_{field}_Z{order}"] = bool(terms)
            checks[f"powers_{field}_Z{order}"] = all(
                tuple(term.get("powers", ())) == expected_powers for term in terms
            )
            checks[f"provenance_{field}_Z{order}"] = all(
                "v3:g5:Z0" in term.get("origin", {}).get("source_ids", [])
                or "v3:g5:Z1" in term.get("origin", {}).get("source_ids", [])
                for term in terms
            )
            got = _formal_value(c, terms, bases)
            target = c.mpf(mp.nstr(reference[f"Z{order}"][field], 85))
            checks[f"nested_integral_{field}_Z{order}"] = bool(
                got.a <= target.a and got.b >= target.b
            )
            if order == 0:
                checks[f"signed_Z0_{field}"] = bool(got.a > 0)
            else:
                checks[f"signed_Z1_{field}"] = bool(got.b < 0)

    checks["chi_Z_branch"] = any(
        term.get("chi_order") == 1 and term.get("source_order") == 0
        for term in response.get("F1", {}).get("Z1", [])
    )
    checks["source_Z1_branch"] = any(
        term.get("chi_order") == 0 and term.get("source_order") == 1
        for term in response.get("F1", {}).get("Z1", [])
    )
    checks["measure_once"] = all(
        term.get("physical_dx_count") == 2
        and term.get("measure_attached_once") is True
        and term.get("coefficient_contains_formal_hb") is False
        for terms in response.get("F1", {}).values()
        for term in terms
    )

    try:
        probe = leaf.endpoint_kernel_masses(
            c,
            hb=hb_interval,
            sc=sc,
            outer_q=(c.mpf(0), c.mpf(1)),
            source_q=(c.mpf(0), c.mpf(1)),
            q_star=1,
        )
        checks["scalar_mpf_endpoint_repair"] = bool(
            probe["F1"].b < c.mpf(1) and probe["d_x_F1"].b < c.mpf(1)
        )
    except (TypeError, ValueError):
        checks["scalar_mpf_endpoint_repair"] = False

    def rejects(**overrides: Any) -> bool:
        candidate = dict(kwargs)
        candidate.update(overrides)
        try:
            leaf.build(c, **candidate)
        except (TypeError, ValueError):
            return True
        return False

    checks["guard_missing_certificate"] = rejects(certificate=None)
    checks["guard_missing_radius_flag"] = rejects(
        certificate={"hb_sc_Z_independent": True},
    )
    checks["guard_missing_chi_Z"] = rejects(chi={"Z0": chi["Z0"]})
    checks["guard_missing_source_Z1"] = rejects(g5={"Z0": source["Z0"]})
    checks["guard_nonendpoint"] = rejects(q_star=c.mpf(99) / 100)

    return {
        "passed": all(checks.values()),
        "check_count": len(checks),
        "checks": checks,
        "synthetic_only": True,
        "reference": "independent 100-digit nested K0 o K0 phase integral at q*=1",
        "hb_interval": "[1/5,2/5]",
        "hb_evaluation": "3/10",
        "sc": "1/1024",
        "source_powers": SOURCE_POWERS,
        "baseline_called": False,
    }


__all__ = ["FIELDS", "SOURCE_POWER_BASIS", "SOURCE_POWERS", "run"]
