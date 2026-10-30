"""Independent synthetic fixture for the finite E51 compensation leaf."""

from __future__ import annotations

from fractions import Fraction
from typing import Any


FIELDS = ("F1", "d_x_F1")


def _q(c: Any, value: Fraction) -> Any:
    return c.mpf(value.numerator) / value.denominator


def _t(c: Any, index: int) -> Any:
    left = _q(c, Fraction(index, 8))
    right = _q(c, Fraction(index + 1, 8))
    return c.make_mpf((left._mpi_[0], right._mpi_[1]))


def _term(c: Any, coefficient: int, powers: tuple[int, int, int, int],
          source: str, *, uncertainty: bool = False) -> dict[str, Any]:
    classes = ["synthetic_uncertainty"] if uncertainty else ["synthetic_signed"]
    return {
        "powers": powers,
        "coefficient": c.mpf(coefficient),
        "source_id": source,
        "source_ids": [source],
        "error_classes": classes,
        "uncertainty": uncertainty,
        "origin": {
            "source_ids": [source], "error_classes": classes,
            "uncertainty": uncertainty,
        },
    }


def _inputs(c: Any) -> tuple[list[dict[str, Any]], list[dict[str, Any]],
                              dict[str, list[dict[str, Any]]]]:
    w0: list[dict[str, Any]] = []
    delta: list[dict[str, Any]] = []
    for cell in range(8):
        t = _t(c, cell)
        w0.append({
            "t": t,
            "F1": {
                "Z0": [_term(c, (-1 if cell % 3 == 0 else 1) * (cell + 1),
                              (0, 1, cell % 2, -1), f"W0:{cell}:Z0")],
                "Z1": [_term(c, (-1 if cell % 2 else 1) * (cell + 2),
                              (1, 0, (cell + 1) % 2, -2), f"W0:{cell}:Z1")],
            },
        })
        delta.append({
            "t": t,
            "F1": {
                "Z0": [_term(c, (-1 if cell % 2 else 1) * (2 + cell),
                              (1, 0, 0, -2), f"Delta:{cell}:Z0")],
                "Z1": [_term(c, (-1 if cell % 3 == 1 else 1) * (3 + cell),
                              (0, 1, 1, -3), f"Delta:{cell}:Z1",
                              uncertainty=(cell == 6))],
            },
        })
    chi = {
        "Z0": [_term(c, 2, (0, 0, 1, 0), "chi:Z0")],
        "Z1": [_term(c, -1, (0, 0, 0, 1), "chi:Z1")],
    }
    return w0, delta, chi


def _endpoint(bits: Any) -> Fraction:
    sign, mantissa, exponent, _ = bits
    value = Fraction(int(mantissa), 1)
    value = value * (1 << exponent) if exponent >= 0 else value / (1 << (-exponent))
    return -value if sign else value


def _contains(interval: Any, expected: Fraction) -> bool:
    return _endpoint(interval._mpi_[0]) <= expected <= _endpoint(interval._mpi_[1])


def run_fixture(build: Any, c: Any) -> dict[str, Any]:
    """Run synthetic nonconstant, signed, Z1, axis, and diagonal checks."""

    w0, delta, chi = _inputs(c)
    result = build(c, w0, delta, chi, include_d_x_F1=True)
    branches = result.get("branches", {})
    checks: dict[str, bool] = {}
    expected_branches = {"G0_P_W0", "G0_P_DeltaW0", "D_P_W0", "D_P_DeltaW0"}
    checks["four_branches"] = set(branches) == expected_branches
    checks["all_target_profiles"] = all(
        len(branches[name].get("target_cells", [])) == 8
        for name in expected_branches
    ) if checks["four_branches"] else False
    checks["axis_and_diagonal"] = all(
        branches[name]["target_cells"][0]["axis_endpoint_continuation"] is True
        and branches[name]["target_cells"][0]["triangular_diagonal_source"] == 0
        and branches[name]["target_cells"][1]["triangular_diagonal_source"] == 1
        for name in expected_branches
    ) if checks["all_target_profiles"] else False

    g0 = branches.get("G0_P_W0", {})
    g0_axis = g0.get("target_cells", [{}])[0].get("by_source", {}).get("0", {})
    g0_axis_mass = g0_axis.get("kernel_masses", {}).get("F1")
    checks["G0_axis_exact_mass"] = (
        g0_axis_mass is not None and _contains(g0_axis_mass, Fraction(1, 128))
    )
    d = branches.get("D_P_DeltaW0", {})
    d_axis = d.get("target_cells", [{}])[0].get("by_source", {}).get("0", {})
    d_axis_mass = d_axis.get("kernel_masses", {}).get("F1")
    checks["D_axis_exact_physical_mass"] = (
        d_axis_mass is not None and _contains(d_axis_mass, Fraction(1, 49152))
    )
    exit_mass = result.get("G0_P_W0_exit_x2_by_source_cell", {}).get("0", {}).get("kernel_masses", {}).get("F1")
    checks["fixed_x2_exact_mass"] = (
        exit_mass is not None and _contains(exit_mass, Fraction(127, 8192))
    )

    p_terms = g0.get("P_q5_by_source_cell", {}).get("1", {}).get(1, [])
    checks["P_Z1_chi_leibniz"] = (
        any("chi_Z*V" in str(term.get("origin")) for term in p_terms)
        and any("chi*V_Z" in str(term.get("origin")) for term in p_terms)
        and all(term.get("chi_level") == "P" for term in p_terms)
    )
    d_terms = d.get("exit_by_source_cell", {}).get("0", {}).get("fields", {}).get("Z1", {}).get("F1", [])
    checks["D_second_chi_level"] = bool(d_terms) and all(
        len(term.get("chi_levels", [])) == 2
        and term.get("chi_level") == "D"
        and term.get("P_order") in (0, 1)
        for term in d_terms
    )
    checks["four_powers_ids"] = all(
        len(term.get("powers", ())) == 4 and term.get("source_ids")
        for name in expected_branches
        for cell in branches[name].get("target_cells", [])
        for source in cell.get("by_source", {}).values()
        for order in ("Z0", "Z1")
        for field in FIELDS
        for term in source.get("fields", {}).get(order, {}).get(field, [])
    )
    checks["signed_and_uncertain_preserved"] = (
        any(term.get("coefficient").a < 0
            for term in d_terms)
        and any(term.get("uncertainty") is True
                for term in branches["G0_P_DeltaW0"].get("target_cells", [])[6]
                .get("by_source", {}).get("6", {}).get("fields", {})
                .get("Z1", {}).get("F1", []))
    )
    checks["finite_only"] = (
        result.get("finite_compensation_only") is True
        and result.get("complete_feedback") is False
        and result.get("complete_Kchi_W1_inlet") is False
        and result.get("raw_C0_consumed") is False
    )
    # Independent continuous fixed-exit integrals for a signed source cell;
    # evaluate all formal source bases at one only for this synthetic check.
    import mpmath as mp
    mp.mp.dps = 80
    x=mp.mpf(2); a=mp.mpf(1)/4; b=mp.mpf(1)/2; tol=mp.mpf('1e-65')
    for name,v0,v1 in (("G0_P_W0",2,-3),("G0_P_DeltaW0",-3,-4),
                       ("D_P_W0",2,-3),("D_P_DeltaW0",-3,-4)):
        for order in (0,1):
            q0=4*v0; q1=-2*v0+4*v1
            def integrand(s):
                k0=s*(1-(s/x)**2)/2
                k1=s*(x*x-s**4/(x*x)-4*s*s*mp.log(x/s))/16
                return ((q0 if order==0 else q1)*k0 if name.startswith('G0_')
                        else (-4*q0 if order==0 else 2*q0-4*q1)*k1)
            ref=mp.quad(integrand,[a,b])
            terms=branches[name]['exit_by_source_cell']['1']['fields']['Z'+str(order)]['F1']
            enclosure=sum((term['coefficient'] for term in terms),c.mpf(0))
            lo=mp.make_mpf(enclosure._mpi_[0]);hi=mp.make_mpf(enclosure._mpi_[1])
            checks[f'continuous_exit_{name}_Z{order}']=bool(lo-tol<=ref<=hi+tol)
    return {
        "passed": all(checks.values()),
        "checks": checks,
        "check_count": len(checks),
        "synthetic_only": True,
        "nonconstant_source": True,
        "signed_delta_source": True,
        "ordinary_Z1_source": True,
        "method": "independent P=+2chi E51 with G0/D target-cell and fixed-x=2 masses",
        "limitations": {
            "fields": list(FIELDS),
            "complete_feedback": False,
            "raw_C0_B1": False,
            "W0_tail_sum": False,
        },
    }


run = run_fixture


__all__ = ["run_fixture", "run"]
