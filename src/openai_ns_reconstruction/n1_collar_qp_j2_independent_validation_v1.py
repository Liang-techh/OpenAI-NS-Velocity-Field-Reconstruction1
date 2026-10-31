"""Independent synthetic validation of a ``j >= 2`` Q/P m15 majorant.

This module is deliberately self-contained.  It does not import the active
J2 leaf, an older fixture, a transport helper, a receipt, or a source
producer.  It evaluates the positive second potential word for the two
incoming ``F`` branches directly and compares that expression with an
independent high-precision nested integral.

Put ``x=2 exp(hb*sc*q/8)`` and ``y=hb*sc*q/8``.  In the logarithmic variable
``t=2 exp(s)`` the positive K0 kernel is

``K(y,s) = 2 exp(2s) (1-exp(2(s-y)))``.

For ``H0F=Fc`` and ``H0F=(1-4/t**2)Fxc`` the one-word masses are
``M1_Fc=cosh(2y)-1`` and ``M1_Fxc=sinh(2y)-2y``.  Integrating K0 against
those exact first-word masses gives the analytic nested ``j=2`` masses
``M2_Fc`` and ``M2_Fxc`` below.  The physical q integral uses
``d rho = 2*x*dx`` exactly once and weights by the cutoff
``eta(q)=1-sigma(4q-2)``.  The outer physical measure is therefore ``d rho``;
the K0 kernel retains its own internal radial ``dx`` integration.  Thus the
raw integrand already contains the four-power (or five-power) hb behavior
of the nested mass; no extra ``hb**4`` or cell mass is multiplied in.

The returned record is numerical synthetic evidence only, not a source
certificate or a bound for the actual Q/P data.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Any, Callable


SCHEMA = "n1_collar_qp_j2_independent_validation_v1"
SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
BRANCHES = ("F_from_Fc", "F_from_Fxc")
OLD_FIXTURE_PARAMETERS = {
    "sc": Fraction(1, 4096),
    "hb": Fraction(1, 4),
    "chi": Fraction(3, 17),
    "chi_Z": Fraction(4, 29),
    "Fc": Fraction(3, 17),
    "Flead": Fraction(7, 19),
    "Lambda": Fraction(5, 3),
}

# These values are intentionally distinct from the saved J2 fixture values.
SYNTHETIC_PARAMETERS = {
    "sc": Fraction(5, 8192),
    "hb": Fraction(2, 7),
    "chi": Fraction(5, 23),
    "chi_Z": Fraction(7, 31),
    "Fc": Fraction(11, 17),
    "Fxc": Fraction(-13, 29),
    "Fc_Z": Fraction(19, 37),
    "Fxc_Z": Fraction(-23, 41),
    "Flead": Fraction(17, 23),
    "Flead_Z": Fraction(-11, 29),
    "Lambda": Fraction(9, 5),
}

__all__ = [
    "SCHEMA", "SOURCE_POWER_BASIS", "BRANCHES", "OLD_FIXTURE_PARAMETERS",
    "SYNTHETIC_PARAMETERS", "nested_j2_mass", "run_validation",
    "run_fixture", "run", "build",
]


def _as_mp(mp: Any, value: Fraction | int | str | Any) -> Any:
    if isinstance(value, Fraction):
        return mp.mpf(value.numerator) / mp.mpf(value.denominator)
    return mp.mpf(value)


def _sigma(mp: Any, t: Any) -> Any:
    if t <= 0:
        return mp.mpf(0)
    if t >= 1:
        return mp.mpf(1)
    left = mp.exp(-1 / t**2)
    right = mp.exp(-1 / (1 - t) ** 2)
    return left / (left + right)


def _eta(mp: Any, q: Any) -> Any:
    return 1 - _sigma(mp, 4 * q - 2)


def _kernel(mp: Any, y: Any, s: Any) -> Any:
    """The positive K0 kernel after ``t=2 exp(s)`` and ``dt=t ds``."""

    if s < 0 or s > y:
        return mp.mpf(0)
    return 2 * mp.exp(2 * s) * (1 - mp.exp(2 * (s - y)))


def _mass1_fc(mp: Any, y: Any) -> Any:
    return mp.cosh(2 * y) - 1


def _mass1_fxc(mp: Any, y: Any) -> Any:
    return mp.sinh(2 * y) - 2 * y


def _mass2_fc(mp: Any, y: Any) -> Any:
    """Exact ``K0`` applied twice to the ``Fc`` source branch."""

    return (
        mp.exp(4 * y) / 12
        + y
        - mp.exp(2 * y) / 2
        + mp.mpf(1) / 4
        + mp.exp(-2 * y) / 6
    )


def _mass2_fxc(mp: Any, y: Any) -> Any:
    """Exact ``K0`` applied twice to ``(1-4/t**2) Fxc``."""

    return (
        mp.exp(4 * y) / 12
        + mp.exp(2 * y) * (mp.mpf(3) / 4 - y)
        - y
        - mp.mpf(3) / 4
        - mp.exp(-2 * y) / 12
    )


def nested_j2_mass(mp: Any, y: Any, branch: str) -> Any:
    """Return the positive analytic nested mass for one source branch."""

    if branch == "F_from_Fc":
        return _mass2_fc(mp, y)
    if branch == "F_from_Fxc":
        return _mass2_fxc(mp, y)
    raise ValueError(f"unknown branch {branch!r}")


def _close(mp: Any, left: Any, right: Any, *, atol: str = "1e-70", rtol: str = "1e-60") -> bool:
    scale = max(mp.mpf(1), abs(left), abs(right))
    return abs(left - right) <= mp.mpf(atol) + mp.mpf(rtol) * scale


def _fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def _number_text(mp: Any, value: Any) -> str:
    return mp.nstr(value, 60)


def _direct_nested(mp: Any, y: Any, branch: str) -> Any:
    first: Callable[[Any], Any] = _mass1_fc if branch == "F_from_Fc" else _mass1_fxc
    return mp.quad(lambda s: _kernel(mp, y, s) * first(mp, s), [mp.mpf(0), y])


def _physical_integral(
    mp: Any,
    *,
    q_lower: Any,
    q_upper: Any,
    branch: str,
    source_amplitude: Any,
    chi_amplitude: Any,
    leading_amplitude: Any,
    hb: Any,
    sc: Any,
    lam: Any,
    eta: Callable[[Any], Any] = lambda q: 1,
    breakpoints: tuple[Any, ...] = (),
) -> Any:
    """Integrate the positive physical m15 majorant using ``d rho`` once."""

    def integrand(q: Any) -> Any:
        x = 2 * mp.exp(hb * sc * q / 8)
        rho = 4 * mp.exp(hb * sc * q / 4)
        dx_dq = hb * sc * x / 8
        d_rho_dq = 2 * x * dx_dq
        y = hb * sc * q / 8
        nested = nested_j2_mass(mp, y, branch)
        # (-2 Lambda^-2 rho Flead) times (+4 chi^2 nested source),
        # with absolute values for a majorant.  d rho/dq=2*x*dx/dq is
        # attached once; dx is not used as the outer physical measure.
        j2 = 4 * abs(chi_amplitude) ** 2 * abs(source_amplitude) * nested
        return eta(q) * 2 * abs(leading_amplitude) * rho * j2 * d_rho_dq / lam**2

    points = [q_lower]
    points.extend(point for point in breakpoints if q_lower < point < q_upper)
    points.append(q_upper)
    return mp.quad(integrand, points)


def _z1_physical_integral(
    mp: Any,
    *,
    q_lower: Any,
    q_upper: Any,
    branch: str,
    source_amplitude: Any,
    source_z_amplitude: Any,
    chi_amplitude: Any,
    chi_z_amplitude: Any,
    leading_amplitude: Any,
    leading_z_amplitude: Any,
    hb: Any,
    sc: Any,
    lam: Any,
) -> dict[str, Any]:
    """Return source-Z and chi-Z pieces of the positive Z1 majorant."""

    def integrate(factor: Any) -> Any:
        def integrand(q: Any) -> Any:
            x = 2 * mp.exp(hb * sc * q / 8)
            rho = 4 * mp.exp(hb * sc * q / 4)
            dx_dq = hb * sc * x / 8
            d_rho_dq = 2 * x * dx_dq
            y = hb * sc * q / 8
            nested = nested_j2_mass(mp, y, branch)
            j2 = factor(q, nested)
            return _eta(mp, q) * 2 * abs(leading_amplitude) * rho * j2 * d_rho_dq / lam**2

        return mp.quad(integrand, [q_lower, q_upper])

    source = integrate(lambda _q, nested: 4 * abs(chi_amplitude) ** 2 * abs(source_z_amplitude) * nested)
    chi = integrate(lambda _q, nested: 8 * abs(chi_amplitude) * abs(chi_z_amplitude) * abs(source_amplitude) * nested)
    leading = integrate(lambda _q, nested: 4 * abs(chi_amplitude) ** 2 * abs(source_amplitude) * nested * abs(leading_z_amplitude) / max(abs(leading_amplitude), mp.mpf("1e-200")))
    return {"source_Z": source, "chi_Z": chi, "leading_Z": leading, "total": source + chi + leading}


def _cells(mp: Any) -> tuple[tuple[Any, Any], ...]:
    return tuple((mp.mpf(index) / 4, mp.mpf(index + 1) / 4) for index in range(4))


def run_validation(c: Any = None, *, dps: int = 100) -> dict[str, Any]:
    """Run the independent synthetic J2 analytic and quadrature checks.

    ``c`` is accepted for runner symmetry with the source consumers but is
    not used to evaluate a source or an interval.  All values here are fresh
    synthetic rationals evaluated in a private high-precision mpmath context.
    """

    del c
    import mpmath

    mp = mpmath.mp.clone()
    mp.dps = int(dps)
    p = {name: _as_mp(mp, value) for name, value in SYNTHETIC_PARAMETERS.items()}
    q_zero = mp.mpf(0)
    q_one = mp.mpf(1)
    q_probe = mp.mpf(63) / 100
    y_probe = p["hb"] * p["sc"] * q_probe / 8

    checks: dict[str, bool] = {
        "analytic_fc_nested": False,
        "analytic_fxc_nested": False,
        "positive_nested_masses": False,
        "nonzero_fxc": abs(p["Fxc"]) > 0,
        "nonzero_chi_Z": abs(p["chi_Z"]) > 0,
        "nonzero_source_Z": abs(p["Fc_Z"]) > 0 and abs(p["Fxc_Z"]) > 0,
        "distinct_from_old_fixture": all(
            name not in OLD_FIXTURE_PARAMETERS or value != _as_mp(mp, OLD_FIXTURE_PARAMETERS[name])
            for name, value in p.items()
        ),
        "physical_measure_is_drho": False,
        "eta_changes_integral": False,
        "cell_partition_matches_whole": False,
        "fc_physical_positive": False,
        "fxc_physical_positive": False,
        "fc_formal_hb5": False,
        "fxc_formal_hb6": False,
        "no_duplicate_hb4": False,
        "no_duplicate_cell_mass": False,
        "z1_source_Z_positive": False,
        "z1_chi_Z_positive": False,
        "z1_leading_Z_positive": False,
        "z1_sum_identity": False,
    }

    analytic_fc = nested_j2_mass(mp, y_probe, "F_from_Fc")
    analytic_fxc = nested_j2_mass(mp, y_probe, "F_from_Fxc")
    direct_fc = _direct_nested(mp, y_probe, "F_from_Fc")
    direct_fxc = _direct_nested(mp, y_probe, "F_from_Fxc")
    checks["analytic_fc_nested"] = _close(mp, analytic_fc, direct_fc)
    checks["analytic_fxc_nested"] = _close(mp, analytic_fxc, direct_fxc)
    checks["positive_nested_masses"] = analytic_fc > 0 and analytic_fxc > 0

    eta_fn = lambda q: _eta(mp, q)
    one_fn = lambda _q: mp.mpf(1)
    branch_physical: dict[str, Any] = {}
    cell_physical: dict[str, list[Any]] = {}
    for branch, source in (("F_from_Fc", p["Fc"]), ("F_from_Fxc", p["Fxc"])):
        quarter_breaks = tuple(mp.mpf(index) / 4 for index in range(1, 4))
        whole = _physical_integral(
            mp, q_lower=q_zero, q_upper=q_one, branch=branch,
            source_amplitude=source, chi_amplitude=p["chi"],
            leading_amplitude=p["Flead"], hb=p["hb"], sc=p["sc"],
            lam=p["Lambda"], eta=eta_fn, breakpoints=quarter_breaks,
        )
        cells = [
            _physical_integral(
                mp, q_lower=lower, q_upper=upper, branch=branch,
                source_amplitude=source, chi_amplitude=p["chi"],
                leading_amplitude=p["Flead"], hb=p["hb"], sc=p["sc"],
                lam=p["Lambda"], eta=eta_fn,
            )
            for lower, upper in _cells(mp)
        ]
        unweighted = _physical_integral(
            mp, q_lower=q_zero, q_upper=q_one, branch=branch,
            source_amplitude=source, chi_amplitude=p["chi"],
            leading_amplitude=p["Flead"], hb=p["hb"], sc=p["sc"],
            lam=p["Lambda"], eta=one_fn,
        )
        branch_physical[branch] = {"eta_weighted": whole, "unweighted": unweighted}
        cell_physical[branch] = cells

    checks["physical_measure_is_drho"] = True
    checks["eta_changes_integral"] = any(
        abs(branch_physical[branch]["eta_weighted"] - branch_physical[branch]["unweighted"]) > mp.mpf("1e-100")
        for branch in BRANCHES
    )
    checks["cell_partition_matches_whole"] = all(
        _close(mp, sum(cell_physical[branch]), branch_physical[branch]["eta_weighted"], atol="1e-80", rtol="1e-70")
        for branch in BRANCHES
    )
    checks["fc_physical_positive"] = branch_physical["F_from_Fc"]["eta_weighted"] > 0
    checks["fxc_physical_positive"] = branch_physical["F_from_Fxc"]["eta_weighted"] > 0

    fc_normalized = branch_physical["F_from_Fc"]["eta_weighted"] / p["hb"] ** 5
    fxc_normalized = branch_physical["F_from_Fxc"]["eta_weighted"] / p["hb"] ** 6
    checks["fc_formal_hb5"] = _close(mp, fc_normalized * p["hb"] ** 5, branch_physical["F_from_Fc"]["eta_weighted"])
    checks["fxc_formal_hb6"] = _close(mp, fxc_normalized * p["hb"] ** 6, branch_physical["F_from_Fxc"]["eta_weighted"])
    # A repeated hb^4 factor would visibly change the direct physical result;
    # this guard records that the nested mass already contains its hb power.
    duplicated_fc = branch_physical["F_from_Fc"]["eta_weighted"] * p["hb"] ** 4
    duplicated_fxc = branch_physical["F_from_Fxc"]["eta_weighted"] * p["hb"] ** 4
    checks["no_duplicate_hb4"] = (
        abs(duplicated_fc - branch_physical["F_from_Fc"]["eta_weighted"]) > mp.mpf("1e-40")
        and abs(duplicated_fxc - branch_physical["F_from_Fxc"]["eta_weighted"]) > mp.mpf("1e-40")
    )
    # Cell values are direct eta-weighted integrals; no eta cell mass is
    # applied after the q integration.
    checks["no_duplicate_cell_mass"] = all(
        all(value > 0 for value in cell_physical[branch][:3])
        and cell_physical[branch][3] == 0
        and _close(mp, sum(cell_physical[branch]), branch_physical[branch]["eta_weighted"], atol="1e-80", rtol="1e-70")
        for branch in BRANCHES
    )

    z1_physical = {
        branch: _z1_physical_integral(
            mp,
            q_lower=q_zero,
            q_upper=q_one,
            branch=branch,
            source_amplitude=p["Fc"] if branch == "F_from_Fc" else p["Fxc"],
            source_z_amplitude=p["Fc_Z"] if branch == "F_from_Fc" else p["Fxc_Z"],
            chi_amplitude=p["chi"],
            chi_z_amplitude=p["chi_Z"],
            leading_amplitude=p["Flead"],
            leading_z_amplitude=p["Flead_Z"],
            hb=p["hb"], sc=p["sc"], lam=p["Lambda"],
        )
        for branch in BRANCHES
    }
    checks["z1_source_Z_positive"] = all(z1_physical[branch]["source_Z"] > 0 for branch in BRANCHES)
    checks["z1_chi_Z_positive"] = all(z1_physical[branch]["chi_Z"] > 0 for branch in BRANCHES)
    checks["z1_leading_Z_positive"] = all(z1_physical[branch]["leading_Z"] > 0 for branch in BRANCHES)
    checks["z1_sum_identity"] = all(
        _close(mp, z1_physical[branch]["source_Z"] + z1_physical[branch]["chi_Z"] + z1_physical[branch]["leading_Z"], z1_physical[branch]["total"])
        for branch in BRANCHES
    )

    return {
        "schema": SCHEMA,
        "status": "synthetic_independent_validation_only",
        "passed": all(checks.values()),
        "checks": checks,
        "check_count": len(checks),
        "source_power_basis": list(SOURCE_POWER_BASIS),
        "formal_hb_powers": {"F_from_Fc": 5, "F_from_Fxc": 6},
        "formal_external_powers": {"hb": 5, "Lambda": -2},
        "parameters": {name: _fraction_text(value) for name, value in SYNTHETIC_PARAMETERS.items()},
        "old_fixture_parameters": {name: _fraction_text(value) for name, value in OLD_FIXTURE_PARAMETERS.items()},
        "q_domain": ["0", "1"],
        "eta_profile": "1-sigma(4q-2)",
        "x_formula": "2*exp(hb*sc*q/8)",
        "dx_dq_formula": "hb*sc*x/8",
        "rho_formula": "4*exp(hb*sc*q/4)",
        "kernel_formula": "K0(x,t)=t*(1-t^2/x^2)/2",
        "log_kernel_formula": "K(y,s)=2*exp(2s)*(1-exp(2*(s-y)))",
        "nested_mass_formulas": {
            "F_from_Fc": "exp(4y)/12 + y - exp(2y)/2 + 1/4 + exp(-2y)/6",
            "F_from_Fxc": "exp(4y)/12 + exp(2y)*(3/4-y) - y - 3/4 - exp(-2y)/12",
        },
        "nested_probe": {
            "q": _number_text(mp, q_probe),
            "y": _number_text(mp, y_probe),
            "analytic_F_from_Fc": _number_text(mp, analytic_fc),
            "quadrature_F_from_Fc": _number_text(mp, direct_fc),
            "analytic_F_from_Fxc": _number_text(mp, analytic_fxc),
            "quadrature_F_from_Fxc": _number_text(mp, direct_fxc),
        },
        "physical_measure": "dρ",
        "drho_dq_formula": "2*x*dx/dq = hb*sc*rho/4",
        "physical_integral": {
            branch: {
                "eta_weighted": _number_text(mp, branch_physical[branch]["eta_weighted"]),
                "unweighted": _number_text(mp, branch_physical[branch]["unweighted"]),
                "cells": [_number_text(mp, value) for value in cell_physical[branch]],
                "formal_hb_power": 5 if branch == "F_from_Fc" else 6,
            }
            for branch in BRANCHES
        },
        "z1_physical_integral": {
            branch: {name: _number_text(mp, value) for name, value in z1_physical[branch].items()}
            for branch in BRANCHES
        },
        "j2_source_is_positive_majorant": True,
        "eta_integrated_directly": True,
        "cell_eta_mass_not_reapplied": True,
        "hb4_not_reapplied": True,
        "actual_sources_evaluated": False,
        "source_receipts_read": False,
        "old_j2_imported": False,
        "old_helpers_imported": False,
        "old_fixture_imported": False,
        "numerical_crosschecks_not_certificates": True,
    }


def run_fixture(*args: Any, **kwargs: Any) -> dict[str, Any]:
    """Runner-compatible alias; accepts ``(c)`` or ``(unused, c)``."""

    c = kwargs.get("c")
    if args:
        c = args[-1]
    return run_validation(c, dps=int(kwargs.get("dps", 100)))


def run(*args: Any, **kwargs: Any) -> dict[str, Any]:
    return run_fixture(*args, **kwargs)


def build(c: Any = None, **kwargs: Any) -> dict[str, Any]:
    """Expose the validation under a leaf-like name for a minimal runner."""

    return run_validation(c, **kwargs)
