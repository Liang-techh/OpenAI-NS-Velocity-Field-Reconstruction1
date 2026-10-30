"""Independent analytic fixture for :mod:`n1_kchi_collar_volterra_v4`.

The fixture is deliberately producer-free: it creates rational interval
boxes and checks the consumer against independent physical-x formulas.  It
contains four cases: chi=0 with nonzero incoming/source/Z1 data, a constant
nonzero chi Picard prefix, a variable chi with nonzero chi_Z, and a local
h5_Z source.  The checks use directed interval containment and rational
arithmetic; they do not use a midpoint or a saved receipt.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Any


FIELDS = ("F", "U", "K", "P", "Fx", "Ux")


def _q(c: Any, value: Fraction | int) -> Any:
    if isinstance(value, int):
        return c.mpf(value)
    return c.mpf(value.numerator) / value.denominator


def _cell(c: Any, left: Fraction, right: Fraction) -> tuple[Any, Any]:
    return _q(c, left), _q(c, right)


def _cells(c: Any) -> list[tuple[Any, Any]]:
    return [_cell(c, Fraction(2) + Fraction(i, 4), Fraction(2) + Fraction(i + 1, 4))
            for i in range(4)]


def _contains(box: Any, value: Any) -> bool:
    value = value if hasattr(value, "a") else value
    return box.a <= value.a and box.b >= value.b


def _mass(c: Any, seam: Fraction, target: Fraction, derivative: bool) -> Any:
    x = _q(c, target)
    a = _q(c, seam)
    b = x
    if derivative:
        return (b ** 4 - a ** 4) / (4 * x ** 3)
    return (b ** 2 - a ** 2) / 4 - (b ** 4 - a ** 4) / (8 * x ** 2)


def _cell_mass(c: Any, left: Fraction, right: Fraction,
               target: Fraction, derivative: bool) -> Any:
    x = _q(c, target)
    a = _q(c, left)
    b = _q(c, right)
    if derivative:
        return (b ** 4 - a ** 4) / (4 * x ** 3)
    return (b ** 2 - a ** 2) / 4 - (b ** 4 - a ** 4) / (8 * x ** 2)


def _incoming(c: Any) -> tuple[dict[str, Any], dict[str, Any]]:
    values = {"F": _q(c, 3), "U": _q(c, 2), "K": _q(c, -1),
              "P": _q(c, 5), "Fx": _q(c, 2), "Ux": _q(c, -3)}
    derivatives = {"F": _q(c, Fraction(1, 2)), "U": _q(c, Fraction(1, 3)),
                   "K": _q(c, Fraction(-1, 4)), "P": _q(c, Fraction(2, 5)),
                   "Fx": _q(c, Fraction(3, 4)), "Ux": _q(c, Fraction(-2, 3))}
    return values, derivatives


def _base_kwargs(c: Any) -> dict[str, Any]:
    values, derivatives = _incoming(c)
    cells = _cells(c)
    zeros = [_q(c, 0) for _ in cells]
    return {
        "x_cells": cells,
        "incoming": values,
        "incoming_Z": derivatives,
        "chi": zeros,
        "chi_Z": zeros,
        "h5": [_q(c, 2) for _ in cells],
        "h5_Z": [_q(c, Fraction(1, 3)) for _ in cells],
    }


def run_fixture(build: Any, c: Any) -> dict[str, Any]:
    """Run independent interval checks; no ancestor producer is called."""

    if hasattr(build, "build"):
        build = build.build
    checks: dict[str, bool] = {}

    # Case 1: chi=0 gives the exact K0 Green source and nonzero incoming jet.
    kw = _base_kwargs(c)
    zero = build(c, **kw, max_order=0)
    missing_z = dict(kw)
    missing_z.pop("incoming_Z")
    try:
        build(c, **missing_z, max_order=0)
    except TypeError:
        checks["missing_incoming_Z_rejected"] = True
    else:
        checks["missing_incoming_Z_rejected"] = False
    partial_direct = dict(kw)
    partial_direct.pop("incoming")
    partial_direct.update({"Fc": _q(c, 1), "Fxc": _q(c, 0)})
    try:
        build(c, **partial_direct, max_order=0)
    except ValueError:
        checks["partial_direct_state_rejected"] = True
    else:
        checks["partial_direct_state_rejected"] = False
    partial_direct_z = dict(kw)
    partial_direct_z.pop("incoming_Z")
    partial_direct_z.update({"Fc_Z": _q(c, 0), "Fxc_Z": _q(c, 0)})
    try:
        build(c, **partial_direct_z, max_order=0)
    except ValueError:
        checks["partial_direct_Z1_rejected"] = True
    else:
        checks["partial_direct_Z1_rejected"] = False
    endpoint = _q(c, Fraction(3))
    seam = Fraction(2)
    values, derivatives = _incoming(c)
    h = values["F"] + values["Fx"] * _q(c, seam) * (1 - (_q(c, seam) / endpoint) ** 2) / 2
    hz = derivatives["F"] + derivatives["Fx"] * _q(c, seam) * (1 - (_q(c, seam) / endpoint) ** 2) / 2
    source = _mass(c, seam, Fraction(3), False) * _q(c, 2)
    source_z = _mass(c, seam, Fraction(3), False) * _q(c, Fraction(1, 3))
    checks["chi_zero_exact_F"] = _contains(zero["state"]["F"], h + source)
    checks["chi_zero_exact_F_Z"] = _contains(zero["state_Z"]["F"], hz + source_z)
    checks["chi_zero_exact_Fx"] = _contains(
        zero["state"]["Fx"], values["Fx"] * (_q(c, seam) / endpoint) ** 3
        + _mass(c, seam, Fraction(3), True) * _q(c, 2)
    )
    checks["h5_Z_once"] = _contains(zero["forcing"]["finite_prefix"]["F_Z"], source_z)

    # Case 1b: chi is identically zero but chi_Z is not.  With N=0 the
    # omitted first word is an A_Z insertion and must remain in the tail.
    chi_z_zero = dict(kw)
    chi_z_zero["chi_Z"] = [_q(c, Fraction(1, 20)) for _ in kw["x_cells"]]
    chi_z_zero["h5"] = [_q(c, 0) for _ in kw["x_cells"]]
    chi_z_zero["h5_Z"] = [_q(c, 0) for _ in kw["x_cells"]]
    chi_z_zero["incoming"] = {"F": _q(c, 1), "U": values["U"], "K": values["K"],
                               "P": values["P"], "Fx": _q(c, 0), "Ux": values["Ux"]}
    chi_z_zero["incoming_Z"] = {key: _q(c, 0) for key in FIELDS}
    chi_z_n0 = build(c, **chi_z_zero, max_order=0)
    expected_chi_z_n0 = _q(c, 0)
    for i in range(4):
        expected_chi_z_n0 = expected_chi_z_n0 - 2 * chi_z_zero["chi_Z"][i] * _cell_mass(
            c, Fraction(2) + Fraction(i, 4), Fraction(2) + Fraction(i + 1, 4),
            Fraction(3), False
        )
    checks["chi_zero_nonzero_chi_Z_N0"] = _contains(
        chi_z_n0["state_Z"]["F"], expected_chi_z_n0
    )

    # Independent K0 transport check for the four states intentionally left
    # outside the Kchi pair.
    d = 2 * c.log(endpoint / _q(c, seam))
    e = c.exp(-d)
    checks["K0_U_transport"] = _contains(zero["state"]["U"], values["U"] + d * values["Ux"])
    checks["K0_K_transport"] = _contains(zero["state"]["K"], e * values["K"] - (1 - e) * values["Ux"])
    checks["K0_P_transport"] = _contains(zero["state"]["P"], values["P"])
    checks["K0_Ux_transport"] = _contains(zero["state"]["Ux"], c.exp(-d / 2) * values["Ux"])

    # A dimensional example with x around 10 makes the K0 value bound larger
    # than its x derivative bound; the tail metadata must retain both rather
    # than silently assuming d_x K0 dominates K0.
    large = dict(chi_z_zero)
    large_count = 32
    large["x_cells"] = [
        _cell(c, Fraction(2) + Fraction(8 * i, large_count),
              Fraction(2) + Fraction(8 * (i + 1), large_count))
        for i in range(large_count)
    ]
    large["chi"] = [_q(c, 0) for _ in range(large_count)]
    large["chi_Z"] = [_q(c, 0) for _ in range(large_count)]
    large["h5"] = [_q(c, 0) for _ in range(large_count)]
    large["h5_Z"] = [_q(c, 0) for _ in range(large_count)]
    large_result = build(c, **large, max_order=0)
    large_facts = large_result["homogeneous"]["tail"]["factorial"]
    checks["K0_sup_exceeds_dxK0_sup_example"] = (
        large_facts["K0_sup"].a > large_facts["d_x_K0_sup"].b
    )

    # Case 2: constant nonzero chi.  The independent first Picard prefix is
    # exact for a constant homogeneous inlet; the consumer's profile box
    # must contain it, while the factorial tail encloses later words.
    nonzero = dict(kw)
    nonzero["chi"] = [_q(c, Fraction(1, 7)) for _ in kw["x_cells"]]
    nonzero["h5"] = [_q(c, 0) for _ in kw["x_cells"]]
    nonzero["h5_Z"] = [_q(c, 0) for _ in kw["x_cells"]]
    nonzero["incoming"] = {"F": _q(c, 1), "U": values["U"], "K": values["K"],
                            "P": values["P"], "Fx": _q(c, 0), "Ux": values["Ux"]}
    nonzero["incoming_Z"] = {key: _q(c, 0) for key in FIELDS}
    first = build(c, **nonzero, max_order=1)
    first_picard = _q(c, 1) - 2 * _q(c, Fraction(1, 7)) * _mass(c, seam, Fraction(3), False)
    checks["constant_nonzero_chi_picard"] = _contains(
        first["homogeneous"]["finite_prefix"]["F"], first_picard
    )
    checks["factorial_tail_explicit"] = (
        first["homogeneous"]["tail"]["symmetric_error"] is True
        and first["homogeneous"]["tail"]["error_only"] is True
        and first["homogeneous"]["tail"]["factorial"]["word_tail_starts_at"] == 2
    )

    # Case 3: variable chi across cells and nonzero chi_Z.  The independent
    # first Z derivative of the homogeneous first word contains chi_Z*F;
    # this catches a silently frozen potential and a dropped Z1 product term.
    variable = dict(nonzero)
    variable["chi"] = [_q(c, Fraction(1 + i, 10)) for i in range(4)]
    variable["chi_Z"] = [_q(c, Fraction(1, 20)) for _ in range(4)]
    variable["incoming_Z"] = {key: _q(c, 0) for key in FIELDS}
    var = build(c, **variable, max_order=1)
    expected_z_first = _q(c, 0)
    for i in range(4):
        expected_z_first = expected_z_first - 2 * variable["chi_Z"][i] * _cell_mass(
            c, Fraction(2) + Fraction(i, 4), Fraction(2) + Fraction(i + 1, 4),
            Fraction(3), False
        )
    checks["variable_chi_cells_retained"] = (
        len(var["chi_cells"]) == 4
        and var["chi_cells"][0]._mpi_ != var["chi_cells"][3]._mpi_
    )
    checks["nonzero_chi_Z_product_rule"] = (
        var["chi_Z_included"] is True
        and var["ordinary_product_rule"]["Z1"] == ((1, 0, 1, "chi_Z * F"), (0, 1, 1, "chi * F_Z"))
        and _contains(var["homogeneous"]["finite_prefix"]["F_Z"], expected_z_first)
    )

    checks["incoming_nonzero_preserved"] = (
        var["incoming"]["F"].a > 0 and var["incoming"]["Fx"].a == 0
        and var["incoming_Z"]["F"].a == 0
    )
    checks["branches_separate"] = (
        var["homogeneous"] is not var["forcing"]
        and var["homogeneous"]["homogeneous_source_separate"] is True
        and var["forcing"]["homogeneous_source_separate"] is True
        and var["core_history_reset"] is False
    )
    checks["physical_dx_once"] = var["source_measure"] == "physical dx = x_right - x_left exactly once"
    checks["directed_no_midpoint"] = (
        var["directed_interval_arithmetic"] is True
        and var["midpoint_values_used"] is False
        and var["ordered_volterra_prefix"] is True
    )
    checks["six_state_contract"] = set(var["state"]) == set(FIELDS) and set(var["state_Z"]) == set(FIELDS)

    return {
        "passed": all(checks.values()),
        "checks": checks,
        "check_count": len(checks),
        "synthetic_only": True,
        "analytic_cases": ["chi=0", "chi=0 with nonzero chi_Z at N=0",
                            "constant nonzero chi first Picard", "variable chi and chi_Z", "nonzero h5_Z"],
        "independent_source_is_cellwise": True,
        "receipt_or_ancestor_read": False,
        "midpoint_values_used": False,
        "limits": {
            "constant_nonzero_case": "checks first polynomial Picard word; later words are covered by factorial tail",
            "chi_Z": "ordinary product rule is checked with nonzero chi_Z and zero incoming Z1 in the first word",
        },
    }


run = run_fixture

__all__ = ["run_fixture", "run"]
