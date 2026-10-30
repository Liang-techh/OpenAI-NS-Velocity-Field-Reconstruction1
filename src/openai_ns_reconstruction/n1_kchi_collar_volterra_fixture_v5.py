"""Self-contained focused fixture for :mod:`n1_kchi_collar_volterra_v5`.

The checks use exact rational physical-x data and directed containment.  They
do not import an ancestor fixture or receipt.  They cover fixed endpoint
validation, row-specific state completeness, the chi=0 analytic source/Z1
transport, and one nonzero-chi Z1/factorial-tail case.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Any


FIELDS = ("F", "U", "K", "P", "Fx", "Ux")


def _q(c: Any, value: Fraction | int) -> Any:
    if isinstance(value, int):
        return c.mpf(value)
    return c.mpf(value.numerator) / value.denominator


def _box(c: Any, lower: Fraction, upper: Fraction) -> Any:
    lo, hi = _q(c, lower), _q(c, upper)
    return c.make_mpf((lo._mpi_[0], hi._mpi_[1]))


def _cell(c: Any, left: Fraction, right: Fraction) -> tuple[Any, Any]:
    return _q(c, left), _q(c, right)


def _cells(c: Any) -> list[tuple[Any, Any]]:
    return [_cell(c, Fraction(2), Fraction(3))]


def _contains(box: Any, value: Any) -> bool:
    return box.a <= value.a and box.b >= value.b


def _mass(
    c: Any,
    seam: Fraction,
    target: Fraction,
    derivative: bool,
) -> Any:
    x = _q(c, target)
    a = _q(c, seam)
    b = x
    if derivative:
        return (b**4 - a**4) / (4 * x**3)
    return (b**2 - a**2) / 4 - (b**4 - a**4) / (8 * x**2)


def _incoming(c: Any) -> tuple[dict[str, Any], dict[str, Any]]:
    values = {
        "F": _q(c, 3),
        "U": _q(c, 2),
        "K": _q(c, -1),
        "P": _q(c, 5),
        "Fx": _q(c, 2),
        "Ux": _q(c, -3),
    }
    derivatives = {
        "F": _q(c, Fraction(1, 2)),
        "U": _q(c, Fraction(1, 3)),
        "K": _q(c, Fraction(-1, 4)),
        "P": _q(c, Fraction(2, 5)),
        "Fx": _q(c, Fraction(3, 4)),
        "Ux": _q(c, Fraction(-2, 3)),
    }
    return values, derivatives


def _base_kwargs(c: Any) -> dict[str, Any]:
    values, derivatives = _incoming(c)
    return {
        "x_cells": _cells(c),
        "incoming": values,
        # Canonical F/U/K/P/Fx/Ux names are valid because this mapping is
        # explicitly the Z1 row.
        "incoming_Z": derivatives,
        "chi": [_q(c, 0)],
        "chi_Z": [_q(c, 0)],
        "h5": [_q(c, 2)],
        "h5_Z": [_q(c, Fraction(1, 3))],
    }


def _zero_z(c: Any) -> dict[str, Any]:
    return {field: _q(c, 0) for field in FIELDS}


def run_fixture(build: Any, c: Any) -> dict[str, Any]:
    """Run focused independent checks for the v5 contract."""

    if hasattr(build, "build"):
        build = build.build

    checks: dict[str, bool] = {}
    kw = _base_kwargs(c)

    # Endpoint boxes whose second-cell dx crosses zero must be rejected by
    # cell validation, before source widths or factorial rates are formed.
    signed_width = dict(kw)
    signed_width["x_cells"] = [
        (_q(c, Fraction(2)), _box(c, Fraction(12, 5), Fraction(5, 2))),
        (
            _box(c, Fraction(5, 2), Fraction(13, 5)),
            _box(c, Fraction(251, 100), Fraction(51, 20)),
        ),
    ]
    signed_width["chi"] = [_q(c, 0), _q(c, 0)]
    signed_width["chi_Z"] = [_q(c, 0), _q(c, 0)]
    signed_width["h5"] = [_q(c, 0), _q(c, 0)]
    signed_width["h5_Z"] = [_q(c, 0), _q(c, 0)]
    try:
        build(c, **signed_width, max_order=0)
    except ValueError:
        checks["signed_width_interval_endpoints_rejected"] = True
    else:
        checks["signed_width_interval_endpoints_rejected"] = False

    values, derivatives = _incoming(c)

    # The ordinary row has five valid Z0 components but its F slot is
    # replaced by a Z1-only alias.  This must not leak across rows.
    ordinary_misplaced = dict(kw)
    ordinary_misplaced["incoming"] = {
        "F_Z": _q(c, 1),
        "U": values["U"],
        "K": values["K"],
        "P": values["P"],
        "Fx": values["Fx"],
        "Ux": values["Ux"],
    }
    try:
        build(c, **ordinary_misplaced, max_order=0)
    except ValueError:
        checks["ordinary_F_Z_alias_with_other_five_rejected"] = True
    else:
        checks["ordinary_F_Z_alias_with_other_five_rejected"] = False

    # The explicit Z1 row has five valid canonical components but omits F.
    # Canonical names are allowed in a complete incoming_Z row, not as a
    # fallback that fills a missing component.
    partial_z1 = dict(kw)
    partial_z1["incoming_Z"] = {
        "U": derivatives["U"],
        "K": derivatives["K"],
        "P": derivatives["P"],
        "Fx": derivatives["Fx"],
        "Ux": derivatives["Ux"],
    }
    try:
        build(c, **partial_z1, max_order=0)
    except ValueError:
        checks["incoming_Z_missing_F_with_other_five_rejected"] = True
    else:
        checks["incoming_Z_missing_F_with_other_five_rejected"] = False

    # chi=0 gives the exact one-cell Green/source transport.  This also
    # proves that a complete canonical incoming_Z row is accepted.
    zero = build(c, **kw, max_order=0)
    seam = Fraction(2)
    endpoint = Fraction(3)
    h = values["F"] + values["Fx"] * _q(c, seam) * (
        1 - (_q(c, seam) / _q(c, endpoint)) ** 2
    ) / 2
    hz = derivatives["F"] + derivatives["Fx"] * _q(c, seam) * (
        1 - (_q(c, seam) / _q(c, endpoint)) ** 2
    ) / 2
    source = _mass(c, seam, endpoint, False) * _q(c, 2)
    source_z = _mass(c, seam, endpoint, False) * _q(c, Fraction(1, 3))
    source_x = _mass(c, seam, endpoint, True) * _q(c, 2)
    source_zx = _mass(c, seam, endpoint, True) * _q(c, Fraction(1, 3))
    checks["chi_zero_source_F"] = _contains(zero["state"]["F"], h + source)
    checks["chi_zero_source_F_Z"] = _contains(zero["state_Z"]["F"], hz + source_z)
    checks["chi_zero_source_Fx"] = _contains(
        zero["state"]["Fx"],
        values["Fx"] * (_q(c, seam) / _q(c, endpoint)) ** 3 + source_x,
    )
    checks["chi_zero_source_Fx_Z"] = _contains(
        zero["state_Z"]["Fx"],
        derivatives["Fx"] * (_q(c, seam) / _q(c, endpoint)) ** 3 + source_zx,
    )
    checks["complete_canonical_incoming_Z_accepted"] = (
        zero["incoming_Z"]["F"].a == derivatives["F"].a
        and zero["incoming_Z"]["Ux"].b == derivatives["Ux"].b
    )

    # Nonzero chi and chi_Z: at max_order=1 the first Z1 word is finite and
    # must contain the chi_Z*F insertion.  The omitted words start at order 2
    # and retain a separate factorial tail.
    nonzero = dict(kw)
    nonzero["incoming"] = {
        "F": _q(c, 1),
        "U": values["U"],
        "K": values["K"],
        "P": values["P"],
        "Fx": _q(c, 0),
        "Ux": values["Ux"],
    }
    nonzero["incoming_Z"] = _zero_z(c)
    nonzero["chi"] = [_q(c, Fraction(1, 7))]
    nonzero["chi_Z"] = [_q(c, Fraction(1, 20))]
    nonzero["h5"] = [_q(c, 0)]
    nonzero["h5_Z"] = [_q(c, 0)]
    first = build(c, **nonzero, max_order=1)
    expected_z_first = -2 * nonzero["chi_Z"][0] * _mass(
        c, seam, endpoint, False
    )
    checks["nonzero_chi_Z_product_rule"] = _contains(
        first["homogeneous"]["finite_prefix"]["F_Z"], expected_z_first
    )
    facts = first["homogeneous"]["tail"]
    checks["factorial_tail_starts_after_word_one"] = (
        facts["factorial"]["word_tail_starts_at"] == 2
        and facts["symmetric_error"] is True
        and facts["error_only"] is True
    )
    checks["factorial_tail_has_Z1_radius"] = (
        facts["radius"]["F"].b > 0 and facts["radius"]["F_Z"].b > 0
    )

    return {
        "passed": all(checks.values()),
        "checks": checks,
        "check_count": len(checks),
        "synthetic_only": True,
        "ancestor_fixture_imported": False,
        "receipt_read": False,
        "midpoint_values_used": False,
        "analytic_cases": [
            "signed-width interval endpoints rejected",
            "ordinary F_Z alias rejected with other five Z0 fields",
            "incoming_Z missing F rejected with other five canonical fields",
            "chi=0 source and canonical Z1 transport",
            "nonzero chi and chi_Z first Z1 word plus factorial tail",
        ],
    }


run = run_fixture

__all__ = ["run_fixture", "run"]
