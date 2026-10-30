"""Independent exact-rational fixture for the low-order core source leaf.

This fixture deliberately does not import the implementation under test.  It
builds the synthetic fields and the original leading operators with SymPy,
then checks a caller supplied ``build`` function through ordinary ``Z0``,
``Z1`` and ``Z2`` values.  The returned source terms keep their formal
``F0_power`` and ``Lambda_power`` labels; the check contracts those labels
with the fixed exact values only after the terms have been returned.

No midpoint, tolerance, persisted artifact, or model field is used here.
"""
from __future__ import annotations

from collections.abc import Mapping, Sequence
from fractions import Fraction


__all__ = ["run_fixture"]


def _fraction(bits):
    """Convert one exact mpmath interval endpoint to a rational."""

    sign, mantissa, exponent, _ = bits
    value = Fraction(int(mantissa), 1)
    value = value * (1 << exponent) if exponent >= 0 else value / (1 << (-exponent))
    return -value if sign else value


def _contains(interval, rational):
    """Test inclusion using the exact binary endpoints, never a midpoint."""

    expected = Fraction(int(rational.p), int(rational.q))
    return _fraction(interval._mpi_[0]) <= expected <= _fraction(interval._mpi_[1])


def _interval(c, value):
    """Make a directed interval from an exact SymPy rational."""

    value = value if hasattr(value, "p") and hasattr(value, "q") else value
    return c.mpf(int(value.p)) / int(value.q)


def _is_interval(value):
    return hasattr(value, "_mpi_") or (hasattr(value, "a") and hasattr(value, "b"))


def _term_parts(term, order):
    """Read one source-factor term without importing the source leaf."""

    if isinstance(term, Mapping):
        f0 = term.get("F0_power", term.get("f0_power"))
        lam = term.get("Lambda_power", term.get("lambda_power"))
        coefficient = term.get("coefficient", term.get("value"))
        if coefficient is None and "coefficients" in term:
            coefficient = term["coefficients"][order]
        if f0 is None or lam is None or coefficient is None:
            raise KeyError("a source term needs F0_power, Lambda_power and coefficient")
        return int(f0), int(lam), coefficient
    if isinstance(term, Sequence) and not isinstance(term, (str, bytes, bytearray)):
        if len(term) == 3:
            return int(term[0]), int(term[1]), term[2]
    raise TypeError("unsupported formal source term")


def _mapping_value(mapping, keys):
    if not isinstance(mapping, Mapping):
        return None
    for key in keys:
        try:
            if key in mapping:
                return mapping[key]
        except (TypeError, ValueError):
            continue
    return None


def _terms(value, order):
    """Normalize the one-cell term container used by a build result."""

    if isinstance(value, Mapping):
        if "coefficient" in value or "value" in value or "coefficients" in value:
            return [value]
        raise TypeError("a source row mapping is not a source term")
    if _is_interval(value):
        raise TypeError("formal powers are required on every source term")
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        values = list(value)
        if len(values) == 3 and not isinstance(values[0], (Mapping, list, tuple)):
            return [values]
        if values and all(isinstance(item, Mapping) or
                          (isinstance(item, Sequence) and not isinstance(item, (str, bytes, bytearray)))
                          for item in values):
            return values
    raise TypeError("unsupported source row")


def _row_terms(result, family, order, row):
    """Read the canonical order -> row -> terms form, with small shape guards."""

    data = result[family]
    order_value = _mapping_value(data, (str(order), order))
    if order_value is not None:
        if family == "g":
            value = _mapping_value(order_value, (str(row), row))
        else:
            value = _mapping_value(order_value, (row, str(row)))
        if value is not None:
            return _terms(value, order)

    # Accept the equivalent row -> order -> terms spelling, which is useful
    # for a 2-D array wrapper without importing an existing leaf.
    row_value = _mapping_value(data, (row, str(row)))
    if row_value is not None:
        value = _mapping_value(row_value, (str(order), order))
        if value is not None:
            return _terms(value, order)

    # A flat spelling such as g4_Z2 or C0[4,1]_Z2 is also unambiguous.
    flat = (
        f"{family}{row}_Z{order}",
        f"{family}[{row}]_Z{order}",
        f"{family}_{row}_Z{order}",
    )
    value = _mapping_value(result, flat)
    if value is not None:
        return _terms(value, order)
    raise KeyError(f"missing {family} row {row} ordinary Z^{order}")


def _recombine(c, terms, anchor, lam, order=0):
    total = c.mpf(0)
    sectors = []
    for term in terms:
        f0_power, lambda_power, coefficient = _term_parts(term, order)
        sectors.append((f0_power, lambda_power))
        total += coefficient * (lam ** lambda_power) * (anchor ** f0_power)
    return total, sectors


def _row_keys(result, family):
    """Return row keys from either order -> row or row -> order output."""

    data = result.get(family)
    if not isinstance(data, Mapping):
        return set()
    for order_key in ("0", 0):
        order_data = data.get(order_key)
        if isinstance(order_data, Mapping):
            return set(order_data)
    return set(data)


def _rational(value, sympy):
    value = sympy.cancel(sympy.simplify(value))
    if not value.is_Rational:
        raise AssertionError(f"fixture expected a rational value, got {value}")
    return value


def _grid(c, expression, r, z, sympy):
    return {
        f"rho{i}_Z{k}": _interval(c, _rational(
            sympy.diff(expression, r, i, z, k).subs({r: 1, z: sympy.Rational(1, 4)}),
            sympy,
        ))
        for i in range(3)
        for k in range(6 - i)
    }


def _axis_grid(c, expression, r, z, sympy):
    return {
        f"rho{i}_Z{k}": _interval(c, _rational(
            sympy.diff(expression, r, i, z, k).subs({r: 0, z: sympy.Rational(1, 4)}),
            sympy,
        ))
        for i in range(3)
        for k in range(6 - i)
    }


def _expected_fields(sympy):
    r, z = sympy.symbols("r z", real=True)
    z0 = sympy.Rational(1, 4)
    delta = sympy.Rational(1, 5)
    lam = sympy.Integer(3)
    j = sympy.Integer(2)
    anchor = sympy.Rational(7, 3)

    d = 1 - z * z
    L = 1 - delta * z * z
    H = z * (1 - delta) / 2 + d * (4 * z + j)
    ell = -lam * L * H / (H * H + (j / 500) ** 2)
    ell_d = [sympy.diff(ell, z, k).subs(z, z0) for k in range(5)]
    relative_derivatives = [sympy.Integer(1)]
    for n in range(5):
        relative_derivatives.append(sum(
            sympy.binomial(n, k) * ell_d[k] * relative_derivatives[n - k]
            for k in range(n + 1)
        ))
    # This is only a local polynomial through the order required by the
    # ordinary Z0/Z1/Z2 checks; it is not an exponential evaluation.
    F0 = anchor * sum(
        relative_derivatives[n] * (z - z0) ** n / sympy.factorial(n)
        for n in range(6)
    )

    M = 1 + r * z**5 / 2 + r**2 * z**3 / 3
    U = M + r * sympy.diff(M, r)
    Phi = 1 + r * z**5 + r**2 * z**3
    F = F0 * Phi

    Q = (2 * z * U - (1 - delta) * z * M - d * sympy.diff(M, z)) / L
    Qr = sympy.diff(Q, r)
    Qrr = sympy.diff(Q, r, 2)
    Qz = sympy.diff(Q, z)
    S = (Q + (1 - delta) * z * Qz / 2 + r * Qr) / L
    S += Q * (Q / 2 + r * Qr)
    S += U * (-2 * z * Q + d * Qz - 2 * z * r * Qr) / L
    Np = -S / 2 + lam * (2 * Qr + r * Qrr)

    def axial(field, weight):
        return (weight * z * field + d * sympy.diff(field, z)
                - 2 * z * r * sympy.diff(field, r)) / L

    Nt = -axial(axial(F, -2 - delta), -3)
    Nz = -axial(axial(U, -1 - delta), -2)
    g = {
        "4": 2 * sympy.sqrt(r) * Np / lam,
        "5": 2 * Nt / lam,
        "6": 2 * Nz / lam - 4 * z * r * Np / (lam**2 * L),
    }

    f = F
    fr = sympy.diff(F, r)
    w = U
    wr = sympy.diff(U, r)
    Hm = z * (1 - delta) / 2 + d * w
    B = Q + (1 - 2 * z * w) / L
    E = z * w - sympy.Rational(1, 2)
    Cf = f + r * fr
    Zf = (f * z * (-2 - delta) + d * sympy.diff(f, z)
          - fr * z * (2 * r)) / L
    Zw = (w * z * (-1 - delta) + d * sympy.diff(w, z)
          - wr * z * (2 * r)) / L

    C0 = {
        "4,1": 4 * sympy.sqrt(r) * f / lam,
        "5,1": 2 * (Q + E * (-2 + delta) / L) / lam,
        "5,2": (2 * Zf - Cf * z * (2 * (-1 + delta)) / L) / lam,
        "5,3": -Cf * z * (2 * (1 + delta)) / L / lam,
        "5,5": B * sympy.sqrt(r) / lam,
        "6,1": -f * z * (8 * r) / L / (lam**2),
        "6,2": 2 * (E * (-1 + delta) / L + Zw
                    - wr * z * ((-1 + delta) * r) / L) / lam,
        "6,3": -wr * z * (2 * (1 + delta) * r) / L / lam,
        "6,4": z * (-4) / L / lam,
        "6,6": B * sympy.sqrt(r) / lam,
    }
    B1 = {
        "5,1": 2 * Hm / L / lam,
        "5,2": -Cf * d * 2 / L / lam,
        "5,3": -Cf * d * 2 / L / lam,
        "6,2": 2 * (Hm - wr * d * r) / L / lam,
        "6,3": -wr * d * (2 * r) / L / lam,
        "6,4": d * 2 / L / lam,
    }
    return {
        "r": r,
        "z": z,
        "z0": z0,
        "delta": delta,
        "Lambda": lam,
        "j": j,
        "anchor": anchor,
        "Phi": Phi,
        "M": M,
        "U": U,
        "g": g,
        "C0": C0,
        "B1": B1,
    }


def run_fixture(build, c):
    """Run independent source and raw-matrix checks against ``build``.

    The caller must provide ``build(c,rho,z,delta,Lambda,j,gridPhi,gridU,
    gridM)``.  It is expected to return ``g``, ``C0`` and ``B1`` in ordinary
    ``Z0/Z1/Z2`` form, with each row represented by formal source terms.
    """

    import sympy as s

    expected = _expected_fields(s)
    r = expected["r"]
    z = expected["z"]
    z0 = expected["z0"]
    delta = expected["delta"]
    lam = expected["Lambda"]
    j = expected["j"]
    anchor = expected["anchor"]

    rho = _interval(c, s.Integer(1))
    z_value = _interval(c, z0)
    delta_value = _interval(c, delta)
    lambda_value = _interval(c, lam)
    j_value = _interval(c, j)
    result = build(
        c,
        rho,
        z_value,
        delta_value,
        lambda_value,
        j_value,
        _grid(c, expected["Phi"], r, z, s),
        _grid(c, expected["U"], r, z, s),
        _grid(c, expected["M"], r, z, s),
    )

    checks = {}
    anchor_value = _interval(c, anchor)
    row_sets = {
        "g": tuple(expected["g"]),
        "C0": tuple(expected["C0"]),
        "B1": tuple(expected["B1"]),
    }
    expected_sectors = {
        "g": {"4": ((0, 0), (0, -1)), "5": ((1, -1),),
              "6": ((0, -1), (0, -2))},
        "C0": {
            "4,1": ((1, -1),), "5,1": ((0, -1),),
            "5,2": ((1, -1),), "5,3": ((1, -1),),
            "5,5": ((0, -1),), "6,1": ((1, -2),),
            "6,2": ((0, -1),), "6,3": ((0, -1),),
            "6,4": ((0, -1),), "6,6": ((0, -1),),
        },
        "B1": {
            "5,1": ((0, -1),), "5,2": ((1, -1),),
            "5,3": ((1, -1),), "6,2": ((0, -1),),
            "6,3": ((0, -1),), "6,4": ((0, -1),),
        },
    }

    for family, rows in row_sets.items():
        for row in rows:
            expression = expected[family][row]
            for order in range(3):
                terms = _row_terms(result, family, order, row)
                total, sectors = _recombine(c, terms, anchor_value, lambda_value, order)
                exact = _rational(s.diff(expression, z, order).subs({r: 1, z: z0}), s)
                checks[f"{family}_{row}_Z{order}"] = _contains(total, exact)
                checks[f"{family}_{row}_Z{order}_formal_sectors"] = (
                    tuple(sorted(sectors)) == tuple(sorted(expected_sectors[family][row]))
                )

    checks["C0_has_exactly_ten_raw_entries"] = set(expected["C0"]) == _row_keys(result, "C0")
    checks["B1_has_exactly_six_raw_entries"] = set(expected["B1"]) == _row_keys(result, "B1")
    c0_61_terms = _row_terms(result, "C0", 0, "6,1")
    _, c0_61_sectors = _recombine(
        c, c0_61_terms, anchor_value, lambda_value, 0
    )
    checks["C0_6_1_is_F0_epsilon_squared"] = (
        tuple(sorted(c0_61_sectors)) == ((1, -2),)
    )
    checks["B1_6_3_keeps_original_L_denominator"] = _contains(
        _recombine(c, _row_terms(result, "B1", 0, "6,3"), anchor_value, lambda_value, 0)[0],
        _rational(expected["B1"]["6,3"].subs({r: 1, z: z0}), s),
    )

    # The source is axis-safe: the synthetic axis grid must not introduce a
    # rho division.  Check the actual g4 row at rho=0 without selecting a
    # value from any interval.
    axis_result = build(
        c,
        c.mpf(0),
        z_value,
        delta_value,
        lambda_value,
        j_value,
        _axis_grid(c, expected["Phi"], r, z, s),
        _axis_grid(c, expected["U"], r, z, s),
        _axis_grid(c, expected["M"], r, z, s),
    )
    axis_total, _ = _recombine(
        c, _row_terms(axis_result, "g", 0, "4"), anchor_value, lambda_value, 0
    )
    checks["g4_axis_is_exactly_zero"] = (
        _fraction(axis_total._mpi_[0]) == 0 and _fraction(axis_total._mpi_[1]) == 0
    )

    return {
        "passed": all(checks.values()),
        "checks": checks,
        "comparison": "independent SymPy ordinary Z0/Z1/Z2 source and raw C0/B1 formulas with exact interval endpoint inclusion",
        "fixture_is_not_actual_field": True,
    }
