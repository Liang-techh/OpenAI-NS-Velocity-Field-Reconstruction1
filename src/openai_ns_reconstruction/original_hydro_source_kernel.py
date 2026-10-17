"""Pure ordinary-Z hydro source arithmetic.

This module contains only local arithmetic for the hydro drive used by the
original source reconstruction.  It deliberately does not load source data,
choose scales, integrate, or publish rows.  A jet is an ordinary derivative
row, so its entries are ``[f, f_Z, f_ZZ, ...]`` rather than Taylor
coefficients.

The fixed-delta entry point returns the ordinary ``Z0..Z2`` hydro drive.  The
sector entry point keeps a formal epsilon polynomial as dictionaries mapping
degrees to ordinary jets.  It returns the numerator and common denominator
separately; no epsilon series division or degree truncation is performed.
All scalar operations are delegated to the supplied mpmath point or
directed-interval values.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

from .original_source_jet_kernel import ordinary_product, ordinary_quotient


def _ordinary_jet(value: Any, order: int, name: str) -> list[Any]:
    """Validate and copy one exact ordinary derivative row."""
    if value is None:
        raise ValueError(f"{name} is required; missing ordinary derivatives")
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise ValueError(f"{name} must be an ordinary derivative sequence")
    if len(value) != order + 1:
        raise ValueError(
            f"{name} must contain ordinary Z0..Z{order}; no zero padding is allowed"
        )
    if any(item is None for item in value):
        raise ValueError(f"{name} contains a missing ordinary derivative")
    return list(value)


def _jet_add(*rows: Sequence[Any]) -> list[Any]:
    if not rows:
        raise ValueError("at least one ordinary jet is required")
    length = len(rows[0])
    if any(len(row) != length for row in rows):
        raise ValueError("ordinary jet orders must match")
    result = list(rows[0])
    for row in rows[1:]:
        result = [left + right for left, right in zip(result, row)]
    return result


def _jet_scale(row: Sequence[Any], scalar: Any) -> list[Any]:
    return [scalar * value for value in row]


def _jet_sub(left: Sequence[Any], right: Sequence[Any]) -> list[Any]:
    return _jet_add(left, _jet_scale(right, -1))


def _jet_divide_integer(row: Sequence[Any], divisor: int) -> list[Any]:
    return [value / divisor for value in row]


def _z_jet(z: Any) -> list[Any]:
    return [z, 1, 0]


def _d_jet(z: Any) -> list[Any]:
    """The ordinary jet of d=1-z^2 through Z2."""
    return [1 - z * z, -2 * z, -2]


def _l_jet(z: Any, delta: Any) -> list[Any]:
    """The ordinary jet of L=1-delta*z^2 through Z2."""
    return [1 - delta * z * z, -2 * delta * z, -2 * delta]


def _shift_derivative(rows: Sequence[Any], name: str) -> list[Any]:
    """Return the first ordinary derivative row, retaining the supplied Z3."""
    checked = _ordinary_jet(rows, 3, name)
    return checked[1:]


def hydro_drive_ordinary(
    z: Any,
    delta: Any,
    v: Sequence[Any],
    m: Sequence[Any],
    a: Sequence[Any],
) -> list[Any]:
    """Compute the hydro drive ordinary derivatives ``Z0..Z2``.

    Parameters
    ----------
    z, delta:
        Scalar mpmath point values or directed interval values.  ``delta`` is
        held fixed while differentiating with respect to ``z``.
    v:
        Exactly three ordinary derivatives ``[v, v_Z, v_ZZ]``.
    m, a:
        Exactly four ordinary derivatives through ``Z3``.  The third entry is
        required because the formula contains ``m_Z`` and ``a_Z`` and the
        result is differentiated through ``Z2``.

    Returns
    -------
    list
        ``[D_hydro, (D_hydro)_Z, (D_hydro)_ZZ]`` as ordinary derivatives.

    The calculation follows

    ``d = 1-z^2``, ``L = 1-delta*z^2``,
    ``W = 1-(1-delta)*z*m-d*m_Z`` and
    ``D_hydro = (-W*v + (1-delta)*(m-z*m_Z)/2
                 + 2*delta*z*a - d*a_Z)/(2*L)``.
    """
    v_row = _ordinary_jet(v, 2, "v")
    m_row = _ordinary_jet(m, 3, "m")
    a_row = _ordinary_jet(a, 3, "a")
    if delta is None:
        raise ValueError("delta is required")

    z_row = _z_jet(z)
    d_row = _d_jet(z)
    m_z = _shift_derivative(m_row, "m")
    a_z = _shift_derivative(a_row, "a")

    z_m = ordinary_product(z_row, m_row[:3])
    z_m_z = ordinary_product(z_row, m_z)
    d_m_z = ordinary_product(d_row, m_z)
    w_row = _jet_sub(
        [1, 0, 0],
        _jet_add(
            _jet_scale(z_m, 1 - delta),
            d_m_z,
        ),
    )

    w_v = ordinary_product(w_row, v_row)
    m_minus_z_m_z = _jet_sub(m_row[:3], z_m_z)
    hydro_numerator = _jet_add(
        _jet_scale(w_v, -1),
        _jet_divide_integer(_jet_scale(m_minus_z_m_z, 1 - delta), 2),
        _jet_scale(ordinary_product(z_row, a_row[:3]), 2 * delta),
        _jet_scale(ordinary_product(d_row, a_z), -1),
    )
    denominator = _jet_scale(_l_jet(z, delta), 2)
    return ordinary_quotient(hydro_numerator, denominator)


def hydro_moment_difference_numerator_ordinary(
    z: Any,
    delta: Any,
    v: Sequence[Any],
    delta_m: Sequence[Any],
    delta_a: Sequence[Any],
) -> list[Any]:
    """Compute the uncancelled linear moment-difference hydro numerator.

    ``v`` is required through ``Z2`` and both ``delta_m`` and ``delta_a``
    through ``Z3``.  The returned row is the ordinary ``Z0..Z2`` numerator

    ``(1-delta)*z*delta_m*v + d*delta_m_Z*v
     + (1-delta)*(delta_m-z*delta_m_Z)/2
     + 2*delta*z*delta_a - d*delta_a_Z``.

    This is a numerator-only identity.  Callers that need the drive divide
    it by the shared ``2*L`` with :func:`ordinary_quotient` after preserving
    any source-proved cancellation.
    """
    v_row = _ordinary_jet(v, 2, "v")
    m_difference = _ordinary_jet(delta_m, 3, "delta_m")
    a_difference = _ordinary_jet(delta_a, 3, "delta_a")
    if delta is None:
        raise ValueError("delta is required")

    z_row = _z_jet(z)
    d_row = _d_jet(z)
    m_z = _shift_derivative(m_difference, "delta_m")
    a_z = _shift_derivative(a_difference, "delta_a")
    z_m = ordinary_product(z_row, m_difference[:3])
    z_m_z = ordinary_product(z_row, m_z)
    return _jet_add(
        _jet_scale(ordinary_product(z_m, v_row), 1 - delta),
        ordinary_product(ordinary_product(d_row, m_z), v_row),
        _jet_divide_integer(_jet_scale(_jet_sub(m_difference[:3], z_m_z), 1 - delta), 2),
        _jet_scale(ordinary_product(z_row, a_difference[:3]), 2 * delta),
        _jet_scale(ordinary_product(d_row, a_z), -1),
    )


def _sector_map(
    value: Any,
    order: int,
    name: str,
) -> dict[int, list[Any]]:
    if value is None:
        raise ValueError(f"{name} sectors are required")
    if not isinstance(value, Mapping) or not value:
        raise ValueError(f"{name} sectors must be a nonempty degree mapping")
    result: dict[int, list[Any]] = {}
    for degree, row in value.items():
        if isinstance(degree, bool) or not isinstance(degree, int) or degree < 0:
            raise ValueError(f"{name} sector degrees must be nonnegative integers")
        result[degree] = _ordinary_jet(row, order, f"{name}[{degree}]")
    return result


def _delta_sector_map(
    delta: Any,
    delta_sectors: Mapping[int, Any] | None,
) -> dict[int, list[Any]]:
    supplied = {0: delta} if delta_sectors is None else delta_sectors
    if not isinstance(supplied, Mapping) or not supplied:
        raise ValueError("delta_sectors must be a nonempty degree mapping")
    result: dict[int, list[Any]] = {}
    for degree, scalar in supplied.items():
        if isinstance(degree, bool) or not isinstance(degree, int) or degree < 0:
            raise ValueError("delta sector degrees must be nonnegative integers")
        if scalar is None:
            raise ValueError(f"delta[{degree}] is missing")
        # Delta is a separately supplied fixed scalar.  It has no ordinary-Z
        # derivatives; represent its scalar sector as a constant jet.
        result[degree] = [scalar, 0, 0]
    return result


def _poly_add(*polynomials: Mapping[int, Sequence[Any]]) -> dict[int, list[Any]]:
    if not polynomials:
        raise ValueError("at least one sector polynomial is required")
    result: dict[int, list[Any]] = {}
    for polynomial in polynomials:
        for degree, row in polynomial.items():
            if degree in result:
                result[degree] = _jet_add(result[degree], row)
            else:
                result[degree] = list(row)
    return result


def _poly_scale(polynomial: Mapping[int, Sequence[Any]], scalar: Any) -> dict[int, list[Any]]:
    return {degree: _jet_scale(row, scalar) for degree, row in polynomial.items()}


def _poly_product(
    left: Mapping[int, Sequence[Any]],
    right: Mapping[int, Sequence[Any]],
) -> dict[int, list[Any]]:
    result: dict[int, list[Any]] = {}
    for left_degree, left_row in left.items():
        for right_degree, right_row in right.items():
            degree = left_degree + right_degree
            product = ordinary_product(left_row, right_row)
            if degree in result:
                result[degree] = _jet_add(result[degree], product)
            else:
                result[degree] = product
    return result


def _poly_sub(
    left: Mapping[int, Sequence[Any]],
    right: Mapping[int, Sequence[Any]],
) -> dict[int, list[Any]]:
    return _poly_add(left, _poly_scale(right, -1))


def _poly_constant(value: Any) -> dict[int, list[Any]]:
    return {0: [value, 0, 0]}


def _poly_trim(rows: Mapping[int, Sequence[Any]], name: str) -> dict[int, list[Any]]:
    result: dict[int, list[Any]] = {}
    for degree, row in rows.items():
        if len(row) == 4:
            result[degree] = _ordinary_jet(row, 3, f"{name}[{degree}]")[:3]
        else:
            result[degree] = _ordinary_jet(row, 2, f"{name}[{degree}]")
    return result


def _poly_shift_derivative(
    rows: Mapping[int, Sequence[Any]],
    name: str,
) -> dict[int, list[Any]]:
    return {
        degree: _shift_derivative(row, f"{name}[{degree}]")
        for degree, row in rows.items()
    }


def hydro_drive_ordinary_epsilon_sectors(
    z: Any,
    delta: Any,
    v_sectors: Mapping[int, Sequence[Any]],
    m_sectors: Mapping[int, Sequence[Any]],
    a_sectors: Mapping[int, Sequence[Any]],
    *,
    delta_sectors: Mapping[int, Any] | None = None,
) -> dict[str, dict[int, list[Any]]]:
    """Return the untruncated epsilon sectors of hydro numerator and 2L.

    ``v_sectors`` maps epsilon degree to ordinary ``Z0..Z2`` rows, while
    ``m_sectors`` and ``a_sectors`` map degree to rows through ``Z3``.  The
    positional ``delta`` is used as the fixed scalar degree-zero sector when
    ``delta_sectors`` is omitted.  To encode a factored identity such as
    ``delta = epsilon * exp(970)``, pass ``delta_sectors={1: exp970}``;
    the positional scalar may then be ``None``.

    The result contains exactly two mappings:

    ``numerator_sectors``
        all degrees of the hydro numerator, including zero-valued sectors
        created by cancellation;
    ``denominator_sectors``
        all degrees of ``2*(1-delta*z^2)``.

    No sector is dropped because it is small or zero, and no formal division
    by the denominator is attempted.
    """
    v_rows = _sector_map(v_sectors, 2, "v")
    m_rows = _sector_map(m_sectors, 3, "m")
    a_rows = _sector_map(a_sectors, 3, "a")
    delta_rows = _delta_sector_map(delta, delta_sectors)

    one = _poly_constant(1)
    z_row = _z_jet(z)
    z_poly = {0: z_row}
    z_squared = _poly_product(z_poly, z_poly)
    d_poly = {0: _d_jet(z)}

    v_poly = _poly_trim(v_rows, "v")
    m_poly = _poly_trim(m_rows, "m")
    a_poly = _poly_trim(a_rows, "a")
    m_z_poly = _poly_shift_derivative(m_rows, "m")
    a_z_poly = _poly_shift_derivative(a_rows, "a")

    one_minus_delta = _poly_sub(one, delta_rows)
    z_m = _poly_product(z_poly, m_poly)
    z_m_z = _poly_product(z_poly, m_z_poly)
    d_m_z = _poly_product(d_poly, m_z_poly)
    w_poly = _poly_sub(
        one,
        _poly_add(
            _poly_product(one_minus_delta, z_m),
            d_m_z,
        ),
    )

    m_minus_z_m_z = _poly_sub(m_poly, z_m_z)
    # Divide only the second term by 2; the common denominator is kept
    # separate below.  No epsilon series division is performed here.
    numerator = _poly_add(
        _poly_scale(_poly_product(w_poly, v_poly), -1),
        {
            degree: _jet_divide_integer(row, 2)
            for degree, row in _poly_product(one_minus_delta, m_minus_z_m_z).items()
        },
        _poly_scale(_poly_product(delta_rows, _poly_product(z_poly, a_poly)), 2),
        _poly_scale(_poly_product(d_poly, a_z_poly), -1),
    )
    denominator = _poly_scale(
        _poly_sub(one, _poly_product(delta_rows, z_squared)),
        2,
    )
    return {
        "numerator_sectors": numerator,
        "denominator_sectors": denominator,
    }


def hydro_moment_difference_ordinary_epsilon_sectors(
    z: Any,
    delta: Any,
    v_sectors: Mapping[int, Sequence[Any]],
    dm_sectors: Mapping[int, Sequence[Any]],
    da_sectors: Mapping[int, Sequence[Any]],
    *,
    delta_sectors: Mapping[int, Any] | None = None,
) -> dict[str, dict[int, list[Any]]]:
    """Return sector polynomials for a linear moment-difference numerator.

    ``v_sectors`` contains rows through ``Z2`` and ``dm_sectors``/
    ``da_sectors`` contain rows through ``Z3``.  The numerator is

    ``(1-delta)*z*dm*v + d*dm_Z*v
     + (1-delta)*(dm-z*dm_Z)/2
     + 2*delta*z*da - d*da_Z``.

    The common denominator ``2*L`` is returned in parallel.  Every epsilon
    degree generated by the supplied sector maps is retained, including
    sectors that cancel to zero.
    """
    v_rows = _sector_map(v_sectors, 2, "v")
    dm_rows = _sector_map(dm_sectors, 3, "dm")
    da_rows = _sector_map(da_sectors, 3, "da")
    delta_rows = _delta_sector_map(delta, delta_sectors)

    one = _poly_constant(1)
    z_poly = {0: _z_jet(z)}
    d_poly = {0: _d_jet(z)}
    z_squared = _poly_product(z_poly, z_poly)
    v_poly = _poly_trim(v_rows, "v")
    dm_poly = _poly_trim(dm_rows, "dm")
    da_poly = _poly_trim(da_rows, "da")
    dm_z_poly = _poly_shift_derivative(dm_rows, "dm")
    da_z_poly = _poly_shift_derivative(da_rows, "da")
    one_minus_delta = _poly_sub(one, delta_rows)

    z_dm = _poly_product(z_poly, dm_poly)
    z_dm_z = _poly_product(z_poly, dm_z_poly)
    dm_minus_z_dm_z = _poly_sub(dm_poly, z_dm_z)
    numerator = _poly_add(
        _poly_product(_poly_product(one_minus_delta, z_dm), v_poly),
        _poly_product(_poly_product(d_poly, dm_z_poly), v_poly),
        {
            degree: _jet_divide_integer(row, 2)
            for degree, row in _poly_product(one_minus_delta, dm_minus_z_dm_z).items()
        },
        _poly_scale(_poly_product(delta_rows, _poly_product(z_poly, da_poly)), 2),
        _poly_scale(_poly_product(d_poly, da_z_poly), -1),
    )
    denominator = _poly_scale(
        _poly_sub(one, _poly_product(delta_rows, z_squared)),
        2,
    )
    return {
        "numerator_sectors": numerator,
        "denominator_sectors": denominator,
    }


def evaluate_ordinary_sectors(
    sectors: Mapping[int, Sequence[Any]],
    epsilon: Any,
) -> list[Any]:
    """Evaluate a complete ordinary-jet sector polynomial at one epsilon."""
    if not isinstance(sectors, Mapping) or not sectors:
        raise ValueError("sector polynomial must be a nonempty degree mapping")
    first_degree = next(iter(sectors))
    first_row = sectors[first_degree]
    if isinstance(first_row, (str, bytes)) or not isinstance(first_row, Sequence):
        raise ValueError(f"sector[{first_degree}] must be an ordinary derivative sequence")
    if len(first_row) < 2:
        raise ValueError(f"sector[{first_degree}] must contain at least ordinary Z0 and Z1")
    order = len(first_row) - 1
    result = _ordinary_jet(first_row, order, f"sector[{first_degree}]")
    result = _jet_scale(result, epsilon ** first_degree)
    for degree, row in sectors.items():
        if degree == first_degree:
            continue
        checked = _ordinary_jet(row, order, f"sector[{degree}]")
        result = _jet_add(result, _jet_scale(checked, epsilon ** degree))
    return result


def hydro_drive_ordinary_epsilon(
    z: Any,
    delta: Any,
    v_sectors: Mapping[int, Sequence[Any]],
    m_sectors: Mapping[int, Sequence[Any]],
    a_sectors: Mapping[int, Sequence[Any]],
    epsilon: Any,
    *,
    delta_sectors: Mapping[int, Any] | None = None,
) -> list[Any]:
    """Evaluate sector numerator and denominator, then use the exact quotient."""
    sectors = hydro_drive_ordinary_epsilon_sectors(
        z,
        delta,
        v_sectors,
        m_sectors,
        a_sectors,
        delta_sectors=delta_sectors,
    )
    numerator = evaluate_ordinary_sectors(sectors["numerator_sectors"], epsilon)
    denominator = evaluate_ordinary_sectors(sectors["denominator_sectors"], epsilon)
    return ordinary_quotient(numerator, denominator)


# Explicit aliases keep the two public concepts discoverable without creating
# alternate arithmetic implementations.
hydro_drive_ordinary_polynomial = hydro_drive_ordinary_epsilon_sectors
hydro_drive_ordinary_sector_decomposition = hydro_drive_ordinary_epsilon_sectors


__all__ = [
    "evaluate_ordinary_sectors",
    "hydro_drive_ordinary",
    "hydro_drive_ordinary_epsilon",
    "hydro_drive_ordinary_epsilon_sectors",
    "hydro_moment_difference_numerator_ordinary",
    "hydro_moment_difference_ordinary_epsilon_sectors",
    "hydro_drive_ordinary_polynomial",
    "hydro_drive_ordinary_sector_decomposition",
]
