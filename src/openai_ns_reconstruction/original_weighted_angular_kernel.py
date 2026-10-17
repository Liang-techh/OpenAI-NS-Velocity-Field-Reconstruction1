"""Pure Laurent--epsilon arithmetic for the original angular comparison.

The kernel works with ordinary ``Z`` derivatives, not Taylor coefficients.  A
sector polynomial is represented by a mapping ``epsilon_degree -> ordinary
Z0..Z2 row``.  The angular gauge term has degree ``-1`` because
``g = -Gprime / epsilon``; the correlated ``delta`` term has degree ``1``
because ``delta = epsilon * delta_scale``.  The two are kept as Laurent
sectors and
are never divided or formally truncated here.

The public entry point only builds the three numerator mode polynomials and
their common ``2*phi*L`` denominator.  It does not load source data, choose
scales, evaluate epsilon, integrate, or perform physical normalization.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from math import comb
from typing import Any


def _ordinary_row(value: Any, order: int, name: str) -> list[Any]:
    """Validate and copy a complete ordinary derivative row."""
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
    return [sum(values) for values in zip(*rows)]


def _jet_scale(row: Sequence[Any], scalar: Any) -> list[Any]:
    return [scalar * value for value in row]


def _jet_divide_integer(row: Sequence[Any], divisor: int) -> list[Any]:
    return [value / divisor for value in row]


def _jet_product(left: Sequence[Any], right: Sequence[Any]) -> list[Any]:
    """Multiply equal-length ordinary rows through Z 3 by Leibniz."""
    if len(left) != len(right):
        raise ValueError("ordinary jet orders must match")
    order = len(left) - 1
    if order not in (1, 2, 3):
        raise ValueError("ordinary rows must contain Z0..Z1, Z0..Z2, or Z0..Z3")
    return [
        sum(comb(n, k) * left[k] * right[n - k] for k in range(n + 1))
        for n in range(order + 1)
    ]


def _jet_shift(row: Sequence[Any]) -> list[Any]:
    if len(row) != 4:
        raise ValueError("a complete ordinary Z0..Z3 row is required for Z differentiation")
    return list(row[1:])


def _trim3(row: Sequence[Any]) -> list[Any]:
    if len(row) != 4:
        raise ValueError("a complete ordinary Z0..Z3 row is required")
    return list(row[:3])


def _z_row(z: Any) -> list[Any]:
    return [z, 1, 0]


def _d_row(z: Any) -> list[Any]:
    return [1 - z * z, -2 * z, -2]


def _poly_constant(value: Any, order: int = 2) -> dict[int, list[Any]]:
    return {0: [value] + [0] * order}


def _poly_add(*polynomials: Mapping[int, Sequence[Any]]) -> dict[int, list[Any]]:
    if not polynomials:
        raise ValueError("at least one epsilon sector polynomial is required")
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


def _poly_divide_integer(
    polynomial: Mapping[int, Sequence[Any]], divisor: int
) -> dict[int, list[Any]]:
    return {
        degree: _jet_divide_integer(row, divisor)
        for degree, row in polynomial.items()
    }


def _poly_product(
    left: Mapping[int, Sequence[Any]],
    right: Mapping[int, Sequence[Any]],
) -> dict[int, list[Any]]:
    result: dict[int, list[Any]] = {}
    for left_degree, left_row in left.items():
        for right_degree, right_row in right.items():
            degree = left_degree + right_degree
            product = _jet_product(left_row, right_row)
            if degree in result:
                result[degree] = _jet_add(result[degree], product)
            else:
                result[degree] = product
    return result


def _poly_trim3(polynomial: Mapping[int, Sequence[Any]]) -> dict[int, list[Any]]:
    return {degree: _trim3(row) for degree, row in polynomial.items()}


def _poly_shift(polynomial: Mapping[int, Sequence[Any]]) -> dict[int, list[Any]]:
    return {degree: _jet_shift(row) for degree, row in polynomial.items()}


def _poly_row(row: Sequence[Any], degree: int) -> dict[int, list[Any]]:
    return {degree: list(row)}


def _angular_polynomial(
    z: Any,
    delta: Mapping[int, Sequence[Any]],
    gauge: Mapping[int, Sequence[Any]],
    h4: Mapping[int, Sequence[Any]],
    k4: Mapping[int, Sequence[Any]],
) -> dict[int, list[Any]]:
    """Build angular sectors for complete H and K rows."""
    one = _poly_constant(1)
    z_poly = _poly_row(_z_row(z), 0)
    d_poly = _poly_row(_d_row(z), 0)
    h3 = _poly_trim3(h4)
    k3 = _poly_trim3(k4)
    h_z = _poly_shift(h4)
    k_z = _poly_shift(k4)

    one_minus_delta = _poly_add(one, _poly_scale(delta, -1))
    h_coefficient = _poly_add(one, _poly_scale(_poly_divide_integer(delta, 2), -1))
    two_delta_minus_one = _poly_add(_poly_scale(delta, 2), _poly_scale(one, -1))

    h_gauge = _poly_add(h_z, _poly_product(gauge, h3))
    k_gauge = _poly_add(k_z, _poly_product(gauge, k3))
    h_transport = _poly_product(_poly_product(one_minus_delta, z_poly), h_gauge)
    k_transport = _poly_product(d_poly, k_gauge)
    k_algebraic = _poly_product(_poly_product(two_delta_minus_one, z_poly), k3)

    return _poly_add(
        _poly_product(h_coefficient, h3),
        _poly_divide_integer(_poly_scale(h_transport, -1), 2),
        _poly_scale(k_transport, -1),
        k_algebraic,
    )


def _w_polynomial(
    z: Any,
    delta: Mapping[int, Sequence[Any]],
    m4: Mapping[int, Sequence[Any]],
) -> dict[int, list[Any]]:
    """Build W = 1 - (1-delta) z m - d m_Z through ordinary Z2."""
    one = _poly_constant(1)
    z_poly = _poly_row(_z_row(z), 0)
    d_poly = _poly_row(_d_row(z), 0)
    m3 = _poly_trim3(m4)
    m_z = _poly_shift(m4)
    one_minus_delta = _poly_add(one, _poly_scale(delta, -1))
    return _poly_add(
        one,
        _poly_scale(
            _poly_add(
                _poly_product(_poly_product(one_minus_delta, z_poly), m3),
                _poly_product(d_poly, m_z),
            ),
            -1,
        ),
    )


def angular_mode_numerator_epsilon_sectors(
    z: Any,
    u0: Sequence[Any],
    phi: Sequence[Any],
    H: Sequence[Any],
    Psi: Sequence[Any],
    Mpsi: Sequence[Any],
    Kpsi: Sequence[Any],
    Gprime: Sequence[Any],
    delta_scale: Any,
) -> dict[str, Any]:
    """Return the three pure Laurent--epsilon angular numerator modes.

    Every input except ``Gprime`` is a raw, undressed complete ordinary row
    through Z3; ``Gprime`` is the ordinary row through Z2 for ``G_Z``.  The
    common ``F0`` amplitude has already been canceled analytically, leaving
    its gauge derivative as ``g = -Gprime/epsilon``.  The returned
    ``numerator_modes`` list contains, in order, the mode 0 baseline, the
    direct m-difference mode 1, and the direct H/K-difference mode 2.  Each
    mode is a mapping from an integer Laurent degree (possibly ``-1``) to an
    ordinary Z0..Z2 row.  ``denominator_sectors`` is the normalized common
    denominator ``2*phi*L = 2*phi*(1-delta*z**2)`` for the returned
    numerators.

    The constructions are:

    * mode 0: ``m = u0 + epsilon*Psi``, ``H = phi``,
      ``K = phi*(u0 + epsilon*Psi)``;
    * mode 1: ``mDiff = epsilon*(Mpsi-Psi)``, with H/K zero;
    * mode 2: ``Hdiff = H-phi`` and
      ``Kdiff = u0*(H-phi) + epsilon*(Kpsi-phi*Psi)``, with W difference
      zero.

    Here ``g = -Gprime/epsilon`` and ``delta = epsilon*delta_scale`` are
    retained as Laurent sectors.  In the source normalization,
    ``delta_scale = exp(970)`` is fixed; it must not be replaced by the
    first-switch ``kappa(t) = 1-sigma(t)``.  No epsilon evaluation, quotient,
    scale selection, source loading, or formal truncation is performed.
    """
    u0_row = _ordinary_row(u0, 3, "u0")
    phi_row = _ordinary_row(phi, 3, "phi")
    h_row = _ordinary_row(H, 3, "H")
    psi_row = _ordinary_row(Psi, 3, "Psi")
    mpsi_row = _ordinary_row(Mpsi, 3, "Mpsi")
    kpsi_row = _ordinary_row(Kpsi, 3, "Kpsi")
    gprime_row = _ordinary_row(Gprime, 2, "Gprime")
    if delta_scale is None:
        raise ValueError("delta_scale is required")

    one3 = _poly_constant(1)
    z_poly = _poly_row(_z_row(z), 0)
    d_poly = _poly_row(_d_row(z), 0)
    delta = _poly_row([delta_scale, 0, 0], 1)
    gauge = _poly_row(_jet_scale(gprime_row, -1), -1)
    phi3 = _poly_row(_trim3(phi_row), 0)
    phi4 = _poly_row(phi_row, 0)
    u04 = _poly_row(u0_row, 0)
    psi4 = _poly_row(psi_row, 1)

    # Mode 0: V=u0+epsilon*Psi, K=phi*V.
    v4 = _poly_add(u04, psi4)
    mode0_m4 = v4
    mode0_h4 = phi4
    mode0_k4 = _poly_product(phi4, v4)
    mode0_w = _w_polynomial(z, delta, mode0_m4)
    mode0_angular = _angular_polynomial(z, delta, gauge, mode0_h4, mode0_k4)
    mode0 = _poly_add(
        _poly_scale(_poly_product(phi3, mode0_w), -2),
        mode0_angular,
    )

    # Mode 1: direct m difference; its H/K difference is exactly zero.
    mdiff4 = _poly_row(_jet_add(mpsi_row, _jet_scale(psi_row, -1)), 1)
    mdiff3 = _poly_trim3(mdiff4)
    mdiff_z = _poly_shift(mdiff4)
    one_minus_delta = _poly_add(one3, _poly_scale(delta, -1))
    mode1_inside = _poly_add(
        _poly_product(_poly_product(one_minus_delta, z_poly), mdiff3),
        _poly_product(d_poly, mdiff_z),
    )
    mode1 = _poly_scale(_poly_product(phi3, mode1_inside), 2)

    # Mode 2: direct H/K difference and zero W difference.
    hdiff4 = _poly_add(_poly_row(h_row, 0), _poly_scale(phi4, -1))
    # Hdiff has degree zero.  Form Kdiff with its degree-zero and degree-one
    # pieces directly so that no independent subtraction of K is introduced.
    u0_hdiff = _poly_product(u04, hdiff4)
    phi_psi = _poly_product(phi4, psi4)
    kpsi_poly = _poly_row(kpsi_row, 1)
    kdiff_eps = _poly_add(kpsi_poly, _poly_scale(phi_psi, -1))
    kdiff4 = _poly_add(u0_hdiff, kdiff_eps)
    mode2 = _angular_polynomial(
        z,
        delta,
        gauge,
        hdiff4,
        kdiff4,
    )

    z_squared = _poly_product(z_poly, z_poly)
    l_factor = _poly_add(one3, _poly_scale(_poly_product(delta, z_squared), -1))
    denominator = _poly_scale(_poly_product(phi3, l_factor), 2)
    return {
        "numerator_modes": [mode0, mode1, mode2],
        "denominator_sectors": denominator,
    }


__all__ = ["angular_mode_numerator_epsilon_sectors"]
