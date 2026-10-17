"""Pure ordinary-``Z`` arithmetic for the original swirl source.

The source inputs in this module are already normalized by the relative
``F0**2`` factor.  In particular, ``b`` and ``c`` are *dressed* ordinary
derivative rows.  This module does not know about the external base
amplitude, scales, radius variables, source data, or integration.  That
separation is deliberate: differentiating a frozen external ``F0`` here
would apply the amplitude gauge twice.

Ordinary rows contain derivatives rather than Taylor coefficients.  Thus a
row through ``Z3`` is ``[f, f_Z, f_ZZ, f_ZZZ]``.  The swirl drive needs the
``Z`` derivatives of ``b`` and ``c`` and therefore requires the complete row
through ``Z3`` even though it returns only ``Z0..Z2``.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

from .original_source_jet_kernel import ordinary_product, ordinary_quotient


def _ordinary_jet(value: Any, order: int, name: str) -> list[Any]:
    """Validate and copy one complete ordinary derivative row."""
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


def _z_jet(z: Any) -> list[Any]:
    return [z, 1, 0]


def _d_jet(z: Any) -> list[Any]:
    """Return the ordinary jet through ``Z2`` of ``d=1-z**2``."""
    return [1 - z * z, -2 * z, -2]


def _l_jet(z: Any, delta: Any) -> list[Any]:
    """Return the ordinary jet through ``Z2`` of ``L=1-delta*z**2``."""
    return [1 - delta * z * z, -2 * delta * z, -2 * delta]


def dress_ordinary_jet(
    relative_amplitude: Sequence[Any],
    raw_moment: Sequence[Any],
) -> list[Any]:
    """Dress a raw ordinary moment by a relative amplitude through ``Z3``.

    Both inputs are complete four-entry rows.  The result is the ordinary
    derivative row of their product, using the binomial product rule,

    ``(a m)^(n) = sum(k=0..n) binom(n,k) a^(k) m^(n-k)``.

    The routine is intentionally generic over mpmath point and directed
    interval values and performs no conversion to Python ``float``.
    """
    amplitude = _ordinary_jet(relative_amplitude, 3, "relative_amplitude")
    moment = _ordinary_jet(raw_moment, 3, "raw_moment")
    # Small fixed binomial table avoids importing a numerical package and
    # keeps every operation delegated to the supplied scalar values.
    binomial = ((1,), (1, 1), (1, 2, 1), (1, 3, 3, 1))
    return [
        sum(
            binomial[order][k] * amplitude[k] * moment[order - k]
            for k in range(order + 1)
        )
        for order in range(4)
    ]


# A concise alias is useful at call sites that already use ``ordinary`` for
# derivative rows.  Both names refer to the one implementation above.
ordinary_dress = dress_ordinary_jet
dress_ordinary = dress_ordinary_jet


def _swirl_numerator_ordinary(
    z: Any,
    delta: Any,
    b: Sequence[Any],
    c: Sequence[Any],
) -> list[Any]:
    """Return the swirl numerator derivatives ``Z0..Z2``."""
    b_row = _ordinary_jet(b, 3, "b")
    c_row = _ordinary_jet(c, 3, "c")
    if z is None:
        raise ValueError("z is required")
    if delta is None:
        raise ValueError("delta is required")

    z_row = _z_jet(z)
    d_row = _d_jet(z)
    b_z = b_row[1:]
    c_z = c_row[1:]

    # N_swirl = -2*delta*z*b + d*b_Z + 2*(1+delta)*z*c - d*c_Z.
    return _jet_add(
        _jet_scale(ordinary_product(z_row, b_row[:3]), -2 * delta),
        ordinary_product(d_row, b_z[:3]),
        _jet_scale(ordinary_product(z_row, c_row[:3]), 2 * (1 + delta)),
        _jet_scale(ordinary_product(d_row, c_z[:3]), -1),
    )


def swirl_numerator_delta_sectors_ordinary(
    z: Any,
    b: Sequence[Any],
    c: Sequence[Any],
) -> dict[int, list[Any]]:
    """Return the ordinary ``Z0..Z2`` numerator sectors in ``delta``.

    For dressed rows ``b`` and ``c`` through ``Z3``, the swirl numerator is

    ``N = N0 + delta*Ndelta``,

    with ``N0 = (1-z**2)*(b_Z-c_Z) + 2*z*c`` and
    ``Ndelta = 2*z*(c-b)``.  The common denominator ``2*(1-delta*z**2)``
    remains the caller's responsibility.  Keeping these pieces separate is
    necessary when a caller combines a factored ``delta`` with other source
    powers before applying scales.
    """
    b_row = _ordinary_jet(b, 3, "b")
    c_row = _ordinary_jet(c, 3, "c")
    if z is None:
        raise ValueError("z is required")

    z_row = _z_jet(z)
    d_row = _d_jet(z)
    b_z = b_row[1:]
    c_z = c_row[1:]
    d_b_z = ordinary_product(d_row, b_z[:3])
    d_c_z = ordinary_product(d_row, c_z[:3])
    z_b = ordinary_product(z_row, b_row[:3])
    z_c = ordinary_product(z_row, c_row[:3])
    return {
        0: _jet_add(
            d_b_z,
            _jet_scale(d_c_z, -1),
            _jet_scale(z_c, 2),
        ),
        1: _jet_add(
            _jet_scale(z_c, 2),
            _jet_scale(z_b, -2),
        ),
    }


def swirl_drive_ordinary(
    z: Any,
    delta: Any,
    b: Sequence[Any],
    c: Sequence[Any],
) -> list[Any]:
    """Compute ordinary ``Z0..Z2`` of the normalized swirl drive.

    The inputs ``b`` and ``c`` must already be dressed relative ``F0**2``
    rows through ``Z3``.  The returned row is the ordinary derivative of

    ``[-2*delta*z*b + (1-z**2)*b_Z + 2*(1+delta)*z*c
       - (1-z**2)*c_Z] / [2*(1-delta*z**2)]``.

    No external amplitude is differentiated or applied here.
    """
    numerator = _swirl_numerator_ordinary(z, delta, b, c)
    denominator = _jet_scale(_l_jet(z, delta), 2)
    return ordinary_quotient(numerator, denominator)


def swirl_moment_modes_ordinary(
    z: Any,
    delta: Any,
    target_square: Sequence[Any],
    dB: Sequence[Any],
    dC: Sequence[Any],
) -> list[list[Any]]:
    """Return the three same-source swirl modal direction numerators.

    Every input is a dressed four-entry ordinary row.  The returned rows are
    ``[base, j1, j2]``, each of length three and each still a numerator (the
    common ``2*(1-delta*z**2)`` denominator is intentionally not divided).
    The directions share one source target:

    * base: ``b = target_square/2``, ``c = target_square``;
    * ``j1``: ``b = 0``, ``c = dC`` where ``dC = C-target_square``;
    * ``j2``: ``b = dB`` where ``dB = B-target_square/2``, ``c = 0``.

    This construction avoids subtracting independently enclosed copies of
    the base target and preserves the source correlation in both directions.
    """
    target = _ordinary_jet(target_square, 3, "target_square")
    delta_b = _ordinary_jet(dB, 3, "dB")
    delta_c = _ordinary_jet(dC, 3, "dC")
    zero = [0, 0, 0, 0]
    base_b = [value / 2 for value in target]
    return [
        _swirl_numerator_ordinary(z, delta, base_b, target),
        _swirl_numerator_ordinary(z, delta, zero, delta_c),
        _swirl_numerator_ordinary(z, delta, delta_b, zero),
    ]


def swirl_moment_mode_numerator_delta_sectors_ordinary(
    z: Any,
    target_square: Sequence[Any],
    dB: Sequence[Any],
    dC: Sequence[Any],
) -> list[dict[int, list[Any]]]:
    """Return delta sectors for the three same-source modal numerators.

    The modes are constructed directly from one shared target row:

    * base: ``b = target_square/2``, ``c = target_square``;
    * ``j1``: ``b = 0``, ``c = dC``;
    * ``j2``: ``b = dB``, ``c = 0``.

    The zero rows are exact inputs.  No independently enclosed base target is
    subtracted to form either direction.  Each returned item has keys ``0``
    and ``1`` for ``N0`` and ``Ndelta``; the common ``2*(1-delta*z**2)``
    denominator is intentionally omitted.
    """
    target = _ordinary_jet(target_square, 3, "target_square")
    delta_b = _ordinary_jet(dB, 3, "dB")
    delta_c = _ordinary_jet(dC, 3, "dC")
    zero = [0, 0, 0, 0]
    base_b = [value / 2 for value in target]
    return [
        swirl_numerator_delta_sectors_ordinary(z, base_b, target),
        swirl_numerator_delta_sectors_ordinary(z, zero, delta_c),
        swirl_numerator_delta_sectors_ordinary(z, delta_b, zero),
    ]


__all__ = [
    "dress_ordinary_jet",
    "ordinary_dress",
    "dress_ordinary",
    "swirl_numerator_delta_sectors_ordinary",
    "swirl_drive_ordinary",
    "swirl_moment_modes_ordinary",
    "swirl_moment_mode_numerator_delta_sectors_ordinary",
]
