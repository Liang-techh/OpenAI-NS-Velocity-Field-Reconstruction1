"""Pure ordinary-``Z`` arithmetic for the native/macro endpoint defect.

The source rows supplied here are the normalized values

``B = U_y / (epsilon * hbB)``.

``native_B`` and ``macro_B`` are ordinary derivative rows, so their first
two entries are ``B`` and ``B_Z``.  The public operator forms the signed
native-minus-macro row *before* applying

``E(B) = 2*Z*B / (1-delta*Z**2)``.

Only the dimensionless ordinary ``Z0/Z1`` result is returned.  The physical
``R0 * epsilon * hbB`` scale, any ``rho2/4`` factor, and the actual-same-
radius statement belong to the caller.  This module does not load source
reports, materialize widths, select interval midpoints, or perform
quadrature.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

from mpmath.libmp import fzero, mpf_cmp


def _is_scalar(value: Any) -> bool:
    """Return whether ``value`` can be used as one ordinary scalar entry."""

    return (
        value is not None
        and not isinstance(value, (str, bytes, Mapping, Sequence))
    )


def _ordinary_row(value: Any, name: str) -> list[Any]:
    """Validate and copy the required ordinary ``Z0/Z1`` source entries."""

    if value is None or isinstance(value, (str, bytes, Mapping)):
        raise ValueError(f"{name} must contain ordinary Z0/Z1 values")
    if not isinstance(value, Sequence) or len(value) < 2:
        raise ValueError(f"{name} must contain ordinary Z0/Z1 values")
    row = list(value[:2])
    if any(not _is_scalar(item) for item in row):
        raise ValueError(f"{name} contains a missing or nested ordinary value")
    return row


def _denominator_crosses_zero(value: Any) -> bool:
    """Check a point or directed interval denominator for zero."""

    if hasattr(value, "_mpi_"):
        lower, upper = value._mpi_
        return mpf_cmp(lower, fzero) <= 0 and mpf_cmp(upper, fzero) >= 0
    return value == 0


def evaluate_seam_defect(
    ctx: Any,
    z: Any,
    delta: Any,
    native_B: Sequence[Any],
    macro_B: Sequence[Any],
    *,
    values_source_proved: Any,
    same_radius_proved: Any,
) -> dict[str, Any]:
    """Return the normalized native-minus-macro endpoint defect.

    ``native_B`` and ``macro_B`` must each provide at least ordinary ``Z0``
    and ``Z1`` entries.  ``values_source_proved`` records the separate proof
    that the native and macro value rows (including ``U`` and ``m`` through
    ``Z2``) are exact source rows; ``same_radius_proved`` records the common
    physical-radius proof.  Both are deliberately exact boolean gates.

    ``ctx`` is retained in the API so point and directed interval callers use
    the same arithmetic contract.  Arithmetic is delegated to the scalar
    values, which makes a ``mpmath`` point context and ``MPIntervalContext``
    behave identically apart from enclosure widening.
    """

    # ``is True`` rejects False, integer 1, interval one, and other merely
    # truthy objects.  A source seam cannot be admitted from a loose flag.
    if values_source_proved is not True:
        raise ValueError("exact native/macro value-source proof is required")
    if same_radius_proved is not True:
        raise ValueError("an exact same-radius proof is required")
    if not _is_scalar(z) or not _is_scalar(delta):
        raise ValueError("z and delta must be scalar ordinary values")

    native = _ordinary_row(native_B, "native_B")
    macro = _ordinary_row(macro_B, "macro_B")

    denominator = 1 - delta * z * z
    if _denominator_crosses_zero(denominator):
        raise ValueError("E(B) denominator crosses zero")

    # Preserve the source correlation: subtract the two signed rows first,
    # then apply the endpoint operator to that one difference.
    source_difference_rows = [
        native[0] - macro[0],
        native[1] - macro[1],
    ]
    difference, difference_z = source_difference_rows
    defect_C0_C1 = [
        2 * z * difference / denominator,
        2 * (difference + z * difference_z) / denominator
        + 4 * delta * z * z * difference / (denominator * denominator),
    ]

    return {
        "defect_C0_C1": defect_C0_C1,
        "source_difference_rows": source_difference_rows,
        "certificate": {
            "passed": True,
            "source_Z_orders": [0, 1],
            "ordinary_Z_orders": [0, 1],
            "values_source_proved": True,
            "same_radius_proved": True,
            "signed_difference_before_operator": True,
            "widths_materialized": False,
            "quadrature_used": False,
            "physical_scale_attached": False,
        },
    }


__all__ = ["evaluate_seam_defect"]
