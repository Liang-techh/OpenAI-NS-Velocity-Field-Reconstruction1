"""Correlated signed density difference for the pulse source formula.

This module is intentionally independent of the older pulse kernels.  It
accepts one frame/source jet and a caller-provided jet difference at the same
source query.  The difference jet is treated as already source-consistent;
the module neither reconstructs the physical jet nor checks source identity.

The source density is

``S = A/L + B + U*C/L``

with

``A = Q + (1-delta)*Z*Q_Z/2 + Q_y``
``B = Q*(Q/2 + Q_y)``
``C = -2*Z*Q + (1-Z**2)*Q_Z - 2*Z*Q_y``
``L = 1 - delta*Z**2``.

Only the correlated difference is formed here.  In particular, this module
does not evaluate ``S(frame + delta) - S(frame)`` and does not evaluate two
separate ``S`` values before subtracting them.  Numeric inputs are left in
their original arithmetic domain, so ``fractions.Fraction``, ``mpmath.mpf``
and other compatible exact or high-precision scalar types are supported.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any


Numeric = Any


@dataclass(frozen=True)
class FieldJets:
    """One frame/source field jet at a fixed source query.

    The seven entries are values, not interval upper bounds.  No conversion
    or finiteness coercion is performed; the caller controls the numeric
    domain and any source-level consistency proof for the difference jet.
    """

    U: Numeric
    U_Z: Numeric
    Q: Numeric
    Q_Z: Numeric
    Q_ZZ: Numeric
    Q_y: Numeric
    Q_yZ: Numeric


@dataclass(frozen=True)
class DensityDelta:
    """Correlated ``Delta S`` and ordinary ``Z1`` result.

    The component fields are retained so a caller can audit the linear,
    quadratic and ``U*C`` pieces without reconstructing either source side.
    """

    delta_S: Numeric
    delta_S_Z: Numeric
    L: Numeric
    L_Z: Numeric
    delta_A: Numeric
    delta_A_Z: Numeric
    delta_B: Numeric
    delta_B_Z: Numeric
    delta_C: Numeric
    delta_C_Z: Numeric
    delta_UC: Numeric
    delta_UC_Z: Numeric

    @property
    def z1(self) -> Numeric:
        """Short alias for the ordinary ``Z1`` correlated density delta."""

        return self.delta_S_Z


def _as_fields(value: FieldJets | dict[str, Numeric], *, label: str) -> FieldJets:
    """Copy a field mapping without coercing its scalar numeric types."""

    if isinstance(value, FieldJets):
        return value
    required = ("U", "U_Z", "Q", "Q_Z", "Q_ZZ", "Q_y", "Q_yZ")
    missing = tuple(name for name in required if name not in value)
    if missing:
        raise ValueError(f"{label} is missing jets: {', '.join(missing)}")
    return FieldJets(**{name: value[name] for name in required})


def source_density_delta_and_Z1(
    frame_jets: FieldJets | dict[str, Numeric],
    delta_jets: FieldJets | dict[str, Numeric],
    *,
    Z: Numeric,
    delta: Numeric,
    delta_z_independent: bool,
) -> DensityDelta:
    """Compute correlated ``Delta S`` and ``Delta S_Z`` directly.

    ``frame_jets`` supplies ``(U, U_Z, Q, Q_Z, Q_ZZ, Q_y, Q_yZ)`` and
    ``delta_jets`` supplies the corresponding signed differences.  The
    physical side is conceptually frame plus delta, but that sum is never
    materialized here.

    The quadratic difference is expanded as required by the source algebra:

    ``Delta B = Delta Q*(Q + Q_y) + Q*Delta Q_y``
    ``           + Delta Q**2/2 + Delta Q*Delta Q_y``.

    The product difference is similarly formed as
    ``Delta(U*C) = Delta U*C + U*Delta C + Delta U*Delta C``.  Their ordinary
    ``Z1`` derivatives are obtained by differentiating those expanded
    expressions, including ``L_Z = -2*delta*Z``.  The caller must explicitly
    provide the source premise ``delta_z_independent=True``: a Z-dependent
    parameter ``delta`` would require additional chain-rule terms.  The
    difference jets themselves may depend on Z and supply their derivatives.
    No source equality or zero claim is inferred from a zero
    difference supplied by the caller.
    """

    if delta_z_independent is not True:
        raise ValueError("delta_z_independent=True is required")

    q = _as_fields(frame_jets, label="frame_jets")
    dq = _as_fields(delta_jets, label="delta_jets")

    # Division by the integer 2 stays in the caller's numeric domain for
    # Fraction, mpmath.mpf and other compatible scalar types.
    h = 1 - delta
    L = 1 - delta * Z * Z
    if L == 0:
        raise ValueError("L = 1 - delta*Z**2 must be nonzero")
    L_Z = -2 * delta * Z

    # Linear source numerator A and its ordinary Z derivative.
    delta_A = dq.Q + h * Z * dq.Q_Z / 2 + dq.Q_y
    delta_A_Z = (
        dq.Q_Z
        + h * (dq.Q_Z + Z * dq.Q_ZZ) / 2
        + dq.Q_yZ
    )

    # B = Q*(Q/2 + Q_y), expanded in the supplied frame/difference basis.
    delta_B = (
        dq.Q * (q.Q + q.Q_y)
        + q.Q * dq.Q_y
        + dq.Q * dq.Q / 2
        + dq.Q * dq.Q_y
    )
    delta_B_Z = (
        dq.Q_Z * (q.Q + q.Q_y)
        + dq.Q * (q.Q_Z + q.Q_yZ)
        + q.Q_Z * dq.Q_y
        + q.Q * dq.Q_yZ
        + dq.Q * dq.Q_Z
        + dq.Q_Z * dq.Q_y
        + dq.Q * dq.Q_yZ
    )

    # C is linear in Q, Q_Z and Q_y.
    delta_C = (
        -2 * Z * dq.Q
        + (1 - Z * Z) * dq.Q_Z
        - 2 * Z * dq.Q_y
    )
    delta_C_Z = (
        -2 * dq.Q
        - 4 * Z * dq.Q_Z
        + (1 - Z * Z) * dq.Q_ZZ
        - 2 * dq.Q_y
        - 2 * Z * dq.Q_yZ
    )

    C = (
        -2 * Z * q.Q
        + (1 - Z * Z) * q.Q_Z
        - 2 * Z * q.Q_y
    )
    C_Z = (
        -2 * q.Q
        - 4 * Z * q.Q_Z
        + (1 - Z * Z) * q.Q_ZZ
        - 2 * q.Q_y
        - 2 * Z * q.Q_yZ
    )

    # Product difference Delta(U*C), still without constructing physical U/C.
    delta_UC = dq.U * C + q.U * delta_C + dq.U * delta_C
    delta_UC_Z = (
        dq.U_Z * C
        + dq.U * C_Z
        + q.U_Z * delta_C
        + q.U * delta_C_Z
        + dq.U_Z * delta_C
        + dq.U * delta_C_Z
    )

    delta_S = delta_A / L + delta_B + delta_UC / L
    delta_S_Z = (
        delta_A_Z / L
        - delta_A * L_Z / (L * L)
        + delta_B_Z
        + delta_UC_Z / L
        - delta_UC * L_Z / (L * L)
    )

    return DensityDelta(
        delta_S=delta_S,
        delta_S_Z=delta_S_Z,
        L=L,
        L_Z=L_Z,
        delta_A=delta_A,
        delta_A_Z=delta_A_Z,
        delta_B=delta_B,
        delta_B_Z=delta_B_Z,
        delta_C=delta_C,
        delta_C_Z=delta_C_Z,
        delta_UC=delta_UC,
        delta_UC_Z=delta_UC_Z,
    )


# A descriptive alias for callers that prefer ``compute`` terminology.
compute_source_density_delta = source_density_delta_and_Z1


__all__ = [
    "DensityDelta",
    "FieldJets",
    "compute_source_density_delta",
    "source_density_delta_and_Z1",
]
