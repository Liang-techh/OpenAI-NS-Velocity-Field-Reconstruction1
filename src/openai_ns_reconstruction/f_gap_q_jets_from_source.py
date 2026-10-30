"""Source-formula algebra adapter for fixed-``xi`` F-gap Q jets.

This module deliberately contains no source-graph imports, node traversal, or
numerical source evaluation.  It implements only the native radial formula

    Q = (2 * Z * U - (1 - delta) * Z * m - (1 - Z**2) * m_Z)
        / (1 - delta * Z**2)

and its ordinary fixed-``xi`` ``Z`` derivatives through the orders needed for
``(U, U_Z, Q, Q_Z, Q_ZZ, Q_y, Q_yZ)``.  The caller supplies actual jet values
or caller-directed interval objects.  Opaque source node ids, ``FunctionRef``
objects, and source-owned claims are intentionally outside this adapter.

The fixed-``xi`` coordinate convention is a caller premise: all rows supplied
below are ordinary derivatives with respect to ``Z`` while ``xi`` is held
fixed.  Endpoint transport, ``xi_Z`` terms, ``R`` ownership, ``Pstar``
scaling, and source capture/namespace authentication are not performed here.

The arithmetic is intentionally generic.  No conversion to ``float``, no
midpoint extraction, and no interval narrowing is performed; ``Fraction`` or
an interval type implementing the ordinary arithmetic operators can be kept
throughout.  The denominator callback is external: it receives the computed
``L`` and must return the built-in singleton ``True``.  The callback result is
only an external gate; this adapter does not claim that it proves nonzero-ness.

The ``Q_y`` and ``Q_yZ`` formulas additionally use the explicit caller
premise ``delta_y = 0`` (and therefore ``(delta_y)_Z = 0``).  The supplied
``Uy1``, ``my1``, and ``my2`` are likewise caller-supplied fixed-``xi``
ordinary-``Z`` derivatives of the corresponding ``y``-jets, e.g.
``my1 = d_Z(d_y m)`` and ``my2 = d_Z^2(d_y m)``.  Neither premise is proved by
this module.

Independent hand check (all inputs written as exact fractions):

    Z = 1/2
    U  = (U0,U1,U2) = (1,2,3)
    m  = (m0,m1,m2,m3) = (2,1,3,4)
    Uy = (Uy0,Uy1) = (5,6)
    my = (my0,my1,my2) = (7,8,9)
    delta = (delta0,delta1,delta2) = (1/2,1/3,1/4)

gives, in return order,

    (1, 2, -2/7, 40/21, 9019/441, -22/7, 266/21).

The canonical hand-check string whose SHA-256 is exposed below is:

    Z=1/2;U=(1,2,3);m=(2,1,3,4);Uy=(5,6);my=(7,8,9);delta=(1/2,1/3,1/4);jets=(1,2,-2/7,40/21,9019/441,-22/7,266/21)

The digest is kept as a checked-in constant so an independent caller can
recompute it without importing or evaluating this module.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Generic, Protocol, TypeVar


Scalar = TypeVar("Scalar")


@dataclass(frozen=True, slots=True)
class FGapFixedXiOrdinaryZInputs(Generic[Scalar]):
    """Caller-provided fixed-``xi`` ordinary-``Z`` jet rows.

    These are values supplied by the caller, not handles extracted from a
    source graph.  The row convention is ordinary derivatives at fixed
    ``xi``, so for example ``m2`` means ``d^2 m / dZ^2`` and is not a divided
    Taylor coefficient.  The ``Q_y`` branch requires the separate caller
    premise ``delta_y = 0``; this packet does not prove that premise.
    """

    U0: Scalar
    U1: Scalar
    U2: Scalar
    m0: Scalar
    m1: Scalar
    m2: Scalar
    m3: Scalar
    Uy0: Scalar
    Uy1: Scalar
    my0: Scalar
    my1: Scalar
    my2: Scalar
    delta0: Scalar
    delta1: Scalar
    delta2: Scalar
    Z: Scalar


@dataclass(frozen=True, slots=True)
class FGapQJets(Generic[Scalar]):
    """Seven fixed-``xi`` ordinary-``Z`` outputs in canonical order."""

    U: Scalar
    U_Z: Scalar
    Q: Scalar
    Q_Z: Scalar
    Q_ZZ: Scalar
    Q_y: Scalar
    Q_yZ: Scalar

    def as_tuple(self) -> tuple[Scalar, Scalar, Scalar, Scalar, Scalar, Scalar, Scalar]:
        """Return ``(U,U_Z,Q,Q_Z,Q_ZZ,Q_y,Q_yZ)`` without changing values."""

        return (self.U, self.U_Z, self.Q, self.Q_Z, self.Q_ZZ, self.Q_y, self.Q_yZ)


class FixedXiOrdinaryZJetProvider(Protocol[Scalar]):
    """Future provider seam for explicitly supplied real jet values.

    A future source-owned provider may implement this method and return actual
    values in ``FGapFixedXiOrdinaryZInputs``.  The caller must separately pass
    ``delta_y_zero_premise=True`` when using the ``Q_y`` branch.  This protocol
    does not authenticate ownership, source namespace, capture, or numerical
    validity; those remain outside this formula adapter.
    """

    def fixed_xi_ordinary_z_inputs(self) -> FGapFixedXiOrdinaryZInputs[Scalar]:
        ...


HAND_CHECK_CANONICAL = (
    "Z=1/2;U=(1,2,3);m=(2,1,3,4);Uy=(5,6);my=(7,8,9);"
    "delta=(1/2,1/3,1/4);jets=(1,2,-2/7,40/21,9019/441,-22/7,266/21)"
)

# SHA-256 of HAND_CHECK_CANONICAL, UTF-8 encoded.  This is a hand-check
# identifier only; it is not a source artifact hash and carries no ownership
# or numerical certification.
HAND_CHECK_SHA256 = "80e3d063e65b2dd827935ee59dba65187d789dbeeb86d994cf8d2a46945579b7"

ADAPTER_KIND = "source_formula_algebra_only"
SOURCE_OWNERSHIP_ASSERTED = False
# This only records that the primitive guard rejects plain integer inputs.  It
# cannot inspect arbitrary objects for graph handles; provenance remains an
# unauthenticated caller responsibility.
OPAQUE_SOURCE_IDS_ACCEPTED = False
OPAQUE_SOURCE_IDS_AUTHENTICATED = False
SOURCE_VALUE_PROVENANCE_CALLER_RESPONSIBLE = True
PSTAR_APPLIED = False


def _reject_non_exact_scalar(name: str, value: Any) -> None:
    """Reject primitive values that could hide a node id or float coercion.

    Exact integral values should be passed as ``Fraction(n, 1)``.  This keeps
    the adapter from silently interpreting an opaque integer node id as a
    source value.  Non-primitive caller-directed interval objects are retained
    and are never converted or inspected for a midpoint.
    """

    if isinstance(value, bool) or isinstance(value, int):
        raise TypeError(
            f"{name} must be Fraction or a caller-directed interval object; "
            "plain integers are not accepted as source values or node ids"
        )
    if isinstance(value, float):
        raise TypeError(f"{name} must not be float; pass Fraction or an interval object")


def _validate_inputs(inputs: FGapFixedXiOrdinaryZInputs[Scalar]) -> None:
    for name in inputs.__slots__:
        _reject_non_exact_scalar(name, getattr(inputs, name))


def reconstruct_f_gap_fixed_xi_q_jets(
    U0: Scalar,
    U1: Scalar,
    U2: Scalar,
    m0: Scalar,
    m1: Scalar,
    m2: Scalar,
    m3: Scalar,
    Uy0: Scalar,
    Uy1: Scalar,
    my0: Scalar,
    my1: Scalar,
    my2: Scalar,
    delta0: Scalar,
    delta1: Scalar,
    delta2: Scalar,
    Z: Scalar,
    *,
    denominator_nonzero_gate: Callable[[Scalar], bool],
    delta_y_zero_premise: bool,
) -> FGapQJets[Scalar]:
    """Rebuild the seven native radial-Q jets at fixed ``xi``.

    ``delta0``, ``delta1``, and ``delta2`` are respectively ``delta``,
    ``delta_Z``, and ``delta_ZZ``.  No ``delta_Z = 0`` simplification is made.
    Likewise, no gap identity such as ``U = 0`` or ``m_y = -m`` is imposed.
    The caller must set ``delta_y_zero_premise`` to the built-in ``True`` to
    acknowledge that the ``Q_y`` and ``Q_yZ`` formulas omit all
    ``delta_y``/``(delta_y)_Z`` terms.  This is a usage premise, not a proof of
    source identity.  ``Uy1``, ``my1``, and ``my2`` must be fixed-``xi``
    ordinary-``Z`` mixed jets supplied by that same caller.

    ``denominator_nonzero_gate`` is a callable receiving the computed ``L``.
    The function proceeds only when its return value ``is True`` (the built-in
    singleton).  It does not compare an interval denominator with zero and does
    not claim that the callback proves nonzero-ness.
    """

    if not callable(denominator_nonzero_gate):
        raise TypeError("denominator_nonzero_gate must be callable")
    if delta_y_zero_premise is not True:
        raise ValueError("delta_y_zero_premise must be the caller's literal True premise")

    inputs = FGapFixedXiOrdinaryZInputs(
        U0=U0,
        U1=U1,
        U2=U2,
        m0=m0,
        m1=m1,
        m2=m2,
        m3=m3,
        Uy0=Uy0,
        Uy1=Uy1,
        my0=my0,
        my1=my1,
        my2=my2,
        delta0=delta0,
        delta1=delta1,
        delta2=delta2,
        Z=Z,
    )
    _validate_inputs(inputs)

    # Native radial-Q ingredients and their ordinary Z derivatives.
    h = 1 - delta0
    h1 = -delta1
    h2 = -delta2
    d = 1 - Z * Z
    d1 = -2 * Z
    d2 = -2
    L = 1 - delta0 * Z * Z
    L1 = -(delta1 * Z * Z + 2 * delta0 * Z)
    L2 = -(delta2 * Z * Z + 4 * delta1 * Z + 2 * delta0)

    gate_result = denominator_nonzero_gate(L)
    if gate_result is not True:
        raise ValueError("denominator_nonzero_gate(L) must return built-in True")

    N = 2 * Z * U0 - h * Z * m0 - d * m1
    N1 = (
        2 * U0
        + 2 * Z * U1
        - h1 * Z * m0
        - h * (m0 + Z * m1)
        - d1 * m1
        - d * m2
    )
    N2 = (
        4 * U1
        + 2 * Z * U2
        - h2 * Z * m0
        - 2 * h1 * (m0 + Z * m1)
        - h * (2 * m1 + Z * m2)
        - d2 * m1
        - 2 * d1 * m2
        - d * m3
    )

    Ny = 2 * Z * Uy0 - h * Z * my0 - d * my1
    Ny1 = (
        2 * Uy0
        + 2 * Z * Uy1
        - h1 * Z * my0
        - h * (my0 + Z * my1)
        - d1 * my1
        - d * my2
    )

    Q = N / L
    Q_Z = N1 / L - N * L1 / (L * L)
    Q_ZZ = (
        N2 / L
        - 2 * N1 * L1 / (L * L)
        - N * L2 / (L * L)
        + 2 * N * L1 * L1 / (L * L * L)
    )
    Q_y = Ny / L
    Q_yZ = Ny1 / L - Ny * L1 / (L * L)

    return FGapQJets(U=U0, U_Z=U1, Q=Q, Q_Z=Q_Z, Q_ZZ=Q_ZZ, Q_y=Q_y, Q_yZ=Q_yZ)


def reconstruct_f_gap_fixed_xi_q_jets_from_provider(
    provider: FixedXiOrdinaryZJetProvider[Scalar],
    *,
    denominator_nonzero_gate: Callable[[Scalar], bool],
    delta_y_zero_premise: bool,
) -> FGapQJets[Scalar]:
    """Apply the same formula to an explicitly supplied future provider.

    The provider is only a value-delivery seam.  This function does not infer
    that it is source-owned and cannot authenticate arbitrary graph handles.
    The caller remains responsible for the ``delta_y = 0`` premise and for
    value provenance.
    """

    inputs = provider.fixed_xi_ordinary_z_inputs()
    if not isinstance(inputs, FGapFixedXiOrdinaryZInputs):
        raise TypeError("provider must return FGapFixedXiOrdinaryZInputs")
    return reconstruct_f_gap_fixed_xi_q_jets(
        inputs.U0,
        inputs.U1,
        inputs.U2,
        inputs.m0,
        inputs.m1,
        inputs.m2,
        inputs.m3,
        inputs.Uy0,
        inputs.Uy1,
        inputs.my0,
        inputs.my1,
        inputs.my2,
        inputs.delta0,
        inputs.delta1,
        inputs.delta2,
        inputs.Z,
        denominator_nonzero_gate=denominator_nonzero_gate,
        delta_y_zero_premise=delta_y_zero_premise,
    )


__all__ = [
    "ADAPTER_KIND",
    "FGapFixedXiOrdinaryZInputs",
    "FGapQJets",
    "FixedXiOrdinaryZJetProvider",
    "HAND_CHECK_CANONICAL",
    "HAND_CHECK_SHA256",
    "OPAQUE_SOURCE_IDS_ACCEPTED",
    "OPAQUE_SOURCE_IDS_AUTHENTICATED",
    "PSTAR_APPLIED",
    "SOURCE_OWNERSHIP_ASSERTED",
    "SOURCE_VALUE_PROVENANCE_CALLER_RESPONSIBLE",
    "reconstruct_f_gap_fixed_xi_q_jets",
    "reconstruct_f_gap_fixed_xi_q_jets_from_provider",
]
