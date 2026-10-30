"""Generic signed ``P-F`` pulse-density algebra.

This module evaluates the algebraic part of the pulse source density and its
ordinary-``Z`` derivative while retaining the arithmetic objects supplied by
the caller.  The intended objects are exact rationals (for example
``fractions.Fraction``) or a directed/interval arithmetic type implementing
the corresponding Python operators.  Built-in ``int`` inputs are upgraded to
``Fraction`` before arithmetic so expressions such as ``q / 2`` cannot create
a binary64 ``float``.  Built-in ``float`` and ``bool`` numeric inputs are
rejected.  Existing non-built-in rational/interval objects are retained
unchanged.  No value is converted to ``float``, ``mpmath.mpf``, a midpoint,
or a raw ``P-F`` subtraction.

The frame and signed-difference jets are supplied independently.  In
particular, ``signed_delta_jets`` must already be the direct signed ``P-F``
jet; this module never constructs it from two raw field packets.  The
denominator proof is also an explicit caller boundary: ``denominator_nonzero``
must return the built-in boolean ``True`` for the computed ``L`` before any
division is attempted.  The callback is not itself proved by this module.

The result is arithmetic infrastructure only.  It does not bind a physical
source, establish F-source identity, perform quadrature, provide a numerical
error certificate, or claim denominator safety beyond the caller's strict
callback assertion.  Endpoint/radius/measure factors are deliberately outside
this module.

For a hand-checkable synthetic example, set ``Z=1/2``, ``delta=0``,
``delta_Z=1/2`` and set all fourteen frame/difference jet components to ``1``.
Then ``h=1``, ``h_Z=-1/2``, ``L=1``, ``L_Z=-1/8`` and the returned values are
``S_F=5/2``, ``S_F_Z=-7/4``, ``dS=3``, ``dS_Z=-177/16``.  This example is an
independent algebra check, not a source or certification claim.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from fractions import Fraction
from typing import Any, Callable


FIELD_NAMES = ("U", "U_Z", "Q", "Q_Z", "Q_ZZ", "Q_y", "Q_yZ")
SIGNED_FIELD_NAMES = (
    "dU",
    "dU_Z",
    "dQ",
    "dQ_Z",
    "dQ_ZZ",
    "dQ_y",
    "dQ_yZ",
)


def _numeric_input(value: Any, label: str) -> Any:
    """Keep the caller's arithmetic domain without an implicit float path.

    Only the exact built-in ``int`` case is normalized, because Python's
    ``int / int`` would otherwise produce a binary64 float in terms such as
    ``q / 2``.  Built-in ``float`` and ``bool`` are rejected at the boundary.
    Non-built-in objects, including Fraction and directed interval types, are
    returned unchanged and must implement the arithmetic used below.
    """

    if isinstance(value, bool):
        raise TypeError(f"{label} must not be bool")
    if type(value) is float:
        raise TypeError(f"{label} must not be a built-in float")
    if type(value) is int:
        return Fraction(value, 1)
    return value


@dataclass(frozen=True)
class FrameJets:
    """Seven caller-owned frame jets for one source-chart point.

    Values are intentionally typed as ``Any``: the arithmetic domain is
    chosen by the caller.  Built-in integers are normalized to exact
    ``Fraction`` values; existing non-built-in arithmetic objects are kept.
    """

    U: Any
    U_Z: Any
    Q: Any
    Q_Z: Any
    Q_ZZ: Any
    Q_y: Any
    Q_yZ: Any

    def __post_init__(self) -> None:
        for name in FIELD_NAMES:
            object.__setattr__(
                self,
                name,
                _numeric_input(getattr(self, name), f"frame_jets.{name}"),
            )


@dataclass(frozen=True)
class SignedDeltaJets:
    """Seven caller-owned direct signed ``P-F`` jets.

    The ``d`` prefix is part of the interface contract.  Constructing these
    values from separate raw P and F inputs is outside this module.
    """

    dU: Any
    dU_Z: Any
    dQ: Any
    dQ_Z: Any
    dQ_ZZ: Any
    dQ_y: Any
    dQ_yZ: Any

    def __post_init__(self) -> None:
        for name in SIGNED_FIELD_NAMES:
            object.__setattr__(
                self,
                name,
                _numeric_input(getattr(self, name), f"signed_delta_jets.{name}"),
            )


@dataclass(frozen=True)
class SignedIntervalAlgebraResult:
    """Unscaled frame/difference density terms and the shared denominator."""

    S_F: Any
    S_F_Z: Any
    dS: Any
    dS_Z: Any
    L: Any

    @property
    def delta_S(self) -> Any:
        """Alias matching callers that spell the signed term ``delta_S``."""

        return self.dS

    @property
    def delta_S_Z(self) -> Any:
        """Alias matching callers that spell the signed term ``delta_S_Z``."""

        return self.dS_Z


# Naming aliases keep the generic interface easy to wire to callers that use
# the field names from ``pulse_signed_scaled_kernel.py``.  These aliases do
# not import that module or add its scalar coercion behavior.
FieldJets = FrameJets
DeltaJets = SignedDeltaJets


def _read_components(
    value: Any,
    names: tuple[str, ...],
    *,
    label: str,
) -> tuple[Any, ...]:
    """Read named components without changing their arithmetic objects."""

    if isinstance(value, Mapping):
        missing = tuple(name for name in names if name not in value)
        if missing:
            raise ValueError(f"{label} is missing: {', '.join(missing)}")
        return tuple(value[name] for name in names)

    missing = tuple(name for name in names if not hasattr(value, name))
    if missing:
        raise TypeError(
            f"{label} must expose components: {', '.join(names)}; "
            f"missing {', '.join(missing)}"
        )
    return tuple(getattr(value, name) for name in names)


def _frame_jets(value: FrameJets | Mapping[str, Any] | Any) -> FrameJets:
    if isinstance(value, FrameJets):
        return value
    return FrameJets(*_read_components(value, FIELD_NAMES, label="frame_jets"))


def _signed_delta_jets(value: SignedDeltaJets | Mapping[str, Any] | Any) -> SignedDeltaJets:
    if isinstance(value, SignedDeltaJets):
        return value
    return SignedDeltaJets(
        *_read_components(value, SIGNED_FIELD_NAMES, label="signed_delta_jets")
    )


@dataclass(frozen=True)
class _FrameTerms:
    A: Any
    A_Z: Any
    B: Any
    B_Z: Any
    C: Any
    C_Z: Any
    numerator: Any
    numerator_Z: Any


@dataclass(frozen=True)
class _DifferenceTerms:
    dA: Any
    dA_Z: Any
    dB: Any
    dB_Z: Any
    dC: Any
    dC_Z: Any
    dUC: Any
    dUC_Z: Any
    numerator: Any
    numerator_Z: Any


def _frame_terms(
    jets: FrameJets,
    *,
    Z: Any,
    h: Any,
    h_Z: Any,
) -> _FrameTerms:
    """Build frame terms without dividing by ``L``."""

    q, q_z, q_zz, q_y, q_yz = (
        jets.Q,
        jets.Q_Z,
        jets.Q_ZZ,
        jets.Q_y,
        jets.Q_yZ,
    )
    A = q + h * Z * q_z / 2 + q_y
    A_Z = q_z + (h_Z * Z * q_z + h * q_z + h * Z * q_zz) / 2 + q_yz
    B = q * (q / 2 + q_y)
    B_Z = q_z * (q / 2 + q_y) + q * (q_z / 2 + q_yz)
    C = -2 * Z * q + (1 - Z * Z) * q_z - 2 * Z * q_y
    C_Z = (
        -2 * q
        - 4 * Z * q_z
        + (1 - Z * Z) * q_zz
        - 2 * q_y
        - 2 * Z * q_yz
    )
    numerator = A + jets.U * C
    numerator_Z = A_Z + jets.U_Z * C + jets.U * C_Z
    return _FrameTerms(A, A_Z, B, B_Z, C, C_Z, numerator, numerator_Z)


def _difference_terms(
    frame: FrameJets,
    signed: SignedDeltaJets,
    frame_terms: _FrameTerms,
    *,
    Z: Any,
    h: Any,
    h_Z: Any,
) -> _DifferenceTerms:
    """Build direct signed-difference terms, including all quadratic terms."""

    dq, dq_z, dq_zz, dq_y, dq_yz = (
        signed.dQ,
        signed.dQ_Z,
        signed.dQ_ZZ,
        signed.dQ_y,
        signed.dQ_yZ,
    )
    q, q_z, q_y, q_yz = (
        frame.Q,
        frame.Q_Z,
        frame.Q_y,
        frame.Q_yZ,
    )

    dA = dq + h * Z * dq_z / 2 + dq_y
    dA_Z = (
        dq_z
        + (h_Z * Z * dq_z + h * dq_z + h * Z * dq_zz) / 2
        + dq_yz
    )

    # Delta B expanded around the frame jet.  The last two terms are the
    # signed quadratic contributions that would be lost by a raw subtraction.
    dB = dq * (q + q_y) + q * dq_y + dq * dq / 2 + dq * dq_y
    dB_Z = (
        dq_z * (q + q_y)
        + dq * (q_z + q_yz)
        + q_z * dq_y
        + q * dq_yz
        + dq * dq_z
        + dq_z * dq_y
        + dq * dq_yz
    )

    dC = -2 * Z * dq + (1 - Z * Z) * dq_z - 2 * Z * dq_y
    dC_Z = (
        -2 * dq
        - 4 * Z * dq_z
        + (1 - Z * Z) * dq_zz
        - 2 * dq_y
        - 2 * Z * dq_yz
    )

    # Delta(U C), including the signed quadratic DeltaU*DeltaC term and its
    # ordinary-Z derivative.
    dUC = frame.U * dC + signed.dU * frame_terms.C + signed.dU * dC
    dUC_Z = (
        frame.U_Z * dC
        + frame.U * dC_Z
        + signed.dU_Z * frame_terms.C
        + signed.dU * frame_terms.C_Z
        + signed.dU_Z * dC
        + signed.dU * dC_Z
    )

    dA_plus_dUC = dA + dUC
    dA_plus_dUC_Z = dA_Z + dUC_Z
    return _DifferenceTerms(
        dA,
        dA_Z,
        dB,
        dB_Z,
        dC,
        dC_Z,
        dUC,
        dUC_Z,
        dA_plus_dUC,
        dA_plus_dUC_Z,
    )


def evaluate_signed_interval_algebra(
    frame_jets: FrameJets | Mapping[str, Any] | Any,
    signed_delta_jets: SignedDeltaJets | Mapping[str, Any] | Any,
    *,
    Z: Any,
    delta: Any,
    delta_Z: Any,
    denominator_nonzero: Callable[[Any], bool],
) -> SignedIntervalAlgebraResult:
    """Evaluate signed density and ordinary-``Z`` algebra in caller arithmetic.

    Parameters
    ----------
    frame_jets:
        Seven frame components ``U, U_Z, Q, Q_Z, Q_ZZ, Q_y, Q_yZ``.
    signed_delta_jets:
        Seven *direct* signed ``P-F`` components ``dU, dU_Z, dQ, dQ_Z,
        dQ_ZZ, dQ_y, dQ_yZ``.  They are never reconstructed here.
    Z, delta, delta_Z:
        Caller-owned arithmetic objects.  Built-in integers are represented
        internally as exact ``Fraction`` values; built-in floats and bools are
        rejected.  Existing interval/rational objects are not coerced.
        ``delta_Z`` is the ordinary-``Z`` derivative used in ``h_Z`` and
        ``L_Z``.
    denominator_nonzero:
        A fail-closed callback receiving the computed ``L``.  Division is
        allowed only when the callback returns the built-in ``True`` object
        (``verdict is True``); ``numpy.bool_`` and other truthy objects are
        rejected.

    Returns
    -------
    SignedIntervalAlgebraResult
        ``S_F``, ``S_F_Z``, direct signed ``dS``, direct signed ``dS_Z``, and
        the shared denominator ``L``.  All values retain the caller's
        arithmetic type.

    Notes
    -----
    This function does not prove the denominator callback, source identity,
    physical units, endpoint ownership, quadrature correctness, or any error
    bound.  It only performs the finite algebra after the caller's explicit
    denominator gate.
    """

    if not callable(denominator_nonzero):
        raise TypeError("denominator_nonzero must be callable")

    frame = _frame_jets(frame_jets)
    signed = _signed_delta_jets(signed_delta_jets)
    Z = _numeric_input(Z, "Z")
    delta = _numeric_input(delta, "delta")
    delta_Z = _numeric_input(delta_Z, "delta_Z")

    h = 1 - delta
    h_Z = -delta_Z
    L = 1 - delta * Z * Z
    L_Z = -delta_Z * Z * Z - 2 * delta * Z

    # Build every numerator before consulting the denominator gate.  This
    # keeps the gate the only permission boundary for division by L.
    frame_terms = _frame_terms(frame, Z=Z, h=h, h_Z=h_Z)
    difference_terms = _difference_terms(
        frame,
        signed,
        frame_terms,
        Z=Z,
        h=h,
        h_Z=h_Z,
    )

    verdict = denominator_nonzero(L)
    if verdict is not True:
        raise ValueError(
            "denominator_nonzero(L) must return the built-in True; "
            "denominator permission was not established"
        )

    S_F = frame_terms.numerator / L + frame_terms.B
    S_F_Z = (
        frame_terms.numerator_Z / L
        + frame_terms.B_Z
        - frame_terms.numerator * L_Z / (L * L)
    )
    dS = difference_terms.numerator / L + difference_terms.dB
    dS_Z = (
        difference_terms.numerator_Z / L
        + difference_terms.dB_Z
        - difference_terms.numerator * L_Z / (L * L)
    )
    return SignedIntervalAlgebraResult(S_F, S_F_Z, dS, dS_Z, L)


# Descriptive aliases for runners using either verb.
compute_signed_interval_algebra = evaluate_signed_interval_algebra
signed_interval_algebra = evaluate_signed_interval_algebra


__all__ = [
    "FIELD_NAMES",
    "SIGNED_FIELD_NAMES",
    "FrameJets",
    "FieldJets",
    "SignedDeltaJets",
    "DeltaJets",
    "SignedIntervalAlgebraResult",
    "evaluate_signed_interval_algebra",
    "compute_signed_interval_algebra",
    "signed_interval_algebra",
]
