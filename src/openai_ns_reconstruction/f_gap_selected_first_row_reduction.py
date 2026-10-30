"""Caller-owned first-row reduction for the F-gap selected coefficients.

The source-level matrix algebra has the following shape.  ``A`` is the first
row of the finite-``mu`` moment matrix and ``D`` is its divided row.  The
positive denominator supplied by the caller is the absolute value of the
matrix determinant, so the signed determinant identity is

``A[0] * D[1] - A[1] * D[0] + positive_determinant = 0``.

If ``C1`` and ``C2`` are the selected coefficients and the first source rows
are represented by ``incoming_rhs`` and ``scaled_rhs``, the two row equations
are checked by the residual interface in :func:`matrix_algebra_residuals`.
The first row then reduces to

``A[0] * C1 + A[1] * C2``
``= -mu * (incoming_factor1 * incoming_m1``
``       + ap * scaled_main_rows1)``.

This module implements only caller-owned ordinary-``Z`` jet arithmetic through
order three.  It does not import a producer, inspect a graph, evaluate a
source function, prove source identity, prove determinant positivity, select a
value, or certify an error.  ``Fraction`` values and non-primitive
caller-directed interval objects are retained by the arithmetic operators.
Plain integers and floats are rejected so an opaque node id or an accidental
binary-rounded value cannot silently become a source input.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import comb
from typing import Any, Generic, Mapping, Sequence, TypeVar


Scalar = TypeVar("Scalar")
SCHEMA = "f_gap_selected_first_row_reduction_v1"
ORDERS = (0, 1, 2, 3)


def _reject_primitive(name: str, value: Any) -> None:
    """Reject primitive values that could hide an id or force rounding."""

    if isinstance(value, bool) or isinstance(value, int):
        raise TypeError(
            f"{name} must be Fraction or a caller-directed interval object; "
            "plain integers are not accepted as source values or node ids"
        )
    if isinstance(value, float):
        raise TypeError(
            f"{name} must not be float; pass Fraction or an interval object"
        )


@dataclass(frozen=True, slots=True)
class OrdinaryZJet3(Generic[Scalar]):
    """Four explicit ordinary derivatives at fixed ``xi``.

    The rows are ordinary derivatives, not divided Taylor coefficients.  No
    omitted row is defaulted to zero.
    """

    rows: tuple[Scalar, Scalar, Scalar, Scalar]

    def __post_init__(self) -> None:
        raw = tuple(self.rows)
        if len(raw) != 4:
            raise ValueError("OrdinaryZJet3 requires exactly four rows")
        for order, value in enumerate(raw):
            _reject_primitive(f"jet row {order}", value)
        object.__setattr__(self, "rows", raw)

    @property
    def z0(self) -> Scalar:
        return self.rows[0]

    @property
    def z1(self) -> Scalar:
        return self.rows[1]

    @property
    def z2(self) -> Scalar:
        return self.rows[2]

    @property
    def z3(self) -> Scalar:
        return self.rows[3]


JetLike = OrdinaryZJet3[Scalar] | Sequence[Scalar]


def make_ordinary_z_jet3(
    z0: Scalar, z1: Scalar, z2: Scalar, z3: Scalar
) -> OrdinaryZJet3[Scalar]:
    """Construct a four-row ordinary-Z jet without implicit defaults."""

    return OrdinaryZJet3((z0, z1, z2, z3))


def _coerce_jet(value: JetLike[Scalar], name: str) -> OrdinaryZJet3[Scalar]:
    if isinstance(value, OrdinaryZJet3):
        return value
    if isinstance(value, (str, bytes)):
        raise TypeError(f"{name} must provide four ordinary-Z rows")
    try:
        rows = tuple(value)
    except TypeError as exc:
        raise TypeError(f"{name} must provide four ordinary-Z rows") from exc
    try:
        return OrdinaryZJet3(rows)  # type: ignore[arg-type]
    except (TypeError, ValueError) as exc:
        raise type(exc)(f"{name}: {exc}") from exc


def _binomial(order: int, split: int) -> Fraction:
    return Fraction(comb(order, split), 1)


def _add(left: Scalar, right: Scalar) -> Scalar:
    return left + right  # type: ignore[operator,return-value]


def _sub(left: Scalar, right: Scalar) -> Scalar:
    return left - right  # type: ignore[operator,return-value]


def _neg(value: Scalar) -> Scalar:
    return -value  # type: ignore[operator,return-value]


def _mul(left: Scalar, right: Scalar) -> Scalar:
    return left * right  # type: ignore[operator,return-value]


def _sum(terms: Sequence[Scalar]) -> Scalar:
    if not terms:
        raise ValueError("at least one term is required")
    total = terms[0]
    for term in terms[1:]:
        total = _add(total, term)
    return total


def _jet_add(*jets: OrdinaryZJet3[Scalar]) -> OrdinaryZJet3[Scalar]:
    if not jets:
        raise ValueError("jet addition requires at least one operand")
    return OrdinaryZJet3(
        tuple(_sum([jet.rows[order] for jet in jets]) for order in ORDERS)  # type: ignore[arg-type]
    )


def _jet_sub(
    left: OrdinaryZJet3[Scalar], right: OrdinaryZJet3[Scalar]
) -> OrdinaryZJet3[Scalar]:
    return OrdinaryZJet3(
        tuple(_sub(left.rows[order], right.rows[order]) for order in ORDERS)  # type: ignore[arg-type]
    )


def _jet_neg(jet: OrdinaryZJet3[Scalar]) -> OrdinaryZJet3[Scalar]:
    return OrdinaryZJet3(tuple(_neg(value) for value in jet.rows))  # type: ignore[arg-type]


def _jet_product(
    left: OrdinaryZJet3[Scalar], right: OrdinaryZJet3[Scalar]
) -> OrdinaryZJet3[Scalar]:
    """Multiply ordinary derivative jets with binomial Leibniz factors."""

    rows: list[Scalar] = []
    for order in ORDERS:
        rows.append(
            _sum(
                [
                    _mul(
                        _mul(_binomial(order, split), left.rows[split]),
                        right.rows[order - split],
                    )
                    for split in range(order + 1)
                ]
            )
        )
    return OrdinaryZJet3(tuple(rows))  # type: ignore[arg-type]


def _jet_scale(jet: OrdinaryZJet3[Scalar], scalar: Scalar) -> OrdinaryZJet3[Scalar]:
    return OrdinaryZJet3(
        tuple(_mul(scalar, value) for value in jet.rows)  # type: ignore[arg-type]
    )


@dataclass(frozen=True, slots=True)
class FGapSelectedFirstRowReduction(Generic[Scalar]):
    """Uncertified caller-owned reduced first-row jet."""

    reduced: OrdinaryZJet3[Scalar]
    incoming_product: OrdinaryZJet3[Scalar]
    selected_product: OrdinaryZJet3[Scalar]
    coordinate: str = "fixed_xi"
    schema: str = SCHEMA
    source_identity_proved: bool = False
    source_values_authenticated: bool = False
    determinant_positivity_proved: bool = False
    error_certified: bool = False

    def as_dict(self) -> dict[str, Any]:
        """Return rows and explicit non-claims without coercing values."""

        return {
            "schema": self.schema,
            "coordinate": self.coordinate,
            "reduced_rows": self.reduced.rows,
            "incoming_product_rows": self.incoming_product.rows,
            "selected_product_rows": self.selected_product.rows,
            "source_identity_proved": self.source_identity_proved,
            "source_values_authenticated": self.source_values_authenticated,
            "determinant_positivity_proved": self.determinant_positivity_proved,
            "error_certified": self.error_certified,
        }


def reduce_f_gap_selected_first_row(
    incoming_factor1: JetLike[Scalar],
    incoming_m1: JetLike[Scalar],
    ap: JetLike[Scalar],
    scaled_main_rows1: JetLike[Scalar],
    mu: JetLike[Scalar],
) -> FGapSelectedFirstRowReduction[Scalar]:
    """Return ``-mu*(incoming_factor1*incoming_m1 + ap*scaled_main_rows1)``.

    All five caller inputs require ordinary-Z rows 0 through 3.  The returned
    jet is an algebraic reduction only; no source owner, source value,
    determinant positivity, numerical interval, or error certificate is
    inferred by this function.
    """

    factor = _coerce_jet(incoming_factor1, "incoming_factor1")
    incoming = _coerce_jet(incoming_m1, "incoming_m1")
    amplitude = _coerce_jet(ap, "ap")
    scaled = _coerce_jet(scaled_main_rows1, "scaled_main_rows1")
    parameter = _coerce_jet(mu, "mu")

    incoming_product = _jet_product(factor, incoming)
    selected_product = _jet_product(amplitude, scaled)
    reduced = _jet_neg(_jet_product(parameter, _jet_add(incoming_product, selected_product)))
    return FGapSelectedFirstRowReduction(
        reduced=reduced,
        incoming_product=incoming_product,
        selected_product=selected_product,
    )


@dataclass(frozen=True, slots=True)
class FGapSelectedFirstRowAlgebraResiduals(Generic[Scalar]):
    """Termwise residual jets for the 2x2 matrix reduction.

    A zero residual is not asserted here.  The caller decides whether its
    supplied rows are same-source rows and how to interpret any resulting
    interval or exact residual.
    """

    determinant_residual: OrdinaryZJet3[Scalar]
    incoming_first_rhs_residual: OrdinaryZJet3[Scalar]
    scaled_first_rhs_residual: OrdinaryZJet3[Scalar]
    matrix_first_row_residual: OrdinaryZJet3[Scalar]
    matrix_divided_row_residual: OrdinaryZJet3[Scalar]
    raw_first_row_reduction_residual: OrdinaryZJet3[Scalar]
    selected_first_row: OrdinaryZJet3[Scalar]
    matrix_first_row_target: OrdinaryZJet3[Scalar]
    raw_first_row_target: OrdinaryZJet3[Scalar]
    schema: str = SCHEMA
    source_identity_proved: bool = False
    source_values_authenticated: bool = False
    determinant_positivity_proved: bool = False
    error_certified: bool = False

    def as_dict(self) -> dict[str, Any]:
        """Return residual rows and explicit non-claims without coercion."""

        return {
            "schema": self.schema,
            "determinant_residual": self.determinant_residual.rows,
            "incoming_first_rhs_residual": self.incoming_first_rhs_residual.rows,
            "scaled_first_rhs_residual": self.scaled_first_rhs_residual.rows,
            "matrix_first_row_residual": self.matrix_first_row_residual.rows,
            "matrix_divided_row_residual": self.matrix_divided_row_residual.rows,
            "raw_first_row_reduction_residual": self.raw_first_row_reduction_residual.rows,
            "selected_first_row": self.selected_first_row.rows,
            "matrix_first_row_target": self.matrix_first_row_target.rows,
            "raw_first_row_target": self.raw_first_row_target.rows,
            "source_identity_proved": self.source_identity_proved,
            "source_values_authenticated": self.source_values_authenticated,
            "determinant_positivity_proved": self.determinant_positivity_proved,
            "error_certified": self.error_certified,
        }


def matrix_algebra_residuals(
    *,
    A: tuple[JetLike[Scalar], JetLike[Scalar]],
    D: tuple[JetLike[Scalar], JetLike[Scalar]],
    positive_determinant: JetLike[Scalar],
    incoming_rhs: tuple[JetLike[Scalar], JetLike[Scalar]],
    scaled_rhs: tuple[JetLike[Scalar], JetLike[Scalar]],
    ap: JetLike[Scalar],
    selected_C1: JetLike[Scalar],
    selected_C2: JetLike[Scalar],
    incoming_factor1: JetLike[Scalar],
    incoming_m1: JetLike[Scalar],
    scaled_main_rows1: JetLike[Scalar],
    mu: JetLike[Scalar],
) -> FGapSelectedFirstRowAlgebraResiduals[Scalar]:
    """Return termwise residuals for the same-source 2x2 algebra.

    The caller supplies the two right-hand-side rows used by the inverse:

    ``incoming_rhs = (-mu*I1, -(I2-I1))`` and
    ``scaled_rhs = (-mu*S1, -(S2-S1))``.

    The first row of ``A*C`` is compared with
    ``incoming_rhs[0] + ap*scaled_rhs[0]``.  Separately, the raw five-input
    reduction is compared with that matrix target.  The positive determinant
    is a caller-provided row and is used only in the signed determinant
    residual; this interface does not prove it positive or nonzero.
    """

    a0, a1 = (_coerce_jet(value, f"A[{index}]") for index, value in enumerate(A))
    d0, d1 = (_coerce_jet(value, f"D[{index}]") for index, value in enumerate(D))
    determinant = _coerce_jet(positive_determinant, "positive_determinant")
    incoming0, incoming1 = (
        _coerce_jet(value, f"incoming_rhs[{index}]")
        for index, value in enumerate(incoming_rhs)
    )
    scaled0, scaled1 = (
        _coerce_jet(value, f"scaled_rhs[{index}]")
        for index, value in enumerate(scaled_rhs)
    )
    amplitude = _coerce_jet(ap, "ap")
    coefficient1 = _coerce_jet(selected_C1, "selected_C1")
    coefficient2 = _coerce_jet(selected_C2, "selected_C2")

    selected_first = _jet_add(_jet_product(a0, coefficient1), _jet_product(a1, coefficient2))
    selected_divided = _jet_add(_jet_product(d0, coefficient1), _jet_product(d1, coefficient2))
    matrix_first_target = _jet_add(incoming0, _jet_product(amplitude, scaled0))
    matrix_divided_target = _jet_add(incoming1, _jet_product(amplitude, scaled1))

    raw_factor = _coerce_jet(incoming_factor1, "incoming_factor1")
    raw_incoming = _coerce_jet(incoming_m1, "incoming_m1")
    raw_scaled = _coerce_jet(scaled_main_rows1, "scaled_main_rows1")
    raw_mu = _coerce_jet(mu, "mu")
    raw_target = _jet_neg(
        _jet_product(
            raw_mu,
            _jet_add(_jet_product(raw_factor, raw_incoming), _jet_product(amplitude, raw_scaled)),
        )
    )

    return FGapSelectedFirstRowAlgebraResiduals(
        determinant_residual=_jet_add(
            _jet_sub(_jet_product(a0, d1), _jet_product(a1, d0)), determinant
        ),
        incoming_first_rhs_residual=_jet_sub(
            incoming0,
            _jet_neg(_jet_product(raw_mu, _jet_product(raw_factor, raw_incoming))),
        ),
        scaled_first_rhs_residual=_jet_sub(
            scaled0,
            _jet_neg(_jet_product(raw_mu, raw_scaled)),
        ),
        matrix_first_row_residual=_jet_sub(selected_first, matrix_first_target),
        matrix_divided_row_residual=_jet_sub(selected_divided, matrix_divided_target),
        raw_first_row_reduction_residual=_jet_sub(matrix_first_target, raw_target),
        selected_first_row=selected_first,
        matrix_first_row_target=matrix_first_target,
        raw_first_row_target=raw_target,
    )


__all__ = [
    "FGapSelectedFirstRowAlgebraResiduals",
    "FGapSelectedFirstRowReduction",
    "OrdinaryZJet3",
    "SCHEMA",
    "make_ordinary_z_jet3",
    "matrix_algebra_residuals",
    "reduce_f_gap_selected_first_row",
]
