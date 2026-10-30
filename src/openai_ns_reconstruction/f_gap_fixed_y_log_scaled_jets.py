"""Conditional fixed-``y`` common-log-scale ordinary-``Z`` F-gap jets.

This module is a small, caller-owned algebra leaf.  It combines the gap mean
source form

``m = A * exp(phi)``,
``A = -Pstar*u*(C1*W11 + C2*W12)``,
``phi = log_end_scale + 13/(2*mu) - 13 - y``

with the radial pressure formula after imposing the caller's gap conditions
``U = U_y = 0``, ``m_y = -m``, and ``delta_y = 0``:

``Q = (-(1-delta)*Z*m - (1-Z**2)*m_Z) / (1-delta*Z**2)``.

The coordinate is explicitly fixed ``y``.  The returned rows are ordinary
derivatives with respect to ``Z`` (not divided Taylor coefficients).  The
common factor ``exp(log_scale)`` is never evaluated: ``scaled_m`` and
``scaled_Q`` are the corresponding rows with that factor removed, and
``scaled_Qy = -scaled_Q`` through order one.  The Bell factors for the
ordinary derivatives of ``exp(phi)`` are included through order three.

For the optional gap contraction ``R = Rp*exp(y)``, the same unmaterialized
rows expose ``log_RQ_scale = log_scale + y``,
``scaled_RQ_y = -scaled_Q0``, and ``scaled_RQ_y_Z = -scaled_Q1``.  Holding
the other supplied rows fixed, the ``+y`` in ``log_RQ_scale`` cancels the
``-y`` in ``log_scale``, so two distinct fixed-``y`` choices have the same
conditional ``R*Q_y/Rp`` value.  Those fields are conditional on fixed ``y``
and ``Rp_Z = 0``; no physical ``Rp`` owner or error-term deletion is
asserted.  Rows for nonzero ``m`` and ``Q`` are retained.

``y`` is a caller premise and is held fixed in the ordinary ``Z`` derivative.
If one substitutes ``xi = mu*y`` into a fixed-``xi`` formula, that does not
give the same derivative convention: fixed ``xi`` differentiates the
``xi/mu`` term while fixed ``y`` differentiates the constant ``y`` term.

All data are either exact :class:`fractions.Fraction` values or opaque,
caller-directed interval objects.  Plain integers and binary floats are
rejected as data.  No source module, graph, producer result, numerical
midpoint, or ``exp`` evaluation is used.  Both the ``mu`` reciprocal and the
radial-Q denominator require caller nonzero gates; these gates are premises,
not proofs.  Source identity, physical-scale certification, and numerical
certification therefore remain explicitly false in the result flags.  The
``delta_Z`` and ``delta_ZZ`` rows are accepted as a generalization; whether
the physical source used a fixed ``delta`` (and whether ``mu`` or other
geometry is source-independent) is a separate, unauthenticated premise.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import comb
from typing import Any, Callable, Generic, Mapping, Sequence, TypeVar


Scalar = TypeVar("Scalar")

SCHEMA = "f_gap_fixed_y_log_scaled_jets_v1"
FIXED_COORDINATE = "fixed_y"
ROW_CONVENTION = "ordinary Z derivatives, not divided Taylor coefficients"
SOURCE_IDENTITY_PROVED = False
PHYSICAL_SCALE_PROVED = False
CERTIFIED = False

_ZERO = Fraction(0, 1)
_ONE = Fraction(1, 1)


def _reject_primitive(name: str, value: Any) -> None:
    """Reject scalar primitives that are not admitted by this algebra leaf."""

    if isinstance(value, bool) or isinstance(value, int):
        raise TypeError(
            f"{name} must be Fraction or a caller-directed interval object; "
            "plain integers are not accepted as data"
        )
    if isinstance(value, float):
        raise TypeError(f"{name} must not be float; pass Fraction or an interval object")


@dataclass(frozen=True, slots=True)
class OrdinaryZJet3(Generic[Scalar]):
    """Four explicit ordinary derivatives in the fixed-``y`` coordinate."""

    rows: tuple[Scalar, Scalar, Scalar, Scalar]

    def __post_init__(self) -> None:
        raw = tuple(self.rows)
        if len(raw) != 4:
            raise ValueError("OrdinaryZJet3 requires exactly four explicit rows")
        for order, value in enumerate(raw):
            _reject_primitive(f"jet row {order}", value)
        object.__setattr__(self, "rows", raw)

    def __iter__(self):
        return iter(self.rows)

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


@dataclass(frozen=True, slots=True)
class OrdinaryZJet2(Generic[Scalar]):
    """Three explicit ordinary derivatives in the fixed-``y`` coordinate."""

    rows: tuple[Scalar, Scalar, Scalar]

    def __post_init__(self) -> None:
        raw = tuple(self.rows)
        if len(raw) != 3:
            raise ValueError("OrdinaryZJet2 requires exactly three explicit rows")
        for order, value in enumerate(raw):
            _reject_primitive(f"jet row {order}", value)
        object.__setattr__(self, "rows", raw)

    def __iter__(self):
        return iter(self.rows)

    @property
    def z0(self) -> Scalar:
        return self.rows[0]

    @property
    def z1(self) -> Scalar:
        return self.rows[1]

    @property
    def z2(self) -> Scalar:
        return self.rows[2]


def make_ordinary_z_jet3(
    z0: Scalar, z1: Scalar, z2: Scalar, z3: Scalar
) -> OrdinaryZJet3[Scalar]:
    """Construct four explicit ordinary-``Z`` rows."""

    return OrdinaryZJet3((z0, z1, z2, z3))


def make_ordinary_z_jet2(z0: Scalar, z1: Scalar, z2: Scalar) -> OrdinaryZJet2[Scalar]:
    """Construct three explicit ordinary-``Z`` rows."""

    return OrdinaryZJet2((z0, z1, z2))


def _coerce_jet3(value: OrdinaryZJet3[Scalar] | Sequence[Scalar], name: str) -> OrdinaryZJet3[Scalar]:
    """Coerce a four-row sequence without interpreting row values."""

    if isinstance(value, OrdinaryZJet3):
        return value
    if isinstance(value, (str, bytes)):
        raise TypeError(f"{name} must provide four ordinary Z rows")
    try:
        rows = tuple(value)
    except TypeError:
        try:
            rows = tuple(getattr(value, "rows"))
        except (AttributeError, TypeError) as exc:
            raise TypeError(f"{name} must provide four ordinary Z rows") from exc
    try:
        return OrdinaryZJet3(rows)  # type: ignore[arg-type]
    except (TypeError, ValueError) as exc:
        raise type(exc)(f"{name}: {exc}") from exc


def _coerce_jet2(value: OrdinaryZJet2[Scalar] | Sequence[Scalar], name: str) -> OrdinaryZJet2[Scalar]:
    """Coerce a three-row sequence without changing its row convention."""

    if isinstance(value, OrdinaryZJet2):
        return value
    if isinstance(value, (str, bytes)):
        raise TypeError(f"{name} must provide three ordinary Z rows")
    try:
        rows = tuple(value)
    except TypeError:
        try:
            rows = tuple(getattr(value, "rows"))
        except (AttributeError, TypeError) as exc:
            raise TypeError(f"{name} must provide three ordinary Z rows") from exc
    try:
        return OrdinaryZJet2(rows)  # type: ignore[arg-type]
    except (TypeError, ValueError) as exc:
        raise type(exc)(f"{name}: {exc}") from exc


def _binomial(n: int, k: int) -> Fraction:
    return Fraction(comb(n, k), 1)


def _add(left: Scalar, right: Scalar) -> Scalar:
    return left + right  # type: ignore[operator,return-value]


def _sub(left: Scalar, right: Scalar) -> Scalar:
    return left - right  # type: ignore[operator,return-value]


def _neg(value: Scalar) -> Scalar:
    return -value  # type: ignore[operator,return-value]


def _mul(left: Scalar, right: Scalar) -> Scalar:
    return left * right  # type: ignore[operator,return-value]


def _div(left: Scalar, right: Scalar) -> Scalar:
    return left / right  # type: ignore[operator,return-value]


def _sum(terms: Sequence[Scalar]) -> Scalar:
    if not terms:
        raise ValueError("at least one ordinary-jet term is required")
    total = terms[0]
    for term in terms[1:]:
        total = _add(total, term)
    return total


def _constant3(value: Scalar) -> OrdinaryZJet3[Scalar]:
    _reject_primitive("fixed scalar", value)
    return OrdinaryZJet3((value, _ZERO, _ZERO, _ZERO))


def _constant2(value: Scalar) -> OrdinaryZJet2[Scalar]:
    _reject_primitive("fixed scalar", value)
    return OrdinaryZJet2((value, _ZERO, _ZERO))


def _jet3_add(*jets: OrdinaryZJet3[Scalar]) -> OrdinaryZJet3[Scalar]:
    if not jets:
        raise ValueError("ordinary-jet addition requires at least one operand")
    return OrdinaryZJet3(
        tuple(_sum([jet.rows[order] for jet in jets]) for order in range(4))  # type: ignore[arg-type]
    )


def _jet3_sub(left: OrdinaryZJet3[Scalar], right: OrdinaryZJet3[Scalar]) -> OrdinaryZJet3[Scalar]:
    return OrdinaryZJet3(
        tuple(_sub(left.rows[order], right.rows[order]) for order in range(4))  # type: ignore[arg-type]
    )


def _jet3_neg(jet: OrdinaryZJet3[Scalar]) -> OrdinaryZJet3[Scalar]:
    return OrdinaryZJet3(tuple(_neg(value) for value in jet.rows))  # type: ignore[arg-type]


def _jet3_product(left: OrdinaryZJet3[Scalar], right: OrdinaryZJet3[Scalar]) -> OrdinaryZJet3[Scalar]:
    """Multiply ordinary rows with the binomial Leibniz factors."""

    rows: list[Scalar] = []
    for order in range(4):
        terms = [
            _mul(_mul(_binomial(order, split), left.rows[split]), right.rows[order - split])
            for split in range(order + 1)
        ]
        rows.append(_sum(terms))
    return OrdinaryZJet3(tuple(rows))  # type: ignore[arg-type]


def _jet3_reciprocal(jet: OrdinaryZJet3[Scalar]) -> OrdinaryZJet3[Scalar]:
    """Differentiate ``jet * reciprocal = 1`` through ordinary order three."""

    rows: list[Scalar] = [_div(_ONE, jet.rows[0])]
    for order in range(1, 4):
        terms = [
            _mul(_mul(_binomial(order, split), jet.rows[split]), rows[order - split])
            for split in range(1, order + 1)
        ]
        rows.append(_div(_neg(_sum(terms)), jet.rows[0]))
    return OrdinaryZJet3(tuple(rows))  # type: ignore[arg-type]


def _jet2_add(*jets: OrdinaryZJet2[Scalar]) -> OrdinaryZJet2[Scalar]:
    if not jets:
        raise ValueError("ordinary-jet addition requires at least one operand")
    return OrdinaryZJet2(
        tuple(_sum([jet.rows[order] for jet in jets]) for order in range(3))  # type: ignore[arg-type]
    )


def _jet2_sub(left: OrdinaryZJet2[Scalar], right: OrdinaryZJet2[Scalar]) -> OrdinaryZJet2[Scalar]:
    return OrdinaryZJet2(
        tuple(_sub(left.rows[order], right.rows[order]) for order in range(3))  # type: ignore[arg-type]
    )


def _jet2_neg(jet: OrdinaryZJet2[Scalar]) -> OrdinaryZJet2[Scalar]:
    return OrdinaryZJet2(tuple(_neg(value) for value in jet.rows))  # type: ignore[arg-type]


def _jet2_product(left: OrdinaryZJet2[Scalar], right: OrdinaryZJet2[Scalar]) -> OrdinaryZJet2[Scalar]:
    """Multiply three ordinary rows with the binomial Leibniz factors."""

    rows: list[Scalar] = []
    for order in range(3):
        terms = [
            _mul(_mul(_binomial(order, split), left.rows[split]), right.rows[order - split])
            for split in range(order + 1)
        ]
        rows.append(_sum(terms))
    return OrdinaryZJet2(tuple(rows))  # type: ignore[arg-type]


def _jet2_reciprocal(jet: OrdinaryZJet2[Scalar]) -> OrdinaryZJet2[Scalar]:
    """Differentiate ``jet * reciprocal = 1`` through ordinary order two."""

    rows: list[Scalar] = [_div(_ONE, jet.rows[0])]
    for order in range(1, 3):
        terms = [
            _mul(_mul(_binomial(order, split), jet.rows[split]), rows[order - split])
            for split in range(1, order + 1)
        ]
        rows.append(_div(_neg(_sum(terms)), jet.rows[0]))
    return OrdinaryZJet2(tuple(rows))  # type: ignore[arg-type]


def _differentiate3(jet: OrdinaryZJet3[Scalar]) -> OrdinaryZJet2[Scalar]:
    """Drop the value row and retain ordinary rows 1 through 3."""

    return OrdinaryZJet2((jet.rows[1], jet.rows[2], jet.rows[3]))


@dataclass(frozen=True, slots=True)
class FGapFixedYLogScaledJets(Generic[Scalar]):
    """Unmaterialized fixed-``y`` mean and radial-Q jet result.

    ``scaled_m`` and ``scaled_Q`` retain ordinary derivative rows.  The
    physical rows are represented by multiplying those rows by
    ``exp(log_scale)`` outside this module.
    """

    log_scale: Scalar
    scaled_m: OrdinaryZJet3[Scalar]
    scaled_Q: OrdinaryZJet2[Scalar]
    scaled_Qy: tuple[Scalar, Scalar]
    phi: OrdinaryZJet3[Scalar]
    amplitude: OrdinaryZJet3[Scalar]
    fixed_y: Scalar
    Z: Scalar
    delta: OrdinaryZJet2[Scalar]
    # Conditional gap identity for (R*Q_y)/Rp with R = Rp*exp(y).  These
    # fields are algebraic data only; the Rp_Z=0 premise is recorded below.
    log_RQ_scale: Scalar
    scaled_RQ_y: Scalar
    scaled_RQ_y_Z: Scalar
    coordinate: str = FIXED_COORDINATE
    schema: str = SCHEMA
    row_convention: str = ROW_CONVENTION
    source_formula_conditioned: bool = True
    fixed_y_premise: bool = True
    gap_U_zero_conditioned: bool = True
    gap_Uy_zero_conditioned: bool = True
    my_identity_conditioned: bool = True
    delta_y_zero_conditioned: bool = True
    Rp_Z_zero_premise: bool = True
    source_parameter_independence_proved: bool = False
    source_mu_independence_proved: bool = False
    source_delta_independence_proved: bool = False
    source_Rp_independence_proved: bool = False
    source_delta_fixed_premise: bool = False
    delta_Z_rows_are_generalization: bool = True
    physical_source_owner_proved: bool = False
    source_Rp_owner_authenticated: bool = False
    error_terms_deletable: bool = False
    rq_y_independent_of_fixed_y: bool = True
    two_distinct_y_rqy_equal_conditioned: bool = True
    two_y_value_equality_proved: bool = False
    nonzero_m_Q_rows_retained: bool = True
    pstar_applied_once: bool = True
    scale_materialized: bool = False
    physical_scale_materialized: bool = False
    source_identity_proved: bool = False
    source_identity: bool = False
    source_values_authenticated: bool = False
    physical_scale_proved: bool = False
    physical_scale_certified: bool = False
    certification: bool = False
    certified: bool = False

    @property
    def scaled_m_rows(self) -> tuple[Scalar, Scalar, Scalar, Scalar]:
        return self.scaled_m.rows

    @property
    def y(self) -> Scalar:
        return self.fixed_y

    @property
    def scaled_m0(self) -> Scalar:
        return self.scaled_m.rows[0]

    @property
    def scaled_m1(self) -> Scalar:
        return self.scaled_m.rows[1]

    @property
    def scaled_m2(self) -> Scalar:
        return self.scaled_m.rows[2]

    @property
    def scaled_m3(self) -> Scalar:
        return self.scaled_m.rows[3]

    @property
    def scaled_Q_rows(self) -> tuple[Scalar, Scalar, Scalar]:
        return self.scaled_Q.rows

    @property
    def scaled_Q0(self) -> Scalar:
        return self.scaled_Q.rows[0]

    @property
    def scaled_Q1(self) -> Scalar:
        return self.scaled_Q.rows[1]

    @property
    def scaled_Q2(self) -> Scalar:
        return self.scaled_Q.rows[2]

    @property
    def scaled_q0(self) -> Scalar:
        return self.scaled_Q.rows[0]

    @property
    def scaled_q1(self) -> Scalar:
        return self.scaled_Q.rows[1]

    @property
    def scaled_q2(self) -> Scalar:
        return self.scaled_Q.rows[2]

    @property
    def scaled_Qy_rows(self) -> tuple[Scalar, Scalar]:
        return self.scaled_Qy

    @property
    def scaled_q_rows(self) -> tuple[Scalar, Scalar, Scalar]:
        return self.scaled_Q.rows

    @property
    def scaled_qy_rows(self) -> tuple[Scalar, Scalar]:
        return self.scaled_Qy

    @property
    def scaled_q(self) -> OrdinaryZJet2[Scalar]:
        """Lower-case alias for callers that use ``q`` naming."""

        return self.scaled_Q

    @property
    def scaled_qy(self) -> tuple[Scalar, Scalar]:
        """Lower-case alias for callers that use ``qy`` naming."""

        return self.scaled_Qy

    @property
    def scaled_Qy0(self) -> Scalar:
        return self.scaled_Qy[0]

    @property
    def scaled_Qy1(self) -> Scalar:
        return self.scaled_Qy[1]

    @property
    def scaled_Q_y0(self) -> Scalar:
        return self.scaled_Qy[0]

    @property
    def scaled_Q_y1(self) -> Scalar:
        return self.scaled_Qy[1]

    @property
    def scaled_Qy_jet(self) -> OrdinaryZJet2[Scalar]:
        """Three-row compatibility view with the conditioned ``Qy=-Q`` rows."""

        return OrdinaryZJet2(
            (self.scaled_Qy[0], self.scaled_Qy[1], _neg(self.scaled_Q.rows[2]))
        )

    @property
    def scaled_RQ_y0(self) -> Scalar:
        return self.scaled_RQ_y

    @property
    def scaled_RQ_y1(self) -> Scalar:
        return self.scaled_RQ_y_Z

    @property
    def scaled_RQ_y_Z1(self) -> Scalar:
        return self.scaled_RQ_y_Z

    @property
    def scaled_RQ_y_rows(self) -> tuple[Scalar, Scalar]:
        return (self.scaled_RQ_y, self.scaled_RQ_y_Z)

    @property
    def scaled_RQy0(self) -> Scalar:
        return self.scaled_RQ_y

    @property
    def scaled_RQy1(self) -> Scalar:
        return self.scaled_RQ_y_Z

    @property
    def rp_Z_zero_premise(self) -> bool:
        """Lower-case alias for the explicit ``Rp_Z=0`` premise flag."""

        return self.Rp_Z_zero_premise

    @property
    def physical_scale(self) -> bool:
        """False until a caller supplies a physical-scale certification."""

        return self.physical_scale_proved

    @property
    def actual_m_statement(self) -> str:
        return "m^(k) = exp(log_scale) * scaled_m[k], k=0..3; exp is not materialized"

    @property
    def actual_Q_statement(self) -> str:
        return "Q^(k) = exp(log_scale) * scaled_Q[k], k=0..2; exp is not materialized"

    @property
    def actual_Qy_statement(self) -> str:
        return "Qy^(k) = exp(log_scale) * scaled_Qy[k], k=0..1; exp is not materialized"

    @property
    def actual_RQy_statement(self) -> str:
        return (
            "R*Qy/Rp = exp(log_RQ_scale) * scaled_RQ_y; "
            "fixed y and Rp_Z=0 are caller premises; the conditional value "
            "is unchanged between two fixed y choices; exp is not materialized"
        )

    @property
    def flags(self) -> Mapping[str, bool]:
        """Explicit nonclaims retained with the algebraic result."""

        return {
            "source_identity": self.source_identity,
            "source_identity_proved": self.source_identity_proved,
            "physical_scale": self.physical_scale,
            "physical_scale_proved": self.physical_scale_proved,
            "physical_scale_certified": self.physical_scale_certified,
            "certification": self.certification,
            "certified": self.certified,
            "source_values_authenticated": self.source_values_authenticated,
            "scale_materialized": self.scale_materialized,
            "physical_scale_materialized": self.physical_scale_materialized,
            "Rp_Z_zero_premise": self.Rp_Z_zero_premise,
            "source_parameter_independence_proved": self.source_parameter_independence_proved,
            "source_mu_independence_proved": self.source_mu_independence_proved,
            "source_delta_independence_proved": self.source_delta_independence_proved,
            "source_Rp_independence_proved": self.source_Rp_independence_proved,
            "source_delta_fixed_premise": self.source_delta_fixed_premise,
            "delta_Z_rows_are_generalization": self.delta_Z_rows_are_generalization,
            "physical_source_owner_proved": self.physical_source_owner_proved,
            "source_Rp_owner_authenticated": self.source_Rp_owner_authenticated,
            "error_terms_deletable": self.error_terms_deletable,
            "rq_y_independent_of_fixed_y": self.rq_y_independent_of_fixed_y,
            "two_distinct_y_rqy_equal_conditioned": self.two_distinct_y_rqy_equal_conditioned,
        }

    def as_tuple(self) -> tuple[Any, ...]:
        """Return the requested values in canonical output order."""

        return (
            self.log_scale,
            *self.scaled_m.rows,
            *self.scaled_Q.rows,
            *self.scaled_Qy,
        )

    def __iter__(self):
        return iter(self.as_tuple())

    def as_dict(self) -> dict[str, Any]:
        """Return rows and flags without evaluating the common exponential."""

        return {
            "schema": self.schema,
            "coordinate": self.coordinate,
            "fixed_y": self.fixed_y,
            "Z": self.Z,
            "delta_rows": self.delta.rows,
            "row_convention": self.row_convention,
            "phi_rows": self.phi.rows,
            "log_scale": self.log_scale,
            "amplitude_rows": self.amplitude.rows,
            "scaled_m_rows": self.scaled_m.rows,
            "scaled_Q_rows": self.scaled_Q.rows,
            "scaled_Qy_rows": self.scaled_Qy,
            "log_RQ_scale": self.log_RQ_scale,
            "scaled_RQ_y": self.scaled_RQ_y,
            "scaled_RQ_y_Z": self.scaled_RQ_y_Z,
            "actual_m": self.actual_m_statement,
            "actual_Q": self.actual_Q_statement,
            "actual_Qy": self.actual_Qy_statement,
            "actual_RQy": self.actual_RQy_statement,
            "source_formula_conditioned": self.source_formula_conditioned,
            "fixed_y_premise": self.fixed_y_premise,
            "gap_U_zero_conditioned": self.gap_U_zero_conditioned,
            "gap_Uy_zero_conditioned": self.gap_Uy_zero_conditioned,
            "my_identity_conditioned": self.my_identity_conditioned,
            "delta_y_zero_conditioned": self.delta_y_zero_conditioned,
            "Rp_Z_zero_premise": self.Rp_Z_zero_premise,
            "source_parameter_independence_proved": self.source_parameter_independence_proved,
            "source_mu_independence_proved": self.source_mu_independence_proved,
            "source_delta_independence_proved": self.source_delta_independence_proved,
            "source_Rp_independence_proved": self.source_Rp_independence_proved,
            "source_delta_fixed_premise": self.source_delta_fixed_premise,
            "delta_Z_rows_are_generalization": self.delta_Z_rows_are_generalization,
            "physical_source_owner_proved": self.physical_source_owner_proved,
            "source_Rp_owner_authenticated": self.source_Rp_owner_authenticated,
            "error_terms_deletable": self.error_terms_deletable,
            "rq_y_independent_of_fixed_y": self.rq_y_independent_of_fixed_y,
            "two_distinct_y_rqy_equal_conditioned": self.two_distinct_y_rqy_equal_conditioned,
            "two_y_value_equality_proved": self.two_y_value_equality_proved,
            "nonzero_m_Q_rows_retained": self.nonzero_m_Q_rows_retained,
            "pstar_applied_once": self.pstar_applied_once,
            "scale_materialized": self.scale_materialized,
            "physical_scale_materialized": self.physical_scale_materialized,
            "source_identity_proved": self.source_identity_proved,
            "source_identity": self.source_identity,
            "source_values_authenticated": self.source_values_authenticated,
            "physical_scale": self.physical_scale,
            "physical_scale_proved": self.physical_scale_proved,
            "physical_scale_certified": self.physical_scale_certified,
            "certification": self.certification,
            "certified": self.certified,
            "flags": dict(self.flags),
        }


def _select_gate(
    primary: Callable[[Scalar], bool] | None,
    secondary: Callable[[Scalar], bool] | None,
    fallback: Callable[[Scalar], bool] | None,
    name: str,
) -> Callable[[Scalar], bool]:
    gate = primary if primary is not None else secondary
    if gate is None:
        gate = fallback
    if gate is None or not callable(gate):
        raise TypeError(f"{name} must be callable")
    return gate


def build_f_gap_fixed_y_log_scaled_jets(
    Pstar: OrdinaryZJet3[Scalar] | Sequence[Scalar],
    u: OrdinaryZJet3[Scalar] | Sequence[Scalar],
    C1: OrdinaryZJet3[Scalar] | Sequence[Scalar],
    C2: OrdinaryZJet3[Scalar] | Sequence[Scalar],
    W11: OrdinaryZJet3[Scalar] | Sequence[Scalar],
    W12: OrdinaryZJet3[Scalar] | Sequence[Scalar],
    log_end_scale: OrdinaryZJet3[Scalar] | Sequence[Scalar],
    mu: OrdinaryZJet3[Scalar] | Sequence[Scalar],
    y: Scalar | None = None,
    Z: Scalar | None = None,
    delta: OrdinaryZJet2[Scalar] | Sequence[Scalar] | None = None,
    delta1: Scalar | None = None,
    delta2: Scalar | None = None,
    *,
    denominator_nonzero_gate: Callable[[Scalar], bool] | None = None,
    q_denominator_nonzero_gate: Callable[[Scalar], bool] | None = None,
    mu_nonzero_gate: Callable[[Scalar], bool] | None = None,
    nonzero_gate: Callable[[Scalar], bool] | None = None,
    fixed_y: Scalar | None = None,
    delta0: Scalar | None = None,
) -> FGapFixedYLogScaledJets[Scalar]:
    """Build conditional fixed-``y`` mean and radial-Q ordinary jets.

    The eight source arguments are explicit ordinary ``Z0`` through ``Z3``
    rows.  ``delta`` supplies ``delta``, ``delta_Z``, and ``delta_ZZ`` (or
    pass the three rows explicitly as ``delta, delta1, delta2``).  The
    scalar ``Z`` is the coordinate value itself, whose coordinate jet is
    ``(Z, 1, 0)``.  The three gap conditions are caller premises: this leaf
    does not independently construct ``U`` or ``U_y`` and sets ``Q_y=-Q``
    through ``Z1``.  Nonzero ``delta_Z``/``delta_ZZ`` rows are a formula
    generalization; a source-side ``delta``-fixed or geometry-independent
    claim is not proved here and is separately flagged false.

    ``y`` is fixed during ordinary ``Z`` differentiation.  In particular,
    ``phi`` is formed as ``log_end_scale + 13/(2*mu) - 13 - y``.  It must not
    be differentiated as a fixed-``xi`` expression after substituting
    ``xi=mu*y``; that fixed-``xi`` convention has different ``Z`` derivatives.

    A gate named ``mu_nonzero_gate`` (or the common ``denominator_nonzero_gate``
    fallback) must return built-in ``True`` on ``mu[0]`` before the reciprocal
    of ``mu`` is formed.  A gate named ``q_denominator_nonzero_gate`` (or the
    same common fallback) must return built-in ``True`` on
    ``L = 1-delta[0]*Z**2`` before the Q quotient is formed.  Passing one
    common gate is allowed when the caller uses the same premise for both
    denominators.  No gate is itself interpreted as an interval proof.

    The result represents every physical row as a common ``exp(log_scale)``
    factor times a returned scaled row.  ``exp`` is intentionally never
    evaluated.  The source identity, physical-scale certification, and
    numerical certification flags are all false.  The optional
    ``log_RQ_scale``/``scaled_RQ_y`` data also record the conditional
    cancellation of ``y`` between ``R = Rp*exp(y)`` and ``Q_y``: with the
    same remaining inputs, two distinct fixed-``y`` choices have the same
    conditional ``R*Q_y/Rp`` value, while nonzero ``m`` and ``Q`` rows are
    retained.
    """

    if fixed_y is not None:
        if y is not None:
            raise TypeError("pass fixed y once, as y or fixed_y, not both")
        y = fixed_y

    # ``None`` is not an admitted scalar; this also gives a clear error for a
    # call that omitted a required positional value.
    if delta0 is not None:
        if delta is not None:
            raise TypeError("pass delta rows once, as delta or delta0/delta1/delta2")
        delta = delta0

    if y is None or Z is None or delta is None:
        raise TypeError("y, Z, and delta are required")

    source_jets = {
        "Pstar": _coerce_jet3(Pstar, "Pstar"),
        "u": _coerce_jet3(u, "u"),
        "C1": _coerce_jet3(C1, "C1"),
        "C2": _coerce_jet3(C2, "C2"),
        "W11": _coerce_jet3(W11, "W11"),
        "W12": _coerce_jet3(W12, "W12"),
        "log_end_scale": _coerce_jet3(log_end_scale, "log_end_scale"),
        "mu": _coerce_jet3(mu, "mu"),
    }
    if delta1 is None and delta2 is None:
        delta_jet = _coerce_jet2(delta, "delta")
    elif delta1 is not None and delta2 is not None:
        # Also accept explicit (delta0, delta1, delta2) rows for callers that
        # do not package the three geometry rows in a jet object.
        delta_jet = OrdinaryZJet2((delta, delta1, delta2))  # type: ignore[arg-type]
    else:
        raise TypeError("delta1 and delta2 must either both be supplied or both omitted")
    _reject_primitive("fixed y", y)
    _reject_primitive("Z", Z)

    # A single common callback is retained as a convenience, while named
    # callbacks can distinguish the two independent caller premises.
    common_gate = denominator_nonzero_gate if denominator_nonzero_gate is not None else nonzero_gate
    mu_fallback = common_gate if common_gate is not None else q_denominator_nonzero_gate
    q_fallback = common_gate if common_gate is not None else mu_nonzero_gate
    mu_gate = _select_gate(mu_nonzero_gate, None, mu_fallback, "mu_nonzero_gate")
    q_gate = _select_gate(q_denominator_nonzero_gate, None, q_fallback, "q_denominator_nonzero_gate")

    # The mu gate is checked before any reciprocal is formed.
    mu0 = source_jets["mu"].rows[0]
    if mu_gate(mu0) is not True:
        raise ValueError("mu_nonzero_gate(mu[0]) must return built-in True")

    coefficient = _jet3_add(
        _jet3_product(source_jets["C1"], source_jets["W11"]),
        _jet3_product(source_jets["C2"], source_jets["W12"]),
    )
    amplitude = _jet3_neg(
        _jet3_product(
            _jet3_product(source_jets["Pstar"], source_jets["u"]),
            coefficient,
        )
    )

    # Fixed y is a constant ordinary-Z jet.  This is intentionally not the
    # fixed-xi numerator (13/2-xi)/mu after an xi=mu*y substitution.
    phi = _jet3_add(
        source_jets["log_end_scale"],
        _jet3_product(_constant3(Fraction(13, 2)), _jet3_reciprocal(source_jets["mu"])),
        _constant3(Fraction(-13, 1)),
        _constant3(_neg(y)),
    )

    # Ordinary Bell factors for exp(phi), after factoring exp(phi[0]).
    phi1, phi2, phi3 = phi.rows[1], phi.rows[2], phi.rows[3]
    bell = OrdinaryZJet3(
        (
            _ONE,
            phi1,
            _add(phi2, _mul(phi1, phi1)),
            _add(
                _add(phi3, _mul(Fraction(3, 1), _mul(phi1, phi2))),
                _mul(_mul(phi1, phi1), phi1),
            ),
        )
    )
    scaled_m = _jet3_product(amplitude, bell)

    # Build Q from scaled m.  The omitted common exponential cancels from the
    # quotient algebra, so Q carries exactly the same external scale.
    z_jet = OrdinaryZJet2((Z, _ONE, _ZERO))
    one_jet2 = _constant2(_ONE)
    z_square = _jet2_product(z_jet, z_jet)
    one_minus_z_square = _jet2_sub(one_jet2, z_square)
    one_minus_delta = _jet2_sub(one_jet2, delta_jet)
    denominator = _jet2_sub(one_jet2, _jet2_product(delta_jet, z_square))
    denominator0 = denominator.rows[0]
    if q_gate(denominator0) is not True:
        raise ValueError("q_denominator_nonzero_gate(1-delta[0]*Z**2) must return built-in True")

    m_value = OrdinaryZJet2((scaled_m.rows[0], scaled_m.rows[1], scaled_m.rows[2]))
    m_z = _differentiate3(scaled_m)
    numerator = _jet2_neg(
        _jet2_add(
            _jet2_product(_jet2_product(one_minus_delta, z_jet), m_value),
            _jet2_product(one_minus_z_square, m_z),
        )
    )
    scaled_Q = _jet2_product(numerator, _jet2_reciprocal(denominator))
    scaled_Qy_jet = _jet2_neg(OrdinaryZJet2((scaled_Q.rows[0], scaled_Q.rows[1], _ZERO)))
    scaled_Qy = (scaled_Qy_jet.rows[0], scaled_Qy_jet.rows[1])
    log_RQ_scale = _add(phi.rows[0], y)
    scaled_RQ_y = scaled_Qy[0]
    scaled_RQ_y_Z = scaled_Qy[1]

    return FGapFixedYLogScaledJets(
        log_scale=phi.rows[0],
        scaled_m=scaled_m,
        scaled_Q=scaled_Q,
        scaled_Qy=scaled_Qy,
        phi=phi,
        amplitude=amplitude,
        fixed_y=y,
        Z=Z,
        delta=delta_jet,
        log_RQ_scale=log_RQ_scale,
        scaled_RQ_y=scaled_RQ_y,
        scaled_RQ_y_Z=scaled_RQ_y_Z,
    )


# A descriptive alias for callers that use ``reconstruct`` terminology.
reconstruct_f_gap_fixed_y_log_scaled_jets = build_f_gap_fixed_y_log_scaled_jets


__all__ = [
    "FIXED_COORDINATE",
    "FGapFixedYLogScaledJets",
    "OrdinaryZJet2",
    "OrdinaryZJet3",
    "ROW_CONVENTION",
    "SCHEMA",
    "SOURCE_IDENTITY_PROVED",
    "PHYSICAL_SCALE_PROVED",
    "CERTIFIED",
    "build_f_gap_fixed_y_log_scaled_jets",
    "make_ordinary_z_jet2",
    "make_ordinary_z_jet3",
    "reconstruct_f_gap_fixed_y_log_scaled_jets",
]
