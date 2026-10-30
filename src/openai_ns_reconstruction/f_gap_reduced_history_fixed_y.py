"""Reduced-history fixed-``y`` F-gap algebra.

This is an independent caller-owned leaf.  It consumes ordinary ``Z`` rows
through order three and performs only exact ``Fraction``/opaque-interval
arithmetic.  The selected first-row relation is used in its reduced form

``sum(A_j C_j) = -mu*(incoming_factor1*incoming_m1
                     + ap*scaled_main_rows1)``

so the gap amplitude is the positive, once-normalized expression

``Pstar*u*mu*(incoming_factor1*incoming_m1
             + ap*scaled_main_rows1)``.

The fixed coordinate is ``y``.  Thus

``phi = log_end_scale + 13/(2*mu) - 13 - y``.

The common exponential is deliberately left outside the returned rows:
``m^(k) = exp(log_scale)*scaled_m[k]`` and likewise for ``Q`` and ``Q_y``.
No producer, source graph, old mathematical component, midpoint, float, or
``exp`` evaluation is used.  The caller must provide both nonzero gates.  All
F-source identity, physical-scale, and certification flags remain false.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import comb
from typing import Any, Callable, Generic, Mapping, Sequence, TypeVar


Scalar = TypeVar("Scalar")

SCHEMA = "f_gap_reduced_history_fixed_y_v1"
FIXED_COORDINATE = "fixed_y"
ROW_CONVENTION = "ordinary Z derivatives, not divided Taylor coefficients"
SOURCE_IDENTITY_PROVED = False
PHYSICAL_SCALE_PROVED = False
CERTIFIED = False
ORDERS3 = (0, 1, 2, 3)
ORDERS2 = (0, 1, 2)
_ZERO = Fraction(0, 1)
_ONE = Fraction(1, 1)


def _reject_primitive(name: str, value: Any) -> None:
    """Reject built-in scalar values other than exact ``Fraction`` values."""

    if value is None:
        raise TypeError(f"{name} is required")
    if isinstance(value, bool) or isinstance(value, int):
        raise TypeError(
            f"{name} must be Fraction or a caller-directed interval object; "
            "plain integers are not accepted"
        )
    if isinstance(value, float):
        raise TypeError(f"{name} must not be float; pass Fraction or an interval object")


@dataclass(frozen=True, slots=True)
class OrdinaryZJet3(Generic[Scalar]):
    """Four explicit ordinary derivatives with respect to ``Z``."""

    rows: tuple[Scalar, Scalar, Scalar, Scalar]

    def __post_init__(self) -> None:
        raw = tuple(self.rows)
        if len(raw) != 4:
            raise ValueError("OrdinaryZJet3 requires exactly four rows")
        for order, value in enumerate(raw):
            _reject_primitive(f"ordinary Z row {order}", value)
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
    """Three explicit ordinary derivatives with respect to ``Z``."""

    rows: tuple[Scalar, Scalar, Scalar]

    def __post_init__(self) -> None:
        raw = tuple(self.rows)
        if len(raw) != 3:
            raise ValueError("OrdinaryZJet2 requires exactly three rows")
        for order, value in enumerate(raw):
            _reject_primitive(f"ordinary Z row {order}", value)
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


Jet3Like = OrdinaryZJet3[Scalar] | Sequence[Scalar]
Jet2Like = OrdinaryZJet2[Scalar] | Sequence[Scalar]


def make_ordinary_z_jet3(
    z0: Scalar, z1: Scalar, z2: Scalar, z3: Scalar
) -> OrdinaryZJet3[Scalar]:
    return OrdinaryZJet3((z0, z1, z2, z3))


def make_ordinary_z_jet2(z0: Scalar, z1: Scalar, z2: Scalar) -> OrdinaryZJet2[Scalar]:
    return OrdinaryZJet2((z0, z1, z2))


def _coerce_jet3(value: Any, name: str) -> OrdinaryZJet3[Any]:
    if isinstance(value, OrdinaryZJet3):
        return value
    if isinstance(value, (str, bytes)):
        raise TypeError(f"{name} must provide four ordinary Z rows")
    if isinstance(value, Mapping) and "rows" in value:
        value = value["rows"]
    try:
        rows = tuple(value)
    except TypeError:
        try:
            rows = tuple(getattr(value, "rows"))
        except (AttributeError, TypeError) as exc:
            raise TypeError(f"{name} must provide four ordinary Z rows") from exc
    return OrdinaryZJet3(rows)


def _coerce_jet2(value: Any, name: str) -> OrdinaryZJet2[Any]:
    if isinstance(value, OrdinaryZJet2):
        return value
    if isinstance(value, (str, bytes)):
        raise TypeError(f"{name} must provide three ordinary Z rows")
    if isinstance(value, Mapping) and "rows" in value:
        value = value["rows"]
    try:
        rows = tuple(value)
    except TypeError:
        try:
            rows = tuple(getattr(value, "rows"))
        except (AttributeError, TypeError) as exc:
            raise TypeError(f"{name} must provide three ordinary Z rows") from exc
    return OrdinaryZJet2(rows)


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


def _div(left: Scalar, right: Scalar) -> Scalar:
    return left / right  # type: ignore[operator,return-value]


def _sum(terms: Sequence[Scalar]) -> Scalar:
    if not terms:
        raise ValueError("at least one term is required")
    total = terms[0]
    for term in terms[1:]:
        total = _add(total, term)
    return total


def _jet3_add(*jets: OrdinaryZJet3[Scalar]) -> OrdinaryZJet3[Scalar]:
    if not jets:
        raise ValueError("ordinary-jet addition requires an operand")
    return OrdinaryZJet3(
        tuple(_sum([jet.rows[order] for jet in jets]) for order in ORDERS3)  # type: ignore[arg-type]
    )


def _jet3_neg(jet: OrdinaryZJet3[Scalar]) -> OrdinaryZJet3[Scalar]:
    return OrdinaryZJet3(tuple(_neg(value) for value in jet.rows))  # type: ignore[arg-type]


def _jet3_scale(jet: OrdinaryZJet3[Scalar], scalar: Scalar) -> OrdinaryZJet3[Scalar]:
    return OrdinaryZJet3(tuple(_mul(scalar, value) for value in jet.rows))  # type: ignore[arg-type]


def _jet3_product(
    left: OrdinaryZJet3[Scalar], right: OrdinaryZJet3[Scalar]
) -> OrdinaryZJet3[Scalar]:
    return OrdinaryZJet3(
        tuple(
            _sum(
                [
                    _mul(
                        _mul(_binomial(order, split), left.rows[split]),
                        right.rows[order - split],
                    )
                    for split in range(order + 1)
                ]
            )
            for order in ORDERS3
        )  # type: ignore[arg-type]
    )


def _jet3_reciprocal(jet: OrdinaryZJet3[Scalar]) -> OrdinaryZJet3[Scalar]:
    rows: list[Scalar] = [_div(_ONE, jet.rows[0])]
    for order in range(1, 4):
        rows.append(
            _div(
                _neg(
                    _sum(
                        [
                            _mul(
                                _mul(_binomial(order, split), jet.rows[split]),
                                rows[order - split],
                            )
                            for split in range(1, order + 1)
                        ]
                    )
                ),
                jet.rows[0],
            )
        )
    return OrdinaryZJet3(tuple(rows))  # type: ignore[arg-type]


def _constant3(value: Scalar) -> OrdinaryZJet3[Scalar]:
    _reject_primitive("constant", value)
    return OrdinaryZJet3((value, _ZERO, _ZERO, _ZERO))


def _jet2_add(*jets: OrdinaryZJet2[Scalar]) -> OrdinaryZJet2[Scalar]:
    if not jets:
        raise ValueError("ordinary-jet addition requires an operand")
    return OrdinaryZJet2(
        tuple(_sum([jet.rows[order] for jet in jets]) for order in ORDERS2)  # type: ignore[arg-type]
    )


def _jet2_neg(jet: OrdinaryZJet2[Scalar]) -> OrdinaryZJet2[Scalar]:
    return OrdinaryZJet2(tuple(_neg(value) for value in jet.rows))  # type: ignore[arg-type]


def _jet2_product(
    left: OrdinaryZJet2[Scalar], right: OrdinaryZJet2[Scalar]
) -> OrdinaryZJet2[Scalar]:
    return OrdinaryZJet2(
        tuple(
            _sum(
                [
                    _mul(
                        _mul(_binomial(order, split), left.rows[split]),
                        right.rows[order - split],
                    )
                    for split in range(order + 1)
                ]
            )
            for order in ORDERS2
        )  # type: ignore[arg-type]
    )


def _jet2_reciprocal(jet: OrdinaryZJet2[Scalar]) -> OrdinaryZJet2[Scalar]:
    rows: list[Scalar] = [_div(_ONE, jet.rows[0])]
    for order in range(1, 3):
        rows.append(
            _div(
                _neg(
                    _sum(
                        [
                            _mul(
                                _mul(_binomial(order, split), jet.rows[split]),
                                rows[order - split],
                            )
                            for split in range(1, order + 1)
                        ]
                    )
                ),
                jet.rows[0],
            )
        )
    return OrdinaryZJet2(tuple(rows))  # type: ignore[arg-type]


def _constant2(value: Scalar) -> OrdinaryZJet2[Scalar]:
    _reject_primitive("constant", value)
    return OrdinaryZJet2((value, _ZERO, _ZERO))


@dataclass(frozen=True, slots=True)
class FGapReducedHistoryFixedY(Generic[Scalar]):
    """Unmaterialized fixed-``y`` reduced-history result."""

    amplitude: OrdinaryZJet3[Scalar]
    phi: OrdinaryZJet3[Scalar]
    scaled_m: OrdinaryZJet3[Scalar]
    log_scale: Scalar
    scaled_Q: OrdinaryZJet2[Scalar]
    scaled_Qy: tuple[Scalar, Scalar]
    fixed_y: Scalar
    Z: Scalar
    delta: OrdinaryZJet2[Scalar]
    reduced_history: OrdinaryZJet3[Scalar]
    source_row_target: OrdinaryZJet3[Scalar]
    incoming_product: OrdinaryZJet3[Scalar]
    selected_product: OrdinaryZJet3[Scalar]
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
    rho_conditions_conditioned: bool = True
    rq_y_independent_of_fixed_y: bool = True
    two_distinct_y_rqy_equal_conditioned: bool = True
    two_y_value_equality_proved: bool = False
    nonzero_m_Q_rows_retained: bool = True
    incoming_history_retained: bool = True
    incoming_nonzero_retained: bool = True
    c_rows_not_subtracted: bool = True
    c1_c2_subtraction_used: bool = False
    pstar_applied_once: bool = True
    pstar_double_applied: bool = False
    scale_materialized: bool = False
    physical_scale_materialized: bool = False
    source_identity: bool = False
    source_identity_proved: bool = False
    source_values_authenticated: bool = False
    physical_scale: bool = False
    physical_scale_proved: bool = False
    physical_scale_certified: bool = False
    certification: bool = False
    certified: bool = False
    error_certified: bool = False
    error_terms_deletable: bool = False
    mu_Z_rows_are_generalization: bool = True
    delta_Z_rows_are_generalization: bool = True
    source_parameter_independence_proved: bool = False
    source_mu_independence_proved: bool = False
    source_delta_independence_proved: bool = False
    source_Rp_independence_proved: bool = False
    physical_source_owner_proved: bool = False
    source_Rp_owner_authenticated: bool = False

    @property
    def scaled_m_rows(self) -> tuple[Scalar, Scalar, Scalar, Scalar]:
        return self.scaled_m.rows

    @property
    def m_rows(self) -> tuple[Scalar, Scalar, Scalar, Scalar]:
        return self.scaled_m.rows

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
    def Q_rows(self) -> tuple[Scalar, Scalar, Scalar]:
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
    def scaled_q_rows(self) -> tuple[Scalar, Scalar, Scalar]:
        return self.scaled_Q.rows

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
    def scaled_q(self) -> OrdinaryZJet2[Scalar]:
        return self.scaled_Q

    @property
    def scaled_Qy_rows(self) -> tuple[Scalar, Scalar]:
        return self.scaled_Qy

    @property
    def scaled_Qy0(self) -> Scalar:
        return self.scaled_Qy[0]

    @property
    def scaled_Qy1(self) -> Scalar:
        return self.scaled_Qy[1]

    @property
    def scaled_qy_rows(self) -> tuple[Scalar, Scalar]:
        return self.scaled_Qy

    @property
    def scaled_qy(self) -> tuple[Scalar, Scalar]:
        return self.scaled_Qy

    @property
    def log_RQy_Rp_scale(self) -> Scalar:
        return self.log_RQ_scale

    @property
    def log_RQy_over_Rp_scale(self) -> Scalar:
        return self.log_RQ_scale

    @property
    def RQy_Rp_log_scale(self) -> Scalar:
        return self.log_RQ_scale

    @property
    def scaled_RQy(self) -> tuple[Scalar, Scalar]:
        return (self.scaled_RQ_y, self.scaled_RQ_y_Z)

    @property
    def scaled_RQy_rows(self) -> tuple[Scalar, Scalar]:
        return self.scaled_RQy

    @property
    def scaled_RQ_y_rows(self) -> tuple[Scalar, Scalar]:
        return self.scaled_RQy

    @property
    def scaled_RQy_over_Rp(self) -> tuple[Scalar, Scalar]:
        return self.scaled_RQy

    @property
    def common_log_scale(self) -> Scalar:
        return self.log_scale

    @property
    def rho_conditions(self) -> Mapping[str, bool]:
        return {
            "fixed_y": self.fixed_y_premise,
            "U_zero": self.gap_U_zero_conditioned,
            "Uy_zero": self.gap_Uy_zero_conditioned,
            "my_equals_minus_m": self.my_identity_conditioned,
            "delta_y_zero": self.delta_y_zero_conditioned,
            "Rp_Z_zero": self.Rp_Z_zero_premise,
        }

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
            "fixed y and Rp_Z=0 are caller premises; exp is not materialized"
        )

    @property
    def flags(self) -> Mapping[str, bool]:
        return {
            "source_identity": self.source_identity,
            "source_identity_proved": self.source_identity_proved,
            "source_values_authenticated": self.source_values_authenticated,
            "physical_scale": self.physical_scale,
            "physical_scale_proved": self.physical_scale_proved,
            "physical_scale_certified": self.physical_scale_certified,
            "certification": self.certification,
            "certified": self.certified,
            "error_certified": self.error_certified,
            "error_terms_deletable": self.error_terms_deletable,
            "scale_materialized": self.scale_materialized,
            "physical_scale_materialized": self.physical_scale_materialized,
            "Rp_Z_zero_premise": self.Rp_Z_zero_premise,
            "rq_y_independent_of_fixed_y": self.rq_y_independent_of_fixed_y,
            "two_distinct_y_rqy_equal_conditioned": self.two_distinct_y_rqy_equal_conditioned,
            "two_y_value_equality_proved": self.two_y_value_equality_proved,
            "nonzero_m_Q_rows_retained": self.nonzero_m_Q_rows_retained,
            "incoming_history_retained": self.incoming_history_retained,
            "incoming_nonzero_retained": self.incoming_nonzero_retained,
            "c_rows_not_subtracted": self.c_rows_not_subtracted,
            "c1_c2_subtraction_used": self.c1_c2_subtraction_used,
            "pstar_applied_once": self.pstar_applied_once,
            "pstar_double_applied": self.pstar_double_applied,
            "mu_Z_rows_are_generalization": self.mu_Z_rows_are_generalization,
            "delta_Z_rows_are_generalization": self.delta_Z_rows_are_generalization,
            "physical_source_owner_proved": self.physical_source_owner_proved,
            "source_Rp_owner_authenticated": self.source_Rp_owner_authenticated,
        }

    def as_tuple(self) -> tuple[Any, ...]:
        return (
            self.log_scale,
            *self.scaled_m.rows,
            *self.scaled_Q.rows,
            *self.scaled_Qy,
            self.log_RQ_scale,
            self.scaled_RQ_y,
            self.scaled_RQ_y_Z,
        )

    def __iter__(self):
        return iter(self.as_tuple())

    def as_dict(self) -> dict[str, Any]:
        return {
            "schema": self.schema,
            "coordinate": self.coordinate,
            "row_convention": self.row_convention,
            "fixed_y": self.fixed_y,
            "Z": self.Z,
            "delta_rows": self.delta.rows,
            "phi_rows": self.phi.rows,
            "amplitude_rows": self.amplitude.rows,
            "reduced_history_rows": self.reduced_history.rows,
            "source_row_target_rows": self.source_row_target.rows,
            "incoming_product_rows": self.incoming_product.rows,
            "selected_product_rows": self.selected_product.rows,
            "log_scale": self.log_scale,
            "common_log_scale": self.common_log_scale,
            "scaled_m_rows": self.scaled_m.rows,
            "scaled_Q_rows": self.scaled_Q.rows,
            "scaled_Qy_rows": self.scaled_Qy,
            "log_RQ_scale": self.log_RQ_scale,
            "log_RQy_Rp_scale": self.log_RQy_Rp_scale,
            "log_RQy_over_Rp_scale": self.log_RQy_over_Rp_scale,
            "RQy_Rp_log_scale": self.RQy_Rp_log_scale,
            "scaled_RQ_y": self.scaled_RQ_y,
            "scaled_RQ_y_Z": self.scaled_RQ_y_Z,
            "scaled_RQy_over_Rp": self.scaled_RQy_over_Rp,
            "actual_m": self.actual_m_statement,
            "actual_Q": self.actual_Q_statement,
            "actual_Qy": self.actual_Qy_statement,
            "actual_RQy": self.actual_RQy_statement,
            "rho_conditions": dict(self.rho_conditions),
            "source_formula_conditioned": self.source_formula_conditioned,
            "fixed_y_premise": self.fixed_y_premise,
            "gap_U_zero_conditioned": self.gap_U_zero_conditioned,
            "gap_Uy_zero_conditioned": self.gap_Uy_zero_conditioned,
            "my_identity_conditioned": self.my_identity_conditioned,
            "delta_y_zero_conditioned": self.delta_y_zero_conditioned,
            "Rp_Z_zero_premise": self.Rp_Z_zero_premise,
            "rq_y_independent_of_fixed_y": self.rq_y_independent_of_fixed_y,
            "two_distinct_y_rqy_equal_conditioned": self.two_distinct_y_rqy_equal_conditioned,
            "two_y_value_equality_proved": self.two_y_value_equality_proved,
            "nonzero_m_Q_rows_retained": self.nonzero_m_Q_rows_retained,
            "incoming_history_retained": self.incoming_history_retained,
            "incoming_nonzero_retained": self.incoming_nonzero_retained,
            "c_rows_not_subtracted": self.c_rows_not_subtracted,
            "c1_c2_subtraction_used": self.c1_c2_subtraction_used,
            "pstar_applied_once": self.pstar_applied_once,
            "pstar_double_applied": self.pstar_double_applied,
            "mu_Z_rows_are_generalization": self.mu_Z_rows_are_generalization,
            "delta_Z_rows_are_generalization": self.delta_Z_rows_are_generalization,
            "source_identity": self.source_identity,
            "source_identity_proved": self.source_identity_proved,
            "source_values_authenticated": self.source_values_authenticated,
            "physical_scale": self.physical_scale,
            "physical_scale_proved": self.physical_scale_proved,
            "physical_scale_certified": self.physical_scale_certified,
            "certification": self.certification,
            "certified": self.certified,
            "error_certified": self.error_certified,
            "error_terms_deletable": self.error_terms_deletable,
            "scale_materialized": self.scale_materialized,
            "physical_scale_materialized": self.physical_scale_materialized,
            "flags": dict(self.flags),
        }


def _resolve_coordinates(
    coordinates: tuple[Any, ...],
    *,
    Z: Any,
    delta0: Any,
    delta1: Any,
    delta2: Any,
    fixed_y: Any,
    y: Any,
    delta: Any,
) -> tuple[Any, Any, Any, Any, Any]:
    """Resolve canonical ``(Z, delta0, delta1, delta2, fixed_y)`` inputs."""

    if coordinates:
        if len(coordinates) == 5:
            if any(item is not None for item in (Z, delta0, delta1, delta2, fixed_y, y, delta)):
                raise TypeError("coordinate rows must be positional or named, not both")
            Z, delta0, delta1, delta2, fixed_y = coordinates
        elif len(coordinates) == 4:
            if any(item is not None for item in (Z, delta0, delta1, delta2, y, delta)):
                raise TypeError("coordinate rows must be positional or named, not both")
            Z, delta0, delta1, delta2 = coordinates
        elif len(coordinates) == 3:
            if any(item is not None for item in (Z, delta0, delta1, delta2, fixed_y, y, delta)):
                raise TypeError("coordinate rows must be positional or named, not both")
            # Compatibility with the earlier fixed-y leaf: (fixed_y, Z, delta rows).
            fixed_y, Z, delta0 = coordinates
        else:
            raise TypeError(
                "coordinates require (Z, delta0, delta1, delta2, fixed_y) "
                "or compatibility (fixed_y, Z, delta_rows)"
            )

    if y is not None:
        if fixed_y is not None:
            raise TypeError("pass fixed_y once, as fixed_y or y")
        fixed_y = y
    if delta is not None:
        if delta0 is not None or delta1 is not None or delta2 is not None:
            raise TypeError("pass delta rows once, as delta or delta0/delta1/delta2")
        delta0 = delta
    if any(item is None for item in (Z, delta0, fixed_y)):
        raise TypeError("Z, delta rows, and fixed_y are required")
    if delta1 is None and delta2 is None:
        delta_jet = _coerce_jet2(delta0, "delta")
    elif delta1 is not None and delta2 is not None:
        delta_jet = OrdinaryZJet2((delta0, delta1, delta2))
    else:
        raise TypeError("delta1 and delta2 must both be supplied or both omitted")
    return Z, delta_jet.rows[0], delta_jet.rows[1], delta_jet.rows[2], fixed_y


def build_f_gap_reduced_history_fixed_y(
    Pstar: Jet3Like[Scalar],
    u: Jet3Like[Scalar],
    mu: Jet3Like[Scalar],
    incoming_factor1: Jet3Like[Scalar],
    incoming_m1: Jet3Like[Scalar],
    ap: Jet3Like[Scalar],
    scaled_main_rows1: Jet3Like[Scalar],
    log_end_scale: Jet3Like[Scalar],
    *coordinates: Any,
    mu_nonzero_gate: Callable[[Scalar], bool],
    q_denominator_nonzero_gate: Callable[[Scalar], bool],
    Z: Scalar | None = None,
    delta0: Scalar | OrdinaryZJet2[Scalar] | Sequence[Scalar] | None = None,
    delta1: Scalar | None = None,
    delta2: Scalar | None = None,
    fixed_y: Scalar | None = None,
    y: Scalar | None = None,
    delta: OrdinaryZJet2[Scalar] | Sequence[Scalar] | None = None,
) -> FGapReducedHistoryFixedY[Scalar]:
    """Build the reduced-history fixed-``y`` ordinary-``Z`` rows.

    The canonical positional coordinate order is ``(Z, delta0, delta1,
    delta2, fixed_y)``.  Named arguments use the same names.  For migration
    only, ``(fixed_y, Z, delta_rows)`` is also accepted.  Both gates are
    mandatory caller premises and must return the built-in ``True`` exactly
    once: the ``mu`` gate is checked before the reciprocal, and the Q gate is
    checked before the radial quotient.
    """

    if not callable(mu_nonzero_gate):
        raise TypeError("mu_nonzero_gate must be callable")
    if not callable(q_denominator_nonzero_gate):
        raise TypeError("q_denominator_nonzero_gate must be callable")

    Z, d0, d1, d2, fixed_y = _resolve_coordinates(
        coordinates,
        Z=Z,
        delta0=delta0,
        delta1=delta1,
        delta2=delta2,
        fixed_y=fixed_y,
        y=y,
        delta=delta,
    )
    _reject_primitive("Z", Z)
    _reject_primitive("fixed_y", fixed_y)

    source = {
        "Pstar": _coerce_jet3(Pstar, "Pstar"),
        "u": _coerce_jet3(u, "u"),
        "mu": _coerce_jet3(mu, "mu"),
        "incoming_factor1": _coerce_jet3(incoming_factor1, "incoming_factor1"),
        "incoming_m1": _coerce_jet3(incoming_m1, "incoming_m1"),
        "ap": _coerce_jet3(ap, "ap"),
        "scaled_main_rows1": _coerce_jet3(scaled_main_rows1, "scaled_main_rows1"),
        "log_end_scale": _coerce_jet3(log_end_scale, "log_end_scale"),
    }
    delta_jet = OrdinaryZJet2((d0, d1, d2))

    mu0 = source["mu"].rows[0]
    if mu_nonzero_gate(mu0) is not True:
        raise ValueError("mu_nonzero_gate(mu[0]) must return built-in True")

    incoming_product = _jet3_product(source["incoming_factor1"], source["incoming_m1"])
    selected_product = _jet3_product(source["ap"], source["scaled_main_rows1"])
    sum_rows = _jet3_add(incoming_product, selected_product)
    reduced_history = _jet3_product(source["mu"], sum_rows)
    source_row_target = _jet3_neg(reduced_history)
    amplitude = _jet3_product(
        _jet3_product(source["Pstar"], source["u"]), reduced_history
    )

    phi = _jet3_add(
        source["log_end_scale"],
        _jet3_scale(_jet3_reciprocal(source["mu"]), Fraction(13, 2)),
        _constant3(Fraction(-13, 1)),
        _constant3(_neg(fixed_y)),
    )
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

    z_jet = OrdinaryZJet2((Z, _ONE, _ZERO))
    z_square = _jet2_product(z_jet, z_jet)
    one = _constant2(_ONE)
    one_minus_delta = _jet2_add(one, _jet2_neg(delta_jet))
    one_minus_z_square = _jet2_add(one, _jet2_neg(z_square))
    denominator = _jet2_add(one, _jet2_neg(_jet2_product(delta_jet, z_square)))
    denominator0 = denominator.rows[0]
    if q_denominator_nonzero_gate(denominator0) is not True:
        raise ValueError(
            "q_denominator_nonzero_gate(1-delta[0]*Z**2) must return built-in True"
        )

    m_value = OrdinaryZJet2((scaled_m.rows[0], scaled_m.rows[1], scaled_m.rows[2]))
    m_z = OrdinaryZJet2((scaled_m.rows[1], scaled_m.rows[2], scaled_m.rows[3]))
    numerator = _jet2_neg(
        _jet2_add(
            _jet2_product(_jet2_product(one_minus_delta, z_jet), m_value),
            _jet2_product(one_minus_z_square, m_z),
        )
    )
    scaled_Q = _jet2_product(numerator, _jet2_reciprocal(denominator))
    scaled_Qy = (_neg(scaled_Q.rows[0]), _neg(scaled_Q.rows[1]))
    log_RQ_scale = _add(phi.rows[0], fixed_y)

    return FGapReducedHistoryFixedY(
        amplitude=amplitude,
        phi=phi,
        scaled_m=scaled_m,
        log_scale=phi.rows[0],
        scaled_Q=scaled_Q,
        scaled_Qy=scaled_Qy,
        fixed_y=fixed_y,
        Z=Z,
        delta=delta_jet,
        reduced_history=reduced_history,
        source_row_target=source_row_target,
        incoming_product=incoming_product,
        selected_product=selected_product,
        log_RQ_scale=log_RQ_scale,
        scaled_RQ_y=scaled_Qy[0],
        scaled_RQ_y_Z=scaled_Qy[1],
    )


reconstruct_f_gap_reduced_history_fixed_y = build_f_gap_reduced_history_fixed_y
build_f_gap_reduced_history_fixed_y_log_scaled_jets = build_f_gap_reduced_history_fixed_y
reconstruct_f_gap_reduced_history_fixed_y_log_scaled_jets = build_f_gap_reduced_history_fixed_y
FGapReducedHistoryFixedYLogScaledJets = FGapReducedHistoryFixedY


@dataclass(frozen=True, slots=True)
class FGapSignedSourceFunctional(Generic[Scalar]):
    """Two unmaterialized signed-gap terms and their ordinary-Z1 audit rows.

    ``linear_physical_Z_mantissa`` and
    ``quadratic_physical_Z_mantissa`` are the physical integral Z1 rows with
    their respective common exponential divided out.  They therefore already
    include the corresponding ``*_log_scale_Z`` contribution.  They are not
    derivatives of the bare mantissas.  A caller must not add
    ``linear_log_scale_Z * linear_mantissa`` or
    ``quadratic_log_scale_Z * quadratic_mantissa`` again: doing so double
    counts the scale derivative.  The ``*_log_scale_Z`` fields remain only as
    separate scale-path audit metadata.
    """

    linear_log_scale: Scalar
    linear_mantissa: Scalar
    linear_log_scale_Z: Scalar
    linear_physical_Z_mantissa: Scalar
    quadratic_log_scale: Scalar
    quadratic_mantissa: Scalar
    quadratic_log_scale_Z: Scalar
    quadratic_physical_Z_mantissa: Scalar
    c: Scalar
    c_Z: Scalar
    gap_decay: Scalar
    gap_decay_Z: Scalar
    coordinate: str = FIXED_COORDINATE
    schema: str = "f_gap_signed_source_functional_v1"
    source_identity_proved: bool = False
    physical_scale_proved: bool = False
    error_certified: bool = False
    certification: bool = False
    scale_materialized: bool = False
    source_Rp_owner_authenticated: bool = False
    source_fixed_R_requires_Rp_Z_zero: bool = True
    Rp_Z_zero_proved: bool = False
    gap_decay_relation_proved: bool = False
    endpoint_same_source_proved: bool = False
    historical_error_removed: bool = False

    @property
    def linear(self) -> tuple[Scalar, Scalar]:
        return (self.linear_log_scale, self.linear_mantissa)

    @property
    def quadratic(self) -> tuple[Scalar, Scalar]:
        return (self.quadratic_log_scale, self.quadratic_mantissa)

    @property
    def flags(self) -> Mapping[str, bool]:
        return {
            "source_identity_proved": self.source_identity_proved,
            "physical_scale_proved": self.physical_scale_proved,
            "error_certified": self.error_certified,
            "certification": self.certification,
            "scale_materialized": self.scale_materialized,
            "source_Rp_owner_authenticated": self.source_Rp_owner_authenticated,
            "source_fixed_R_requires_Rp_Z_zero": self.source_fixed_R_requires_Rp_Z_zero,
            "Rp_Z_zero_proved": self.Rp_Z_zero_proved,
            "gap_decay_relation_proved": self.gap_decay_relation_proved,
            "endpoint_same_source_proved": self.endpoint_same_source_proved,
            "historical_error_removed": self.historical_error_removed,
        }

    def as_dict(self) -> dict[str, Any]:
        return {
            "schema": self.schema,
            "coordinate": self.coordinate,
            "linear_log_scale": self.linear_log_scale,
            "linear_mantissa": self.linear_mantissa,
            "linear_log_scale_Z": self.linear_log_scale_Z,
            "linear_physical_Z_mantissa": self.linear_physical_Z_mantissa,
            "quadratic_log_scale": self.quadratic_log_scale,
            "quadratic_mantissa": self.quadratic_mantissa,
            "quadratic_log_scale_Z": self.quadratic_log_scale_Z,
            "quadratic_physical_Z_mantissa": self.quadratic_physical_Z_mantissa,
            "c": self.c,
            "c_Z": self.c_Z,
            "gap_decay": self.gap_decay,
            "gap_decay_Z": self.gap_decay_Z,
            "source_fixed_R_requires_Rp_Z_zero": self.source_fixed_R_requires_Rp_Z_zero,
            "Rp_Z_zero_proved": self.Rp_Z_zero_proved,
            "gap_decay_relation_proved": self.gap_decay_relation_proved,
            "endpoint_same_source_proved": self.endpoint_same_source_proved,
            "historical_error_removed": self.historical_error_removed,
            "flags": dict(self.flags),
        }


def analytic_gap_signed_source_functional(
    scaled_Q: Jet2Like[Scalar],
    Z: Scalar,
    delta0: Scalar,
    delta1: Scalar,
    K: Jet3Like[Scalar],
    logRp: Scalar,
    logRp_Z: Scalar,
    y_lower: Scalar,
    y_upper: Scalar,
    y_lower_Z: Scalar,
    y_upper_Z: Scalar,
    gap_decay: Scalar,
    *,
    q_denominator_nonzero_gate: Callable[[Scalar], bool],
) -> FGapSignedSourceFunctional[Scalar]:
    """Return the conditional signed gap functional without materializing ``exp``.

    ``K`` is the caller-supplied ordinary-Z jet for
    ``log_end_scale + 13/(2*mu) - 13``.  ``gap_decay`` is the caller-supplied
    directed value ``1-exp(-(y_upper-y_lower))``; it is not recomputed here.
    With ``qj = scaled_Q[j]`` and ``L = 1-delta*Z**2``, the two terms are

    ``exp(2*logRp + K + y_upper) * c*q1*gap_decay`` and
    ``exp(2*logRp + 2*K) * (-q0**2*(y_upper-y_lower)/4)``,

    where ``c=(1-delta)*Z/(4*L)``.  The returned
    ``linear_physical_Z_mantissa`` and
    ``quadratic_physical_Z_mantissa`` are the physical integral Z1 rows after
    division by their displayed common exponential.  They already include
    the displayed ``*_log_scale_Z`` terms, so a caller must not add those
    scale terms again.  The caller supplies ``gap_decay`` as a directed
    premise; this function neither validates the endpoint relation nor
    treats an interval upper bound or midpoint as its value.
    ``logRp_Z`` may be nonzero.  The fixed-R requirement is recorded, but
    ``Rp_Z_zero_proved`` remains false and no physical fixed-R equivalence is
    asserted.
    """

    if not callable(q_denominator_nonzero_gate):
        raise TypeError("q_denominator_nonzero_gate must be callable")
    q = _coerce_jet2(scaled_Q, "scaled_Q")
    k = _coerce_jet3(K, "K")
    for name, value in (
        ("Z", Z),
        ("delta0", delta0),
        ("delta1", delta1),
        ("logRp", logRp),
        ("logRp_Z", logRp_Z),
        ("y_lower", y_lower),
        ("y_upper", y_upper),
        ("y_lower_Z", y_lower_Z),
        ("y_upper_Z", y_upper_Z),
        ("gap_decay", gap_decay),
    ):
        _reject_primitive(name, value)
    delta = OrdinaryZJet2((delta0, delta1, _ZERO))
    z = OrdinaryZJet2((Z, _ONE, _ZERO))
    z2 = _jet2_product(z, z)
    denominator = _jet2_add(_constant2(_ONE), _jet2_neg(_jet2_product(delta, z2)))
    L, L_Z = denominator.rows[0], denominator.rows[1]
    if q_denominator_nonzero_gate(L) is not True:
        raise ValueError(
            "q_denominator_nonzero_gate(1-delta[0]*Z**2) must return built-in True"
        )

    one_minus_delta = OrdinaryZJet2((_sub(_ONE, delta0), _neg(delta1), _ZERO))
    numerator = _jet2_product(one_minus_delta, z)
    c = _div(numerator.rows[0], _mul(Fraction(4, 1), L))
    numerator_Z = numerator.rows[1]
    c_Z = _div(
        _sub(_mul(numerator_Z, L), _mul(numerator.rows[0], L_Z)),
        _mul(Fraction(4, 1), _mul(L, L)),
    )

    delta_y = _sub(y_upper, y_lower)
    delta_y_Z = _sub(y_upper_Z, y_lower_Z)
    gap_decay_Z = _mul(_sub(_ONE, gap_decay), delta_y_Z)
    q0, q1, q2 = q.rows
    linear_mantissa = _mul(_mul(c, q1), gap_decay)
    linear_physical_Z_mantissa = _sum(
        (
            _mul(_mul(c_Z, q1), gap_decay),
            _mul(
                c,
                _mul(
                    _add(q2, _mul(_add(_mul(Fraction(2, 1), logRp_Z), y_upper_Z), q1)),
                    gap_decay,
                ),
            ),
            _mul(_mul(c, q1), gap_decay_Z),
        )
    )
    quadratic_mantissa = _div(_neg(_mul(_mul(q0, q0), delta_y)), Fraction(4, 1))
    quadratic_physical_Z_mantissa = _div(
        _neg(
            _sum(
                (
                    _mul(_mul(Fraction(2, 1), logRp_Z), _mul(_mul(q0, q0), delta_y)),
                    _mul(_mul(Fraction(2, 1), _mul(q0, q1)), delta_y),
                    _mul(_mul(q0, q0), delta_y_Z),
                )
            )
        ),
        Fraction(4, 1),
    )
    linear_log_scale = _sum((_mul(Fraction(2, 1), logRp), k.rows[0], y_upper))
    linear_log_scale_Z = _sum((_mul(Fraction(2, 1), logRp_Z), k.rows[1], y_upper_Z))
    quadratic_log_scale = _sum((_mul(Fraction(2, 1), logRp), _mul(Fraction(2, 1), k.rows[0])))
    quadratic_log_scale_Z = _sum(
        (_mul(Fraction(2, 1), logRp_Z), _mul(Fraction(2, 1), k.rows[1]))
    )
    return FGapSignedSourceFunctional(
        linear_log_scale=linear_log_scale,
        linear_mantissa=linear_mantissa,
        linear_log_scale_Z=linear_log_scale_Z,
        linear_physical_Z_mantissa=linear_physical_Z_mantissa,
        quadratic_log_scale=quadratic_log_scale,
        quadratic_mantissa=quadratic_mantissa,
        quadratic_log_scale_Z=quadratic_log_scale_Z,
        quadratic_physical_Z_mantissa=quadratic_physical_Z_mantissa,
        c=c,
        c_Z=c_Z,
        gap_decay=gap_decay,
        gap_decay_Z=gap_decay_Z,
    )


build_gap_signed_source_functional = analytic_gap_signed_source_functional
build_f_gap_signed_source_functional = analytic_gap_signed_source_functional
build_f_gap_analytic_signed_source_functional = analytic_gap_signed_source_functional


def hand_check_denominator_nonzero(value: Fraction) -> bool:
    return value != _ZERO


def hand_check_fixture() -> Mapping[str, Any]:
    """Return an exact, independent hand-check packet (the builder is not run)."""

    f = Fraction
    constants = (f(1), f(0), f(0), f(0))
    inputs = {
        "Pstar": constants,
        "u": constants,
        "mu": (f(1, 2), f(0), f(0), f(0)),
        "incoming_factor1": constants,
        "incoming_m1": constants,
        "ap": constants,
        "scaled_main_rows1": constants,
        "log_end_scale": (f(23), f(0), f(0), f(0)),
        "Z": f(1, 2),
        "delta0": f(1, 4),
        "delta1": f(0),
        "delta2": f(0),
        "fixed_y": f(22),
        "mu_nonzero_gate": hand_check_denominator_nonzero,
        "q_denominator_nonzero_gate": hand_check_denominator_nonzero,
    }
    return {
        "inputs": inputs,
        "expected": {
            "amplitude": (f(1), f(0), f(0), f(0)),
            "phi": (f(1), f(0), f(0), f(0)),
            "scaled_m": (f(1), f(0), f(0), f(0)),
            "scaled_Q": (f(-2, 5), f(-68, 75), f(-784, 1125)),
            "scaled_Qy": (f(2, 5), f(68, 75)),
            "log_scale": f(1),
            "log_RQ_scale": f(23),
            "scaled_RQ_y": f(2, 5),
            "scaled_RQ_y_Z": f(68, 75),
            "reduced_history": (f(1), f(0), f(0), f(0)),
            "source_row_target": (f(-1), f(0), f(0), f(0)),
        },
        "unused_constant_jet": constants,
        "source_identity_proved": False,
        "physical_scale_proved": False,
        "error_certified": False,
        "coordinate": FIXED_COORDINATE,
    }


__all__ = [
    "CERTIFIED",
    "FIXED_COORDINATE",
    "FGapReducedHistoryFixedY",
    "FGapReducedHistoryFixedYLogScaledJets",
    "FGapSignedSourceFunctional",
    "OrdinaryZJet2",
    "OrdinaryZJet3",
    "PHYSICAL_SCALE_PROVED",
    "ROW_CONVENTION",
    "SCHEMA",
    "SOURCE_IDENTITY_PROVED",
    "analytic_gap_signed_source_functional",
    "build_f_gap_analytic_signed_source_functional",
    "build_f_gap_reduced_history_fixed_y",
    "build_f_gap_reduced_history_fixed_y_log_scaled_jets",
    "build_f_gap_signed_source_functional",
    "build_gap_signed_source_functional",
    "hand_check_denominator_nonzero",
    "hand_check_fixture",
    "make_ordinary_z_jet2",
    "make_ordinary_z_jet3",
    "reconstruct_f_gap_reduced_history_fixed_y",
    "reconstruct_f_gap_reduced_history_fixed_y_log_scaled_jets",
]
