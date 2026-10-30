"""Fixed-physical-radius transport for a scaled gap radial ``Q`` field.

This leaf is deliberately independent of every source producer, graph, old
mathematical module, interval evaluator, and test fixture.  The caller owns
the source meaning of all inputs.  It supplies four ordinary fixed-``y``
rows of the common-scale-factored mean ``m`` and the three ordinary rows of
``delta``.  The rows are kept as exact arithmetic or as caller-directed
interval objects; no midpoint, binary float, source value, or exponential is
created here.

The radial derivative is reconstructed first:

``m_R = m_Z|y + r*m``,
``Q_R = -(h*Z*m + d*m_R)/L``,
``h = 1-delta``, ``d = 1-Z**2``, ``L = 1-delta*Z**2``.

The returned ``fixed_y_scaled_Q_R`` rows are the ordinary fixed-``y`` rows
of this radial field, with the same external common scale as the input.  The
physical-radius rows then use the caller premise ``Q_y = -Q``:

``Q_Z|R  = Q_Z|y + r*Q``,
``Q_ZZ|R = Q_ZZ|y + 2*r*Q_Z|y + (r_Z+r**2)*Q``.

The ``r``, ``r_Z``, and ``r_ZZ`` values are recorded as caller premises only.
No endpoint or error term is removed, and every source/physical/error proof
flag remains false.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import comb
from collections.abc import Mapping, Sequence
from typing import Any, Callable, Generic, TypeVar


Scalar = TypeVar("Scalar")

SCHEMA = "gap_fixed_radius_transport_v2"
ROW_CONVENTION = "ordinary Z derivatives, not divided Taylor coefficients"
FIXED_COORDINATE = "fixed_y"
SOURCE_IDENTITY_PROVED = False
PHYSICAL_SCALE_PROVED = False
ERROR_CERTIFIED = False
CERTIFIED = False

_ZERO = Fraction(0, 1)
_ONE = Fraction(1, 1)
_TWO = Fraction(2, 1)
_MISSING = object()


class GapFixedRadiusTransportInputError(TypeError):
    """Raised when a caller input is not an admitted exact/interval value."""


def _is_midpoint(value: Any) -> bool:
    """Reject common explicit midpoint wrappers without inspecting intervals."""

    type_name = type(value).__name__.lower()
    if "midpoint" in type_name:
        return True
    for attribute in ("is_midpoint", "_is_midpoint"):
        try:
            if getattr(value, attribute, False) is True:
                return True
        except Exception:  # pragma: no cover - opaque caller object
            return False
    if isinstance(value, Mapping) and "midpoint" in value:
        return True
    return False


def _reject_scalar(name: str, value: Any) -> None:
    """Reject primitive point substitutions while retaining caller arithmetic."""

    if value is _MISSING or value is None:
        raise GapFixedRadiusTransportInputError(f"{name} is required")
    if isinstance(value, bool) or isinstance(value, int):
        raise GapFixedRadiusTransportInputError(
            f"{name} must be Fraction or a caller-directed interval object; "
            "plain integers are not admitted"
        )
    if isinstance(value, float):
        raise GapFixedRadiusTransportInputError(
            f"{name} must not be binary float; pass Fraction or an interval object"
        )
    if isinstance(value, (str, bytes, bytearray, list, tuple, set, frozenset, dict)):
        raise GapFixedRadiusTransportInputError(
            f"{name} must be one scalar, not a container or text value"
        )
    if _is_midpoint(value):
        raise GapFixedRadiusTransportInputError(
            f"{name} must not be a midpoint or midpoint-selected value"
        )


def _binomial(order: int, split: int) -> Fraction:
    return Fraction(comb(order, split), 1)


def _add(left: Scalar, right: Scalar) -> Scalar:
    return left + right  # type: ignore[operator,return-value]


def _sub(left: Scalar, right: Scalar) -> Scalar:
    return left - right  # type: ignore[operator,return-value]


def _neg(value: Scalar) -> Scalar:
    return -value  # type: ignore[operator,return-value]


def _mul(left: Scalar, right: Scalar) -> Scalar:
    try:
        return left * right  # type: ignore[operator,return-value]
    except TypeError:
        # Some directed interval implementations expose scalar multiplication
        # only through the reflected operation.  No coercion or midpoint is
        # introduced by this fallback.
        return right * left  # type: ignore[operator,return-value]


def _div(left: Scalar, right: Scalar) -> Scalar:
    if isinstance(left, Fraction) and isinstance(right, Fraction):
        return left / right  # type: ignore[return-value]
    return left / right  # type: ignore[operator,return-value]


def _sum(terms: Sequence[Scalar]) -> Scalar:
    if not terms:
        raise ValueError("at least one ordinary-jet term is required")
    total = terms[0]
    for term in terms[1:]:
        total = _add(total, term)
    return total


@dataclass(frozen=True, slots=True)
class OrdinaryZJet3(Generic[Scalar]):
    """Four explicit ordinary fixed-``y`` derivatives in ``Z``."""

    rows: tuple[Scalar, Scalar, Scalar, Scalar]

    def __post_init__(self) -> None:
        raw = tuple(self.rows)
        if len(raw) != 4:
            raise ValueError("OrdinaryZJet3 requires exactly four rows")
        for order, value in enumerate(raw):
            _reject_scalar(f"ordinary Z row {order}", value)
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
    """Three explicit ordinary fixed-``y`` derivatives in ``Z``."""

    rows: tuple[Scalar, Scalar, Scalar]

    def __post_init__(self) -> None:
        raw = tuple(self.rows)
        if len(raw) != 3:
            raise ValueError("OrdinaryZJet2 requires exactly three rows")
        for order, value in enumerate(raw):
            _reject_scalar(f"ordinary Z row {order}", value)
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
    """Build four explicit ordinary rows without changing their convention."""

    return OrdinaryZJet3((z0, z1, z2, z3))


def make_ordinary_z_jet2(
    z0: Scalar, z1: Scalar, z2: Scalar
) -> OrdinaryZJet2[Scalar]:
    """Build three explicit ordinary rows without changing their convention."""

    return OrdinaryZJet2((z0, z1, z2))


def _coerce_jet3(value: Any, name: str) -> OrdinaryZJet3[Any]:
    if isinstance(value, OrdinaryZJet3):
        return value
    if isinstance(value, (str, bytes)):
        raise GapFixedRadiusTransportInputError(
            f"{name} must provide four ordinary Z rows"
        )
    if isinstance(value, Mapping) and "rows" in value:
        value = value["rows"]
    try:
        rows = tuple(value)
    except (TypeError, ValueError) as exc:
        try:
            rows = tuple(getattr(value, "rows"))
        except (AttributeError, TypeError, ValueError) as row_exc:
            raise GapFixedRadiusTransportInputError(
                f"{name} must provide four ordinary Z rows"
            ) from row_exc
    try:
        return OrdinaryZJet3(rows)
    except (TypeError, ValueError) as exc:
        raise type(exc)(f"{name}: {exc}") from exc


def _coerce_jet2(value: Any, name: str) -> OrdinaryZJet2[Any]:
    if isinstance(value, OrdinaryZJet2):
        return value
    if isinstance(value, (str, bytes)):
        raise GapFixedRadiusTransportInputError(
            f"{name} must provide three ordinary Z rows"
        )
    if isinstance(value, Mapping) and "rows" in value:
        value = value["rows"]
    try:
        rows = tuple(value)
    except (TypeError, ValueError) as exc:
        try:
            rows = tuple(getattr(value, "rows"))
        except (AttributeError, TypeError, ValueError) as row_exc:
            raise GapFixedRadiusTransportInputError(
                f"{name} must provide three ordinary Z rows"
            ) from row_exc
    try:
        return OrdinaryZJet2(rows)
    except (TypeError, ValueError) as exc:
        raise type(exc)(f"{name}: {exc}") from exc


def _constant2(value: Scalar) -> OrdinaryZJet2[Scalar]:
    _reject_scalar("fixed scalar", value)
    return OrdinaryZJet2((value, _ZERO, _ZERO))  # type: ignore[arg-type]


def _jet2_add(*jets: OrdinaryZJet2[Scalar]) -> OrdinaryZJet2[Scalar]:
    if not jets:
        raise ValueError("ordinary-jet addition requires an operand")
    return OrdinaryZJet2(
        tuple(_sum([jet.rows[order] for jet in jets]) for order in range(3))  # type: ignore[arg-type]
    )


def _jet2_sub(
    left: OrdinaryZJet2[Scalar], right: OrdinaryZJet2[Scalar]
) -> OrdinaryZJet2[Scalar]:
    return OrdinaryZJet2(
        tuple(_sub(left.rows[order], right.rows[order]) for order in range(3))  # type: ignore[arg-type]
    )


def _jet2_neg(jet: OrdinaryZJet2[Scalar]) -> OrdinaryZJet2[Scalar]:
    return OrdinaryZJet2(tuple(_neg(value) for value in jet.rows))  # type: ignore[arg-type]


def _jet2_product(
    left: OrdinaryZJet2[Scalar], right: OrdinaryZJet2[Scalar]
) -> OrdinaryZJet2[Scalar]:
    """Multiply ordinary rows with exact binomial Leibniz factors."""

    rows: list[Scalar] = []
    for order in range(3):
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
    return OrdinaryZJet2(tuple(rows))  # type: ignore[arg-type]


def _jet2_reciprocal(jet: OrdinaryZJet2[Scalar]) -> OrdinaryZJet2[Scalar]:
    """Differentiate ``jet * reciprocal = 1`` through ordinary order two."""

    rows: list[Scalar] = [_div(_ONE, jet.rows[0])]  # type: ignore[arg-type]
    for order in range(1, 3):
        terms = [
            _mul(
                _mul(_binomial(order, split), jet.rows[split]),
                rows[order - split],
            )
            for split in range(1, order + 1)
        ]
        rows.append(_div(_neg(_sum(terms)), jet.rows[0]))
    return OrdinaryZJet2(tuple(rows))  # type: ignore[arg-type]


@dataclass(frozen=True, slots=True)
class UnmaterializedCommonScale(Generic[Scalar]):
    """A scale token whose exponential is intentionally never evaluated."""

    log_scale: Scalar | None = None
    statement: str = "exp(common_log_scale)"
    materialized: bool = False


@dataclass(frozen=True, slots=True)
class GapFixedRadiusTransport(Generic[Scalar]):
    """Scaled fixed-radius ``Q`` rows and their caller-premise metadata."""

    common_scale: Any
    scaled_m: OrdinaryZJet3[Scalar]
    scaled_m_R: OrdinaryZJet2[Scalar]
    fixed_y_scaled_Q_R: OrdinaryZJet2[Scalar]
    Q_R: Scalar
    Q_RZ: Scalar
    Q_RZZ: Scalar
    fixed_y_Q_RZ: Scalar
    fixed_y_Q_RZZ: Scalar
    r: Scalar
    r_Z: Scalar
    r_ZZ: Scalar
    Z: Scalar
    delta: OrdinaryZJet2[Scalar]
    h: OrdinaryZJet2[Scalar]
    d: OrdinaryZJet2[Scalar]
    L: OrdinaryZJet2[Scalar]
    numerator: OrdinaryZJet2[Scalar]
    denominator_nonzero_gate_supplied: bool = True
    Q_y_equals_minus_Q_caller_premise: bool = True
    m_y_equals_minus_m_caller_premise: bool = True
    fixed_y_rows_are_scaled: bool = True
    r_is_logRp_Z_caller_premise: bool = True
    r_Z_is_derivative_caller_premise: bool = True
    r_ZZ_is_derivative_caller_premise: bool = True
    common_scale_materialized: bool = False
    source_identity_proved: bool = False
    physical_scale_proved: bool = False
    error_certified: bool = False
    certified: bool = False
    endpoint_terms_cancelled: bool = False
    errors_removed: bool = False

    @property
    def fixed_y_Q_R(self) -> OrdinaryZJet2[Scalar]:
        return self.fixed_y_scaled_Q_R

    @property
    def schema(self) -> str:
        return SCHEMA

    @property
    def coordinate(self) -> str:
        return FIXED_COORDINATE

    @property
    def row_convention(self) -> str:
        return ROW_CONVENTION

    @property
    def fixed_y_scaled_Q_RZ(self) -> Scalar:
        return self.fixed_y_scaled_Q_R.rows[1]

    @property
    def fixed_y_scaled_Q_RZZ(self) -> Scalar:
        return self.fixed_y_scaled_Q_R.rows[2]

    @property
    def common_log_scale(self) -> Any:
        """Expose a supplied log-scale token without evaluating ``exp``."""

        if isinstance(self.common_scale, UnmaterializedCommonScale):
            return self.common_scale.log_scale
        return self.common_scale

    @property
    def scaled_Q_R(self) -> OrdinaryZJet2[Scalar]:
        return self.fixed_y_scaled_Q_R

    @property
    def scaled_Q_R_rows(self) -> tuple[Scalar, Scalar, Scalar]:
        return self.fixed_y_scaled_Q_R.rows

    @property
    def Q_R_rows(self) -> tuple[Scalar, Scalar, Scalar]:
        return (self.Q_R, self.Q_RZ, self.Q_RZZ)

    @property
    def fixed_radius_scaled_Q_R_rows(self) -> tuple[Scalar, Scalar, Scalar]:
        """The returned fixed-physical-``R`` rows, before the common scale."""

        return self.Q_R_rows

    @property
    def scaled_Q_RZ(self) -> Scalar:
        return self.Q_RZ

    @property
    def scaled_Q_RZZ(self) -> Scalar:
        return self.Q_RZZ

    @property
    def actual_Q_statement(self) -> str:
        return (
            "(Q_R, Q_Z|R, Q_ZZ|R) = exp(common_log_scale) * "
            "(Q_R, Q_RZ, Q_RZZ); exp is not materialized"
        )

    @property
    def actual_Q_RZ_statement(self) -> str:
        return (
            "Q_Z|R^(k) uses the same exp(common_log_scale) factor; "
            "the fixed-radius rows are returned scaled"
        )

    @property
    def source_evidence_false(self) -> bool:
        return not self.source_identity_proved

    @property
    def physical_evidence_false(self) -> bool:
        return not self.physical_scale_proved

    @property
    def error_evidence_false(self) -> bool:
        return not self.error_certified

    def as_dict(self) -> dict[str, Any]:
        """Return a lossless metadata view without coercing caller scalars."""

        return {
            "schema": SCHEMA,
            "common_scale": self.common_scale,
            "common_log_scale": self.common_log_scale,
            "common_scale_materialized": self.common_scale_materialized,
            "row_convention": ROW_CONVENTION,
            "coordinate": FIXED_COORDINATE,
            "scaled_m_rows": self.scaled_m.rows,
            "scaled_m_R_rows": self.scaled_m_R.rows,
            "fixed_y_scaled_Q_R_rows": self.fixed_y_scaled_Q_R.rows,
            "fixed_y_scaled_Q_RZ": self.fixed_y_scaled_Q_RZ,
            "fixed_y_scaled_Q_RZZ": self.fixed_y_scaled_Q_RZZ,
            "Q_R": self.Q_R,
            "Q_RZ": self.Q_RZ,
            "Q_RZZ": self.Q_RZZ,
            "fixed_radius_scaled_Q_R_rows": self.fixed_radius_scaled_Q_R_rows,
            "fixed_y_Q_RZ": self.fixed_y_Q_RZ,
            "fixed_y_Q_RZZ": self.fixed_y_Q_RZZ,
            "r": self.r,
            "r_Z": self.r_Z,
            "r_ZZ": self.r_ZZ,
            "Z": self.Z,
            "delta_rows": self.delta.rows,
            "L_rows": self.L.rows,
            "actual_Q": self.actual_Q_statement,
            "actual_Q_RZ": self.actual_Q_RZ_statement,
            "premises": {
                "L_nonzero_gate_supplied": self.denominator_nonzero_gate_supplied,
                "Q_y_equals_minus_Q": self.Q_y_equals_minus_Q_caller_premise,
                "m_y_equals_minus_m": self.m_y_equals_minus_m_caller_premise,
                "r_is_logRp_Z": self.r_is_logRp_Z_caller_premise,
                "r_Z_is_derivative": self.r_Z_is_derivative_caller_premise,
                "r_ZZ_is_derivative": self.r_ZZ_is_derivative_caller_premise,
            },
            "source_identity_proved": self.source_identity_proved,
            "source_evidence_false": self.source_evidence_false,
            "physical_scale_proved": self.physical_scale_proved,
            "physical_evidence_false": self.physical_evidence_false,
            "error_certified": self.error_certified,
            "error_evidence_false": self.error_evidence_false,
            "certified": self.certified,
            "endpoint_terms_cancelled": self.endpoint_terms_cancelled,
            "errors_removed": self.errors_removed,
        }


def _resolve_scale(
    common_scale: Any,
    log_scale: Any,
) -> Any:
    if common_scale is not _MISSING and log_scale is not _MISSING:
        raise TypeError("pass common_scale or log_scale, not both")
    if log_scale is not _MISSING:
        _reject_scalar("log_scale", log_scale)
        return UnmaterializedCommonScale(log_scale=log_scale)
    if common_scale is _MISSING:
        return UnmaterializedCommonScale()
    # A supplied token is retained exactly.  It is metadata only and never
    # participates in arithmetic; the returned rows remain scale-factored.
    if isinstance(common_scale, (bool, int, float, Fraction)):
        raise GapFixedRadiusTransportInputError(
            "common_scale must remain an unmaterialized caller token"
        )
    if _is_midpoint(common_scale):
        raise GapFixedRadiusTransportInputError(
            "common_scale must not be a midpoint-selected value"
        )
    return common_scale


def build_gap_fixed_radius_transport(
    scaled_m: Jet3Like[Scalar],
    Z: Scalar,
    delta0: Scalar,
    delta1: Scalar,
    delta2: Scalar,
    r: Scalar,
    r_Z: Scalar,
    r_ZZ: Scalar,
    *,
    denominator_nonzero_gate: Callable[[Scalar], bool],
    common_scale: Any = _MISSING,
    log_scale: Any = _MISSING,
) -> GapFixedRadiusTransport[Scalar]:
    """Build fixed-physical-``R`` scaled ``Q`` rows through ordinary ``Z2``.

    ``scaled_m`` contains the four ordinary fixed-``y`` rows ``m0`` through
    ``m3``.  The first radial row is reconstructed as ``m1 + r*m0``.  The
    denominator gate is mandatory and is called exactly once on ``L[0]``
    before a reciprocal is formed.  ``common_scale`` (or ``log_scale``) is
    only carried as an opaque token; no exponential is evaluated.
    """

    if not callable(denominator_nonzero_gate):
        raise TypeError("denominator_nonzero_gate must be callable")

    m = _coerce_jet3(scaled_m, "scaled_m")
    for name, value in (
        ("Z", Z),
        ("delta0", delta0),
        ("delta1", delta1),
        ("delta2", delta2),
        ("r", r),
        ("r_Z", r_Z),
        ("r_ZZ", r_ZZ),
    ):
        _reject_scalar(name, value)

    delta = OrdinaryZJet2((delta0, delta1, delta2))
    z_jet = OrdinaryZJet2((Z, _ONE, _ZERO))
    r_jet = OrdinaryZJet2((r, r_Z, r_ZZ))
    m_value = OrdinaryZJet2((m.rows[0], m.rows[1], m.rows[2]))
    m_fixed_y_Z = OrdinaryZJet2((m.rows[1], m.rows[2], m.rows[3]))

    # m_R = m_Z|y + r*m, with the product using ordinary binomial rows.
    m_R = _jet2_add(m_fixed_y_Z, _jet2_product(r_jet, m_value))

    one = _constant2(_ONE)
    z_square = _jet2_product(z_jet, z_jet)
    h = _jet2_sub(one, delta)
    d = _jet2_sub(one, z_square)
    L = _jet2_sub(one, _jet2_product(delta, z_square))

    # This is the only gate invocation, and it precedes the reciprocal.
    if denominator_nonzero_gate(L.rows[0]) is not True:
        raise ValueError(
            "denominator_nonzero_gate(1-delta[0]*Z**2) must return built-in True"
        )

    numerator = _jet2_neg(
        _jet2_add(
            _jet2_product(_jet2_product(h, z_jet), m_value),
            _jet2_product(d, m_R),
        )
    )
    fixed_y_scaled_Q_R = _jet2_product(numerator, _jet2_reciprocal(L))

    q0, q1, q2 = fixed_y_scaled_Q_R.rows
    Q_R = q0
    Q_RZ = _add(q1, _mul(r, q0))
    Q_RZZ = _add(
        _add(q2, _mul(_mul(_TWO, r), q1)),
        _mul(_add(r_Z, _mul(r, r)), q0),
    )

    return GapFixedRadiusTransport(
        common_scale=_resolve_scale(common_scale, log_scale),
        scaled_m=m,
        scaled_m_R=m_R,
        fixed_y_scaled_Q_R=fixed_y_scaled_Q_R,
        Q_R=Q_R,
        Q_RZ=Q_RZ,
        Q_RZZ=Q_RZZ,
        fixed_y_Q_RZ=q1,
        fixed_y_Q_RZZ=q2,
        r=r,
        r_Z=r_Z,
        r_ZZ=r_ZZ,
        Z=Z,
        delta=delta,
        h=h,
        d=d,
        L=L,
        numerator=numerator,
    )


def build_gap_fixed_radius_transport_from_rows(
    m0: Scalar,
    m1: Scalar,
    m2: Scalar,
    m3: Scalar,
    Z: Scalar,
    delta0: Scalar,
    delta1: Scalar,
    delta2: Scalar,
    r: Scalar,
    r_Z: Scalar,
    r_ZZ: Scalar,
    *,
    denominator_nonzero_gate: Callable[[Scalar], bool],
    common_scale: Any = _MISSING,
    log_scale: Any = _MISSING,
) -> GapFixedRadiusTransport[Scalar]:
    """Row-wise convenience interface for callers without a jet wrapper."""

    return build_gap_fixed_radius_transport(
        OrdinaryZJet3((m0, m1, m2, m3)),
        Z,
        delta0,
        delta1,
        delta2,
        r,
        r_Z,
        r_ZZ,
        denominator_nonzero_gate=denominator_nonzero_gate,
        common_scale=common_scale,
        log_scale=log_scale,
    )


# Descriptive aliases keep the formula interface discoverable without adding
# another implementation or importing any prior mathematical component.
reconstruct_gap_fixed_radius_transport = build_gap_fixed_radius_transport
transport_gap_fixed_radius = build_gap_fixed_radius_transport
compute_gap_fixed_radius_transport = build_gap_fixed_radius_transport


__all__ = [
    "CERTIFIED",
    "ERROR_CERTIFIED",
    "FIXED_COORDINATE",
    "GapFixedRadiusTransport",
    "GapFixedRadiusTransportInputError",
    "OrdinaryZJet2",
    "OrdinaryZJet3",
    "PHYSICAL_SCALE_PROVED",
    "ROW_CONVENTION",
    "SCHEMA",
    "SOURCE_IDENTITY_PROVED",
    "UnmaterializedCommonScale",
    "build_gap_fixed_radius_transport",
    "build_gap_fixed_radius_transport_from_rows",
    "compute_gap_fixed_radius_transport",
    "make_ordinary_z_jet2",
    "make_ordinary_z_jet3",
    "reconstruct_gap_fixed_radius_transport",
    "transport_gap_fixed_radius",
]


