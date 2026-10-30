"""Caller-owned closed-form signed ``P-F`` gap differences.

This leaf contains only the algebra needed after a caller has supplied the
gap closed-form ports.  It never imports a source producer, evaluates an
exponential, forms a raw ``P-F`` subtraction, or combines terms carrying
different logarithmic scales.

For each of the four independent rows
``linear``, ``quadratic``, ``endpoint_upper`` and ``endpoint_lower`` the
caller supplies

``LF``
    The frame log-scale token.  It is retained exactly and is never
    exponentiated or otherwise materialised.
``aF``
    The frame mantissa normalised with that row's ``LF``; it is not the
    physical frame value.
``delta_a``
    The direct caller-supplied correlated mantissa difference ``aP - aF``
    between the separately normalised ``aP`` (at ``LP``) and ``aF`` (at
    ``LF``).  It is a mantissa premise, not the physical ``P-F`` value and not
    a subtraction performed by this module; ``(1+t)`` supplies the scale ratio.
``aF_physZ``
    The complete frame physical ordinary-``Z1`` mantissa for that row,
    including the row's scale and any moving effects already supplied by the
    caller.
``delta_physZ``
    The direct caller-supplied signed physical ordinary-``Z1`` difference.
``t``
    The caller's value of ``expm1(LP - LF)``.  It is retained as supplied; no
    exponential or ``expm1`` routine is called here.

The closed form is applied independently to every row:

``signed_value_mantissa = (1 + t) * delta_a + t * aF``
``signed_physicalZ_mantissa = (1 + t) * delta_physZ + t * aF_physZ``

``aF_physZ`` is already a complete physical ``Z1`` row.  The leaf therefore
does not add ``LF_Z`` or ``t_Z`` terms.  The returned four components keep
their own untouched ``LF`` token and carry the m15 orientation
``linear+``, ``quadratic+``, ``upper-``, ``lower+``.  They are not summed: the
four scales need not agree, and retaining both endpoint terms is part of the
  contract.  ``ordinary Z1`` here names only the derivative-row convention.
  ``Pi1``, pressure ``S1``, ``P1`` and ``m15`` are distinct roles that this leaf
  does not identify.  The preserved role distinction is ``P1 != m15``.

The companion gap-decay algebra accepts, independently for the same four
rows, caller values

``hF = exp(-wF)``, ``d = 1 - exp(-(wP-wF))``,
``etaF = wF_Z`` and ``deta = etaP - etaF``.

Without computing any of those quantities, it returns

``delta_g   = hF * d``
``delta_g_Z = hF * ((1 - d) * deta - d * etaF)``.

Both APIs preserve ``Fraction`` and caller-owned interval/symbolic arithmetic
objects, reject missing values, built-in integers and floats, and reject
obvious midpoint-selected wrappers.  All source ownership, F identity,
interval certification, quadrature/error certification, and physical-scale
claims remain false.  Direct-difference identities and the supplied
``expm1``/gap-decay ports are caller premises only.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from fractions import Fraction
from typing import Any


SCHEMA = "gap_closed_form_signed_difference_v1"
GAP_DECAY_SCHEMA = "gap_closed_form_gap_decay_difference_v1"

ROW_NAMES = ("linear", "quadratic", "endpoint_upper", "endpoint_lower")
ROW_LABELS = ("linear+", "quadratic+", "upper-", "lower+")
ROW_SIGNS = ("+", "+", "-", "+")
ROW_CONVENTION = (
    "ordinary Z1 derivative rows only; Pi1, pressure S1, P1 and m15 roles "
    "are not identified here"
)
P1_IS_M15 = False
P1_NOT_M15 = True

F_OWNERSHIP = False
SOURCE_IDENTITY_PROVED = False
INTERVAL_CERTIFIED = False
PHYSICAL_SCALE_PROVED = False
QUADRATURE_ERROR_CERTIFIED = False
ERROR_CERTIFIED = False

_MISSING = object()
_ONE = Fraction(1, 1)


class GapClosedFormSignedDifferenceInputError(TypeError):
    """Raised when a caller-owned scalar or row is missing or unsafe."""


# The shorter name is convenient for callers that only need the boundary
# contract.  It is an alias, not a second implementation.
SignedDifferenceInputError = GapClosedFormSignedDifferenceInputError


def _is_midpoint(value: Any) -> bool:
    """Reject common explicit midpoint wrappers without inspecting them."""

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


def _require_scalar(value: Any, name: str) -> Any:
    """Keep caller arithmetic unchanged while rejecting point shortcuts."""

    if value is _MISSING:
        raise GapClosedFormSignedDifferenceInputError(
            f"missing required caller input: {name}"
        )
    if value is None:
        raise GapClosedFormSignedDifferenceInputError(
            f"caller input {name} cannot be None"
        )
    if isinstance(value, bool):
        raise GapClosedFormSignedDifferenceInputError(
            f"caller input {name} cannot be bool"
        )
    if isinstance(value, (int, float, complex)):
        raise GapClosedFormSignedDifferenceInputError(
            f"caller input {name} must be Fraction or a caller interval; "
            "plain int, float and complex values are rejected"
        )
    if isinstance(
        value,
        (str, bytes, bytearray, list, tuple, set, frozenset, dict),
    ):
        raise GapClosedFormSignedDifferenceInputError(
            f"caller input {name} must be one scalar, not a container or text value"
        )
    if _is_midpoint(value):
        raise GapClosedFormSignedDifferenceInputError(
            f"caller input {name} must not be a midpoint-selected value"
        )
    # Fraction, directed interval and symbolic scalar objects are retained.
    # No coercion or point extraction is performed.
    return value


def _pick(value: Any, names: tuple[str, ...], label: str) -> Any:
    """Read one named field from a mapping/object without coercion."""

    if isinstance(value, Mapping):
        for name in names:
            if name in value:
                return value[name]
    else:
        for name in names:
            try:
                return getattr(value, name)
            except AttributeError:
                continue
    raise GapClosedFormSignedDifferenceInputError(
        f"{label} is missing; expected one of {', '.join(names)}"
    )


@dataclass(frozen=True, slots=True)
class GapClosedFormSignedDifferenceInput:
    """Six caller-owned ports for one closed-form signed row."""

    LF: Any
    aF: Any
    delta_a: Any
    aF_physZ: Any
    delta_physZ: Any
    t: Any

    def __post_init__(self) -> None:
        for name in ("LF", "aF", "delta_a", "aF_physZ", "delta_physZ", "t"):
            _require_scalar(getattr(self, name), f"row.{name}")

    @property
    def frame_logscale(self) -> Any:
        """Alias for the untouched frame log-scale token."""

        return self.LF

    @property
    def frame_mantissa(self) -> Any:
        return self.aF

    @property
    def direct_delta_a(self) -> Any:
        return self.delta_a

    @property
    def frame_physZ(self) -> Any:
        return self.aF_physZ

    @property
    def direct_delta_physZ(self) -> Any:
        return self.delta_physZ

    @property
    def expm1_gap(self) -> Any:
        """The caller's supplied ``expm1(LP-LF)`` port."""

        return self.t


# Compact aliases make the six-port contract discoverable without changing it.
SignedGapRow = GapClosedFormSignedDifferenceInput
GapSignedDifferenceRow = GapClosedFormSignedDifferenceInput


def _coerce_signed_row(value: Any, label: str) -> GapClosedFormSignedDifferenceInput:
    if isinstance(value, GapClosedFormSignedDifferenceInput):
        return value
    if isinstance(value, Mapping) or not isinstance(value, (str, bytes, bytearray)):
        if isinstance(value, Sequence) and not isinstance(value, Mapping):
            raw = tuple(value)
            if len(raw) != 6:
                raise GapClosedFormSignedDifferenceInputError(
                    f"{label} must provide six ports: LF, aF, delta_a, "
                    "aF_physZ, delta_physZ and t"
                )
            return GapClosedFormSignedDifferenceInput(*raw)
        return GapClosedFormSignedDifferenceInput(
            LF=_pick(value, ("LF", "LF_j", "frame_logscale", "log_scale"), label + ".LF"),
            aF=_pick(value, ("aF", "frame_mantissa", "aF_mantissa"), label + ".aF"),
            delta_a=_pick(
                value,
                ("delta_a", "direct_delta_a", "aP_minus_aF"),
                label + ".delta_a",
            ),
            aF_physZ=_pick(
                value,
                ("aF_physZ", "aF_phys_Z", "frame_physZ", "frame_phys_z"),
                label + ".aF_physZ",
            ),
            delta_physZ=_pick(
                value,
                ("delta_physZ", "delta_phys_Z", "direct_delta_physZ"),
                label + ".delta_physZ",
            ),
            t=_pick(value, ("t", "expm1_gap", "expm1_delta_logscale"), label + ".t"),
        )
    raise GapClosedFormSignedDifferenceInputError(
        f"{label} must be a six-port row, mapping or row object"
    )


def _resolve_rows(
    rows: Mapping[str, Any] | Sequence[Any] | None,
    *,
    linear: Any,
    quadratic: Any,
    endpoint_upper: Any,
    endpoint_lower: Any,
    label: str,
) -> tuple[GapClosedFormSignedDifferenceInput, ...]:
    supplied_keywords = (linear, quadratic, endpoint_upper, endpoint_lower)
    if rows is not None and any(value is not _MISSING for value in supplied_keywords):
        raise GapClosedFormSignedDifferenceInputError(
            f"{label}: pass rows or named row arguments, not both"
        )
    if rows is None:
        raw_rows = supplied_keywords
    elif isinstance(rows, Mapping):
        extras = set(rows) - set(ROW_NAMES)
        missing = [name for name in ROW_NAMES if name not in rows]
        if missing:
            raise GapClosedFormSignedDifferenceInputError(
                f"{label} is missing rows: {', '.join(missing)}"
            )
        if extras:
            raise GapClosedFormSignedDifferenceInputError(
                f"{label} has unexpected rows: {', '.join(sorted(map(str, extras)))}"
            )
        raw_rows = tuple(rows[name] for name in ROW_NAMES)
    else:
        if isinstance(rows, (str, bytes, bytearray)):
            raise GapClosedFormSignedDifferenceInputError(
                f"{label} must provide exactly four named rows"
            )
        try:
            raw_rows = tuple(rows)
        except (TypeError, ValueError) as error:
            raise GapClosedFormSignedDifferenceInputError(
                f"{label} must provide exactly four named rows"
            ) from error
        if len(raw_rows) != 4:
            raise GapClosedFormSignedDifferenceInputError(
                f"{label} requires exactly four rows: {', '.join(ROW_NAMES)}"
            )
    if any(value is _MISSING for value in raw_rows):
        missing = [
            name for name, value in zip(ROW_NAMES, raw_rows) if value is _MISSING
        ]
        raise GapClosedFormSignedDifferenceInputError(
            f"{label} is missing rows: {', '.join(missing)}"
        )
    return tuple(
        _coerce_signed_row(value, f"{label}.{name}")
        for name, value in zip(ROW_NAMES, raw_rows)
    )


def _add(left: Any, right: Any) -> Any:
    try:
        return left + right
    except TypeError:
        return right + left


def _mul(left: Any, right: Any) -> Any:
    try:
        return left * right
    except TypeError:
        return right * left


def _sub(left: Any, right: Any) -> Any:
    try:
        return left - right
    except TypeError:
        return -(right - left)


def _neg(value: Any) -> Any:
    return -value


def _closed_form_values(
    row: GapClosedFormSignedDifferenceInput, label: str
) -> tuple[Any, Any]:
    """Apply the two supplied-port identities without constructing P or F."""

    try:
        one_plus_t = _add(_ONE, row.t)
        value = _add(
            _mul(one_plus_t, row.delta_a),
            _mul(row.t, row.aF),
        )
        physical_z = _add(
            _mul(one_plus_t, row.delta_physZ),
            _mul(row.t, row.aF_physZ),
        )
    except (ArithmeticError, OverflowError, TypeError, ValueError) as error:
        raise GapClosedFormSignedDifferenceInputError(
            f"caller-owned arithmetic failed for {label}: {error}"
        ) from error
    return value, physical_z


@dataclass(frozen=True, slots=True)
class SignedGapComponent:
    """One unmerged signed row and its m15 orientation."""

    name: str
    sign: str
    LF: Any
    signed_value_mantissa: Any
    signed_physicalZ_mantissa: Any
    contribution_mantissa: Any
    contribution_physicalZ_mantissa: Any
    source_ownership: bool = False
    source_identity_proved: bool = False
    interval_certified: bool = False
    quadrature_error_certified: bool = False
    error_certified: bool = False

    @property
    def label(self) -> str:
        return self.name

    @property
    def frame_logscale(self) -> Any:
        return self.LF

    @property
    def value_mantissa(self) -> Any:
        """Raw direct signed formula before the m15 orientation sign."""

        return self.signed_value_mantissa

    @property
    def physicalZ_mantissa(self) -> Any:
        return self.signed_physicalZ_mantissa

    @property
    def signed_mantissa(self) -> Any:
        return self.signed_value_mantissa

    @property
    def signed_physZ_mantissa(self) -> Any:
        return self.signed_physicalZ_mantissa

    @property
    def aPhysZ_mantissa(self) -> Any:
        """Alias for the returned direct signed physical ``Z1`` mantissa."""

        return self.signed_physicalZ_mantissa

    @property
    def F_ownership(self) -> bool:
        return False

    @property
    def physical_scale_proved(self) -> bool:
        return False

    @property
    def certified(self) -> bool:
        return False

    @property
    def source_evidence_false(self) -> bool:
        return True

    @property
    def interval_evidence_false(self) -> bool:
        return True

    @property
    def m15_contribution(self) -> Any:
        return self.contribution_mantissa

    @property
    def m15_contribution_physZ(self) -> Any:
        return self.contribution_physicalZ_mantissa


@dataclass(frozen=True, slots=True)
class GapClosedFormSignedDifferenceResult:
    """Four separate scaled components of the signed m15 decomposition."""

    linear_plus: SignedGapComponent
    quadratic_plus: SignedGapComponent
    upper_minus: SignedGapComponent
    lower_plus: SignedGapComponent
    schema: str = SCHEMA
    row_convention: str = ROW_CONVENTION
    ordinary_Z1: bool = True
    P1_is_m15: bool = P1_IS_M15
    P1_not_m15: bool = P1_NOT_M15
    scales_materialized: bool = False
    frame_logscale_materialized: bool = False
    t_expm1_caller_premise: bool = True
    exponential_evaluated: bool = False
    direct_difference_identity_caller_premise: bool = True
    direct_difference_identity_proved: bool = False
    source_ownership: bool = False
    F_ownership: bool = False
    source_identity_proved: bool = False
    same_source_proved: bool = False
    physical_scale_proved: bool = False
    interval_certified: bool = False
    certified: bool = False
    quadrature_error_certified: bool = False
    error_certified: bool = False
    endpoint_terms_retained: bool = True
    terms_not_combined: bool = True
    LF_Z_added: bool = False
    t_Z_added: bool = False

    @property
    def linear(self) -> SignedGapComponent:
        return self.linear_plus

    @property
    def quadratic(self) -> SignedGapComponent:
        return self.quadratic_plus

    @property
    def endpoint_upper(self) -> SignedGapComponent:
        return self.upper_minus

    @property
    def endpoint_lower(self) -> SignedGapComponent:
        return self.lower_plus

    @property
    def terms(self) -> tuple[SignedGapComponent, ...]:
        return (
            self.linear_plus,
            self.quadratic_plus,
            self.upper_minus,
            self.lower_plus,
        )

    def __iter__(self):
        return iter(self.terms)

    def as_dict(self) -> dict[str, Any]:
        """Return metadata and caller objects without point conversion."""

        return {
            "schema": self.schema,
            "row_convention": self.row_convention,
            "ordinary_Z1": self.ordinary_Z1,
            "P1_is_m15": self.P1_is_m15,
            "P1_not_m15": self.P1_not_m15,
            "scales_materialized": self.scales_materialized,
            "frame_logscale_materialized": self.frame_logscale_materialized,
            "t_expm1_caller_premise": self.t_expm1_caller_premise,
            "exponential_evaluated": self.exponential_evaluated,
            "direct_difference_identity_caller_premise": self.direct_difference_identity_caller_premise,
            "direct_difference_identity_proved": self.direct_difference_identity_proved,
            "source_ownership": self.source_ownership,
            "F_ownership": self.F_ownership,
            "source_identity_proved": self.source_identity_proved,
            "same_source_proved": self.same_source_proved,
            "physical_scale_proved": self.physical_scale_proved,
            "interval_certified": self.interval_certified,
            "certified": self.certified,
            "quadrature_error_certified": self.quadrature_error_certified,
            "error_certified": self.error_certified,
            "endpoint_terms_retained": self.endpoint_terms_retained,
            "terms_not_combined": self.terms_not_combined,
            "LF_Z_added": self.LF_Z_added,
            "t_Z_added": self.t_Z_added,
            "terms": {
                term.name: {
                    "sign": term.sign,
                    "LF": term.LF,
                    "signed_value_mantissa": term.signed_value_mantissa,
                    "signed_physicalZ_mantissa": term.signed_physicalZ_mantissa,
                    "contribution_mantissa": term.contribution_mantissa,
                    "contribution_physicalZ_mantissa": term.contribution_physicalZ_mantissa,
                    "source_ownership": term.source_ownership,
                    "source_identity_proved": term.source_identity_proved,
                    "interval_certified": term.interval_certified,
                    "quadrature_error_certified": term.quadrature_error_certified,
                    "error_certified": term.error_certified,
                }
                for term in self.terms
            },
        }


def build_gap_closed_form_signed_difference(
    rows: Mapping[str, Any] | Sequence[Any] | None = None,
    *,
    linear: Any = _MISSING,
    quadratic: Any = _MISSING,
    endpoint_upper: Any = _MISSING,
    endpoint_lower: Any = _MISSING,
) -> GapClosedFormSignedDifferenceResult:
    """Build four independent closed-form signed ``P-F`` components.

    ``rows`` may be a mapping with exactly the four names in :data:`ROW_NAMES`,
    a sequence in that order, or omitted in favour of four named arguments.
    Each row is a :class:`GapClosedFormSignedDifferenceInput`, a six-item
    sequence ``(LF, aF, delta_a, aF_physZ, delta_physZ, t)``, or an object/
    mapping exposing those fields.  No result component is summed with another.
    """

    resolved = _resolve_rows(
        rows,
        linear=linear,
        quadratic=quadratic,
        endpoint_upper=endpoint_upper,
        endpoint_lower=endpoint_lower,
        label="signed rows",
    )
    components: list[SignedGapComponent] = []
    for row_name, row_label, sign, row in zip(
        ROW_NAMES, ROW_LABELS, ROW_SIGNS, resolved
    ):
        value, physical_z = _closed_form_values(row, row_name)
        if sign == "-":
            contribution = _neg(value)
            contribution_physical_z = _neg(physical_z)
        else:
            contribution = value
            contribution_physical_z = physical_z
        components.append(
            SignedGapComponent(
                name=row_label,
                sign=sign,
                LF=row.LF,
                signed_value_mantissa=value,
                signed_physicalZ_mantissa=physical_z,
                contribution_mantissa=contribution,
                contribution_physicalZ_mantissa=contribution_physical_z,
            )
        )
    return GapClosedFormSignedDifferenceResult(*components)


# Descriptive aliases do not introduce another formula or import a legacy
# implementation.
evaluate_gap_closed_form_signed_difference = build_gap_closed_form_signed_difference
gap_closed_form_signed_difference = build_gap_closed_form_signed_difference
compute_gap_closed_form_signed_difference = build_gap_closed_form_signed_difference


@dataclass(frozen=True, slots=True)
class CorrelatedGapDecayDifferenceInput:
    """Four caller-owned ports for one correlated gap-decay difference row."""

    hF: Any
    d: Any
    etaF: Any
    deta: Any

    def __post_init__(self) -> None:
        for name in ("hF", "d", "etaF", "deta"):
            _require_scalar(getattr(self, name), f"gap row.{name}")

    @property
    def h_F(self) -> Any:
        return self.hF

    @property
    def eta_F(self) -> Any:
        return self.etaF

    @property
    def delta_eta(self) -> Any:
        return self.deta


CorrelatedGapRow = CorrelatedGapDecayDifferenceInput
GapDecayDifferenceRow = CorrelatedGapDecayDifferenceInput


def _coerce_gap_row(value: Any, label: str) -> CorrelatedGapDecayDifferenceInput:
    if isinstance(value, CorrelatedGapDecayDifferenceInput):
        return value
    if isinstance(value, Mapping) or not isinstance(value, (str, bytes, bytearray)):
        if isinstance(value, Sequence) and not isinstance(value, Mapping):
            raw = tuple(value)
            if len(raw) != 4:
                raise GapClosedFormSignedDifferenceInputError(
                    f"{label} must provide hF, d, etaF and deta"
                )
            return CorrelatedGapDecayDifferenceInput(*raw)
        return CorrelatedGapDecayDifferenceInput(
            hF=_pick(value, ("hF", "h_F", "frame_decay"), label + ".hF"),
            d=_pick(value, ("d", "delta_g_factor", "gap_delta"), label + ".d"),
            etaF=_pick(value, ("etaF", "eta_F", "frame_eta"), label + ".etaF"),
            deta=_pick(
                value,
                ("deta", "delta_eta", "direct_delta_eta"),
                label + ".deta",
            ),
        )
    raise GapClosedFormSignedDifferenceInputError(
        f"{label} must be a four-port gap-decay row, mapping or row object"
    )


def _correlated_gap_values(
    row: CorrelatedGapDecayDifferenceInput, label: str
) -> tuple[Any, Any]:
    """Apply the caller-supplied ``hF,d,etaF,deta`` gap identities."""

    try:
        one_minus_d = _sub(_ONE, row.d)
        delta_g = _mul(row.hF, row.d)
        delta_g_z = _mul(
            row.hF,
            _sub(
                _mul(one_minus_d, row.deta),
                _mul(row.d, row.etaF),
            ),
        )
    except (ArithmeticError, OverflowError, TypeError, ValueError) as error:
        raise GapClosedFormSignedDifferenceInputError(
            f"caller-owned arithmetic failed for {label}: {error}"
        ) from error
    return delta_g, delta_g_z


@dataclass(frozen=True, slots=True)
class CorrelatedGapDecayComponent:
    """One unmerged correlated gap-decay difference and ordinary Z1 row."""

    name: str
    sign: str
    signed_gap_mantissa: Any
    signed_gap_physZ_mantissa: Any
    contribution_mantissa: Any
    contribution_physZ_mantissa: Any
    source_ownership: bool = False
    source_identity_proved: bool = False
    interval_certified: bool = False
    quadrature_error_certified: bool = False
    error_certified: bool = False

    @property
    def delta_g(self) -> Any:
        return self.signed_gap_mantissa

    @property
    def delta_g_Z(self) -> Any:
        return self.signed_gap_physZ_mantissa

    @property
    def signed_physZ_mantissa(self) -> Any:
        return self.signed_gap_physZ_mantissa

    @property
    def delta_g_physZ(self) -> Any:
        return self.signed_gap_physZ_mantissa

    @property
    def F_ownership(self) -> bool:
        return False

    @property
    def physical_scale_proved(self) -> bool:
        return False

    @property
    def certified(self) -> bool:
        return False

    @property
    def source_evidence_false(self) -> bool:
        return True

    @property
    def interval_evidence_false(self) -> bool:
        return True

    @property
    def m15_contribution(self) -> Any:
        return self.contribution_mantissa

    @property
    def m15_contribution_physZ(self) -> Any:
        return self.contribution_physZ_mantissa


@dataclass(frozen=True, slots=True)
class CorrelatedGapDecayDifferenceResult:
    """Four separate correlated gap-decay rows with no cross-scale sum."""

    linear_plus: CorrelatedGapDecayComponent
    quadratic_plus: CorrelatedGapDecayComponent
    upper_minus: CorrelatedGapDecayComponent
    lower_plus: CorrelatedGapDecayComponent
    schema: str = GAP_DECAY_SCHEMA
    row_convention: str = ROW_CONVENTION
    ordinary_Z1: bool = True
    P1_is_m15: bool = P1_IS_M15
    P1_not_m15: bool = P1_NOT_M15
    direct_gap_ports_caller_premise: bool = True
    direct_gap_ports_proved: bool = False
    exponential_evaluated: bool = False
    source_ownership: bool = False
    F_ownership: bool = False
    source_identity_proved: bool = False
    same_source_proved: bool = False
    physical_scale_proved: bool = False
    interval_certified: bool = False
    certified: bool = False
    quadrature_error_certified: bool = False
    error_certified: bool = False
    endpoint_terms_retained: bool = True
    terms_not_combined: bool = True

    @property
    def linear(self) -> CorrelatedGapDecayComponent:
        return self.linear_plus

    @property
    def quadratic(self) -> CorrelatedGapDecayComponent:
        return self.quadratic_plus

    @property
    def endpoint_upper(self) -> CorrelatedGapDecayComponent:
        return self.upper_minus

    @property
    def endpoint_lower(self) -> CorrelatedGapDecayComponent:
        return self.lower_plus

    @property
    def terms(self) -> tuple[CorrelatedGapDecayComponent, ...]:
        return (
            self.linear_plus,
            self.quadratic_plus,
            self.upper_minus,
            self.lower_plus,
        )

    def __iter__(self):
        return iter(self.terms)

    def as_dict(self) -> dict[str, Any]:
        return {
            "schema": self.schema,
            "row_convention": self.row_convention,
            "ordinary_Z1": self.ordinary_Z1,
            "P1_is_m15": self.P1_is_m15,
            "P1_not_m15": self.P1_not_m15,
            "direct_gap_ports_caller_premise": self.direct_gap_ports_caller_premise,
            "direct_gap_ports_proved": self.direct_gap_ports_proved,
            "exponential_evaluated": self.exponential_evaluated,
            "source_ownership": self.source_ownership,
            "F_ownership": self.F_ownership,
            "source_identity_proved": self.source_identity_proved,
            "same_source_proved": self.same_source_proved,
            "physical_scale_proved": self.physical_scale_proved,
            "interval_certified": self.interval_certified,
            "certified": self.certified,
            "quadrature_error_certified": self.quadrature_error_certified,
            "error_certified": self.error_certified,
            "endpoint_terms_retained": self.endpoint_terms_retained,
            "terms_not_combined": self.terms_not_combined,
            "terms": {
                term.name: {
                    "sign": term.sign,
                    "signed_gap_mantissa": term.signed_gap_mantissa,
                    "signed_gap_physZ_mantissa": term.signed_gap_physZ_mantissa,
                    "contribution_mantissa": term.contribution_mantissa,
                    "contribution_physZ_mantissa": term.contribution_physZ_mantissa,
                    "source_ownership": term.source_ownership,
                    "source_identity_proved": term.source_identity_proved,
                    "interval_certified": term.interval_certified,
                    "quadrature_error_certified": term.quadrature_error_certified,
                    "error_certified": term.error_certified,
                }
                for term in self.terms
            },
        }


def build_correlated_gap_decay_difference(
    rows: Mapping[str, Any] | Sequence[Any] | None = None,
    *,
    linear: Any = _MISSING,
    quadratic: Any = _MISSING,
    endpoint_upper: Any = _MISSING,
    endpoint_lower: Any = _MISSING,
) -> CorrelatedGapDecayDifferenceResult:
    """Build four independent ``Delta g``/ordinary-``Z1`` rows.

    The four rows are required and are never summed.  Each row is a
    :class:`CorrelatedGapDecayDifferenceInput`, a four-item sequence
    ``(hF, d, etaF, deta)``, or a mapping/object exposing those fields.
    """

    resolved_raw: tuple[Any, ...]
    supplied_keywords = (linear, quadratic, endpoint_upper, endpoint_lower)
    if rows is not None and any(value is not _MISSING for value in supplied_keywords):
        raise GapClosedFormSignedDifferenceInputError(
            "gap-decay rows: pass rows or named row arguments, not both"
        )
    if rows is None:
        resolved_raw = supplied_keywords
    elif isinstance(rows, Mapping):
        extras = set(rows) - set(ROW_NAMES)
        missing = [name for name in ROW_NAMES if name not in rows]
        if missing:
            raise GapClosedFormSignedDifferenceInputError(
                "gap-decay rows are missing: " + ", ".join(missing)
            )
        if extras:
            raise GapClosedFormSignedDifferenceInputError(
                "gap-decay rows have unexpected names: "
                + ", ".join(sorted(map(str, extras)))
            )
        resolved_raw = tuple(rows[name] for name in ROW_NAMES)
    else:
        if isinstance(rows, (str, bytes, bytearray)):
            raise GapClosedFormSignedDifferenceInputError(
                "gap-decay rows must provide exactly four named rows"
            )
        try:
            resolved_raw = tuple(rows)
        except (TypeError, ValueError) as error:
            raise GapClosedFormSignedDifferenceInputError(
                "gap-decay rows must provide exactly four named rows"
            ) from error
        if len(resolved_raw) != 4:
            raise GapClosedFormSignedDifferenceInputError(
                "gap-decay rows require exactly four rows: " + ", ".join(ROW_NAMES)
            )
    if any(value is _MISSING for value in resolved_raw):
        missing = [
            name
            for name, value in zip(ROW_NAMES, resolved_raw)
            if value is _MISSING
        ]
        raise GapClosedFormSignedDifferenceInputError(
            "gap-decay rows are missing: " + ", ".join(missing)
        )

    components: list[CorrelatedGapDecayComponent] = []
    for row_name, row_label, sign, raw in zip(
        ROW_NAMES, ROW_LABELS, ROW_SIGNS, resolved_raw
    ):
        row = _coerce_gap_row(raw, "gap-decay rows." + row_name)
        delta_g, delta_g_z = _correlated_gap_values(row, row_name)
        if sign == "-":
            contribution = _neg(delta_g)
            contribution_z = _neg(delta_g_z)
        else:
            contribution = delta_g
            contribution_z = delta_g_z
        components.append(
            CorrelatedGapDecayComponent(
                name=row_label,
                sign=sign,
                signed_gap_mantissa=delta_g,
                signed_gap_physZ_mantissa=delta_g_z,
                contribution_mantissa=contribution,
                contribution_physZ_mantissa=contribution_z,
            )
        )
    return CorrelatedGapDecayDifferenceResult(*components)


evaluate_correlated_gap_decay_difference = build_correlated_gap_decay_difference
correlated_gap_decay_difference = build_correlated_gap_decay_difference
compute_correlated_gap_decay_difference = build_correlated_gap_decay_difference


def make_gap_closed_form_signed_difference_input(
    LF: Any,
    aF: Any,
    delta_a: Any,
    aF_physZ: Any,
    delta_physZ: Any,
    t: Any,
) -> GapClosedFormSignedDifferenceInput:
    """Construct one six-port row without changing caller arithmetic."""

    return GapClosedFormSignedDifferenceInput(
        LF=LF,
        aF=aF,
        delta_a=delta_a,
        aF_physZ=aF_physZ,
        delta_physZ=delta_physZ,
        t=t,
    )


make_signed_gap_row = make_gap_closed_form_signed_difference_input


def make_correlated_gap_decay_difference_input(
    hF: Any,
    d: Any,
    etaF: Any,
    deta: Any,
) -> CorrelatedGapDecayDifferenceInput:
    """Construct one four-port gap-decay row without evaluating exponentials."""

    return CorrelatedGapDecayDifferenceInput(hF=hF, d=d, etaF=etaF, deta=deta)


make_correlated_gap_row = make_correlated_gap_decay_difference_input


__all__ = [
    "ERROR_CERTIFIED",
    "F_OWNERSHIP",
    "GAP_DECAY_SCHEMA",
    "INTERVAL_CERTIFIED",
    "PHYSICAL_SCALE_PROVED",
    "QUADRATURE_ERROR_CERTIFIED",
    "ROW_CONVENTION",
    "ROW_LABELS",
    "ROW_NAMES",
    "ROW_SIGNS",
    "P1_IS_M15",
    "P1_NOT_M15",
    "SCHEMA",
    "SOURCE_IDENTITY_PROVED",
    "CorrelatedGapDecayComponent",
    "CorrelatedGapDecayDifferenceInput",
    "CorrelatedGapDecayDifferenceResult",
    "CorrelatedGapRow",
    "GapClosedFormSignedDifferenceInput",
    "GapClosedFormSignedDifferenceInputError",
    "GapClosedFormSignedDifferenceResult",
    "GapDecayDifferenceRow",
    "GapSignedDifferenceRow",
    "SignedDifferenceInputError",
    "SignedGapComponent",
    "SignedGapRow",
    "build_correlated_gap_decay_difference",
    "build_gap_closed_form_signed_difference",
    "compute_correlated_gap_decay_difference",
    "compute_gap_closed_form_signed_difference",
    "correlated_gap_decay_difference",
    "evaluate_correlated_gap_decay_difference",
    "evaluate_gap_closed_form_signed_difference",
    "gap_closed_form_signed_difference",
    "make_correlated_gap_decay_difference_input",
    "make_correlated_gap_row",
    "make_gap_closed_form_signed_difference_input",
    "make_signed_gap_row",
]
